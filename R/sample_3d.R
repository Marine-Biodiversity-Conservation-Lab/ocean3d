# Sampling SpatVector3D lines into points (#105): the build step SPEC.md
# plans for tracks ("Extraction on the three steps", 1. Build). Points follow
# WGS84 geodesics between vertices, with depth and time interpolated along
# each segment; extract_to_track() samples its tracks with sample_3d().

#' Sample 3D lines into points
#'
#' Turns each line of a [SpatVector3D-class], such as a dive track or a
#' vessel route, into points along it: its vertices, and with a spacing or
#' step, points in between. Each point has the depth and time interpolated
#' along the segment it lies on. The points can then go into the
#' intersect step ([intersect_3d()], [mask()]) or be buffered with
#' [buffer_3d()].
#'
#' Points between vertices follow the WGS84 geodesic (the shortest path on
#' the ellipsoid), so a segment across the antimeridian goes the short way.
#' Depth and time change linearly with distance along each segment.
#' Distances are horizontal, not three-dimensional.
#'
#' A line with a time at some vertices but not others is an error: its time
#' between those vertices isn't known. Give every vertex a time, or none.
#'
#' @param x A [SpatVector3D-class] of lines, from [as_vector3d()].
#' @param sampling How to place the points:
#'   * `"vertices"` keeps only the vertices.
#'   * `"balanced"` spaces points evenly along each segment, as many as the
#'     strictest step limit needs. With calendar steps, it keeps the step
#'     dates and spaces points evenly between them.
#'   * `"merged_steps"` places the evenly spaced points that the distance
#'     and depth limits need, and separately a point at each time step, and
#'     keeps both sets, so the points can be unevenly spaced.
#'
#'   `NULL` (default) gives `"balanced"` when a step is given and
#'   `"vertices"` otherwise.
#' @param spacing_km Largest horizontal distance between points, in
#'   kilometres.
#' @param depth_step_m Largest change in depth between points, in metres.
#' @param time_step,time_unit Largest change in time between points, and its
#'   unit: seconds, minutes, hours, days, or whole calendar months or years.
#'   Needs a time at every vertex. If the lines have dates without a time of
#'   day, the step must be whole days, months or years. Calendar steps count
#'   from the first vertex of each segment, so a segment that starts on
#'   31 January with monthly steps gets 29 February, 31 March, 30 April.
#' @param max_samples Largest total number of points, across all lines.
#'   `sample_3d()` stops with an error if the steps need more.
#'
#' @return A [SpatVector3D-class] of points, one per sample, in order along
#'   each line. Each point's feature ID is the line's ID and the sample
#'   number (`"A_3"`). Columns: `track_id` (the line's ID), `sample_id`
#'   (1, 2, ... along the line), `distance_km` (distance from the line's
#'   start), `support_km` and `support_seconds` (half the length and
#'   duration of the segments either side, as weights for summaries along
#'   the track; `NA` without times), `is_interpolated` (`FALSE` for the
#'   line's own vertices), and the line's columns that have one value along
#'   it. If the lines have dates without a time of day, each point's time is
#'   the UTC date it falls on.
#' @seealso [buffer_3d()] to turn the points into polygons;
#'   [extract_to_track()], which samples tracks with `sample_3d()` and
#'   summarises the values along them.
#' @examples
#' dive <- data.frame(tag = "A", lon = c(150, 150.1, 150.2), lat = -20,
#'                    depth = c(0, 120, 15),
#'                    time = as.POSIXct(c("2020-01-01 00:00",
#'                                        "2020-01-01 02:00",
#'                                        "2020-01-01 03:00"), tz = "UTC"))
#' track <- as_vector3d(dive, type = "lines", id = "tag", time = "time")
#'
#' # A point at least every 2 km, with depth and time along each segment.
#' samples <- sample_3d(track, spacing_km = 2)
#' as.data.frame(samples)[c("track_id", "sample_id", "distance_km", "depth",
#'                          "time")]
#'
#' # Buffer the samples to follow the dive more closely than one buffer
#' # around the whole line.
#' buffer_3d(samples, radius_km = 1, depth_window = 10)
#' @export
sample_3d <- function(x, sampling = NULL, spacing_km = NULL,
                      depth_step_m = NULL, time_step = NULL,
                      time_unit = c("hours", "minutes", "seconds", "days",
                                    "months", "years"),
                      max_samples = 100000L) {
  if (!methods::is(x, "SpatVector3D")) {
    .abort_vector3d(
      "`x` must be SpatVector3D lines, from as_vector3d(type = \"lines\")."
    )
  }
  .check_geomtype(
    x, "lines",
    "Points are already samples, and polygons have no path to follow."
  )
  if (!is.null(sampling)) {
    sampling <- match.arg(sampling, c("vertices", "balanced", "merged_steps"))
  }
  time_unit <- match.arg(time_unit)
  sampling <- .track_sampling(sampling, spacing_km, depth_step_m, time_step,
                              time_unit, max_samples)
  # Vertices in feature order, then in order along each line, as
  # buffer_3d() reads them: the sampler takes each line's rows in order.
  v <- x@vertices
  ids <- unique(v$feature_id)
  v <- v[order(match(v$feature_id, ids), v$vertex), , drop = FALSE]
  .check_vertex_times(v, if (!is.null(time_step)) "`time_step`")
  if (!is.null(time_step) && !time_unit %in% c("months", "years")) {
    .check_whole_days(x, time_step * .time_unit_seconds(time_unit),
                      "`time_step`")
  }
  timed <- !all(is.na(v$time))

  # The points' own columns can't also come from the lines.
  extra <- .feature_attributes(v)
  generated <- c("track_id", "sample_id", "distance_km", "support_km",
                 "support_seconds", "is_interpolated")
  clash <- intersect(names(extra), generated)
  if (length(clash)) {
    .abort_vector3d(glue::glue(
      "The lines have columns {toString(clash)}, which sample_3d() adds to ",
      "each point. Rename them first."
    ))
  }

  options <- list(sampling = sampling, spacing_km = spacing_km,
                  depth_step_m = depth_step_m, time_step = time_step,
                  time_unit = time_unit)
  # The lines share one budget of `max_samples` points.
  pieces <- list()
  remaining <- max_samples
  for (id in ids) {
    on_line <- v[v$feature_id == id, , drop = FALSE]
    piece <- .track_sample_one(
      data.frame(lon = on_line$lon, lat = on_line$lat, depth = on_line$depth,
                 time = as.numeric(on_line$time)),
      options, remaining
    )
    remaining <- remaining - nrow(piece)
    pieces[[length(pieces) + 1L]] <- cbind(track_id = id, piece,
                                           stringsAsFactors = FALSE)
  }
  s <- do.call(rbind, pieces)
  s$time <- .POSIXct(s$time, tz = "UTC")

  time <- s$time
  if (identical(x@time_resolution, "daily")) {
    time <- .time_comparison_values(time, "daily")
  }
  points <- data.frame(
    feature_id = paste(s$track_id, s$sample_id, sep = "_"),
    vertex = 1L,
    lon = .normalise_longitude(s$lon),
    lat = s$lat,
    depth = s$depth,
    time = if (timed) time else .POSIXct(rep(NA_real_, nrow(s)), tz = "UTC"),
    track_id = s$track_id,
    sample_id = s$sample_id,
    distance_km = s$distance_km,
    support_km = s$support_km,
    support_seconds = if (timed) s$support_seconds else NA_real_,
    is_interpolated = s$is_interpolated,
    stringsAsFactors = FALSE
  )
  extra <- extra[match(points$track_id, ids), , drop = FALSE]
  .new_vector3d("points", x@time_resolution,
                vertices = .add_attributes(points, extra))
}


# Seconds in one elapsed (non-calendar) time unit.
.time_unit_seconds <- function(unit) {
  switch(unit, days = 86400, hours = 3600, minutes = 60, seconds = 1)
}


# The one producer of `ocean3d_error_track_sampling`: step limits, a sample
# budget or a track sample_3d() can't sample. Raised by .track_sampling()
# and the sampler.
.abort_track_sampling <- function(message) {
  rlang::abort(message, class = "ocean3d_error_track_sampling", call = NULL)
}


# The sampling mode, checked with its step limits and `max_samples`, for
# sample_3d(). NULL picks "balanced" when a step is given and "vertices"
# otherwise. Steps must be positive, calendar steps whole months or years,
# and `max_samples` a positive whole number; vertices sampling takes no
# steps, and the stepped modes need one.
.track_sampling <- function(sampling, spacing_km, depth_step_m, time_step,
                            time_unit, max_samples) {
  positive <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0
  }
  for (name in c("spacing_km", "depth_step_m", "time_step")) {
    value <- get(name)
    if (!is.null(value) && !positive(value)) {
      .abort_track_sampling(glue::glue(
        "`{name}` must be one finite positive number."
      ))
    }
  }
  if (!positive(max_samples) || max_samples != floor(max_samples)) {
    .abort_track_sampling("`max_samples` must be a positive integer.")
  }
  has_steps <- !is.null(spacing_km) || !is.null(depth_step_m) ||
    !is.null(time_step)
  if (is.null(sampling)) sampling <- if (has_steps) "balanced" else "vertices"
  if (time_unit %in% c("months", "years") && !is.null(time_step) &&
        time_step != floor(time_step)) {
    .abort_track_sampling("Calendar steps require whole months or years.")
  }
  if (sampling %in% c("balanced", "merged_steps") && !has_steps) {
    .abort_track_sampling(
      "Balanced or merged_steps sampling needs an explicit spacing or step."
    )
  }
  if (sampling == "vertices" && has_steps) {
    .abort_track_sampling("Step limits cannot be used with vertices sampling.")
  }
  sampling
}


.track_geodesic <- function(lon, lat) {
  radians <- pi / 180
  xyz <- cbind(
    cos(lat * radians) * cos(lon * radians),
    cos(lat * radians) * sin(lon * radians),
    sin(lat * radians)
  )
  a <- xyz[-nrow(xyz), , drop = FALSE]
  b <- xyz[-1L, , drop = FALSE]
  cross <- cbind(
    a[, 2L] * b[, 3L] - a[, 3L] * b[, 2L],
    a[, 3L] * b[, 1L] - a[, 1L] * b[, 3L],
    a[, 1L] * b[, 2L] - a[, 2L] * b[, 1L]
  )
  angle <- atan2(sqrt(rowSums(cross^2)), rowSums(a * b))
  if (any(pi - angle < 1e-08)) {
    .abort_track_sampling(
      "Antipodal endpoints need an intermediate waypoint."
    )
  }
  # Each segment is a straight line from the origin of an azimuthal
  # equidistant projection centred on its start, i.e. the WGS84 geodesic.
  crs <- paste0(
    "+proj=aeqd +lat_0=", lat[-length(lat)], " +lon_0=", lon[-length(lon)],
    " +datum=WGS84 +units=m +no_defs"
  )
  xy <- lapply(
    seq_along(crs), function(i) {
      sf::sf_project("EPSG:4326", crs[i], cbind(lon[i + 1L], lat[i + 1L]))
    }
  )
  lengths <- vapply(xy, function(x) sqrt(sum(x^2)) / 1000, numeric(1))
  if (any(!is.finite(lengths))) {
    .abort_track_sampling(
      "A geodesic could not be constructed; add a waypoint."
    )
  }
  list(
    length_km = lengths, crs = crs, xy = xy,
    longitude = lon, latitude = lat
  )
}


# `value` shifted by whole calendar `months` in UTC, with the day clamped to
# the last day of a shorter month: 31 January + 1 month is 29 February.
.shift_temporal_months <- function(value, months) {
  parts <- as.POSIXlt(value, tz = "UTC")
  index <- (parts$year + 1900L) * 12L + parts$mon + months
  first <- as.Date(sprintf("%04d-%02d-01", index %/% 12L, index %% 12L + 1L))
  next_index <- index + 1L
  next_first <- as.Date(sprintf(
    "%04d-%02d-01", next_index %/% 12L, next_index %% 12L + 1L
  ))
  last_day <- as.integer(next_first - first)
  as.POSIXct(first + min(parts$mday, last_day) - 1L, tz = "UTC") +
    parts$hour * 3600 + parts$min * 60 + parts$sec
}


.track_time_breaks <- function(start, end, step, unit, budget) {
  if (is.null(step) || !is.finite(start) || !is.finite(end) || end == start) {
    return(c(0, 1))
  }
  # .track_sampling() has checked `step`, including whole calendar steps.
  if (unit %in% c("months", "years")) {
    # Break k is k steps from the start, not one step from the previous
    # break: end-of-month clamping (31 January + 1 month = 29 February)
    # would otherwise carry forward and every later break would drift (#146).
    first <- as.POSIXct(start, origin = "1970-01-01", tz = "UTC")
    values <- start
    k <- 1
    repeat {
      months <- if (unit == "years") 12 * k * step else k * step
      following <- as.numeric(.shift_temporal_months(first, months))
      if (length(values) >= budget) {
        .abort_track_sampling("time_step exceeds max_samples.")
      }
      values <- c(values, min(end, following))
      if (following >= end) break
      k <- k + 1
    }
    return((values - start) / (end - start))
  }
  seconds <- step * .time_unit_seconds(unit)
  n <- max(1, ceiling((end - start) / seconds))
  if (!is.finite(n) || n + 1 > budget) {
    .abort_track_sampling("time_step exceeds max_samples.")
  }
  seq(0, 1, length.out = n + 1)
}


# The samples of one line, from its vertices `v` (lon, lat, depth and `time`
# in seconds, NA without times) in order, with `options` from sample_3d()
# and at most `budget` points. Each segment is cut at fractions of its
# length: evenly for the distance and depth steps, at the time steps, or
# both for "merged_steps"; depth and time change linearly along it. Each
# point's support is half the length and duration of the cuts either side.
.track_sample_one <- function(v, options, budget) {
  arc <- .track_geodesic(v$lon, v$lat)
  depths <- v$depth
  times <- v$time
  dt <- diff(times)
  if (any(dt < 0, na.rm = TRUE)) {
    .abort_track_sampling("Dates must follow traversal order.")
  }
  pieces <- list()
  for (i in seq_along(dt)) {
    cuts <- c(0, 1)
    ends <- c(i, i + 1L)
    if (options$sampling != "vertices") {
      count <- if (is.null(options$spacing_km)) {
        1
      } else {
        max(1, ceiling(arc$length_km[i] / options$spacing_km))
      }
      if (!is.null(options$depth_step_m) && all(is.finite(depths[ends]))) {
        count <- max(
          count, ceiling(abs(diff(depths[ends])) / options$depth_step_m)
        )
      }
      if (!is.finite(count) || count + 1 > budget) {
        .abort_track_sampling("Track sampling exceeds max_samples.")
      }
      time_cuts <- .track_time_breaks(
        times[i], times[i + 1L], options$time_step,
        options$time_unit, budget
      )
      if (options$sampling == "balanced") {
        if (!is.null(options$time_step) &&
              options$time_unit %in% c("months", "years")) {
          # Preserve calendar anchors; subdivide each span only as needed.
          cuts <- 0
          for (j in seq_len(length(time_cuts) - 1L)) {
            first <- time_cuts[j]
            last <- time_cuts[j + 1L]
            horizontal <- if (is.null(options$spacing_km)) {
              1
            } else {
              arc$length_km[i] * (last - first) / options$spacing_km
            }
            vertical <- if (is.null(options$depth_step_m) ||
                              any(!is.finite(depths[ends]))) {
              1
            } else {
              abs(diff(depths[ends])) * (last - first) / options$depth_step_m
            }
            n <- max(1, ceiling(horizontal), ceiling(vertical))
            if (length(cuts) + n > budget) {
              .abort_track_sampling("Track sampling exceeds max_samples.")
            }
            cuts <- c(cuts, seq(first, last, length.out = n + 1L)[-1L])
          }
        } else {
          count <- max(count, length(time_cuts) - 1L)
          if (count + 1L > budget) {
            .abort_track_sampling("Track sampling exceeds max_samples.")
          }
          cuts <- seq(0, 1, length.out = count + 1L)
        }
      } else {
        cuts <- sort(unique(c(seq(0, 1, length.out = count + 1L), time_cuts)))
      }
    }
    fractions <- utils::head(cuts, -1L)
    width <- diff(cuts)
    used <- sum(vapply(pieces, nrow, integer(1)))
    if (length(fractions) + used + 1 > budget) {
      .abort_track_sampling("Track sampling exceeds max_samples.")
    }
    # Linear along the segment, with the vertex's own value at its start.
    along <- function(x) {
      values <- x[i] + fractions * (x[i + 1L] - x[i])
      values[fractions == 0] <- x[i]
      values
    }
    # Points along the segment's straight line in its projection; the
    # vertex keeps its own coordinates.
    xy <- sf::sf_project(arc$crs[i], "EPSG:4326",
                         outer(fractions, as.numeric(arc$xy[[i]])))
    xy[fractions == 0, ] <- c(arc$longitude[i], arc$latitude[i])
    pieces[[i]] <- data.frame(
      lon = xy[, 1L], lat = xy[, 2L], depth = along(depths),
      time = along(times), is_interpolated = fractions != 0,
      distance_km = sum(utils::head(arc$length_km, i - 1L)) +
        fractions * arc$length_km[i],
      support_km = width * arc$length_km[i], support_seconds = width * dt[i]
    )
  }
  n <- length(arc$longitude)
  last <- data.frame(
    lon = arc$longitude[n], lat = arc$latitude[n], depth = depths[n],
    time = times[n], is_interpolated = FALSE,
    distance_km = sum(arc$length_km), support_km = 0, support_seconds = 0
  )
  result <- do.call(rbind, c(pieces, list(last)))
  if (nrow(result) > budget) {
    .abort_track_sampling("Track sampling exceeds max_samples.")
  }
  # Each point stands for half the cut either side of it.
  for (name in c("support_km", "support_seconds")) {
    intervals <- utils::head(result[[name]], -1L)
    result[[name]] <- (c(0, intervals) + c(intervals, 0)) / 2
  }
  result$sample_id <- seq_len(nrow(result))
  result
}
