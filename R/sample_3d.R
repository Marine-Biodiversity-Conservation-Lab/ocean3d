# Sampling SpatVector3D lines into points (#105): the build step SPEC.md
# plans for tracks ("Extraction on the three steps", 1. Build). The sampling
# itself is extract_to_track()'s (R/extract_to_track.R): WGS84 geodesics
# between vertices, with depth and time interpolated along each segment.

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
#' @param sampling How to place the points, as in [extract_to_track()]:
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
#'   [extract_to_track()], which samples tracks the same way.
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

  vertices <- data.frame(track_id = v$feature_id, lon = v$lon, lat = v$lat,
                         depth = v$depth, time = v$time,
                         stringsAsFactors = FALSE)
  options <- list(
    sampling = sampling, spacing_km = spacing_km, depth_step_m = depth_step_m,
    time_step = time_step, time_unit = time_unit, interpolate_depth = TRUE,
    interpolate_time = timed, land_mask = NULL, land_policy = "skip"
  )
  s <- .track_sample_all(
    vertices, options, "lon", "lat", "depth", "time",
    time_values = if (timed) as.numeric(v$time), depth_used = TRUE,
    time_used = timed, max_samples = max_samples
  )$samples

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


# Shift times by `amount` units: elapsed seconds, or UTC calendar months and
# years, clamped to the end of a shorter month. Used by sample_3d() and the
# track functions.
.buffer_shift_time <- function(values, amount, unit) {
  # Months/years are UTC calendar shifts, with end-of-month clamping.
  # Reuse the scalar calendar helper already used by point diagnostics.
  if (unit %in% c("months", "years")) {
    months <- if (unit == "years") amount * 12 else amount
    shifted <- values
    for (i in which(!is.na(values))) {
      shifted[i] <- tryCatch(
        .shift_temporal_months(values[i], months),
        error = function(e) {
          stop(
            paste0("The requested time buffer exceeds the supported ",
                   "calendar range. "),
            "Check `time_before`, `time_after` and `time_unit`.",
            call. = FALSE
          )
        }
      )
    }
  } else {
    shifted <- values + amount * .time_unit_seconds(unit)
  }

  if (any(!is.na(values) & !is.finite(as.numeric(shifted)))) {
    stop(
      paste0("The requested time buffer exceeds the supported ",
             "calendar range. "),
      "Check `time_before`, `time_after` and `time_unit`.",
      call. = FALSE
    )
  }
  shifted
}
