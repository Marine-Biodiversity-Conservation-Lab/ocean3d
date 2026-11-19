# Buffering SpatVector3D points and lines into polygons (#114): the 3D
# counterpart of sf::st_buffer(), with a depth window and time window.

#' Buffer 3D points or lines into polygons
#'
#' Draws a polygon around each point or line of a [SpatVector3D-class], with
#' a depth range and an optional time interval, as `sf::st_buffer()` does in
#' 2D. A detection at 35 m with `radius_km = 10` and `depth_window = 25`
#' becomes a circle of 10 km radius from 10 to 60 m deep.
#'
#' The outline is a circle (for a line, the area within `radius_km` of it)
#' measured on the WGS84 ellipsoid, made with `terra::buffer()`. Outlines
#' that cross the 180th meridian (the antimeridian) are split there into
#' two parts.
#'
#' Each polygon has one depth range and one time interval, so a line's
#' polygon covers the line's whole depth and time range. Depth runs from
#' `depth_window` above the shallowest vertex to `depth_window` below the
#' deepest, and time from `time_window` before the first time to
#' `time_window` after the last. For a dive from the surface to 300 m, the
#' polygon spans 0 m to 300 m plus `depth_window` along the whole track. To
#' follow the dive more closely, split the line into shorter lines or points
#' and buffer those.
#'
#' @param x A [SpatVector3D-class] of points or lines, from [as_vector3d()].
#' @param radius_km Horizontal radius in kilometres: one positive number,
#'   below 10000.
#' @param depth_window Depth in metres added above and below each feature's
#'   depths: one number, 0 or more. The top stops at the surface (0 m).
#'   With `0`, the polygon spans just the feature's own depth range.
#' @param time_window Optional time added before and after each feature's
#'   times: a number of days, or a `difftime` (e.g. `as.difftime(6, units =
#'   "hours")`). Every vertex needs a time. If `x` has dates without hours,
#'   the window must be whole days, and the interval covers whole UTC days as
#'   in [as_vector3d()]. With `NULL` (default), each polygon keeps its
#'   feature's own times: a point's time, a line's first to last time, or no
#'   time. A line with times at some vertices but not others is an error
#'   either way.
#'
#' @return A [SpatVector3D-class] of polygons, one per feature of `x`, with
#'   the same IDs. For points, the other columns of `x` are kept. For lines,
#'   columns with one value along each line (such as a tag name) are kept, and
#'   columns that change from vertex to vertex (such as speed) are dropped.
#' @seealso [as_vector3d()], [geomtype()]
#' @examples
#' obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(35, 120),
#'                   date = as.Date(c("2020-01-01", "2020-01-02")))
#' detections <- as_vector3d(obs, time = "date")
#'
#' # 10 km circles, 25 m above and below each detection, a day either side.
#' buffers <- buffer_3d(detections, radius_km = 10, depth_window = 25,
#'                      time_window = 1)
#' as.data.frame(buffers)[c("feature_id", "depth_min", "depth_max",
#'                          "time_start", "time_end")]
#'
#' # A corridor 5 km either side of a dive track.
#' dive <- data.frame(tag = "A", lon = c(150, 150.1, 150.2), lat = -20,
#'                    depth = c(0, 120, 15))
#' track <- as_vector3d(dive, type = "lines", id = "tag")
#' buffer_3d(track, radius_km = 5, depth_window = 10)
#' @export
buffer_3d <- function(x, radius_km, depth_window, time_window = NULL) {
  if (!methods::is(x, "SpatVector3D") || x@type == "polygons") {
    .abort_vector3d(paste0(
      "`x` must be SpatVector3D points or lines, from as_vector3d(). ",
      "Polygons already have a depth range."
    ))
  }
  if (length(x) == 0L) {
    .abort_vector3d("`x` has no points or lines to buffer.")
  }
  if (!is.numeric(radius_km) || length(radius_km) != 1L ||
        !is.finite(radius_km) || radius_km <= 0 || radius_km >= 10000) {
    .abort_vector3d(
      "`radius_km` must be one positive number of kilometres, below 10000."
    )
  }
  if (!is.numeric(depth_window) || length(depth_window) != 1L ||
        !is.finite(depth_window) || depth_window < 0) {
    .abort_vector3d("`depth_window` must be one number of metres, 0 or more.")
  }

  # Vertices in feature order, then in order along each line: the geometry
  # below pairs them up row by row.
  v <- x@vertices
  ids <- unique(v$feature_id)
  v <- v[order(match(v$feature_id, ids), v$vertex), , drop = FALSE]
  feature <- factor(v$feature_id, levels = ids)
  feature_min <- function(values) {
    vapply(split(values, feature), min, numeric(1), USE.NAMES = FALSE)
  }
  feature_max <- function(values) {
    vapply(split(values, feature), max, numeric(1), USE.NAMES = FALSE)
  }

  # Each feature's earliest and latest time, or NA for a feature without
  # times. A line timed at only some vertices has no clear time range.
  timed <- !is.na(v$time)
  timed_count <- tabulate(feature[timed], nbins = length(ids))
  vertex_count <- tabulate(feature, nbins = length(ids))
  partial <- ids[timed_count > 0L & timed_count < vertex_count]
  if (length(partial)) {
    noun <- if (length(partial) == 1L) "Line" else "Lines"
    verb <- if (length(partial) == 1L) "has" else "have"
    .abort_vector3d(glue::glue(
      "{noun} {toString(partial)} of `x` {verb} a time at some vertices but ",
      "not others. Give every vertex a time, or none."
    ))
  }
  seconds_utc <- as.numeric(v$time)
  start <- feature_min(ifelse(timed, seconds_utc, Inf))
  end <- feature_max(ifelse(timed, seconds_utc, -Inf))
  start[is.infinite(start)] <- NA
  end[is.infinite(end)] <- NA
  start <- .POSIXct(start, tz = "UTC")
  end <- .POSIXct(end, tz = "UTC")
  if (!is.null(time_window)) {
    seconds <- if (inherits(time_window, "difftime")) {
      as.numeric(time_window, units = "secs")
    } else if (is.numeric(time_window)) {
      time_window * 86400
    }
    if (length(seconds) != 1L || !is.finite(seconds) || seconds < 0) {
      .abort_vector3d(
        "`time_window` must be one number of days, 0 or more, or a difftime."
      )
    }
    untimed <- ids[timed_count == 0L]
    if (length(untimed)) {
      noun <- if (length(untimed) == 1L) "Feature" else "Features"
      verb <- if (length(untimed) == 1L) "has" else "have"
      .abort_vector3d(glue::glue(
        "`time_window` needs a time at every vertex. {noun} ",
        "{toString(untimed)} of `x` {verb} no time; give `time` in ",
        "as_vector3d()."
      ))
    }
    if (identical(x@time_resolution, "daily")) {
      days <- seconds / 86400
      if (abs(days - round(days)) > 1e-6) {
        .abort_vector3d(paste0(
          "`x` has dates without hours, so `time_window` must be whole ",
          "days. Give the times with hours in as_vector3d() for a finer ",
          "window."
        ))
      }
      seconds <- round(days) * 86400
    }
    start <- start - seconds
    end <- end + seconds
  }

  # A line that crosses the antimeridian (179.9 to -179.9) is unwrapped so
  # each step is the short way round; terra would otherwise buffer it the
  # long way. Points are wrapped into [-180, 180).
  shapes <- if (x@type == "points") {
    terra::vect(cbind(.normalise_longitude(v$lon), v$lat), crs = "EPSG:4326")
  } else {
    lon <- unlist(lapply(split(v$lon, feature), .unwrap_longitude),
                  use.names = FALSE)
    terra::vect(cbind(object = as.integer(feature), part = 1L, x = lon,
                      y = v$lat), type = "lines", crs = "EPSG:4326")
  }
  outlines <- terra::buffer(shapes, radius_km * 1000, quadsegs = 180L)

  features <- data.frame(
    feature_id = ids,
    depth_min = pmax(0, feature_min(v$depth) - depth_window),
    depth_max = feature_max(v$depth) + depth_window,
    time_start = start,
    time_end = end,
    stringsAsFactors = FALSE
  )
  features$geometry <- sf::st_geometry(sf::st_as_sf(outlines))
  # Other columns: kept if they have one value per feature.
  extra <- v[setdiff(names(v), c("feature_id", "vertex", "lon", "lat",
                                 "depth", "time"))]
  constant <- vapply(extra, function(column) {
    all(tapply(column, feature, function(values) {
      length(unique(values)) == 1L
    }))
  }, logical(1))
  extra <- extra[!duplicated(feature), constant, drop = FALSE]
  .new_vector3d("polygons", x@time_resolution,
                features = .add_attributes(features, extra))
}
