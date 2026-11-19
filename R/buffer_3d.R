# Buffering SpatVector3D points and lines into polygons (#114): the 3D
# counterpart of sf::st_buffer(), with a depth window and time window.

#' Buffer 3D points or lines into polygons
#'
#' Draws a polygon around each point or line of a [SpatVector3D-class], with
#' a depth range and an optional time interval, as `sf::st_buffer()` does in
#' 2D. A detection at 35 m buffered by 10 km and 25 m becomes a 10 km circle
#' from 10 to 60 m.
#'
#' Each point or line becomes one polygon with the same ID. The outline is a
#' geodesic circle (or, for a line, the area within `radius_km` of it) on the
#' WGS84 ellipsoid, as in [extract_to_buffer()]; outlines that cross the
#' antimeridian are split there.
#'
#' A polygon has one depth range and one time interval, so a line's polygon
#' covers its whole depth and time range: from `depth_window` above its
#' shallowest vertex to `depth_window` below its deepest, and from
#' `time_window` before its first time to `time_window` after its last. For
#' a dive that goes from the surface to 300 m, that is the whole water column
#' along the track. To follow the dive more closely, split the line into
#' shorter lines or points first and buffer those.
#'
#' @param x A [SpatVector3D-class] of points or lines, from [as_vector3d()].
#' @param radius_km Horizontal radius in kilometres: one positive number,
#'   below 10000.
#' @param depth_window Depth in metres added above and below each feature's
#'   depths: one non-negative number. The top stops at the surface (0 m).
#'   With `0`, the polygon spans just the feature's own depth range.
#' @param time_window Optional time added before and after each feature's
#'   times: a number of days, or a `difftime` (e.g. `as.difftime(6, units =
#'   "hours")`). Needs a time at every vertex. With daily times, it must be
#'   whole days, and the interval covers whole UTC days as in
#'   [as_vector3d()]. `NULL` (default) keeps each feature's own times: a
#'   point's polygon is at the point's time, a line's polygon spans its first
#'   to last time, and features without time give polygons without time.
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
  if (!is.numeric(radius_km) || length(radius_km) != 1L ||
        !is.finite(radius_km) || radius_km <= 0 || radius_km >= 10000) {
    .abort_vector3d(
      "`radius_km` must be one positive number of kilometres, below 10000."
    )
  }
  if (!is.numeric(depth_window) || length(depth_window) != 1L ||
        !is.finite(depth_window) || depth_window < 0) {
    .abort_vector3d(
      "`depth_window` must be one non-negative number of metres."
    )
  }

  v <- x@vertices
  ids <- unique(v$feature_id)
  feature <- factor(v$feature_id, levels = ids)
  first <- function(values) unname(sapply(split(values, feature), min))
  last <- function(values) unname(sapply(split(values, feature), max))

  # Each feature's first and last time, NA for a feature without times.
  timed <- !is.na(v$time)
  seconds_utc <- as.numeric(v$time)
  start <- first(ifelse(timed, seconds_utc, Inf))
  end <- last(ifelse(timed, seconds_utc, -Inf))
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
      .abort_vector3d(paste0(
        "`time_window` must be one non-negative number of days, or a ",
        "difftime."
      ))
    }
    if (!all(timed)) {
      missing <- unique(v$feature_id[!timed])
      which <- if (length(missing) == 1L) "Feature" else "Features"
      verb <- if (length(missing) == 1L) "has" else "have"
      .abort_vector3d(glue::glue(
        "`time_window` needs a time at every vertex. {which} ",
        "{toString(missing)} of `x` {verb} none; give `time` in as_vector3d()."
      ))
    }
    if (identical(x@time_resolution, "daily") && seconds %% 86400 != 0) {
      .abort_vector3d(paste0(
        "`x` has dates without hours, so `time_window` must be whole days. ",
        "Give the times with hours in as_vector3d() for a finer window."
      ))
    }
    start <- start - seconds
    end <- end + seconds
  }

  # A line that crosses the antimeridian (179.9 to -179.9) is unwrapped so
  # each step is the short way round; terra would otherwise buffer it the
  # long way. Points are wrapped into [-180, 180) as in extract_to_buffer().
  lon <- if (x@type == "points") {
    .normalise_longitude(v$lon)
  } else {
    unlist(lapply(split(v$lon, feature), function(values) {
      values[1L] + cumsum(c(0, (diff(values) + 180) %% 360 - 180))
    }), use.names = FALSE)
  }
  shapes <- if (x@type == "points") {
    terra::vect(cbind(lon, v$lat), crs = "EPSG:4326")
  } else {
    terra::vect(cbind(object = as.integer(feature), part = 1L, x = lon,
                      y = v$lat), type = "lines", crs = "EPSG:4326")
  }
  outlines <- terra::buffer(shapes, radius_km * 1000, quadsegs = 180L)

  features <- data.frame(
    feature_id = ids,
    depth_min = pmax(0, first(v$depth) - depth_window),
    depth_max = last(v$depth) + depth_window,
    time_start = start,
    time_end = end,
    stringsAsFactors = FALSE
  )
  features$geometry <- sf::st_geometry(sf::st_as_sf(outlines))
  extra <- v[setdiff(names(v), c("feature_id", "vertex", "lon", "lat",
                                 "depth", "time"))]
  constant <- vapply(extra, function(column) {
    all(tapply(column, feature, function(values) {
      length(unique(values)) == 1L
    }))
  }, logical(1))
  extra <- extra[!duplicated(feature), constant, drop = FALSE]
  .new_vector3d("polygons", x@time_resolution,
                features = cbind(features, extra))
}
