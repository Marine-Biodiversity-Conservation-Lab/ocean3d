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
#' [buffer_3d()], so a track is handled as a set of points or buffers.
#'
#' Points between vertices follow the WGS84 geodesic (the shortest path on
#' the ellipsoid), so a segment across the antimeridian goes the short way.
#' Depth and time change linearly with distance along each segment.
#' Distances are horizontal, not three-dimensional.
#'
#' @param x A [SpatVector3D-class] of lines, from [as_vector3d()].
#' @param sampling How to place the points, as in [extract_to_track()]:
#'   `"vertices"` keeps only the vertices; `"balanced"` divides each segment
#'   evenly so that every step limit below holds; `"merged_steps"` combines
#'   the distance and depth cuts with the time cuts. `NULL` (default) gives
#'   `"balanced"` when a step is given and `"vertices"` otherwise.
#' @param spacing_km Largest horizontal distance between points, in
#'   kilometres.
#' @param depth_step_m Largest change in depth between points, in metres.
#' @param time_step,time_unit Largest change in time between points, and its
#'   unit: seconds, minutes, hours, days, or whole calendar months or years.
#'   Needs a time at every vertex. Calendar steps count from the first vertex
#'   of each segment, so a route that starts on 31 January with monthly steps
#'   gets 29 February, 31 March, 30 April.
#' @param max_samples Most points to make, in total.
#'
#' @return A [SpatVector3D-class] of points, one per sample, in order along
#'   each line. Each point's feature ID is the line's ID and the sample
#'   number (`"A_3"`). Columns: `track_id` (the line's ID), `sample_id`
#'   (1, 2, ... along the line), `distance_km` (distance from the line's
#'   start), `support_km` and `support_seconds` (half the length and
#'   duration of the segments either side, as weights for summaries along
#'   the track; `NA` without times), `is_interpolated` (`FALSE` for the
#'   line's own vertices), and the line's columns that have one value along
#'   it. With dates and no hours, each point's time is the UTC day it falls
#'   on.
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
  if (!methods::is(x, "SpatVector3D") || x@type != "lines") {
    .abort_vector3d(paste0(
      "`x` must be SpatVector3D lines, from as_vector3d(type = \"lines\"). ",
      "Points are already samples, and polygons have no path to follow."
    ))
  }
  if (!is.null(sampling)) {
    sampling <- match.arg(sampling, c("vertices", "balanced", "merged_steps"))
  }
  time_unit <- match.arg(time_unit)
  sampling <- .track_sampling(sampling, spacing_km, depth_step_m, time_step,
                              time_unit)
  .track_positive(max_samples, "max_samples", optional = FALSE)
  v <- x@vertices
  ids <- unique(v$feature_id)
  v <- v[order(match(v$feature_id, ids), v$vertex), , drop = FALSE]
  timed <- !all(is.na(v$time))
  if (!is.null(time_step) && !timed) {
    .abort_vector3d(
      "`time_step` needs a time at each vertex; give `time` in as_vector3d()."
    )
  }

  vertices <- data.frame(track_id = v$feature_id, lon = v$lon, lat = v$lat,
                         depth = v$depth, time = v$time,
                         stringsAsFactors = FALSE)
  options <- list(
    sampling = sampling, spacing_km = spacing_km, depth_step_m = depth_step_m,
    time_step = time_step, time_unit = time_unit, interpolate_depth = TRUE,
    interpolate_time = timed, land_mask = NULL, land_policy = "skip"
  )
  time_values <- if (timed) as.numeric(v$time)
  remaining <- max_samples
  pieces <- lapply(ids, function(id) {
    sampled <- .track_sample_one(
      vertices, which(vertices$track_id == id), options, "lon", "lat",
      "depth", "time", time_values, depth_used = TRUE, time_used = timed,
      budget = remaining
    )$samples
    remaining <<- remaining - nrow(sampled)
    sampled
  })
  s <- do.call(rbind, pieces)

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
  extra <- .feature_attributes(v)
  extra <- extra[match(points$track_id, ids), , drop = FALSE]
  .new_vector3d("points", x@time_resolution,
                vertices = .add_attributes(points, extra))
}
