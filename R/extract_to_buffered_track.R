#-------------------------------------------------------------------------------
# Environmental buffers along ordered tracks
#-------------------------------------------------------------------------------

#' Extract environmental summaries in buffers along tracks
#'
#' Prepare a trajectory using the same sampling and marine-mask rules as
#' [extract_to_track()], then apply [extract_to_buffer()] at each position.
#'
#' @inheritParams extract_to_track
#' @inheritParams extract_to_buffer
#' @param buffer_weighting Weighting within each buffer: weighted or equal.
#' @param depth_step_m Optional maximum depth change between inserted samples,
#'   in metres. Changing depths also require interpolate_depth = TRUE.
#' @param method Point-axis depth selection: 2d, nearest, surface or bottom.
#'   NULL chooses per variable. The all method is unavailable for buffers.
#' @param min_valid_fraction Minimum valid support within each buffer. NULL
#'   uses zero for most_similar and 0.8 for other statistics.
#' @param track_min_valid_fraction Coverage threshold for the track summary.
#'   min_valid_fraction controls each individual buffer instead.
#' @param interpolate_target Allow linear interpolation of target columns
#'   between vertices. This permission is independent of depth/time
#'   interpolation.
#'   Constant targets can be retained without it. Original NA targets stay NA.
#' @param buffer_time_unit Units of time_before and time_after around each
#'   sample.
#'   time_unit controls trajectory sampling steps instead.
#'
#' @returns A list containing samples, summary and settings. Each sample has
#'   buffer diagnostics, including all winning locations for most_similar.
#'   Summaries describe the sequence of buffer results, not a union of buffers.
#'   Overlapping buffers can reuse the same environmental records.
#'   Summary counts n_observed and n_assumed distinguish original vertices
#'   from inserted positions; n_samples is their total, including missing
#'   buffer results. settings$sampling_counts contains these counts once per
#'   track, also when statistics = NULL.
#'   Track summary variation_status follows extract_to_track(): not_requested,
#'   available, insufficient_samples, or not_available. It describes SD/SE
#'   requested through statistics along the track, not the local buffer
#'   statistic. Consult status for coverage and support problems.
#'   most_similar adds the buffer engine's _match_lon, _match_lat, _match_depth
#'   and _match_time columns to samples. Ambiguous fields and fields on invalid
#'   connections are NA; original route coordinates are retained.
#'
#' @details No environmental interpolation or automatic marine routing is done.
#'   most_similar selects by absolute difference in the variable's units; it
#'   does not estimate an animal's route or replace the original coordinates.
#'   Equal weighting in cells mode excludes zero-support original vertices.
#' @seealso [extract_to_track()], [extract_to_buffer()]
#' @examples
#' nc_file <- system.file("extdata", "example_2d.nc", package = "ocean3d")
#' track <- data.frame(lon = c(0, 1), lat = 40,
#'   date = as.Date("2020-01-01"))
#' result <- extract_to_buffered_track(track, nc_file, var = "temp",
#'   horizontal_radius_km = 1, infer_bounds = TRUE, quiet = TRUE)
#' result$summary
#' @export
extract_to_buffered_track <- function(
  data, nc, var = NULL, track_col = NULL,
  sampling = c("vertices", "balanced", "merged_steps", "cells"),
  spacing_km = NULL,
  interpolate_depth = FALSE, interpolate_time = FALSE,
  interpolate_target = FALSE,
  depth_step_m = NULL, time_step = NULL, time_unit = c(
    "hours", "minutes",
    "seconds", "days", "months", "years"
  ), horizontal_radius_km = 0,
  depth_radius_m = 0, lon_cell_radius = NULL, lat_cell_radius = NULL,
  depth_layer_radius = NULL, time_before = 0, time_after = 0,
  buffer_time_unit = c(
    "days",
    "hours", "minutes", "seconds", "months", "years"
  ), statistic = "mean",
  buffer_weighting = c("weighted", "equal"), min_valid_fraction = NULL,
  track_min_valid_fraction = 0.8, infer_bounds = FALSE, method = NULL,
  statistics = "mean", weighting = c("distance", "duration", "equal"),
  target_col = NULL, max_target_difference = NULL, target_ties = c(
    "NA",
    "lower", "higher"
  ), lon_col = "lon", lat_col = "lat", depth_col = "depth",
  date_col = "date", output_col = NULL, time_resolution = c(
    "auto",
    "datetime", "daily", "monthly", "yearly"
  ), max_samples = 100000L,
  max_segment_km = NULL, max_gap_hours = NULL, z_as_depth = FALSE,
  distance_model = c("WGS84", "sphere"), land_mask = NULL, marine_values = NULL,
  mask_type = c("land", "marine"), mask_domain = NULL, land_policy = c(
    "skip",
    "error"
  ), mask_step_km = 0.5, quiet = FALSE, ...
) {
  .track_flag(interpolate_target, "interpolate_target")
  if (identical(method, "all")) {
    stop("Use one buffer depth method; 'all' is unavailable.", call. = FALSE)
  }
  time_unit <- match.arg(time_unit)
  buffer_time_unit <- match.arg(buffer_time_unit)
  buffer_weighting <- match.arg(buffer_weighting)
  weighting <- match.arg(weighting)
  time_resolution <- match.arg(time_resolution)
  distance_model <- match.arg(distance_model)
  mask_type <- match.arg(mask_type)
  land_policy <- match.arg(land_policy)
  target_ties <- match.arg(target_ties)
  preparation <- list(
    data = data, nc = nc, var = var, track_col = track_col,
    spacing_km = spacing_km, interpolate_depth = interpolate_depth,
    interpolate_time = interpolate_time, depth_step_m = depth_step_m,
    time_step = time_step, time_unit = time_unit, method = method,
    statistics = statistics, weighting = weighting,
    min_valid_fraction = track_min_valid_fraction,
    lon_col = lon_col, lat_col = lat_col, depth_col = depth_col,
    date_col = date_col, output_col = output_col,
    time_resolution = time_resolution,
    max_samples = max_samples, max_segment_km = max_segment_km,
    max_gap_hours = max_gap_hours,
    z_as_depth = z_as_depth, distance_model = distance_model,
    land_mask = land_mask,
    marine_values = marine_values, mask_type = mask_type,
    mask_domain = mask_domain,
    land_policy = land_policy, mask_step_km = mask_step_km, quiet = quiet
  )
  if (!missing(sampling)) {
    preparation$sampling <- match.arg(sampling)
  }
  dots <- list(...)
  prepare_names <- c(
    "file_col", "lon_dim", "lat_dim", "depth_dim",
    "time_dim", "time_match", "max_time_diff", "spatial_match",
    "max_distance_km"
  )
  if (
      length(dots) && (
                       is.null(names(dots)) || any(
                         !names(dots) %in%
                           prepare_names
                       ))) {
    stop("Use supported named buffer selection options in ... .")
  }
  context <- do.call(.prepare_track, c(preparation, dots))
  samples <- context$samples
  variables <- unique(
    vapply(
      context$metadata, function(x) x$variable,
      character(1)
    )
  )
  if (identical(statistic, "most_similar")) {
    targets <- .target_columns(context$vertices, variables, target_col)
    for (column in unique(unname(targets))) {
      values <- context$vertices[[column]]
      samples[[column]] <- vapply(
        seq_len(nrow(samples)), function(i) {
          if (!is.na(samples$source_row[i])) {
            return(values[samples$source_row[i]])
          }
          rows <- which(context$vertices$track_id == samples$track_id[i])
          .track_interpolate(
            values[rows], samples$segment_id[i],
            samples$segment_fraction[i], interpolate_target,
            "interpolate_target",
            TRUE
          )
        }, numeric(1)
      )
    }
  }
  resolution <- if (is.null(context$prepared)) {
    time_resolution
  } else {
    context$prepared$resolution
  }
  if (!is.null(context$prepared)) {
    .warn_time_zone(context$prepared, paste0("`data$", date_col, "`"))
  }
  output <- .buffer_output_names(
    variables, output_col, names(samples),
    TRUE, locations = statistic == "most_similar"
  )
  samples <- withCallingHandlers(do.call(
    extract_to_buffer, c(
      list(
        data = samples, nc = nc,
        var = variables, horizontal_radius_km = horizontal_radius_km,
        depth_radius_m = depth_radius_m, lon_cell_radius = lon_cell_radius,
        lat_cell_radius = lat_cell_radius,
        depth_layer_radius = depth_layer_radius,
        time_before = time_before, time_after = time_after,
        time_unit = buffer_time_unit,
        statistic = statistic, weighting = buffer_weighting,
        min_valid_fraction = min_valid_fraction,
        infer_bounds = infer_bounds, method = method, target_col = target_col,
        max_target_difference = max_target_difference,
        target_ties = target_ties,
        lon_col = lon_col, lat_col = lat_col, depth_col = depth_col,
        date_col = date_col, output_col = output_col,
        time_resolution = resolution,
        return_match_info = TRUE, quiet = quiet
      ), dots
    )
  ), ocean3d_target_location_warning = function(warning) {
    invokeRestart("muffleWarning")
  })
  samples <- .track_mask_invalid(samples, unname(output))
  if (statistic == "most_similar") {
    for (column in unname(output)) {
      name <- paste0(column, "_match_info")
      invalid <- vapply(
        samples[[name]], function(x) identical(x$status, "invalid_connection"),
        logical(1)
      )
      selected_locations <- .buffer_target_locations(
        samples[[name]], samples[[column]], column, eligible = !invalid
      )
      samples[[name]] <- selected_locations$diagnostics
      samples[names(selected_locations$columns)] <- selected_locations$columns
    }
  }
  summary <- .track_summaries(
    samples, unname(output), statistics,
    weighting, track_min_valid_fraction
  )
  settings <- .track_settings(
    context, samples, statistics = statistics,
    statistic = statistic,
    buffer_weighting = buffer_weighting,
    track_weighting = weighting,
    track_min_valid_fraction = track_min_valid_fraction,
    buffer_min_valid_fraction = if (is.null(min_valid_fraction)) {
      if (statistic == "most_similar") 0 else 0.8
    } else {
      min_valid_fraction
    },
    interpolate_target = interpolate_target,
    mask_checked = !is.null(land_mask),
    buffer = list(
      horizontal_radius_km = horizontal_radius_km,
      depth_radius_m = depth_radius_m, lon_cell_radius = lon_cell_radius,
      lat_cell_radius = lat_cell_radius,
      depth_layer_radius = depth_layer_radius,
      time_before = time_before, time_after = time_after,
      time_unit = buffer_time_unit, infer_bounds = infer_bounds,
      method = method
    ),
    target = list(
      columns = target_col, ties = target_ties,
      max_difference = max_target_difference,
      interpolate = interpolate_target
    )
  )
  list(samples = samples, summary = summary, settings = settings)
}
