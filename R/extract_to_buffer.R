# -----------------------------------------------------------------------------
# Buffer extraction from netCDF files
# -----------------------------------------------------------------------------

.validate_buffer_options <- function(
  horizontal_radius_km = 0,
  depth_radius_m = 0,
  time_before = 0,
  time_after = 0,
  time_unit = c("days", "hours", "minutes", "seconds", "months", "years"),
  statistic = "mean",
  weighting = c("weighted", "equal"),
  min_valid_fraction = 0.8,
  infer_bounds = FALSE,
  method = NULL,
  return_match_info = TRUE,
  quiet = FALSE
) {
  sizes <- list(
    horizontal_radius_km = horizontal_radius_km,
    depth_radius_m = depth_radius_m,
    time_before = time_before,
    time_after = time_after
  )

  for (name in names(sizes)) {
    value <- sizes[[name]]
    if (!is.numeric(value) || length(value) != 1L ||
          is.na(value) || !is.finite(value) || value < 0) {
      stop(
        "`", name, "` must be one finite, non-negative number.",
        call. = FALSE
      )
    }
  }

  if (!any(unlist(sizes) > 0)) {
    stop(
      "At least one buffer window must be positive. ",
      "Use extract_to_point() for point selection.",
      call. = FALSE
    )
  }

  time_unit <- match.arg(time_unit)
  if (time_unit %in% c("months", "years")) {
    for (name in c("time_before", "time_after")) {
      value <- sizes[[name]]
      if (value != floor(value)) {
        stop(
          "`", name, "` must be a whole number when `time_unit` is '",
          time_unit, "'. Use days or hours for fractional durations.",
          call. = FALSE
        )
      }
    }
  }

  allowed_statistics <- c("mean", "sd", "se", "min", "max")
  if (!is.character(statistic) || length(statistic) != 1L ||
        is.na(statistic) || !statistic %in% allowed_statistics) {
    stop("`statistic` must be one of: ",
         paste(allowed_statistics, collapse = ", "), ".", call. = FALSE)
  }
  weighting <- match.arg(weighting)

  if (!is.numeric(min_valid_fraction) ||
        length(min_valid_fraction) != 1L ||
        is.na(min_valid_fraction) || !is.finite(min_valid_fraction) ||
        min_valid_fraction < 0 || min_valid_fraction > 1) {
    stop(
      "`min_valid_fraction` must be one number between 0 and 1.",
      call. = FALSE
    )
  }

  switches <- list(
    infer_bounds = infer_bounds,
    return_match_info = return_match_info,
    quiet = quiet
  )
  for (name in names(switches)) {
    value <- switches[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }

  if (depth_radius_m > 0) {
    if (!is.null(method)) {
      stop(
        "`depth_radius_m` defines a depth interval. ",
        "Leave `method = NULL` instead of selecting one depth layer.",
        call. = FALSE
      )
    }
    resolved_method <- "interval"
  } else {
    resolved_method <- if (is.null(method)) {
      "nearest"
    } else {
      match.arg(method, c("nearest", "surface", "bottom", "2d"))
    }
  }

  c(
    sizes,
    list(
      time_unit = time_unit,
      statistic = statistic,
      weighting = weighting,
      min_valid_fraction = min_valid_fraction,
      infer_bounds = infer_bounds,
      method = resolved_method,
      return_match_info = return_match_info,
      quiet = quiet
    )
  )
}


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
            "The requested time buffer exceeds the supported calendar range. ",
            "Check `time_before`, `time_after` and `time_unit`.",
            call. = FALSE
          )
        }
      )
    }
  } else {
    seconds <- switch(
      unit, days = 86400, hours = 3600, minutes = 60, seconds = 1
    )
    shifted <- values + amount * seconds
  }

  if (any(!is.na(values) & !is.finite(as.numeric(shifted)))) {
    stop(
      "The requested time buffer exceeds the supported calendar range. ",
      "Check `time_before`, `time_after` and `time_unit`.",
      call. = FALSE
    )
  }
  shifted
}


.buffer_observation_windows <- function(
  data,
  options,
  depth_col = "depth",
  date_col = "date",
  id_col = NULL
) {
  # Input has already been normalised to a table by the public entry point.
  .check_required_cols(data, id_col)

  n <- nrow(data)
  windows <- data.frame(
    row = seq_len(n),
    depth_lower_m = rep(NA_real_, n),
    depth_upper_m = rep(NA_real_, n),
    time_lower_utc = as.POSIXct(rep(NA_character_, n), tz = "UTC"),
    time_upper_utc = as.POSIXct(rep(NA_character_, n), tz = "UTC")
  )
  prepared_time <- NULL

  if (options$depth_radius_m > 0) {
    if (is.null(depth_col)) {
      stop(
        "A depth buffer requires `depth_col` naming an observation column.",
        call. = FALSE
      )
    }
    .check_required_cols(data, depth_col)
    depth <- data[[depth_col]]

    if (!is.numeric(depth)) {
      stop(
        "`data$", depth_col, "` must contain depths in metres, ",
        "positive downward. Explicit NA values are allowed.",
        call. = FALSE
      )
    }

    bad <- which(!is.na(depth) & (!is.finite(depth) | depth < 0))
    if (length(bad)) {
      stop(
        "`data$", depth_col, "` contains invalid observation depths. ",
        "Use finite, non-negative metres or NA. Affected rows: ",
        length(bad), "; first rows: ", toString(head(bad, 5L)), ".",
        call. = FALSE
      )
    }

    windows$depth_lower_m <- depth - options$depth_radius_m
    windows$depth_upper_m <- depth + options$depth_radius_m
    if (any(!is.na(depth) & !is.finite(windows$depth_upper_m))) {
      stop("The requested depth buffer exceeds the numeric range.",
           call. = FALSE)
    }
  }

  if (options$time_before > 0 || options$time_after > 0) {
    if (is.null(date_col)) {
      stop(
        "A time buffer requires `date_col` naming an observation column.",
        call. = FALSE
      )
    }
    .check_required_cols(data, date_col)

    prepared_time <- .prepare_time_input(
      data[[date_col]],
      input_name = paste0("`data$", date_col, "`"),
      ids = if (is.null(id_col)) NULL else data[[id_col]]
    )

    windows$time_lower_utc <- .buffer_shift_time(
      prepared_time$values, -options$time_before, options$time_unit
    )
    windows$time_upper_utc <- .buffer_shift_time(
      prepared_time$values, options$time_after, options$time_unit
    )
  }

  # The public entry point will emit the timezone warning once per extraction.
  list(windows = windows, prepared_time = prepared_time)
}


.buffer_read_bounds <- function(nc, coordinate, n) {
  name <- .read_time_attribute(nc, coordinate, "bounds")
  if (is.null(name)) return(NULL)

  if (!is.character(name) || length(name) != 1L ||
        is.na(name) || !nzchar(name) || !name %in% names(nc$var)) {
    stop("Invalid or missing bounds variable for `", coordinate, "`.",
         call. = FALSE)
  }

  dimension <- if (coordinate %in% names(nc$dim)) {
    coordinate
  } else {
    dims <- nc$var[[coordinate]]$dim
    if (length(dims) != 1L) {
      stop("Buffer coordinates must be one-dimensional.", call. = FALSE)
    }
    dims[[1L]]$name
  }
  bounds_var <- nc$var[[name]]
  axes <- vapply(bounds_var$dim, function(x) x$name, character(1))
  coordinate_axis <- which(axes == dimension)
  raw <- ncdf4::ncvar_get(nc, name, collapse_degen = FALSE)

  if (length(dim(raw)) != 2L || length(coordinate_axis) != 1L) {
    stop("Unsupported bounds dimensions for `", coordinate, "`.",
         call. = FALSE)
  }
  if (coordinate_axis == 1L) raw <- t(raw)
  if (!identical(dim(raw), c(2L, as.integer(n))) ||
        !is.numeric(raw) || any(!is.finite(raw))) {
    stop("Bounds for `", coordinate,
         "` must contain two finite endpoints per coordinate.",
         call. = FALSE)
  }

  units <- .read_time_attribute(nc, coordinate, "units")
  bounds_units <- .read_time_attribute(nc, name, "units")
  if (!is.null(bounds_units) && !identical(bounds_units, units)) {
    stop("Bounds for `", coordinate, "` use different units. ",
         "Convert bounds and coordinate to the same units first.",
         call. = FALSE)
  }

  list(name = name, values = raw)
}


.buffer_axis_intervals <- function(centres, bounds, axis, coordinate) {
  if (!length(centres) || any(!is.finite(centres))) {
    stop("Buffer coordinates must contain finite values.", call. = FALSE)
  }
  steps <- diff(centres)
  if (length(steps) && !(all(steps > 0) || all(steps < 0))) {
    stop("Coordinate `", coordinate,
         "` must be strictly increasing or decreasing, without duplicates.",
         call. = FALSE)
  }

  clipped_at_surface <- rep(FALSE, length(centres))
  if (is.null(bounds)) {
    if (length(centres) < 2L) {
      stop("Cannot infer bounds from a single coordinate in `", coordinate,
           "`. Supply explicit bounds; one centre has no known extent.",
           call. = FALSE)
    }
    index <- order(centres)
    sorted <- centres[index]
    gaps <- diff(sorted)
    edges <- c(
      sorted[1L] - gaps[1L] / 2,
      head(sorted, -1L) + gaps / 2,
      tail(sorted, 1L) + tail(gaps, 1L) / 2
    )
    lower <- upper <- numeric(length(centres))
    lower[index] <- head(edges, -1L)
    upper[index] <- tail(edges, -1L)
    if (axis == "depth") {
      clipped_at_surface <- lower < 0
      lower <- pmax(lower, 0)
    }
    basis <- "inferred_midpoints"
  } else {
    lower <- pmin(bounds[1L, ], bounds[2L, ])
    upper <- pmax(bounds[1L, ], bounds[2L, ])
    basis <- "explicit_bounds"
  }

  if (any(!is.finite(lower)) || any(!is.finite(upper)) ||
        any(upper <= lower) || any(!is.finite(upper - lower)) ||
        any(centres < lower | centres > upper)) {
    stop("Invalid bounds for `", coordinate,
         "`: intervals need positive widths and must contain their centres.",
         call. = FALSE)
  }
  index <- order(lower)
  if (length(index) > 1L &&
        any(tail(lower[index], -1L) < head(cummax(upper[index]), -1L))) {
    stop("Bounds overlap in `", coordinate,
         "`. Resolve overlapping support before buffer aggregation.",
         call. = FALSE)
  }

  # Keep original record order so values and intervals retain their alignment.
  list(
    centres = centres, lower = lower, upper = upper,
    source_index = seq_along(centres), basis = basis,
    clipped_at_surface = clipped_at_surface
  )
}


.buffer_read_axis <- function(
  nc,
  coordinate,
  axis = c("depth", "time"),
  infer_bounds = FALSE
) {
  axis <- match.arg(axis)
  if (!is.logical(infer_bounds) || length(infer_bounds) != 1L ||
        is.na(infer_bounds)) {
    stop("`infer_bounds` must be TRUE or FALSE.", call. = FALSE)
  }
  raw <- as.vector(.get_dim_values(nc, coordinate))
  bounds <- .buffer_read_bounds(nc, coordinate, length(raw))
  units <- .read_time_attribute(nc, coordinate, "units")
  calendar <- NULL
  direction_basis <- NULL

  if (axis == "depth") {
    depth <- .read_depth_axis(nc, coordinate, raw)
    centres <- depth$depth_m
    direction_basis <- depth$direction_basis

    if (!is.null(bounds)) {
      positive <- .read_time_attribute(nc, coordinate, "positive")
      if (is.null(positive)) {
        positive <- if (all(raw >= 0)) "down" else "up"
      }
      bounds_positive <- .read_time_attribute(nc, bounds$name, "positive")
      if (!is.null(bounds_positive) &&
            !identical(tolower(trimws(bounds_positive)),
                       tolower(trimws(positive)))) {
        stop("Depth bounds and coordinate have conflicting directions.",
             call. = FALSE)
      }
      # Unique endpoints avoid treating shared edges as duplicate depths.
      endpoints <- unique(as.vector(bounds$values))
      converted <- .prepare_depth_axis(endpoints, units, positive)$depth_m
      bounds$values <- matrix(
        converted[match(as.vector(bounds$values), endpoints)], nrow = 2L
      )
    }
  } else {
    if (!is.null(.read_time_attribute(nc, coordinate, "climatology"))) {
      stop("Climatological coordinates require separate buffer handling. ",
           "Supply an ordinary dated time axis instead.", call. = FALSE)
    }
    centres <- as.numeric(.decode_netcdf_time(nc, coordinate))
    calendar <- .read_time_attribute(nc, coordinate, "calendar")
    if (is.null(calendar)) calendar <- "standard"
    if (!is.null(bounds)) {
      bounds_calendar <- .read_time_attribute(nc, bounds$name, "calendar")
      if (!is.null(bounds_calendar) &&
            !identical(tolower(bounds_calendar), tolower(calendar))) {
        stop("Time bounds and coordinate use different calendars.",
             call. = FALSE)
      }
      bounds$values <- matrix(
        as.numeric(.decode_netcdf_time(
          nc, coordinate, offsets = as.vector(bounds$values)
        )),
        nrow = 2L
      )
    }
  }

  if (is.null(bounds) && !infer_bounds) {
    stop("Coordinate `", coordinate, "` has no explicit bounds. ",
         "Use a dataset with bounds or explicitly set `infer_bounds = TRUE` ",
         "to estimate midpoint intervals.", call. = FALSE)
  }

  result <- .buffer_axis_intervals(
    centres,
    if (is.null(bounds)) NULL else bounds$values,
    axis,
    coordinate
  )
  c(result, list(
    coordinate = coordinate, axis = axis,
    units = if (axis == "depth") "m" else "UTC seconds since 1970-01-01",
    source_units = units,
    bounds_name = if (is.null(bounds)) NULL else bounds$name,
    calendar = calendar, direction_basis = direction_basis
  ))
}


# Summarise non-overlapping depth or time intervals.
# Time coordinates are supplied as UTC seconds; depth coordinates as metres.
.buffer_interval_mean <- function(
  values,
  lower,
  upper,
  requested,
  weighting = c("weighted", "equal"),
  min_valid_fraction = 0.8
) {
  weighting <- match.arg(weighting)

  if (!is.numeric(requested) || length(requested) != 2L ||
        any(!is.finite(requested)) ||
        requested[2L] <= requested[1L] ||
        !is.finite(diff(requested))) {
    stop(
      "`requested` must contain two finite, increasing interval limits.",
      call. = FALSE
    )
  }

  if (!is.numeric(values) || !is.numeric(lower) || !is.numeric(upper) ||
        length(values) != length(lower) ||
        length(values) != length(upper) ||
        any(!is.finite(lower)) || any(!is.finite(upper)) ||
        any(upper <= lower)) {
    stop(
      "Values and interval limits must have matching lengths, ",
      "with finite, increasing limits for every interval.",
      call. = FALSE
    )
  }

  if (!is.numeric(min_valid_fraction) ||
        length(min_valid_fraction) != 1L ||
        is.na(min_valid_fraction) ||
        !is.finite(min_valid_fraction) ||
        min_valid_fraction < 0 || min_valid_fraction > 1) {
    stop(
      "`min_valid_fraction` must be one number between 0 and 1.",
      call. = FALSE
    )
  }

  widths <- pmax(
    0,
    pmin(upper, requested[2L]) - pmax(lower, requested[1L])
  )
  selected <- which(widths > 0)

  # Overlap is checked only within the requested interval.
  if (length(selected) > 1L) {
    clipped_lower <- pmax(lower[selected], requested[1L])
    clipped_upper <- pmin(upper[selected], requested[2L])
    order <- order(clipped_lower)
    previous_end <- head(cummax(clipped_upper[order]), -1L)

    if (any(tail(clipped_lower[order], -1L) < previous_end)) {
      stop(
        "Selected intervals overlap. Resolve overlapping support ",
        "before aggregation; no source is chosen automatically.",
        call. = FALSE
      )
    }
  }

  valid <- selected[is.finite(values[selected])]
  requested_width <- diff(requested)
  coverage_fraction <- min(1, sum(widths[selected]) / requested_width)
  valid_fraction <- min(1, sum(widths[valid]) / requested_width)

  status <- if (!length(selected)) {
    "empty"
  } else if (!length(valid)) {
    "all_missing"
  } else if (valid_fraction < min_valid_fraction) {
    "insufficient_coverage"
  } else if (valid_fraction < 1) {
    "partial"
  } else {
    "complete"
  }

  value <- NA_real_
  if (status %in% c("complete", "partial")) {
    weights <- if (weighting == "weighted") {
      widths[valid]
    } else {
      rep(1, length(valid))
    }
    value <- sum(values[valid] * (weights / sum(weights)))
  }

  list(
    value = value,
    coverage_fraction = coverage_fraction,
    valid_fraction = valid_fraction,
    n_intersecting = length(selected),
    n_valid = length(valid),
    contributing_rows = valid,
    status = status
  )
}


# Select one physical layer for a time-buffer sample.
.buffer_point_value <- function(nc, var, dims, horizontal, time_index,
                                depth_axis, method, target_depth) {
  if (method == "2d") {
    value <- .extract_array(
      nc, var, dims, horizontal$lon_idx, horizontal$lat_idx,
      time_idx = time_index
    )
    return(list(value = as.numeric(value)[1L], depth_m = NA_real_,
                depth_index = if (is.null(dims$depth)) NA_integer_ else 1L))
  }

  if (method == "surface") {
    index <- depth_axis$surface_index
    value <- .extract_array(
      nc, var, dims, horizontal$lon_idx, horizontal$lat_idx,
      time_idx = time_index, depth_start = index, depth_count = 1L
    )
    return(list(value = as.numeric(value)[1L],
                depth_m = depth_axis$depth_m[index], depth_index = index))
  }

  profile <- as.numeric(.extract_array(
    nc, var, dims, horizontal$lon_idx, horizontal$lat_idx,
    time_idx = time_index, depth_start = 1L,
    depth_count = length(depth_axis$depth_m)
  ))
  valid <- which(is.finite(profile))
  if (!length(valid)) {
    return(list(value = NA_real_, depth_m = NA_real_,
                depth_index = NA_integer_))
  }

  ordered <- valid[order(depth_axis$depth_m[valid])]
  index <- if (method == "bottom") {
    tail(ordered, 1L)
  } else {
    ordered[.nearest_index(depth_axis$depth_m[ordered], target_depth)]
  }
  list(value = profile[index], depth_m = depth_axis$depth_m[index],
       depth_index = index)
}


.buffer_rejected_result <- function(status) {
  list(
    value = NA_real_, coverage_fraction = NA_real_, valid_fraction = NA_real_,
    n_intersecting = 0L, n_valid = 0L, contributing_rows = integer(),
    depth_indices = integer(), time_indices = integer(),
    selected_depth_m = numeric(), status = status
  )
}


# Read a single variable from one source, with depth and/or time windows.
# Multi-source merging and horizontal buffers are added in later blocks.
.buffer_extract_one_source <- function(
  data,
  source,
  var,
  options,
  lon_col = "lon",
  lat_col = "lat",
  depth_col = "depth",
  date_col = "date",
  id_col = NULL,
  lon_dim = NULL,
  lat_dim = NULL,
  depth_dim = NULL,
  time_dim = NULL,
  time_match = c("nearest", "exact"),
  time_resolution = c("auto", "datetime", "daily", "monthly", "yearly"),
  max_time_diff = NULL,
  spatial_match = c("within", "nearest"),
  max_distance_km = NULL
) {
  depth_buffer <- options$depth_radius_m > 0
  time_buffer <- options$time_before > 0 || options$time_after > 0
  if (options$horizontal_radius_km > 0 || !(depth_buffer || time_buffer)) {
    stop("This internal reader supports depth and time buffers. ",
         "Horizontal buffers are not implemented yet.",
         call. = FALSE)
  }
  time_match <- match.arg(time_match)
  time_resolution <- match.arg(time_resolution)
  spatial_match <- .validate_spatial_matching(spatial_match, max_distance_km)
  if (time_buffer && !is.null(max_time_diff)) {
    stop("`max_time_diff` applies to point time selection. ",
         "A time buffer uses `time_before` and `time_after` instead.",
         call. = FALSE)
  }
  if (!is.null(max_time_diff) &&
        (!is.numeric(max_time_diff) || length(max_time_diff) != 1L ||
         is.na(max_time_diff) || !is.finite(max_time_diff) ||
         max_time_diff < 0)) {
    stop("`max_time_diff` must be NULL or one finite, non-negative number.",
         call. = FALSE)
  }
  .check_required_cols(data, c(lon_col, lat_col, id_col))
  if (!is.numeric(data[[lon_col]]) || !is.numeric(data[[lat_col]])) {
    stop("Observation longitude and latitude columns must be numeric.",
         call. = FALSE)
  }

  opened <- .open_source(source)
  on.exit(.close_source(opened), add = TRUE)
  nc <- opened$nc
  dims <- .detect_dims(
    nc, var, lon_dim = lon_dim, lat_dim = lat_dim,
    depth_dim = depth_dim, time_dim = time_dim
  )
  if (depth_buffer && is.null(dims$depth)) {
    stop("The requested depth buffer requires a depth dimension in `",
         var, "`.", call. = FALSE)
  }
  if (time_buffer && is.null(dims$time)) {
    stop("The requested time buffer requires a time dimension in `",
         var, "`.", call. = FALSE)
  }
  if (time_buffer && options$method != "2d" && is.null(dims$depth)) {
    stop("This variable has no depth dimension. Use `method = '2d'`.",
         call. = FALSE)
  }
  if (time_buffer && options$method == "2d" && !is.null(dims$depth) &&
        length(.get_dim_values(nc, dims$depth)) > 1L) {
    stop("This variable has multiple depth layers. Select surface, bottom ",
         "or nearest instead of `method = '2d'`.", call. = FALSE)
  }
  if (!is.null(dims$time)) {
    .require_observation_time(data, date_col, var, nc$filename)
  }
  if (time_buffer && options$method == "nearest") {
    if (is.null(depth_col)) {
      stop("Nearest-depth selection requires `depth_col`.", call. = FALSE)
    }
    .check_required_cols(data, depth_col)
    target <- data[[depth_col]]
    if (!is.numeric(target) ||
          any(!is.na(target) & (!is.finite(target) | target < 0))) {
      stop("Observation depths must be finite, non-negative metres or NA.",
           call. = FALSE)
    }
  }

  windows <- .buffer_observation_windows(
    data, options, depth_col, date_col, id_col
  )
  prepared_time <- windows$prepared_time
  if (!is.null(dims$time) && is.null(prepared_time)) {
    prepared_time <- .prepare_time_input(
      data[[date_col]], time_resolution,
      input_name = paste0("`data$", date_col, "`"),
      ids = if (is.null(id_col)) NULL else data[[id_col]]
    )
  }

  lon <- .get_dim_values(nc, dims$lon)
  lat <- .get_dim_values(nc, dims$lat)
  lon_axis <- .prepare_spatial_axis(nc, dims$lon, lon, longitude = TRUE)
  lat_axis <- .prepare_spatial_axis(nc, dims$lat, lat)
  depth_support <- if (depth_buffer) {
    .buffer_read_axis(nc, dims$depth, "depth", options$infer_bounds)
  } else {
    NULL
  }
  time_support <- if (time_buffer) {
    .buffer_read_axis(nc, dims$time, "time", options$infer_bounds)
  } else {
    NULL
  }
  depth_axis <- if (time_buffer && !depth_buffer && options$method != "2d") {
    .read_depth_axis(nc, dims$depth, .get_dim_values(nc, dims$depth))
  } else {
    NULL
  }
  timestamps <- if (!is.null(dims$time) && !time_buffer) {
    .decode_netcdf_time(nc, dims$time)
  } else {
    NULL
  }
  if (!is.null(timestamps) && any(!is.na(prepared_time$values))) {
    .check_time_resolution(timestamps, prepared_time$resolution)
  }

  results <- vector("list", nrow(data))
  for (i in seq_len(nrow(data))) {
    horizontal <- .match_horizontal_cell(
      lon_axis, lat_axis, lon, lat,
      data[[lon_col]][i], data[[lat_col]][i],
      spatial_match, max_distance_km
    )
    result <- .buffer_rejected_result("spatial_rejected")
    if (horizontal$accepted) {
      if (!is.null(prepared_time) && is.na(prepared_time$values[i])) {
        result <- .buffer_rejected_result("missing_time")
      } else if ((depth_buffer || options$method == "nearest") &&
                   is.na(data[[depth_col]][i])) {
        result <- .buffer_rejected_result("missing_depth")
      } else if (depth_buffer && time_buffer) {
        requested_depth <- c(
          windows$windows$depth_lower_m[i],
          windows$windows$depth_upper_m[i]
        )
        requested_time <- as.numeric(c(
          windows$windows$time_lower_utc[i],
          windows$windows$time_upper_utc[i]
        ))

        depth_hit <- which(
          depth_support$upper > requested_depth[1L] &
            depth_support$lower < requested_depth[2L]
        )
        time_hit <- which(
          time_support$upper > requested_time[1L] &
            time_support$lower < requested_time[2L]
        )

        depth_widths <- pmin(
          depth_support$upper[depth_hit], requested_depth[2L]
        ) - pmax(depth_support$lower[depth_hit], requested_depth[1L])

        time_widths <- pmin(
          time_support$upper[time_hit], requested_time[2L]
        ) - pmax(time_support$lower[time_hit], requested_time[1L])

        values <- matrix(
          NA_real_, nrow = length(depth_hit), ncol = length(time_hit)
        )
        if (length(depth_hit) && length(time_hit)) {
          for (j in seq_along(time_hit)) {
            profile <- as.numeric(.extract_array(
              nc, var, dims, horizontal$lon_idx, horizontal$lat_idx,
              time_idx = time_hit[j],
              depth_start = min(depth_hit),
              depth_count = max(depth_hit) - min(depth_hit) + 1L
            ))
            values[, j] <- profile[depth_hit - min(depth_hit) + 1L]
          }
        }

        result <- .buffer_joint_mean(
          values, depth_widths, time_widths,
          diff(requested_depth), diff(requested_time),
          options$weighting, options$min_valid_fraction
        )
        result$depth_indices <-
          depth_hit[result$contributing_cells$depth]
        result$time_indices <-
          time_hit[result$contributing_cells$time]
        result$selected_depth_m <-
          depth_support$centres[result$depth_indices]
      } else if (depth_buffer) {
        time_index <- if (is.null(timestamps)) {
          NULL
        } else {
          .match_time_index(
            timestamps, prepared_time$values[i], time_match,
            max_time_diff, prepared_time$resolution
          )
        }
        if (!is.null(time_index) && is.na(time_index)) {
          result <- .buffer_rejected_result("time_rejected")
        } else {
          requested <- c(
            windows$windows$depth_lower_m[i], windows$windows$depth_upper_m[i]
          )
          hit <- which(depth_support$upper > requested[1L] &
                         depth_support$lower < requested[2L])
          values <- numeric()
          if (length(hit)) {
            # Read a contiguous subset, then retain only intersecting layers.
            profile <- as.numeric(.extract_array(
              nc, var, dims, horizontal$lon_idx, horizontal$lat_idx,
              time_idx = time_index, depth_start = min(hit),
              depth_count = max(hit) - min(hit) + 1L
            ))
            values <- profile[hit - min(hit) + 1L]
          }
          result <- .buffer_interval_mean(
            values, depth_support$lower[hit], depth_support$upper[hit],
            requested, options$weighting, options$min_valid_fraction
          )
          result$depth_indices <- hit[result$contributing_rows]
          result$time_indices <- if (result$n_valid) {
            if (is.null(time_index)) integer() else time_index
          } else {
            integer()
          }
          result$selected_depth_m <-
            depth_support$centres[result$depth_indices]
        }
      } else {
        requested <- as.numeric(c(
          windows$windows$time_lower_utc[i], windows$windows$time_upper_utc[i]
        ))
        hit <- which(time_support$upper > requested[1L] &
                       time_support$lower < requested[2L])
        samples <- lapply(hit, function(index) {
          .buffer_point_value(
            nc, var, dims, horizontal, index, depth_axis, options$method,
            if (options$method == "nearest") data[[depth_col]][i] else NA_real_
          )
        })
        values <- vapply(samples, function(x) x$value, numeric(1))
        depths <- vapply(samples, function(x) x$depth_m, numeric(1))
        result <- .buffer_interval_mean(
          values, time_support$lower[hit], time_support$upper[hit],
          requested, options$weighting, options$min_valid_fraction
        )
        result$time_indices <- hit[result$contributing_rows]
        depth_indices <- vapply(samples, function(x) x$depth_index, integer(1))
        result$depth_indices <- depth_indices[result$contributing_rows]
        result$selected_depth_m <- depths[result$contributing_rows]
      }
    }
    result$horizontal <- horizontal
    results[[i]] <- result
  }

  list(
    observations = data,
    values = vapply(results, function(x) x$value, numeric(1)),
    diagnostics = results,
    prepared_time = prepared_time,
    windows = windows$windows,
    metadata = list(
      variable = var, source_file = nc$filename,
      depth_support = depth_support, time_support = time_support
    )
  )
}


# Combine intersecting depth thicknesses and time durations.
# Widths must come from validated, non-overlapping axis bounds.
.buffer_joint_mean <- function(
  values,
  depth_widths,
  time_widths,
  requested_depth_width,
  requested_time_width,
  weighting = c("weighted", "equal"),
  min_valid_fraction = 0.8
) {
  weighting <- match.arg(weighting)

  for (widths in list(depth_widths, time_widths)) {
    if (!is.numeric(widths) || any(!is.finite(widths)) ||
          any(widths < 0)) {
      stop("Intersecting widths must be finite, non-negative numbers.",
           call. = FALSE)
    }
  }
  for (width in list(requested_depth_width, requested_time_width)) {
    if (!is.numeric(width) || length(width) != 1L ||
          !is.finite(width) || width <= 0) {
      stop("Requested windows must have finite, positive widths.",
           call. = FALSE)
    }
  }
  if (!is.matrix(values) || !is.numeric(values) ||
        !identical(dim(values), c(length(depth_widths), length(time_widths)))) {
    stop("The value matrix must match the depth and time intervals.",
         call. = FALSE)
  }
  if (!is.numeric(min_valid_fraction) ||
        length(min_valid_fraction) != 1L ||
        is.na(min_valid_fraction) || !is.finite(min_valid_fraction) ||
        min_valid_fraction < 0 || min_valid_fraction > 1) {
    stop("`min_valid_fraction` must be one number between 0 and 1.",
         call. = FALSE)
  }

  # Normalise each axis separately to avoid multiplying large physical widths.
  support <- outer(
    depth_widths / requested_depth_width,
    time_widths / requested_time_width
  )
  selected <- which(support > 0)
  usable <- support > 0 & is.finite(values)
  valid <- which(usable)

  coverage_fraction <- min(1, sum(support[selected]))
  valid_fraction <- min(1, sum(support[valid]))

  status <- if (!length(selected)) {
    "empty"
  } else if (!length(valid)) {
    "all_missing"
  } else if (valid_fraction < min_valid_fraction) {
    "insufficient_coverage"
  } else if (valid_fraction < 1) {
    "partial"
  } else {
    "complete"
  }

  value <- NA_real_
  if (status %in% c("complete", "partial")) {
    weights <- if (weighting == "weighted") {
      support[valid]
    } else {
      rep(1, length(valid))
    }
    value <- sum(values[valid] * (weights / sum(weights)))
  }

  cells <- which(usable, arr.ind = TRUE)
  list(
    value = value,
    coverage_fraction = coverage_fraction,
    valid_fraction = valid_fraction,
    n_intersecting = length(selected),
    n_valid = length(valid),
    contributing_rows = valid,
    contributing_cells = data.frame(
      depth = cells[, 1L],
      time = cells[, 2L]
    ),
    status = status
  )
}


# Build a distance-based buffer and retain an equal-area representation.
.buffer_horizontal_geometry <- function(lon, lat, horizontal_radius_km) {
  for (value in list(lon, lat, horizontal_radius_km)) {
    if (!is.numeric(value) || length(value) != 1L ||
          !is.finite(value)) {
      stop("Buffer coordinates and radius must be finite scalar numbers.",
           call. = FALSE)
    }
  }
  if (lat < -90 || lat > 90) {
    stop("Observation latitude must be between -90 and 90 degrees.",
         call. = FALSE)
  }
  if (horizontal_radius_km <= 0 || horizontal_radius_km >= 10000) {
    stop("This geometry reader requires a positive radius below 10000 km.",
         call. = FALSE)
  }

  lon <- .normalise_longitude(lon)
  distance_crs <- sprintf(
    paste0(
      "+proj=aeqd +lat_0=%.15g +lon_0=%.15g ",
      "+datum=WGS84 +units=m +no_defs"
    ),
    lat, lon
  )
  area_crs <- sprintf(
    paste0(
      "+proj=laea +lat_0=%.15g +lon_0=%.15g ",
      "+datum=WGS84 +units=m +no_defs"
    ),
    lat, lon
  )

  centre <- sf::st_sfc(
    sf::st_point(c(0, 0)), crs = distance_crs
  )
  circle <- sf::st_buffer(
    centre, dist = horizontal_radius_km * 1000, nQuadSegs = 180L
  )
  geometry <- sf::st_transform(circle, area_crs)
  area_m2 <- as.numeric(sf::st_area(geometry))

  if (!isTRUE(all(sf::st_is_valid(geometry))) ||
        any(sf::st_is_empty(geometry)) ||
        !is.finite(area_m2) || area_m2 <= 0) {
    stop("The requested horizontal buffer could not be constructed. ",
         "Check its coordinates and radius.", call. = FALSE)
  }

  list(
    geometry = geometry,
    area_m2 = area_m2,
    lon = lon,
    lat = lat,
    horizontal_radius_km = horizontal_radius_km,
    distance_crs = sf::st_crs(distance_crs),
    geometry_method = "WGS84 distance buffer, equal-area intersections",
    boundary_segments = 720L
  )
}


# Cells must already have valid boundaries in the buffer's equal-area CRS.
.buffer_horizontal_support <- function(buffer, cells) {
  if (!inherits(cells, c("sf", "sfc"))) {
    stop("Cells must be supplied as sf or sfc polygons.", call. = FALSE)
  }
  cells <- sf::st_geometry(cells)
  if (is.na(sf::st_crs(cells)) ||
        !isTRUE(sf::st_crs(cells) == sf::st_crs(buffer$geometry))) {
    stop("Cells and buffer must use the same projected equal-area CRS.",
         call. = FALSE)
  }
  if (!all(sf::st_geometry_type(cells) %in% c("POLYGON", "MULTIPOLYGON")) ||
        !isTRUE(all(sf::st_is_valid(cells))) ||
        any(sf::st_is_empty(cells))) {
    stop("Cell boundaries must be valid, non-empty polygons.",
         call. = FALSE)
  }

  empty <- data.frame(
    cell_index = integer(), area_m2 = numeric(), fraction = numeric()
  )
  clipped <- sf::st_intersection(cells, buffer$geometry)
  if (!length(clipped)) return(empty)

  areas <- as.numeric(sf::st_area(clipped))
  if (any(!is.finite(areas))) {
    stop("Cell intersections produced non-finite areas.", call. = FALSE)
  }
  keep <- which(areas > 0)
  if (!length(keep)) return(empty)

  index <- attr(clipped, "idx")[keep, 1L]
  clipped <- clipped[keep]

  # Shared edges are allowed; overlapping interiors would double-count area.
  if (length(clipped) > 1L) {
    overlaps <- sf::st_relate(
      clipped, pattern = "2********", sparse = TRUE
    )
    if (any(lengths(overlaps) > 1L)) {
      stop("Selected cell interiors overlap within the buffer. ",
           "Resolve overlapping support before aggregation.",
           call. = FALSE)
    }
  }

  totals <- rowsum(areas[keep], index, reorder = TRUE)
  area_m2 <- as.numeric(totals[, 1L])
  data.frame(
    cell_index = as.integer(rownames(totals)),
    area_m2 = area_m2,
    fraction = area_m2 / buffer$area_m2
  )
}


# Validate the public request without changing the earlier interval helpers.
.buffer_validate_request <- function(horizontal_radius_km, depth_radius_m,
                                     lon_cell_radius, lat_cell_radius,
                                     depth_layer_radius, ...) {
  counts <- list(lon_cell_radius = lon_cell_radius,
                 lat_cell_radius = lat_cell_radius,
                 depth_layer_radius = depth_layer_radius)
  for (name in names(counts)) {
    value <- counts[[name]]
    if (!is.null(value) &&
          (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
           value < 0 || value != floor(value))) {
      stop("`", name, "` must be NULL or one non-negative whole number.",
           call. = FALSE)
    }
  }
  horizontal_index <- !is.null(lon_cell_radius) || !is.null(lat_cell_radius)
  depth_index <- !is.null(depth_layer_radius)
  for (value in list(horizontal_radius_km, depth_radius_m)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0) {
      stop("Metric radii must be finite, non-negative scalar numbers.", call. = FALSE)
    }
  }
  if (horizontal_index && horizontal_radius_km != 0) {
    stop("Choose horizontal kilometres or cell radii, not both.", call. = FALSE)
  }
  if (depth_index && depth_radius_m != 0) {
    stop("Choose depth metres or a layer radius, not both.", call. = FALSE)
  }
  options <- .validate_buffer_options(
    horizontal_radius_km = if (horizontal_index) 1 else horizontal_radius_km,
    depth_radius_m = if (depth_index) 1 else depth_radius_m, ...
  )
  if (!any(c(horizontal_radius_km, depth_radius_m, options$time_before,
             options$time_after, unlist(counts)) > 0)) {
    stop("At least one buffer window must be positive. Use extract_to_point().",
         call. = FALSE)
  }
  options$horizontal_radius_km <- horizontal_radius_km
  options$depth_radius_m <- depth_radius_m
  c(options, counts, list(horizontal_index = horizontal_index,
                         depth_index = depth_index))
}


.buffer_index_neighbours <- function(order, seed, radius, periodic = FALSE) {
  position <- match(seed, order)
  if (2 * radius + 1 > length(order)) return(NULL)
  indices <- seq.int(position - radius, position + radius)
  if (periodic) {
    indices <- (indices - 1) %% length(order) + 1L
  } else if (any(indices < 1 | indices > length(order))) {
    return(NULL)
  }
  order[indices]
}


.buffer_spatial_axis <- function(nc, coordinate, longitude, infer_bounds) {
  axis <- .prepare_spatial_axis(
    nc, coordinate, .get_dim_values(nc, coordinate), longitude
  )
  if (axis$basis != "explicit_bounds" && !infer_bounds) {
    stop("Coordinate `", coordinate, "` needs cell bounds for aggregation. ",
         "Supply bounds or explicitly set `infer_bounds = TRUE`.",
         call. = FALSE)
  }
  if (any(axis$upper <= axis$lower)) {
    stop("A singleton spatial coordinate needs explicit cell bounds.",
         call. = FALSE)
  }
  .buffer_axis_intervals(
    axis$centres, rbind(axis$lower, axis$upper), "time", coordinate
  )
  if (longitude &&
        (anyDuplicated(.normalise_longitude(axis$centres)) ||
         sum(axis$upper - axis$lower) > 360 + 1e-8)) {
    stop("Longitude cells contain duplicate or overlapping periodic support.",
         call. = FALSE)
  }
  axis
}


# Subdivide parallel and meridian edges before projecting them.
# Longitude strips avoid folded rings at the antimeridian and poles.
.buffer_project_rectangle <- function(xmin, xmax, ymin, ymax, crs) {
  breaks <- seq(xmin, xmax, length.out = max(2L, ceiling((xmax - xmin) / 90) + 1L))
  pieces <- lapply(seq_len(length(breaks) - 1L), function(i) {
    x <- seq(breaks[i], breaks[i + 1L],
             length.out = max(2L, ceiling((breaks[i + 1L] - breaks[i]) / 0.01) + 1L))
    y <- seq(ymin, ymax, length.out = max(2L, ceiling((ymax - ymin) / 0.01) + 1L))
    ring <- rbind(cbind(x, ymin), cbind(tail(x, 1L), y),
                  cbind(rev(x), ymax), cbind(x[1L], rev(y)))
    ring[, 1L] <- .normalise_longitude(ring[, 1L])
    points <- sf::st_sfc(lapply(seq_len(nrow(ring)), function(j) {
      sf::st_point(ring[j, ])
    }), crs = 4326)
    projected <- sf::st_coordinates(sf::st_transform(points, crs))[, 1:2, drop = FALSE]
    if (any(!is.finite(projected))) {
      stop("Cell geometry cannot be projected for this buffer. ",
           "Use a smaller window or a supported regional grid.", call. = FALSE)
    }
    keep <- c(TRUE, rowSums(abs(projected[-1L, , drop = FALSE] -
                               projected[-nrow(projected), , drop = FALSE])) > 1e-8)
    projected <- projected[keep, , drop = FALSE]
    if (!all(projected[1L, ] == projected[nrow(projected), ])) {
      projected <- rbind(projected, projected[1L, ])
    }
    sf::st_polygon(list(projected))
  })
  geometry <- sf::st_union(sf::st_sfc(pieces, crs = crs))
  if (!isTRUE(all(sf::st_is_valid(geometry)))) {
    stop("Projected cell boundaries are invalid. Check the grid bounds.",
         call. = FALSE)
  }
  geometry[[1L]]
}


.buffer_cell_polygons <- function(meta, buffer, lon_indices = NULL,
                                  lat_indices = NULL, clip = TRUE) {
  x <- meta$lon_axis
  y <- meta$lat_axis
  if (is.null(lon_indices)) lon_indices <- seq_along(x$centres)
  if (is.null(lat_indices)) lat_indices <- seq_along(y$centres)
  # A conservative spherical cap bounds all WGS84 distances used here.
  angle <- buffer$horizontal_radius_km * 1000 / 6300000
  latitude <- if (clip) {
    pmax(-90, pmin(90, buffer$lat + c(-1, 1) * angle * 180 / pi))
  } else {
    c(-90, 90)
  }
  longitude_width <- if (!clip || abs(buffer$lat) * pi / 180 + angle >= pi / 2) {
    180
  } else {
    asin(sin(angle) / cos(buffer$lat * pi / 180)) * 180 / pi
  }
  longitude <- buffer$lon + c(-1, 1) * longitude_width
  lat_indices <- lat_indices[y$upper[lat_indices] > latitude[1L] &
                              y$lower[lat_indices] < latitude[2L]]
  mapping <- list()
  polygons <- list()
  for (i in lon_indices) {
    shift <- 360 * round((buffer$lon - x$centres[i]) / 360)
    intervals <- lapply(c(-360, 0, 360), function(extra) {
      c(max(x$lower[i] + shift + extra, longitude[1L]),
        min(x$upper[i] + shift + extra, longitude[2L]))
    })
    intervals <- Filter(function(z) z[2L] > z[1L], intervals)
    if (!length(intervals)) next
    for (j in lat_indices) {
      pieces <- lapply(intervals, function(z) {
        .buffer_project_rectangle(
          z[1L], z[2L], max(y$lower[j], latitude[1L]),
          min(y$upper[j], latitude[2L]), sf::st_crs(buffer$geometry)
        )
      })
      polygons[[length(polygons) + 1L]] <-
        sf::st_union(sf::st_sfc(pieces, crs = sf::st_crs(buffer$geometry)))[[1L]]
      mapping[[length(mapping) + 1L]] <- c(i, j)
    }
  }
  indices <- if (length(mapping)) do.call(rbind, mapping) else matrix(integer(), 0, 2)
  list(
    cells = sf::st_sfc(polygons, crs = sf::st_crs(buffer$geometry)),
    lon = indices[, 1L], lat = indices[, 2L]
  )
}


.buffer_detect_variable <- function(nc) {
  coordinates <- names(Filter(
    function(x) isTRUE(x$dimvarid$id >= 0), nc$dim
  ))
  metadata <- unlist(lapply(unique(c(coordinates, names(nc$var))), function(name) {
    unlist(lapply(c("bounds", "climatology", "coordinates"), function(attribute) {
      value <- .read_time_attribute(nc, name, attribute)
      if (is.character(value)) unlist(strsplit(value, "[[:space:]]+")) else NULL
    }))
  }))
  variables <- setdiff(names(nc$var), c(names(nc$dim), metadata))
  if (length(variables) != 1L) {
    stop("Automatic buffer variable detection requires one data variable. ",
         "Specify `var` explicitly; candidates: ", toString(variables), ".",
         call. = FALSE)
  }
  variables
}


.buffer_source_metadata <- function(nc, var, source_index, options,
                                    method, dimension_args) {
  dims <- do.call(.detect_dims, c(list(nc = nc, var = var), dimension_args))
  spatial <- options$horizontal_radius_km > 0 || options$horizontal_index
  vertical <- options$depth_radius_m > 0 || options$depth_index
  temporal <- options$time_before > 0 || options$time_after > 0
  if (vertical && is.null(dims$depth)) {
    stop("Variable `", var, "` needs a depth axis for the requested buffer.",
         call. = FALSE)
  }
  if (temporal && is.null(dims$time)) {
    stop("Variable `", var, "` needs a time axis for the requested buffer.",
         call. = FALSE)
  }
  selected_method <- if (vertical) {
    "interval"
  } else if (is.null(method)) {
    if (is.null(dims$depth)) "2d" else "nearest"
  } else {
    method
  }
  if (selected_method != "2d" && is.null(dims$depth)) {
    stop("Variable `", var, "` has no depth axis; use method = '2d'.",
         call. = FALSE)
  }
  raw_depth <- .get_dim_values(nc, dims$depth)
  if (selected_method == "2d" && length(raw_depth) > 1L) {
    stop("Variable `", var, "` has multiple layers; choose a 3D method.",
         call. = FALSE)
  }
  lon <- .get_dim_values(nc, dims$lon)
  lat <- .get_dim_values(nc, dims$lat)
  lon_axis <- if (spatial) {
    .buffer_spatial_axis(nc, dims$lon, TRUE, options$infer_bounds)
  } else {
    .prepare_spatial_axis(nc, dims$lon, lon, TRUE)
  }
  lat_axis <- if (spatial) {
    .buffer_spatial_axis(nc, dims$lat, FALSE, options$infer_bounds)
  } else {
    .prepare_spatial_axis(nc, dims$lat, lat)
  }
  if (!is.null(dims$time) &&
        !is.null(.read_time_attribute(nc, dims$time, "climatology"))) {
    stop("Climatological coordinates require separate buffer handling.",
         call. = FALSE)
  }
  list(
    nc = nc, var = var, source_index = source_index, dims = dims,
    source_file = nc$filename,
    variable_units = .read_time_attribute(nc, var, "units"),
    lon = lon, lat = lat, lon_axis = lon_axis, lat_axis = lat_axis,
    method = selected_method,
    depth_axis = if (!is.null(dims$depth) && selected_method != "2d") {
      .read_depth_axis(nc, dims$depth, raw_depth)
    } else NULL,
    depth_support = if (vertical) {
      .buffer_read_axis(nc, dims$depth, "depth", options$infer_bounds)
    } else NULL,
    time_values = as.numeric(.decode_netcdf_time(nc, dims$time)),
    time_support = if (temporal) {
      .buffer_read_axis(nc, dims$time, "time", options$infer_bounds)
    } else NULL
  )
}


.buffer_space_query <- function(meta, lon, lat, options, buffer,
                                 spatial_match, max_distance_km) {
  spatial <- options$horizontal_radius_km > 0 || options$horizontal_index
  if (!spatial) {
    match <- .match_horizontal_cell(
      meta$lon_axis, meta$lat_axis, meta$lon, meta$lat, lon, lat,
      spatial_match, max_distance_km
    )
    if (!match$accepted) return(list(status = "spatial_rejected"))
    return(list(status = "ready", cells = data.frame(
      lon = match$lon_idx, lat = match$lat_idx, fraction = 1,
      distance_km = match$distance_km, key = "point"
    ), geometry = NULL, buffer = NULL))
  }
  if (options$horizontal_index) {
    match <- .match_horizontal_cell(
      meta$lon_axis, meta$lat_axis, meta$lon, meta$lat, lon, lat,
      spatial_match, max_distance_km
    )
    if (!match$accepted) return(list(status = "spatial_rejected"))
    x_radius <- if (is.null(options$lon_cell_radius)) 0 else options$lon_cell_radius
    y_radius <- if (is.null(options$lat_cell_radius)) 0 else options$lat_cell_radius
    periodic <- abs(sum(meta$lon_axis$upper - meta$lon_axis$lower) - 360) < 1e-8
    x <- .buffer_index_neighbours(
      order(meta$lon_axis$centres), match$lon_idx, x_radius, periodic
    )
    y <- .buffer_index_neighbours(order(meta$lat), match$lat_idx, y_radius)
    if (is.null(x) || is.null(y)) return(list(status = "incomplete_index_window"))
    grid <- .buffer_cell_polygons(meta, buffer, x, y, clip = FALSE)
    buffer$geometry <- sf::st_union(grid$cells)
    buffer$area_m2 <- as.numeric(sf::st_area(buffer$geometry))
    buffer$geometry_method <- "Whole-cell index neighbourhood"
  } else {
    grid <- .buffer_cell_polygons(meta, buffer)
  }
  support <- .buffer_horizontal_support(buffer, grid$cells)
  cells <- data.frame(
    lon = grid$lon[support$cell_index], lat = grid$lat[support$cell_index],
    fraction = support$fraction, distance_km = rep(0, nrow(support))
  )
  cells$key <- vapply(seq_len(nrow(cells)), function(i) {
    x <- cells$lon[i]
    y <- cells$lat[i]
    paste(format(c(.normalise_longitude(meta$lon_axis$lower[x]),
                   meta$lon_axis$upper[x] - meta$lon_axis$lower[x],
                   meta$lat_axis$lower[y], meta$lat_axis$upper[y]),
                 digits = 17, scientific = TRUE), collapse = ":")
  }, character(1))
  geometry <- if (nrow(cells)) {
    clipped <- sf::st_intersection(grid$cells[support$cell_index], buffer$geometry)
    clipped[match(seq_len(nrow(cells)), attr(clipped, "idx")[, 1L])]
  } else {
    sf::st_sfc(crs = sf::st_crs(buffer$geometry))
  }
  list(status = "ready", cells = cells, geometry = geometry, buffer = buffer)
}


.buffer_depth_query <- function(meta, target, options) {
  active <- options$depth_radius_m > 0 || options$depth_index
  if (!active) {
    index <- if (is.null(meta$dims$depth)) NA_integer_ else {
      if (meta$method == "2d") 1L else seq_along(meta$depth_axis$depth_m)
    }
    return(list(status = "ready", index = index, fraction = rep(1, length(index)),
                lower = rep(NA_real_, length(index)),
                upper = rep(NA_real_, length(index))))
  }
  support <- meta$depth_support
  if (options$depth_index) {
    sorted <- order(support$centres)
    seed <- sorted[.nearest_index(support$centres[sorted], target)]
    index <- .buffer_index_neighbours(sorted, seed, options$depth_layer_radius)
    if (is.null(index)) return(list(status = "incomplete_index_window"))
    lower <- support$lower[index]
    upper <- support$upper[index]
    span <- sum(upper - lower)
  } else {
    request <- target + c(-1, 1) * options$depth_radius_m
    index <- which(support$upper > request[1L] & support$lower < request[2L])
    lower <- pmax(support$lower[index], request[1L])
    upper <- pmin(support$upper[index], request[2L])
    span <- diff(request)
  }
  if (!is.finite(span) || span <= 0) {
    stop("The requested depth buffer has no finite positive width.", call. = FALSE)
  }
  list(status = "ready", index = index, lower = lower, upper = upper,
       fraction = (upper - lower) / span)
}


.buffer_time_query <- function(meta, window, selected, temporal) {
  if (temporal) {
    support <- meta$time_support
    span <- diff(window)
    if (!is.finite(span) || span <= 0) {
      stop("The requested time buffer has no finite positive width.", call. = FALSE)
    }
    index <- which(support$upper > window[1L] & support$lower < window[2L])
    lower <- pmax(support$lower[index], window[1L])
    upper <- pmin(support$upper[index], window[2L])
    return(list(index = index, lower = lower, upper = upper,
                fraction = (upper - lower) / span))
  }
  if (is.null(meta$dims$time)) {
    return(list(index = NA_integer_, lower = NA_real_, upper = NA_real_, fraction = 1))
  }
  index <- which(meta$time_values == selected)
  list(index = index, lower = rep(selected, length(index)),
       upper = rep(selected, length(index)), fraction = rep(1, length(index)))
}


.buffer_empty_samples <- function() {
  data.frame(
    source_index = integer(), lon_index = integer(), lat_index = integer(),
    depth_index = integer(), time_index = integer(), value = numeric(),
    area_fraction = numeric(), depth_fraction = numeric(), time_fraction = numeric(),
    depth_m = numeric(), timestamp = numeric(), depth_lower = numeric(),
    depth_upper = numeric(), time_lower = numeric(), time_upper = numeric(),
    distance_km = numeric(), cell_key = character(), geometry = I(list())
  )
}


# Retain individual values and support fractions until the final statistic.
.buffer_collect_source_samples <- function(meta, space, depth, time) {
  if (!nrow(space$cells) || !length(depth$index) || !length(time$index)) {
    return(.buffer_empty_samples())
  }
  pieces <- list()
  for (i in seq_len(nrow(space$cells))) {
    cell <- space$cells[i, ]
    for (j in seq_along(time$index)) {
      z <- depth$index
      t <- time$index[j]
      start <- if (is.null(meta$dims$depth)) NULL else min(z)
      count <- if (is.null(start)) NULL else max(z) - start + 1L
      raw <- as.numeric(.extract_array(
        meta$nc, meta$var, meta$dims, cell$lon, cell$lat,
        time_idx = if (is.na(t)) NULL else t,
        depth_start = start, depth_count = count
      ))
      values <- if (is.null(start)) raw[1L] else raw[z - start + 1L]
      physical_depth <- if (is.null(meta$depth_axis)) {
        rep(NA_real_, length(z))
      } else meta$depth_axis$depth_m[z]
      geometry <- if (is.null(space$geometry)) NULL else space$geometry[[i]]
      pieces[[length(pieces) + 1L]] <- data.frame(
        source_index = meta$source_index,
        lon_index = cell$lon, lat_index = cell$lat,
        depth_index = z, time_index = t, value = values,
        area_fraction = cell$fraction, depth_fraction = depth$fraction,
        time_fraction = time$fraction[j], depth_m = physical_depth,
        timestamp = if (is.na(t)) NA_real_ else meta$time_values[t],
        depth_lower = depth$lower, depth_upper = depth$upper,
        time_lower = time$lower[j], time_upper = time$upper[j],
        distance_km = cell$distance_km, cell_key = cell$key,
        geometry = I(rep(list(geometry), length(z)))
      )
    }
  }
  do.call(rbind, pieces)
}


.buffer_sample_keys <- function(samples, depth = FALSE) {
  parts <- list(samples$time_lower, samples$time_upper)
  if (depth) parts <- c(parts, list(samples$depth_lower, samples$depth_upper))
  do.call(paste, c(lapply(parts, function(x) {
    format(x, digits = 17, scientific = TRUE, trim = TRUE)
  }), sep = ":"))
}


.buffer_select_point_samples <- function(samples, method, target_depth,
                                        spatial, vertical) {
  if (!nrow(samples)) return(samples)
  if (!spatial) {
    groups <- split(seq_len(nrow(samples)), .buffer_sample_keys(samples, vertical))
    keep <- unlist(lapply(groups, function(rows) {
      rows[samples$distance_km[rows] == min(samples$distance_km[rows])]
    }), use.names = FALSE)
    samples <- samples[sort(keep), , drop = FALSE]
  }
  if (!vertical && method != "2d") {
    groups <- split(seq_len(nrow(samples)), paste(
      samples$cell_key, .buffer_sample_keys(samples), sep = ":"
    ))
    keep <- unlist(lapply(groups, function(rows) {
      eligible <- rows[is.finite(samples$value[rows])]
      if (method == "surface" || !length(eligible)) eligible <- rows
      selected <- if (method == "surface") {
        min(samples$depth_m[eligible])
      } else if (method == "bottom") {
        max(samples$depth_m[eligible])
      } else {
        ordered <- eligible[order(samples$depth_m[eligible])]
        samples$depth_m[ordered[.nearest_index(samples$depth_m[ordered], target_depth)]]
      }
      # Keep duplicate support even if one copy is NA: overlap remains an error.
      rows[samples$depth_m[rows] == selected]
    }), use.names = FALSE)
    samples <- samples[sort(keep), , drop = FALSE]
  }
  samples
}


.buffer_record_overlap <- function(a, b, vertical, temporal) {
  for (i in seq.int(1L, nrow(a), by = 128L)) {
    x <- a[i:min(nrow(a), i + 127L), , drop = FALSE]
    for (j in seq.int(1L, nrow(b), by = 128L)) {
      y <- b[j:min(nrow(b), j + 127L), , drop = FALSE]
      hit <- matrix(TRUE, nrow(x), nrow(y))
      if (vertical) {
        hit <- hit & outer(x$depth_lower, y$depth_upper, "<") &
          outer(x$depth_upper, y$depth_lower, ">")
      }
      if (temporal) {
        hit <- hit & outer(x$time_lower, y$time_upper, "<") &
          outer(x$time_upper, y$time_lower, ">")
      }
      if (any(hit)) return(TRUE)
    }
  }
  FALSE
}


.buffer_check_source_overlap <- function(samples, metadata, spatial,
                                         vertical, temporal) {
  sources <- unique(samples$source_index)
  if (length(sources) < 2L) return(invisible(NULL))
  pairs <- utils::combn(sources, 2L)
  for (k in seq_len(ncol(pairs))) {
    a <- samples[samples$source_index == pairs[1L, k], , drop = FALSE]
    b <- samples[samples$source_index == pairs[2L, k], , drop = FALSE]
    a_cells <- !duplicated(paste(a$lon_index, a$lat_index))
    b_cells <- !duplicated(paste(b$lon_index, b$lat_index))
    a_ids <- paste(a$lon_index, a$lat_index)
    b_ids <- paste(b$lon_index, b$lat_index)
    a_unique <- which(a_cells)
    b_unique <- which(b_cells)
    relations <- if (spatial) {
      # sfg objects do not retain a CRS; use the shared projected query CRS.
      crs <- attr(samples, "buffer_crs")
      sf::st_relate(
        sf::st_sfc(a$geometry[a_unique], crs = crs),
        sf::st_sfc(b$geometry[b_unique], crs = crs),
        pattern = "2********", sparse = TRUE
      )
    } else {
      rep(list(seq_along(b_unique)), length(a_unique))
    }
    for (i in seq_along(relations)) {
      for (j in relations[[i]]) {
        x <- a[a_ids == a_ids[a_unique[i]], , drop = FALSE]
        y <- b[b_ids == b_ids[b_unique[j]], , drop = FALSE]
        if (.buffer_record_overlap(x, y, vertical, temporal)) {
          files <- vapply(pairs[, k], function(index) {
            metadata[[which(vapply(metadata, function(m) m$source_index == index,
                                  logical(1)))[1L]]]$source_file
          }, character(1))
          stop("Source files overlap within the requested buffer: ",
               paste(basename(files), collapse = ", "), ". ",
               "Choose complementary sources; overlapping support is not averaged.",
               call. = FALSE)
        }
      }
    }
  }
  invisible(NULL)
}


.buffer_statistic <- function(values, weights, statistic) {
  if (statistic == "min") return(min(values))
  if (statistic == "max") return(max(values))

  weights <- weights / sum(weights)
  average <- sum(values * weights)
  if (statistic == "mean") return(average)
  if (length(values) < 2L) return(NA_real_)

  correction <- sum(weights * (1 - weights))
  if (correction <= 0) return(NA_real_)
  deviation <- sqrt(sum(weights * (values - average)^2) / correction)
  if (statistic == "sd") return(deviation)

  effective_n <- 1 / sum(weights^2)
  deviation / sqrt(effective_n)
}


.buffer_summarise_samples <- function(samples, weighting, min_valid_fraction,
                                      statistic = "mean") {
  support <- with(samples, area_fraction * depth_fraction * time_fraction)
  selected <- support > 0
  valid <- selected & is.finite(samples$value)
  coverage <- min(1, sum(support[selected]))
  fraction <- min(1, sum(support[valid]))
  # Geometric overlays and products can differ from unity by machine precision.
  if (abs(coverage - 1) < 1e-10) coverage <- 1
  if (abs(fraction - 1) < 1e-10) fraction <- 1
  meets_threshold <- fraction >= min_valid_fraction ||
    abs(fraction - min_valid_fraction) < 1e-10
  status <- if (!any(selected)) "empty" else if (!any(valid)) {
    "all_missing"
  } else if (!meets_threshold) "insufficient_coverage" else {
    if (fraction < 1) "partial" else "complete"
  }
  value <- NA_real_
  effective_n <- NA_real_
  if (any(valid)) {
    weights <- if (weighting == "weighted") support[valid] else rep(1, sum(valid))
    effective_n <- 1 / sum((weights / sum(weights))^2)
  }
  if (any(valid) && meets_threshold) {
    value <- .buffer_statistic(samples$value[valid], weights, statistic)
    if (is.na(value) && statistic %in% c("sd", "se")) {
      status <- "insufficient_variation_samples"
    }
  }
  list(value = value, coverage_fraction = coverage, valid_fraction = fraction,
       n_intersecting = sum(selected), n_valid = sum(valid), status = status,
       contributing_rows = which(valid), effective_n = effective_n)
}


.buffer_check_index_grids <- function(metadata, options) {
  signature <- function(meta) {
    axis <- function(x, longitude = FALSE) {
      centres <- if (longitude) .normalise_longitude(x$centres) else x$centres
      order <- order(centres)
      lower <- if (longitude) .normalise_longitude(x$lower) else x$lower
      upper <- if (longitude) x$upper - x$lower else x$upper
      unname(cbind(centres, lower, upper)[order, , drop = FALSE])
    }
    result <- list()
    if (options$horizontal_index) {
      result$lon <- axis(meta$lon_axis, TRUE)
      result$lat <- axis(meta$lat_axis)
    }
    if (options$depth_index) result$depth <- axis(meta$depth_support)
    result
  }
  if (length(metadata) > 1L && (options$horizontal_index || options$depth_index)) {
    reference <- signature(metadata[[1L]])
    for (meta in metadata[-1L]) {
      if (!isTRUE(all.equal(reference, signature(meta), tolerance = 1e-12))) {
        stop("Index buffers require matching grid/layer coordinates and bounds ",
             "across sources for the same variable. Use metric buffers instead.",
             call. = FALSE)
      }
    }
  }
  invisible(NULL)
}


.buffer_output_names <- function(variables, output_col, input_names, match_info) {
  outputs <- if (is.null(output_col)) {
    variables
  } else if (length(variables) == 1L && length(output_col) == 1L) {
    unname(output_col)
  } else {
    if (is.null(names(output_col)) || anyDuplicated(names(output_col)) ||
          !setequal(names(output_col), variables)) {
      stop("For multiple variables, `output_col` must name each variable once.",
           call. = FALSE)
    }
    unname(output_col[variables])
  }
  if (!is.character(outputs) || anyNA(outputs) || any(!nzchar(outputs)) ||
        anyDuplicated(outputs)) {
    stop("Output column names must be distinct, non-empty strings.", call. = FALSE)
  }
  generated <- c(outputs, if (match_info) paste0(outputs, "_match_info"))
  if (anyDuplicated(generated) || any(generated %in% input_names)) {
    stop("Output or diagnostic names would overwrite existing columns. ",
         "Choose different `output_col` names.", call. = FALSE)
  }
  stats::setNames(outputs, variables)
}


.buffer_diagnostic <- function(result, samples, metadata, options, requested) {
  if (is.null(result$effective_n)) result$effective_n <- NA_real_
  valid <- result$contributing_rows
  used <- samples[valid, , drop = FALSE]
  range_or_null <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) range(x) else NULL
  }
  indices <- unique(samples$source_index)
  files <- vapply(indices, function(index) {
    metadata[[which(vapply(metadata, function(m) m$source_index == index,
                          logical(1)))[1L]]]$source_file
  }, character(1))
  bounds <- lapply(metadata, function(m) list(
    file = m$source_file,
    longitude = m$lon_axis$basis, latitude = m$lat_axis$basis,
    depth = if (is.null(m$depth_support)) NULL else m$depth_support$basis,
    time = if (is.null(m$time_support)) NULL else m$time_support$basis
  ))
  result$value <- NULL
  result$contributing_rows <- NULL
  c(result, list(
    requested = requested, statistic = options$statistic,
    weighting = options$weighting,
    min_valid_fraction = options$min_valid_fraction,
    variable = metadata[[1L]]$var,
    sources_intersecting = files,
    sources_used = vapply(unique(used$source_index), function(index) {
      metadata[[which(vapply(metadata, function(m) m$source_index == index,
                            logical(1)))[1L]]]$source_file
    }, character(1)),
    bounds = bounds,
    depth_range_m = range_or_null(used$depth_m),
    depth_support_range_m = range_or_null(c(used$depth_lower, used$depth_upper)),
    time_range_utc = if (any(is.finite(used$timestamp))) {
      as.POSIXct(range_or_null(used$timestamp), origin = "1970-01-01", tz = "UTC")
    } else NULL,
    time_support_range_utc = if (any(is.finite(used$time_lower))) {
      as.POSIXct(range_or_null(c(used$time_lower, used$time_upper)),
                 origin = "1970-01-01", tz = "UTC")
    } else NULL,
    geometry_method = if (options$horizontal_index) {
      "Whole-cell neighbourhood; cell edges at most 0.01 degrees"
    } else if (options$horizontal_radius_km > 0) {
      "WGS84; 720-segment circle; cell edges at most 0.01 degrees"
    } else "Point-selected horizontal cell"
  ))
}


.show_buffer_report <- function(diagnostics, data, variable, output, id_col) {
  affected <- which(vapply(diagnostics, function(x) x$status != "complete", logical(1)))
  if (!length(affected)) return(invisible(NULL))
  reasons <- c(
    empty = "No source covers the requested buffer.",
    all_missing = "The buffer is covered, but all selected values are missing.",
    insufficient_coverage = "Valid coverage is below the requested threshold; result: NA.",
    insufficient_variation_samples = paste(
      "SD or SE requires at least two valid samples with usable weights; result: NA."
    ),
    partial = "The result uses partial coverage and meets the requested threshold.",
    missing_time = "The observation date is NA; result: NA.",
    missing_depth = "The observation depth is NA; result: NA.",
    invalid_coordinates = "Observation coordinates are missing or invalid; result: NA.",
    incomplete_index_window = paste(
      "The index window extends beyond the grid; its full support is unknown.",
      "Result: NA. Reduce the window or use wider grid coverage."
    ),
    time_rejected = "No timestamp satisfies the point-time selection; result: NA."
  )
  examples <- vapply(head(affected, 5L), function(i) {
    info <- diagnostics[[i]]
    identifier <- if (is.null(id_col)) "" else paste0(" (id: ", data[[id_col]][i], ")")
    coverage <- if (is.finite(info$coverage_fraction)) {
      paste0(" Coverage: ", .selection_number(100 * info$coverage_fraction),
             "%; valid: ", .selection_number(100 * info$valid_fraction),
             "%; required: ", .selection_number(100 * info$min_valid_fraction), "%. ")
    } else " "
    paste0("Row ", i, identifier, ": ", reasons[[info$status]], coverage)
  }, character(1))
  message(structure(list(
    message = paste0(
      "Buffer report: ", length(affected), " observation(s) need attention.\n",
      "netCDF variable: ", variable, " \u00b7 output column: ", output, "\n",
      paste(examples, collapse = "\n"),
      if (length(affected) > 5L) paste0("\nShowing 5 of ", length(affected), " cases.") else "",
      "\nCheck coordinates, source coverage and missing values. ",
      "Change the window or threshold only if appropriate for your analysis.\n"
    ), call = NULL
  ), class = c("ocean3d_buffer_report", "message", "condition")))
  invisible(NULL)
}


#' Summarise environmental values around observations
#'
#' Calculate one summary per observation and variable within horizontal, depth
#' and/or temporal windows. Combine complementary netCDF sources before
#' assessing coverage. Observations and input order are retained.
#'
#' @param data Observation table, matrix, named list or sf POINT object.
#' @param nc netCDF path, paths, opened ncdf4 connection, list, or file table.
#' @param var Variable name or vector of names. NULL detects one data variable
#'   per file, excluding referenced bounds and coordinate variables.
#' @param horizontal_radius_km Horizontal geodesic radius in kilometres.
#'   Positive radii below 10000 km are supported.
#' @param lon_cell_radius,lat_cell_radius Optional integer radii in grid cells.
#'   A radius of one on both axes requests a 3 by 3 neighbourhood. Mutually
#'   exclusive with a positive horizontal_radius_km.
#' @param depth_radius_m Symmetric depth radius in metres, positive downward.
#' @param depth_layer_radius Optional integer radius in physical depth order.
#'   Mutually exclusive with a positive depth_radius_m.
#' @param time_before,time_after Non-negative temporal distances before and
#'   after the observation. Both use time_unit.
#' @param time_unit One of days, hours, minutes, seconds, months or years.
#'   Days are 24 elapsed hours. Months and years require whole numbers and
#'   shift the UTC calendar, clamping to the destination month's last day.
#' @param statistic One of "mean", "sd", "se", "min" or "max".
#'   SD uses the sample variance correction. SE is SD divided by the square
#'   root of the effective sample size. Min and max ignore weight magnitudes.
#' @param weighting "weighted" uses intersecting area, thickness and duration;
#'   "equal" gives equal weight to valid samples. Coverage uses physical
#'   support in both modes. A point-selected axis contributes a factor of one.
#' @param min_valid_fraction Minimum valid fraction of the entire requested
#'   support, between zero and one. Empty or entirely missing buffers return NA.
#' @param infer_bounds Allow missing bounds to be estimated using neighbouring
#'   coordinate midpoints, with half-gap extrapolation at outer edges. Inferred
#'   depth support is clipped at zero metres. A singleton needs explicit bounds.
#' @param method Depth selection on an unbuffered vertical axis: nearest,
#'   surface, bottom or 2d. NULL uses 2d without a depth axis and nearest with
#'   one. Leave NULL when requesting a depth interval or layer neighbourhood.
#' @param lon_col,lat_col,depth_col,date_col Observation coordinate columns.
#'   Dates are required only for temporal variables. Date and ISO date-only
#'   text centre time buffers at midnight UTC; POSIXt preserves instants.
#' @param id_col Optional observation identifier column used in diagnostics.
#' @param file_col File-path column when nc is a table.
#' @param lon_dim,lat_dim,depth_dim,time_dim Optional netCDF dimension overrides.
#' @param time_match,time_resolution,max_time_diff Point-time selection options
#'   on an unbuffered time axis; see [extract_to_point()]. max_time_diff is in
#'   days and cannot be combined with a temporal window.
#' @param spatial_match,max_distance_km Point-cell or index-seed selection
#'   options; see [extract_to_point()]. They do not expand a kilometre buffer.
#' @param output_col Output name for one variable, or a named character vector
#'   mapping multiple variable names to distinct output columns.
#' @param return_match_info Add a list column named <output>_match_info per
#'   variable, containing coverage, source, bounds and selection diagnostics.
#' @param quiet Hide explanatory summaries. Time-zone warnings remain active.
#'
#' @returns The normalised observation table with one numeric summary column
#'   per variable and, optionally, a diagnostic list column per output. sf
#'   inputs retain their geometry. Existing columns are never overwritten.
#'
#' @section Support and limitations:
#' For normalised weights p, corrected variance is sum(p * (x - mean)^2)
#' divided by sum(p * (1 - p)); effective sample size is 1 / sum(p^2).
#' Equal weighting reproduces ordinary sample SD and SE. SD and SE require
#' at least two valid samples; otherwise they return NA with a diagnostic.
#' SE assumes independent samples and does not correct spatial or temporal
#' autocorrelation. It is not an automatic estimate of environmental uncertainty.
#' All statistics, including min and max, require the requested valid coverage.
#'
#' At least one window must be positive. Bounds describe represented support;
#' inferred intervals do not establish that instantaneous values are interval
#' means. The calculation treats each sample as constant over its assigned
#' support. No interpolation or nearest fallback expands an empty buffer.
#'
#' Geographic grids must have separable one-dimensional longitude and latitude
#' axes. Circles use 720 segments, with WGS84 distance and equal-area projections;
#' cell edges are subdivided at at most 0.01 degrees before projection. This is
#' a polygonal approximation, independent of the caller's sf_use_s2 setting.
#' Unsupported projections, calendars, climatologies or transformed vertical
#' coordinates require explicit preprocessing.
#'
#' Complementary sources are combined before applying min_valid_fraction.
#' Overlapping selected support is an error, including duplicate values or
#' a missing copy. Sources for one variable must use matching value units.
#' Point-time selection is resolved across eligible sources before aggregation.
#' Index neighbourhoods require matching grid/layer coordinates and bounds
#' across sources. If requested indices extend beyond available coverage, the
#' result is NA because the full requested physical support is unknown.
#'
#' @examples
#' \dontrun{
#' observations <- data.frame(lon = 0, lat = 38, depth = 60,
#'                            date = as.Date("2020-01-15"))
#' extract_to_buffer(observations, "temperature.nc", var = "thetao",
#'                   horizontal_radius_km = 5, method = "surface")
#' extract_to_buffer(observations, "temperature.nc", var = "thetao",
#'                   depth_radius_m = 10, time_before = 3)
#' }
#' @seealso [extract_to_point()], [extract_to_area()]
#' @export
extract_to_buffer <- function(
  data, nc, var = NULL,
  horizontal_radius_km = 0,
  lon_cell_radius = NULL, lat_cell_radius = NULL,
  depth_radius_m = 0, depth_layer_radius = NULL,
  time_before = 0, time_after = 0,
  time_unit = c("days", "hours", "minutes", "seconds", "months", "years"),
  statistic = "mean", weighting = c("weighted", "equal"),
  min_valid_fraction = 0.8, infer_bounds = FALSE, method = NULL,
  lon_col = "lon", lat_col = "lat", depth_col = "depth", date_col = "date",
  id_col = NULL, file_col = "file",
  lon_dim = NULL, lat_dim = NULL, depth_dim = NULL, time_dim = NULL,
  time_match = c("nearest", "exact"),
  time_resolution = c("auto", "datetime", "daily", "monthly", "yearly"),
  max_time_diff = NULL,
  spatial_match = c("within", "nearest"), max_distance_km = NULL,
  output_col = NULL, return_match_info = TRUE, quiet = FALSE
) {
  options <- .buffer_validate_request(
    horizontal_radius_km, depth_radius_m,
    lon_cell_radius, lat_cell_radius, depth_layer_radius,
    time_before = time_before, time_after = time_after, time_unit = time_unit,
    statistic = statistic, weighting = weighting,
    min_valid_fraction = min_valid_fraction, infer_bounds = infer_bounds,
    method = method, return_match_info = return_match_info, quiet = quiet
  )
  spatial <- horizontal_radius_km > 0 || options$horizontal_index
  vertical <- depth_radius_m > 0 || options$depth_index
  temporal <- time_before > 0 || time_after > 0
  time_match <- match.arg(time_match)
  time_resolution <- match.arg(time_resolution)
  spatial_match <- .validate_spatial_matching(spatial_match, max_distance_km)
  if (horizontal_radius_km >= 10000) {
    stop("Horizontal radii must be below 10000 km.", call. = FALSE)
  }
  if (horizontal_radius_km > 0 && spatial_match != "within") {
    stop("Kilometre buffers use their requested footprint; remove nearest matching.",
         call. = FALSE)
  }
  if (!is.null(max_time_diff) &&
        (temporal || !is.numeric(max_time_diff) || length(max_time_diff) != 1L ||
         !is.finite(max_time_diff) || max_time_diff < 0)) {
    stop("`max_time_diff` needs a finite non-negative value and an unbuffered time axis.",
         call. = FALSE)
  }
  column_names <- list(lon_col = lon_col, lat_col = lat_col, file_col = file_col)
  optional_names <- list(depth_col = depth_col, date_col = date_col, id_col = id_col,
                         lon_dim = lon_dim, lat_dim = lat_dim,
                         depth_dim = depth_dim, time_dim = time_dim)
  column_names <- c(column_names, Filter(Negate(is.null), optional_names))
  for (name in names(column_names)) {
    value <- column_names[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("`", name, "` must name one column or dimension.", call. = FALSE)
    }
  }
  if (!is.null(var) && (!is.character(var) || !length(var) || anyNA(var) ||
                       any(!nzchar(var)) || anyDuplicated(var))) {
    stop("`var` must be NULL or distinct, non-empty variable names.", call. = FALSE)
  }
  data <- .normalize_points_input(data, lon_col, lat_col)
  if (anyDuplicated(names(data))) stop("Observation column names must be unique.", call. = FALSE)
  .check_required_cols(data, c(lon_col, lat_col, id_col))
  if (!is.numeric(data[[lon_col]]) || !is.numeric(data[[lat_col]])) {
    stop("Observation longitude and latitude must be numeric.", call. = FALSE)
  }

  sources <- .as_netcdf_sources(nc, file_col)
  if (!length(sources)) stop("At least one netCDF source is required.", call. = FALSE)
  opened <- list()
  on.exit(lapply(opened, .close_source), add = TRUE)
  metadata <- list()
  dimension_args <- list(lon_dim = lon_dim, lat_dim = lat_dim,
                          depth_dim = depth_dim, time_dim = time_dim)
  for (i in seq_along(sources)) {
    opened[[i]] <- .open_source(sources[[i]])
    handle <- opened[[i]]$nc
    variables <- if (is.null(var)) .buffer_detect_variable(handle) else {
      intersect(var, names(handle$var))
    }
    if (!length(variables)) {
      stop("No requested variable is present in ", basename(handle$filename), ".",
           call. = FALSE)
    }
    for (variable in variables) {
      metadata[[length(metadata) + 1L]] <- .buffer_source_metadata(
        handle, variable, i, options, method, dimension_args
      )
    }
  }
  variables <- unique(vapply(metadata, function(x) x$var, character(1)))
  if (!is.null(var) && !setequal(var, variables)) {
    stop("Requested variables were not found: ", toString(setdiff(var, variables)),
         call. = FALSE)
  }
  outputs <- .buffer_output_names(variables, output_col, names(data), return_match_info)
  has_time <- any(vapply(metadata, function(x) !is.null(x$dims$time), logical(1)))
  needs_depth <- vertical || any(vapply(metadata, function(x) x$method == "nearest", logical(1)))
  depth_values <- rep(NA_real_, nrow(data))
  window_data <- data
  if (needs_depth) {
    if (is.null(depth_col)) stop("The requested selection requires `depth_col`.", call. = FALSE)
    .check_required_cols(data, depth_col)
    depth_values <- data[[depth_col]]
    if (all(is.na(depth_values))) depth_values <- rep(NA_real_, nrow(data))
    if (!is.numeric(depth_values) ||
          any(!is.na(depth_values) & (!is.finite(depth_values) | depth_values < 0))) {
      stop("Observation depths must be finite, non-negative metres or NA.", call. = FALSE)
    }
    window_data[[depth_col]] <- depth_values
  }
  prepared <- NULL
  if (has_time) {
    for (meta in metadata) {
      if (!is.null(meta$dims$time)) {
        .require_observation_time(data, date_col, meta$var, meta$source_file)
      }
    }
    prepared <- .prepare_time_input(
      data[[date_col]], if (temporal) "auto" else time_resolution,
      input_name = paste0("`data$", date_col, "`"),
      ids = if (is.null(id_col)) NULL else data[[id_col]]
    )
  }
  windows <- .buffer_observation_windows(window_data, options, depth_col, date_col, id_col)$windows
  out <- data
  reports <- list()
  for (variable in variables) {
    group <- Filter(function(x) x$var == variable, metadata)
    .buffer_check_index_grids(group, options)
    methods <- unique(vapply(group, function(x) x$method, character(1)))
    dynamic <- vapply(group, function(x) !is.null(x$dims$time), logical(1))
    if (length(methods) != 1L || (any(dynamic) && !all(dynamic))) {
      stop("Sources for `", variable, "` have incompatible static/depth dimensions.", call. = FALSE)
    }
    for (meta in group) {
      if (!identical(meta$variable_units, group[[1L]]$variable_units)) {
        stop("Sources for `", variable, "` use different value units. Convert them first.",
             call. = FALSE)
      }
    }
    values <- rep(NA_real_, nrow(data))
    diagnostics <- vector("list", nrow(data))
    for (row in seq_len(nrow(data))) {
      lon <- data[[lon_col]][row]
      lat <- data[[lat_col]][row]
      status <- if (!is.finite(lon) || !is.finite(lat) || lat < -90 || lat > 90) {
        "invalid_coordinates"
      } else if (any(dynamic) && is.na(prepared$values[row])) {
        "missing_time"
      } else if ((vertical || methods == "nearest") && is.na(depth_values[row])) {
        "missing_depth"
      } else NULL
      samples <- .buffer_empty_samples()
      request <- list(
        longitude = lon, latitude = lat,
        horizontal_radius_km = horizontal_radius_km,
        lon_cell_radius = lon_cell_radius, lat_cell_radius = lat_cell_radius,
        depth_radius_m = depth_radius_m, depth_layer_radius = depth_layer_radius,
        depth = depth_values[row], time_before = time_before, time_after = time_after,
        date_utc = if (any(dynamic)) prepared$values[row] else NULL,
        time_unit = options$time_unit,
        time_window_utc = if (temporal) c(windows$time_lower_utc[row], windows$time_upper_utc[row]) else NULL
      )
      if (is.null(status)) {
        buffer <- if (spatial) .buffer_horizontal_geometry(
          lon, lat, if (options$horizontal_index) 1 else horizontal_radius_km
        ) else NULL
        queries <- lapply(group, function(meta) {
          space <- .buffer_space_query(meta, lon, lat, options, buffer,
                                        spatial_match, max_distance_km)
          depth <- .buffer_depth_query(meta, depth_values[row], options)
          list(meta = meta, space = space, depth = depth)
        })
        if (any(vapply(queries, function(x) {
          x$space$status == "incomplete_index_window" || x$depth$status == "incomplete_index_window"
        }, logical(1)))) {
          status <- "incomplete_index_window"
        } else {
          eligible <- Filter(function(x) {
            x$space$status == "ready" && nrow(x$space$cells) > 0 && length(x$depth$index) > 0
          }, queries)
          selected <- NULL
          if (any(dynamic) && !temporal && length(eligible)) {
            available <- sort(unique(unlist(lapply(eligible, function(x) x$meta$time_values))))
            times <- as.POSIXct(available, origin = "1970-01-01", tz = "UTC")
            index <- .match_time_index(times, prepared$values[row], time_match,
                                       max_time_diff, prepared$resolution)
            if (is.na(index)) status <- "time_rejected" else selected <- available[index]
          }
          if (is.null(status)) {
            pieces <- lapply(eligible, function(query) {
              time <- .buffer_time_query(query$meta, if (temporal) {
                as.numeric(c(windows$time_lower_utc[row], windows$time_upper_utc[row]))
              } else NULL, selected, temporal)
              .buffer_collect_source_samples(query$meta, query$space, query$depth, time)
            })
            if (length(pieces)) samples <- do.call(rbind, pieces)
            samples <- .buffer_select_point_samples(samples, methods, depth_values[row], spatial, vertical)
            if (spatial) attr(samples, "buffer_crs") <- sf::st_crs(buffer$geometry)
            .buffer_check_source_overlap(samples, group, spatial, vertical, temporal)
          }
        }
      }
      result <- if (is.null(status)) {
        .buffer_summarise_samples(samples, options$weighting, min_valid_fraction,
                                  options$statistic)
      } else .buffer_rejected_result(status)
      values[row] <- result$value
      diagnostics[[row]] <- .buffer_diagnostic(result, samples, group, options, request)
    }
    output <- outputs[[variable]]
    out[[output]] <- values
    if (return_match_info) out[[paste0(output, "_match_info")]] <- diagnostics
    reports[[variable]] <- diagnostics
  }
  if (!is.null(prepared)) .warn_time_zone(prepared, paste0("`data$", date_col, "`"))
  if (!quiet) {
    for (variable in variables) {
      .show_buffer_report(reports[[variable]], data, variable, outputs[[variable]], id_col)
    }
  }
  out
}
