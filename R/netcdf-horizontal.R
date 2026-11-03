# Longitude and latitude axes and horizontal cell matching. Coordinates are
# checked and given cell edges (from CF bounds or inferred from the centres),
# longitudes are compared in either convention (0-360 or -180-180), and each
# observation is matched to the cell that contains it or to the nearest cell
# by great-circle distance. Shared by point, buffer and track extraction.

.validate_spatial_matching <- function(spatial_match, max_distance_km) {
  spatial_match <- match.arg(spatial_match, c("within", "nearest"))
  if (spatial_match == "within") {
    if (!is.null(max_distance_km)) {
      stop(
        "`max_distance_km` requires spatial_match = 'nearest'. ",
        "Use 'nearest' to allow nearby cells outside the grid, ",
        "or remove the distance argument.",
        call. = FALSE
      )
    }
    return(spatial_match)
  }
  valid_distance <- is.numeric(max_distance_km) &&
    length(max_distance_km) == 1L && !is.na(max_distance_km) &&
    is.finite(max_distance_km) && max_distance_km >= 0
  if (!valid_distance) {
    stop(
      "Nearest spatial matching requires one finite, ",
      "non-negative max_distance_km value. ",
      "For example, use max_distance_km = 5 for five kilometres.",
      call. = FALSE
    )
  }
  spatial_match
}

.normalise_longitude <- function(values) {
  (values + 180) %% 360 - 180
}

.containing_axis_indices <- function(axis, target) {
  if (length(target) != 1L || !is.finite(target)) return(integer())
  if (axis$longitude) {
    targets <- target + 360 * round((axis$centres - target) / 360)
  } else {
    if (target < -90 || target > 90) return(integer())
    targets <- rep(target, length(axis$centres))
  }
  # A singleton axis without bounds has zero width. Allow single-precision
  # rounding of the stored coordinate, but no real offset.
  tol <- if (identical(axis$basis, "singleton_coordinate")) {
    pmax(1e-8, abs(axis$centres) * 2^-23)
  } else {
    1e-8
  }
  which(targets >= axis$lower - tol & targets <= axis$upper + tol)
}

# Great-circle distance on a sphere with the s2 mean Earth radius, so values
# agree with sf::st_distance() under s2, independently of sf settings.
.great_circle_km <- function(lon1, lat1, lon2, lat2,
                             radius_km = .earth_radius_km) {
  p <- pi / 180
  h <- sin((lat2 - lat1) * p / 2)^2 +
    cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * radius_km * asin(pmin(1, sqrt(h)))
}

# Closest cell centre among the given longitude and latitude indices.
.nearest_horizontal_cell <- function(lon_values, lat_values,
                                     lon_indices, lat_indices,
                                     lon_target, lat_target) {
  if (!length(lon_indices) || !length(lat_indices)) return(NULL)
  lon_idx <- rep(lon_indices, times = length(lat_indices))
  lat_idx <- rep(lat_indices, each = length(lon_indices))
  distances <- .great_circle_km(lon_target, lat_target,
                                lon_values[lon_idx], lat_values[lat_idx])
  tied <- which(distances <= min(distances) + .distance_tie_km)
  position <- (lat_idx[tied] - 1) * length(lon_values) + lon_idx[tied]
  winner <- tied[which.min(position)]
  list(
    lon_idx = as.integer(lon_idx[winner]),
    lat_idx = as.integer(lat_idx[winner]),
    selected_lon = lon_values[lon_idx[winner]],
    selected_lat = lat_values[lat_idx[winner]],
    distance_km = distances[winner]
  )
}

# Closest cell centre over the whole rectilinear grid. On every latitude row
# the distance grows with the wrapped longitude difference, so only the
# longitude column(s) closest to the target can hold the nearest centre.
# This is exact and avoids computing a distance to every grid cell.
.nearest_horizontal_cells <- function(lon_values, lat_values,
                                      lon_targets, lat_targets) {
  result <- vector("list", length(lon_targets))
  if (!length(lon_values) || !length(lat_values)) return(result)
  valid <- which(is.finite(lon_targets) & is.finite(lat_targets) &
                   lat_targets >= -90 & lat_targets <= 90)
  for (k in valid) {
    lon_gap <- abs(.normalise_longitude(lon_values - lon_targets[k]))
    columns <- which(lon_gap <= min(lon_gap) + 1e-9)
    result[[k]] <- .nearest_horizontal_cell(
      lon_values, lat_values, columns, seq_along(lat_values),
      lon_targets[k], lat_targets[k]
    )
  }
  result
}

# Match every observation to a horizontal cell. Results are stored as one
# vector per field (not one list per row) to keep memory proportional to the
# number of observations.
.match_horizontal_cells <- function(lon_axis, lat_axis,
                                    lon_values, lat_values,
                                    lon_targets, lat_targets,
                                    spatial_match, max_distance_km) {
  n <- length(lon_targets)
  matches <- list(
    accepted = rep(FALSE, n), reason = rep("invalid_coordinates", n),
    lon_idx = rep(NA_integer_, n), lat_idx = rep(NA_integer_, n),
    selected_lon = rep(NA_real_, n), selected_lat = rep(NA_real_, n),
    distance_km = rep(NA_real_, n), inside_coverage = rep(FALSE, n),
    inside_selected_cell = rep(FALSE, n),
    requested_lon = lon_targets, requested_lat = lat_targets,
    spatial_match = spatial_match, max_distance_km = max_distance_km
  )
  nearest <- if (spatial_match == "nearest") {
    .nearest_horizontal_cells(lon_values, lat_values, lon_targets, lat_targets)
  }
  for (i in seq_len(n)) {
    lon_target <- lon_targets[i]
    lat_target <- lat_targets[i]
    if (!is.finite(lon_target) || !is.finite(lat_target) ||
          lat_target < -90 || lat_target > 90) {
      next
    }
    containing_lon <- .containing_axis_indices(lon_axis, lon_target)
    containing_lat <- .containing_axis_indices(lat_axis, lat_target)
    candidate <- if (spatial_match == "nearest") {
      nearest[[i]]
    } else {
      .nearest_horizontal_cell(lon_values, lat_values, containing_lon,
                               containing_lat, lon_target, lat_target)
    }
    if (is.null(candidate)) {
      matches$reason[i] <- "outside_grid"
      next
    }
    accepted <- spatial_match == "within" ||
      candidate$distance_km <= max_distance_km
    matches$accepted[i] <- accepted
    matches$reason[i] <- if (accepted) "matched" else "distance_limit_exceeded"
    matches$lon_idx[i] <- candidate$lon_idx
    matches$lat_idx[i] <- candidate$lat_idx
    matches$selected_lon[i] <- candidate$selected_lon
    matches$selected_lat[i] <- candidate$selected_lat
    matches$distance_km[i] <- candidate$distance_km
    matches$inside_coverage[i] <- length(containing_lon) > 0L &&
      length(containing_lat) > 0L
    matches$inside_selected_cell[i] <-
      candidate$lon_idx %in% containing_lon &&
      candidate$lat_idx %in% containing_lat
  }
  matches
}

# One observation's horizontal match, for reports.
.horizontal_match_row <- function(matches, i) {
  settings <- c("spatial_match", "max_distance_km")
  row <- lapply(matches[setdiff(names(matches), settings)], `[[`, i)
  row$spatial_match <- matches$spatial_match
  row["max_distance_km"] <- list(matches$max_distance_km)
  row
}

# Longitudes as one continuous monotonic axis. An axis that is already
# monotonic is left as it is, however large its steps (#95). Otherwise the
# axis may cross the file's longitude seam (e.g. 350, 10, 30): take every step
# the same way around the globe, east or west, and keep the direction that
# covers less than a full circle. At most one direction can, for three or more
# centres. If neither does, the values are returned unchanged and the
# monotonic check in .prepare_spatial_axis() rejects them.
.unwrap_longitude <- function(values) {
  steps <- diff(values)
  if (all(steps > 0) || all(steps < 0)) return(values)
  east <- steps %% 360
  west <- east - 360
  if (any(east == 0)) return(values)
  if (sum(east) < 360) return(values[1L] + c(0, cumsum(east)))
  if (sum(west) > -360) return(values[1L] + c(0, cumsum(west)))
  values
}

.prepare_spatial_axis <- function(nc, coordinate, values, longitude = FALSE) {
  if (!is.numeric(values) || !length(values) ||
        any(!is.finite(values))) {
    stop("Spatial coordinates must be finite numbers.", call. = FALSE)
  }
  if (!longitude && any(values < -90 | values > 90)) {
    stop("Latitude coordinates are outside [-90, 90].", call. = FALSE)
  }
  centres <- if (longitude) .unwrap_longitude(values) else values
  steps <- diff(centres)
  if (length(steps) && !(all(steps > 0) || all(steps < 0))) {
    seam <- if (longitude) {
      " Longitudes may cross the 0/360 or -180/180 seam once."
    } else {
      ""
    }
    rlang::abort(glue::glue(
      "Spatial coordinates for `{coordinate}` must form a monotonic axis ",
      "(as CF requires).{seam} ",
      "Sort the file along this dimension before extraction."
    ), class = "ocean3d_error_spatial_axis", call = NULL)
  }
  explicit <- .read_bounds(nc, coordinate, length(values))$values
  if (!is.null(explicit)) {
    intervals <- explicit
    if (longitude) {
      for (j in seq_along(values)) {
        width <- abs(diff(explicit[, j]))
        if (width >= 360) {
          intervals[, j] <- centres[j] + c(-180, 180)
        } else {
          offsets <- (explicit[, j] - values[j] + 180) %% 360 - 180
          intervals[, j] <- centres[j] + offsets
        }
      }
    }
    lower <- apply(intervals, 2L, min)
    upper <- apply(intervals, 2L, max)
    if (any(centres < lower - 1e-8 | centres > upper + 1e-8)) {
      warning(
        "Spatial bounds for `", coordinate, "` do not contain their cell ",
        "centres; inferring cell edges from the coordinate centres instead.",
        call. = FALSE
      )
      explicit <- NULL
    }
  }
  if (!is.null(explicit)) {
    basis <- "explicit_bounds"
  } else if (length(centres) == 1L) {
    lower <- upper <- centres
    basis <- "singleton_coordinate"
  } else {
    order <- order(centres)
    sorted <- centres[order]
    gaps <- diff(sorted)
    edges <- c(
      sorted[1L] - gaps[1L] / 2,
      head(sorted, -1L) + gaps / 2,
      tail(sorted, 1L) + tail(gaps, 1L) / 2
    )
    lower <- upper <- numeric(length(centres))
    lower[order] <- head(edges, -1L)
    upper[order] <- tail(edges, -1L)
    regular <- all(
      abs(gaps - gaps[1L]) <= max(1e-8, gaps[1L] * 1e-8)
    )
    basis <- if (regular) "inferred_regular" else "inferred_irregular"
  }
  if (!longitude) {
    lower <- pmax(lower, -90)
    upper <- pmin(upper, 90)
  }
  list(
    centres = centres, lower = lower, upper = upper,
    longitude = longitude, basis = basis
  )
}
