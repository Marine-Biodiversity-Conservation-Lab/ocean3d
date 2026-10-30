#-------------------------------------------------------------------------------
# Environmental extraction along ordered tracks
#-------------------------------------------------------------------------------

# Internal helpers
# -------------------------------------------------------------

.track_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
}


.track_positive <- function(x, name, optional = TRUE) {
  if (optional && is.null(x)) {
    return(invisible(NULL))
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be one finite positive number.", call. = FALSE)
  }
  invisible(NULL)
}


.track_input <- function(
  data, track_col, lon_col, lat_col, depth_col,
  z_as_depth
) {
  if (inherits(data, "sf")) {
    if (is.na(sf::st_crs(data))) {
      stop("Track geometries need a CRS.", call. = FALSE)
    }
    types <- as.character(sf::st_geometry_type(data))
    if (any(types != "LINESTRING") || any(sf::st_is_empty(data))) {
      stop(
        paste0(
          "Use non-empty LINESTRING features; split multipart ",
          "lines explicitly."
        ),
        call. = FALSE
      )
    }
    if (!nrow(data)) {
      stop("Supply at least one track.", call. = FALSE)
    }
    lines <- sf::st_geometry(sf::st_transform(data, 4326))
    attributes <- sf::st_drop_geometry(data)
    ids <- if (is.null(track_col)) {
      as.character(seq_len(nrow(data)))
    } else {
      .check_required_cols(attributes, track_col)
      as.character(attributes[[track_col]])
    }
    if (anyDuplicated(ids)) {
      stop("Each LINESTRING needs a distinct track ID.", call. = FALSE)
    }
    rows <- lapply(
      seq_along(lines), function(i) {
        coordinates <- unclass(lines[[i]])
        if (inherits(lines[[i]], "XYM") || inherits(lines[[i]], "XYZM")) {
          stop(
            "Measured (M) lines are unsupported; supply time explicitly.",
            call. = FALSE
          )
        }
        n <- nrow(coordinates)
        if (n < 2L) {
          stop("Each track needs at least two vertices.", call. = FALSE)
        }
        row <- attributes[rep(i, n), , drop = FALSE]
        for (name in names(attributes)) {
          if (is.list(attributes[[name]])) {
            values <- attributes[[name]][[i]]
            if (!length(values) %in% c(1L, n)) {
              stop(
                "List column `", name,
                "` must have one value or one per vertex ",
                "in feature ", i, ".",
                call. = FALSE
              )
            }
            row[[name]] <- rep(values, length.out = n)
          }
        }
        if (ncol(coordinates) == 3L) {
          if (!z_as_depth) {
            stop(
              paste0(
                "XYZ lines require `z_as_depth = TRUE`: Z must be ",
                "metres positive down."
              ),
              call. = FALSE
            )
          }
          if (depth_col %in% names(attributes)) {
            stop(
              "Supply depth either as Z or in `depth_col`, not both.",
              call. = FALSE
            )
          }
          row[[depth_col]] <- coordinates[, 3L]
        }
        row[[lon_col]] <- coordinates[, 1L]
        row[[lat_col]] <- coordinates[, 2L]
        row$track_id <- ids[i]
        row
      }
    )
    result <- do.call(rbind, rows)
  } else {
    result <- .normalize_points_input(data, lon_col, lat_col)
    .check_required_cols(result, c(lon_col, lat_col, track_col))
    if (!nrow(result)) {
      stop("Supply at least one track.", call. = FALSE)
    }
    result$track_id <- if (is.null(track_col)) {
      "1"
    } else {
      as.character(result[[track_col]])
    }
  }
  if (anyDuplicated(names(result))) {
    stop("Input column names must be unique.", call. = FALSE)
  }
  if (anyNA(result$track_id) || any(!nzchar(result$track_id))) {
    stop("Track IDs must be non-missing and non-empty.", call. = FALSE)
  }
  for (name in c(lon_col, lat_col)) {
    if (!is.numeric(result[[name]]) || any(!is.finite(result[[name]]))) {
      stop(
        "`", name,
        "` must contain finite coordinates. Split tracks at missing ",
        "locations; unknown gaps are not joined.",
        call. = FALSE
      )
    }
  }
  if (any(abs(result[[lat_col]]) > 90)) {
    stop("Latitudes must be within [-90, 90].", call. = FALSE)
  }
  if (any(table(result$track_id) < 2L)) {
    stop(
      "Each track needs at least two vertices, in traversal order.",
      call. = FALSE
    )
  }
  rownames(result) <- NULL
  result
}


.track_sources <- function(nc, var, method, file_col, dimension_args) {
  sources <- .as_netcdf_sources(nc, file_col)
  if (!length(sources)) {
    stop("Supply at least one netCDF source.", call. = FALSE)
  }
  opened <- list()
  on.exit(lapply(opened, .close_source), add = TRUE)
  metadata <- list()
  for (i in seq_along(sources)) {
    opened[[i]] <- .open_source(sources[[i]])
    handle <- opened[[i]]$nc
    variables <- if (is.null(var)) {
      .detect_netcdf_var(handle)
    } else {
      intersect(var, names(handle$var))
    }
    if (!length(variables)) {
      stop(
        "No requested variable in ", basename(handle$filename),
        ".",
        call. = FALSE
      )
    }
    for (variable in variables) {
      dims <- do.call(
        .detect_dims, c(
          list(nc = handle, var = variable),
          dimension_args
        )
      )
      selected_method <- if (is.null(method)) {
        if (is.null(dims$depth)) {
          "2d"
        } else {
          "nearest"
        }
      } else {
        method
      }
      metadata[[length(metadata) + 1L]] <- list(
        variable = variable,
        source = sources[[i]], file = handle$filename, method = selected_method,
        temporal = !is.null(dims$time), longitude = .get_dim_values(
          handle,
          dims$lon
        ), latitude = .get_dim_values(handle, dims$lat),
        lon_axis = .prepare_spatial_axis(
          handle, dims$lon, .get_dim_values(
            handle,
            dims$lon
          ),
          longitude = TRUE
        ), lat_axis = .prepare_spatial_axis(
          handle,
          dims$lat, .get_dim_values(handle, dims$lat)
        ), depth_m = if (is.null(dims$depth)) {
          NULL
        } else {
          .read_depth_axis(
            handle, dims$depth, .get_dim_values(
              handle,
              dims$depth
            )
          )$depth_m
        }, time_values = if (is.null(dims$time)) {
          numeric()
        } else {
          as.numeric(.decode_netcdf_time(handle, dims$time))
        }, units = .read_time_attribute(handle, variable, "units")
      )
    }
  }
  found <- unique(vapply(metadata, function(x) x$variable, character(1)))
  if (!is.null(var) && !setequal(found, var)) {
    stop(
      "Variables not found: ", toString(setdiff(var, found)),
      ".",
      call. = FALSE
    )
  }
  for (variable in found) {
    group <- Filter(
      function(x) {
        x$variable == variable
      }, metadata
    )
    compatible <- vapply(
      group, function(x) {
        identical(x$method, group[[1L]]$method) && identical(
          x$temporal,
          group[[1L]]$temporal
        ) && identical(x$units, group[[1L]]$units)
      }, logical(1)
    )
    if (!all(compatible)) {
      stop(
        "Sources for `", variable,
        "` have incompatible dimensions or value units.",
        call. = FALSE
      )
    }
  }
  metadata
}


.track_interpolate <- function(x, i, fractions, enabled, name, used) {
  a <- x[i]
  b <- x[i + 1L]
  if (
      any(fractions > 0 & fractions < 1) && used && !enabled && !is.na(a) &&
        !is.na(b) && a != b) {
    stop(
      "`", name, paste0(
        "` changes between vertices. Enable its interpolation ",
        "explicitly "
      ), "or use sampling = 'vertices'.",
      call. = FALSE
    )
  }
  values <- rep(NA_real_, length(fractions))
  if (enabled || (!is.na(a) && !is.na(b) && a == b)) {
    values <- a + fractions * (b - a)
  }
  values[fractions == 0] <- a
  values[fractions == 1] <- b
  values
}


.track_sampling_counts <- function(samples) {
  rows <- lapply(unique(samples$track_id), function(id) {
    selected <- samples[samples$track_id == id, , drop = FALSE]
    data.frame(
      track_id = id, n_observed = sum(!selected$is_interpolated),
      n_assumed = sum(selected$is_interpolated), n_samples = nrow(selected)
    )
  })
  do.call(rbind, rows)
}


.track_summaries <- function(
  samples, columns, statistics, weighting,
  threshold
) {
  if (!length(statistics)) {
    return(NULL)
  }
  rows <- list()
  for (id in unique(samples$track_id)) {
    for (column in columns) {
      selected <- samples[samples$track_id == id, , drop = FALSE]
      cell_intervals <- "sample_kind" %in% names(selected) &&
        any(selected$sample_kind == "cell_interval")
      weights <- switch(weighting,
        distance = selected$support_km,
        duration = selected$support_seconds,
        equal = if (cell_intervals) {
          as.numeric(selected$sample_kind == "cell_interval")
        } else {
          rep(1, nrow(selected))
        }
      )
      usable <- is.finite(weights) & weights > 0
      available <- switch(weighting,
        distance = selected$available_km,
        duration = selected$available_seconds,
        equal = weights *
          ifelse(
            selected$support_km > 0,
            selected$available_km / selected$support_km,
            as.numeric(
              !selected$is_interpolated | selected$segment_status %in%
                c("marine", "not_checked")
            )
          )
      )
      if (is.null(available)) {
        available <- weights
      }
      valid <- usable & is.finite(selected[[column]]) & is.finite(available) &
        available > 0
      total <- sum(weights[usable])
      # Missing duration support makes the full denominator
      # unknowable.
      unknown <- any(!is.finite(weights))
      fraction <- if (unknown || total <= 0) {
        NA_real_
      } else {
        sum(available[valid]) / total
      }
      accepted <- is.finite(fraction) && any(valid) && fraction +
        1e-10 >= threshold
      status <- if (!is.finite(fraction)) {
        "unknown_or_zero_support"
      } else if (!any(valid)) {
        "all_missing"
      } else if (!accepted) {
        "insufficient_coverage"
      } else if (fraction < 1 - 1e-10) {
        "partial"
      } else {
        "complete"
      }
      row <- data.frame(
        track_id = id, output_column = column,
        weighting = weighting, n_samples = nrow(selected),
        n_observed = sum(!selected$is_interpolated),
        n_assumed = sum(selected$is_interpolated), n_valid = sum(valid),
        valid_fraction = fraction, total_distance_km = sum(selected$support_km),
        total_duration_seconds = if (anyNA(selected$support_seconds)) {
          NA_real_
        } else {
          sum(selected$support_seconds)
        }, effective_n = if (any(valid)) {
          1 / sum((available[valid] / sum(available[valid]))^2)
        } else {
          NA_real_
        }, status = status, stringsAsFactors = FALSE
      )
      for (statistic in statistics) {
        row[[statistic]] <- if (accepted) {
          .buffer_statistic(
            selected[[column]][valid], available[valid],
            statistic
          )
        } else {
          NA_real_
        }
      }
      variation <- intersect(statistics, c("sd", "se"))
      row$variation_status <- if (!length(variation)) {
        "not_requested"
      } else if (accepted && sum(valid) < 2L) {
        "insufficient_samples"
      } else if (accepted && all(is.finite(unlist(row[variation])))) {
        "available"
      } else {
        "not_available"
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  do.call(rbind, rows)
}


# Geodesics, calendar steps and grid intersections for tracks
# -------------------

.track_check_cell_sources <- function(metadata) {
  variables <- unique(vapply(metadata, function(m) m$variable, character(1)))
  for (variable in variables) {
    sources <- Filter(function(m) m$variable == variable, metadata)
    if (length(sources) < 2L) next
    for (i in seq_len(length(sources) - 1L)) {
      for (j in seq.int(i + 1L, length(sources))) {
        a <- sources[[i]]
        b <- sources[[j]]
        shift <- 360 * round((mean(a$longitude) - mean(b$longitude)) / 360)
        x <- min(max(a$lon_axis$upper), max(b$lon_axis$upper + shift)) -
          max(min(a$lon_axis$lower), min(b$lon_axis$lower + shift))
        y <- min(max(a$lat_axis$upper), max(b$lat_axis$upper)) -
          max(min(a$lat_axis$lower), min(b$lat_axis$lower))
        same_lon <- isTRUE(all.equal(
          sort(a$longitude), sort(b$longitude + shift), tolerance = 1e-10
        ))
        same_lat <- isTRUE(all.equal(
          sort(a$latitude), sort(b$latitude), tolerance = 1e-10
        ))
        same <- same_lon && same_lat
        singleton_x <- length(a$longitude) == 1L ||
          length(b$longitude) == 1L
        singleton_y <- length(a$latitude) == 1L ||
          length(b$latitude) == 1L
        shared_x <- x > 1e-9 || (x >= -1e-9 && singleton_x)
        shared_y <- y > 1e-9 || (y >= -1e-9 && singleton_y)
        if (shared_x && shared_y && !same) {
          stop("Cells mode cannot combine overlapping different grids. ",
               "Use balanced sampling or align the sources explicitly.",
               call. = FALSE)
        }
      }
    }
  }
  invisible(NULL)
}


.track_geodesic <- function(lon, lat, model = "WGS84") {
  radians <- pi / 180
  xyz <- cbind(
    cos(lat * radians) * cos(lon * radians), cos(lat * radians) *
      sin(lon * radians), sin(lat * radians)
  )
  dot <- rowSums(xyz[-nrow(xyz), , drop = FALSE] * xyz[-1L, , drop = FALSE])
  previous <- xyz[-nrow(xyz), , drop = FALSE]
  following <- xyz[-1L, , drop = FALSE]
  cross <- cbind(
    previous[, 2L] * following[, 3L] - previous[, 3L] *
      following[, 2L], previous[, 3L] * following[, 1L] - previous[
      ,
      1L
    ] * following[, 3L], previous[, 1L] * following[, 2L] - previous[
      ,
      2L
    ] * following[, 1L]
  )
  angle <- atan2(sqrt(rowSums(cross^2)), dot)
  if (any(pi - angle < 1e-08)) {
    stop("Antipodal endpoints need an intermediate waypoint.", call. = FALSE)
  }
  crs <- vapply(
    seq_len(length(lon) - 1L), function(i) {
      paste0("+proj=aeqd +lat_0=", lat[i], " +lon_0=", lon[i],
             " +datum=WGS84 +units=m +no_defs")
    }, character(1)
  )
  xy <- lapply(
    seq_along(crs), function(i) {
      sf::sf_project(
        "EPSG:4326", crs[i], cbind(
          lon[i + 1L], lat[i +
                             1L]
        )
      )
    }
  )
  lengths <- if (model == "sphere") {
    angle * 6371.0088
  } else {
    vapply(xy, function(x) sqrt(sum(x^2)) / 1000, numeric(1))
  }
  if (any(!is.finite(lengths))) {
    stop(
      "A geodesic could not be constructed; add a waypoint.",
      call. = FALSE
    )
  }
  list(
    xyz = xyz, angle = angle, length_km = lengths, crs = crs, xy = xy,
    longitude = lon, latitude = lat, model = model
  )
}


.track_geodesic_positions <- function(arc, segment, fractions) {
  if (arc$model == "WGS84") {
    xy <- outer(fractions, as.numeric(arc$xy[[segment]]))
    result <- sf::sf_project(arc$crs[segment], "EPSG:4326", xy)
  } else {
    angle <- arc$angle[segment]
    a <- arc$xyz[segment, ]
    b <- arc$xyz[segment + 1L, ]
    xyz <- if (angle < 1e-12) {
      matrix(rep(a, each = length(fractions)), ncol = 3L)
    } else {
      outer(sin((1 - fractions) * angle) / sin(angle), a) + outer(
        sin(
          fractions *
            angle
        ) / sin(angle), b
      )
    }
    result <- cbind(
      atan2(xyz[, 2L], xyz[, 1L]) * 180 / pi, atan2(
        xyz[
          ,
          3L
        ], sqrt(xyz[, 1L]^2 + xyz[, 2L]^2)
      ) * 180 / pi
    )
  }
  # Retain supplied endpoint coordinates, including stationary
  # segments.
  result[fractions == 0, ] <- c(arc$longitude[segment], arc$latitude[segment])
  result[fractions == 1, ] <- c(
    arc$longitude[segment + 1L], arc$latitude[segment +
                                                1L]
  )
  colnames(result) <- c("lon", "lat")
  result
}


.track_time_breaks <- function(start, end, step, unit, budget) {
  if (
      is.null(step) || !is.finite(start) || !is.finite(end) || end ==
        start) {
    return(c(0, 1))
  }
  .track_positive(step, "time_step", optional = FALSE)
  if (unit %in% c("months", "years")) {
    if (step != floor(step)) {
      stop(
        "Calendar time steps must be whole months or years.",
        call. = FALSE
      )
    }
    current <- as.POSIXct(start, origin = "1970-01-01", tz = "UTC")
    values <- start
    while (as.numeric(current) < end) {
      following <- .buffer_shift_time(current, step, unit)
      if (following <= current) {
        stop("Calendar step did not advance time.")
      }
      if (length(values) >= budget) {
        stop("time_step exceeds max_samples.")
      }
      values <- c(values, min(end, as.numeric(following)))
      current <- following
    }
    return((values - start) / (end - start))
  }
  seconds <- step * switch(unit,
    seconds = 1,
    minutes = 60,
    hours = 3600,
    days = 86400
  )
  n <- max(1, ceiling((end - start) / seconds))
  if (!is.finite(n) || n + 1 > budget) {
    stop("time_step exceeds max_samples.")
  }
  seq(0, 1, length.out = n + 1)
}


.track_midpoint_boundaries <- function(values, lower, upper, budget) {
  values <- sort(unique(values[is.finite(values)]))
  if (length(values) < 2L) {
    return(numeric())
  }
  steps <- diff(values)
  if (all(abs(steps - steps[1L]) < max(abs(steps[1L]) * 1e-10, 1e-12))) {
    first <- max(min(values), lower)
    last <- min(max(values), upper)
    if (last <= first) {
      return(numeric())
    }
    half <- steps[1L] / 2
    a <- max(0, floor((first - values[1L]) / half))
    b <- ceiling((last - values[1L]) / half)
    if (b - a > budget) {
      stop("Grid intersections exceed max_samples.")
    }
    return(values[1L] + seq(a, b) * half)
  }
  if (length(values)^2 > 2e+06) {
    stop(
      "Too many irregular coordinate pairs for cells mode. ",
      "Use a shorter source window or balanced sampling.",
      call. = FALSE
    )
  }
  mids <- unique(as.numeric(outer(values, values, "+") / 2))
  mids[mids > lower & mids < upper]
}


.track_grid_breaks <- function(
  arc, segment, metadata, depths, times,
  resolution, matching, budget
) {
  position <- function(f) .track_geodesic_positions(arc, segment, f)
  reference <- arc$longitude[segment]
  coordinate <- function(f, axis) {
    value <- position(f)[, axis]
    if (axis == 1L) {
      value <- reference + (value - reference + 180) %% 360 - 180
    }
    value
  }
  cuts <- c(0, 1)
  # Latitude may turn within one geodesic; bracket each monotonic
  # portion.
  turns <- c(
    0, 1, stats::optimize(
      function(f) coordinate(f, 2L), c(0, 1),
      tol = 1e-12
    )$minimum, stats::optimize(
      function(f) -coordinate(f, 2L),
      c(0, 1),
      tol = 1e-12
    )$minimum
  )
  turns <- sort(unique(round(turns, 12)))
  for (axis in 1:2) {
    brackets <- if (axis == 1L) {
      c(0, 1)
    } else {
      turns
    }
    extrema <- coordinate(brackets, axis)
    if (axis == 2L && any(abs(extrema) > 89.999999)) {
      stop("Cells mode needs waypoints away from a pole.", call. = FALSE)
    }
    edges <- unique(
      unlist(
        lapply(
          metadata, function(m) {
            a <- if (axis == 1L) {
              m$lon_axis
            } else {
              m$lat_axis
            }
            edges <- c(a$lower, a$upper)
            order <- order(a$centres)
            overlapping <- length(order) > 1L && any(
              tail(a$lower[order], -1L) < head(a$upper[order], -1L) - 1e-9
            )
            if (overlapping) {
              centres <- a$centres
              if (axis == 1L) {
                centres <- reference + (centres - reference + 180) %% 360 - 180
              }
              edges <- c(edges, .track_midpoint_boundaries(
                centres, min(extrema), max(extrema), budget
              ))
            }
            edges
          }
        )
      )
    )
    if (axis == 1L) {
      edges <- reference + (edges - reference + 180) %% 360 - 180
    }
    edges <- edges[is.finite(edges) & edges > min(extrema) & edges <
                     max(extrema)]
    if (length(edges) > budget) {
      stop("Grid intersections exceed max_samples.")
    }
    for (j in seq_len(length(brackets) - 1L)) {
      interval <- brackets[c(j, j + 1L)]
      at_ends <- coordinate(interval, axis)
      for (edge in edges[edges > min(at_ends) & edges < max(at_ends)]) {
        cuts <- c(
          cuts, stats::uniroot(
            function(f) {
              coordinate(f, axis) -
                edge
            }, interval,
            tol = 1e-12
          )$root
        )
      }
    }
  }
  if (all(is.finite(depths)) && diff(depths) != 0) {
    levels <- unique(
      unlist(
        lapply(
          metadata, function(m) {
            if (m$method %in% c("nearest", "all")) {
              m$depth_m
            } else {
              NULL
            }
          }
        )
      )
    )
    edges <- .track_midpoint_boundaries(
      levels, min(depths), max(depths),
      budget
    )
    cuts <- c(cuts, (edges - depths[1L]) / diff(depths))
  }
  if (all(is.finite(times)) && diff(times) > 0) {
    if (resolution %in% c("daily", "monthly", "yearly")) {
      step <- switch(resolution,
        daily = "day",
        monthly = "month",
        yearly = "year"
      )
      first <- .time_comparison_values(
        as.POSIXct(times[1L], origin = "1970-01-01", tz = "UTC"),
        resolution
      )
      seconds_min <- switch(resolution,
        daily = 86400,
        monthly = 28 *
          86400,
        yearly = 365 * 86400
      )
      if (diff(times) / seconds_min > budget) {
        stop("Time grid exceeds max_samples.")
      }
      edges <- as.numeric(
        seq(
          first, as.POSIXct(
            times[2L] + seconds_min,
            origin = "1970-01-01", tz = "UTC"
          ),
          by = step
        )
      )
    } else {
      values <- unique(
        unlist(
          lapply(
            metadata, function(m) {
              v <- sort(m$time_values)
              c(
                v[v >= times[1L] & v <= times[2L]], tail(
                  v[v < times[1L]],
                  1L
                ), head(v[v > times[2L]], 1L)
              )
            }
          )
        )
      )
      edges <- if (identical(matching$time_match, "exact")) {
        values
      } else {
        .track_midpoint_boundaries(
          values, times[1L], times[2L],
          budget
        )
      }
      if (!is.null(matching$max_time_diff)) {
        edges <- c(
          edges, values - matching$max_time_diff * 86400,
          values + matching$max_time_diff * 86400
        )
      }
    }
    cuts <- c(cuts, (edges - times[1L]) / diff(times))
  }
  cuts <- sort(
    unique(
      c(
        0, 1, cuts[is.finite(cuts) & cuts > 0 & cuts <
                     1]
      )
    )
  )
  if (length(cuts) > budget + 1) {
    stop("Grid intersections exceed max_samples.")
  }
  cuts
}


.track_segment_geometry <- function(arc, segment, step_km = 0.5) {
  # Geometry predicates inspect the complete connecting polyline,
  # not samples.
  n <- max(1, ceiling(arc$length_km[segment] / step_km))
  if (n > 1e+05) {
    stop("Mask geometry exceeds the segment safety limit.")
  }
  xy <- .track_geodesic_positions(
    arc, segment, seq(
      0, 1,
      length.out = n +
        1L
    )
  )
  xy[, 1L] <- .normalise_longitude(xy[, 1L])
  sf::st_sfc(sf::st_linestring(unname(xy)), crs = 4326)
}


.track_land_check <- function(
  arc, segment, mask, marine_values, mask_type,
  domain, step_km
) {
  # .prepare_track() validates `mask` and `domain` and unions polygons once.
  if (is.null(mask)) {
    return("not_checked")
  }
  line <- .track_segment_geometry(arc, segment, step_km)
  if (inherits(mask, "SpatRaster")) {
    wrapped <- sf::st_wrap_dateline(line, quiet = TRUE)
    route <- terra::project(terra::vect(wrapped), terra::crs(mask))
    # Raster extents are rectangles in native coordinates, not geodesic rings;
    # a line lies inside a rectangle exactly when its bounding box does.
    route_box <- as.vector(terra::ext(route))
    mask_box <- as.vector(terra::ext(mask))
    inside <- route_box[["xmin"]] >= mask_box[["xmin"]] &&
      route_box[["xmax"]] <= mask_box[["xmax"]] &&
      route_box[["ymin"]] >= mask_box[["ymin"]] &&
      route_box[["ymax"]] <= mask_box[["ymax"]]
    if (!inside) {
      return("unknown_mask_coverage")
    }
    cells <- terra::extract(
      mask, route, cells = TRUE, touches = TRUE, raw = TRUE
    )
    values <- cells[, 2L]
    if (!length(values) || any(!is.finite(values))) {
      return("unknown_mask_coverage")
    }
    if (any(!values %in% marine_values)) {
      "land_crossing"
    } else {
      "marine"
    }
  } else {
    # Polygon predicates run in the segment's local azimuthal projection.
    crs <- arc$crs[segment]
    local <- sf::st_transform(line, crs)
    geometry <- sf::st_transform(mask, crs)
    if (mask_type == "land" && length(sf::st_intersects(local,
                                                        geometry)[[1L]])) {
      return("land_crossing")
    }
    if (!is.null(domain)) {
      known <- sf::st_transform(domain, crs)
      if (!length(sf::st_covered_by(local, known)[[1L]])) {
        return("unknown_mask_coverage")
      }
    } else if (mask_type == "land") {
      return("unknown_mask_coverage")
    }
    if (mask_type == "marine") {
      inside <- length(sf::st_contains_properly(geometry, local)[[1L]]) >
        0
      if (!inside) {
        if (is.null(domain)) return("unknown_mask_coverage")
        return("land_crossing")
      }
    }
    "marine"
  }
}


# Track sampling and support
# ---------------------------------------------------

.track_sample_one <- function(
  vertices, rows, options, lon_col, lat_col,
  depth_col, date_col, time_values, depth_used, time_used, budget
) {
  arc <- .track_geodesic(
    vertices[[lon_col]][rows], vertices[[lat_col]][rows],
    options$distance_model
  )
  times <- if (is.null(time_values)) {
    rep(NA_real_, length(rows))
  } else {
    time_values[rows]
  }
  depths <- if (depth_col %in% names(vertices) &&
                  is.numeric(vertices[[depth_col]])) {
    vertices[[depth_col]][rows]
  } else {
    rep(NA_real_, length(rows))
  }
  dt <- diff(times)
  if (any(dt < 0, na.rm = TRUE)) {
    stop("Dates must follow traversal order.")
  }
  if (!is.null(options$max_segment_km) &&
        any(arc$length_km > options$max_segment_km)) {
    stop(
      paste0(
        "A segment exceeds max_segment_km; split the track or ",
        "add waypoints."
      )
    )
  }
  if (
      !is.null(options$max_gap_hours) && any(
        dt > options$max_gap_hours *
          3600,
        na.rm = TRUE
      )) {
    stop("A time gap exceeds max_gap_hours; split the track.")
  }
  land <- vapply(
    seq_along(dt), function(i) {
      .track_land_check(
        arc, i, options$land_mask, options$marine_values,
        options$mask_type, options$mask_domain, options$mask_step_km
      )
    }, character(1)
  )
  good <- land %in% c("marine", "not_checked")
  if (options$land_policy == "error" && any(!good)) {
    invalid <- which(!good)
    shown <- head(invalid, 5L)
    reasons <- c(
      land_crossing = "crosses or touches land",
      unknown_mask_coverage = "mask coverage is unknown"
    )
    examples <- paste0("segment ", shown, ": ", reasons[land[shown]])
    stop(
      "Track '", vertices$track_id[rows[1L]], "': ", length(invalid),
      " invalid marine segment(s). ", paste(examples, collapse = "; "),
      call. = FALSE
    )
  }
  pieces <- list()
  for (i in seq_along(dt)) {
    cuts <- c(0, 1)
    if (options$sampling != "vertices") {
      count <- if (is.null(options$spacing_km)) {
        1
      } else {
        max(1, ceiling(arc$length_km[i] / options$spacing_km))
      }
      if (
          !is.null(options$depth_step_m) && all(
            is.finite(
              depths[c(
                i,
                i + 1L
              )]
            )
          )) {
        count <- max(
          count, ceiling(
            abs(
              diff(
                depths[c(
                  i, i +
                    1L
                )]
              )
            ) / options$depth_step_m
          )
        )
      }
      if (!is.finite(count) || count + 1 > budget) {
        stop("Track sampling exceeds max_samples.")
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
                              any(!is.finite(depths[c(i, i + 1L)]))) {
              1
            } else {
              abs(diff(depths[c(i, i + 1L)])) * (last - first) /
                options$depth_step_m
            }
            n <- max(1, ceiling(horizontal), ceiling(vertical))
            if (length(cuts) + n > budget) {
              stop("Track sampling exceeds max_samples.")
            }
            cuts <- c(cuts, seq(first, last, length.out = n + 1L)[-1L])
          }
        } else {
          count <- max(count, length(time_cuts) - 1L)
          if (count + 1L > budget) {
            stop("Track sampling exceeds max_samples.")
          }
          cuts <- seq(0, 1, length.out = count + 1L)
        }
      } else {
        cuts <- sort(unique(c(seq(0, 1, length.out = count + 1L), time_cuts)))
      }
      if (options$sampling == "cells") {
        cuts <- sort(
          unique(
            c(
              cuts, .track_grid_breaks(
                arc, i,
                options$metadata, depths[c(i, i + 1L)], times[c(
                  i,
                  i + 1L
                )], options$resolution, options$matching, budget
              )
            )
          )
        )
      }
    }
    a <- head(cuts, -1L)
    b <- tail(cuts, -1L)
    fractions <- if (options$sampling == "cells") {
      (a + b) / 2
    } else {
      a
    }
    if (
        length(fractions) + sum(vapply(pieces, nrow, integer(1))) +
          1 > budget) {
      stop("Track sampling exceeds max_samples.")
    }
    xy <- .track_geodesic_positions(arc, i, fractions)
    original <- fractions == 0
    row <- vertices[rep(rows[i], length(fractions)), , drop = FALSE]
    for (name in setdiff(names(row), "track_id")) row[[name]][!original] <- NA
    row[[lon_col]] <- xy[, 1L]
    row[[lat_col]] <- xy[, 2L]
    if (depth_used || options$interpolate_depth) {
      row[[depth_col]] <- .track_interpolate(
        depths, i, fractions,
        options$interpolate_depth, "interpolate_depth", depth_used
      )
    }
    if (!is.null(time_values)) {
      row[[date_col]] <- as.POSIXct(
        .track_interpolate(
          times, i,
          fractions, options$interpolate_time, "interpolate_time",
          time_used
        ),
        origin = "1970-01-01", tz = "UTC"
      )
    }
    row$source_row <- ifelse(original, rows[i], NA_integer_)
    row$segment_id <- i
    row$segment_fraction <- fractions
    row$segment_start_fraction <- a
    row$segment_end_fraction <- b
    row$is_interpolated <- !original
    row$distance_km <- sum(head(arc$length_km, i - 1L)) + fractions *
      arc$length_km[i]
    row$segment_status <- land[i]
    row$support_km <- (b - a) * arc$length_km[i]
    row$support_seconds <- (b - a) * dt[i]
    row$available_km <- row$support_km * good[i]
    row$available_seconds <- row$support_seconds * good[i]
    pieces[[i]] <- row
  }
  # Keep original observations as zero-support rows in cells mode.
  if (options$sampling == "cells") {
    last <- vertices[rows, , drop = FALSE]
    if (!is.null(time_values)) {
      last[[date_col]] <- as.POSIXct(
        times,
        origin = "1970-01-01",
        tz = "UTC"
      )
    }
    last$source_row <- rows
    last$segment_id <- pmin(seq_along(rows), length(dt))
    last$segment_fraction <- c(rep(0, length(dt)), 1)
    last$segment_start_fraction <- last$segment_fraction
    last$segment_end_fraction <- last$segment_fraction
    last$is_interpolated <- FALSE
    last$distance_km <- c(0, cumsum(arc$length_km))
    last$segment_status <- "observation"
    last$support_km <- last$support_seconds <- 0
    last$available_km <- last$available_seconds <- 0
  } else {
    last <- tail(pieces[[length(pieces)]], 1L)
    last[intersect(names(vertices), names(last))] <- vertices[tail(
      rows,
      1L
    ), , drop = FALSE]
    if (!is.null(time_values)) {
      last[[date_col]] <- as.POSIXct(
        tail(times, 1L),
        origin = "1970-01-01",
        tz = "UTC"
      )
    }
    last$source_row <- tail(rows, 1L)
    last$segment_fraction <- 1
    last$segment_start_fraction <- last$segment_end_fraction <- 1
    last$is_interpolated <- FALSE
    last$distance_km <- sum(arc$length_km)
    last$support_km <- last$support_seconds <- 0
    last$available_km <- last$available_seconds <- 0
  }
  result <- do.call(rbind, c(pieces, list(last)))
  if (nrow(result) > budget) {
    stop("Track sampling exceeds max_samples.")
  }
  if (options$sampling != "cells") {
    for (
         pair in list(
           c("support_km", "available_km"), c(
             "support_seconds",
             "available_seconds"
           )
         )) {
      for (name in pair) {
        intervals <- head(result[[name]], -1L)
        result[[name]] <- (c(0, intervals) + c(intervals, 0)) / 2
      }
    }
  }
  result$sample_id <- seq_len(nrow(result))
  result$sample_kind <- if (options$sampling == "cells") {
    ifelse(result$is_interpolated, "cell_interval", "observation")
  } else {
    ifelse(result$is_interpolated, options$sampling, "observation")
  }
  if (options$sampling == "cells") {
    result <- result[order(result$distance_km), , drop = FALSE]
  }
  result
}


.track_collapse_cells <- function(
  samples, columns, context, lon_col,
  lat_col, depth_col, date_col
) {
  pieces <- list()
  for (id in unique(samples$track_id)) {
    rows <- which(samples$track_id == id)
    vertices <- context$vertices[context$vertices$track_id == id, ,
      drop = FALSE
    ]
    arc <- .track_geodesic(
      vertices[[lon_col]], vertices[[lat_col]],
      context$options$distance_model
    )
    for (segment in unique(samples$segment_id[rows])) {
      selected <- rows[samples$segment_id[rows] == segment &
                         samples$sample_kind[rows] ==
                           "cell_interval"]
      if (!length(selected)) {
        next
      }
      keys <- vapply(
        selected, function(i) {
          info <- lapply(
            columns, function(column) {
              x <- samples[[paste0(column, "_match_info")]][[i]]
              x[c(
                "status", "source_file", "lon_index", "lat_index",
                "depth_m", "time_index"
              )]
            }
          )
          values <- unname(unlist(samples[i, columns, drop = FALSE]))
          paste(
            utils::capture.output(
              dput(
                list(
                  info, samples$segment_status[i],
                  values
                )
              )
            ),
            collapse = ""
          )
        }, character(1)
      )
      groups <- cumsum(c(TRUE, keys[-1L] != head(keys, -1L)))
      for (group in unique(groups)) {
        indices <- selected[groups == group]
        row <- samples[indices[1L], , drop = FALSE]
        for (
             name in c(
               "support_km", "support_seconds", "available_km",
               "available_seconds"
             )) {
          row[[name]] <- sum(samples[[name]][indices])
        }
        row$segment_start_fraction <-
          min(samples$segment_start_fraction[indices])
        row$segment_end_fraction <- max(samples$segment_end_fraction[indices])
        fraction <- (row$segment_start_fraction + row$segment_end_fraction) / 2
        row$segment_fraction <- fraction
        xy <- .track_geodesic_positions(arc, segment, fraction)
        row[[lon_col]] <- xy[1L, 1L]
        row[[lat_col]] <- xy[1L, 2L]
        row$distance_km <- sum(
          head(
            arc$length_km, segment -
              1L
          )
        ) + arc$length_km[segment] * fraction
        if (
            context$options$interpolate_depth && depth_col %in%
              names(vertices)) {
          row[[depth_col]] <- .track_interpolate(
            vertices[[depth_col]],
            segment, fraction, TRUE, "depth", TRUE
          )
        }
        if (!is.null(context$prepared)) {
          times <- as.numeric(
            context$prepared$values[context$vertices$track_id ==
                                      id]
          )
          row[[date_col]] <- as.POSIXct(
            .track_interpolate(
              times,
              segment, fraction, context$options$interpolate_time,
              "time", TRUE
            ),
            origin = "1970-01-01", tz = "UTC"
          )
        }
        for (column in columns) {
          name <- paste0(column, "_match_info")
          info <- row[[name]][[1L]]
          if (identical(info$status, "selected")) {
            matching_sources <- Filter(
              function(m) {
                identical(m$file, info$source_file) && identical(
                  m$variable,
                  info$variable
                )
              }, context$metadata
            )
            if (length(matching_sources)) {
              m <- matching_sources[[1L]]
              horizontal <- .horizontal_match_row(.match_horizontal_cells(
                m$lon_axis,
                m$lat_axis, m$longitude, m$latitude, row[[lon_col]],
                row[[lat_col]], "within", NULL
              ), 1L)
              info$distance_km <- horizontal$distance_km
            }
            if (!is.null(context$prepared)) {
              info$time_displacement_days <- abs(
                as.numeric(
                  difftime(
                    .time_comparison_values(
                      info$timestamp_utc,
                      context$options$resolution
                    ), .time_comparison_values(
                      row[[date_col]],
                      context$options$resolution
                    ),
                    units = "days"
                  )
                )
              )
            }
            row[[name]][[1L]] <- info
          }
        }
        pieces[[length(pieces) + 1L]] <- row
      }
    }
  }
  observations <- samples[samples$sample_kind == "observation", , drop = FALSE]
  result <- do.call(rbind, c(pieces, list(observations)))
  result <- result[order(
    match(result$track_id, unique(samples$track_id)),
    result$distance_km
  ), , drop = FALSE]
  result$sample_id <- stats::ave(
    seq_len(nrow(result)), result$track_id, FUN = seq_along
  )
  rownames(result) <- NULL
  result
}


# Public function
# --------------------------------------------------------------

#' Extract environmental values along tracks and lines
#'
#' Sample ordered trajectories and reuse [extract_to_point()] for environmental
#' selection. Return sampled positions and optional summaries per track.
#'
#' @param data Ordered vertex table, or sf LINESTRING features with a CRS.
#'   Tables may be data frames, matrices or named lists. Rows are traversed in
#'   their supplied order within each track; they are never sorted by date.
#'   sf attributes may be scalar per line or list columns with one value per
#'   vertex.
#' @param nc,var netCDF sources and selected variable names. See
#'   [extract_to_buffer()].
#' @param track_col Optional track identifier column. Without it, a table is one
#'   track and each sf feature is a separate track. IDs are returned as text.
#' @param sampling 'vertices' (default) uses observations only; 'balanced'
#'   uses a shared subdivision satisfying all explicit step limits;
#'   'merged_steps' combines distance/depth cuts with temporal cuts;
#'   'cells' splits each connection at
#'   grid boundaries and samples interval midpoints, retaining original
#'   vertices.
#'   Repeated visits are retained. Cells mode currently uses spatial within
#'   matching.
#' @param spacing_km Maximum horizontal spacing in kilometres when inserting
#'   samples.
#'   Supply spacing or another step explicitly for balanced or merged_steps.
#' @param interpolate_depth,interpolate_time Allow linear interpolation of depth
#'   and elapsed UTC time between vertices. FALSE by default. Constant endpoint
#'   values are retained; missing endpoints produce missing intermediate values.
#' @param depth_step_m,time_step_seconds Optional maximum changes between
#'   samples. Require the corresponding interpolation switch. They supplement
#'   horizontal spacing, including for stationary or vertically moving tracks.
#' @param time_step,time_unit Optional temporal step and units: seconds,
#'   minutes,
#'   hours, elapsed days, whole calendar months or years. Balanced sampling uses
#'   one equal subdivision per original segment for elapsed-time limits.
#'   Calendar steps advance from each preceding step in UTC, clamping to the
#'   destination month's end. Balanced sampling preserves these calendar anchors
#'   and subdivides each calendar span to satisfy distance/depth limits.
#' @param distance_model WGS84 ellipsoidal geodesics (default), or sphere.
#' @param land_mask Optional sf polygon or single-layer terra SpatRaster mask.
#'   With no mask, marine connectivity is not checked.
#'   Geographic rasters must use longitudes from -180 to 180; rotate or
#'   reproject other conventions explicitly. Antimeridian connections are split.
#' @param marine_values Numeric raster codes representing marine cells,
#'   including the underlying codes of categorical rasters.
#' @param mask_type Meaning of polygon interiors: land or marine.
#' @param mask_domain Optional polygon domain for polygon masks. A land
#'   polygon alone cannot establish that its exterior is known sea.
#'   For rasters, coverage is defined by their extent and missing cells.
#' @param land_policy Skip invalid connections (default) or raise an error.
#'   Original observations remain; invalid support stays in the denominator.
#' @param mask_step_km Maximum segment length of the polyline used for mask
#'   predicates. This is an approximation to the geodesic; coastline accuracy
#'   also depends on the supplied mask. Touching land counts as an intersection.
#' @param method Depth method passed to [extract_to_point()], including 'all'.
#'   NULL uses nearest for variables with depth and 2d for variables without it.
#' @param statistics NULL for samples only, or distinct names from mean, sd, se,
#'   min and max. Summaries describe sampled environmental selections.
#' @param weighting 'distance', 'duration' or 'equal'. Distance and duration
#'   assign half each adjacent segment's horizontal length or elapsed seconds
#'   to its endpoint samples. Equal gives each sample the same weight.
#' @param min_valid_fraction Minimum valid fraction of total requested track
#'   support for a summary, from zero to one. Missing values stay in its
#'   denominator.
#' @param lon_col,lat_col,depth_col,date_col Vertex coordinate columns. Dates
#'   are
#'   required for selected temporal variables, duration weights or time-gap
#'   limits.
#' @param output_col Output name or named variable-to-column mapping.
#' @param time_resolution Point-time precision; see [extract_to_point()]. With
#'   auto, the original observations determine precision, before interpolation.
#' @param max_samples Maximum total sampled positions, a positive integer.
#' @param max_segment_km,max_gap_hours Optional limits that reject long segments
#'   or known elapsed time gaps. Split disconnected trajectories explicitly.
#' @param z_as_depth Interpret sf XYZ geometry as metres positive downwards.
#'   FALSE by default because geometric elevation may use a different
#'   convention.
#' @param quiet Hide explanatory reports; scientific and time-zone warnings
#'   remain.
#' @param ... Point-selection options: file_col, dimension overrides,
#'   time_match,
#'   max_time_diff, custom_time_interval, spatial_match, max_distance_km, and
#'   time_report, spatial_report, depth_report, verbose.
#'
#' @returns A list with samples (coordinates, original vertex attributes,
#'   source-row references, sample support and extracted values), summary (one
#'   row
#'   per track and output, or NULL), and settings. Tracks follow first
#'   appearance order. Summary counts n_observed and n_assumed distinguish
#'   original vertices from inserted positions; n_samples is their total,
#'   including missing results and invalid connections. The same counts are
#'   available once per track in settings$sampling_counts, even when
#'   statistics = NULL. Cell-interval midpoints count as assumed positions.
#'   variation_status is not_requested when neither SD nor SE is requested,
#'   available when the requested variation statistics are finite, or
#'   insufficient_samples when accepted support has fewer than two valid
#'   samples. not_available covers rejected support or non-finite statistics;
#'   consult status for coverage and support problems.
#'   Samples follow traversal order. For vertices-only sampling the global
#'   input row order is retained. Non-coordinate attributes of inserted samples
#'   are NA. source_row indexes the expanded vertex table for sf lines.
#'
#' @section Interpretation and limits:
#' Positions follow WGS84 shortest geodesics by default. Distances are
#' horizontal, not three-dimensional travel distances. Antipodal endpoints need
#' a waypoint. Interpolated coordinates describe assumptions, not observations.
#'
#' Balanced sampling approximates line integration. Cells mode finds numerical
#' intersections with rectilinear spatial cell boundaries and selection changes
#' in depth and time. It retains separate segment visits, including returns to
#' the same cell. Zero-support original vertices remain available for
#'   inspection.
#' Cells mode is not a volume-weighted integration of physical averaging
#'   periods;
#' the temporal rules still compare timestamp labels. Custom temporal precision
#' and spatial nearest matching are currently unsupported in cells mode.
#' Overlapping sources with different horizontal grids require explicit
#' alignment or another sampling mode.
#' Environmental values are selected using extract_to_point(), never
#'   interpolated.
#' Missing coordinates are rejected and multipart lines must be split
#'   explicitly.
#'
#' SD and SE use the same reliability-weight correction and effective sample
#' size as extract_to_buffer(). Adjacent samples are usually correlated: SE is
#' descriptive under an independence assumption, not an uncertainty estimate.
#' Densifying a track does not create independent environmental observations.
#' Distance-weighted summaries of stationary tracks are NA; duration weighting
#' requires known time support and returns NA if any adjacent duration is
#'   unknown.
#' settings records options and candidate files. Actual selected records are
#' stored in the per-output match_info columns in samples.
#'
#' @examples
#' nc_file <- system.file("extdata", "example_2d.nc", package = "ocean3d")
#' track <- data.frame(
#'   lon = c(0, 1), lat = c(40, 40),
#'   date = as.Date(c("2020-01-01", "2020-01-01"))
#' )
#' result <- extract_to_track(track, nc_file,
#'   var = "temp", method = "2d",
#'   spacing_km = 20, quiet = TRUE
#' )
#' result$summary
#' @seealso [extract_to_point()], [extract_to_buffer()]
#' @export
extract_to_track <- function(
  data, nc, var = NULL, track_col = NULL,
  sampling = c("vertices", "balanced", "merged_steps", "cells"),
  spacing_km = NULL,
  interpolate_depth = FALSE, interpolate_time = FALSE, depth_step_m = NULL,
  time_step_seconds = NULL, method = NULL, statistics = "mean", weighting = c(
    "distance",
    "duration", "equal"
  ), min_valid_fraction = 0.8, lon_col = "lon",
  lat_col = "lat", depth_col = "depth", date_col = "date", output_col = NULL,
  time_resolution = c(
    "auto", "datetime", "daily", "monthly", "yearly",
    "custom"
  ), max_samples = 100000L, max_segment_km = NULL, max_gap_hours = NULL,
  z_as_depth = FALSE, quiet = FALSE, time_step = NULL, time_unit = c(
    "hours",
    "minutes", "seconds", "days", "months", "years"
  ), distance_model = c(
    "WGS84",
    "sphere"
  ), land_mask = NULL, marine_values = NULL, mask_type = c(
    "land",
    "marine"
  ), mask_domain = NULL, land_policy = c("skip", "error"),
  mask_step_km = 0.5, ...
) {
  call <- match.call()
  call[[1L]] <- quote(.prepare_track)
  context <- eval(call, list(.prepare_track = .prepare_track), parent.frame())
  samples <- context$samples
  metadata <- context$metadata
  prepared <- context$prepared
  options <- context$options
  sampling <- options$sampling
  weighting <- match.arg(weighting)
  time_resolution <- match.arg(time_resolution)
  dots <- list(...)
  variables <- unique(vapply(metadata, function(x) x$variable, character(1)))
  resolution <- time_resolution
  if (time_resolution == "auto" && !is.null(prepared)) {
    resolution <- prepared$resolution
  }
  prefixes <- .buffer_output_names(
    variables, output_col, names(samples),
    FALSE
  )
  output_columns <- character()
  for (variable in variables) {
    group <- Filter(
      function(x) {
        x$variable == variable
      }, metadata
    )
    selected_method <- group[[1L]]$method
    columns <- .output_column_names(selected_method, variable,
                                    prefixes[[variable]])
    if (any(columns %in% names(samples))) {
      stop(
        paste0(
          "Extraction output would overwrite an existing column; ",
          "change `output_col`."
        ),
        call. = FALSE
      )
    }
    selected_sources <- lapply(group, function(x) x$source$source)
    arguments <- c(
      list(
        data = samples, nc = selected_sources, var = variable,
        method = selected_method, lon_col = lon_col, lat_col = lat_col,
        depth_col = depth_col, date_col = date_col,
        output_col = prefixes[[variable]],
        return_match_info = TRUE, time_resolution = resolution
      ), dots
    )
    if (quiet) {
      arguments$time_report <- FALSE
      arguments$spatial_report <- FALSE
      arguments$depth_report <- FALSE
    }
    samples <- do.call(extract_to_point, arguments)
    output_columns <- c(output_columns, columns)
  }
  if (!is.null(prepared)) {
    .warn_time_zone(prepared, paste0("`data$", date_col, "`"))
  }
  if (sampling == "cells") {
    samples <- .track_collapse_cells(
      samples, output_columns, context,
      lon_col, lat_col, depth_col, date_col
    )
  }
  invalid <- samples$is_interpolated & !samples$segment_status %in%
    c("marine", "not_checked")
  for (column in output_columns) {
    samples[[column]][invalid] <- NA_real_
    name <- paste0(column, "_match_info")
    for (i in which(invalid)) {
      candidate <- samples[[name]][[i]]
      samples[[name]][[i]] <- list(
        status = "invalid_connection", source_file = NA_character_,
        track_status = samples$segment_status[i],
        environmental_candidate = candidate
      )
    }
  }
  summary <- .track_summaries(
    samples, output_columns, statistics,
    weighting, min_valid_fraction
  )
  if (!quiet && !is.null(summary) && any(summary$status != "complete")) {
    message(
      "Track report: ", sum(summary$status != "complete"),
      paste0(
        " track/output summary row(s) need attention. See ",
        "`result$summary` for "
      ), paste0(
        "coverage and status; NA samples remain part of the ",
        "requested support."
      )
    )
  }
  list(
    samples = samples, summary = summary, settings = c(
      options[setdiff(
        names(options),
        c("metadata", "land_mask", "mask_domain")
      )], list(
        weighting = weighting,
        sampling_counts = .track_sampling_counts(samples),
        statistics = statistics, min_valid_fraction = min_valid_fraction,
        time_resolution = if (time_resolution == "auto" && !is.null(prepared)) {
          prepared$resolution
        } else {
          time_resolution
        }, candidate_sources = unique(
          vapply(
            metadata, function(x) x$file,
            character(1)
          )
        ),
        integration = if (sampling == "cells") "Grid interval support" else
          "Endpoint half-interval support (approximate)"
      )
    )
  )
}


.prepare_track <- function(
  data, nc, var = NULL, track_col = NULL, sampling = c(
    "vertices",
    "balanced", "merged_steps", "cells"
  ), spacing_km = NULL, interpolate_depth = FALSE,
  interpolate_time = FALSE, depth_step_m = NULL, time_step_seconds = NULL,
  method = NULL, statistics = "mean", weighting = c(
    "distance", "duration",
    "equal"
  ), min_valid_fraction = 0.8, lon_col = "lon", lat_col = "lat",
  depth_col = "depth", date_col = "date", output_col = NULL,
  time_resolution = c(
    "auto",
    "datetime", "daily", "monthly", "yearly", "custom"
  ), max_samples = 100000L,
  max_segment_km = NULL, max_gap_hours = NULL, z_as_depth = FALSE,
  quiet = FALSE, time_step = NULL, time_unit = c(
    "hours", "minutes",
    "seconds", "days", "months", "years"
  ), distance_model = c(
    "WGS84",
    "sphere"
  ), land_mask = NULL, marine_values = NULL, mask_type = c(
    "land",
    "marine"
  ), mask_domain = NULL, land_policy = c("skip", "error"),
  mask_step_km = 0.5, ...
) {
  if (
      missing(sampling) && any(
        !vapply(
          list(
            spacing_km, depth_step_m,
            time_step_seconds, time_step
          ), is.null, logical(1)
        )
      )) {
    sampling <- "balanced"
  }
  sampling <- match.arg(sampling)
  time_unit <- match.arg(time_unit)
  time_resolution <- match.arg(time_resolution)
  spatial_policy <- list(...)$spatial_match
  if (!is.null(spatial_policy)) {
    spatial_policy <- match.arg(spatial_policy, c("within", "nearest"))
  }
  distance_model <- match.arg(distance_model)
  mask_type <- match.arg(mask_type)
  land_policy <- match.arg(land_policy)
  .track_positive(mask_step_km, "mask_step_km", FALSE)
  .track_positive(time_step, "time_step")
  if (!is.null(time_step_seconds)) {
    if (!is.null(time_step)) {
      stop("Use time_step or time_step_seconds, not both.")
    }
    time_step <- time_step_seconds
    time_unit <- "seconds"
  }
  if (
      time_unit %in% c("months", "years") && !is.null(time_step) &&
        time_step != floor(time_step)) {
    stop("Calendar steps require whole months or years.")
  }
  if (!is.null(time_step) && !interpolate_time) {
    stop("time_step requires interpolate_time.")
  }
  if (
      sampling %in% c("balanced", "merged_steps") && all(
        vapply(
          list(
            spacing_km, depth_step_m,
            time_step
          ), is.null, logical(1)
        )
      )) {
    stop("Balanced or merged_steps sampling needs an explicit spacing or step.")
  }
  if (
      sampling == "cells" && (
                              identical(spatial_policy, "nearest") ||
                                identical(time_resolution, "custom"))) {
    stop(
      paste0(
        "Cells mode requires spatial within matching and ",
        "standard time precision."
      )
    )
  }
  weighting <- match.arg(weighting)
  time_resolution <- match.arg(time_resolution)
  if (!is.null(method)) {
    method <- match.arg(
      method, c(
        "2d", "nearest", "surface", "bottom",
        "all"
      )
    )
  }
  for (
       name in c(
         "interpolate_depth", "interpolate_time", "z_as_depth",
         "quiet"
       )) {
    .track_flag(get(name), name)
  }
  for (
       name in c(
         "spacing_km", "depth_step_m", "time_step_seconds",
         "max_segment_km", "max_gap_hours"
       )) {
    .track_positive(get(name), name)
  }
  .track_positive(max_samples, "max_samples", optional = FALSE)
  if (max_samples != floor(max_samples)) {
    stop("`max_samples` must be an integer.", call. = FALSE)
  }
  if (
      !is.numeric(min_valid_fraction) || length(min_valid_fraction) !=
        1L || !is.finite(min_valid_fraction) || min_valid_fraction <
        0 || min_valid_fraction > 1) {
    stop("`min_valid_fraction` must be between zero and one.", call. = FALSE)
  }
  supported_statistics <- c("mean", "sd", "se", "min", "max")
  bad_statistics <- !is.character(statistics) || !length(statistics) ||
    anyNA(statistics) || anyDuplicated(statistics) ||
    any(!statistics %in% supported_statistics)
  if (!is.null(statistics) && bad_statistics) {
    stop(
      paste0(
        "`statistics` must be NULL or distinct supported ",
        "statistic names."
      ),
      call. = FALSE
    )
  }
  if (
      !is.null(var) && (
                        !is.character(var) || !length(var) || anyNA(var) ||
                          any(!nzchar(var)) || anyDuplicated(var))) {
    stop("`var` must be NULL or distinct variable names.", call. = FALSE)
  }
  for (name in c("lon_col", "lat_col", "depth_col", "date_col", "track_col")) {
    value <- get(name)
    if (is.null(value) && name == "track_col") {
      next
    }
    if (
        !is.character(value) || length(value) != 1L || is.na(value) ||
          !nzchar(value)) {
      stop("`", name, "` must name one column.", call. = FALSE)
    }
  }
  coordinates <- c(lon_col, lat_col, depth_col, date_col)
  reserved <- c(
    "source_row", "segment_id", "segment_fraction", "is_interpolated",
    "distance_km", "sample_id", "support_km", "support_seconds",
    "available_km", "available_seconds", "segment_status",
    "segment_start_fraction",
    "segment_end_fraction", "sample_kind"
  )
  if (
      anyDuplicated(coordinates) || any(
        coordinates %in% c(
          reserved,
          "track_id"
        )
      ) || (!is.null(track_col) && track_col %in% coordinates)) {
    stop(
      "Coordinate, track and generated column names must be distinct.",
      call. = FALSE
    )
  }
  input_names <- if (is.matrix(data)) {
    colnames(data)
  } else {
    names(data)
  }
  if (anyDuplicated(input_names)) {
    stop("Input column names must be unique.", call. = FALSE)
  }
  if (
      any(reserved %in% input_names) || (
                                         "track_id" %in% input_names &&
                                           !identical(track_col, "track_id"))) {
    stop(
      paste0(
        "Input contains reserved track output columns. Rename ",
        "them first."
      ),
      call. = FALSE
    )
  }
  invalid_depth_step <- !is.null(depth_step_m) && !interpolate_depth
  invalid_time_step <- !is.null(time_step_seconds) && !interpolate_time
  if (invalid_depth_step || invalid_time_step) {
    stop(
      "Step limits require the corresponding interpolation switch.",
      call. = FALSE
    )
  }
  steps <- list(spacing_km, depth_step_m, time_step_seconds, time_step)
  has_steps <- any(!vapply(steps, is.null, logical(1)))
  if (sampling == "vertices" && has_steps) {
    stop(
      "Step limits cannot be used with vertices sampling.",
      call. = FALSE
    )
  }
  dots <- list(...)
  if (!is.null(dots$time_match)) {
    dots$time_match <- match.arg(dots$time_match, c("nearest", "exact"))
  }
  allowed <- c(
    "file_col", "lon_dim", "lat_dim", "depth_dim", "time_dim",
    "time_match", "max_time_diff", "custom_time_interval", "spatial_match",
    "max_distance_km", "time_report", "spatial_report", "depth_report",
    "verbose"
  )
  bad_dots <- is.null(names(dots)) || any(!nzchar(names(dots))) ||
    anyDuplicated(names(dots)) || any(!names(dots) %in% allowed)
  if (length(dots) && bad_dots) {
    stop(
      "Use named supported point-selection options in `...`.",
      call. = FALSE
    )
  }
  dimension_args <- dots[intersect(
    names(dots), c(
      "lon_dim", "lat_dim",
      "depth_dim", "time_dim"
    )
  )]
  metadata <- .track_sources(
    nc, var, method, if (is.null(dots$file_col)) {
      "file"
    } else {
      dots$file_col
    }, dimension_args
  )
  if (sampling == "cells") .track_check_cell_sources(metadata)
  time_used <- any(vapply(metadata, function(x) x$temporal, logical(1)))
  depth_used <- any(
    vapply(
      metadata, function(x) {
        x$method %in% c("nearest", "all")
      }, logical(1)
    )
  )
  vertices <- .track_input(
    data, track_col, lon_col, lat_col, depth_col,
    z_as_depth
  )
  if (depth_used || interpolate_depth) {
    .check_required_cols(vertices, depth_col)
    if (
        !is.numeric(vertices[[depth_col]]) || any(
          !is.na(vertices[[depth_col]]) &
            (
             !is.finite(vertices[[depth_col]]) | vertices[[depth_col]] <
               0)
        )) {
      stop("Depths must be non-negative metres or NA.", call. = FALSE)
    }
  }
  prepared <- NULL
  time_values <- NULL
  parse_time <- time_used || interpolate_time || weighting == "duration" ||
    !is.null(max_gap_hours)
  if (parse_time) {
    .check_required_cols(vertices, date_col)
    prepared <- .prepare_time_input(
      vertices[[date_col]], if (
                                time_resolution ==
                                  "custom") {
        "datetime"
      } else {
        time_resolution
      },
      input_name = paste0("`data$", date_col, "`"), ids = vertices$track_id
    )
    time_values <- as.numeric(prepared$values)
  }

  if (inherits(land_mask, "SpatRaster")) {
    if (!is.null(mask_domain)) {
      stop(
        "mask_domain applies to polygon masks; use raster NA cells.",
        call. = FALSE
      )
    }
    if (terra::nlyr(land_mask) != 1L || !nzchar(terra::crs(land_mask))) {
      stop("A marine raster needs one layer and a CRS.", call. = FALSE)
    }
    if (
        is.null(marine_values) || !is.numeric(marine_values) ||
          !length(marine_values) ||
          any(!is.finite(marine_values))) {
      stop(
        "Specify marine_values explicitly for a raster mask.",
        call. = FALSE
      )
    }
    outside_longitudes <- terra::xmin(land_mask) < -180 - 1e-8 ||
      terra::xmax(land_mask) > 180 + 1e-8
    if (terra::is.lonlat(land_mask) && outside_longitudes) {
      stop(
        "Rotate or reproject geographic masks to -180..180 longitudes.",
        call. = FALSE
      )
    }
  } else if (!is.null(land_mask)) {
    polygon_types <- c("POLYGON", "MULTIPOLYGON")
    if (
        !inherits(land_mask, c("sf", "sfc")) || is.na(sf::st_crs(land_mask)) ||
          any(!sf::st_is_valid(land_mask)) ||
          any(!sf::st_geometry_type(land_mask) %in% polygon_types)) {
      stop("A polygon mask needs valid polygons with a CRS.", call. = FALSE)
    }
    land_mask <- sf::st_union(sf::st_geometry(land_mask))
    if (!is.null(mask_domain)) {
      if (!inherits(mask_domain, c("sf", "sfc")) ||
            is.na(sf::st_crs(mask_domain)) ||
            any(!sf::st_geometry_type(mask_domain) %in% polygon_types) ||
            !all(sf::st_is_valid(mask_domain) %in% TRUE)) {
        stop("mask_domain needs polygons with a CRS.", call. = FALSE)
      }
      mask_domain <- sf::st_union(sf::st_geometry(mask_domain))
    }
  }
  options <- list(
    sampling = sampling, spacing_km = spacing_km,
    interpolate_depth = interpolate_depth,
    interpolate_time = interpolate_time, depth_step_m = depth_step_m,
    time_step_seconds = time_step_seconds, time_step = time_step,
    time_unit = time_unit, distance_model = distance_model,
    land_mask = land_mask,
    marine_values = marine_values, mask_type = mask_type,
    mask_domain = mask_domain,
    land_policy = land_policy, mask_step_km = mask_step_km, metadata = metadata,
    resolution = if (is.null(prepared)) {
      time_resolution
    } else {
      prepared$resolution
    },
    matching = dots, max_segment_km = max_segment_km,
    max_gap_hours = max_gap_hours
  )
  pieces <- list()
  remaining <- max_samples
  for (id in unique(vertices$track_id)) {
    piece <- .track_sample_one(
      vertices, which(
        vertices$track_id ==
          id
      ), options, lon_col, lat_col, depth_col, date_col, time_values,
      depth_used, time_used || weighting == "duration", remaining
    )
    remaining <- remaining - nrow(piece)
    pieces[[length(pieces) + 1L]] <- piece
  }
  samples <- do.call(rbind, pieces)
  if (sampling == "vertices") {
    samples <- samples[order(samples$source_row), , drop = FALSE]
  }
  rownames(samples) <- NULL
  list(
    samples = samples, metadata = metadata, prepared = prepared,
    options = options, vertices = vertices
  )
}
