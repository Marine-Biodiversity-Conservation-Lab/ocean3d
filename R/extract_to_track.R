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


.track_input <- function(data, track_col, lon_col, lat_col, depth_col,
                         z_as_depth) {
  if (inherits(data, "sf")) {
    if (is.na(sf::st_crs(data))) {
      stop("Track geometries need a CRS.", call. = FALSE)
    }
    types <- as.character(sf::st_geometry_type(data))
    if (any(types != "LINESTRING") || any(sf::st_is_empty(data))) {
      stop(
        "Use non-empty LINESTRING features; split multipart lines explicitly.",
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
    coordinates <- sf::st_coordinates(lines)
    if ("M" %in% colnames(coordinates)) {
      stop(
        "Measured (M) lines are unsupported; supply time explicitly.",
        call. = FALSE
      )
    }
    n <- tabulate(coordinates[, "L1"], length(lines))
    if (any(n < 2L)) {
      stop("Each track needs at least two vertices.", call. = FALSE)
    }
    result <- attributes[rep(seq_along(lines), n), , drop = FALSE]
    for (name in names(attributes)) {
      if (is.list(attributes[[name]])) {
        values <- attributes[[name]]
        bad <- which(lengths(values) != 1L & lengths(values) != n)
        if (length(bad)) {
          stop(
            "List column `", name, "` must have one value or one per ",
            "vertex in feature ", bad[1L], ".",
            call. = FALSE
          )
        }
        result[[name]] <- do.call(
          c, unname(Map(function(x, k) rep(x, length.out = k), values, n))
        )
      }
    }
    if ("Z" %in% colnames(coordinates)) {
      if (!z_as_depth) {
        stop(
          "XYZ lines require `z_as_depth = TRUE`: Z must be metres positive ",
          "down.",
          call. = FALSE
        )
      }
      if (depth_col %in% names(attributes)) {
        stop(
          "Supply depth either as Z or in `depth_col`, not both.",
          call. = FALSE
        )
      }
      result[[depth_col]] <- unname(coordinates[, "Z"])
    }
    result[[lon_col]] <- unname(coordinates[, "X"])
    result[[lat_col]] <- unname(coordinates[, "Y"])
    result$track_id <- rep(ids, n)
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


.track_interpolate <- function(x, i, fractions, enabled, name, used) {
  a <- x[i]
  b <- x[i + 1L]
  changes <- !is.na(a) && !is.na(b) && a != b
  if (any(fractions > 0 & fractions < 1) && used && !enabled && changes) {
    stop(
      "`", name, "` changes between vertices. Enable its interpolation ",
      "explicitly or use sampling = 'vertices'.",
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


# Mark inserted samples on connections that failed the marine-mask check.
# Their values become NA; each match_info record keeps its environmental
# selection, moves its status to `environmental_status`, and gains
# `status = "invalid_connection"` and the connection's `track_status`.
.track_mask_invalid <- function(samples, columns) {
  invalid <- samples$is_interpolated & !samples$segment_status %in%
    c("marine", "not_checked")
  for (column in columns) {
    samples[[column]][invalid] <- NA_real_
    name <- paste0(column, "_match_info")
    for (i in which(invalid)) {
      info <- samples[[name]][[i]]
      info$environmental_status <- info$status
      info$status <- "invalid_connection"
      info$track_status <- samples$segment_status[i]
      samples[[name]][[i]] <- info
    }
  }
  samples
}


# Settings shared by extract_to_track() and extract_to_buffered_track(),
# followed by the caller's own entries in `...`.
.track_settings <- function(context, samples, ..., statistics) {
  options <- context$options
  c(
    options[setdiff(names(options), c("metadata", "land_mask", "mask_domain"))],
    list(
      sampling_counts = .track_sampling_counts(samples),
      statistics = statistics,
      candidate_sources = unique(
        vapply(context$metadata, function(x) x$file, character(1))
      )
    ),
    list(...)
  )
}


# Depth and UTC time of positions at `fractions` along `segment`, following
# the interpolation permissions in `options`. `track` is one element of the
# per-track list returned by .track_sample_one().
.track_interpolate_row <- function(row, track, segment, fractions, options,
                                   depth_col, date_col) {
  if (track$depth_used || options$interpolate_depth) {
    row[[depth_col]] <- .track_interpolate(
      track$depths, segment, fractions, options$interpolate_depth,
      "interpolate_depth", track$depth_used
    )
  }
  if (!is.null(track$times)) {
    row[[date_col]] <- as.POSIXct(
      .track_interpolate(
        track$times, segment, fractions, options$interpolate_time,
        "interpolate_time", track$time_used
      ),
      origin = "1970-01-01", tz = "UTC"
    )
  }
  row
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


.track_summaries <- function(samples, columns, statistics, weighting,
                             threshold) {
  if (!length(statistics)) {
    return(NULL)
  }
  counts <- .track_sampling_counts(samples)
  rows <- list()
  for (id in unique(samples$track_id)) {
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
    usable <- is.finite(weights) & weights > 0
    total <- sum(weights[usable])
    # Missing duration support makes the full denominator unknowable.
    known <- !any(!is.finite(weights)) && total > 0
    support <- ifelse(is.finite(available) & available > 0, available, 0)
    support <- support / if (known) total else max(1, sum(support))
    track_counts <- counts[counts$track_id == id, , drop = FALSE]
    for (column in columns) {
      # Each sample's share of the track's support, in the buffer engine's
      # terms: available support stands in for the area fraction.
      frame <- data.frame(
        area_fraction = support, depth_fraction = 1, time_fraction = 1,
        value = selected[[column]]
      )
      summarise <- function(statistic) {
        result <- .buffer_summarise_samples(
          frame, "weighted", threshold, statistic
        )
        if (!known) {
          result$status <- "empty"
          result$valid_fraction <- NA_real_
          result$value <- NA_real_
        }
        result
      }
      overall <- summarise("mean")
      row <- data.frame(
        track_id = id, output_column = column, weighting = weighting,
        n_samples = track_counts$n_samples,
        n_observed = track_counts$n_observed,
        n_assumed = track_counts$n_assumed, n_valid = overall$n_valid,
        valid_fraction = overall$valid_fraction,
        total_distance_km = sum(selected$support_km),
        total_duration_seconds = if (anyNA(selected$support_seconds)) {
          NA_real_
        } else {
          sum(selected$support_seconds)
        },
        effective_n = overall$effective_n, status = overall$status,
        stringsAsFactors = FALSE
      )
      variation <- list()
      for (statistic in statistics) {
        result <- if (statistic == "mean") overall else summarise(statistic)
        row[[statistic]] <- result$value
        if (statistic %in% c("sd", "se")) {
          variation[[statistic]] <- result
        }
      }
      accepted <- overall$status %in% c("partial", "complete")
      row$variation_status <- if (!length(variation)) {
        "not_requested"
      } else if (any(vapply(variation, function(x) {
        x$status == "insufficient_variation_samples"
      }, logical(1)))) {
        "insufficient_variation_samples"
      } else if (accepted && all(is.finite(unlist(row[names(variation)])))) {
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
    stop("Antipodal endpoints need an intermediate waypoint.", call. = FALSE)
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
    stop(
      "A geodesic could not be constructed; add a waypoint.",
      call. = FALSE
    )
  }
  list(
    length_km = lengths, crs = crs, xy = xy,
    longitude = lon, latitude = lat
  )
}


.track_geodesic_positions <- function(arc, segment, fractions) {
  xy <- outer(fractions, as.numeric(arc$xy[[segment]]))
  result <- sf::sf_project(arc$crs[segment], "EPSG:4326", xy)
  # Retain supplied endpoint coordinates, including stationary
  # segments.
  result[fractions == 0, ] <- c(arc$longitude[segment], arc$latitude[segment])
  result[fractions == 1, ] <- c(
    arc$longitude[segment + 1L], arc$latitude[segment + 1L]
  )
  colnames(result) <- c("lon", "lat")
  result
}


.track_time_breaks <- function(start, end, step, unit, budget) {
  if (is.null(step) || !is.finite(start) || !is.finite(end) || end == start) {
    return(c(0, 1))
  }
  # .prepare_track() has validated `step`, including whole calendar steps.
  if (unit %in% c("months", "years")) {
    current <- as.POSIXct(start, origin = "1970-01-01", tz = "UTC")
    values <- start
    while (as.numeric(current) < end) {
      following <- .buffer_shift_time(current, step, unit)
      if (length(values) >= budget) {
        stop("time_step exceeds max_samples.")
      }
      values <- c(values, min(end, as.numeric(following)))
      current <- following
    }
    return((values - start) / (end - start))
  }
  seconds <- step * .time_unit_seconds(unit)
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


.track_grid_breaks <- function(arc, segment, metadata, depths, times,
                               resolution, matching, budget) {
  reference <- arc$longitude[segment]
  wrap <- function(lon) reference + (lon - reference + 180) %% 360 - 180
  coordinate <- function(f, axis) {
    value <- .track_geodesic_positions(arc, segment, f)[, axis]
    if (axis == 1L) wrap(value) else value
  }
  # Distinct values of `field(m)` across all source descriptions.
  axis_values <- function(field) {
    unique(unlist(lapply(metadata, field)))
  }
  cuts <- c(0, 1)
  # Latitude may turn within one geodesic; bracket each monotonic portion.
  latitude <- function(f) coordinate(f, 2L)
  turns <- c(
    0, 1,
    stats::optimize(latitude, c(0, 1), tol = 1e-12)$minimum,
    stats::optimize(latitude, c(0, 1), tol = 1e-12, maximum = TRUE)$maximum
  )
  turns <- sort(unique(round(turns, 12)))
  for (axis in 1:2) {
    brackets <- if (axis == 1L) c(0, 1) else turns
    extrema <- coordinate(brackets, axis)
    if (axis == 2L && any(abs(extrema) > 89.999999)) {
      stop("Cells mode needs waypoints away from a pole.", call. = FALSE)
    }
    edges <- axis_values(function(m) {
      a <- if (axis == 1L) m$lon_axis else m$lat_axis
      edges <- c(a$lower, a$upper)
      order <- order(a$centres)
      overlapping <- length(order) > 1L && any(
        tail(a$lower[order], -1L) < head(a$upper[order], -1L) - 1e-9
      )
      if (overlapping) {
        centres <- if (axis == 1L) wrap(a$centres) else a$centres
        edges <- c(edges, .track_midpoint_boundaries(
          centres, min(extrema), max(extrema), budget
        ))
      }
      edges
    })
    if (axis == 1L) {
      edges <- wrap(edges)
    }
    edges <- edges[
      is.finite(edges) & edges > min(extrema) & edges < max(extrema)
    ]
    if (length(edges) > budget) {
      stop("Grid intersections exceed max_samples.")
    }
    for (j in seq_len(length(brackets) - 1L)) {
      interval <- brackets[c(j, j + 1L)]
      at_ends <- coordinate(interval, axis)
      for (edge in edges[edges > min(at_ends) & edges < max(at_ends)]) {
        root <- stats::uniroot(
          function(f) coordinate(f, axis) - edge, interval, tol = 1e-12
        )$root
        cuts <- c(cuts, root)
      }
    }
  }
  if (all(is.finite(depths)) && diff(depths) != 0) {
    levels <- axis_values(function(m) {
      if (m$method %in% c("nearest", "all")) m$depth_m
    })
    edges <- .track_midpoint_boundaries(
      levels, min(depths), max(depths), budget
    )
    cuts <- c(cuts, (edges - depths[1L]) / diff(depths))
  }
  if (all(is.finite(times)) && diff(times) > 0) {
    if (resolution %in% c("daily", "monthly", "yearly")) {
      step <- switch(resolution,
        daily = "day", monthly = "month", yearly = "year"
      )
      first <- .time_comparison_values(
        as.POSIXct(times[1L], origin = "1970-01-01", tz = "UTC"), resolution
      )
      seconds_min <- switch(resolution,
        daily = 86400, monthly = 28 * 86400, yearly = 365 * 86400
      )
      if (diff(times) / seconds_min > budget) {
        stop("Time grid exceeds max_samples.")
      }
      last <- as.POSIXct(
        times[2L] + seconds_min, origin = "1970-01-01", tz = "UTC"
      )
      edges <- as.numeric(seq(first, last, by = step))
    } else {
      values <- axis_values(function(m) {
        v <- sort(m$time_values)
        c(
          v[v >= times[1L] & v <= times[2L]],
          tail(v[v < times[1L]], 1L),
          head(v[v > times[2L]], 1L)
        )
      })
      edges <- if (identical(matching$time_match, "exact")) {
        values
      } else {
        .track_midpoint_boundaries(values, times[1L], times[2L], budget)
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
  cuts <- sort(unique(c(0, 1, cuts[is.finite(cuts) & cuts > 0 & cuts < 1])))
  if (length(cuts) > budget + 1) {
    stop("Grid intersections exceed max_samples.")
  }
  cuts
}


.track_land_check <- function(
  arc, segment, mask, marine_values, mask_type,
  domain, step_km
) {
  # .prepare_track() validates `mask` and `domain` and unions polygons once.
  if (is.null(mask)) {
    return("not_checked")
  }
  # Mask predicates inspect the whole connecting polyline, not the samples.
  # It is built with the same geodesic positions as the samples.
  n <- max(1, ceiling(arc$length_km[segment] / step_km))
  if (n > 1e+05) {
    stop("Mask geometry exceeds the segment safety limit.")
  }
  xy <- .track_geodesic_positions(arc, segment, seq(0, 1, length.out = n + 1L))
  xy[, 1L] <- .normalise_longitude(xy[, 1L])
  line <- sf::st_sfc(sf::st_linestring(unname(xy)), crs = 4326)
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
  arc <- .track_geodesic(vertices[[lon_col]][rows], vertices[[lat_col]][rows])
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
  track <- list(
    arc = arc, depths = depths,
    times = if (is.null(time_values)) NULL else times,
    depth_used = depth_used, time_used = time_used
  )
  dt <- diff(times)
  if (any(dt < 0, na.rm = TRUE)) {
    stop("Dates must follow traversal order.")
  }
  if (!is.null(options$max_segment_km) &&
        any(arc$length_km > options$max_segment_km)) {
    stop("A segment exceeds max_segment_km; split the track or add waypoints.")
  }
  if (!is.null(options$max_gap_hours) &&
        any(dt > options$max_gap_hours * 3600, na.rm = TRUE)) {
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
                              any(!is.finite(depths[ends]))) {
              1
            } else {
              abs(diff(depths[ends])) * (last - first) / options$depth_step_m
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
        breaks <- .track_grid_breaks(
          arc, i, options$metadata, depths[ends], times[ends],
          options$resolution, options$matching, budget
        )
        cuts <- sort(unique(c(cuts, breaks)))
      }
    }
    a <- head(cuts, -1L)
    b <- tail(cuts, -1L)
    fractions <- if (options$sampling == "cells") (a + b) / 2 else a
    used <- sum(vapply(pieces, nrow, integer(1)))
    if (length(fractions) + used + 1 > budget) {
      stop("Track sampling exceeds max_samples.")
    }
    xy <- .track_geodesic_positions(arc, i, fractions)
    original <- fractions == 0
    row <- vertices[rep(rows[i], length(fractions)), , drop = FALSE]
    for (name in setdiff(names(row), "track_id")) row[[name]][!original] <- NA
    row[[lon_col]] <- xy[, 1L]
    row[[lat_col]] <- xy[, 2L]
    row <- .track_interpolate_row(
      row, track, i, fractions, options, depth_col, date_col
    )
    row$source_row <- ifelse(original, rows[i], NA_integer_)
    row$segment_id <- i
    row$segment_fraction <- fractions
    row$segment_start_fraction <- a
    row$segment_end_fraction <- b
    row$is_interpolated <- !original
    row$distance_km <- sum(head(arc$length_km, i - 1L)) +
      fractions * arc$length_km[i]
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
      last[[date_col]] <- as.POSIXct(times, origin = "1970-01-01", tz = "UTC")
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
    last[intersect(names(vertices), names(last))] <-
      vertices[tail(rows, 1L), , drop = FALSE]
    if (!is.null(time_values)) {
      last[[date_col]] <- as.POSIXct(
        tail(times, 1L), origin = "1970-01-01", tz = "UTC"
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
    for (name in c("support_km", "available_km", "support_seconds",
                   "available_seconds")) {
      intervals <- head(result[[name]], -1L)
      result[[name]] <- (c(0, intervals) + c(intervals, 0)) / 2
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
  list(samples = result, track = track)
}


.track_collapse_cells <- function(samples, columns, context, lon_col,
                                  lat_col, depth_col, date_col) {
  pieces <- list()
  for (id in unique(samples$track_id)) {
    rows <- which(samples$track_id == id)
    track <- context$tracks[[id]]
    arc <- track$arc
    for (segment in unique(samples$segment_id[rows])) {
      selected <- rows[samples$segment_id[rows] == segment &
                         samples$sample_kind[rows] == "cell_interval"]
      if (!length(selected)) {
        next
      }
      fields <- c(
        "status", "source_file", "lon_index", "lat_index", "depth_m",
        "time_index"
      )
      keys <- vapply(selected, function(i) {
        info <- lapply(columns, function(column) {
          samples[[paste0(column, "_match_info")]][[i]][fields]
        })
        values <- unname(unlist(samples[i, columns, drop = FALSE]))
        key <- list(info, samples$segment_status[i], values)
        paste(utils::capture.output(dput(key)), collapse = "")
      }, character(1))
      groups <- cumsum(c(TRUE, keys[-1L] != head(keys, -1L)))
      for (group in unique(groups)) {
        indices <- selected[groups == group]
        row <- samples[indices[1L], , drop = FALSE]
        for (name in c("support_km", "support_seconds", "available_km",
                       "available_seconds")) {
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
        row$distance_km <- sum(head(arc$length_km, segment - 1L)) +
          arc$length_km[segment] * fraction
        row <- .track_interpolate_row(
          row, track, segment, fraction, context$options, depth_col, date_col
        )
        for (column in columns) {
          name <- paste0(column, "_match_info")
          info <- row[[name]][[1L]]
          if (identical(info$status, "selected")) {
            matching_sources <- Filter(function(m) {
              identical(m$file, info$source_file) &&
                identical(m$variable, info$variable)
            }, context$metadata)
            if (length(matching_sources)) {
              m <- matching_sources[[1L]]
              cells <- .match_horizontal_cells(
                m$lon_axis, m$lat_axis, m$longitude, m$latitude,
                row[[lon_col]], row[[lat_col]], "within", NULL
              )
              info$distance_km <- .horizontal_match_row(cells, 1L)$distance_km
            }
            if (!is.null(context$prepared)) {
              resolution <- context$options$resolution
              info$time_displacement_days <- abs(as.numeric(difftime(
                .time_comparison_values(info$timestamp_utc, resolution),
                .time_comparison_values(row[[date_col]], resolution),
                units = "days"
              )))
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
#' @param sampling NULL (default) uses 'balanced' when spacing_km,
#'   depth_step_m, time_step_seconds or time_step is supplied, and
#'   'vertices' otherwise. 'vertices' uses observations only; 'balanced'
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
#'   insufficient_variation_samples when accepted support has fewer than two
#'   valid samples. not_available covers rejected support or non-finite
#'   statistics; consult status for coverage and support problems. status
#'   follows extract_to_buffer(): empty (no known support), all_missing,
#'   insufficient_coverage, partial or complete.
#'   Samples follow traversal order. For vertices-only sampling the global
#'   input row order is retained. Non-coordinate attributes of inserted samples
#'   are NA. source_row indexes the expanded vertex table for sf lines.
#'
#' @section Interpretation and limits:
#' Positions follow WGS84 shortest geodesics. Distances are
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
  data, nc, var = NULL, track_col = NULL, sampling = NULL, spacing_km = NULL,
  interpolate_depth = FALSE, interpolate_time = FALSE, depth_step_m = NULL,
  time_step_seconds = NULL, method = NULL, statistics = "mean",
  weighting = c("distance", "duration", "equal"), min_valid_fraction = 0.8,
  lon_col = "lon", lat_col = "lat", depth_col = "depth", date_col = "date",
  output_col = NULL,
  time_resolution = c(
    "auto", "datetime", "daily", "monthly", "yearly", "custom"
  ),
  max_samples = 100000L, max_segment_km = NULL, max_gap_hours = NULL,
  z_as_depth = FALSE, quiet = FALSE, time_step = NULL,
  time_unit = c("hours", "minutes", "seconds", "days", "months", "years"),
  land_mask = NULL, marine_values = NULL, mask_type = c("land", "marine"),
  mask_domain = NULL, land_policy = c("skip", "error"), mask_step_km = 0.5,
  ...
) {
  if (!is.null(sampling)) {
    sampling <- match.arg(
      sampling, c("vertices", "balanced", "merged_steps", "cells")
    )
  }
  weighting <- match.arg(weighting)
  time_resolution <- match.arg(time_resolution)
  time_unit <- match.arg(time_unit)
  mask_type <- match.arg(mask_type)
  land_policy <- match.arg(land_policy)
  dots <- list(...)
  track_args <- mget(setdiff(names(formals(extract_to_track)), "..."))
  context <- .prepare_track(track_args, dots)
  samples <- context$samples
  metadata <- context$metadata
  prepared <- context$prepared
  options <- context$options
  sampling <- options$sampling
  variables <- unique(vapply(metadata, function(x) x$variable, character(1)))
  resolution <- if (time_resolution == "auto" && !is.null(prepared)) {
    prepared$resolution
  } else {
    time_resolution
  }
  prefixes <- .buffer_output_names(variables, output_col, names(samples), FALSE)
  output_columns <- character()
  for (variable in variables) {
    group <- Filter(function(x) x$variable == variable, metadata)
    selected_method <- group[[1L]]$method
    columns <- .output_column_names(
      selected_method, variable, prefixes[[variable]]
    )
    if (any(columns %in% names(samples))) {
      stop(
        "Extraction output would overwrite an existing column; ",
        "change `output_col`.",
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
      samples, output_columns, context, lon_col, lat_col, depth_col, date_col
    )
  }
  samples <- .track_mask_invalid(samples, output_columns)
  summary <- .track_summaries(
    samples, output_columns, statistics, weighting, min_valid_fraction
  )
  if (!quiet && !is.null(summary) && any(summary$status != "complete")) {
    message(
      "Track report: ", sum(summary$status != "complete"),
      " track/output summary row(s) need attention. See `result$summary` ",
      "for coverage and status; NA samples remain part of the requested ",
      "support."
    )
  }
  list(
    samples = samples, summary = summary,
    settings = .track_settings(
      context, samples, statistics = statistics,
      weighting = weighting, min_valid_fraction = min_valid_fraction,
      time_resolution = resolution,
      integration = if (sampling == "cells") "Grid interval support" else
        "Endpoint half-interval support (approximate)"
    )
  )
}


# Validate the arguments of extract_to_track() or extract_to_buffered_track()
# (a named list, enumerations already matched) and sample the tracks.
# `sampling = NULL` chooses balanced sampling when a step is supplied and
# vertices otherwise. `dots` holds the point-selection options.
.prepare_track <- function(args, dots) {
  steps <- list(
    args$spacing_km, args$depth_step_m, args$time_step_seconds, args$time_step
  )
  has_steps <- any(!vapply(steps, is.null, logical(1)))
  if (is.null(args$sampling)) {
    args$sampling <- if (has_steps) "balanced" else "vertices"
  }
  spatial_policy <- dots$spatial_match
  if (!is.null(spatial_policy)) {
    spatial_policy <- match.arg(spatial_policy, c("within", "nearest"))
  }
  .track_positive(args$mask_step_km, "mask_step_km", FALSE)
  .track_positive(args$time_step, "time_step")
  if (!is.null(args$time_step_seconds)) {
    if (!is.null(args$time_step)) {
      stop("Use time_step or time_step_seconds, not both.")
    }
    args$time_step <- args$time_step_seconds
    args$time_unit <- "seconds"
  }
  if (args$time_unit %in% c("months", "years") && !is.null(args$time_step) &&
        args$time_step != floor(args$time_step)) {
    stop("Calendar steps require whole months or years.")
  }
  if (!is.null(args$time_step) && !args$interpolate_time) {
    stop("time_step requires interpolate_time.")
  }
  no_steps <- is.null(args$spacing_km) && is.null(args$depth_step_m) &&
    is.null(args$time_step)
  if (args$sampling %in% c("balanced", "merged_steps") && no_steps) {
    stop("Balanced or merged_steps sampling needs an explicit spacing or step.")
  }
  custom_matching <- identical(spatial_policy, "nearest") ||
    identical(args$time_resolution, "custom")
  if (args$sampling == "cells" && custom_matching) {
    stop(
      "Cells mode requires spatial within matching and standard time precision."
    )
  }
  if (!is.null(args$method)) {
    args$method <- match.arg(
      args$method, c("2d", "nearest", "surface", "bottom", "all")
    )
  }
  for (name in c("interpolate_depth", "interpolate_time", "z_as_depth",
                 "quiet")) {
    .track_flag(args[[name]], name)
  }
  for (name in c("spacing_km", "depth_step_m", "time_step_seconds",
                 "max_segment_km", "max_gap_hours")) {
    .track_positive(args[[name]], name)
  }
  .track_positive(args$max_samples, "max_samples", optional = FALSE)
  if (args$max_samples != floor(args$max_samples)) {
    stop("`max_samples` must be an integer.", call. = FALSE)
  }
  fraction <- args$min_valid_fraction
  if (!is.numeric(fraction) || length(fraction) != 1L ||
        !is.finite(fraction) || fraction < 0 || fraction > 1) {
    stop("`min_valid_fraction` must be between zero and one.", call. = FALSE)
  }
  supported_statistics <- c("mean", "sd", "se", "min", "max")
  bad_statistics <- !is.character(args$statistics) ||
    !length(args$statistics) || anyNA(args$statistics) ||
    anyDuplicated(args$statistics) ||
    any(!args$statistics %in% supported_statistics)
  if (!is.null(args$statistics) && bad_statistics) {
    stop(
      "`statistics` must be NULL or distinct supported statistic names.",
      call. = FALSE
    )
  }
  bad_var <- !is.character(args$var) || !length(args$var) ||
    anyNA(args$var) || any(!nzchar(args$var)) || anyDuplicated(args$var)
  if (!is.null(args$var) && bad_var) {
    stop("`var` must be NULL or distinct variable names.", call. = FALSE)
  }
  for (name in c("lon_col", "lat_col", "depth_col", "date_col", "track_col")) {
    value <- args[[name]]
    if (is.null(value) && name == "track_col") {
      next
    }
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
          !nzchar(value)) {
      stop("`", name, "` must name one column.", call. = FALSE)
    }
  }
  coordinates <- c(args$lon_col, args$lat_col, args$depth_col, args$date_col)
  reserved <- c(
    "source_row", "segment_id", "segment_fraction", "is_interpolated",
    "distance_km", "sample_id", "support_km", "support_seconds",
    "available_km", "available_seconds", "segment_status",
    "segment_start_fraction", "segment_end_fraction", "sample_kind"
  )
  if (anyDuplicated(coordinates) ||
        any(coordinates %in% c(reserved, "track_id")) ||
        (!is.null(args$track_col) && args$track_col %in% coordinates)) {
    stop(
      "Coordinate, track and generated column names must be distinct.",
      call. = FALSE
    )
  }
  input_names <- if (is.matrix(args$data)) {
    colnames(args$data)
  } else {
    names(args$data)
  }
  if (anyDuplicated(input_names)) {
    stop("Input column names must be unique.", call. = FALSE)
  }
  if (any(reserved %in% input_names) ||
        ("track_id" %in% input_names &&
           !identical(args$track_col, "track_id"))) {
    stop(
      "Input contains reserved track output columns. Rename them first.",
      call. = FALSE
    )
  }
  # time_step_seconds without interpolate_time is rejected with time_step.
  if (!is.null(args$depth_step_m) && !args$interpolate_depth) {
    stop(
      "Step limits require the corresponding interpolation switch.",
      call. = FALSE
    )
  }
  if (args$sampling == "vertices" && has_steps) {
    stop("Step limits cannot be used with vertices sampling.", call. = FALSE)
  }
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
  dimension_args <- dots[
    intersect(names(dots), c("lon_dim", "lat_dim", "depth_dim", "time_dim"))
  ]
  # Tracks keep point extraction's support for climatological time axes.
  described <- .describe_netcdf_sources(
    args$nc, args$var,
    if (is.null(dots$file_col)) "file" else dots$file_col, dimension_args,
    function(handle, variable, dims, i) {
      lon <- .get_dim_values(handle, dims$lon)
      lat <- .get_dim_values(handle, dims$lat)
      list(
        variable = variable, source_index = i, file = handle$filename,
        method = if (!is.null(args$method)) {
          args$method
        } else if (is.null(dims$depth)) {
          "2d"
        } else {
          "nearest"
        },
        temporal = !is.null(dims$time), longitude = lon, latitude = lat,
        lon_axis = .prepare_spatial_axis(
          handle, dims$lon, lon, longitude = TRUE
        ),
        lat_axis = .prepare_spatial_axis(handle, dims$lat, lat),
        depth_m = if (!is.null(dims$depth)) {
          .read_depth_axis(
            handle, dims$depth, .get_dim_values(handle, dims$depth)
          )$depth_m
        },
        time_values = if (is.null(dims$time)) {
          numeric()
        } else {
          as.numeric(.decode_netcdf_time(handle, dims$time))
        }
      )
    }
  )
  lapply(described$opened, .close_source)
  metadata <- lapply(described$metadata, function(m) {
    m$source <- described$sources[[m$source_index]]
    m
  })
  if (args$sampling == "cells") .track_check_cell_sources(metadata)
  time_used <- any(vapply(metadata, function(x) x$temporal, logical(1)))
  depth_used <- any(vapply(
    metadata, function(x) x$method %in% c("nearest", "all"), logical(1)
  ))
  vertices <- .track_input(
    args$data, args$track_col, args$lon_col, args$lat_col, args$depth_col,
    args$z_as_depth
  )
  if (depth_used || args$interpolate_depth) {
    .check_required_cols(vertices, args$depth_col)
    depths <- vertices[[args$depth_col]]
    if (!is.numeric(depths) ||
          any(!is.na(depths) & (!is.finite(depths) | depths < 0))) {
      stop("Depths must be non-negative metres or NA.", call. = FALSE)
    }
  }
  prepared <- NULL
  time_values <- NULL
  parse_time <- time_used || args$interpolate_time ||
    args$weighting == "duration" || !is.null(args$max_gap_hours)
  for (m in Filter(function(x) x$temporal, metadata)) {
    .require_observation_time(vertices, args$date_col, m$variable, m$file)
  }
  if (parse_time) {
    .check_required_cols(vertices, args$date_col)
    prepared <- .prepare_time_input(
      vertices[[args$date_col]],
      if (args$time_resolution == "custom") "datetime" else
        args$time_resolution,
      input_name = paste0("`data$", args$date_col, "`"),
      ids = vertices$track_id
    )
    time_values <- as.numeric(prepared$values)
  }
  if (inherits(args$land_mask, "SpatRaster")) {
    if (!is.null(args$mask_domain)) {
      stop(
        "mask_domain applies to polygon masks; use raster NA cells.",
        call. = FALSE
      )
    }
    if (terra::nlyr(args$land_mask) != 1L ||
          !nzchar(terra::crs(args$land_mask))) {
      stop("A marine raster needs one layer and a CRS.", call. = FALSE)
    }
    values <- args$marine_values
    if (!is.numeric(values) || !length(values) || any(!is.finite(values))) {
      stop("Specify marine_values explicitly for a raster mask.", call. = FALSE)
    }
    outside_longitudes <- terra::xmin(args$land_mask) < -180 - 1e-8 ||
      terra::xmax(args$land_mask) > 180 + 1e-8
    if (terra::is.lonlat(args$land_mask) && outside_longitudes) {
      stop(
        "Rotate or reproject geographic masks to -180..180 longitudes.",
        call. = FALSE
      )
    }
  } else if (!is.null(args$land_mask)) {
    polygon_types <- c("POLYGON", "MULTIPOLYGON")
    if (!inherits(args$land_mask, c("sf", "sfc")) ||
          is.na(sf::st_crs(args$land_mask)) ||
          any(!sf::st_is_valid(args$land_mask)) ||
          any(!sf::st_geometry_type(args$land_mask) %in% polygon_types)) {
      stop("A polygon mask needs valid polygons with a CRS.", call. = FALSE)
    }
    args$land_mask <- sf::st_union(sf::st_geometry(args$land_mask))
    if (!is.null(args$mask_domain)) {
      if (!inherits(args$mask_domain, c("sf", "sfc")) ||
            is.na(sf::st_crs(args$mask_domain)) ||
            any(!sf::st_geometry_type(args$mask_domain) %in% polygon_types) ||
            !all(sf::st_is_valid(args$mask_domain) %in% TRUE)) {
        stop("mask_domain needs polygons with a CRS.", call. = FALSE)
      }
      args$mask_domain <- sf::st_union(sf::st_geometry(args$mask_domain))
    }
  }
  options <- c(
    args[c(
      "sampling", "spacing_km", "interpolate_depth", "interpolate_time",
      "depth_step_m", "time_step_seconds", "time_step", "time_unit",
      "land_mask", "marine_values", "mask_type", "mask_domain",
      "land_policy", "mask_step_km"
    )],
    list(
      metadata = metadata,
      resolution = if (is.null(prepared)) args$time_resolution else
        prepared$resolution,
      matching = dots
    ),
    args[c("max_segment_km", "max_gap_hours")]
  )
  pieces <- list()
  tracks <- list()
  remaining <- args$max_samples
  for (id in unique(vertices$track_id)) {
    sampled <- .track_sample_one(
      vertices, which(vertices$track_id == id), options,
      args$lon_col, args$lat_col, args$depth_col, args$date_col, time_values,
      depth_used, time_used || args$weighting == "duration", remaining
    )
    remaining <- remaining - nrow(sampled$samples)
    pieces[[length(pieces) + 1L]] <- sampled$samples
    tracks[[id]] <- sampled$track
  }
  samples <- do.call(rbind, pieces)
  if (args$sampling == "vertices") {
    samples <- samples[order(samples$source_row), , drop = FALSE]
  }
  rownames(samples) <- NULL
  list(
    samples = samples, metadata = metadata, prepared = prepared,
    options = options, vertices = vertices, tracks = tracks
  )
}
