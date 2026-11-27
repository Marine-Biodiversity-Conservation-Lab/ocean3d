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


.track_input <- function(data, track_col, lon_col, lat_col) {
  # Same inputs as extract_to_point(): tables, matrices, named lists and sf
  # POINT features. Rows are vertices in traversal order.
  result <- .normalize_points_input(data, lon_col, lat_col)
  if (inherits(result, "sf")) {
    result <- sf::st_drop_geometry(result)
  }
  .check_required_cols(result, c(lon_col, lat_col, track_col))
  if (!nrow(result)) {
    stop("Supply at least one track.", call. = FALSE)
  }
  result$track_id <- if (is.null(track_col)) {
    "1"
  } else {
    as.character(result[[track_col]])
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


# The sampling mode, checked with its step limits and `max_samples`, for
# extract_to_track() and sample_3d(). NULL picks "balanced" when a step is
# given and "vertices" otherwise. Steps must be positive, calendar steps
# whole months or years, and `max_samples` a positive whole number;
# vertices sampling takes no steps, and the stepped modes need one. The one
# producer of `ocean3d_error_track_sampling`.
.track_sampling <- function(sampling, spacing_km, depth_step_m, time_step,
                            time_unit, max_samples) {
  fail <- function(message) {
    rlang::abort(message, class = "ocean3d_error_track_sampling", call = NULL)
  }
  positive <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0
  }
  for (name in c("spacing_km", "depth_step_m", "time_step")) {
    value <- get(name)
    if (!is.null(value) && !positive(value)) {
      fail(glue::glue("`{name}` must be one finite positive number."))
    }
  }
  if (!positive(max_samples) || max_samples != floor(max_samples)) {
    fail("`max_samples` must be a positive integer.")
  }
  has_steps <- !is.null(spacing_km) || !is.null(depth_step_m) ||
    !is.null(time_step)
  if (is.null(sampling)) sampling <- if (has_steps) "balanced" else "vertices"
  if (time_unit %in% c("months", "years") && !is.null(time_step) &&
        time_step != floor(time_step)) {
    fail("Calendar steps require whole months or years.")
  }
  if (sampling %in% c("balanced", "merged_steps") && !has_steps) {
    fail("Balanced or merged_steps sampling needs an explicit spacing or step.")
  }
  if (sampling == "vertices" && has_steps) {
    fail("Step limits cannot be used with vertices sampling.")
  }
  sampling
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
    # Break k is k steps from the start, not one step from the previous
    # break: end-of-month clamping (31 January + 1 month = 29 February)
    # would otherwise carry forward and every later break would drift (#146).
    first <- as.POSIXct(start, origin = "1970-01-01", tz = "UTC")
    values <- start
    k <- 1
    repeat {
      following <- as.numeric(.buffer_shift_time(first, k * step, unit))
      if (length(values) >= budget) {
        stop("time_step exceeds max_samples.")
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


# Every track of `vertices` sampled with .track_sample_one(), sharing one
# budget of `max_samples` points, for extract_to_track() and sample_3d().
# Returns the samples, in input row order for vertices sampling and in
# track order otherwise, and the per-track list .track_collapse_cells()
# reads.
.track_sample_all <- function(vertices, options, lon_col, lat_col, depth_col,
                              date_col, time_values, depth_used, time_used,
                              max_samples) {
  pieces <- list()
  tracks <- list()
  remaining <- max_samples
  for (id in unique(vertices$track_id)) {
    sampled <- .track_sample_one(
      vertices, which(vertices$track_id == id), options, lon_col, lat_col,
      depth_col, date_col, time_values, depth_used, time_used, remaining
    )
    remaining <- remaining - nrow(sampled$samples)
    pieces[[length(pieces) + 1L]] <- sampled$samples
    tracks[[id]] <- sampled$track
  }
  samples <- do.call(rbind, pieces)
  if (options$sampling == "vertices") {
    samples <- samples[order(samples$source_row), , drop = FALSE]
  }
  rownames(samples) <- NULL
  list(samples = samples, tracks = tracks)
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


# -----------------------------------------------------------------------------
# Track extraction: a convenience wrapper over the three steps (#116)
# -----------------------------------------------------------------------------
#
# The vertex table becomes SpatVector3D lines, sampled into points along each
# track with sample_3d(), and with a radius buffered with buffer_3d()
# (build); each sample's value is read with extract(), or inside its buffer
# with crop(), as_voxel() (rast() for a variable without depth) and mask()
# (intersect); summarise_3d() gives a buffer's mean, and the samples' values
# are summarised per track, weighted by the distance or time each sample
# stands for (summarise). The wrapper calls only exported functions of those
# steps.

#' Summarise environmental values along tracks
#'
#' Samples points along each track, such as a dive or a vessel route, reads
#' the value at each sample, and summarises the samples of each track. It is
#' a convenience wrapper over the three steps of an analysis, and the
#' [workflow article](../articles/workflow-values-in-shapes.html) shows the
#' same steps written out:
#'
#' 1. **Build:** the vertices of each track become a line with
#'    [as_vector3d()], and [sample_3d()] places points along it, with depth
#'    and time interpolated between vertices. With `horizontal_radius_km`,
#'    [buffer_3d()] turns each sample into a buffer.
#' 2. **Intersect:** without a radius, [extract()][extract-SpatProxy] reads
#'    the value at each sample, in the cell that holds it and at the depth
#'    level `method` picks. With a radius, the data for the sample's time
#'    step are read with [crop()] and [as_voxel()] (or [rast()][rast-SpatProxy]
#'    for a variable without depth), and [mask()] keeps the values inside its
#'    buffer: the cells the circle touches, and the levels whose slab of
#'    water overlaps the depth window. [summarise_3d()] gives their mean,
#'    weighted by cell volume.
#' 3. **Summarise:** the statistics of each track's sample values, each
#'    sample weighted by `weighting`.
#'
#' Each sample is matched to the time step nearest its time, or equal to it
#' with `time_match = "exact"`, as [extract()][extract-SpatProxy] matches
#' points. A sample gets `NA` when it is outside the grid, on a cell or level
#' without a value, or has no time step that matches. One warning counts
#' these samples, by reason. The track summaries leave them out, and `n` is
#' the number of samples with a value.
#'
#' With `weighting = "distance"` each sample stands for half the distance to
#' the samples either side (`support_km`), so a stretch sampled more densely
#' doesn't count for more; with `"duration"` for half the time
#' (`support_seconds`). The weighted standard deviation is the one
#' [summarise_3d()] uses. A track whose samples all have zero weight, such
#' as a track that doesn't move with distance weights, gets equal weights.
#' Neighbouring samples are usually correlated, so the standard deviation
#' describes their spread, not the uncertainty of the mean.
#'
#' @param data Vertices of the tracks, in the order they were travelled: a
#'   data frame, matrix, named list or sf POINT object. sf points are
#'   transformed to longitude and latitude. To use sf lines, convert them
#'   first, e.g. with
#'   `as.data.frame(sf::st_coordinates(sf::st_transform(lines, 4326)))` and
#'   `lon_col = "X", lat_col = "Y", track_col = "L1"`.
#' @param nc The data: a netCDF file path, several paths for a series of
#'   files on one grid, an open [ncdf4::nc_open()] connection, or a
#'   [SpatProxy-class] from [as_proxy()].
#' @param var Name of the variable, for a netCDF source. `NULL` detects the
#'   only data variable of the file.
#' @param track_col Column with track IDs. `NULL` (default) reads all of
#'   `data` as one track, with ID `"1"`. Each track needs at least two
#'   vertices.
#' @param sampling,spacing_km,depth_step_m,time_step,time_unit,max_samples
#'   Where the samples go along each track, as in [sample_3d()]: only the
#'   vertices (default without a step), or with a spacing in kilometres or
#'   a depth or time step, points in between as well.
#' @param horizontal_radius_km Radius of a circle around each sample, in
#'   kilometres, below 10000 as in [buffer_3d()]. `0` (default) reads the
#'   value at the sample itself.
#' @param depth_radius_m Distance in metres above and below each sample's
#'   depth, for the buffers. Needed with a radius for a variable with depth,
#'   and used only with a radius.
#' @param method For a variable with depth and no radius, the depth level
#'   each sample takes, as in [extract()][extract-SpatProxy]: `"nearest"`
#'   (default; the level nearest the sample's depth that has a value),
#'   `"nearest_level"` (the nearest level, even if it has no value),
#'   `"bottom"` (the deepest level with a value) or `"surface"`. Leave
#'   `NULL` with a radius; for a variable without depth, leave `NULL` or use
#'   `"2d"`.
#' @param statistic One or more of `"mean"` (default), `"sd"`, `"min"` and
#'   `"max"`, of the samples of each track.
#' @param weighting How each sample is weighted in the track summaries:
#'   `"distance"` (default) by `support_km`, `"duration"` by
#'   `support_seconds` (needs dates), or `"equal"`.
#' @param lon_col,lat_col,depth_col,date_col Vertex columns. The depth is
#'   needed for a variable with depth, unless `method` is `"bottom"` or
#'   `"surface"`, and for `depth_step_m`. The date is needed for data with
#'   dated time steps, for `time_step` and for duration weights. Dates and
#'   date-times are read as in [as_vector3d()].
#' @param lon_dim,lat_dim,depth_dim,time_dim Dimension names for a netCDF
#'   source, as in [as_proxy()].
#' @param time_match,max_time_diff How a sample's time is matched to a time
#'   step, as in [as_voxel()]: the nearest step (default) or an exact one,
#'   and the largest gap in days for a nearest match.
#' @param output_col Name of the value column in `samples`. `NULL` uses the
#'   variable's name.
#'
#' @returns A list of two data frames:
#'   * `samples`: one row per sample, in order along each track, with
#'     `track_id`, `sample_id`, the position (in `lon_col` and `lat_col`),
#'     the depth and date when used (in `depth_col` and `date_col`),
#'     `distance_km` from the track's start, `support_km`,
#'     `support_seconds`, `is_interpolated` (`FALSE` for the track's own
#'     vertices), and the value in `output_col`. With a radius, also
#'     `<output>_n`, the number of values in the sample's buffer.
#'   * `summary`: one row per track, with `track_id`, a column per
#'     statistic, and `n`, the number of samples with a value.
#'
#' @seealso [sample_3d()], [extract_to_point()], [extract_to_buffer()], and
#'   the [workflow article](../articles/workflow-values-in-shapes.html).
#' @examples
#' # Daily temperature on a 10 x 10 grid of 1 degree cells at four depths.
#' path <- tempfile(fileext = ".nc")
#' lon <- ncdf4::ncdim_def("lon", "degrees_east", seq(0.5, 9.5))
#' lat <- ncdf4::ncdim_def("lat", "degrees_north", seq(0.5, 9.5))
#' depth <- ncdf4::ncdim_def("depth", "m", c(5, 15, 30, 60))
#' time <- ncdf4::ncdim_def("time", "days since 2020-01-01", 0:2)
#' temp <- ncdf4::ncvar_def("temp", "degree_Celsius",
#'                          list(lon, lat, depth, time), -9999)
#' nc <- ncdf4::nc_create(path, temp)
#' ncdf4::ncvar_put(nc, temp, seq_len(10 * 10 * 4 * 3) / 10)
#' ncdf4::nc_close(nc)
#'
#' # A dive from the surface to 50 m over two days.
#' dive <- data.frame(lon = c(2.5, 6.5), lat = c(3.5, 5.5), depth = c(0, 50),
#'                    date = as.Date(c("2020-01-01", "2020-01-03")))
#'
#' # The value at a sample at least every 50 km, and the mean along the dive.
#' track <- extract_to_track(dive, path, spacing_km = 50)
#' track$summary
#'
#' # The mean within 30 km and 10 m above and below each sample instead.
#' extract_to_track(dive, path, spacing_km = 50, horizontal_radius_km = 30,
#'                  depth_radius_m = 10)$summary
#' @export
extract_to_track <- function(
  data, nc, var = NULL, track_col = NULL,
  sampling = NULL, spacing_km = NULL, depth_step_m = NULL, time_step = NULL,
  time_unit = c("hours", "minutes", "seconds", "days", "months", "years"),
  max_samples = 100000L,
  horizontal_radius_km = 0, depth_radius_m = 0, method = NULL,
  statistic = "mean", weighting = c("distance", "duration", "equal"),
  lon_col = "lon", lat_col = "lat", depth_col = "depth", date_col = "date",
  lon_dim = NULL, lat_dim = NULL, depth_dim = NULL, time_dim = NULL,
  time_match = c("nearest", "exact"), max_time_diff = NULL,
  output_col = NULL
) {
  time_unit <- match.arg(time_unit)
  statistic <- match.arg(statistic, c("mean", "sd", "min", "max"),
                         several.ok = TRUE)
  weighting <- match.arg(weighting)
  time_match <- match.arg(time_match)
  stop_track <- function(message) {
    rlang::abort(message, class = "ocean3d_error_track_input", call = NULL)
  }
  for (amount in list(horizontal_radius_km, depth_radius_m)) {
    if (!is.numeric(amount) || length(amount) != 1L || !is.finite(amount) ||
          amount < 0) {
      stop_track(paste(
        "`horizontal_radius_km` and `depth_radius_m` must each be one",
        "number, 0 or more."
      ))
    }
  }
  buffered <- horizontal_radius_km > 0
  if (!buffered && depth_radius_m > 0) {
    stop_track("`depth_radius_m` needs a `horizontal_radius_km` above 0.")
  }

  proxy <- if (methods::is(nc, "SpatProxy")) {
    nc
  } else {
    as_proxy(nc, var = var, lon_dim = lon_dim, lat_dim = lat_dim,
             depth_dim = depth_dim, time_dim = time_dim)
  }
  var_name <- proxy@variable$name
  has_depth <- !is.null(proxy@axes$depth)
  if (!has_depth) {
    if (!is.null(method) && !identical(method, "2d")) {
      stop_track(glue::glue(
        "Variable `{var_name}` has no depth, so `method` doesn't apply; ",
        "leave it NULL."
      ))
    }
    method <- NULL
  } else if (buffered) {
    if (!is.null(method)) stop_track("Leave `method` NULL with a radius.")
    if (depth_radius_m == 0) {
      stop_track(glue::glue(
        "Variable `{var_name}` has depth, so buffers need a ",
        "`depth_radius_m` above 0."
      ))
    }
  } else {
    if (is.null(method)) method <- "nearest"
    rules <- c("nearest", "nearest_level", "bottom", "surface")
    if (!rlang::is_string(method) || !method %in% rules) {
      stop_track(glue::glue("`method` must be one of {toString(rules)}."))
    }
  }
  # Dated time steps, even one, are matched to each sample's time; a
  # climatology gives its steps to every sample.
  time_axis <- proxy@axes$time
  timed <- !is.null(time_axis) && !isTRUE(time_axis$climatology)
  needs_depth <- (has_depth && (buffered ||
                                  method %in% c("nearest", "nearest_level"))) ||
    !is.null(depth_step_m)
  needs_time <- timed || weighting == "duration" || !is.null(time_step)

  # One plain table of the vertices, whatever form `data` takes; sf points
  # give their coordinates in WGS84.
  table <- if (inherits(data, "sf")) {
    xy <- sf::st_coordinates(sf::st_transform(sf::st_geometry(data), 4326))
    out_table <- sf::st_drop_geometry(data)
    out_table[[lon_col]] <- xy[, "X"]
    out_table[[lat_col]] <- xy[, "Y"]
    out_table
  } else {
    as.data.frame(data, stringsAsFactors = FALSE)
  }
  wanted <- c(lon_col, lat_col, track_col, if (needs_depth) depth_col,
              if (needs_time) date_col)
  missing_cols <- setdiff(wanted, names(table))
  if (length(missing_cols)) {
    stop_track(glue::glue(
      "`data` needs the column(s) {toString(missing_cols)}."
    ))
  }
  output <- if (is.null(output_col)) var_name else output_col
  if (!rlang::is_string(output) || !nzchar(output)) {
    stop_track("`output_col` must be one non-empty name.")
  }

  # Build: each track becomes a line, sampled into points.
  work <- data.frame(
    track = if (is.null(track_col)) "1" else as.character(table[[track_col]]),
    lon = table[[lon_col]], lat = table[[lat_col]],
    depth = if (needs_depth) table[[depth_col]] else 0,
    stringsAsFactors = FALSE
  )
  if (needs_time) work$time <- table[[date_col]]
  lines <- as_vector3d(work, type = "lines", id = "track",
                       time = if (needs_time) "time")
  samples <- sample_3d(lines, sampling = sampling, spacing_km = spacing_km,
                       depth_step_m = depth_step_m, time_step = time_step,
                       time_unit = time_unit, max_samples = max_samples)
  out <- as.data.frame(samples)
  if (weighting == "duration" && anyNA(out$support_seconds)) {
    stop_track(paste(
      "`weighting = \"duration\"` needs a date at every vertex. Give the",
      "dates, or use \"distance\" or \"equal\"."
    ))
  }

  # Intersect: one value per sample. The time step each sample takes comes
  # from extract(), which matches all the samples at once.
  rule <- if (!is.null(method)) list(method = method)
  found <- do.call(terra::extract, c(list(proxy, samples), rule, list(
    time_match = time_match, max_time_diff = max_time_diff
  )))
  no_step <- if (timed) is.na(found$time_step) else rep(FALSE, nrow(out))
  value <- found$value
  if (buffered) {
    # For a variable without depth the depth window doesn't matter, but
    # mask() keeps no cell of a plain raster for a buffer with no thickness,
    # so it gets one metre.
    circles <- buffer_3d(samples, radius_km = horizontal_radius_km,
                         depth_window = if (has_depth) depth_radius_m else 1)
    value <- rep(NA_real_, nrow(out))
    n <- rep(0L, nrow(out))
    for (i in which(!no_step)) {
      circle <- circles[i]
      # Cropped first: an error raised while S4 evaluates an argument for
      # dispatch loses its class. A buffer off the grid gets NA.
      source <- tryCatch(crop(proxy, circle),
                         ocean3d_error_proxy_crop = function(e) NULL)
      if (is.null(source)) next
      when <- if (timed) found$time_step[i]
      layers <- if (has_depth) {
        as_voxel(source, time = when, time_match = "exact")
      } else {
        terra::rast(source, time = when, time_match = "exact")
      }
      # The buffers keep the sample's time, which the voxel already matches,
      # so mask()'s message that it drops polygon times is muted.
      kept <- suppressMessages(mask(layers, circle))
      summary <- summarise_3d(kept, statistic = "mean")
      value[i] <- summary$mean
      n[i] <- summary$n
    }
  }

  # One warning for every sample that got NA.
  missing_value <- is.na(value) & !no_step
  reasons <- c(
    if (any(no_step)) glue::glue(
      "{sum(no_step)} with no time step that matches"
    ),
    if (any(missing_value)) glue::glue(
      "{sum(missing_value)} outside the grid or without a value"
    )
  )
  if (length(reasons)) {
    rlang::warn(glue::glue(
      "{sum(is.na(value))} of {nrow(out)} samples get NA: ",
      "{glue::glue_collapse(reasons, sep = ' and ')}."
    ), class = "ocean3d_warning_track_samples")
  }

  # Summarise: the statistics of each track's samples, weighted by the
  # distance or time each stands for, as summarise_3d() weights values.
  weights <- switch(weighting, distance = out$support_km,
                    duration = out$support_seconds, equal = rep(1, nrow(out)))
  ids <- unique(out$track_id)
  summary <- do.call(rbind, lapply(ids, function(id) {
    keep <- out$track_id == id & !is.na(value)
    v <- value[keep]
    w <- weights[keep]
    if (sum(w) <= 0) w <- rep(1, length(v))
    w <- w / sum(w)
    average <- sum(w * v)
    correction <- 1 - sum(w^2)
    stats <- list(
      mean = if (length(v)) average else NA_real_,
      sd = if (length(v) > 1L && correction > 0) {
        sqrt(sum(w * (v - average)^2) / correction)
      } else {
        NA_real_
      },
      min = if (length(v)) min(v) else NA_real_,
      max = if (length(v)) max(v) else NA_real_
    )
    data.frame(track_id = id, stats[statistic], n = length(v),
               stringsAsFactors = FALSE)
  }))
  rownames(summary) <- NULL

  # The samples, track and sample number first, with the caller's column
  # names for position, depth and date, and the columns not used dropped.
  first <- c("track_id", "sample_id")
  out <- out[c(first, setdiff(names(out), c(first, "feature_id", "vertex",
                                            if (!needs_depth) "depth",
                                            if (!needs_time) "time")))]
  renamed <- c(lon = lon_col, lat = lat_col, depth = depth_col,
               time = date_col)
  at <- match(names(renamed), names(out))
  names(out)[at[!is.na(at)]] <- renamed[!is.na(at)]
  added <- c(output, if (buffered) paste0(output, "_n"))
  if (any(added %in% names(out))) {
    stop_track(glue::glue(
      "The samples already have a column named ",
      "{toString(intersect(added, names(out)))}. Choose another `output_col`."
    ))
  }
  out[[output]] <- value
  if (buffered) out[[paste0(output, "_n")]] <- n
  list(samples = out, summary = summary)
}

# Open every netCDF source, select the requested (or detected) variables and
# describe each source/variable pair with `describe(nc, var, dims,
# source_index)`. Checks that each variable uses one method, is temporal in
# all of its sources or none, and has identical value units. Returns the
# records, the sources and the open handles; the caller closes the handles
# (they are closed here if a check fails). Used by the track functions.
.describe_netcdf_sources <- function(nc, var, file_col, dimension_args,
                                     describe) {
  sources <- .as_netcdf_sources(nc, file_col)
  if (!length(sources)) {
    stop("At least one netCDF source is required.", call. = FALSE)
  }
  opened <- list()
  complete <- FALSE
  on.exit(if (!complete) lapply(opened, .close_source), add = TRUE)
  metadata <- list()
  checks <- list()
  for (i in seq_along(sources)) {
    opened[[i]] <- .open_source(sources[[i]])
    handle <- opened[[i]]$nc
    variables <- if (is.null(var)) {
      .detect_netcdf_var(handle)
    } else {
      intersect(var, names(handle$var))
    }
    if (!length(variables)) {
      file <- basename(handle$filename)
      stop(glue::glue("No requested variable is present in {file}."),
           call. = FALSE)
    }
    for (variable in variables) {
      dims <- do.call(.detect_dims, c(list(nc = handle, var = variable),
                                      dimension_args))
      record <- describe(handle, variable, dims, i)
      metadata[[length(metadata) + 1L]] <- record
      checks[[length(checks) + 1L]] <- list(
        variable = variable, method = record$method,
        temporal = !is.null(dims$time),
        units = .read_attribute(handle, variable, "units")
      )
    }
  }
  found <- unique(vapply(checks, function(x) x$variable, character(1)))
  if (!is.null(var) && !setequal(var, found)) {
    missing <- toString(setdiff(var, found))
    stop(glue::glue("Requested variables were not found: {missing}."),
         call. = FALSE)
  }
  for (variable in found) {
    group <- Filter(function(x) x$variable == variable, checks)
    methods <- unique(vapply(group, function(x) x$method, character(1)))
    temporal <- vapply(group, function(x) x$temporal, logical(1))
    if (length(methods) != 1L || (any(temporal) && !all(temporal))) {
      stop(glue::glue(
        "Sources for `{variable}` have incompatible static/depth dimensions."
      ), call. = FALSE)
    }
    units <- lapply(group, function(x) x$units)
    if (!all(vapply(units, identical, logical(1), units[[1L]]))) {
      stop(glue::glue(
        "Sources for `{variable}` use different value units. ",
        "Convert them first."
      ), call. = FALSE)
    }
  }
  complete <- TRUE
  list(metadata = metadata, sources = sources, opened = opened)
}




# One summary of `samples` (a value and its area, depth and time fractions
# per row), with the share of the requested support that has values.
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
  status <- if (!any(selected)) {
    "empty"
  } else if (!any(valid)) {
    "all_missing"
  } else if (!meets_threshold) {
    "insufficient_coverage"
  } else {
    if (fraction < 1) "partial" else "complete"
  }
  value <- NA_real_
  effective_n <- NA_real_
  if (any(valid)) {
    weights <- if (weighting == "weighted") support[valid] else rep(1,
                                                                    sum(valid))
    effective_n <- 1 / sum((weights / sum(weights))^2)
  }
  if (any(valid) && meets_threshold) {
    # The weighted statistics are shared with summarise_3d().
    value <- .weighted_statistic(samples$value[valid], weights, statistic)
    if (is.na(value) && statistic %in% c("sd", "se")) {
      status <- "insufficient_variation_samples"
    }
  }
  list(value = value, coverage_fraction = coverage, valid_fraction = fraction,
       n_intersecting = sum(selected), n_valid = sum(valid), status = status,
       contributing_rows = which(valid), effective_n = effective_n)
}


# The output column for each variable, checked against the input's names.
.buffer_output_names <- function(variables, output_col, input_names) {
  outputs <- if (is.null(output_col)) {
    variables
  } else if (length(variables) == 1L && length(output_col) == 1L) {
    unname(output_col)
  } else {
    if (is.null(names(output_col)) || anyDuplicated(names(output_col)) ||
          !setequal(names(output_col), variables)) {
      stop(paste0("For multiple variables, `output_col` must name each ",
                  "variable once."),
           call. = FALSE)
    }
    unname(output_col[variables])
  }
  if (!is.character(outputs) || anyNA(outputs) || any(!nzchar(outputs)) ||
        anyDuplicated(outputs)) {
    stop("Output column names must be distinct, non-empty strings.",
         call. = FALSE)
  }
  if (any(outputs %in% input_names)) {
    stop("Output names would overwrite existing columns. ",
         "Choose different `output_col` names.", call. = FALSE)
  }
  stats::setNames(outputs, variables)
}


# Shift times by `amount` units: elapsed seconds, or UTC calendar months and
# years, clamped to the end of a shorter month. Used by the track functions.
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
