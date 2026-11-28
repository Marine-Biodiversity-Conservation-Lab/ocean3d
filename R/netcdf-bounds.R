# CF cell bounds: the bounds variable named by a coordinate's `bounds`
# attribute is read and checked against the coordinate (shape, finite values,
# matching units). Without bounds, edges are inferred halfway between
# neighbouring centres. .read_interval_axis() builds and checks the cell
# intervals of every axis (longitude, latitude, depth and time) for the core
# object (as_proxy()).

# Cell edges halfway between neighbouring centres, with the outer edges half
# a step beyond the end centres. `centres` (at least two) may be in any
# order; `lower` and `upper` keep that order. Used by .axis_intervals().
.midpoint_edges <- function(centres) {
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
  list(lower = lower, upper = upper)
}

# Units can be spelled differently and still mean the same thing: angles
# (degrees_east, degree_east, degrees_E, ...) and lengths (m, metres, ...).
.same_units <- function(a, b) {
  a <- tolower(trimws(a))
  b <- tolower(trimws(b))
  if (identical(a, b)) return(TRUE)
  angle <- "^degrees?(_?(east|north|e|n))?$"
  if (grepl(angle, a) && grepl(angle, b)) return(TRUE)
  isTRUE(.length_unit_factor(a) == .length_unit_factor(b))
}

# Unusable or inconsistent cell bounds (#131): the only producer of
# `ocean3d_error_bounds`, for .read_bounds() and .axis_intervals(). The
# message names the file, so a bad file can be found in a multi-file call.
.abort_bounds <- function(coordinate, file, problem) {
  where <- if (is.null(file)) "" else glue::glue(" in {file}")
  rlang::abort(glue::glue("Bounds for `{coordinate}`{where} {problem}."),
               class = "ocean3d_error_bounds", call = NULL)
}

# Reads the CF bounds variable of a coordinate. Returns list(name, values),
# where values is a 2 x n matrix of cell edges, or NULL when the coordinate
# declares no bounds. Unusable bounds (missing variable, wrong shape,
# non-finite values, different units from the coordinate) are an error:
# declared edges are used as given, never replaced by inferred ones.
.read_bounds <- function(nc, coordinate, n) {
  name <- .read_attribute(nc, coordinate, "bounds")
  if (is.null(name)) return(NULL)
  unusable <- function(problem) {
    .abort_bounds(coordinate, basename(nc$filename), problem)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !nzchar(name)) {
    return(unusable("have an invalid bounds reference"))
  }
  # In a netCDF-4 group, CF names are relative to the coordinate's group.
  if (!name %in% names(nc$var) && grepl("/", coordinate, fixed = TRUE)) {
    grouped <- file.path(dirname(coordinate), name)
    if (grouped %in% names(nc$var)) name <- grouped
  }
  if (!name %in% names(nc$var)) {
    return(unusable(paste0(
      "were not found (missing bounds variable `", name, "`)"
    )))
  }

  # A coordinate variable that isn't itself a dimension must be 1D; its
  # bounds then follow that dimension.
  dimension <- if (coordinate %in% names(nc$dim)) {
    coordinate
  } else {
    dims <- nc$var[[coordinate]]$dim
    if (length(dims) != 1L) {
      stop("Coordinate `", coordinate, "` must be one-dimensional.",
           call. = FALSE)
    }
    dims[[1L]]$name
  }
  axes <- vapply(nc$var[[name]]$dim, function(x) x$name, character(1))
  coordinate_axis <- which(axes == dimension)
  raw <- ncdf4::ncvar_get(nc, name, collapse_degen = FALSE)
  if (length(dim(raw)) != 2L || length(coordinate_axis) != 1L) {
    return(unusable("have unsupported dimensions"))
  }
  if (coordinate_axis == 1L) raw <- t(raw)
  if (!identical(dim(raw), c(2L, as.integer(n))) || !is.numeric(raw) ||
        any(!is.finite(raw))) {
    return(unusable("must contain two finite endpoints per coordinate"))
  }

  coordinate_units <- .read_attribute(nc, coordinate, "units")
  bounds_units <- .read_attribute(nc, name, "units")
  if (!is.null(bounds_units) && !is.null(coordinate_units) &&
        !.same_units(bounds_units, coordinate_units)) {
    unusable(glue::glue(
      "and the coordinate use different units (`{bounds_units}` and ",
      "`{coordinate_units}`). Convert bounds and coordinate to the same ",
      "units first"
    ))
  }
  list(name = name, values = raw)
}

# Stops with `ocean3d_error_axis_order` unless `values` run strictly one way,
# so that each value is one cell. Both callers reject non-finite values
# first (.axis_intervals() and .prepare_depth_axis()); an NA would still count
# as out of order rather than fail inside `if`. The only
# producer of that class, so every stored axis follows the same rule. Called
# by .axis_intervals() (every axis read by .read_interval_axis()) and by
# .read_depth_axis().
.check_axis_order <- function(values, coordinate, file = NULL) {
  steps <- diff(values)
  if (length(steps) &&
        !(isTRUE(all(steps > 0)) || isTRUE(all(steps < 0)))) {
    where <- if (is.null(file)) "" else glue::glue(" in {file}")
    rlang::abort(glue::glue(
      "Coordinate `{coordinate}`{where} must be strictly increasing or ",
      "decreasing, without duplicates. Sort the file along this dimension ",
      "first."
    ), class = "ocean3d_error_axis_order", call = NULL)
  }
  invisible(values)
}

# Shared edges and outer edges can differ by rounding, up to single precision
# for values stored as float32.
.interval_tolerance <- function(lower, upper) {
  max(1e-8, max(abs(c(lower, upper)), na.rm = TRUE) * 2^-23)
}

# Why cell intervals `lower` to `upper` don't fit their `centres`, as the end
# of a sentence starting "Bounds for ...", or NULL when they do: each needs a
# positive width (or zero, when `zero_width`), must contain its centre and
# mustn't overlap another, allowing for rounding. Returns the problem rather
# than stopping, so the SpatVoxel validity method can report it too; readers
# stop through .axis_intervals().
.interval_problem <- function(centres, lower, upper, zero_width = FALSE) {
  if (any(!is.finite(lower)) || any(!is.finite(upper)) ||
        any(!is.finite(upper - lower))) {
    return("need positive widths and must contain their centres")
  }
  tol <- .interval_tolerance(lower, upper)
  if ((!zero_width && any(upper <= lower)) ||
        any(centres < lower - tol | centres > upper + tol)) {
    return("need positive widths and must contain their centres")
  }
  index <- order(lower)
  if (length(index) > 1L &&
        any(tail(lower[index], -1L) < head(cummax(upper[index]), -1L) - tol)) {
    return("overlap. Resolve overlapping cells first")
  }
  NULL
}

# Cell intervals of one axis: `centres` with `bounds` (a 2 x n matrix of
# edges, or NULL to infer them halfway between centres, cut at the surface for
# depth). Checks that the centres run one way without duplicates and that the
# intervals have positive widths, contain their centres and don't overlap
# (also modulo 360 when `periodic`, for longitude cells a horizontal buffer
# could count twice), allowing for float32 rounding. Bad bounds stop with
# `ocean3d_error_bounds`. A single centre without bounds is an error, unless
# `allow_singleton`, when it becomes a zero-width cell. Values keep their
# storage order. Called by .read_interval_axis() for every axis, and by the
# core object's time axis (as_proxy()).
.axis_intervals <- function(centres, bounds, axis, coordinate,
                            allow_singleton = FALSE, periodic = FALSE,
                            file = NULL) {
  if (!length(centres) || any(!is.finite(centres))) {
    stop("Coordinates must contain finite values.", call. = FALSE)
  }
  .check_axis_order(centres, coordinate, file)

  clipped_at_surface <- rep(FALSE, length(centres))
  if (is.null(bounds) && length(centres) < 2L) {
    if (!allow_singleton) {
      stop("Cannot infer bounds from a single coordinate in `", coordinate,
           "`. Supply explicit bounds; one centre has no known extent.",
           call. = FALSE)
    }
    # One cell without bounds has no known extent: a zero-width cell.
    lower <- upper <- centres
    basis <- "singleton_coordinate"
  } else if (is.null(bounds)) {
    edges <- .midpoint_edges(centres)
    lower <- edges$lower
    upper <- edges$upper
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

  problem <- .interval_problem(centres, lower, upper,
                               zero_width = basis == "singleton_coordinate")
  if (!is.null(problem)) .abort_bounds(coordinate, file, problem)
  if (periodic && length(centres) > 1L) {
    # Longitude cells must not overlap modulo 360 either (e.g. 350 to 370
    # and 0 to 20), or a cell would be counted twice.
    tol <- .interval_tolerance(lower, upper)
    start <- lower %% 360
    end <- start + upper - lower
    sorted <- order(start)
    if (any(start[sorted][-1L] < head(cummax(end[sorted]), -1L) - tol) ||
          max(end) - 360 > min(start) + tol) {
      .abort_bounds(coordinate, file, paste(
        "overlap around the globe (longitude cells modulo 360).",
        "Resolve overlapping cells first"
      ))
    }
  }

  # Keep original record order so values and intervals retain their alignment.
  list(
    centres = centres, lower = lower, upper = upper,
    source_index = seq_along(centres), basis = basis,
    clipped_at_surface = clipped_at_surface
  )
}

# Reads one coordinate with its cell intervals, through .axis_intervals():
# longitude and latitude in degrees, depth in metres positive down (with CF
# bounds converted the same way), time in UTC seconds. Longitudes are
# unwrapped so neighbouring centres differ by less than 180 degrees, and each
# cell's bounds are wrapped around its centre (a cell 360 degrees or wider
# becomes centre +/- 180); latitude edges are clamped to +/-90. `periodic`
# also rejects longitude cells that overlap modulo 360. Without CF bounds,
# edges are inferred only when `infer_bounds`, and a single coordinate is a
# zero-width cell only when `allow_singleton`; both default to TRUE for
# longitude and latitude (the defaults are evaluated after match.arg()).
# Depth order is checked by .read_depth_axis() before .axis_intervals()
# checks it again. Used by the core object (as_proxy()).
.read_interval_axis <- function(
  nc,
  coordinate,
  axis = c("depth", "time", "lon", "lat"),
  infer_bounds = axis %in% c("lon", "lat"),
  allow_singleton = axis %in% c("lon", "lat"),
  periodic = FALSE
) {
  axis <- match.arg(axis)
  if (!is.logical(infer_bounds) || length(infer_bounds) != 1L ||
        is.na(infer_bounds)) {
    stop("`infer_bounds` must be TRUE or FALSE.", call. = FALSE)
  }
  raw <- as.vector(.get_dim_values(nc, coordinate))
  if (axis %in% c("lon", "lat")) {
    problem <- if (!is.numeric(raw) || !length(raw) || any(!is.finite(raw))) {
      "must be finite numbers"
    } else if (axis == "lat" && any(raw < -90 | raw > 90)) {
      "must be within [-90, 90]"
    }
    if (!is.null(problem)) {
      rlang::abort(glue::glue(
        "Coordinate `{coordinate}` in {basename(nc$filename)} {problem}."
      ), class = "ocean3d_error_spatial_coordinates", call = NULL)
    }
  }
  bounds <- .read_bounds(nc, coordinate, length(raw))
  units <- .read_attribute(nc, coordinate, "units")
  calendar <- NULL
  direction_basis <- NULL

  if (axis == "lon") {
    centres <- .unwrap_longitude(raw)
    if (!is.null(bounds)) {
      # Each column is one cell: its edges as offsets from the stored centre,
      # wrapped to [-180, 180) and moved to the unwrapped centre.
      edges <- bounds$values
      whole <- abs(edges[1L, ] - edges[2L, ]) >= 360
      edges <- rep(centres, each = 2L) +
        .normalise_longitude(edges - rep(raw, each = 2L))
      edges[, whole] <- rbind(centres[whole] - 180, centres[whole] + 180)
      bounds$values <- edges
    }
  } else if (axis == "lat") {
    centres <- raw
  } else if (axis == "depth") {
    depth <- .read_depth_axis(nc, coordinate, raw)
    centres <- depth$depth_m
    direction_basis <- depth$direction_basis

    if (!is.null(bounds)) {
      bounds$values <- .depth_bounds_m(nc, coordinate, raw, bounds)
    }
  } else {
    centres <- as.numeric(.decode_netcdf_time(nc, coordinate))
    calendar <- .read_attribute(nc, coordinate, "calendar")
    if (is.null(calendar)) calendar <- "standard"
    if (!is.null(bounds)) {
      .check_time_bounds_calendar(nc, coordinate, bounds)
      bounds$values <- .decode_netcdf_time(
        nc, coordinate, bounds = bounds$values
      )
    }
  }

  if (is.null(bounds) && !infer_bounds) {
    stop("Coordinate `", coordinate, "` has no explicit bounds. ",
         paste0("Use a dataset with bounds or explicitly set ",
                "`infer_bounds = TRUE` "),
         "to estimate midpoint intervals.", call. = FALSE)
  }

  result <- .axis_intervals(
    centres,
    if (is.null(bounds)) NULL else bounds$values,
    axis,
    coordinate,
    allow_singleton,
    periodic = axis == "lon" && periodic,
    file = basename(nc$filename)
  )
  if (axis == "lat") {
    result$lower <- pmax(result$lower, -90)
    result$upper <- pmin(result$upper, 90)
  }
  c(result, list(
    coordinate = coordinate, axis = axis, longitude = axis == "lon",
    units = switch(axis, depth = "m", time = "UTC seconds since 1970-01-01",
                   lon = "degrees_east", lat = "degrees_north"),
    source_units = units,
    bounds_name = if (is.null(bounds)) NULL else bounds$name,
    calendar = calendar, direction_basis = direction_basis
  ))
}
