# CF cell bounds: the bounds variable named by a coordinate's `bounds`
# attribute is read and checked against the coordinate (shape, finite values,
# matching units). Used for longitude and latitude edges in point extraction
# and for the cell edges buffer extraction weights by. Without bounds, edges
# are inferred halfway between neighbouring centres. Also builds and checks
# the cell intervals of depth and time axes, for buffer extraction and the
# core object.

# Cell edges halfway between neighbouring centres, with the outer edges half
# a step beyond the end centres. `centres` (at least two) may be in any
# order; `lower` and `upper` keep that order. `regular` is TRUE when the
# steps are equal. Shared by the longitude and latitude axes, buffer
# extraction's depth and time intervals, and the core object's depth axis.
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
  regular <- all(abs(gaps - gaps[1L]) <= max(1e-8, gaps[1L] * 1e-8))
  list(lower = lower, upper = upper, regular = regular)
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

# Reads the CF bounds variable of a coordinate. Returns list(name, values),
# where values is a 2 x n matrix of cell edges, or NULL when the coordinate
# declares no bounds. Unusable bounds (missing variable, wrong shape,
# non-finite values, different units from the coordinate) are an error:
# declared edges are used as given, never replaced by inferred ones.
.read_bounds <- function(nc, coordinate, n) {
  name <- .read_time_attribute(nc, coordinate, "bounds")
  if (is.null(name)) return(NULL)
  unusable <- function(problem) {
    stop("Bounds for `", coordinate, "` ", problem, ".", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !nzchar(name)) {
    return(unusable("have an invalid bounds reference"))
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

  coordinate_units <- .read_time_attribute(nc, coordinate, "units")
  bounds_units <- .read_time_attribute(nc, name, "units")
  if (!is.null(bounds_units) && !is.null(coordinate_units) &&
        !.same_units(bounds_units, coordinate_units)) {
    stop(
      "Bounds for `", coordinate, "` and the coordinate use different ",
      "units (`", bounds_units, "` and `", coordinate_units, "`). ",
      "Convert bounds and coordinate to the same units first.",
      call. = FALSE
    )
  }
  list(name = name, values = raw)
}

# Cell intervals of one axis: `centres` with `bounds` (a 2 x n matrix of
# edges, or NULL to infer them halfway between centres, cut at the surface for
# depth). Checks that the centres run one way without duplicates and that the
# intervals have positive widths, contain their centres and don't overlap
# (also modulo 360 when `periodic`, for longitude). Bad bounds stop with
# `ocean3d_error_bounds`. A single centre without bounds is an error, unless
# `allow_singleton`, when it becomes a zero-width cell. Values keep their
# storage order. Shared by every axis reader: longitude and latitude
# (.prepare_spatial_axis()), buffer extraction's depth and time, and the core
# object's depth and time
# and the core object's depth axis (as_proxy()).
.axis_intervals <- function(centres, bounds, axis, coordinate,
                            allow_singleton = FALSE, periodic = FALSE) {
  if (!length(centres) || any(!is.finite(centres))) {
    stop("Coordinates must contain finite values.", call. = FALSE)
  }
  steps <- diff(centres)
  if (length(steps) && !(all(steps > 0) || all(steps < 0))) {
    rlang::abort(glue::glue(
      "Coordinate `{coordinate}` must be strictly increasing or ",
      "decreasing, without duplicates."
    ), class = "ocean3d_error_axis_order", call = NULL)
  }

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

  if (any(!is.finite(lower)) || any(!is.finite(upper)) ||
        (basis != "singleton_coordinate" && any(upper <= lower)) ||
        any(!is.finite(upper - lower)) ||
        any(centres < lower - 1e-8 | centres > upper + 1e-8)) {
    rlang::abort(glue::glue(
      "Invalid bounds for `{coordinate}`: intervals need positive widths ",
      "and must contain their centres."
    ), class = "ocean3d_error_bounds", call = NULL)
  }
  index <- order(lower)
  if (length(index) > 1L &&
        any(tail(lower[index], -1L) < head(cummax(upper[index]), -1L))) {
    rlang::abort(glue::glue(
      "Bounds overlap in `{coordinate}`. Resolve overlapping cells first."
    ), class = "ocean3d_error_bounds", call = NULL)
  }
  if (periodic && length(index) > 1L) {
    # Longitude cells must not overlap modulo 360 either (e.g. 350 to 370
    # and 0 to 20), or a cell would be counted twice.
    start <- lower %% 360
    end <- start + upper - lower
    sorted <- order(start)
    if (any(start[sorted][-1L] < head(cummax(end[sorted]), -1L) - 1e-8) ||
          max(end) - 360 > min(start) + 1e-8) {
      rlang::abort(glue::glue(
        "Longitude cells in `{coordinate}` overlap around the globe. ",
        "Resolve overlapping cells first."
      ), class = "ocean3d_error_bounds", call = NULL)
    }
  }

  # Keep original record order so values and intervals retain their alignment.
  list(
    centres = centres, lower = lower, upper = upper,
    source_index = seq_along(centres), basis = basis,
    clipped_at_surface = clipped_at_surface
  )
}

# Reads a depth or time coordinate with its cell intervals, through
# .axis_intervals(): depth in metres positive down (with CF bounds converted
# the same way), time in UTC seconds. Without CF bounds, edges are inferred
# only when `infer_bounds = TRUE`. Used by buffer extraction and the core
# object's depth axis.
.read_interval_axis <- function(
  nc,
  coordinate,
  axis = c("depth", "time"),
  infer_bounds = FALSE,
  allow_singleton = FALSE
) {
  axis <- match.arg(axis)
  if (!is.logical(infer_bounds) || length(infer_bounds) != 1L ||
        is.na(infer_bounds)) {
    stop("`infer_bounds` must be TRUE or FALSE.", call. = FALSE)
  }
  raw <- as.vector(.get_dim_values(nc, coordinate))
  bounds <- .read_bounds(nc, coordinate, length(raw))
  units <- .read_time_attribute(nc, coordinate, "units")
  calendar <- NULL
  direction_basis <- NULL

  if (axis == "depth") {
    depth <- .read_depth_axis(nc, coordinate, raw)
    centres <- depth$depth_m
    direction_basis <- depth$direction_basis

    if (!is.null(bounds)) {
      bounds$values <- .depth_bounds_m(nc, coordinate, raw, bounds)
    }
  } else {
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
    allow_singleton
  )
  c(result, list(
    coordinate = coordinate, axis = axis,
    units = if (axis == "depth") "m" else "UTC seconds since 1970-01-01",
    source_units = units,
    bounds_name = if (is.null(bounds)) NULL else bounds$name,
    calendar = calendar, direction_basis = direction_basis
  ))
}
