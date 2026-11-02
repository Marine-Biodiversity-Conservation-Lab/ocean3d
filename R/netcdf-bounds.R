# CF cell bounds: the bounds variable named by a coordinate's `bounds`
# attribute is read and checked against the coordinate (shape, finite values,
# matching units). Used for longitude and latitude edges in point extraction
# and for the cell edges buffer extraction weights by.

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

# Reads the CF bounds variable of a coordinate, for point and buffer
# extraction. Returns list(name, values), where values is a 2 x n matrix of
# cell edges, or NULL when the coordinate declares no bounds.
#
# Unusable bounds (missing variable, wrong shape, non-finite values) are an
# error when `strict = TRUE`. Buffer extraction uses that, because its weights
# depend on the declared edges. Otherwise they give a warning and NULL, and
# point extraction infers cell edges from the centres. Bounds in different
# units from their coordinate are always an error.
.read_bounds <- function(nc, coordinate, n, strict = FALSE) {
  name <- .read_time_attribute(nc, coordinate, "bounds")
  if (is.null(name)) return(NULL)
  unusable <- function(problem) {
    if (strict) {
      stop("Bounds for `", coordinate, "` ", problem, ".", call. = FALSE)
    }
    warning(
      "Bounds for `", coordinate, "` ", problem,
      "; inferring cell edges from the coordinate centres instead.",
      call. = FALSE
    )
    NULL
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
