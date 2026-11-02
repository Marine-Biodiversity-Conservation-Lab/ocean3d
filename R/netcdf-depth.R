# The depth axis of a netCDF variable, as depths in metres increasing
# downward: units are converted to metres, the `positive` direction is read
# (or inferred for heights and altitudes), and formula-based or pressure
# coordinates are rejected. Shared by point, buffer and track extraction.

# Metres per unit for the supported length units, or NA.
.length_unit_factor <- function(units) {
  switch(
    tolower(trimws(units)),
    "m" = 1,
    "meter" = 1,
    "meters" = 1,
    "metre" = 1,
    "metres" = 1,
    "km" = 1000,
    "kilometer" = 1000,
    "kilometers" = 1000,
    "kilometre" = 1000,
    "kilometres" = 1000,
    "cm" = 0.01,
    "centimeter" = 0.01,
    "centimeters" = 0.01,
    "centimetre" = 0.01,
    "centimetres" = 0.01,
    NA_real_
  )
}

# Interpret physical depths without changing their positions in the source.
.prepare_depth_axis <- function(values, units, positive = NULL) {
  if (!is.numeric(values) || !length(values) ||
        any(!is.finite(values))) {
    stop("Depth coordinates must contain finite numeric values.", call. = FALSE)
  }
  if (!is.character(units) || length(units) != 1L || is.na(units)) {
    stop(
      "Depth units cannot be converted to metres. ",
      "Supply a coordinate with explicit supported length units. ",
      "For a file with a single depth level, extract2d() reads it ",
      "without depth metadata.",
      call. = FALSE
    )
  }
  factor <- .length_unit_factor(units)
  if (!is.finite(factor)) {
    stop(
      "Depth units cannot be converted to metres. ",
      "Check the vertical coordinate metadata; supported units are ",
      "m, meter(s), metre(s), km, kilometer(s), kilometre(s), cm, ",
      "centimeter(s) and centimetre(s). For a file with a single depth ",
      "level, extract2d() reads it without depth metadata.",
      call. = FALSE
    )
  }
  if (is.null(positive)) {
    positive <- if (all(values >= 0)) {
      "down"
    } else if (all(values <= 0)) {
      "up"
    } else {
      stop(
        "Depth direction is ambiguous. ",
        "The coordinate contains mixed signs and has no ",
        "usable positive attribute.",
        call. = FALSE
      )
    }
    direction_basis <- "inferred"
  } else {
    if (!is.character(positive) || length(positive) != 1L || is.na(positive)) {
      stop(
        "The depth positive attribute must be 'up' or 'down'.",
        call. = FALSE
      )
    }
    positive <- tolower(trimws(positive))
    direction_basis <- "metadata"
  }
  if (!positive %in% c("up", "down")) {
    stop("The depth positive attribute must be 'up' or 'down'.", call. = FALSE)
  }
  depth_m <- as.numeric(values) * factor
  if (positive == "up") depth_m <- -depth_m
  if (any(depth_m < 0)) {
    stop(
      "The vertical coordinates do not describe supported ",
      "non-negative depths below the reference surface. ",
      "Check the sign and vertical reference.",
      call. = FALSE
    )
  }
  if (anyDuplicated(depth_m)) {
    stop(
      "The depth coordinate contains duplicate physical depths.",
      call. = FALSE
    )
  }
  list(
    depth_m = depth_m,
    surface_index = which.min(depth_m),
    direction_basis = direction_basis
  )
}

.read_depth_axis <- function(nc, coordinate, values) {
  if (!.has_coordinate_variable(nc, coordinate)) {
    stop(
      "The depth dimension `", coordinate, "` has no coordinate variable, ",
      "so its physical depths and units are unknown. Add a coordinate ",
      "variable with depth values and units to the file before extraction. ",
      "For a file with a single depth level, extract2d() reads it without ",
      "depth metadata.",
      call. = FALSE
    )
  }
  standard_name <- .normalise_axis_attribute(
    .read_time_attribute(nc, coordinate, "standard_name")
  )
  formula_terms <- .read_time_attribute(nc, coordinate, "formula_terms")
  if (!is.null(formula_terms) ||
        grepl("^(ocean|atmosphere)_.*coordinate$", standard_name) ||
        grepl("pressure", standard_name, fixed = TRUE)) {
    stop(
      "This vertical coordinate requires a physical transformation ",
      "before it can be interpreted as depth in metres. ",
      "Convert sigma, formula-based or pressure coordinates explicitly ",
      "before point extraction.",
      call. = FALSE
    )
  }
  positive <- .read_time_attribute(nc, coordinate, "positive")
  # Heights and altitudes increase upward; CF requires `positive` on them,
  # but do not read a non-compliant height axis as depths below the surface.
  if (is.null(positive) && standard_name %in% c("height", "altitude")) {
    positive <- "up"
  }
  .prepare_depth_axis(
    values,
    units = .read_time_attribute(nc, coordinate, "units"),
    positive = positive
  )
}
