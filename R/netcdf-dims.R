# Which dimension of a netCDF variable is longitude, latitude, depth and time.
# Dimensions are matched from CF metadata (`axis`, `standard_name`, `units`),
# then by name. Coordinate values and attributes are read for both coordinate
# variables and plain dimensions, and two-dimensional (curvilinear) coordinates
# are rejected. Shared by point, buffer and track extraction.

.get_dim_values <- function(nc, dim_name) {
  if (is.null(dim_name)) {
    return(NULL)
  }

  if (dim_name %in% names(nc$dim)) {
    return(nc$dim[[dim_name]]$vals)
  }

  # Some files store coordinate values as variables rather than only as
  # dimension values. This fallback covers those files.
  if (dim_name %in% names(nc$var)) {
    return(ncdf4::ncvar_get(nc, dim_name))
  }

  stop(
    "Could not read dimension or coordinate variable: ", dim_name,
    call. = FALSE
  )
}

.has_coordinate_variable <- function(nc, coordinate) {
  if (coordinate %in% names(nc$dim)) {
    return(isTRUE(nc$dim[[coordinate]]$create_dimvar))
  }
  coordinate %in% names(nc$var)
}

.read_time_attribute <- function(nc, coordinate, attribute) {
  # ncdf4 prints warnings to stdout for dimensions without a coordinate
  # variable, so return early instead of asking it for their attributes.
  if (!.has_coordinate_variable(nc, coordinate)) return(NULL)
  att <- ncdf4::ncatt_get(nc, coordinate, attribute)

  if (!isTRUE(att$hasatt)) {
    return(NULL)
  }

  att$value
}

# Read a coordinate attribute, including dimension metadata.
.axis_attribute <- function(nc, dim_name, attribute) {
  if (dim_name %in% names(nc$var)) {
    att <- ncdf4::ncatt_get(nc, dim_name, attribute)
    if (isTRUE(att$hasatt)) return(att$value)
  }
  if (dim_name %in% names(nc$dim) && attribute == "units") {
    return(nc$dim[[dim_name]]$units)
  }
  NULL
}

.normalise_axis_attribute <- function(x) {
  if (!is.character(x) || length(x) != 1L ||
        is.na(x)) {
    return("")
  }
  tolower(trimws(x))
}

# Detect CF axis metadata associated with a dimension.
.axis_from_metadata <- function(nc, dim_name) {
  axis <- toupper(.normalise_axis_attribute(
    .axis_attribute(nc, dim_name, "axis")
  ))
  standard_name <- .normalise_axis_attribute(
    .axis_attribute(nc, dim_name, "standard_name")
  )
  units <- .normalise_axis_attribute(
    .axis_attribute(nc, dim_name, "units")
  )
  lon_units <- c(
                 "degrees_east", "degree_east",
                 "degrees_e", "degree_e")

  lat_units <- c(
    "degrees_north", "degree_north",
    "degrees_n", "degree_n"
  )
  if (standard_name == "longitude" || units %in% lon_units) {
    return("lon")
  }
  if (standard_name == "latitude" || units %in% lat_units) {
    return("lat")
  }
  if (standard_name == "time" ||
    grepl(
      "^(milliseconds?|seconds?|minutes?|hours?|days?) since ",
      units
    )) {
    return("time")
  }
  if (standard_name %in% c("depth", "height", "altitude")) {
    return("depth")
  }

  if (axis == "X") return("lon")
  if (axis == "Y") return("lat")
  if (axis == "Z") return("depth")
  if (axis == "T") return("time")

  NULL
}

# Find one dimension belonging to the selected variable.
.resolve_variable_axis <- function(
  nc, var_dims, explicit, candidates, axis_name,
  required = FALSE) {

  if (!is.null(explicit)) {
    if (!is.character(explicit) ||
          length(explicit) != 1L ||
          is.na(explicit) ||
          !nzchar(trimws(explicit))) {
      stop(
        "`", axis_name,
        "_dim` must be a single non-empty dimension name.",
        call. = FALSE
      )
    }

    if (!explicit %in% var_dims) {
      stop(
        "Explicit `", axis_name, "_dim = '", explicit,
        "'` is not a dimension of the selected variable. ",
        "Available dimensions: ", toString(var_dims),
        call. = FALSE
      )
    }

    return(explicit)
  }

  metadata_matches <- var_dims[
    vapply(
      var_dims,
      function(d) {
        identical(.axis_from_metadata(nc, d), axis_name)
      },
      logical(1)
    )
  ]

  if (length(metadata_matches) > 1L) {
    stop(
      "Ambiguous ", axis_name, " axis: ",
      toString(metadata_matches),
      ". Specify `", axis_name, "_dim` explicitly.",
      call. = FALSE
    )
  }

  if (length(metadata_matches) == 1L) {
    return(metadata_matches)
  }

  name_matches <- candidates[candidates %in% var_dims]

  if (length(name_matches) > 1L) {
    stop(
      "Ambiguous ", axis_name, " dimension names: ",
      toString(name_matches),
      ". Specify `", axis_name, "_dim` explicitly.",
      call. = FALSE
    )
  }

  if (length(name_matches) == 1L) {
    return(name_matches)
  }

  if (isTRUE(required)) {
    stop(
      "Could not identify the ", axis_name,
      " dimension of the selected variable. ",
      "Available dimensions: ", toString(var_dims),
      ". Specify `", axis_name, "_dim` if necessary.",
      call. = FALSE
    )
  }

  NULL
}

# Reject grids that need unsupported 2D coordinate lookup.
.check_auxiliary_coordinates <- function(nc, var, var_dims) {
  coordinates <- .axis_attribute(nc, var, "coordinates")

  if (is.null(coordinates)) return(invisible(TRUE))

  coordinates <- strsplit(
    trimws(as.character(coordinates)[1]),
    "\\s+"
  )[[1]]

  for (coordinate in coordinates) {
    if (!coordinate %in% names(nc$var)) next

    standard_name <- .normalise_axis_attribute(
      .axis_attribute(nc, coordinate, "standard_name")
    )

    if (!standard_name %in% c("longitude", "latitude")) next

    coord_dims <- vapply(
      nc$var[[coordinate]]$dim,
      function(d) d$name,
      character(1)
    )

    if (length(coord_dims) > 1L &&
          all(coord_dims %in% var_dims)) {
      stop(
        "Variable `", var, "` uses a two-dimensional ",
        "geographic coordinate `", coordinate, "`. ",
        "Curvilinear longitude/latitude grids are not ",
        "yet supported by extract_to_point().",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.detect_dims <- function(nc,
                         var,
                         lon_dim = NULL,
                         lat_dim = NULL,
                         depth_dim = NULL,
                         time_dim = NULL) {

  if (!var %in% names(nc$var)) {
    stop(
      "Variable `", var, "` not found in netCDF file. ",
      "Available variables are: ",
      toString(names(nc$var)),
      call. = FALSE
    )
  }

  variable_dims <- nc$var[[var]]$dim

  var_dims <- vapply(
    variable_dims,
    function(d) d$name,
    character(1)
  )

  dim_lengths <- vapply(
    variable_dims,
    function(d) as.integer(d$len),
    integer(1)
  )

  .check_auxiliary_coordinates(nc, var, var_dims)

  lon_name <- .resolve_variable_axis(
    nc, var_dims, lon_dim,
    c("lon", "longitude", "x", "nav_lon"),
    "lon", required = TRUE
  )

  lat_name <- .resolve_variable_axis(
    nc, var_dims, lat_dim,
    c("lat", "latitude", "y", "nav_lat"),
    "lat", required = TRUE
  )

  depth_name <- .resolve_variable_axis(
    nc, var_dims, depth_dim,
    c("depth", "deptht", "lev", "level", "z"),
    "depth"
  )

  time_name <- .resolve_variable_axis(
    nc, var_dims, time_dim,
    c("time", "time_counter", "valid_time", "t"),
    "time"
  )

  assigned <- unlist(
    list(lon_name, lat_name, depth_name, time_name),
    use.names = FALSE
  )

  if (anyDuplicated(assigned)) {
    stop(
      "A netCDF dimension was assigned to multiple axes. ",
      "Check the dimension names and CF metadata.",
      call. = FALSE
    )
  }

  unknown <- setdiff(var_dims, assigned)

  unknown_lengths <- dim_lengths[
    match(unknown, var_dims)
  ]

  unsupported <- unknown[unknown_lengths > 1L]

  if (length(unsupported) > 0L) {
    stop(
      "Unrecognised non-singleton dimension(s) ",
      "in variable `", var, "`: ",
      toString(unsupported), ". ",
      "Specify the appropriate *_dim argument or ",
      "use a dataset with supported dimensions.",
      call. = FALSE
    )
  }

  singleton_dims <- unknown[unknown_lengths == 1L]

  projected_units <- c(
    "m", "meter", "meters", "metre", "metres",
    "km", "kilometer", "kilometers",
    "kilometre", "kilometres"
  )

  for (axis_info in list(
    list(name = lon_name, axis = "longitude"),
    list(name = lat_name, axis = "latitude")
  )) {
    units <- .normalise_axis_attribute(
      .axis_attribute(nc, axis_info$name, "units")
    )

    if (units %in% projected_units) {
      stop(
        "Dimension `", axis_info$name,
        "` uses projected coordinates (", units, "). ",
        "Geographic longitude/latitude extraction ",
        "requires coordinates in degrees. ",
        "Projected grids are not yet supported.",
        call. = FALSE
      )
    }
  }

  list(
    lon = lon_name,
    lat = lat_name,
    depth = depth_name,
    time = time_name,
    var_dims = var_dims,
    singleton_dims = singleton_dims
  )
}
