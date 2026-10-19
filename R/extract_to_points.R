#-------------------------------------------------------------------------------
# Point-based extraction from netCDF files
#-------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------
# These functions are intentionally small and focused. This makes the code
# easier to test and easier to maintain inside a package.
# -----------------------------------------------------------------------------

.check_required_cols <- function(x, cols, x_name = "data") {
  cols <- stats::na.omit(cols)
  missing_cols <- setdiff(cols, names(x))

  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required column(s) in `%s`: %s",
        x_name,
        toString(missing_cols)
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.is_ncdf4 <- function(x) {
  inherits(x, "ncdf4")
}


.normalize_points_input <- function(x, lon_col = "lon", lat_col = "lat") {
  # data.frame subclasses such as tibbles already use the internal format.
  if (is.data.frame(x) && !inherits(x, "sf")) {
    return(x)
  }

  # sf objects can provide longitude and latitude directly from POINT geometry.
  if (inherits(x, "sf")) {
    geom_type <- unique(as.character(sf::st_geometry_type(x)))

    if (length(geom_type) != 1L || geom_type != "POINT") {
      stop("`sf` input must contain POINT geometries.", call. = FALSE)
    }

    coords_source <- x

    if (!isTRUE(sf::st_is_longlat(coords_source))) {
      if (is.na(sf::st_crs(coords_source))) {
        stop(
          "`sf` input must have a geographic CRS or a CRS that can be ",
          "transformed to longitude/latitude.",
          call. = FALSE
        )
      }

      coords_source <- sf::st_transform(coords_source, 4326)
    }

    coords <- sf::st_coordinates(coords_source)
    x[[lon_col]] <- coords[, 1]
    x[[lat_col]] <- coords[, 2]

    return(x)
  }

  # Matrices are converted to a data frame while preserving column names.
  if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop("Matrix input must have column names.", call. = FALSE)
    }

    return(as.data.frame(x))
  }

  # Named lists are converted to a data frame.
  if (is.list(x)) {
    if (is.null(names(x)) || any(names(x) == "")) {
      stop("List input must be named.", call. = FALSE)
    }

    return(as.data.frame(x))
  }

  stop(
    "`data` must be a data frame, tibble, sf object, matrix, or named list.",
    call. = FALSE
  )
}


.recycle_point_inputs <- function(lon, lat, depth = NULL, date = NULL) {
  inputs <- list(lon = lon, lat = lat, depth = depth, date = date)
  lengths <- vapply(
    inputs[!vapply(inputs, is.null, logical(1))],
    length,
    integer(1)
  )
  n <- max(lengths)

  incompatible <- lengths != 1L & lengths != n
  if (any(incompatible)) {
    stop(
      "Direct coordinate inputs must have compatible lengths: each must ",
      "have length 1 or the same maximum length.",
      call. = FALSE
    )
  }

  recycle <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 1L && n > 1L) rep(x, n) else x
  }

  lapply(inputs, recycle)
}


.as_netcdf_sources <- function(nc, file_col = "file") {
  # This helper standardises all accepted input formats into the same
  # internal representation. Each source keeps the original object plus a type
  # flag telling the extraction code whether the file needs to be opened and
  # closed.

  if (is.character(nc)) {
    sources <- lapply(nc, function(x) list(source = x, type = "path"))
    return(sources)
  }

  # ncdf4 objects and data frames are lists too, so both are matched before
  # the list branch below.
  if (.is_ncdf4(nc)) {
    return(list(list(source = nc, type = "connection")))
  }

  if (is.data.frame(nc)) {
    .check_required_cols(nc, file_col, "nc")
    files <- as.character(nc[[file_col]])
    sources <- lapply(files, function(x) list(source = x, type = "path"))
    return(sources)
  }

  # Users often pass nested lists when they group files by variable or by year.
  # Each element is resolved on its own, recursing into nested lists, and the
  # results are joined in order.
  if (is.list(nc)) {
    sources <- lapply(nc, .as_netcdf_sources, file_col = file_col)
    return(Reduce(c, sources, list()))
  }

  stop(
    "Unsupported `nc` input. Use a file path, vector of file paths, ",
    "opened ncdf4 object, list, or data frame with a file column.",
    call. = FALSE
  )
}

.open_source <- function(source) {
  if (identical(source$type, "path")) {
    if (!file.exists(source$source)) {
      stop("netCDF file does not exist: ", source$source, call. = FALSE)
    }

    return(list(nc = ncdf4::nc_open(source$source), close = TRUE))
  }

  if (identical(source$type, "connection")) {
    return(list(nc = source$source, close = FALSE))
  }

  stop("Unknown internal netCDF source type.", call. = FALSE)
}

.close_source <- function(opened) {
  if (isTRUE(opened$close)) {
    try(ncdf4::nc_close(opened$nc), silent = TRUE)
  }

  invisible(TRUE)
}


.detect_netcdf_var <- function(nc) {
  # When `var = NULL`, the package tries to detect the variable automatically.
  # This is useful for workflows where each netCDF file contains one main
  # environmental variable, for example one file for NPPV and another file for
  # salinity.
  #
  # The detection is deliberately conservative. If more than one variable is
  # present, the function stops and asks the user to provide `var` manually.
  # This avoids extracting the wrong variable by accident.

  vars <- names(nc$var)

  if (length(vars) == 1L) {
    return(vars[1])
  }

  stop(
    "Automatic variable detection requires exactly one variable in each ",
    "netCDF file. ",
    "Variables found: ",
    toString(vars),
    ". Please specify `var` manually.",
    call. = FALSE
  )
}

.resolve_source_var <- function(source, var) {
  # If the user provides `var`, the same variable name is used for all sources.
  # If `var = NULL`, each source is opened briefly and its variable is detected.

  if (!is.null(var)) {
    return(var)
  }

  opened <- .open_source(source)
  on.exit(.close_source(opened), add = TRUE)

  .detect_netcdf_var(opened$nc)
}

.resolve_output_prefix <- function(method, var, output_prefix, auto_var) {
  # If the user provides `output_prefix` or `output_col`, that name is
  # respected. If `var = NULL`, the output name is generated from the detected
  # variable and the extraction method, so multi-variable calls create clear
  # column names:
  #   bottom  -> seabottom_nppv, seabottom_so
  #   surface -> surface_nppv, surface_so
  #   nearest -> nearest_nppv, nearest_so
  #   all     -> surface_*, nearest_* and seabottom_*
  #
  # If the user provided `var`, the previous behaviour is kept.

  if (!is.null(output_prefix)) {
    return(output_prefix)
  }

  if (!isTRUE(auto_var)) {
    return(NULL)
  }

  if (identical(method, "surface")) {
    return(paste0("surface_", var))
  }

  if (identical(method, "bottom")) {
    return(paste0("seabottom_", var))
  }

  if (identical(method, "nearest")) {
    return(paste0("nearest_", var))
  }

  var
}

.match_name <- function(available,
                        candidates,
                        required = TRUE,
                        what = "dimension") {
  hit <- candidates[candidates %in% available]

  if (length(hit) > 0) {
    return(hit[1])
  }

  if (isTRUE(required)) {
    stop(
      sprintf(
        "Could not find %s. Tried: %s. Available names are: %s",
        what,
        toString(candidates),
        toString(available)
      ),
      call. = FALSE
    )
  }

  NULL
}

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

.read_time_attribute <- function(nc, coordinate, attribute) {
  att <- ncdf4::ncatt_get(nc, coordinate, attribute)

  if (!isTRUE(att$hasatt)) {
    return(NULL)
  }

  att$value
}


.decode_netcdf_time <- function(nc, time_dim) {
  if (is.null(time_dim)) {
    return(NULL)
  }

  units <- .read_time_attribute(nc, time_dim, "units")

  if (!is.character(units) || length(units) != 1L ||
        is.na(units) || !nzchar(trimws(units))) {
    stop(
      "Time coordinate `", time_dim,
      "` has no usable `units` attribute. ",
      "Expected units such as 'hours since 2020-01-01 00:00:00'. ",
      "Check the file metadata or use a corrected dataset.",
      call. = FALSE
    )
  }

  calendar <- .read_time_attribute(nc, time_dim, "calendar")

  if (is.null(calendar)) {
    calendar <- "standard"
  }

  supported <- c("standard", "gregorian", "proleptic_gregorian")

  if (!is.character(calendar) || length(calendar) != 1L ||
        is.na(calendar) ||
        !tolower(calendar) %in% supported) {
    stop(
      "Time coordinate `", time_dim,
      "` uses an unsupported calendar: ", toString(calendar), ". ",
      "This reader currently supports standard, gregorian and ",
      "proleptic_gregorian calendars.",
      call. = FALSE
    )
  }

  pattern <- paste0(
    "^\\s*",
    "(milliseconds?|seconds?|minutes?|hours?|days?)",
    "\\s+since\\s+(.+?)\\s*$"
  )

  parts <- regmatches(
    units,
    regexec(pattern, units, ignore.case = TRUE, perl = TRUE)
  )[[1]]

  if (length(parts) != 3L) {
    stop(
      "Unsupported time units for `", time_dim, "`: '", units, "'. ",
      "Expected milliseconds, seconds, minutes, hours or days ",
      "followed by 'since' and a reference date-time. ",
      "Months, years and climatologies require separate handling.",
      call. = FALSE
    )
  }

  unit <- sub("s$", "", tolower(parts[2]))

  seconds_per_unit <- switch(
    unit,
    millisecond = 0.001,
    second = 1,
    minute = 60,
    hour = 3600,
    day = 86400
  )

  origin_text <- trimws(parts[3])
  origin_text <- sub("T", " ", origin_text, fixed = TRUE)
  origin_text <- sub(
    "\\s*(UTC|Z)$", "",
    origin_text,
    ignore.case = TRUE
  )
  origin_text <- trimws(origin_text)

  # This first reader supports origins written in UTC without
  # numeric timezone offsets. Reject other forms explicitly.
  date_pattern <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}$"
  datetime_pattern <- paste0(
    "^[0-9]{4}-[0-9]{2}-[0-9]{2} ",
    "[0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]+)?$"
  )

  if (grepl(date_pattern, origin_text)) {
    origin_text <- paste(origin_text, "00:00:00")
  } else if (!grepl(datetime_pattern, origin_text)) {
    stop(
      "Unsupported reference date-time in units '", units, "'. ",
      "Expected YYYY-MM-DD or YYYY-MM-DD HH:MM:SS in UTC, ",
      "optionally followed by Z or UTC. ",
      "Numeric timezone offsets are not supported by this reader yet.",
      call. = FALSE
    )
  }

  origin <- as.POSIXct(
    strptime(
      origin_text,
      format = "%Y-%m-%d %H:%M:%OS",
      tz = "UTC"
    )
  )

  # Round-trip validation catches invalid dates and times that
  # a platform might otherwise normalise.
  if (is.na(origin) ||
        format(origin, "%Y-%m-%d %H:%M:%S", tz = "UTC") !=
          substr(origin_text, 1L, 19L)) {
    stop(
      "Invalid reference date-time in units '", units, "'.",
      call. = FALSE
    )
  }

  values <- as.vector(.get_dim_values(nc, time_dim))

  if (!is.numeric(values) || any(!is.finite(values))) {
    stop(
      "Time coordinate `", time_dim,
      "` must contain finite numeric offsets.",
      call. = FALSE
    )
  }

  timestamps <- origin + values * seconds_per_unit

  # R's ordinary date-time arithmetic does not implement the
  # historical mixed Julian/Gregorian calendar.
  if (tolower(calendar) %in% c("standard", "gregorian")) {
    transition <- as.POSIXct(
      "1582-10-15 00:00:00",
      tz = "UTC"
    )

    if (origin < transition || any(timestamps < transition)) {
      stop(
        "Historical standard/gregorian dates before 1582-10-15 ",
        "require a calendar-aware decoder.",
        call. = FALSE
      )
    }
  }

  timestamps
}


.nearest_index <- function(values, target) {
  if (length(values) == 0 || all(is.na(values)) || is.na(target)) {
    return(NA_integer_)
  }

  which.min(abs(values - target))
}

.prepare_time_input <- function(
  x,
  time_resolution = c(
    "auto", "datetime", "daily", "monthly", "yearly"
  )) {
  time_resolution <- match.arg(time_resolution)

  if (inherits(x, "POSIXt")) {
    values <- as.POSIXct(x, tz = "UTC")
    has_time <- TRUE
  } else if (inherits(x, "Date")) {
    values <- as.POSIXct(x, tz = "UTC")
    has_time <- FALSE
  } else if (is.character(x)) {
    date_pattern <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}$"
    datetime_pattern <- paste0(
      "^[0-9]{4}-[0-9]{2}-[0-9]{2}[ T]",
      "[0-9]{2}:[0-9]{2}:[0-9]{2}Z?$"
    )
    present <- !is.na(x)
    date_only <- present & grepl(date_pattern, x)
    datetime <- present & grepl(datetime_pattern, x)
    bad_format <- which(present & !date_only & !datetime)

    if (length(bad_format) > 0L) {
      shown <- head(bad_format, 10L)
      stop(
        "Time values could not be converted. Invalid format at ",
        "row/position(s): ", toString(shown), ". Values: ",
        toString(shQuote(x[shown])), ". ",
        "Use YYYY-MM-DD or YYYY-MM-DD HH:MM:SS in UTC. ",
        "Date and POSIXct objects are also accepted.",
        call. = FALSE
      )
    }

    if (time_resolution == "auto" &&
          any(date_only) && any(datetime)) {
      stop(
        "Time input mixes dates with and without hours. ",
        "Set time_resolution explicitly or use a consistent format.",
        call. = FALSE
      )
    }

    text <- sub("T", " ", x, fixed = TRUE)
    text <- sub("Z$", "", text)
    text[date_only] <- paste(text[date_only], "00:00:00")
    values <- as.POSIXct(
      strptime(text, "%Y-%m-%d %H:%M:%S", tz = "UTC")
    )
    round_trip <- format(
      values, "%Y-%m-%d %H:%M:%S", tz = "UTC"
    )
    invalid <- which(
      present & (is.na(values) | round_trip != text)
    )

    if (length(invalid) > 0L) {
      shown <- head(invalid, 10L)
      stop(
        "Time values could not be converted. Invalid date/time at ",
        "row/position(s): ", toString(shown), ". Values: ",
        toString(shQuote(x[shown])), ".",
        call. = FALSE
      )
    }

    has_time <- any(datetime)
    if (time_resolution == "datetime" && any(date_only)) {
      stop(
        "time_resolution = 'datetime' requires hours for every ",
        "non-missing observation. Supply POSIXct values or ",
        "YYYY-MM-DD HH:MM:SS text.",
        call. = FALSE
      )
    }
  } else {
    stop(
      "Time values could not be converted. ",
      "Use Date, POSIXct or ISO-formatted character values.",
      call. = FALSE
    )
  }

  if (time_resolution == "auto") {
    time_resolution <- if (has_time) "datetime" else "daily"
  }

  if (time_resolution == "datetime" &&
        !has_time && any(!is.na(values))) {
    stop(
      "time_resolution = 'datetime' requires observation hours. ",
      "Supply POSIXct values or YYYY-MM-DD HH:MM:SS text.",
      call. = FALSE
    )
  }

  list(values = values, resolution = time_resolution)
}


.time_comparison_values <- function(values, resolution) {
  if (inherits(values, "Date")) {
    values <- as.POSIXct(values, tz = "UTC")
  }
  if (resolution == "datetime") return(values)

  pattern <- switch(
    resolution,
    daily = "%Y-%m-%d",
    monthly = "%Y-%m-01",
    yearly = "%Y-01-01"
  )
  as.POSIXct(
    as.Date(format(values, pattern, tz = "UTC")),
    tz = "UTC"
  )
}


.check_time_resolution <- function(values, resolution) {
  keys <- .time_comparison_values(values, resolution)
  valid <- !is.na(keys)

  if (anyDuplicated(keys[valid])) {
    stop(
      "Multiple dataset timestamps map to the same '",
      resolution, "' comparison value. ",
      "Selecting the first would discard information. ",
      "Use a finer time_resolution, select the intended records, ",
      "or aggregate the data explicitly before extraction.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


.match_time_index <- function(
  values,
  target,
  time_match = c("nearest", "exact"),
  max_time_diff = NULL,
  time_resolution = c(
    "auto", "datetime", "daily", "monthly", "yearly"
  )) {
  time_match <- match.arg(time_match)
  time_resolution <- match.arg(time_resolution)

  if (length(values) == 0L ||
        all(is.na(values)) ||
        is.na(target)) {
    return(NA_integer_)
  }

  prepared <- .prepare_time_input(target, time_resolution)
  resolution <- prepared$resolution
  .check_time_resolution(values, resolution)

  available <- .time_comparison_values(values, resolution)
  requested <- .time_comparison_values(
    prepared$values, resolution
  )

  if (time_match == "exact") {
    hit <- which(available == requested)
    return(
      if (length(hit) > 0L) as.integer(hit[1L]) else NA_integer_
    )
  }

  differences <- abs(as.numeric(
    difftime(available, requested, units = "days")
  ))
  idx <- which.min(differences)

  if (!is.null(max_time_diff) &&
        differences[idx] > max_time_diff) {
    return(NA_integer_)
  }
  as.integer(idx)
}

.prepare_lon <- function(lon_values, lon_target) {
  # Many global netCDF files store longitude from 0 to 360, whereas observation
  # tables often use -180 to 180. This conversion is applied only when the file
  # clearly uses 0 to 360 and the observation longitude is negative.
  if (is.na(lon_target)) {
    return(lon_target)
  }

  if (min(lon_values, na.rm = TRUE) >= 0 &&
        max(lon_values, na.rm = TRUE) > 180 &&
        lon_target < 0) {
    lon_target <- lon_target + 360
  }

  lon_target
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


.build_start_count <- function(var_dims, index_list,
                               singleton_dims = character(0)) {
  # Build start/count in the exact dimension order of the variable.
  # Every dimension must be explicitly indexed, except known singletons.
  
  start <- integer(length(var_dims))
  count <- integer(length(var_dims))
  
  for (i in seq_along(var_dims)) {
    dim_name <- var_dims[i]
    idx <- index_list[[dim_name]]
    
    if (!is.null(idx)) {
      start[i] <- as.integer(idx$start)
      count[i] <- as.integer(idx$count)
    } else if (dim_name %in% singleton_dims) {
      start[i] <- 1L
      count[i] <- 1L
    } else {
      stop(
        "No extraction index was assigned to dimension `",
        dim_name, "`. Refusing to read the full dimension ",
        "because this could return an incorrect value.",
        call. = FALSE
      )
    }
  }
  
  list(start = start, count = count)
}


.extract_array <- function(nc,
                           var,
                           dims,
                           lon_idx,
                           lat_idx,
                           time_idx = NULL,
                           depth_start = NULL,
                           depth_count = NULL) {
  index_list <- list()
  index_list[[dims$lon]] <- list(start = lon_idx, count = 1L)
  index_list[[dims$lat]] <- list(start = lat_idx, count = 1L)

  if (!is.null(dims$time) && !is.null(time_idx)) {
    index_list[[dims$time]] <- list(start = time_idx, count = 1L)
  }

  if (!is.null(dims$depth)) {
    if (!is.null(depth_start)) {
      # Explicit depth selection for 3D extraction methods.
      index_list[[dims$depth]] <- list(
        start = depth_start,
        count = depth_count
      )
    } else if (nc$dim[[dims$depth]]$len == 1L) {
      # A single depth level can be selected without ambiguity.
      index_list[[dims$depth]] <- list(
        start = 1L,
        count = 1L
      )
    } else {
      stop(
        "Variable `", var, "` has multiple depth levels ",
        "but no depth index was provided.",
        call. = FALSE
      )
    }
  }
  
  sc <- .build_start_count(dims$var_dims, index_list, singleton_dims = dims$singleton_dims)
  
  ncdf4::ncvar_get(
    nc = nc,
    varid = var,
    start = sc$start,
    count = sc$count
  )
}

.extract_one_file <- function(data,
                              source,
                              var,
                              method,
                              lon_col,
                              lat_col,
                              date_col,
                              depth_col,
                              id_col,
                              lon_dim,
                              lat_dim,
                              depth_dim,
                              time_dim,
                              time_match,
                              max_time_diff,
                              output_prefix,
                              verbose,
                              time_resolution = "auto") {
  opened <- .open_source(source)
  on.exit(.close_source(opened), add = TRUE)
  nc <- opened$nc

  dims <- .detect_dims(
    nc = nc,
    var = var,
    lon_dim = lon_dim,
    lat_dim = lat_dim,
    depth_dim = depth_dim,
    time_dim = time_dim
  )

  if (method != "2d" && is.null(dims$depth)) {
    stop(
      "Method `", method, "` requires a depth dimension, but none was found.",
      call. = FALSE
    )
  }
  
  # A 2D extraction cannot choose between multiple depth levels.
  if (method == "2d" && !is.null(dims$depth)) {
    if (length(.get_dim_values(nc, dims$depth)) > 1L) {
      stop(
        "Variable `", var, "` has multiple depth levels. ",
        "Use method = 'surface' or another 3D method ",
        "(bottom, nearest, all).",
        call. = FALSE
      )
    }
  }
  
  lon_values <- .get_dim_values(nc, dims$lon)
  lat_values <- .get_dim_values(nc, dims$lat)
  depth_values <- .get_dim_values(nc, dims$depth)
  time_values <- .decode_netcdf_time(nc, dims$time)
  resolved_time_resolution <- time_resolution

  if (!is.null(dims$time) && !is.null(date_col)) {
    prepared <- .prepare_time_input(
      data[[date_col]], time_resolution
    )
    date_values <- prepared$values
    resolved_time_resolution <- prepared$resolution
    .check_time_resolution(
      time_values, resolved_time_resolution
    )
  } else {
    date_values <- rep(
      as.POSIXct(NA_character_, tz = "UTC"),
      nrow(data)
    )
  }

  prefix <- if (is.null(output_prefix)) var else output_prefix

  if (method == "all") {
    # `all` now returns three summary columns, not one column per depth layer.
    out <- data
    # nearest valid depth summary column.
    out[[paste0("nearest_", prefix)]] <- NA_real_
    # surface layer summary column.
    out[[paste0("surface_", prefix)]] <- NA_real_
    # deepest valid layer summary column.
    out[[paste0("seabottom_", prefix)]] <- NA_real_
  } else {
    out <- data
    out[[prefix]] <- NA_real_
  }

  for (i in seq_len(nrow(data))) {
    if (isTRUE(verbose) && i %% 100 == 0) {
      message("Extracting row ", i, " of ", nrow(data))
    }

    lon_target <- .prepare_lon(lon_values, data[[lon_col]][i])
    lat_target <- data[[lat_col]][i]

    lon_idx <- .nearest_index(lon_values, lon_target)
    lat_idx <- .nearest_index(lat_values, lat_target)

    if (is.na(lon_idx) || is.na(lat_idx)) next

    time_idx <- NULL
    if (!is.null(time_values)) {
      time_idx <- .match_time_index(
        values = time_values,
        target = date_values[i],
        time_match = time_match,
        max_time_diff = max_time_diff,
        time_resolution = resolved_time_resolution
      )
      if (is.na(time_idx)) next
    }

    if (method == "2d") {
      val <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx
      )
      out[[prefix]][i] <- as.numeric(val)[1]
    }

    if (method == "surface") {
      val <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = 1L
      )
      out[[prefix]][i] <- as.numeric(val)[1]
    }

    if (method == "nearest") {
      depth_target <- data[[depth_col]][i]

      # read the full profile before selecting the nearest depth.
      # this allows the function to ignore depth layers that are NA,
      # which usually correspond to layers below the seabed or outside
      # the valid water column for that grid cell.
      profile <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = length(depth_values)
      )
      profile <- as.numeric(profile)
      valid <- which(!is.na(profile))

      if (length(valid) > 0 && !is.na(depth_target)) {
        depth_idx <- valid[which.min(abs(depth_values[valid] - depth_target))]
        out[[prefix]][i] <- profile[depth_idx]
      }
    }

    if (method == "bottom") {
      profile <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = length(depth_values)
      )
      profile <- as.numeric(profile)
      valid <- which(!is.na(profile))
      if (length(valid) > 0) {
        out[[prefix]][i] <- profile[utils::tail(valid, 1)]
      }
    }

    if (method == "all") {
      # `all` now extracts one vertical profile and summarises it into
      # nearest valid depth, surface and bottom columns. This mirrors
      # applying extract3d_nearest(), extract3d_surface() and
      # extract3d_bottom() separately, but avoids reading the same profile
      # three times.
      profile <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = length(depth_values)
      )
      profile <- as.numeric(profile)
      valid <- which(!is.na(profile))

      out[[paste0("surface_", prefix)]][i] <- profile[1]

      if (length(valid) > 0) {
        out[[paste0("seabottom_", prefix)]][i] <- profile[utils::tail(valid, 1)]

        depth_target <- data[[depth_col]][i]
        if (!is.na(depth_target)) {
          depth_idx <- valid[which.min(abs(depth_values[valid] - depth_target))]
          out[[paste0("nearest_", prefix)]][i] <- profile[depth_idx]
        }
      }
    }
  }

  out
}

.merge_source_outputs <- function(outputs, original_names) {
  # When several files are provided, the function extracts from each file and
  # then combines the results. For ordinary single-variable extraction, the
  # first non-missing value across files is used. This is useful when files
  # represent different time periods or spatial tiles.
  if (length(outputs) == 1L) {
    return(outputs[[1]])
  }

  out <- outputs[[1]][original_names]
  new_cols <- unique(unlist(lapply(
    outputs,
    function(x) setdiff(names(x), original_names)
  )))

  for (col in new_cols) {
    values <- rep(NA_real_, nrow(out))

    for (x in outputs) {
      if (!col %in% names(x)) next
      replace <- is.na(values) & !is.na(x[[col]])
      values[replace] <- x[[col]][replace]
    }

    out[[col]] <- values
  }

  out
}

# -----------------------------------------------------------------------------
# Main extraction function
# -----------------------------------------------------------------------------

#' Extract values from netCDFs to one or more point observations
#'
#' This is the general point-extraction interface. It can extract environmental
#' values for one or multiple observations. The user chooses the type of
#' extraction with the `method` argument. The more specific functions
#' `extract2d()`, `extract3d_surface()`, `extract3d_bottom()`,
#' `extract3d_nearest()` and `extract3d_all()` are wrappers around this
#' function.
#'
#' @param data Optional object containing observation points. Can be a data
#'   frame, tibble, `sf` object with POINT geometries, matrix with column names,
#'   or named list. If `NULL`, coordinates can be supplied directly through
#'   `lon`, `lat`, and optionally `depth` and `date`. For `sf` objects,
#'   longitude and latitude are obtained from the POINT geometries; projected
#'   geometries are transformed internally to longitude and latitude when a CRS
#'   is available.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of
#'   file paths, opened `ncdf4` object, list of opened objects, or a data frame
#'   with a file column.
#' @param lon Optional numeric longitude value or vector for direct point
#'   extraction when `data = NULL`.
#' @param lat Optional numeric latitude value or vector for direct point
#'   extraction when `data = NULL`.
#' @param depth Optional numeric observation depth value or vector for direct
#'   point extraction. Required when `data = NULL` and `method = "nearest"` or
#'   `method = "all"`.
#' @param date Optional observation date or vector of dates for direct point
#'   extraction. Used when `data = NULL` and the netCDF variable has a time
#'   dimension. Accepts Date and POSIXt objects, or character values in
#'   "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS" format (UTC). Direct input
#'   vectors may have length 1 or a common compatible length.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param method Extraction method. Options are `"2d"`, `"surface"`, `"bottom"`,
#'   `"nearest"` and `"all"`. With `method = "all"`, the function returns
#'   nearest, surface and bottom outputs together.
#' @param lon_col Name of the longitude column in `data`.
#' @param lat_col Name of the latitude column in `data`.
#' @param date_col Name of the date column in `data`. Required only when the
#'   netCDF variable has a time dimension.
#' @param depth_col Name of the observation-depth column in `data`. Required
#'   only when `method = "nearest"` or `method = "all"`.
#' @param id_col Optional observation identifier column. It is not required for
#'   extraction, but it is checked when provided to help users detect mistakes.
#' @param file_col Column containing file paths when `nc` is a data frame.
#' @param lon_dim Optional name of the longitude dimension in the netCDF file.
#' @param lat_dim Optional name of the latitude dimension in the netCDF file.
#' @param depth_dim Optional name of the depth dimension in the netCDF file.
#' @param time_dim Optional name of the time dimension in the netCDF file.
#' @param time_match Character. How observation dates are matched to netCDF time
#'   values. `"nearest"` (default) selects the closest available date; `"exact"`
#'   requires an exact date match and returns `NA` when no match exists.
#' @param max_time_diff Optional non-negative numeric value giving the maximum
#'   allowed difference, in days, between an observation date and the nearest
#'   netCDF date when `time_match = "nearest"`. If `NULL` (default), no maximum
#'   temporal distance is imposed. Ignored when `time_match = "exact"`.
#' @param output_prefix Optional name used for the extracted output column. If
#'   omitted, `var` is used.
#' @param output_col Optional alias for `output_prefix`. Use this when you want
#'   to define the exact name of the output column directly. If both
#'   `output_col` and `output_prefix` are provided, `output_col` is used. With
#'   `method = "all"`, this name is used as the suffix in `nearest_*`,
#'   `surface_*` and `seabottom_*` columns.
#' @param time_resolution Precision used to compare observation times with
#'   dataset timestamps. "auto" uses daily comparison for Date objects
#'   and date-only text, and datetime comparison for POSIXt objects and
#'   date-time text. Explicit alternatives are "datetime", "daily",
#'   "monthly" and "yearly". Calendar comparisons use UTC.
#'   These options compare timestamp labels; they do not infer averaging
#'   intervals or calculate temporal summaries. Multiple dataset records
#'   mapping to the same comparison value cause an error.
#'   With nearest matching, monthly and yearly comparisons use the start
#'   of the labelled month or year. max_time_diff remains in days and
#'   applies to these comparison values. Use exact matching to require
#'   the same labelled day, month or year.
#' @param verbose Logical. If `TRUE`, progress messages are printed every 100
#'   rows.
#'
#' @details
#' Numeric netCDF time coordinates must provide units expressed as
#' milliseconds, seconds, minutes, hours or days since a reference
#' date-time. Singular and plural unit names are accepted.
#'
#' Reference date-times are interpreted in UTC. A trailing "Z" or
#' "UTC" is accepted; numeric timezone offsets such as "+02:00"
#' are not currently supported.
#'
#' Supported calendars are "standard", "gregorian" and
#' "proleptic_gregorian". The first two are supported only for
#' reference dates and timestamps on or after 1582-10-15.
#' Unsupported units, calendars or missing time units produce an error.
#'
#' Temporal resolution controls comparison of timestamp labels.
#' It does not identify averaging intervals, calculate temporal means,
#' or interpret climatological periods.
#'
#' @return If `data` is provided, the input observations with extracted values
#'   appended. Data frames, tibbles and `sf` objects retain their structure;
#'   matrices and named lists are returned as data frames. If `data = NULL`,
#'   direct extraction returns a numeric value for one point and one output, a
#'   numeric vector for multiple points and one output, a named numeric vector
#'   for one point and multiple outputs (for example, `method = "all"`), or a
#'   data frame for multiple points and multiple outputs.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' # Extract the value closest to a single observation
#' extract_to_point(
#'   nc = nc_file,
#'   lon = 0,
#'   lat = 40,
#'   depth = 60,
#'   date = "2020-01-03",
#'   var = "temp",
#'   method = "nearest"
#' )
#'
#' # Extract values for multiple observations
#' extract_to_point(
#'   nc = nc_file,
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   depth = c(60, 100),
#'   date = "2020-01-03",
#'   var = "temp",
#'   method = "nearest"
#' )
#'
#' # Match a complete observation timestamp.
#' # The example file has a reading at midnight on 2020-01-03.
#' extract_to_point(
#'   nc = nc_file,
#'   lon = 0,
#'   lat = 40,
#'   depth = 60,
#'   date = as.POSIXct("2020-01-03 00:15:00", tz = "UTC"),
#'   var = "temp",
#'   method = "nearest",
#'   time_resolution = "datetime",
#'   max_time_diff = 30 / 1440
#' )
#'
#' # Compare monthly labels when the file has one record per month.
#' # Here the timestamp labels are illustrative.
#' monthly_labels <- as.Date(c("2020-01-15", "2020-02-15"))
#' observation_date <- as.Date("2020-02-28")
#' format(monthly_labels, "%Y-%m") ==
#'   format(observation_date, "%Y-%m")
#'
#' # For extraction from such a file, use:
#' # time_resolution = "monthly", time_match = "exact".
#' # Multiple records in the same month cause an error;
#' # they are not averaged automatically.
#'
#' @export

extract_to_point <- function(data = NULL,
                             nc,
                             lon = NULL,
                             lat = NULL,
                             depth = NULL,
                             date = NULL,
                             var = NULL,
                             method = c(
                               "2d", "surface", "bottom", "nearest", "all"
                             ),
                             lon_col = "lon",
                             lat_col = "lat",
                             date_col = "date",
                             depth_col = "depth",
                             id_col = NULL,
                             file_col = "file",
                             lon_dim = NULL,
                             lat_dim = NULL,
                             depth_dim = NULL,
                             time_dim = NULL,
                             time_match = c("nearest", "exact"),
                             max_time_diff = NULL,
                             output_prefix = NULL,
                             output_col = NULL,
                             verbose = FALSE,
                             time_resolution = c(
                               "auto", "datetime", "daily",
                               "monthly", "yearly"
                             )) {
  method <- match.arg(method)
  time_match <- match.arg(time_match)
  time_resolution <- match.arg(time_resolution)

  if (!is.null(max_time_diff) &&
        (!is.numeric(max_time_diff) || length(max_time_diff) != 1L ||
           is.na(max_time_diff) || max_time_diff < 0)) {
    stop("`max_time_diff` must be `NULL` or one non-negative numeric value ",
         "in days.",
         call. = FALSE)
  }

  direct_input <- is.null(data)

  if (direct_input) {
    if (!is.numeric(lon) || length(lon) == 0L || anyNA(lon)) {
      stop("`lon` must contain numeric values.", call. = FALSE)
    }

    if (!is.numeric(lat) || length(lat) == 0L || anyNA(lat)) {
      stop("`lat` must contain numeric values.", call. = FALSE)
    }

    if (method %in% c("nearest", "all") &&
          (!is.numeric(depth) || length(depth) == 0L || anyNA(depth))) {
      stop(
        "`depth` must contain numeric values for `nearest` and `all` ",
        "extraction.",
        call. = FALSE
      )
    }

    direct_values <- .recycle_point_inputs(
      lon = lon,
      lat = lat,
      depth = depth,
      date = date
    )

    data <- data.frame(
      lon = direct_values$lon,
      lat = direct_values$lat
    )

    lon_col <- "lon"
    lat_col <- "lat"

    if (!is.null(direct_values$depth)) {
      data$depth <- direct_values$depth
      depth_col <- "depth"
    }

    if (!is.null(direct_values$date)) {
      # Validate without replacing the original precision.
      .prepare_time_input(direct_values$date, time_resolution)
      data$date <- direct_values$date
      date_col <- "date"
    } else {
      date_col <- NULL
    }
  } else {
    data <- .normalize_points_input(
      data,
      lon_col = lon_col,
      lat_col = lat_col
    )
  }

  if (!is.null(var) &&
        (!is.character(var) || length(var) != 1L || is.na(var))) {
    stop(
      "`var` must be `NULL` or a single character value naming one netCDF ",
      "variable.",
      call. = FALSE
    )
  }

  # allow users to define the output column name with `output_col`.
  # `output_prefix` is kept for backward compatibility with previous code.
  if (!is.null(output_col)) {
    if (!is.character(output_col) ||
          length(output_col) != 1L ||
          is.na(output_col)) {
      stop(
        "`output_col` must be a single character value naming the output ",
        "column.",
        call. = FALSE
      )
    }
    output_prefix <- output_col
  }

  required_cols <- c(lon_col, lat_col)
  if (!is.null(date_col)) required_cols <- c(required_cols, date_col)
  # `all` needs depth to calculate nearest-depth output.
  if (method %in% c("nearest", "all")) {
    required_cols <- c(required_cols, depth_col)
  }
  if (!is.null(id_col)) required_cols <- c(required_cols, id_col)
  .check_required_cols(data, required_cols, "data")

  sources <- .as_netcdf_sources(nc, file_col = file_col)

  auto_var <- is.null(var)

  outputs <- lapply(sources, function(src) {
    source_var <- .resolve_source_var(src, var)

    source_output_prefix <- .resolve_output_prefix(
      method = method,
      var = source_var,
      output_prefix = output_prefix,
      auto_var = auto_var
    )

    .extract_one_file(
      data = data,
      source = src,
      var = source_var,
      method = method,
      lon_col = lon_col,
      lat_col = lat_col,
      date_col = date_col,
      depth_col = depth_col,
      id_col = id_col,
      lon_dim = lon_dim,
      lat_dim = lat_dim,
      depth_dim = depth_dim,
      time_dim = time_dim,
      time_match = time_match,
      max_time_diff = max_time_diff,
      output_prefix = source_output_prefix,
      verbose = verbose,
      time_resolution = time_resolution
    )
  })

  out <- .merge_source_outputs(outputs, original_names = names(data))

  if (direct_input) {
    value_cols <- setdiff(names(out), names(data))

    if (length(value_cols) == 1L) {
      return(out[[value_cols]])
    }

    if (nrow(out) == 1L) {
      return(unlist(out[1, value_cols], use.names = TRUE))
    }

    return(out[value_cols])
  }

  out
}

# -----------------------------------------------------------------------------
# Wrapper functions
# -----------------------------------------------------------------------------
# These functions keep explicit names for the most common extraction modes.
# They are easier to discover and easier to document, while the internal
# workflow stays centralised in `extract_to_point()`.
# -----------------------------------------------------------------------------

#' Extract values from a 2D netCDF variable
#'
#' Use this function when the variable has longitude and latitude dimensions,
#' and optionally a time dimension, but no depth dimension.
#'
#' @param data Object containing observation points. Can be a data frame,
#'   tibble, `sf` object with POINT geometries, matrix with column names, or
#'   named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of
#'   file paths, opened `ncdf4` object, list of opened objects, or a data frame
#'   with a file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with extracted values appended. Data frames,
#'   tibbles and `sf` objects retain their structure; matrices and named lists
#'   are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_2d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract2d(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export
extract2d <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "2d", ...)
}

#' Extract the surface layer from a 3D netCDF variable
#'
#' The surface layer is defined as the first available depth layer in the
#' netCDF file. This is usually the shallowest layer.
#'
#' @param data Object containing observation points. Can be a data frame,
#'   tibble, `sf` object with POINT geometries, matrix with column names, or
#'   named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of
#'   file paths, opened `ncdf4` object, list of opened objects, or a data frame
#'   with a file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with extracted values appended. Data frames,
#'   tibbles and `sf` objects retain their structure; matrices and named lists
#'   are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract3d_surface(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export
extract3d_surface <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "surface", ...)
}

#' Extract the bottom available layer from a 3D netCDF variable
#'
#' For each observation, this function extracts the full vertical profile at
#' the nearest longitude, latitude and time cell, then returns the deepest
#' non-missing value. This is useful when the bathymetry of the netCDF grid
#' means that deeper layers are missing over shallow areas.
#'
#' @param data Object containing observation points. Can be a data frame,
#'   tibble, `sf` object with POINT geometries, matrix with column names, or
#'   named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of
#'   file paths, opened `ncdf4` object, list of opened objects, or a data frame
#'   with a file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with extracted values appended. Data frames,
#'   tibbles and `sf` objects retain their structure; matrices and named lists
#'   are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract3d_bottom(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export
extract3d_bottom <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "bottom", ...)
}

#' Extract the nearest valid depth layer from a 3D netCDF variable
#'
#' For each observation, this function reads the full vertical profile and
#' finds the valid, non-missing netCDF depth layer closest to the observation
#' depth. This avoids returning `NA` when the geometrically closest depth layer
#' is invalid at that grid cell, for example below the seabed.
#'
#' @param data Object containing observation points. Can be a data frame,
#'   tibble, `sf` object with POINT geometries, matrix with column names, or
#'   named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of
#'   file paths, opened `ncdf4` object, list of opened objects, or a data frame
#'   with a file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with extracted values appended. Data frames,
#'   tibbles and `sf` objects retain their structure; matrices and named lists
#'   are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   depth = c(60, 100),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract3d_nearest(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export

extract3d_nearest <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "nearest", ...)
}

#' Extract nearest, surface and bottom values from a 3D netCDF variable
#'
#' For each observation, this function extracts the full vertical profile at
#' the nearest longitude, latitude and time cell. It then returns three summary
#' columns: nearest valid depth, surface and bottom.
#'
#' @param data Object containing observation points. Can be a data frame,
#'   tibble, `sf` object with POINT geometries, matrix with column names, or
#'   named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of
#'   file paths, opened `ncdf4` object, list of opened objects, or a data frame
#'   with a file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with three extracted columns: `nearest_*`,
#'   `surface_*` and `seabottom_*`. Data frames, tibbles and `sf` objects retain
#'   their structure; matrices and named lists are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   depth = c(60, 100),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract3d_all(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export

extract3d_all <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "all", ...)
}
