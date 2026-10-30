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


.decode_netcdf_time <- function(nc, time_dim, offsets = NULL) {
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

  values <- if (is.null(offsets)) {
    as.vector(.get_dim_values(nc, time_dim))
  } else {
    as.vector(offsets)
  }

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


# Read temporal metadata while the connection is open. Bounds describe
# represented intervals; their widths are not sampling cadences.
.read_temporal_metadata <- function(nc, time_dim, var, timestamps) {
  if (is.null(time_dim)) return(list(static = TRUE, variable = var))
  bounds_name <- .read_time_attribute(nc, time_dim, "bounds")
  climatology <- .read_time_attribute(nc, time_dim, "climatology")
  bounds_variable <- which(vapply(nc$var, function(x) {
    identical(x$name, bounds_name)
  }, logical(1)))
  bounds <- NULL
  bounds_status <- "absent"
  if (!is.null(bounds_name)) {
    bounds_status <- "unavailable"
    if (is.character(bounds_name) && length(bounds_name) == 1L &&
          !is.na(bounds_name) && length(bounds_variable) == 1L) {
      bounds <- tryCatch({
        raw <- ncdf4::ncvar_get(nc, bounds_name, collapse_degen = FALSE)
        axes <- vapply(
          nc$var[[bounds_variable]]$dim,
          function(x) x$name, character(1)
        )
        time_axis <- which(axes == time_dim)
        if (length(dim(raw)) != 2L || length(time_axis) != 1L) {
          stop("Unrecognised temporal bounds axes")
        }
        if (time_axis == 1L) raw <- t(raw)
        if (!identical(dim(raw), c(2L, length(timestamps)))) {
          stop("Unrecognised temporal bounds shape")
        }
        bound_units <- .read_time_attribute(nc, bounds_name, "units")
        if (!is.null(bound_units) &&
          !identical(
            bound_units, .read_time_attribute(nc, time_dim, "units")
          )) {
          stop("Temporal bounds use different units")
        }
        decoded <- .decode_netcdf_time(nc, time_dim, as.vector(raw))
        intervals <- matrix(as.numeric(decoded), nrow = 2L)
        if (any(intervals[1L, ] > intervals[2L, ])) {
          stop("Temporal bounds are not ordered")
        }
        intervals
      }, error = function(e) NULL)
      if (!is.null(bounds)) bounds_status <- "read"
    }
  }
  calendar <- .read_time_attribute(nc, time_dim, "calendar")
  if (is.null(calendar)) calendar <- "standard"
  list(
    static = FALSE, variable = var,
    units = .read_time_attribute(nc, time_dim, "units"),
    calendar = tolower(calendar),
    variable_units = .read_time_attribute(nc, var, "units"),
    cell_methods = .read_time_attribute(nc, var, "cell_methods"),
    bounds_name = bounds_name, bounds = bounds,
    bounds_status = bounds_status, climatology = climatology,
    timestamps = timestamps,
    cadence = .infer_temporal_cadence(timestamps)
  )
}

.infer_temporal_cadence <- function(timestamps, basis = "source") {
  values <- sort(unique(as.numeric(timestamps[!is.na(timestamps)])))
  result <- list(type = "unknown", seconds = NA_real_, months = NA_integer_,
                 basis = basis, n = length(values))
  if (length(values) < 2L) return(result)
  dates <- as.POSIXct(values, origin = "1970-01-01", tz = "UTC")
  parts <- as.POSIXlt(dates, tz = "UTC")
  month_index <- (parts$year + 1900L) * 12L + parts$mon
  month_steps <- diff(month_index)
  clock <- format(dates, "%H:%M:%OS6", tz = "UTC")
  day <- as.integer(format(dates, "%d", tz = "UTC"))
  next_month <- as.Date(format(as.Date(dates) + 32 - day, "%Y-%m-01"))
  month_end <- as.Date(dates) == next_month - 1
  calendar_regular <- all(month_steps > 0) &&
    length(unique(month_steps)) == 1L && length(unique(clock)) == 1L &&
    (length(unique(day)) == 1L || all(month_end))
  if (calendar_regular) {
    result$months <- as.integer(month_steps[1L])
    result$type <- if (result$months %% 12L == 0L) "yearly" else "monthly"
    return(result)
  }
  steps <- diff(values)
  tolerance <- max(1e-6, abs(steps[1L]) * 1e-9)
  if (all(abs(steps - steps[1L]) <= tolerance)) {
    result$type <- if (abs(steps[1L] - 86400) <= tolerance) "daily" else "fixed"
    result$seconds <- steps[1L]
  } else {
    result$type <- "irregular"
  }
  result
}

# Calendar shifts clamp the day to the final day of the destination month.
.shift_temporal_months <- function(value, months) {
  parts <- as.POSIXlt(value, tz = "UTC")
  index <- (parts$year + 1900L) * 12L + parts$mon + months
  first <- as.Date(sprintf("%04d-%02d-01", index %/% 12L, index %% 12L + 1L))
  next_index <- index + 1L
  next_first <- as.Date(sprintf(
    "%04d-%02d-01", next_index %/% 12L, next_index %% 12L + 1L
  ))
  last_day <- as.integer(next_first - first)
  as.POSIXct(first + min(parts$mday, last_day) - 1L, tz = "UTC") +
    parts$hour * 3600 + parts$min * 60 + parts$sec
}

.diagnose_temporal_selection <- function(provenance, requested, resolution,
                                         metadata, custom_seconds = NULL) {
  records <- list()
  for (column in names(provenance)) {
    selected <- provenance[[column]]
    for (i in seq_len(nrow(selected))) {
      source <- selected$source_index[i]
      timestamp <- selected$timestamp[i]

      if (is.na(source) || is.na(timestamp) || is.na(requested[i])) next
      info <- metadata[[source]]
      if (isTRUE(info$static)) next
      cadence <- info$cadence
      # A singleton may borrow an observed cadence from compatible sources
      # of the same variable. Never pool different calendars or units.
      if (identical(cadence$type, "unknown") && is.null(info$climatology)) {
        compatible <- vapply(metadata, function(x) {
          !isTRUE(x$static) && is.null(x$climatology) &&
            identical(x$variable, info$variable) &&
            identical(x$calendar, info$calendar) &&
            identical(x$variable_units, info$variable_units) &&
            identical(x$cell_methods, info$cell_methods)
        }, logical(1))
        combined <- unlist(lapply(metadata[compatible], function(x) {
          as.numeric(x$timestamps)
        }), use.names = FALSE)
        cadence <- .infer_temporal_cadence(
          as.POSIXct(combined, origin = "1970-01-01", tz = "UTC"),
          basis = "source_collection"
        )
      }
      delta <- as.numeric(difftime(timestamp, requested[i], units = "secs"))
      exceeds_cadence <- FALSE
      if (cadence$type %in% c("daily", "fixed")) {
        exceeds_cadence <- abs(delta) > cadence$seconds + 1e-6
      } else if (cadence$type %in% c("monthly", "yearly")) {
        earliest <- .shift_temporal_months(timestamp, -cadence$months)
        latest <- .shift_temporal_months(timestamp, cadence$months)
        exceeds_cadence <- requested[i] < earliest || requested[i] > latest
      }
      reason <- NA_character_
      if (resolution == "custom") {
        if (abs(delta) > custom_seconds) {
          reason <- "custom_interval_exceeded"
        }
      } else if (resolution %in% c("daily", "monthly", "yearly")) {
        pattern <- switch(resolution, daily = "%Y-%m-%d",
                          monthly = "%Y-%m", yearly = "%Y")
        if (format(timestamp, pattern, tz = "UTC") !=
              format(requested[i], pattern, tz = "UTC")) {
          reason <- "period_changed"
        }
      } else if (delta != 0) {
        if (!is.null(info$climatology) ||
              cadence$type %in% c("unknown", "irregular")) {
          reason <- "cadence_unknown"
        } else if (cadence$type == "daily" &&
                     format(timestamp, "%Y-%m-%d", tz = "UTC") !=
                       format(requested[i], "%Y-%m-%d", tz = "UTC")) {
          reason <- "period_changed"
        } else if (exceeds_cadence) {
          reason <- "cadence_exceeded"
        }
      }
      records[[length(records) + 1L]] <- data.frame(
        row = i, column = column, source_index = source,
        variable = selected$variable[i],
        requested = requested[i], selected = timestamp,
        difference_seconds = delta, resolution = resolution,
        cadence = cadence$type, cadence_basis = cadence$basis,
        cadence_seconds = cadence$seconds, cadence_months = cadence$months,
        custom_seconds = if (is.null(custom_seconds)) {
          NA_real_
        } else {
          custom_seconds
        },
        cadence_exceeded = exceeds_cadence,
        source_file = if (is.null(info$source_file)) {
          NA_character_
        } else {
          info$source_file
        },
        status = if (delta == 0) "matched" else "shifted",
        reason = reason, stringsAsFactors = FALSE
      )
    }
  }
  if (!length(records)) return(data.frame())
  do.call(rbind, records)
}

.validate_custom_time_interval <- function(resolution, interval) {
  if (resolution != "custom") {
    if (!is.null(interval)) {
      stop(
        "`custom_time_interval` requires time_resolution = 'custom'. ",
        "Set that resolution or remove custom_time_interval.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (!inherits(interval, "difftime") || length(interval) != 1L) {
    stop(
      "`custom_time_interval` must be one difftime value. ",
      "For example, use as.difftime(3, units = 'hours').",
      call. = FALSE
    )
  }
  seconds <- as.numeric(interval, units = "secs")
  if (is.na(seconds) || !is.finite(seconds) || seconds <= 0) {
    stop(
      "`custom_time_interval` must be finite and greater than zero. ",
      "Use a positive duration, such as as.difftime(3, units = 'hours').",
      call. = FALSE
    )
  }
  seconds
}

.temporal_duration_text <- function(seconds) {
  seconds <- abs(seconds)
  if (!is.finite(seconds)) return("an unknown interval")
  divisor <- if (seconds >= 3600) 3600 else if (seconds >= 60) 60 else 1
  unit <- if (divisor == 3600) {
    "hour"
  } else if (divisor == 60) {
    "minute"
  } else {
    "second"
  }
  amount <- seconds / divisor
  paste0(
    format(signif(amount, 6), trim = TRUE, scientific = FALSE), " ",
    unit, if (amount == 1) "" else "s"
  )
}

.temporal_cadence_text <- function(record) {
  type <- record$cadence
  if (type == "daily") return("daily (one observed step is 24 hours)")
  if (type == "fixed") {
    return(paste0("one observed step every ",
                  .temporal_duration_text(record$cadence_seconds)))
  }
  if (type %in% c("monthly", "yearly")) {
    months <- record$cadence_months
    number <- if (type == "yearly") months / 12 else months
    unit <- if (type == "yearly") "year" else "month"
    return(paste0(number, " calendar ", unit, if (number == 1) "" else "s"))
  }
  if (type == "irregular") return("irregular; no reliable sampling interval")
  "unknown; not enough information to infer a sampling interval"
}

# Combine depth-summary columns only when they describe the same selection.
.temporal_example_rows <- function(records) {
  keys <- setdiff(names(records), "column")
  examples <- unique(records[keys])
  examples$column <- vapply(seq_len(nrow(examples)), function(i) {
    same <- rep(TRUE, nrow(records))
    for (key in keys) {
      value <- examples[[key]][i]
      matches <- if (is.na(value)) {
        is.na(records[[key]])
      } else {
        !is.na(records[[key]]) & records[[key]] == value
      }
      same <- same & matches
    }
    paste(unique(records$column[same]), collapse = ", ")
  }, character(1))
  head(examples, 3L)
}

.format_temporal_example <- function(record) {
  reason <- record$reason
  status <- record$status
  delta <- record$difference_seconds
  location <- paste0(
    "Row ", record$row, "; variable: ", record$variable,
    "; output: ", record$column, "."
  )
  requested <- if (is.na(record$requested)) {
    "unavailable"
  } else {
    format(record$requested, "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
  }
  if (status == "no_value") {
    detail <- record[["detail"]]
    if (is.null(detail) || is.na(detail)) {
      detail <- paste0(
        "No source supplied a valid value under the requested settings.",
        " The result is NA; no selected timestamp or file can be assigned.",
        "\nCheck the observation, source coverage and matching limits."
      )
    } else {
      detail <- paste0(detail, "\nResult: NA.")
    }
    return(paste0(
      "No value was extracted.\n", location,
      "\nObservation: ", requested, "\n", detail
    ))
  }
  if (status == "not_assessable") {
    return(paste0(
      "The temporal comparison could not be completed.\n", location,
      "\nObservation: ", requested,
      "\nA requested or selected timestamp is unavailable."
    ))
  }
  title <- "A value with a different timestamp was selected."
  if (!is.na(reason)) {
    title <- switch(reason,
      period_changed = "A value from a different calendar period was selected.",
      cadence_exceeded = paste0(
        "The selected value is more than one sampling interval away."
      ),
      custom_interval_exceeded = paste0(
        "The selected value exceeds your custom time interval."
      ),
      cadence_unknown = paste0(
        "A different timestamp was selected, but its cadence is unknown."
      )
    )
  }
  explanation <- ""
  if (identical(reason, "period_changed")) {
    period <- if (record$resolution == "datetime") "day" else switch(
      record$resolution, daily = "day", monthly = "month", yearly = "year"
    )
    explanation <- if (record$resolution == "datetime") {
      paste0(
        "The dataset has a daily cadence. The selected value belongs to a ",
        "different day. This warning is based on the calendar day, ",
        "not only the number of elapsed hours."
      )
    } else {
      paste0(
        "The requested comparison uses ", record$resolution,
        " precision. The selected value belongs to a different ", period,
        ", so a warning is required even if the elapsed time is short."
      )
    }
    if (isTRUE(record$cadence_exceeded)) {
      explanation <- paste0(
        explanation,
        " It is also more than one observed sampling interval away."
      )
    }
  } else if (identical(reason, "custom_interval_exceeded")) {
    explanation <- paste0(
      "Your custom warning interval is ",
      .temporal_duration_text(record$custom_seconds),
      ". The selected value is further away than this interval.",
      " This warning does not reject the extracted value."
    )
  } else if (identical(reason, "cadence_exceeded")) {
    explanation <- paste0(
      "The time difference exceeds one observed sampling interval.",
      " The extracted value is returned, but it is further away in time ",
      "than the detected cadence."
    )
  } else if (identical(reason, "cadence_unknown")) {
    explanation <- paste0(
      "The available timestamps do not establish a reliable sampling interval.",
      " The time difference cannot be assessed against a detected cadence.",
      " Use time_resolution = 'custom' with custom_time_interval ",
      "if you know the appropriate warning interval."
    )
  } else if (record$resolution == "custom") {
    explanation <- paste0(
      "The difference is within your custom warning interval of ",
      .temporal_duration_text(record$custom_seconds),
      ". No temporal warning is required."
    )
  } else if (record$resolution %in% c("daily", "monthly", "yearly")) {
    explanation <- paste0(
      "The selected value remains within the requested ",
      record$resolution, " calendar period. No temporal warning is required."
    )
  } else {
    explanation <- paste0(
      "The difference is within one observed sampling interval.",
      " No temporal warning is required."
    )
  }
  file <- if (is.na(record$source_file)) {
    paste0("source ", record$source_index)
  } else {
    basename(gsub("\\", "/", record$source_file, fixed = TRUE))
  }
  direction <- if (delta > 0) {
    "later"
  } else if (delta < 0) {
    "earlier"
  } else {
    "unchanged"
  }
  paste0(
    title, "\n", location,
    "\nObservation: ", requested,
    "\nSelected data: ",
    format(record$selected, "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
    "\nDifference: ", .temporal_duration_text(delta), " ", direction, ".",
    "\nComparison resolution: ", record$resolution, ".",
    "\nDetected cadence: ", .temporal_cadence_text(record), ".",
    "\n", explanation, "\nFile: ", file
  )
}

# Formatting affects messages only, never matching or extracted values.
.selection_number <- function(value, digits = 2L, small = FALSE) {
  if (!is.finite(value)) return("unavailable")
  if (small && value > 0 && value < 10^-digits) {
    return(paste0("<", formatC(10^-digits, format = "f", digits = digits)))
  }
  rounded <- round(value, digits)
  if (rounded == 0) rounded <- 0
  text <- formatC(rounded, format = "f", digits = digits)
  sub("\\.$", "", sub("0+$", "", text))
}

.selection_point <- function(lon, lat) {
  paste0("(", .selection_number(lon, 5L), ", ",
         .selection_number(lat, 5L), ")")
}

# Depth axes stored as float32 read back slightly off their nominal values,
# so report comparisons allow a small relative tolerance.
.depth_tolerance <- function(requested) {
  1e-6 * max(1, abs(requested))
}

# Use the source's physical grid, including layers with missing variable values.
.depth_selection_reference <- function(depth_m, requested, chosen) {
  layers <- sort(as.numeric(depth_m))
  tol <- .depth_tolerance(requested)
  difference <- abs(chosen - requested)
  closest <- min(abs(layers - requested))
  outside <- requested < min(layers) - tol || requested > max(layers) + tol
  reference <- list(
    range = range(layers), outside = outside, difference = difference,
    closer_missing = difference > closest + tol,
    lower = NA_real_, upper = NA_real_, spacing = NA_real_,
    basis = "single_layer", status = "normal"
  )
  if (length(layers) == 1L) {
    if (difference > tol) reference$status <- "spacing_unknown"
    return(reference)
  }
  exact <- which(abs(layers - requested) <= tol)[1L]
  if (!is.na(exact)) {
    neighbours <- c(exact - 1L, exact + 1L)
    neighbours <- neighbours[neighbours >= 1L & neighbours <= length(layers)]
    nearest <- neighbours[which.min(abs(layers[neighbours] - requested))]
    pair <- sort(c(layers[exact], layers[nearest]))
    reference$basis <- "adjacent"
  } else if (outside) {
    pair <- if (requested < min(layers)) head(layers, 2L) else tail(layers, 2L)
    reference$basis <- "nearest_edge"
  } else {
    lower <- max(which(layers < requested))
    pair <- layers[c(lower, lower + 1L)]
    reference$basis <- "surrounding"
  }
  reference$lower <- pair[1L]
  reference$upper <- pair[2L]
  reference$spacing <- diff(pair)
  if (difference > reference$spacing + tol) {
    reference$status <- "distant_depth"
  }
  # Missing values explain the displacement only if a closer layer would
  # itself have been within the spacing.
  reference$closer_missing <- reference$closer_missing &&
    closest <= reference$spacing + tol
  reference
}

.format_depth_reference <- function(reference) {
  if (reference$basis == "single_layer") {
    return("Only one depth layer; vertical spacing cannot be assessed.")
  }
  label <- switch(reference$basis,
                  surrounding = "Surrounding layers: ",
                  adjacent = "Requested layer and nearest neighbour: ",
                  nearest_edge = "Nearest edge layers: ")
  comparison <- .selection_comparison(reference$difference, reference$spacing)
  text <- paste0(
    label, .selection_number(reference$lower, small = TRUE), "\u2013",
    .selection_number(reference$upper, small = TRUE), " m (",
    comparison$limit, " m apart)."
  )
  if (reference$outside) {
    text <- paste0(
      "Requested depth is outside the layer range (",
      .selection_number(reference$range[1L], small = TRUE), "\u2013",
      .selection_number(reference$range[2L], small = TRUE), " m).\n", text
    )
  }
  if (reference$status == "distant_depth") {
    text <- paste0(text, "\nSelection exceeds this spacing",
                   if (reference$closer_missing) {
                     " because closer layers have no valid value."
                   } else {
                     "."
                   })
  }
  text
}

# Why one source gave no value for one observation and output column.
.source_rejection <- function(candidate, i, column, data, method, depth_col) {
  horizontal <- candidate$spatial_matches
  covered <- isTRUE(horizontal$accepted[i])
  status <- candidate$extraction_status[i]
  if (!covered) {
    category <- "spatial"
    reason <- switch(horizontal$reason[i],
      outside_grid = "Point outside cell coverage (within matching).",
      distance_limit_exceeded = {
        comparison <- .selection_comparison(
          horizontal$distance_km[i], horizontal$max_distance_km
        )
        paste0("Nearest cell: ", comparison$value,
               " km away; allowed: ", comparison$limit, " km.")
      },
      "Observation coordinates are missing or invalid."
    )
  } else if (status == "date_missing") {
    category <- "temporal"
    reason <- paste0(
      "No observation date was supplied, but this source has a time ",
      "dimension; supply `date` or `date_col`."
    )
  } else if (status == "time_rejected") {
    category <- "temporal"
    reason <- "No timestamp satisfied the requested temporal settings."
  } else {
    category <- if (method == "2d") "spatial" else "depth"
    missing_depth <- method %in% c("nearest", "all") &&
      (method == "nearest" || column == candidate$columns[1L]) &&
      is.na(data[[depth_col]][i])
    reason <- if (missing_depth) {
      "Requested depth is missing; no nearest layer can be selected."
    } else {
      "No valid value in the selected cell/layer for this output."
    }
  }
  list(category = category, covered = covered, reason = reason,
       file = basename(candidate$metadata$source_file))
}

# Show each distinct reason once, with a few of the files that gave it.
.group_source_reasons <- function(reasons, files, max_files = 3L) {
  vapply(unique(reasons), function(reason) {
    hit <- unique(files[reasons == reason])
    paste0(
      reason, "\nFile", if (length(hit) > 1L) "s", ": ",
      paste(head(hit, max_files), collapse = ", "),
      if (length(hit) > max_files) {
        paste0(" and ", length(hit) - max_files, " more")
      }
    )
  }, character(1), USE.NAMES = FALSE)
}

# Collect report records. Normal selections are collected only for reports
# in full (TRUE) mode. Missing outputs are explained once: spatial and depth
# causes here, temporal causes in the time report (returned as `temporal`).
.selection_report_records <- function(outputs, provenance, data, method,
                                      depth_col, spatial_report = TRUE,
                                      depth_report = TRUE) {
  records <- list()
  explained <- list()
  temporal <- list()
  candidates <- lapply(outputs, attr, which = ".source_candidate",
                       exact = TRUE)
  wanted <- function(mode, normal) {
    !identical(mode, FALSE) && (isTRUE(mode) || !normal)
  }
  for (column in names(provenance)) {
    selected <- provenance[[column]]
    explained[[column]] <- rep(FALSE, nrow(data))
    temporal[[column]] <- rep(NA_character_, nrow(data))
    producers <- candidates[vapply(candidates, function(x) {
      column %in% x$columns
    }, logical(1))]
    for (i in seq_len(nrow(data))) {
      source <- selected$source_index[i]
      if (is.na(source)) {
        rejections <- lapply(producers, .source_rejection,
                             i = i, column = column, data = data,
                             method = method, depth_col = depth_col)
        categories <- vapply(rejections, `[[`, "", "category")
        covered <- vapply(rejections, `[[`, NA, "covered")
        # Sources whose grid never covers the point do not decide the cause.
        relevant <- if (any(covered)) categories[covered] else categories
        category <- if ("spatial" %in% relevant) {
          "spatial"
        } else if ("depth" %in% relevant) {
          "depth"
        } else {
          "temporal"
        }
        reasons <- .group_source_reasons(
          vapply(rejections, `[[`, "", "reason"),
          vapply(rejections, `[[`, "", "file")
        )
        if (category == "temporal") {
          temporal[[column]][i] <- paste(reasons, collapse = "\n")
          next
        }
        explained[[column]][i] <- TRUE
        mode <- if (category == "spatial") spatial_report else depth_report
        if (!wanted(mode, normal = FALSE)) next
        records[[length(records) + 1L]] <- list(
          row = i, column = column,
          variable = paste(unique(vapply(producers, `[[`, "", "variable")),
                           collapse = ", "),
          category = category, status = "missing", reasons = reasons
        )
        next
      }
      candidate <- candidates[[source]]
      base <- list(row = i, column = column, variable = candidate$variable,
                   source_index = source,
                   source_file = candidate$metadata$source_file)
      matches <- candidate$spatial_matches
      status <- if (isTRUE(matches$inside_selected_cell[i])) {
        "normal"
      } else if (isTRUE(matches$inside_coverage[i])) {
        "outside_cell"
      } else {
        "outside_coverage"
      }
      if (wanted(spatial_report, normal = status == "normal")) {
        records[[length(records) + 1L]] <- c(base, list(
          category = "spatial", status = status,
          horizontal = .horizontal_match_row(matches, i)
        ))
      }
      if (method == "2d" || identical(depth_report, FALSE)) next
      chosen <- candidate$depth_selected[[column]][i]
      if (!is.finite(chosen)) next
      kind <- if (method == "all") {
        c("nearest", "surface", "bottom")[match(column, candidate$columns)]
      } else {
        method
      }
      requested <- if (kind == "nearest") data[[depth_col]][i] else NA_real_
      reference <- if (kind == "nearest") {
        .depth_selection_reference(
          candidate$depth_axis$depth_m, requested, chosen
        )
      }
      status <- if (is.null(reference)) "normal" else reference$status
      if (!wanted(depth_report, normal = status == "normal")) next
      records[[length(records) + 1L]] <- c(base, list(
        category = "depth", status = status,
        kind = kind, requested = requested, chosen = chosen,
        depth_reference = reference
      ))
    }
  }
  list(records = records, explained = explained, temporal = temporal)
}

.format_selection_example <- function(record) {
  location <- paste0(
    "Row ", record$row, " \u00b7 netCDF variable: ", record$variable,
    " \u00b7 output column: ", record$column
  )
  if (record$status == "missing") {
    return(paste(c(location, record$reasons, "Result: NA."), collapse = "\n"))
  }
  file <- paste0("File: ", basename(record$source_file))
  if (record$category == "spatial") {
    horizontal <- record$horizontal
    label <- if (record$status == "normal") {
      "containing cell centre: "
    } else {
      "nearest cell centre: "
    }
    distance <- paste0("Distance: ",
                       .selection_number(horizontal$distance_km, small = TRUE),
                       " km")
    if (!is.null(horizontal$max_distance_km)) {
      comparison <- .selection_comparison(
        horizontal$distance_km, horizontal$max_distance_km
      )
      distance <- paste0("Distance: ", comparison$value,
                         " km; allowed: ", comparison$limit, " km")
    }
    return(paste0(
      location, "\nRequested: ",
      .selection_point(horizontal$requested_lon, horizontal$requested_lat),
      " \u2192 ", label,
      .selection_point(horizontal$selected_lon, horizontal$selected_lat),
      "\n", distance, ".\n", file
    ))
  }
  chosen <- .selection_number(record$chosen, small = TRUE)
  detail <- if (record$kind == "nearest") {
    delta <- record$chosen - record$requested
    comparison <- .selection_comparison(
      abs(delta), record$depth_reference$spacing
    )
    paste0(
      "Requested: ", .selection_number(record$requested, small = TRUE),
      " m \u2192 selected: ", chosen, " m",
      if (abs(delta) <= .depth_tolerance(record$requested)) {
        " (exact match)."
      } else {
        paste0(" (", comparison$value, " m ",
               if (delta < 0) "shallower)." else "deeper).")
      }, "\nNo vertical interpolation."
    )
  } else if (record$kind == "surface") {
    paste0("Shallowest layer selected: ", chosen, " m.")
  } else {
    paste0("Deepest layer with a valid value: ", chosen,
           " m. This does not establish the seabed depth.")
  }
  if (record$kind == "nearest" &&
        abs(record$chosen - record$requested) >
          .depth_tolerance(record$requested)) {
    detail <- paste(detail, .format_depth_reference(record$depth_reference),
                    sep = "\n")
  }
  paste(location, detail, file, sep = "\n")
}

.selection_count <- function(count, singular, plural = paste0(singular, "s")) {
  paste(count, if (count == 1L) singular else plural)
}

.show_selection_reports <- function(outputs, provenance, data, method,
                                    depth_col, spatial_report, depth_report) {
  collected <- .selection_report_records(outputs, provenance, data, method,
                                         depth_col, spatial_report,
                                         depth_report)
  modes <- list(spatial = spatial_report, depth = depth_report)
  categories <- vapply(collected$records, `[[`, "", "category")
  for (category in names(modes)) {
    mode <- modes[[category]]
    records <- collected$records[categories == category]
    if (!length(records)) next
    statuses <- vapply(records, function(x) x$status, character(1))
    # Stable ordering puts missing results first, then deviations, then normal.
    priority <- ifelse(statuses == "missing", 0L,
                       ifelse(statuses == "normal", 2L, 1L))
    records <- records[order(priority)]
    statuses <- statuses[order(priority)]
    count <- length(unique(vapply(records, function(x) x$row, integer(1))))
    heading <- switch(
      category, spatial = "Spatial selection: ", depth = "Depth selection: "
    )
    action <- if (all(statuses == "outside_coverage")) {
      " outside cell coverage used the nearest cell"
    } else if (all(statuses == "outside_cell")) {
      " used a nearest cell that does not contain the requested point"
    } else if (all(statuses == "distant_depth")) {
      " used a depth more than one local layer spacing away"
    } else if (all(statuses == "spacing_unknown")) {
      " used a different depth; vertical spacing is unavailable"
    } else if (all(statuses == "missing")) {
      " had no value extracted"
    } else if (all(statuses == "normal")) {
      " matched the requested selection"
    } else {
      " with selection details"
    }
    examples <- head(records, 5L)
    detail <- vapply(examples, .format_selection_example, character(1))
    omitted <- length(records) - length(examples)
    message(structure(list(
      message = paste0(
        heading, .selection_count(count, "observation"), action, ".\n",
        .selection_count(length(records), "output value"),
        if (isTRUE(mode)) " reported" else " affected",
        if (any(statuses == "missing")) {
          paste0("; ", sum(statuses == "missing"), " missing")
        } else {
          ""
        }, ".",
        if (omitted > 0L) {
          paste0(" Showing ", length(examples), " of ", length(records),
                 " cases.")
        } else {
          ""
        }, "\n\n", paste(detail, collapse = "\n\n"),
        if (omitted > 0L) paste0("\n", omitted, " further examples omitted.")
        else "", "\n"
      ), call = NULL, records = records
    ), class = c(paste0("ocean3d_", category, "_report"), "message",
                 "condition")))
  }
  collected[c("explained", "temporal")]
}


.warn_temporal_selection <- function(diagnostics) {
  if (!nrow(diagnostics)) return(invisible(NULL))
  flagged <- diagnostics[!is.na(diagnostics$reason), , drop = FALSE]
  if (!nrow(flagged)) return(invisible(NULL))
  examples <- .temporal_example_rows(flagged)
  detail <- vapply(seq_len(nrow(examples)), function(i) {
    .format_temporal_example(examples[i, , drop = FALSE])
  }, character(1))
  warning(structure(list(
    message = paste0(
      "Temporal warning: ", length(unique(flagged$row)),
      " observation(s) need attention.\n\n",
      paste(detail, collapse = "\n\n")
    ),
    call = NULL, diagnostics = flagged
  ), class = c("ocean3d_time_warning", "warning", "condition")))
  invisible(NULL)
}

# `missing` comes from .show_selection_reports(): `explained` marks missing
# outputs already reported as spatial or depth misses, and `temporal` holds
# the per-source reasons for missing outputs with a temporal cause.
.show_time_report <- function(provenance, requested, metadata,
                              diagnostics, missing = NULL) {
  rows <- list()
  diagnostic_keys <- if (nrow(diagnostics)) {
    paste(diagnostics$row, diagnostics$column, sep = ":")
  } else {
    character()
  }
  for (column in names(provenance)) {
    selected <- provenance[[column]]
    explained <- missing$explained[[column]]
    if (is.null(explained)) explained <- rep(FALSE, nrow(selected))
    temporal <- missing$temporal[[column]]
    if (is.null(temporal)) temporal <- rep(NA_character_, nrow(selected))
    for (i in seq_len(nrow(selected))) {
      source_index <- selected$source_index[i]
      if (is.na(source_index) && explained[i]) next
      hit <- match(paste(i, column, sep = ":"), diagnostic_keys)
      if (!is.na(hit)) {
        row <- diagnostics[hit, , drop = FALSE]
        row$detail <- NA_character_
        rows[[length(rows) + 1L]] <- row
        next
      }
      if (!is.na(source_index) && isTRUE(metadata[[source_index]]$static)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        row = i, column = column, source_index = source_index,
        variable = if (is.na(selected$variable[i])) {
          column
        } else {
          selected$variable[i]
        },
        requested = requested[i], selected = selected$timestamp[i],
        difference_seconds = NA_real_, resolution = NA_character_,
        cadence = "unknown", cadence_basis = NA_character_,
        cadence_seconds = NA_real_, cadence_months = NA_integer_,
        custom_seconds = NA_real_, cadence_exceeded = FALSE,
        source_file = NA_character_,
        status = if (is.na(source_index)) "no_value" else "not_assessable",
        reason = NA_character_,
        detail = if (is.na(source_index)) temporal[i] else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(invisible(NULL))
  report <- do.call(rbind, rows)
  details <- report$detail
  report$detail <- NULL
  # Warning examples already explain flagged selections.
  # Do not print them twice.
  information <- report$status != "matched" & is.na(report$reason)
  if (!any(information)) return(invisible(report))
  shown <- report[information, , drop = FALSE]
  shown$detail <- details[information]
  examples <- .temporal_example_rows(shown)
  detail <- vapply(seq_len(nrow(examples)), function(i) {
    .format_temporal_example(examples[i, , drop = FALSE])
  }, character(1))
  message(structure(list(
    message = paste0(
      "Temporal report: ", length(unique(shown$row)),
      " observation(s) have additional information.\n\n",
      paste(detail, collapse = "\n\n"), "\n"
    ),
    call = NULL, report = report
  ), class = c("ocean3d_time_report", "message", "condition")))
  invisible(report)
}


.output_column_names <- function(method, var, prefix = NULL) {
  prefix <- if (is.null(prefix)) var else prefix
  if (method == "all") {
    return(paste0(c("nearest_", "surface_", "seabottom_"), prefix))
  }
  prefix
}

.source_output_columns <- function(output, original_names) {
  candidate <- attr(output, ".source_candidate", exact = TRUE)
  if (!is.null(candidate$columns)) return(candidate$columns)
  # Compatibility with existing internal merge tests.
  setdiff(names(output), original_names)
}

.plan_source_outputs <- function(sources, var, method, output_prefix,
                                 protected_names) {
  if (!length(sources)) {
    stop("At least one netCDF source is required.", call. = FALSE)
  }
  plan <- lapply(sources, function(source) {
    variable <- .resolve_source_var(source, var)
    prefix <- .resolve_output_prefix(
      method = method, var = variable, output_prefix = output_prefix,
      auto_var = is.null(var)
    )
    list(
      variable = variable, prefix = prefix,
      columns = .output_column_names(method, variable, prefix)
    )
  })
  mapping <- do.call(rbind, lapply(plan, function(item) {
    data.frame(
      column = item$columns, variable = item$variable,
      stringsAsFactors = FALSE
    )
  }))
  protected <- intersect(mapping$column, protected_names)
  if (length(protected)) {
    stop(
      "Output columns would overwrite extraction inputs: ",
      toString(unique(protected)), ". Choose another output name.",
      call. = FALSE
    )
  }
  variables_by_column <- split(mapping$variable, mapping$column)
  collisions <- names(Filter(
    function(x) length(unique(x)) > 1L, variables_by_column
  ))
  if (length(collisions)) {
    stop(
      "Different variables cannot share these output columns: ",
      toString(collisions), ". Use separate output names.",
      call. = FALSE
    )
  }
  plan
}

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
.great_circle_km <- function(lon1, lat1, lon2, lat2, radius_km = 6371.01) {
  p <- pi / 180
  h <- sin((lat2 - lat1) * p / 2)^2 +
    cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * radius_km * asin(pmin(1, sqrt(h)))
}

# Distances that differ by less than this are treated as ties, so rounding
# noise cannot choose between cells that are equally close (for example a
# point on a shared cell edge). Ties keep the first cell in storage order.
.distance_tie_km <- 1e-9

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

.selection_comparison <- function(value, limit) {
  shown <- c(.selection_number(value, small = TRUE),
             .selection_number(limit, small = TRUE))
  if (is.finite(value) && is.finite(limit) && value != limit &&
        shown[1L] == shown[2L]) {
    for (digits in 3:15) {
      shown <- c(.selection_number(value, digits, small = TRUE),
                 .selection_number(limit, digits, small = TRUE))
      if (shown[1L] != shown[2L]) break
    }
    if (shown[1L] == shown[2L]) {
      shown <- format(c(value, limit), digits = 17L, scientific = FALSE,
                      trim = TRUE)
    }
  }
  list(value = shown[1L], limit = shown[2L])
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
  units <- tolower(trimws(units))
  factor <- switch(
    units,
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

# Select layers by physical depth, independently of storage order. Returns
# the selected values and the physical depth (m) of each selected layer.
.summarise_depth_profile <- function(profile, depth_m, target_depth) {
  index <- c(nearest = NA_integer_, surface = which.min(depth_m),
             seabottom = NA_integer_)
  valid <- which(!is.na(profile))
  if (length(valid)) {
    index["seabottom"] <- valid[which.max(depth_m[valid])]
    if (!is.na(target_depth)) {
      # A tie chooses the shallower valid layer.
      valid <- valid[order(depth_m[valid])]
      index["nearest"] <- valid[which.min(abs(depth_m[valid] - target_depth))]
    }
  }
  list(
    value = stats::setNames(profile[index], names(index)),
    depth_m = stats::setNames(depth_m[index], names(index))
  )
}


.unwrap_longitude <- function(values) {
  if (length(values) < 2L) return(values)
  original_steps <- diff(values)
  steps <- (original_steps + 180) %% 360 - 180
  ambiguous <- abs(abs(steps) - 180) < 1e-8
  steps[ambiguous] <- original_steps[ambiguous]
  values[1L] + c(0, cumsum(steps))
}

# Angle units differ in spelling (degrees_east, degree_east, degrees_E, ...).
.same_spatial_units <- function(a, b) {
  a <- tolower(trimws(a))
  b <- tolower(trimws(b))
  if (identical(a, b)) return(TRUE)
  angle <- "^degrees?(_?(east|north|e|n))?$"
  grepl(angle, a) && grepl(angle, b)
}

# Returns a 2 x n matrix of cell edges, or NULL when the file has no usable
# bounds. Unusable bounds are reported and edges are then inferred from the
# coordinate centres, as when no bounds are declared.
.read_spatial_bounds <- function(nc, coordinate, values) {
  name <- .read_time_attribute(nc, coordinate, "bounds")
  if (is.null(name)) return(NULL)
  ignore <- function(problem) {
    warning(
      "Spatial bounds for `", coordinate, "` ", problem,
      "; inferring cell edges from the coordinate centres instead.",
      call. = FALSE
    )
    NULL
  }
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    return(ignore("have an invalid bounds reference"))
  }
  if (!name %in% names(nc$var)) {
    return(ignore(paste0("were not found (variable `", name, "`)")))
  }
  variable <- nc$var[[name]]
  axes <- vapply(variable$dim, function(x) x$name, character(1))
  coordinate_axis <- which(axes == coordinate)
  raw <- ncdf4::ncvar_get(nc, name, collapse_degen = FALSE)
  if (length(dim(raw)) != 2L || length(coordinate_axis) != 1L) {
    return(ignore("have unsupported dimensions"))
  }
  if (coordinate_axis == 1L) raw <- t(raw)
  if (!identical(dim(raw), c(2L, length(values))) ||
        any(!is.finite(raw))) {
    return(ignore("contain missing or invalid values"))
  }
  coordinate_units <- .read_time_attribute(nc, coordinate, "units")
  bounds_units <- .read_time_attribute(nc, name, "units")
  if (!is.null(bounds_units) && !is.null(coordinate_units) &&
        !.same_spatial_units(bounds_units, coordinate_units)) {
    stop(
      "Spatial bounds and coordinates use different units (`",
      bounds_units, "` and `", coordinate_units, "`).",
      call. = FALSE
    )
  }
  raw
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
    stop(
      "Spatial coordinates must form a monotonic axis (as CF requires). ",
      "Sort the file along this dimension before extraction.",
      call. = FALSE
    )
  }
  explicit <- .read_spatial_bounds(nc, coordinate, values)
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

.prepare_time_input <- function(
  x,
  time_resolution = c(
    "auto", "datetime", "daily", "monthly", "yearly"
  ),
  input_name = "time input",
  ids = NULL
) {
  time_resolution <- match.arg(time_resolution)
  present <- !is.na(x)
  date_only <- rep(FALSE, length(x))
  unzoned <- rep(FALSE, length(x))

  fail <- function(bad, reason) {
    shown <- head(bad, 5L)
    original <- encodeString(as.character(x[shown]), quote = "\"")

    identifiers <- if (is.null(ids)) {
      rep("", length(shown))
    } else {
      paste0(
        " (id: ",
        encodeString(as.character(ids[shown]), quote = "\""),
        ")"
      )
    }

    examples <- paste0(
      "row/position ", shown, identifiers, ": ", original
    )
    omitted <- length(bad) - length(shown)

    stop(
      reason,
      "\nInput: ", input_name,
      "\nAffected observations: ", length(bad), ".",
      "\n", paste(examples, collapse = "\n"),
      if (omitted > 0L) {
        paste0("\nShowing ", length(shown), " of ", length(bad), " cases.")
      } else {
        ""
      },
      "\nUse Date, POSIXct/POSIXlt, YYYY-MM-DD, ",
      "YYYY-MM-DD HH:MM or YYYY-MM-DD HH:MM:SS.",
      "\nText timestamps are interpreted in UTC. ",
      "For another time zone, supply a POSIXct/POSIXlt object.",
      call. = FALSE
    )
  }

  if (inherits(x, "POSIXt")) {
    # Keep the instant; all calendar comparisons below explicitly use UTC.
    values <- as.POSIXct(x)
    attr(values, "tzone") <- "UTC"
    has_time <- TRUE
  } else if (inherits(x, "Date")) {
    values <- as.POSIXct(x, tz = "UTC")
    has_time <- FALSE
  } else if (is.atomic(x) && all(is.na(x))) {
    values <- as.POSIXct(rep(NA_character_, length(x)), tz = "UTC")
    has_time <- FALSE
  } else if (is.character(x)) {
    date_pattern <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}"
    clock_pattern <- paste0(date_pattern, "[ T][0-9]{2}:[0-9]{2}")
    date_only <- present & grepl(paste0(date_pattern, "$"), x)
    minutes <- present & grepl(paste0(clock_pattern, "Z?$"), x)
    seconds <- present & grepl(paste0(clock_pattern, ":[0-9]{2}Z?$"), x)
    datetime <- minutes | seconds
    unzoned <- datetime & !grepl("Z$", x)
    recognised <- date_only | datetime

    text <- rep(NA_character_, length(x))
    text[recognised] <- sub(
      "T", " ", sub("Z$", "", x[recognised]), fixed = TRUE
    )
    text[minutes] <- paste0(text[minutes], ":00")
    text[date_only] <- paste0(text[date_only], " 00:00:00")
    values <- as.POSIXct(strptime(text, "%Y-%m-%d %H:%M:%S", tz = "UTC"))
    round_trip <- format(values, "%Y-%m-%d %H:%M:%S", tz = "UTC")
    bad <- which(present & (!recognised | is.na(values) | round_trip != text))
    if (length(bad)) {
      fail(bad, "Time values could not be converted. Invalid date/time.")
    }

    if (time_resolution == "auto" && any(date_only) && any(datetime)) {
      fail(which(present), paste0(
        "Time input mixes dates with and without hours. ",
        "Set time_resolution explicitly or use a consistent format."
      ))
    }
    has_time <- any(datetime)
  } else {
    fail(which(present),
         "Time values could not be converted. Unsupported input type.")
  }

  bad <- which(present & !is.finite(as.numeric(values)))
  if (length(bad)) {
    fail(bad, "Time values could not be converted. Invalid date/time.")
  }
  if (time_resolution == "auto") {
    time_resolution <- if (has_time) "datetime" else "daily"
  }
  if (time_resolution == "datetime" && any(date_only)) {
    fail(which(date_only), paste0(
      "time_resolution = 'datetime' requires hours for every ",
      "non-missing observation."
    ))
  }
  if (time_resolution == "datetime" && !has_time && any(present)) {
    fail(which(present),
         "time_resolution = 'datetime' requires observation hours.")
  }
  list(values = values, resolution = time_resolution,
       unzoned_rows = which(unzoned))
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


# A source may offer a temporal candidate without implying a coverage interval.
# Matching and tolerance are handled by .match_time_index().
.temporal_source_eligible <- function(time_values, target, resolution) {
  if (is.null(time_values)) return(TRUE)
  if (length(time_values) == 0L || is.na(target)) return(FALSE)

  available <- .time_comparison_values(time_values, resolution)
  requested <- .time_comparison_values(target, resolution)
  !is.na(requested) && any(!is.na(available))
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

  sc <- .build_start_count(
    dims$var_dims, index_list,
    singleton_dims = dims$singleton_dims
  )

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
                              time_resolution = "auto",
                              multi_source = FALSE,
                              spatial_match = "within",
                              max_distance_km = NULL,
                              horizontal_cache = NULL,
                              direct_input = FALSE) {
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

  if (!is.null(dims$time)) {
    named_column <- is.character(date_col) &&
      length(date_col) == 1L &&
      !is.na(date_col) &&
      nzchar(date_col)

    if (!(named_column && date_col %in% names(data))) {
      correction <- if (direct_input) {
        "Supply the `date` argument."
      } else if (named_column) {
        paste0(
          "Column `", date_col, "` is missing from `data`. ",
          "Supply that column or set `date_col` to an existing date column."
        )
      } else {
        "Set `date_col` to an existing date column in `data`."
      }

      stop(
        "netCDF variable `", var, "` has a time dimension.\n",
        "File: ", basename(nc$filename), "\n",
        correction,
        "\nAn explicit NA date is allowed; a missing date input is not.",
        call. = FALSE
      )
    }
  }

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
  spatial_axes <- list(
    lon = .prepare_spatial_axis(
      nc, dims$lon, lon_values, longitude = TRUE
    ),
    lat = .prepare_spatial_axis(
      nc, dims$lat, lat_values, longitude = FALSE
    )
  )
  depth_values <- .get_dim_values(nc, dims$depth)
  depth_axis <- if (method == "2d") {
    NULL
  } else {
    .read_depth_axis(nc, dims$depth, depth_values)
  }
  time_values <- .decode_netcdf_time(nc, dims$time)
  temporal_metadata <- .read_temporal_metadata(nc, dims$time, var, time_values)
  temporal_metadata$source_file <- if (source$type == "path") {
    source$source
  } else {
    nc$filename
  }
  resolved_time_resolution <- time_resolution

  if (!is.null(dims$time) && !is.null(date_col)) {
    prepared <- .prepare_time_input(
      data[[date_col]], time_resolution,
      input_name = if (direct_input) {
        "`date`"
      } else {
        paste0("`data$", date_col, "`")
      },
      ids = if (is.null(id_col)) NULL else data[[id_col]]
    )
    date_values <- prepared$values
    resolved_time_resolution <- prepared$resolution
    if (any(!is.na(date_values))) {
      .check_time_resolution(
        time_values, resolved_time_resolution
      )
    }
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

  horizontal_matches <- .cached_horizontal_matches(
    horizontal_cache,
    key = list(lon_values, lat_values, spatial_axes),
    compute = function() {
      .match_horizontal_cells(
        lon_axis = spatial_axes$lon, lat_axis = spatial_axes$lat,
        lon_values = lon_values, lat_values = lat_values,
        lon_targets = data[[lon_col]], lat_targets = data[[lat_col]],
        spatial_match = spatial_match, max_distance_km = max_distance_km
      )
    }
  )
  accepted <- horizontal_matches$accepted
  extraction_status <- ifelse(accepted, "time_rejected", "spatial_rejected")
  if (!is.null(time_values)) {
    extraction_status[accepted & is.na(date_values)] <- "date_missing"
  }
  depth_selected <- lapply(.output_column_names(method, var, output_prefix),
                           function(x) rep(NA_real_, nrow(data)))
  names(depth_selected) <- .output_column_names(method, var, output_prefix)
  candidate_distance <- rep(Inf, nrow(data))
  candidate_time <- rep(NA_real_, nrow(data))
  candidate_index <- rep(NA_integer_, nrow(data))

  for (i in seq_len(nrow(data))) {
    if (isTRUE(verbose) && i %% 100 == 0) {
      message("Extracting row ", i, " of ", nrow(data))
    }
    if (!accepted[i]) next
    lon_idx <- horizontal_matches$lon_idx[i]
    lat_idx <- horizontal_matches$lat_idx[i]

    time_idx <- NULL
    if (!is.null(time_values)) {

      # Single-timestamp sources can offer nearest candidates too.
      # Compare their distances across sources during merging.
      if (isTRUE(multi_source) &&
        !.temporal_source_eligible(
          time_values = time_values,
          target = date_values[i],
          resolution = resolved_time_resolution
        )) {
        next
      }

      time_idx <- .match_time_index(
        values = time_values,
        target = date_values[i],
        time_match = time_match,
        max_time_diff = max_time_diff,
        time_resolution = resolved_time_resolution
      )

      if (is.na(time_idx)) next

      available <- .time_comparison_values(
        time_values[time_idx], resolved_time_resolution
      )
      requested <- .time_comparison_values(
        date_values[i], resolved_time_resolution
      )
      candidate_distance[i] <- abs(as.numeric(
        difftime(available, requested, units = "days")
      ))
      candidate_time[i] <- as.numeric(time_values[time_idx])
      candidate_index[i] <- time_idx
    }

    if (is.null(time_values)) {
      # Static sources tie on time. The merge then keeps input order, or the
      # geographically closer cell under spatial_match = "nearest".
      candidate_distance[i] <- 0
    }

    extraction_status[i] <- "value_missing"
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
        depth_start = depth_axis$surface_index,
        depth_count = 1L
      )
      out[[prefix]][i] <- as.numeric(val)[1L]
      depth_selected[[prefix]][i] <-
        depth_axis$depth_m[depth_axis$surface_index]
    }

    if (method %in% c("nearest", "bottom", "all")) {
      # Read once in the original file order; select by physical depth.
      profile <- as.numeric(.extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = length(depth_values)
      ))
      depth_target <- if (method == "bottom") NA_real_ else data[[depth_col]][i]
      summary <- .summarise_depth_profile(
        profile, depth_axis$depth_m, depth_target
      )
      kinds <- switch(method, all = c("nearest", "surface", "seabottom"),
                      bottom = "seabottom", nearest = "nearest")
      columns <- if (method == "all") paste0(kinds, "_", prefix) else prefix
      for (k in seq_along(kinds)) {
        out[[columns[k]]][i] <- summary$value[[kinds[k]]]
        depth_selected[[columns[k]]][i] <- summary$depth_m[[kinds[k]]]
      }
    }
  }
  attr(out, ".source_candidate") <- list(
    distance_days = candidate_distance,
    distance_km = ifelse(accepted, horizontal_matches$distance_km, Inf),
    selected_lon = ifelse(accepted, horizontal_matches$selected_lon, NA_real_),
    selected_lat = ifelse(accepted, horizontal_matches$selected_lat, NA_real_),
    spatial_matches = horizontal_matches,
    depth_axis = depth_axis,
    depth_selected = depth_selected,
    extraction_status = extraction_status,
    timestamp = candidate_time,
    time_index = candidate_index,
    variable = var,
    metadata = temporal_metadata,
    columns = .output_column_names(method, var, output_prefix)
  )
  out
}

# Sources that share one grid (for example a stack of daily files) give the
# same horizontal matches, so compute them once per distinct grid.
.cached_horizontal_matches <- function(cache, key, compute) {
  if (is.null(cache)) return(compute())
  for (entry in cache$entries) {
    if (identical(entry$key, key)) return(entry$matches)
  }
  matches <- compute()
  cache$entries[[length(cache$entries) + 1L]] <- list(
    key = key, matches = matches
  )
  matches
}

.merge_source_outputs <- function(
  outputs,
  original_names,
  return_provenance = FALSE,
  spatial_match = "within"
) {
  if (length(outputs) == 1L) {
    out <- outputs[[1L]]
    candidate <- attr(out, ".source_candidate", exact = TRUE)
    output_columns <- .source_output_columns(out, original_names)
    attr(out, ".source_candidate") <- NULL
    if (isTRUE(return_provenance)) {
      provenance <- list()
      for (col in output_columns) {
        valid <- !is.na(out[[col]]) & is.finite(candidate$distance_days)
        provenance[[col]] <- data.frame(
          source_index = ifelse(valid, 1L, NA_integer_),
          variable = ifelse(valid, candidate$variable, NA_character_),
          time_index = ifelse(valid, candidate$time_index, NA_integer_),
          timestamp = as.POSIXct(
            ifelse(valid, candidate$timestamp, NA_real_),
            origin = "1970-01-01", tz = "UTC"
          ),
          distance_days = ifelse(valid, candidate$distance_days, NA_real_),
          distance_km = ifelse(valid, candidate$distance_km, NA_real_),
          selected_lon = ifelse(valid, candidate$selected_lon, NA_real_),
          selected_lat = ifelse(valid, candidate$selected_lat, NA_real_),
          stringsAsFactors = FALSE
        )
      }
      return(list(data = out, provenance = provenance))
    }
    return(out)
  }

  out <- outputs[[1L]][original_names]
  attr(out, ".source_candidate") <- NULL
  new_cols <- unique(unlist(lapply(
    outputs, .source_output_columns, original_names = original_names
  )))
  provenance <- list()

  for (col in new_cols) {
    values <- rep(NA_real_, nrow(out))
    best_distance <- rep(Inf, nrow(out))
    best_spatial_distance <- rep(Inf, nrow(out))
    selected_lon <- rep(NA_real_, nrow(out))
    selected_lat <- rep(NA_real_, nrow(out))
    selected_source <- rep(NA_integer_, nrow(out))
    selected_time <- rep(NA_real_, nrow(out))
    selected_index <- rep(NA_integer_, nrow(out))
    selected_variable <- rep(NA_character_, nrow(out))

    for (j in seq_along(outputs)) {
      x <- outputs[[j]]
      produced_columns <- .source_output_columns(x, original_names)
      if (!col %in% produced_columns) next
      candidate <- attr(x, ".source_candidate", exact = TRUE)
      if (is.null(candidate)) {
        stop("Missing internal source-candidate metadata.", call. = FALSE)
      }
      distance <- candidate$distance_days
      spatial_distance <- candidate$distance_km
      better_time <- distance < best_distance
      same_time <- distance == best_distance
      # With "within" matching every eligible cell contains the point, so
      # input order breaks temporal ties. With "nearest" matching the closer
      # cell centre wins first; the tolerance keeps float rounding of
      # otherwise identical grids from overriding input order.
      better_space <- spatial_match == "nearest" &
        spatial_distance < best_spatial_distance - 1e-3
      replace <- !is.na(x[[col]]) & is.finite(distance) &
        is.finite(spatial_distance) &
        (is.na(values) | better_time | (same_time & better_space))
      values[replace] <- x[[col]][replace]
      best_distance[replace] <- distance[replace]
      best_spatial_distance[replace] <- spatial_distance[replace]
      selected_lon[replace] <- candidate$selected_lon[replace]
      selected_lat[replace] <- candidate$selected_lat[replace]
      selected_source[replace] <- j
      selected_time[replace] <- candidate$timestamp[replace]
      selected_index[replace] <- candidate$time_index[replace]
      selected_variable[replace] <- candidate$variable
    }

    out[[col]] <- values
    provenance[[col]] <- data.frame(
      source_index = selected_source,
      variable = selected_variable,
      time_index = selected_index,
      timestamp = as.POSIXct(
        selected_time, origin = "1970-01-01", tz = "UTC"
      ),
      distance_days = ifelse(
        is.na(selected_source), NA_real_, best_distance
      ),
      distance_km = ifelse(
        is.na(selected_source), NA_real_, best_spatial_distance
      ),
      selected_lon = selected_lon,
      selected_lat = selected_lat,
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(return_provenance)) {
    return(list(data = out, provenance = provenance))
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
#'   "YYYY-MM-DD", "YYYY-MM-DD HH:MM" or "YYYY-MM-DD HH:MM:SS" format.
#'   Text timestamps are UTC; POSIXt objects preserve their actual instant.
#'   Explicit NA values are retained. Direct input vectors may have length 1
#'   or a common compatible length.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param method Extraction method. Options are `"2d"`, `"surface"`, `"bottom"`,
#'   `"nearest"` and `"all"`. With `method = "all"`, the function returns
#'   nearest, surface and bottom outputs together.
#' @param lon_col Name of the longitude column in `data`.
#' @param lat_col Name of the latitude column in `data`.
#' @param date_col Name of the date column in `data`. Required only when a
#'   selected netCDF variable has a time dimension. A missing column is an
#'   error for temporal variables; an explicit missing date retains the row
#'   and returns NA for its temporal output. Static variables do not require
#'   a date column and do not parse unused date values.
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
#' @section Required observation dates:
#' Date requirements follow each selected variable's dimensions, including
#' when a file contains both static and temporal variables. Static variables
#' can be extracted without date or date_col input. Temporal variables require
#' date for direct input, or an existing column named by date_col for tables;
#' date_col = NULL does not bypass this requirement. With multiple sources,
#' every selected temporal variable requires a date input. Explicit NA entries,
#' including NA without a date class, retain their observations; temporal
#' outputs are NA while independent static outputs can still be extracted.
#' An entirely missing date vector does not trigger calendar-label ambiguity
#' checks because no observation has a date to compare.
#' @section Observation time zones and formats:
#' Date inputs have daily precision. POSIXct and POSIXlt inputs preserve the
#' instant represented, including their time-zone offset; matching and calendar
#' labels use UTC. For example, 2020-01-02 00:30 in Europe/Madrid is
#' 2020-01-01 23:30 UTC, so daily matching uses January 1.
#' Character input must use the ISO formats listed for date. A space or T may
#' separate date and time, and a trailing Z is accepted for date-times. Minutes
#' without seconds are interpreted with zero seconds. Text has no inferred
#' local zone. Date-times without Z generate one warning per extraction,
#' even when time_report is FALSE, because they are interpreted as UTC.
#' Date-only text and POSIXt inputs do not generate this warning. Numeric
#' offsets and zone names in text are rejected; supply a
#' POSIXt object constructed with the intended zone instead. Resolve ambiguous
#' or nonexistent local times when constructing that object; this function
#' cannot recover the intended local time from an already converted instant.
#' Conversion errors identify the input, affected observations, original values
#' and identifiers from id_col when available. At most five cases are shown.
#' Original NA values are not conversion failures. For recycled direct inputs,
#' positions refer to the resulting observation rows.
#' @section Physical depth coordinates:
#' Observation depths must be in metres, positive downward (0 at the
#' surface); negative observation depths are rejected, and depths in other
#' units must be converted first. Dataset coordinates are converted from m,
#' meter(s), metre(s), km, kilometer(s), kilometre(s), cm, centimeter(s) or
#' centimetre(s). The positive attribute determines direction when supplied;
#' height and altitude axes without it are read as positive-up. Otherwise
#' uniformly non-negative values imply positive-down and non-positive values
#' imply positive-up. Mixed signs without direction metadata are rejected.
#' The depth coordinate must have units, including for surface and bottom
#' extraction; a file with a single depth level can be read with extract2d().
#' Surface means the physically shallowest layer, even if its value is missing.
#' Bottom means the deepest non-missing layer. Nearest selects the valid layer
#' closest to the observation depth; ties select the shallower layer.
#' Selection does not depend on the storage order of the depth coordinate.
#' Duplicate physical depths, unsupported units, coordinates above the assumed
#' reference surface, and sigma, formula-based or pressure coordinates require
#' explicit correction or conversion before extraction.
#' @section Multiple sources:
#' Spatial matching uses the same rules for one or multiple sources.
#' With spatial_match = "within", only cells containing the point are eligible.
#' With spatial_match = "nearest", the requested point may be outside grid
#' coverage. The nearest cell centre must be within max_distance_km of it.
#' Explicit spatial cell bounds are preferred. Otherwise, or when the bounds
#' variable is missing or invalid (with a warning), edges are inferred from
#' adjacent centres; singleton axes without bounds require coincidence, up to
#' single-precision rounding of the stored coordinate.
#' Longitude comparisons support equivalent coordinates across the antimeridian.
#' A point on an edge shared by two cells uses the first cell in storage order.
#' Among eligible sources the closest temporal candidate with a valid value
#' wins. With "within" matching, equal temporal distances are resolved by
#' input order, so list the preferred source first. With "nearest" matching,
#' they are resolved by geographic distance to the cell centre (differences
#' under 1 m are ignored), then by input order. Existing output columns are
#' replaced, including by NA. Output names cannot overwrite extraction inputs
#' or combine different variables in the same column.
#' @section Custom temporal intervals:
#' With time_resolution = "custom", observations require hours as for datetime.
#' For example, custom_time_interval = as.difftime(3, units = "hours") warns
#' above three hours of displacement. Exact matching requires the same
#' timestamp.
#' Auto matching remains unchanged and uses the observed cadence for
#' diagnostics.
#' @section Spatial and depth reports:
#' The default NULL mode reports only missing outputs, selected cells that do
#' not contain the requested point, and unusually distant depth selections.
#' TRUE requests a full report, including normal and exact selections; FALSE
#' hides the corresponding report. Returned values are identical in all modes.
#' Spatial reports show the requested point, selected cell centre, distance
#' and any explicit distance limit. With within matching, a point may be
#' inside a cell without coinciding with its centre.
#' Depth reports describe physical layers in metres, positive down from the
#' reference surface. Surface means the shallowest layer; bottom means the
#' deepest layer with a valid value, not an independently measured seabed.
#' Values are selected without vertical interpolation. Normal nearest-depth
#' approximations are silent in automatic mode. A difference strictly greater
#' than the local layer spacing is reported. Spacing uses the two physical
#' layers surrounding the requested depth, including layers with missing
#' values. At an exact layer, its nearest neighbour defines the spacing; outside
#' the range, the two nearest edge layers are used. Singleton axes have no
#' spacing to assess a nonzero displacement. Messages identify the reference
#' layers and only attribute a displacement to missing values when closer
#' physical layers must lack valid values under nearest-valid selection.
#' Missing outputs are explained once, with rejection reasons for each
#' relevant source; a repeated reason is shown once with up to three of its
#' files. The cause is taken from the sources whose grid covers the point,
#' unless none does. Temporal causes (including a missing observation date)
#' are reported by time_report; spatial and depth causes by their
#' corresponding switches, so FALSE hides them. Mixed causes are grouped under
#' spatial first, then depth, then time, with all reasons shown.
#' Reports count distinct observations and output values separately. At most
#' five cases are shown, prioritising missing values; truncation and the count
#' omitted are stated. Coordinates use up to five decimals; metres/kilometres
#' normally use up to two. Extra decimals distinguish values close to a limit
#' when rounding would otherwise make them appear equal. Calculations retain
#' full precision. Small positive distances use
#' <0.01 rather than zero. Depth comparisons in reports allow a relative
#' tolerance of 1e-6, so single-precision depth axes match their nominal
#' depths.
#'
#' @section Temporal diagnostics:
#' Diagnostics are emitted after source selection, for values actually returned.
#' Daily, monthly and yearly comparisons warn on a change of calendar period
#' in UTC. Datetime comparisons warn when displacement exceeds one observed
#' cadence step, or when a nonzero displacement has no reliable cadence.
#' With a detected daily cadence, datetime comparison also warns when the
#' selected value falls on a different UTC day, even below 24 hours.
#' Explanatory warnings include the requested and selected timestamps, their
#' difference, the comparison resolution, cadence, row, variable and file name.
#' Reports explain non-warning displacements and missing outputs without
#' repeating examples already explained by a temporal warning.
#' Cadence inference does not change matching or impose a maximum time gap.
#' Singleton sources may use the observed sequence of compatible files for the
#' same variable. Time bounds describe intervals, not sampling frequency.
#' @param spatial_match Horizontal selection rule. "within" (default) requires
#'   the point to lie inside a cell. "nearest" selects the geographically
#'   closest cell centre, including outside the grid, up to max_distance_km.
#'   The rule applies equally to one or multiple sources. A selected missing
#'   value does not trigger a search for another cell in the same source.
#' @param max_distance_km Required finite non-negative distance in kilometres
#'   when spatial_match = "nearest". Distance is the great-circle distance
#'   from the observation to the cell centre on a sphere of radius 6371.01 km
#'   (the radius sf::st_distance() uses with s2), and does not depend on sf
#'   settings. A cell exactly at the limit is accepted. Must be NULL with
#'   "within".
#' @param custom_time_interval Positive difftime value used as the warning
#'   interval when time_resolution = "custom". Matching uses full timestamps.
#'   This interval does not reject candidates; use max_time_diff for that.
#' @param spatial_report NULL (default) reports missing spatial outputs and
#'   selections of cells that do not contain the requested point. TRUE also
#'   reports normal and exact selections; FALSE hides the spatial report.
#' @param depth_report NULL (default) reports missing depth outputs, selections
#'   more than one local layer spacing away, and displaced singleton layers
#'   whose spacing cannot be assessed. TRUE also reports normal surface, bottom
#'   and exact depth selections; FALSE hides the depth report. Depths are in
#'   metres, positive down. Report modes do not change returned values.
#' @param time_report Logical. If TRUE, the default, display a temporal summary
#'   when returned values are displaced in time, or outputs are missing for a
#'   temporal reason (no timestamp within the matching settings, or no
#'   observation date for a source with a time dimension). Outputs missing
#'   because of location or depth are reported by spatial_report and
#'   depth_report instead. The return format is unchanged. FALSE hides the
#'   summary but does not suppress temporal warnings.
#' @param time_resolution Precision used to compare observation times with
#'   dataset timestamps. "auto" uses daily comparison for Date objects
#'   and date-only text, and datetime comparison for POSIXt objects and
#'   date-time text. Explicit alternatives are "datetime", "daily",
#'   "monthly", "yearly" and "custom". Calendar comparisons use UTC.
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
#' # Horizontal matching: this point is outside the example grid's cell edges.
#' # The default within rule returns NA; nearest accepts a centre within 60 km.
#' extract_to_point(
#'   nc = nc_file, lon = -0.6, lat = 40, date = "2020-01-03",
#'   var = "temp", method = "surface", spatial_match = "within"
#' )
#' extract_to_point(
#'   nc = nc_file, lon = -0.6, lat = 40, date = "2020-01-03",
#'   var = "temp", method = "surface", spatial_match = "nearest",
#'   max_distance_km = 60
#' )
#'
#' # Report modes do not change the extracted value.
#' # NULL is automatic: the normal nearest-depth approximation stays silent.
#' automatic <- extract_to_point(
#'   nc = nc_file, lon = 0, lat = 40, depth = 60, date = "2020-01-03",
#'   var = "temp", method = "nearest",
#'   spatial_report = NULL, depth_report = NULL
#' )
#' # TRUE explains the selected cell and depth, including normal selections.
#' detailed <- extract_to_point(
#'   nc = nc_file, lon = 0, lat = 40, depth = 60, date = "2020-01-03",
#'   var = "temp", method = "nearest",
#'   spatial_report = TRUE, depth_report = TRUE
#' )
#' # FALSE hides informational reports; temporal warnings remain active.
#' quiet <- extract_to_point(
#'   nc = nc_file, lon = 0, lat = 40, depth = 60, date = "2020-01-03",
#'   var = "temp", method = "nearest",
#'   spatial_report = FALSE, depth_report = FALSE, time_report = FALSE
#' )
#' identical(automatic, detailed)
#' identical(automatic, quiet)
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
                               "monthly", "yearly", "custom"
                             ),
                             custom_time_interval = NULL,
                             time_report = TRUE,
                             spatial_match = c("within", "nearest"),
                             max_distance_km = NULL,
                             spatial_report = NULL,
                             depth_report = NULL) {
  for (option in c("spatial_report", "depth_report")) {
    value <- get(option)
    if (!is.null(value) &&
          (!is.logical(value) || length(value) != 1L || is.na(value))) {
      stop("`", option, "` must be NULL, TRUE or FALSE.", call. = FALSE)
    }
  }
  spatial_match <- .validate_spatial_matching(spatial_match, max_distance_km)
  method <- match.arg(method)
  time_match <- match.arg(time_match)
  time_resolution <- match.arg(time_resolution)
  custom_seconds <- .validate_custom_time_interval(
    time_resolution, custom_time_interval
  )
  if (!is.logical(time_report) || length(time_report) != 1L ||
        is.na(time_report)) {
    stop(
      "`time_report` must be TRUE or FALSE. ",
      "Use TRUE to show explanatory reports or FALSE to hide them; ",
      "temporal warnings remain active.", call. = FALSE
    )
  }
  comparison_resolution <- if (time_resolution == "custom") {
    "datetime"
  } else {
    time_resolution
  }

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

  # `all` needs depth to calculate nearest-depth output.
  if (method %in% c("nearest", "all")) {
    required_cols <- c(required_cols, depth_col)
  }
  if (!is.null(id_col)) required_cols <- c(required_cols, id_col)
  .check_required_cols(data, required_cols, "data")
  if (method %in% c("nearest", "all")) {
    requested_depth <- data[[depth_col]]
    if (!is.numeric(requested_depth) ||
          any(requested_depth < 0, na.rm = TRUE)) {
      stop(
        "Observation depths must be numeric metres, positive downward ",
        "(0 at the surface). Convert heights (negative values) or km/cm ",
        "depths before extraction.",
        call. = FALSE
      )
    }
  }

  sources <- .as_netcdf_sources(nc, file_col = file_col)

  protected_names <- intersect(
    c(lon_col, lat_col, date_col, depth_col, id_col), names(data)
  )
  if (inherits(data, "sf")) {
    protected_names <- unique(c(
      protected_names, attr(data, "sf_column")
    ))
  }
  source_plan <- .plan_source_outputs(
    sources = sources, var = var, method = method,
    output_prefix = output_prefix, protected_names = protected_names
  )
  horizontal_cache <- new.env(parent = emptyenv())
  horizontal_cache$entries <- list()
  outputs <- lapply(seq_along(sources), function(j) {
    src <- sources[[j]]
    source_var <- source_plan[[j]]$variable
    source_output_prefix <- source_plan[[j]]$prefix

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
      time_resolution = comparison_resolution,
      multi_source = length(sources) > 1L,
      spatial_match = spatial_match,
      max_distance_km = max_distance_km,
      horizontal_cache = horizontal_cache,
      direct_input = direct_input
    )
  })

  merged <- .merge_source_outputs(
    outputs,
    original_names = names(data),
    return_provenance = TRUE,
    spatial_match = spatial_match
  )
  out <- merged$data
  source_provenance <- merged$provenance
  missing <- .show_selection_reports(
    outputs, source_provenance, data, method, depth_col,
    spatial_report, depth_report
  )
  temporal_metadata <- lapply(outputs, function(x) {
    attr(x, ".source_candidate", exact = TRUE)$metadata
  })
  temporal_diagnostics <- data.frame()
  has_temporal_source <- any(vapply(
    temporal_metadata, function(x) !isTRUE(x$static), logical(1)
  ))
  if (has_temporal_source &&
        !is.null(date_col)) {
    input_name <- if (direct_input) {
      "`date`"
    } else {
      paste0("`data$", date_col, "`")
    }
    prepared <- .prepare_time_input(
      data[[date_col]], comparison_resolution,
      input_name = input_name,
      ids = if (is.null(id_col)) NULL else data[[id_col]]
    )
    n_unzoned <- length(prepared$unzoned_rows)
    if (n_unzoned) {
      warning(structure(
        list(message = paste0(
          "Time zone warning: ", n_unzoned,
          if (n_unzoned == 1L) " observation has" else " observations have",
          " times without a time zone.\nInput: ", input_name,
          ". These times are being interpreted as UTC.\n",
          "If they represent local time, supply POSIXct/POSIXlt values ",
          "constructed with the correct time zone.\n",
          "If they are UTC, add Z: \"2020-01-02T19:54Z\"."
        ), call = NULL),
        class = c("ocean3d_timezone_warning", "warning", "condition")
      ))
    }
    diagnostic_resolution <- if (time_resolution == "custom") {
      "custom"
    } else {
      prepared$resolution
    }
    temporal_diagnostics <- .diagnose_temporal_selection(
      source_provenance, prepared$values, diagnostic_resolution,
      temporal_metadata, custom_seconds = custom_seconds
    )
    .warn_temporal_selection(temporal_diagnostics)
    if (isTRUE(time_report)) {
      .show_time_report(
        source_provenance, prepared$values, temporal_metadata,
        temporal_diagnostics, missing
      )
    }
  }

  if (direct_input) {
    value_cols <- unique(unlist(lapply(
      source_plan, function(x) x$columns
    )))

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
#' The surface layer is defined as the physically shallowest depth layer in the
#' netCDF file, whatever the storage order of the depth coordinate.
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
#' For each observation, this function reads the vertical profile at the cell
#' and timestamp selected under spatial_match and time_match. Spatial matching
#' defaults to "within", requiring the point to lie inside a cell; "nearest"
#' requires an explicit max_distance_km limit. The deepest layer with a valid
#' value is returned. This layer does not establish the actual seabed depth.
#' See [extract_to_point()] for matching rules and report options.
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
#' For each observation, this function reads the vertical profile at the cell
#' and timestamp selected under spatial_match and time_match. Spatial matching
#' defaults to "within"; "nearest" requires an explicit max_distance_km limit.
#' It returns three summaries: the valid layer nearest the requested depth,
#' the physically shallowest layer, and the deepest layer with a valid value.
#' See [extract_to_point()] for depth conventions and report options.
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
