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
    return(paste0(
      "No value was extracted.\n", location,
      "\nObservation: ", requested,
      "\nNo source supplied a valid value under the requested settings.",
      " The result is NA; no selected timestamp or file can be assigned.",
      "\nCheck the observation, source coverage and matching limits."
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

.show_time_report <- function(provenance, requested, metadata,
                              diagnostics, sources) {
  rows <- list()
  for (column in names(provenance)) {
    selected <- provenance[[column]]
    for (i in seq_len(nrow(selected))) {
      source_index <- selected$source_index[i]
      timestamp <- selected$timestamp[i]
      hit <- if (nrow(diagnostics)) {
        which(diagnostics$row == i & diagnostics$column == column)
      } else {
        integer()
      }
      if (length(hit)) {
        rows[[length(rows) + 1L]] <- diagnostics[hit[1L], , drop = FALSE]
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
        requested = requested[i], selected = timestamp,
        difference_seconds = NA_real_, resolution = NA_character_,
        cadence = "unknown", cadence_basis = NA_character_,
        cadence_seconds = NA_real_, cadence_months = NA_integer_,
        custom_seconds = NA_real_, cadence_exceeded = FALSE,
        source_file = NA_character_,
        status = if (is.na(source_index)) "no_value" else "not_assessable",
        reason = NA_character_, stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(invisible(NULL))
  report <- do.call(rbind, rows)
  interesting <- report[report$status != "matched", , drop = FALSE]
  if (!nrow(interesting)) return(invisible(report))
  # Warning examples already explain flagged selections.
  # Do not print them twice.
  information <- interesting[is.na(interesting$reason), , drop = FALSE]
  if (!nrow(information)) return(invisible(report))
  examples <- .temporal_example_rows(information)
  detail <- vapply(seq_len(nrow(examples)), function(i) {
    .format_temporal_example(examples[i, , drop = FALSE])
  }, character(1))
  message(structure(list(
    message = paste0(
      "Temporal report: ", length(unique(information$row)),
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

.unwrap_longitude <- function(values) {
  if (length(values) < 2L) return(values)
  original_steps <- diff(values)
  steps <- (original_steps + 180) %% 360 - 180
  ambiguous <- abs(abs(steps) - 180) < 1e-8
  steps[ambiguous] <- original_steps[ambiguous]
  values[1L] + c(0, cumsum(steps))
}

.read_spatial_bounds <- function(nc, coordinate, values) {
  name <- .read_time_attribute(nc, coordinate, "bounds")
  if (is.null(name)) return(NULL)
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("Invalid spatial bounds reference.", call. = FALSE)
  }
  index <- which(vapply(
    nc$var, function(x) identical(x$name, name), logical(1)
  ))
  if (length(index) != 1L) {
    stop("Spatial bounds variable was not found: ", name, call. = FALSE)
  }
  variable <- nc$var[[index]]
  axes <- vapply(variable$dim, function(x) x$name, character(1))
  coordinate_axis <- which(axes == coordinate)
  raw <- ncdf4::ncvar_get(nc, name, collapse_degen = FALSE)
  if (length(dim(raw)) != 2L || length(coordinate_axis) != 1L) {
    stop("Unsupported spatial bounds dimensions.", call. = FALSE)
  }
  if (coordinate_axis == 1L) raw <- t(raw)
  if (!identical(dim(raw), c(2L, length(values))) ||
        any(!is.finite(raw))) {
    stop("Invalid spatial bounds values.", call. = FALSE)
  }
  coordinate_units <- .read_time_attribute(nc, coordinate, "units")
  bounds_units <- .read_time_attribute(nc, name, "units")
  if (!is.null(bounds_units) &&
        !identical(bounds_units, coordinate_units)) {
    stop(
      "Spatial bounds and coordinates use different units.",
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
    stop("Spatial coordinates must form a monotonic axis.", call. = FALSE)
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
      stop(
        "Spatial bounds do not contain their cell centres.",
        call. = FALSE
      )
    }
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

.spatial_axis_index <- function(axis, target) {
  if (length(target) != 1L || !is.finite(target)) return(NA_integer_)
  if (axis$longitude) {
    targets <- target + 360 * round((axis$centres - target) / 360)
  } else {
    if (target < -90 || target > 90) return(NA_integer_)
    targets <- rep(target, length(axis$centres))
  }
  contained <- which(
    targets >= axis$lower - 1e-8 & targets <= axis$upper + 1e-8
  )
  if (!length(contained)) return(NA_integer_)
  distances <- abs(axis$centres[contained] - targets[contained])
  as.integer(contained[which.min(distances)])
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


# A source may offer a temporal candidate without implying a coverage interval.
# Matching and tolerance are handled by .match_time_index().
.temporal_source_eligible <- function(time_values, target, resolution) {
  if (is.null(time_values)) return(TRUE)
  if (length(time_values) == 0L || is.na(target)) return(FALSE)

  available <- .time_comparison_values(time_values, resolution)
  requested <- .time_comparison_values(target, resolution)
  !is.na(requested) && any(!is.na(available))
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
                              multi_source = FALSE) {
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
  spatial_axes <- NULL
  if (isTRUE(multi_source)) {
    spatial_axes <- list(
      lon = .prepare_spatial_axis(
        nc, dims$lon, lon_values, longitude = TRUE
      ),
      lat = .prepare_spatial_axis(
        nc, dims$lat, lat_values, longitude = FALSE
      )
    )
  }
  depth_values <- .get_dim_values(nc, dims$depth)
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

  candidate_distance <- rep(Inf, nrow(data))
  candidate_time <- rep(NA_real_, nrow(data))
  candidate_index <- rep(NA_integer_, nrow(data))

  for (i in seq_len(nrow(data))) {
    if (isTRUE(verbose) && i %% 100 == 0) {
      message("Extracting row ", i, " of ", nrow(data))
    }

    if (isTRUE(multi_source)) {
      lon_idx <- .spatial_axis_index(
        spatial_axes$lon, data[[lon_col]][i]
      )
      lat_idx <- .spatial_axis_index(
        spatial_axes$lat, data[[lat_col]][i]
      )
    } else {
      lon_target <- .prepare_lon(lon_values, data[[lon_col]][i])
      lat_target <- data[[lat_col]][i]
      lon_idx <- .nearest_index(lon_values, lon_target)
      lat_idx <- .nearest_index(lat_values, lat_target)
    }

    if (is.na(lon_idx) || is.na(lat_idx)) next


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
      # Static sources retain first-valid-source priority.
      candidate_distance[i] <- 0
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

  attr(out, ".source_candidate") <- list(
    distance_days = candidate_distance,
    timestamp = candidate_time,
    time_index = candidate_index,
    variable = var,
    metadata = temporal_metadata,
    columns = .output_column_names(method, var, output_prefix)
  )
  out
}

.merge_source_outputs <- function(
  outputs,
  original_names,
  return_provenance = FALSE
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
      # Strict inequality preserves input priority only for temporal ties.
      replace <- !is.na(x[[col]]) & is.finite(distance) &
        (is.na(values) | distance < best_distance)
      values[replace] <- x[[col]][replace]
      best_distance[replace] <- distance[replace]
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
#' @section Multiple sources:
#' With multiple sources, only cells containing the observation are eligible.
#' Explicit spatial cell bounds are preferred. Otherwise, edges are inferred
#' from adjacent centres; singleton axes without bounds require coincidence.
#' Longitude comparisons support equivalent coordinates across the antimeridian.
#' Among eligible sources the closest temporal candidate with a valid value
#' wins, with input order breaking ties. Existing output columns are replaced,
#' including by NA. Output names cannot overwrite extraction inputs or combine
#' different variables in the same column.
#' @section Custom temporal intervals:
#' With time_resolution = "custom", observations require hours as for datetime.
#' For example, custom_time_interval = as.difftime(3, units = "hours") warns
#' above three hours of displacement. Exact matching requires the same
#' timestamp.
#' Auto matching remains unchanged and uses the observed cadence for
#' diagnostics.
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
#' @param custom_time_interval Positive difftime value used as the warning
#'   interval when time_resolution = "custom". Matching uses full timestamps.
#'   This interval does not reject candidates; use max_time_diff for that.
#' @param time_report Logical. If TRUE, the default, display a temporal summary
#'   when returned values are displaced in time or outputs are missing.
#'   The return format is unchanged. FALSE hides the summary but does not
#'   suppress temporal warnings.
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
                             time_report = TRUE) {
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
      # Validate without replacing the original precision.
      .prepare_time_input(direct_values$date, comparison_resolution)
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
      multi_source = length(sources) > 1L
    )
  })

  merged <- .merge_source_outputs(
    outputs,
    original_names = names(data),
    return_provenance = TRUE
  )
  out <- merged$data
  source_provenance <- merged$provenance
  temporal_metadata <- lapply(outputs, function(x) {
    attr(x, ".source_candidate", exact = TRUE)$metadata
  })
  temporal_diagnostics <- data.frame()
  has_temporal_source <- any(vapply(
    temporal_metadata, function(x) !isTRUE(x$static), logical(1)
  ))
  if (has_temporal_source &&
        !is.null(date_col)) {
    prepared <- .prepare_time_input(data[[date_col]], comparison_resolution)
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
        temporal_diagnostics, sources
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
