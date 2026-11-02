# The time axis of a netCDF variable and matching observations to it. CF time
# units are decoded to UTC date-times, observation times are checked and given
# a precision (date or date-time), and each observation is matched to the
# nearest or an exactly equal time step. Shared by point, buffer and track
# extraction.

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

    examples <- glue::glue("row/position {shown}{identifiers}: {original}")
    if (length(bad) > length(shown)) {
      examples <- c(
        examples,
        glue::glue("Showing {length(shown)} of {length(bad)} cases.")
      )
    }
    examples <- paste(examples, collapse = "\n")

    rlang::abort(glue::glue(
      "{reason}\n",
      "Input: {input_name}\n",
      "Affected observations: {length(bad)}.\n",
      "{examples}\n",
      "Use Date, POSIXct/POSIXlt, YYYY-MM-DD, ",
      "YYYY-MM-DD HH:MM or YYYY-MM-DD HH:MM:SS.\n",
      "Text timestamps are interpreted in UTC. ",
      "For another time zone, supply a POSIXct/POSIXlt object."
    ), call = NULL)
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
