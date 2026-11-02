# The time axis of a netCDF variable and matching observations to it. CF time
# units are decoded to UTC date-times, observation times are checked and given
# a precision (date or date-time), and each observation is matched to the
# nearest or an exactly equal time step. Shared by point, buffer and track
# extraction.

# Calendars the package decodes. Files on other CF calendars (noleap,
# 360_day, all_leap, julian) stop here until output columns and voxels can
# show dates that R's date-times can't hold (SPEC.md, "Time decoding"). The
# core object (#113) is to use the same check, so support can be added in one
# place.
.check_time_calendar <- function(calendar, time_dim) {
  if (is.null(calendar)) return("standard")
  supported <- c("standard", "gregorian", "proleptic_gregorian")
  if (is.character(calendar) && length(calendar) == 1L &&
        !is.na(calendar) && tolower(calendar) %in% supported) {
    return(tolower(calendar))
  }
  rlang::abort(glue::glue(
    "Time coordinate `{time_dim}` uses an unsupported calendar: ",
    "{toString(calendar)}.\n",
    "ocean3d currently supports the standard, gregorian and ",
    "proleptic_gregorian calendars. Convert the file to one of these ",
    "first, for example with CDO."
  ), class = "ocean3d_error_calendar", call = NULL)
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

  calendar <- .check_time_calendar(
    .read_time_attribute(nc, time_dim, "calendar"), time_dim
  )

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

  # CFtime reads the units and reference date-time (UTC or Z suffixes,
  # unpadded dates, missing seconds, years before 1000) and converts offsets
  # to date-times without the rounding error of origin + offset * seconds.
  axis <- withCallingHandlers(
    tryCatch(
      CFtime::CFtime(units, calendar, values),
      error = function(e) {
        rlang::abort(glue::glue(
          "Unsupported time units for `{time_dim}`: '{units}'.\n",
          "Expected a unit such as days or hours followed by 'since' and ",
          "a reference date-time ({conditionMessage(e)})."
        ), class = "ocean3d_error_time_units", call = NULL)
      }
    ),
    # Repeated or unsorted time steps are dealt with when observations are
    # matched to the axis.
    warning = function(w) {
      if (grepl("monotonic", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  origin <- axis$cal$origin
  if (!identical(origin$tz, "+0000")) {
    rlang::abort(glue::glue(
      "Unsupported reference date-time in units '{units}': time zone ",
      "offset {origin$tz}.\n",
      "Numeric time zone offsets are not supported yet; the reference ",
      "date-time must be in UTC (optionally followed by Z or UTC)."
    ), class = "ocean3d_error_time_units", call = NULL)
  }

  # Months and years have no fixed length, so a series of them is ambiguous
  # (CF defines a month as 1/12 of a year; producers often mean calendar
  # months). A single step needs no conversion.
  if (axis$unit %in% c("months", "years") && length(unique(values)) > 1L) {
    rlang::abort(glue::glue(
      "Unsupported time units for `{time_dim}`: '{units}'.\n",
      "A time axis of several months or years is ambiguous; months, years ",
      "and climatologies require separate handling."
    ), class = "ocean3d_error_time_units", call = NULL)
  }

  axis$as_timestamp(asPOSIX = TRUE)
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
