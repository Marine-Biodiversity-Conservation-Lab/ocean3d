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

# Decodes a CF time axis to UTC date-times. With `bounds` (a 2 x n matrix of
# offsets in the axis units, one column per time step), returns those bounds
# instead, as a 2 x n matrix of UTC seconds since 1970. CFtime holds the
# bounds with the axis, so shared edges of back-to-back intervals
# (0, 31, 31, 60) aren't read as a repeated time axis, which CFtime warns
# about (#113). Bounds must contain their time step; either edge may come
# first, and each column keeps its order.
.decode_netcdf_time <- function(nc, time_dim, bounds = NULL) {
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

  values <- as.vector(.get_dim_values(nc, time_dim))
  offsets <- c(values, bounds)

  if (!is.numeric(offsets) || any(!is.finite(offsets))) {
    stop(
      "Time coordinate `", time_dim,
      "` must contain finite numeric offsets.",
      call. = FALSE
    )
  }

  # CFtime reads the units and reference date-time (UTC or Z suffixes,
  # unpadded dates, missing seconds, years before 1000) and converts offsets
  # to date-times without the rounding error of origin + offset * seconds.
  axis <- tryCatch(
    CFtime::CFtime(units, calendar, values),
    error = function(e) {
      rlang::abort(glue::glue(
        "Unsupported time units for `{time_dim}`: '{units}'.\n",
        "Expected a unit such as days or hours followed by 'since' and ",
        "a reference date-time ({conditionMessage(e)})."
      ), class = "ocean3d_error_time_units", call = NULL)
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
  if (axis$unit %in% c("months", "years") && length(unique(offsets)) > 1L) {
    rlang::abort(glue::glue(
      "Unsupported time units for `{time_dim}`: '{units}'.\n",
      "A time axis of several months or years is ambiguous; months, years ",
      "and climatologies require separate handling."
    ), class = "ocean3d_error_time_units", call = NULL)
  }

  if (is.null(bounds)) return(axis$as_timestamp(asPOSIX = TRUE))

  # CFtime only warns about bounds it can't use, and then drops them.
  reversed <- bounds[1L, ] > bounds[2L, ]
  withCallingHandlers(
    axis$set_bounds(rbind(
      pmin(bounds[1L, ], bounds[2L, ]), pmax(bounds[1L, ], bounds[2L, ])
    )),
    warning = function(w) {
      rlang::abort(glue::glue(
        "Invalid bounds for `{time_dim}`: each interval must contain its ",
        "time step ({conditionMessage(w)})."
      ), class = "ocean3d_error_time_bounds", call = NULL)
    }
  )

  # Bounds come back only as text: dates, or date-times with seconds that
  # may have a fraction.
  stamps <- axis$get_bounds(format = "timestamp")
  format <- if (any(grepl("T", stamps, fixed = TRUE))) {
    "%Y-%m-%dT%H:%M:%OS"
  } else {
    "%Y-%m-%d"
  }
  seconds <- matrix(
    as.numeric(as.POSIXct(strptime(stamps, format, tz = "UTC"))), nrow = 2L
  )
  seconds[, reversed] <- seconds[2:1, reversed]
  seconds
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

  if (resolution == "daily") {
    # The start of each UTC day. Arithmetic gives the same keys as formatting
    # to text and back, about 600 times faster (#89).
    return(.POSIXct(floor(as.numeric(values) / 86400) * 86400, tz = "UTC"))
  }

  pattern <- switch(
    resolution,
    monthly = "%Y-%m-01",
    yearly = "%Y-01-01"
  )
  as.POSIXct(
    as.Date(format(values, pattern, tz = "UTC")),
    tz = "UTC"
  )
}

# Match observation times to a time axis. `targets` are already prepared
# (.prepare_time_input()), and `resolution` is their resolved comparison
# resolution. The keys for the axis and for the observations are each
# computed once (#89). A time step repeated at that resolution is only an
# error when an observation selects it (#91). Returns one index per target,
# NA where nothing matches.
.match_time_indices <- function(values, targets,
                                time_match = c("nearest", "exact"),
                                max_time_diff = NULL, resolution,
                                file = NULL) {
  time_match <- match.arg(time_match)
  out <- rep(NA_integer_, length(targets))
  if (!length(values) || all(is.na(values))) return(out)

  available <- as.numeric(.time_comparison_values(values, resolution))
  requested <- as.numeric(.time_comparison_values(targets, resolution))
  repeated <- available %in% available[duplicated(available) &
                                         !is.na(available)]

  for (i in which(!is.na(requested))) {
    if (time_match == "exact") {
      idx <- match(requested[i], available)
    } else {
      gap_days <- abs(available - requested[i]) / 86400
      # An observation exactly halfway between two steps takes the one on
      # its own UTC day, so a midnight observation between noon-stamped days
      # and a noon observation between midnight-stamped days both keep their
      # day (#71). Otherwise the earlier step wins, as before.
      closest <- which(gap_days == min(gap_days, na.rm = TRUE))
      same_day <- closest[floor(as.numeric(values[closest]) / 86400) ==
                            floor(as.numeric(targets[i]) / 86400)]
      idx <- if (length(same_day) == 1L) same_day else closest[1L]
      if (!is.null(max_time_diff) && gap_days[idx] > max_time_diff) {
        idx <- NA_integer_
      }
    }
    if (!is.na(idx) && repeated[idx]) {
      key <- format(.POSIXct(available[idx], tz = "UTC"),
                    "%Y-%m-%d %H:%M:%S UTC")
      count <- sum(available == available[idx], na.rm = TRUE)
      where <- if (is.null(file)) "" else glue::glue(" in {file}")
      rlang::abort(glue::glue(
        "Multiple dataset timestamps map to the same '{resolution}' ",
        "comparison value: {key} appears {count} times{where}.\n",
        "Selecting the first would discard information. Use a finer ",
        "time_resolution, select the intended records, or aggregate the ",
        "data explicitly before extraction."
      ), class = "ocean3d_error_repeated_time", call = NULL)
    }
    out[i] <- as.integer(idx)
  }
  out
}

# One observation at a time, for callers that match row by row.
.match_time_index <- function(
  values,
  target,
  time_match = c("nearest", "exact"),
  max_time_diff = NULL,
  time_resolution = c(
    "auto", "datetime", "daily", "monthly", "yearly"
  ),
  file = NULL) {
  time_match <- match.arg(time_match)
  time_resolution <- match.arg(time_resolution)

  if (length(values) == 0L ||
        all(is.na(values)) ||
        is.na(target)) {
    return(NA_integer_)
  }

  prepared <- .prepare_time_input(target, time_resolution)
  .match_time_indices(
    values, prepared$values, time_match, max_time_diff, prepared$resolution,
    file = file
  )
}

# A source may offer a temporal candidate without implying a coverage interval.
# Matching and tolerance are handled by .match_time_indices().
.temporal_source_eligible <- function(time_values, target, resolution) {
  if (is.null(time_values)) return(TRUE)
  !is.na(target) && any(!is.na(time_values))
}
