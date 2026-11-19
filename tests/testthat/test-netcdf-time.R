# CF time decoding and matching (#88-#92, #71, #111, #113), tested on the
# shared helpers directly: the core object (#113) reuses them, while the
# extraction functions that call them may be retired (#116; AGENTS.md,
# "Tests"). Each file is built with make_grid_file() (helper-netcdf-grid.R)
# without a depth axis; time bounds are 2 x n matrices in the axis units.

decode_time <- function(f, bounds = NULL, readunlim = TRUE) {
  nc <- ncdf4::nc_open(f, readunlim = readunlim)
  on.exit(ncdf4::nc_close(nc))
  .decode_netcdf_time(nc, "time", bounds = bounds)
}

# Index of the time step each observation time is matched to.
match_times <- function(axis, times, ...) {
  prepared <- .prepare_time_input(times)
  .match_time_indices(axis, prepared$values,
                      resolution = prepared$resolution, ...)
}

utc <- function(x) as.POSIXct(x, tz = "UTC")

test_that("a single-step months axis is read and a longer one rejected (#88)", {
  one <- make_grid_file(depth = NULL, time = 0,
                        time_units = "months since 1955-01-01 00:00:00")
  several <- make_grid_file(depth = NULL, time = c(0, 1, 2),
                            time_units = "months since 2020-01-01")
  on.exit(unlink(c(one, several)), add = TRUE)

  expect_equal(decode_time(one), utc("1955-01-01"))
  err <- expect_error(decode_time(several), class = "ocean3d_error_time_units")
  expect_match(conditionMessage(err), "months since 2020-01-01", fixed = TRUE)
  expect_match(conditionMessage(err),
               "several months or years is ambiguous", fixed = TRUE)
})

test_that("reference times with a numeric time zone offset are rejected", {
  f <- make_grid_file(depth = NULL, time = c(0, 24),
                      time_units = "hours since 2020-01-01 05:00:00+05:00")
  on.exit(unlink(f), add = TRUE)

  err <- expect_error(decode_time(f), class = "ocean3d_error_time_units")
  expect_match(conditionMessage(err), "time zone offset +0500", fixed = TRUE)
})

test_that("fractional-day timestamps decode and match exactly (#90)", {
  base <- as.numeric(as.Date("2020-06-15"))
  f <- make_grid_file(depth = NULL, time = base + (0:143) / 144,
                      time_units = "days since 1970-01-01")
  on.exit(unlink(f), add = TRUE)
  times <- utc("2020-06-15") + (0:143) * 600

  axis <- decode_time(f)
  expect_equal(axis, times)
  expect_identical(match_times(axis, times, time_match = "exact"), 1:144)
})

test_that("calendars other than the standard ones stop with a clear error", {
  for (calendar in c("noleap", "360_day", "all_leap", "julian")) {
    f <- make_grid_file(depth = NULL, time = c(0, 1),
                        time_units = "days since 2000-01-01",
                        calendar = calendar)
    on.exit(unlink(f), add = TRUE)
    err <- expect_error(decode_time(f), class = "ocean3d_error_calendar",
                        label = calendar)
    expect_match(conditionMessage(err),
                 paste("unsupported calendar:", calendar), fixed = TRUE)
  }
})

test_that("the standard calendars are accepted under each name", {
  expected <- utc(c("2020-01-01", "2020-01-03", "2020-01-06"))
  for (calendar in c("standard", "gregorian", "proleptic_gregorian",
                     "Gregorian")) {
    f <- make_grid_file(depth = NULL, time = c(0, 2, 5), calendar = calendar)
    on.exit(unlink(f), add = TRUE)
    expect_equal(decode_time(f), expected, label = calendar)
  }
})

test_that("a readunlim = FALSE connection reads the time axis (#92)", {
  f <- make_grid_file(depth = NULL, time = c(0, 2, 5), time_unlim = TRUE)
  on.exit(unlink(f), add = TRUE)

  expect_equal(decode_time(f, readunlim = FALSE),
               utc(c("2020-01-01", "2020-01-03", "2020-01-06")))
})

test_that("many dates match a long daily axis in one pass (#89)", {
  f <- make_grid_file(depth = NULL, time = 0:3649,
                      time_units = "days since 2010-01-01")
  on.exit(unlink(f), add = TRUE)
  set.seed(1)
  dates <- as.Date("2010-01-01") + sample(0:3649, 500, replace = TRUE)

  expect_identical(match_times(decode_time(f), dates),
                   as.integer(dates - as.Date("2010-01-01")) + 1L)
})

test_that("a repeated time step is only an error if it is selected (#91)", {
  f <- make_grid_file(depth = NULL, time = c(0, 1, 1, 2))
  on.exit(unlink(f), add = TRUE)

  # CFtime warns about the repeated offset; that warning is not silenced.
  expect_warning(axis <- decode_time(f),
                 "Offsets not monotonically increasing", fixed = TRUE)
  expect_identical(match_times(axis, as.Date("2020-01-01")), 1L)
  expect_identical(match_times(axis, as.Date("2020-01-03")), 4L)

  err <- expect_error(
    match_times(axis, as.Date("2020-01-02"), file = basename(f)),
    class = "ocean3d_error_repeated_time"
  )
  expect_match(conditionMessage(err),
               "Multiple dataset timestamps", fixed = TRUE)
  expect_match(conditionMessage(err),
               "2020-01-02 00:00:00 UTC appears 2 times", fixed = TRUE)
  expect_match(conditionMessage(err), basename(f), fixed = TRUE)
})

test_that("midnight date-times keep their day on noon-stamped files (#71)", {
  f <- make_grid_file(depth = NULL, time = c(613620, 613644, 613668),
                      time_units = "hours since 1950-01-01")
  on.exit(unlink(f), add = TRUE)
  times <- utc(c("2020-01-01", "2020-01-02", "2020-01-03"))

  axis <- decode_time(f)
  expect_equal(axis, times + 12 * 3600)
  expect_identical(match_times(axis, times), 1:3)
  # Dates compare by day, so a noon time step is no time away.
  for (time_match in c("nearest", "exact")) {
    expect_identical(match_times(axis, as.Date(times), time_match = time_match,
                                 max_time_diff = 0), 1:3)
  }
})

# Monthly means stamped mid-month, with back-to-back bounds: January and
# February 2020 share the 1 February edge (offsets 0, 31, 31, 60).
monthly_bounds <- matrix(c(0, 31, 31, 60), nrow = 2L)
month_edges <- matrix(as.numeric(utc(
  c("2020-01-01", "2020-02-01", "2020-02-01", "2020-03-01")
)), nrow = 2L)

read_metadata <- function(f) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc))
  .read_temporal_metadata(nc, "time", "temp", .decode_netcdf_time(nc, "time"))
}

test_that("back-to-back time bounds decode without a warning (#113)", {
  f <- make_grid_file(depth = NULL, time = c(15, 45),
                      time_bounds = monthly_bounds)
  on.exit(unlink(f), add = TRUE)

  expect_no_warning(edges <- decode_time(f, bounds = monthly_bounds))
  expect_equal(edges, month_edges)

  # The bounds reach the time metadata as the month edges, even when warnings
  # are errors.
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  expect_no_error(read_metadata(f))
})

test_that("time bounds may list the later edge first (#113)", {
  f <- make_grid_file(depth = NULL, time = c(15, 45),
                      time_bounds = monthly_bounds[2:1, ])
  on.exit(unlink(f), add = TRUE)

  expect_no_warning(edges <- decode_time(f, bounds = monthly_bounds[2:1, ]))
  expect_equal(edges, month_edges[2:1, ])
})

test_that("time bounds that miss their time step are rejected (#113)", {
  # January's bounds end on 10 January, before its 16 January time step.
  bad <- matrix(c(0, 10, 31, 60), nrow = 2L)
  f <- make_grid_file(depth = NULL, time = c(15, 45), time_bounds = bad)
  on.exit(unlink(f), add = TRUE)

  err <- expect_error(decode_time(f, bounds = bad),
                      class = "ocean3d_error_time_bounds")
  expect_match(conditionMessage(err), "Invalid bounds for `time`", fixed = TRUE)
  expect_match(conditionMessage(err),
               "must contain its time step", fixed = TRUE)

  # Point extraction doesn't read time bounds, so they don't stop it.
  expect_no_error(read_metadata(f))
})

test_that("the time axis builder returns a CFtime object with its bounds", {
  f <- make_grid_file(depth = NULL, time = c(15, 45),
                      time_bounds = monthly_bounds[2:1, ])
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  plain <- .netcdf_time_axis(nc, "time")
  expect_s3_class(plain, "CFTime")
  expect_identical(plain$offsets, c(15, 45))
  expect_null(plain$get_bounds())

  # Bounds that list the later edge first are stored lower edge first.
  bounded <- .netcdf_time_axis(nc, "time", bounds = monthly_bounds[2:1, ])
  expect_equal(bounded$get_bounds(), monthly_bounds)
})


test_that("climatological time axes are recognised by attribute or year 0", {
  check <- function(f) {
    nc <- ncdf4::nc_open(f)
    on.exit(ncdf4::nc_close(nc))
    .check_climatology(nc, "time", basename(f))
  }
  dated <- make_grid_file(depth = NULL, time = c(0, 1))
  woa <- make_grid_file(depth = NULL, time = c(0.5, 1.5),
                        time_units = "months since 0000-01-01 00:00:00")
  flagged <- make_grid_file(depth = NULL, time = c(0, 1), time_atts = list(
    climatology = "climatology_bounds"
  ))
  on.exit(unlink(c(dated, woa, flagged)), add = TRUE)

  expect_null(check(dated))
  # Days since year 0 (MATLAB datenums) are dates, not a climatology.
  datenum <- make_grid_file(depth = NULL, time = c(737791, 737792),
                            time_units = "days since 0000-01-01 00:00:00",
                            calendar = "proleptic_gregorian")
  on.exit(unlink(datenum), add = TRUE)
  expect_null(check(datenum))
  err <- expect_error(check(woa), class = "ocean3d_error_time_climatology")
  expect_match(conditionMessage(err), basename(woa), fixed = TRUE)
  expect_match(conditionMessage(err), "months since 0000-01-01", fixed = TRUE)
  expect_error(check(flagged),
               class = "ocean3d_error_time_climatology")
})

test_that("a single monthly step with bounds is read", {
  # A one-month mean: the two bound edges aren't extra time steps.
  f <- make_grid_file(depth = NULL, time = 0,
                      time_units = "months since 2020-01-01",
                      time_bounds = rbind(0, 1))
  on.exit(unlink(f), add = TRUE)
  expect_equal(
    extract_to_point(nc = f, lon = 1.5, lat = 10.5, var = "temp",
                     date = as.Date("2020-01-15"), method = "2d"),
    2
  )
  expect_identical(as_proxy(f)@axes$time$edges, "bounds")
})

# Observation times (.prepare_time_input()) ------------------------------------

test_that("time validation identifies malformed values not inventing hours", {
  expect_error(
    .prepare_time_input(c("2020-01-02", "06/01/2020")),
    "row/position"
  )
  expect_error(
    .prepare_time_input("2020-02-30"),
    "Invalid date/time"
  )
  expect_error(
    .prepare_time_input(
      c("2020-01-02", "2020-01-02 19:54:00"),
      "datetime"
    ),
    "requires hours for every"
  )
  expect_error(
    .prepare_time_input(as.Date("2020-01-02"), "datetime"),
    "requires observation hours"
  )
  parsed <- .prepare_time_input(
    c("2020-01-02", NA_character_),
    "daily"
  )
  expect_true(is.na(parsed$values[2]))
})

test_that("minute and second ISO text share datetime precision", {
  x <- c("2020-01-02 19:54", "2020-01-02T19:54Z",
         "2020-01-02 19:54:30", NA_character_)
  parsed <- .prepare_time_input(x)
  expect_identical(parsed$resolution, "datetime")
  expect_equal(format(parsed$values, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
               c("2020-01-02 19:54:00", "2020-01-02 19:54:00",
                 "2020-01-02 19:54:30", NA_character_))
  expect_identical(.prepare_time_input("2020-01-02")$resolution, "daily")
  expect_error(.prepare_time_input(c("2020-01-02", "2020-01-02 19:54")),
               "mixes dates", fixed = TRUE)
})

test_that("explicit missing vectors are not conversion failures", {
  for (x in list(NA, c(NA, NA), NA_real_, NA_character_,
                  as.Date(NA_character_))) {
    for (resolution in c("auto", "datetime", "daily", "monthly", "yearly")) {
      parsed <- .prepare_time_input(x, resolution)
      expect_length(parsed$values, length(x))
      expect_true(all(is.na(parsed$values)))
    }
  }
  mixed <- .prepare_time_input(c("2020-01-02 19:54", NA_character_))
  expect_false(is.na(mixed$values[1L]))
  expect_true(is.na(mixed$values[2L]))
  expect_error(.prepare_time_input(c(1, NA_real_)), "Unsupported input type",
               fixed = TRUE)
})

test_that("time errors identify values and IDs while bounding examples", {
  x <- c("2020-01-02", NA_character_, "06/01/2020", "2020-02-30",
         rep("2020-01-02junk", 5))
  error <- tryCatch(.prepare_time_input(
    x, input_name = "`data$observed_at`", ids = paste0("station-", 1:9)
  ), error = identity)
  expect_s3_class(error, "error")
  text <- conditionMessage(error)
  expect_match(text, "Input: `data$observed_at`", fixed = TRUE)
  expect_match(text, "Affected observations: 7", fixed = TRUE)
  expect_match(text, "row/position 3 (id: \"station-3\"): \"06/01/2020\"",
               fixed = TRUE)
  expect_match(text, "2020-02-30", fixed = TRUE)
  expect_match(text, "Showing 5 of 7 cases", fixed = TRUE)
  expect_false(grepl("row/position 2", text, fixed = TRUE))
  expect_false(grepl("station-8", text, fixed = TRUE))
  expect_match(text, "YYYY-MM-DD HH:MM", fixed = TRUE)
  single <- tryCatch(.prepare_time_input("2020-01-02junk"), error = identity)
  expect_false(grepl("Showing", conditionMessage(single), fixed = TRUE))
})

test_that("POSIX time zones preserve instants across UTC day boundaries", {
  local <- as.POSIXct(c("2020-01-02 00:30:00", "2020-07-02 00:30:00"),
                      tz = "Europe/Madrid")
  for (x in list(local, as.POSIXlt(local))) {
    parsed <- .prepare_time_input(x)
    expect_identical(parsed$resolution, "datetime")
    expect_equal(format(parsed$values, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                 c("2020-01-01 23:30:00", "2020-07-01 22:30:00"))
    expect_equal(as.numeric(parsed$values), as.numeric(local))
  }
  transition <- as.POSIXct(c("2020-03-29 01:30:00", "2020-03-29 03:30:00"),
                           tz = "Europe/Madrid")
  parsed <- .prepare_time_input(as.POSIXlt(transition))
  expect_equal(format(parsed$values, "%H:%M:%S", tz = "UTC"),
               c("00:30:00", "01:30:00"))
  expect_equal(as.numeric(diff(parsed$values), units = "hours"), 1)
  text <- .prepare_time_input("2020-01-02 00:30")
  expect_equal(format(text$values, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
               "2020-01-02 00:30:00")
  for (bad in c("2020-01-02 00:30+01:00", "2020-01-02 00:30 Europe/Madrid")) {
    expect_error(.prepare_time_input(bad), "supply a POSIXct/POSIXlt object",
                 fixed = TRUE)
  }
})

test_that("only text clock times without a zone are unzoned", {
  x <- c("2020-01-01 00:00", "2020-01-01T00:00Z", "2020-01-01", NA_character_)
  expect_identical(.prepare_time_input(x, "daily")$unzoned_rows, 1L)
  for (x in list("2020-01-01", as.Date("2020-01-01"),
                  "2020-01-01T12:00Z", "2020-01-01 12:00:00Z",
                  as.POSIXct("2020-01-01 12:00:00", tz = "Europe/Madrid"),
                  as.POSIXlt("2020-01-01 12:00:00", tz = "Europe/Madrid"),
                  NA)) {
    expect_length(.prepare_time_input(x)$unzoned_rows, 0L)
  }
})

# Matching one observation time (.match_time_index()) --------------------------

test_that("temporal matching handles empty and entirely missing time axes", {
  for (values in list(as.Date(character()), as.Date(c(NA, NA)))) {
    for (method in c("nearest", "exact")) {
      expect_identical(
        .match_time_index(
          values, as.Date("2020-01-03"), time_match = method
        ),
        NA_integer_
      )
    }
  }
})

test_that("daily comparison keeps late observations on the same day", {
  observations <- as.POSIXct(
    c(
      "2020-01-02 00:01:00",
      "2020-01-02 12:00:00",
      "2020-01-02 19:54:00",
      "2020-01-02 23:59:00",
      "2020-01-03 00:00:00"
    ),
    tz = "UTC"
  )

  for (hour in c("00:00:00", "12:00:00")) {
    available <- as.POSIXct(
      paste(
        c("2020-01-01", "2020-01-02", "2020-01-03"),
        hour
      ),
      tz = "UTC"
    )

    for (matching in c("nearest", "exact")) {
      indices <- vapply(seq_along(observations), function(i) {
        .match_time_index(
          available, observations[i],
          time_match = matching,
          time_resolution = "daily"
        )
      }, integer(1))

      expect_equal(indices, c(2L, 2L, 2L, 2L, 3L))
    }
  }
})

test_that("datetime matching preserves minutes and time limits", {
  available <- as.POSIXct(
    c("2020-01-02 19:00:00", "2020-01-02 20:00:00"),
    tz = "UTC"
  )
  target <- as.POSIXct("2020-01-02 19:54:00", tz = "UTC")

  expect_identical(.match_time_index(available, target), 2L)
  expect_identical(
    .match_time_index(available, target, time_match = "exact"),
    NA_integer_
  )
  expect_identical(
    .match_time_index(
      available, target, max_time_diff = 5 / 1440
    ),
    NA_integer_
  )
  expect_identical(
    .match_time_index(
      available, target, max_time_diff = 7 / 1440
    ),
    2L
  )
})

test_that("monthly and yearly selection compares calendar labels", {
  months <- as.POSIXct(
    c("2020-02-15", "2020-03-15", "2020-04-15"),
    tz = "UTC"
  )

  expect_identical(
    .match_time_index(
      months, as.Date("2020-03-31"),
      time_match = "exact",
      time_resolution = "monthly"
    ),
    2L
  )

  years <- as.POSIXct(
    c("2019-07-01", "2020-07-01", "2021-07-01"),
    tz = "UTC"
  )

  expect_identical(
    .match_time_index(
      years, as.Date("2020-12-31"),
      time_match = "exact",
      time_resolution = "yearly"
    ),
    2L
  )
})

test_that("coarse selection rejects multiple records in one period", {
  available <- as.POSIXct(
    c("2020-01-02 19:00:00", "2020-01-02 20:00:00"),
    tz = "UTC"
  )

  expect_error(
    .match_time_index(available, as.Date("2020-01-02")),
    "Multiple dataset timestamps", class = "ocean3d_error_repeated_time"
  )
  for (resolution in c("monthly", "yearly")) {
    expect_error(
      .match_time_index(
        available, as.Date("2020-01-02"),
        time_resolution = resolution
      ),
      "Multiple dataset timestamps", class = "ocean3d_error_repeated_time"
    )
  }
})

test_that("temporal eligibility needs a target and a non-missing time axis", {
  axis <- as.POSIXct(c(NA_character_, "2020-01-01"), tz = "UTC")
  target <- as.POSIXct("2020-01-10", tz = "UTC")
  expect_true(.temporal_source_eligible(axis, target, "daily"))
  # Static sources have no time axis to restrict them.
  expect_true(.temporal_source_eligible(NULL, target, "daily"))
  expect_false(.temporal_source_eligible(axis, target[NA], "daily"))
  expect_false(.temporal_source_eligible(axis[1], target, "daily"))
  expect_false(.temporal_source_eligible(axis[0], target, "daily"))
})

# Time metadata for point extraction (.read_temporal_metadata()) ---------------

test_that("metadata keeps bounds separate from sampling cadence", {
  f <- make_grid_file(depth = NULL, time = 0)
  on.exit(unlink(f), add = TRUE)
  set_atts(f, "temp", cell_methods = "time: mean")

  info <- read_metadata(f)
  expect_identical(info$cell_methods, "time: mean")
  expect_identical(info$cadence$type, "unknown")
  expect_identical(info$calendar, "standard")
  expect_match(info$units, "days since")
})

test_that("explicit bounds describe intervals without inventing cadence", {
  f <- make_grid_file(depth = NULL, time = 0,
                      time_bounds = matrix(c(-0.5, 0.5), nrow = 2L))
  on.exit(unlink(f), add = TRUE)
  expect_identical(read_metadata(f)$cadence$type, "unknown")
})
