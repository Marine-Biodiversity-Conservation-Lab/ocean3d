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
  expect_error(decode_time(several), class = "ocean3d_error_time_units")
  err <- tryCatch(decode_time(several), error = conditionMessage)
  expect_match(err, "months since 2020-01-01", fixed = TRUE)
  expect_match(err, "several months or years is ambiguous", fixed = TRUE)
})

test_that("reference times with a numeric time zone offset are rejected", {
  f <- make_grid_file(depth = NULL, time = c(0, 24),
                      time_units = "hours since 2020-01-01 05:00:00+05:00")
  on.exit(unlink(f), add = TRUE)

  expect_error(decode_time(f), class = "ocean3d_error_time_units")
  err <- tryCatch(decode_time(f), error = conditionMessage)
  expect_match(err, "time zone offset +0500", fixed = TRUE)
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
    expect_error(decode_time(f), class = "ocean3d_error_calendar",
                 label = calendar)
    err <- tryCatch(decode_time(f), error = conditionMessage)
    expect_match(err, paste("unsupported calendar:", calendar), fixed = TRUE)
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

  expect_error(match_times(axis, as.Date("2020-01-02"), file = basename(f)),
               class = "ocean3d_error_repeated_time")
  err <- tryCatch(
    match_times(axis, as.Date("2020-01-02"), file = basename(f)),
    error = conditionMessage
  )
  expect_match(err, "Multiple dataset timestamps", fixed = TRUE)
  expect_match(err, "2020-01-02 00:00:00 UTC appears 2 times", fixed = TRUE)
  expect_match(err, basename(f), fixed = TRUE)
})

test_that("midnight date-times keep their day on noon-stamped files (#71)", {
  f <- make_grid_file(depth = NULL, time = c(613620, 613644, 613668),
                      time_units = "hours since 1950-01-01")
  on.exit(unlink(f), add = TRUE)
  times <- utc(c("2020-01-01", "2020-01-02", "2020-01-03"))

  axis <- decode_time(f)
  expect_equal(axis, times + 12 * 3600)
  expect_identical(match_times(axis, times), 1:3)
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

  expect_error(decode_time(f, bounds = bad),
               class = "ocean3d_error_time_bounds")
  err <- tryCatch(decode_time(f, bounds = bad), error = conditionMessage)
  expect_match(err, "Invalid bounds for `time`", fixed = TRUE)
  expect_match(err, "must contain its time step", fixed = TRUE)

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
