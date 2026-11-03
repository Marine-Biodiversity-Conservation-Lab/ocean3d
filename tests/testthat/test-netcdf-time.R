# CF time decoding and matching (#88-#92, #71, #111, #113), tested on the
# shared helpers directly: the core object (#113) reuses them, while the
# extraction functions that call them may be retired (#116; AGENTS.md,
# "Tests"). Each file has a time axis and an optional 2 x n matrix of time
# bounds in the axis units.

make_time_file <- function(time_vals, units, calendar = NULL, unlim = FALSE,
                           bounds = NULL) {
  f <- tempfile(fileext = ".nc")
  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  tm <- ncdf4::ncdim_def("time", units, time_vals, unlim = unlim)
  vars <- list(ncdf4::ncvar_def(
    "temp", "1", list(lon, lat, tm), -9999, prec = "double"
  ))
  if (!is.null(bounds)) {
    nv <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
    vars <- c(vars, list(ncdf4::ncvar_def(
      "time_bounds", units, list(nv, tm), prec = "double"
    )))
  }
  nc <- ncdf4::nc_create(f, vars)
  nt <- length(time_vals)
  ncdf4::ncvar_put(nc, "temp", array(0, c(2, 2, nt)), count = c(2, 2, nt))
  if (!is.null(calendar)) ncdf4::ncatt_put(nc, "time", "calendar", calendar)
  if (!is.null(bounds)) {
    ncdf4::ncvar_put(nc, "time_bounds", bounds)
    ncdf4::ncatt_put(nc, "time", "bounds", "time_bounds")
  }
  ncdf4::nc_close(nc)
  f
}

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
  one <- make_time_file(0, "months since 1955-01-01 00:00:00")
  several <- make_time_file(c(0, 1, 2), "months since 2020-01-01")
  on.exit(unlink(c(one, several)), add = TRUE)

  expect_equal(decode_time(one), utc("1955-01-01"))
  expect_error(decode_time(several), class = "ocean3d_error_time_units")
  err <- tryCatch(decode_time(several), error = conditionMessage)
  expect_match(err, "months since 2020-01-01", fixed = TRUE)
  expect_match(err, "several months or years is ambiguous", fixed = TRUE)
})

test_that("reference times with a numeric time zone offset are rejected", {
  f <- make_time_file(c(0, 24), "hours since 2020-01-01 05:00:00+05:00")
  on.exit(unlink(f), add = TRUE)

  expect_error(decode_time(f), class = "ocean3d_error_time_units")
  err <- tryCatch(decode_time(f), error = conditionMessage)
  expect_match(err, "time zone offset +0500", fixed = TRUE)
})

test_that("fractional-day timestamps decode and match exactly (#90)", {
  base <- as.numeric(as.Date("2020-06-15"))
  f <- make_time_file(base + (0:143) / 144, "days since 1970-01-01")
  on.exit(unlink(f), add = TRUE)
  times <- utc("2020-06-15") + (0:143) * 600

  axis <- decode_time(f)
  expect_equal(axis, times)
  expect_identical(match_times(axis, times, time_match = "exact"), 1:144)
})

test_that("calendars other than the standard ones stop with a clear error", {
  for (calendar in c("noleap", "360_day", "all_leap", "julian")) {
    f <- make_time_file(c(0, 1), "days since 2000-01-01", calendar)
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
    f <- make_time_file(c(0, 2, 5), "days since 2020-01-01", calendar)
    on.exit(unlink(f), add = TRUE)
    expect_equal(decode_time(f), expected, label = calendar)
  }
})

test_that("a readunlim = FALSE connection reads the time axis (#92)", {
  f <- make_time_file(c(0, 2, 5), "days since 2020-01-01", unlim = TRUE)
  on.exit(unlink(f), add = TRUE)

  expect_equal(decode_time(f, readunlim = FALSE),
               utc(c("2020-01-01", "2020-01-03", "2020-01-06")))
})

test_that("many dates match a long daily axis in one pass (#89)", {
  f <- make_time_file(0:3649, "days since 2010-01-01")
  on.exit(unlink(f), add = TRUE)
  set.seed(1)
  dates <- as.Date("2010-01-01") + sample(0:3649, 500, replace = TRUE)

  expect_identical(match_times(decode_time(f), dates),
                   as.integer(dates - as.Date("2010-01-01")) + 1L)
})

test_that("a repeated time step is only an error if it is selected (#91)", {
  f <- make_time_file(c(0, 1, 1, 2), "days since 2020-01-01")
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
  f <- make_time_file(c(613620, 613644, 613668), "hours since 1950-01-01")
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
  f <- make_time_file(c(15, 45), "days since 2020-01-01",
                      bounds = monthly_bounds)
  on.exit(unlink(f), add = TRUE)

  expect_no_warning(edges <- decode_time(f, bounds = monthly_bounds))
  expect_equal(edges, month_edges)

  # The bounds reach the time metadata as the month edges, even when warnings
  # are errors.
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  info <- read_metadata(f)
  expect_identical(info$bounds_status, "read")
  expect_equal(info$bounds, month_edges)
})

test_that("time bounds may list the later edge first (#113)", {
  f <- make_time_file(c(15, 45), "days since 2020-01-01",
                      bounds = monthly_bounds[2:1, ])
  on.exit(unlink(f), add = TRUE)

  expect_no_warning(edges <- decode_time(f, bounds = monthly_bounds[2:1, ]))
  expect_equal(edges, month_edges[2:1, ])
})

test_that("time bounds that miss their time step are rejected (#113)", {
  # January's bounds end on 10 January, before its 16 January time step.
  bad <- matrix(c(0, 10, 31, 60), nrow = 2L)
  f <- make_time_file(c(15, 45), "days since 2020-01-01", bounds = bad)
  on.exit(unlink(f), add = TRUE)

  expect_error(decode_time(f, bounds = bad),
               class = "ocean3d_error_time_bounds")
  err <- tryCatch(decode_time(f, bounds = bad), error = conditionMessage)
  expect_match(err, "Invalid bounds for `time`", fixed = TRUE)
  expect_match(err, "must contain its time step", fixed = TRUE)

  # The time metadata marks them unavailable rather than stopping, as decided
  # on 3 November 2026 (#123).
  info <- read_metadata(f)
  expect_identical(info$bounds_status, "unavailable")
  expect_null(info$bounds)
})

test_that("the time axis builder returns a CFtime object with its bounds", {
  f <- make_time_file(c(15, 45), "days since 2020-01-01",
                      bounds = monthly_bounds[2:1, ])
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
