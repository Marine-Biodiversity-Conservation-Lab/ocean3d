# CF time decoding through the exported extraction functions (#88, #90, #92,
# #111, #113). Each file is a 2 x 2 grid with a time axis; `temp` at time step
# t equals t, so each extracted value is the index of the selected time step.
# `bounds` is an optional 2 x n matrix of time bounds in the axis units.

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
  ncdf4::ncvar_put(nc, "temp", array(rep(seq_len(nt), each = 4), c(2, 2, nt)),
                   count = c(2, 2, nt))
  if (!is.null(calendar)) ncdf4::ncatt_put(nc, "time", "calendar", calendar)
  if (!is.null(bounds)) {
    ncdf4::ncvar_put(nc, "time_bounds", bounds)
    ncdf4::ncatt_put(nc, "time", "bounds", "time_bounds")
  }
  ncdf4::nc_close(nc)
  f
}

at_date <- function(f, date) {
  out <- extract2d(data.frame(lon = 0, lat = 40, date = date), f,
                   var = "temp", time_report = FALSE)
  out$temp
}

test_that("a single-step months axis is read and a longer one rejected (#88)", {
  one <- make_time_file(0, "months since 1955-01-01 00:00:00")
  several <- make_time_file(c(0, 1, 2), "months since 2020-01-01")
  on.exit(unlink(c(one, several)), add = TRUE)

  # The only time step is 1955, so the temporal warning is expected.
  expect_warning(
    value <- at_date(one, as.Date("2020-01-03")),
    class = "ocean3d_time_warning"
  )
  expect_equal(value, 1)
  expect_error(
    at_date(several, as.Date("2020-01-03")),
    class = "ocean3d_error_time_units"
  )
  err <- tryCatch(at_date(several, as.Date("2020-01-03")),
                  error = function(e) conditionMessage(e))
  expect_match(err, "months since 2020-01-01", fixed = TRUE)
  expect_match(err, "several months or years is ambiguous", fixed = TRUE)
})

test_that("reference times with a numeric time zone offset are rejected", {
  f <- make_time_file(c(0, 24), "hours since 2020-01-01 05:00:00+05:00")
  on.exit(unlink(f), add = TRUE)

  expect_error(at_date(f, as.Date("2020-01-02")),
               class = "ocean3d_error_time_units")
  err <- tryCatch(at_date(f, as.Date("2020-01-02")),
                  error = function(e) conditionMessage(e))
  expect_match(err, "time zone offset +0500", fixed = TRUE)
})

test_that("fractional-day timestamps match exactly (#90)", {
  base <- as.numeric(as.Date("2020-06-15"))
  f <- make_time_file(base + (0:143) / 144, "days since 1970-01-01")
  on.exit(unlink(f), add = TRUE)
  times <- as.POSIXct("2020-06-15", tz = "UTC") + (0:143) * 600

  out <- extract2d(
    data.frame(lon = 0, lat = 40, date = times), f, var = "temp",
    time_match = "exact", time_resolution = "datetime", time_report = FALSE
  )
  expect_equal(out$temp, as.numeric(1:144))
})

test_that("calendars other than the standard ones stop with a clear error", {
  for (calendar in c("noleap", "360_day", "all_leap", "julian")) {
    f <- make_time_file(c(0, 1), "days since 2000-01-01", calendar)
    on.exit(unlink(f), add = TRUE)
    expect_error(at_date(f, as.Date("2000-01-02")),
                 class = "ocean3d_error_calendar", label = calendar)
    err <- tryCatch(at_date(f, as.Date("2000-01-02")),
                    error = function(e) conditionMessage(e))
    expect_match(err, paste("unsupported calendar:", calendar), fixed = TRUE)
  }
})

test_that("the standard calendars are accepted under each name", {
  for (calendar in c("standard", "gregorian", "proleptic_gregorian",
                     "Gregorian")) {
    f <- make_time_file(c(0, 2, 5), "days since 2020-01-01", calendar)
    on.exit(unlink(f), add = TRUE)
    expect_equal(at_date(f, as.Date("2020-01-03")), 2, label = calendar)
  }
})

test_that("a readunlim = FALSE connection reads the time axis (#92)", {
  f <- make_time_file(c(0, 2, 5), "days since 2020-01-01", unlim = TRUE)
  on.exit(unlink(f), add = TRUE)
  con <- ncdf4::nc_open(f, readunlim = FALSE)
  on.exit(ncdf4::nc_close(con), add = TRUE)

  expect_equal(at_date(f, as.Date("2020-01-03")), 2)
  expect_equal(at_date(con, as.Date("2020-01-03")), 2)
})

test_that("many dates match a long daily axis in one pass (#89)", {
  f <- make_time_file(0:3649, "days since 2010-01-01")
  on.exit(unlink(f), add = TRUE)
  set.seed(1)
  dates <- as.Date("2010-01-01") + sample(0:3649, 500, replace = TRUE)

  out <- extract2d(data.frame(lon = 0, lat = 40, date = dates), f,
                   var = "temp", time_report = FALSE)
  expect_equal(out$temp, as.numeric(dates - as.Date("2010-01-01")) + 1)
})

test_that("a repeated time step only stops extraction if requested (#91)", {
  f <- make_time_file(c(0, 1, 1, 2), "days since 2020-01-01")
  on.exit(unlink(f), add = TRUE)
  # CFtime warns about the repeated offset; that warning is not silenced.
  repeat_warning <- "Offsets not monotonically increasing"

  expect_warning(first <- at_date(f, as.Date("2020-01-01")),
                 repeat_warning, fixed = TRUE)
  expect_equal(first, 1)
  expect_warning(last <- at_date(f, as.Date("2020-01-03")),
                 repeat_warning, fixed = TRUE)
  expect_equal(last, 4)

  expect_error(suppressWarnings(at_date(f, as.Date("2020-01-02"))),
               class = "ocean3d_error_repeated_time")
  err <- tryCatch(suppressWarnings(at_date(f, as.Date("2020-01-02"))),
                  error = function(e) conditionMessage(e))
  expect_match(err, "Multiple dataset timestamps", fixed = TRUE)
  expect_match(err, "2020-01-02 00:00:00 UTC appears 2 times", fixed = TRUE)
  expect_match(err, basename(f), fixed = TRUE)
})

test_that("midnight date-times keep their day on noon-stamped files (#71)", {
  f <- make_time_file(c(613620, 613644, 613668), "hours since 1950-01-01")
  on.exit(unlink(f), add = TRUE)
  times <- as.POSIXct(c("2020-01-01", "2020-01-02", "2020-01-03"), tz = "UTC")

  out <- extract2d(data.frame(lon = 0, lat = 40, date = times), f,
                   var = "temp", time_report = FALSE)
  expect_equal(out$temp, c(1, 2, 3))
})

# Monthly means stamped mid-month, with back-to-back bounds: January and
# February 2020 share the 1 February edge (offsets 0, 31, 31, 60).
monthly_bounds <- matrix(c(0, 31, 31, 60), nrow = 2L)

buffer_at <- function(f, date) {
  out <- extract_to_buffer(
    data.frame(lon = 0, lat = 40, date = date), f, var = "temp",
    time_before = 10, time_after = 10, time_unit = "days",
    infer_bounds = TRUE, output_col = "temp", quiet = TRUE
  )
  out$temp
}

test_that("back-to-back time bounds decode without a warning (#113)", {
  f <- make_time_file(c(15, 45), "days since 2020-01-01",
                      bounds = monthly_bounds)
  on.exit(unlink(f), add = TRUE)

  expect_no_warning(point <- at_date(f, as.Date("2020-02-15")))
  expect_equal(point, 2)
  expect_no_warning(point <- extract_to_point(
    data.frame(lon = 0, lat = 40, date = as.Date("2020-01-16")), f,
    var = "temp", time_report = FALSE
  ))
  expect_equal(point$temp, 1)

  # 21 January to 10 February: 11 days of January (1) and 9 of February (2).
  expect_no_warning(value <- buffer_at(f, "2020-01-31"))
  expect_equal(value, (11 * 1 + 9 * 2) / 20)

  # The bounds reach the time report as the month edges, even when warnings
  # are errors.
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  info <- .read_temporal_metadata(
    nc, "time", "temp", .decode_netcdf_time(nc, "time")
  )
  expect_identical(info$bounds_status, "read")
  edges <- as.numeric(as.POSIXct(
    c("2020-01-01", "2020-02-01", "2020-02-01", "2020-03-01"), tz = "UTC"
  ))
  expect_equal(info$bounds, matrix(edges, nrow = 2L))
})

test_that("time bounds may list the later edge first (#113)", {
  f <- make_time_file(c(15, 45), "days since 2020-01-01",
                      bounds = monthly_bounds[2:1, ])
  on.exit(unlink(f), add = TRUE)

  expect_no_warning(value <- buffer_at(f, "2020-01-31"))
  expect_equal(value, (11 * 1 + 9 * 2) / 20)
})

test_that("time bounds that miss their time step stop buffer extraction", {
  # January's bounds end on 10 January, before its 16 January time step.
  f <- make_time_file(c(15, 45), "days since 2020-01-01",
                      bounds = matrix(c(0, 10, 31, 60), nrow = 2L))
  on.exit(unlink(f), add = TRUE)

  expect_error(buffer_at(f, "2020-01-31"),
               class = "ocean3d_error_time_bounds")
  err <- tryCatch(buffer_at(f, "2020-01-31"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "Invalid bounds for `time`", fixed = TRUE)
  expect_match(err, "must contain its time step", fixed = TRUE)

  # Point extraction doesn't select values by time bounds, so it carries on.
  expect_no_warning(point <- at_date(f, as.Date("2020-02-15")))
  expect_equal(point, 2)
})
