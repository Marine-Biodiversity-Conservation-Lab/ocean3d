# Helpers ----------------------------------------------------------------------

.buffer_test_rectangle <- function(xmin, xmax, ymin, ymax) {
  sf::st_polygon(list(rbind(
    c(xmin, ymin),
    c(xmax, ymin),
    c(xmax, ymax),
    c(xmin, ymax),
    c(xmin, ymin)
  )))
}


.make_buffer_grid_nc <- function(
  path, cell_values = matrix(c(0, 10, 20, 30), 2L),
  lon = c(-0.01, 0.01), lat = c(-0.01, 0.01),
  lon_bounds = matrix(c(-0.02, 0, 0, 0.02), 2L),
  lat_bounds = matrix(c(-0.02, 0, 0, 0.02), 2L),
  var = "temp", spatial_bounds = TRUE,
  field = matrix(0), depth = NULL, time = NULL,
  depth_bounds = NULL, time_bounds = NULL
) {
  x <- ncdf4::ncdim_def("lon", "degrees_east", lon)
  y <- ncdf4::ncdim_def("lat", "degrees_north", lat)
  nv <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
  z <- if (is.null(depth)) NULL else ncdf4::ncdim_def("depth", "m", depth)
  t <- if (is.null(time)) NULL else ncdf4::ncdim_def(
    "time", "hours since 2020-01-01 00:00:00", time
  )
  dims <- list()
  if (!is.null(t)) dims <- c(dims, list(t))
  dims <- c(dims, list(y))
  if (!is.null(z)) dims <- c(dims, list(z))
  dims <- c(dims, list(x))
  variables <- list(ncdf4::ncvar_def(var, "degree_Celsius", dims, -9999, prec = "double"))
  bounds <- list()
  if (spatial_bounds) bounds <- list(lon = lon_bounds, lat = lat_bounds)
  if (!is.null(depth_bounds)) bounds$depth <- depth_bounds
  if (!is.null(time_bounds)) bounds$time <- time_bounds
  dimensions <- list(lon = x, lat = y, depth = z, time = t)
  for (name in names(bounds)) {
    variables <- c(variables, list(ncdf4::ncvar_def(
      paste0(name, "_bounds"), dimensions[[name]]$units,
      list(nv, dimensions[[name]]), prec = "double"
    )))
  }
  handle <- ncdf4::nc_create(path, variables)
  on.exit(ncdf4::nc_close(handle))
  nz <- if (is.null(z)) 1L else length(depth)
  nt <- if (is.null(t)) 1L else length(time)
  stopifnot(identical(dim(field), c(nz, nt)))
  values <- array(NA_real_, c(nt, length(lat), nz, length(lon)))
  for (i in seq_along(lon)) for (j in seq_along(lat)) {
    for (k in seq_len(nz)) for (l in seq_len(nt)) {
      values[l, j, k, i] <- field[k, l] + cell_values[i, j]
    }
  }
  shape <- vapply(dims, function(d) length(d$vals), integer(1))
  ncdf4::ncvar_put(handle, var, array(values, shape))
  for (name in names(bounds)) {
    ncdf4::ncvar_put(handle, paste0(name, "_bounds"), bounds[[name]])
    ncdf4::ncatt_put(handle, name, "bounds", paste0(name, "_bounds"))
  }
  if (!is.null(z)) ncdf4::ncatt_put(handle, "depth", "positive", "down")
  invisible(path)
}


.make_buffer_axis_nc <- function(
  path,
  centres,
  bounds = NULL,
  axis = "depth",
  units = if (axis == "depth") "m" else "hours since 2020-01-01 00:00:00",
  attributes = list(),
  bounds_units = units,
  bounds_reference = "axis_bounds",
  coordinate_first = FALSE,
  bounds_attributes = list()
) {
  dimension <- ncdf4::ncdim_def(axis, units, centres)
  value <- ncdf4::ncvar_def("value", "1", list(dimension), prec = "double")
  variables <- list(value)
  if (!is.null(bounds)) {
    vertices <- ncdf4::ncdim_def(
      "nv", "", seq_len(nrow(bounds)), create_dimvar = FALSE
    )
    dimensions <- if (coordinate_first) {
      list(dimension, vertices)
    } else {
      list(vertices, dimension)
    }
    variables <- c(variables, list(ncdf4::ncvar_def(
      "axis_bounds", bounds_units, dimensions, missval = -9999,
      prec = "double"
    )))
  }
  nc <- ncdf4::nc_create(path, variables)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_put(nc, "value", seq_along(centres))
  if (!is.null(bounds)) {
    ncdf4::ncvar_put(
      nc, "axis_bounds", if (coordinate_first) t(bounds) else bounds
    )
  }
  if (!is.null(bounds) || bounds_reference != "axis_bounds") {
    ncdf4::ncatt_put(nc, axis, "bounds", bounds_reference)
  }
  for (name in names(attributes)) {
    ncdf4::ncatt_put(nc, axis, name, attributes[[name]])
  }
  for (name in names(bounds_attributes)) {
    ncdf4::ncatt_put(nc, "axis_bounds", name, bounds_attributes[[name]])
  }
  invisible(path)
}


.read_buffer_axis_nc <- function(path, coordinate = "depth", ...) {
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  .buffer_read_axis(nc, coordinate, ...)
}


.make_buffer_source_nc <- function(
  path,
  field,
  depth = NULL,
  time = NULL,
  depth_bounds = NULL,
  time_bounds = NULL
) {
  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  vertices <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
  depth_dim <- if (is.null(depth)) {
    NULL
  } else {
    ncdf4::ncdim_def("depth", "m", depth)
  }
  time_dim <- if (is.null(time)) {
    NULL
  } else {
    ncdf4::ncdim_def("time", "hours since 2020-01-01 00:00:00", time)
  }
  n_depth <- if (is.null(depth)) 1L else length(depth)
  n_time <- if (is.null(time)) 1L else length(time)
  stopifnot(is.matrix(field), identical(dim(field), c(n_depth, n_time)))

  # Deliberately use a different order from lon/lat/depth/time.
  dims <- list()
  if (!is.null(time_dim)) dims <- c(dims, list(time_dim))
  dims <- c(dims, list(lat))
  if (!is.null(depth_dim)) dims <- c(dims, list(depth_dim))
  dims <- c(dims, list(lon))
  value <- ncdf4::ncvar_def(
    "temp", "degree_Celsius", dims, missval = -9999, prec = "double"
  )
  variables <- list(value)
  if (!is.null(depth_bounds)) {
    variables <- c(variables, list(ncdf4::ncvar_def(
      "depth_bounds", "m", list(vertices, depth_dim), prec = "double"
    )))
  }
  if (!is.null(time_bounds)) {
    variables <- c(variables, list(ncdf4::ncvar_def(
      "time_bounds", "hours since 2020-01-01 00:00:00",
      list(vertices, time_dim), prec = "double"
    )))
  }
  nc <- ncdf4::nc_create(path, variables)
  on.exit(ncdf4::nc_close(nc))
  values <- array(NA_real_, c(n_time, 2L, n_depth, 2L))
  for (x in 1:2) {
    for (y in 1:2) {
      for (z in seq_len(n_depth)) {
        for (t in seq_len(n_time)) {
          values[t, y, z, x] <- field[z, t] +
            c(0, 100)[x] + c(0, 200)[y]
        }
      }
    }
  }
  shape <- vapply(dims, function(x) length(x$vals), integer(1))
  ncdf4::ncvar_put(nc, "temp", array(values, shape))
  if (!is.null(depth_dim)) ncdf4::ncatt_put(nc, "depth", "positive", "down")
  if (!is.null(depth_bounds)) {
    ncdf4::ncvar_put(nc, "depth_bounds", depth_bounds)
    ncdf4::ncatt_put(nc, "depth", "bounds", "depth_bounds")
  }
  if (!is.null(time_bounds)) {
    ncdf4::ncvar_put(nc, "time_bounds", time_bounds)
    ncdf4::ncatt_put(nc, "time", "bounds", "time_bounds")
  }
  invisible(path)
}


.read_buffer_source_nc <- function(source, data, options, ...) {
  .buffer_extract_one_source(
    data, .as_netcdf_sources(source)[[1L]], "temp", options, ...
  )
}


# Request validation -----------------------------------------------------------

test_that("buffers require a positive window and valid scalar sizes", {
  expect_error(.validate_buffer_options(), "extract_to_point", fixed = TRUE)

  sizes <- c(
    "horizontal_radius_km", "depth_radius_m", "time_before", "time_after"
  )
  for (name in sizes) {
    for (value in list(-1, NA_real_, Inf, "5", c(1, 2), numeric(), NULL)) {
      args <- list(horizontal_radius_km = 5)
      args[name] <- list(value)
      expect_error(
        do.call(.validate_buffer_options, args), name, fixed = TRUE
      )
    }
  }

  options <- .validate_buffer_options(time_before = 0.5)
  expect_equal(options$time_before, 0.5)
  expect_identical(options$time_unit, "days")
  expect_identical(options$weighting, "weighted")
  expect_false(options$infer_bounds)
})


test_that("statistics, coverage thresholds and switches are validated", {
  expect_error(
    .validate_buffer_options(horizontal_radius_km = 5, statistic = "sum"),
    "`statistic` must be one of", fixed = TRUE
  )
  expect_error(
    .validate_buffer_options(horizontal_radius_km = 5, weighting = "other")
  )

  for (value in list(-0.1, 1.1, NA_real_, Inf, c(0.5, 0.8), NULL)) {
    expect_error(
      .validate_buffer_options(
        horizontal_radius_km = 5, min_valid_fraction = value
      ),
      "min_valid_fraction", fixed = TRUE
    )
  }

  for (name in c("infer_bounds", "return_match_info", "quiet")) {
    for (value in list(NA, 1, c(TRUE, FALSE), NULL)) {
      args <- list(horizontal_radius_km = 5)
      args[name] <- list(value)
      expect_error(
        do.call(.validate_buffer_options, args), name, fixed = TRUE
      )
    }
  }

  for (threshold in c(0, 0.8, 1)) {
    options <- .validate_buffer_options(
      horizontal_radius_km = 5,
      min_valid_fraction = threshold,
      weighting = "equal"
    )
    expect_equal(options$min_valid_fraction, threshold)
    expect_identical(options$weighting, "equal")
  }
})


test_that("calendar sizes are whole numbers; fixed durations allow fractions", {
  for (unit in c("months", "years")) {
    for (name in c("time_before", "time_after")) {
      args <- list(time_unit = unit)
      args[[name]] <- 0.5
      expect_error(
        do.call(.validate_buffer_options, args),
        paste0("`", name, "` must be a whole number"), fixed = TRUE
      )
    }
    options <- .validate_buffer_options(time_before = 1, time_unit = unit)
    expect_identical(options$time_unit, unit)
  }
  for (unit in c("days", "hours", "minutes", "seconds")) {
    options <- .validate_buffer_options(time_after = 0.5, time_unit = unit)
    expect_equal(options$time_after, 0.5)
  }
  expect_error(.validate_buffer_options(time_before = 1, time_unit = "weeks"))
})


test_that("depth intervals and single-layer methods are distinct", {
  options <- .validate_buffer_options(depth_radius_m = 10)
  expect_identical(options$method, "interval")

  for (method in c("nearest", "surface", "bottom", "2d")) {
    expect_error(
      .validate_buffer_options(depth_radius_m = 10, method = method),
      "Leave `method = NULL`", fixed = TRUE
    )
    options <- .validate_buffer_options(
      horizontal_radius_km = 5, method = method
    )
    expect_identical(options$method, method)
  }
  options <- .validate_buffer_options(horizontal_radius_km = 5)
  expect_identical(options$method, "nearest")
  expect_error(
    .validate_buffer_options(horizontal_radius_km = 5, method = "other")
  )
})


# Observation windows ---------------------------------------------------------

test_that("depth windows preserve requested support and missing rows", {
  data <- data.frame(depth = c(60, 5, NA_real_))
  original <- data
  options <- .validate_buffer_options(depth_radius_m = 10)
  result <- .buffer_observation_windows(data, options)

  expect_equal(result$windows$row, 1:3)
  expect_equal(result$windows$depth_lower_m, c(50, -5, NA_real_))
  expect_equal(result$windows$depth_upper_m, c(70, 15, NA_real_))
  expect_true(all(is.na(result$windows$time_lower_utc)))
  expect_null(result$prepared_time)
  expect_identical(data, original)

  expect_error(
    .buffer_observation_windows(data.frame(lon = 0), options),
    "depth", fixed = TRUE
  )
  expect_error(
    .buffer_observation_windows(data, options, depth_col = NULL),
    "depth_col", fixed = TRUE
  )
  for (depth in list("60", -1, Inf)) {
    expect_error(
      .buffer_observation_windows(data.frame(depth = depth), options),
      "depth", fixed = TRUE
    )
  }
  custom <- .buffer_observation_windows(
    data.frame(z = 60), options, depth_col = "z"
  )
  expect_equal(custom$windows$depth_lower_m, 50)
  expect_equal(custom$windows$depth_upper_m, 70)
})


test_that("time windows use elapsed hours across daylight saving changes", {
  data <- data.frame(
    observed_at = as.POSIXct(
      c("2020-03-29 03:30:00", NA_character_), tz = "Europe/Madrid"
    )
  )
  original <- data
  for (unit in c("days", "hours")) {
    options <- .validate_buffer_options(
      time_before = if (unit == "days") 0.5 else 12,
      time_after = if (unit == "days") 0.25 else 6,
      time_unit = unit
    )
    result <- .buffer_observation_windows(
      data, options, date_col = "observed_at"
    )
    expect_equal(
      format(result$windows$time_lower_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
      c("2020-03-28 13:30", NA_character_)
    )
    expect_equal(
      format(result$windows$time_upper_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
      c("2020-03-29 07:30", NA_character_)
    )
    expect_true(all(is.na(result$windows$depth_lower_m)))
    expect_identical(data, original)
  }
})


test_that("minutes and seconds preserve sub-minute window boundaries", {
  data <- data.frame(date = "2020-01-01T12:00:00Z")
  for (unit in c("minutes", "seconds")) {
    options <- .validate_buffer_options(
      time_before = if (unit == "minutes") 0.5 else 30,
      time_after = if (unit == "minutes") 1.5 else 90,
      time_unit = unit
    )
    result <- .buffer_observation_windows(data, options)
    expect_equal(
      format(result$windows$time_lower_utc, "%H:%M:%S", tz = "UTC"),
      "11:59:30"
    )
    expect_equal(
      format(result$windows$time_upper_utc, "%H:%M:%S", tz = "UTC"),
      "12:01:30"
    )
  }
})


test_that("calendar months clamp month ends and preserve UTC clock time", {
  data <- data.frame(date = c(
    "2020-03-31T12:34:56Z", "2021-03-31T12:34:56Z",
    "2020-01-31T12:34:56Z", NA_character_
  ))
  options <- .validate_buffer_options(
    time_before = 1, time_after = 1, time_unit = "months"
  )
  result <- .buffer_observation_windows(data, options)
  expect_equal(
    format(result$windows$time_lower_utc, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    c("2020-02-29 12:34:56", "2021-02-28 12:34:56",
      "2019-12-31 12:34:56", NA_character_)
  )
  expect_equal(
    format(result$windows$time_upper_utc, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    c("2020-04-30 12:34:56", "2021-04-30 12:34:56",
      "2020-02-29 12:34:56", NA_character_)
  )
})


test_that("calendar years handle leap days and preserve unchanged endpoints", {
  options <- .validate_buffer_options(time_before = 1, time_unit = "years")
  result <- .buffer_observation_windows(
    data.frame(date = "2020-02-29T19:30:00Z"), options
  )
  expect_equal(
    format(result$windows$time_lower_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
    "2019-02-28 19:30"
  )
  expect_equal(
    format(result$windows$time_upper_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
    "2020-02-29 19:30"
  )
  options <- .validate_buffer_options(time_after = 1, time_unit = "years")
  result <- .buffer_observation_windows(
    data.frame(date = "2020-02-29T19:30:00Z"), options
  )
  expect_equal(
    format(result$windows$time_upper_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
    "2021-02-28 19:30"
  )
})


test_that("Date centres use midnight and unused coordinates are ignored", {
  options <- .validate_buffer_options(time_before = 3)
  result <- .buffer_observation_windows(
    data.frame(date = as.Date("2020-01-04"), depth = "unused"), options
  )
  expect_equal(
    format(result$windows$time_lower_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
    "2020-01-01 00:00"
  )
  expect_equal(
    format(result$windows$time_upper_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
    "2020-01-04 00:00"
  )
  expect_error(
    .buffer_observation_windows(data.frame(lon = 0), options),
    "date", fixed = TRUE
  )
  expect_error(
    .buffer_observation_windows(
      data.frame(date = NA), options, date_col = NULL
    ),
    "date_col", fixed = TRUE
  )

  options <- .validate_buffer_options(horizontal_radius_km = 5)
  result <- .buffer_observation_windows(
    data.frame(date = "unused", depth = "unused"), options
  )
  expect_equal(nrow(result$windows), 1L)
  expect_null(result$prepared_time)
})


test_that("combined windows retain timezone information for a later warning", {
  options <- .validate_buffer_options(
    horizontal_radius_km = 5, depth_radius_m = 10, time_before = 1
  )
  result <- .buffer_observation_windows(
    data.frame(depth = 60, date = "2020-01-02 12:00"), options
  )
  expect_equal(result$windows$depth_lower_m, 50)
  expect_equal(result$windows$depth_upper_m, 70)
  expect_equal(
    format(result$windows$time_lower_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
    "2020-01-01 12:00"
  )
  expect_equal(
    format(result$windows$time_upper_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
    "2020-01-02 12:00"
  )
  expect_identical(result$prepared_time$unzoned_rows, 1L)
})


test_that("missing and empty rows are preserved across window types", {
  for (unit in c("days", "months", "years")) {
    options <- .validate_buffer_options(
      depth_radius_m = 10, time_before = 1, time_unit = unit
    )
    result <- .buffer_observation_windows(
      data.frame(depth = c(NA_real_, NA_real_), date = c(NA, NA)), options
    )
    expect_equal(result$windows$row, 1:2)
    expect_true(all(is.na(result$windows$depth_lower_m)))
    expect_true(all(is.na(result$windows$time_lower_utc)))
    expect_true(all(is.na(result$windows$time_upper_utc)))

    empty <- .buffer_observation_windows(
      data.frame(depth = numeric(), date = character()), options
    )
    expect_equal(nrow(empty$windows), 0L)
  }
})


test_that("temporal conversion errors retain column and observation context", {
  options <- .validate_buffer_options(time_before = 1)
  data <- data.frame(sample = c("a", "b"),
                     observed_at = c("2020-01-02", "06/01/2020"))
  error <- tryCatch(
    .buffer_observation_windows(
      data, options, date_col = "observed_at", id_col = "sample"
    ),
    error = identity
  )
  expect_s3_class(error, "error")
  expect_match(
    conditionMessage(error), "Input: `data$observed_at`", fixed = TRUE
  )
  expect_match(
    conditionMessage(error), "row/position 2 (id: \"b\")", fixed = TRUE
  )
  expect_match(conditionMessage(error), "06/01/2020", fixed = TRUE)
})


# Interval aggregation --------------------------------------------------------

test_that("interval means distinguish weighted and equal averaging", {
  weighted <- .buffer_interval_mean(
    values = c(10, 20),
    lower = c(0, 10),
    upper = c(10, 30),
    requested = c(5, 25)
  )
  equal <- .buffer_interval_mean(
    values = c(10, 20),
    lower = c(0, 10),
    upper = c(10, 30),
    requested = c(5, 25),
    weighting = "equal"
  )

  expect_equal(weighted$value, 17.5)
  expect_equal(equal$value, 15)
  expect_equal(weighted$coverage_fraction, 1)
  expect_equal(weighted$valid_fraction, 1)
  expect_equal(weighted$n_intersecting, 2)
  expect_equal(weighted$n_valid, 2)
  expect_identical(weighted$contributing_rows, 1:2)
  expect_identical(weighted$status, "complete")
})


test_that("coverage refers to the full request and accepts the threshold", {
  for (weighting in c("weighted", "equal")) {
    result <- .buffer_interval_mean(
      values = 10, lower = 0, upper = 10,
      requested = c(0, 20), weighting = weighting
    )
    expect_true(is.na(result$value))
    expect_equal(result$coverage_fraction, 0.5)
    expect_equal(result$valid_fraction, 0.5)
    expect_identical(result$status, "insufficient_coverage")

    accepted <- .buffer_interval_mean(
      values = 10, lower = 0, upper = 10,
      requested = c(0, 20), weighting = weighting,
      min_valid_fraction = 0.5
    )
    expect_equal(accepted$value, 10)
    expect_identical(accepted$status, "partial")
  }
})


test_that("coverage and valid support are reported separately", {
  result <- .buffer_interval_mean(
    values = c(NA_real_, 20),
    lower = c(0, 10), upper = c(10, 20),
    requested = c(0, 20), min_valid_fraction = 0.5
  )

  expect_equal(result$value, 20)
  expect_equal(result$coverage_fraction, 1)
  expect_equal(result$valid_fraction, 0.5)
  expect_equal(result$n_intersecting, 2)
  expect_equal(result$n_valid, 1)
  expect_identical(result$contributing_rows, 2L)
  expect_identical(result$status, "partial")
})


test_that("coverage thresholds distinguish below, equal and above", {
  for (threshold in c(0.7999, 0.8, 0.8001)) {
    result <- .buffer_interval_mean(
      values = c(10, 20),
      lower = c(0, 6), upper = c(4, 10),
      requested = c(0, 10),
      min_valid_fraction = threshold
    )
    expect_equal(result$coverage_fraction, 0.8)
    expect_equal(result$valid_fraction, 0.8)

    if (threshold <= 0.8) {
      expect_equal(result$value, 15)
      expect_identical(result$status, "partial")
    } else {
      expect_true(is.na(result$value))
      expect_identical(result$status, "insufficient_coverage")
    }
  }
})


test_that("empty and all-missing intervals never produce a mean", {
  empty <- .buffer_interval_mean(
    values = c(10, 20),
    lower = c(0, 20), upper = c(10, 30),
    requested = c(10, 20),
    min_valid_fraction = 0
  )
  expect_true(is.na(empty$value))
  expect_identical(empty$status, "empty")
  expect_equal(empty$n_intersecting, 0)
  expect_equal(empty$coverage_fraction, 0)

  missing <- .buffer_interval_mean(
    values = c(NA_real_, NaN, Inf),
    lower = c(0, 10, 20), upper = c(10, 20, 30),
    requested = c(0, 30),
    min_valid_fraction = 0
  )
  expect_true(is.na(missing$value))
  expect_identical(missing$status, "all_missing")
  expect_equal(missing$coverage_fraction, 1)
  expect_equal(missing$valid_fraction, 0)
})


# Horizontal buffer geometry ---------------------------------------------------

test_that("horizontal buffers use kilometres and preserve radial distances", {
  skip_if_not_installed("sf")
  buffer <- .buffer_horizontal_geometry(0, 0, 1)

  # A 1 km circle has approximately pi square kilometres.
  expect_equal(buffer$area_m2 / (pi * 1000^2), 1, tolerance = 1e-4)
  expect_identical(buffer$boundary_segments, 720L)

  vertices <- sf::st_coordinates(sf::st_transform(
    buffer$geometry, buffer$distance_crs
  ))
  distances <- sqrt(vertices[, 1L]^2 + vertices[, 2L]^2)
  expect_equal(
    distances, rep(1000, length(distances)), tolerance = 1e-6
  )
})


test_that("cell support distinguishes halves and quarters of the buffer", {
  skip_if_not_installed("sf")
  buffer <- .buffer_horizontal_geometry(0, 0, 1)
  cells <- sf::st_sfc(
    .buffer_test_rectangle(-2000, 0, -2000, 2000),
    .buffer_test_rectangle(0, 2000, 0, 2000),
    .buffer_test_rectangle(0, 2000, -2000, 0),
    crs = sf::st_crs(buffer$geometry)
  )

  support <- .buffer_horizontal_support(buffer, cells)
  expect_identical(support$cell_index, 1:3)
  expect_equal(support$fraction, c(0.5, 0.25, 0.25), tolerance = 1e-6)
  expect_equal(sum(support$fraction), 1, tolerance = 1e-6)
  expect_true(all(support$area_m2 > 0))
})


test_that("partial and empty horizontal support are kept separate", {
  skip_if_not_installed("sf")
  buffer <- .buffer_horizontal_geometry(0, 0, 1)
  cells <- sf::st_sfc(
    .buffer_test_rectangle(-2000, 0, -2000, 2000),
    .buffer_test_rectangle(3000, 4000, 3000, 4000),
    crs = sf::st_crs(buffer$geometry)
  )

  partial <- .buffer_horizontal_support(buffer, cells)
  expect_identical(partial$cell_index, 1L)
  expect_equal(partial$fraction, 0.5, tolerance = 1e-6)

  empty <- .buffer_horizontal_support(buffer, cells[2L])
  expect_equal(nrow(empty), 0L)

  edge <- as.numeric(sf::st_bbox(buffer$geometry)["xmax"])
  touching <- sf::st_sfc(
    .buffer_test_rectangle(edge, edge + 1000, -2000, 2000),
    crs = sf::st_crs(buffer$geometry)
  )
  expect_equal(
    nrow(.buffer_horizontal_support(buffer, touching)), 0L
  )
})


test_that("overlapping cell interiors are rejected", {
  skip_if_not_installed("sf")
  buffer <- .buffer_horizontal_geometry(0, 0, 1)
  cells <- sf::st_sfc(
    .buffer_test_rectangle(-2000, 2000, -2000, 2000),
    .buffer_test_rectangle(0, 2000, -2000, 2000),
    crs = sf::st_crs(buffer$geometry)
  )

  expect_error(
    .buffer_horizontal_support(buffer, cells),
    "cell interiors overlap", fixed = TRUE
  )
})


test_that("local buffers support antimeridian and high-latitude centres", {
  skip_if_not_installed("sf")
  for (point in list(c(179.999, 0), c(-179.999, 0), c(0, 89.9))) {
    buffer <- .buffer_horizontal_geometry(point[1L], point[2L], 1)
    expect_true(all(sf::st_is_valid(buffer$geometry)))
    expect_equal(buffer$area_m2 / (pi * 1000^2), 1, tolerance = 1e-4)
  }
})


test_that("horizontal geometry does not depend on or change the s2 setting", {
  skip_if_not_installed("sf")
  original <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(original)))

  areas <- numeric()
  for (engine in c(TRUE, FALSE)) {
    suppressMessages(sf::sf_use_s2(engine))
    buffer <- .buffer_horizontal_geometry(0, 38, 5)
    expect_identical(sf::sf_use_s2(), engine)
    areas <- c(areas, buffer$area_m2)
  }
  expect_equal(areas[1L], areas[2L])
})


test_that("invalid horizontal requests and incompatible CRS are rejected", {
  skip_if_not_installed("sf")
  for (radius in c(-1, 0, NA_real_, Inf, 10000)) {
    expect_error(.buffer_horizontal_geometry(0, 0, radius))
  }
  expect_error(.buffer_horizontal_geometry(0, 91, 1), "latitude")

  buffer <- .buffer_horizontal_geometry(0, 0, 1)
  cells <- sf::st_sfc(
    .buffer_test_rectangle(-1, 1, -1, 1), crs = 4326
  )
  expect_error(
    .buffer_horizontal_support(buffer, cells),
    "same projected equal-area CRS", fixed = TRUE
  )
})


# Public extraction and complementary sources ---------------------------------

test_that("the public function averages horizontal cells and preserves inputs", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(path)
  data <- data.frame(id = "a", lon = 0, lat = 0, date = "unused")
  out <- extract_to_buffer(data, path, horizontal_radius_km = 1, quiet = TRUE)
  expect_named(out, c(names(data), "temp", "temp_match_info"))
  expect_identical(out[names(data)], data)
  expect_equal(out$temp, 15, tolerance = 1e-6)
  expect_equal(out$temp_match_info[[1L]]$coverage_fraction, 1)
  expect_equal(out$temp_match_info[[1L]]$n_valid, 4)
  expect_identical(out$temp_match_info[[1L]]$status, "complete")
  simple <- extract_to_buffer(
    data, path, horizontal_radius_km = 1, return_match_info = FALSE, quiet = TRUE
  )
  expect_named(simple, c(names(data), "temp"))
})


test_that("the public function combines area, depth and time before averaging", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(
    path, cell_values = matrix(c(0, 100, 200, 300), 2L),
    field = matrix(c(0, 10, 20, 40), 2L), depth = c(5, 20), time = c(5, 20),
    depth_bounds = matrix(c(0, 10, 10, 30), 2L),
    time_bounds = matrix(c(0, 10, 10, 30), 2L)
  )
  data <- data.frame(lon = 0, lat = 0, depth = 15, date = "2020-01-01T15:00Z")
  for (weighting in c("weighted", "equal")) {
    out <- extract_to_buffer(
      data, path, var = "temp", horizontal_radius_km = 1,
      depth_radius_m = 10, time_before = 10, time_after = 10, time_unit = "hours",
      weighting = weighting, quiet = TRUE
    )
    expect_equal(out$temp, if (weighting == "weighted") 178.125 else 167.5,
                 tolerance = 1e-6)
    expect_equal(out$temp_match_info[[1L]]$n_valid, 16)
    expect_equal(out$temp_match_info[[1L]]$valid_fraction, 1)
  }
})


test_that("temporal files complete coverage together in either input order", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  .make_buffer_source_nc(paths[1L], matrix(10), time = 5,
                         time_bounds = matrix(c(0, 10), 2L))
  .make_buffer_source_nc(paths[2L], matrix(20), time = 20,
                         time_bounds = matrix(c(10, 30), 2L))
  data <- data.frame(lon = 0, lat = 40, date = "2020-01-01T15:00Z")
  for (files in list(paths, rev(paths))) {
    out <- extract_to_buffer(
      data, files, var = "temp", time_before = 10, time_after = 10,
      time_unit = "hours", quiet = TRUE
    )
    expect_equal(out$temp, 17.5)
    expect_equal(out$temp_match_info[[1L]]$coverage_fraction, 1)
    expect_length(out$temp_match_info[[1L]]$sources_used, 2)
    expect_identical(out$temp_match_info[[1L]]$status, "complete")
  }
})


test_that("spatial tiles complete a circle without double-counting shared edges", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  .make_buffer_grid_nc(
    paths[1L], matrix(10, 1L, 2L), lon = -0.01,
    lon_bounds = matrix(c(-0.02, 0), 2L)
  )
  .make_buffer_grid_nc(
    paths[2L], matrix(30, 1L, 2L), lon = 0.01,
    lon_bounds = matrix(c(0, 0.02), 2L)
  )
  for (files in list(paths, rev(paths))) {
    out <- extract_to_buffer(
      data.frame(lon = 0, lat = 0), files, var = "temp",
      horizontal_radius_km = 1, quiet = TRUE
    )
    expect_equal(out$temp, 20, tolerance = 1e-6)
    expect_equal(out$temp_match_info[[1L]]$valid_fraction, 1)
  }
})


test_that("overlapping temporal and spatial sources are rejected explicitly", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  .make_buffer_source_nc(paths[1L], matrix(10), time = 5,
                         time_bounds = matrix(c(0, 20), 2L))
  .make_buffer_source_nc(paths[2L], matrix(20), time = 20,
                         time_bounds = matrix(c(10, 30), 2L))
  error <- tryCatch(extract_to_buffer(
    data.frame(lon = 0, lat = 40, date = "2020-01-01T15:00Z"),
    paths, var = "temp", time_before = 10, time_after = 10,
    time_unit = "hours", quiet = TRUE
  ), error = identity)
  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), "Source files overlap", fixed = TRUE)
  for (path in paths) expect_match(conditionMessage(error), basename(path), fixed = TRUE)
  .make_buffer_grid_nc(paths[1L])
  expect_error(extract_to_buffer(
    data.frame(lon = 0, lat = 0), c(paths[1L], paths[1L]),
    var = "temp", horizontal_radius_km = 1, quiet = TRUE
  ), "Source files overlap", fixed = TRUE)
})


test_that("point time is resolved globally before a vertical mean", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  bounds <- matrix(c(0, 10, 10, 30), 2L)
  .make_buffer_source_nc(
    paths[1L], matrix(c(10, 20, 50, 60), 2L), depth = c(5, 20),
    time = c(0, 48), depth_bounds = bounds
  )
  .make_buffer_source_nc(
    paths[2L], matrix(c(30, 40), 2L), depth = c(5, 20),
    time = 24, depth_bounds = bounds
  )
  data <- data.frame(lon = 0, lat = 40, depth = 15, date = as.Date("2020-01-02"))
  for (files in list(paths, rev(paths))) {
    out <- extract_to_buffer(data, files, var = "temp", depth_radius_m = 10, quiet = TRUE)
    expect_equal(out$temp, 37.5)
    expect_length(out$temp_match_info[[1L]]$sources_used, 1)
  }
})


test_that("point depth is selected across source layers without merging means", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  .make_buffer_source_nc(paths[1L], matrix(10), depth = 5, time = 12,
                         time_bounds = matrix(c(0, 24), 2L))
  .make_buffer_source_nc(paths[2L], matrix(20), depth = 20, time = 12,
                         time_bounds = matrix(c(0, 24), 2L))
  data <- data.frame(lon = 0, lat = 40, depth = 18, date = "2020-01-01T12:00Z")
  for (files in list(paths, rev(paths))) {
    out <- extract_to_buffer(data, files, var = "temp", time_before = 12,
                             time_after = 12, time_unit = "hours", quiet = TRUE)
    expect_equal(out$temp, 20)
    expect_equal(out$temp_match_info[[1L]]$depth_range_m, c(20, 20))
  }
})


test_that("index layer windows retain physical weights and reject grid edges", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(30, 20, 10), 3L), depth = c(35, 15, 5),
    depth_bounds = matrix(c(20, 50, 10, 20, 0, 10), 2L)
  )
  data <- data.frame(lon = 0, lat = 40, depth = c(15, 5))
  for (weighting in c("weighted", "equal")) {
    out <- extract_to_buffer(
      data, path, depth_layer_radius = 1, weighting = weighting,
      min_valid_fraction = 0, quiet = TRUE
    )
    expect_equal(out$temp, c(if (weighting == "weighted") 24 else 20, NA_real_))
    expect_identical(out$temp_match_info[[2L]]$status, "incomplete_index_window")
  }
})


test_that("index cell windows average a full neighbourhood and retain edge NAs", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(
    path, matrix(c(0, 10, 20), 3L), lon = c(-0.01, 0.01, 0.03), lat = 0,
    lon_bounds = matrix(c(-0.02, 0, 0, 0.02, 0.02, 0.04), 2L),
    lat_bounds = matrix(c(-0.01, 0.01), 2L)
  )
  out <- extract_to_buffer(
    data.frame(lon = c(0.01, -0.01), lat = 0), path,
    lon_cell_radius = 1, lat_cell_radius = 0, quiet = TRUE
  )
  expect_equal(out$temp, c(10, NA_real_), tolerance = 1e-5)
  expect_identical(out$temp_match_info[[2L]]$status, "incomplete_index_window")
})


test_that("actual antimeridian cells intersect the requested circle", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(
    path, lon = c(179.99, 180.01),
    lon_bounds = matrix(c(179.98, 180, 180, 180.02), 2L)
  )
  out <- extract_to_buffer(data.frame(lon = c(180, -180), lat = 0), path,
                           horizontal_radius_km = 1, quiet = TRUE)
  expect_equal(out$temp, c(15, 15), tolerance = 1e-5)
  expect_equal(out$temp_match_info[[1L]]$valid_fraction, 1)
})


test_that("multiple variables have distinct outputs and collision-safe diagnostics", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  .make_buffer_grid_nc(paths[1L])
  .make_buffer_grid_nc(paths[2L], matrix(c(2, 4, 6, 8), 2L), var = "oxygen")
  data <- data.frame(lon = 0, lat = 0)
  out <- extract_to_buffer(
    data, paths, var = c("temp", "oxygen"), horizontal_radius_km = 1,
    output_col = c(temp = "temperature", oxygen = "oxygen_mean"), quiet = TRUE
  )
  expect_named(out, c("lon", "lat", "temperature", "temperature_match_info",
                       "oxygen_mean", "oxygen_mean_match_info"))
  expect_equal(out$temperature, 15, tolerance = 1e-6)
  expect_equal(out$oxygen_mean, 5, tolerance = 1e-6)
  data$temp <- 999
  expect_error(extract_to_buffer(data, paths[1L], horizontal_radius_km = 1),
                "overwrite", fixed = TRUE)
})


test_that("coverage reports explain partial results and quiet suppresses only reports", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(path, matrix(10, 1L, 2L), lon = 0.01,
                       lon_bounds = matrix(c(0, 0.02), 2L))
  data <- data.frame(lon = 0, lat = 0)
  expect_message(out <- extract_to_buffer(
    data, path, horizontal_radius_km = 1, min_valid_fraction = 0.5
  ), "partial coverage", fixed = TRUE)
  expect_equal(out$temp, 10)
  expect_equal(out$temp_match_info[[1L]]$valid_fraction, 0.5, tolerance = 1e-6)
  expect_no_message(rejected <- extract_to_buffer(
    data, path, horizontal_radius_km = 1, quiet = TRUE
  ))
  expect_true(is.na(rejected$temp))

  .make_buffer_source_nc(path, matrix(10), time = 12,
                         time_bounds = matrix(c(0, 24), 2L))
  expect_warning(timed <- extract_to_buffer(
    data.frame(lon = 0, lat = 40, date = "2020-01-01 12:00"), path,
    time_before = 12, time_after = 12, time_unit = "hours", quiet = TRUE
  ), class = "ocean3d_timezone_warning")
  expect_equal(timed$temp, 10)
})


test_that("spatial bounds inference is explicit and statistics are validated", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(path, spatial_bounds = FALSE)
  data <- data.frame(lon = 0, lat = 0)
  expect_error(extract_to_buffer(data, path, horizontal_radius_km = 1),
                "infer_bounds = TRUE", fixed = TRUE)
  out <- extract_to_buffer(data, path, horizontal_radius_km = 1,
                           infer_bounds = TRUE, quiet = TRUE)
  expect_equal(out$temp, 15, tolerance = 1e-6)
  expect_error(extract_to_buffer(data, path, horizontal_radius_km = 1,
                                 statistic = "sum"),
                "`statistic` must be one of", fixed = TRUE)
  expect_error(extract_to_buffer(data, path, horizontal_radius_km = 1,
                                 lon_cell_radius = 1), "not both", fixed = TRUE)
})


test_that("public extraction retains sf geometries and caller-owned connections", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  .make_buffer_grid_nc(path)
  handle <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(handle), add = TRUE, after = FALSE)
  data <- sf::st_as_sf(data.frame(id = "a", x = 0, y = 0), coords = c("x", "y"), crs = 4326)
  out <- extract_to_buffer(data, handle, horizontal_radius_km = 1, quiet = TRUE)
  expect_s3_class(out, "sf")
  expect_identical(sf::st_geometry(out), sf::st_geometry(data))
  expect_equal(out$temp, 15, tolerance = 1e-6)
  expect_equal(as.vector(ncdf4::ncvar_get(handle, "temp")), c(0, 20, 10, 30))
})


test_that("interval ordering does not change the mean", {
  result <- .buffer_interval_mean(
    values = c(20, 10),
    lower = c(10, 0), upper = c(30, 10),
    requested = c(5, 25)
  )
  expect_equal(result$value, 17.5)
  expect_identical(result$status, "complete")
})


test_that("overlapping support and invalid intervals are rejected", {
  expect_error(
    .buffer_interval_mean(
      values = c(10, 20),
      lower = c(0, 5), upper = c(10, 15),
      requested = c(0, 15)
    ),
    "Selected intervals overlap", fixed = TRUE
  )
  expect_error(
    .buffer_interval_mean(10, 0, 10, requested = c(5, 5)),
    "requested", fixed = TRUE
  )
  expect_error(
    .buffer_interval_mean(c(10, 20), 0, 10, requested = c(0, 10)),
    "matching lengths", fixed = TRUE
  )
  expect_error(
    .buffer_interval_mean(10, 10, 0, requested = c(0, 10)),
    "interval limits", fixed = TRUE
  )
})


# Axis bounds ------------------------------------------------------------------

test_that("explicit depth bounds feed the known weighted and equal means", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_axis_nc(
    path, c(5, 20), matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  support <- .read_buffer_axis_nc(path)
  expect_equal(support$centres, c(5, 20))
  expect_equal(support$lower, c(0, 10))
  expect_equal(support$upper, c(10, 30))
  expect_identical(support$basis, "explicit_bounds")
  expect_identical(support$bounds_name, "axis_bounds")
  expect_identical(support$units, "m")
  expect_false(any(support$clipped_at_surface))

  weighted <- .buffer_interval_mean(
    c(10, 20), support$lower, support$upper, c(5, 25)
  )
  equal <- .buffer_interval_mean(
    c(10, 20), support$lower, support$upper, c(5, 25), weighting = "equal"
  )
  expect_equal(weighted$value, 17.5)
  expect_equal(equal$value, 15)
  expect_equal(weighted$valid_fraction, 1)
})


test_that("reversed positive-up depths and shared bounds convert to metres", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_axis_nc(
    path, c(-2000, -500),
    matrix(c(-3000, -1000, -1000, 0), nrow = 2L),
    units = "cm", attributes = list(positive = "up")
  )
  support <- .read_buffer_axis_nc(path)
  expect_equal(support$centres, c(20, 5))
  expect_equal(support$lower, c(10, 0))
  expect_equal(support$upper, c(30, 10))
  expect_identical(support$source_index, 1:2)
  expect_identical(support$direction_basis, "metadata")
  result <- .buffer_interval_mean(
    c(20, 10), support$lower, support$upper, c(5, 25)
  )
  expect_equal(result$value, 17.5)
})


test_that("both declared bounds dimension orders are supported", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  for (i in seq_along(paths)) {
    .make_buffer_axis_nc(
      paths[i], c(5, 20, 40),
      matrix(c(0, 10, 10, 30, 30, 50), nrow = 2L),
      coordinate_first = i == 2L
    )
    support <- .read_buffer_axis_nc(paths[i])
    expect_equal(support$lower, c(0, 10, 30))
    expect_equal(support$upper, c(10, 30, 50))
  }
})


test_that("missing bounds require explicit permission to infer", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_axis_nc(path, c(0, 10, 30))
  expect_error(.read_buffer_axis_nc(path), "infer_bounds = TRUE", fixed = TRUE)
  support <- .read_buffer_axis_nc(path, infer_bounds = TRUE)
  expect_equal(support$lower, c(0, 5, 20))
  expect_equal(support$upper, c(5, 20, 40))
  expect_identical(support$basis, "inferred_midpoints")
  expect_identical(support$clipped_at_surface, c(TRUE, FALSE, FALSE))

  # The request is not clipped at the surface: only half has support.
  result <- .buffer_interval_mean(
    c(10, 20, 30), support$lower, support$upper, c(-5, 5)
  )
  expect_equal(result$coverage_fraction, 0.5)
  expect_equal(result$valid_fraction, 0.5)
  expect_true(is.na(result$value))
})


test_that("a singleton needs explicit bounds and keeps its extent", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  .make_buffer_axis_nc(paths[1], 10)
  expect_error(
    .read_buffer_axis_nc(paths[1], infer_bounds = TRUE),
    "single coordinate", fixed = TRUE
  )
  .make_buffer_axis_nc(paths[2], 10, matrix(c(5, 15), nrow = 2L))
  support <- .read_buffer_axis_nc(paths[2])
  expect_equal(support$lower, 5)
  expect_equal(support$upper, 15)
})


test_that("time bounds decode to UTC and allow start-labelled intervals", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_axis_nc(
    path, c(0, 24), matrix(c(0, 24, 24, 48), nrow = 2L),
    axis = "time", attributes = list(calendar = "gregorian")
  )
  support <- .read_buffer_axis_nc(path, "time", axis = "time")
  expect_equal(support$lower, as.numeric(as.POSIXct(
    c("2020-01-01", "2020-01-02"), tz = "UTC"
  )))
  expect_equal(support$upper, as.numeric(as.POSIXct(
    c("2020-01-02", "2020-01-03"), tz = "UTC"
  )))
  expect_identical(support$calendar, "gregorian")
  expect_identical(support$basis, "explicit_bounds")
})


test_that("irregular time inference uses neighbouring midpoints", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_axis_nc(path, c(0, 2, 8), axis = "time")
  support <- .read_buffer_axis_nc(
    path, "time", axis = "time", infer_bounds = TRUE
  )
  origin <- as.numeric(as.POSIXct("2020-01-01", tz = "UTC"))
  expect_equal(support$lower, origin + c(-3600, 3600, 18000))
  expect_equal(support$upper, origin + c(3600, 18000, 39600))
  expect_identical(support$basis, "inferred_midpoints")
})


test_that("invalid declared bounds never fall back to inference", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  cases <- list(
    list(bounds = matrix(c(0, 10, 5, 30), nrow = 2L), pattern = "overlap"),
    list(bounds = matrix(c(0, 4, 10, 30), nrow = 2L), pattern = "centres"),
    list(bounds = matrix(c(5, 5, 10, 30), nrow = 2L),
         pattern = "positive widths"),
    list(bounds = matrix(c(NA, 10, 10, 30), nrow = 2L), pattern = "finite"),
    list(bounds = matrix(0:5, nrow = 3L), pattern = "two finite endpoints")
  )
  for (case in cases) {
    .make_buffer_axis_nc(path, c(5, 20), case$bounds)
    expect_error(
      .read_buffer_axis_nc(path, infer_bounds = TRUE),
      case$pattern, fixed = TRUE
    )
  }
  .make_buffer_axis_nc(path, c(5, 20), bounds_reference = "missing_bounds")
  expect_error(
    .read_buffer_axis_nc(path, infer_bounds = TRUE),
    "missing bounds variable", fixed = TRUE
  )
})


test_that("incompatible units and directions are reported explicitly", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  bounds <- matrix(c(0, 10, 10, 30), nrow = 2L)
  .make_buffer_axis_nc(path, c(5, 20), bounds, bounds_units = "cm")
  expect_error(.read_buffer_axis_nc(path), "different units", fixed = TRUE)
  .make_buffer_axis_nc(
    path, c(5, 20), bounds, attributes = list(positive = "down"),
    bounds_attributes = list(positive = "up")
  )
  expect_error(
    .read_buffer_axis_nc(path), "conflicting directions", fixed = TRUE
  )
})


test_that("unsupported calendars and climatologies are not inferred", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  for (attributes in list(
    list(calendar = "360_day"), list(climatology = "climatology_bounds")
  )) {
    .make_buffer_axis_nc(path, c(0, 24), axis = "time", attributes = attributes)
    expect_error(
      .read_buffer_axis_nc(path, "time", axis = "time", infer_bounds = TRUE),
      if (is.null(attributes$calendar)) {
        "Climatological"
      } else {
        "unsupported calendar"
      },
      fixed = TRUE
    )
  }
  .make_buffer_axis_nc(
    path, c(0, 24), matrix(c(0, 24, 24, 48), nrow = 2L),
    axis = "time", bounds_attributes = list(calendar = "360_day")
  )
  expect_error(
    .read_buffer_axis_nc(path, "time", axis = "time"),
    "different calendars", fixed = TRUE
  )
})


test_that("duplicates and unordered coordinates are rejected", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  for (centres in list(c(10, 10), c(10, 30, 20))) {
    .make_buffer_axis_nc(path, centres)
    expect_error(.read_buffer_axis_nc(path, infer_bounds = TRUE))
  }
})


test_that("inferred depth intervals retain reversed source order", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_axis_nc(path, c(30, 10, 0))
  support <- .read_buffer_axis_nc(path, infer_bounds = TRUE)
  expect_equal(support$lower, c(20, 5, 0))
  expect_equal(support$upper, c(40, 20, 5))
  expect_identical(support$source_index, 1:3)
})


test_that("kilometres convert to metres and transformed axes are rejected", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_axis_nc(
    path, c(0.005, 0.02),
    matrix(c(0, 0.01, 0.01, 0.03), nrow = 2L), units = "km"
  )
  support <- .read_buffer_axis_nc(path)
  expect_equal(support$centres, c(5, 20))
  expect_equal(support$lower, c(0, 10))
  expect_equal(support$upper, c(10, 30))
  for (attributes in list(
    list(standard_name = "sea_water_pressure"),
    list(formula_terms = "sigma: sigma eta: eta depth: depth")
  )) {
    .make_buffer_axis_nc(path, c(5, 20), attributes = attributes)
    expect_error(
      .read_buffer_axis_nc(path, infer_bounds = TRUE),
      "physical transformation", fixed = TRUE
    )
  }
})


test_that("reading an axis leaves a caller-owned connection usable", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  .make_buffer_axis_nc(
    path, c(5, 20), matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE, after = FALSE)
  support <- .buffer_read_axis(nc, "depth")
  expect_equal(support$centres, c(5, 20))
  values <- ncdf4::ncvar_get(nc, "value")
  expect_identical(dim(values), 2L)
  expect_equal(as.vector(values), c(1, 2))
})


# Source extraction ------------------------------------------------------------

test_that("vertical extraction reads known cells and preserves observations", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(10, 20), nrow = 2L), depth = c(5, 20),
    depth_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(
    id = c("a", "b"), lon = c(0, 1), lat = c(40, 41),
    depth = c(15, 15), date = "unused"
  )
  original <- data
  options <- .validate_buffer_options(depth_radius_m = 10)
  result <- .read_buffer_source_nc(path, data, options, id_col = "id")
  expect_identical(result$observations, original)
  expect_identical(data, original)
  expect_equal(result$values, c(17.5, 317.5))
  expect_identical(result$diagnostics[[1L]]$status, "complete")
  expect_identical(result$diagnostics[[1L]]$depth_indices, 1:2)
  expect_equal(result$diagnostics[[1L]]$n_valid, 2)
  expect_identical(result$metadata$depth_support$basis, "explicit_bounds")
  expect_null(result$prepared_time)

  equal <- .read_buffer_source_nc(path, data, .validate_buffer_options(
    depth_radius_m = 10, weighting = "equal"
  ))
  expect_equal(equal$values, c(15, 315))
})


test_that("vertical extraction retains reversed layer alignment", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(20, 10), nrow = 2L), depth = c(20, 5),
    depth_bounds = matrix(c(10, 30, 0, 10), nrow = 2L)
  )
  result <- .read_buffer_source_nc(
    path, data.frame(lon = 0, lat = 40, depth = 15),
    .validate_buffer_options(depth_radius_m = 10)
  )
  expect_equal(result$values, 17.5)
  expect_equal(result$diagnostics[[1L]]$selected_depth_m, c(20, 5))
})


test_that("a vertical buffer uses point time matching without time bounds", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(10, 20, 30, 50), nrow = 2L),
    depth = c(5, 20), time = c(0, 24),
    depth_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(lon = 0, lat = 40, depth = 15,
                     date = as.Date("2020-01-02"))
  options <- .validate_buffer_options(depth_radius_m = 10)
  result <- .read_buffer_source_nc(path, data, options, time_match = "exact")
  expect_equal(result$values, 45)
  expect_identical(result$diagnostics[[1L]]$time_indices, 2L)
  expect_identical(result$observations$date, data$date)

  data$date <- as.Date("2020-01-03")
  rejected <- .read_buffer_source_nc(
    path, data, options, time_match = "exact"
  )
  expect_true(is.na(rejected$values))
  expect_identical(rejected$diagnostics[[1L]]$status, "time_rejected")
  limited <- .read_buffer_source_nc(path, data, options, max_time_diff = 0.5)
  expect_true(is.na(limited$values))
  expect_identical(limited$diagnostics[[1L]]$status, "time_rejected")
})


test_that("inferred support retains full-request coverage", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(10, 20), nrow = 2L), depth = c(5, 20)
  )
  data <- data.frame(lon = 0, lat = 40, depth = 15)
  expect_error(
    .read_buffer_source_nc(
      path, data, .validate_buffer_options(depth_radius_m = 10)
    ),
    "infer_bounds = TRUE", fixed = TRUE
  )
  result <- .read_buffer_source_nc(path, data, .validate_buffer_options(
    depth_radius_m = 10, infer_bounds = TRUE
  ))
  expect_equal(result$values, 16.25)
  expect_identical(result$metadata$depth_support$basis, "inferred_midpoints")

  data$depth <- 0
  partial <- .read_buffer_source_nc(path, data, .validate_buffer_options(
    depth_radius_m = 5, infer_bounds = TRUE
  ))
  expect_equal(partial$diagnostics[[1L]]$coverage_fraction, 0.5)
  expect_true(is.na(partial$values))
  expect_identical(partial$diagnostics[[1L]]$status, "insufficient_coverage")
})


test_that("time extraction weights only the intersecting record durations", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(10, 20), nrow = 1L), time = c(5, 20),
    time_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(lon = 0, lat = 40, date = "2020-01-01T15:00Z")
  for (weighting in c("weighted", "equal")) {
    result <- .read_buffer_source_nc(path, data, .validate_buffer_options(
      time_before = 10, time_after = 10, time_unit = "hours",
      method = "2d", weighting = weighting
    ))
    expect_equal(result$values, if (weighting == "weighted") 17.5 else 15)
    expect_equal(result$diagnostics[[1L]]$coverage_fraction, 1)
    expect_identical(result$diagnostics[[1L]]$time_indices, 1:2)
    expect_identical(result$observations, data)
    expect_length(result$prepared_time$unzoned_rows, 0L)
  }
})


test_that("equivalent local and UTC times produce the same buffer", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(10, 20), nrow = 1L), time = c(5, 20),
    time_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  options <- .validate_buffer_options(
    time_before = 10, time_after = 10, time_unit = "hours", method = "2d"
  )
  utc <- as.POSIXct("2020-01-01 15:00:00", tz = "UTC")
  madrid <- as.POSIXct("2020-01-01 16:00:00", tz = "Europe/Madrid")
  for (date in list(utc, madrid, as.POSIXlt(madrid))) {
    data <- data.frame(lon = 0, lat = 40)
    data$date <- date
    result <- .read_buffer_source_nc(path, data, options)
    expect_equal(result$values, 17.5)
    expect_equal(
      format(result$windows$time_lower_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
      "2020-01-01 05:00"
    )
    expect_equal(
      format(result$windows$time_upper_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
      "2020-01-02 01:00"
    )
  }
})


test_that("time buffers select surface, bottom and nearest physical layers", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(100, 10, 200, 20), nrow = 2L),
    depth = c(20, 5), time = c(5, 20),
    time_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(lon = 0, lat = 40, depth = 12.5,
                     date = "2020-01-01T15:00Z")
  for (method in c("surface", "bottom", "nearest")) {
    result <- .read_buffer_source_nc(path, data, .validate_buffer_options(
      time_before = 10, time_after = 10, time_unit = "hours", method = method
    ))
    expect_equal(result$values, if (method == "bottom") 175 else 17.5)
    expect_equal(
      result$diagnostics[[1L]]$selected_depth_m,
      if (method == "bottom") c(20, 20) else c(5, 5)
    )
  }
})


test_that("a missing surface is not replaced while nearest uses valid layers", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(20, NA_real_, 30, 10), nrow = 2L),
    depth = c(20, 5), time = c(5, 20),
    time_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(lon = 0, lat = 40, depth = 5,
                     date = "2020-01-01T15:00Z")
  nearest <- .read_buffer_source_nc(path, data, .validate_buffer_options(
    time_before = 10, time_after = 10, time_unit = "hours"
  ))
  expect_equal(nearest$values, 12.5)
  expect_equal(nearest$diagnostics[[1L]]$selected_depth_m, c(20, 5))

  surface <- .read_buffer_source_nc(path, data, .validate_buffer_options(
    time_before = 10, time_after = 10, time_unit = "hours", method = "surface"
  ))
  expect_true(is.na(surface$values))
  expect_equal(surface$diagnostics[[1L]]$coverage_fraction, 1)
  expect_equal(surface$diagnostics[[1L]]$valid_fraction, 0.75)
  expect_identical(surface$diagnostics[[1L]]$status, "insufficient_coverage")

  partial <- .read_buffer_source_nc(path, data, .validate_buffer_options(
    time_before = 10, time_after = 10, time_unit = "hours",
    method = "surface", min_valid_fraction = 0.75
  ))
  expect_equal(partial$values, 10)
  expect_identical(partial$diagnostics[[1L]]$time_indices, 2L)
  expect_identical(partial$diagnostics[[1L]]$status, "partial")
})


test_that("empty buffers and missing coordinates preserve all rows", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(10, 20), nrow = 2L), depth = c(5, 20),
    depth_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(
    id = c("ok", "empty", "missing", "outside"),
    lon = c(0, 0, 0, 100), lat = 40, depth = c(15, 60, NA_real_, 15)
  )
  result <- .read_buffer_source_nc(
    path, data, .validate_buffer_options(depth_radius_m = 10), id_col = "id"
  )
  expect_identical(result$observations, data)
  expect_equal(result$values, c(17.5, NA_real_, NA_real_, NA_real_))
  expect_identical(
    vapply(result$diagnostics, function(x) x$status, character(1)),
    c("complete", "empty", "missing_depth", "spatial_rejected")
  )
  expect_equal(result$diagnostics[[2L]]$coverage_fraction, 0)
  expect_true(is.na(result$diagnostics[[3L]]$coverage_fraction))
})


test_that("a temporal singleton with bounds can fill a requested window", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(10), time = 12, time_bounds = matrix(c(0, 24), nrow = 2L)
  )
  options <- .validate_buffer_options(
    time_before = 12, time_after = 12, time_unit = "hours", method = "2d"
  )
  data <- data.frame(lon = 0, lat = 40, date = "2020-01-01T12:00Z")
  result <- .read_buffer_source_nc(path, data, options)
  expect_equal(result$values, 10)
  expect_equal(result$diagnostics[[1L]]$valid_fraction, 1)
  expect_equal(result$diagnostics[[1L]]$n_intersecting, 1)

  data$date <- NA_character_
  missing <- .read_buffer_source_nc(path, data, options)
  expect_identical(missing$observations, data)
  expect_true(is.na(missing$values))
  expect_identical(missing$diagnostics[[1L]]$status, "missing_time")
  data$date <- "2020-01-10T12:00Z"
  empty <- .read_buffer_source_nc(path, data, options)
  expect_true(is.na(empty$values))
  expect_identical(empty$diagnostics[[1L]]$status, "empty")
})


test_that("missing axes and horizontal buffers fail explicitly", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(path, matrix(10))
  data <- data.frame(lon = 0, lat = 40, depth = 15,
                     date = "2020-01-01T12:00Z")
  expect_error(
    .read_buffer_source_nc(
      path, data, .validate_buffer_options(depth_radius_m = 10)
    ),
    "requires a depth dimension", fixed = TRUE
  )
  expect_error(
    .read_buffer_source_nc(path, data, .validate_buffer_options(
      time_before = 1, method = "2d"
    )),
    "requires a time dimension", fixed = TRUE
  )
  expect_error(
    .read_buffer_source_nc(path, data, .validate_buffer_options(
      depth_radius_m = 10, time_before = 1
    )),
    "requires a depth dimension", fixed = TRUE
  )
  expect_error(
    .read_buffer_source_nc(path, data, .validate_buffer_options(
      horizontal_radius_km = 5, time_before = 1
    )),
    "not implemented yet", fixed = TRUE
  )
})


test_that("custom columns and missing observation dates keep their meaning", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(10), time = 12, time_bounds = matrix(c(0, 24), nrow = 2L)
  )
  options <- .validate_buffer_options(
    time_before = 12, time_after = 12, time_unit = "hours", method = "2d"
  )
  data <- data.frame(
    x = c(0, 0), y = c(40, 40), sample = c("a", "b"),
    observed_at = c("2020-01-01T12:00Z", NA_character_)
  )
  result <- .read_buffer_source_nc(
    path, data, options, lon_col = "x", lat_col = "y",
    date_col = "observed_at", id_col = "sample"
  )
  expect_equal(result$values, c(10, NA_real_))
  expect_identical(result$observations, data)
  data$observed_at <- NULL
  expect_error(
    .read_buffer_source_nc(
      path, data, options, lon_col = "x", lat_col = "y",
      date_col = "observed_at", id_col = "sample"
    ),
    "Column `observed_at` is missing", fixed = TRUE
  )
})


test_that("reading source values keeps caller-owned connections open", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  .make_buffer_source_nc(
    path, matrix(c(10, 20), nrow = 2L), depth = c(5, 20),
    depth_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE, after = FALSE)
  result <- .read_buffer_source_nc(
    nc, data.frame(lon = 0, lat = 40, depth = 15),
    .validate_buffer_options(depth_radius_m = 10)
  )
  expect_equal(result$values, 17.5)
  expect_equal(as.vector(ncdf4::ncvar_get(nc, "temp"))[1L], 10)
})


# Combined depth and time buffers ----------------------------------------------

test_that("combined extraction uses joint depth and time weights", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(0, 10, 20, 40), nrow = 2L),
    depth = c(5, 20), time = c(5, 20),
    depth_bounds = matrix(c(0, 10, 10, 30), nrow = 2L),
    time_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(
    id = "a", lon = 0, lat = 40, depth = 15,
    date = "2020-01-01T15:00Z"
  )

  for (weighting in c("weighted", "equal")) {
    result <- .read_buffer_source_nc(path, data, .validate_buffer_options(
      depth_radius_m = 10,
      time_before = 10, time_after = 10, time_unit = "hours",
      weighting = weighting
    ))
    expect_equal(result$values, if (weighting == "weighted") 28.125 else 17.5)
    expect_identical(result$observations, data)
    expect_equal(result$diagnostics[[1L]]$coverage_fraction, 1)
    expect_equal(result$diagnostics[[1L]]$valid_fraction, 1)
    expect_equal(result$diagnostics[[1L]]$n_intersecting, 4)
    expect_identical(
      result$diagnostics[[1L]]$depth_indices, c(1L, 2L, 1L, 2L)
    )
    expect_identical(
      result$diagnostics[[1L]]$time_indices, c(1L, 1L, 2L, 2L)
    )
  }
})


test_that("missing joint samples do not produce a mean of separate means", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(0, 10, 20, NA_real_), nrow = 2L),
    depth = c(5, 20), time = c(5, 20),
    depth_bounds = matrix(c(0, 10, 10, 30), nrow = 2L),
    time_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(
    lon = 0, lat = 40, depth = 15, date = "2020-01-01T15:00Z"
  )

  # The missing sample occupies 9/16 of the requested support.
  # Valid samples have relative weights 1, 3 and 3.
  for (weighting in c("weighted", "equal")) {
    for (threshold in c(0.4374, 0.4375, 0.4376)) {
      result <- .read_buffer_source_nc(path, data, .validate_buffer_options(
        depth_radius_m = 10,
        time_before = 10, time_after = 10, time_unit = "hours",
        weighting = weighting, min_valid_fraction = threshold
      ))
      expect_equal(result$diagnostics[[1L]]$coverage_fraction, 1)
      expect_equal(result$diagnostics[[1L]]$valid_fraction, 0.4375)
      expect_equal(result$diagnostics[[1L]]$n_valid, 3)

      if (threshold <= 0.4375) {
        expect_equal(result$values, if (weighting == "weighted") 90 / 7 else 10)
        expect_identical(result$diagnostics[[1L]]$status, "partial")
      } else {
        expect_true(is.na(result$values))
        expect_identical(
          result$diagnostics[[1L]]$status, "insufficient_coverage"
        )
      }
    }
  }
})


test_that("partial coverage on both axes uses the complete requested buffer", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(0, 10, 20, 40), nrow = 2L),
    depth = c(5, 20), time = c(5, 20),
    depth_bounds = matrix(c(0, 10, 10, 30), nrow = 2L),
    time_bounds = matrix(c(0, 10, 10, 30), nrow = 2L)
  )
  data <- data.frame(
    lon = 0, lat = 40, depth = 10, date = "2020-01-01T10:00Z"
  )

  # Each axis covers 30 of the requested 40 units: joint coverage is 9/16.
  for (threshold in c(0.8, 0.5625)) {
    result <- .read_buffer_source_nc(path, data, .validate_buffer_options(
      depth_radius_m = 20,
      time_before = 20, time_after = 20, time_unit = "hours",
      min_valid_fraction = threshold
    ))
    expect_equal(result$windows$depth_lower_m, -10)
    expect_equal(result$diagnostics[[1L]]$coverage_fraction, 0.5625)
    expect_equal(result$diagnostics[[1L]]$valid_fraction, 0.5625)

    if (threshold == 0.8) {
      expect_true(is.na(result$values))
    } else {
      expect_equal(result$values, 220 / 9)
      expect_identical(result$diagnostics[[1L]]$status, "partial")
    }
  }
})


test_that("combined extraction preserves reversed axes and rejected rows", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(
    path, matrix(c(40, 20, 10, 0), nrow = 2L),
    depth = c(20, 5), time = c(20, 5),
    depth_bounds = matrix(c(10, 30, 0, 10), nrow = 2L),
    time_bounds = matrix(c(10, 30, 0, 10), nrow = 2L)
  )
  data <- data.frame(
    lon = 0, lat = 40, depth = c(15, 40, NA_real_, 15),
    date = c(rep("2020-01-01T15:00Z", 3), NA_character_)
  )
  result <- .read_buffer_source_nc(path, data, .validate_buffer_options(
    depth_radius_m = 10,
    time_before = 10, time_after = 10, time_unit = "hours"
  ))
  expect_identical(result$observations, data)
  expect_equal(result$values, c(28.125, NA_real_, NA_real_, NA_real_))
  expect_identical(
    vapply(result$diagnostics, function(x) x$status, character(1)),
    c("complete", "empty", "missing_depth", "missing_time")
  )
  expect_equal(result$diagnostics[[2L]]$coverage_fraction, 0)
})


test_that("buffer statistics use known equal and unequal weight results", {
  values <- c(10, 20)
  statistics <- c("mean", "sd", "se", "min", "max")
  expected <- c(17.5, sqrt(50), sqrt(31.25), 10, 20)
  for (i in seq_along(statistics)) {
    expect_equal(.buffer_statistic(values, c(1, 3), statistics[i]), expected[i])
    expect_equal(.buffer_statistic(values, c(10, 30), statistics[i]), expected[i])
  }
  expect_equal(.buffer_statistic(c(2, 4, 6), c(1, 1, 1), "sd"), 2)
  expect_equal(.buffer_statistic(c(2, 4, 6), c(1, 1, 1), "se"), 2 / sqrt(3))
  for (statistic in c("sd", "se")) {
    expect_equal(.buffer_statistic(c(7, 7, 7), c(1, 2, 3), statistic), 0)
    expect_true(is.na(.buffer_statistic(7, 1, statistic)))
  }
  for (value in list(NULL, NA_character_, 1, c("min", "max"), "SD", "s")) {
    expect_error(.validate_buffer_options(depth_radius_m = 1, statistic = value),
                 "`statistic` must be one of", fixed = TRUE)
  }
})


test_that("all public statistics combine files before aggregation in either order", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths))
  .make_buffer_source_nc(paths[1L], matrix(10), time = 5,
                         time_bounds = matrix(c(0, 10), 2L))
  .make_buffer_source_nc(paths[2L], matrix(20), time = 20,
                         time_bounds = matrix(c(10, 30), 2L))
  data <- data.frame(lon = 0, lat = 40, date = "2020-01-01T15:00Z")
  statistics <- c("mean", "sd", "se", "min", "max")
  for (files in list(paths, rev(paths))) for (weighting in c("weighted", "equal")) {
    expected <- if (weighting == "weighted") {
      c(17.5, sqrt(50), sqrt(31.25), 10, 20)
    } else c(15, sqrt(50), 5, 10, 20)
    for (i in seq_along(statistics)) {
      out <- extract_to_buffer(
        data, files, var = "temp", time_before = 10, time_after = 10,
        time_unit = "hours", statistic = statistics[i], weighting = weighting,
        output_col = "summary", quiet = TRUE
      )
      expect_named(out, c(names(data), "summary", "summary_match_info"))
      expect_equal(out$summary, expected[i])
      info <- out$summary_match_info[[1L]]
      expect_identical(info$statistic, statistics[i])
      expect_equal(info$effective_n, if (weighting == "weighted") 1.6 else 2)
      expect_identical(info$status, "complete")
      expect_length(info$sources_used, 2)
    }
  }
})


test_that("all statistics respect missing values and physical coverage", {
  samples <- data.frame(
    value = c(10, 20, NA_real_), area_fraction = c(0.2, 0.6, 0.2),
    depth_fraction = 1, time_fraction = 1
  )
  statistics <- c("mean", "sd", "se", "min", "max")
  expected <- c(17.5, sqrt(50), sqrt(31.25), 10, 20)
  for (i in seq_along(statistics)) {
    out <- .buffer_summarise_samples(samples, "weighted", 0.8, statistics[i])
    expect_equal(out$value, expected[i])
    expect_equal(out$valid_fraction, 0.8)
    expect_equal(out$effective_n, 1.6)
    expect_identical(out$status, "partial")
    rejected <- .buffer_summarise_samples(samples, "weighted", 0.9, statistics[i])
    expect_true(is.na(rejected$value))
    expect_identical(rejected$status, "insufficient_coverage")
    missing <- samples
    missing$value[] <- NA_real_
    rejected <- .buffer_summarise_samples(missing, "weighted", 0, statistics[i])
    expect_true(is.na(rejected$value))
    expect_identical(rejected$status, "all_missing")
    empty <- .buffer_summarise_samples(samples[FALSE, ], "equal", 0, statistics[i])
    expect_true(is.na(empty$value))
    expect_identical(empty$status, "empty")
  }
})


test_that("public SD and SE explain an insufficient number of valid samples", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_source_nc(path, matrix(10), time = 5,
                         time_bounds = matrix(c(0, 10), 2L))
  data <- data.frame(lon = 0, lat = 40, date = "2020-01-01T05:00Z")
  for (statistic in c("sd", "se")) {
    expect_message(out <- extract_to_buffer(
      data, path, time_before = 5, time_after = 5, time_unit = "hours",
      statistic = statistic
    ), "at least two valid samples", fixed = TRUE)
    expect_true(is.na(out$temp))
    expect_identical(out$temp_match_info[[1L]]$status,
                     "insufficient_variation_samples")
    expect_equal(out$temp_match_info[[1L]]$effective_n, 1)
    expect_message(extract_to_buffer(
      data, path, time_before = 5, time_after = 5, time_unit = "hours",
      statistic = statistic, quiet = TRUE
    ), NA)
  }
  for (statistic in c("mean", "min", "max")) {
    out <- extract_to_buffer(data, path, time_before = 5, time_after = 5,
                              time_unit = "hours", statistic = statistic)
    expect_equal(out$temp, 10)
    expect_identical(out$temp_match_info[[1L]]$status, "complete")
  }
})


test_that("automatic detection ignores dimensions without coordinate variables", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  .make_buffer_grid_nc(path)
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_equal(nc$dim$nv$dimvarid$id, -1)
  expect_output(variable <- .buffer_detect_variable(nc), NA)
  expect_identical(variable, "temp")
})


test_that("joint aggregation retains singleton and all-missing cases", {
  single <- .buffer_joint_mean(
    matrix(42), 10, 24, 10, 24
  )
  expect_equal(single$value, 42)
  expect_equal(single$valid_fraction, 1)
  expect_equal(single$n_intersecting, 1)

  missing <- .buffer_joint_mean(
    matrix(NA_real_, nrow = 2L, ncol = 2L),
    c(5, 15), c(5, 15), 20, 20,
    min_valid_fraction = 0
  )
  expect_true(is.na(missing$value))
  expect_identical(missing$status, "all_missing")
  expect_equal(missing$coverage_fraction, 1)
  expect_equal(missing$valid_fraction, 0)
})
