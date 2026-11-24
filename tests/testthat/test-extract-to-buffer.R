# Helpers ----------------------------------------------------------------------

# make_grid_file() with the buffer tests' defaults: each cell holds its
# `cell_values` entry (lon x lat) plus `field` (depth x time). Dimensions are
# stored in a different order from lon, lat, depth, time on purpose. Other
# arguments, such as depth_bounds, time_bounds and var, go to make_grid_file(),
# except `positive`, `var_units`, `time_units` and `dim_order`, which this
# wrapper sets (bounds go through `lon_bounds`, `lat_bounds` and
# `spatial_bounds`).
.make_buffer_grid_nc <- function(
  path,
  cell_values = matrix(c(0, 10, 20, 30), 2L),
  lon = c(-0.01, 0.01),
  lat = c(-0.01, 0.01),
  lon_bounds = matrix(c(-0.02, 0, 0, 0.02), 2L),
  lat_bounds = matrix(c(-0.02, 0, 0, 0.02), 2L),
  spatial_bounds = TRUE,
  field = matrix(0),
  depth = NULL,
  time = NULL,
  ...
) {
  stopifnot(identical(dim(field), c(max(length(depth), 1L),
                                    max(length(time), 1L))))
  make_grid_file( # nolint: object_usage_linter.
    lon, lat, depth, time, path = path, values = outer(cell_values, field, "+"),
    lon_bounds = if (spatial_bounds) lon_bounds,
    lat_bounds = if (spatial_bounds) lat_bounds,
    positive = if (!is.null(depth)) "down", var_units = "degree_Celsius",
    time_units = "hours since 2020-01-01 00:00:00",
    dim_order = c("time", "lat", "depth", "lon"), ...
  )
}

# The cell values of the source files that cover lon 0 to 1 and lat 40 to 41.
.buffer_source_cells <- matrix(c(0, 100, 200, 300), 2L)


# Request validation -----------------------------------------------------------

test_that("most_similar exports locations without diagnostics", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path, depth = c(0, 20), time = c(0, 24))
  data <- data.frame(
    lon = 1, lat = 0, depth = 20, date = as.Date("2020-01-01"), target = 24030
  )
  run <- function(data, ...) {
    extract_to_buffer(
      data, path, var = "temp", lon_cell_radius = 1, lat_cell_radius = 0,
      depth_radius_m = 20, time_after = 1, time_unit = "days",
      statistic = "most_similar", target_col = "target", output_col = "picked",
      infer_bounds = TRUE, return_match_info = FALSE, quiet = TRUE, ...
    )
  }
  out <- run(data)
  expect_equal(out$picked, 24030)
  expect_equal(out$picked_match_lon, 1)
  expect_equal(out$picked_match_lat, 0)
  expect_equal(out$picked_match_depth, 20)
  expect_equal(out$picked_match_time, .track_clock(24))
  expect_identical(attr(out$picked_match_time, "tzone"), "UTC")
  expect_false("picked_match_info" %in% names(out))
  expect_equal(out[names(data)], data)
  missing <- run(transform(data, target = NA_real_))
  expect_true(all(is.na(missing[grep("picked", names(missing))])))
  rejected <- run(transform(data, target = 24031), max_target_difference = 0)
  expect_true(all(is.na(rejected[grep("picked", names(rejected))])))
})


test_that("buffer extraction stops on a climatological time axis", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path, depth = c(0, 20), time = c(0, 24),
                 time_atts = list(climatology = "climatology_bounds"))
  data <- data.frame(lon = 1, lat = 0, depth = 20,
                     date = as.Date("2020-01-01"))
  # Without a time buffer too: the file's time axis can't be decoded as dates.
  err <- expect_error(
    extract_to_buffer(data, path, var = "temp", depth_radius_m = 20,
                      infer_bounds = TRUE, quiet = TRUE),
    class = "ocean3d_error_time_climatology"
  )
  expect_match(conditionMessage(err), basename(path), fixed = TRUE)
})


test_that("identical winners expose ambiguity in only nonunique axes", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path, depth = c(0, 20), time = c(0, 24), values = 10)
  data <- data.frame(
    lon = 1, lat = 0, depth = 20, date = as.Date("2020-01-01"), target = 10
  )
  expect_warning(out <- extract_to_buffer(
    data, path, var = "temp", lon_cell_radius = 0, lat_cell_radius = 0,
    depth_radius_m = 20, time_after = 1, statistic = "most_similar",
    target_col = "target", infer_bounds = TRUE, quiet = TRUE
  ), class = "ocean3d_target_location_warning")
  expect_equal(out$temp, 10)
  expect_equal(out$temp_match_lon, 1)
  expect_equal(out$temp_match_lat, 0)
  expect_true(is.na(out$temp_match_depth))
  expect_true(is.na(out$temp_match_time))
  expect_equal(nrow(out$temp_match_info[[1]]$target_match$matches), 4)
  expect_identical(out$temp_match_info[[1]]$target_match$location_status,
                   "ambiguous_location")
})


test_that("location warnings are bounded and independent of quiet or reports", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path, values = 10)
  data <- data.frame(lon = rep(1, 12), lat = 0, target = 10)
  captured <- list()
  out <- withCallingHandlers(extract_to_buffer(
    data, path, var = "temp", lon_cell_radius = 1, lat_cell_radius = 0,
    statistic = "most_similar", target_col = "target", target_ties = "higher",
    infer_bounds = TRUE, return_match_info = FALSE, quiet = TRUE
  ), ocean3d_target_location_warning = function(warning) {
    captured[[length(captured) + 1L]] <<- warning
    invokeRestart("muffleWarning")
  })
  expect_length(captured, 1)
  expect_identical(captured[[1]]$rows, 1:12)
  expect_match(captured[[1]]$message, "Rows: 1, 2, 3, 4, 5; ...", fixed = TRUE)
  expect_match(captured[[1]]$message, "return_match_info = TRUE", fixed = TRUE)
  expect_match(captured[[1]]$message, "narrower buffer", fixed = TRUE)
  expect_match(captured[[1]]$message, "target_ties does not choose",
               fixed = TRUE)
  expect_equal(out$temp, rep(10, 12))
  expect_true(all(is.na(out$temp_match_lon)))
  expect_true(all(is.na(out$temp_match_lat)))
  expect_equal(out$lon, data$lon)
})


test_that("location columns protect input names and retain sf geometry", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path, extra_vars = list(static = 7))
  data <- data.frame(lon = 1, lat = 0, target = 12, static_target = 7)
  for (suffix in c("_match_lon", "_match_lat", "_match_depth", "_match_time")) {
    collision <- data
    collision[[paste0("temp", suffix)]] <- "existing"
    expect_error(extract_to_buffer(
      collision, path, var = "temp", horizontal_radius_km = 1,
      statistic = "most_similar", target_col = "target", infer_bounds = TRUE
    ), "overwrite")
  }
  out <- extract_to_buffer(
    data, path, var = c("temp", "static"), horizontal_radius_km = 1,
    statistic = "most_similar", target_col = c(
      temp = "target", static = "static_target"
    ), output_col = c(temp = "temperature", static = "identifier"),
    infer_bounds = TRUE, quiet = TRUE
  )
  expect_equal(out$temperature_match_lon, 1)
  expect_equal(out$identifier_match_lon, 1)
  expect_error(extract_to_buffer(
    data, path, var = c("temp", "static"), horizontal_radius_km = 1,
    statistic = "most_similar", target_col = c(
      temp = "target", static = "static_target"
    ), output_col = c(temp = "picked", static = "picked_match_lon"),
    infer_bounds = TRUE
  ), "overwrite")
  points <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326,
                         remove = FALSE)
  selected <- extract_to_buffer(
    points, path, var = "temp", horizontal_radius_km = 1,
    statistic = "most_similar", target_col = "target", infer_bounds = TRUE,
    quiet = TRUE
  )
  expect_s3_class(selected, "sf")
  expect_identical(sf::st_geometry(selected), sf::st_geometry(points))
  expect_equal(selected$temp_match_lon, 1)
})


test_that(
  "buffers require a positive window and valid scalar sizes",
  {
    expect_error(.validate_buffer_options(), "extract_to_point", fixed = TRUE)

    sizes <- c(
      "horizontal_radius_km",
      "depth_radius_m",
      "time_before",
      "time_after"
    )
    for (name in sizes) {
      for (value in list(-1, NA_real_, Inf, "5", c(1, 2), numeric(), NULL)) {
        args <- list(horizontal_radius_km = 5)
        args[name] <- list(value)
        expect_error(
          do.call(.validate_buffer_options, args),
          name,
          fixed = TRUE
        )
      }
    }

    options <- .validate_buffer_options(time_before = 0.5)
    expect_equal(options$time_before, 0.5)
    expect_identical(options$time_unit, "days")
    expect_identical(options$weighting, "weighted")
    expect_false(options$infer_bounds)
  }
)


test_that(
  "statistics, coverage thresholds and switches are validated",
  {
    expect_error(
      .validate_buffer_options(horizontal_radius_km = 5, statistic = "sum"),
      "`statistic` must be one of",
      fixed = TRUE
    )
    expect_error(
      .validate_buffer_options(horizontal_radius_km = 5, weighting = "other")
    )

    for (value in list(-0.1, 1.1, NA_real_, Inf, c(0.5, 0.8), NULL)) {
      expect_error(
        .validate_buffer_options(
          horizontal_radius_km = 5,
          min_valid_fraction = value
        ),
        "min_valid_fraction",
        fixed = TRUE
      )
    }

    for (name in c("infer_bounds", "return_match_info", "quiet")) {
      for (value in list(NA, 1, c(TRUE, FALSE), NULL)) {
        args <- list(horizontal_radius_km = 5)
        args[name] <- list(value)
        expect_error(
          do.call(.validate_buffer_options, args),
          name,
          fixed = TRUE
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
  }
)


test_that(
  "calendar sizes are whole numbers; f ... ns allow fractions",
  {
    for (unit in c("months", "years")) {
      for (name in c("time_before", "time_after")) {
        args <- list(time_unit = unit)
        args[[name]] <- 0.5
        expect_error(
          do.call(.validate_buffer_options, args),
          paste0("`", name, "` must be a whole number"),
          fixed = TRUE
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
  }
)


test_that(
  "depth intervals and single-layer methods are distinct",
  {
    options <- .validate_buffer_options(depth_radius_m = 10)
    expect_identical(options$method, "interval")

    for (method in c("nearest", "surface", "bottom", "2d")) {
      expect_error(
        .validate_buffer_options(depth_radius_m = 10, method = method),
        "Leave `method = NULL`",
        fixed = TRUE
      )
      options <- .validate_buffer_options(
        horizontal_radius_km = 5,
        method = method
      )
      expect_identical(options$method, method)
    }
    options <- .validate_buffer_options(horizontal_radius_km = 5)
    expect_identical(options$method, "nearest")
    expect_error(
      .validate_buffer_options(horizontal_radius_km = 5, method = "other")
    )
  }
)


# Observation windows ---------------------------------------------------------

test_that(
  "depth windows preserve requested support and missing rows",
  {
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
      "depth",
      fixed = TRUE
    )
    expect_error(
      .buffer_observation_windows(data, options, depth_col = NULL),
      "depth_col",
      fixed = TRUE
    )
    for (depth in list("60", -1, Inf)) {
      expect_error(
        .buffer_observation_windows(data.frame(depth = depth), options),
        "depth",
        fixed = TRUE
      )
    }
    custom <- .buffer_observation_windows(
      data.frame(z = 60),
      options,
      depth_col = "z"
    )
    expect_equal(custom$windows$depth_lower_m, 50)
    expect_equal(custom$windows$depth_upper_m, 70)
  }
)


test_that(
  "time windows use elapsed hours acro ... ght saving changes",
  {
    data <- data.frame(
      observed_at = as.POSIXct(
        c("2020-03-29 03:30:00", NA_character_),
        tz = "Europe/Madrid"
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
        data,
        options,
        date_col = "observed_at"
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
  }
)


test_that(
  "minutes and seconds preserve sub-minute window boundaries",
  {
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
  }
)


test_that(".buffer_shift_time() clamps calendar shifts to the month end", {
  shift <- function(x, amount, unit) {
    shifted <- .buffer_shift_time(as.POSIXct(x, tz = "UTC"), amount, unit)
    format(shifted, "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }
  ends <- c("2020-01-31 12:34:56", "2020-03-31 12:34:56",
            "2021-03-31 12:34:56", NA)
  expect_equal(shift(ends, 1, "months"),
               c("2020-02-29 12:34:56", "2020-04-30 12:34:56",
                 "2021-04-30 12:34:56", NA))
  expect_equal(shift(ends, -1, "months"),
               c("2019-12-31 12:34:56", "2020-02-29 12:34:56",
                 "2021-02-28 12:34:56", NA))
  # Years are 12 months, so a leap day moves to 28 February.
  leap_day <- "2020-02-29 19:30:00"
  expect_equal(shift(leap_day, 1, "years"), "2021-02-28 19:30:00")
  expect_equal(shift(leap_day, -1, "years"), "2019-02-28 19:30:00")
})

test_that(
  "a one-sided calendar window leaves the other endpoint unchanged",
  {
    options <- .validate_buffer_options(time_before = 1, time_unit = "years")
    result <- .buffer_observation_windows(
      data.frame(date = "2020-02-29T19:30:00Z"),
      options
    )
    expect_equal(
      format(result$windows$time_upper_utc, "%Y-%m-%d %H:%M", tz = "UTC"),
      "2020-02-29 19:30"
    )
  }
)


test_that(
  "Date centres use midnight and unuse ... inates are ignored",
  {
    options <- .validate_buffer_options(time_before = 3)
    result <- .buffer_observation_windows(
      data.frame(date = as.Date("2020-01-04"), depth = "unused"),
      options
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
      "date",
      fixed = TRUE
    )
    expect_error(
      .buffer_observation_windows(
        data.frame(date = NA),
        options,
        date_col = NULL
      ),
      "date_col",
      fixed = TRUE
    )

    options <- .validate_buffer_options(horizontal_radius_km = 5)
    result <- .buffer_observation_windows(
      data.frame(date = "unused", depth = "unused"),
      options
    )
    expect_equal(nrow(result$windows), 1L)
    expect_null(result$prepared_time)
  }
)


test_that(
  "combined windows retain timezone in ... or a later warning",
  {
    options <- .validate_buffer_options(
      horizontal_radius_km = 5,
      depth_radius_m = 10,
      time_before = 1
    )
    result <- .buffer_observation_windows(
      data.frame(depth = 60, date = "2020-01-02 12:00"),
      options
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
  }
)


test_that(
  "missing and empty rows are preserved across window types",
  {
    for (unit in c("days", "months", "years")) {
      options <- .validate_buffer_options(
        depth_radius_m = 10,
        time_before = 1,
        time_unit = unit
      )
      result <- .buffer_observation_windows(
        data.frame(depth = c(NA_real_, NA_real_), date = c(NA, NA)),
        options
      )
      expect_equal(result$windows$row, 1:2)
      expect_true(all(is.na(result$windows$depth_lower_m)))
      expect_true(all(is.na(result$windows$time_lower_utc)))
      expect_true(all(is.na(result$windows$time_upper_utc)))

      empty <- .buffer_observation_windows(
        data.frame(depth = numeric(), date = character()),
        options
      )
      expect_equal(nrow(empty$windows), 0L)
    }
  }
)


test_that(
  "temporal conversion errors retain c ... bservation context",
  {
    options <- .validate_buffer_options(time_before = 1)
    data <- data.frame(
      sample = c("a", "b"),
      observed_at = c("2020-01-02", "06/01/2020")
    )
    error <- tryCatch(
      .buffer_observation_windows(
        data,
        options,
        date_col = "observed_at",
        id_col = "sample"
      ),
      error = identity
    )
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      "Input: `data$observed_at`",
      fixed = TRUE
    )
    expect_match(
      conditionMessage(error),
      "row/position 2 (id: \"b\")",
      fixed = TRUE
    )
    expect_match(conditionMessage(error), "06/01/2020", fixed = TRUE)
  }
)


# Horizontal buffer geometry ---------------------------------------------------

# Degrees of longitude per kilometre along the WGS84 equator.
.km_to_equator_degrees <- function(km) km / (6378.137 * pi / 180)


test_that("horizontal buffers use kilometres and geodesic distances", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(
    path, matrix(c(0, 10), 2L, 1L), lat = 0,
    lat_bounds = matrix(c(-0.02, 0.02), 2L)
  )
  # The meridian at lon 0 cuts a 1 km circle 0.5 km from its centre; the
  # eastern cell holds the circular segment beyond that chord.
  segment <- (acos(0.5) - 0.5 * sqrt(0.75)) / pi
  out <- extract_to_buffer(
    data.frame(lon = -.km_to_equator_degrees(0.5), lat = 0), path,
    horizontal_radius_km = 1, quiet = TRUE
  )
  expect_equal(out$temp, 10 * segment, tolerance = 1e-4)
  expect_equal(out$temp_match_info[[1L]]$coverage_fraction, 1)
})


test_that("circles containing a pole cover the cells around it", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  # Four 90-degree cells meet at the north pole.
  .make_buffer_grid_nc(
    path, matrix(c(0, 10, 20, 30), 4L, 1L), lon = c(-135, -45, 45, 135),
    lat = 89.95, lon_bounds = matrix(c(-180, -90, -90, 0, 0, 90, 90, 180), 2L),
    lat_bounds = matrix(c(89.9, 90), 2L)
  )
  out <- extract_to_buffer(
    data.frame(lon = c(0, 45), lat = c(90, 89.99)), path,
    horizontal_radius_km = 5, quiet = TRUE
  )
  info <- out$temp_match_info
  # A circle centred on the pole covers the four cells equally.
  expect_equal(out$temp[1L], 15, tolerance = 1e-6)
  for (i in 1:2) {
    expect_equal(info[[i]]$coverage_fraction, 1, tolerance = 1e-6)
    expect_equal(info[[i]]$n_valid, 4)
  }
  # Off the pole, the cell under the centre (value 20) covers most of it.
  expect_gt(out$temp[2L], 15)
})


test_that("longitude cells that overlap modulo 360 are rejected", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  # Increasing and non-overlapping as stored, but [340, 360] is [-20, 0].
  .make_buffer_grid_nc(
    path, matrix(1:4, 4L, 1L), lon = c(0, 120, 240, 350), lat = 0,
    lon_bounds = matrix(c(-15, 15, 105, 135, 225, 255, 340, 360), 2L),
    lat_bounds = matrix(c(-1, 1), 2L)
  )
  err <- expect_error(
    extract_to_buffer(data.frame(lon = 0, lat = 0), path,
                      horizontal_radius_km = 1, quiet = TRUE),
    class = "ocean3d_error_bounds"
  )
  expect_match(conditionMessage(err), "overlap around the globe")
})

test_that("invalid horizontal requests are rejected", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(path)
  data <- data.frame(lon = 0, lat = 0)
  for (radius in c(-1, 0, NA_real_, Inf, 10000)) {
    expect_error(extract_to_buffer(data, path, horizontal_radius_km = radius))
  }
  out <- extract_to_buffer(
    data.frame(lon = 0, lat = 91), path, horizontal_radius_km = 1, quiet = TRUE
  )
  expect_identical(out$temp_match_info[[1L]]$status, "invalid_coordinates")
  expect_error(
    extract_to_buffer(data, path, horizontal_radius_km = 1,
                      lon_cell_radius = 1),
    "not both", fixed = TRUE
  )
  .make_buffer_grid_nc(path, spatial_bounds = FALSE)
  expect_error(
    extract_to_buffer(data, path, horizontal_radius_km = 1),
    "infer_bounds = TRUE", fixed = TRUE
  )
})


# Public extraction and complementary sources ---------------------------------

test_that("cell support distinguishes halves and quarters of the buffer", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(path, matrix(c(0, 0, 0, 40), 2L))
  quarters <- extract_to_buffer(
    data.frame(lon = 0, lat = 0), path, horizontal_radius_km = 1, quiet = TRUE
  )
  expect_equal(quarters$temp, 10, tolerance = 1e-6)

  .make_buffer_grid_nc(
    path, matrix(c(0, 40), 2L, 1L), lat = 0,
    lat_bounds = matrix(c(-0.02, 0.02), 2L)
  )
  halves <- extract_to_buffer(
    data.frame(lon = 0, lat = 0), path, horizontal_radius_km = 1, quiet = TRUE
  )
  expect_equal(halves$temp, 20, tolerance = 1e-6)
  expect_equal(halves$temp_match_info[[1L]]$n_valid, 2)
})


test_that("partial and empty horizontal support are kept separate", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(
    path, matrix(10, 1L, 1L), lon = 0.01, lat = 0,
    lon_bounds = matrix(c(0, 0.02), 2L), lat_bounds = matrix(c(-0.02, 0.02), 2L)
  )
  # The second circle ends 5 m west of the only cell, close enough for the
  # cell to be a candidate.
  expect_no_warning(out <- extract_to_buffer(
    data.frame(lon = c(0, -.km_to_equator_degrees(1.005)), lat = 0), path,
    horizontal_radius_km = 1, min_valid_fraction = 0.5, quiet = TRUE
  ))
  info <- out$temp_match_info
  expect_equal(out$temp, c(10, NA_real_))
  expect_equal(info[[1L]]$coverage_fraction, 0.5, tolerance = 1e-6)
  expect_identical(info[[1L]]$status, "partial")
  expect_identical(info[[2L]]$status, "empty")
  expect_equal(info[[2L]]$n_intersecting, 0)
})


test_that("horizontal geometry does not depend on or change the s2 setting", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_buffer_grid_nc(path)
  original <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(original)), add = TRUE)

  values <- numeric()
  for (engine in c(TRUE, FALSE)) {
    suppressMessages(sf::sf_use_s2(engine))
    out <- extract_to_buffer(
      data.frame(lon = 0.004, lat = 0.003), path,
      horizontal_radius_km = 1, quiet = TRUE
    )
    expect_identical(sf::sf_use_s2(), engine)
    values <- c(values, out$temp)
  }
  expect_equal(values[1L], values[2L])
})


test_that(
  "horizontal means preserve inputs",
  {
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
      data,
      path,
      horizontal_radius_km = 1,
      return_match_info = FALSE,
      quiet = TRUE
    )
    expect_named(simple, c(names(data), "temp"))
  }
)


test_that(
  "joint means combine area, depth & time",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_buffer_grid_nc(
      path,
      cell_values = matrix(c(0, 100, 200, 300), 2L),
      field = matrix(c(0, 10, 20, 40), 2L),
      depth = c(5, 20),
      time = c(5, 20),
      depth_bounds = matrix(c(0, 10, 10, 30), 2L),
      time_bounds = matrix(c(0, 10, 10, 30), 2L)
    )
    data <- data.frame(lon = 0, lat = 0, depth = 15, date = "2020-01-01T15:00Z")
    for (weighting in c("weighted", "equal")) {
      out <- extract_to_buffer(
        data,
        path,
        var = "temp",
        horizontal_radius_km = 1,
        depth_radius_m = 10,
        time_before = 10,
        time_after = 10,
        time_unit = "hours",
        weighting = weighting,
        quiet = TRUE
      )
      expect_equal(
        out$temp,
        if (weighting == "weighted") 178.125 else 167.5,
        tolerance = 1e-6
      )
      expect_equal(out$temp_match_info[[1L]]$n_valid, 16)
      expect_equal(out$temp_match_info[[1L]]$valid_fraction, 1)
    }
  }
)


test_that(
  "spatial tiles share edges without double-counting",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    .make_buffer_grid_nc(
      paths[1L],
      matrix(10, 1L, 2L),
      lon = -0.01,
      lon_bounds = matrix(c(-0.02, 0), 2L)
    )
    .make_buffer_grid_nc(
      paths[2L],
      matrix(30, 1L, 2L),
      lon = 0.01,
      lon_bounds = matrix(c(0, 0.02), 2L)
    )
    for (files in list(paths, rev(paths))) {
      out <- extract_to_buffer(
        data.frame(lon = 0, lat = 0),
        files,
        var = "temp",
        horizontal_radius_km = 1,
        quiet = TRUE
      )
      expect_equal(out$temp, 20, tolerance = 1e-6)
      expect_equal(out$temp_match_info[[1L]]$valid_fraction, 1)
    }
  }
)


test_that(
  "overlapping temporal and spatial so ... ejected explicitly",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    .make_buffer_grid_nc(
      paths[1L], .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(10),
      time = 5,
      time_bounds = matrix(c(0, 20), 2L)
    )
    .make_buffer_grid_nc(
      paths[2L], .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(20),
      time = 20,
      time_bounds = matrix(c(10, 30), 2L)
    )
    error <- tryCatch(
      extract_to_buffer(
        data.frame(lon = 0, lat = 40, date = "2020-01-01T15:00Z"),
        paths,
        var = "temp",
        time_before = 10,
        time_after = 10,
        time_unit = "hours",
        quiet = TRUE
      ),
      error = identity
    )
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), "Source files overlap", fixed = TRUE)
    for (path in paths) expect_match(
      conditionMessage(error),
      basename(path),
      fixed = TRUE
    )
    .make_buffer_grid_nc(paths[1L])
    expect_error(
      extract_to_buffer(
        data.frame(lon = 0, lat = 0),
        c(paths[1L], paths[1L]),
        var = "temp",
        horizontal_radius_km = 1,
        quiet = TRUE
      ),
      "Source files overlap",
      fixed = TRUE
    )
  }
)


test_that(
  "point time is resolved globally before a vertical mean",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    bounds <- matrix(c(0, 10, 10, 30), 2L)
    .make_buffer_grid_nc(
      paths[1L], .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(c(10, 20, 50, 60), 2L),
      depth = c(5, 20),
      time = c(0, 48),
      depth_bounds = bounds
    )
    .make_buffer_grid_nc(
      paths[2L], .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(c(30, 40), 2L),
      depth = c(5, 20),
      time = 24,
      depth_bounds = bounds
    )
    data <- data.frame(
      lon = 0,
      lat = 40,
      depth = 15,
      date = as.Date("2020-01-02")
    )
    for (files in list(paths, rev(paths))) {
      out <- extract_to_buffer(
        data,
        files,
        var = "temp",
        depth_radius_m = 10,
        quiet = TRUE
      )
      expect_equal(out$temp, 37.5)
      expect_length(out$temp_match_info[[1L]]$sources_used, 1)
    }
  }
)


test_that(
  "point depth selection combines source layers",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    .make_buffer_grid_nc(
      paths[1L], .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(10),
      depth = 5,
      time = 12,
      time_bounds = matrix(c(0, 24), 2L)
    )
    .make_buffer_grid_nc(
      paths[2L], .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(20),
      depth = 20,
      time = 12,
      time_bounds = matrix(c(0, 24), 2L)
    )
    data <- data.frame(
      lon = 0,
      lat = 40,
      depth = 18,
      date = "2020-01-01T12:00Z"
    )
    for (files in list(paths, rev(paths))) {
      out <- extract_to_buffer(
        data,
        files,
        var = "temp",
        time_before = 12,
        time_after = 12,
        time_unit = "hours",
        quiet = TRUE
      )
      expect_equal(out$temp, 20)
      expect_equal(out$temp_match_info[[1L]]$depth_range_m, c(20, 20))
    }
  }
)


test_that(
  "index layer windows retain physical ...  reject grid edges",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_buffer_grid_nc(
      path, .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(c(30, 20, 10), 3L),
      depth = c(35, 15, 5),
      depth_bounds = matrix(c(20, 50, 10, 20, 0, 10), 2L)
    )
    data <- data.frame(lon = 0, lat = 40, depth = c(15, 5))
    for (weighting in c("weighted", "equal")) {
      out <- extract_to_buffer(
        data,
        path,
        depth_layer_radius = 1,
        weighting = weighting,
        min_valid_fraction = 0,
        quiet = TRUE
      )
      expect_equal(
        out$temp,
        c(
          if (
            weighting == "weighted"
          ) 24 else 20,
          NA_real_
        )
      )
      expect_identical(
        out$temp_match_info[[2L]]$status,
        "incomplete_index_window"
      )
    }
  }
)


test_that(
  "index windows retain edge NAs",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_buffer_grid_nc(
      path,
      matrix(c(0, 10, 20), 3L),
      lon = c(-0.01, 0.01, 0.03),
      lat = 0,
      lon_bounds = matrix(c(-0.02, 0, 0, 0.02, 0.02, 0.04), 2L),
      lat_bounds = matrix(c(-0.01, 0.01), 2L)
    )
    out <- extract_to_buffer(
      data.frame(lon = c(0.01, -0.01), lat = 0),
      path,
      lon_cell_radius = 1,
      lat_cell_radius = 0,
      quiet = TRUE
    )
    expect_equal(out$temp, c(10, NA_real_), tolerance = 1e-5)
    expect_identical(
      out$temp_match_info[[2L]]$status,
      "incomplete_index_window"
    )
  }
)


test_that(
  "actual antimeridian cells intersect the requested circle",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_buffer_grid_nc(
      path,
      lon = c(179.99, 180.01),
      lon_bounds = matrix(c(179.98, 180, 180, 180.02), 2L)
    )
    out <- extract_to_buffer(
      data.frame(lon = c(180, -180), lat = 0),
      path,
      horizontal_radius_km = 1,
      quiet = TRUE
    )
    expect_equal(out$temp, c(15, 15), tolerance = 1e-5)
    expect_equal(out$temp_match_info[[1L]]$valid_fraction, 1)
  }
)


test_that(
  "multiple outputs have collision-safe diagnostics",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    .make_buffer_grid_nc(paths[1L])
    .make_buffer_grid_nc(paths[2L], matrix(c(2, 4, 6, 8), 2L), var = "oxygen")
    data <- data.frame(lon = 0, lat = 0)
    out <- extract_to_buffer(
      data,
      paths,
      var = c("temp", "oxygen"),
      horizontal_radius_km = 1,
      output_col = c(temp = "temperature", oxygen = "oxygen_mean"),
      quiet = TRUE
    )
    expect_named(
      out,
      c(
        "lon",
        "lat",
        "temperature",
        "temperature_match_info",
        "oxygen_mean",
        "oxygen_mean_match_info"
      )
    )
    expect_equal(out$temperature, 15, tolerance = 1e-6)
    expect_equal(out$oxygen_mean, 5, tolerance = 1e-6)
    data$temp <- 999
    expect_error(
      extract_to_buffer(data, paths[1L], horizontal_radius_km = 1),
      "overwrite",
      fixed = TRUE
    )
  }
)


test_that(
  "coverage reports explain partial results",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_buffer_grid_nc(
      path,
      matrix(10, 1L, 2L),
      lon = 0.01,
      lon_bounds = matrix(c(0, 0.02), 2L)
    )
    data <- data.frame(lon = 0, lat = 0)
    expect_message(
      out <- extract_to_buffer(
        data,
        path,
        horizontal_radius_km = 1,
        min_valid_fraction = 0.5
      ),
      "partial coverage",
      fixed = TRUE
    )
    expect_equal(out$temp, 10)
    expect_equal(
      out$temp_match_info[[1L]]$valid_fraction,
      0.5,
      tolerance = 1e-6
    )
    expect_no_message(
      rejected <- extract_to_buffer(
        data,
        path,
        horizontal_radius_km = 1,
        quiet = TRUE
      )
    )
    expect_true(is.na(rejected$temp))
  }
)


test_that(
  "unzoned text times warn even when quiet",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_buffer_grid_nc(
      path, .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(10),
      time = 12,
      time_bounds = matrix(c(0, 24), 2L)
    )
    expect_warning(
      timed <- extract_to_buffer(
        data.frame(lon = 0, lat = 40, date = "2020-01-01 12:00"),
        path,
        time_before = 12,
        time_after = 12,
        time_unit = "hours",
        quiet = TRUE
      ),
      class = "ocean3d_timezone_warning"
    )
    expect_equal(timed$temp, 10)
  }
)


test_that("a caller-owned connection is read and left open", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  .make_buffer_grid_nc(path)
  handle <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(handle), add = TRUE, after = FALSE)
  out <- extract_to_buffer(
    data.frame(lon = 0, lat = 0), handle, horizontal_radius_km = 1,
    quiet = TRUE
  )
  expect_equal(out$temp, 15, tolerance = 1e-6)
  expect_equal(as.vector(ncdf4::ncvar_get(handle, "temp")), c(0, 20, 10, 30))
})

# Axis bounds ------------------------------------------------------------------

# Source extraction ------------------------------------------------------------

test_that("buffer SD and SE of identical values are zero", {
  for (statistic in c("sd", "se")) {
    expect_equal(.weighted_statistic(c(7, 7, 7), c(1, 2, 3), statistic), 0)
  }
})

test_that(
  "statistics combine files in either order",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    .make_buffer_grid_nc(
      paths[1L], .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(10),
      time = 5,
      time_bounds = matrix(c(0, 10), 2L)
    )
    .make_buffer_grid_nc(
      paths[2L], .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(20),
      time = 20,
      time_bounds = matrix(c(10, 30), 2L)
    )
    data <- data.frame(lon = 0, lat = 40, date = "2020-01-01T15:00Z")
    statistics <- c("mean", "sd", "se", "min", "max")
    for (files in list(paths, rev(paths))) for (
      weighting in c(
        "weighted",
        "equal"
      )
    ) {
      expected <- if (weighting == "weighted") {
        c(17.5, sqrt(50), sqrt(31.25), 10, 20)
      } else {
        c(15, sqrt(50), 5, 10, 20)
      }
      for (i in seq_along(statistics)) {
        out <- extract_to_buffer(
          data,
          files,
          var = "temp",
          time_before = 10,
          time_after = 10,
          time_unit = "hours",
          statistic = statistics[i],
          weighting = weighting,
          output_col = "summary",
          quiet = TRUE
        )
        expect_named(out, c(names(data), "summary", "summary_match_info"))
        expect_equal(out$summary, expected[i])
        info <- out$summary_match_info[[1L]]
        expect_identical(info$statistic, statistics[i])
        expect_equal(info$effective_n, if (weighting == "weighted") 1.6 else 2)
        expect_identical(info$status, "complete")
        expect_equal(info$coverage_fraction, 1)
        expect_length(info$sources_used, 2)
      }
    }
  }
)


test_that(
  "all statistics respect missing valu ...  physical coverage",
  {
    samples <- data.frame(
      value = c(10, 20, NA_real_),
      area_fraction = c(0.2, 0.6, 0.2),
      depth_fraction = 1,
      time_fraction = 1
    )
    statistics <- c("mean", "sd", "se", "min", "max")
    expected <- c(17.5, sqrt(50), sqrt(31.25), 10, 20)
    for (i in seq_along(statistics)) {
      out <- .buffer_summarise_samples(samples, "weighted", 0.8, statistics[i])
      expect_equal(out$value, expected[i])
      expect_equal(out$valid_fraction, 0.8)
      expect_equal(out$effective_n, 1.6)
      expect_identical(out$status, "partial")
      rejected <- .buffer_summarise_samples(
        samples,
        "weighted",
        0.9,
        statistics[i]
      )
      expect_true(is.na(rejected$value))
      expect_identical(rejected$status, "insufficient_coverage")
      missing <- samples
      missing$value[] <- NA_real_
      rejected <- .buffer_summarise_samples(
        missing,
        "weighted",
        0,
        statistics[i]
      )
      expect_true(is.na(rejected$value))
      expect_identical(rejected$status, "all_missing")
      empty <- .buffer_summarise_samples(
        samples[FALSE, ],
        "equal",
        0,
        statistics[i]
      )
      expect_true(is.na(empty$value))
      expect_identical(empty$status, "empty")
    }
  }
)


test_that(
  "public SD and SE explain an insuffi ... r of valid samples",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_buffer_grid_nc(
      path, .buffer_source_cells,
      lon = c(0, 1), lat = c(40, 41), spatial_bounds = FALSE,
      field = matrix(10),
      time = 5,
      time_bounds = matrix(c(0, 10), 2L)
    )
    data <- data.frame(lon = 0, lat = 40, date = "2020-01-01T05:00Z")
    for (statistic in c("sd", "se")) {
      expect_message(
        out <- extract_to_buffer(
          data,
          path,
          time_before = 5,
          time_after = 5,
          time_unit = "hours",
          statistic = statistic
        ),
        "at least two valid samples",
        fixed = TRUE
      )
      expect_true(is.na(out$temp))
      expect_identical(
        out$temp_match_info[[1L]]$status,
        "insufficient_variation_samples"
      )
      expect_equal(out$temp_match_info[[1L]]$effective_n, 1)
      expect_message(
        extract_to_buffer(
          data,
          path,
          time_before = 5,
          time_after = 5,
          time_unit = "hours",
          statistic = statistic,
          quiet = TRUE
        ),
        NA
      )
    }
    for (statistic in c("mean", "min", "max")) {
      out <- extract_to_buffer(
        data,
        path,
        time_before = 5,
        time_after = 5,
        time_unit = "hours",
        statistic = statistic
      )
      expect_equal(out$temp, 10)
      expect_identical(out$temp_match_info[[1L]]$status, "complete")
    }
  }
)


test_that(
  "most_similar selects values and only resolves genuine ties",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    data <- data.frame(lon = 1, lat = 0, target = 12)
    run <- function(data = data, ...) {
      extract_to_buffer(
        data,
        path,
        var = "temp",
        lon_cell_radius = 1,
        lat_cell_radius = 0,
        infer_bounds = TRUE,
        statistic = "most_similar",
        target_col = "target",
        quiet = TRUE,
        ...
      )
    }
    out <- run(data)
    expect_equal(out$temp, 10)
    expect_equal(
      out$temp_match_info[[1L]]$target_match$absolute_difference,
      2
    )
    expect_equal(
      out$temp_match_info[[1L]]$target_match$matches$longitude,
      1
    )
    expect_equal(
      run(
        transform(data, target = 15),
        target_ties = "lower"
      )$temp,
      10
    )
    expect_equal(
      run(
        transform(data, target = 15),
        target_ties = "higher"
      )$temp,
      20
    )
    expect_warning(
      tied <- run(transform(data, target = 15)),
      class = "ocean3d_target_location_warning"
    )
    expect_true(is.na(tied$temp))
    expect_equal(tied$temp_match_info[[1L]]$status, "ambiguous_target")
    expect_equal(
      nrow(tied$temp_match_info[[1L]]$target_match$ties),
      2
    )
    expect_true(is.na(run(data, max_target_difference = 1)$temp))
    expect_equal(run(data, target_ties = "higher")$temp, 10)
    expect_true(is.na(run(transform(data, target = NA_real_))$temp))
    expect_error(run(transform(data, target = Inf)), "finite")
    expect_error(run(transform(data, target = "12")), "numeric")
  }
)


test_that(
  "most_similar accepts partial valid  ...  retain 80 percent",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, missing = TRUE)
    data <- data.frame(lon = 1, lat = 0, target = 11)
    args <- list(
      data = data,
      nc = path,
      var = "temp",
      lon_cell_radius = 1,
      lat_cell_radius = 0,
      infer_bounds = TRUE,
      quiet = TRUE
    )
    out <- do.call(
      extract_to_buffer,
      c(
        args,
        list(
          statistic = "most_similar",
          target_col = "target"
        )
      )
    )
    expect_equal(out$temp, 10)
    expect_equal(
      out$temp_match_info[[1L]]$valid_fraction,
      2 / 3,
      tolerance = 1e-05
    )
    expect_true(is.na(do.call(extract_to_buffer, args)$temp))
    expect_true(
      is.na(
        do.call(
          extract_to_buffer,
          c(
            args,
            list(
              statistic = "most_similar",
              target_col = "target",
              min_valid_fraction = 0.8
            )
          )
        )$temp
      )
    )
  }
)


test_that(
  "target selection across disjoint so ... nt of source order",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    .make_track_nc(paths[1L], lon = 0:2)
    .make_track_nc(paths[2L], lon = 3:5)
    data <- data.frame(lon = 2.5, lat = 0, target = 25)
    for (sources in list(paths, rev(paths))) {
      args <- list(
        data = data,
        nc = sources,
        var = "temp",
        horizontal_radius_km = 100,
        infer_bounds = TRUE,
        statistic = "most_similar",
        target_col = "target",
        quiet = TRUE
      )
      expect_warning(
        out <- do.call(extract_to_buffer, args),
        class = "ocean3d_target_location_warning"
      )
      expect_true(is.na(out$temp))
      expect_equal(
        sort(out$temp_match_info[[1L]]$target_match$ties$value),
        c(20, 30)
      )
      expect_equal(
        do.call(
          extract_to_buffer,
          c(
            args,
            list(target_ties = "lower")
          )
        )$temp,
        20
      )
      high <- do.call(
        extract_to_buffer,
        c(
          args,
          list(target_ties = "higher")
        )
      )
      expect_equal(high$temp, 30)
      files_used <- high$temp_match_info[[1L]]$target_match$matches$source_file
      expect_equal(
        normalizePath(files_used),
        normalizePath(paths[2L])
      )
    }
  }
)


test_that(
  "identical winners retain all locati ...  buffers remain NA",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, values = 10)
    data <- data.frame(lon = 1, lat = 0, target = 10)
    args <- list(
      data = data,
      nc = path,
      var = "temp",
      lon_cell_radius = 1,
      lat_cell_radius = 0,
      infer_bounds = TRUE,
      statistic = "most_similar",
      target_col = "target",
      quiet = TRUE
    )
    expect_warning(
      out <- do.call(extract_to_buffer, args),
      class = "ocean3d_target_location_warning"
    )
    expect_equal(out$temp, 10)
    expect_equal(
      out$temp_match_info[[1L]]$target_match$matches$longitude,
      0:2
    )
    handle <- ncdf4::nc_open(path, write = TRUE)
    ncdf4::ncvar_put(handle, "temp", matrix(NA_real_, 3, 2))
    ncdf4::nc_close(handle)
    expect_true(is.na(do.call(extract_to_buffer, args)$temp))
  }
)


test_that(
  "target mappings allow several varia ... uous configuration",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, extra_vars = list(static = 7))
    data <- data.frame(
      lon = 1,
      lat = 0,
      target_temp = 12,
      target_static = 6
    )
    out <- extract_to_buffer(
      data,
      path,
      var = c(
        "temp",
        "static"
      ),
      statistic = "most_similar",
      target_col = c(
        temp = "target_temp",
        static = "target_static"
      ),
      horizontal_radius_km = 1,
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_equal(out$temp, 10)
    expect_equal(out$static, 7)
    expect_error(
      extract_to_buffer(data, path, target_ties = "higher"),
      "require"
    )
    expect_error(
      extract_to_buffer(
        data,
        path,
        var = "temp",
        statistic = "most_similar",
        target_col = "absent",
        horizontal_radius_km = 1,
        infer_bounds = TRUE
      ),
      "absent"
    )
  }
)

test_that(
  "a target selects an actual candidate without averaging",
  {
    candidates <- data.frame(
      candidate_id = c("a", "b", "c"),
      value = c(
        12,
        15,
        20
      ),
      support_fraction = c(0.2, 0.3, 0.5),
      file = c(
        "x.nc",
        "y.nc",
        "z.nc"
      )
    )
    out <- .select_target_value(candidates, 16)
    expect_equal(out$value, 15)
    expect_equal(out$absolute_difference, 1)
    expect_identical(out$status, "matched")
    expect_identical(candidates$file[out$selected_rows], "y.nc")
    expect_equal(out$valid_fraction, 1)
    expect_equal(.select_target_value(candidates, 12)$value, 12)
  }
)

test_that(
  "file or candidate order never resolves a tie silently",
  {
    candidates <- data.frame(
      candidate_id = c("a", "b"),
      value = c(
        14,
        18
      ),
      support_fraction = 0.5
    )
    for (data in list(candidates, candidates[2:1, ])) {
      ambiguous <- .select_target_value(data, 16)
      expect_true(is.na(ambiguous$value))
      expect_identical(ambiguous$status, "ambiguous_target")
      expect_equal(ambiguous$n_tied_candidates, 2)
      expect_equal(
        .select_target_value(data, 16, ties = "lower")$value,
        14
      )
      expect_equal(
        .select_target_value(data, 16, ties = "higher")$value,
        18
      )
    }
  }
)

test_that(
  "the same winning value can have several honest locations",
  {
    candidates <- data.frame(
      candidate_id = c("a", "b", "c"),
      value = c(
        15,
        15,
        30
      ),
      support_fraction = c(0.2, 0.3, 0.5)
    )
    out <- .select_target_value(candidates, 16)
    expect_equal(out$value, 15)
    expect_equal(out$selected_rows, 1:2)
    expect_identical(out$status, "matched_multiple_locations")
    expect_equal(out$n_tied_candidates, 2)
    expect_true(is.na(.select_target_value(candidates, NA_real_)$value))
    expect_identical(
      .select_target_value(candidates, NA_real_)$status,
      "missing_target"
    )
  }
)

test_that(
  "coverage includes missing data and zero support cannot win",
  {
    candidates <- data.frame(
      candidate_id = c("a", "b", "outside"),
      value = c(20, NA, 16),
      support_fraction = c(0.79, 0.21, 0)
    )
    rejected <- .select_target_value(candidates, 16)
    expect_identical(rejected$status, "insufficient_coverage")
    expect_true(is.na(rejected$value))
    accepted <- .select_target_value(
      candidates,
      16,
      min_valid_fraction = 0.79
    )
    expect_equal(accepted$value, 20)
    expect_equal(accepted$valid_fraction, 0.79)
    expect_equal(accepted$coverage_fraction, 1)
    expect_equal(accepted$n_valid, 1)
    candidates$value[1] <- Inf
    expect_identical(
      .select_target_value(candidates, 16)$status,
      "all_missing"
    )
    candidates$support_fraction <- 0
    expect_identical(
      .select_target_value(candidates, 16)$status,
      "empty"
    )
    expect_identical(
      .select_target_value(
        candidates[FALSE, ],
        16
      )$status,
      "empty"
    )
  }
)

test_that(
  "environmental difference limits retain diagnostics",
  {
    candidates <- data.frame(
      candidate_id = "a",
      value = 15,
      support_fraction = 1
    )
    rejected <- .select_target_value(
      candidates,
      16,
      max_target_difference = 0.5
    )
    expect_true(is.na(rejected$value))
    expect_equal(rejected$closest_available_value, 15)
    expect_equal(rejected$absolute_difference, 1)
    expect_identical(rejected$status, "outside_target_limit")
    accepted <- .select_target_value(candidates, 16, max_target_difference = 1)
    expect_equal(accepted$value, 15)
  }
)

test_that(
  "weights do not override proximity once coverage is valid",
  {
    candidates <- data.frame(
      candidate_id = c("a", "b"),
      value = c(
        15.9,
        30
      ),
      support_fraction = c(0.01, 0.99)
    )
    expect_equal(.select_target_value(candidates, 16)$value, 15.9)
    candidates$support_fraction <- c(0.99, 0.01)
    expect_equal(.select_target_value(candidates, 16)$value, 15.9)
  }
)

test_that(
  "numerical ties are handled without  ...  exact target hits",
  {
    candidates <- data.frame(
      candidate_id = c("a", "b"),
      value = c(
        10.1,
        10.3
      ),
      support_fraction = 0.5
    )
    out <- .select_target_value(candidates, 10.2)
    expect_identical(out$status, "ambiguous_target")
    candidates$value <- c(1e+12, 1e+12 + 0.001)
    expect_equal(
      .select_target_value(candidates, 1e+12)$selected_rows,
      1
    )
    candidates$value <- c(1e-20, 3e-20)
    expect_equal(
      .select_target_value(candidates, 1.1e-20)$value,
      1e-20
    )
    candidates$value <- c(1, 1.2)
    expect_equal(
      .select_target_value(
        candidates,
        0,
        tie_tolerance = 0
      )$value,
      1
    )
  }
)

test_that(
  "invalid requests fail before any selection",
  {
    candidates <- data.frame(
      candidate_id = "a",
      value = 15,
      support_fraction = 1
    )
    for (bad in list("16", Inf, c(1, 2), numeric())) {
      expect_error(.select_target_value(candidates, bad), "Target")
    }
    for (bad in list(-1, 2, NA, "half")) {
      expect_error(
        .select_target_value(candidates, 16, min_valid_fraction = bad),
        "min_valid_fraction"
      )
    }
    for (bad in list(-1, Inf, NA, c(1, 2))) {
      expect_error(
        .select_target_value(
          candidates,
          16,
          max_target_difference = bad
        ),
        "max_target_difference"
      )
    }
    expect_error(.select_target_value(candidates[FALSE], 16), "Candidates")
    doubled <- rbind(candidates, candidates)
    expect_error(.select_target_value(doubled, 16), "IDs")
    candidates$support_fraction <- -1
    expect_error(.select_target_value(candidates, 16), "non-negative")
    candidates$support_fraction <- 1.1
    expect_error(.select_target_value(candidates, 16), "overlap")
  }
)

test_that(
  "tie tolerance cannot override the e ... l difference limit",
  {
    candidates <- data.frame(
      candidate_id = c("a", "b"),
      value = c(
        1,
        1.2
      ),
      support_fraction = 0.5
    )
    out <- .select_target_value(
      candidates,
      0,
      max_target_difference = 1.1,
      tie_tolerance = 0.3,
      ties = "higher"
    )
    expect_equal(out$value, 1)
    expect_lte(out$absolute_difference, 1.1)
    expect_equal(out$n_tied_candidates, 1)
    expect_error(
      .select_target_value(
        candidates,
        0,
        tie_tolerance = -1
      ),
      "tie_tolerance"
    )
  }
)

test_that(
  "unrepresentable differences produce an explicit status",
  {
    candidates <- data.frame(
      candidate_id = "a",
      value = 1e+308,
      support_fraction = 1
    )
    out <- .select_target_value(candidates, -1e+308)
    expect_true(is.na(out$value))
    expect_identical(out$status, "difference_overflow")
  }
)
test_that(
  "similarity respects point-selected depth and time",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, depth = c(0, 20), time = c(0, 24))
    data <- data.frame(
      lon = 1,
      lat = 0,
      depth = 0,
      date = as.Date("2020-01-01"),
      target = 30
    )
    point_axes <- extract_to_buffer(
      data,
      path,
      var = "temp",
      horizontal_radius_km = 1,
      statistic = "most_similar",
      target_col = "target",
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_equal(point_axes$temp, 10)
    vertical <- extract_to_buffer(
      data,
      path,
      var = "temp",
      horizontal_radius_km = 1,
      depth_radius_m = 20,
      statistic = "most_similar",
      target_col = "target",
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_equal(vertical$temp, 30)
    data$target <- 24010
    same_time <- extract_to_buffer(
      data,
      path,
      var = "temp",
      horizontal_radius_km = 1,
      statistic = "most_similar",
      target_col = "target",
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_equal(same_time$temp, 10)
    window <- extract_to_buffer(
      data,
      path,
      var = "temp",
      horizontal_radius_km = 1,
      time_after = 1,
      time_unit = "days",
      statistic = "most_similar",
      target_col = "target",
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_equal(window$temp, 24010)
  }
)
