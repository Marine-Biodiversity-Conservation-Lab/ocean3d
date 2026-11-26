# extract_to_point() and its extract2d()/extract3d_*() family: a wrapper over
# as_vector3d(), as_proxy() and extract() (#116). Its own behaviour is tested
# here: input forms, output naming, method and date checks. How values are
# picked (cells, depth levels, time steps, series, the antimeridian) is
# extract()'s, tested in test-proxy-extract.R. Synthetic files come from
# make_grid_file() (helper-netcdf-grid.R).

# lon 0.5, 1.5, 2.5 x lat 10.5, 11.5 x depth 5, 15, 40, 80 m x days 0-4, with
# the deepest level empty in the first cell and the 15 m level empty in the
# second. The expected values below are the old extract_to_point()'s, which
# the wrapper reproduces.
point_grid_values <- function() {
  set.seed(1)
  vals <- array(round(stats::runif(120) * 30, 2), c(3, 2, 4, 5))
  vals[1, 1, 4, ] <- NA
  vals[2, 1, 2, ] <- NA
  vals
}
point_grid_file <- function(depth = c(5, 15, 40, 80), ...) {
  vals <- point_grid_values()
  if (identical(depth, c(80, 40, 15, 5))) vals <- vals[, , 4:1, ]
  if (is.null(depth)) vals <- vals[, , 1, ]
  make_grid_file(c(0.5, 1.5, 2.5), c(10.5, 11.5), depth = depth, time = 0:4,
                 values = vals, ...)
}
point_obs <- function() {
  data.frame(lon = c(0.2, 1.7, 2.9, 0.6, 1.1, 3.4),
             lat = c(10.2, 10.9, 11.8, 10.4, 10.0, 11),
             depth = c(8, 14, 70, 75, 30, 10),
             date = as.Date("2020-01-01") + c(0, 1, 4, 2, 3, 1))
}

test_that("extract_to_point() picks each method's level at each time", {
  f <- point_grid_file()
  on.exit(unlink(f), add = TRUE)
  obs <- point_obs()
  # The 6th observation is off the grid. Surface keeps the shallowest level
  # even where it is empty; bottom and nearest skip empty levels.
  vals <- point_grid_values()
  expect_equal(extract_to_point(obs, f, var = "temp", method = "surface")$temp,
               c(vals[1, 1, 1, 1], vals[2, 1, 1, 2], vals[3, 2, 1, 5],
                 vals[1, 1, 1, 3], vals[2, 1, 1, 4], NA))
  expect_equal(extract_to_point(obs, f, var = "temp", method = "bottom")$temp,
               c(vals[1, 1, 3, 1], vals[2, 1, 4, 2], vals[3, 2, 4, 5],
                 vals[1, 1, 3, 3], vals[2, 1, 4, 4], NA))
  # Nearest to 14 m in the second cell skips the empty 15 m level.
  expect_equal(extract_to_point(obs, f, var = "temp", method = "nearest")$temp,
               c(vals[1, 1, 1, 1], vals[2, 1, 1, 2], vals[3, 2, 4, 5],
                 vals[1, 1, 3, 3], vals[2, 1, 3, 4], NA))

  # The same values whichever way the depth axis is stored.
  g <- point_grid_file(depth = c(80, 40, 15, 5))
  on.exit(unlink(g), add = TRUE)
  for (m in c("surface", "bottom", "nearest")) {
    expect_identical(extract_to_point(obs, g, var = "temp", method = m)$temp,
                     extract_to_point(obs, f, var = "temp", method = m)$temp,
                     label = m)
  }
})

test_that("extract_to_point() matches times as extract() does", {
  f <- point_grid_file()
  on.exit(unlink(f), add = TRUE)
  vals <- point_grid_values()
  one <- function(date, ...) {
    extract_to_point(nc = f, lon = 0.5, lat = 10.5, date = date,
                     var = "temp", method = "surface", ...)
  }
  # A date by UTC day; a date-time to the nearest step; exact needs the day.
  expect_equal(one(as.Date("2020-01-03")), vals[1, 1, 1, 3])
  expect_equal(one(as.POSIXct("2020-01-03 20:00", tz = "UTC")),
               vals[1, 1, 1, 4])
  expect_equal(one(as.Date("2020-01-09")), vals[1, 1, 1, 5])
  expect_true(is.na(one(as.Date("2020-01-09"), time_match = "exact")))
  expect_true(is.na(one(as.Date("2020-01-09"), max_time_diff = 2)))
  # A series of two files on one grid joins its time steps.
  s1 <- make_grid_file(c(0.5, 1.5, 2.5), c(10.5, 11.5), depth = NULL,
                       time = 0:2, values = vals[, , 1, 1:3])
  s2 <- make_grid_file(c(0.5, 1.5, 2.5), c(10.5, 11.5), depth = NULL,
                       time = 3:4, values = vals[, , 1, 4:5])
  on.exit(unlink(c(s1, s2)), add = TRUE)
  obs <- point_obs()
  expect_equal(extract_to_point(obs, c(s1, s2), var = "temp")$temp,
               extract_to_point(obs, list(s1, s2), var = "temp")$temp)
  expect_equal(extract_to_point(obs, c(s1, s2), var = "temp")$temp[3],
               vals[3, 2, 1, 5])
})

test_that("extract_to_point() takes tables, sf, matrices, lists and vectors", {
  f <- point_grid_file(depth = NULL)
  on.exit(unlink(f), add = TRUE)
  obs <- point_obs()
  want <- extract_to_point(obs, f, var = "temp")$temp

  pts <- sf::st_as_sf(obs, coords = c("lon", "lat"), crs = 4326)
  projected <- sf::st_transform(pts, 3857)
  cases <- list(
    list(label = "tibble", data = tibble::as_tibble(obs), class = "tbl_df"),
    list(label = "sf", data = pts, class = "sf"),
    list(label = "projected sf", data = projected, class = "sf"),
    list(label = "named list", data = as.list(obs), class = "data.frame")
  )
  for (case in cases) {
    out <- extract_to_point(case$data, f, var = "temp")
    expect_s3_class(out, case$class)
    expect_equal(out$temp, want, label = case$label)
  }
  mat <- as.matrix(obs[c("lon", "lat")])
  static <- make_grid_file(c(0.5, 1.5, 2.5), c(10.5, 11.5), depth = NULL,
                           time = NULL)
  on.exit(unlink(static), add = TRUE)
  expect_s3_class(extract_to_point(mat, static), "data.frame")

  # Direct input: one value per point, length-1 inputs recycled.
  expect_equal(extract_to_point(nc = f, lon = obs$lon, lat = obs$lat,
                                date = obs$date, var = "temp"), want)
  expect_equal(extract_to_point(nc = f, lon = c(0.5, 1.5), lat = 10.5,
                                date = as.Date("2020-01-01"), var = "temp"),
               point_grid_values()[1:2, 1, 1, 1])
  expect_error(extract_to_point(nc = f, lon = 1:3, lat = 1:2,
                                date = as.Date("2020-01-01")),
               "length 1 or the same length",
               class = "ocean3d_error_point_input")
})

test_that("observations without coordinates or depth get NA", {
  f <- point_grid_file()
  on.exit(unlink(f), add = TRUE)
  obs <- point_obs()[1:3, ]
  obs$lon[2] <- NA
  obs$depth[3] <- NA
  out <- extract_to_point(obs, f, var = "temp", method = "nearest")
  expect_true(is.na(out$temp[2]))
  expect_true(is.na(out$temp[3]))
  expect_false(is.na(out$temp[1]))
  # Depth is ignored by the other methods.
  expect_false(is.na(extract_to_point(obs, f, var = "temp",
                                      method = "bottom")$temp[3]))
  # A missing date gives NA, with extract()'s warning.
  obs$date[1] <- NA
  expect_warning(out <- extract_to_point(obs, f, var = "temp",
                                         method = "bottom"),
                 class = "ocean3d_warning_point_time")
  expect_true(is.na(out$temp[1]))
})

test_that("extract_to_point() names its output column", {
  f <- point_grid_file()
  on.exit(unlink(f), add = TRUE)
  obs <- point_obs()
  cases <- list(
    list(args = list(var = "temp", method = "nearest"), col = "temp"),
    list(args = list(method = "nearest"), col = "nearest_temp"),
    list(args = list(method = "surface"), col = "surface_temp"),
    list(args = list(method = "bottom"), col = "seabottom_temp"),
    list(args = list(method = "bottom", output_prefix = "b"), col = "b"),
    list(args = list(method = "bottom", output_prefix = "b",
                     output_col = "c"), col = "c")
  )
  for (case in cases) {
    out <- do.call(extract_to_point, c(list(obs, f), case$args))
    expect_true(case$col %in% names(out), label = case$col)
  }
  # The family: one method each.
  expect_identical(extract3d_surface(obs, f)$surface_temp,
                   extract_to_point(obs, f, method = "surface")$surface_temp)
  expect_identical(extract3d_bottom(obs, f, var = "temp")$temp,
                   extract_to_point(obs, f, var = "temp",
                                    method = "bottom")$temp)
  expect_identical(extract3d_nearest(obs, f)$nearest_temp,
                   extract_to_point(obs, f, method = "nearest")$nearest_temp)
  obs$temp <- 1
  expect_error(extract_to_point(obs, f, var = "temp", method = "bottom"),
               "already has a column", class = "ocean3d_error_point_input")
})

test_that("a single depth level reads with method = \"2d\"", {
  vals <- point_grid_values()
  f <- make_grid_file(c(0.5, 1.5, 2.5), c(10.5, 11.5), depth = 0.5,
                      time = 0:4, values = vals[, , 1, ])
  on.exit(unlink(f), add = TRUE)
  expect_equal(extract_to_point(nc = f, lon = 0.5, lat = 10.5,
                                date = as.Date("2020-01-02"), var = "temp"),
               vals[1, 1, 1, 2])
})

test_that("extract_to_point() rejects input it can't use", {
  f <- point_grid_file()
  f2 <- point_grid_file(depth = NULL)
  on.exit(unlink(c(f, f2)), add = TRUE)
  obs <- point_obs()
  cases <- list(
    list(label = "2d on several levels", nc = f, data = obs,
         args = list(method = "2d"), pattern = "several depth levels"),
    list(label = "a depth method on 2D", nc = f2, data = obs,
         args = list(method = "bottom"), pattern = "has no depth axis"),
    list(label = "a table of files", data = obs,
         nc = data.frame(file = f), args = list(),
         pattern = "table of files"),
    list(label = "no date column", nc = f, data = obs[c("lon", "lat")],
         args = list(method = "surface"), pattern = "has no column date"),
    list(label = "no date", nc = f, data = NULL,
         args = list(lon = 0.5, lat = 10.5, method = "surface"),
         pattern = "need a date"),
    list(label = "no depth column for nearest", nc = f,
         data = obs[c("lon", "lat", "date")],
         args = list(method = "nearest"), pattern = "has no column depth"),
    list(label = "missing id column", nc = f, data = obs,
         args = list(method = "surface", id_col = "tag"),
         pattern = "has no column tag")
  )
  for (case in cases) {
    expect_error(
      do.call(extract_to_point, c(list(data = case$data, nc = case$nc),
                                  case$args)),
      case$pattern, class = "ocean3d_error_point_input", label = case$label
    )
  }
  # Negative depths are rejected when the points are built.
  bad <- obs
  bad$depth[1] <- -5
  expect_error(extract_to_point(bad, f, method = "nearest"),
               "positive metres", class = "ocean3d_error_vector3d")
})

test_that("an open connection is read and left open", {
  f <- point_grid_file(depth = NULL)
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE, after = FALSE)
  out <- extract_to_point(point_obs(), nc, var = "temp")
  expect_equal(out$temp, extract_to_point(point_obs(), f, var = "temp")$temp)
  expect_no_error(ncdf4::ncatt_get(nc, 0))
})

# Shared internals ------------------------------------------------------------

# The tests below call internals that extract_to_buffer() and the legacy
# point extraction behind extract_to_track() still use (netCDF sources,
# horizontal matching, depth axes, time diagnostics and reports). The legacy
# ones go when extract_to_track() is rebuilt on sample_3d() and extract()
# (#116).

without_selection_reports <- function(expr) {
  quiet <- function(condition) invokeRestart("muffleMessage")
  withCallingHandlers(
    expr, ocean3d_time_report = quiet, ocean3d_spatial_report = quiet,
    ocean3d_depth_report = quiet
  )
}

make_test_netcdf <- function(path, var = "temp", with_depth = TRUE,
                             with_time = TRUE, bottom_na = FALSE,
                             time_units = "days since 2020-01-01",
                             lon_vals = c(0, 1)) {
  depth <- if (with_depth) c(0, 50, 100)
  time <- if (with_time) c(0, 2, 5)
  shape <- c(length(lon_vals), 2L, max(length(depth), 1L),
             max(length(time), 1L))
  values <- array(seq_len(prod(shape)), shape)
  if (with_depth && bottom_na) values[1, 1, 3, ] <- NA_real_
  make_grid_file( # nolint: object_usage_linter.
    lon_vals, c(40, 41), depth, time, time_units = time_units, var = var,
    var_units = "1", values = values, path = path
  )
}

make_obs <- function(date = as.Date("2020-01-03"), depth = 60) {
  data.frame(lon = 0, lat = 40, date = date, depth = depth)
}

capture_temporal_warning <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(expr, ocean3d_time_warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

make_b_diagnostic_input <- function(requested, selected, axis = selected,
                                    resolution = "datetime") {
  utc <- function(x) as.POSIXct(x, tz = "UTC")
  infer_cadence <- getFromNamespace(".infer_temporal_cadence", "ocean3d")
  diagnose <- getFromNamespace(".diagnose_temporal_selection", "ocean3d")
  metadata <- list(list(
    static = FALSE, variable = "temp", calendar = "standard",
    variable_units = "degree_Celsius", cell_methods = NULL,
    climatology = NULL, timestamps = utc(axis),
    cadence = infer_cadence(utc(axis))
  ))
  provenance <- list(temp = data.frame(
    source_index = rep(1L, length(selected)), variable = "temp",
    timestamp = utc(selected)
  ))
  diagnose(provenance, utc(requested), resolution, metadata)
}

make_horizontal_test_nc <- function(path, lon, lat, values, days = NULL) {
  make_grid_file( # nolint: object_usage_linter.
    lon, lat, depth = NULL, time = days, values = values,
    var_units = "degree_Celsius", path = path
  )
}

make_depth_axis_test_nc <- function(path, depths, values,
                                    units = "m", positive = NULL,
                                    depth_name = "depth",
                                    standard_name = NULL,
                                    formula_terms = NULL) {
  # One cell at lon 0, lat 40 with no time axis. NULL attributes are left out.
  make_grid_file( # nolint: object_usage_linter.
    lon = 0, lat = 40, depth = depths, time = NULL,
    values = values, depth_units = units,
    positive = positive, depth_name = depth_name,
    depth_atts = Filter(Negate(is.null), list(
      standard_name = standard_name, formula_terms = formula_terms
    )),
    var_units = "degree_Celsius", path = path
  )
}

test_that("nested netCDF lists keep every source, in order", {
  # a data frame inside a list used to replace the sources collected before it
  files_df <- data.frame(file = c("b.nc", "c.nc"))
  sources <- .as_netcdf_sources(list("a.nc", list(files_df, list("d.nc"))))

  expect_identical(vapply(sources, `[[`, character(1), "source"),
                   c("a.nc", "b.nc", "c.nc", "d.nc"))
})

test_that("source handling rejects unsupported inputs and missing paths", {
  expect_error(.as_netcdf_sources(42), "Unsupported")
  expect_error(.as_netcdf_sources(list(42)), "Unsupported")
  expect_error(
    .open_source(list(type = "path", source = tempfile())),
    "netCDF file does not exist",
    class = "ocean3d_error_netcdf_missing"
  )
})

test_that("depth wrappers name columns and extract values", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)
  cases <- list(
    list(fun = extract3d_surface, col = "surface_temp", value = 13),
    list(fun = extract3d_nearest, col = "nearest_temp", value = 17),
    list(fun = extract3d_bottom, col = "seabottom_temp", value = 21)
  )
  x <- make_obs()
  for (case in cases) {
    out <- case$fun(x, f, var = NULL)
    expect_named(out, c(names(x), case$col))
    expect_equal(out[[case$col]], case$value)
    expect_equal(out[names(x)], x)
  }
})

test_that("merging records provenance separately for each output column", {
  stamp <- as.POSIXct(c("2020-01-01", "2020-01-11"), tz = "UTC")
  a <- data.frame(id = 1:2, nearest_temp = c(1, NA), surface_temp = c(NA, 2))
  b <- data.frame(id = 1:2, nearest_temp = c(11, 12), surface_temp = c(13, 14))
  attr(a, ".source_candidate") <- list(
    distance_days = c(9, 1), timestamp = as.numeric(stamp),
    distance_km = c(0, 0), selected_lon = c(0, 0), selected_lat = c(40, 40),
    time_index = c(1L, 2L), variable = "temp"
  )
  attr(b, ".source_candidate") <- list(
    distance_days = c(1, 9), timestamp = as.numeric(rev(stamp)),
    distance_km = c(0, 0), selected_lon = c(0, 0), selected_lat = c(40, 40),
    time_index = c(2L, 1L), variable = "temp"
  )
  merged <- .merge_source_outputs(list(a, b), "id", return_provenance = TRUE)
  expect_equal(merged$data$nearest_temp, c(11, 12))
  expect_equal(merged$data$surface_temp, c(13, 2))
  expect_identical(merged$provenance$nearest_temp$source_index, c(2L, 2L))
  expect_identical(merged$provenance$surface_temp$source_index, c(2L, 1L))
  expect_equal(merged$provenance$nearest_temp$timestamp, rev(stamp))
  expect_identical(merged$provenance$nearest_temp$time_index, c(2L, 1L))
  expect_equal(merged$provenance$nearest_temp$distance_days, c(1, 9))
  expect_identical(merged$provenance$nearest_temp$variable, rep("temp", 2))
  expect_null(attr(merged$data, ".source_candidate", exact = TRUE))
  expect_equal(.merge_source_outputs(list(a, b), "id"), merged$data)
})

test_that("cadence inference distinguishes fixed and calendar sequences", {
  utc <- function(x) as.POSIXct(x, tz = "UTC")
  expect_identical(.infer_temporal_cadence(utc(c(
    "2020-01-01 00:00:00", "2020-01-01 03:00:00",
    "2020-01-01 06:00:00"
  )))$seconds, 10800)
  expect_identical(.infer_temporal_cadence(utc(c(
    "2020-01-01", "2020-01-02", "2020-01-03"
  )))$type, "daily")
  monthly <- .infer_temporal_cadence(utc(c(
    "2020-01-31", "2020-02-29", "2020-03-31"
  )))
  expect_identical(monthly$type, "monthly")
  expect_identical(monthly$months, 1L)
  expect_identical(.infer_temporal_cadence(utc(c(
    "2019-01-01", "2020-01-01", "2021-01-01"
  )))$type, "yearly")
  expect_identical(.infer_temporal_cadence(utc(c(
    "2020-01-01", "2020-01-03", "2020-01-08"
  )))$type, "irregular")
  expect_identical(.infer_temporal_cadence(utc("2020-01-01"))$type, "unknown")
})

test_that("calendar period changes warn even below 24 hours", {
  daily <- make_b_diagnostic_input(
    "2020-10-19 19:23:00", "2020-10-20 00:00:00", resolution = "daily"
  )
  expect_identical(daily$reason, "period_changed")
  expect_lt(abs(daily$difference_seconds), 86400)
  expect_warning(
    .warn_temporal_selection(daily), class = "ocean3d_time_warning"
  )
  expect_identical(make_b_diagnostic_input(
    "2020-01-31 23:00:00", "2020-02-01 00:00:00", resolution = "monthly"
  )$reason, "period_changed")
  expect_identical(make_b_diagnostic_input(
    "2020-12-31 23:00:00", "2021-01-01 00:00:00", resolution = "yearly"
  )$reason, "period_changed")
  expect_true(is.na(make_b_diagnostic_input(
    "2020-10-19 01:00:00", "2020-10-19 23:00:00", resolution = "daily"
  )$reason))
})

test_that("datetime diagnostics use cadence and calendar limits", {
  axis <- c("2020-01-01 00:00:00", "2020-01-01 03:00:00",
            "2020-01-01 06:00:00")
  expect_true(is.na(make_b_diagnostic_input(
    "2020-01-01 09:00:00", "2020-01-01 06:00:00", axis
  )$reason))
  expect_identical(make_b_diagnostic_input(
    "2020-01-01 09:00:01", "2020-01-01 06:00:00", axis
  )$reason, "cadence_exceeded")
  months <- c("2020-01-31", "2020-02-29", "2020-03-31")
  expect_true(is.na(make_b_diagnostic_input(
    "2020-02-29", "2020-01-31", months
  )$reason))
  expect_identical(make_b_diagnostic_input(
    "2020-03-01", "2020-01-31", months
  )$reason, "cadence_exceeded")
})

test_that("unknown and irregular cadences do not invent thresholds", {
  expect_identical(make_b_diagnostic_input(
    "2020-01-02", "2020-01-01"
  )$reason, "cadence_unknown")
  expect_identical(make_b_diagnostic_input(
    "2020-01-09", "2020-01-08", c("2020-01-01", "2020-01-03", "2020-01-08")
  )$reason, "cadence_unknown")
  expect_true(is.na(make_b_diagnostic_input(
    "2020-01-01", "2020-01-01"
  )$reason))
})

test_that("custom intervals require explicit positive durations", {
  interval <- as.difftime(3, units = "hours")
  expect_equal(
    .validate_custom_time_interval("custom", interval), 10800
  )
  expect_error(.validate_custom_time_interval("custom", NULL), "difftime")
  expect_error(.validate_custom_time_interval("custom", 3), "difftime")
  expect_error(
    .validate_custom_time_interval(
      "custom", as.difftime(0, units = "hours")
    ),
    "greater than zero"
  )
  expect_error(
    .validate_custom_time_interval(
      "custom", as.difftime(Inf, units = "hours")
    ),
    "finite"
  )
  expect_error(.validate_custom_time_interval("daily", interval), "requires")
})

test_that("daily datetime warnings explain day changes below 24 hours", {
  axis <- c("2020-01-01", "2020-01-02", "2020-01-03")
  result <- make_b_diagnostic_input(
    "2020-01-02 12:00:00", "2020-01-03 00:00:00", axis
  )
  expect_identical(result$reason, "period_changed")
  expect_false(result$cadence_exceeded)
  captured <- capture_temporal_warning(.warn_temporal_selection(result))
  text <- conditionMessage(captured$warnings[[1L]])
  expect_match(text, "12 hours later", fixed = TRUE)
  expect_match(text, "different day", fixed = TRUE)
  expect_match(text, "Comparison resolution: datetime", fixed = TRUE)
  expect_match(text, "Detected cadence: daily", fixed = TRUE)
  expect_match(text, "Row 1; variable: temp", fixed = TRUE)
  expect_match(text, "Observation: 2020-01-02 12:00:00 UTC", fixed = TRUE)
  expect_match(text, "Selected data: 2020-01-03 00:00:00 UTC", fixed = TRUE)
  late <- make_b_diagnostic_input(
    "2020-01-05 12:00:00", "2020-01-03 00:00:00", axis
  )
  expect_identical(late$reason, "period_changed")
  expect_true(late$cadence_exceeded)
  captured <- capture_temporal_warning(.warn_temporal_selection(late))
  text <- conditionMessage(captured$warnings[[1L]])
  expect_match(text, "60 hours earlier", fixed = TRUE)
  expect_match(text, "also more than one observed", fixed = TRUE)
  within_day <- make_b_diagnostic_input(
    "2020-01-02 12:00:00", "2020-01-02 00:00:00", axis
  )
  expect_true(is.na(within_day$reason))
})

test_that("unknown cadence warnings explain how to set an interval", {
  result <- make_b_diagnostic_input("2020-01-02", "2020-01-01")
  captured <- capture_temporal_warning(.warn_temporal_selection(result))
  text <- conditionMessage(captured$warnings[[1L]])
  expect_match(
    text, "do not establish a reliable sampling interval", fixed = TRUE
  )
  expect_match(text, "custom_time_interval", fixed = TRUE)
})

test_that("nearest spatial matching requires an explicit valid distance", {
  expect_identical(.validate_spatial_matching("within", NULL), "within")
  expect_identical(.validate_spatial_matching("nearest", 0), "nearest")
  expect_identical(.validate_spatial_matching("nearest", 5), "nearest")
  expect_error(.validate_spatial_matching("within", 5), "requires")
  for (limit in list(NULL, NA_real_, Inf, -1, c(1, 2), "5")) {
    expect_error(
      .validate_spatial_matching("nearest", limit),
      "finite"
    )
  }
})

test_that("formula-based and pressure axes require explicit conversion", {
  paths <- replicate(3, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_depth_axis_test_nc(
    paths[1], c(0, 0.5, 1), c(11, 22, 33), units = "1",
    standard_name = "ocean_sigma_coordinate"
  )
  make_depth_axis_test_nc(
    paths[2], c(0, 50, 100), c(11, 22, 33),
    formula_terms = "sigma: sigma eta: eta depth: bathymetry"
  )
  make_depth_axis_test_nc(
    paths[3], c(0, 50, 100), c(11, 22, 33), units = "dbar",
    standard_name = "sea_water_pressure"
  )
  for (path in paths) {
    expect_error(extract3d_surface(
      data.frame(lon = 0, lat = 40), path, var = "temp", date_col = NULL
    ), "physical transformation")
  }
})

test_that("invalid file depth metadata fail before vertical extraction", {
  paths <- replicate(4, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_depth_axis_test_nc(paths[1], c(-50, 0, 50), c(11, 22, 33))
  make_depth_axis_test_nc(
    paths[2], c(0, -50, -100), c(11, 22, 33), positive = "down"
  )
  make_depth_axis_test_nc(paths[3], c(0, 50, 100), c(11, 22, 33), units = "")
  make_depth_axis_test_nc(paths[4], c(0, 50, 50), c(11, 22, 33))
  patterns <- c("ambiguous", "reference surface", "units", "duplicate")
  for (i in seq_along(paths)) {
    expect_error(extract3d_surface(
      data.frame(lon = 0, lat = 40), paths[i], var = "temp", date_col = NULL
    ), patterns[i])
  }
})

test_that("report numbers are short without hiding small differences", {
  expect_equal(.selection_number(1.80399482019311e-7, 5L), "0")
  expect_equal(.selection_number(-1e-7, 5L), "0")
  expect_equal(.selection_number(38.0208320617676, 5L), "38.02083")
  expect_equal(.selection_number(1.01823663711548), "1.02")
  expect_equal(.selection_number(20), "20")
  expect_equal(.selection_number(0), "0")
  expect_equal(.selection_number(0.001, small = TRUE), "<0.01")
  expect_equal(.selection_number(0.01, small = TRUE), "0.01")
})

test_that("local spacing uses surrounding layers instead of a global mean", {
  reference <- .depth_selection_reference(c(150, 67, 55, 0), 60, 150)
  expect_equal(c(reference$lower, reference$upper, reference$spacing),
               c(55, 67, 12))
  expect_identical(reference$basis, "surrounding")
  expect_identical(reference$status, "distant_depth")
  expect_true(reference$closer_missing)
  exact <- .depth_selection_reference(c(0, 55, 67, 150), 55, 150)
  expect_identical(exact$basis, "adjacent")
  expect_equal(c(exact$lower, exact$upper, exact$spacing), c(55, 67, 12))
  tied <- .depth_selection_reference(c(0, 50, 100), 50, 0)
  expect_equal(tied$spacing, 50)
  expect_identical(tied$status, "normal")
  expect_equal(.depth_selection_reference(c(0, 50, 100), 0, 100)$spacing, 50)
  expect_equal(.depth_selection_reference(c(0, 50, 100), 100, 0)$spacing, 50)
})

test_that("spacing outside the range uses edge layers and a strict threshold", {
  for (requested in c(-100, 200)) {
    chosen <- if (requested < 0) 0 else 100
    reference <- .depth_selection_reference(c(100, 0, 50), requested, chosen)
    expect_true(reference$outside)
    expect_equal(reference$spacing, 50)
    expect_identical(reference$basis, "nearest_edge")
    expect_identical(reference$status, "distant_depth")
    expect_false(reference$closer_missing)
    text <- .format_depth_reference(reference)
    expect_match(text, "outside the layer range", fixed = TRUE)
    expect_false(grepl("no valid value", text, fixed = TRUE))
  }
  boundary <- .depth_selection_reference(c(0, 50, 100), 150, 100)
  expect_identical(boundary$status, "normal")
  above <- .depth_selection_reference(c(0, 50, 100), 150.001, 100)
  expect_identical(above$status, "distant_depth")
  singleton <- .depth_selection_reference(50, 60, 50)
  expect_identical(singleton$status, "spacing_unknown")
  expect_true(is.na(singleton$spacing))
  expect_identical(.depth_selection_reference(50, 50, 50)$status, "normal")
})

test_that("comparison precision distinguishes close values", {
  expect_identical(.selection_comparison(1.114, 1.111),
                   list(value = "1.114", limit = "1.111"))
  expect_identical(.selection_comparison(1.111, 1.114),
                   list(value = "1.111", limit = "1.114"))
  expect_identical(.selection_comparison(1.111, 1.111),
                   list(value = "1.11", limit = "1.11"))
  expect_identical(.selection_comparison(12.634, 20),
                   list(value = "12.63", limit = "20"))
  tiny <- .selection_comparison(0.00002, 0.00001)
  expect_identical(tiny, list(value = "0.00002", limit = "0.00001"))
  close <- .selection_comparison(1 + .Machine$double.eps, 1)
  expect_false(identical(close$value, close$limit))
})

test_that("static variables use default date_col without a date column", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_test_netcdf(path, with_depth = FALSE, with_time = FALSE)
  x <- data.frame(id = c("a", "b"), lon = c(0, 1), lat = c(40, 41))
  out <- extract2d(x, path, var = "temp")
  expect_named(out, c("id", "lon", "lat", "temp"))
  expect_equal(out$temp, c(1, 4))
  expect_equal(out[names(x)], x)
  expect_equal(extract_to_point(
    nc = path, lon = 0, lat = 40, var = "temp", method = "2d"
  ), 1)
})

test_that("an existing typed NA date preserves observations & missing values", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_test_netcdf(path, with_depth = FALSE)
  x <- data.frame(id = c("a", "b"), lon = c(0, 1), lat = c(40, 41),
                  date = as.Date(c("2020-01-03", NA_character_)))
  # A missing date gives NA, with extract()'s warning.
  expect_warning(out <- extract2d(x, path, var = "temp"),
                 "has no time", class = "ocean3d_warning_point_time")
  expect_named(out, c("id", "lon", "lat", "date", "temp"))
  expect_equal(out[names(x)], x)
  expect_equal(out$temp, c(5, NA_real_))
  expect_warning(value <- extract_to_point(
    nc = path, lon = 0, lat = 40, date = as.Date(NA_character_),
    var = "temp", method = "2d"
  ), "has no time", class = "ocean3d_warning_point_time")
  expect_length(value, 1L)
  expect_true(is.na(value))
})

test_that("static inputs do not parse unused date values", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_test_netcdf(path, with_depth = FALSE, with_time = FALSE)
  x <- data.frame(lon = 0, lat = 40, date = "not-a-date")
  out <- extract2d(x, path, var = "temp")
  expect_named(out, c("lon", "lat", "date", "temp"))
  expect_equal(out$temp, 1)
  expect_equal(out$date, "not-a-date")
  expect_equal(extract_to_point(
    nc = path, lon = 0, lat = 40, date = "not-a-date",
    var = "temp", method = "2d"
  ), 1)
})

test_that("custom date columns and scalar date recycling keep known values", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_test_netcdf(path, with_depth = FALSE)
  x <- data.frame(lon = c(0, 1), lat = c(40, 41),
                  observed_at = as.Date(c("2020-01-01", "2020-01-03")))
  out <- extract2d(x, path, var = "temp", date_col = "observed_at")
  expect_named(out, c("lon", "lat", "observed_at", "temp"))
  expect_equal(out$observed_at, x$observed_at)
  expect_equal(out$temp, c(1, 8))
  expect_equal(extract_to_point(
    nc = path, lon = c(0, 1), lat = c(40, 41),
    date = as.Date("2020-01-03"), var = "temp", method = "2d"
  ), c(5, 8))
})

test_that("explicit zones and date-only inputs need no timezone warning", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_test_netcdf(path, with_depth = FALSE, with_time = FALSE)
  expect_no_warning(value <- extract_to_point(
    nc = path, lon = 0, lat = 40, date = "2020-01-01 12:00",
    var = "temp", method = "2d"
  ))
  expect_equal(value, 1)
})

test_that("height axes without a positive attribute are not depths", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(2, 10, 20), c(11, 22, 33),
                          standard_name = "height")
  expect_error(
    extract3d_surface(data.frame(lon = 0, lat = 40), path, var = "temp",
                      date_col = NULL, depth_dim = "depth"),
    "reference surface"
  )
})

test_that("a depth dimension without coordinate variable fails clearly", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  lon <- ncdf4::ncdim_def("lon", "degrees_east", 0)
  lat <- ncdf4::ncdim_def("lat", "degrees_north", 40)
  depth <- ncdf4::ncdim_def("depth", "", 1:3, create_dimvar = FALSE)
  variable <- ncdf4::ncvar_def("temp", "degree_Celsius",
                               list(lon, lat, depth), -9999)
  nc <- ncdf4::nc_create(path, variable)
  ncdf4::ncvar_put(nc, "temp", array(c(11, 22, 33), c(1, 1, 3)))
  ncdf4::nc_close(nc)
  printed <- capture.output(error <- tryCatch(
    extract3d_surface(data.frame(lon = 0, lat = 40), path, var = "temp",
                      date_col = NULL),
    error = identity
  ))
  expect_match(conditionMessage(error), "has no coordinate variable")
  expect_match(conditionMessage(error), "extract2d()", fixed = TRUE)
  expect_length(printed, 0L)
})

test_that("unusable spatial bounds stop extraction", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_horizontal_test_nc(path, c(0, 1), c(40, 41), 1:4)
  nc <- ncdf4::nc_open(path, write = TRUE)
  ncdf4::ncatt_put(nc, "lon", "bounds", "lon_bnds")
  ncdf4::nc_close(nc)
  x <- data.frame(lon = c(0, 1, 5), lat = 40)
  err <- expect_error(
    extract2d(x, path, var = "temp", date_col = NULL),
    class = "ocean3d_error_bounds"
  )
  # The message names the file, to find it in a multi-file call.
  expected <- paste0("Bounds for `lon` in ", basename(path), " were not found")
  expect_match(conditionMessage(err), expected, fixed = TRUE)
})

test_that("depth references tolerate float32 axes and blame missing layers", {
  layers <- c(0, 50, 100)
  far <- .depth_selection_reference(layers, 200, 50)
  expect_identical(far$status, "distant_depth")
  expect_false(far$closer_missing)
  near <- .depth_selection_reference(layers, 120, 50)
  expect_true(near$closer_missing)
  float32 <- c(0.49402499198913574, 1.54137504100799561, 3.81949496269226074)
  exact <- .depth_selection_reference(float32, 1.541375, float32[2])
  expect_identical(exact$basis, "adjacent")
  expect_identical(exact$status, "normal")
  deepest <- .depth_selection_reference(float32, 3.819495, float32[2])
  expect_false(deepest$outside)
})

test_that("bad spatial bounds stop point extraction and the core object", {
  check <- function(lon_bounds, pattern) {
    f <- make_grid_file(depth = NULL, time = NULL, lon_bounds = lon_bounds)
    on.exit(unlink(f))
    err <- expect_error(
      extract_to_point(nc = f, lon = 1.5, lat = 10.5, var = "temp",
                       method = "2d"),
      class = "ocean3d_error_bounds"
    )
    expect_match(conditionMessage(err), pattern)
    expect_error(as_proxy(f), class = "ocean3d_error_bounds")
  }
  # Centres 0.5, 1.5, 2.5: the first cell's bounds miss its centre.
  check(rbind(c(1, 1, 2), c(1.2, 2, 3)), "must contain their centres")
  # The first two cells overlap between 0.8 and 1.
  check(rbind(c(0, 0.8, 2), c(1, 2, 3)), "Bounds for `lon` in .* overlap\\.")
})

test_that("bounds names inside a netCDF-4 group are relative to the group", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  lon <- ncdf4::ncdim_def("grp/lon", "degrees_east", c(0.5, 1.5))
  lat <- ncdf4::ncdim_def("grp/lat", "degrees_north", c(10.5, 11.5))
  nv <- ncdf4::ncdim_def("grp/nv", "", 1:2, create_dimvar = FALSE)
  nc <- ncdf4::nc_create(f, list(
    ncdf4::ncvar_def("grp/temp", "degC", list(lon, lat), -9999),
    ncdf4::ncvar_def("grp/lon_bnds", "degrees_east", list(nv, lon))
  ), force_v4 = TRUE)
  ncdf4::ncvar_put(nc, "grp/temp", matrix(1:4, 2))
  ncdf4::ncvar_put(nc, "grp/lon_bnds", rbind(c(0, 1), c(1, 2)))
  ncdf4::ncatt_put(nc, "grp/lon", "bounds", "lon_bnds")
  ncdf4::nc_close(nc)

  expect_equal(
    extract_to_point(nc = f, var = "grp/temp", lon = 1.5, lat = 10.5,
                     method = "2d"),
    2
  )
  expect_identical(as_proxy(f, var = "grp/temp")@axes$lon$edges, "bounds")
})

