# extract_to_point(): a wrapper over as_vector3d(), as_proxy() and extract()
# (#116). The first tests cover its own behaviour: input forms, output naming,
# method and date checks. How values are picked (cells, depth levels, time
# steps, series, the antimeridian) is extract()'s, tested in
# test-proxy-extract.R. The rest test shared internals (see below). Synthetic
# files come from make_grid_file() (helper-netcdf-grid.R).

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
  make_grid_file( # nolint: object_usage_linter.
    c(0.5, 1.5, 2.5), c(10.5, 11.5), depth = depth, time = 0:4,
    values = vals, ...
  )
}
point_obs <- function() {
  data.frame(lon = c(0.2, 1.7, 2.9, 0.6, 1.1, 3.4),
             lat = c(10.2, 10.9, 11.8, 10.4, 10.1, 11),
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
  # A SpatProxy is used as it is.
  expect_equal(extract_to_point(obs, as_proxy(f))$temp, want)

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
  obs$temp <- 1
  expect_error(extract_to_point(obs, f, var = "temp", method = "bottom"),
               "already has a column", class = "ocean3d_error_point_input")
})

test_that("a single dated time step is matched like several", {
  # One file per day: the old extract_to_point() gave 1, 2 with the nearest
  # rule, 1, NA with "exact" or max_time_diff = 1, and needed a date.
  f <- make_grid_file(c(0.5, 1.5), c(10.5, 11.5), depth = NULL, time = 0,
                      values = 1:4)
  on.exit(unlink(f), add = TRUE)
  obs <- data.frame(lon = c(0.5, 1.5), lat = 10.5,
                    date = as.Date(c("2020-01-01", "2020-03-01")))
  expect_equal(extract_to_point(obs, f, var = "temp")$temp, c(1, 2))
  expect_equal(extract_to_point(obs, f, var = "temp",
                                time_match = "exact")$temp, c(1, NA))
  expect_equal(extract_to_point(obs, f, var = "temp",
                                max_time_diff = 1)$temp, c(1, NA))
  expect_error(extract_to_point(obs[c("lon", "lat")], f, var = "temp"),
               "has no column date", class = "ocean3d_error_point_input")
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
         pattern = "has no column tag"),
    # A matrix with a date column holds its coordinates as text.
    list(label = "text coordinates", nc = f, data = as.matrix(obs),
         args = list(method = "surface"),
         pattern = "Columns `lon`, `lat` must be numeric"),
    list(label = "sf without a CRS", nc = f,
         data = sf::st_as_sf(obs, coords = c("lon", "lat")),
         args = list(method = "surface"), pattern = "has no CRS")
  )
  for (case in cases) {
    expect_error(
      do.call(extract_to_point, c(list(data = case$data, nc = case$nc),
                                  case$args)),
      case$pattern, class = "ocean3d_error_point_input", label = case$label
    )
  }
  # An unreadable date is reported by its column, without advice about
  # arguments extract_to_point() doesn't have.
  bad <- obs
  bad$date <- c("2020-01-01", "2020-01-02 10:00", rep("2020-01-03", 4))
  err <- expect_error(extract_to_point(bad, f, method = "surface"),
                      "Input: column `date`",
                      class = "ocean3d_error_time_input")
  expect_no_match(conditionMessage(err), "time_resolution")
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

# The tests below call internals that the proxy, extract_to_buffer() and
# extract_to_track() still share: netCDF sources, horizontal matching, depth
# and time axes, bounds and report numbers. Tests that reach only the legacy
# point extraction are in test-extract-to-track.R.

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

test_that("depth methods name columns and extract values", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)
  cases <- list(
    list(method = "surface", col = "surface_temp", value = 13),
    list(method = "nearest", col = "nearest_temp", value = 17),
    list(method = "bottom", col = "seabottom_temp", value = 21)
  )
  x <- make_obs()
  for (case in cases) {
    out <- extract_to_point(x, f, var = NULL, method = case$method)
    expect_named(out, c(names(x), case$col))
    expect_equal(out[[case$col]], case$value)
    expect_equal(out[names(x)], x)
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
    expect_error(extract_to_point(
      data.frame(lon = 0, lat = 40), path, var = "temp", method = "surface",
      date_col = NULL
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
    expect_error(extract_to_point(
      data.frame(lon = 0, lat = 40), paths[i], var = "temp",
      method = "surface", date_col = NULL
    ), patterns[i])
  }
})

test_that("static variables use default date_col without a date column", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_test_netcdf(path, with_depth = FALSE, with_time = FALSE)
  x <- data.frame(id = c("a", "b"), lon = c(0, 1), lat = c(40, 41))
  out <- extract_to_point(x, path, var = "temp")
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
  expect_warning(out <- extract_to_point(x, path, var = "temp"),
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
  out <- extract_to_point(x, path, var = "temp")
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
  out <- extract_to_point(x, path, var = "temp", date_col = "observed_at")
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
    extract_to_point(data.frame(lon = 0, lat = 40), path, var = "temp",
                     method = "surface", date_col = NULL,
                     depth_dim = "depth"),
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
    extract_to_point(data.frame(lon = 0, lat = 40), path, var = "temp",
                     method = "surface", date_col = NULL),
    error = identity
  ))
  expect_match(conditionMessage(error), "has no coordinate variable")
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
    extract_to_point(x, path, var = "temp", date_col = NULL),
    class = "ocean3d_error_bounds"
  )
  # The message names the file, to find it in a multi-file call.
  expected <- paste0("Bounds for `lon` in ", basename(path), " were not found")
  expect_match(conditionMessage(err), expected, fixed = TRUE)
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
