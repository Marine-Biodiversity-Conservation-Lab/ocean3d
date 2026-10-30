# Silence only known informational reports in value/selection tests.
# Report assertions keep their own captures; other conditions remain visible.
without_selection_reports <- function(expr) {
  quiet <- function(condition) invokeRestart("muffleMessage")
  withCallingHandlers(
    expr, ocean3d_time_report = quiet, ocean3d_spatial_report = quiet,
    ocean3d_depth_report = quiet
  )
}

# Verify the expected B diagnostic while preserving the extraction value.
# Other warning classes remain visible to testthat.
expect_period_change_warning <- function(expr) {
  captured <- list()
  value <- withCallingHandlers(expr, ocean3d_time_warning = function(w) {
    captured[[length(captured) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  testthat::expect_length(captured, 1L)
  for (w in captured) {
    testthat::expect_s3_class(w, "ocean3d_time_warning")
    testthat::expect_gt(nrow(w$diagnostics), 0L)
    testthat::expect_true(all(w$diagnostics$reason == "period_changed"))
  }
  value
}

# Helpers ----------------------------------------------------------------------

make_test_netcdf <- function(path, var = "temp", with_depth = TRUE,
                             with_time = TRUE, bottom_na = FALSE,
                             time_units = "days since 2020-01-01",
                             lon_vals = c(0, 1), coordinate_var = FALSE) {
  lon <- ncdf4::ncdim_def("lon", "degrees_east", vals = lon_vals)
  lat <- ncdf4::ncdim_def("lat", "degrees_north", vals = c(40, 41))
  dims <- list(lon, lat)

  if (with_depth) {
    depth <- ncdf4::ncdim_def("depth", "m", vals = c(0, 50, 100))
    dims <- c(dims, list(depth))
  }

  if (with_time) {
    time <- ncdf4::ncdim_def(
      "time", if (is.null(time_units)) "" else time_units,
      vals = if (is.null(time_units)) 1:3 else c(0, 2, 5),
      unlim = FALSE, create_dimvar = !is.null(time_units)
    )
    dims <- c(dims, list(time))
  }

  v <- ncdf4::ncvar_def(var, "1", dim = dims, missval = -9999, prec = "double")
  vars <- list(v)
  if (coordinate_var) {
    vars <- c(vars, list(ncdf4::ncvar_def(
      "longitude", "degrees_east", list(lon), prec = "double"
    )))
  }
  nc <- ncdf4::nc_create(path, vars)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  shape <- vapply(dims, function(x) length(x$vals), integer(1))
  values <- array(seq_len(prod(shape)), dim = shape)

  if (with_depth && bottom_na) {
    if (with_time) values[1, 1, 3, ] <- NA_real_
    else values[1, 1, 3] <- NA_real_
  }

  if (coordinate_var) ncdf4::ncvar_put(nc, "longitude", lon_vals + 10)
  ncdf4::ncvar_put(nc, var, values)
  invisible(path)
}

make_obs <- function(date = as.Date("2020-01-03"), depth = 60) {
  data.frame(lon = 0, lat = 40, date = date, depth = depth)
}

# Generic extraction -----------------------------------------------------------

test_that("extract_to_point extracts from a 2D netCDF variable", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)

  out <- extract_to_point(
    data = data.frame(lon = 0, lat = 40),
    nc = f, var = "temp", method = "2d",
    date_col = NULL,
    time_report = FALSE
  )

  expect_true("temp" %in% names(out))
  expect_identical(out$temp, 1)
})

test_that("extract_to_point detects a single variable automatically", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, var = "oxygen")

  out <- extract3d_surface(make_obs(), f, var = NULL,
    time_report = FALSE
  )

  expect_true("surface_oxygen" %in% names(out))
})

test_that("extract_to_point() errors when var = NULL matches two variables", {
  # automatic variable detection is ambiguous with more than one candidate
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  v1 <- ncdf4::ncvar_def("temp", "1", list(lon, lat), -9999)
  v2 <- ncdf4::ncvar_def("salt", "1", list(lon, lat), -9999)
  nc <- ncdf4::nc_create(f, list(v1, v2))
  ncdf4::ncvar_put(nc, "temp", matrix(1:4, 2, 2))
  ncdf4::ncvar_put(nc, "salt", matrix(5:8, 2, 2))
  ncdf4::nc_close(nc)

  expect_error(
    extract_to_point(
      data.frame(lon = 0, lat = 40),
      f, var = NULL, method = "2d", date_col = NULL,
      time_report = FALSE
    ),
    "exactly one variable"
  )
})


# Depth methods ----------------------------------------------------------------

test_that("3D wrappers return surface, nearest, bottom and all summaries", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, bottom_na = TRUE)
  x <- make_obs()

  surface <- extract3d_surface(x, f, var = "temp",
    time_report = FALSE
  )
  nearest <- extract3d_nearest(x, f, var = "temp",
    time_report = FALSE
  )
  bottom <- extract3d_bottom(x, f, var = "temp",
    time_report = FALSE
  )
  all <- extract3d_all(x, f, var = "temp",
    time_report = FALSE
  )

  expect_identical(surface$temp, all$surface_temp)
  expect_identical(nearest$temp, all$nearest_temp)
  expect_identical(bottom$temp, all$seabottom_temp)
  expect_false(is.na(bottom$temp))
})

test_that("nearest depth ignores missing layers below the valid water column", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, bottom_na = TRUE)

  out <- extract3d_nearest(make_obs(depth = 100), f, var = "temp",
    time_report = FALSE
  )

  expect_false(is.na(out$temp))
})

# Time matching ----------------------------------------------------------------

test_that("time_match = nearest uses the closest available netCDF date", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  out <- expect_period_change_warning(extract3d_surface(
    make_obs(date = as.Date("2020-01-02")),
    f, var = "temp", time_match = "nearest",
    time_report = FALSE
  ))

  expect_false(is.na(out$temp))
})

test_that("time_match = exact returns NA when the date is absent", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  out <- extract3d_surface(
    make_obs(date = as.Date("2020-01-02")),
    f, var = "temp", time_match = "exact",
    time_report = FALSE
  )

  expect_true(is.na(out$temp))
})

test_that("max_time_diff limits nearest temporal matching", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  accepted <- expect_period_change_warning(extract3d_surface(
    make_obs(date = as.Date("2020-01-04")),
    f, var = "temp", time_match = "nearest", max_time_diff = 1,
    time_report = FALSE
  ))
  rejected <- extract3d_surface(
    make_obs(date = as.Date("2020-01-10")),
    f, var = "temp", time_match = "nearest", max_time_diff = 1,
    time_report = FALSE
  )

  expect_false(is.na(accepted$temp))
  expect_true(is.na(rejected$temp))
})

test_that("max_time_diff rejects invalid values", {
  x <- data.frame(lon = 0, lat = 40)
  for (bad in list(-1, "1", c(1, 2), numeric(), NA_real_)) {
    expect_error(
      extract_to_point(
        x, "unused.nc", var = "temp", method = "2d",
        date_col = NULL, max_time_diff = bad,
        time_report = FALSE
      ),
      "non-negative numeric value"
    )
  }
})

test_that("multiple netCDF files add one output column per detected variable", {
  skip_if_not_installed("ncdf4")
  f1 <- tempfile(fileext = ".nc")
  f2 <- tempfile(fileext = ".nc")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  make_test_netcdf(f1, var = "temp")
  make_test_netcdf(f2, var = "salt")

  out <- extract3d_bottom(make_obs(), c(f1, f2), var = NULL,
    time_report = FALSE
  )

  expect_true(all(c("seabottom_temp", "seabottom_salt") %in% names(out)))
})

test_that("netCDF input can be a list, opened connection, or data frame", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  from_list <- extract3d_surface(make_obs(), list(f), var = "temp",
    time_report = FALSE
  )

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  from_connection <- extract3d_surface(make_obs(), nc, var = "temp",
    time_report = FALSE
  )

  from_df <- extract3d_surface(
    make_obs(), data.frame(file = f), var = "temp", file_col = "file",
    time_report = FALSE
  )

  expect_identical(from_list$temp, from_connection$temp)
  expect_identical(from_list$temp, from_df$temp)
})

test_that("nested netCDF lists keep every source, in order", {
  # a data frame inside a list used to replace the sources collected before it
  files_df <- data.frame(file = c("b.nc", "c.nc"))
  sources <- .as_netcdf_sources(list("a.nc", list(files_df, list("d.nc"))))

  expect_identical(vapply(sources, `[[`, character(1), "source"),
                   c("a.nc", "b.nc", "c.nc", "d.nc"))
})

test_that("extract_to_point reports missing required observation columns", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  expect_error(
    extract3d_nearest(
      data.frame(lon = 0, lat = 40, date = as.Date("2020-01-03")),
      f, var = "temp",
      time_report = FALSE
    ),
    "Missing required column"
  )
})

# Direct point extraction ------------------------------------------------------

test_that("extract_to_point extracts a single value from coordinates", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  out <- extract_to_point(
    nc = f,
    lon = 0,
    lat = 40,
    depth = 60,
    var = "temp",
    method = "nearest",
    date_col = NULL,
    time_report = FALSE
  )

  expect_length(out, 1)
  expect_identical(as.numeric(out), 5)
})

test_that("extract_to_point requires a valid longitude", {
  expect_error(
    extract_to_point(
      nc = "unused.nc",
      lon = NULL,
      lat = 40,
      method = "2d",
      date_col = NULL,
      time_report = FALSE
    ),
    "`lon` must contain numeric values"
  )
})

test_that("extract_to_point requires a valid latitude", {
  expect_error(
    extract_to_point(
      nc = "unused.nc",
      lon = 0,
      lat = NULL,
      method = "2d",
      date_col = NULL,
      time_report = FALSE
    ),
    "`lat` must contain numeric values"
  )
})

test_that("extract_to_point requires depth for nearest extraction", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  expect_error(
    extract_to_point(
      nc = f,
      lon = 0,
      lat = 40,
      method = "nearest",
      date_col = NULL,
      time_report = FALSE
    ),
    "`depth` must contain numeric values"
  )
})

test_that("extract_to_point rejects invalid dates", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  expect_error(
    extract_to_point(
      nc = f,
      lon = 0,
      lat = 40,
      depth = 60,
      date = "not-a-date",
      var = "temp",
      method = "nearest",
      time_report = FALSE
    ),
    "Time values could not be converted"
  )
})

test_that("extract_to_point extracts 2d values without depth", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)

  out <- extract_to_point(
    nc = f,
    lon = 0,
    lat = 40,
    var = "temp",
    method = "2d",
    date_col = NULL,
    time_report = FALSE
  )

  expect_length(out, 1)
  expect_identical(as.numeric(out), 1)
})


# Input formats ---------------------------------------------------------------
test_that("extract_to_point accepts a matrix", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  x <- matrix(
    c(0, 40, 60),
    nrow = 1,
    dimnames = list(NULL, c("lon", "lat", "depth"))
  )

  out <- extract_to_point(
    data = x,
    nc = f,
    var = "temp",
    method = "nearest",
    date_col = NULL,
    time_report = FALSE
  )

  expect_identical(as.numeric(out$temp), 5)
})

test_that("extract_to_point accepts a named list", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  x <- list(
    lon = 0,
    lat = 40,
    depth = 60
  )

  out <- extract_to_point(
    data = x,
    nc = f,
    var = "temp",
    method = "nearest",
    date_col = NULL,
    time_report = FALSE
  )

  expect_identical(as.numeric(out$temp), 5)
})

test_that("extract_to_point accepts sf point geometries", {
  skip_if_not_installed("ncdf4")
  skip_if_not_installed("sf")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  x <- sf::st_sf(
    depth = 60,
    geometry = sf::st_sfc(sf::st_point(c(0, 40)), crs = 4326)
  )

  out <- extract_to_point(
    data = x,
    nc = f,
    var = "temp",
    method = "nearest",
    date_col = NULL,
    time_report = FALSE
  )

  expect_identical(as.numeric(out$temp), 5)
})

test_that("extract_to_point transforms projected sf coordinates", {
  skip_if_not_installed("ncdf4")
  skip_if_not_installed("sf")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  x <- sf::st_sf(
    depth = 60,
    geometry = sf::st_sfc(sf::st_point(c(0, 40)), crs = 4326)
  )
  x <- sf::st_transform(x, 3857)

  out <- extract_to_point(
    data = x,
    nc = f,
    var = "temp",
    method = "nearest",
    date_col = NULL,
    time_report = FALSE
  )

  expect_s3_class(out, "sf")
  expect_identical(as.numeric(out$temp), 5)
})

test_that("extract_to_point accepts tibbles", {
  skip_if_not_installed("ncdf4")
  skip_if_not_installed("tibble")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  x <- tibble::tibble(
    lon = 0,
    lat = 40,
    depth = 60
  )

  out <- extract_to_point(
    data = x,
    nc = f,
    var = "temp",
    method = "nearest",
    date_col = NULL,
    time_report = FALSE
  )

  expect_s3_class(out, "tbl_df")
  expect_identical(as.numeric(out$temp), 5)
})

test_that("extract_to_point accepts coordinate vectors", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  out <- extract_to_point(
    nc = f,
    lon = c(0, 1),
    lat = c(40, 41),
    depth = c(60, 100),
    var = "temp",
    method = "nearest",
    time_report = FALSE
  )

  expect_length(out, 2)
  expect_identical(as.numeric(out), c(5, 12))
})

test_that("direct coordinate vectors must have compatible lengths", {
  expect_error(
    extract_to_point(
      nc = "unused.nc",
      lon = c(0, 1),
      lat = 40,
      depth = c(50, 100, 150),
      method = "nearest",
      time_report = FALSE
    ),
    "compatible lengths"
  )
})

test_that("extract_to_point preserves data frame input behaviour", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  x <- data.frame(
    lon = c(0, 1),
    lat = c(40, 41),
    depth = c(60, 100)
  )

  out <- extract_to_point(
    data = x,
    nc = f,
    var = "temp",
    method = "nearest",
    date_col = NULL,
    time_report = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 2L)
  expect_identical(as.numeric(out$temp), c(5, 12))
})

test_that("direct inputs recycle scalar values", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)

  out <- extract_to_point(
    nc = f,
    lon = c(0, 1),
    lat = 40,
    depth = c(60, 100),
    var = "temp",
    method = "nearest",
    time_report = FALSE
  )

  expect_length(out, 2)
  expect_identical(as.numeric(out), c(5, 10))
})

# Input validation -------------------------------------------------------------

test_that("extract_to_point rejects sf geometries other than POINT", {
  skip_if_not_installed("sf")

  x <- sf::st_sf(
    depth = 60,
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 40, 1, 41), ncol = 2, byrow = TRUE)),
      crs = 4326
    )
  )

  expect_error(
    extract_to_point(
      data = x,
      nc = "unused.nc",
      method = "nearest",
      date_col = NULL,
      time_report = FALSE
    ),
    "POINT geometries"
  )
})

test_that("extract_to_point rejects projected sf input without a CRS", {
  skip_if_not_installed("sf")

  x <- sf::st_sf(
    depth = 60,
    geometry = sf::st_sfc(sf::st_point(c(0, 40)))
  )

  expect_error(
    extract_to_point(
      data = x,
      nc = "unused.nc",
      method = "nearest",
      date_col = NULL,
      time_report = FALSE
    ),
    "geographic CRS"
  )
})

test_that("extract_to_point rejects matrices without column names", {
  x <- matrix(c(0, 40, 60), nrow = 1)

  expect_error(
    extract_to_point(
      data = x,
      nc = "unused.nc",
      method = "nearest",
      date_col = NULL,
      time_report = FALSE
    ),
    "Matrix input must have column names"
  )
})

test_that("extract_to_point rejects unnamed lists", {
  x <- list(0, 40, 60)

  expect_error(
    extract_to_point(
      data = x,
      nc = "unused.nc",
      method = "nearest",
      date_col = NULL,
      time_report = FALSE
    ),
    "List input must be named"
  )
})

test_that("extract_to_point rejects unsupported data input types", {
  expect_error(
    extract_to_point(
      data = factor("a"),
      nc = "unused.nc",
      method = "2d",
      date_col = NULL,
      time_report = FALSE
    ),
    "`data` must be a data frame"
  )
})


# Public entry points and netCDF helpers ---------------------------------------

test_that("extract2d preserves observations and extracts known cells", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  x <- data.frame(id = c("a", "b"), lon = c(0, 1), lat = c(40, 41))
  out <- extract2d(x, f, var = "temp", date_col = NULL,
    time_report = FALSE
  )
  expect_equal(out[names(x)], x)
  expect_equal(out$temp, c(1, 4))
})

test_that("variable and output column names must be single strings", {
  x <- data.frame(lon = 0, lat = 40)
  for (bad in list(42, c("a", "b"), character(), NA_character_)) {
    expect_error(extract_to_point(
      x, "unused.nc", var = bad, method = "2d", date_col = NULL,
      time_report = FALSE
    ), "var.*must be")
    expect_error(extract_to_point(
      x, "unused.nc", var = "temp", output_col = bad,
      method = "2d", date_col = NULL,
      time_report = FALSE
    ), "output_col.*must be")
  }
})

test_that("depth methods reject files without a depth axis", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  for (method in c("surface", "bottom", "nearest", "all")) {
    expect_error(extract_to_point(
      make_obs(), f, var = "temp", method = method, date_col = NULL,
      time_report = FALSE
    ), "requires a depth dimension")
  }
})

test_that("a single direct observation returns named depth summaries", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)
  out <- extract_to_point(nc = f, lon = 0, lat = 40, depth = 60,
    var = "temp", method = "all",
    time_report = FALSE
  )
  expect_equal(out, c(nearest_temp = 5, surface_temp = 1, seabottom_temp = 9))
})

test_that("verbose extraction reports progress without changing results", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  x <- data.frame(lon = rep(0, 100), lat = rep(40, 100))
  expect_message(
    out <- extract2d(x, f, var = "temp", date_col = NULL, verbose = TRUE,
      time_report = FALSE
    ),
    "Extracting row 100 of 100"
  )
  expect_equal(out$temp, rep(1, 100))
  expect_silent(extract2d(x, f, var = "temp", date_col = NULL,
    time_report = FALSE
  ))
})

test_that("source handling rejects unsupported inputs and missing paths", {
  expect_error(.as_netcdf_sources(42), "Unsupported")
  expect_error(.as_netcdf_sources(list(42)), "Unsupported")
  expect_error(
    .open_source(list(type = "path", source = tempfile())),
    "netCDF file does not exist"
  )
})

test_that("extraction leaves caller-owned netCDF connections open", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  opened <- .open_source(list(type = "connection", source = nc))
  expect_identical(opened$nc, nc)
  expect_false(opened$close)
  out <- extract2d(data.frame(lon = 0, lat = 40), nc,
    var = "temp", date_col = NULL,
    time_report = FALSE
  )
  expect_equal(out$temp, 1)
  expect_equal(ncdf4::ncvar_get(nc, "temp"), matrix(1:4, 2, 2))
})

test_that("dimension helpers read coordinates and report missing names", {
  expect_error(
    .match_name(c("x", "y"), "lon"),
    "Available names are: x, y"
  )
  expect_null(.match_name("x", "lon", required = FALSE))
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, coordinate_var = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_equal(as.vector(.get_dim_values(nc, "longitude")), c(10, 11))
  expect_null(.get_dim_values(nc, NULL))
  expect_error(.get_dim_values(nc, "missing"),
               "Could not read dimension or coordinate variable")
  expect_error(.detect_dims(nc, "missing"),
               "Available variables are: temp, longitude")

})

test_that("longitude conversion and missing coordinates are handled", {
  without_selection_reports({
    f <- tempfile(fileext = ".nc")
    on.exit(unlink(f), add = TRUE)
    make_test_netcdf(f, with_depth = FALSE, with_time = FALSE,
                     lon_vals = c(0, 359))
    x <- data.frame(lon = c(-1, NA, 0), lat = c(40, 40, NA))
    out <- extract2d(x, f, var = "temp", date_col = NULL,
      time_report = FALSE
    )
    expect_equal(out$temp, c(2, NA, NA))
    expect_equal(nrow(out), 3L)
  })
})


test_that("time without usable units gives an informative error", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  # Keep the time coordinate variable, but give it empty units.
  make_test_netcdf(f, time_units = "")

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_error(
    .decode_netcdf_time(nc, "time"),
    "no usable.*units"
  )
})

test_that("output prefixes respect automatic and explicit naming", {
  for (method in c("surface", "bottom", "nearest")) {
    prefix <- switch(method, surface = "surface_", bottom = "seabottom_",
                     nearest = "nearest_")
    expect_equal(
      .resolve_output_prefix(method, "temp", NULL, TRUE),
      paste0(prefix, "temp")
    )
    expect_null(.resolve_output_prefix(method, "temp", NULL, FALSE))
    expect_equal(
      .resolve_output_prefix(method, "temp", "custom", TRUE),
      "custom"
    )
  }
})

test_that("multiple files fill missing values & preserve first-file priority", {
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  for (f in paths) make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  nc <- ncdf4::nc_open(paths[1], write = TRUE)
  tryCatch(
    ncdf4::ncvar_put(nc, "temp", matrix(c(NA, 2, 3, 4), 2, 2)),
    finally = ncdf4::nc_close(nc)
  )
  nc <- ncdf4::nc_open(paths[2], write = TRUE)
  tryCatch(
    ncdf4::ncvar_put(nc, "temp", matrix(11:14, 2, 2)),
    finally = ncdf4::nc_close(nc)
  )
  x <- data.frame(lon = c(0, 1), lat = c(40, 41))
  out <- extract2d(x, paths, var = "temp", date_col = NULL,
    time_report = FALSE
  )
  expect_equal(out$temp, c(11, 4))
  reverse <- extract2d(x, rev(paths), var = "temp", date_col = NULL,
    time_report = FALSE
  )
  expect_equal(reverse$temp, c(11, 14))
  expect_equal(out[names(x)], x)
})

test_that("entirely missing profiles return missing depth summaries", {
  without_selection_reports({
    f <- tempfile(fileext = ".nc")
    on.exit(unlink(f), add = TRUE)
    make_test_netcdf(f, with_time = FALSE)
    nc <- ncdf4::nc_open(f, write = TRUE)
    tryCatch(
      ncdf4::ncvar_put(nc, "temp", rep(NA_real_, 3),
                       start = c(1, 1, 1), count = c(1, 1, 3)),
      finally = ncdf4::nc_close(nc)
    )
    for (method in c("nearest", "bottom", "all")) {
      out <- extract_to_point(make_obs(), f, var = "temp", method = method,
        date_col = NULL,
        time_report = FALSE
      )
      result_cols <- if (method == "all") {
        c("nearest_temp", "surface_temp", "seabottom_temp")
      } else {
        "temp"
      }
      expect_named(out, c(names(make_obs()), result_cols))
      expect_true(all(is.na(out[result_cols])))
      expect_equal(nrow(out), 1L)
    }
  })
})

test_that("missing target depth preserves surface and bottom summaries", {
  without_selection_reports({
    f <- tempfile(fileext = ".nc")
    on.exit(unlink(f), add = TRUE)
    make_test_netcdf(f, with_time = FALSE)
    x <- make_obs(depth = NA_real_)
    nearest <- extract3d_nearest(x, f, var = "temp", date_col = NULL,
      time_report = FALSE
    )
    all <- extract3d_all(x, f, var = "temp", date_col = NULL,
      time_report = FALSE
    )
    expect_true(is.na(nearest$temp))
    expect_true(is.na(all$nearest_temp))
    expect_equal(all$surface_temp, 1)
    expect_equal(all$seabottom_temp, 9)
  })
})

test_that("output_col overrides output_prefix in public extraction", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  x <- data.frame(lon = 0, lat = 40)
  for (var in list(NULL, "temp")) {
    out <- extract2d(x, f, var = var, date_col = NULL,
      output_prefix = "ignored", output_col = "chosen",
      time_report = FALSE
    )
    expect_named(out, c("lon", "lat", "chosen"))
    expect_equal(out$chosen, 1)
  }
})

test_that("nested sources retain paths and caller-owned connections", {
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  for (f in paths) {
    make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  }

  first <- ncdf4::nc_open(paths[1], write = TRUE)
  tryCatch(
    ncdf4::ncvar_put(
      first, "temp", matrix(c(NA, 2, 3, 4), 2, 2)
    ),
    finally = ncdf4::nc_close(first)
  )
  second <- ncdf4::nc_open(paths[2], write = TRUE)
  tryCatch(
    ncdf4::ncvar_put(second, "temp", matrix(11:14, 2, 2)),
    finally = ncdf4::nc_close(second)
  )

  nc <- ncdf4::nc_open(paths[2])
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  x <- data.frame(lon = c(0, 1), lat = c(40, 41))
  nested <- list(
    year = list(path = paths[1], opened = list(connection = nc))
  )
  out <- extract2d(x, nested, var = "temp", date_col = NULL,
    time_report = FALSE
  )

  expect_equal(out$temp, c(11, 4))
  expect_equal(out[names(x)], x)
  expect_equal(ncdf4::ncvar_get(nc, "temp"), matrix(11:14, 2, 2))
})

test_that("direct vectors return all depth summaries", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)
  out <- extract_to_point(
    nc = f, lon = c(0, 1), lat = c(40, 41), depth = c(60, 100),
    var = "temp", method = "all",
    time_report = FALSE
  )
  expect_equal(out, data.frame(
    nearest_temp = c(5, 12), surface_temp = c(1, 4),
    seabottom_temp = c(9, 12)
  ))
})


test_that("id_col must exist and is preserved when supplied", {
  x <- data.frame(lon = 0, lat = 40)
  expect_error(extract2d(
    x, "unused.nc", var = "temp", date_col = NULL, id_col = "sample",
    time_report = FALSE
  ), "Missing required column.*sample")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  x$sample <- "station-a"
  out <- extract2d(x, f, var = "temp", date_col = NULL, id_col = "sample",
    time_report = FALSE
  )
  expect_equal(out$sample, "station-a")
  expect_equal(out$temp, 1)
})
test_that("direct coordinates and a valid date select the correct time slice", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)
  for (date in list("2020-01-03", as.Date("2020-01-03"))) {
    out <- extract_to_point(
      nc = f, lon = 0, lat = 40, depth = 60, date = date,
      var = "temp", method = "nearest", time_match = "exact",
      time_report = FALSE
    )
    expect_equal(out, 17)
  }
})

test_that("exact matching selects values and skips missing dates", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)
  x <- data.frame(
    lon = c(0, 0, 0), lat = c(40, 40, 40),
    date = as.Date(c("2020-01-01", "2020-01-03", NA))
  )
  exact <- extract3d_surface(x, f, var = "temp", time_match = "exact",
    time_report = FALSE
  )
  nearest <- extract3d_surface(x, f, var = "temp", time_match = "nearest",
    time_report = FALSE
  )
  expect_equal(exact$temp, c(1, 13, NA))
  expect_equal(nearest$temp, c(1, 13, NA))
  expect_equal(exact[names(x)], x)
})

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

test_that("extract2d automatically names detected variable with time match", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, var = "oxygen", with_depth = FALSE)
  x <- data.frame(lon = c(0, 1), lat = c(40, 41),
                  date = as.Date(c("2020-01-01", "2020-01-03")))
  out <- extract2d(x, f, time_match = "exact",
    time_report = FALSE
  )
  expect_named(out, c(names(x), "oxygen"))
  expect_equal(out$oxygen, c(1, 8))
})

test_that("nested data frames use a custom file column", {
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_test_netcdf(paths[1], var = "temp")
  make_test_netcdf(paths[2], var = "salt")
  sources <- list(year = list(files = data.frame(nc_path = paths)))
  x <- make_obs()
  out <- extract3d_surface(x, sources, file_col = "nc_path",
    time_report = FALSE
  )
  expect_equal(out[names(x)], x)
  expect_equal(out$surface_temp, 13)
  expect_equal(out$surface_salt, 13)
  expect_equal(out, extract3d_surface(x, paths,
    time_report = FALSE
  ))
  expect_error(
    extract3d_surface(x, list(data.frame(wrong = paths)),
      file_col = "nc_path",
      time_report = FALSE
    ),
    "Missing required column.*nc_path"
  )
})

make_time_decoder_fixture <- function(
  path,
  values,
  units = "hours since 2020-01-01 00:00:00",
  calendar = "standard") {
  time <- ncdf4::ncdim_def(
    "time",
    units,
    vals = values
  )

  variable <- ncdf4::ncvar_def(
    "example",
    "1",
    dim = list(time),
    prec = "double"
  )

  nc <- ncdf4::nc_create(path, variable)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  ncdf4::ncatt_put(nc, "time", "calendar", calendar)
  ncdf4::ncvar_put(nc, "example", seq_along(values))

  invisible(path)
}


test_that("time decoding preserves noon across multiple days", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  make_time_decoder_fixture(f, values = c(12, 36, 132))

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expected <- as.POSIXct(
    c(
      "2020-01-01 12:00:00",
      "2020-01-02 12:00:00",
      "2020-01-06 12:00:00"
    ),
    tz = "UTC"
  )

  expect_equal(.decode_netcdf_time(nc, "time"), expected)
})


test_that("time decoding preserves the hour in the reference origin", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  make_time_decoder_fixture(
    f,
    values = c(0, 1, 24),
    units = "hours since 2020-01-01 06:30:00"
  )

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expected <- as.POSIXct(
    c(
      "2020-01-01 06:30:00",
      "2020-01-01 07:30:00",
      "2020-01-02 06:30:00"
    ),
    tz = "UTC"
  )

  expect_equal(.decode_netcdf_time(nc, "time"), expected)
})


test_that("supported time units decode to independently known instants", {
  cases <- list(
    list(unit = "milliseconds", values = c(0, 3600000)),
    list(unit = "second", values = c(0, 3600)),
    list(unit = "minutes", values = c(0, 60)),
    list(unit = "hour", values = c(0, 1)),
    list(unit = "days", values = c(0, 1 / 24))
  )

  paths <- character()
  on.exit(unlink(paths), add = TRUE)

  expected <- as.POSIXct(
    c("2020-01-01 00:00:00", "2020-01-01 01:00:00"),
    tz = "UTC"
  )

  for (case in cases) {
    f <- tempfile(fileext = ".nc")
    paths <- c(paths, f)

    make_time_decoder_fixture(
      f,
      values = case$values,
      units = paste(case$unit, "since 2020-01-01")
    )

    nc <- ncdf4::nc_open(f)

    tryCatch(
      expect_equal(.decode_netcdf_time(nc, "time"), expected),
      finally = ncdf4::nc_close(nc)
    )
  }
})


test_that("time decoding rejects empty units instead of inventing an origin", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  make_time_decoder_fixture(
    f,
    values = c(0, 2, 5),
    units = ""
  )

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_error(
    .decode_netcdf_time(nc, "time"),
    "no usable.*units"
  )
})


test_that("time decoding rejects unsupported calendars", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  make_time_decoder_fixture(
    f,
    values = c(0, 24),
    calendar = "360_day"
  )

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_error(
    .decode_netcdf_time(nc, "time"),
    "unsupported calendar"
  )
})

# Timestamp-label comparison precision -----------------------------------------

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

test_that("auto preserves observation precision", {
  daily <- .prepare_time_input("2020-01-02")
  hourly <- .prepare_time_input("2020-01-02 19:54:00")

  expect_identical(daily$resolution, "daily")
  expect_identical(hourly$resolution, "datetime")

  expect_error(
    .prepare_time_input(
      c("2020-01-02", "2020-01-02 19:54:00")
    ),
    "mixes dates"
  )
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
    "Multiple dataset timestamps"
  )
  for (resolution in c("monthly", "yearly")) {
    expect_error(
      .match_time_index(
        available, as.Date("2020-01-02"),
        time_resolution = resolution
      ),
      "Multiple dataset timestamps"
    )
  }
})

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

test_that("public extraction preserves hours and supports daily selection", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  make_test_netcdf(
    f,
    with_depth = FALSE,
    time_units = "hours since 2020-01-01"
  )

  nc <- ncdf4::nc_open(f, write = TRUE)
  tryCatch(
    ncdf4::ncvar_put(nc, "time", c(12, 36, 60)),
    finally = ncdf4::nc_close(nc)
  )

  target <- as.POSIXct("2020-01-02 19:54:00", tz = "UTC")

  daily <- extract_to_point(
    nc = f, lon = 0, lat = 40, date = target,
    var = "temp", method = "2d",
    time_resolution = "daily",
    time_match = "exact",
    time_report = FALSE
  )
  expect_equal(daily, 5)

  table_input <- extract2d(
    data.frame(lon = 0, lat = 40, date = target),
    f, var = "temp",
    time_resolution = "daily",
    time_match = "exact",
    time_report = FALSE
  )
  expect_equal(table_input$temp, 5)

  nc <- ncdf4::nc_open(f, write = TRUE)
  tryCatch(
    ncdf4::ncvar_put(nc, "time", c(18, 19, 20)),
    finally = ncdf4::nc_close(nc)
  )

  hourly <- extract_to_point(
    nc = f, lon = 0, lat = 40,
    date = as.POSIXct("2020-01-01 19:54:00", tz = "UTC"),
    var = "temp", method = "2d",
    time_report = FALSE
  )
  expect_equal(hourly, 9)

  expect_error(
    extract_to_point(
      nc = f, lon = 0, lat = 40, date = "2020-01-01",
      var = "temp", method = "2d",
      time_report = FALSE
    ),
    "Multiple dataset timestamps"
  )
})

test_that("noon-stamped dates select correct slices with zero tolerance", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(
    f, with_depth = FALSE,
    time_units = "hours since 2020-01-01"
  )

  nc <- ncdf4::nc_open(f, write = TRUE)
  tryCatch(
    ncdf4::ncvar_put(nc, "time", c(12, 36, 132)),
    finally = ncdf4::nc_close(nc)
  )

  x <- data.frame(
    lon = rep(0, 3), lat = rep(40, 3),
    date = as.Date(c("2020-01-01", "2020-01-02", "2020-01-06"))
  )
  for (matching in c("nearest", "exact")) {
    out <- extract2d(
      x, f, var = "temp", time_match = matching,
      max_time_diff = 0,
      time_report = FALSE
    )
    expect_equal(out$temp, c(1, 5, 9))
    expect_equal(out[names(x)], x)
  }
})

test_that("every depth wrapper forwards output naming arguments", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)
  x <- make_obs()
  wrappers <- list(
    surface = extract3d_surface,
    nearest = extract3d_nearest,
    bottom = extract3d_bottom,
    all = extract3d_all
  )

  for (method in names(wrappers)) {
    out <- wrappers[[method]](
      x, f, var = "temp", output_col = "chosen"
    )
    if (method == "all") {
      expected_cols <- c(
        "nearest_chosen", "surface_chosen", "seabottom_chosen"
      )
      expect_named(out, c(names(x), expected_cols))
      expect_equal(
        unname(unlist(out[expected_cols])),
        c(17, 13, 21)
      )
    } else {
      expect_named(out, c(names(x), "chosen"))
      expected <- switch(method, surface = 13, nearest = 17, bottom = 21)
      expect_equal(out$chosen, expected)
    }
  }
})


# Issue #72: variable-scoped netCDF dimensions -------------------------------

test_that("dimension detection uses only the selected variable", {
  skip_if_not_installed("ncdf4")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  other_lon <- ncdf4::ncdim_def(
    "longitude", "degrees_east", c(10, 11)
  )
  time <- ncdf4::ncdim_def(
    "time", "days since 2020-01-01", c(0, 1)
  )

  static <- ncdf4::ncvar_def(
    "static", "1", list(lon, lat), -9999
  )
  dynamic <- ncdf4::ncvar_def(
    "dynamic", "1", list(other_lon, lat, time), -9999
  )

  nc <- ncdf4::nc_create(f, list(static, dynamic))
  ncdf4::ncvar_put(nc, "static", matrix(1:4, 2, 2))
  ncdf4::ncvar_put(
    nc, "dynamic", array(1:8, c(2, 2, 2))
  )
  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  dims <- .detect_dims(nc, "static")

  expect_identical(dims$lon, "lon")
  expect_identical(dims$lat, "lat")
  expect_null(dims$time)
  expect_identical(dims$var_dims, c("lon", "lat"))
})


test_that("valid_time is recognised as a time dimension", {
  skip_if_not_installed("ncdf4")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  valid_time <- ncdf4::ncdim_def(
    "valid_time", "days since 2020-01-01", c(0, 1)
  )

  v <- ncdf4::ncvar_def(
    "temp", "1", list(lon, lat, valid_time), -9999
  )

  nc <- ncdf4::nc_create(f, v)
  ncdf4::ncvar_put(nc, "temp", array(1:8, c(2, 2, 2)))
  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  dims <- .detect_dims(nc, "temp")

  expect_identical(dims$time, "valid_time")
})


test_that("unknown non-singleton dimensions are rejected", {
  skip_if_not_installed("ncdf4")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  member <- ncdf4::ncdim_def("member", "1", 1:3)

  v <- ncdf4::ncvar_def(
    "temp", "1", list(lon, lat, member), -9999
  )

  nc <- ncdf4::nc_create(f, v)
  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_error(
    .detect_dims(nc, "temp"),
    "Unrecognised non-singleton dimension"
  )
})


test_that("unknown singleton dimensions are identified", {
  skip_if_not_installed("ncdf4")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  member <- ncdf4::ncdim_def("member", "1", 1)

  v <- ncdf4::ncvar_def(
    "temp", "1", list(lon, lat, member), -9999
  )

  nc <- ncdf4::nc_create(f, v)
  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  dims <- .detect_dims(nc, "temp")

  expect_identical(dims$singleton_dims, "member")

  sc <- .build_start_count(
    dims$var_dims,
    list(
      lon = list(start = 1L, count = 1L),
      lat = list(start = 1L, count = 1L)
    ),
    singleton_dims = dims$singleton_dims
  )

  expect_identical(sc$start, c(1L, 1L, 1L))
  expect_identical(sc$count, c(1L, 1L, 1L))
})


test_that("explicit dimensions are validated", {
  skip_if_not_installed("ncdf4")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  make_test_netcdf(
    f, with_depth = FALSE, with_time = FALSE
  )

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_identical(
    .detect_dims(nc, "temp", lon_dim = "lon")$lon,
    "lon"
  )

  expect_error(
    .detect_dims(nc, "temp", lon_dim = "longitude"),
    "not a dimension of the selected variable"
  )
})


test_that("projected coordinates are rejected", {
  skip_if_not_installed("ncdf4")

  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  x <- ncdf4::ncdim_def("x", "m", c(1000, 2000))
  y <- ncdf4::ncdim_def("y", "m", c(3000, 4000))

  v <- ncdf4::ncvar_def(
    "temp", "1", list(x, y), -9999
  )

  nc <- ncdf4::nc_create(f, v)
  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_error(
    .detect_dims(nc, "temp"),
    "projected coordinates"
  )
})


test_that("2d extraction rejects variables with multiple depth levels", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  make_test_netcdf(
    f,
    with_depth = TRUE,
    with_time = FALSE
  )

  x <- data.frame(lon = 0, lat = 40)

  expect_error(
    extract2d(x, f, var = "temp", date_col = NULL,
      time_report = FALSE
    ),
    "multiple depth levels"
  )
})


test_that("2d extraction accepts a single depth level", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  depth <- ncdf4::ncdim_def("depth", "m", 5)

  v <- ncdf4::ncvar_def(
    "temp", "degree_Celsius",
    list(lon, lat, depth), -9999
  )

  nc <- ncdf4::nc_create(f, v)
  ncdf4::ncvar_put(
    nc, "temp",
    array(c(11, 12, 13, 14), dim = c(2, 2, 1))
  )
  ncdf4::nc_close(nc)

  x <- data.frame(lon = 0, lat = 40)

  out <- extract2d(
    x, f, var = "temp", date_col = NULL,
    time_report = FALSE
  )

  expect_equal(out$temp, 11)
})


test_that("unassigned dimensions cannot be read silently", {
  expect_error(
    .build_start_count(
      c("lon", "lat", "member"),
      list(
        lon = list(start = 1L, count = 1L),
        lat = list(start = 1L, count = 1L)
      )
    ),
    "No extraction index was assigned"
  )
})


# additional dimension-detection regression tests ----------------

test_that("CF metadata identifies axes with non-standard names", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  dim_x <- ncdf4::ncdim_def(
    "horizontal_a", "degrees_east", c(0, 1)
  )
  dim_y <- ncdf4::ncdim_def(
    "horizontal_b", "degrees_north", c(40, 41)
  )
  dim_t <- ncdf4::ncdim_def(
    "forecast_step", "days since 2020-01-01", c(0, 1)
  )

  v <- ncdf4::ncvar_def(
    "temp", "degree_Celsius",
    list(dim_x, dim_y, dim_t), -9999
  )

  nc <- ncdf4::nc_create(f, v)
  ncdf4::ncvar_put(
    nc, "temp", array(1:8, dim = c(2, 2, 2))
  )
  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  dims <- .detect_dims(nc, "temp")

  expect_identical(dims$lon, "horizontal_a")
  expect_identical(dims$lat, "horizontal_b")
  expect_identical(dims$time, "forecast_step")
})


test_that("ambiguous longitude dimensions produce an informative error", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def(
    "lon", "degrees_east", c(0, 1)
  )
  longitude <- ncdf4::ncdim_def(
    "longitude", "degrees_east", c(10, 11)
  )
  lat <- ncdf4::ncdim_def(
    "lat", "degrees_north", c(40, 41)
  )

  v <- ncdf4::ncvar_def(
    "temp", "degree_Celsius",
    list(lon, longitude, lat), -9999
  )

  nc <- ncdf4::nc_create(f, v)
  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_error(
    .detect_dims(nc, "temp"),
    "Ambiguous lon axis"
  )
})


test_that("2D auxiliary geographic coordinates are rejected", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  x <- ncdf4::ncdim_def("x", "1", c(1, 2))
  y <- ncdf4::ncdim_def("y", "1", c(1, 2))

  temp <- ncdf4::ncvar_def(
    "temp", "degree_Celsius",
    list(x, y), -9999
  )
  longitude <- ncdf4::ncvar_def(
    "longitude", "degrees_east",
    list(x, y), -9999
  )
  latitude <- ncdf4::ncvar_def(
    "latitude", "degrees_north",
    list(x, y), -9999
  )

  nc <- ncdf4::nc_create(
    f, list(temp, longitude, latitude)
  )

  ncdf4::ncatt_put(
    nc, "temp", "coordinates", "longitude latitude"
  )
  ncdf4::ncatt_put(
    nc, "longitude", "standard_name", "longitude"
  )
  ncdf4::ncatt_put(
    nc, "latitude", "standard_name", "latitude"
  )

  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  expect_error(
    .detect_dims(nc, "temp"),
    "Curvilinear longitude/latitude grids are not"
  )
})


test_that("static extraction ignores unrelated temporal dimensions", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def(
    "lon", "degrees_east", c(0, 1)
  )
  lat <- ncdf4::ncdim_def(
    "lat", "degrees_north", c(40, 41)
  )
  time <- ncdf4::ncdim_def(
    "time", "days since 2020-01-01", c(0, 1)
  )

  static <- ncdf4::ncvar_def(
    "static", "1", list(lon, lat), -9999
  )
  dynamic <- ncdf4::ncvar_def(
    "dynamic", "1", list(lon, lat, time), -9999
  )

  nc <- ncdf4::nc_create(f, list(static, dynamic))

  ncdf4::ncvar_put(
    nc, "static", matrix(c(11, 12, 13, 14), 2, 2)
  )
  ncdf4::ncvar_put(
    nc, "dynamic", array(1:8, c(2, 2, 2))
  )

  ncdf4::nc_close(nc)

  x <- data.frame(
    lon = c(0, 1),
    lat = c(40, 41)
  )

  out <- extract2d(
    x, f, var = "static", date_col = NULL,
    time_report = FALSE
  )

  expect_equal(out$static, c(11, 14))
  expect_equal(nrow(out), 2L)
})


# -------------------------------------------------------------------------
# Issue #73: Temporal source eligibility
# -------------------------------------------------------------------------

test_that("temporal candidates are not restricted to datetime ranges", {
  timestamps <- as.POSIXct(
    c("2019-01-01 00:00:00", "2019-12-31 00:00:00"),
    tz = "UTC"
  )

  expect_true(
    .temporal_source_eligible(
      timestamps,
      as.POSIXct("2019-06-01 12:00:00", tz = "UTC"),
      "datetime"
    )
  )

  expect_true(
    .temporal_source_eligible(
      timestamps,
      as.POSIXct("2020-07-01 00:00:00", tz = "UTC"),
      "datetime"
    )
  )

  expect_true(
    .temporal_source_eligible(
      timestamps,
      as.POSIXct("2019-12-31 00:00:00", tz = "UTC"),
      "datetime"
    )
  )
})


test_that("temporal candidates may lie outside the requested day", {
  timestamps <- as.POSIXct(
    c("2020-01-01 12:00:00", "2020-01-03 12:00:00"),
    tz = "UTC"
  )

  expect_true(
    .temporal_source_eligible(
      timestamps,
      as.POSIXct("2020-01-03 19:54:00", tz = "UTC"),
      "daily"
    )
  )

  expect_true(
    .temporal_source_eligible(
      timestamps,
      as.POSIXct("2020-01-04 00:01:00", tz = "UTC"),
      "daily"
    )
  )
})


test_that("single-timestamp sources can offer nearest candidates", {
  timestamp <- as.POSIXct(
    "2020-07-01 12:00:00",
    tz = "UTC"
  )

  expect_true(
    .temporal_source_eligible(
      timestamp,
      as.POSIXct("2020-07-01 19:54:00", tz = "UTC"),
      "daily"
    )
  )

  expect_true(
    .temporal_source_eligible(
      timestamp,
      as.POSIXct("2020-07-02 00:00:00", tz = "UTC"),
      "daily"
    )
  )

  expect_true(
    .temporal_source_eligible(
      timestamp,
      as.POSIXct("2020-07-01 19:54:00", tz = "UTC"),
      "datetime"
    )
  )
})


test_that("static sources have no temporal restriction", {
  expect_true(
    .temporal_source_eligible(
      NULL,
      as.POSIXct("2020-07-01", tz = "UTC"),
      "daily"
    )
  )
})


test_that("temporal eligibility handles irregular and missing inputs", {
  timestamps <- as.POSIXct(
    c("2020-01-10", "2020-01-01", "2020-01-04"),
    tz = "UTC"
  )

  expect_true(
    .temporal_source_eligible(
      timestamps,
      as.POSIXct("2020-01-05", tz = "UTC"),
      "daily"
    )
  )

  expect_false(
    .temporal_source_eligible(
      timestamps,
      as.POSIXct(NA_character_, tz = "UTC"),
      "daily"
    )
  )

  expect_false(
    .temporal_source_eligible(
      as.POSIXct(character(0), tz = "UTC"),
      as.POSIXct("2020-01-05", tz = "UTC"),
      "daily"
    )
  )
})

test_that("multi-source extraction selects the correct temporal file", {

  make_temporal_nc <- function(path, year, value) {
    lon <- ncdf4::ncdim_def("lon", "degrees_east", vals = 0)
    lat <- ncdf4::ncdim_def("lat", "degrees_north", vals = 40)

    time <- ncdf4::ncdim_def(
      "time",
      paste0("days since ", year, "-01-01 00:00:00"),
      vals = c(0, 364)
    )

    variable <- ncdf4::ncvar_def(
      "temp", "degree_Celsius",
      list(lon, lat, time),
      missval = -9999,
      prec = "double"
    )

    nc <- ncdf4::nc_create(path, variable)
    on.exit(ncdf4::nc_close(nc), add = TRUE)

    ncdf4::ncvar_put(
      nc, "temp",
      array(value, dim = c(1, 1, 2))
    )
  }

  file_2019 <- tempfile(fileext = ".nc")
  file_2020 <- tempfile(fileext = ".nc")

  on.exit(
    unlink(c(file_2019, file_2020)),
    add = TRUE
  )

  make_temporal_nc(file_2019, 2019, 1931)
  make_temporal_nc(file_2020, 2020, 2007)

  observations <- data.frame(
    id = c("obs_2019", "obs_2020"),
    lon = c(0, 0),
    lat = c(40, 40),
    date = as.Date(c("2019-07-01", "2020-07-01"))
  )

  forward <- expect_period_change_warning(extract2d(
    data = observations,
    nc = c(file_2019, file_2020),
    var = "temp",
    time_match = "nearest",
    time_report = FALSE
  ))

  reverse <- expect_period_change_warning(extract2d(
    data = observations,
    nc = c(file_2020, file_2019),
    var = "temp",
    time_match = "nearest",
    time_report = FALSE
  ))

  expect_equal(forward$temp, c(1931, 2007))
  expect_equal(reverse$temp, c(1931, 2007))

  expect_identical(forward$id, observations$id)
  expect_identical(reverse$id, observations$id)
  expect_identical(forward$date, observations$date)
  expect_identical(reverse$date, observations$date)

  # A single source retains the previous nearest behaviour.
  single <- expect_period_change_warning(extract2d(
    data = observations[2, , drop = FALSE],
    nc = file_2019,
    var = "temp",
    time_match = "nearest",
    time_report = FALSE
  ))

  expect_equal(single$temp, 1931)
})


# Issue #73 A2: compare temporal candidates across sources -----------------

make_a2_temporal_nc <- function(path, days, values, with_depth = FALSE) {
  lon <- ncdf4::ncdim_def("lon", "degrees_east", 0)
  lat <- ncdf4::ncdim_def("lat", "degrees_north", 40)
  time <- ncdf4::ncdim_def("time", "days since 2020-01-01", days)
  dims <- list(lon, lat)
  if (with_depth) {
    dims <- c(dims, list(ncdf4::ncdim_def("depth", "m", c(0, 50, 100))))
  }
  dims <- c(dims, list(time))
  variable <- ncdf4::ncvar_def(
    "temp", "degree_Celsius", dims, -9999, prec = "double"
  )
  nc <- ncdf4::nc_create(path, variable)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  shape <- vapply(dims, function(x) length(x$vals), integer(1))
  ncdf4::ncvar_put(nc, "temp", array(values, dim = shape))
}

test_that("nearest compares singleton and multi-record files in either order", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_a2_temporal_nc(paths[1], c(0, 4), c(1, 5))
  make_a2_temporal_nc(paths[2], 10, 11)
  x <- data.frame(
    id = c("gap", "outside", "missing"), lon = 0, lat = 40,
    date = as.Date(c("2020-01-10", "2020-01-20", NA_character_))
  )
  for (files in list(paths, rev(paths))) {
    out <- expect_period_change_warning(extract2d(x, files, var = "temp",
      time_report = FALSE
    ))
    expect_equal(out$temp, c(11, 11, NA_real_))
    expect_equal(out[names(x)], x)
    expect_null(attr(out, ".source_candidate", exact = TRUE))
    expect_equal(
      extract2d(x, files, var = "temp", time_match = "exact",
        time_report = FALSE
      )$temp,
      rep(NA_real_, 3)
    )
    expect_equal(
      extract2d(x, files, var = "temp", max_time_diff = 0.5,
        time_report = FALSE
      )$temp,
      rep(NA_real_, 3)
    )
    expect_equal(
      expect_period_change_warning(
        extract2d(x, files, var = "temp", max_time_diff = 1,
          time_report = FALSE
        )
      )$temp,
      c(11, NA_real_, NA_real_)
    )
    expect_equal(
      expect_period_change_warning(
        extract2d(
          data = NULL, nc = files, lon = 0, lat = 40,
          date = as.Date("2020-01-10"), var = "temp",
          time_report = FALSE
        )
      ),
      11
    )
  }
})

test_that("singleton candidates respect exact matching and temporal ties", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_a2_temporal_nc(paths[1], 0, 1)
  make_a2_temporal_nc(paths[2], 10, 11)
  x <- data.frame(lon = 0, lat = 40, date = as.Date("2020-01-10"))
  for (files in list(paths, rev(paths))) {
    expect_equal(
      expect_period_change_warning(
        extract2d(x, files, var = "temp",
          time_report = FALSE
        )
      )$temp, 11
    )
    x$date <- as.Date("2020-01-11")
    expect_equal(
      extract2d(x, files, var = "temp", time_match = "exact",
        time_report = FALSE
      )$temp, 11
    )
    x$date <- as.Date("2020-01-10")
  }
  x$date <- as.Date("2020-01-06")
  expect_equal(
    expect_period_change_warning(
      extract2d(x, paths, var = "temp",
        time_report = FALSE
      )
    )$temp, 1
  )
  expect_equal(
    expect_period_change_warning(
      extract2d(x, rev(paths), var = "temp",
        time_report = FALSE
      )
    )$temp, 11
  )
})

test_that("missing nearest values fall back within matching tolerance", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    paths <- replicate(2, tempfile(fileext = ".nc"))
    on.exit(unlink(paths), add = TRUE)
    make_a2_temporal_nc(paths[1], 8, NA_real_)
    make_a2_temporal_nc(paths[2], 0, 1)
    x <- data.frame(lon = 0, lat = 40, date = as.Date("2020-01-10"))
    for (files in list(paths, rev(paths))) {
      expect_equal(
        expect_period_change_warning(
          extract2d(x, files, var = "temp",
            time_report = FALSE
          )
        )$temp, 1
      )
      expect_true(is.na(
        extract2d(x, files, var = "temp", max_time_diff = 2,
          time_report = FALSE
        )$temp
      ))
    }
  })
})

test_that("3D wrappers preserve vertical summaries across files", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_a2_temporal_nc(paths[1], 0, c(1, 2, 3), with_depth = TRUE)
  make_a2_temporal_nc(paths[2], 10, c(11, 12, 13), with_depth = TRUE)
  x <- data.frame(lon = 0, lat = 40, depth = 60,
                  date = as.Date("2020-01-10"))
  for (files in list(paths, rev(paths))) {
    expect_equal(
      expect_period_change_warning(
        extract3d_surface(x, files, var = "temp",
          time_report = FALSE
        )
      )$temp, 11
    )
    expect_equal(
      expect_period_change_warning(
        extract3d_nearest(x, files, var = "temp",
          time_report = FALSE
        )
      )$temp, 12
    )
    expect_equal(
      expect_period_change_warning(
        extract3d_bottom(x, files, var = "temp",
          time_report = FALSE
        )
      )$temp, 13
    )
    out <- expect_period_change_warning(extract3d_all(x, files, var = "temp",
      time_report = FALSE
    ))
    expect_equal(out$nearest_temp, 12)
    expect_equal(out$surface_temp, 11)
    expect_equal(out$seabottom_temp, 13)
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

test_that("temporal candidates handle partially and entirely missing axes", {
  timestamps <- as.POSIXct(c(NA_character_, "2020-01-01"), tz = "UTC")
  target <- as.POSIXct("2020-01-10", tz = "UTC")
  expect_true(.temporal_source_eligible(timestamps, target, "daily"))
  expect_false(.temporal_source_eligible(timestamps[1], target, "daily"))
})


# Issue #73 B: diagnostics for the selected temporal values.
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

test_that("single source provenance supports warnings and missing results", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 1, 11)
  x <- data.frame(lon = 0, lat = 40, date = as.Date("2020-01-01"))
  result <- capture_temporal_warning(extract2d(x, path, var = "temp",
    time_report = FALSE
  ))
  expect_equal(result$value$temp, 11)
  expect_length(result$warnings, 1L)
  expect_identical(result$warnings[[1L]]$diagnostics$source_index, 1L)
  expect_identical(result$warnings[[1L]]$diagnostics$reason, "period_changed")
  expect_null(attr(result$value, ".source_candidate", exact = TRUE))
  missing <- capture_temporal_warning(
    extract2d(x, path, var = "temp", time_match = "exact",
      time_report = FALSE
    )
  )
  expect_true(is.na(missing$value$temp))
  expect_length(missing$warnings, 0L)
  direct <- capture_temporal_warning(extract2d(
    data = NULL, nc = path, lon = 0, lat = 40, date = x$date, var = "temp",
    time_report = FALSE
  ))
  expect_equal(direct$value, 11)
  expect_length(direct$warnings, 1L)
})

test_that("singleton collections infer cadence for selected values", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(3, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  for (j in seq_along(paths)) make_a2_temporal_nc(paths[j], j - 1, j)
  x <- data.frame(lon = 0, lat = 40, date = as.POSIXct(
    "2020-01-02 12:00:00", tz = "UTC"
  ))
  for (files in list(paths, rev(paths))) {
    result <- capture_temporal_warning(extract2d(x, files, var = "temp",
      time_report = FALSE
    ))
    if (identical(files, paths)) {
      expect_length(result$warnings, 0L)
    } else {
      expect_length(result$warnings, 1L)
      expect_identical(
        result$warnings[[1L]]$diagnostics$reason, "period_changed"
      )
    }
  }
  x$date <- as.POSIXct("2020-01-05 12:00:00", tz = "UTC")
  result <- capture_temporal_warning(extract2d(x, paths, var = "temp",
    time_report = FALSE
  ))
  expect_equal(result$value$temp, 3)
  expect_identical(result$warnings[[1L]]$diagnostics$cadence_basis,
                   "source_collection")
  expect_identical(result$warnings[[1L]]$diagnostics$reason, "period_changed")
})

test_that("all summaries share one warning without duplicating examples", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 1, c(11, 12, 13), with_depth = TRUE)
  x <- data.frame(lon = 0, lat = 40, depth = 60, date = as.Date("2020-01-01"))
  result <- capture_temporal_warning(extract3d_all(x, path, var = "temp",
    time_report = FALSE
  ))
  expect_length(result$warnings, 1L)
  expect_equal(nrow(result$warnings[[1L]]$diagnostics), 3)
  message <- conditionMessage(result$warnings[[1L]])
  expect_match(message, "1 observation(s)", fixed = TRUE)
  matches <- gregexpr("Row 1; variable: temp", message, fixed = TRUE)
  expect_equal(lengths(regmatches(message, matches)), 1L)
})

test_that("metadata keeps bounds separate from sampling cadence", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 0, 1)
  nc <- ncdf4::nc_open(path, write = TRUE)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  ncdf4::ncatt_put(nc, "temp", "cell_methods", "time: mean")
  ncdf4::ncatt_put(nc, "time", "bounds", "missing_bounds")
  info <- .read_temporal_metadata(
    nc, "time", "temp", .decode_netcdf_time(nc, "time")
  )
  expect_identical(info$cell_methods, "time: mean")
  expect_identical(info$bounds_status, "unavailable")
  expect_identical(info$cadence$type, "unknown")
  expect_identical(info$calendar, "standard")
  expect_match(info$units, "days since")
})

test_that("explicit bounds describe intervals without inventing cadence", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 0, 1)
  nc <- ncdf4::nc_open(path, write = TRUE)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  vertex <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
  bounds <- ncdf4::ncvar_def(
    "time_bounds", "days since 2020-01-01",
    list(vertex, nc$dim$time), prec = "double"
  )
  nc <- ncdf4::ncvar_add(nc, bounds)
  ncdf4::ncvar_put(nc, "time_bounds", matrix(c(-0.5, 0.5), nrow = 2L))
  ncdf4::ncatt_put(nc, "time", "bounds", "time_bounds")
  info <- .read_temporal_metadata(
    nc, "time", "temp", .decode_netcdf_time(nc, "time")
  )
  expect_identical(info$bounds_status, "read")
  expect_equal(info$bounds[2, 1] - info$bounds[1, 1], 86400)
  expect_identical(info$cadence$type, "unknown")
})

test_that("static variables do not produce temporal diagnostics", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_test_netcdf(path, with_depth = FALSE, with_time = FALSE)
  result <- capture_temporal_warning(extract2d(
    data.frame(lon = 0, lat = 40), path, var = "temp", date_col = NULL,
    time_report = FALSE
  ))
  expect_equal(result$value$temp, 1)
  expect_length(result$warnings, 0L)
})


# Issue #73 C1: custom intervals and visible temporal reports.
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

test_that("custom warns above its interval without rejecting values", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 0, 11)
  x <- data.frame(
    lon = 0, lat = 40,
    date = as.POSIXct("2020-01-01 03:00:00", tz = "UTC")
  )
  at_limit <- capture_temporal_warning(extract2d(
    x, path, var = "temp", time_resolution = "custom",
    custom_time_interval = as.difftime(3, units = "hours"),
    time_report = FALSE
  ))
  expect_equal(at_limit$value$temp, 11)
  expect_length(at_limit$warnings, 0L)
  x$date <- x$date + 1
  above_limit <- capture_temporal_warning(extract2d(
    x, path, var = "temp", time_resolution = "custom",
    custom_time_interval = as.difftime(3, units = "hours"),
    time_report = FALSE
  ))
  expect_equal(above_limit$value$temp, 11)
  expect_length(above_limit$warnings, 1L)
  expect_identical(
    above_limit$warnings[[1L]]$diagnostics$reason,
    "custom_interval_exceeded"
  )
})

test_that("custom exact matching still requires the same timestamp", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 0, 11)
  result <- extract2d(
    data = NULL, nc = path, lon = 0, lat = 40,
    date = as.POSIXct("2020-01-01 01:00:00", tz = "UTC"),
    var = "temp", time_match = "exact", time_resolution = "custom",
    custom_time_interval = as.difftime(3, units = "hours"),
    time_report = FALSE
  )
  expect_true(is.na(result))
})

test_that("time reports are visible without changing results", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 0, 11)
  x <- data.frame(
    lon = 0, lat = 40,
    date = as.POSIXct("2020-01-01 01:00:00", tz = "UTC")
  )
  expect_message(
    visible <- extract2d(
      x, path, var = "temp", time_resolution = "custom",
      custom_time_interval = as.difftime(3, units = "hours")
    ),
    class = "ocean3d_time_report"
  )
  expect_message(
    hidden <- extract2d(
      x, path, var = "temp", time_resolution = "custom",
      custom_time_interval = as.difftime(3, units = "hours"),
      time_report = FALSE
    ),
    NA
  )
  expect_identical(visible, hidden)
  x$date <- as.POSIXct("2020-01-01", tz = "UTC")
  expect_message(extract2d(x, path, var = "temp"), NA)
})

test_that("hiding the report preserves temporal warnings", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 0, 11)
  x <- data.frame(lon = 0, lat = 40, date = as.Date("2020-01-02"))
  expect_message(
    expect_warning(
      extract2d(x, path, var = "temp", time_report = FALSE),
      class = "ocean3d_time_warning"
    ),
    NA
  )
})

# C2: explicit output columns and preflight variable validation.
test_that("existing output columns are replaced across all sources", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_a2_temporal_nc(paths[1], 0, 11)
  make_a2_temporal_nc(paths[2], 1, 22)
  x <- data.frame(
    lon = 0, lat = 40,
    date = as.Date(c("2020-01-01", "2020-01-02")),
    temp = c(999, 999)
  )
  for (files in list(paths, rev(paths))) {
    out <- extract2d(
      x, files, var = "temp", time_match = "exact", time_report = FALSE
    )
    expect_equal(out$temp, c(11, 22))
    expect_identical(out$date, x$date)
  }
  x$date <- as.Date(c("2020-01-03", "2020-01-04"))
  out <- extract2d(
    x, paths, var = "temp", time_match = "exact", time_report = FALSE
  )
  expect_true(all(is.na(out$temp)))
})

test_that("different variables cannot share an output name", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_test_netcdf(
    paths[1], var = "temp", with_depth = FALSE, with_time = FALSE
  )
  make_test_netcdf(
    paths[2], var = "salt", with_depth = FALSE, with_time = FALSE
  )
  x <- data.frame(lon = 0, lat = 40)
  for (argument in c("output_col", "output_prefix")) {
    args <- list(data = x, nc = paths, date_col = NULL, time_report = FALSE)
    args[[argument]] <- "value"
    expect_error(do.call(extract2d, args), "Different variables")
  }
  out <- extract2d(x, paths, date_col = NULL, time_report = FALSE)
  expect_true(all(c("temp", "salt") %in% names(out)))
})

test_that("outputs cannot overwrite extraction input columns", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_test_netcdf(path, with_depth = FALSE, with_time = FALSE)
  expect_error(
    extract2d(
      data.frame(lon = 0, lat = 40), path, var = "temp",
      date_col = NULL, output_col = "lon",
      time_report = FALSE
    ),
    "overwrite extraction inputs"
  )
})

# C3: spatial coverage, cell edges and longitude conventions.
test_that("spatial indices respect edges and longitude conventions", {
  axis <- list(
    centres = c(0, 1), lower = c(-0.5, 0.5), upper = c(0.5, 1.5),
    longitude = FALSE, basis = "inferred_regular"
  )
  expect_identical(.containing_axis_indices(axis, -0.5), 1L)
  expect_identical(.containing_axis_indices(axis, 1.5), 2L)
  expect_length(.containing_axis_indices(axis, 1.6), 0L)
  axis <- list(
    centres = c(179, 181), lower = c(178, 180), upper = c(180, 182),
    longitude = TRUE, basis = "inferred_regular"
  )
  expect_identical(.containing_axis_indices(axis, -179), 2L)
  expect_identical(.containing_axis_indices(axis, 181), 2L)
  expect_length(.containing_axis_indices(axis, 0), 0L)
})

test_that("singleton axes do not invent spatial coverage", {
  axis <- list(
    centres = 40, lower = 40, upper = 40,
    longitude = FALSE, basis = "singleton_coordinate"
  )
  expect_identical(.containing_axis_indices(axis, 40), 1L)
  expect_length(.containing_axis_indices(axis, 41), 0L)
})

test_that("spatial eligibility selects the correct tile in either order", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    paths <- replicate(2, tempfile(fileext = ".nc"))
    on.exit(unlink(paths), add = TRUE)
    make_test_netcdf(
      paths[1], with_depth = FALSE, with_time = FALSE,
      lon_vals = c(-10, -5)
    )
    make_test_netcdf(
      paths[2], with_depth = FALSE, with_time = FALSE,
      lon_vals = c(0, 5)
    )
    nc <- ncdf4::nc_open(paths[1], write = TRUE)
    tryCatch(
      ncdf4::ncvar_put(nc, "temp", matrix(18, 2, 2)),
      finally = ncdf4::nc_close(nc)
    )
    nc <- ncdf4::nc_open(paths[2], write = TRUE)
    tryCatch(
      ncdf4::ncvar_put(nc, "temp", matrix(24, 2, 2)),
      finally = ncdf4::nc_close(nc)
    )
    x <- data.frame(lon = c(2, -7, 30), lat = c(40, 40, 40))
    for (files in list(paths, rev(paths))) {
      out <- extract2d(
        x, files, var = "temp", date_col = NULL, time_report = FALSE
      )
      expect_equal(out$temp, c(24, 18, NA_real_))
    }
  })
})

make_spatial_axis_nc <- function(path, centres, bounds = NULL) {
  lon <- ncdf4::ncdim_def("lon", "degrees_east", centres)
  lat <- ncdf4::ncdim_def("lat", "degrees_north", 40)
  variables <- list(ncdf4::ncvar_def(
    "temp", "degree_Celsius", list(lon, lat), -9999
  ))
  if (!is.null(bounds)) {
    vertex <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
    variables[[2L]] <- ncdf4::ncvar_def(
      "lon_bounds", "degrees_east", list(vertex, lon), prec = "double"
    )
  }
  nc <- ncdf4::nc_create(path, variables)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  ncdf4::ncvar_put(nc, "temp", matrix(1, length(centres), 1))
  if (!is.null(bounds)) {
    ncdf4::ncvar_put(nc, "lon_bounds", bounds)
    ncdf4::ncatt_put(nc, "lon", "bounds", "lon_bounds")
  }
}

test_that("explicit spatial bounds override inferred edges", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_spatial_axis_nc(
    path, centres = c(0, 2),
    bounds = matrix(c(-0.2, 0.2, 1.8, 2.2), nrow = 2L)
  )
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  axis <- .prepare_spatial_axis(nc, "lon", c(0, 2), longitude = TRUE)
  expect_identical(axis$basis, "explicit_bounds")
  expect_identical(.containing_axis_indices(axis, 0.1), 1L)
  expect_length(.containing_axis_indices(axis, 1), 0L)
})

test_that("irregular descending axes preserve original indices", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_spatial_axis_nc(path, centres = c(5, 2, 0))
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  axis <- .prepare_spatial_axis(nc, "lon", c(5, 2, 0), longitude = TRUE)
  expect_identical(axis$basis, "inferred_irregular")
  expect_identical(.containing_axis_indices(axis, 5), 1L)
  expect_identical(.containing_axis_indices(axis, 0), 3L)
  expect_length(.containing_axis_indices(axis, 20), 0L)
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

test_that("explanatory warnings keep full paths internally", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, c(0, 1, 2), c(11, 12, 13))
  x <- data.frame(
    lon = 0, lat = 40,
    date = as.POSIXct("2020-01-02 12:00:00", tz = "UTC")
  )
  expect_message(
    capture_temporal_warning(extract2d(
      x, path, var = "temp", time_report = TRUE
    )),
    "No temporal warning is required", fixed = TRUE
  )
  x$date <- as.POSIXct("2020-01-05 12:00:00", tz = "UTC")
  captured <- capture_temporal_warning(extract2d(
    x, path, var = "temp", time_report = FALSE
  ))
  text <- conditionMessage(captured$warnings[[1L]])
  expect_match(text, basename(path), fixed = TRUE)
  expect_false(grepl(dirname(path), text, fixed = TRUE))
  expect_identical(captured$warnings[[1L]]$diagnostics$source_file, path)
})

test_that("missing output reports explain what is unavailable", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 0, 11)
  expect_message(
    out <- extract2d(
      data = NULL, nc = path, lon = 0, lat = 40,
      date = as.Date("2020-01-02"), var = "temp",
      time_match = "exact", time_report = TRUE
    ),
    "No timestamp satisfied the requested temporal settings", fixed = TRUE
  )
  expect_true(is.na(out))
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


# Issue #74 A: geographic selection and explicit horizontal limits.
make_horizontal_test_nc <- function(path, lon, lat, values, days = NULL) {
  dims <- list(
    ncdf4::ncdim_def("lon", "degrees_east", lon),
    ncdf4::ncdim_def("lat", "degrees_north", lat)
  )
  if (!is.null(days)) {
    dims <- c(dims, list(ncdf4::ncdim_def(
      "time", "days since 2020-01-01", days
    )))
  }
  variable <- ncdf4::ncvar_def(
    "temp", "degree_Celsius", dims, -9999, prec = "double"
  )
  nc <- ncdf4::nc_create(path, variable)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  shape <- vapply(dims, function(x) length(x$vals), integer(1))
  ncdf4::ncvar_put(nc, "temp", array(values, dim = shape))
}

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

test_that("within rejects distant coordinates with one or several files", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    paths <- replicate(2, tempfile(fileext = ".nc"))
    on.exit(unlink(paths), add = TRUE)
    for (path in paths) make_horizontal_test_nc(path, c(0, 1), c(40, 41), 11)
    x <- data.frame(lon = c(1.5, 1.6, 20), lat = 40)
    for (files in list(paths[1], paths, rev(paths))) {
      out <- extract2d(x, files, var = "temp", date_col = NULL)
      expect_equal(out$temp, c(11, NA_real_, NA_real_))
    }
  })
})

test_that("nearest accepts nearby points outside the grid up to its limit", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path), add = TRUE)
    make_horizontal_test_nc(path, 0, 0, 11)
    x <- data.frame(lon = 0.01, lat = 0)
    out <- extract2d(
      x, path, var = "temp", date_col = NULL,
      spatial_match = "nearest", max_distance_km = 2
    )
    expect_equal(out$temp, 11)
    out <- extract2d(
      x, path, var = "temp", date_col = NULL,
      spatial_match = "nearest", max_distance_km = 0.5
    )
    expect_true(is.na(out$temp))
    expect_equal(extract2d(
      data = NULL, nc = path, lon = 0, lat = 0, var = "temp",
      spatial_match = "nearest", max_distance_km = 0
    ), 11)
  })
})

test_that("a cell exactly at the kilometre limit is accepted", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path), add = TRUE)
    make_horizontal_test_nc(path, 0, 0, 11)
    boundary <- .great_circle_km(0.01, 0, 0, 0)
    x <- data.frame(lon = 0.01, lat = 0)
    expect_equal(extract2d(
      x, path, var = "temp", date_col = NULL,
      spatial_match = "nearest", max_distance_km = boundary
    )$temp, 11)
    expect_true(is.na(extract2d(
      x, path, var = "temp", date_col = NULL,
      spatial_match = "nearest", max_distance_km = boundary - 1e-6
    )$temp))
  })
})

test_that("geographic proximity is not independent degree subtraction", {
  # At high latitude the geographically closest centre is not necessarily
  # the one with the smallest latitude difference.
  result <- .nearest_horizontal_cell(60, c(80, 85), 1L, 1:2, 0, 80)
  expect_identical(result$lat_idx, 2L)
  expect_identical(result$lon_idx, 1L)
  expect_gt(result$distance_km, 0)
  whole <- .nearest_horizontal_cells(60, c(80, 85), 0, 80)[[1L]]
  expect_equal(result, whole)
})

test_that("longitudes use both conventions and wrap around the antimeridian", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_horizontal_test_nc(paths[1], c(-160, 170), 40, c(11, 22))
  expect_equal(extract2d(
    data.frame(lon = 200, lat = 40), paths[1],
    var = "temp", date_col = NULL
  )$temp, 11)
  make_horizontal_test_nc(paths[2], c(0, 359), 40, c(11, 22))
  for (longitude in c(-0.1, 359.9)) {
    expect_equal(extract2d(
      data.frame(lon = longitude, lat = 40), paths[2],
      var = "temp", date_col = NULL
    )$temp, 11)
  }
  selected <- .nearest_horizontal_cell(
    c(179, -170), 40, 1:2, 1L, -179.9, 40
  )
  expect_identical(selected$lon_idx, 1L)
  expect_equal(.normalise_longitude(c(-160, 200)), c(-160, -160))
})

test_that("equal temporal candidates prefer the spatially closer source", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    paths <- replicate(2, tempfile(fileext = ".nc"))
    on.exit(unlink(paths), add = TRUE)
    make_horizontal_test_nc(paths[1], 0, 0, 11, days = 0)
    make_horizontal_test_nc(paths[2], 0.02, 0, 22, days = 0)
    x <- data.frame(lon = 0.015, lat = 0, date = as.Date("2020-01-01"))
    for (files in list(paths, rev(paths))) {
      out <- extract2d(
        x, files, var = "temp", spatial_match = "nearest",
        max_distance_km = 5, time_report = FALSE
      )
      expect_equal(out$temp, 22)
    }
  })
})

test_that("temporal proximity remains the first selection criterion", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    paths <- replicate(2, tempfile(fileext = ".nc"))
    on.exit(unlink(paths), add = TRUE)
    make_horizontal_test_nc(paths[1], 0, 0, 11, days = 0)
    make_horizontal_test_nc(paths[2], 0.02, 0, 22, days = 1)
    x <- data.frame(lon = 0, lat = 0, date = as.Date("2020-01-02"))
    for (files in list(paths, rev(paths))) {
      expect_equal(extract2d(
        x, files, var = "temp", spatial_match = "nearest",
        max_distance_km = 5, time_report = FALSE
      )$temp, 22)
    }
  })
})

test_that("input order breaks ties in both time and geographic distance", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_horizontal_test_nc(paths[1], 0, 0, 11, days = 0)
  make_horizontal_test_nc(paths[2], 0, 0, 22, days = 0)
  x <- data.frame(lon = 0, lat = 0, date = as.Date("2020-01-01"))
  expect_equal(extract2d(x, paths, var = "temp")$temp, 11)
  expect_equal(extract2d(x, rev(paths), var = "temp")$temp, 22)
})

test_that("nearest does not search past a missing cell within the same file", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    paths <- replicate(2, tempfile(fileext = ".nc"))
    on.exit(unlink(paths), add = TRUE)
    make_horizontal_test_nc(paths[1], c(0, 0.02), 0, c(NA_real_, 22))
    x <- data.frame(lon = 0, lat = 0)
    expect_true(is.na(extract2d(
      x, paths[1], var = "temp", date_col = NULL,
      spatial_match = "nearest", max_distance_km = 5
    )$temp))
    # Another eligible file may still supply its selected valid candidate.
    make_horizontal_test_nc(paths[2], 0.01, 0, 33)
    for (files in list(paths, rev(paths))) {
      expect_equal(extract2d(
        x, files, var = "temp", date_col = NULL,
        spatial_match = "nearest", max_distance_km = 5
      )$temp, 33)
    }
  })
})

match_one_cell <- function(lon_axis, lat_axis, lon_values, lat_values,
                           lon_target, lat_target, spatial_match,
                           max_distance_km) {
  match_cells <- getFromNamespace(".match_horizontal_cells", "ocean3d")
  match_row <- getFromNamespace(".horizontal_match_row", "ocean3d")
  match_row(match_cells(
    lon_axis, lat_axis, lon_values, lat_values, lon_target, lat_target,
    spatial_match, max_distance_km
  ), 1L)
}

test_that("horizontal selection records accepted and rejected candidates", {
  inside <- list(
    centres = 0, lower = 0, upper = 0, longitude = FALSE
  )
  longitude <- inside
  longitude$longitude <- TRUE
  accepted <- match_one_cell(
    longitude, inside, 0, 0, 0.01, 0, "nearest", 2
  )
  expect_true(accepted$accepted)
  expect_identical(accepted$reason, "matched")
  expect_equal(accepted$selected_lon, 0)
  expect_equal(accepted$selected_lat, 0)
  expect_gt(accepted$distance_km, 1)
  rejected <- match_one_cell(
    longitude, inside, 0, 0, 0.01, 0, "nearest", 0.5
  )
  expect_false(rejected$accepted)
  expect_identical(rejected$reason, "distance_limit_exceeded")
  expect_gt(rejected$distance_km, 1)
  outside <- match_one_cell(
    longitude, inside, 0, 0, 0.01, 0, "within", NULL
  )
  expect_identical(outside$reason, "outside_grid")
  invalid <- match_one_cell(
    longitude, inside, 0, 0, NA_real_, 0, "nearest", 5
  )
  expect_identical(invalid$reason, "invalid_coordinates")
})

test_that("nearest horizontal selection also works through 3D wrappers", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path), add = TRUE)
    make_test_netcdf(path, with_time = FALSE)
    x <- data.frame(lon = -0.6, lat = 40, depth = 60)
    out <- extract3d_all(
      x, path, var = "temp", date_col = NULL,
      spatial_match = "nearest", max_distance_km = 100
    )
    expect_equal(out$surface_temp, 1)
    expect_equal(out$nearest_temp, 5)
    expect_equal(out$seabottom_temp, 9)
    expect_true(is.na(extract3d_surface(
      x, path, var = "temp", date_col = NULL
    )$temp))
  })
})


# Issue #74 B: physical depths independently of coordinate storage order.
make_depth_axis_test_nc <- function(path, depths, values,
                                    units = "m", positive = NULL,
                                    depth_name = "depth",
                                    standard_name = NULL,
                                    formula_terms = NULL) {
  lon <- ncdf4::ncdim_def("lon", "degrees_east", 0)
  lat <- ncdf4::ncdim_def("lat", "degrees_north", 40)
  depth <- ncdf4::ncdim_def(depth_name, units, depths)
  variable <- ncdf4::ncvar_def(
    "temp", "degree_Celsius", list(lon, lat, depth),
    -9999, prec = "double"
  )
  nc <- ncdf4::nc_create(path, variable)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  ncdf4::ncvar_put(nc, "temp", array(values, c(1, 1, length(depths))))
  attributes <- list(
    positive = positive, standard_name = standard_name,
    formula_terms = formula_terms
  )
  for (name in names(attributes)) {
    if (!is.null(attributes[[name]])) {
      ncdf4::ncatt_put(nc, depth_name, name, attributes[[name]])
    }
  }
}

test_that("physical depth conversion retains original coordinate positions", {
  descending <- .prepare_depth_axis(c(100, 50, 0), "m", "down")
  expect_equal(descending$depth_m, c(100, 50, 0))
  expect_identical(descending$surface_index, 3L)
  expect_identical(descending$direction_basis, "metadata")
  upward <- .prepare_depth_axis(c(0, -50, -100), "m", " UP ")
  expect_equal(upward$depth_m, c(0, 50, 100))
  expect_identical(upward$surface_index, 1L)
  expect_identical(upward$direction_basis, "metadata")
  for (unit in c("m", "meter", "meters", "metre", "metres")) {
    expect_equal(.prepare_depth_axis(c(0, 50), unit)$depth_m, c(0, 50))
  }
  expect_equal(.prepare_depth_axis(c(0, 0.05), "km")$depth_m, c(0, 50))
  expect_equal(.prepare_depth_axis(c(0, 5000), "cm")$depth_m, c(0, 50))
  inferred <- .prepare_depth_axis(c(-100, -50, 0), "m")
  expect_equal(inferred$depth_m, c(100, 50, 0))
  expect_identical(inferred$direction_basis, "inferred")
})

test_that("ambiguous or unsupported depth coordinates are rejected", {
  expect_error(.prepare_depth_axis(c(-50, 0, 50), "m"), "ambiguous")
  expect_error(.prepare_depth_axis(c(0, 50), "dbar"), "converted to metres")
  expect_error(.prepare_depth_axis(c(0, 50), NULL), "explicit supported")
  expect_error(.prepare_depth_axis(c(0, 50), NA_character_), "units")
  expect_error(.prepare_depth_axis(c(0, 50), "m", "sideways"), "up.*down")
  expect_error(.prepare_depth_axis(c(0, 50), "m", NA_character_), "up.*down")
  expect_error(.prepare_depth_axis(c(0, NA_real_), "m"), "finite")
  expect_error(.prepare_depth_axis(numeric(), "m"), "finite")
  expect_error(.prepare_depth_axis(c(0, Inf), "m"), "finite")
  expect_error(.prepare_depth_axis(c(0, 50, 50), "m"), "duplicate")
  expect_error(.prepare_depth_axis(c(0, -50), "m", "down"), "reference surface")
  expect_error(.prepare_depth_axis(c(0, 50), "m", "up"), "reference surface")
})

test_that("depth methods select physical layers for every storage order", {
  skip_if_not_installed("ncdf4")
  cases <- list(
    list(depth = c(0, 50, 100), values = c(11, 22, 33), positive = "down"),
    list(depth = c(100, 50, 0), values = c(33, 22, 11), positive = "down"),
    list(depth = c(50, 0, 100), values = c(22, 11, 33), positive = "down"),
    list(depth = c(0, -50, -100), values = c(11, 22, 33), positive = "up"),
    list(depth = c(-100, -50, 0), values = c(33, 22, 11), positive = "up")
  )
  paths <- replicate(length(cases), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  x <- data.frame(lon = 0, lat = 40, depth = 60)
  for (i in seq_along(cases)) {
    item <- cases[[i]]
    make_depth_axis_test_nc(
      paths[i], item$depth, item$values,
      positive = item$positive, depth_name = "z"
    )
    out <- extract3d_all(x, paths[i], var = "temp", date_col = NULL)
    expect_equal(out$surface_temp, 11)
    expect_equal(out$nearest_temp, 22)
    expect_equal(out$seabottom_temp, 33)
    expect_equal(extract3d_surface(
      x, paths[i], var = "temp", date_col = NULL
    )$temp, out$surface_temp)
    expect_equal(extract3d_nearest(
      x, paths[i], var = "temp", date_col = NULL
    )$temp, out$nearest_temp)
    expect_equal(extract3d_bottom(
      x, paths[i], var = "temp", date_col = NULL
    )$temp, out$seabottom_temp)
  }
})

test_that("file depth units convert to observation metres", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_depth_axis_test_nc(paths[1], c(0.1, 0.05, 0), c(33, 22, 11), "km")
  make_depth_axis_test_nc(paths[2], c(0, -5000, -10000), c(11, 22, 33), "cm")
  x <- data.frame(lon = 0, lat = 40, depth = 60)
  for (path in paths) {
    out <- extract3d_all(x, path, var = "temp", date_col = NULL)
    expect_equal(out$surface_temp, 11)
    expect_equal(out$nearest_temp, 22)
    expect_equal(out$seabottom_temp, 33)
  }
})

test_that("nearest depth ties choose the shallower layer in both orders", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_depth_axis_test_nc(paths[1], c(0, 50, 100), c(11, 22, 33))
  make_depth_axis_test_nc(paths[2], c(100, 50, 0), c(33, 22, 11))
  x <- data.frame(lon = 0, lat = 40, depth = 75)
  for (path in paths) {
    expect_equal(extract3d_nearest(
      x, path, var = "temp", date_col = NULL
    )$temp, 22)
  }
})

test_that("surface preserves a missing shallowest layer in descending files", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path), add = TRUE)
    make_depth_axis_test_nc(path, c(100, 50, 0), c(33, 22, NA_real_))
    x <- data.frame(lon = 0, lat = 40, depth = 60)
    out <- extract3d_all(x, path, var = "temp", date_col = NULL)
    expect_true(is.na(out$surface_temp))
    expect_equal(out$nearest_temp, 22)
    expect_equal(out$seabottom_temp, 33)
    expect_true(is.na(extract3d_surface(
      x, path, var = "temp", date_col = NULL
    )$temp))
  })
})

test_that("bottom and nearest ignore missing layers by physical depth", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(100, 50, 0), c(NA_real_, 22, 11))
  x <- data.frame(lon = 0, lat = 40, depth = 90)
  out <- extract3d_all(x, path, var = "temp", date_col = NULL)
  expect_equal(out$surface_temp, 11)
  expect_equal(out$nearest_temp, 22)
  expect_equal(out$seabottom_temp, 22)
})

test_that("missing requested depths preserve the other physical summaries", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path), add = TRUE)
    make_depth_axis_test_nc(path, c(100, 50, 0), c(33, 22, 11))
    out <- extract3d_all(
      data.frame(lon = 0, lat = 40, depth = NA_real_),
      path, var = "temp", date_col = NULL
    )
    expect_true(is.na(out$nearest_temp))
    expect_equal(out$surface_temp, 11)
    expect_equal(out$seabottom_temp, 33)
    summary <- .summarise_depth_profile(
      rep(NA_real_, 3), c(100, 50, 0), 60
    )
    expect_true(all(is.na(summary$value)))
  })
})

test_that("depth metadata are read from the dimension coordinate variable", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(0, -50, -100), c(11, 22, 33), positive = "up")
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  axis <- .read_depth_axis(nc, "depth", nc$dim$depth$vals)
  expect_identical(axis$direction_basis, "metadata")
  expect_equal(axis$depth_m, c(0, 50, 100))
})

test_that("formula-based and pressure axes require explicit conversion", {
  skip_if_not_installed("ncdf4")
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
  skip_if_not_installed("ncdf4")
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

test_that("physical depth summaries also support direct point input", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(-100, -50, 0), c(33, 22, 11), positive = "up")
  out <- extract3d_all(
    data = NULL, nc = path, lon = 0, lat = 40, depth = 60, var = "temp"
  )
  expect_equal(
    out, c(nearest_temp = 22, surface_temp = 11, seabottom_temp = 33)
  )
})

# Issue #74 C: explanatory reports preserve values and source attribution.
capture_selection_reports <- function(expr) {
  reports <- list()
  value <- withCallingHandlers(expr, message = function(condition) {
    reports[[length(reports) + 1L]] <<- condition
    invokeRestart("muffleMessage")
  })
  list(value = value, reports = reports,
       text = paste(vapply(reports, conditionMessage, character(1)),
                    collapse = "\n"))
}

test_that("spatial reports explain centres and explicit distance limits", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_horizontal_test_nc(path, c(0, 1), c(40, 41), 1:4)
  x <- data.frame(lon = 0.1, lat = 40)
  within <- capture_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL, spatial_report = TRUE
  ))
  expect_match(within$text, "matched the requested selection", fixed = TRUE)
  expect_match(within$text, "containing cell centre: (0, 40)", fixed = TRUE)
  nearby <- capture_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL, spatial_report = TRUE,
    spatial_match = "nearest", max_distance_km = 20
  ))
  expect_match(nearby$text, "allowed: 20 km", fixed = TRUE)
  rejected <- capture_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL, spatial_report = TRUE,
    spatial_match = "nearest", max_distance_km = 1
  ))
  expect_true(is.na(rejected$value$temp))
  expect_match(rejected$text, "allowed: 1 km", fixed = TRUE)
  expect_match(rejected$text, "Result: NA", fixed = TRUE)
  silent <- capture_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL,
    spatial_report = FALSE, depth_report = FALSE
  ))
  expect_length(silent$reports, 0L)
  expect_equal(silent$value, within$value)
})

test_that("depth reports describe physical layers without interpolation", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(100, 50, 0), c(33, 22, 11))
  x <- data.frame(lon = 0, lat = 40, depth = 60)
  result <- capture_selection_reports(extract3d_all(
    x, path, var = "temp", date_col = NULL, depth_report = TRUE
  ))
  expect_match(result$text, "Requested: 60 m", fixed = TRUE)
  expect_match(result$text, "10 m shallower", fixed = TRUE)
  expect_match(result$text, "selected: 50 m", fixed = TRUE)
  expect_match(result$text, "Shallowest layer selected", fixed = TRUE)
  expect_match(result$text, "does not establish the seabed", fixed = TRUE)
  expect_match(result$text, "No vertical interpolation", fixed = TRUE)
  quiet <- capture_selection_reports(extract3d_all(
    x, path, var = "temp", date_col = NULL, depth_report = FALSE
  ))
  expect_length(quiet$reports, 0L)
  expect_equal(quiet$value, result$value)
})

test_that("missing layer reports retain rejected sources without a winner", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  for (path in paths) {
    make_depth_axis_test_nc(path, c(0, 50, 100), c(11, 22, 33))
  }
  result <- capture_selection_reports(extract3d_nearest(
    data.frame(lon = 0, lat = 40, depth = NA_real_), paths,
    var = "temp", date_col = NULL
  ))
  expect_true(is.na(result$value$temp))
  expect_match(result$text, "Requested depth is missing", fixed = TRUE)
  for (path in paths) expect_match(result$text, basename(path), fixed = TRUE)
  expect_false(grepl("Selected depth:", result$text, fixed = TRUE))
})

test_that("reports follow the winning source separately for each summary", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_depth_axis_test_nc(paths[1], c(0, 50, 100), c(NA, 22, 33))
  make_depth_axis_test_nc(paths[2], c(5, 60, 120), c(44, 55, 66))
  result <- capture_selection_reports(extract3d_all(
    data.frame(lon = 0, lat = 40, depth = 60), paths,
    var = "temp", date_col = NULL, depth_report = TRUE
  ))
  expect_equal(result$value$surface_temp, 44)
  expect_equal(result$value$nearest_temp, 22)
  expect_match(result$text, "Shallowest layer selected: 5 m", fixed = TRUE)
  expect_match(result$text, "selected: 50 m", fixed = TRUE)
  for (path in paths) expect_match(result$text, basename(path), fixed = TRUE)
})

test_that("selection report switches require scalar logical values", {
  for (bad in list(NA, 1, "yes", logical(), c(TRUE, FALSE))) {
    expect_error(extract_to_point(spatial_report = bad), "spatial_report")
    expect_error(extract_to_point(depth_report = bad), "depth_report")
  }
})

test_that("missing reports follow their cause and are never repeated", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_a2_temporal_nc(path, 0, 11)
  x <- data.frame(lon = 0, lat = 40, date = as.Date("2020-01-02"))
  temporal <- capture_selection_reports(extract2d(
    x, path, var = "temp", time_match = "exact"
  ))
  expect_length(temporal$reports, 1L)
  expect_s3_class(temporal$reports[[1L]], "ocean3d_time_report")
  expect_match(temporal$text, "No timestamp satisfied", fixed = TRUE)
  x$lon <- 20
  spatial <- capture_selection_reports(extract2d(x, path, var = "temp"))
  expect_length(spatial$reports, 1L)
  expect_s3_class(spatial$reports[[1L]], "ocean3d_spatial_report")
  expect_match(spatial$text, "outside cell coverage", fixed = TRUE)
  hidden <- capture_selection_reports(extract2d(
    x, path, var = "temp", spatial_report = FALSE
  ))
  expect_length(hidden$reports, 0L)
  expect_equal(hidden$value, spatial$value)
})

test_that("missing requested depths obey only the depth report switch", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(0, 50, 100), c(11, 22, 33))
  x <- data.frame(lon = 0, lat = 40, depth = NA_real_)
  shown <- capture_selection_reports(extract3d_nearest(
    x, path, var = "temp", date_col = NULL, spatial_report = FALSE
  ))
  expect_length(shown$reports, 1L)
  expect_s3_class(shown$reports[[1L]], "ocean3d_depth_report")
  hidden <- capture_selection_reports(extract3d_nearest(
    x, path, var = "temp", date_col = NULL, depth_report = FALSE
  ))
  expect_length(hidden$reports, 0L)
  expect_equal(shown$value, hidden$value)
})

test_that("automatic reports stay silent for normal selections", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_horizontal_test_nc(path, c(0, 1), c(40, 41), 1:4)
  x <- data.frame(lon = c(0, 0.1), lat = 40)
  automatic <- capture_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL
  ))
  expect_length(automatic$reports, 0L)
  nearest <- capture_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL,
    spatial_match = "nearest", max_distance_km = 20
  ))
  expect_length(nearest$reports, 0L)
  full <- capture_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL, spatial_report = TRUE
  ))
  expect_length(full$reports, 1L)
  expect_match(full$text, "2 observations", fixed = TRUE)
  expect_match(full$text, "netCDF variable: temp", fixed = TRUE)
  expect_match(full$text, "output column: temp", fixed = TRUE)
  expect_match(full$text, "Distance: 0 km", fixed = TRUE)
  expect_equal(automatic$value, full$value)
  expect_equal(automatic$value, nearest$value)
})

test_that("automatic depth reports distinguish normal and different layers", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(100, 50, 0), c(33, 22, 11))
  x <- data.frame(lon = 0, lat = 40, depth = 50)
  automatic <- capture_selection_reports(extract3d_all(
    x, path, var = "temp", date_col = NULL
  ))
  expect_length(automatic$reports, 0L)
  full <- capture_selection_reports(extract3d_all(
    x, path, var = "temp", date_col = NULL, depth_report = TRUE
  ))
  expect_length(full$reports, 1L)
  expect_match(full$text, "exact match", fixed = TRUE)
  expect_match(full$text, "3 output values reported", fixed = TRUE)
  expect_equal(automatic$value, full$value)
  x$depth <- 60
  different <- capture_selection_reports(extract3d_all(
    x, path, var = "temp", date_col = NULL
  ))
  expect_length(different$reports, 0L)
  hidden <- capture_selection_reports(extract3d_all(
    x, path, var = "temp", date_col = NULL, depth_report = FALSE
  ))
  expect_length(hidden$reports, 0L)
  expect_equal(hidden$value, different$value)
})

test_that("automatic spatial reports explain accepted external selections", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_horizontal_test_nc(path, c(0, 1), c(40, 41), 1:4)
  result <- capture_selection_reports(extract2d(
    data.frame(lon = 1.6, lat = 40), path, var = "temp",
    output_col = "temperature", date_col = NULL,
    spatial_match = "nearest", max_distance_km = 100
  ))
  expect_length(result$reports, 1L)
  expect_match(result$text, "outside cell coverage used the nearest cell",
               fixed = TRUE)
  expect_match(result$text, "netCDF variable: temp", fixed = TRUE)
  expect_match(result$text, "output column: temperature", fixed = TRUE)
  expect_match(result$text, "nearest cell centre: (1, 40)", fixed = TRUE)
})

test_that("cell containment is separate from overall coverage", {
  lon <- list(centres = c(0, 1), lower = c(-0.01, -0.5),
              upper = c(0.01, 1.5), longitude = TRUE)
  lat <- list(centres = 0, lower = 0, upper = 0, longitude = FALSE)
  result <- match_one_cell(lon, lat, c(0, 1), 0, 0.1, 0,
                           "nearest", 100)
  expect_true(result$accepted)
  expect_true(result$inside_coverage)
  expect_false(result$inside_selected_cell)
  expect_equal(result$lon_idx, 1L)
  lon <- list(centres = c(179, 181), lower = c(178, 180),
              upper = c(180, 182), longitude = TRUE)
  wrapped <- match_one_cell(lon, lat, c(179, -179), 0, 181, 0,
                            "nearest", 100)
  expect_true(wrapped$inside_selected_cell)
  expect_true(wrapped$inside_coverage)
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

test_that("large reports count observations and prioritise missing results", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(0, 50, 100), c(11, 22, 33))
  x <- data.frame(lon = 0, lat = 40, depth = c(NA, rep(60, 6)))
  result <- capture_selection_reports(extract3d_all(
    x, path, var = "temp", date_col = NULL, depth_report = TRUE
  ))
  expect_length(result$reports, 1L)
  expect_match(result$text, "7 observations", fixed = TRUE)
  expect_match(
    result$text, "21 output values reported; 1 missing", fixed = TRUE
  )
  expect_match(result$text, "Showing 5 of 21 cases", fixed = TRUE)
  expect_match(result$text, "16 further examples omitted", fixed = TRUE)
  expect_identical(result$reports[[1L]]$records[[1L]]$status, "missing")
  expect_length(result$reports[[1L]]$records, 21L)
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

test_that("normal nearest layers are silent unless a report is requested", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(50, 61.52, 75), c(11, 22, 33))
  x <- data.frame(lon = 0, lat = 40, depth = 60)
  automatic <- capture_selection_reports(extract3d_nearest(
    x, path, var = "temp", date_col = NULL
  ))
  expect_length(automatic$reports, 0L)
  expect_equal(automatic$value$temp, 22)
  full <- capture_selection_reports(extract3d_nearest(
    x, path, var = "temp", date_col = NULL, depth_report = TRUE
  ))
  expect_equal(automatic$value, full$value)
  expect_match(full$text, "1.52 m deeper", fixed = TRUE)
  expect_match(full$text, "11.52 m apart", fixed = TRUE)
  expect_match(full$text, "Surrounding layers:", fixed = TRUE)
  expect_false(grepl("Showing", full$text, fixed = TRUE))
})

test_that("distant layers explain spacing and proven missing neighbours", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_depth_axis_test_nc(paths[1], c(55, 67, 150), c(NA, NA, 33))
  make_depth_axis_test_nc(paths[2], c(-15000, -6700, -5500), c(33, NA, NA),
                          units = "cm", positive = "up")
  x <- data.frame(lon = 0, lat = 40, depth = c(60, 55))
  for (path in paths) {
    automatic <- capture_selection_reports(extract3d_nearest(
      x, path, var = "temp", date_col = NULL
    ))
    expect_equal(automatic$value$temp, c(33, 33))
    expect_match(
      automatic$text, "more than one local layer spacing", fixed = TRUE
    )
    expect_match(automatic$text, "12 m apart", fixed = TRUE)
    expect_match(
      automatic$text, "closer layers have no valid value", fixed = TRUE
    )
    expect_match(automatic$text, "Requested layer and nearest neighbour:",
                 fixed = TRUE)
    expect_false(grepl("Showing", automatic$text, fixed = TRUE))
    hidden <- capture_selection_reports(extract3d_nearest(
      x, path, var = "temp", date_col = NULL, depth_report = FALSE
    ))
    expect_length(hidden$reports, 0L)
    expect_equal(hidden$value, automatic$value)
  }
})

test_that("out of range and singleton reports do not invent missing layers", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_depth_axis_test_nc(paths[1], c(0, 50, 100), c(11, 22, 33))
  make_depth_axis_test_nc(paths[2], 50, 22)
  x <- data.frame(lon = 0, lat = 40, depth = 200)
  outside <- capture_selection_reports(extract3d_nearest(
    x, paths[1], var = "temp", date_col = NULL
  ))
  expect_match(outside$text, "outside the layer range", fixed = TRUE)
  expect_match(outside$text, "Nearest edge layers:", fixed = TRUE)
  expect_match(outside$text, "50 m apart", fixed = TRUE)
  expect_false(grepl("no valid value", outside$text, fixed = TRUE))
  singleton <- capture_selection_reports(extract3d_nearest(
    x, paths[2], var = "temp", date_col = NULL
  ))
  expect_match(singleton$text, "Only one depth layer", fixed = TRUE)
  expect_match(singleton$text, "spacing cannot be assessed", fixed = TRUE)
  expect_false(grepl("Showing", singleton$text, fixed = TRUE))
})

test_that("depth reference follows the source that supplies the value", {
  skip_if_not_installed("ncdf4")
  paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_depth_axis_test_nc(paths[1], c(0, 50, 100), rep(NA_real_, 3))
  make_depth_axis_test_nc(paths[2], c(55, 67, 150), c(NA, NA, 33))
  for (sources in list(paths, rev(paths))) {
    result <- capture_selection_reports(extract3d_nearest(
      data.frame(lon = 0, lat = 40, depth = 60), sources,
      var = "temp", date_col = NULL
    ))
    expect_equal(result$value$temp, 33)
    expect_match(result$text, "12 m apart", fixed = TRUE)
    expect_match(result$text, basename(paths[2]), fixed = TRUE)
    expect_false(grepl(basename(paths[1]), result$text, fixed = TRUE))
  }
})

test_that("quiet value tests preserve unexpected conditions", {
  expect_message(without_selection_reports(message("other message")),
                 "other message", fixed = TRUE)
  expect_warning(without_selection_reports(warning("other warning")),
                 "other warning", fixed = TRUE)
  expect_error(without_selection_reports(stop("other error")),
               "other error", fixed = TRUE)
})

test_that("great-circle distances match sf under s2", {
  original <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(original)), add = TRUE)
  suppressMessages(sf::sf_use_s2(TRUE))
  a <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326)
  points <- rbind(c(1, 1), c(179.5, -60), c(-179.9, 89), c(0.01, 0))
  for (k in seq_len(nrow(points))) {
    b <- sf::st_sfc(sf::st_point(points[k, ]), crs = 4326)
    expect_equal(
      .great_circle_km(0, 0, points[k, 1], points[k, 2]),
      as.numeric(sf::st_distance(a, b)) / 1000,
      tolerance = 1e-12
    )
  }
})

test_that("matching does not depend on the sf s2 setting", {
  skip_if_not_installed("ncdf4")
  original <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(original)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_horizontal_test_nc(path, c(0, 1), c(40, 41), 1:4)
  x <- data.frame(lon = c(0.1, 1.4), lat = 40)
  expect_equal(without_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL
  ))$temp, c(1, 2))
  expect_equal(without_selection_reports(extract2d(
    x, path, var = "temp", date_col = NULL,
    spatial_match = "nearest", max_distance_km = 100
  ))$temp, c(1, 2))
})

test_that("whole-grid nearest search matches a search over every cell", {
  set.seed(42)
  lon <- c(seq(-179.5, 179.5, by = 7), 180)
  lat <- sort(c(-90, runif(20, -89, 89), 90))
  x <- c(runif(60, -360, 360), 180, -180, 179.99, 0, 10)
  y <- c(runif(60, -90, 90), 0, 89.9, -89.9, 90, -90)
  fast <- .nearest_horizontal_cells(lon, lat, x, y)
  for (k in seq_along(x)) {
    slow <- .nearest_horizontal_cell(lon, lat, seq_along(lon),
                                     seq_along(lat), x[k], y[k])
    expect_equal(fast[[k]]$distance_km, slow$distance_km, tolerance = 1e-12)
    # Away from the poles the selected cell is the same too.
    if (abs(y[k]) < 89) {
      expect_identical(fast[[k]][c("lon_idx", "lat_idx")],
                       slow[c("lon_idx", "lat_idx")])
    }
  }
  invalid <- .nearest_horizontal_cells(lon, lat, c(NA, 0, Inf), c(0, 91, 0))
  expect_true(all(vapply(invalid, is.null, logical(1))))
  expect_length(.nearest_horizontal_cells(lon, lat, numeric(), numeric()), 0L)
})

test_that("exact distance ties keep the first cell in storage order", {
  batch <- .nearest_horizontal_cells(c(0, 0), 0, c(0, 0), c(0, 0))
  expect_identical(batch[[1L]]$lon_idx, 1L)
  expect_identical(batch[[2L]]$lon_idx, 1L)
  expect_equal(batch[[1L]]$distance_km, 0)
  # A point on a shared edge is equally far from both centres; rounding
  # noise must not decide, so the first (west/south) cell is kept.
  lon <- seq(-4.5, 4.5, 1)
  lat <- seq(35.5, 44.5, 1)
  for (y in lat) {
    edge <- .nearest_horizontal_cell(lon, lat, 4:5, which(lat == y), -1, y)
    expect_identical(edge$lon_idx, 4L)
    nearest <- .nearest_horizontal_cells(lon, lat, -1, y)[[1L]]
    expect_identical(nearest$lon_idx, 4L)
  }
  # A point just off the edge still gets the truly closer cell.
  off <- .nearest_horizontal_cell(lon, lat, 4:5, 1L, -1 + 1e-6, 35.5)
  expect_identical(off$lon_idx, 5L)
})

test_that("nearest candidates preserve matching limits and containment", {
  lon <- list(centres = c(0, 1), lower = c(-0.5, 0.5),
              upper = c(0.5, 1.5), longitude = TRUE)
  lat <- list(centres = 40, lower = 39, upper = 41, longitude = FALSE)
  x <- c(0, 0.1, 1.6, NA)
  y <- rep(40, 4)
  nearest <- .nearest_horizontal_cells(c(0, 1), 40, x, y)
  distance <- nearest[[3L]]$distance_km
  for (limit in c(0, 1, 100, distance, distance - 1e-6)) {
    matches <- .match_horizontal_cells(lon, lat, c(0, 1), 40, x, y,
                                       "nearest", limit)
    for (i in 1:3) {
      expect_identical(matches$lon_idx[i], nearest[[i]]$lon_idx)
      expect_identical(matches$accepted[i],
                       nearest[[i]]$distance_km <= limit)
    }
    expect_identical(matches$reason[4L], "invalid_coordinates")
  }
  expect_false(matches$inside_coverage[3L])
  expect_true(matches$inside_selected_cell[2L])
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

test_that("spatial rejection text explains nearly identical distances", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_horizontal_test_nc(path, 0, 0, 11)
  distance <- .nearest_horizontal_cells(0, 0, 0.01, 0)[[1L]]$distance_km
  limit <- distance - 1e-6
  result <- capture_selection_reports(extract2d(
    data.frame(lon = 0.01, lat = 0), path, var = "temp", date_col = NULL,
    spatial_match = "nearest", max_distance_km = limit
  ))
  expect_true(is.na(result$value$temp))
  comparison <- .selection_comparison(distance, limit)
  expect_false(identical(comparison$value, comparison$limit))
  expect_match(result$text, paste0("Nearest cell: ", comparison$value),
               fixed = TRUE)
  expect_match(result$text, paste0("allowed: ", comparison$limit), fixed = TRUE)
})

test_that("depth threshold text distinguishes a small excess", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(0, 50, 100.001), c(NA, NA, 33))
  result <- capture_selection_reports(extract3d_nearest(
    data.frame(lon = 0, lat = 40, depth = 50), path,
    var = "temp", date_col = NULL
  ))
  expect_match(result$text, "50.001 m deeper", fixed = TRUE)
  expect_match(result$text, "50 m apart", fixed = TRUE)
  expect_match(result$text, "exceeds this spacing", fixed = TRUE)
})

# PR #81 review fixes ----------------------------------------------------------

test_that("observation depths must be metres, positive downward", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_depth_axis_test_nc(path, c(0, -50, -100), c(10, 20, 30),
                          positive = "up")
  expect_error(
    extract3d_nearest(data.frame(lon = 0, lat = 40, depth = -50), path,
                      var = "temp", date_col = NULL),
    "positive downward"
  )
  expect_error(
    extract3d_all(data.frame(lon = 0, lat = 40, depth = "50"), path,
                  var = "temp", date_col = NULL),
    "positive downward"
  )
  expect_equal(without_selection_reports(extract3d_nearest(
    data.frame(lon = 0, lat = 40, depth = c(50, NA)), path,
    var = "temp", date_col = NULL
  ))$temp, c(20, NA))
})

test_that("full spellings of km and cm depth units are converted", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  for (units in c("centimeters", "centimetres", "Centimeter")) {
    expect_equal(.prepare_depth_axis(c(0, 500), units)$depth_m, c(0, 5))
  }
  for (units in c("kilometers", "kilometres", "kilometer")) {
    expect_equal(.prepare_depth_axis(c(0, 0.5), units)$depth_m, c(0, 500))
  }
  make_depth_axis_test_nc(path, c(500, 1500, 2500), c(11, 22, 33),
                          units = "centimeters", positive = "down")
  expect_equal(without_selection_reports(extract3d_nearest(
    data.frame(lon = 0, lat = 40, depth = 14), path,
    var = "temp", date_col = NULL
  ))$temp, 22)
})

test_that("height axes without a positive attribute are not depths", {
  skip_if_not_installed("ncdf4")
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
  skip_if_not_installed("ncdf4")
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

test_that("single-precision station coordinates match their own station", {
  station <- as.numeric(c(-122.40000152587891, 36.79999923706055))
  lon <- list(centres = station[1], lower = station[1], upper = station[1],
              longitude = TRUE, basis = "singleton_coordinate")
  lat <- list(centres = station[2], lower = station[2], upper = station[2],
              longitude = FALSE, basis = "singleton_coordinate")
  expect_identical(.containing_axis_indices(lon, -122.4), 1L)
  expect_identical(.containing_axis_indices(lat, 36.8), 1L)
  expect_length(.containing_axis_indices(lon, -122.401), 0L)
})

test_that("unusable spatial bounds fall back to inferred edges", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_horizontal_test_nc(path, c(0, 1), c(40, 41), 1:4)
  nc <- ncdf4::nc_open(path, write = TRUE)
  ncdf4::ncatt_put(nc, "lon", "bounds", "lon_bnds")
  ncdf4::nc_close(nc)
  x <- data.frame(lon = c(0, 1, 5), lat = 40)
  expect_warning(
    out <- without_selection_reports(extract2d(
      x, path, var = "temp", date_col = NULL
    )),
    "were not found"
  )
  expect_equal(out$temp, c(1, 2, NA))
  expect_true(.same_spatial_units("degree_east", "degrees_east"))
  expect_true(.same_spatial_units("degrees_E", "degrees"))
  expect_false(.same_spatial_units("degrees_east", "radians"))
})

test_that("within matching keeps input order between overlapping sources", {
  without_selection_reports({
    skip_if_not_installed("ncdf4")
    paths <- replicate(2, tempfile(fileext = ".nc"))
    on.exit(unlink(paths), add = TRUE)
    fine <- seq(0.05, 0.95, by = 0.1)
    make_horizontal_test_nc(paths[1], fine, fine + 40, rep(100, 100))
    make_horizontal_test_nc(paths[2], c(0, 1), c(40, 41), rep(1, 4))
    x <- data.frame(lon = 0.01, lat = 40.01)
    expect_equal(extract2d(x, paths, var = "temp", date_col = NULL)$temp, 100)
    expect_equal(extract2d(x, rev(paths), var = "temp",
                           date_col = NULL)$temp, 1)
    # Nearest matching prefers the closer centre, whatever the order.
    for (files in list(paths, rev(paths))) {
      expect_equal(extract2d(
        x, files, var = "temp", date_col = NULL,
        spatial_match = "nearest", max_distance_km = 50
      )$temp, 1)
    }
  })
})

test_that("off-site tiles do not turn a temporal miss into a spatial one", {
  skip_if_not_installed("ncdf4")
  paths <- replicate(2, tempfile(fileext = ".nc"))
  on.exit(unlink(paths), add = TRUE)
  make_horizontal_test_nc(paths[1], c(10, 11), c(40, 41), 1:8, days = c(0, 4))
  make_horizontal_test_nc(paths[2], 0, 40, 5, days = 0)
  x <- data.frame(lon = 0, lat = 40, date = as.Date("2020-01-05"))
  result <- capture_selection_reports(extract2d(
    x, paths, var = "temp", max_time_diff = 1, spatial_report = FALSE
  ))
  expect_true(is.na(result$value$temp))
  expect_length(result$reports, 1L)
  expect_s3_class(result$reports[[1L]], "ocean3d_time_report")
  expect_match(result$text, "Observation: 2020-01-05", fixed = TRUE)
  expect_match(result$text, "No timestamp satisfied", fixed = TRUE)
  expect_match(result$text, "Point outside cell coverage", fixed = TRUE)
  expect_identical(result$reports[[1L]]$report$status, "no_value")
})

test_that("repeated missing reasons are grouped across many files", {
  skip_if_not_installed("ncdf4")
  paths <- vapply(1:8, function(k) tempfile(fileext = ".nc"), character(1))
  on.exit(unlink(paths), add = TRUE)
  for (k in seq_along(paths)) {
    make_horizontal_test_nc(paths[k], 0, 40, k, days = k - 1L)
  }
  x <- data.frame(lon = 0, lat = 40, date = as.Date("2020-03-01"))
  result <- capture_selection_reports(extract2d(
    x, paths, var = "temp", max_time_diff = 1
  ))
  expect_length(result$reports, 1L)
  expect_identical(lengths(regmatches(
    result$text, gregexpr("No timestamp satisfied", result$text)
  )), 1L)
  expect_match(result$text, "and 5 more", fixed = TRUE)
})

test_that("a missing date gets its own temporal reason", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_horizontal_test_nc(path, 0, 40, 1:3, days = 0:2)
  result <- capture_selection_reports(extract2d(
    data.frame(lon = 0, lat = 40), path, var = "temp", date_col = NULL
  ))
  expect_true(is.na(result$value$temp))
  expect_length(result$reports, 1L)
  expect_match(result$text, "No observation date was supplied", fixed = TRUE)
  quiet <- capture_selection_reports(extract2d(
    data.frame(lon = 0, lat = 40), path, var = "temp", date_col = NULL,
    time_report = FALSE
  ))
  expect_length(quiet$reports, 0L)
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
