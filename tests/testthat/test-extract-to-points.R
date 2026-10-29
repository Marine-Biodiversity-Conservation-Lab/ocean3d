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
  expect_equal(.prepare_lon(c(0, 359), -1), 359)
  expect_true(is.na(.prepare_lon(c(0, 359), NA_real_)))
  expect_identical(.nearest_index(numeric(), 0), NA_integer_)
  expect_identical(.nearest_index(c(NA_real_, NA_real_), 0),
                   NA_integer_)
  expect_identical(.nearest_index(c(0, 1), NA_real_), NA_integer_)
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

test_that("missing target depth preserves surface and bottom summaries", {
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
    time_index = c(1L, 2L), variable = "temp"
  )
  attr(b, ".source_candidate") <- list(
    distance_days = c(1, 9), timestamp = as.numeric(rev(stamp)),
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
  expect_identical(.spatial_axis_index(axis, -0.5), 1L)
  expect_identical(.spatial_axis_index(axis, 1.5), 2L)
  expect_true(is.na(.spatial_axis_index(axis, 1.6)))
  axis <- list(
    centres = c(179, 181), lower = c(178, 180), upper = c(180, 182),
    longitude = TRUE, basis = "inferred_regular"
  )
  expect_identical(.spatial_axis_index(axis, -179), 2L)
  expect_identical(.spatial_axis_index(axis, 181), 2L)
  expect_true(is.na(.spatial_axis_index(axis, 0)))
})

test_that("singleton axes do not invent spatial coverage", {
  axis <- list(
    centres = 40, lower = 40, upper = 40,
    longitude = FALSE, basis = "singleton_coordinate"
  )
  expect_identical(.spatial_axis_index(axis, 40), 1L)
  expect_true(is.na(.spatial_axis_index(axis, 41)))
})

test_that("spatial eligibility selects the correct tile in either order", {
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
  expect_identical(.spatial_axis_index(axis, 0.1), 1L)
  expect_true(is.na(.spatial_axis_index(axis, 1)))
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
  expect_identical(.spatial_axis_index(axis, 5), 1L)
  expect_identical(.spatial_axis_index(axis, 0), 3L)
  expect_true(is.na(.spatial_axis_index(axis, 20)))
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
    "No source supplied a valid value", fixed = TRUE
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
