# Helpers ----------------------------------------------------------------------

make_test_netcdf <- function(path, var = "temp", with_depth = TRUE,
                             with_time = TRUE, bottom_na = FALSE) {
  lon <- ncdf4::ncdim_def("lon", "degrees_east", vals = c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", vals = c(40, 41))
  dims <- list(lon, lat)

  if (with_depth) {
    depth <- ncdf4::ncdim_def("depth", "m", vals = c(0, 50, 100))
    dims <- c(dims, list(depth))
  }

  if (with_time) {
    time <- ncdf4::ncdim_def(
      "time", "days since 2020-01-01",
      vals = c(0, 2, 5), unlim = FALSE
    )
    dims <- c(dims, list(time))
  }

  v <- ncdf4::ncvar_def(var, "1", dim = dims, missval = -9999, prec = "double")
  nc <- ncdf4::nc_create(path, v)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  shape <- vapply(dims, function(x) length(x$vals), integer(1))
  values <- array(seq_len(prod(shape)), dim = shape)

  if (with_depth && bottom_na) {
    if (with_time) values[1, 1, 3, ] <- NA_real_
    else values[1, 1, 3] <- NA_real_
  }

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
    date_col = NULL
  )

  expect_true("temp" %in% names(out))
  expect_identical(out$temp, 1)
})

test_that("extract_to_point detects a single variable automatically", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, var = "oxygen")

  out <- extract3d_surface(make_obs(), f, var = NULL)

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
      f, var = NULL, method = "2d", date_col = NULL
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

  surface <- extract3d_surface(x, f, var = "temp")
  nearest <- extract3d_nearest(x, f, var = "temp")
  bottom <- extract3d_bottom(x, f, var = "temp")
  all <- extract3d_all(x, f, var = "temp")

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

  out <- extract3d_nearest(make_obs(depth = 100), f, var = "temp")

  expect_false(is.na(out$temp))
})

# Time matching ----------------------------------------------------------------

test_that("time_match = nearest uses the closest available netCDF date", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  out <- extract3d_surface(
    make_obs(date = as.Date("2020-01-02")),
    f, var = "temp", time_match = "nearest"
  )

  expect_false(is.na(out$temp))
})

test_that("time_match = exact returns NA when the date is absent", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  out <- extract3d_surface(
    make_obs(date = as.Date("2020-01-02")),
    f, var = "temp", time_match = "exact"
  )

  expect_true(is.na(out$temp))
})

test_that("max_time_diff limits nearest temporal matching", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  accepted <- extract3d_surface(
    make_obs(date = as.Date("2020-01-04")),
    f, var = "temp", time_match = "nearest", max_time_diff = 1
  )
  rejected <- extract3d_surface(
    make_obs(date = as.Date("2020-01-10")),
    f, var = "temp", time_match = "nearest", max_time_diff = 1
  )

  expect_false(is.na(accepted$temp))
  expect_true(is.na(rejected$temp))
})

test_that("max_time_diff rejects invalid values", {
  expect_error(
    extract_to_point(
      data.frame(lon = 0, lat = 40),
      nc = "unused.nc", var = "temp", method = "2d",
      date_col = NULL, max_time_diff = -1
    ),
    "non-negative numeric value"
  )
})

# Input forms and validation ---------------------------------------------------

test_that("multiple netCDF files add one output column per detected variable", {
  skip_if_not_installed("ncdf4")
  f1 <- tempfile(fileext = ".nc")
  f2 <- tempfile(fileext = ".nc")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  make_test_netcdf(f1, var = "temp")
  make_test_netcdf(f2, var = "salt")

  out <- extract3d_bottom(make_obs(), c(f1, f2), var = NULL)

  expect_true(all(c("seabottom_temp", "seabottom_salt") %in% names(out)))
})

test_that("netCDF input can be a list, opened connection, or data frame", {
  skip_if_not_installed("ncdf4")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f)

  from_list <- extract3d_surface(make_obs(), list(f), var = "temp")

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  from_connection <- extract3d_surface(make_obs(), nc, var = "temp")

  from_df <- extract3d_surface(
    make_obs(), data.frame(file = f), var = "temp", file_col = "file"
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
      f, var = "temp"
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
    date_col = NULL
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
      date_col = NULL
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
      date_col = NULL
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
      date_col = NULL
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
      method = "nearest"
    ),
    "`date` could not be converted to a valid date"
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
    date_col = NULL
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
    date_col = NULL
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
    date_col = NULL
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
    date_col = NULL
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
    date_col = NULL
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
    date_col = NULL
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
    method = "nearest"
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
      method = "nearest"
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
    date_col = NULL
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
    method = "nearest"
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
      date_col = NULL
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
      date_col = NULL
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
      date_col = NULL
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
      date_col = NULL
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
      date_col = NULL
    ),
    "`data` must be a data frame"
  )
})

