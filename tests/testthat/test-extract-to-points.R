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
  x <- data.frame(lon = 0, lat = 40)
  for (bad in list(-1, "1", c(1, 2), numeric(), NA_real_)) {
    expect_error(
      extract_to_point(
        x, "unused.nc", var = "temp", method = "2d",
        date_col = NULL, max_time_diff = bad
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


# Public entry points and netCDF helpers ---------------------------------------

test_that("extract2d preserves observations and extracts known cells", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  x <- data.frame(id = c("a", "b"), lon = c(0, 1), lat = c(40, 41))
  out <- extract2d(x, f, var = "temp", date_col = NULL)
  expect_equal(out[names(x)], x)
  expect_equal(out$temp, c(1, 4))
})

test_that("variable and output column names must be single strings", {
  x <- data.frame(lon = 0, lat = 40)
  for (bad in list(42, c("a", "b"), character(), NA_character_)) {
    expect_error(extract_to_point(
      x, "unused.nc", var = bad, method = "2d", date_col = NULL
    ), "var.*must be")
    expect_error(extract_to_point(
      x, "unused.nc", var = "temp", output_col = bad,
      method = "2d", date_col = NULL
    ), "output_col.*must be")
  }
})

test_that("depth methods reject files without a depth axis", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  for (method in c("surface", "bottom", "nearest", "all")) {
    expect_error(extract_to_point(
      make_obs(), f, var = "temp", method = method, date_col = NULL
    ), "requires a depth dimension")
  }
})

test_that("a single direct observation returns named depth summaries", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_time = FALSE)
  out <- extract_to_point(nc = f, lon = 0, lat = 40, depth = 60,
                          var = "temp", method = "all")
  expect_equal(out, c(nearest_temp = 5, surface_temp = 1, seabottom_temp = 9))
})

test_that("verbose extraction reports progress without changing results", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  x <- data.frame(lon = rep(0, 100), lat = rep(40, 100))
  expect_message(
    out <- extract2d(x, f, var = "temp", date_col = NULL, verbose = TRUE),
    "Extracting row 100 of 100"
  )
  expect_equal(out$temp, rep(1, 100))
  expect_silent(extract2d(x, f, var = "temp", date_col = NULL))
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
                   var = "temp", date_col = NULL)
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
  out <- extract2d(x, f, var = "temp", date_col = NULL)
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
  out <- extract2d(x, paths, var = "temp", date_col = NULL)
  expect_equal(out$temp, c(11, 4))
  reverse <- extract2d(x, rev(paths), var = "temp", date_col = NULL)
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
                            date_col = NULL)
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
  nearest <- extract3d_nearest(x, f, var = "temp", date_col = NULL)
  all <- extract3d_all(x, f, var = "temp", date_col = NULL)
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
                     output_prefix = "ignored", output_col = "chosen")
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
  out <- extract2d(x, nested, var = "temp", date_col = NULL)

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
    var = "temp", method = "all"
  )
  expect_equal(out, data.frame(
    nearest_temp = c(5, 12), surface_temp = c(1, 4),
    seabottom_temp = c(9, 12)
  ))
})


test_that("id_col must exist and is preserved when supplied", {
  x <- data.frame(lon = 0, lat = 40)
  expect_error(extract2d(
    x, "unused.nc", var = "temp", date_col = NULL, id_col = "sample"
  ), "Missing required column.*sample")
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  x$sample <- "station-a"
  out <- extract2d(x, f, var = "temp", date_col = NULL, id_col = "sample")
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
      var = "temp", method = "nearest", time_match = "exact"
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
  exact <- extract3d_surface(x, f, var = "temp", time_match = "exact")
  nearest <- extract3d_surface(x, f, var = "temp", time_match = "nearest")
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
  out <- extract2d(x, f, time_match = "exact")
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
  out <- extract3d_surface(x, sources, file_col = "nc_path")
  expect_equal(out[names(x)], x)
  expect_equal(out$surface_temp, 13)
  expect_equal(out$surface_salt, 13)
  expect_equal(out, extract3d_surface(x, paths))
  expect_error(
    extract3d_surface(x, list(data.frame(wrong = paths)),
                      file_col = "nc_path"),
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
    time_match = "exact"
  )
  expect_equal(daily, 5)

  table_input <- extract2d(
    data.frame(lon = 0, lat = 40, date = target),
    f, var = "temp",
    time_resolution = "daily",
    time_match = "exact"
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
    var = "temp", method = "2d"
  )
  expect_equal(hourly, 9)

  expect_error(
    extract_to_point(
      nc = f, lon = 0, lat = 40, date = "2020-01-01",
      var = "temp", method = "2d"
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
      max_time_diff = 0
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
    extract2d(x, f, var = "temp", date_col = NULL),
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
    x, f, var = "temp", date_col = NULL
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
    x, f, var = "static", date_col = NULL
  )
  
  expect_equal(out$static, c(11, 14))
  expect_equal(nrow(out), 2L)
})
