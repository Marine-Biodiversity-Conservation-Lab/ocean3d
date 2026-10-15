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

  if (coordinate_var) ncdf4::ncvar_put(nc, "longitude", lon_vals)
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


# Public entry points and netCDF helpers ----------------------------------------

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
  expect_error(ocean3d:::.as_netcdf_sources(42), "Unsupported")
  expect_error(ocean3d:::.as_netcdf_sources(list(42)), "Unsupported")
  expect_error(
    ocean3d:::.open_source(list(type = "path", source = tempfile())),
    "netCDF file does not exist"
  )
})

test_that("extraction leaves caller-owned netCDF connections open", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  opened <- ocean3d:::.open_source(list(type = "connection", source = nc))
  expect_identical(opened$nc, nc)
  expect_false(opened$close)
  out <- extract2d(data.frame(lon = 0, lat = 40), nc,
                   var = "temp", date_col = NULL)
  expect_equal(out$temp, 1)
  expect_equal(ncdf4::ncvar_get(nc, "temp"), matrix(1:4, 2, 2))
})

test_that("dimension helpers read coordinate variables and report missing names", {
  expect_error(
    ocean3d:::.match_name(c("x", "y"), "lon"),
    "Available names are: x, y"
  )
  expect_null(ocean3d:::.match_name("x", "lon", required = FALSE))
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, coordinate_var = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_equal(as.vector(ocean3d:::.get_dim_values(nc, "longitude")), c(0, 1))
  expect_null(ocean3d:::.get_dim_values(nc, NULL))
  expect_error(ocean3d:::.get_dim_values(nc, "missing"),
               "Could not read dimension or coordinate variable")
  expect_error(ocean3d:::.detect_dims(nc, "missing"),
               "Available variables are: temp, longitude")
  expect_equal(
    ocean3d:::.build_start_count(
      c("lon", "lat"), list(lon = list(start = 2L, count = 1L))
    ), list(start = c(2L, 1L), count = c(1L, -1L))
  )
})

test_that("longitude conversion and missing coordinates are handled", {
  expect_equal(ocean3d:::.prepare_lon(c(0, 359), -1), 359)
  expect_true(is.na(ocean3d:::.prepare_lon(c(0, 359), NA_real_)))
  expect_identical(ocean3d:::.nearest_index(numeric(), 0), NA_integer_)
  expect_identical(ocean3d:::.nearest_index(c(NA_real_, NA_real_), 0),
                   NA_integer_)
  expect_identical(ocean3d:::.nearest_index(c(0, 1), NA_real_), NA_integer_)
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE,
                   lon_vals = c(0, 359))
  x <- data.frame(lon = c(-1, NA, 0), lat = c(40, 40, NA))
  out <- extract2d(x, f, var = "temp", date_col = NULL)
  expect_equal(out$temp, c(2, NA, NA))
  expect_equal(nrow(out), 3L)
})

test_that("time helpers convert seconds minutes and hours to days", {
  paths <- character()
  on.exit(unlink(paths), add = TRUE)
  for (unit in c("seconds", "minutes", "hours")) {
    f <- tempfile(fileext = ".nc")
    paths <- c(paths, f)
    make_test_netcdf(f, time_units = paste(unit, "since 2020-01-01"))
    nc <- ncdf4::nc_open(f)
    tryCatch({
      per_day <- switch(unit, seconds = 86400, minutes = 1440, hours = 24)
      info <- ocean3d:::.get_time_origin(nc, "time")
      expect_equal(info$origin, "2020-01-01")
      expect_equal(info$multiplier, 1 / per_day)
      expect_equal(
        c(ocean3d:::.convert_time_values(nc, "time")),
        as.Date(c(0, 2, 5) / per_day, origin = "2020-01-01")
      )
      expect_null(ocean3d:::.convert_time_values(nc, NULL))
    }, finally = ncdf4::nc_close(nc))
  }
})

test_that("time without units falls back to days since the Unix epoch", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, time_units = NULL)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_null(ocean3d:::.get_time_origin(nc, "time"))
  expect_equal(c(ocean3d:::.convert_time_values(nc, "time")),
               as.Date(1:3, origin = "1970-01-01"))
})

test_that("output prefixes respect automatic and explicit naming", {
  for (method in c("surface", "bottom", "nearest")) {
    prefix <- switch(method, surface = "surface_", bottom = "seabottom_",
                     nearest = "nearest_")
    expect_equal(
      ocean3d:::.resolve_output_prefix(method, "temp", NULL, TRUE),
      paste0(prefix, "temp")
    )
    expect_null(ocean3d:::.resolve_output_prefix(method, "temp", NULL, FALSE))
    expect_equal(
      ocean3d:::.resolve_output_prefix(method, "temp", "custom", TRUE),
      "custom"
    )
  }
})
test_that("multiple files fill missing values and preserve first-file priority", {
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
    result_cols <- setdiff(names(out), names(make_obs()))
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
    expect_equal(names(out), c("lon", "lat", "chosen"))
    expect_equal(out$chosen, 1)
  }
})

test_that("named nested sources accept paths and caller-owned connections", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, with_depth = FALSE, with_time = FALSE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  x <- data.frame(lon = c(0, 1), lat = c(40, 41))
  expected <- extract2d(x, f, var = "temp", date_col = NULL)
  nested <- list(year = list(path = f, opened = list(connection = nc)))
  out <- extract2d(x, nested, var = "temp", date_col = NULL)
  expect_equal(out, expected)
  expect_equal(ncdf4::ncvar_get(nc, "temp"), matrix(1:4, 2, 2))
})

test_that("direct coordinate vectors return all depth summaries as a data frame", {
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

test_that("max_time_diff rejects non-scalar non-numeric and missing limits", {
  x <- data.frame(lon = 0, lat = 40)
  for (bad in list("1", c(1, 2), numeric(), NA_real_)) {
    expect_error(extract_to_point(
      x, "unused.nc", var = "temp", method = "2d",
      date_col = NULL, max_time_diff = bad
    ), "non-negative numeric value")
  }
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

test_that("exact temporal matching selects known values and skips missing dates", {
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
        ocean3d:::.match_time_index(
          values, as.Date("2020-01-03"), time_match = method
        ),
        NA_integer_
      )
    }
  }
})

test_that("public depth wrappers automatically name columns and extract values", {
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
    expect_equal(names(out), c(names(x), case$col))
    expect_equal(out[[case$col]], case$value)
    expect_equal(out[names(x)], x)
  }
})

test_that("extract2d automatically names the detected variable with time matching", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  make_test_netcdf(f, var = "oxygen", with_depth = FALSE)
  x <- data.frame(lon = c(0, 1), lat = c(40, 41),
                  date = as.Date(c("2020-01-01", "2020-01-03")))
  out <- extract2d(x, f, time_match = "exact")
  expect_equal(names(out), c(names(x), "oxygen"))
  expect_equal(out$oxygen, c(1, 8))
})

test_that("nested data frame sources read every file using a custom file column", {
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