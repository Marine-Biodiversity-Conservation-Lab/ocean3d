# Axis detection from CF metadata (#93), tested on the shared helpers
# directly: the core object (#113) reuses `.detect_dims()`, while the
# extraction functions that call it may be retired (#116; AGENTS.md, "Tests").
# The files use ordinary dimension coordinates made with ncdim_def(), which
# ncdf4 lists in `nc$dim` rather than `nc$var`.

# A 2 x 2 file whose dimensions are named `x` and `y` and labelled the other
# way round by their `standard_name` and `axis` attributes. `group` puts the
# variable and its dimensions in a netCDF-4 group.
make_swapped_axes_file <- function(group = NULL) {
  name <- function(x) if (is.null(group)) x else paste0(group, "/", x)
  f <- tempfile(fileext = ".nc")
  x <- ncdf4::ncdim_def(name("x"), "degrees", c(10, 20))
  y <- ncdf4::ncdim_def(name("y"), "degrees", c(30, 40))
  nc <- ncdf4::nc_create(
    f, ncdf4::ncvar_def(name("v"), "1", list(x, y), -999),
    force_v4 = !is.null(group)
  )
  ncdf4::ncvar_put(nc, name("v"), matrix(1:4, 2))
  ncdf4::ncatt_put(nc, name("x"), "standard_name", "latitude")
  ncdf4::ncatt_put(nc, name("x"), "axis", "Y")
  ncdf4::ncatt_put(nc, name("y"), "standard_name", "longitude")
  ncdf4::ncatt_put(nc, name("y"), "axis", "X")
  ncdf4::nc_close(nc)
  f
}

# A NEMO-style file: index dimensions `x` and `y`, and 2D `nav_lon` /
# `nav_lat` coordinates named by the variable's `coordinates` attribute.
# `standard_name = FALSE` leaves the coordinates with units only. The
# `coordinates` attribute holds names relative to the variable's group.
make_curvilinear_file <- function(standard_name = FALSE, group = NULL) {
  name <- function(x) if (is.null(group)) x else paste0(group, "/", x)
  f <- tempfile(fileext = ".nc")
  x <- ncdf4::ncdim_def(name("x"), "", 1:2)
  y <- ncdf4::ncdim_def(name("y"), "", 1:2)
  nc <- ncdf4::nc_create(f, list(
    ncdf4::ncvar_def(name("v"), "1", list(x, y), -999),
    ncdf4::ncvar_def(name("nav_lon"), "degrees_east", list(x, y), -999),
    ncdf4::ncvar_def(name("nav_lat"), "degrees_north", list(x, y), -999)
  ), force_v4 = !is.null(group))
  ncdf4::ncvar_put(nc, name("v"), matrix(1:4, 2))
  ncdf4::ncvar_put(nc, name("nav_lon"), matrix(c(10, 20, 10, 20), 2))
  ncdf4::ncvar_put(nc, name("nav_lat"), matrix(c(40, 40, 50, 50), 2))
  ncdf4::ncatt_put(nc, name("v"), "coordinates", "nav_lon nav_lat")
  if (standard_name) {
    ncdf4::ncatt_put(nc, name("nav_lon"), "standard_name", "longitude")
    ncdf4::ncatt_put(nc, name("nav_lat"), "standard_name", "latitude")
  }
  ncdf4::nc_close(nc)
  f
}

detect <- function(f, var = "v", ...) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc))
  .detect_dims(nc, var, ...)
}

test_that("axis and standard_name on dimension coordinates set the axes", {
  f <- make_swapped_axes_file()
  on.exit(unlink(f), add = TRUE)

  dims <- detect(f)
  expect_identical(dims$lon, "y")
  expect_identical(dims$lat, "x")
})

test_that("axis attributes are read on dimensions inside a group", {
  f <- make_swapped_axes_file(group = "g")
  on.exit(unlink(f), add = TRUE)

  dims <- detect(f, var = "g/v")
  expect_identical(dims$lon, "g/y")
  expect_identical(dims$lat, "g/x")
})

test_that("an axis attribute alone identifies a dimension", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  a <- ncdf4::ncdim_def("a", "1", c(0, 1))
  b <- ncdf4::ncdim_def("b", "1", c(40, 41))
  nc <- ncdf4::nc_create(f, ncdf4::ncvar_def("v", "1", list(a, b), -999))
  ncdf4::ncatt_put(nc, "a", "axis", "X")
  ncdf4::ncatt_put(nc, "b", "axis", "Y")
  ncdf4::nc_close(nc)

  dims <- detect(f)
  expect_identical(dims$lon, "a")
  expect_identical(dims$lat, "b")
})

test_that("an altitude dimension with axis = 'Z' is the depth axis", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  alt <- ncdf4::ncdim_def("altitude", "m", c(0, -10))
  nc <- ncdf4::nc_create(
    f, ncdf4::ncvar_def("v", "1", list(lon, lat, alt), -999)
  )
  ncdf4::ncvar_put(nc, "v", array(1:8, c(2, 2, 2)))
  ncdf4::ncatt_put(nc, "altitude", "standard_name", "altitude")
  ncdf4::ncatt_put(nc, "altitude", "axis", "Z")
  ncdf4::nc_close(nc)

  expect_identical(detect(f)$depth, "altitude")
})

test_that("2D geographic coordinates are rejected by units or standard_name", {
  for (standard_name in c(FALSE, TRUE)) {
    f <- make_curvilinear_file(standard_name = standard_name)
    on.exit(unlink(f), add = TRUE)

    err <- expect_error(detect(f), class = "ocean3d_error_curvilinear_grid")
    expect_match(conditionMessage(err),
                 "two-dimensional geographic coordinate `nav_lon`",
                 fixed = TRUE)
    expect_match(conditionMessage(err),
                 "regular longitude/latitude grid", fixed = TRUE)
  }
})

test_that("2D geographic coordinates inside a group are rejected", {
  f <- make_curvilinear_file(group = "g")
  on.exit(unlink(f), add = TRUE)

  expect_error(detect(f, var = "g/v"),
               class = "ocean3d_error_curvilinear_grid")
})

test_that("an explicit dimension does not bypass the curvilinear check", {
  f <- make_curvilinear_file()
  on.exit(unlink(f), add = TRUE)

  expect_error(
    detect(f, lon_dim = "x", lat_dim = "y"),
    class = "ocean3d_error_curvilinear_grid"
  )
})

# Dimension detection and reading (.detect_dims(), .get_dim_values()) ----------

test_that("dimension helpers read coordinates and report missing names", {
  f <- make_grid_file(lon = c(0, 1))
  on.exit(unlink(f), add = TRUE)
  # Add a `longitude` coordinate variable on the lon dimension.
  nc <- ncdf4::nc_open(f, write = TRUE)
  nc <- ncdf4::ncvar_add(nc, ncdf4::ncvar_def(
    "longitude", "degrees_east", list(nc$dim$lon), prec = "double"
  ))
  ncdf4::ncvar_put(nc, "longitude", c(10, 11))
  ncdf4::nc_close(nc)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_equal(as.vector(.get_dim_values(nc, "longitude")), c(10, 11))
  expect_null(.get_dim_values(nc, NULL))
  expect_error(.get_dim_values(nc, "missing"),
               "Could not read dimension or coordinate variable")
  expect_error(.detect_dims(nc, "missing"),
               "Available variables are: temp, longitude")

})

test_that("dimension detection uses only the selected variable", {
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
  f <- make_grid_file(depth = NULL, time = NULL)
  on.exit(unlink(f), add = TRUE)

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
