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
  # Altitude is positive up, so -10 m is 10 m deep.
  out <- extract_to_point(
    nc = f, var = "v", lon = 0, lat = 40, depth = 10, method = "nearest"
  )
  expect_identical(as.numeric(out), 5)
})

test_that("point extraction reads the cell CF metadata points to", {
  f <- make_swapped_axes_file()
  on.exit(unlink(f), add = TRUE)

  # lon = 40 is the second `y`, lat = 10 the first `x`: v[1, 2] = 3.
  out <- extract_to_point(nc = f, var = "v", lon = 40, lat = 10)
  expect_identical(as.numeric(out), 3)
})

test_that("2D geographic coordinates are rejected by units or standard_name", {
  for (standard_name in c(FALSE, TRUE)) {
    f <- make_curvilinear_file(standard_name = standard_name)
    on.exit(unlink(f), add = TRUE)

    expect_error(detect(f), class = "ocean3d_error_curvilinear_grid")
    err <- tryCatch(detect(f), error = conditionMessage)
    expect_match(err, "two-dimensional geographic coordinate `nav_lon`",
                 fixed = TRUE)
    expect_match(err, "regular longitude/latitude grid", fixed = TRUE)
  }
})

test_that("2D geographic coordinates inside a group are rejected", {
  f <- make_curvilinear_file(group = "g")
  on.exit(unlink(f), add = TRUE)

  expect_error(detect(f, var = "g/v"),
               class = "ocean3d_error_curvilinear_grid")
})

test_that("point extraction rejects a units-only curvilinear grid", {
  f <- make_curvilinear_file()
  on.exit(unlink(f), add = TRUE)

  expect_error(
    extract_to_point(nc = f, var = "v", lon = 10, lat = 40),
    class = "ocean3d_error_curvilinear_grid"
  )
})

test_that("an explicit dimension does not bypass the curvilinear check", {
  f <- make_curvilinear_file()
  on.exit(unlink(f), add = TRUE)

  expect_error(
    detect(f, lon_dim = "x", lat_dim = "y"),
    class = "ocean3d_error_curvilinear_grid"
  )
})
