# Longitude axes with steps over 180 degrees (#95). The shared helpers are
# tested directly, because the core object (#113) reuses them, and point
# extraction checks the result end to end (AGENTS.md, "Tests"). Each file has
# a longitude axis, one latitude row at 40 N, the value 10 * i in longitude
# cell i, and optional 2 x n longitude bounds.

make_lon_file <- function(lons, bounds = NULL) {
  f <- tempfile(fileext = ".nc")
  lon <- ncdf4::ncdim_def("lon", "degrees_east", lons)
  lat <- ncdf4::ncdim_def("lat", "degrees_north", 40)
  vars <- list(ncdf4::ncvar_def(
    "temp", "1", list(lon, lat), -9999, prec = "double"
  ))
  if (!is.null(bounds)) {
    nv <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
    vars <- c(vars, list(ncdf4::ncvar_def(
      "lon_bounds", "degrees_east", list(nv, lon), prec = "double"
    )))
  }
  nc <- ncdf4::nc_create(f, vars)
  ncdf4::ncvar_put(nc, "temp", matrix(10 * seq_along(lons), ncol = 1L))
  if (!is.null(bounds)) {
    ncdf4::ncvar_put(nc, "lon_bounds", bounds)
    ncdf4::ncatt_put(nc, "lon", "bounds", "lon_bounds")
  }
  ncdf4::nc_close(nc)
  f
}

prepare_lon_axis <- function(f) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc))
  .prepare_spatial_axis(nc, "lon", .get_dim_values(nc, "lon"),
                        longitude = TRUE)
}

point_value <- function(files, lon) {
  extract_to_point(nc = files, lon = lon, lat = 40, var = "temp",
                   method = "2d")
}

test_that("monotonic longitude axes are left as they are", {
  expect_identical(.unwrap_longitude(c(0, 200, 300)), c(0, 200, 300))
  expect_identical(.unwrap_longitude(c(300, 200, 0)), c(300, 200, 0))
  expect_identical(.unwrap_longitude(c(0, 190, 200)), c(0, 190, 200))
  expect_identical(.unwrap_longitude(c(-170, 20, 170)), c(-170, 20, 170))
  expect_identical(.unwrap_longitude(c(0, 200)), c(0, 200))
  expect_identical(.unwrap_longitude(5), 5)
})

test_that("longitude axes crossing the seam are unwrapped", {
  expect_equal(.unwrap_longitude(c(350, 10, 30)), c(350, 370, 390))
  expect_equal(.unwrap_longitude(c(10, -10, 190)), c(10, -10, -170))
  # Steps over 180 degrees are kept when the axis crosses the seam.
  expect_equal(.unwrap_longitude(c(100, 300, 50)), c(100, 300, 410))
  expect_equal(.unwrap_longitude(c(50, 300, 100)), c(50, -60, -260))
})

test_that("longitude axes that cannot be made monotonic are rejected", {
  # Wraps more than once around the globe, in either direction.
  f <- make_lon_file(c(0, 170, 340, 150))
  on.exit(unlink(f))
  expect_identical(.unwrap_longitude(c(0, 170, 340, 150)),
                   c(0, 170, 340, 150))
  err <- expect_error(prepare_lon_axis(f),
                      class = "ocean3d_error_spatial_axis")
  expect_match(conditionMessage(err), "`lon` must form a monotonic axis",
               fixed = TRUE)
  expect_match(conditionMessage(err), "seam", fixed = TRUE)
  # A repeated longitude is not a valid axis either.
  expect_identical(.unwrap_longitude(c(0, 360, 10)), c(0, 360, 10))
})

test_that("coarse longitude cells use their explicit bounds", {
  f <- make_lon_file(
    c(0, 200, 300),
    bounds = matrix(c(-40, 100, 100, 250, 250, 310), nrow = 2L)
  )
  on.exit(unlink(f))
  axis <- prepare_lon_axis(f)
  expect_identical(axis$basis, "explicit_bounds")
  expect_identical(as.vector(axis$centres), c(0, 200, 300))
  expect_identical(.containing_axis_indices(axis, 120), 2L)
  expect_identical(.containing_axis_indices(axis, 260), 3L)
  expect_identical(.containing_axis_indices(axis, -30), 1L)
  expect_length(.containing_axis_indices(axis, 315), 0L)
})

test_that("point extraction reads coarse monotonic longitude axes", {
  axes <- list(c(0, 200, 300), c(300, 200, 0), c(0, 190, 200),
               c(-170, 20, 170))
  for (lons in axes) {
    f <- make_lon_file(lons)
    on.exit(unlink(f), add = TRUE)
    expect_identical(point_value(f, lons[3]), 30)
    expect_identical(point_value(f, lons[2]), 20)
    # The same file passed twice is read the same way.
    expect_identical(point_value(c(f, f), lons[3]), 30)
  }
})

test_that("point extraction uses explicit bounds of coarse longitude cells", {
  f <- make_lon_file(
    c(0, 200, 300),
    bounds = matrix(c(-40, 100, 100, 250, 250, 310), nrow = 2L)
  )
  on.exit(unlink(f))
  expect_identical(point_value(f, 260), 30)
  expect_identical(point_value(f, -30), 10)
  expect_message(
    expect_identical(point_value(f, 315), NA_real_),
    class = "ocean3d_spatial_report"
  )
})
