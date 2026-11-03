# The core object (#113): as_proxy() on synthetic netCDF files built by
# make_grid_file() (helper-netcdf-grid.R).

proxy_of <- function(...) {
  f <- make_grid_file(...)
  on.exit(unlink(f))
  as_proxy(f, var = "temp")
}

proxy_error <- function(...) {
  f <- make_grid_file(...)
  on.exit(unlink(f))
  list(
    file = basename(f),
    error = tryCatch(as_proxy(f), error = identity)
  )
}

test_that("an ordinary lon, lat, depth and time grid is described", {
  f <- make_grid_file()
  on.exit(unlink(f), add = TRUE)
  x <- as_proxy(f)

  expect_s4_class(x, "SpatProxy")
  expect_identical(x@source, list(type = "netcdf", path = normalizePath(f)))
  expect_identical(x@variable$name, "temp")
  expect_identical(x@variable$dims, c("lon", "lat", "depth", "time"))
  expect_identical(x@variable$dim_lengths,
                   c(lon = 3L, lat = 2L, depth = 2L, time = 2L))
  expect_identical(x@crs, "EPSG:4326")

  expect_equal(x@axes$lon$centres, c(0.5, 1.5, 2.5))
  expect_equal(x@axes$lon$lower, c(0, 1, 2))
  expect_identical(x@axes$lon$edges, "inferred")
  expect_equal(x@axes$lat$upper, c(11, 12))
  expect_equal(x@axes$depth$centres, c(5, 15))
  # The surface edge is cut at 0 m.
  expect_equal(x@axes$depth$lower, c(0, 10))
  expect_equal(x@axes$depth$upper, c(10, 20))
  expect_identical(x@axes$depth$storage, "shallow_first")

  time <- x@axes$time$cftime
  expect_s3_class(time, "CFTime")
  expect_identical(time$offsets, c(0, 1))
  expect_identical(time$cal$name, "standard")
})

test_that("explicit bounds are used before inferred edges", {
  x <- proxy_of(
    lon_bounds = rbind(c(0, 1, 2.2), c(1, 2.2, 3)),
    lat_bounds = rbind(c(10, 11), c(11, 12)),
    depth_bounds = rbind(c(0, 8), c(8, 30)),
    time_bounds = rbind(c(-0.5, 0.5), c(0.5, 1.5))
  )

  expect_identical(x@axes$lon$edges, "bounds")
  expect_equal(x@axes$lon$upper, c(1, 2.2, 3))
  expect_identical(x@axes$depth$edges, "bounds")
  expect_equal(x@axes$depth$lower, c(0, 8))
  expect_equal(x@axes$depth$upper, c(8, 30))
  expect_identical(x@axes$time$edges, "bounds")
  expect_equal(x@axes$time$cftime$get_bounds(),
               rbind(c(-0.5, 0.5), c(0.5, 1.5)))

  # Without bounds, the time axis has none: they aren't inferred.
  plain <- proxy_of()
  expect_identical(plain@axes$time$edges, "none")
  expect_null(plain@axes$time$cftime$get_bounds())
})

test_that("the longitude convention is recorded", {
  expect_identical(proxy_of(lon = c(90, 270))@axes$lon$convention, "0_360")
  expect_identical(proxy_of(lon = c(-90, 90))@axes$lon$convention,
                   "-180_180")
})

test_that("depth is metres positive down in its storage order", {
  reversed <- proxy_of(depth = c(15, 5))
  expect_equal(reversed@axes$depth$centres, c(15, 5))
  expect_equal(reversed@axes$depth$lower, c(10, 0))
  expect_identical(reversed@axes$depth$storage, "deep_first")

  centimetres <- proxy_of(depth = c(500, 1500), depth_units = "cm")
  expect_equal(centimetres@axes$depth$centres, c(5, 15))

  heights <- proxy_of(depth = c(-5, -15), positive = "up",
                      depth_bounds = rbind(c(0, -10), c(-10, -20)))
  expect_equal(heights@axes$depth$centres, c(5, 15))
  expect_equal(heights@axes$depth$upper, c(10, 20))

  f <- make_grid_file(depth = c(5, 15, 10))
  on.exit(unlink(f), add = TRUE)
  expect_error(as_proxy(f), class = "ocean3d_error_depth_order")
})

test_that("a grid without depth or time has no such axis", {
  x <- proxy_of(depth = NULL, time = NULL)
  expect_null(x@axes$depth)
  expect_null(x@axes$time)
  expect_identical(x@variable$dims, c("lon", "lat"))

  surface <- proxy_of(depth = 0.5, time = NULL)
  expect_identical(surface@axes$depth$edges, "none")
})

test_that("packing and chunk sizes are recorded", {
  x <- proxy_of(prec = "short", chunksizes = c(3, 1, 1, 1),
                scale = 0.01, offset = 20)
  packing <- x@variable$packing
  expect_identical(packing$storage, "short")
  expect_equal(packing$scale_factor, 0.01)
  expect_equal(packing$add_offset, 20)
  expect_equal(packing$fill_value, -32767)
  expect_null(packing$missing_value)
  expect_identical(x@variable$chunksizes,
                   c(lon = 3L, lat = 1L, depth = 1L, time = 1L))

  expect_null(proxy_of()@variable$chunksizes)
})

test_that("a time axis that isn't increasing stops and names the file", {
  for (time in list(c(2, 1, 0), c(0, 2, 1))) {
    out <- proxy_error(time = time)
    expect_s3_class(out$error, "ocean3d_error_time_order")
    expect_match(conditionMessage(out$error), out$file, fixed = TRUE)
    expect_match(conditionMessage(out$error), "sorttimestamp", fixed = TRUE)
  }
  backwards <- proxy_error(time = c(2, 1, 0))
  expect_match(conditionMessage(backwards$error), "runs backwards")

  # Sorting can't fix repeated time steps, so they have their own error.
  repeated <- proxy_error(time = c(0, 1, 1))
  expect_s3_class(repeated$error, "ocean3d_error_time_repeated")
  expect_match(conditionMessage(repeated$error), repeated$file, fixed = TRUE)
  expect_match(conditionMessage(repeated$error), "repeats a time step")
})

test_that("unsupported calendars and climatologies stop", {
  calendar <- proxy_error(calendar = "noleap")
  expect_s3_class(calendar$error, "ocean3d_error_calendar")
  expect_match(conditionMessage(calendar$error), "noleap", fixed = TRUE)

  # WOA-style monthly climatology, which CFtime can't read.
  woa <- proxy_error(time = c(0.5, 1.5),
                     time_units = "months since 0000-01-01 00:00:00")
  expect_s3_class(woa$error, "ocean3d_error_time_climatology")
  expect_match(conditionMessage(woa$error), woa$file, fixed = TRUE)

  flagged <- proxy_error(time_atts = list(climatology = "climatology_bnds"))
  expect_s3_class(flagged$error, "ocean3d_error_time_climatology")
})

test_that("the CRS comes from a geographic grid mapping", {
  add_mapping <- function(f, wkt) {
    nc <- ncdf4::nc_open(f, write = TRUE)
    nc <- ncdf4::ncvar_add(nc, ncdf4::ncvar_def("crs", "", list(),
                                                prec = "integer"))
    ncdf4::ncatt_put(nc, "crs", "crs_wkt", wkt)
    ncdf4::ncatt_put(nc, "temp", "grid_mapping", "crs")
    ncdf4::nc_close(nc)
  }
  f <- make_grid_file()
  g <- make_grid_file()
  on.exit(unlink(c(f, g)), add = TRUE)

  etrs89 <- terra::crs("EPSG:4258")
  add_mapping(f, etrs89)
  expect_identical(as_proxy(f, var = "temp")@crs, etrs89)

  add_mapping(g, terra::crs("EPSG:3857"))
  expect_error(as_proxy(g, var = "temp"), class = "ocean3d_error_grid_mapping")
})

test_that("an open connection is read and left open", {
  f <- make_grid_file()
  nc <- ncdf4::nc_open(f)
  on.exit({
    ncdf4::nc_close(nc)
    unlink(f)
  }, add = TRUE)

  x <- as_proxy(nc)
  expect_identical(x@source$path, normalizePath(f))
  expect_identical(nc$var$temp$name, "temp")
  expect_no_error(ncdf4::ncvar_get(nc, "temp"))

  expect_error(as_proxy(c(f, f)), class = "ocean3d_error_proxy_input")
})

test_that("a SpatProxy survives saveRDS() and checks its axes", {
  x <- proxy_of(time_bounds = rbind(c(-0.5, 0.5), c(0.5, 1.5)))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path)
  y <- readRDS(path)
  expect_equal(y@axes$time$cftime$get_bounds(),
               x@axes$time$cftime$get_bounds())
  expect_true(methods::validObject(y))

  y@axes$lon$lower <- y@axes$lon$lower[-1L]
  expect_error(methods::validObject(y), "lon axis needs one centre")
})

test_that("the print shows each axis", {
  x <- proxy_of(prec = "short", chunksizes = c(3, 1, 1, 1),
                time_bounds = rbind(c(-0.5, 0.5), c(0.5, 1.5)))
  printed <- capture.output(print(x))

  # The source line holds the temporary file's random name.
  expect_match(printed[2L], "^source     : file.*\\.nc$")
  expect_identical(printed[-2L], c(
    "class      : SpatProxy (experimental)",
    "variable   : temp (degC), stored as short, packed, chunks 3 x 1 x 1 x 1",
    "dimensions : 3 x 2 x 2 x 2 (lon x lat x depth x time)",
    "lon        : 0.5 to 2.5, 3 cells, -180 to 180, edges inferred",
    "lat        : 10.5 to 11.5, 2 cells, edges inferred",
    "depth      : 5 to 15 m, 2 levels, shallow first, edges inferred",
    paste("time       : 2020-01-01 to 2020-01-02, 2 steps,",
          "standard calendar, with bounds"),
    "crs        : EPSG:4326"
  ))

  bare <- capture.output(print(proxy_of(depth = NULL, time = NULL)))
  expect_identical(bare[3:4], c(
    "variable   : temp (degC), stored as double",
    "dimensions : 3 x 2 (lon x lat)"
  ))
  expect_identical(bare[7:8], c("depth      : none", "time       : none"))
})
