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
  err <- expect_error(as_proxy(f), class = "ocean3d_error_axis_order")
  expect_match(conditionMessage(err), "`depth` must be strictly increasing")
})

test_that("depth bounds that don't contain their centres are rejected", {
  f <- make_grid_file(depth = c(5, 15),
                      depth_bounds = rbind(c(20, 30), c(30, 40)))
  on.exit(unlink(f), add = TRUE)
  expect_error(as_proxy(f), "must contain their centres")
})

test_that("a grid without depth or time has no such axis", {
  x <- proxy_of(depth = NULL, time = NULL)
  expect_null(x@axes$depth)
  expect_null(x@axes$time)
  expect_identical(x@variable$dims, c("lon", "lat"))

  surface <- proxy_of(depth = 0.5, time = NULL)
  expect_identical(surface@axes$depth$edges, "none")
})

test_that("chunk sizes are recorded", {
  x <- proxy_of(prec = "short", chunksizes = c(3, 1, 1, 1),
                scale = 0.01, offset = 20)
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
  expect_s3_class(repeated$error, "ocean3d_error_repeated_time")
  expect_match(conditionMessage(repeated$error), repeated$file, fixed = TRUE)
  expect_match(conditionMessage(repeated$error),
               "Time step 1 (days since 2020-01-01) appears 2 times",
               fixed = TRUE)
})

test_that("bad time units, offsets or bounds calendars get their own error", {
  no_units <- make_grid_file(time = c(0, 1, 1), time_units = "")
  non_finite <- make_grid_file(time = c(3, 2, NaN, 1))
  calendars <- make_grid_file(time = c(0.5, 1.5),
                              time_bounds = rbind(c(0, 1), c(1, 2)))
  on.exit(unlink(c(no_units, non_finite, calendars)), add = TRUE)
  expect_error(as_proxy(no_units), "has no usable `units` attribute")
  expect_error(as_proxy(non_finite), "must contain finite numeric offsets")
  set_atts(calendars, "time_bnds", calendar = "noleap")
  expect_error(as_proxy(calendars), "use different calendars")
})

test_that("contiguous netCDF-4 storage has no chunk sizes", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0.5, 1.5))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(10.5, 11.5))
  nc <- ncdf4::nc_create(f, ncdf4::ncvar_def("temp", "degC", list(lon, lat)),
                         force_v4 = TRUE)
  ncdf4::nc_close(nc)
  expect_null(as_proxy(f)@variable$chunksizes)
})

test_that("a connection whose relative path can't be found stops", {
  dir <- tempfile("proxy-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  file.copy(make_grid_file(), file.path(dir, "grid.nc"))
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  nc <- ncdf4::nc_open("grid.nc")
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_identical(as_proxy(nc)@source$path,
                   normalizePath(file.path(dir, "grid.nc")))
  # After the working directory changes, the relative path is ambiguous.
  setwd(old)
  expect_error(as_proxy(nc), class = "ocean3d_error_proxy_source")
})

test_that("longitude and latitude need coordinate variables", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  lon <- ncdf4::ncdim_def("lon", "", 1:3, create_dimvar = FALSE)
  lat <- ncdf4::ncdim_def("lat", "", 1:2, create_dimvar = FALSE)
  nc <- ncdf4::nc_create(f, ncdf4::ncvar_def("temp", "degC", list(lon, lat)))
  ncdf4::nc_close(nc)
  err <- expect_error(as_proxy(f), class = "ocean3d_error_coordinate_variable")
  expect_match(conditionMessage(err), "`lon` has no coordinate variable")
})

test_that("rotated-pole grid mappings are rejected", {
  add_mapping <- function(f, name = NULL, wkt = NULL) {
    nc <- ncdf4::nc_open(f, write = TRUE)
    nc <- ncdf4::ncvar_add(nc, ncdf4::ncvar_def("crs", "", list(),
                                                prec = "integer"))
    if (!is.null(name)) ncdf4::ncatt_put(nc, "crs", "grid_mapping_name", name)
    if (!is.null(wkt)) ncdf4::ncatt_put(nc, "crs", "crs_wkt", wkt)
    ncdf4::ncatt_put(nc, "temp", "grid_mapping", "crs")
    ncdf4::nc_close(nc)
  }
  named <- make_grid_file()
  derived <- make_grid_file()
  on.exit(unlink(c(named, derived)), add = TRUE)
  # A rotated mapping name, even with a geographic WKT.
  add_mapping(named, "rotated_latitude_longitude", terra::crs("EPSG:4326"))
  expect_error(as_proxy(named), class = "ocean3d_error_grid_mapping")
  # A rotated-pole WKT, which terra::is.lonlat() calls lon/lat.
  rotated <- terra::crs(paste("+proj=ob_tran +o_proj=longlat +o_lon_p=0",
                              "+o_lat_p=37.5 +lon_0=357.5 +datum=WGS84"))
  add_mapping(derived, wkt = rotated)
  expect_error(as_proxy(derived), class = "ocean3d_error_grid_mapping")
})

test_that("overlapping time bounds stop", {
  f <- make_grid_file(time = c(0.5, 1.5),
                      time_bounds = rbind(c(0, 0.8), c(1.2, 2)))
  on.exit(unlink(f), add = TRUE)
  err <- expect_error(as_proxy(f), class = "ocean3d_error_bounds")
  expect_match(conditionMessage(err), "Bounds for `time` in .* overlap\\.")
})

test_that("unsupported calendars and climatologies stop", {
  calendar <- proxy_error(calendar = "noleap")
  expect_s3_class(calendar$error, "ocean3d_error_calendar")

  # WOA-style monthly climatology, which CFtime can't read. The calendar and
  # climatology checks are tested directly in test-netcdf-time.R.
  woa <- proxy_error(time = c(0.5, 1.5),
                     time_units = "months since 0000-01-01 00:00:00")
  expect_s3_class(woa$error, "ocean3d_error_time_climatology")
})

test_that("the CRS comes from a geographic grid mapping", {
  add_mapping <- function(f, wkt, grid_mapping = "crs") {
    nc <- ncdf4::nc_open(f, write = TRUE)
    nc <- ncdf4::ncvar_add(nc, ncdf4::ncvar_def("crs", "", list(),
                                                prec = "integer"))
    ncdf4::ncatt_put(nc, "crs", "crs_wkt", wkt)
    ncdf4::ncatt_put(nc, "temp", "grid_mapping", grid_mapping)
    ncdf4::nc_close(nc)
  }
  f <- make_grid_file()
  g <- make_grid_file()
  h <- make_grid_file(depth = NULL, time = NULL)
  on.exit(unlink(c(f, g, h)), add = TRUE)

  # The grid-mapping variable isn't data, so `var` is detected.
  etrs89 <- terra::crs("EPSG:4258")
  add_mapping(f, etrs89)
  expect_identical(as_proxy(f)@crs, etrs89)

  add_mapping(g, terra::crs("EPSG:3857"))
  expect_error(as_proxy(g), class = "ocean3d_error_grid_mapping")

  # Point extraction detects the variable the same way, also with CF's
  # extended form of the attribute.
  add_mapping(h, etrs89, grid_mapping = "crs: lat lon")
  expect_equal(
    extract_to_point(nc = h, lon = 1.5, lat = 10.5, method = "2d"), 2
  )
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

test_that("unlimited, empty and unsaved time axes", {
  unlimited_file <- function(times) {
    f <- tempfile(fileext = ".nc")
    dims <- list(
      ncdf4::ncdim_def("lon", "degrees_east", c(0.5, 1.5, 2.5)),
      ncdf4::ncdim_def("lat", "degrees_north", c(10.5, 11.5)),
      ncdf4::ncdim_def("time", "days since 2020-01-01", times, unlim = TRUE)
    )
    nc <- ncdf4::nc_create(f, ncdf4::ncvar_def("temp", "degC", dims, -9999))
    if (length(times)) ncdf4::ncvar_put(nc, "temp", seq_len(6 * length(times)))
    ncdf4::nc_close(nc)
    f
  }
  f <- unlimited_file(c(0, 1))
  empty <- unlimited_file(numeric(0))
  on.exit(unlink(c(f, empty)), add = TRUE)

  # A classic file with an unlimited dimension has no chunk sizes.
  x <- as_proxy(f)
  expect_null(x@variable$chunksizes)
  expect_identical(dim(.read_proxy_window(x)), c(3L, 2L, 1L, 2L))

  err <- expect_error(as_proxy(empty), class = "ocean3d_error_empty_dimension")
  expect_match(conditionMessage(err), "no values along time")

  # Values written through a connection but not yet saved are read.
  nc <- ncdf4::nc_open(f, write = TRUE)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  ncdf4::ncvar_put(nc, "temp", rep(100, 6), start = c(1, 1, 1),
                   count = c(3, 2, 1))
  y <- as_proxy(nc)
  expect_identical(as.vector(.read_proxy_window(y, list(time = c(1, 1)))),
                   rep(100, 6))
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

  bad <- y
  bad@axes$lon$lower <- bad@axes$lon$lower[-1L]
  expect_error(methods::validObject(bad), "lon axis needs a dimension")

  # Problems give a message, not an R error from inside the check.
  bad <- y
  bad@axes$lat$dim <- "nope"
  expect_error(methods::validObject(bad), "lat axis needs a dimension")
  bad <- y
  bad@variable$dim_lengths[["time"]] <- 3L
  expect_error(methods::validObject(bad), "one offset per time step")
})

test_that("the print shows each axis", {
  x <- proxy_of(prec = "short", chunksizes = c(3, 1, 1, 1),
                time_bounds = rbind(c(-0.5, 0.5), c(0.5, 1.5)))
  printed <- capture.output(print(x))

  # The source line holds the temporary file's random name.
  expect_match(printed[2L], "^source     : file.*\\.nc$")
  expect_identical(printed[-2L], c(
    "class      : SpatProxy (experimental)",
    "variable   : temp (degC), chunks 3 x 1 x 1 x 1",
    "dimensions : 3 x 2 x 2 x 2 (lon x lat x depth x time)",
    "lon        : 0.5 to 2.5, 3 cells, edges inferred",
    "lat        : 10.5 to 11.5, 2 cells, edges inferred",
    "depth      : 5 to 15 m, 2 levels, shallow first, edges inferred",
    paste("time       : 2020-01-01 to 2020-01-02, 2 steps,",
          "standard calendar, with bounds"),
    "crs        : EPSG:4326"
  ))

  bare <- capture.output(print(proxy_of(depth = NULL, time = NULL)))
  expect_identical(bare[3:4], c(
    "variable   : temp (degC)",
    "dimensions : 3 x 2 (lon x lat)"
  ))
  expect_identical(bare[7:8], c("depth      : none", "time       : none"))
})

# The index-window reader. The file must outlive the proxy, so each test
# removes it at the end.
proxy_and_file <- function(...) {
  f <- make_grid_file(...)
  list(proxy = as_proxy(f, var = "temp"), file = f)
}

stored_values <- function(f) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_get(nc, "temp", collapse_degen = FALSE)
}

test_that("a window reads the same values as the whole variable", {
  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)
  full <- stored_values(g$file)

  expect_identical(.read_proxy_window(g$proxy), full)
  window <- .read_proxy_window(
    g$proxy, list(lon = c(2, 3), depth = c(2, 2), time = c(1, 1))
  )
  expect_identical(window, full[2:3, , 2, 1, drop = FALSE])
})

test_that("values come back in lon, lat, depth, time order", {
  # A storage order whose reordering isn't its own inverse.
  g <- proxy_and_file(dim_order = c("depth", "lon", "time", "lat"))
  on.exit(unlink(g$file), add = TRUE)
  fixed <- aperm(stored_values(g$file), c(2, 4, 1, 3))

  expect_identical(dim(.read_proxy_window(g$proxy)), c(3L, 2L, 2L, 2L))
  expect_identical(.read_proxy_window(g$proxy), fixed)
  expect_identical(
    .read_proxy_window(g$proxy, list(lat = c(2, 2), time = c(2, 2))),
    fixed[, 2, , 2, drop = FALSE]
  )

  # A deep-first depth axis is read in storage order: index 1 is deepest.
  deep <- proxy_and_file(depth = c(15, 5))
  on.exit(unlink(deep$file), add = TRUE)
  expect_identical(
    .read_proxy_window(deep$proxy, list(depth = c(1, 1))),
    stored_values(deep$file)[, , 1, , drop = FALSE]
  )
})

test_that("missing axes and other length-1 dimensions have length 1", {
  flat <- proxy_and_file(depth = NULL, time = NULL)
  on.exit(unlink(flat$file), add = TRUE)
  expect_identical(dim(.read_proxy_window(flat$proxy)), c(3L, 2L, 1L, 1L))
  no_depth <- proxy_and_file(depth = NULL)
  on.exit(unlink(no_depth$file), add = TRUE)
  expect_identical(dim(.read_proxy_window(no_depth$proxy)), c(3L, 2L, 1L, 2L))

  # A length-1 `band` dimension between lon and lat.
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)
  dims <- list(
    ncdf4::ncdim_def("lon", "degrees_east", c(0.5, 1.5, 2.5)),
    ncdf4::ncdim_def("band", "", 1L, create_dimvar = FALSE),
    ncdf4::ncdim_def("lat", "degrees_north", c(10.5, 11.5))
  )
  nc <- ncdf4::nc_create(f, ncdf4::ncvar_def("temp", "degC", dims, -9999))
  ncdf4::ncvar_put(nc, "temp", array(1:6, c(3, 1, 2)))
  ncdf4::nc_close(nc)
  expect_identical(
    .read_proxy_window(as_proxy(f), list(lat = c(2, 2))),
    array(c(4, 5, 6), c(3, 1, 1, 1))
  )
})

test_that("packed values are unpacked and fill values are NA", {
  g <- proxy_and_file(prec = "short", chunksizes = c(1, 1, 1, 1))
  on.exit(unlink(g$file), add = TRUE)
  nc <- ncdf4::nc_open(g$file, write = TRUE)
  ncdf4::ncvar_put(nc, "temp", NA, start = c(2, 1, 1, 1),
                   count = c(1, 1, 1, 1))
  ncdf4::nc_close(nc)

  values <- .read_proxy_window(g$proxy, list(depth = c(1, 1), time = c(1, 1)))
  expect_equal(as.vector(values), c(1, NA, 3, 4, 5, 6) * 0.1)
})

test_that("bad read windows and changed files stop", {
  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)
  bad <- list(
    list(lon = c(0, 2)), list(lon = c(2, 4)), list(lon = c(3, 2)),
    list(lon = 2), list(lon = c(1.5, 2)), list(lon = c(NA, 2)),
    list(lon = c("1", "2"))
  )
  for (window in bad) {
    err <- expect_error(.read_proxy_window(g$proxy, window),
                        class = "ocean3d_error_read_window")
    expect_match(conditionMessage(err), "`lon` must be .* from 1 to 3")
  }
  # Not a list of named axes, or an axis named twice.
  unnamed <- list(
    list(x = c(1, 1)), list(c(1, 1)), c(lon = 2),
    list(lon = c(1, 1), lon = c(2, 3))
  )
  for (window in unnamed) {
    err <- expect_error(.read_proxy_window(g$proxy, window),
                        class = "ocean3d_error_read_window")
    expect_match(conditionMessage(err), "list naming each axis at most once")
  }

  # An axis the variable doesn't have has length 1.
  flat <- proxy_and_file(depth = NULL, time = NULL)
  on.exit(unlink(flat$file), add = TRUE)
  expect_identical(dim(.read_proxy_window(flat$proxy, list(depth = c(1, 1)))),
                   c(3L, 2L, 1L, 1L))
  err <- expect_error(.read_proxy_window(flat$proxy, list(depth = c(2, 5))),
                      class = "ocean3d_error_read_window")
  expect_match(conditionMessage(err), "`depth` must be .* from 1 to 1")
})

test_that("a file that no longer matches the proxy stops", {
  changed <- function(g, msg) {
    err <- expect_error(.read_proxy_window(g$proxy),
                        class = "ocean3d_error_proxy_changed")
    expect_match(conditionMessage(err), msg)
    expect_match(conditionMessage(err), "Build the SpatProxy again")
  }
  rewrite <- function(g, ...) {
    other <- make_grid_file(...)
    file.copy(other, g$file, overwrite = TRUE)
    unlink(other)
  }
  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)

  # Another grid, or the same sizes with depth and time swapped.
  rewrite(g, lon = c(0.5, 1.5))
  changed(g, "no longer has the dimensions")
  rewrite(g)
  expect_no_error(.read_proxy_window(g$proxy))
  rewrite(g, dim_order = c("lon", "lat", "time", "depth"))
  changed(g, "no longer has the dimensions")

  # The variable or the file is gone.
  renamed <- g
  renamed$proxy@variable$name <- "salt"
  changed(renamed, "no longer has the variable `salt`")
  unlink(g$file)
  changed(g, "Can't find")
})

test_that("a proxy of a URL source is read", {
  store <- file.path(tempfile(), "grid.zarr")
  dir.create(dirname(store))
  on.exit(unlink(dirname(store), recursive = TRUE), add = TRUE)
  url <- paste0("file://", store, "#mode=nczarr,file")
  dims <- list(
    ncdf4::ncdim_def("lon", "degrees_east", c(0.5, 1.5, 2.5)),
    ncdf4::ncdim_def("lat", "degrees_north", c(10.5, 11.5))
  )
  made <- tryCatch({
    nc <- ncdf4::nc_create(url, ncdf4::ncvar_def("temp", "degC", dims, -9999),
                           force_v4 = TRUE)
    ncdf4::ncvar_put(nc, "temp", array(1:6, c(3, 2)))
    ncdf4::nc_close(nc)
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(made, "the netCDF library can't write NCZarr stores")

  nc <- ncdf4::nc_open(url)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  x <- as_proxy(nc)
  expect_identical(x@source$path, url)
  expect_identical(.read_proxy_window(x, list(lat = c(2, 2))),
                   array(c(4, 5, 6), c(3, 1, 1, 1)))
})
