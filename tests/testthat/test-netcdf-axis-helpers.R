# Shared axis helpers used by point, track and buffer extraction and the core
# object (#113), tested directly (AGENTS.md, "Tests").

test_that("depth bounds are converted to metres positive down", {
  # Centres at 5 m and 15 m, stored as heights in centimetres.
  f <- make_grid_file(
    depth = c(-500, -1500), depth_units = "cm", positive = "up",
    depth_bounds = rbind(c(0, -1000), c(-1000, -2000)), time = NULL
  )
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  bounds <- .read_bounds(nc, "depth", 2L)
  expect_equal(
    .depth_bounds_m(nc, "depth", c(-500, -1500), bounds),
    rbind(c(0, 10), c(10, 20))
  )
})

test_that("midpoint edges keep the order of the centres", {
  regular <- .midpoint_edges(c(30, 20, 10))
  expect_equal(regular$lower, c(25, 15, 5))
  expect_equal(regular$upper, c(35, 25, 15))

  irregular <- .midpoint_edges(c(5, 15, 35))
  expect_equal(irregular$lower, c(0, 10, 25))
  expect_equal(irregular$upper, c(10, 25, 45))
})

test_that("bounds checks allow float32 rounding on global grids", {
  as_float32 <- function(x) {
    readBin(writeBin(x, raw(), size = 4L), "double", length(x), size = 4L)
  }
  # GLORYS-style 1/12 degree global longitudes stored as float32: the
  # inferred outer edges span slightly more than 360 degrees.
  lon <- as_float32(-180 + (0:4319) / 12)
  f <- make_grid_file(lon = lon, depth = NULL, time = NULL)
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_no_error(.read_interval_axis(nc, "lon", "lon", periodic = TRUE))
})

test_that("a repeated first meridian is only rejected around the globe", {
  f <- make_grid_file(lon = 0:360, depth = NULL, time = NULL)
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_no_error(.read_interval_axis(nc, "lon", "lon"))
  expect_error(.read_interval_axis(nc, "lon", "lon", periodic = TRUE),
               class = "ocean3d_error_bounds")
})

test_that("longitude cells that overlap modulo 360 are rejected", {
  # Increasing and non-overlapping as stored, but [340, 360] is [-20, 0].
  lon <- c(0, 120, 240, 350)
  f <- make_grid_file(
    lon = lon, lat = 0, depth = NULL, time = NULL,
    lon_bounds = matrix(c(-15, 15, 105, 135, 225, 255, 340, 360), 2L),
    lat_bounds = matrix(c(-1, 1), 2L)
  )
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_no_error(.read_interval_axis(nc, "lon", "lon"))
  err <- expect_error(
    .read_interval_axis(nc, "lon", "lon", periodic = TRUE),
    class = "ocean3d_error_bounds"
  )
  expect_match(conditionMessage(err), "overlap around the globe")
})

# Horizontal axes --------------------------------------------------------------

test_that("irregular descending axes preserve original indices", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_grid_file(lon = c(5, 2, 0), lat = 40, depth = NULL, time = NULL,
                 path = path)
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  axis <- .read_interval_axis(nc, "lon", "lon")
  expect_identical(axis$basis, "inferred_midpoints")
  # Each cell keeps its place in the file, with edges halfway between.
  expect_equal(axis$centres, c(5, 2, 0))
  expect_equal(axis$lower, c(3.5, 1, -1))
  expect_equal(axis$upper, c(6.5, 3.5, 1))
})

test_that("spatial axes check their values and wrap whole-globe cells", {
  bad_lat <- make_grid_file(lon = 0, lat = c(80, 95), depth = NULL,
                            time = NULL)
  missing_lon <- make_grid_file(lon = c(0, NaN), depth = NULL, time = NULL)
  # One longitude cell 360 degrees wide, stored from 0 to 360.
  globe <- make_grid_file(lon = 10, lat = 0, depth = NULL, time = NULL,
                          lon_bounds = matrix(c(0, 360), 2L),
                          lat_bounds = matrix(c(-1, 1), 2L))
  on.exit(unlink(c(bad_lat, missing_lon, globe)), add = TRUE)
  read <- function(f, axis, ...) {
    nc <- ncdf4::nc_open(f)
    on.exit(ncdf4::nc_close(nc))
    .read_interval_axis(nc, axis, axis, ...)
  }
  err <- expect_error(read(bad_lat, "lat"),
                      class = "ocean3d_error_spatial_coordinates")
  expect_match(conditionMessage(err), glue::glue(
    "`lat` in {basename(bad_lat)} must be within [-90, 90]"
  ), fixed = TRUE)
  err <- expect_error(read(missing_lon, "lon"),
                      class = "ocean3d_error_spatial_coordinates")
  expect_match(conditionMessage(err), "`lon` in .* must be finite numbers")
  expect_error(read(globe, "lon", infer_bounds = NA), "TRUE or FALSE")
  axis <- read(globe, "lon")
  expect_equal(c(axis$lower, axis$upper), c(-170, 190))
  expect_identical(axis$basis, "explicit_bounds")
})

test_that("latitude edges past the poles are clamped to +/-90", {
  inferred <- make_grid_file(lon = 0, lat = c(89, 90), depth = NULL,
                             time = NULL)
  explicit <- make_grid_file(lon = 0, lat = c(-89.5, 89.5), depth = NULL,
                             time = NULL, lat_bounds = cbind(c(-91, 0),
                                                             c(0, 91)))
  on.exit(unlink(c(inferred, explicit)), add = TRUE)
  for (f in c(inferred, explicit)) {
    nc <- ncdf4::nc_open(f)
    axis <- .read_interval_axis(nc, "lat", "lat")
    ncdf4::nc_close(nc)
    expect_identical(c(min(axis$lower), max(axis$upper)),
                     c(if (f == inferred) 88.5 else -90, 90))
  }
})

test_that("longitude and latitude stored out of order are rejected", {
  f <- make_grid_file(lon = c(2, 0, 5), lat = c(41, 40, 42), depth = NULL,
                      time = NULL)
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  for (axis in c("lon", "lat")) {
    err <- expect_error(.read_interval_axis(nc, axis, axis),
                        class = "ocean3d_error_axis_order")
    expect_match(conditionMessage(err), glue::glue(
      "`{axis}` in {basename(f)} must be strictly increasing"
    ))
  }
})

test_that("longitudes wrap around the antimeridian", {
  expect_equal(.normalise_longitude(c(-160, 200)), c(-160, -160))
})

# Depth axes (.prepare_depth_axis()) -------------------------------------------

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
  for (unit in c("centimeters", "centimetres", "Centimeter")) {
    expect_equal(.prepare_depth_axis(c(0, 500), unit)$depth_m, c(0, 5))
  }
  for (unit in c("kilometers", "kilometres", "kilometer")) {
    expect_equal(.prepare_depth_axis(c(0, 0.5), unit)$depth_m, c(0, 500))
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
  values <- "ocean3d_error_depth_values"
  err <- expect_error(.prepare_depth_axis(c(0, 50, 50), "m"), class = values)
  expect_match(conditionMessage(err), "duplicate depth(s) 50.", fixed = TRUE)
  expect_no_match(conditionMessage(err), "sign")
  for (args in list(list(c(0, -50), "m", "down"), list(c(0, 50), "m", "up"))) {
    err <- expect_error(do.call(.prepare_depth_axis, args), class = values)
    expect_match(conditionMessage(err),
                 "reference surface must be positive metres")
    expect_match(conditionMessage(err), "Check its sign")
  }
})

test_that("a stored depth axis must run one way", {
  f <- make_grid_file(depth = c(50, 0, 100), time = NULL)
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  err <- expect_error(.read_depth_axis(nc, "depth", c(50, 0, 100)),
                      class = "ocean3d_error_axis_order")
  expect_match(conditionMessage(err), glue::glue(
    "`depth` in {basename(f)} must be strictly increasing or decreasing"
  ))
  # Duplicates are caught first, with their own message.
  err <- expect_error(.read_depth_axis(nc, "depth", c(0, 50, 50)),
                      class = "ocean3d_error_depth_values")
  expect_match(conditionMessage(err), "duplicate depth")
  expect_equal(.read_depth_axis(nc, "depth", c(100, 50, 0))$depth_m,
               c(100, 50, 0))
})

# Bounds -----------------------------------------------------------------------

test_that("unit spellings are compared by meaning", {
  expect_true(.same_units("degree_east", "degrees_east"))
  expect_true(.same_units("degrees_E", "degrees"))
  expect_false(.same_units("degrees_east", "radians"))
  expect_true(.same_units("m", "metres"))
  expect_false(.same_units("m", "cm"))
  expect_false(.same_units("days since 2020-01-01", "hours since 2020-01-01"))
})

test_that("shared bounds edges that differ by rounding are accepted", {
  # 0.1 degree cells whose shared edges are stored once each.
  lon <- seq(100, 140, by = 0.1)
  edges <- seq(99.95, 140.05, by = 0.1)
  f <- make_grid_file(lon = lon, depth = NULL, time = NULL,
                      lon_bounds = rbind(head(edges, -1L), tail(edges, -1L)))
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  axis <- .read_interval_axis(nc, "lon", "lon")
  expect_identical(axis$basis, "explicit_bounds")
  expect_equal(c(axis$lower[201], axis$upper[201]), c(119.95, 120.05))
})

# Interval axes (.read_interval_axis()) ----------------------------------------

# Reads one depth or time axis of a file as intervals.
read_axis <- function(f, coordinate = "depth", ...) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc))
  .read_interval_axis(nc, coordinate, ...) # nolint: object_usage_linter.
}

# A file with one variable on a depth axis whose bounds variable is stored
# `bounds` as given, which make_grid_file() can't do for bounds stored
# coordinate first or with more than two vertices.
raw_bounds_file <- function(centres, bounds, coordinate_first = FALSE) {
  f <- tempfile(fileext = ".nc")
  depth <- ncdf4::ncdim_def("depth", "m", centres)
  nv <- ncdf4::ncdim_def("nv", "", seq_len(if (coordinate_first) ncol(bounds)
                                           else nrow(bounds)),
                         create_dimvar = FALSE)
  dims <- if (coordinate_first) list(depth, nv) else list(nv, depth)
  nc <- ncdf4::nc_create(f, list(
    ncdf4::ncvar_def("value", "1", list(depth), prec = "double"),
    ncdf4::ncvar_def("depth_bnds", "m", dims, prec = "double")
  ))
  ncdf4::ncvar_put(nc, "depth_bnds", bounds)
  ncdf4::ncatt_put(nc, "depth", "bounds", "depth_bnds")
  ncdf4::nc_close(nc)
  f
}

test_that("both declared bounds dimension orders are supported", {
  bounds <- matrix(c(0, 10, 10, 30, 30, 50), nrow = 2L)
  files <- c(make_grid_file(depth = c(5, 20, 40), time = NULL,
                            depth_bounds = bounds),
             raw_bounds_file(c(5, 20, 40), t(bounds), coordinate_first = TRUE))
  on.exit(unlink(files), add = TRUE)
  for (f in files) {
    support <- read_axis(f)
    expect_equal(support$lower, c(0, 10, 30))
    expect_equal(support$upper, c(10, 30, 50))
  }
})

test_that("missing bounds require explicit permission to infer", {
  f <- make_grid_file(depth = c(0, 10, 30), time = NULL)
  on.exit(unlink(f), add = TRUE)
  expect_error(read_axis(f), "infer_bounds = TRUE", fixed = TRUE)
  support <- read_axis(f, infer_bounds = TRUE)
  expect_equal(support$lower, c(0, 5, 20))
  expect_equal(support$upper, c(5, 20, 40))
  expect_identical(support$basis, "inferred_midpoints")
  expect_identical(support$clipped_at_surface, c(TRUE, FALSE, FALSE))
})

test_that("a singleton needs explicit bounds and keeps its extent", {
  bare <- make_grid_file(depth = 10, time = NULL)
  bounded <- make_grid_file(depth = 10, time = NULL,
                            depth_bounds = matrix(c(5, 15), nrow = 2L))
  on.exit(unlink(c(bare, bounded)), add = TRUE)
  expect_error(read_axis(bare, infer_bounds = TRUE), "single coordinate",
               fixed = TRUE)
  support <- read_axis(bounded)
  expect_equal(support$lower, 5)
  expect_equal(support$upper, 15)
})

test_that("time bounds decode to UTC and allow labelled intervals", {
  f <- make_grid_file(depth = NULL, time = c(0, 24),
                      time_units = "hours since 2020-01-01 00:00:00",
                      calendar = "gregorian",
                      time_bounds = matrix(c(0, 24, 24, 48), nrow = 2L))
  on.exit(unlink(f), add = TRUE)
  support <- read_axis(f, "time", axis = "time")
  expect_equal(support$lower,
               as.numeric(as.POSIXct(c("2020-01-01", "2020-01-02"),
                                     tz = "UTC")))
  expect_equal(support$upper,
               as.numeric(as.POSIXct(c("2020-01-02", "2020-01-03"),
                                     tz = "UTC")))
  expect_identical(support$calendar, "gregorian")
  expect_identical(support$basis, "explicit_bounds")
})

test_that("irregular time inference uses neighbouring midpoints", {
  f <- make_grid_file(depth = NULL, time = c(0, 2, 8),
                      time_units = "hours since 2020-01-01 00:00:00")
  on.exit(unlink(f), add = TRUE)
  support <- read_axis(f, "time", axis = "time", infer_bounds = TRUE)
  origin <- as.numeric(as.POSIXct("2020-01-01", tz = "UTC"))
  expect_equal(support$lower, origin + c(-3600, 3600, 18000))
  expect_equal(support$upper, origin + c(3600, 18000, 39600))
  expect_identical(support$basis, "inferred_midpoints")
})

test_that("invalid declared bounds never fall back to inference", {
  cases <- list(
    list(bounds = matrix(c(0, 10, 5, 30), nrow = 2L), pattern = "overlap"),
    list(
      bounds = matrix(c(5, 5, 10, 30), nrow = 2L),
      pattern = "positive widths"
    ),
    list(bounds = matrix(c(NA, 10, 10, 30), nrow = 2L), pattern = "finite")
  )
  for (case in cases) {
    f <- make_grid_file(depth = c(5, 20), time = NULL,
                        depth_bounds = case$bounds)
    on.exit(unlink(f), add = TRUE)
    expect_error(read_axis(f, infer_bounds = TRUE), case$pattern,
                 fixed = TRUE)
  }
  three <- raw_bounds_file(c(5, 20), matrix(0:5, nrow = 3L))
  missing <- make_grid_file(depth = c(5, 20), time = NULL,
                            depth_atts = list(bounds = "missing_bounds"))
  on.exit(unlink(c(three, missing)), add = TRUE)
  expect_error(read_axis(three, infer_bounds = TRUE), "two finite endpoints",
               fixed = TRUE)
  expect_error(read_axis(missing, infer_bounds = TRUE),
               "missing bounds variable", fixed = TRUE)
})

test_that("incompatible units and directions are reported explicitly", {
  bounds <- matrix(c(0, 10, 10, 30), nrow = 2L)
  f <- make_grid_file(depth = c(5, 20), time = NULL, depth_bounds = bounds,
                      positive = "down")
  on.exit(unlink(f), add = TRUE)
  set_atts(f, "depth_bnds", units = "cm")
  expect_error(read_axis(f), "different units", fixed = TRUE)
  set_atts(f, "depth_bnds", units = "m", positive = "up")
  expect_error(read_axis(f), "conflicting directions", fixed = TRUE)
})

test_that("kilometres convert to metres and transformed axes are rejected", {
  f <- make_grid_file(depth = c(0.005, 0.02), depth_units = "km", time = NULL,
                      depth_bounds = matrix(c(0, 0.01, 0.01, 0.03), nrow = 2L))
  on.exit(unlink(f), add = TRUE)
  support <- read_axis(f)
  expect_equal(support$centres, c(5, 20))
  expect_equal(support$lower, c(0, 10))
  expect_equal(support$upper, c(10, 30))
  for (atts in list(
    list(standard_name = "sea_water_pressure"),
    list(formula_terms = "sigma: sigma eta: eta depth: depth")
  )) {
    g <- make_grid_file(depth = c(5, 20), time = NULL, depth_atts = atts)
    on.exit(unlink(g), add = TRUE)
    expect_error(read_axis(g, infer_bounds = TRUE), "physical transformation",
                 fixed = TRUE)
  }
})
