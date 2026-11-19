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
  expect_true(regular$regular)

  irregular <- .midpoint_edges(c(5, 15, 35))
  expect_equal(irregular$lower, c(0, 10, 25))
  expect_equal(irregular$upper, c(10, 25, 45))
  expect_false(irregular$regular)
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
  expect_no_error(.prepare_spatial_axis(nc, "lon", lon, TRUE, periodic = TRUE))
})

test_that("a repeated first meridian is only rejected around the globe", {
  f <- make_grid_file(lon = 0:360, depth = NULL, time = NULL)
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_no_error(.prepare_spatial_axis(nc, "lon", 0:360, TRUE))
  expect_error(.prepare_spatial_axis(nc, "lon", 0:360, TRUE, periodic = TRUE),
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
  expect_no_error(.prepare_spatial_axis(nc, "lon", lon, TRUE))
  err <- expect_error(
    .prepare_spatial_axis(nc, "lon", lon, TRUE, periodic = TRUE),
    class = "ocean3d_error_bounds"
  )
  expect_match(conditionMessage(err), "overlap around the globe")
})

# Horizontal axes and cell matching --------------------------------------------

test_that("spatial indices respect edges and longitude conventions", {
  axis <- list(
    centres = c(0, 1), lower = c(-0.5, 0.5), upper = c(0.5, 1.5),
    longitude = FALSE, basis = "inferred_regular"
  )
  expect_identical(.containing_axis_indices(axis, -0.5), 1L)
  expect_identical(.containing_axis_indices(axis, 1.5), 2L)
  expect_length(.containing_axis_indices(axis, 1.6), 0L)
  axis <- list(
    centres = c(179, 181), lower = c(178, 180), upper = c(180, 182),
    longitude = TRUE, basis = "inferred_regular"
  )
  expect_identical(.containing_axis_indices(axis, -179), 2L)
  expect_identical(.containing_axis_indices(axis, 181), 2L)
  expect_length(.containing_axis_indices(axis, 0), 0L)
})

test_that("singleton axes do not invent spatial coverage", {
  axis <- list(
    centres = 40, lower = 40, upper = 40,
    longitude = FALSE, basis = "singleton_coordinate"
  )
  expect_identical(.containing_axis_indices(axis, 40), 1L)
  expect_length(.containing_axis_indices(axis, 41), 0L)
})

test_that("single-precision station coordinates match their own station", {
  station <- as.numeric(c(-122.40000152587891, 36.79999923706055))
  lon <- list(centres = station[1], lower = station[1], upper = station[1],
              longitude = TRUE, basis = "singleton_coordinate")
  lat <- list(centres = station[2], lower = station[2], upper = station[2],
              longitude = FALSE, basis = "singleton_coordinate")
  expect_identical(.containing_axis_indices(lon, -122.4), 1L)
  expect_identical(.containing_axis_indices(lat, 36.8), 1L)
  expect_length(.containing_axis_indices(lon, -122.401), 0L)
})

test_that("irregular descending axes preserve original indices", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path), add = TRUE)
  make_grid_file(lon = c(5, 2, 0), lat = 40, depth = NULL, time = NULL,
                 path = path)
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  axis <- .prepare_spatial_axis(nc, "lon", c(5, 2, 0), longitude = TRUE)
  expect_identical(axis$basis, "inferred_irregular")
  expect_identical(.containing_axis_indices(axis, 5), 1L)
  expect_identical(.containing_axis_indices(axis, 0), 3L)
  expect_length(.containing_axis_indices(axis, 20), 0L)
})

test_that("longitude and latitude stored out of order are rejected", {
  f <- make_grid_file(lon = c(2, 0, 5), lat = c(41, 40, 42), depth = NULL,
                      time = NULL)
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  for (case in list(list("lon", c(2, 0, 5), TRUE),
                    list("lat", c(41, 40, 42), FALSE))) {
    err <- expect_error(
      .prepare_spatial_axis(nc, case[[1]], case[[2]], longitude = case[[3]]),
      class = "ocean3d_error_axis_order"
    )
    expect_match(conditionMessage(err), glue::glue(
      "`{case[[1]]}` in {basename(f)} must be strictly increasing"
    ))
  }
})

test_that("geographic proximity is not independent degree subtraction", {
  # At high latitude the geographically closest centre is not necessarily
  # the one with the smallest latitude difference.
  result <- .nearest_horizontal_cell(60, c(80, 85), 1L, 1:2, 0, 80)
  expect_identical(result$lat_idx, 2L)
  expect_identical(result$lon_idx, 1L)
  expect_gt(result$distance_km, 0)
  whole <- .nearest_horizontal_cells(60, c(80, 85), 0, 80)[[1L]]
  expect_equal(result, whole)
})

test_that("nearest cells and longitudes wrap around the antimeridian", {
  selected <- .nearest_horizontal_cell(
    c(179, -170), 40, 1:2, 1L, -179.9, 40
  )
  expect_identical(selected$lon_idx, 1L)
  expect_equal(.normalise_longitude(c(-160, 200)), c(-160, -160))
})

match_one_cell <- function(lon_axis, lat_axis, lon_values, lat_values,
                           lon_target, lat_target, spatial_match,
                           max_distance_km) {
  match_cells <- getFromNamespace(".match_horizontal_cells", "ocean3d")
  match_row <- getFromNamespace(".horizontal_match_row", "ocean3d")
  match_row(match_cells(
    lon_axis, lat_axis, lon_values, lat_values, lon_target, lat_target,
    spatial_match, max_distance_km
  ), 1L)
}

test_that("horizontal selection records accepted and rejected candidates", {
  inside <- list(
    centres = 0, lower = 0, upper = 0, longitude = FALSE
  )
  longitude <- inside
  longitude$longitude <- TRUE
  accepted <- match_one_cell(
    longitude, inside, 0, 0, 0.01, 0, "nearest", 2
  )
  expect_true(accepted$accepted)
  expect_identical(accepted$reason, "matched")
  expect_equal(accepted$selected_lon, 0)
  expect_equal(accepted$selected_lat, 0)
  expect_gt(accepted$distance_km, 1)
  rejected <- match_one_cell(
    longitude, inside, 0, 0, 0.01, 0, "nearest", 0.5
  )
  expect_false(rejected$accepted)
  expect_identical(rejected$reason, "distance_limit_exceeded")
  expect_gt(rejected$distance_km, 1)
  outside <- match_one_cell(
    longitude, inside, 0, 0, 0.01, 0, "within", NULL
  )
  expect_identical(outside$reason, "outside_grid")
  invalid <- match_one_cell(
    longitude, inside, 0, 0, NA_real_, 0, "nearest", 5
  )
  expect_identical(invalid$reason, "invalid_coordinates")
})

test_that("cell containment is separate from overall coverage", {
  lon <- list(centres = c(0, 1), lower = c(-0.01, -0.5),
              upper = c(0.01, 1.5), longitude = TRUE)
  lat <- list(centres = 0, lower = 0, upper = 0, longitude = FALSE)
  result <- match_one_cell(lon, lat, c(0, 1), 0, 0.1, 0,
                           "nearest", 100)
  expect_true(result$accepted)
  expect_true(result$inside_coverage)
  expect_false(result$inside_selected_cell)
  expect_equal(result$lon_idx, 1L)
  lon <- list(centres = c(179, 181), lower = c(178, 180),
              upper = c(180, 182), longitude = TRUE)
  wrapped <- match_one_cell(lon, lat, c(179, -179), 0, 181, 0,
                            "nearest", 100)
  expect_true(wrapped$inside_selected_cell)
  expect_true(wrapped$inside_coverage)
})

test_that("great-circle distances match sf under s2", {
  original <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(original)), add = TRUE)
  suppressMessages(sf::sf_use_s2(TRUE))
  a <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326)
  points <- rbind(c(1, 1), c(179.5, -60), c(-179.9, 89), c(0.01, 0))
  for (k in seq_len(nrow(points))) {
    b <- sf::st_sfc(sf::st_point(points[k, ]), crs = 4326)
    expect_equal(
      .great_circle_km(0, 0, points[k, 1], points[k, 2]),
      as.numeric(sf::st_distance(a, b)) / 1000,
      tolerance = 1e-12
    )
  }
})

test_that("matching does not depend on the sf s2 setting", {
  original <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(original)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))
  lon <- list(centres = c(0, 1), lower = c(-0.5, 0.5), upper = c(0.5, 1.5),
              longitude = TRUE)
  lat <- list(centres = c(40, 41), lower = c(39.5, 40.5),
              upper = c(40.5, 41.5), longitude = FALSE)
  for (spatial_match in c("within", "nearest")) {
    matches <- .match_horizontal_cells(lon, lat, c(0, 1), c(40, 41),
                                       c(0.1, 1.4), c(40, 40),
                                       spatial_match, 100)
    expect_identical(matches$lon_idx, 1:2)
    expect_identical(matches$lat_idx, c(1L, 1L))
  }
})

test_that("whole-grid nearest search matches a search over every cell", {
  set.seed(42)
  lon <- c(seq(-179.5, 179.5, by = 7), 180)
  lat <- sort(c(-90, runif(20, -89, 89), 90))
  x <- c(runif(60, -360, 360), 180, -180, 179.99, 0, 10)
  y <- c(runif(60, -90, 90), 0, 89.9, -89.9, 90, -90)
  fast <- .nearest_horizontal_cells(lon, lat, x, y)
  for (k in seq_along(x)) {
    slow <- .nearest_horizontal_cell(lon, lat, seq_along(lon),
                                     seq_along(lat), x[k], y[k])
    expect_equal(fast[[k]]$distance_km, slow$distance_km, tolerance = 1e-12)
    # Away from the poles the selected cell is the same too.
    if (abs(y[k]) < 89) {
      expect_identical(fast[[k]][c("lon_idx", "lat_idx")],
                       slow[c("lon_idx", "lat_idx")])
    }
  }
  invalid <- .nearest_horizontal_cells(lon, lat, c(NA, 0, Inf), c(0, 91, 0))
  expect_true(all(vapply(invalid, is.null, logical(1))))
  expect_length(.nearest_horizontal_cells(lon, lat, numeric(), numeric()), 0L)
})

test_that("exact distance ties keep the first cell in storage order", {
  batch <- .nearest_horizontal_cells(c(0, 0), 0, c(0, 0), c(0, 0))
  expect_identical(batch[[1L]]$lon_idx, 1L)
  expect_identical(batch[[2L]]$lon_idx, 1L)
  expect_equal(batch[[1L]]$distance_km, 0)
  # A point on a shared edge is equally far from both centres; rounding
  # noise must not decide, so the first (west/south) cell is kept.
  lon <- seq(-4.5, 4.5, 1)
  lat <- seq(35.5, 44.5, 1)
  for (y in lat) {
    edge <- .nearest_horizontal_cell(lon, lat, 4:5, which(lat == y), -1, y)
    expect_identical(edge$lon_idx, 4L)
    nearest <- .nearest_horizontal_cells(lon, lat, -1, y)[[1L]]
    expect_identical(nearest$lon_idx, 4L)
  }
  # A point just off the edge still gets the truly closer cell.
  off <- .nearest_horizontal_cell(lon, lat, 4:5, 1L, -1 + 1e-6, 35.5)
  expect_identical(off$lon_idx, 5L)
})

test_that("nearest candidates preserve matching limits and containment", {
  lon <- list(centres = c(0, 1), lower = c(-0.5, 0.5),
              upper = c(0.5, 1.5), longitude = TRUE)
  lat <- list(centres = 40, lower = 39, upper = 41, longitude = FALSE)
  x <- c(0, 0.1, 1.6, NA)
  y <- rep(40, 4)
  nearest <- .nearest_horizontal_cells(c(0, 1), 40, x, y)
  distance <- nearest[[3L]]$distance_km
  for (limit in c(0, 1, 100, distance, distance - 1e-6)) {
    matches <- .match_horizontal_cells(lon, lat, c(0, 1), 40, x, y,
                                       "nearest", limit)
    for (i in 1:3) {
      expect_identical(matches$lon_idx[i], nearest[[i]]$lon_idx)
      expect_identical(matches$accepted[i],
                       nearest[[i]]$distance_km <= limit)
    }
    expect_identical(matches$reason[4L], "invalid_coordinates")
  }
  expect_false(matches$inside_coverage[3L])
  expect_true(matches$inside_selected_cell[2L])
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
  expect_error(.prepare_depth_axis(c(0, 50, 50), "m"), "duplicate")
  expect_error(.prepare_depth_axis(c(0, -50), "m", "down"), "reference surface")
  expect_error(.prepare_depth_axis(c(0, 50), "m", "up"), "reference surface")
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
  expect_error(.read_depth_axis(nc, "depth", c(0, 50, 50)), "duplicate")
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
  axis <- .prepare_spatial_axis(nc, "lon", lon, longitude = TRUE)
  expect_identical(axis$basis, "explicit_bounds")
  expect_identical(.containing_axis_indices(axis, 120), 201L)
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
