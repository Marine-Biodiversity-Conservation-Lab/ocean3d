# extract() on a SpatProxy (#163): values at many points at once, which must
# match the per-shape loop of #115 point by point. Synthetic files from
# make_grid_file() (helper-netcdf-grid.R).

# The value the per-shape loop gives each point: crop, read its time step,
# mask with `method`, and take the one value kept (or NA).
loop_values <- function(p, pts, method = "nearest") {
  n <- length(pts)
  vapply(seq_len(n), function(i) {
    s <- pts[i]
    v <- tryCatch(as_voxel(crop(p, s), time = s)[[1L]],
                  ocean3d_error_proxy_crop = function(e) NULL)
    if (is.null(v)) return(NA_real_)
    kept <- terra::values(mask(v, s, method = method))
    kept <- kept[!is.na(kept)]
    if (length(kept)) kept else NA_real_
  }, numeric(1))
}

grid_proxy <- function(..., values = NULL) {
  args <- list(...)
  if (is.null(values)) {
    set.seed(1)
    # An axis left out (NULL) counts as one.
    n <- prod(pmax(1L, lengths(args[c("lon", "lat", "depth", "time")])))
    values <- stats::runif(n)
  }
  path <- do.call(make_grid_file, c(args, list(values = values)))
  as_proxy(path)
}

timed_points <- function(n, lon, lat, depth = c(0, 120), days = 4) {
  set.seed(2)
  as_vector3d(data.frame(
    lon = stats::runif(n, lon[1L], lon[2L]),
    lat = stats::runif(n, lat[1L], lat[2L]),
    depth = stats::runif(n, depth[1L], depth[2L]),
    t = as.POSIXct("2020-01-01", tz = "UTC") + stats::runif(n, 0, days * 86400)
  ), time = "t")
}

test_that("extract() gives the per-shape loop's values for every method", {
  p <- grid_proxy(lon = seq(0.5, 9.5, 1), lat = seq(10.5, 14.5, 1),
                  depth = c(0, 10, 50, 100), time = 0:4)
  pts <- timed_points(40, c(0, 10), c(10, 15))
  for (method in c("nearest", "nearest_level", "bottom", "surface")) {
    out <- extract(p, pts, method = method)
    expect_identical(out$value, loop_values(p, pts, method), label = method)
  }
  out <- extract(p, pts)
  expect_named(out, c("feature_id", "value", "time_step", "depth_level"))
  expect_identical(out$feature_id, as.data.frame(pts)$feature_id)
  # Each point's step and level, as the loop reads them.
  i <- 3L
  v <- as_voxel(crop(p, pts[i]), time = pts[i])[[1L]]
  expect_identical(out$time_step[i], terra::time(v)[1L])
  expect_true(out$depth_level[i] %in% depths(v))
})

test_that("extract() reads levels with no value as `method` says", {
  # Every second level is empty, so the nearest level and the nearest level
  # with a value differ.
  values <- rep(c(1, NA), each = 10 * 5)
  p <- grid_proxy(lon = seq(0.5, 9.5, 1), lat = seq(10.5, 14.5, 1),
                  depth = c(0, 10, 50, 100), time = NULL,
                  values = rep(values, 2))
  pts <- as_vector3d(data.frame(lon = 2.2, lat = 12.2, depth = 12))
  # By default, the nearest level with a value, as extract_to_point() picks
  # it; "nearest_level" takes the empty 10 m level.
  expect_identical(extract(p, pts)$depth_level, 0)
  expect_true(is.na(extract(p, pts, method = "nearest_level")$value))
  expect_identical(extract(p, pts, method = "bottom")$depth_level, 50)
})

test_that("extract() matches a series of files and a 2D variable", {
  first <- make_grid_file(lon = seq(0.5, 4.5, 1), lat = c(10.5, 11.5),
                          depth = c(0, 50), time = 0:1,
                          values = stats::runif(5 * 2 * 2 * 2))
  second <- make_grid_file(lon = seq(0.5, 4.5, 1), lat = c(10.5, 11.5),
                           depth = c(0, 50), time = 2:3,
                           values = stats::runif(5 * 2 * 2 * 2))
  series <- as_proxy(c(first, second))
  pts <- timed_points(20, c(0, 5), c(10, 12), depth = c(0, 60), days = 3)
  expect_identical(extract(series, pts)$value, loop_values(series, pts))

  flat <- grid_proxy(lon = seq(0.5, 4.5, 1), lat = c(10.5, 11.5),
                     depth = NULL, time = 0:3)
  out <- extract(flat, pts)
  n <- length(pts)
  by_rast <- vapply(seq_len(n), function(i) {
    r <- terra::rast(crop(flat, pts[i]), time = pts[i])[[1L]]
    kept <- terra::values(mask(r, pts[i]))
    kept[!is.na(kept)]
  }, numeric(1))
  expect_identical(out$value, by_rast)
  expect_true(all(is.na(out$depth_level)))
  # A cell of a 2D variable holds one value, so `method` doesn't apply.
  expect_error(extract(flat, pts, method = "bottom"),
               "method", class = "ocean3d_error_unused_argument")
})

test_that("extract() reads points either side of the antimeridian", {
  p <- grid_proxy(lon = seq(-179.5, 179.5, 1), lat = c(-0.5, 0.5),
                  depth = c(0, 50), time = 0:1)
  pts <- as_vector3d(data.frame(
    lon = c(179.7, -179.7, 179.2), lat = c(0.2, -0.2, 0.4), depth = 10,
    t = as.POSIXct(c("2020-01-01", "2020-01-02", "2020-01-02"), tz = "UTC")
  ), time = "t")
  expect_identical(extract(p, pts)$value, loop_values(p, pts))

  # A 0 to 360 grid with points given from -180 to 180.
  east <- grid_proxy(lon = seq(0.5, 359.5, 1), lat = c(-0.5, 0.5),
                     depth = c(0, 50), time = 0:1)
  west <- as_vector3d(data.frame(
    lon = c(-170.3, 10.4, -0.6), lat = c(0.2, -0.2, 0.4), depth = 10,
    t = as.POSIXct(c("2020-01-01", "2020-01-02", "2020-01-02"), tz = "UTC")
  ), time = "t")
  expect_identical(extract(east, west)$value, loop_values(east, west))
})

test_that("points off the grid or the time axis give NA", {
  p <- grid_proxy(lon = seq(0.5, 4.5, 1), lat = c(10.5, 11.5),
                  depth = c(0, 50), time = 0:2)
  pts <- as_vector3d(data.frame(
    lon = c(1.2, 40, 41, 2.5), lat = c(10.6, 30, 30, 11.2), depth = 5,
    t = as.POSIXct(c("2020-01-01", "2020-01-02", "2020-01-02", "2020-02-01"),
                   tz = "UTC")
  ), time = "t")
  # Points 2 and 3 are off the grid and read together; point 4 is a month
  # after the last step.
  out <- extract(p, pts, max_time_diff = 2)
  expect_false(is.na(out$value[1L]))
  expect_true(all(is.na(out$value[2:4])))
  expect_true(is.na(out$time_step[4L]))

  # A point without a time on a source with several steps.
  untimed <- as_vector3d(data.frame(lon = 1.2, lat = 10.6, depth = 5))
  expect_warning(out <- extract(p, untimed), "1 of 1 point has no time",
                 class = "ocean3d_warning_point_time")
  expect_true(is.na(out$value))
  half_timed <- as_vector3d(data.frame(
    lon = c(1.2, 2.2), lat = 10.6, depth = 5,
    t = as.POSIXct(c("2020-01-01", NA), tz = "UTC")
  ), time = "t")
  expect_warning(out <- extract(p, half_timed), "1 of 2 points has no time",
                 class = "ocean3d_warning_point_time")
  expect_identical(is.na(out$value), c(FALSE, TRUE))
  # A one-step climatology reads its step for every point, whatever its
  # time.
  one <- grid_proxy(lon = seq(0.5, 4.5, 1), lat = c(10.5, 11.5),
                    depth = c(0, 50), time = 6,
                    time_units = "months since 0000-01-01 00:00:00")
  expect_false(is.na(extract(one, untimed)$value))
})

test_that("a small read budget splits the reads without changing values", {
  p <- grid_proxy(lon = seq(0.5, 9.5, 1), lat = seq(10.5, 14.5, 1),
                  depth = c(0, 10, 50, 100), time = 0:4)
  pts <- timed_points(40, c(0, 10), c(10, 15))
  whole <- extract(p, pts)
  rlang::local_options(ocean3d.max_read_mb = 1e-5)
  expect_identical(extract(p, pts), whole)
})

test_that("a point on a cell edge takes the cell its own crop gives", {
  # On a 1 degree grid, (1, 11) lies on the edge between cells. Alone or
  # with a point far to the north-east, it takes the same cell as in the
  # per-shape loop.
  p <- grid_proxy(lon = seq(0.5, 4.5, 1), lat = seq(10.5, 12.5, 1),
                  depth = 0, time = NULL,
                  values = seq_len(5 * 3))
  alone <- as_vector3d(data.frame(lon = 1, lat = 11, depth = 0))
  paired <- as_vector3d(data.frame(lon = c(1, 4.5), lat = c(11, 12.5),
                                   depth = 0))
  expect_identical(extract(p, alone)$value, loop_values(p, alone))
  expect_identical(extract(p, paired)$value[1L], extract(p, alone)$value)

  # Points on whole degrees, at the normal and a tiny budget.
  grid <- expand.grid(lon = 0:5, lat = 10:13)
  on_edges <- as_vector3d(cbind(grid, depth = 0))
  want <- loop_values(p, on_edges)
  expect_identical(extract(p, on_edges)$value, want)
  rlang::local_options(ocean3d.max_read_mb = 1e-5)
  expect_identical(extract(p, on_edges)$value, want)
})

test_that("extract() finds points on a projected grid", {
  # A voxel on 1 km LAEA cells near (0, 0), and points given in lon/lat.
  r <- terra::rast(nrows = 4, ncols = 5, xmin = 0, xmax = 5000, ymin = 0,
                   ymax = 4000, nlyrs = 2,
                   crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
  terra::values(r) <- seq_len(40)
  p <- as_proxy(as_voxel(r, depths = c(0, 50)))
  metres <- cbind(c(1200, 3700, 4400, 600), c(500, 2600, 3900, 3100))
  ll <- terra::crds(terra::project(terra::vect(metres, crs = terra::crs(r)),
                                   "EPSG:4326"))
  pts <- as_vector3d(data.frame(lon = ll[, 1], lat = ll[, 2],
                                depth = c(5, 45, 5, 30)))
  want <- loop_values(p, pts)
  expect_false(anyNA(want))
  expect_identical(extract(p, pts)$value, want)
  rlang::local_options(ocean3d.max_read_mb = 1e-5)
  expect_identical(extract(p, pts)$value, want)
})

test_that("a point at 180 takes the same cell alone and grouped", {
  # On a -180 to 180 grid, 180 is the grid's eastern edge; alone, a point
  # there is read across the seam and takes the westernmost column.
  p <- grid_proxy(lon = seq(-179.5, 179.5, 1), lat = c(-0.5, 0.5),
                  depth = 0, time = NULL,
                  values = seq_len(360 * 2))
  alone <- as_vector3d(data.frame(lon = 180, lat = 0.2, depth = 0))
  paired <- as_vector3d(data.frame(lon = c(180, 179.7), lat = c(0.2, 0.3),
                                   depth = 0))
  want <- loop_values(p, alone)
  expect_identical(extract(p, alone)$value, want)
  expect_identical(extract(p, paired)$value[1L], want)
  rlang::local_options(ocean3d.max_read_mb = 1e-5)
  expect_identical(extract(p, paired)$value[1L], want)
})

test_that("points far apart are read as small windows, not the box between", {
  p <- grid_proxy(lon = seq(0.5, 199.5, 1), lat = seq(-59.5, 59.5, 1),
                  depth = c(0, 50), time = NULL)
  far <- as_vector3d(data.frame(lon = c(5.2, 190.7), lat = c(-50.3, 52.6),
                                depth = 10))
  read_cells <- integer(0)
  reader <- .read_proxy_window
  testthat::local_mocked_bindings(.read_proxy_window = function(x, ...) {
    sizes <- x@variable$dim_lengths
    read_cells <<- c(read_cells, prod(sizes[c("lon", "lat")]))
    reader(x, ...)
  })
  out <- extract(p, far)
  # Two reads of a few cells each, not one read of the 186 x 104 box.
  expect_length(read_cells, 2L)
  expect_true(all(read_cells <= 9L))
  expect_identical(out$value, loop_values(p, far))
})

test_that("bind = TRUE returns the points with their values added", {
  p <- grid_proxy(lon = seq(0.5, 4.5, 1), lat = c(10.5, 11.5),
                  depth = c(0, 50), time = 0:2)
  pts <- as_vector3d(data.frame(
    id = c("c", "a", "b"), tag = c("T1", "T2", "T3"),
    lon = c(3.2, 1.4, 40), lat = c(11.2, 10.6, 30), depth = c(40, 5, 5),
    t = as.POSIXct(c("2020-01-02", "2020-01-01", "2020-01-03"), tz = "UTC")
  ), time = "t", id = "id")
  short <- extract(p, pts)
  bound <- extract(p, pts, bind = TRUE)
  expect_s4_class(bound, "SpatVector3D")
  table <- as.data.frame(bound)
  # The same features in the same order, with their own columns.
  expect_identical(table$feature_id, c("c", "a", "b"))
  expect_identical(table$tag, c("T1", "T2", "T3"))
  expect_identical(table$depth, c(40, 5, 5))
  expect_identical(table$time, as.data.frame(pts)$time)
  # The values, the step used and the level used, as the short form gives.
  expect_identical(table$value, short$value)
  expect_identical(table$time_step, short$time_step)
  expect_identical(table$depth_level, short$depth_level)

  # A column named like an added one stops.
  clash <- as_vector3d(data.frame(lon = 1.4, lat = 10.6, depth = 5,
                                  value = 9, t = as.POSIXct("2020-01-01",
                                                            tz = "UTC")),
                       time = "t")
  expect_error(extract(p, clash, bind = TRUE), "value",
               class = "ocean3d_error_vector3d")
  expect_error(extract(p, pts, bind = NA), "`bind` must be TRUE or FALSE",
               class = "ocean3d_error_extract_bind")
})

test_that("extract() rejects input it can't use", {
  p <- grid_proxy(lon = seq(0.5, 4.5, 1), lat = c(10.5, 11.5),
                  depth = c(0, 50), time = NULL)
  line <- as_vector3d(data.frame(id = "a", lon = c(1, 2), lat = 11, depth = 5),
                      type = "lines", id = "id")
  pts <- as_vector3d(data.frame(lon = 1, lat = 11, depth = 5))
  expect_error(extract(p, line), "points",
               class = "ocean3d_error_vector3d_type")
  expect_error(extract(p, pts, method = "deepest"), "must be one of",
               class = "ocean3d_error_mask_method")
  expect_error(extract(p, pts, max_time_diff = -1), "max_time_diff",
               class = "ocean3d_error_voxel_input")
  expect_error(extract(p, pts, bounds = "top"), "bounds",
               class = "ocean3d_error_unused_argument")
})

test_that("a single dated time step is matched to each point's time", {
  f <- make_grid_file(lon = c(0.5, 1.5), lat = c(10.5, 11.5), depth = NULL,
                      time = 0, values = 1:4)
  on.exit(unlink(f), add = TRUE)
  p <- as_proxy(f)
  pts <- as_vector3d(data.frame(
    lon = c(0.5, 1.5), lat = 10.5, depth = 0,
    date = as.Date(c("2020-01-01", "2020-03-01"))
  ), time = "date")
  expect_equal(extract(p, pts)$value, c(1, 2))
  expect_equal(extract(p, pts, time_match = "exact")$value, c(1, NA))
  expect_equal(extract(p, pts, max_time_diff = 1)$value, c(1, NA))
  untimed <- as_vector3d(data.frame(lon = 0.5, lat = 10.5, depth = 0))
  expect_warning(out <- extract(p, untimed), "a dated time step",
                 class = "ocean3d_warning_point_time")
  expect_true(is.na(out$value))
})

test_that("a single unbounded longitude or latitude reads points on it", {
  # A point time series (one cell, three levels) and a transect along one
  # latitude: a point on the coordinate reads it, a point off it gets NA.
  f <- make_grid_file(lon = 150.25, lat = -20.5, depth = c(5, 15, 40),
                      time = 0:2, values = 1:9)
  g <- make_grid_file(lon = c(0.5, 1.5, 2.5), lat = 40.25, depth = NULL,
                      time = 0:1, values = 1:6)
  on.exit(unlink(c(f, g)), add = TRUE)
  series <- as_vector3d(data.frame(
    lon = c(150.25, 150.25, 150.3), lat = -20.5, depth = c(10, 30, 10),
    date = as.Date("2020-01-01") + c(0, 2, 1)
  ), time = "date")
  p <- as_proxy(f)
  expect_equal(extract(p, series, method = "surface")$value, c(1, 7, NA))
  expect_equal(extract(p, series, method = "bottom")$value, c(3, 9, NA))
  expect_equal(extract(p, series)$value, c(1, 9, NA))
  transect <- as_vector3d(data.frame(
    lon = c(0.5, 1.7, 2.5, 1), lat = c(40.25, 40.25, 40.25, 40.3), depth = 0,
    date = as.Date("2020-01-02")
  ), time = "date")
  expect_equal(extract(as_proxy(g), transect)$value, c(4, 5, 6, NA))
  # A voxel needs the cell's width, so as_voxel() still stops.
  expect_error(as_voxel(p, time = 1), "single cell without bounds",
               class = "ocean3d_error_irregular_grid")
})
