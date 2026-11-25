# extract() on a SpatProxy (#163): values at many points at once, which must
# match the per-shape loop of #115 point by point. Synthetic files from
# make_grid_file() (helper-netcdf-grid.R).

# The value the per-shape loop gives each point: crop, read its time step,
# mask with `method`, and take the one value kept (or NA).
loop_values <- function(p, pts, method = "nearest") {
  vapply(seq_len(length(pts)), function(i) {
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
  expect_named(out, c("feature_id", "value", "time", "depth"))
  expect_identical(out$feature_id, as.data.frame(pts)$feature_id)
  # Each point's step and level, as the loop reads them.
  i <- 3L
  v <- as_voxel(crop(p, pts[i]), time = pts[i])[[1L]]
  expect_identical(out$time[i], terra::time(v)[1L])
  expect_true(out$depth[i] %in% depths(v))
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
  expect_identical(extract(p, pts)$depth, 0)
  expect_true(is.na(extract(p, pts, method = "nearest_level")$value))
  expect_identical(extract(p, pts, method = "bottom")$depth, 50)
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
  by_rast <- vapply(seq_len(length(pts)), function(i) {
    r <- terra::rast(crop(flat, pts[i]), time = pts[i])[[1L]]
    kept <- terra::values(mask(r, pts[i]))
    kept[!is.na(kept)]
  }, numeric(1))
  expect_identical(out$value, by_rast)
  expect_true(all(is.na(out$depth)))
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
  expect_true(is.na(out$time[4L]))

  # A point without a time on a source with several steps.
  untimed <- as_vector3d(data.frame(lon = 1.2, lat = 10.6, depth = 5))
  expect_warning(out <- extract(p, untimed), "1 of 1 point has no time",
                 class = "ocean3d_warning_point_time")
  expect_true(is.na(out$value))
  # A source with one step reads it for every point, whatever its time.
  one <- grid_proxy(lon = seq(0.5, 4.5, 1), lat = c(10.5, 11.5),
                    depth = c(0, 50), time = 0)
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
