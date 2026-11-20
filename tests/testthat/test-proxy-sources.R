# Sources of the core object (#113): a series of netCDF files and a SpatVoxel
# already in memory, besides one file. Synthetic files from make_grid_file()
# (helper-netcdf-grid.R).

# One time step of one file, in fixed lon, lat, depth, time order.
read_step <- function(f, step) {
  .read_proxy_window(as_proxy(f, var = "temp"), list(time = c(step, step)))
}

test_that("a series joins its files' time axes in time order", {
  a <- make_grid_file(time = c(0, 1))
  b <- make_grid_file(time = c(2, 3))
  # In other units: 2020-01-05 and 2020-01-06.
  h <- make_grid_file(time = c(96, 120), time_units = "hours since 2020-01-01")
  on.exit(unlink(c(a, b, h)), add = TRUE)

  # Given out of order; ordered by their first date.
  s <- as_proxy(c(b, h, a), var = "temp")
  expect_identical(s@source$type, "series")
  expect_identical(s@source$path, normalizePath(c(a, b, h)))
  expect_identical(s@source$steps, c(2L, 2L, 2L))
  expect_identical(s@variable$dim_lengths[["time"]], 6L)
  expect_equal(s@axes$time$cftime$as_timestamp(asPOSIX = TRUE),
               seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "day",
                   length.out = 6))
  expect_match(capture.output(print(s))[2], "3 files")

  # A window across files reads each step from its own file.
  w <- .read_proxy_window(s, list(time = c(2, 5)))
  expect_identical(dim(w), c(3L, 2L, 2L, 4L))
  expect_identical(w[, , , 1, drop = FALSE], read_step(a, 2))
  expect_identical(w[, , , 2, drop = FALSE], read_step(b, 1))
  expect_identical(w[, , , 4, drop = FALSE], read_step(h, 1))

  # A voxel from a series, by date.
  v <- as_voxel(s, time = "2020-01-05")
  expect_equal(terra::time(v)[1], as.POSIXct("2020-01-05", tz = "UTC"))

  # Bounds are kept only when every file has them.
  ab <- make_grid_file(time = c(0, 1),
                       time_bounds = rbind(c(-0.5, 0.5), c(0.5, 1.5)))
  on.exit(unlink(ab), add = TRUE)
  expect_identical(as_proxy(c(ab, b), var = "temp")@axes$time$edges, "none")
})

test_that("a series needs one variable on one grid, each step once", {
  a <- make_grid_file(time = c(0, 1))
  overlap <- make_grid_file(time = c(1, 2))
  grid <- make_grid_file(time = c(5, 6), lon = c(0.5, 1.5))
  depth <- make_grid_file(time = c(5, 6), depth = c(5, 25))
  untimed <- make_grid_file(time = NULL)
  clim <- make_grid_file(time = 6,
                         time_units = "months since 0000-01-01 00:00:00")
  on.exit(unlink(c(a, overlap, grid, depth, untimed, clim)), add = TRUE)

  cases <- list(
    list(files = c(a, a), pattern = "more than once"),
    list(files = c(a, overlap), pattern = "repeats or overlaps"),
    list(files = c(a, grid), pattern = "different grid"),
    list(files = c(a, depth), pattern = "different grid"),
    list(files = c(a, untimed), pattern = "no dated time axis"),
    list(files = c(a, clim), pattern = "no dated time axis")
  )
  for (case in cases) {
    err <- expect_no_warning(expect_error(
      as_proxy(case$files, var = "temp"), class = "ocean3d_error_proxy_series"
    ))
    expect_match(conditionMessage(err), case$pattern)
  }
})

test_that("a voxel in memory is a source like a file", {
  f <- make_grid_file(depth_bounds = rbind(c(0, 10), c(10, 20)))
  on.exit(unlink(f), add = TRUE)
  file_proxy <- as_proxy(f, var = "temp")
  v <- as_voxel(file_proxy, time = 2)

  p <- as_proxy(v)
  expect_s4_class(p, "SpatProxy")
  expect_identical(p@source$type, "voxel")
  expect_match(capture.output(print(p))[2], "SpatVoxel in memory")
  # The same cells, depths, depth edges and date.
  expect_equal(sort(p@axes$lon$centres), file_proxy@axes$lon$centres)
  expect_equal(sort(p@axes$lat$centres), file_proxy@axes$lat$centres)
  expect_equal(p@axes$depth$lower, c(0, 10))
  expect_equal(p@axes$depth$upper, c(10, 20))
  expect_equal(p@axes$time$cftime$as_timestamp(asPOSIX = TRUE),
               as.POSIXct("2020-01-02", tz = "UTC"))

  # Its values, in its own storage order (terra's rows run north to south).
  from_voxel <- .read_proxy_window(p)
  from_file <- read_step(f, 2)
  expect_identical(from_voxel[, 2:1, , , drop = FALSE], from_file)

  # as_voxel() gives the voxel back.
  expect_identical(as_voxel(p), v)
  expect_error(as_voxel(p, time = 2), class = "ocean3d_error_proxy_time")

  # A voxel without a date or stored edges: no time axis, and the slabs of
  # its default `bounds`.
  plain <- make_voxel(1, c(0, 100, 300))
  q <- as_proxy(plain)
  expect_null(q@axes$time)
  expect_equal(q@axes$depth$lower, c(0, 50, 200))
  expect_equal(q@axes$depth$upper, c(50, 200, 300))
})
