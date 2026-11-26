# extract_to_buffer() is a wrapper over as_vector3d(), buffer_3d(), crop(),
# as_voxel() or rast(), mask() and summarise_3d() (#116), whose own tests
# cover their rules. These tests cover what the wrapper adds: how each
# observation becomes a buffer, which time steps it reads, how a time window
# is pooled, and the output columns.

# A 10 x 10 grid of 1 degree cells at depths 5, 15, 30 and 60 m on three
# daily steps (1 to 3 January 2020). Each value is its position in storage
# order, lon fastest, then lat, depth and time, so it can be worked out by
# hand: cell (i, j), level k, step t holds i + 10 (j - 1) + 100 (k - 1) +
# 400 (t - 1).
buffer_file <- function(depth = c(5, 15, 30, 60), ...) {
  make_grid_file(lon = seq(0.5, 9.5, 1), lat = seq(0.5, 9.5, 1),
                 depth = depth, time = c(0, 1, 2), ...)
}
at <- function(i, j, k, t, levels = 4L) {
  i + 10 * (j - 1) + 10 * 10 * (k - 1) + 10 * 10 * levels * (t - 1)
}
# An observation in cell (6, 6) at 20 m on 2 January.
obs <- function(...) {
  data.frame(lon = 5.2, lat = 5.3, depth = 20,
             date = as.Date("2020-01-02"), ...)
}

test_that("without a buffer, the observation's cell and method's level count", {
  path <- buffer_file()
  on.exit(unlink(path))
  # 20 m: the nearest level is 15 m (level 2); the bottom is 60 m, the
  # surface 5 m.
  cases <- list(
    list(method = NULL, level = 2L),
    list(method = "nearest", level = 2L),
    list(method = "nearest_level", level = 2L),
    list(method = "bottom", level = 4L),
    list(method = "surface", level = 1L)
  )
  for (case in cases) {
    out <- extract_to_buffer(obs(), path, method = case$method)
    label <- if (is.null(case$method)) "NULL" else case$method
    expect_equal(out$temp, at(6, 6, case$level, 2), label = label)
    expect_identical(out$temp_n, 1L, label = label)
  }
})

test_that("a depth window takes the levels whose slabs it overlaps", {
  path <- buffer_file()
  on.exit(unlink(path))
  # 10 to 30 m overlaps the slabs of 15 m (10 to 22.5) and 30 m (22.5 to 45).
  values <- at(6, 6, 2:3, 2)
  weighted <- extract_to_buffer(obs(), path, depth_radius_m = 10)
  expect_equal(weighted$temp, sum(values * c(12.5, 22.5)) / 35)
  equal <- extract_to_buffer(obs(), path, depth_radius_m = 10,
                             weighting = "equal")
  expect_equal(equal$temp, mean(values))
  expect_identical(equal$temp_n, 2L)
  expect_equal(
    extract_to_buffer(obs(), path, depth_radius_m = 10, statistic = "max")$temp,
    max(values)
  )
})

test_that("a circle takes every cell it touches", {
  path <- buffer_file()
  on.exit(unlink(path))
  # 60 km from (5.2, 5.3) reaches into the cells west and south of (6, 6).
  out <- extract_to_buffer(obs(), path, horizontal_radius_km = 60,
                           method = "nearest_level", weighting = "equal")
  cells <- rbind(c(5, 5), c(6, 5), c(5, 6), c(6, 6))
  expect_identical(out$temp_n, nrow(cells))
  expect_equal(out$temp, mean(at(cells[, 1], cells[, 2], 2, 2)))
})

test_that("a time window pools the values of every step inside it", {
  path <- buffer_file()
  on.exit(unlink(path))
  pooled <- extract_to_buffer(obs(), path, time_before = 1, time_after = 1,
                              weighting = "equal")
  expect_identical(pooled$temp_n, 3L)
  expect_equal(pooled$temp, mean(at(6, 6, 2, 1:3)))
  # A day after only: steps 2 and 3, with the whole of the last day.
  after <- extract_to_buffer(obs(), path, time_after = 1, statistic = "min")
  expect_equal(after$temp, at(6, 6, 2, 2))
  expect_identical(after$temp_n, 2L)
  # Hours on date-times.
  timed <- obs()
  timed$date <- as.POSIXct("2020-01-02 06:00", tz = "UTC")
  hours <- extract_to_buffer(timed, path, time_before = 12, time_after = 12,
                             time_unit = "hours", weighting = "equal")
  expect_equal(hours$temp, at(6, 6, 2, 2))
})

test_that("a variable without depth reads one layer per step", {
  path <- buffer_file(depth = NULL)
  on.exit(unlink(path))
  out <- extract_to_buffer(obs(), path)
  expect_equal(out$temp, at(6, 6, 1, 2, levels = 1L))
  pooled <- extract_to_buffer(obs(), path, time_before = 1,
                              weighting = "equal", method = "2d")
  expect_equal(pooled$temp, mean(at(6, 6, 1, 1:2, levels = 1L)))
  expect_error(extract_to_buffer(obs(), path, method = "bottom"),
               "has no depth", class = "ocean3d_error_buffer_input")
})

test_that("observations it can't place get NA, with one warning", {
  path <- buffer_file()
  on.exit(unlink(path))
  data <- data.frame(lon = c(5.2, 50, 5.2, 5.2), lat = 5.3,
                     depth = c(20, 20, NA, 20),
                     date = as.Date(c("2020-01-02", "2020-01-02",
                                      "2020-01-02", NA)))
  expect_warning(out <- extract_to_buffer(data, path),
                 "2 of 4 observations", class = "ocean3d_warning_buffer_rows")
  expect_equal(out$temp, c(at(6, 6, 2, 2), NA, NA, NA))
  expect_identical(out$temp_n, c(1L, 0L, 0L, 0L))
  # A window with no step in it.
  late <- obs()
  late$date <- as.Date("2020-02-01")
  expect_identical(extract_to_buffer(late, path, time_before = 1)$temp_n, 0L)
})

test_that("a series, a proxy and sf input give the same values", {
  paths <- c(
    make_grid_file(lon = seq(0.5, 9.5, 1), lat = seq(0.5, 9.5, 1),
                   depth = c(5, 15), time = 0:1, values = 1),
    make_grid_file(lon = seq(0.5, 9.5, 1), lat = seq(0.5, 9.5, 1),
                   depth = c(5, 15), time = 2, values = 3)
  )
  on.exit(unlink(paths))
  # Steps 1 and 2 hold 1, step 3 holds 3.
  pooled <- extract_to_buffer(obs(), paths, time_after = 1,
                              weighting = "equal")
  expect_equal(pooled$temp, 2)
  expect_equal(
    extract_to_buffer(obs(), as_proxy(paths), time_after = 1,
                      weighting = "equal")$temp,
    2
  )
  points <- sf::st_as_sf(obs(), coords = c("lon", "lat"), crs = 4326)
  from_sf <- extract_to_buffer(points, paths, time_after = 1,
                               weighting = "equal")
  expect_s3_class(from_sf, "sf")
  expect_equal(from_sf$temp, 2)
})

test_that("the output column is named and never overwrites", {
  path <- buffer_file()
  on.exit(unlink(path))
  out <- extract_to_buffer(obs(), path, output_col = "thetao")
  expect_named(out, c(names(obs()), "thetao", "thetao_n"))
  expect_error(extract_to_buffer(obs(temp = 1), path),
               "already has a column named temp",
               class = "ocean3d_error_buffer_input")
})

test_that("extract_to_buffer() rejects requests it can't serve", {
  path <- buffer_file()
  on.exit(unlink(path))
  cases <- list(
    list(args = list(horizontal_radius_km = -1), pattern = "0 or more"),
    list(args = list(time_before = 1, max_time_diff = 1),
         pattern = "without a time window"),
    list(args = list(depth_radius_m = 10, method = "bottom"),
         pattern = "Leave `method` NULL"),
    list(args = list(method = "deepest"), pattern = "must be one of"),
    list(args = list(time_before = 6, time_unit = "hours"),
         pattern = "must\\s+be in days"),
    list(args = list(depth_col = "z"), pattern = "needs the column")
  )
  for (case in cases) {
    expect_error(do.call(extract_to_buffer, c(list(obs(), path), case$args)),
                 case$pattern, class = "ocean3d_error_buffer_input")
  }
  # Dropped options are no longer arguments.
  expect_error(extract_to_buffer(obs(), path, statistic = "most_similar"))
  expect_error(extract_to_buffer(obs(), path, min_valid_fraction = 0.5))
})
