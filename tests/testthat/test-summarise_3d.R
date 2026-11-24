# summarise_3d() (#116): statistics of a voxel's values across its cells and
# depth levels. make_grid()'s LAEA cells are 1 km x 1 km, so volume weights
# are the slab thicknesses.

# A 3 x 3 voxel on two levels, 0 and 100 m, whose slabs are 50 m each under
# "midpoint": cell 1 is NA at 0 m, cell 9 at both levels.
two_levels <- function() {
  make_voxel(cbind(c(NA, 2:8, NA), c(11:18, NA)), c(0, 100), varname = "temp")
}

test_that("statistics run over every value, weighted by volume or equally", {
  # One cell at seven levels, slabs 5, 10, 10, 10, 35, 80 and 50 m thick.
  r <- terra::rast(make_grid(1, 1), nlyrs = 7)
  terra::values(r) <- t(c(20, 19, 18, 17, 16, 8, 4))
  v <- as_voxel(r, depths = c(0, 10, 20, 30, 40, 100, 200), varname = "temp")
  values <- c(20, 19, 18, 17, 16, 8, 4)
  thickness <- c(5, 10, 10, 10, 35, 80, 50)

  weighted <- summarise_3d(v, statistic = c("mean", "sd", "min", "max",
                                            "median"))
  expect_identical(names(weighted),
                   c("time", "mean", "sd", "min", "max", "median", "n"))
  expect_equal(weighted$mean, sum(values * thickness) / sum(thickness))
  w <- thickness / sum(thickness)
  expect_equal(weighted$sd,
               sqrt(sum(w * (values - weighted$mean)^2) / sum(w * (1 - w))))
  expect_identical(weighted[c("min", "max")], data.frame(min = 4, max = 20))
  # Half the water (100 m of 200 m) is reached within the 8 degree level.
  expect_identical(weighted$median, 8)
  expect_identical(weighted$n, 7L)
  expect_true(is.na(weighted$time))

  # Equal weights agree with base R.
  equal <- summarise_3d(v, statistic = c("mean", "sd", "median"),
                        weighting = "equal")
  expect_equal(equal$mean, mean(values))
  expect_equal(equal$sd, stats::sd(values))
  expect_equal(equal$median, stats::median(values))

  # The columns follow the order asked for.
  expect_identical(names(summarise_3d(v, statistic = c("max", "mean"))),
                   c("time", "max", "mean", "n"))
})

test_that("the median averages the middle pair under equal weights", {
  v <- make_voxel(cbind(c(1, 2, 3, 4, rep(NA, 5)), rep(NA, 9)), c(0, 100))
  expect_identical(summarise_3d(v, "median", weighting = "equal")$median, 2.5)
  expect_identical(summarise_3d(v, "median")$median, 2.5)
})

test_that("values and weights pair up across cells and levels", {
  v <- two_levels()
  # 50 m slabs everywhere, so the weighted mean equals the equal one here.
  out <- summarise_3d(v)
  expect_equal(out$mean, mean(c(2:8, 11:18)))
  expect_identical(out$n, 15L)
  expect_equal(out$mean, summarise_3d(v, weighting = "equal")$mean)

  # `bounds` moves the weights: under "top" the 100 m level has no height, so
  # only the 0 m level counts.
  expect_equal(summarise_3d(v, bounds = "top")$mean, mean(2:8))
  # Stored edges win over `bounds`: 0 m stands for 10 m and 100 m for 90 m.
  edged <- as_voxel(v, depth_bounds = cbind(c(0, 10), c(10, 100)))
  expect_equal(summarise_3d(edged)$mean,
               (10 * sum(2:8) + 90 * sum(11:18)) / (10 * 7 + 90 * 8))
})

test_that("by = 'depth' gives one row per level and 'time' one per voxel", {
  v <- two_levels()
  by_depth <- summarise_3d(v, statistic = c("mean", "max"), by = "depth")
  expect_identical(names(by_depth),
                   c("time", "depth", "mean", "max", "n"))
  expect_identical(by_depth$depth, c(0, 100))
  expect_equal(by_depth$mean, c(mean(2:8), mean(11:18)))
  expect_identical(by_depth$n, c(7L, 8L))
  expect_identical(summarise_3d(v, by = "time"), summarise_3d(v))
})

test_that("a list of voxels gives one row per time step, with its time", {
  v <- two_levels()
  later <- v + 1
  terra::time(v) <- rep(as.Date("2020-01-01"), 2)
  terra::time(later) <- rep(as.POSIXct("2020-01-02 06:00", tz = "UTC"), 2)

  out <- summarise_3d(list(v, later))
  expect_equal(out$time, as.POSIXct(c("2020-01-01 00:00", "2020-01-02 06:00"),
                                    tz = "UTC"))
  expect_equal(out$mean, summarise_3d(v)$mean + c(0, 1))

  # Steps aren't pooled with by = "depth": two rows per step.
  by_depth <- summarise_3d(list(v, later), by = "depth")
  expect_identical(by_depth$depth, c(0, 100, 0, 100))
  expect_identical(format(by_depth$time, "%d"), c("01", "01", "02", "02"))
})

test_that("a voxel with no values gives NA with n = 0", {
  out <- summarise_3d(make_voxel(NA, c(0, 100)), c("mean", "sd", "median"))
  expect_identical(out[c("mean", "sd", "median", "n")],
                   data.frame(mean = NA_real_, sd = NA_real_,
                              median = NA_real_, n = 0L))
  # A single value has no standard deviation.
  one <- make_voxel(cbind(c(5, rep(NA, 8)), rep(NA, 9)), c(0, 100))
  expect_identical(summarise_3d(one, "sd")$sd, NA_real_)
})

test_that("values with no volume are weighted equally", {
  # A single level without edges stands for a zero-thickness slab, so every
  # volume weight is 0; the values still get their equal-weighted summary.
  one_level <- make_voxel(cbind(c(1, 2, 3, 4, rep(NA, 5))), 50)
  out <- summarise_3d(one_level, c("mean", "sd", "median"))
  expect_equal(out[c("mean", "sd", "median", "n")],
               data.frame(mean = 2.5, sd = stats::sd(1:4), median = 2.5,
                          n = 4L))
  # Under "top" the deepest level has no height either.
  by_depth <- summarise_3d(two_levels(), c("mean", "sd"), by = "depth",
                           bounds = "top")
  expect_equal(by_depth$mean, c(mean(2:8), mean(11:18)))
  expect_equal(by_depth$sd, c(stats::sd(2:8), stats::sd(11:18)))
  # The weighted statistic itself falls back to equal weights.
  expect_identical(.weighted_statistic(c(1, 3), c(0, 0), "mean"), 2)
})

test_that("a value with no volume doesn't join the median's middle pair", {
  # Under "top" the 200 m level has no height: weights 100, 100 and 0 m, so
  # the median is halfway between the two values with weight, 1 and 3.
  r <- terra::rast(make_grid(1, 1), nlyrs = 3)
  terra::values(r) <- t(c(1, 3, 2))
  v <- as_voxel(r, depths = c(0, 100, 200))
  expect_identical(summarise_3d(v, "median", bounds = "top")$median, 2)
})

test_that("volumes and weights stay doubles when terra writes to disk", {
  # terra writes intermediates to temporary files for large rasters or when
  # memory is short, in 4-byte floats by default; the answers mustn't change.
  e <- make_range_rast(rep(0, 9), rep(200, 9))
  v <- envelope_to_voxel(e, depths = c(0, 100, 200))
  field <- two_levels()
  in_memory <- list(volume(e), volume(v),
                    volume(v, per_cell = TRUE, per_level = TRUE),
                    summarise_3d(field, c("mean", "sd")))
  old <- terra::terraOptions(print = FALSE)$todisk
  terra::terraOptions(todisk = TRUE)
  on.exit(terra::terraOptions(todisk = old), add = TRUE)
  on_disk <- list(volume(e), volume(v),
                  volume(v, per_cell = TRUE, per_level = TRUE),
                  summarise_3d(field, c("mean", "sd")))
  expect_identical(on_disk[[1]], in_memory[[1]])
  expect_identical(on_disk[[2]], in_memory[[2]])
  expect_identical(on_disk[[2]], on_disk[[1]])
  expect_identical(terra::values(on_disk[[3]]), terra::values(in_memory[[3]]))
  expect_identical(on_disk[[4]], in_memory[[4]])
})

test_that("a voxel whose time is a year or month stops, as in as_proxy()", {
  v <- two_levels()
  terra::time(v, tstep = "years") <- rep(2020L, 2)
  expect_error(summarise_3d(v), "is a year",
               class = "ocean3d_error_voxel_time")
  expect_error(as_proxy(v), class = "ocean3d_error_voxel_time")
})

test_that("summarise_3d() rejects input it can't summarise", {
  v <- two_levels()
  # One object is checked as every 3D input is.
  expect_error(summarise_3d(make_range_rast(rep(0, 9), rep(100, 9))),
               "got SpatEnvelope", class = "ocean3d_error_not_3d")
  expect_error(summarise_3d(methods::as(v, "SpatRaster")), "got SpatRaster",
               class = "ocean3d_error_not_3d")
  cases <- list(
    list(x = list(v, 1), args = list(), pattern = "nothing else"),
    list(x = list(), args = list(), pattern = "one or more SpatVoxels"),
    list(x = v, args = list(by = "cell"), pattern = "`by` must be NULL")
  )
  for (case in cases) {
    expect_error(do.call(summarise_3d, c(list(case$x), case$args)),
                 case$pattern, class = "ocean3d_error_summarise_input")
  }
  expect_error(summarise_3d(v, statistic = "mode"), "should be one of")
  # A tagged voxel that is no longer valid is rejected as elsewhere, in a
  # list too.
  expect_error(summarise_3d(v[[c(1, 1)]]), class = "ocean3d_error_invalid_3d")
  expect_error(summarise_3d(list(v, v[[c(1, 1)]])),
               class = "ocean3d_error_invalid_3d")
})
