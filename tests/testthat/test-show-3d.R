# show() for SpatVoxel and SpatEnvelope (#61): the class and its 3D
# information first, then terra's description of the grid.

shown <- function(x) utils::capture.output(print(x))

test_that("a voxel prints its class, variable, depths, edges and time", {
  v <- as_voxel(make_multidepth_rast(), depths = c(0, 100, 200, 300),
                varname = "temp")
  out <- shown(v)
  expect_match(out[1L], "^class +: SpatVoxel$")
  expect_match(out, "^variable +: temp$", all = FALSE)
  # terra sets no variable name for a raster built in R.
  expect_match(shown(as_voxel(make_multidepth_rast())), "^variable +: none$",
               all = FALSE)
  expect_match(out, "^depth \\(m\\) +: 0, 100, 200, 300 \\(4 levels\\)$",
               all = FALSE)
  expect_match(out, "^depth edges : none; `bounds` applies$", all = FALSE)
  expect_match(out, "^time +: none$", all = FALSE)
  # terra's lines follow, except the ones the 3D lines replace.
  expect_match(out, "^size +: 2, 2, 4", all = FALSE)
  expect_match(out, "^names +: temp_depth=0", all = FALSE)
  expect_no_match(out, "SpatRaster")
  # terra's own `depth` and `time` lines are replaced, not repeated.
  expect_no_match(out, "^depth +:")
  expect_identical(sum(grepl("^time +:", out)), 1L)

  # Edges, a time step, and many levels shown as a range.
  edged <- as_voxel(v, depth_bounds = cbind(c(0, 50, 150, 250),
                                            c(50, 150, 250, 350)))
  terra::time(edged) <- rep(as.POSIXct("2020-01-02 06:00", tz = "UTC"), 4)
  out <- shown(edged)
  expect_match(out, "^depth edges : from depth_bounds$", all = FALSE)
  expect_match(out, "^time +: 2020-01-02 06:00:00 UTC$", all = FALSE)

  many <- make_voxel(matrix(1, nrow = 9, ncol = 7),
                     depths = c(0, 10, 20, 30, 40, 50, 60))
  expect_match(shown(many), "^depth \\(m\\) +: 0 to 60 \\(7 levels\\)$",
               all = FALSE)
  expect_match(shown(make_voxel(matrix(1, 9, 1), depths = 50)),
               "^depth \\(m\\) +: 50 \\(1 level\\)$", all = FALSE)
})

test_that("an envelope prints its class, depth range and occupied cells", {
  e <- make_range_rast(c(0, 10, NA, 20, 0, NA, 0, 0, 0),
                       c(100, 50, NA, 300, 100, NA, 100, 100, 100))
  out <- shown(e)
  expect_match(out[1L], "^class +: SpatEnvelope$")
  expect_match(out, "^depth \\(m\\) +: 0 to 300$", all = FALSE)
  expect_match(out, "^cells +: 7 of 9 occupied$", all = FALSE)
  expect_match(out, "^names +: depth_min, depth_max$", all = FALSE)
  expect_no_match(out, "SpatRaster")

  # An empty envelope has no depth range.
  empty <- make_range_rast(rep(NA_real_, 9), rep(NA_real_, 9))
  out <- shown(empty)
  expect_match(out, "^depth \\(m\\) +: none$", all = FALSE)
  expect_match(out, "^cells +: 0 of 9 occupied$", all = FALSE)
})

test_that("printing a tagged object that is no longer valid doesn't stop", {
  v <- as_voxel(make_multidepth_rast(depths = c(0, 100, 200, 300)))
  # `[[` keeps the class tag; two layers at one depth aren't a voxel.
  out <- shown(v[[c(1, 1)]])
  expect_match(out[1L], "^class +: SpatVoxel$")
  expect_match(out[2L], "^status +: not a valid SpatVoxel: .*distinct depth")
  expect_no_match(out, "^depth \\(m\\)")

  e <- make_range_rast(rep(0, 9), rep(100, 9))
  names(e) <- c("a", "b")
  out <- shown(e)
  expect_match(out[2L], "^status +: not a valid SpatEnvelope: layers must")
  expect_no_match(out, "^cells")
})
