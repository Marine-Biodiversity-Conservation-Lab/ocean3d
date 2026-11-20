test_that("depths() reads the depth axis of a SpatVoxel", {
  v <- as_voxel(make_multidepth_rast(depths = c(0, 100, 200, 300)))
  expect_identical(depths(v), c(0, 100, 200, 300))
})

test_that("depths() reads terra::depth() of a plain SpatRaster", {
  r <- make_multidepth_rast(depths = c(0, 50, 100, 500))
  expect_false(methods::is(r, "SpatVoxel"))
  expect_identical(depths(r), c(0, 50, 100, 500))
  # A depth set to 0 m counts as set.
  surface <- terra::rast(nrows = 2, ncols = 2, vals = 1)
  terra::depth(surface) <- 0
  expect_identical(depths(surface), 0)
})

test_that("depths() never reads layer names", {
  # Names in the `{variable}_depth={value}` pattern but no terra::depth().
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2, vals = 1)
  names(r) <- c("t_depth=0", "t_depth=100")
  err <- expect_error(depths(r), class = "ocean3d_error_no_depths")
  expect_match(conditionMessage(err), "no depths in terra::depth\\(\\)")
  expect_match(conditionMessage(err), "Layer names aren't read as depths")

  # Renamed layers keep their depths.
  v <- as_voxel(make_multidepth_rast(depths = c(0, 100, 200, 300)))
  names(v) <- c("t_depth=9", "a", "b", "c")
  expect_identical(depths(v), c(0, 100, 200, 300))
})

test_that("depths() rejects input that isn't a raster", {
  for (x in list(1:3, NULL, "t_an_depth=0")) {
    err <- expect_error(depths(x), class = "ocean3d_error_depths_input")
    expect_match(conditionMessage(err), "must be a SpatVoxel or a SpatRaster")
  }
})

test_that("depth_bounds() reads a voxel's edges for its current layers", {
  edges <- cbind(c(0, 1, 2), c(1, 2, 4))
  r <- terra::rast(make_grid(), nlyrs = 3, vals = 1)
  # Given deepest first: the table comes back shallow to deep.
  v <- as_voxel(r, depths = c(3, 0.5, 1.5), depth_bounds = edges[c(3, 1, 2), ])

  b <- depth_bounds(v)
  expect_s3_class(b, "data.frame")
  expect_identical(names(b), c("depth", "lower", "upper"))
  expect_equal(b$depth, depths(v))
  expect_equal(b$lower, edges[, 1])
  expect_equal(b$upper, edges[, 2])

  # A subset gives the rows for the layers it kept, and the table goes back
  # into as_voxel() as it is.
  expect_equal(depth_bounds(v[[2:3]])$depth, c(1.5, 3))
  again <- as_voxel(r, depths = c(0.5, 1.5, 3), depth_bounds = b)
  expect_equal(depth_bounds(again), b)
})

test_that("depth_bounds() is NULL without edges for every layer", {
  plain <- as_voxel(terra::rast(make_grid(), nlyrs = 2, vals = 1),
                    depths = c(0.5, 1.5))
  expect_null(depth_bounds(plain))

  top <- as_voxel(terra::rast(make_grid(), nlyrs = 1, vals = 1), depths = 0.2,
                  depth_bounds = cbind(0, 0.3))
  expect_null(depth_bounds(c(top, plain)))
})

test_that("depth_bounds() needs a SpatVoxel", {
  r <- terra::rast(make_grid(), nlyrs = 2, vals = 1)
  names(r) <- c("t_depth=0", "t_depth=100")
  err <- expect_error(depth_bounds(r), class = "ocean3d_error_not_voxel")
  expect_match(conditionMessage(err), "must be a SpatVoxel; got SpatRaster")
})
