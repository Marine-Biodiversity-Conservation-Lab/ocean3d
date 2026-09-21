test_that("depths() reads the depth axis of a SpatVoxel", {
  v <- as_voxel(make_multidepth_rast(depths = c(0, 100, 200, 300)))
  expect_identical(depths(v), c(0, 100, 200, 300))
})

test_that("depths() works on a plain SpatRaster with conforming names", {
  r <- make_multidepth_rast(depths = c(0, 50, 100, 500))
  expect_false(methods::is(r, "SpatVoxel"))
  expect_identical(depths(r), c(0, 50, 100, 500))
})

test_that("depths() accepts a character vector of layer names", {
  expect_identical(depths(c("t_an_depth=0", "t_an_depth=1500")), c(0, 1500))
})

test_that("depths() recovers depths from terra::global() row names", {
  v <- as_voxel(make_multidepth_rast(depths = c(0, 100, 200, 300)))
  per_depth <- terra::global(v, "mean", na.rm = TRUE)
  expect_identical(depths(rownames(per_depth)), c(0, 100, 200, 300))
})

test_that("depths() returns NA for names that do not follow the convention", {
  expect_identical(depths(c("t_an_depth=0", "nope")), c(0, NA))
})

test_that("depths() errors when no name follows the convention", {
  expect_error(depths("nope"), "'\\{variable\\}_depth=\\{value\\}'")

  r <- terra::rast(nrows = 2, ncols = 2)
  terra::values(r) <- 1:4
  names(r) <- "not_a_depth_layer"
  expect_error(depths(r), "'\\{variable\\}_depth=\\{value\\}'")
})

test_that("depths() parses decimal and negative values in names", {
  # Negative depths are rejected by as_voxel(), but the parser itself is the
  # shared primitive and must report what the name actually says.
  expect_identical(depths(c("x_depth=2.5", "x_depth=-5")), c(2.5, -5))
})

test_that("depths() rejects input that is neither a raster nor character", {
  expect_error(depths(1:3), "SpatVoxel")
  expect_error(depths(NULL), "SpatVoxel")
})
