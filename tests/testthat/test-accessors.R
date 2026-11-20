test_that("depths() reads depths from voxels, rasters and layer names", {
  v <- as_voxel(make_multidepth_rast(depths = c(0, 100, 200, 300)))
  r <- make_multidepth_rast(depths = c(0, 50, 100, 500))
  expect_false(methods::is(r, "SpatVoxel"))

  cases <- list(
    list(label = "SpatVoxel", x = v, depths = c(0, 100, 200, 300)),
    list(label = "plain SpatRaster with conforming names", x = r,
         depths = c(0, 50, 100, 500)),
    list(label = "layer names", x = c("t_an_depth=0", "t_an_depth=1500"),
         depths = c(0, 1500)),
    list(label = "terra::global() row names",
         x = rownames(terra::global(v, "mean", na.rm = TRUE)),
         depths = c(0, 100, 200, 300)),
    list(label = "a name that doesn't follow the convention",
         x = c("t_an_depth=0", "nope"), depths = c(0, NA)),
    # Negative depths are rejected by as_voxel(), but the parser itself is
    # the shared primitive and must report what the name actually says.
    list(label = "decimal and negative values",
         x = c("x_depth=2.5", "x_depth=-5"), depths = c(2.5, -5))
  )
  for (case in cases) {
    expect_identical(depths(case$x), case$depths, label = case$label)
  }
})

test_that("depths() errors on input it can't read depths from", {
  r <- terra::rast(nrows = 2, ncols = 2)
  terra::values(r) <- 1:4
  names(r) <- "not_a_depth_layer"
  convention <- "'\\{variable\\}_depth=\\{value\\}'"
  cases <- list(
    list(label = "no name follows the convention", x = "nope",
         pattern = convention),
    list(label = "no layer name follows the convention", x = r,
         pattern = convention),
    list(label = "integer vector", x = 1:3, pattern = "SpatVoxel"),
    list(label = "NULL", x = NULL, pattern = "SpatVoxel")
  )
  for (case in cases) {
    expect_error(depths(case$x), case$pattern, label = case$label)
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
