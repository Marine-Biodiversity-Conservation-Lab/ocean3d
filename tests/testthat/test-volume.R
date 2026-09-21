# Helpers are defined locally so this does not depend on definitions that live
# in another test file. The LAEA projection makes every cell 1 km x 1 km, so the
# expected volumes below are hand-computable.
make_grid <- function(ncol = 3, nrow = 3) {
  terra::rast(nrows = nrow, ncols = ncol,
              xmin = 0, xmax = ncol * 1000,
              ymin = 0, ymax = nrow * 1000,
              crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
}

# Helper: a SpatEnvelope with per-cell depth limits (NA marks absence).
make_range_rast <- function(depth_min_vals, depth_max_vals, ncol = 3, nrow = 3) {
  dmin <- make_grid(ncol, nrow)
  terra::values(dmin) <- depth_min_vals
  names(dmin) <- "depth_min"

  dmax <- terra::rast(dmin)
  terra::values(dmax) <- depth_max_vals
  names(dmax) <- "depth_max"

  as_envelope(c(dmin, dmax))
}

# Helper: a SpatVoxel from a cells x depths matrix of values (NA = unoccupied).
make_voxel <- function(vals, depths, ncol = 3, nrow = 3) {
  r <- terra::rast(make_grid(ncol, nrow), nlyrs = length(depths))
  terra::values(r) <- vals
  as_voxel(r, depths = depths, varname = "presence")
}

test_that("calc_volume_overlap() overlap volume never exceeds volume of either input range", {
  # Range A: covers all 9 cells, 0-100m depth
  a <- make_range_rast(
    depth_min_vals = rep(0, 9),
    depth_max_vals = rep(100, 9)
  )

  # Range B: covers only the center 3 cells, 50-200m depth
  b <- make_range_rast(
    depth_min_vals = c(NA, NA, NA, 50, 50, 50, NA, NA, NA),
    depth_max_vals = c(NA, NA, NA, 200, 200, 200, NA, NA, NA)
  )

  result <- calc_volume_overlap(a, b)
  vol_layers <- result[[c("volume_a", "volume_b", "volume_overlap")]]
  totals <- terra::global(vol_layers, "sum", na.rm = TRUE)$sum
  names(totals) <- c("volume_a", "volume_b", "volume_overlap")

  expect_lte(totals[["volume_overlap"]], totals[["volume_a"]])
  expect_lte(totals[["volume_overlap"]] <= totals[["volume_b"]])
  expect_gt(totals[["volume_overlap"]] > 0)
})

test_that("calc_volume_overlap() overlap is NA where ranges don't spatially overlap", {
  # A: top 3 cells
  a <- make_range_rast(
    depth_min_vals = c(0, 0, 0, NA, NA, NA, NA, NA, NA),
    depth_max_vals = c(100, 100, 100, NA, NA, NA, NA, NA, NA)
  )

  # B: bottom 3 cells
  b <- make_range_rast(
    depth_min_vals = c(NA, NA, NA, NA, NA, NA, 0, 0, 0),
    depth_max_vals = c(NA, NA, NA, NA, NA, NA, 100, 100, 100)
  )

  result <- calc_volume_overlap(a, b)

  # No cell has both ranges present, so every overlap cell is NA -- distinct
  # from the 0 a cell gets when both are present but vertically disjoint.
  expect_true(all(is.na(terra::values(result[["volume_overlap"]]))))
  expect_true(all(is.na(terra::values(result[["depth_min_overlap"]]))))
})

test_that("calc_volume_overlap() overlap is zero when depth ranges don't overlap", {
  # Same cells, but A is 0-50m and B is 100-200m
  a <- make_range_rast(
    depth_min_vals = rep(0, 9),
    depth_max_vals = rep(50, 9)
  )

  b <- make_range_rast(
    depth_min_vals = rep(100, 9),
    depth_max_vals = rep(200, 9)
  )

  result <- calc_volume_overlap(a, b)
  overlap_total <- terra::global(result[["volume_overlap"]], "sum", na.rm = TRUE)$sum

  expect_equal(overlap_total, 0)
  # Both present but vertically disjoint: volume is 0, the depth limits are NA.
  expect_true(all(terra::values(result[["volume_overlap"]]) == 0))
  expect_true(all(is.na(terra::values(result[["depth_max_overlap"]]))))
})

test_that("calc_volume_overlap() overlap is zero for adjacent depth ranges, as voxels too", {
  # A ends exactly where B begins. The envelopes share no water, and neither
  # do their voxels: the boundary level belongs to B alone.
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(100, 9), rep(200, 9))
  depths <- c(0, 100, 200)

  env_overlap <- terra::values(calc_volume_overlap(a, b)[["volume_overlap"]])
  expect_true(all(env_overlap == 0))

  va <- envelope_to_voxel(a, depths)
  vb <- envelope_to_voxel(b, depths)
  vox_overlap <- terra::values(calc_volume_overlap(va, vb)[["volume_overlap"]])
  expect_true(all(vox_overlap == 0))

  # the mixed path promotes the envelope and must agree
  mixed_overlap <- terra::values(calc_volume_overlap(a, vb)[["volume_overlap"]])
  expect_true(all(mixed_overlap == 0))
})

test_that("calc_volume_overlap() full overlap when ranges are identical", {
  a <- make_range_rast(
    depth_min_vals = rep(0, 9),
    depth_max_vals = rep(100, 9)
  )

  result <- calc_volume_overlap(a, a)
  vol_layers <- result[[c("volume_a", "volume_b", "volume_overlap")]]
  totals <- terra::global(vol_layers, "sum", na.rm = TRUE)$sum
  names(totals) <- c("volume_a", "volume_b", "volume_overlap")

  expect_equal(totals[["volume_overlap"]], totals[["volume_a"]])
  expect_equal(totals[["volume_overlap"]], totals[["volume_b"]])
})

test_that("volume() returns correct value for uniform grid", {
  # 9 cells, each 1km x 1km, 100m depth = 9 * 1 * 0.1 = 0.9 km³
  r <- make_range_rast(
    depth_min_vals = rep(0, 9),
    depth_max_vals = rep(100, 9)
  )

  vol <- volume(r)

  # Cell area depends on projection; just check it's positive and reasonable
  expect_gt(vol, 0)
  # With 1km cells: 9 cells * 1 km² * 0.1 km depth = 0.9 km³
  expect_equal(vol, 0.9, tolerance = 0.01)
})

test_that("volume() of a vect_to_envelope() output is positive and finite", {
  skip_if_not_installed("sf")

  grid <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5,
                      ymin = 0, ymax = 5, vals = NA, crs = "EPSG:4326")
  seafloor <- grid
  terra::values(seafloor) <- 500
  poly <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(1, 1), c(4, 1), c(4, 4), c(1, 4), c(1, 1)))),
    crs = "EPSG:4326"
  ))

  out <- vect_to_envelope(poly, grid, depth_min = 0,
                          depth_max = list(200, seafloor))

  expect_s4_class(out, "SpatEnvelope")
  v <- volume(out)
  expect_true(is.finite(v))
  expect_gt(v, 0)
})

# ---- SpatVoxel volume -------------------------------------------------------

test_that("volume() sums slab thicknesses over occupied voxels", {
  # 9 cells, occupied at 0 m and 100 m but not at 300 m.
  v <- make_voxel(
    vals = cbind(rep(1, 9), rep(1, 9), rep(NA, 9)),
    depths = c(0, 100, 300)
  )

  # "top": the level names the top of its slab, so 0 m stands for 0-100 m and
  # 100 m for 100-300 m => 300 m of occupied water per cell.
  expect_equal(volume(v), 9 * 0.3, tolerance = 0.01)

  # "midpoint": edges fall halfway to each neighbour, so 0 m stands for 0-50 m
  # and 100 m for 50-200 m => 200 m per cell.
  expect_equal(volume(v, bounds = "midpoint"), 9 * 0.2, tolerance = 0.01)
})

test_that("volume() excludes interior gaps that an envelope would fill", {
  # Occupied at 0, 200 and 300 m, with a gap at 100 m.
  v <- make_voxel(
    vals = cbind(rep(1, 9), rep(NA, 9), rep(1, 9), rep(1, 9)),
    depths = c(0, 100, 200, 300)
  )

  # Slabs are 100, 100, 100 and 0 m, so the occupied levels give 100 + 100 + 0.
  expect_equal(volume(v), 9 * 0.2, tolerance = 0.01)

  # voxel_to_envelope() is lossy in exactly this way: it spans 0-300 m solid.
  expect_equal(volume(voxel_to_envelope(v)), 9 * 0.3, tolerance = 0.01)
  expect_lt(volume(v), volume(voxel_to_envelope(v)))
})

test_that("volume() round-trips through envelope_to_voxel() on aligned levels", {
  e <- make_range_rast(
    depth_min_vals = rep(0, 9),
    depth_max_vals = rep(200, 9)
  )

  # Under "top" the levels are 0, 100 and 200 m with slabs 100, 100 and 0 m,
  # so the discretization neither loses nor invents any water here.
  v <- envelope_to_voxel(e, depths = c(0, 100, 200))
  expect_equal(volume(v), volume(e))
})

test_that("volume() honours a custom occupancy predicate", {
  v <- make_voxel(
    vals = cbind(rep(0.9, 9), rep(0.1, 9), rep(0.9, 9)),
    depths = c(0, 100, 200)
  )

  # Default !is.na() counts all three levels (100 + 100 + 0 m). A threshold
  # drops the middle one, leaving only the 100 m slab the surface level names.
  expect_equal(volume(v), 9 * 0.2, tolerance = 0.01)
  expect_equal(volume(occupied(v, function(x) x > 0.5)), 9 * 0.1,
               tolerance = 0.01)
  expect_equal(volume(occupied(v, function(x) x > 0.99)), 0)
})

# ---- SpatVoxel and mixed overlap --------------------------------------------

test_that("calc_volume_overlap() on voxels reports occupied volume, not spanned volume", {
  depths <- c(0, 100, 200, 300)

  # A is solid over the top three levels; B has a gap at 100 m.
  a <- make_voxel(cbind(rep(1, 9), rep(1, 9), rep(1, 9), rep(NA, 9)), depths)
  b <- make_voxel(cbind(rep(1, 9), rep(NA, 9), rep(1, 9), rep(NA, 9)), depths)

  result <- calc_volume_overlap(a, b)
  vals <- terra::values(result)

  # Depth layers are summary bounds: B spans 0-200 m despite the gap.
  expect_true(all(vals[, "depth_min_b"] == 0))
  expect_true(all(vals[, "depth_max_b"] == 200))

  # Volumes are not: B occupies 100 + 100 m, A occupies 100 + 100 + 100 m,
  # and the overlap is the two levels they share.
  expect_equal(terra::global(result[["volume_a"]], "sum", na.rm = TRUE)$sum,
               9 * 0.3, tolerance = 0.01)
  expect_equal(terra::global(result[["volume_b"]], "sum", na.rm = TRUE)$sum,
               9 * 0.2, tolerance = 0.01)
  expect_equal(terra::global(result[["volume_overlap"]], "sum", na.rm = TRUE)$sum,
               9 * 0.2, tolerance = 0.01)
})

test_that("calc_volume_overlap() overlap never exceeds either voxel's volume", {
  depths <- c(0, 50, 100, 150)
  a <- make_voxel(cbind(rep(1, 9), rep(1, 9), rep(NA, 9), rep(NA, 9)), depths)
  b <- make_voxel(
    cbind(c(NA, NA, NA, 1, 1, 1, NA, NA, NA),
          c(NA, NA, NA, 1, 1, 1, NA, NA, NA),
          c(NA, NA, NA, 1, 1, 1, NA, NA, NA),
          rep(NA, 9)),
    depths
  )

  result <- calc_volume_overlap(a, b)
  totals <- terra::global(
    result[[c("volume_a", "volume_b", "volume_overlap")]], "sum", na.rm = TRUE
  )$sum

  expect_lte(totals[3], totals[1])
  expect_lte(totals[3], totals[2])
  expect_gt(totals[3], 0)
})

test_that("calc_volume_overlap() promotes an envelope onto the voxel's levels", {
  depths <- c(0, 50, 100, 150, 200)

  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(50, 9), rep(200, 9))
  av <- envelope_to_voxel(a, depths = depths)
  bv <- envelope_to_voxel(b, depths = depths)

  all_voxel <- terra::values(calc_volume_overlap(av, bv))

  # Either argument may be the envelope; both give the all-voxel answer.
  expect_equal(terra::values(calc_volume_overlap(a, bv)), all_voxel)
  expect_equal(terra::values(calc_volume_overlap(av, b)), all_voxel)
})

test_that("calc_volume_overlap() rejects voxels on different depth levels", {
  a <- make_voxel(cbind(rep(1, 9), rep(1, 9)), depths = c(0, 100))
  b <- make_voxel(cbind(rep(1, 9), rep(1, 9)), depths = c(0, 200))

  expect_error(calc_volume_overlap(a, b), "same depth levels")
})

# ---- output class -----------------------------------------------------------

test_that("calc_volume_overlap() returns a plain SpatRaster", {
  # terra keeps the input's class on `[[` and `c()`, so the nine-layer stack
  # used to come back tagged SpatEnvelope while failing that class's rules.
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(50, 9), rep(200, 9))
  depths <- c(0, 50, 100, 150, 200)
  va <- envelope_to_voxel(a, depths)
  vb <- envelope_to_voxel(b, depths)

  for (out in list(calc_volume_overlap(a, b), calc_volume_overlap(va, vb),
                   calc_volume_overlap(a, vb), calc_volume_overlap(va, b))) {
    expect_identical(class(out)[[1]], "SpatRaster")
    expect_equal(terra::nlyr(out), 9L)
  }
})

# ---- input contract ---------------------------------------------------------

test_that("volume functions reject bare SpatRasters", {
  dmin <- make_grid()
  terra::values(dmin) <- rep(0, 9)
  names(dmin) <- "depth_min"
  dmax <- terra::rast(dmin)
  terra::values(dmax) <- rep(100, 9)
  names(dmax) <- "depth_max"

  # Carrying the right layer names is not the same as being a SpatEnvelope.
  bare <- c(dmin, dmax)
  expect_s4_class(bare, "SpatRaster")
  expect_false(methods::is(bare, "SpatEnvelope"))

  expect_error(volume(bare), "not a bare SpatRaster")
  expect_error(calc_volume_overlap(bare, bare), "not a bare SpatRaster")

  # Mixed with a real envelope, in either position.
  e <- as_envelope(bare)
  expect_error(calc_volume_overlap(e, bare), "not a bare SpatRaster")
  expect_error(calc_volume_overlap(bare, e), "not a bare SpatRaster")
})

test_that("volume functions reject inputs on different grids", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(0, 16), rep(100, 16), ncol = 4, nrow = 4)

  expect_error(calc_volume_overlap(a, b), "same grid")
})

test_that("voxel-only arguments are an error for an envelope", {
  e <- make_range_rast(rep(0, 9), rep(100, 9))

  expect_error(volume(e, bounds = "midpoint"), "unused argument")
  expect_error(calc_volume_overlap(e, e, bounds = "midpoint"), "unused argument")
})
