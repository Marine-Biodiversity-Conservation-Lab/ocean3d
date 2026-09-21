# Helpers are defined locally so this file does not depend on another test
# file. The LAEA projection makes every cell 1 km x 1 km, so volumes below are
# hand-computable.
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
make_voxel <- function(vals, depths, ncol = 3, nrow = 3, varname = "presence") {
  r <- terra::rast(make_grid(ncol, nrow), nlyrs = length(depths))
  terra::values(r) <- vals
  as_voxel(r, depths = depths, varname = varname)
}

# Helper: a single-layer footprint on the grid (NA = outside).
make_footprint <- function(vals, ncol = 3, nrow = 3) {
  terra::setValues(make_grid(ncol, nrow), vals)
}

# Helper: a polygon covering the centres of the left two columns of a 3x3
# grid, in the grid's CRS.
make_left_polygon <- function() {
  terra::vect("POLYGON ((0 0, 2000 0, 2000 3000, 0 3000, 0 0))",
              crs = terra::crs(make_grid()))
}

# ---- intersect_3d: envelopes ------------------------------------------------

test_that("intersect_3d() on envelopes is the per-cell interval intersection", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(50, 9), rep(200, 9))

  out <- intersect_3d(a, b)

  expect_s4_class(out, "SpatEnvelope")
  expect_identical(class(out)[[1]], "SpatEnvelope")
  expect_true(methods::validObject(out))
  expect_identical(terra::values(out)[, "depth_min"], rep(50, 9))
  expect_identical(terra::values(out)[, "depth_max"], rep(100, 9))
})

test_that("intersect_3d() on envelopes is empty where the domains do not meet", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))

  # Touching only: a shared edge is no shared water.
  touching <- make_range_rast(rep(100, 9), rep(200, 9))
  expect_true(all(is.na(terra::values(intersect_3d(a, touching)))))

  # Vertically disjoint.
  deeper <- make_range_rast(rep(300, 9), rep(400, 9))
  expect_true(all(is.na(terra::values(intersect_3d(a, deeper)))))

  # Spatially disjoint: A in the top row, B in the bottom row.
  top <- make_range_rast(c(0, 0, 0, rep(NA, 6)), c(100, 100, 100, rep(NA, 6)))
  bottom <- make_range_rast(c(rep(NA, 6), 0, 0, 0), c(rep(NA, 6), 100, 100, 100))
  expect_true(all(is.na(terra::values(intersect_3d(top, bottom)))))
})

test_that("intersect_3d() on envelopes keeps only cells present in both", {
  a <- make_range_rast(c(0, 0, NA, 0, 0, 0, 0, 0, 0),
                       c(100, 100, NA, 100, 100, 100, 100, 100, 100))
  b <- make_range_rast(c(NA, 50, 50, 50, 50, 50, 50, 50, 50),
                       c(NA, 200, 200, 200, 200, 200, 200, 200, 200))

  out <- intersect_3d(a, b)
  expect_identical(terra::values(out)[, "depth_min"], c(NA, 50, NA, rep(50, 6)))
  expect_identical(terra::values(out)[, "depth_max"], c(NA, 100, NA, rep(100, 6)))
})

test_that("intersect_3d() is symmetric", {
  a <- make_range_rast(c(0, 0, NA, 0, 0, 0, 0, 0, 0),
                       c(100, 100, NA, 100, 100, 100, 100, 100, 100))
  b <- make_range_rast(rep(50, 9), c(200, 60, 200, 200, 200, 200, 200, 200, 200))
  expect_identical(terra::values(intersect_3d(a, b)), terra::values(intersect_3d(b, a)))
})

test_that("intersect_3d() rejects voxel-only arguments for a pair of envelopes", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  expect_error(intersect_3d(a, a, bounds = "top"), "unused argument")
})

# ---- intersect_3d: voxels and mixed input -----------------------------------

test_that("intersect_3d() on voxels is the co-occupied set of levels", {
  depths <- c(0, 100, 200)
  # Cell 1: occupied at all levels in A, at 0 and 200 in B (interior gap).
  # Cell 2: A only. Cell 3: neither.
  a <- make_voxel(cbind(c(1, 1, NA), c(1, 1, NA), c(1, 1, NA)), depths, ncol = 3, nrow = 1)
  b <- make_voxel(cbind(c(1, NA, NA), c(NA, NA, NA), c(1, NA, NA)), depths, ncol = 3, nrow = 1)

  out <- intersect_3d(a, b)

  expect_s4_class(out, "SpatVoxel")
  expect_identical(class(out)[[1]], "SpatVoxel")
  expect_identical(names(out), paste0("presence_depth=", depths))
  expect_identical(depths(out), depths)
  expect_identical(unname(terra::values(out)),
               unname(cbind(c(1, NA, NA), c(NA, NA, NA), c(1, NA, NA))))
})

test_that("intersect_3d() on voxels intersects presence and requires shared depths", {
  depths <- c(0, 100)
  a <- make_voxel(cbind(c(5, 5, 5), c(5, 5, 5)), depths, ncol = 3, nrow = 1, varname = "t")
  b <- make_voxel(cbind(c(1, 10, NA), c(10, 1, NA)), depths, ncol = 3, nrow = 1, varname = "t")

  out <- intersect_3d(occupied(a, function(v) v > 3),
                      occupied(b, function(v) v > 3))
  expect_identical(unname(terra::values(out)),
               unname(cbind(c(NA, 1, NA), c(1, NA, NA))))

  c3 <- make_voxel(cbind(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1)), c(0, 100, 200), ncol = 3, nrow = 1)
  expect_error(intersect_3d(a, c3), "same depth levels")
})

test_that("intersect_3d() on mixed input promotes the envelope by hand's rule", {
  depths <- c(0, 50, 100, 150, 200)
  e <- make_range_rast(rep(0, 9), rep(100, 9))
  v <- make_voxel(matrix(1, nrow = 9, ncol = 5), depths)

  by_hand <- intersect_3d(envelope_to_voxel(e, depths), v)
  expect_identical(terra::values(intersect_3d(e, v)), terra::values(by_hand))
  expect_identical(terra::values(intersect_3d(v, e)), terra::values(by_hand))
  expect_identical(names(intersect_3d(e, v)), paste0("presence_depth=", depths))

  # `bounds` reaches the promotion.
  mid <- intersect_3d(envelope_to_voxel(e, depths, bounds = "midpoint"), v)
  expect_identical(terra::values(intersect_3d(e, v, bounds = "midpoint")),
               terra::values(mid))
})

# ---- intersect_3d: 2D input -------------------------------------------------

test_that("intersect_3d() with a footprint raster restricts an envelope horizontally", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  fp <- make_footprint(c(1, NA, 0, 7, NA, NA, 1, 1, 1))

  out <- intersect_3d(a, fp)
  expect_identical(class(out)[[1]], "SpatEnvelope")
  expect_identical(terra::values(out)[, "depth_min"], c(0, NA, 0, 0, NA, NA, 0, 0, 0))
  expect_identical(terra::values(out)[, "depth_max"], c(100, NA, 100, 100, NA, NA, 100, 100, 100))

  # Mirrored order gives the same object.
  expect_identical(terra::values(intersect_3d(fp, a)), terra::values(out))
})

test_that("intersect_3d() with a footprint raster restricts a voxel horizontally", {
  depths <- c(0, 100)
  v <- make_voxel(cbind(c(3, NA, 3), c(NA, 3, 3)), depths, ncol = 3, nrow = 1, varname = "t")
  fp <- make_footprint(c(1, 1, NA), ncol = 3, nrow = 1)

  out <- intersect_3d(v, fp)
  expect_identical(class(out)[[1]], "SpatVoxel")
  expect_identical(names(out), paste0("presence_depth=", depths))
  # Presence, not the field values, and the third cell is gone at every level.
  expect_identical(unname(terra::values(out)),
               unname(cbind(c(1, NA, NA), c(NA, 1, NA))))
  expect_identical(terra::values(intersect_3d(fp, v)), terra::values(out))
})

test_that("intersect_3d() accepts polygons as SpatVector, sf and sfc", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  poly <- make_left_polygon()
  expected_min <- c(0, 0, NA, 0, 0, NA, 0, 0, NA)

  expect_identical(terra::values(intersect_3d(a, poly))[, "depth_min"], expected_min)
  expect_identical(terra::values(intersect_3d(poly, a))[, "depth_min"], expected_min)

  poly_sf <- sf::st_as_sf(poly)
  expect_identical(terra::values(intersect_3d(a, poly_sf))[, "depth_min"], expected_min)
  expect_identical(terra::values(intersect_3d(a, sf::st_geometry(poly_sf)))[, "depth_min"],
               expected_min)
})

test_that("intersect_3d() projects polygons in another CRS onto the grid", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  poly_ll <- terra::project(make_left_polygon(), "EPSG:4326")

  expect_identical(terra::values(intersect_3d(a, poly_ll))[, "depth_min"],
               terra::values(intersect_3d(a, make_left_polygon()))[, "depth_min"])
})

test_that("intersect_3d() rejects input it cannot read as a footprint", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  fp <- make_footprint(rep(1, 9))

  expect_error(intersect_3d(fp, fp), "neither `x` nor `y`")
  expect_error(intersect_3d(a, "hello"), "SpatVector, sf, sfc")
  expect_error(intersect_3d(a, c(fp, fp)), "single-layer footprint")
  expect_error(intersect_3d(a, make_footprint(rep(1, 4), ncol = 2, nrow = 2)),
               "same grid")
})

test_that("intersect_3d() requires a shared grid for two 3D inputs", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(0, 4), rep(100, 4), ncol = 2, nrow = 2)
  expect_error(intersect_3d(a, b), "same grid")
})

# ---- intersect_3d: composition with volume ----------------------------------

test_that("volume(intersect_3d(a, b)) equals the overlap volume calc_volume_overlap() reports", {
  depths <- c(0, 50, 100, 150, 200)
  e1 <- make_range_rast(c(0, 0, NA, 0, 0, 0, 0, 0, 0),
                        c(100, 100, NA, 100, 100, 100, 120, 100, 100))
  e2 <- make_range_rast(rep(50, 9), c(200, 60, 200, 200, 200, 200, 200, 200, 200))
  v1 <- envelope_to_voxel(e1, depths)
  v2 <- envelope_to_voxel(e2, depths)

  overlap_sum <- function(x, y) {
    terra::global(calc_volume_overlap(x, y)[["volume_overlap"]], "sum",
                  na.rm = TRUE)[[1]]
  }

  expect_identical(volume(intersect_3d(e1, e2)), overlap_sum(e1, e2))
  expect_identical(volume(intersect_3d(v1, v2)), overlap_sum(v1, v2))
  expect_identical(volume(intersect_3d(e1, v2)), overlap_sum(e1, v2))
  expect_identical(volume(intersect_3d(v1, e2)), overlap_sum(v1, e2))
})

# ---- intersects_3d ----------------------------------------------------------

# Four cells covering the whole truth table: 1 = both present, depths meet;
# 2 = both present, depths disjoint; 3 = A only; 4 = neither.
truth_envelopes <- function() {
  list(
    a = make_range_rast(c(0, 0, 0, NA), c(100, 100, 100, NA), ncol = 4, nrow = 1),
    b = make_range_rast(c(50, 300, NA, NA), c(200, 400, NA, NA), ncol = 4, nrow = 1)
  )
}

test_that("intersects_3d() on envelopes gives the tri-state truth table", {
  e <- truth_envelopes()
  out <- intersects_3d(e$a, e$b)

  expect_identical(class(out)[[1]], "SpatRaster")
  expect_true(terra::is.bool(out))
  expect_identical(names(out), "intersects")
  expect_identical(terra::values(out)[, 1], c(TRUE, FALSE, FALSE, NA))
  expect_identical(terra::values(intersects_3d(e$b, e$a))[, 1], c(TRUE, FALSE, FALSE, NA))
})

test_that("intersects_3d() treats touching intervals as not intersecting", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(100, 9), rep(200, 9))
  expect_identical(terra::values(intersects_3d(a, b))[, 1], rep(FALSE, 9))
})

test_that("adjacent depth ranges do not intersect as voxels either", {
  # The envelope rule, touching is not overlapping, must survive promotion to
  # depth levels: [0, 100] and [100, 200] share only the boundary at 100 m,
  # so their voxels occupy no common level and every 3D verb agrees with the
  # envelope answer, whichever combination of types it is asked about.
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(100, 9), rep(200, 9))
  depths <- c(0, 100, 200)
  va <- envelope_to_voxel(a, depths)
  vb <- envelope_to_voxel(b, depths)

  # voxel vs voxel
  expect_identical(terra::values(intersects_3d(va, vb))[, 1], rep(FALSE, 9))
  expect_identical(terra::values(intersects_3d(vb, va))[, 1], rep(FALSE, 9))
  expect_true(all(is.na(terra::values(intersect_3d(va, vb)))))

  # envelope vs voxel, both orders
  expect_identical(terra::values(intersects_3d(a, vb))[, 1], rep(FALSE, 9))
  expect_identical(terra::values(intersects_3d(vb, a))[, 1], rep(FALSE, 9))
  expect_true(all(is.na(terra::values(intersect_3d(a, vb)))))
  expect_true(all(is.na(terra::values(intersect_3d(va, b)))))

  # and the same under the midpoint convention, with levels whose slab edge
  # falls at 100 m (halfway between 50 and 150)
  depths_mid <- c(0, 50, 150, 200)
  va_mid <- envelope_to_voxel(a, depths_mid, bounds = "midpoint")
  vb_mid <- envelope_to_voxel(b, depths_mid, bounds = "midpoint")
  expect_identical(terra::values(intersects_3d(va_mid, vb_mid))[, 1], rep(FALSE, 9))
  expect_identical(terra::values(intersects_3d(a, vb_mid, bounds = "midpoint"))[, 1], rep(FALSE, 9))

  # overlapping by any amount is still an intersection, so this is the
  # boundary itself and not a wider gap
  c_ <- make_range_rast(rep(99, 9), rep(200, 9))
  expect_identical(terra::values(intersects_3d(va, envelope_to_voxel(c_, depths)))[, 1], rep(TRUE, 9))
})

test_that("intersects_3d() on voxels gives the tri-state truth table", {
  depths <- c(0, 100, 200)
  a <- make_voxel(cbind(c(1, 1, 1, NA), c(1, NA, 1, NA), c(NA, NA, 1, NA)),
                  depths, ncol = 4, nrow = 1)
  b <- make_voxel(cbind(c(NA, NA, NA, NA), c(1, NA, NA, NA), c(1, 1, NA, NA)),
                  depths, ncol = 4, nrow = 1)

  out <- intersects_3d(a, b)
  expect_true(terra::is.bool(out))
  expect_identical(class(out)[[1]], "SpatRaster")
  expect_identical(terra::values(out)[, 1], c(TRUE, FALSE, FALSE, NA))
})

test_that("intersects_3d() on mixed input equals promoting the envelope by hand", {
  depths <- c(0, 50, 100, 150, 200, 300, 400)
  e <- truth_envelopes()$a
  v <- envelope_to_voxel(truth_envelopes()$b, depths)

  by_hand <- intersects_3d(envelope_to_voxel(e, depths), v)
  expect_identical(terra::values(intersects_3d(e, v))[, 1], terra::values(by_hand)[, 1])
  expect_identical(terra::values(intersects_3d(v, e))[, 1], terra::values(by_hand)[, 1])
  expect_identical(terra::values(by_hand)[, 1], c(TRUE, FALSE, FALSE, NA))
})

test_that("intersects_3d() is TRUE exactly where intersect_3d() is non-empty", {
  e <- truth_envelopes()
  shared <- intersect_3d(e$a, e$b)
  # nolint start: scalar_in_linter
  # use %in% instead of ==, because NA == TRUE returns NA, compared to NA %in% TRUE returns FALSE
  expect_identical(terra::values(intersects_3d(e$a, e$b))[, 1] %in% TRUE,
               !is.na(terra::values(shared)[, "depth_min"]))
  # nolint end

  depths <- c(0, 50, 100, 150, 200, 300, 400)
  va <- envelope_to_voxel(e$a, depths)
  vb <- envelope_to_voxel(e$b, depths)
  shared_v <- intersect_3d(va, vb)
  # nolint start: scalar_in_linter
  # use %in% instead of ==, because NA == TRUE returns NA, compared to NA %in% TRUE returns FALSE
  expect_identical(terra::values(intersects_3d(va, vb))[, 1] %in% TRUE,
               rowSums(!is.na(terra::values(shared_v))) > 0)
  # nolint end
})

test_that("intersects_3d() is TRUE exactly where calc_volume_overlap() finds overlap", {
  a <- make_range_rast(c(0, 0, 0, 0, 0, 0, NA, NA, NA),
                       c(100, 100, 100, 40, 40, 40, NA, NA, NA))
  b <- make_range_rast(c(50, 50, 50, 50, 50, 50, 0, 0, 0),
                       c(200, 200, 200, 200, 200, 200, 100, 100, 100))

  from_volume <- terra::values(calc_volume_overlap(a, b))[, "depth_min_overlap"]
  # nolint start: scalar_in_linter
  # use %in% instead of ==, because NA == TRUE returns NA, compared to NA %in% TRUE returns FALSE
  expect_identical(terra::values(intersects_3d(a, b))[, 1] %in% TRUE, !is.na(from_volume))
  # nolint end
  # ...and FALSE, not NA, where both are present but disjoint.
  expect_identical(terra::values(intersects_3d(a, b))[, 1], c(rep(TRUE, 3), rep(FALSE, 6)))
})

test_that("intersects_3d() with 2D input tests presence on both sides", {
  a <- make_range_rast(c(0, 0, NA, NA), c(100, 100, NA, NA), ncol = 4, nrow = 1)
  fp <- make_footprint(c(1, NA, 1, NA), ncol = 4, nrow = 1)

  expect_identical(terra::values(intersects_3d(a, fp))[, 1], c(TRUE, FALSE, FALSE, NA))
  expect_identical(terra::values(intersects_3d(fp, a))[, 1], c(TRUE, FALSE, FALSE, NA))

  v <- envelope_to_voxel(a, c(0, 100))
  expect_identical(terra::values(intersects_3d(v, fp))[, 1], c(TRUE, FALSE, FALSE, NA))
  expect_identical(terra::values(intersects_3d(fp, v))[, 1], c(TRUE, FALSE, FALSE, NA))

  # Polygons: the left two columns of the 3x3 grid.
  a9 <- make_range_rast(c(0, 0, 0, NA, NA, NA, 0, 0, 0),
                        c(100, 100, 100, NA, NA, NA, 100, 100, 100))
  expect_identical(terra::values(intersects_3d(a9, make_left_polygon()))[, 1],
               c(TRUE, TRUE, FALSE, FALSE, FALSE, NA, TRUE, TRUE, FALSE))
})

test_that("a stack of intersects_3d() layers sums to a richness map", {
  e <- truth_envelopes()
  stack <- c(intersects_3d(e$a, e$b), intersects_3d(e$a, e$a))
  expect_identical(terra::values(sum(stack, na.rm = TRUE))[, 1], c(2, 1, 1, NA))
})

test_that("intersects_3d() rejects input without a 3D object or on another grid", {
  fp <- make_footprint(rep(1, 9))
  expect_error(intersects_3d(fp, fp), "neither `x` nor `y`")

  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(0, 4), rep(100, 4), ncol = 2, nrow = 2)
  expect_error(intersects_3d(a, b), "same grid")
  expect_error(intersects_3d(a, a, bounds = "top"), "unused argument")
})

# ---- mask -------------------------------------------------------------------

# A temperature field on five levels over a 3x3 grid, every cell a distinct
# value so a misplaced mask would show.
make_field <- function(depths = c(0, 50, 100, 150, 200)) {
  vals <- matrix(seq_len(9 * length(depths)), nrow = 9)
  make_voxel(vals, depths, varname = "temp")
}

test_that("mask(voxel, envelope) equals the hand-written idiom exactly", {
  field <- make_field()
  e <- make_range_rast(c(0, 0, NA, 60, 60, 60, 0, 0, 0),
                       c(100, 30, NA, 120, 500, 120, 200, 200, 200))

  by_hand <- terra::mask(field, envelope_to_voxel(e, depths(field)))
  out <- mask(field, e)

  expect_identical(class(out)[[1]], "SpatVoxel")
  expect_identical(names(out), names(field))
  expect_identical(terra::values(out), terra::values(by_hand))

  # `bounds` reaches the promotion.
  by_hand_mid <- terra::mask(field, envelope_to_voxel(e, depths(field),
                                                      bounds = "midpoint"))
  expect_identical(terra::values(mask(field, e, bounds = "midpoint")),
               terra::values(by_hand_mid))
})

test_that("mask(voxel, envelope) consults the depth axis", {
  # Before: terra's method masked layer by layer and left the voxel unchanged.
  field <- make_field()
  # [0, 75] runs through the 0 and 50 m slabs and stops short of 100.
  shallow <- make_range_rast(rep(0, 9), rep(75, 9))

  out <- mask(field, shallow)
  expect_false(anyNA(terra::values(out[[c("temp_depth=0", "temp_depth=50")]])))
  expect_true(all(is.na(terra::values(out[[c("temp_depth=100", "temp_depth=150",
                                               "temp_depth=200")]]))))
})

test_that("mask(voxel, voxel) requires the same depth levels", {
  # Before: a 3-level mask was applied by position to the first 3 of 5 layers.
  field <- make_field()
  other <- envelope_to_voxel(make_range_rast(rep(0, 9), rep(100, 9)),
                             depths = c(0, 100, 200))
  expect_error(mask(field, other), "same depth levels")

  same <- envelope_to_voxel(make_range_rast(rep(0, 9), c(rep(100, 8), NA)),
                            depths = depths(field))
  out <- mask(field, same)
  expect_identical(class(out)[[1]], "SpatVoxel")
  expect_identical(terra::values(out),
               terra::values(terra::mask(methods::as(field, "SpatRaster"),
                                         methods::as(same, "SpatRaster"))))

  # `fun` decides where the mask is occupied.
  cold <- make_field()
  out_fun <- mask(field, occupied(cold, function(v) v > 20))
  expect_identical(is.na(terra::values(out_fun)), terra::values(cold) <= 20)
})

test_that("mask(envelope, envelope) keeps whole intervals where the domains meet", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(c(50, 300, NA, 50, 50, 50, 50, 50, 50),
                       c(200, 400, NA, 200, 200, 200, 200, 200, 200))

  out <- mask(a, b)
  expect_identical(class(out)[[1]], "SpatEnvelope")
  expect_identical(terra::values(out)[, "depth_min"], c(0, NA, NA, rep(0, 6)))
  expect_identical(terra::values(out)[, "depth_max"], c(100, NA, NA, rep(100, 6)))

  # Before: two envelopes 100 m apart came back unchanged.
  deeper <- make_range_rast(rep(300, 9), rep(400, 9))
  expect_true(all(is.na(terra::values(mask(a, deeper)))))
})

test_that("mask(envelope, voxel) keeps cells where the voxel reaches the envelope", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  v <- envelope_to_voxel(make_range_rast(c(50, 300, NA, rep(50, 6)),
                                         c(200, 400, NA, rep(200, 6))),
                         depths = c(0, 50, 100, 150, 200, 300, 400))
  out <- mask(a, v)
  expect_identical(class(out)[[1]], "SpatEnvelope")
  expect_identical(terra::values(out)[, "depth_min"], c(0, NA, NA, rep(0, 6)))
})

test_that("mask(plain raster, voxel) keeps cells occupied at any depth", {
  # Before: terra read the voxel's first layer, so a cell occupied only at
  # depth was dropped.
  depths <- c(0, 100)
  v <- make_voxel(cbind(c(NA, 1, NA), c(1, NA, NA)), depths, ncol = 3, nrow = 1)
  fp <- make_footprint(c(7, 8, 9), ncol = 3, nrow = 1)

  out <- mask(fp, v)
  expect_identical(class(out)[[1]], "SpatRaster")
  expect_identical(terra::values(out)[, 1], c(7, 8, NA))
  expect_identical(terra::values(mask(fp, occupied(v, function(x) x > 5)))[, 1], rep(NA_real_, 3))

  e <- make_range_rast(c(0, NA, 0), c(100, NA, 100), ncol = 3, nrow = 1)
  expect_identical(terra::values(mask(fp, e))[, 1], c(7, NA, 9))
})

test_that("mask() passes terra's own arguments through", {
  field <- make_field()
  e <- make_range_rast(rep(0, 9), c(rep(100, 8), NA))
  inside <- mask(field, e)
  outside <- mask(field, e, inverse = TRUE)
  expect_identical(is.na(terra::values(outside)), !is.na(terra::values(inside)))
})

test_that("mask() on a 3D object by 2D input is still terra's, and plain rasters are untouched", {
  field <- make_field()
  fp <- make_footprint(c(1, NA, 1, 1, 1, 1, 1, 1, NA))

  by_raster <- mask(field, fp)
  expect_identical(class(by_raster)[[1]], "SpatVoxel")
  expect_true(all(is.na(terra::values(by_raster)[c(2, 9), ])))
  expect_false(anyNA(terra::values(by_raster)[-c(2, 9), ]))

  by_polygon <- mask(field, make_left_polygon())
  expect_identical(class(by_polygon)[[1]], "SpatVoxel")
  expect_true(all(is.na(terra::values(by_polygon)[c(3, 6, 9), ])))

  plain <- methods::as(field, "SpatRaster")
  expect_identical(terra::values(mask(plain, fp)),
               terra::values(terra::mask(plain, fp)))
  expect_identical(class(mask(plain, fp))[[1]], "SpatRaster")
})

test_that("mask() names the offending argument for a grid mismatch", {
  field <- make_field()
  e <- make_range_rast(rep(0, 4), rep(100, 4), ncol = 2, nrow = 2)
  expect_error(mask(field, e), "`x` and `mask` must be on the same grid")
})

test_that("the query verbs reject a tagged object that is no longer valid", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  # `[[` keeps the SpatEnvelope tag on a single layer.
  broken <- a[["depth_min"]]
  expect_s4_class(broken, "SpatEnvelope")

  expect_error(intersect_3d(a, broken), "`y` is tagged SpatEnvelope")
  expect_error(intersects_3d(broken, a), "`x` is tagged SpatEnvelope")
  expect_error(mask(make_footprint(rep(1, 9)), broken), "`mask` is tagged")
  expect_error(mask(make_field(), broken), "as_envelope\\(\\)")
})

test_that("constructors accept tagged input as a footprint or template", {
  field <- make_field()
  # A voxel layer used as a footprint is still tagged SpatVoxel.
  e <- as_envelope(field[[1]], depth_min = 0, depth_max = 100)
  expect_identical(class(e)[[1]], "SpatEnvelope")
  expect_identical(terra::values(e)[, "depth_min"], rep(0, 9))

  poly <- sf::st_as_sf(make_left_polygon())
  e2 <- vect_to_envelope(poly, field, depth_min = 0, depth_max = 100)
  expect_identical(class(e2)[[1]], "SpatEnvelope")
  expect_identical(terra::values(e2)[, "depth_max"], c(100, 100, NA, 100, 100, NA, 100, 100, NA))
})

