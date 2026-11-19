# make_grid(), make_range_rast(), make_voxel(), make_field() and
# make_left_polygon() are in helper-example_objects.R; make_grid(vals = ) is a
# single-layer footprint.
# lintr checks this file on its own and can't see them, so the helpers below
# that call them carry object_usage nolints.

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

test_that("intersect_3d() on envelopes is empty where domains do not meet", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))

  # Touching only: a shared edge is no shared water.
  touching <- make_range_rast(rep(100, 9), rep(200, 9))
  expect_true(all(is.na(terra::values(intersect_3d(a, touching)))))

  # Vertically disjoint.
  deeper <- make_range_rast(rep(300, 9), rep(400, 9))
  expect_true(all(is.na(terra::values(intersect_3d(a, deeper)))))

  # Spatially disjoint: A in the top row, B in the bottom row.
  top <- make_range_rast(c(0, 0, 0, rep(NA, 6)),
                         c(100, 100, 100, rep(NA, 6)))
  bottom <- make_range_rast(c(rep(NA, 6), 0, 0, 0),
                            c(rep(NA, 6), 100, 100, 100))
  expect_true(all(is.na(terra::values(intersect_3d(top, bottom)))))
})

test_that("intersect_3d() keeps only cells present in both, either way round", {
  # A is absent from cell 3
  a <- make_range_rast(c(0, 0, NA, 0, 0, 0, 0, 0, 0),
                       c(100, 100, NA, 100, 100, 100, 100, 100, 100))

  cases <- list(
    list(label = "B absent from cell 1",
         b = make_range_rast(c(NA, 50, 50, 50, 50, 50, 50, 50, 50),
                             c(NA, 200, 200, 200, 200, 200, 200, 200, 200)),
         min = c(NA, 50, NA, rep(50, 6)), max = c(NA, 100, NA, rep(100, 6))),
    list(label = "B shallower than A in cell 2",
         b = make_range_rast(rep(50, 9),
                             c(200, 60, 200, 200, 200, 200, 200, 200, 200)),
         min = c(50, 50, NA, rep(50, 6)), max = c(100, 60, NA, rep(100, 6)))
  )
  for (case in cases) {
    out <- intersect_3d(a, case$b)

    expect_identical(terra::values(out)[, "depth_min"], case$min,
                     label = case$label)
    expect_identical(terra::values(out)[, "depth_max"], case$max,
                     label = case$label)
    # symmetric
    expect_identical(terra::values(intersect_3d(case$b, a)),
                     terra::values(out), label = case$label)
  }
})

# ---- intersect_3d: voxels and mixed input -----------------------------------

test_that("intersect_3d() on voxels is the co-occupied set of levels", {
  depths <- c(0, 100, 200)
  # Cell 1: occupied at all levels in A, at 0 and 200 in B (interior gap).
  # Cell 2: A only. Cell 3: neither.
  a <- make_voxel(cbind(c(1, 1, NA), c(1, 1, NA), c(1, 1, NA)), depths,
                  ncol = 3, nrow = 1)
  b <- make_voxel(cbind(c(1, NA, NA), c(NA, NA, NA), c(1, NA, NA)), depths,
                  ncol = 3, nrow = 1)

  out <- intersect_3d(a, b)

  expect_s4_class(out, "SpatVoxel")
  expect_identical(class(out)[[1]], "SpatVoxel")
  expect_named(out, paste0("presence_depth=", depths))
  expect_identical(depths(out), depths)
  expect_identical(unname(terra::values(out)),
               unname(cbind(c(1, NA, NA), c(NA, NA, NA), c(1, NA, NA))))
})

test_that("intersect_3d() on mixed input promotes the envelope to a voxel", {
  # the envelope is promoted to the voxel's levels, as done by hand below.
  depths <- c(0, 50, 100, 150, 200)
  e <- make_range_rast(rep(0, 9), rep(100, 9))
  v <- make_voxel(matrix(1, nrow = 9, ncol = 5), depths)

  by_hand <- intersect_3d(envelope_to_voxel(e, depths), v)
  expect_identical(terra::values(intersect_3d(e, v)), terra::values(by_hand))
  expect_identical(terra::values(intersect_3d(v, e)), terra::values(by_hand))
  expect_named(intersect_3d(e, v), paste0("presence_depth=", depths))

  # `bounds` reaches the promotion.
  mid <- intersect_3d(envelope_to_voxel(e, depths, bounds = "midpoint"), v)
  expect_identical(terra::values(intersect_3d(e, v, bounds = "midpoint")),
               terra::values(mid))
})

# ---- intersect_3d: 2D input -------------------------------------------------

test_that("intersect_3d() restricts an envelope to a footprint raster", {
  # the footprint restricts cells horizontally; depths are left as they are.
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  fp <- make_grid(vals = c(1, NA, 0, 7, NA, NA, 1, 1, 1))

  out <- intersect_3d(a, fp)
  expect_identical(class(out)[[1]], "SpatEnvelope")
  expect_identical(terra::values(out)[, "depth_min"],
                   c(0, NA, 0, 0, NA, NA, 0, 0, 0))
  expect_identical(terra::values(out)[, "depth_max"],
                   c(100, NA, 100, 100, NA, NA, 100, 100, 100))

  # Mirrored order gives the same object.
  expect_identical(terra::values(intersect_3d(fp, a)), terra::values(out))
})

test_that("intersect_3d() restricts a voxel to a footprint raster", {
  # the footprint restricts cells horizontally, at every depth level.
  depths <- c(0, 100)
  v <- make_voxel(cbind(c(3, NA, 3), c(NA, 3, 3)), depths,
                  ncol = 3, nrow = 1, varname = "t")
  fp <- make_grid(vals = c(1, 1, NA), ncol = 3, nrow = 1)

  out <- intersect_3d(v, fp)
  expect_identical(class(out)[[1]], "SpatVoxel")
  expect_named(out, paste0("presence_depth=", depths))
  # Presence, not the field values, and the third cell is gone at every level.
  expect_identical(unname(terra::values(out)),
               unname(cbind(c(1, NA, NA), c(NA, 1, NA))))
  expect_identical(terra::values(intersect_3d(fp, v)), terra::values(out))
})

test_that("intersect_3d() accepts polygons as SpatVector, sf and sfc", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  poly <- make_left_polygon()
  poly_sf <- sf::st_as_sf(poly)
  expected_min <- c(0, 0, NA, 0, 0, NA, 0, 0, NA)

  cases <- list(
    "SpatVector" = list(a, poly),
    "SpatVector first" = list(poly, a),
    "sf" = list(a, poly_sf),
    "sfc" = list(a, sf::st_geometry(poly_sf)),
    # polygons in another CRS are projected onto the grid
    "SpatVector in EPSG:4326" = list(a, terra::project(poly, "EPSG:4326"))
  )
  for (form in names(cases)) {
    out <- do.call(intersect_3d, cases[[form]])
    expect_identical(terra::values(out)[, "depth_min"], expected_min,
                     label = form)
  }
})

# ---- intersect_3d: composition with volume ----------------------------------

test_that("volume(intersect_3d(a, b)) equals summed calc_volume_overlap()", {
  depths <- c(0, 50, 100, 150, 200)
  e1 <- make_range_rast(c(0, 0, NA, 0, 0, 0, 0, 0, 0),
                        c(100, 100, NA, 100, 100, 100, 120, 100, 100))
  e2 <- make_range_rast(rep(50, 9),
                        c(200, 60, 200, 200, 200, 200, 200, 200, 200))
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
# nolint start: object_usage_linter.
truth_envelopes <- function() {
  list(
    a = make_range_rast(c(0, 0, 0, NA), c(100, 100, 100, NA),
                        ncol = 4, nrow = 1),
    b = make_range_rast(c(50, 300, NA, NA), c(200, 400, NA, NA),
                        ncol = 4, nrow = 1)
  )
}
# nolint end

test_that("intersects_3d() on envelopes gives the tri-state truth table", {
  e <- truth_envelopes()

  cases <- list(
    list(label = "truth table", a = e$a, b = e$b,
         expected = c(TRUE, FALSE, FALSE, NA)),
    # touching intervals share an edge and no water
    list(label = "touching intervals",
         a = make_range_rast(rep(0, 9), rep(100, 9)),
         b = make_range_rast(rep(100, 9), rep(200, 9)),
         expected = rep(FALSE, 9))
  )
  for (case in cases) {
    out <- intersects_3d(case$a, case$b)

    expect_identical(class(out)[[1]], "SpatRaster", label = case$label)
    expect_true(terra::is.bool(out), label = case$label)
    expect_named(out, "intersects", label = case$label)
    expect_identical(terra::values(out)[, 1], case$expected,
                     label = case$label)
    expect_identical(terra::values(intersects_3d(case$b, case$a))[, 1],
                     case$expected, label = case$label)
  }
})

test_that("adjacent depth ranges do not intersect as voxels either", {
  # The envelope rule, touching is not overlapping, must survive promotion to
  # depth levels: [0, 100] and [100, 200] share only the boundary at 100 m,
  # so their voxels occupy no common level and every 3D verb agrees with the
  # envelope answer, whichever combination of types it is asked about.
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(100, 9), rep(200, 9))
  depths <- c(0, 100, 200)
  va <- envelope_to_voxel(a, depths, bounds = "top")
  vb <- envelope_to_voxel(b, depths, bounds = "top")

  # voxel vs voxel
  expect_identical(terra::values(intersects_3d(va, vb))[, 1], rep(FALSE, 9))
  expect_identical(terra::values(intersects_3d(vb, va))[, 1], rep(FALSE, 9))
  expect_true(all(is.na(terra::values(intersect_3d(va, vb)))))

  # envelope vs voxel, both orders
  expect_identical(terra::values(intersects_3d(a, vb, bounds = "top"))[, 1],
                   rep(FALSE, 9))
  expect_identical(terra::values(intersects_3d(vb, a, bounds = "top"))[, 1],
                   rep(FALSE, 9))
  expect_true(all(is.na(terra::values(intersect_3d(a, vb, bounds = "top")))))
  expect_true(all(is.na(terra::values(intersect_3d(va, b, bounds = "top")))))

  # and the same under the midpoint convention, with levels whose slab edge
  # falls at 100 m (halfway between 50 and 150)
  depths_mid <- c(0, 50, 150, 200)
  va_mid <- envelope_to_voxel(a, depths_mid, bounds = "midpoint")
  vb_mid <- envelope_to_voxel(b, depths_mid, bounds = "midpoint")
  expect_identical(terra::values(intersects_3d(va_mid, vb_mid))[, 1],
                   rep(FALSE, 9))
  expect_identical(
    terra::values(intersects_3d(a, vb_mid, bounds = "midpoint"))[, 1],
    rep(FALSE, 9)
  )

  # overlapping by any amount is still an intersection, so this is the
  # boundary itself and not a wider gap
  c_ <- make_range_rast(rep(99, 9), rep(200, 9))
  expect_identical(
    terra::values(intersects_3d(
      va, envelope_to_voxel(c_, depths, bounds = "top")
    ))[, 1],
    rep(TRUE, 9)
  )
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

test_that("intersects_3d() on mixed input promotes the envelope to a voxel", {
  # the envelope is promoted to the voxel's levels, as done by hand below.
  depths <- c(0, 50, 100, 150, 200, 300, 400)
  e <- truth_envelopes()$a
  v <- envelope_to_voxel(truth_envelopes()$b, depths)

  by_hand <- intersects_3d(envelope_to_voxel(e, depths), v)
  expect_identical(terra::values(intersects_3d(e, v))[, 1],
                   terra::values(by_hand)[, 1])
  expect_identical(terra::values(intersects_3d(v, e))[, 1],
                   terra::values(by_hand)[, 1])
  expect_identical(terra::values(by_hand)[, 1], c(TRUE, FALSE, FALSE, NA))
})

test_that("intersects_3d() is TRUE exactly where intersect_3d() is non-empty", {
  e <- truth_envelopes()
  shared <- intersect_3d(e$a, e$b)
  # nolint start: scalar_in_linter
  # use %in% instead of ==, because NA == TRUE returns NA, compared to
  # NA %in% TRUE returns FALSE
  expect_identical(terra::values(intersects_3d(e$a, e$b))[, 1] %in% TRUE,
               !is.na(terra::values(shared)[, "depth_min"]))
  # nolint end

  depths <- c(0, 50, 100, 150, 200, 300, 400)
  va <- envelope_to_voxel(e$a, depths)
  vb <- envelope_to_voxel(e$b, depths)
  shared_v <- intersect_3d(va, vb)
  # nolint start: scalar_in_linter
  # use %in% instead of ==, because NA == TRUE returns NA, compared to
  # NA %in% TRUE returns FALSE
  expect_identical(terra::values(intersects_3d(va, vb))[, 1] %in% TRUE,
               rowSums(!is.na(terra::values(shared_v))) > 0)
  # nolint end
})

test_that("intersects_3d() with 2D input tests presence on both sides", {
  a <- make_range_rast(c(0, 0, NA, NA), c(100, 100, NA, NA),
                       ncol = 4, nrow = 1)
  fp <- make_grid(vals = c(1, NA, 1, NA), ncol = 4, nrow = 1)

  expect_identical(terra::values(intersects_3d(a, fp))[, 1],
                   c(TRUE, FALSE, FALSE, NA))
  expect_identical(terra::values(intersects_3d(fp, a))[, 1],
                   c(TRUE, FALSE, FALSE, NA))

  v <- envelope_to_voxel(a, c(0, 100))
  expect_identical(terra::values(intersects_3d(v, fp))[, 1],
                   c(TRUE, FALSE, FALSE, NA))
  expect_identical(terra::values(intersects_3d(fp, v))[, 1],
                   c(TRUE, FALSE, FALSE, NA))

  # Polygons: the left two columns of the 3x3 grid.
  a9 <- make_range_rast(c(0, 0, 0, NA, NA, NA, 0, 0, 0),
                        c(100, 100, 100, NA, NA, NA, 100, 100, 100))
  expect_identical(terra::values(intersects_3d(a9, make_left_polygon()))[, 1],
               c(TRUE, TRUE, FALSE, FALSE, FALSE, NA, TRUE, TRUE, FALSE))
})

# ---- mask -------------------------------------------------------------------

test_that("mask(voxel, envelope) equals the hand-written idiom exactly", {
  field <- make_field()
  e <- make_range_rast(c(0, 0, NA, 60, 60, 60, 0, 0, 0),
                       c(100, 30, NA, 120, 500, 120, 200, 200, 200))

  by_hand <- terra::mask(field, envelope_to_voxel(e, depths(field)))
  out <- mask(field, e)

  expect_identical(class(out)[[1]], "SpatVoxel")
  expect_named(out, names(field))
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

test_that("mask(envelope, envelope) keeps whole intervals where domains meet", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(c(50, 300, NA, 50, 50, 50, 50, 50, 50),
                       c(200, 400, NA, 200, 200, 200, 200, 200, 200))

  out <- mask(a, b)
  expect_identical(class(out)[[1]], "SpatEnvelope")
  expect_identical(terra::values(out)[, "depth_min"], c(0, NA, NA, rep(0, 6)))
  expect_identical(terra::values(out)[, "depth_max"],
                   c(100, NA, NA, rep(100, 6)))

  # Before: two envelopes 100 m apart came back unchanged.
  deeper <- make_range_rast(rep(300, 9), rep(400, 9))
  expect_true(all(is.na(terra::values(mask(a, deeper)))))
})

test_that("mask(envelope, voxel) keeps cells where voxel meets envelope", {
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
  fp <- make_grid(vals = c(7, 8, 9), ncol = 3, nrow = 1)

  out <- mask(fp, v)
  expect_identical(class(out)[[1]], "SpatRaster")
  expect_identical(terra::values(out)[, 1], c(7, 8, NA))
  expect_identical(
    terra::values(mask(fp, occupied(v, function(x) x > 5)))[, 1],
    rep(NA_real_, 3)
  )

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

test_that("mask() by 2D input applies terra's mask to 3D and plain rasters", {
  # a 2D mask needs no depth handling: a 3D object keeps its class, and a
  # plain raster gets exactly what terra returns.
  field <- make_field()
  fp <- make_grid(vals = c(1, NA, 1, 1, 1, 1, 1, 1, NA))

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

# ---- input errors -----------------------------------------------------------

test_that("the 3D verbs reject input they cannot compare", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  # an envelope on another grid
  a_2x2 <- make_range_rast(rep(0, 4), rep(100, 4), ncol = 2, nrow = 2)
  fp <- make_grid(vals = rep(1, 9))
  small <- make_grid(vals = rep(1, 4), ncol = 2, nrow = 2)
  t2 <- make_voxel(cbind(c(5, 5, 5), c(5, 5, 5)), c(0, 100),
                   ncol = 3, nrow = 1, varname = "t")
  p3 <- make_voxel(cbind(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1)), c(0, 100, 200),
                   ncol = 3, nrow = 1)
  field <- make_field()

  cases <- list(
    # voxel-only arguments for envelopes
    list(call = function() intersect_3d(a, a, bounds = "top"),
         pattern = "unused argument"),
    list(call = function() intersects_3d(a, a, bounds = "top"),
         pattern = "unused argument"),
    # voxels on different depth levels
    list(call = function() intersect_3d(t2, p3), pattern = "same depth levels"),
    # inputs on different grids
    list(call = function() intersect_3d(a, a_2x2), pattern = "same grid"),
    list(call = function() intersects_3d(a, a_2x2), pattern = "same grid"),
    list(call = function() intersect_3d(a, small), pattern = "same grid"),
    # mask() names the offending argument
    list(call = function() mask(field, a_2x2),
         pattern = "`x` and `mask` must be on the same grid"),
    # input that cannot be read as a footprint
    list(call = function() intersect_3d(fp, fp),
         pattern = "neither `x` nor `y`"),
    list(call = function() intersects_3d(fp, fp),
         pattern = "neither `x` nor `y`"),
    list(call = function() intersect_3d(a, "hello"),
         pattern = "SpatVector, sf, sfc"),
    list(call = function() intersect_3d(a, c(fp, fp)),
         pattern = "single-layer footprint")
  )
  for (case in cases) {
    expect_error(case$call(), case$pattern,
                 label = paste(deparse(body(case$call)), collapse = " "))
  }
})

test_that("the query verbs reject a tagged object that is no longer valid", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  # `[[` keeps the SpatEnvelope tag on a single layer.
  broken <- a[["depth_min"]]
  expect_s4_class(broken, "SpatEnvelope")

  expect_error(intersect_3d(a, broken), "`y` is tagged SpatEnvelope")
  expect_error(intersects_3d(broken, a), "`x` is tagged SpatEnvelope")
  expect_error(mask(make_grid(vals = rep(1, 9)), broken), "`mask` is tagged")
  expect_error(mask(make_field(), broken), "as_envelope\\(\\)")
})

# ---- terra::intersect() guard -----------------------------------------------

test_that("terra::intersect() on a 3D object errors and names replacements", {
  # Before: intersect(a, b) on two envelopes 100 m apart returned a tagged
  # SpatEnvelope of TRUE/TRUE that passed validObject().
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(200, 9), rep(300, 9))
  v <- envelope_to_voxel(a, c(0, 50, 100))
  fp <- make_grid(vals = rep(1, 9))

  expect_error(terra::intersect(a, b), "intersects_3d\\(\\)")
  expect_error(terra::intersect(v, v), "intersect_3d\\(\\)")
  expect_error(terra::intersect(a, v), "ignores the depth axis")
  expect_error(terra::intersect(v, a), "ignores the depth axis")
  expect_error(terra::intersect(a, fp), "ignores the depth axis")
  expect_error(terra::intersect(fp, v), "ignores the depth axis")
})

test_that("terra's own intersect() methods are untouched", {
  fp1 <- make_grid(vals = c(1, 2, NA, NA), ncol = 4, nrow = 1)
  fp2 <- make_grid(vals = c(9, NA, 7, NA), ncol = 4, nrow = 1)
  expect_identical(terra::values(terra::intersect(fp1, fp2))[, 1],
                   c(TRUE, FALSE, FALSE, NA))

  a <- make_range_rast(rep(0, 9), rep(100, 9))
  ex <- terra::intersect(a, terra::ext(0, 1000, 0, 1000))
  expect_s4_class(ex, "SpatExtent")

  poly <- make_left_polygon()
  expect_s4_class(terra::intersect(poly, poly), "SpatVector")

  expect_identical(intersect(1:3, 2:5), 2:3)
})

# ---- depth edges (#138) -----------------------------------------------------

test_that("voxel results keep the depth edges of their inputs", {
  depths <- c(0.5, 1.5, 3)
  edges <- cbind(c(0, 1.4, 2), c(1.4, 2, 4))
  field <- as_voxel(make_voxel(seq_len(27), depths), depth_bounds = edges)
  plain <- make_voxel(1, depths)

  expect_equal(edges_of(intersect_3d(field, plain)), edges)
  expect_equal(edges_of(intersect_3d(plain, field)), edges)
  expect_equal(edges_of(mask(field, plain)), edges)
  # Masking by a voxel that carries the edges keeps them too, so the result
  # stands for the same water as the intersection.
  expect_equal(edges_of(mask(plain, field)), edges)
  expect_equal(volume(mask(plain, field)), volume(intersect_3d(plain, field)))

  # A 2D footprint keeps the voxel's edges, and so its volume.
  fp <- terra::rast(make_grid(), vals = 1)
  expect_equal(edges_of(intersect_3d(field, fp)), edges)
  expect_equal(volume(intersect_3d(field, fp)), volume(field))

  # Masking by an envelope uses the edges: [1.1, 1.3] is in the 0.5 m
  # level's slab (0-1.4 m), not the 1.5 m level's as under "midpoint".
  e <- make_range_rast(rep(1.1, 9), rep(1.3, 9))
  masked <- mask(field, e)
  expect_equal(edges_of(masked), edges)
  expect_equal(terra::global(!is.na(masked), "sum")[, 1], c(9, 0, 0))
  expect_equal(
    terra::global(!is.na(mask(plain, e, bounds = "midpoint")), "sum")[, 1],
    c(0, 9, 0)
  )
})

test_that("voxel pairs with different depth edges stop", {
  depths <- c(0.5, 1.5, 3)
  a <- as_voxel(make_voxel(1, depths), depth_bounds = cbind(c(0, 1, 2),
                                                            c(1, 2, 4)))
  b <- as_voxel(make_voxel(1, depths), depth_bounds = cbind(c(0, 1, 2),
                                                            c(1, 2, 5)))
  for (call in list(quote(intersect_3d(a, b)), quote(intersects_3d(a, b)),
                    quote(mask(a, b)))) {
    err <- expect_error(eval(call), class = "ocean3d_error_depth_edges_differ")
    expect_match(conditionMessage(err), "different depth edges")
  }
})
