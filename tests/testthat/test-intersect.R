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
  other_crs <- make_grid(vals = rep(1, 9))
  terra::crs(other_crs) <- "EPSG:3857"
  mismatch <- "ocean3d_error_grid_mismatch"
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
    list(call = function() intersect_3d(a, a_2x2),
         pattern = "`y` does not align with `x`", class = mismatch),
    list(call = function() intersects_3d(a, a_2x2),
         pattern = "`y` does not align with `x`", class = mismatch),
    list(call = function() intersect_3d(a, small),
         pattern = "`y` does not align with `x`", class = mismatch),
    list(call = function() intersect_3d(a, other_crs),
         pattern = "`y` has a different CRS than `x`", class = mismatch),
    # a 2D object given first is named as the caller passed it
    list(call = function() intersect_3d(small, a),
         pattern = "`x` does not align with `y`", class = mismatch),
    # mask() names the offending argument
    list(call = function() mask(field, a_2x2),
         pattern = "`mask` does not align with `x`", class = mismatch),
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
    expect_error(case$call(), case$pattern, class = case$class,
                 label = paste(deparse(body(case$call)), collapse = " "))
  }
})

test_that("a bad 2D input is named as the caller passed it", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  fp <- make_grid(vals = rep(1, 9))
  for (query in list(intersect_3d, intersects_3d)) {
    expect_error(query(a, c(fp, fp)), "`y` must be a single-layer footprint",
                 class = "ocean3d_error_footprint_layers")
    expect_error(query(c(fp, fp), a), "`x` must be a single-layer footprint",
                 class = "ocean3d_error_footprint_layers")
    expect_error(query("hello", a), "`x` must be a SpatEnvelope.*got character",
                 class = "ocean3d_error_footprint_type")
  }
  # The 3D side is checked before the pair is swapped.
  broken <- a[["depth_min"]]
  expect_error(intersects_3d(fp, broken), "`y` is tagged SpatEnvelope")
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
  # An envelope placed on the voxel takes its edges.
  e_full <- make_range_rast(rep(0, 9), rep(4, 9))
  expect_equal(edges_of(intersect_3d(e_full, field)), edges)
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

# SpatVector3D polygons (#115) ----------------------------------------------

test_that("intersect_3d() of polygons and an envelope compares depths", {
  for (result in list(intersect_3d(polygons_0_100(), envelope_50_200()),
                      intersect_3d(envelope_50_200(), polygons_0_100()))) {
    expect_s4_class(result, "SpatEnvelope")
    vals <- terra::values(result)
    expect_identical(sum(inside()), 9L)
    expect_true(all(vals[inside(), "depth_min"] == 50))
    expect_true(all(vals[inside(), "depth_max"] == 100))
    expect_true(all(is.na(vals[!inside(), ])))
  }
  # The overlap is 9 cells x 50 m, half the polygons' 100 m. The tolerance
  # allows for the polygons' geodesic edges (see test-volume.R).
  overlap <- volume(intersect_3d(polygons_0_100(), envelope_50_200()))
  expect_equal(overlap, volume(polygons_0_100()) / 2, tolerance = 1e-3)
})

test_that("intersect_3d() of polygons and a voxel uses the voxel's levels", {
  voxel <- envelope_to_voxel(envelope_50_200(), depths = c(0, 50, 100, 150))
  polygons <- polygons_0_100()
  on_grid <- vector3d_to_envelope(polygons, voxel)
  expect_identical(
    terra::values(intersect_3d(polygons, voxel, bounds = "top")),
    terra::values(intersect_3d(on_grid, voxel, bounds = "top"))
  )
  expect_identical(terra::values(intersect_3d(voxel, polygons)),
                   terra::values(intersect_3d(on_grid, voxel)))
  expect_s4_class(intersect_3d(voxel, polygons), "SpatVoxel")
})

test_that("queries with a voxel keep the depth gap between polygons", {
  # Two bands over the polygon's cells, 0-50 m and 200-300 m, on levels
  # whose midpoint slabs are 0-50, 50-150, 150-250 and 250-300 m: the 100 m
  # level lies between them.
  bands <- rbind(make_polygon(), make_polygon())
  bands$top <- c(0, 200)
  bands$bottom <- c(50, 300)
  bands <- as_vector3d(bands, depth_min = "top", depth_max = "bottom")
  layers <- lapply(1:4, function(i) make_template(i))
  voxel <- as_voxel(terra::rast(layers), depths = c(0, 100, 200, 300),
                    varname = "temp")

  on_levels <- vector3d_to_voxel(bands, voxel)
  for (shared in list(intersect_3d(bands, voxel),
                      intersect_3d(voxel, bands))) {
    expect_identical(terra::values(shared), terra::values(on_levels))
  }
  masked <- terra::values(mask(voxel, bands))
  expect_true(all(masked[inside(), c(1, 3, 4)] == rep(c(1, 3, 4), each = 9)))
  expect_true(all(is.na(masked[inside(), 2])))
  # The envelope route fills the gap, so a voxel masked by the envelope keeps
  # the 100 m level.
  filled <- terra::values(mask(voxel, vector3d_to_envelope(bands, voxel)))
  expect_true(all(filled[inside(), 2] == 2))

  # `bounds` reaches the placement: under "top" the slabs are 0-100,
  # 100-200, 200-300 and nothing, so the 300 m level is empty and the 100 m
  # level is covered by neither band.
  top <- terra::values(intersect_3d(bands, voxel, bounds = "top"))
  expect_true(all(top[inside(), c(1, 3)] == 1))
  expect_true(all(is.na(top[inside(), c(2, 4)])))
  expect_identical(terra::values(intersects_3d(bands, voxel, bounds = "top")),
                   terra::values(intersects_3d(voxel, bands, bounds = "top")))
  expect_true(all(terra::values(mask(voxel, bands, bounds = "top"))[
    inside(), c(2, 4)
  ] |> is.na()))

  # Without a voxel there are no levels to place, so `bounds` is an error,
  # as for two envelopes.
  for (x in list(envelope_50_200(), make_template(1))) {
    expect_error(mask(x, bands, bounds = "top"), "unused argument",
                 class = "ocean3d_error_unused_argument")
  }
})

test_that("intersects_3d() of polygons tests depth overlap per cell", {
  deep <- as_envelope(make_template(1), depth_min = 300, depth_max = 400)
  cases <- list(
    list(label = "overlapping depths", y = envelope_50_200(), inside = TRUE),
    list(label = "disjoint depths", y = deep, inside = FALSE)
  )
  for (case in cases) {
    for (answer in list(intersects_3d(polygons_0_100(), case$y),
                        intersects_3d(case$y, polygons_0_100()))) {
      vals <- terra::values(answer)[, 1]
      expect_true(all(vals[inside()] == case$inside), label = case$label)
      # Outside the polygons only the envelope is present.
      expect_false(any(vals[!inside()]), label = case$label)
    }
  }
})

test_that("polygons that miss the grid share no space", {
  far <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(20, 20), c(21, 20), c(21, 21), c(20, 21), c(20, 20)
  ))), crs = 4326))
  far <- as_vector3d(far, depth_min = 0, depth_max = 100)
  shared <- intersect_3d(far, envelope_50_200())
  expect_s4_class(shared, "SpatEnvelope")
  expect_true(all(is.na(terra::values(shared))))
  expect_false(any(terra::values(intersects_3d(far, envelope_50_200()))))
})

test_that("mask() by polygons keeps values where they reach", {
  # Levels 0, 100 and 200 m stand for 0-50, 50-150 and 150-200 m, so the
  # polygons' 0-100 m reaches the first two.
  layers <- lapply(1:3, function(i) make_template(i))
  voxel <- as_voxel(terra::rast(layers), depths = c(0, 100, 200),
                    varname = "temp")
  masked <- mask(voxel, polygons_0_100())
  expect_s4_class(masked, "SpatVoxel")
  vals <- terra::values(masked)
  expect_true(all(vals[inside(), 1:2] == rep(1:2, each = 9)))
  expect_true(all(is.na(vals[inside(), 3])))
  expect_true(all(is.na(vals[!inside(), ])))

  # An envelope is kept whole where the depths overlap; a plain raster where
  # the polygons are.
  expect_identical(
    terra::values(mask(envelope_50_200(), polygons_0_100())),
    terra::values(mask(envelope_50_200(),
                       vector3d_to_envelope(polygons_0_100(), make_template())))
  )
  plain <- terra::values(mask(make_template(7), polygons_0_100()))[, 1]
  expect_identical(!is.na(plain), inside())
})

test_that("an empty polygon set shares no space with a grid", {
  far <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(20, 20), c(21, 20), c(21, 21), c(20, 21), c(20, 20)
  ))), crs = 4326))
  far <- as_vector3d(far, depth_min = 0, depth_max = 100)
  # intersect_3d() of two polygon sets that share nothing has no polygons.
  none <- intersect_3d(polygons_0_100(), far)
  expect_identical(length(none), 0L)

  shared <- intersect_3d(none, envelope_50_200())
  expect_s4_class(shared, "SpatEnvelope")
  expect_true(all(is.na(terra::values(shared))))
  expect_false(any(terra::values(intersects_3d(envelope_50_200(), none))))
  masked <- mask(envelope_50_200(), none)
  expect_true(all(is.na(terra::values(masked))))
})

test_that("comparing polygons with a grid drops their times with a message", {
  dated <- sf::st_sf(day = as.Date("2020-01-01"),
                     geometry = sf::st_geometry(make_polygon()))
  dated <- as_vector3d(dated, time = "day", depth_min = 0, depth_max = 100)
  expect_message(intersect_3d(dated, envelope_50_200()),
                 class = "ocean3d_message_time_dropped")
})

test_that("lines aren't compared with a grid or polygons", {
  line <- as_vector3d(data.frame(id = "a", lon = c(1, 2), lat = c(1, 2),
                                 depth = c(1, 2)), type = "lines", id = "id")
  voxel <- envelope_to_voxel(envelope_50_200(), depths = c(0, 100))
  for (call in list(function() intersect_3d(line, envelope_50_200()),
                    function() intersects_3d(envelope_50_200(), line),
                    function() mask(envelope_50_200(), line),
                    function() intersect_3d(voxel, line),
                    function() mask(voxel, line),
                    function() mask(make_template(1), line),
                    function() intersect_3d(line, polygons_0_100()),
                    function() intersects_3d(polygons_0_100(), line))) {
    expect_error(call(), "This needs polygons or points, but these are lines",
                 class = "ocean3d_error_vector3d_type")
    expect_error(call(), "buffer_3d() turns lines into polygons",
                 class = "ocean3d_error_vector3d_type", fixed = TRUE)
  }
})

# SpatVector3D points (#115) ------------------------------------------------

test_that("a point covers its cell and the level nearest its depth", {
  voxel <- voxel_1_2_3()
  # Cell 13 (the middle) at 40 m: nearest level 0 m. Cell 23 at 160 m:
  # nearest level 200 m. The third point is off the grid.
  points <- points_at(lon = c(2.5, 2.5, 20), lat = c(2.5, 0.5, 2.5),
                      depth = c(40, 160, 10))
  for (shared in list(intersect_3d(points, voxel),
                      intersect_3d(voxel, points))) {
    expect_s4_class(shared, "SpatVoxel")
    expect_identical(depths(shared), c(0, 100, 200))
    vals <- terra::values(shared)
    expect_identical(which(!is.na(vals[, 1])), 13L)
    expect_identical(which(!is.na(vals[, 2])), integer(0))
    expect_identical(which(!is.na(vals[, 3])), 23L)
    expect_true(all(vals[c(13, 23), c(1, 3)][c(1, 4)] == 1))
  }
  # A point halfway between two levels takes the shallower one, and the
  # result carries the voxel's depth edges.
  edged <- as_voxel(voxel, depth_bounds = cbind(c(0, 50, 150),
                                                c(50, 150, 250)))
  halfway <- intersect_3d(points_at(2.5, 2.5, 50), edged)
  expect_identical(which(!is.na(terra::values(halfway)[, 1])), 13L)
  expect_true(all(is.na(terra::values(halfway)[, 2:3])))
  expect_identical(depth_bounds(halfway), depth_bounds(edged))

  # intersects_3d(): TRUE in the point's cell, FALSE where only the voxel is.
  answer <- terra::values(intersects_3d(points, voxel))[, 1]
  expect_identical(which(answer), c(13L, 23L))
  expect_false(any(answer[-c(13, 23)]))
  expect_identical(terra::values(intersects_3d(voxel, points)),
                   terra::values(intersects_3d(points, voxel)))
})

test_that("mask(voxel, points) keeps the value at each point's level only", {
  voxel <- voxel_1_2_3()
  points <- points_at(lon = c(2.5, 2.5), lat = c(2.5, 0.5),
                      depth = c(40, 160))
  masked <- mask(voxel, points)
  expect_s4_class(masked, "SpatVoxel")
  vals <- terra::values(masked)
  expect_equal(unname(vals[13, ]), c(1, NA, NA))
  expect_equal(unname(vals[23, ]), c(NA, NA, 3))
  expect_true(all(is.na(vals[-c(13, 23), ])))

  # With "nearest_level", the level comes from the point's depth, not from
  # where the voxel has values: a point at a level that is NA there masks to
  # nothing.
  gappy <- voxel
  gappy[[1]] <- terra::setValues(make_template(), NA)
  gappy <- as_voxel(gappy)
  expect_true(all(is.na(terra::values(
    mask(gappy, points_at(2.5, 2.5, 40), method = "nearest_level")
  ))))

  # Points take the nearest level, so `bounds` is an error, not ignored.
  for (call in list(function() intersect_3d(points, voxel, bounds = "top"),
                    function() mask(voxel, points, bounds = "top"))) {
    expect_error(call(), "points take a depth level of their own",
                 class = "ocean3d_error_unused_argument")
  }
  # A grid without a CRS can't place them.
  bare <- voxel
  terra::crs(bare) <- ""
  expect_error(mask(bare, points), "has no CRS",
               class = "ocean3d_error_envelope_input")

  # Points off the grid share nothing, and so do no points at all.
  far <- points_at(20, 20, 40)
  expect_true(all(is.na(terra::values(intersect_3d(far, voxel)))))
  expect_true(all(is.na(terra::values(mask(voxel, far)))))
  expect_false(any(terra::values(intersects_3d(far, voxel))))
})

test_that("mask(voxel, points) picks each point's level with `method`", {
  # One column of values at 0, 100, 200 and 300 m in cell 13: NA at 100 m
  # and 300 m, so the nearest level and the bottom read the values.
  vals <- matrix(NA_real_, 25, 4)
  vals[13, ] <- c(10, NA, 30, NA)
  voxel <- make_voxel(vals, c(0, 100, 200, 300), varname = "temp",
                      grid = make_template())
  kept <- function(depth, method) {
    out <- terra::values(mask(voxel, points_at(2.5, 2.5, depth),
                              method = method))
    out[!is.na(out)]
  }
  cases <- list(
    # The nearest level is NA here, so it keeps nothing.
    list(depth = 110, method = "nearest_level", want = numeric(0)),
    # The nearest level with a value: 200 m is nearer than 0 m.
    list(depth = 110, method = "nearest", want = 30),
    # Halfway between two levels with values, the shallower wins.
    list(depth = 100, method = "nearest", want = 10),
    # The deepest level with a value, whatever the point's depth.
    list(depth = 0, method = "bottom", want = 30),
    # The shallowest level, whatever the point's depth.
    list(depth = 290, method = "surface", want = 10)
  )
  for (case in cases) {
    expect_equal(kept(case$depth, case$method), case$want,
                 info = paste(case$method, case$depth))
  }
  # Only the chosen cell and level keep a value.
  out <- terra::values(mask(voxel, points_at(2.5, 2.5, 0), method = "bottom"))
  expect_equal(which(!is.na(out)), 13 + 25 * 2)
  # Without `method`, a point takes the nearest level with a value.
  expect_equal(terra::values(mask(voxel, points_at(2.5, 2.5, 110))),
               terra::values(mask(voxel, points_at(2.5, 2.5, 110),
                                  method = "nearest")))

  # A cell with no value at any level, and a point off the grid, keep nothing.
  for (method in c("nearest", "bottom", "surface")) {
    empty <- mask(voxel, points_at(c(0.5, 20), c(0.5, 20), c(0, 0)),
                  method = method)
    expect_true(all(is.na(terra::values(empty))), info = method)
  }

  # `method` applies only to points on a voxel.
  band <- as_vector3d(make_polygon(), depth_min = 0, depth_max = 100)
  envelope <- as_envelope(make_template(1), depth_min = 0, depth_max = 100)
  for (call in list(function() mask(voxel, band, method = "bottom"),
                    function() {
                      mask(envelope, points_at(2.5, 2.5, 50), method = "bottom")
                    })) {
    expect_error(call(), "applies only to a variable with depth",
                 class = "ocean3d_error_unused_argument")
  }
  expect_error(mask(voxel, points_at(2.5, 2.5, 0), method = "deepest"),
               "must be one of", class = "ocean3d_error_mask_method")

  # intersect_3d() and intersects_3d() read no values, so `method` stays an
  # unused argument there.
  point <- points_at(2.5, 2.5, 0)
  for (call in list(function() intersect_3d(point, voxel, method = "bottom"),
                    function() {
                      intersects_3d(point, voxel, method = "bottom")
                    })) {
    expect_error(call(), "mask\\(\\) of a SpatVoxel by points",
                 class = "ocean3d_error_unused_argument")
  }
})

test_that("many points pick their levels once per cell and depth", {
  vals <- matrix(NA_real_, 25, 4)
  vals[13, ] <- c(10, NA, 30, NA)
  voxel <- make_voxel(vals, c(0, 100, 200, 300), varname = "temp",
                      grid = make_template())
  # 1e5 points at three depths in two cells.
  n <- 1e5
  many <- points_at(lon = rep(c(2.5, 0.5), length.out = n), lat = 2.5,
                    depth = rep(c(0, 110, 290), length.out = n))
  calls <- 0L
  helper <- .summarise_depth_profile
  testthat::local_mocked_bindings(.summarise_depth_profile = function(...) {
    calls <<- calls + 1L
    helper(...)
  })
  for (method in c("nearest_level", "nearest", "bottom")) {
    calls <- 0L
    out <- terra::values(mask(voxel, many, method = method))
    # One call per distinct profile and depth, not one per point.
    expect_lte(calls, 6L)
    # The same levels as three single points in cell 13.
    one <- vapply(c(0, 110, 290), function(z) {
      v <- terra::values(mask(voxel, points_at(2.5, 2.5, z), method = method))
      sum((!is.na(v[13, ])) * seq_len(4))
    }, numeric(1))
    expect_equal(sort(unique(which(!is.na(out[13, ])))),
                 sort(unique(one[one > 0])), info = method)
  }
})

test_that("mask(voxel, points, method) reads a voxel from a netCDF file", {
  # Two cells, four levels, with empty levels in each column, as near the
  # seafloor. The expected values are the levels extract_to_point() picks
  # with the same names: the nearest level with a value, the deepest level
  # with a value, and the shallowest level.
  values <- array(NA_real_, c(2, 2, 4))
  values[1, 1, ] <- c(1, 2, NA, 4)
  values[2, 1, ] <- c(5, NA, 7, NA)
  path <- make_grid_file(lon = c(0.5, 1.5), lat = c(10.5, 11.5),
                         depth = c(0, 10, 20, 30), time = NULL,
                         values = values)
  on.exit(unlink(path), add = TRUE)
  voxel <- as_voxel(as_proxy(path))
  obs <- data.frame(lon = c(0.5, 0.5, 1.5, 1.5), lat = 10.5,
                    depth = c(12, 30, 0, 16))
  expected <- list(nearest = c(2, 4, 5, 7), bottom = c(4, 4, 7, 7),
                   surface = c(1, 1, 5, 5))
  for (method in names(expected)) {
    got <- vapply(seq_len(nrow(obs)), function(i) {
      v <- terra::values(mask(voxel, as_vector3d(obs[i, ]), method = method))
      v <- v[!is.na(v)]
      if (length(v)) v else NA_real_
    }, numeric(1))
    expect_equal(got, expected[[method]], info = method)
  }
})

test_that("points land on a grid in the other longitude convention", {
  east <- terra::rast(nrows = 2, ncols = 2, xmin = 178, xmax = 182,
                      ymin = 0, ymax = 2, crs = "EPSG:4326")
  voxel <- as_voxel(terra::rast(list(terra::setValues(east, 1:4),
                                     terra::setValues(east, 5:8))),
                    depths = c(0, 100), varname = "temp")
  # Longitude -179 is 181 on this grid: the second column.
  shared <- intersect_3d(points_at(c(-179, 178.5), c(1.5, 0.5), c(0, 100)),
                         voxel)
  expect_identical(terra::values(shared)[, 1], c(NA, 1, NA, NA))
  expect_identical(terra::values(shared)[, 2], c(NA, NA, 1, NA))
})

test_that("points and an envelope share the depths the points reach", {
  e <- envelope_50_200()
  # Two points in cell 13, at the envelope's top and inside it; one in cell
  # 21 below the envelope; one off the grid.
  points <- points_at(lon = c(2.5, 2.5, 0.5, 20), lat = c(2.5, 2.5, 0.5, 2.5),
                      depth = c(50, 120, 300, 10))
  for (shared in list(intersect_3d(points, e), intersect_3d(e, points))) {
    expect_s4_class(shared, "SpatEnvelope")
    vals <- terra::values(shared)
    expect_identical(vals[13, ], c(depth_min = 50, depth_max = 120))
    expect_true(all(is.na(vals[-13, ])))
  }
  # A lone point is a zero-thickness interval: points have no volume.
  lone <- intersect_3d(points_at(2.5, 2.5, 200), e)
  expect_identical(terra::values(lone)[13, ],
                   c(depth_min = 200, depth_max = 200))
  expect_identical(volume(lone), 0)
  # Just outside the interval shares nothing.
  expect_true(all(is.na(terra::values(intersect_3d(points_at(2.5, 2.5, 201),
                                                   e)))))
  # A point in a cell the envelope doesn't reach shares nothing.
  holed <- as_envelope(terra::setValues(make_template(), c(NA, rep(1, 24))),
                       depth_min = 50, depth_max = 200)
  expect_true(all(is.na(terra::values(intersect_3d(points_at(0.5, 4.5, 100),
                                                   holed)))))

  # intersects_3d(): TRUE where a point lies in the interval, FALSE where
  # the depths miss or only one side is present, NA where neither is.
  answer <- terra::values(intersects_3d(e, points))[, 1]
  expect_true(answer[13])
  expect_false(answer[21])
  expect_false(answer[1])
  holed_answer <- terra::values(intersects_3d(holed, points_at(0.5, 4.5,
                                                               100)))[, 1]
  expect_false(holed_answer[1])
  expect_true(all(holed_answer[-1] == FALSE))
})

test_that("mask() by points keeps cells the points reach", {
  e <- envelope_50_200()
  points <- points_at(lon = c(2.5, 0.5), lat = c(2.5, 0.5),
                      depth = c(120, 300))
  masked <- mask(e, points)
  expect_s4_class(masked, "SpatEnvelope")
  vals <- terra::values(masked)
  # The envelope's own interval is kept whole where a point lies in it.
  expect_identical(vals[13, ], c(depth_min = 50, depth_max = 200))
  expect_true(all(is.na(vals[-13, ])))
  # A plain raster keeps the cells that hold a point, whatever their depth.
  plain <- terra::values(mask(make_template(7), points))[, 1]
  expect_identical(which(!is.na(plain)), c(13L, 21L))
  expect_identical(unname(plain[c(13, 21)]), c(7, 7))
  # terra's own arguments pass through.
  outside <- terra::values(mask(make_template(7), points, inverse = TRUE))
  expect_true(all(is.na(outside[c(13, 21), 1])))
  expect_identical(sum(!is.na(outside)), 23L)
})

test_that("points inside polygons are matched by outline, depth and time", {
  polygons <- polygons_0_100()
  # In the outline at the bottom depth (inside), in the outline but too
  # deep, and outside the outline.
  points <- as_vector3d(data.frame(
    id = c("p1", "p2", "p3"), lon = c(2.5, 2.5, 0.5), lat = c(2.5, 2.5, 0.5),
    depth = c(100, 120, 10), tag = c("a", "b", "c")
  ), id = "id")
  shared <- intersect_3d(points, polygons)
  expect_s4_class(shared, "SpatVector3D")
  expect_identical(geomtype(shared), "points")
  f <- as.data.frame(shared)
  expect_identical(f$feature_id, "1")
  expect_identical(f$feature_id_x, "p1")
  expect_identical(f$feature_id_y, "1")
  expect_identical(f$depth, 100)
  expect_null(f$tag)
  # The other way round names the pair the other way round.
  swapped <- as.data.frame(intersect_3d(polygons, points))
  expect_identical(swapped$feature_id_x, "1")
  expect_identical(swapped$feature_id_y, "p1")
  expect_identical(swapped[c("lon", "lat", "depth")],
                   f[c("lon", "lat", "depth")])

  # intersects_3d() lists the same pairs from each side.
  expect_identical(intersects_3d(points, polygons),
                   list(p1 = "1", p2 = character(0), p3 = character(0)))
  expect_identical(intersects_3d(polygons, points), list(`1` = "p1"))

  # A point in two polygons gives one point per pair; none inside gives
  # no points.
  two <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(
    sf::st_geometry(make_polygon())[[1]], sf::st_geometry(make_polygon())[[1]],
    crs = 4326
  )), depth_min = 0, depth_max = 100)
  pairs <- as.data.frame(intersect_3d(points, two))
  expect_identical(pairs$feature_id, c("1", "2"))
  expect_identical(pairs$feature_id_y, c("1", "2"))
  expect_identical(length(intersect_3d(points_at(0.5, 0.5, 10), polygons)),
                   0L)
})

test_that("points match polygon times with both ends included", {
  dated <- sf::st_sf(start = as.Date("2020-01-01"),
                     end = as.Date("2020-01-15"),
                     geometry = sf::st_geometry(make_polygon()))
  dated <- as_vector3d(dated, time = c("start", "end"), depth_min = 0,
                       depth_max = 100)
  # On the first day, on the last day (dates include it), the day after,
  # and without a time.
  points <- as_vector3d(data.frame(
    id = c("first", "last", "after", "untimed"), lon = 2.5, lat = 2.5,
    depth = 50,
    day = as.Date(c("2020-01-01", "2020-01-15", "2020-01-16", NA))
  ), id = "id", time = "day")
  inside <- as.data.frame(intersect_3d(points, dated))
  expect_identical(inside$feature_id_x, c("first", "last", "untimed"))
  expect_identical(inside$time[1:2],
                   as.POSIXct(c("2020-01-01", "2020-01-15"), tz = "UTC"))
  expect_identical(intersect_3d(points, dated)@time_resolution, "daily")

  # A polygon without a time takes points at any time.
  expect_identical(length(intersect_3d(points, polygons_0_100())), 4L)

  # With hours, the ends are exact instants.
  hourly <- sf::st_sf(start = as.POSIXct("2020-01-05 06:00", tz = "UTC"),
                      end = as.POSIXct("2020-01-05 12:00", tz = "UTC"),
                      geometry = sf::st_geometry(make_polygon()))
  hourly <- as_vector3d(hourly, time = c("start", "end"), depth_min = 0,
                        depth_max = 100)
  timed <- as_vector3d(data.frame(
    id = c("at_end", "past_end"), lon = 2.5, lat = 2.5, depth = 50,
    when = as.POSIXct(c("2020-01-05 12:00", "2020-01-05 12:01"), tz = "UTC")
  ), id = "id", time = "when")
  expect_identical(as.data.frame(intersect_3d(timed, hourly))$feature_id_x,
                   "at_end")
  # Dates on one side and hours on the other can't be compared.
  expect_error(intersect_3d(points, hourly), "dates and the other has hours",
               class = "ocean3d_error_time_resolution_differ")
})

test_that("two sets of points share no space", {
  points <- points_at(2.5, 2.5, 10)
  for (call in list(function() intersect_3d(points, points),
                    function() intersects_3d(points, points_at(1, 1, 1)))) {
    expect_error(call(), "Two sets of points share no space",
                 class = "ocean3d_error_query_input")
  }
})

test_that("intersect_3d() of two polygon sets shares outline, depth and time", {
  square <- function(x0, x1) {
    sf::st_polygon(list(rbind(c(x0, 0), c(x1, 0), c(x1, 1), c(x0, 1),
                              c(x0, 0))))
  }
  # a1 and a2: 0-100 m, 1 to 15 January.
  a <- sf::st_sf(id = c("a1", "a2"),
                 start = as.Date("2020-01-01"), end = as.Date("2020-01-15"),
                 geometry = sf::st_sfc(square(0, 2), square(5, 6), crs = 4326))
  a <- as_vector3d(a, id = "id", time = c("start", "end"),
                   depth_min = 0, depth_max = 100)
  # b1 overlaps a1 in outline and depth; b2 overlaps a1 in outline but only
  # touches its depth range; b3 overlaps a2. No times.
  b <- sf::st_sf(id = c("b1", "b2", "b3"), top = c(50, 100, 0),
                 geometry = sf::st_sfc(square(1, 3), square(0, 1),
                                       square(5.5, 7), crs = 4326))
  b <- as_vector3d(b, id = "id", depth_min = "top", depth_max = 200)

  shared <- intersect_3d(a, b)
  expect_s4_class(shared, "SpatVector3D")
  f <- as.data.frame(shared)
  expect_identical(f$feature_id_x, c("a1", "a2"))
  expect_identical(f$feature_id_y, c("b1", "b3"))
  expect_identical(f$depth_min, c(50, 0))
  expect_identical(f$depth_max, c(100, 100))
  # b has no times, so the pairs keep a's interval.
  expect_identical(f$time_start, a@features$time_start)
  expect_identical(shared@time_resolution, "daily")
  # a1 and b1 share a 1 x 1 degree square over 50-100 m.
  one_square <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(square(1, 2),
                                                            crs = 4326)),
                            depth_min = 50, depth_max = 100)
  expect_equal(volume(intersect_3d(a, b)),
               volume(one_square) + volume(as_vector3d(
                 sf::st_sf(geometry = sf::st_sfc(square(5.5, 6), crs = 4326)),
                 depth_min = 0, depth_max = 100
               )),
               # s2 puts the new corners where the geodesic edges cross, a
               # little off the parallels.
               tolerance = 1e-4)
})

test_that("intersects_3d() of polygon sets lists the pairs that share space", {
  square <- function(x0, x1) {
    sf::st_polygon(list(rbind(c(x0, 0), c(x1, 0), c(x1, 1), c(x0, 1),
                              c(x0, 0))))
  }
  a <- sf::st_sf(id = c("a1", "a2", "a3"),
                 geometry = sf::st_sfc(square(0, 2), square(5, 6),
                                       square(10, 11), crs = 4326))
  a <- as_vector3d(a, id = "id", depth_min = 0, depth_max = 100)
  # b1 and b3 overlap a1 in outline and depth; b2 overlaps a1 but only
  # touches its depth range; b4 only shares an edge with a2; nothing meets
  # a3.
  b <- sf::st_sf(id = c("b3", "b2", "b1", "b4"), top = c(0, 100, 50, 0),
                 geometry = sf::st_sfc(square(1, 3), square(0, 1),
                                       square(1.5, 2), square(6, 7),
                                       crs = 4326))
  b <- as_vector3d(b, id = "id", depth_min = "top", depth_max = 200)

  expect_identical(intersects_3d(a, b),
                   list(a1 = c("b3", "b1"), a2 = character(0),
                        a3 = character(0)))
  # The other way round lists the same pairs from b's side.
  expect_identical(intersects_3d(b, a),
                   list(b3 = "a1", b2 = character(0), b1 = "a1",
                        b4 = character(0)))
})

test_that("intersects_3d() of polygon sets skips outlines that only touch", {
  rect <- function(x0, x1, y0 = 0, y1 = 1) {
    sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1),
                              c(x0, y0))))
  }
  # a1: 0-100 m, 1 to 15 January; a2: 0-100 m, no time.
  a <- sf::st_sf(id = c("a1", "a2"),
                 start = as.Date(c("2020-01-01", NA)),
                 end = as.Date(c("2020-01-15", NA)),
                 geometry = sf::st_sfc(rect(0, 1), rect(10, 11), crs = 4326))
  a <- as_vector3d(a, id = "id", time = c("start", "end"),
                   depth_min = 0, depth_max = 100)
  cases <- list(
    list(label = "overlap", geometry = rect(0.5, 2), top = 0,
         start = NA, end = NA, a1 = TRUE),
    list(label = "whole shared edge", geometry = rect(1, 2), top = 0,
         start = NA, end = NA, a1 = FALSE),
    list(label = "part of an edge", geometry = rect(1, 2, 0.5, 0.7), top = 0,
         start = NA, end = NA, a1 = FALSE),
    # a1's corner lies on b1's edge, and the two edges overlap from there.
    list(label = "edges that overlap in part", geometry = rect(1, 2, 0.5, 3),
         top = 0, start = NA, end = NA, a1 = FALSE),
    list(label = "corner", geometry = rect(1, 2, 1, 2), top = 0,
         start = NA, end = NA, a1 = FALSE),
    list(label = "touching depths", geometry = rect(0.5, 2), top = 100,
         start = NA, end = NA, a1 = FALSE),
    # Dates include their end day, so sharing 15 January is an overlap.
    list(label = "one shared day", geometry = rect(0.5, 2), top = 0,
         start = "2020-01-15", end = "2020-01-31", a1 = TRUE),
    list(label = "next day", geometry = rect(0.5, 2), top = 0,
         start = "2020-01-16", end = "2020-01-31", a1 = FALSE)
  )
  old <- suppressMessages(sf::sf_use_s2())
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  for (s2 in c(TRUE, FALSE)) {
    suppressMessages(sf::sf_use_s2(s2))
    for (case in cases) {
      label <- glue::glue("{case$label}, s2 = {s2}")
      # b2 overlaps a2, which has no time, so it is listed whatever b's
      # times are.
      b <- sf::st_sf(id = c("b1", "b2"), top = c(case$top, 0),
                     start = as.Date(case$start), end = as.Date(case$end),
                     geometry = sf::st_sfc(case$geometry, rect(10.5, 12),
                                           crs = 4326))
      b <- as_vector3d(b, id = "id", time = c("start", "end"),
                       depth_min = "top", depth_max = 200)
      expected <- list(a1 = if (case$a1) "b1" else character(0), a2 = "b2")
      expect_identical(suppressMessages(intersects_3d(a, b)), expected,
                       label = label)
    }
  }
})

test_that("polygon sets with edges at computed coordinates share area", {
  # Corners at x0 / 3 lie on the slanted edge from (0, 0) to (3, 1) as near
  # as floating point allows. On top of the edge, the two share only a line;
  # a nanodegree above it, a sliver with an area. On the sphere, that sliver
  # comes with a line, as a collection, and it is the only pair.
  pg <- function(m) sf::st_polygon(list(rbind(m, m[1, ])))
  a <- as_vector3d(sf::st_sf(id = "a1", geometry = sf::st_sfc(
    pg(rbind(c(0, 0), c(3, 1), c(0, 1))), crs = 4326
  )), id = "id", depth_min = 0, depth_max = 100)
  box <- function(x0, lift) {
    p <- pg(rbind(c(x0, x0 / 3 + lift), c(x0 + 0.1, (x0 + 0.1) / 3 + lift),
                  c(x0 + 0.1, x0 / 3 - 0.5), c(x0, x0 / 3 - 0.5)))
    as_vector3d(sf::st_sf(id = "b1", geometry = sf::st_sfc(p, crs = 4326)),
                id = "id", depth_min = 0, depth_max = 100)
  }
  cases <- list(
    list(label = "on the edge", b = box(1.1, 0), n = 0L),
    list(label = "on the edge, at the corner", b = box(2.9, 0), n = 0L),
    list(label = "a sliver above it", b = box(2.9, 1e-9), n = 1L)
  )
  old <- suppressMessages(sf::sf_use_s2())
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  for (s2 in c(TRUE, FALSE)) {
    suppressMessages(sf::sf_use_s2(s2))
    for (case in cases) {
      label <- glue::glue("{case$label}, s2 = {s2}")
      shared <- suppressMessages(intersect_3d(a, case$b))
      expect_identical(length(shared), case$n, label = label)
      expect_identical(suppressMessages(intersects_3d(a, case$b)),
                       list(a1 = rep("b1", case$n)), label = label)
    }
  }
})

test_that("intersect_3d() of polygon sets takes an empty set", {
  far <- function(x0) {
    sf::st_polygon(list(rbind(c(x0, 0), c(x0 + 1, 0), c(x0 + 1, 1),
                              c(x0, 1), c(x0, 0))))
  }
  a <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(far(0), crs = 4326)),
                   depth_min = 0, depth_max = 100)
  b <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(far(5), crs = 4326)),
                   depth_min = 0, depth_max = 100)
  none <- suppressMessages(intersect_3d(a, b))
  old <- suppressMessages(sf::sf_use_s2())
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  for (s2 in c(TRUE, FALSE)) {
    suppressMessages(sf::sf_use_s2(s2))
    for (pair in list(list(none, a), list(a, none))) {
      shared <- suppressMessages(intersect_3d(pair[[1]], pair[[2]]))
      expect_s4_class(shared, "SpatVector3D")
      expect_identical(length(shared), 0L, label = glue::glue("s2 = {s2}"))
    }
    expect_identical(suppressMessages(intersects_3d(a, none)),
                     list(`1` = character(0)))
  }
})

test_that("intersect_3d() of polygon sets gives multipolygons on the sphere", {
  # With sf's spherical geometry, when every pair that meets shares a
  # polygon and one shares a multipolygon, every pair gets a multipolygon,
  # as sf::st_intersection() of two sets gives.
  rect <- function(x0, x1) {
    rbind(c(x0, 0), c(x1, 0), c(x1, 1), c(x0, 1), c(x0, 0))
  }
  two_parts <- sf::st_multipolygon(list(list(rect(0, 1)), list(rect(3, 4))))
  x <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(two_parts, crs = 4326)),
                   depth_min = 0, depth_max = 100)
  y <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(rect(0.5, 3.5))), sf::st_polygon(list(rect(0.2, 0.4))),
    crs = 4326
  )), depth_min = 0, depth_max = 100)
  old <- suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  geometry <- as.data.frame(intersect_3d(x, y))$geometry
  expect_identical(as.character(sf::st_geometry_type(geometry)),
                   c("MULTIPOLYGON", "MULTIPOLYGON"))
})

test_that("intersect_3d() of polygon sets says once that it is planar", {
  # With planar geometry, sf says that it treats longitude and latitude as
  # planar once per call.
  square <- function(x0) {
    sf::st_polygon(list(rbind(c(x0, 0), c(x0 + 2, 0), c(x0 + 2, 1),
                              c(x0, 1), c(x0, 0))))
  }
  x <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(square(0), crs = 4326)),
                   depth_min = 0, depth_max = 100)
  y <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(square(1), crs = 4326)),
                   depth_min = 0, depth_max = 100)
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  said <- character(0)
  withCallingHandlers(intersect_3d(x, y), message = function(m) {
    said <<- c(said, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  expect_identical(sum(grepl("planar", said)), 1L)
})

test_that("intersect_3d() of polygon sets keeps the polygon part of a mix", {
  # With planar geometry, a two-part polygon that overlaps a square in one
  # part and only touches it in the other gives a polygon and a line.
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  square <- function(x0, x1) {
    rbind(c(x0, 0), c(x1, 0), c(x1, 1), c(x0, 1), c(x0, 0))
  }
  two_parts <- sf::st_multipolygon(list(list(square(0, 1)),
                                        list(square(3, 4))))
  x <- as_vector3d(sf::st_sf(geometry = sf::st_sfc(two_parts, crs = 4326)),
                   depth_min = 0, depth_max = 100)
  y <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(square(0.5, 3))),
                                       crs = 4326))
  y <- as_vector3d(y, depth_min = 0, depth_max = 100)

  shared <- suppressMessages(intersect_3d(x, y))
  expect_identical(length(shared), 1L)
  geometry <- as.data.frame(shared)$geometry
  expect_identical(as.character(sf::st_geometry_type(geometry)), "POLYGON")
  expect_identical(unname(sf::st_bbox(geometry)[c("xmin", "xmax")]),
                   c(0.5, 1))
})

test_that("intersect_3d() of polygon sets compares time intervals", {
  on_days <- function(start, end) {
    p <- sf::st_sf(start = as.Date(start), end = as.Date(end),
                   geometry = sf::st_geometry(make_polygon()))
    as_vector3d(p, time = c("start", "end"), depth_min = 0, depth_max = 100)
  }
  january <- on_days("2020-01-01", "2020-01-15")
  cases <- list(
    list(label = "overlapping days", y = on_days("2020-01-10", "2020-01-20"),
         n = 1L, start = "2020-01-10", end = "2020-01-15"),
    # Dates include their end day, so sharing 15 January is an overlap.
    list(label = "one shared day", y = on_days("2020-01-15", "2020-01-31"),
         n = 1L, start = "2020-01-15", end = "2020-01-15"),
    list(label = "next day", y = on_days("2020-01-16", "2020-01-31"),
         n = 0L)
  )
  for (case in cases) {
    shared <- intersect_3d(january, case$y)
    expect_identical(length(shared), case$n, label = case$label)
    if (case$n > 0L) {
      f <- as.data.frame(shared)
      expect_identical(format(f$time_start), case$start, label = case$label)
      expect_identical(format(f$time_end), case$end, label = case$label)
    }
  }
  # Sharing nothing is an empty result, with no volume.
  expect_identical(volume(intersect_3d(january,
                                       on_days("2021-01-01", "2021-01-02"))),
                   0)

  # With hours, windows that meet at one instant overlap at that instant,
  # so a polygon at a single time can meet one whose window ends then.
  at_hours <- function(start, end) {
    p <- sf::st_sf(start = as.POSIXct(start, tz = "UTC"),
                   end = as.POSIXct(end, tz = "UTC"),
                   geometry = sf::st_geometry(make_polygon()))
    as_vector3d(p, time = c("start", "end"), depth_min = 0, depth_max = 100)
  }
  meeting <- intersect_3d(at_hours("2020-01-05 06:00", "2020-01-05 12:00"),
                          at_hours("2020-01-05 12:00", "2020-01-05 18:00"))
  expect_identical(length(meeting), 1L)
  expect_identical(as.data.frame(meeting)$time_start,
                   as.data.frame(meeting)$time_end)

  hourly <- sf::st_sf(time = as.POSIXct("2020-01-05 06:00", tz = "UTC"),
                      geometry = sf::st_geometry(make_polygon()))
  hourly <- as_vector3d(hourly, time = "time", depth_min = 0, depth_max = 10)
  expect_error(intersect_3d(january, hourly), "dates and the other has hours",
               class = "ocean3d_error_time_resolution_differ")
})

test_that("intersect_3d() of a SpatVector3D rejects 2D input", {
  # A 2D object isn't a second set of polygons.
  flat <- sf::st_geometry(make_polygon())
  points <- as_vector3d(data.frame(lon = 2.5, lat = 2.5, depth = 10))
  for (call in list(function() intersect_3d(polygons_0_100(), flat),
                    function() intersects_3d(flat, polygons_0_100()),
                    function() intersect_3d(points, flat))) {
    err <- expect_error(call(), class = "ocean3d_error_query_input")
    expect_match(conditionMessage(err),
                 "SpatVector3D polygons or points can be compared with")
    expect_match(conditionMessage(err), "not with sfc_POLYGON")
  }
})
