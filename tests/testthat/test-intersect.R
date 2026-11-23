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

test_that("points and lines aren't compared with a grid yet", {
  points <- as_vector3d(data.frame(lon = 2.5, lat = 2.5, depth = 10))
  for (call in list(function() intersect_3d(points, envelope_50_200()),
                    function() intersects_3d(envelope_50_200(), points),
                    function() mask(envelope_50_200(), points))) {
    expect_error(call(), "This needs polygons, but these are points",
                 class = "ocean3d_error_vector3d_type")
    expect_error(call(), "buffer_3d() turns points into polygons",
                 class = "ocean3d_error_vector3d_type", fixed = TRUE)
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

test_that("intersect_3d() of polygon sets rejects points and lines", {
  points <- as_vector3d(data.frame(lon = 2.5, lat = 2.5, depth = 10))
  expect_error(intersect_3d(polygons_0_100(), points),
               "This needs polygons, but these are points",
               class = "ocean3d_error_vector3d_type")
  # A 2D object isn't a second set of polygons.
  flat <- sf::st_geometry(make_polygon())
  for (call in list(function() intersect_3d(polygons_0_100(), flat),
                    function() intersects_3d(flat, polygons_0_100()))) {
    err <- expect_error(call(), class = "ocean3d_error_query_input")
    expect_match(conditionMessage(err),
                 "SpatVector3D polygons can be compared with")
    expect_match(conditionMessage(err), "not with sfc_POLYGON")
  }
})
