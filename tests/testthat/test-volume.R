test_that("calc_volume_overlap() totals match hand-computed volumes", {
  # Every cell is 1 km x 1 km, so each total is cells x thickness in km.
  # `tolerance` allows for the cell areas terra computes; `same` rows also
  # check that the overlap equals both volumes exactly.
  depths <- c(0, 50, 100, 150)
  centre <- c(NA, NA, NA, 1, 1, 1, NA, NA, NA)
  # A covers all 9 cells, 0-100 m
  a <- make_range_rast(rep(0, 9), rep(100, 9))

  cases <- list(
    # B covers only the centre 3 cells, 50-200 m: they share 50-100 m there
    list(label = "envelopes",
         a = a, b = make_range_rast(50 * centre, 200 * centre),
         expected = c(0.9, 0.45, 0.15), tolerance = 0.01),
    # full overlap when the ranges are identical
    list(label = "identical envelopes", a = a, b = a,
         expected = c(0.9, 0.9, 0.9), tolerance = 0.01, same = TRUE),
    # A occupies 0 and 50 m (0-100 m) in every cell; B occupies 0, 50 and
    # 100 m (0-150 m) in the centre 3 cells
    list(label = "voxels",
         a = make_voxel(cbind(rep(1, 9), rep(1, 9), rep(NA, 9), rep(NA, 9)),
                        depths),
         b = make_voxel(cbind(centre, centre, centre, rep(NA, 9)), depths),
         expected = c(0.9, 0.45, 0.3), tolerance = 0.01)
  )
  for (case in cases) {
    result <- calc_volume_overlap(case$a, case$b)
    totals <- terra::global(
      result[[c("volume_a", "volume_b", "volume_overlap")]], "sum",
      na.rm = TRUE
    )$sum

    expect_identical(totals, case$expected, tolerance = case$tolerance,
                     label = case$label)
    if (isTRUE(case$same)) {
      expect_identical(totals[[3]], totals[[1]], label = case$label)
      expect_identical(totals[[3]], totals[[2]], label = case$label)
    }
  }
})

test_that("calc_volume_overlap() is NA where ranges share no cells", {
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

test_that("calc_volume_overlap() overlap is 0 for disjoint depth ranges", {
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
  overlap_total <- terra::global(result[["volume_overlap"]], "sum",
                                 na.rm = TRUE)$sum

  expect_identical(overlap_total, 0)
  # Both present but vertically disjoint: volume is 0, the depth limits are NA.
  expect_true(all(terra::values(result[["volume_overlap"]]) == 0))
  expect_true(all(is.na(terra::values(result[["depth_max_overlap"]]))))
})

test_that("volume() sums the occupied water in km^3", {
  # Every cell is 1 km x 1 km, so a volume is 9 cells x thickness in km.
  # Voxel levels stand for slabs: under "top" a level names the top of its
  # slab, under "midpoint" the edges fall halfway to each neighbour.
  # `tolerance` allows for the cell areas terra computes; rows without one are
  # exact.
  # 9 cells, occupied at 0 m and 100 m but not at 300 m
  two_levels <- make_voxel(cbind(rep(1, 9), rep(1, 9), rep(NA, 9)),
                           depths = c(0, 100, 300))
  # occupied at 0, 200 and 300 m, with a gap at 100 m
  gap <- make_voxel(cbind(rep(1, 9), rep(NA, 9), rep(1, 9), rep(1, 9)),
                    depths = c(0, 100, 200, 300))
  e200 <- make_range_rast(rep(0, 9), rep(200, 9))
  scores <- make_voxel(cbind(rep(0.9, 9), rep(0.1, 9), rep(0.9, 9)),
                       depths = c(0, 100, 200))

  cases <- list(
    list(label = "envelope", x = make_range_rast(rep(0, 9), rep(100, 9)),
         expected = 0.9, tolerance = 0.01),
    # 0 m stands for 0-100 m and 100 m for 100-300 m: 300 m per cell
    list(label = "voxel, top", x = two_levels, expected = 9 * 0.3,
         tolerance = 0.01),
    # 0 m stands for 0-50 m and 100 m for 50-200 m: 200 m per cell
    list(label = "voxel, midpoint", x = two_levels, bounds = "midpoint",
         expected = 9 * 0.2, tolerance = 0.01),
    # slabs are 100, 100, 100 and 0 m, so the occupied levels give 100 +
    # 100 + 0 m
    list(label = "voxel with an interior gap", x = gap, expected = 9 * 0.2,
         tolerance = 0.01),
    # voxel_to_envelope() is lossy in exactly this way: it spans 0-300 m solid
    list(label = "the gap filled by voxel_to_envelope()",
         x = voxel_to_envelope(gap), expected = 9 * 0.3, tolerance = 0.01),
    # Under "top" the levels are 0, 100 and 200 m with slabs 100, 100 and 0 m,
    # so the discretization neither loses nor invents any water here.
    list(label = "envelope_to_voxel() on aligned levels",
         x = envelope_to_voxel(e200, depths = c(0, 100, 200)),
         expected = volume(e200)),
    # Default !is.na() counts all three levels (100 + 100 + 0 m). A threshold
    # drops the middle one, leaving only the 100 m slab the surface level
    # names.
    list(label = "default occupancy", x = scores, expected = 9 * 0.2,
         tolerance = 0.01),
    list(label = "occupied where x > 0.5",
         x = occupied(scores, function(x) x > 0.5), expected = 9 * 0.1,
         tolerance = 0.01),
    list(label = "occupied nowhere",
         x = occupied(scores, function(x) x > 0.99), expected = 0)
  )
  for (case in cases) {
    vol <- if (is.null(case$bounds)) {
      volume(case$x)
    } else {
      volume(case$x, bounds = case$bounds)
    }
    expect_identical(vol, case$expected, tolerance = case$tolerance,
                     label = case$label)
  }
})

# ---- SpatVoxel and mixed overlap --------------------------------------------

test_that("calc_volume_overlap() uses occupied, not spanned, voxel volume", {
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
  expect_identical(terra::global(result[["volume_a"]], "sum", na.rm = TRUE)$sum,
               9 * 0.3, tolerance = 0.01)
  expect_identical(terra::global(result[["volume_b"]], "sum", na.rm = TRUE)$sum,
               9 * 0.2, tolerance = 0.01)
  expect_identical(
    terra::global(result[["volume_overlap"]], "sum", na.rm = TRUE)$sum,
    9 * 0.2, tolerance = 0.01
  )
})

test_that("calc_volume_overlap() returns a SpatRaster, promoting envelopes", {
  # Either argument may be the envelope; both give the all-voxel answer. terra
  # keeps the input's class on `[[` and `c()`, so the nine-layer stack used to
  # come back tagged SpatEnvelope while failing that class's rules.
  depths <- c(0, 50, 100, 150, 200)

  a <- make_range_rast(rep(0, 9), rep(100, 9))
  b <- make_range_rast(rep(50, 9), rep(200, 9))
  av <- envelope_to_voxel(a, depths = depths)
  bv <- envelope_to_voxel(b, depths = depths)

  all_voxel <- calc_volume_overlap(av, bv)
  cases <- list(
    "two voxels" = all_voxel,
    "two envelopes" = calc_volume_overlap(a, b),
    "envelope and voxel" = calc_volume_overlap(a, bv),
    "voxel and envelope" = calc_volume_overlap(av, b)
  )
  for (combo in names(cases)) {
    out <- cases[[combo]]
    expect_identical(class(out)[[1]], "SpatRaster", label = combo)
    expect_identical(terra::nlyr(out), 9, label = combo)
  }
  expect_identical(terra::values(cases[["envelope and voxel"]]),
                   terra::values(all_voxel))
  expect_identical(terra::values(cases[["voxel and envelope"]]),
                   terra::values(all_voxel))
})

# ---- input errors -----------------------------------------------------------

test_that("volume functions reject input they cannot measure", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  # an envelope on another grid
  a_2x2 <- make_range_rast(rep(0, 4), rep(100, 4), ncol = 2, nrow = 2)
  v100 <- make_voxel(cbind(rep(1, 9), rep(1, 9)), depths = c(0, 100))
  v200 <- make_voxel(cbind(rep(1, 9), rep(1, 9)), depths = c(0, 200))

  # Carrying the right layer names is not the same as being a SpatEnvelope.
  bare <- methods::as(a, "SpatRaster")
  expect_s4_class(bare, "SpatRaster")
  expect_false(methods::is(bare, "SpatEnvelope"))

  cases <- list(
    # voxel-only arguments for envelopes
    list(call = function() volume(a, bounds = "midpoint"),
         pattern = "unused argument"),
    list(call = function() calc_volume_overlap(a, a, bounds = "midpoint"),
         pattern = "unused argument"),
    list(call = function() calc_volume_overlap(v100, v200),
         pattern = "same depth levels"),
    list(call = function() calc_volume_overlap(a, a_2x2),
         pattern = "same grid"),
    # bare SpatRasters, alone or with a real envelope in either position
    list(call = function() volume(bare), pattern = "not a bare SpatRaster"),
    list(call = function() calc_volume_overlap(bare, bare),
         pattern = "not a bare SpatRaster"),
    list(call = function() calc_volume_overlap(a, bare),
         pattern = "not a bare SpatRaster"),
    list(call = function() calc_volume_overlap(bare, a),
         pattern = "not a bare SpatRaster")
  )
  for (case in cases) {
    expect_error(case$call(), case$pattern,
                 label = paste(deparse(body(case$call)), collapse = " "))
  }
})
