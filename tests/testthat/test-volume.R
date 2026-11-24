test_that("volume() sums the occupied water in km^3", {
  # Every cell is 1 km x 1 km, so a volume is 9 cells x thickness in km.
  # Voxel levels stand for slabs: under "top" a level names the top of its
  # slab, under "midpoint", the default, the edges fall halfway to each
  # neighbour.
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
    list(label = "voxel, top", x = two_levels, bounds = "top",
         expected = 9 * 0.3, tolerance = 0.01),
    # 0 m stands for 0-50 m and 100 m for 50-200 m: 200 m per cell
    list(label = "voxel, midpoint", x = two_levels, bounds = "midpoint",
         expected = 9 * 0.2, tolerance = 0.01),
    list(label = "voxel, default", x = two_levels,
         expected = volume(two_levels, bounds = "midpoint")),
    # slabs are 50, 100, 100 and 50 m, so the occupied levels give 50 +
    # 100 + 50 m
    list(label = "voxel with an interior gap", x = gap, expected = 9 * 0.2,
         tolerance = 0.01),
    # voxel_to_envelope() is lossy in exactly this way: it spans 0-300 m solid
    list(label = "the gap filled by voxel_to_envelope()",
         x = voxel_to_envelope(gap), expected = 9 * 0.3, tolerance = 0.01),
    # The levels 0, 100 and 200 m stand for 0-50, 50-150 and 150-200 m, which
    # tile 0-200 m, so the discretization neither loses nor invents any water
    # here.
    list(label = "envelope_to_voxel() on aligned levels",
         x = envelope_to_voxel(e200, depths = c(0, 100, 200)),
         expected = volume(e200)),
    # Default !is.na() counts all three levels (50 + 100 + 50 m). A threshold
    # drops the middle one, leaving the 50 m slabs of the outer levels.
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

# ---- per_cell ---------------------------------------------------------------

test_that("volume(per_cell = TRUE) sums to the total", {
  depths <- c(0, 100, 200, 300)
  gap <- make_voxel(cbind(rep(1, 9), rep(NA, 9), rep(1, 9), rep(1, 9)),
                    depths)
  edged <- as_voxel(make_voxel(1, c(0.5, 1.5, 3)),
                    depth_bounds = cbind(c(0, 1, 2), c(1, 2, 4)))
  cases <- list(
    "envelope" = list(x = make_range_rast(c(0:7, NA), c(100:107, NA))),
    "voxel, midpoint" = list(x = gap, bounds = "midpoint"),
    "voxel, top" = list(x = gap, bounds = "top"),
    "voxel with depth edges" = list(x = edged)
  )
  for (label in names(cases)) {
    case <- cases[[label]]
    args <- if (is.null(case$bounds)) list() else list(bounds = case$bounds)
    total <- do.call(volume, c(list(case$x), args))
    per_cell <- do.call(volume, c(list(case$x, per_cell = TRUE), args))

    expect_identical(class(per_cell)[[1]], "SpatRaster", label = label)
    expect_identical(names(per_cell), "volume", label = label)
    expect_true(terra::compareGeom(per_cell, case$x), label = label)
    expect_equal(terra::global(per_cell, "sum", na.rm = TRUE)[[1]], total,
                 label = label)
  }
})

test_that("volume(per_cell = TRUE) is NA outside the domain", {
  # Cell 1 is outside the envelope; cell 2 has its depths but zero thickness.
  e <- make_range_rast(c(NA, 50, rep(0, 7)), c(NA, 50, rep(100, 7)))
  vals <- terra::values(volume(e, per_cell = TRUE))[, 1]
  expect_true(is.na(vals[1]))
  expect_identical(vals[2], 0)
  expect_true(all(vals[3:9] > 0))

  # Cell 1 has no occupied level: NA, as outside an envelope, not 0. Cell 2
  # is occupied only where the value fails the occupancy rule.
  v <- make_voxel(cbind(c(NA, 0.1, rep(1, 7)), c(NA, 0.1, rep(1, 7))),
                  depths = c(0, 100))
  vals <- terra::values(volume(occupied(v, function(x) x > 0.5),
                               per_cell = TRUE))[, 1]
  expect_true(all(is.na(vals[1:2])))
  expect_true(all(vals[3:9] > 0))

  # A voxel occupied nowhere is NA everywhere.
  nowhere <- occupied(v, function(x) x > 2)
  expect_true(all(is.na(terra::values(volume(nowhere, per_cell = TRUE)))))
})

# ---- per_level (#116) -------------------------------------------------------

test_that("volume(per_level = TRUE) gives each level its slab times the cell", {
  # Uneven levels: under "midpoint" the slabs are 25, 75, 100 and 50 m thick,
  # and every cell is 1 km x 1 km.
  depths <- c(0, 50, 150, 250)
  v <- make_voxel(cbind(rep(1, 9), c(NA, rep(1, 8)), rep(1, 9), rep(NA, 9)),
                  depths)
  per_level <- volume(v, per_cell = TRUE, per_level = TRUE)

  expect_s4_class(per_level, "SpatVoxel")
  expect_identical(depths(per_level), depths)
  expect_identical(names(per_level), paste0("volume_depth=", depths))
  vals <- unname(terra::values(per_level, mat = TRUE))
  expect_equal(vals[1, ], c(0.025, NA, 0.1, NA), tolerance = 0.01)
  expect_equal(vals[2, 2], 0.075, tolerance = 0.01)
  # Unoccupied levels are NA, as everywhere a shape is absent.
  expect_true(all(is.na(vals[, 4])))

  # The levels sum to the per-cell volume, and keep the voxel's edges.
  expect_equal(rowSums(vals, na.rm = TRUE),
               terra::values(volume(v, per_cell = TRUE))[, 1])
  edged <- as_voxel(v, depth_bounds = cbind(c(0, 40, 100, 200),
                                            c(40, 100, 200, 300)))
  per_level <- volume(edged, per_cell = TRUE, per_level = TRUE)
  expect_equal(depth_bounds(per_level), depth_bounds(edged))
  expect_equal(unname(terra::values(per_level)[2, 1]), 0.04, tolerance = 0.01)
})

test_that("per_level needs per_cell and a voxel", {
  v <- make_voxel(1, c(0, 100))
  cases <- list(
    list(x = v, args = list(per_level = TRUE), pattern = "needs `per_cell"),
    list(x = v, args = list(per_cell = TRUE, per_level = NA),
         pattern = "TRUE or FALSE"),
    list(x = make_range_rast(rep(0, 9), rep(100, 9)),
         args = list(per_cell = TRUE, per_level = TRUE),
         pattern = "Envelopes have no depth levels"),
    list(x = as_vector3d(make_polygon(), depth_min = 0, depth_max = 100),
         args = list(per_cell = TRUE, per_level = TRUE),
         pattern = "Polygons have no depth levels")
  )
  for (case in cases) {
    expect_error(do.call(volume, c(list(case$x), case$args)), case$pattern,
                 class = "ocean3d_error_volume_per_level")
  }
})

# ---- measuring an overlap ---------------------------------------------------

test_that("volume(intersect_3d()) matches hand-computed overlaps", {
  # Every cell is 1 km x 1 km, so each total is cells x thickness in km.
  # `tolerance` allows for the cell areas terra computes. Each case gives the
  # volumes of a, b and their overlap.
  depths <- c(0, 50, 100, 150)
  centre <- c(NA, NA, NA, 1, 1, 1, NA, NA, NA)
  # A covers all 9 cells, 0-100 m
  a <- make_range_rast(rep(0, 9), rep(100, 9))

  cases <- list(
    # B covers only the centre 3 cells, 50-200 m: they share 50-100 m there
    "envelopes" = list(a = a, b = make_range_rast(50 * centre, 200 * centre),
                       expected = c(0.9, 0.45, 0.15)),
    "identical envelopes" = list(a = a, b = a, expected = c(0.9, 0.9, 0.9)),
    # By default the levels stand for 0-25, 25-75, 75-125 and 125-150 m. A
    # occupies 0 and 50 m (75 m) in every cell; B occupies 0, 50 and 100 m
    # (125 m) in the centre 3 cells, and they share 75 m there.
    "voxels" = list(
      a = make_voxel(cbind(rep(1, 9), rep(1, 9), rep(NA, 9), rep(NA, 9)),
                     depths),
      b = make_voxel(cbind(centre, centre, centre, rep(NA, 9)), depths),
      expected = c(0.675, 0.375, 0.225)
    ),
    # Under "top" each of the upper three levels stands for 100 m. A is solid
    # over them; B has a gap at 100 m, so it occupies 200 m, not the 300 m
    # it spans, and that is all they share.
    "voxels with a gap" = list(
      a = make_voxel(cbind(rep(1, 9), rep(1, 9), rep(1, 9), rep(NA, 9)),
                     c(0, 100, 200, 300)),
      b = make_voxel(cbind(rep(1, 9), rep(NA, 9), rep(1, 9), rep(NA, 9)),
                     c(0, 100, 200, 300)),
      bounds = "top", expected = c(2.7, 1.8, 1.8)
    )
  )
  for (label in names(cases)) {
    case <- cases[[label]]
    args <- if (is.null(case$bounds)) list() else list(bounds = case$bounds)
    vol <- function(x, ...) do.call(volume, c(list(x, ...), args))
    ov <- intersect_3d(case$a, case$b)
    totals <- c(vol(case$a), vol(case$b), vol(ov))
    expect_equal(totals, case$expected, tolerance = 0.01, label = label)

    # The overlap never exceeds either input, in total or in any cell.
    expect_lte(totals[[3]], min(totals[1:2]), label = label)
    cells <- terra::values(c(vol(case$a, per_cell = TRUE),
                             vol(case$b, per_cell = TRUE),
                             vol(ov, per_cell = TRUE)))
    shared <- !is.na(cells[, 3])
    expect_true(all(cells[shared, 3] <= cells[shared, 1] + 1e-12 &
                      cells[shared, 3] <= cells[shared, 2] + 1e-12),
                label = label)
  }
})

# ---- input errors -----------------------------------------------------------

test_that("volume() rejects input it cannot measure", {
  a <- make_range_rast(rep(0, 9), rep(100, 9))
  v100 <- make_voxel(cbind(rep(1, 9), rep(1, 9)), depths = c(0, 100))

  # Carrying the right layer names is not the same as being a SpatEnvelope.
  bare <- methods::as(a, "SpatRaster")
  expect_false(methods::is(bare, "SpatEnvelope"))

  for (x in list(a, v100)) {
    err <- expect_error(volume(x, per_cell = NA),
                        class = "ocean3d_error_per_cell")
    expect_match(conditionMessage(err), "`per_cell` must be TRUE or FALSE")
  }

  for (call in list(function() volume(bare),
                    function() volume(bare, per_cell = TRUE))) {
    err <- expect_error(call(), class = "ocean3d_error_not_3d")
    expect_match(conditionMessage(err),
                 "must be a SpatVoxel or SpatEnvelope; got SpatRaster")
  }

  cases <- list(
    # voxel-only arguments for envelopes
    list(call = function() volume(a, bounds = "midpoint"),
         pattern = "unused argument")
  )
  for (case in cases) {
    expect_error(case$call(), case$pattern,
                 label = paste(deparse(body(case$call)), collapse = " "))
  }
})

# ---- depth edges (#138) -----------------------------------------------------

test_that("volume() uses a voxel's depth edges before `bounds`", {
  # Levels at 0.5, 1.5 and 3 m. Under "midpoint" the column is 0.5-3 m
  # (2.5 m); the stored edges make it 0-4 m.
  depths <- c(0.5, 1.5, 3)
  plain <- make_voxel(1, depths)
  edged <- as_voxel(plain, depth_bounds = cbind(c(0, 1, 2), c(1, 2, 4)))

  expect_equal(volume(edged), volume(plain, bounds = "midpoint") * 4 / 2.5)
  expect_equal(volume(edged, bounds = "top"),
               volume(edged, bounds = "midpoint"))

  # A layer subset keeps the edges of the layers it kept: 1-2 and 2-4 m.
  expect_equal(volume(edged[[2:3]]), volume(edged) * 3 / 4)
})


# SpatVector3D polygons (#115) ----------------------------------------------

test_that("volume() of polygons is area times depth range", {
  # make_polygon() covers 3 x 3 one-degree cells of make_template(), so its
  # volume is close to that of its envelope on the grid. Not equal: the
  # polygon's edges are geodesics, which bow away from the parallels that
  # bound the cells, by about 0.02% here.
  polygons <- as_vector3d(make_polygon(), depth_min = 0, depth_max = 100)
  envelope <- vector3d_to_envelope(polygons, make_template())
  expect_equal(volume(polygons), volume(envelope), tolerance = 1e-3)

  # Two polygons with the same outline, 100 m and 10 m thick: the total is
  # the sum of each, counting their overlap twice.
  two <- sf::st_sf(top = c(0, 50), bottom = c(100, 60),
                   geometry = rep(sf::st_geometry(make_polygon()), 2))
  two <- as_vector3d(two, depth_min = "top", depth_max = "bottom")
  expect_equal(volume(two), 1.1 * volume(polygons))
})

test_that("volume() of polygons ignores time and allows no thickness", {
  geometry <- sf::st_geometry(make_polygon())
  dated <- sf::st_sf(start = as.Date("2020-01-01"),
                     end = as.Date("2020-12-31"), geometry = geometry)
  timed <- as_vector3d(dated, time = c("start", "end"),
                       depth_min = 0, depth_max = 100)
  untimed <- as_vector3d(make_polygon(), depth_min = 0, depth_max = 100)
  expect_identical(volume(timed), volume(untimed))

  flat <- as_vector3d(make_polygon(), depth_min = 20, depth_max = 20)
  expect_identical(volume(flat), 0)
})

test_that("volume() of polygons rejects per_cell and voxel arguments", {
  polygons <- as_vector3d(make_polygon(), depth_min = 0, depth_max = 100)
  expect_error(volume(polygons, per_cell = TRUE), "Polygons have no cells",
               class = "ocean3d_error_volume_per_cell")
  expect_error(volume(polygons, bounds = "top"), "unused argument",
               class = "ocean3d_error_unused_argument")
})

# ---- SpatProxy: volume through time -----------------------------------------

test_that("volume() of a proxy gives one row per time step", {
  # Values 1..24 over 3 x 2 cells, 2 depths, 2 steps; a cutoff leaves some
  # cells unoccupied at some levels.
  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)
  p <- g$proxy
  warm <- function(v) v > 8

  out <- volume(p)
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("time", "volume_km3"))
  expect_equal(out$time, p@axes$time$cftime$as_timestamp(asPOSIX = TRUE))
  expect_equal(out$volume_km3,
               c(volume(as_voxel(p, time = 1)), volume(as_voxel(p, time = 2))))

  with_fun <- volume(p, fun = warm)
  expect_equal(with_fun$volume_km3, vapply(1:2, function(i) {
    volume(occupied(as_voxel(p, time = i), warm))
  }, numeric(1)))
  expect_lt(with_fun$volume_km3[1L], out$volume_km3[1L])

  # Per cell: one layer per step, dated.
  cells <- volume(p, per_cell = TRUE, fun = warm)
  expect_s4_class(cells, "SpatRaster")
  expect_false(methods::is(cells, "SpatVoxel"))
  expect_identical(names(cells), c("volume", "volume"))
  expect_equal(terra::time(cells), out$time)
  for (i in 1:2) {
    expect_equal(
      terra::values(cells[[i]]),
      terra::values(volume(occupied(as_voxel(p, time = i), warm),
                           per_cell = TRUE))
    )
  }
})

test_that("volume() of a proxy without dates gives one undated row", {
  flat <- proxy_and_file(time = NULL)
  on.exit(unlink(flat$file), add = TRUE)
  out <- volume(flat$proxy)
  expect_identical(nrow(out), 1L)
  expect_s3_class(out$time, "POSIXct")
  expect_true(is.na(out$time))
  expect_equal(out$volume_km3, volume(as_voxel(flat$proxy)))
  expect_false(terra::timeInfo(volume(flat$proxy, per_cell = TRUE))$time)

  clim <- proxy_and_file(time = 6,
                         time_units = "months since 0000-01-01 00:00:00")
  on.exit(unlink(clim$file), add = TRUE)
  expect_true(is.na(volume(clim$proxy)$time))
})

test_that("reading a step in groups of depth layers gives the same volume", {
  # Four levels, with a gap in the data so some groups have no occupied level
  # in some cells.
  values <- seq_len(3 * 2 * 4 * 2)
  values[c(1, 7, 13, 19, 2, 8)] <- NA
  g <- proxy_and_file(depth = c(5, 15, 25, 35), values = values)
  on.exit(unlink(g$file), add = TRUE)
  whole <- volume(g$proxy)
  whole_cells <- volume(g$proxy, per_cell = TRUE)

  # A budget below one layer's size reads one layer at a time.
  rlang::local_options(ocean3d.max_read_mb = 1e-9)
  expect_equal(volume(g$proxy), whole)
  expect_equal(terra::values(volume(g$proxy, per_cell = TRUE)),
               terra::values(whole_cells))
  # A voxel in memory without stored edges keeps the proxy's halfway slabs
  # in every group, so one level at a time gives the whole volume.
  v <- as_voxel(terra::rast(nrows = 2, ncols = 2, nlyrs = 4, vals = 1),
                depths = c(5, 15, 25, 35))
  rlang::local_options(ocean3d.max_read_mb = 1024)
  whole <- volume(as_proxy(v))$volume_km3
  whole_cells <- terra::values(volume(as_proxy(v), per_cell = TRUE))
  expect_gt(whole, 0)
  rlang::local_options(ocean3d.max_read_mb = 1e-9)
  expect_equal(volume(as_proxy(v))$volume_km3, whole)
  expect_equal(terra::values(volume(as_proxy(v), per_cell = TRUE)),
               whole_cells)
})

test_that("a proxy cropped across the antimeridian is measured in groups", {
  # A -180 to 180 grid cropped to the cells either side of 180 degrees: the
  # depth groups go through the same reader, which joins the two ends.
  globe <- make_grid_file(lon = seq(-179.5, 179.5, by = 1),
                          lat = c(10.5, 11.5), depth = c(5, 15, 25))
  on.exit(unlink(globe), add = TRUE)
  small <- terra::crop(as_proxy(globe, var = "temp"),
                       terra::ext(179.2, 180.8, 10, 12))
  expect_equal(small@axes$lon$centres, c(179.5, 180.5))
  whole <- volume(small, per_cell = TRUE)
  expect_equal(terra::values(whole[[1]]),
               terra::values(volume(as_voxel(small, time = 1),
                                    per_cell = TRUE)))
  rlang::local_options(ocean3d.max_read_mb = 1e-9)
  expect_equal(terra::values(volume(small, per_cell = TRUE)),
               terra::values(whole))
  expect_equal(volume(small)$volume_km3,
               vapply(1:2, function(i) volume(as_voxel(small, time = i)), 0))
})

test_that("volume() of a series measures every file's steps", {
  f1 <- make_grid_file(time = c(0, 1))
  f2 <- make_grid_file(time = c(2, 3), values = 2 * seq_len(24))
  on.exit(unlink(c(f1, f2)), add = TRUE)
  p <- as_proxy(c(f1, f2), var = "temp")
  out <- volume(p, fun = function(v) v > 30)
  expect_identical(nrow(out), 4L)
  expect_equal(out$volume_km3, vapply(1:4, function(i) {
    volume(occupied(as_voxel(p, time = i), function(v) v > 30))
  }, numeric(1)))
  expect_equal(out$volume_km3[1:2], c(0, 0))
})

test_that("volume() of a proxy rejects what it can't measure", {
  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)
  expect_error(volume(g$proxy, bounds = "top"), "unused argument",
               class = "ocean3d_error_unused_argument")
  expect_error(volume(g$proxy, per_cell = TRUE, per_level = TRUE),
               "measured one time step at a time",
               class = "ocean3d_error_volume_per_level")
  expect_error(volume(g$proxy, fun = "warm"), "must be NULL or a function",
               class = "ocean3d_error_volume_fun")
  for (budget in list(NA, "512", -1, c(1, 2))) {
    rlang::local_options(ocean3d.max_read_mb = budget)
    expect_error(volume(g$proxy), "`ocean3d.max_read_mb` must be one",
                 class = "ocean3d_error_read_budget")
  }
  rlang::local_options(ocean3d.max_read_mb = NULL)
  no_depth <- proxy_and_file(depth = NULL)
  on.exit(unlink(no_depth$file), add = TRUE)
  expect_error(volume(no_depth$proxy), "no depth axis",
               class = "ocean3d_error_proxy_no_depth")
})
