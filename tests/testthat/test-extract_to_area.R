# Helper: synthetic multi-depth voxel with standard layer naming. Covers a 4x4
# grid in lon/lat with values that make each cell / layer distinguishable.
# Kept here: no shared helper (helper-example_objects.R) has this extent,
# grid or values.
make_area_voxel <- function(depths = c(0, 50, 100, 500, 1000),
                            variable = "tan",
                            ncol = 4, nrow = 4) {
  layers <- lapply(seq_along(depths), function(i) {
    r <- terra::rast(
      nrows = nrow, ncols = ncol,
      xmin = -10, xmax = 10, ymin = -10, ymax = 10,
      crs = "EPSG:4326"
    )
    terra::values(r) <- seq_len(ncol * nrow) + depths[i] * 0.01
    r
  })
  as_voxel(terra::rast(layers), depths = depths, varname = variable)
}

# Unlike make_polygon(), an sfc square centred on the raster.
make_area_polygon <- function() {
  # Square covering the centre of the raster.
  sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(-5, -5), c(5, -5), c(5, 5), c(-5, 5), c(-5, -5)
    ))),
    crs = "EPSG:4326"
  )
}

test_that("extract_to_area() keeps the levels whose slabs the band reaches", {
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())

  # Midpoint slabs of 0, 50, 100, 500 and 1000 m: 0-25, 25-75, 75-300,
  # 300-750 and 750-1000 m. `kept` are the levels with values inside the area.
  cases <- list(
    list(label = "40 to 600 m", min = 40, max = 600, kept = c(50, 100, 500)),
    list(label = "no bounds", min = NULL, max = NULL,
         kept = c(0, 50, 100, 500, 1000)),
    list(label = "min only", min = 90, max = NULL, kept = c(100, 500, 1000)),
    list(label = "max only", min = NULL, max = 90, kept = c(0, 50, 100)),
    # A band that only touches slab edges keeps only the slab it fills.
    list(label = "25 to 75 m", min = 25, max = 75, kept = 50)
  )
  for (case in cases) {
    out <- extract_to_area(a, v, min_depth = case$min, max_depth = case$max)

    # Every layer comes back, with its depth; levels outside the band are NA.
    expect_identical(depths(out), depths(v), label = case$label)
    present <- terra::global(out, "notNA")[, 1] > 0
    expect_identical(depths(out)[present], case$kept, label = case$label)
  }

  # Under bounds = "top" the 40 m top falls in the 0-50 m slab of the 0 m
  # level, and the band ends inside the 500-1000 m slab.
  top <- extract_to_area(a, v, min_depth = 40, max_depth = 600, bounds = "top")
  present <- terra::global(top, "notNA")[, 1] > 0
  expect_identical(depths(top)[present], c(0, 50, 100, 500))
})

test_that("extract_to_area() without a band keeps every level", {
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())
  # Under either `bounds`, including "top", where the deepest level has no
  # thickness.
  for (bounds in c("midpoint", "top")) {
    out <- extract_to_area(a, v, bounds = bounds)
    expect_true(all(terra::global(out, "notNA")[, 1] == 4), label = bounds)
  }
  # With only min_depth, the band runs to the deepest level, which "top"
  # gives no thickness.
  top <- extract_to_area(a, v, min_depth = 0, bounds = "top")
  expect_identical(terra::global(top, "notNA")[, 1], c(4, 4, 4, 4, 0))
})

test_that("extract_to_area() keeps the values inside the area", {
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())
  out <- extract_to_area(a, v, min_depth = 40, max_depth = 600)
  # The same values as masking by hand.
  shape <- as_vector3d(a, depth_min = 40, depth_max = 600)
  expect_identical(terra::values(out),
                   terra::values(mask(terra::crop(v, shape), shape)))
})

test_that("extract_to_area() gives a voxel and its proxy the same cells", {
  v <- make_area_voxel()
  # A square whose edges fall inside cells: every cell it touches counts.
  square <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-6, -6), c(6, -6), c(6, 6), c(-6, 6), c(-6, -6)
  ))), crs = 4326))
  from_voxel <- extract_to_area(square, v, min_depth = 40, max_depth = 600)
  from_proxy <- extract_to_area(square, as_proxy(v), 40, 600)
  expect_identical(terra::values(from_voxel), terra::values(from_proxy))
  expect_identical(terra::global(from_voxel, "notNA")[, 1],
                   c(0, 16, 16, 16, 0))
})

test_that("extract_to_area() finds an area on a 0-360 voxel", {
  r <- terra::rast(nrows = 2, ncols = 4, xmin = 170, xmax = 190, ymin = -5,
                   ymax = 5, nlyrs = 2, crs = "EPSG:4326", vals = 1)
  v <- as_voxel(r, depths = c(0, 100))
  area <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-179, -4), c(-176, -4), c(-176, -1), c(-179, -1), c(-179, -4)
  ))), crs = 4326))
  out <- extract_to_area(area, v)
  expect_identical(terra::global(out, "notNA")[, 1], c(1, 1))
})

test_that("extract_to_area() warns about a band that reaches no level", {
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())
  for (band in list(c(50, 50), c(1200, 1500), c(1200, NA))) {
    max_depth <- if (is.na(band[2])) NULL else band[2]
    expect_warning(
      out <- extract_to_area(a, v, min_depth = band[1], max_depth = max_depth),
      "overlaps none of the voxel's depth levels",
      class = "ocean3d_warning_area_band"
    )
    expect_true(all(terra::global(out, "notNA")[, 1] == 0))
  }
})

test_that("extract_to_area() doesn't warn when the band reaches empty levels", {
  # Levels 0, 50 and 100 m with no data at 50 m: a 30-70 m band reaches only
  # that level, so the answer is a voxel with no values, without a warning.
  r <- terra::rast(nrows = 4, ncols = 4, xmin = -10, xmax = 10, ymin = -10,
                   ymax = 10, nlyrs = 3, crs = "EPSG:4326")
  terra::values(r) <- cbind(rep(1, 16), NA, rep(1, 16))
  v <- as_voxel(r, depths = c(0, 50, 100))
  a <- sf::st_sf(geometry = make_area_polygon())
  expect_no_warning(out <- extract_to_area(a, v, 30, 70))
  expect_true(all(terra::global(out, "notNA")[, 1] == 0))
  # A band below the voxel still overlaps no level, so it warns.
  expect_warning(extract_to_area(a, v, 200, 300),
                 class = "ocean3d_warning_area_band")
})

test_that("extract_to_area() uses only the area's outline", {
  v <- make_area_voxel()
  poly <- make_area_polygon()
  # Columns named like the ones SpatVector3D uses don't get in the way.
  a <- sf::st_sf(depth_max = 7, feature_id = "reef", geometry = poly)
  expect_identical(
    terra::values(extract_to_area(a, v, 40, 600)),
    terra::values(extract_to_area(sf::st_sf(geometry = poly), v, 40, 600))
  )
})

test_that("extract_to_area() reads a SpatProxy", {
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())
  out <- extract_to_area(a, as_proxy(v), min_depth = 40, max_depth = 600)
  expect_identical(terra::values(out),
                   terra::values(extract_to_area(a, v, 40, 600)))
})

test_that("extract_to_area() returns a SpatVoxel", {
  out <- extract_to_area(sf::st_sf(geometry = make_area_polygon()),
                         make_area_voxel())
  expect_s4_class(out, "SpatVoxel")
  expect_true(methods::validObject(out))
})

test_that("extract_to_area() crops to sf, sfc and SpatVector areas", {
  v <- make_area_voxel()
  poly <- make_area_polygon()
  from_sfc <- extract_to_area(poly, v, min_depth = 0, max_depth = 100)

  # Rows with `same_extent` give the same extent as the sfc polygon. An area
  # in another CRS is reprojected onto the voxel's.
  cases <- list(
    list(label = "sfc", area = poly, same_extent = TRUE),
    list(label = "sf", area = sf::st_sf(geometry = poly), same_extent = TRUE),
    list(label = "SpatVector", area = terra::vect(poly), same_extent = TRUE),
    list(label = "sfc in EPSG:3857",
         area = sf::st_transform(poly, "EPSG:3857"), same_extent = FALSE)
  )
  for (case in cases) {
    out <- extract_to_area(case$area, v, min_depth = 0, max_depth = 100)

    expect_true(terra::same.crs(out, v), label = case$label)
    # The extent shrinks from -10..10 to roughly -5..5.
    e <- terra::ext(out)
    expect_lt(e[2] - e[1], 20, label = case$label)
    if (case$same_extent) {
      expect_identical(as.vector(e), as.vector(terra::ext(from_sfc)),
                       label = case$label)
    }
  }
})

test_that("extract_to_area() rejects input it can't use", {
  r <- terra::rast(nrows = 4, ncols = 4)
  terra::values(r) <- 1:16
  cases <- list(
    list(label = "rast_3d that is not a SpatVoxel",
         area = sf::st_sf(geometry = make_area_polygon()), rast_3d = r,
         pattern = "SpatVoxel or a SpatProxy"),
    list(label = "area that is not vector geometry", area = "nope",
         rast_3d = make_area_voxel(), pattern = "sf, sfc, or SpatVector"),
    list(label = "area without a CRS",
         area = terra::vect("POLYGON ((-5 -5, 5 -5, 5 5, -5 5, -5 -5))"),
         rast_3d = make_area_voxel(), pattern = "has no CRS")
  )
  for (case in cases) {
    expect_error(extract_to_area(case$area, case$rast_3d), case$pattern,
                 class = "ocean3d_error_area_input", label = case$label)
  }
  # A band whose top is deeper than its bottom.
  expect_error(
    extract_to_area(sf::st_sf(geometry = make_area_polygon()),
                    make_area_voxel(), min_depth = 600, max_depth = 40),
    "depth_min must not be deeper than depth_max",
    class = "ocean3d_error_vector3d"
  )
})
