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

test_that("extract_to_area() keeps the values inside the area", {
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())
  out <- extract_to_area(a, v, min_depth = 40, max_depth = 600)
  # The same values as masking by hand.
  shape <- as_vector3d(a, depth_min = 40, depth_max = 600)
  expect_identical(terra::values(out),
                   terra::values(mask(terra::crop(v, shape), shape)))
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
  from_sfc <- extract_to_area(poly, v, min_depth = 0, max_depth = 0)

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
    out <- extract_to_area(case$area, v, min_depth = 0, max_depth = 0)

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
         rast_3d = make_area_voxel(), pattern = "sf, sfc, or SpatVector")
  )
  for (case in cases) {
    expect_error(extract_to_area(case$area, case$rast_3d), case$pattern,
                 class = "ocean3d_error_area_input", label = case$label)
  }
})
