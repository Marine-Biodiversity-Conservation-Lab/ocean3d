# Helper: synthetic multi-depth voxel with standard layer naming. Covers a 4x4
# grid in lon/lat with values that make each cell / layer distinguishable.
# Defined locally so this file does not depend on definitions in another one.
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

make_area_polygon <- function() {
  # Square covering the centre of the raster.
  sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(-5, -5), c(5, -5), c(5, 5), c(-5, 5), c(-5, -5)
    ))),
    crs = "EPSG:4326"
  )
}

test_that("extract_to_area() selects the correct depth range", {
  skip_if_not_installed("sf")
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())

  out <- extract_to_area(a, v, min_depth = 40, max_depth = 600)

  # depths 50, 100, 500 are the nearest standard layers bracketing 40-600
  expect_identical(terra::nlyr(out), 3)
  expect_true(all(grepl("depth=(50|100|500)$", names(out))))
})

test_that("extract_to_area() returns a SpatVoxel", {
  skip_if_not_installed("sf")
  out <- extract_to_area(sf::st_sf(geometry = make_area_polygon()),
                         make_area_voxel())
  expect_s4_class(out, "SpatVoxel")
  expect_true(methods::validObject(out))
})

test_that("extract_to_area() crops to area", {
  skip_if_not_installed("sf")
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())

  out <- extract_to_area(a, v, min_depth = 0, max_depth = 0)

  # Extent should shrink from -10..10 to roughly -5..5
  e <- terra::ext(out)
  expect_lt(e[2] - e[1], 20)
})

test_that("extract_to_area() keeps every layer when depth bounds are NULL", {
  skip_if_not_installed("sf")
  v <- make_area_voxel()
  out <- extract_to_area(sf::st_sf(geometry = make_area_polygon()), v)

  expect_identical(terra::nlyr(out), terra::nlyr(v))
  expect_identical(depths(out), depths(v))
})

test_that("extract_to_area() honours a single open depth bound", {
  skip_if_not_installed("sf")
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())

  # min only: from the nearest layer to 90 m down to the deepest layer
  expect_identical(depths(extract_to_area(a, v, min_depth = 90)), c(100, 500, 1000))
  # max only: from the shallowest layer down to the nearest layer to 90 m
  expect_identical(depths(extract_to_area(a, v, max_depth = 90)), c(0, 50, 100))
})

test_that("extract_to_area() accepts sf, sfc and SpatVector areas", {
  skip_if_not_installed("sf")
  v <- make_area_voxel()
  poly <- make_area_polygon()

  from_sfc <- extract_to_area(poly, v, min_depth = 0, max_depth = 0)
  from_sf <- extract_to_area(sf::st_sf(geometry = poly), v,
                             min_depth = 0, max_depth = 0)
  from_vect <- extract_to_area(terra::vect(poly), v,
                               min_depth = 0, max_depth = 0)

  expect_identical(as.vector(terra::ext(from_sfc)), as.vector(terra::ext(from_sf)))
  expect_identical(as.vector(terra::ext(from_sfc)), as.vector(terra::ext(from_vect)))
})

test_that("extract_to_area() reprojects an area in another CRS", {
  skip_if_not_installed("sf")
  v <- make_area_voxel()
  a <- sf::st_transform(make_area_polygon(), "EPSG:3857")

  out <- extract_to_area(a, v, min_depth = 0, max_depth = 0)

  expect_true(terra::same.crs(out, v))
  expect_lt(terra::ext(out)[2] - terra::ext(out)[1], 20)
})

test_that("extract_to_area() rejects a rast_3d that is not a SpatVoxel", {
  r <- terra::rast(nrows = 4, ncols = 4)
  terra::values(r) <- 1:16
  names(r) <- "not_a_depth_layer"
  # The voxel is checked before `area`, so a NULL area still surfaces this.
  expect_error(extract_to_area(NULL, r), "SpatVoxel")
})

test_that("extract_to_area() rejects an area that is not vector geometry", {
  expect_error(extract_to_area("nope", make_area_voxel()), "sf, sfc, or SpatVector")
})
