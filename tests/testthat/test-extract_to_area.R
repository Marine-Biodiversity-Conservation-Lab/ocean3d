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

test_that("extract_to_area() selects the correct depth range", {
  v <- make_area_voxel()
  a <- sf::st_sf(geometry = make_area_polygon())

  # `depths` are the layers kept from v's 0, 50, 100, 500 and 1000 m.
  cases <- list(
    # 50, 100 and 500 m are the nearest standard layers bracketing 40-600 m.
    list(label = "40 to 600 m", min = 40, max = 600, depths = c(50, 100, 500)),
    list(label = "no bounds", min = NULL, max = NULL,
         depths = c(0, 50, 100, 500, 1000)),
    # From the nearest layer to 90 m down to the deepest layer.
    list(label = "min only", min = 90, max = NULL, depths = c(100, 500, 1000)),
    # From the shallowest layer down to the nearest layer to 90 m.
    list(label = "max only", min = NULL, max = 90, depths = c(0, 50, 100))
  )
  for (case in cases) {
    out <- extract_to_area(a, v, min_depth = case$min, max_depth = case$max)

    expect_identical(terra::nlyr(out), as.numeric(length(case$depths)),
                     label = case$label)
    # The crop to the polygon keeps the depths in terra::depth().
    expect_identical(depths(out), case$depths, label = case$label)
    expect_identical(terra::depth(out), case$depths, label = case$label)
    expect_identical(names(out), paste0("tan_depth=", case$depths),
                     label = case$label)
  }
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
  names(r) <- "not_a_depth_layer"
  cases <- list(
    # The voxel is checked before `area`, so a NULL area still surfaces this.
    list(label = "rast_3d that is not a SpatVoxel", area = NULL, rast_3d = r,
         pattern = "SpatVoxel"),
    list(label = "area that is not vector geometry", area = "nope",
         rast_3d = make_area_voxel(), pattern = "sf, sfc, or SpatVector")
  )
  for (case in cases) {
    expect_error(extract_to_area(case$area, case$rast_3d), case$pattern,
                 label = case$label)
  }
})
