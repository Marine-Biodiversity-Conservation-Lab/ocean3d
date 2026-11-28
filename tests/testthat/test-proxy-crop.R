# Cropping a SpatProxy (#113): crop() cuts a proxy down to an area without
# reading values, and as_voxel() then reads only that window. Synthetic files
# from make_grid_file() (helper-netcdf-grid.R).

grid_lon <- seq(0.5, 9.5, by = 1)
grid_lat <- seq(10.5, 14.5, by = 1)

test_that("a cropped proxy reads only the cells of its area", {
  f <- make_grid_file(lon = grid_lon, lat = grid_lat)
  on.exit(unlink(f), add = TRUE)
  p <- as_proxy(f, var = "temp")
  area <- terra::ext(2, 5, 11, 13)
  small <- terra::crop(p, area)

  expect_s4_class(small, "SpatProxy")
  expect_equal(small@axes$lon$centres, c(2.5, 3.5, 4.5))
  expect_equal(small@axes$lat$centres, c(11.5, 12.5))
  expect_identical(small@variable$dim_lengths[c("lon", "lat")],
                   c(lon = 3L, lat = 2L))
  # The same voxel as cropping the whole one.
  v <- as_voxel(small, time = 2)
  whole <- terra::crop(as_voxel(p, time = 2), area, snap = "out")
  expect_equal(terra::values(v), terra::values(whole))
  expect_equal(as.vector(terra::ext(v)), as.vector(terra::ext(whole)))
  expect_identical(depths(v), depths(whole))

  # A read of the cropped proxy is that part of the file, and cropping again
  # moves further into it.
  full <- .read_proxy_window(p, list(time = c(1, 1)))
  expect_identical(.read_proxy_window(small, list(time = c(1, 1))),
                   full[3:5, 2:3, , , drop = FALSE])
  smaller <- terra::crop(small, terra::ext(3, 4, 11, 12))
  expect_equal(smaller@axes$lon$centres, 3.5)
  expect_identical(.read_proxy_window(smaller, list(time = c(1, 1))),
                   full[4, 2, , , drop = FALSE])
})

test_that("crop() takes rasters, vectors and sf, in any CRS", {
  f <- make_grid_file(lon = grid_lon, lat = grid_lat)
  on.exit(unlink(f), add = TRUE)
  p <- as_proxy(f, var = "temp")
  want <- terra::crop(p, terra::ext(2, 5, 11, 13))@axes

  square <- terra::vect("POLYGON ((2 11, 5 11, 5 13, 2 13, 2 11))",
                        crs = "EPSG:4326")
  areas <- list(
    vector = square,
    sf = sf::st_as_sf(square),
    raster = terra::rast(square, ncols = 3, nrows = 2),
    projected = terra::project(square, "EPSG:3857")
  )
  for (kind in names(areas)) {
    got <- terra::crop(p, areas[[kind]])@axes
    expect_equal(got$lon$centres, want$lon$centres, label = kind)
    expect_equal(got$lat$centres, want$lat$centres, label = kind)
  }

  # In a polar projection, the top edge of a box bends north between its
  # corners: the corners alone reach 77 N, the edge 80.8 N.
  arctic <- make_grid_file(lon = seq(-89.5, -0.5, by = 1),
                           lat = seq(65.5, 84.5, by = 1), depth = 5, time = 0)
  on.exit(unlink(arctic), add = TRUE)
  polar <- terra::as.polygons(terra::ext(-1e6, 1e6, -2e6, -1e6),
                              crs = "EPSG:3413")
  rows <- terra::crop(as_proxy(arctic, var = "temp"), polar)@axes$lat$centres
  expect_equal(max(rows), 80.5)

  # A box around the North Pole: its edges reach only 80.8 N, but the area
  # includes every longitude up to the pole.
  north <- make_grid_file(lon = seq(-179.5, 179.5, by = 5),
                          lat = seq(70.5, 89.5, by = 1), depth = 5, time = 0)
  on.exit(unlink(north), add = TRUE)
  around <- terra::as.polygons(terra::ext(-1e6, 1e6, -1e6, 1e6),
                               crs = "EPSG:3413")
  cropped <- terra::crop(as_proxy(north, var = "temp"), around)@axes
  expect_equal(max(cropped$lat$centres), 89.5)
  expect_length(cropped$lon$centres, 72)

  # A polar box across the antimeridian, without the pole: about 174 to 186
  # degrees east. On a 0 to 360 grid it is one window; on a -180 to 180 grid
  # it is read from both ends of the grid, rather than taking every
  # longitude.
  across <- terra::as.polygons(terra::ext(-1.1e6, -0.9e6, 0.9e6, 1.1e6),
                               crs = "EPSG:3413")
  east <- make_grid_file(lon = seq(0.5, 359.5, by = 1),
                         lat = seq(70.5, 84.5, by = 1), depth = 5, time = 0)
  west <- make_grid_file(lon = seq(-179.5, 179.5, by = 1),
                         lat = seq(70.5, 84.5, by = 1), depth = 5, time = 0)
  on.exit(unlink(c(east, west)), add = TRUE)
  lon <- terra::crop(as_proxy(east, var = "temp"), across)@axes$lon$centres
  expect_equal(range(lon), c(174.5, 185.5))
  both_ends <- terra::crop(as_proxy(west, var = "temp"), across)
  expect_equal(both_ends@axes$lon$centres, lon)

  # A box covering every longitude (Web Mercator's full width) doesn't cross
  # the antimeridian: it keeps every longitude of a -180 to 180 grid.
  equator <- make_grid_file(lon = seq(-177.5, 177.5, by = 5),
                            lat = seq(-4.5, 4.5, by = 1), depth = 5, time = 0)
  on.exit(unlink(equator), add = TRUE)
  wide <- terra::as.polygons(
    terra::ext(-20037508, 20037508, -1e6, 1e6), crs = "EPSG:3857"
  )
  expect_length(terra::crop(as_proxy(equator, var = "temp"),
                            wide)@axes$lon$centres, 72)
})

test_that("crop() of a voxel in map units doesn't wrap longitudes", {
  # 100 m cells from 0 to 1000 m: x = 100 to 200 is one column, not also the
  # columns 360 m away.
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0,
                   ymax = 1000, crs = "EPSG:3857", vals = 1)
  v <- as_voxel(r, depths = 0)
  got <- terra::crop(as_proxy(v), terra::ext(100, 200, 0, 1000))
  expect_equal(got@axes$lon$centres, 150)
})

test_that("crop() matches longitudes across conventions", {
  # A 0 to 360 grid, and an area given from -180 to 180.
  east <- make_grid_file(lon = seq(180.5, 189.5, by = 1), lat = grid_lat)
  on.exit(unlink(east), add = TRUE)
  p <- as_proxy(east, var = "temp")
  expect_equal(terra::crop(p, terra::ext(-178, -175, 11, 13))@axes$lon$centres,
               c(182.5, 183.5, 184.5))

  cases <- list(
    list(x = p, area = terra::ext(50, 60, 11, 13),
         pattern = "doesn't overlap the grid's longitudes"),
    list(x = p, area = terra::ext(182, 184, 40, 50),
         pattern = "doesn't overlap the grid's latitudes")
  )
  for (case in cases) {
    err <- expect_error(terra::crop(case$x, case$area),
                        class = "ocean3d_error_proxy_crop")
    expect_match(conditionMessage(err), case$pattern)
  }
})

test_that("series and voxel sources crop too", {
  a <- make_grid_file(lon = grid_lon, lat = grid_lat, time = c(0, 1))
  b <- make_grid_file(lon = grid_lon, lat = grid_lat, time = c(2, 3))
  on.exit(unlink(c(a, b)), add = TRUE)
  area <- terra::ext(2, 5, 11, 13)

  # A series reads each file's window at the same offset.
  s <- as_proxy(c(a, b), var = "temp")
  expect_equal(terra::values(as_voxel(terra::crop(s, area), time = 3)),
               terra::values(terra::crop(as_voxel(s, time = 3), area,
                                         snap = "out")))

  # A voxel in memory is cropped as a raster.
  v <- as_voxel(as_proxy(a, var = "temp"), time = 1)
  cropped <- terra::crop(as_proxy(v), area)
  expect_identical(cropped@source$type, "voxel")
  expect_equal(terra::values(as_voxel(cropped)),
               terra::values(terra::crop(v, area, snap = "out")))
  # An area in another CRS is projected first, as for a file.
  square <- terra::vect("POLYGON ((2 11, 5 11, 5 13, 2 13, 2 11))",
                        crs = "EPSG:4326")
  projected <- terra::crop(as_proxy(v), terra::project(square, "EPSG:3857"))
  expect_equal(projected@axes$lon$centres, cropped@axes$lon$centres)
  expect_equal(projected@axes$lat$centres, cropped@axes$lat$centres)
})

test_that("crop() takes a SpatVector3D, so a buffer reads only its cells", {
  f <- make_grid_file(lon = grid_lon, lat = grid_lat)
  on.exit(unlink(f), add = TRUE)
  p <- as_proxy(f, var = "temp")
  detection <- as_vector3d(data.frame(lon = 4.5, lat = 12.5, depth = 10))
  # A circle of about 0.5 degrees radius reaches the neighbouring cells only.
  buffer <- buffer_3d(detection, radius_km = 60, depth_window = 5)
  cropped <- terra::crop(p, buffer)
  expect_equal(cropped@axes$lon$centres, c(3.5, 4.5, 5.5))
  expect_equal(cropped@axes$lat$centres, c(11.5, 12.5, 13.5))
  # The same cells as cropping with the buffer's extent.
  expect_equal(terra::values(as_voxel(cropped, time = 1)),
               terra::values(terra::crop(as_voxel(p, time = 1),
                                         terra::ext(buffer), snap = "out")))
  # A point is an area of no width: the cell it falls in, or both cells when
  # it sits on their shared edge (grid cells run 4 to 5 and 5 to 6 here).
  expect_equal(terra::crop(p, detection)@axes$lon$centres, 4.5)
  on_edge <- as_vector3d(data.frame(lon = 5, lat = 12.5, depth = 10))
  expect_equal(terra::crop(p, on_edge)@axes$lon$centres, c(4.5, 5.5))
  expect_equal(terra::crop(p, on_edge)@axes$lat$centres, 12.5)
  # Off the grid, a point takes no cell.
  expect_error(terra::crop(p, as_vector3d(data.frame(lon = 20, lat = 12.5,
                                                     depth = 10))),
               "doesn't overlap", class = "ocean3d_error_proxy_crop")
})

test_that("crop() reads across the antimeridian from both ends of the grid", {
  # The voxel's value at each cell, against the file's value at the cell's
  # own longitude, which may be 360 degrees away on the grid.
  expect_cells_from_file <- function(v, p) {
    stored <- .read_proxy_window(p, list(time = c(1, 1)))
    lon <- p@axes$lon$centres
    xy <- terra::xyFromCell(v, seq_len(terra::ncell(v)))
    i <- match(.normalise_longitude(xy[, 1]), .normalise_longitude(lon))
    j <- match(xy[, 2], p@axes$lat$centres)
    expect_false(anyNA(i))
    want <- vapply(seq_along(i), function(k) stored[i[k], j[k], 1L, 1L], 0)
    expect_equal(terra::values(v)[, 1], want)
  }
  detection <- as_vector3d(data.frame(lon = 179.8, lat = 12.5, depth = 10))
  buffer <- buffer_3d(detection, radius_km = 60, depth_window = 5)

  # A -180 to 180 grid: the cells at 179.5 and -179.5 (as 180.5), and no
  # others, with the voxel extent running on past 180.
  globe <- make_grid_file(lon = seq(-179.5, 179.5, by = 1), lat = grid_lat)
  on.exit(unlink(globe), add = TRUE)
  p <- as_proxy(globe, var = "temp")
  small <- terra::crop(p, buffer)
  expect_equal(small@axes$lon$centres, c(179.5, 180.5))
  v <- as_voxel(small, time = 1)
  expect_equal(as.vector(terra::ext(v))[c("xmin", "xmax")],
               c(xmin = 179, xmax = 181))
  expect_cells_from_file(v, p)
  # mask() puts the buffer on the voxel's grid past 180 too: it keeps cells
  # on both sides of the antimeridian.
  kept <- !is.na(terra::values(mask(v, buffer))[, 1])
  x <- terra::xFromCell(v, seq_len(terra::ncell(v)))
  expect_true(any(kept[x > 180]) && any(kept[x < 180]))
  # Cropping again stays on the joined axis.
  again <- terra::crop(small, terra::ext(179.6, 180.4, 12, 13))
  expect_equal(again@axes$lon$centres, c(179.5, 180.5))

  # A 0 to 360 grid reads the same cells as one window, and an area across
  # its own edge (0 degrees) from both ends.
  east <- make_grid_file(lon = seq(0.5, 359.5, by = 1), lat = grid_lat)
  on.exit(unlink(east), add = TRUE)
  q <- as_proxy(east, var = "temp")
  expect_equal(terra::crop(q, buffer)@axes$lon$centres, c(179.5, 180.5))
  near_zero <- terra::crop(q, terra::ext(-1.5, 1.5, 11, 13))
  expect_equal(near_zero@axes$lon$centres, c(-1.5, -0.5, 0.5, 1.5))
  expect_cells_from_file(as_voxel(near_zero, time = 1), q)

  # A series reads both runs from each file.
  a <- make_grid_file(lon = seq(-179.5, 179.5, by = 1), lat = grid_lat,
                      time = 0)
  b <- make_grid_file(lon = seq(-179.5, 179.5, by = 1), lat = grid_lat,
                      time = 1, values = 1e5 + seq_len(360 * 5 * 2))
  on.exit(unlink(c(a, b)), add = TRUE)
  s <- as_proxy(c(a, b), var = "temp")
  cut <- terra::crop(s, buffer)
  rows <- match(cut@axes$lat$centres, grid_lat)
  expect_identical(.read_proxy_window(cut),
                   .read_proxy_window(s)[c(360, 1), rows, , , drop = FALSE])

  # A voxel in memory is cropped on each side and joined.
  voxel <- as_voxel(p, time = 1)
  from_voxel <- as_voxel(terra::crop(as_proxy(voxel), buffer))
  expect_equal(terra::values(from_voxel), terra::values(v))
  expect_equal(as.vector(terra::ext(from_voxel)), as.vector(terra::ext(v)))
  expect_identical(depth_bounds(from_voxel), depth_bounds(v))
  expect_equal(terra::time(from_voxel), terra::time(v))
})

test_that("a WOA-style file reads through as_proxy(), crop(), as_voxel()", {
  # WOA's depth bounds and one-step climatology (`months since` year 0).
  f <- make_grid_file(lon = grid_lon, lat = grid_lat, depth = c(0, 5, 10),
                      depth_bounds = rbind(c(0, 2.5, 7.5), c(2.5, 7.5, 12.5)),
                      time = 6,
                      time_units = "months since 0000-01-01 00:00:00",
                      var = "t_an")
  on.exit(unlink(f), add = TRUE)
  v <- as_voxel(terra::crop(as_proxy(f, var = "t_an"),
                            terra::ext(2, 5, 11, 13)))
  expect_s4_class(v, "SpatVoxel")
  expect_identical(dim(v), c(2, 3, 3))
  expect_equal(depth_bounds(v)$lower, c(0, 2.5, 7.5))
  expect_equal(depth_bounds(v)$upper, c(2.5, 7.5, 12.5))
  expect_identical(unique(terra::varnames(v)), "t_an")
  expect_true(terra::same.crs(v, "EPSG:4326"))
  expect_false(terra::timeInfo(v)$time)
})

test_that("as_proxy() picks the `ma` field out of a WOA-style file", {
  # A monthly WOA file holds each field as its own variable; `t_ma` is the
  # month minus the annual climatology.
  n <- length(grid_lon) * length(grid_lat) * 3
  f <- make_grid_file(lon = grid_lon, lat = grid_lat, depth = c(0, 5, 10),
                      depth_bounds = rbind(c(0, 2.5, 7.5), c(2.5, 7.5, 12.5)),
                      time = 0.5,
                      time_units = "months since 0000-01-01 00:00:00",
                      var = "t_an", extra_vars = list(t_ma = -seq_len(n)))
  on.exit(unlink(f), add = TRUE)

  area <- terra::ext(2, 5, 11, 13)
  an <- as_voxel(terra::crop(as_proxy(f, var = "t_an"), area))
  v <- as_voxel(terra::crop(as_proxy(f, var = "t_ma"), area))
  expect_identical(unique(terra::varnames(v)), "t_ma")
  expect_identical(terra::nlyr(v), 3)
  expect_equal(terra::values(v), -terra::values(an), ignore_attr = TRUE)
  expect_equal(depth_bounds(v), depth_bounds(an))
})
