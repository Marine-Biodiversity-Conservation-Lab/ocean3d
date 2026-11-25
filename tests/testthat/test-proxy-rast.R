# rast() on a SpatProxy (#116): time steps of a variable without depth read
# as plain one-layer SpatRasters. Synthetic files from make_grid_file()
# (helper-netcdf-grid.R); its values are 1, 2, 3, ... in lon, lat, depth,
# time storage order.

test_that("a 2D variable reads as one plain layer per step", {
  g <- proxy_and_file(depth = NULL, time = c(0, 1, 2))
  on.exit(unlink(g$file), add = TRUE)
  nc <- ncdf4::nc_open(g$file)
  stored <- ncdf4::ncvar_get(nc, "temp")
  ncdf4::nc_close(nc)

  r <- terra::rast(g$proxy, time = 2)
  expect_false(methods::is(r, "SpatVoxel"))
  expect_equal(terra::nlyr(r), 1)
  expect_identical(names(r), "temp")
  expect_equal(terra::time(r), as.POSIXct("2020-01-02", tz = "UTC"))
  # Each cell centre holds the stored value of step 2.
  for (i in 1:3) {
    for (j in 1:2) {
      xy <- cbind(c(0.5, 1.5, 2.5)[i], c(10.5, 11.5)[j])
      expect_equal(terra::extract(r, xy)[[1]], stored[i, j, 2])
    }
  }

  # A date picks the same step, and a window gives a list.
  expect_identical(terra::values(terra::rast(g$proxy, time = "2020-01-02")),
                   terra::values(r))
  w <- terra::rast(g$proxy, time = c("2020-01-01", "2020-01-02"))
  expect_length(w, 2L)
  expect_equal(do.call(c, lapply(w, terra::time)),
               as.POSIXct(c("2020-01-01", "2020-01-02"), tz = "UTC"))
})

test_that("a variable with depth is read with as_voxel(), not rast()", {
  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)
  expect_error(terra::rast(g$proxy, time = 1), "Read it with\\s+as_voxel\\(\\)",
               class = "ocean3d_error_proxy_has_depth")
  # A voxel in memory has depth too.
  v <- as_voxel(g$proxy, time = 1)
  expect_error(terra::rast(as_proxy(v)), "has a depth axis",
               class = "ocean3d_error_proxy_has_depth")
})

test_that("a cropped proxy reads only its cells", {
  g <- proxy_and_file(depth = NULL)
  on.exit(unlink(g$file), add = TRUE)
  cropped <- terra::rast(terra::crop(g$proxy, terra::ext(1, 2, 10, 11)),
                         time = 1)
  expect_equal(dim(cropped), c(1, 1, 1))
  whole <- terra::rast(g$proxy, time = 1)
  expect_equal(unname(terra::values(cropped)[1, 1]),
               terra::extract(whole, cbind(1.5, 10.5))[[1]])
})

test_that("a 2D step works with mask() and summarise_3d()", {
  g <- proxy_and_file(depth = NULL)
  on.exit(unlink(g$file), add = TRUE)
  r <- terra::rast(g$proxy, time = 1)

  # A point keeps its own cell; its depth doesn't apply to a 2D layer.
  pt <- as_vector3d(data.frame(lon = 1.5, lat = 10.5, depth = 300))
  kept <- terra::values(mask(r, pt))[, 1]
  expect_equal(which(!is.na(kept)), terra::cellFromXY(r, cbind(1.5, 10.5)))
  # A 2D layer has one value per cell, so there is no depth level to choose.
  expect_error(mask(r, pt, method = "bottom"),
               "`method` applies only when mask\\(\\) keeps a SpatVoxel",
               class = "ocean3d_error_unused_argument")

  # A polygon keeps the cells it touches, with no message about depth.
  square <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0.2, 10.2), c(1.8, 10.2), c(1.8, 10.8), c(0.2, 10.8), c(0.2, 10.2)
  ))), crs = 4326))
  poly <- as_vector3d(square, depth_min = 0, depth_max = 50)
  masked <- expect_no_message(mask(r, poly))
  expect_identical(sum(!is.na(terra::values(masked))), 2L)

  # Weighted by cell area, which shrinks with latitude.
  area <- terra::values(terra::cellSize(r, unit = "km"))[, 1]
  vals <- terra::values(r)[, 1]
  out <- summarise_3d(r, statistic = c("mean", "min"))
  expect_equal(out$mean, sum(vals * area) / sum(area))
  expect_identical(out$n, 6L)
  expect_equal(out$time, as.POSIXct("2020-01-01", tz = "UTC"))
  expect_equal(summarise_3d(r, weighting = "equal")$mean, mean(vals))
  expect_identical(nrow(summarise_3d(list(r, r + 1))), 2L)

  expect_error(summarise_3d(r, by = "depth"), "has a single level",
               class = "ocean3d_error_summarise_input")
})

test_that("rast() and as_voxel() say how to read each kind of variable", {
  g <- proxy_and_file(depth = NULL)
  on.exit(unlink(g$file), add = TRUE)
  expect_error(as_voxel(g$proxy, time = 1), "terra::rast\\(\\) on the proxy",
               class = "ocean3d_error_proxy_no_depth")
  expect_error(terra::rast(g$proxy, time = 1, depths = 1), "takes only",
               class = "ocean3d_error_voxel_input")
  expect_error(terra::rast(g$proxy, time = 1, max_time_diff = -1),
               "`max_time_diff` must be one number",
               class = "ocean3d_error_voxel_input")
  expect_error(terra::rast(g$proxy), "has 2 time steps",
               class = "ocean3d_error_proxy_time")
})
