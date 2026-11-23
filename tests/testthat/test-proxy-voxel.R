# The voxel builder (#113): as_voxel() on a SpatProxy reads one time step
# into a SpatVoxel. Synthetic files from make_grid_file()
# (helper-netcdf-grid.R).

# The value a voxel holds at each cell centre and depth, against the stored
# array read in fixed lon, lat, depth order for time step `step`.
expect_voxel_matches_file <- function(v, p, step = 1L) {
  stored <- .read_proxy_window(p, list(time = c(step, step)))
  for (i in seq_along(p@axes$lon$centres)) {
    for (j in seq_along(p@axes$lat$centres)) {
      xy <- cbind(p@axes$lon$centres[i], p@axes$lat$centres[j])
      got <- unlist(terra::extract(v, xy))
      want <- stored[i, j, order(p@axes$depth$centres), 1L]
      expect_equal(unname(got), want)
    }
  }
}

test_that("a time step of a proxy becomes a voxel with the file's values", {
  layouts <- list(
    list(),
    list(dim_order = c("depth", "lon", "time", "lat")),
    list(lat = c(11.5, 10.5)),
    list(depth = c(15, 5))
  )
  for (args in layouts) {
    g <- do.call(proxy_and_file, args)
    on.exit(unlink(g$file), add = TRUE)
    v <- as_voxel(g$proxy, time = 2)
    expect_s4_class(v, "SpatVoxel")
    expect_identical(depths(v), c(5, 15))
    expect_voxel_matches_file(v, g$proxy, step = 2L)
  }

  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)
  v <- as_voxel(g$proxy, time = 1)
  # The extent comes from the cell edges, and the CRS is always set.
  expect_equal(as.vector(terra::ext(v)),
               c(xmin = 0, xmax = 3, ymin = 10, ymax = 12))
  expect_true(terra::same.crs(v, "EPSG:4326"))
  expect_identical(names(v), c("temp_depth=5", "temp_depth=15"))
  expect_identical(terra::units(v), c("degC", "degC"))
  expect_equal(terra::time(v), rep(as.POSIXct("2020-01-01", tz = "UTC"), 2))

  # Packed values are unpacked.
  packed <- proxy_and_file(prec = "short")
  on.exit(unlink(packed$file), add = TRUE)
  expect_equal(terra::values(as_voxel(packed$proxy, time = 1)),
               terra::values(v) * 0.1)
})

test_that("a voxel keeps the proxy's depth edges", {
  # WOA-style: explicit depth bounds.
  woa <- proxy_and_file(depth = c(0, 5, 10),
                        depth_bounds = rbind(c(0, 2.5, 7.5),
                                             c(2.5, 7.5, 12.5)))
  on.exit(unlink(woa$file), add = TRUE)
  b <- depth_bounds(as_voxel(woa$proxy, time = 1))
  expect_equal(b$lower, c(0, 2.5, 7.5))
  expect_equal(b$upper, c(2.5, 7.5, 12.5))

  # Copernicus-style: levels at cell centres without bounds. The voxel gets
  # the edges the proxy inferred, from the surface to half a step below the
  # deepest level, not stopping at the end levels.
  cmems <- proxy_and_file(depth = c(0.494, 1.541, 2.646))
  on.exit(unlink(cmems$file), add = TRUE)
  v <- as_voxel(cmems$proxy, time = 1)
  expect_equal(depth_bounds(v)$lower, cmems$proxy@axes$depth$lower)
  expect_equal(depth_bounds(v)$upper, cmems$proxy@axes$depth$upper)
  expect_equal(depth_bounds(v)$lower[1], 0)
  expect_gt(max(depth_bounds(v)$upper), 2.646)
  # The voxel has edges, so `bounds` doesn't change its volume, and a subset
  # measures its share of the whole column.
  expect_equal(volume(v, bounds = "top"), volume(v))
  expect_equal(volume(v[[1:2]]) + volume(v[[3]]), volume(v))

  # A single level without bounds has no edges to store, so `bounds` applies
  # (note on #113).
  single <- proxy_and_file(depth = 50)
  on.exit(unlink(single$file), add = TRUE)
  one <- as_voxel(single$proxy, time = 1)
  expect_identical(single$proxy@axes$depth$edges, "none")
  expect_null(depth_bounds(one))
  expect_identical(depths(one), 50)
})

test_that("`time` picks a step by number or date", {
  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)
  expect_identical(terra::values(as_voxel(g$proxy, time = "2020-01-02")),
                   terra::values(as_voxel(g$proxy, time = 2)))

  errors <- list(
    list(time = NULL, pattern = "has 2 time steps. Choose one"),
    list(time = 3, pattern = "a step from 1 to 2"),
    list(time = "2021-06-01", pattern = "No time step matches")
  )
  for (case in errors) {
    err <- expect_error(as_voxel(g$proxy, time = case$time),
                        class = "ocean3d_error_proxy_time")
    expect_match(conditionMessage(err), case$pattern)
  }

  # No time axis: the only step, and `time` is an error.
  flat <- proxy_and_file(time = NULL)
  on.exit(unlink(flat$file), add = TRUE)
  expect_s4_class(as_voxel(flat$proxy), "SpatVoxel")
  expect_error(as_voxel(flat$proxy, time = 1),
               class = "ocean3d_error_proxy_time")

  # A one-step climatology has no date to set.
  clim <- proxy_and_file(time = 6,
                         time_units = "months since 0000-01-01 00:00:00")
  on.exit(unlink(clim$file), add = TRUE)
  expect_false(terra::timeInfo(as_voxel(clim$proxy))$time)
})

test_that("a proxy that can't be a voxel stops", {
  no_depth <- proxy_and_file(depth = NULL)
  on.exit(unlink(no_depth$file), add = TRUE)
  expect_error(as_voxel(no_depth$proxy, time = 1),
               class = "ocean3d_error_proxy_no_depth")

  # Cells that differ in width by more than float32 rounding.
  bumpy <- proxy_and_file(lon = c(0.5, 1.5, 2.7))
  on.exit(unlink(bumpy$file), add = TRUE)
  err <- expect_error(as_voxel(bumpy$proxy, time = 1),
                      class = "ocean3d_error_irregular_grid")
  expect_match(conditionMessage(err), "longitude axis isn't regular")

  # float32 rounding is regular enough, including the drift that builds up
  # along a 1/12 degree axis (about 0.02% of a cell on Copernicus Marine).
  rounded <- proxy_and_file(lon = c(0.5, 1.5000001, 2.4999999))
  on.exit(unlink(rounded$file), add = TRUE)
  expect_s4_class(as_voxel(rounded$proxy, time = 1), "SpatVoxel")
  drift <- proxy_and_file(lon = c(1, 3, 5) / 24 + c(0, 0, 1.5e-5))
  on.exit(unlink(drift$file), add = TRUE)
  expect_s4_class(as_voxel(drift$proxy, time = 1), "SpatVoxel")

  # `depths` and `depth_bounds` come from the proxy; `time` needs one.
  g <- proxy_and_file()
  on.exit(unlink(g$file), add = TRUE)
  expect_error(as_voxel(g$proxy, depths = c(1, 2)),
               class = "ocean3d_error_voxel_input")
  expect_error(as_voxel(terra::rast(), time = 1),
               class = "ocean3d_error_voxel_input")
})
