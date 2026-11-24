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
    list(time = 3, pattern = "step number from 1 to 2"),
    list(time = "2021-06-01", time_match = "exact",
         pattern = "No time step matches"),
    list(time = "2021-06-01", max_time_diff = 30,
         pattern = "No time step is within 30 days"),
    list(time = c("2020-01-01", "2020-01-02", "2020-01-03"),
         pattern = "one date or date-time, a window of two"),
    list(time = NA, pattern = "must not be NA"),
    list(time = c("2020-01-01", NA), pattern = "must not be NA"),
    list(time = 1.5, pattern = "one whole step number from 1 to 2"),
    list(time = c(1, 2), pattern = "one whole step number from 1 to 2"),
    list(time = c("2020-01-01", "2020-01-02 06:00"),
         pattern = "couldn't be read as dates or date-times")
  )
  for (case in errors) {
    expect_error(
      as_voxel(g$proxy, time = case$time,
               time_match = case$time_match %||% "nearest",
               max_time_diff = case$max_time_diff),
      case$pattern, class = "ocean3d_error_proxy_time"
    )
  }
  # Several steps on one day: a date is ambiguous, a date-time isn't.
  hourly <- proxy_and_file(time = c(0, 6, 12),
                           time_units = "hours since 2020-01-01")
  on.exit(unlink(hourly$file), add = TRUE)
  expect_error(as_voxel(hourly$proxy, time = as.Date("2020-01-01")),
               "Give a date-time", class = "ocean3d_error_proxy_time")
  expect_s4_class(as_voxel(hourly$proxy, time = "2020-01-01 06:00Z"),
                  "SpatVoxel")
  # Text without a time zone is read as UTC, with the extraction warning.
  expect_warning(as_voxel(hourly$proxy, time = "2020-01-01 06:00"),
                 class = "ocean3d_timezone_warning")
  expect_error(as_voxel(g$proxy, time = 1, max_time_diff = -1),
               class = "ocean3d_error_voxel_input")

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

# Time matching in as_voxel() (#115): a date picks the nearest or the exact
# step, a window or a polygon gives a list of voxels, one per step.
test_that("a date reads the nearest step, or the exact one on request", {
  # Steps at midnight on 1, 2 and 3 January.
  g <- proxy_and_file(time = c(0, 1, 2))
  on.exit(unlink(g$file), add = TRUE)
  step_of <- function(...) {
    v <- as_voxel(g$proxy, ...)
    expect_s4_class(v, "SpatVoxel")
    format(terra::time(v)[1L], "%Y-%m-%d")
  }
  # Halfway between two steps the earlier wins; before the first, the first.
  expect_identical(step_of(time = "2020-01-02 12:00Z"), "2020-01-02")
  expect_identical(step_of(time = "2020-01-02 13:00Z"), "2020-01-03")
  expect_identical(step_of(time = "2019-12-01"), "2020-01-01")
  expect_identical(step_of(time = "2020-01-03", time_match = "exact"),
                   "2020-01-03")
  expect_identical(step_of(time = "2020-01-02", max_time_diff = 0),
                   "2020-01-02")
  # Date-times with no step at their instant are exact misses.
  expect_error(
    as_voxel(g$proxy, time = "2020-01-02 06:00Z", time_match = "exact"),
    class = "ocean3d_error_proxy_time"
  )
})

test_that("a time window gives a list of voxels, one per step inside it", {
  g <- proxy_and_file(time = c(0, 1, 2))
  on.exit(unlink(g$file), add = TRUE)
  three <- as_voxel(g$proxy, time = c("2020-01-01", "2020-01-03"))
  expect_type(three, "list")
  expect_null(names(three))
  expect_length(three, 3L)
  for (i in 1:3) {
    expect_s4_class(three[[i]], "SpatVoxel")
    expect_voxel_matches_file(three[[i]], g$proxy, step = i)
  }
  expect_equal(vapply(three, function(v) terra::time(v)[1L], numeric(1)),
               as.numeric(as.POSIXct(c("2020-01-01", "2020-01-02",
                                       "2020-01-03"), tz = "UTC")))
  # A window inside one step is still a list. With dates, the end day counts
  # whole; with date-times, the ends are exact.
  one <- as_voxel(g$proxy, time = as.Date(c("2020-01-02", "2020-01-02")))
  expect_type(one, "list")
  expect_length(one, 1L)
  expect_identical(terra::values(one[[1L]]),
                   terra::values(as_voxel(g$proxy, time = 2)))
  at <- as.POSIXct(c("2020-01-01 12:00", "2020-01-03 00:00"), tz = "UTC")
  expect_length(as_voxel(g$proxy, time = at), 2L)

  errors <- list(
    list(time = as.POSIXct(c("2020-01-01 01:00", "2020-01-01 12:00"),
                           tz = "UTC"),
         pattern = "No time step falls in the window"),
    list(time = c("2020-01-03", "2020-01-01"),
         pattern = "start must not be after its end")
  )
  for (case in errors) {
    err <- expect_error(as_voxel(g$proxy, time = case$time),
                        class = "ocean3d_error_proxy_time")
    expect_match(conditionMessage(err), case$pattern)
  }

  # A window over a series reads a step from each file.
  a <- make_grid_file(time = c(0, 1))
  b <- make_grid_file(time = c(2, 3), values = 100 + seq_len(24))
  on.exit(unlink(c(a, b)), add = TRUE)
  series <- as_proxy(c(a, b), var = "temp")
  both <- as_voxel(series, time = c("2020-01-02", "2020-01-03"))
  expect_length(both, 2L)
  expect_identical(terra::values(both[[1L]]),
                   terra::values(as_voxel(as_proxy(a, var = "temp"),
                                          time = 2)))
  expect_identical(terra::values(both[[2L]]),
                   terra::values(as_voxel(as_proxy(b, var = "temp"),
                                          time = 1)))

  # A voxel in memory is its own only step.
  dated <- make_voxel(1, c(0, 10))
  terra::time(dated) <- rep(as.POSIXct("2020-01-02", tz = "UTC"), 2)
  covered <- as_voxel(as_proxy(dated), time = c("2020-01-01", "2020-01-05"))
  expect_type(covered, "list")
  expect_length(covered, 1L)
  expect_error(as_voxel(as_proxy(dated), time = c("2020-02-01", "2020-02-05")),
               class = "ocean3d_error_proxy_time")
})

test_that("a SpatVector3D feature reads its steps, always as a list", {
  g <- proxy_and_file(time = c(0, 1, 2))
  on.exit(unlink(g$file), add = TRUE)
  point <- as_vector3d(data.frame(lon = 1, lat = 11, depth = 5,
                                  date = as.Date("2020-01-03")),
                       time = "date")
  # A point is a list of one, so a loop over shapes gets one kind of result.
  v <- as_voxel(g$proxy, time = point)
  expect_type(v, "list")
  expect_length(v, 1L)
  expect_equal(terra::time(v[[1L]])[1L], as.POSIXct("2020-01-03", tz = "UTC"))
  # A buffer a day either side of the point covers two of the three steps.
  buffer <- buffer_3d(point, radius_km = 10, depth_window = 5,
                      time_window = 1)
  expect_length(as_voxel(g$proxy, time = buffer), 2L)
  # A polygon at a single date covers that day: a list of one.
  expect_length(as_voxel(g$proxy, time = buffer_3d(point, 10, 5)), 1L)

  # On a one-step source, the shape's time, or lack of one, is ignored:
  # dated observations read a World Ocean Atlas climatology.
  untimed <- as_vector3d(data.frame(lon = 1, lat = 11, depth = 5))
  one <- proxy_and_file(time = 0)
  flat <- proxy_and_file(time = NULL)
  clim <- proxy_and_file(time = 6,
                         time_units = "months since 0000-01-01 00:00:00")
  on.exit(unlink(c(one$file, flat$file, clim$file)), add = TRUE)
  dated <- make_voxel(1, c(0, 10))
  terra::time(dated) <- rep(as.POSIXct("2019-06-01", tz = "UTC"), 2)
  sources <- list(one$proxy, flat$proxy, clim$proxy, as_proxy(dated))
  for (p in sources) {
    for (shape in list(point, buffer, untimed)) {
      out <- as_voxel(p, time = shape)
      expect_type(out, "list")
      expect_length(out, 1L)
      expect_s4_class(out[[1L]], "SpatVoxel")
    }
  }
  # A plain date on a climatology is still an error.
  expect_error(as_voxel(clim$proxy, time = "2020-01-03"), "climatology",
               class = "ocean3d_error_proxy_time")
  expect_error(as_voxel(g$proxy, time = untimed), "shape has no time",
               class = "ocean3d_error_proxy_time")

  # One feature at a time, and no lines.
  two_points <- as_vector3d(data.frame(lon = c(1, 2), lat = 11, depth = 5))
  expect_error(as_voxel(g$proxy, time = two_points),
               "one feature; this one has 2",
               class = "ocean3d_error_proxy_time")
  line <- as_vector3d(data.frame(id = "a", lon = c(1, 2), lat = 11,
                                 depth = 5),
                      type = "lines", id = "id")
  expect_error(as_voxel(g$proxy, time = line),
               "needs polygons or points, but these are lines",
               class = "ocean3d_error_vector3d_type")

  # The time arguments belong to a proxy.
  expect_error(as_voxel(make_voxel(1, c(0, 10)), time_match = "exact"),
               class = "ocean3d_error_voxel_input")
})
