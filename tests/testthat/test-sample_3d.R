# sample_3d(): SpatVector3D lines sampled into points (#105).

dive <- function() {
  as_vector3d(
    data.frame(tag = "A", vessel = "V1", speed = c(1, 2, 3),
               lon = c(150, 150.1, 150.2), lat = -20, depth = c(0, 120, 15),
               time = as.POSIXct(c("2020-01-01 00:00", "2020-01-01 02:00",
                                   "2020-01-01 03:00"), tz = "UTC")),
    type = "lines", id = "tag", time = "time"
  )
}

test_that("vertices sampling keeps each line's vertices as points", {
  s <- sample_3d(dive())
  expect_s4_class(s, "SpatVector3D")
  expect_identical(geomtype(s), "points")
  v <- as.data.frame(s)
  expect_identical(v$feature_id, c("A_1", "A_2", "A_3"))
  expect_identical(v$track_id, rep("A", 3))
  expect_identical(v$sample_id, 1:3)
  expect_equal(v$depth, c(0, 120, 15))
  expect_equal(v$time, dive()@vertices$time)
  expect_false(any(v$is_interpolated))
  expect_equal(v$distance_km[1], 0)
  # A column with one value along the line is kept; one that changes isn't.
  expect_identical(v$vessel, rep("V1", 3))
  expect_false("speed" %in% names(v))
})

test_that("a spacing inserts points with depth and time along the segment", {
  v <- as.data.frame(sample_3d(dive(), spacing_km = 2))
  expect_true(all(diff(v$distance_km) <= 2 + 1e-9))
  expect_identical(v$sample_id, seq_len(nrow(v)))
  # The second vertex is at 120 m and 02:00; points before it are on the
  # first segment, so depth and time grow in proportion to distance.
  first <- v$distance_km <= v$distance_km[v$depth == 120]
  share <- v$distance_km[first] / max(v$distance_km[first])
  expect_equal(v$depth[first], 120 * share)
  expect_equal(as.numeric(v$time[first] - v$time[1], units = "hours"),
               2 * share)
  expect_identical(sum(!v$is_interpolated), 3L)
  # Support weights add up to the line's length and duration.
  expect_equal(sum(v$support_km), max(v$distance_km))
  expect_equal(sum(v$support_seconds), 3 * 3600)
})

test_that("depth and time steps limit the change between points", {
  v <- as.data.frame(sample_3d(dive(), depth_step_m = 10))
  expect_true(all(abs(diff(v$depth)) <= 10 + 1e-9))
  v <- as.data.frame(sample_3d(dive(), time_step = 30, time_unit = "minutes"))
  expect_true(all(diff(as.numeric(v$time)) <= 1800 + 1e-6))
})

test_that("monthly steps count from the start and don't drift (#146)", {
  route <- as_vector3d(
    data.frame(id = "C", lon = c(0, 2), lat = 0, depth = 5,
               date = as.Date(c("2020-01-31", "2020-06-30"))),
    type = "lines", id = "id", time = "date"
  )
  s <- sample_3d(route, time_step = 1, time_unit = "months")
  expect_identical(
    format(as.data.frame(s)$time, "%Y-%m-%d"),
    c("2020-01-31", "2020-02-29", "2020-03-31", "2020-04-30", "2020-05-31",
      "2020-06-30")
  )
  # Dates without hours stay dates: each point is on its UTC day.
  v <- as.data.frame(sample_3d(route, spacing_km = 50))
  expect_true(all(as.numeric(v$time) %% 86400 == 0))
  expect_identical(s@time_resolution, "daily")
})

test_that("a line across the antimeridian goes the short way", {
  line <- as_vector3d(
    data.frame(id = "B", lon = c(179.5, -179.5), lat = 0, depth = c(0, 10)),
    type = "lines", id = "id"
  )
  v <- as.data.frame(sample_3d(line, spacing_km = 30))
  expect_true(all(abs(v$lon) >= 179.5 - 1e-9))
  expect_lt(max(v$distance_km), 112)
  expect_equal(v$depth, seq(0, 10, length.out = nrow(v)))
  expect_true(all(is.na(v$time)))
  expect_true(all(is.na(v$support_seconds)))
})

test_that("several lines keep their own IDs, in order", {
  two <- as_vector3d(
    data.frame(id = rep(c("a", "b"), each = 2), lon = c(0, 1, 10, 11),
               lat = 0, depth = c(5, 15, 20, 30)),
    type = "lines", id = "id"
  )
  v <- as.data.frame(sample_3d(two, spacing_km = 50))
  expect_identical(unique(v$track_id), c("a", "b"))
  expect_false(anyDuplicated(v$feature_id) > 0)
  expect_equal(v$distance_km[v$track_id == "b"][1], 0)
})

test_that("the samples can be buffered into polygons", {
  b <- buffer_3d(sample_3d(dive(), spacing_km = 5), radius_km = 1,
                 depth_window = 10)
  expect_identical(geomtype(b), "polygons")
  expect_identical(length(b), length(sample_3d(dive(), spacing_km = 5)))
})

test_that("sample_3d() rejects input it can't sample", {
  pts <- as_vector3d(data.frame(lon = 150, lat = -20, depth = 10))
  untimed <- as_vector3d(
    data.frame(id = "B", lon = c(0, 1), lat = 0, depth = 5),
    type = "lines", id = "id"
  )
  expect_error(sample_3d(pts), "must be SpatVector3D lines",
               class = "ocean3d_error_vector3d")
  expect_error(sample_3d(untimed, time_step = 1), "needs a time at each",
               class = "ocean3d_error_vector3d")
  expect_error(sample_3d(dive(), sampling = "vertices", spacing_km = 1),
               "cannot be used with vertices")
  expect_error(sample_3d(dive(), sampling = "balanced"), "needs an explicit")
  expect_error(sample_3d(dive(), spacing_km = -1), "positive")
  expect_error(sample_3d(dive(), time_step = 1.5, time_unit = "months"),
               "whole months")
  expect_error(sample_3d(dive(), spacing_km = 0.01, max_samples = 10),
               "max_samples")
})
