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

test_that("calendar shifts clamp to the end of a shorter month", {
  shift <- function(x, months) {
    format(.shift_temporal_months(as.POSIXct(x, tz = "UTC"), months),
           "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }
  ends <- c("2020-01-31 12:34:56", "2020-03-31 12:34:56",
            "2021-03-31 12:34:56")
  expect_equal(vapply(ends, shift, "", months = 1), c(
    "2020-02-29 12:34:56", "2020-04-30 12:34:56", "2021-04-30 12:34:56"
  ), ignore_attr = TRUE)
  expect_equal(vapply(ends, shift, "", months = -1), c(
    "2019-12-31 12:34:56", "2020-02-29 12:34:56", "2021-02-28 12:34:56"
  ), ignore_attr = TRUE)
  # Years are 12 months, so a leap day moves to 28 February.
  expect_equal(shift("2020-02-29 19:30:00", 12), "2021-02-28 19:30:00")
  expect_equal(shift("2020-02-29 19:30:00", -12), "2019-02-28 19:30:00")
})

test_that("balanced sampling meets every step limit with one even spacing", {
  # 2 degrees along the equator (222 km), 150 m down and 24 hours: the
  # hourly step needs the most points, 24 segments, and so sets the spacing
  # for distance and depth too.
  line <- as_vector3d(
    data.frame(id = "A", lon = c(0, 2), lat = 0, depth = c(30, 180),
               time = as.POSIXct(c("2020-01-01 00:00", "2020-01-02 00:00"),
                                 tz = "UTC")),
    type = "lines", id = "id", time = "time"
  )
  v <- as.data.frame(sample_3d(line, spacing_km = 30, depth_step_m = 10,
                               time_step = 1))
  expect_equal(nrow(v), 25)
  expect_equal(v$depth, seq(30, 180, length.out = 25))
  expect_equal(v$lon, seq(0, 2, length.out = 25))
  expect_equal(as.numeric(v$time - v$time[1], units = "hours"), 0:24)
  # The same limits with merged_steps: the 15 depth steps and 24 hourly
  # steps are kept as separate cuts.
  merged <- as.data.frame(sample_3d(line, sampling = "merged_steps",
                                    spacing_km = 30, depth_step_m = 10,
                                    time_step = 1))
  expect_equal(nrow(merged), 37)
  expect_lte(max(diff(merged$depth)), 10)
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

test_that("merged_steps keeps distance points and time-step points", {
  # Two hours over one 5.2 km segment: a 2 km spacing alone gives points
  # every 40 minutes, and hourly steps add a point at 01:00.
  line <- as_vector3d(
    data.frame(id = "A", lon = c(150, 150.05), lat = -20, depth = 0,
               time = as.POSIXct(c("2020-01-01 00:00", "2020-01-01 02:00"),
                                 tz = "UTC")),
    type = "lines", id = "id", time = "time"
  )
  even <- as.data.frame(sample_3d(line, spacing_km = 2))
  merged <- as.data.frame(sample_3d(line, sampling = "merged_steps",
                                    spacing_km = 2, time_step = 1))
  hour <- as.POSIXct("2020-01-01 01:00", tz = "UTC")
  expect_false(any(abs(as.numeric(even$time - hour)) < 1e-6))
  expect_true(any(abs(as.numeric(merged$time - hour)) < 1e-6))
  expect_true(all(round(even$distance_km, 9) %in%
                    round(merged$distance_km, 9)))
})

test_that("sample_3d() rejects input it can't sample", {
  pts <- as_vector3d(data.frame(lon = 150, lat = -20, depth = 10))
  untimed <- as_vector3d(
    data.frame(id = "B", lon = c(0, 1), lat = 0, depth = 5),
    type = "lines", id = "id"
  )
  partly <- dive()
  partly@vertices$time[2] <- NA
  # One timed line and one without times.
  mixed <- as_vector3d(
    data.frame(id = rep(c("a", "b"), each = 2), lon = c(0, 1, 5, 6), lat = 0,
               depth = 5,
               time = as.POSIXct(c("2020-01-01", "2020-01-02", NA, NA),
                                 tz = "UTC")),
    type = "lines", id = "id", time = "time"
  )
  daily <- as_vector3d(
    data.frame(id = "D", lon = c(0, 1), lat = 0, depth = 5,
               date = as.Date(c("2020-01-01", "2020-01-05"))),
    type = "lines", id = "id", time = "date"
  )
  clashing <- as_vector3d(
    data.frame(seg = "S", track_id = "T1", lon = c(0, 1), lat = 0, depth = 5),
    type = "lines", id = "seg"
  )
  cases <- list(
    list(quote(sample_3d(data.frame(lon = 1))), "must be SpatVector3D lines",
         "ocean3d_error_vector3d"),
    list(quote(sample_3d(pts)), "Points are already samples",
         "ocean3d_error_vector3d_type"),
    list(quote(sample_3d(untimed, time_step = 1)),
         "`time_step` needs a time at every vertex", "ocean3d_error_vector3d"),
    list(quote(sample_3d(mixed, time_step = 1)),
         "Feature b of `x` has no time", "ocean3d_error_vector3d"),
    list(quote(sample_3d(partly)),
         "Line A of `x` has a time at some vertices", "ocean3d_error_vector3d"),
    list(quote(sample_3d(daily, time_step = 6)), "must be whole days",
         "ocean3d_error_vector3d"),
    list(quote(sample_3d(clashing)),
         "columns track_id, which sample_3d() adds", "ocean3d_error_vector3d"),
    list(quote(sample_3d(dive(), sampling = "vertices", spacing_km = 1)),
         "cannot be used with vertices", "ocean3d_error_track_sampling"),
    list(quote(sample_3d(dive(), sampling = "balanced")), "needs an explicit",
         "ocean3d_error_track_sampling"),
    list(quote(sample_3d(dive(), spacing_km = -1)),
         "`spacing_km` must be one finite positive number",
         "ocean3d_error_track_sampling"),
    list(quote(sample_3d(dive(), time_step = 1.5, time_unit = "months")),
         "whole months", "ocean3d_error_track_sampling"),
    list(quote(sample_3d(dive(), max_samples = 10.5)), "positive integer",
         "ocean3d_error_track_sampling")
  )
  for (case in cases) {
    expect_error(eval(case[[1]]), case[[2]], class = case[[3]], fixed = TRUE,
                 label = deparse(case[[1]]))
  }
  # A whole number of days on dates is fine.
  expect_s4_class(sample_3d(daily, time_step = 1, time_unit = "days"),
                  "SpatVector3D")
  expect_error(sample_3d(dive(), spacing_km = 0.01, max_samples = 10),
               "Track sampling exceeds max_samples",
               class = "ocean3d_error_track_sampling")
})
