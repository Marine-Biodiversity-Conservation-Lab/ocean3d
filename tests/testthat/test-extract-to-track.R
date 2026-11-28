# extract_to_track() is a wrapper over as_vector3d(), sample_3d(), extract(),
# buffer_3d(), crop(), as_voxel() or rast(), mask() and summarise_3d() (#116),
# whose own tests cover their rules. These tests cover what the wrapper adds:
# how the vertex table becomes tracks, the value each sample takes, the
# summaries per track and their weights, and the output.

# A 10 x 10 grid of 1 degree cells at depths 5, 15, 30 and 60 m on three
# daily steps (1 to 3 January 2020). Each value is its position in storage
# order, lon fastest, then lat, depth and time, so it can be worked out by
# hand: cell (i, j), level k, step t holds i + 10 (j - 1) + 100 (k - 1) +
# 400 (t - 1).
track_file <- function(depth = c(5, 15, 30, 60), time = c(0, 1, 2), ...) {
  make_grid_file(lon = seq(0.5, 9.5, 1), lat = seq(0.5, 9.5, 1),
                 depth = depth, time = time, ...)
}
at <- function(i, j, k, t, levels = 4L) {
  i + 10 * (j - 1) + 10 * 10 * (k - 1) + 10 * 10 * levels * (t - 1)
}

test_that("each vertex takes the value in its cell, level and time step", {
  path <- track_file()
  on.exit(unlink(path))
  # Cell (3, 4) at 20 m on 1 January, then cell (6, 4) at 50 m on 3 January.
  dive <- data.frame(lon = c(2.2, 5.2), lat = 3.3, depth = c(20, 50),
                     date = as.Date(c("2020-01-01", "2020-01-03")))
  out <- extract_to_track(dive, path)
  expect_named(out, c("samples", "summary"))
  expect_named(out$samples, c(
    "track_id", "sample_id", "lon", "lat", "depth", "date", "distance_km",
    "support_km", "support_seconds", "is_interpolated", "temp"
  ))
  expect_identical(out$samples$track_id, c("1", "1"))
  expect_false(any(out$samples$is_interpolated))
  # 20 m is nearest 15 m (level 2), and 50 m nearest 60 m (level 4).
  values <- c(at(3, 4, 2, 1), at(6, 4, 4, 3))
  expect_equal(out$samples$temp, values)
  # Two vertices stand for half the track each.
  expect_equal(out$summary,
               data.frame(track_id = "1", mean = mean(values), n = 2L))

  rules <- list(nearest_level = c(2, 4), bottom = c(4, 4), surface = c(1, 1))
  for (rule in names(rules)) {
    k <- rules[[rule]]
    expect_equal(extract_to_track(dive, path, method = rule)$samples$temp,
                 c(at(3, 4, k[1], 1), at(6, 4, k[2], 3)), label = rule)
  }
  # The bottom and the surface don't need a depth.
  flat <- dive[c("lon", "lat", "date")]
  expect_equal(extract_to_track(flat, path, method = "bottom")$samples$temp,
               c(at(3, 4, 4, 1), at(6, 4, 4, 3)))
})

test_that("track summaries weight each sample by distance, time or equally", {
  path <- track_file()
  on.exit(unlink(path))
  # Three vertices along 0.5 N, at the centres of cells 1, 2 and 5, one hour
  # apart: the first segment is a third as long as the second.
  route <- data.frame(lon = c(0.5, 1.5, 4.5), lat = 0.5, depth = 5,
                      date = as.POSIXct("2020-01-01", tz = "UTC") +
                        c(0, 3600, 7200))
  values <- c(1, 2, 5)
  run <- function(...) {
    extract_to_track(route, path, statistic = c("mean", "sd", "min", "max"),
                     ...)$summary
  }
  # Half of each segment either side: 1, 4 and 3 eighths of the length.
  distance <- run()
  expect_equal(distance$mean, sum(values * c(1, 4, 3)) / 8, tolerance = 1e-4)
  w <- c(1, 4, 3) / 8
  expect_equal(distance$sd,
               sqrt(sum(w * (values - sum(w * values))^2) / (1 - sum(w^2))),
               tolerance = 1e-4)
  expect_equal(c(distance$min, distance$max), c(1, 5))
  expect_identical(distance$n, 3L)
  # Equal halves of the two hours: 1, 2 and 1 quarters.
  expect_equal(run(weighting = "duration")$mean, (1 + 4 + 5) / 4)
  equal <- run(weighting = "equal")
  expect_equal(c(equal$mean, equal$sd), c(mean(values), stats::sd(values)))

  # A track that doesn't move has no length, so its samples count equally.
  still <- data.frame(lon = 5.2, lat = 5.3, depth = 20,
                      date = as.Date(c("2020-01-01", "2020-01-03")))
  expect_equal(extract_to_track(still, path)$summary$mean,
               mean(c(at(6, 6, 2, 1), at(6, 6, 2, 3))))
})

test_that("samples between vertices take the values the steps give", {
  path <- track_file()
  on.exit(unlink(path))
  dive <- data.frame(lon = c(0.6, 4.4), lat = c(0.7, 3.2), depth = c(5, 55),
                     date = as.POSIXct(c("2020-01-01 00:00",
                                         "2020-01-02 18:00"), tz = "UTC"))
  out <- extract_to_track(dive, path, spacing_km = 40, depth_step_m = 10)
  # The same steps written out.
  line <- as_vector3d(cbind(id = "1", dive), type = "lines", id = "id",
                      time = "date")
  samples <- sample_3d(line, spacing_km = 40, depth_step_m = 10)
  found <- terra::extract(as_proxy(path), samples)
  expect_equal(out$samples$temp, found$value)
  expect_equal(out$samples$depth, as.data.frame(samples)$depth)
  expect_gt(nrow(out$samples), 2)
  weights <- as.data.frame(samples)$support_km
  expect_equal(out$summary$mean, sum(found$value * weights) / sum(weights))
})

test_that("a radius summarises the values in each sample's buffer", {
  path <- track_file()
  on.exit(unlink(path))
  # Two vertices in cell (6, 6) at 20 m on 2 January. 10 to 30 m overlaps
  # the slabs of 15 m (10 to 22.5) and 30 m (22.5 to 45).
  pair <- data.frame(lon = c(5.2, 5.4), lat = 5.3, depth = 20,
                     date = as.Date("2020-01-02"))
  out <- extract_to_track(pair, path, horizontal_radius_km = 10,
                          depth_radius_m = 10)
  expected <- sum(at(6, 6, 2:3, 2) * c(12.5, 22.5)) / 35
  expect_equal(out$samples$temp, rep(expected, 2))
  expect_identical(out$samples$temp_n, c(2L, 2L))
  expect_equal(out$summary$mean, expected)

  # A variable without depth: the cell the circle touches.
  flat <- make_grid_file(lon = seq(0.5, 9.5, 1), lat = seq(0.5, 9.5, 1),
                         depth = NULL, time = c(0, 1, 2))
  on.exit(unlink(flat), add = TRUE)
  out <- extract_to_track(pair[c("lon", "lat", "date")], flat,
                          horizontal_radius_km = 10)
  expect_equal(out$samples$temp, rep(6 + 50 + 100, 2))
  expect_identical(out$samples$temp_n, c(1L, 1L))
})

test_that("samples without a value get NA, counted in one warning", {
  path <- track_file()
  on.exit(unlink(path))
  # Track a ends off the grid; track b is a month after the last step.
  route <- data.frame(
    trip = c("a", "a", "b", "b"), lon = c(5.2, 12, 5.2, 5.4), lat = 5.3,
    depth = 20,
    date = as.Date(c("2020-01-02", "2020-01-02", "2020-02-02", "2020-02-02"))
  )
  for (radius in c(0, 10)) {
    expect_warning(
      out <- extract_to_track(route, path, track_col = "trip",
                              max_time_diff = 1,
                              horizontal_radius_km = radius,
                              depth_radius_m = if (radius) 10 else 0),
      paste("3 of 4 samples get NA: 2 with no time step that matches and 1",
            "outside the grid or without a value"),
      class = "ocean3d_warning_track_samples", fixed = TRUE
    )
    expect_identical(is.na(out$samples$temp), c(FALSE, TRUE, TRUE, TRUE))
    expect_identical(out$summary$n, c(1L, 0L))
    expect_true(is.na(out$summary$mean[2]))
  }
})

test_that("tracks keep their IDs and order, whatever form the vertices take", {
  path <- track_file()
  on.exit(unlink(path))
  route <- data.frame(trip = c(2, 2, 1, 1), lon = c(0.6, 1.6, 2.6, 3.6),
                      lat = 0.5, depth = 5)
  undated <- track_file(time = NULL)
  on.exit(unlink(undated), add = TRUE)
  expected <- extract_to_track(route, undated, track_col = "trip")
  expect_identical(expected$summary$track_id, c("2", "1"))
  expect_equal(expected$summary$mean, c(1.5, 3.5))
  expect_equal(expected$samples$temp, 1:4)

  points <- sf::st_as_sf(route, coords = c("lon", "lat"), crs = 4326)
  inputs <- list(matrix = as.matrix(route), list = as.list(route),
                 sf = points, proxy = NULL)
  for (form in names(inputs)) {
    out <- if (form == "proxy") {
      extract_to_track(route, as_proxy(undated), track_col = "trip")
    } else {
      extract_to_track(inputs[[form]], undated, track_col = "trip")
    }
    expect_equal(out$summary, expected$summary, label = form)
    expect_equal(out$samples$temp, 1:4, label = form)
  }
  # A dated file needs the dates; one without them doesn't.
  expect_error(extract_to_track(route, path, track_col = "trip"),
               "Missing required column(s) in `data`: date",
               class = "ocean3d_error_missing_columns", fixed = TRUE)
  # Lines are given as their vertices.
  line <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(cbind(c(0.6, 1.6), 0.5)), crs = 4326
  ))
  expect_error(extract_to_track(line, undated), "Convert lines",
               class = "ocean3d_error_points_input")
})

test_that("dates come back as dates, and undated samples are counted", {
  path <- track_file()
  on.exit(unlink(path))
  dive <- data.frame(trip = c("a", "a", "b", "b"), lon = c(2.2, 5.2, 2.2, 5.2),
                     lat = 3.3, depth = 20,
                     date = as.Date(c("2020-01-01", "2020-01-03", NA, NA)))
  expect_warning(
    out <- extract_to_track(dive, path, track_col = "trip"),
    "2 of 4 samples get NA: 2 without a date.",
    class = "ocean3d_warning_track_samples", fixed = TRUE
  )
  expect_s3_class(out$samples$date, "Date")
  expect_equal(out$samples$date[1:2], dive$date[1:2])
  expect_identical(out$summary$n, c(2L, 0L))
})

test_that("a buffer needs cells with a width", {
  # One longitude without bounds, as in a transect: the value at each
  # sample can be read, but a buffer can't be placed.
  transect <- make_grid_file(lon = 5, lat = seq(0.5, 9.5), depth = c(5, 15),
                             time = NULL)
  on.exit(unlink(transect))
  along <- data.frame(lon = 5, lat = c(2.2, 6.3), depth = 10)
  # Cells 3 and 7 along the latitudes; 10 m is halfway between the levels,
  # so the shallower one.
  expect_equal(extract_to_track(along, transect)$samples$temp, c(3, 7))
  expect_error(
    extract_to_track(along, transect, horizontal_radius_km = 20,
                     depth_radius_m = 5),
    "single longitude without cell bounds", class = "ocean3d_error_track_input"
  )
})

test_that("extract_to_track() rejects input it can't use", {
  path <- track_file()
  flat <- track_file(depth = NULL)
  undated <- track_file(time = NULL)
  on.exit(unlink(c(path, flat, undated)))
  dive <- data.frame(lon = c(2.2, 5.2), lat = 3.3, depth = c(20, 50),
                     date = as.Date(c("2020-01-01", "2020-01-03")))
  cases <- list(
    list(quote(extract_to_track(dive, path, horizontal_radius_km = -1)),
         "must each be one number, 0 or more"),
    list(quote(extract_to_track(dive, path, depth_radius_m = 10)),
         "`depth_radius_m` needs a `horizontal_radius_km` above 0"),
    list(quote(extract_to_track(dive, path, horizontal_radius_km = 10)),
         "buffers need a `depth_radius_m` above 0"),
    list(quote(extract_to_track(dive, path, horizontal_radius_km = 10,
                                depth_radius_m = 5, method = "bottom")),
         "Leave `method` NULL with a radius"),
    list(quote(extract_to_track(dive, path, method = "all")),
         "`method` must be one of nearest, nearest_level, bottom, surface"),
    list(quote(extract_to_track(dive, flat, method = "surface")),
         "has no depth, so `method` doesn't apply"),
    list(quote(extract_to_track(dive, path, output_col = "")),
         "must be one non-empty name"),
    list(quote(extract_to_track(dive, path, output_col = "distance_km")),
         "already have a column named distance_km"),
    list(quote(extract_to_track(stats::setNames(dive, c("x_n", "lat", "depth",
                                                        "date")),
                                path, output_col = "x", lon_col = "x_n",
                                horizontal_radius_km = 10,
                                depth_radius_m = 5)),
         "already have a column named x_n")
  )
  for (case in cases) {
    expect_error(eval(case[[1]]), case[[2]],
                 class = "ocean3d_error_track_input", fixed = TRUE,
                 label = deparse(case[[1]]))
  }
  untimed <- dive
  untimed$date <- NA
  expect_error(
    extract_to_track(untimed, undated, weighting = "duration"),
    "needs a date at every vertex", class = "ocean3d_error_track_input"
  )
  for (case in list(list(dive["lon"], "lat, depth, date"),
                    list(dive[c("lon", "lat", "depth")], "date"))) {
    expect_error(
      extract_to_track(case[[1]], if (ncol(case[[1]]) == 1) path else undated,
                       weighting = if (ncol(case[[1]]) == 1) "distance" else
                         "duration"),
      paste("Missing required column(s) in `data`:", case[[2]]),
      class = "ocean3d_error_missing_columns", fixed = TRUE
    )
  }
  # The build and sampling steps check their own input.
  expect_error(extract_to_track(dive[1, ], path),
               class = "ocean3d_error_vector3d")
  expect_error(extract_to_track(dive, path, sampling = "cells"),
               "should be one of")
  expect_error(extract_to_track(dive, path, statistic = "se"),
               "should be one of")
  expect_error(extract_to_track(dive, path, quiet = TRUE), "unused argument")
})
