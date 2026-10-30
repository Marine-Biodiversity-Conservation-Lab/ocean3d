test_that("buffered target locations follow the selected record", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path)
  route <- data.frame(lon = c(0, 2), lat = 0, target = 10)
  out <- extract_to_buffered_track(
    route, path, spacing_km = 112, horizontal_radius_km = 120,
    statistic = "most_similar", target_col = "target", infer_bounds = TRUE,
    quiet = TRUE
  )
  expect_equal(out$samples$temp, c(10, 10, 10))
  expect_equal(out$samples$temp_match_lon, c(1, 1, 1))
  expect_equal(out$samples$temp_match_lat, c(0, 0, 0))
  expect_equal(out$samples$lon, c(0, 1, 2))
  expect_true(all(is.na(out$samples$temp_match_depth)))
  expect_true(all(is.na(out$samples$temp_match_time)))
  expect_error(extract_to_buffered_track(
    transform(route, temp_match_lon = 9), path, horizontal_radius_km = 1,
    statistic = "most_similar", target_col = "target", infer_bounds = TRUE
  ), "overwrite")
})


test_that("buffered warnings exclude invalid connections and clear locations", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path)
  handle <- ncdf4::nc_open(path, write = TRUE)
  ncdf4::ncvar_put(handle, "temp", matrix(10, 3, 2))
  ncdf4::nc_close(handle)
  route <- data.frame(lon = c(0, 2), lat = 0, target = 10)
  island <- .track_test_polygon(0.9, 1.1, -0.2, 0.2)
  domain <- .track_test_polygon(-1, 3, -1, 2)
  captured <- list()
  out <- withCallingHandlers(extract_to_buffered_track(
    route, path, spacing_km = 112, horizontal_radius_km = 120,
    statistic = "most_similar", target_col = "target", infer_bounds = TRUE,
    land_mask = island, mask_domain = domain, quiet = TRUE
  ), ocean3d_target_location_warning = function(warning) {
    captured[[length(captured) + 1L]] <<- warning
    invokeRestart("muffleWarning")
  })
  expect_length(captured, 1)
  expect_identical(captured[[1]]$rows, c(1L, 3L))
  expect_equal(out$samples$temp, c(10, NA, 10))
  expect_true(is.na(out$samples$temp_match_lon[2]))
  expect_true(is.na(out$samples$temp_match_lat[2]))
  expect_true(is.na(out$samples$temp_match_depth[2]))
  expect_true(is.na(out$samples$temp_match_time[2]))
  expect_identical(out$samples$temp_match_info[[2]]$status,
                   "invalid_connection")
})


test_that("balanced sampling shares one subdivision for all limits", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path)
  route <- data.frame(
    lon = c(0, 2), lat = 0, depth = c(30, 180),
    date = .track_clock(c(0, 24))
  )
  arguments <- list(
    data = route, nc = path, spacing_km = 30, depth_step_m = 10,
    time_step = 1, time_unit = "hours", interpolate_depth = TRUE,
    interpolate_time = TRUE, quiet = TRUE
  )
  balanced <- do.call(extract_to_track, arguments)
  expect_identical(balanced$settings$sampling, "balanced")
  expect_equal(nrow(balanced$samples), 25)
  expect_equal(balanced$samples$depth, seq(30, 180, length.out = 25))
  expect_equal(balanced$samples$date, .track_clock(0:24))
  expect_equal(balanced$samples$lon, seq(0, 2, length.out = 25))
  expect_lte(max(diff(balanced$samples$distance_km)), 30)
  expect_equal(balanced$samples$source_row[c(1, 25)], c(1, 2))
  expect_true(all(is.na(balanced$samples$source_row[2:24])))
  expect_equal(balanced$summary$n_observed, 2)
  expect_equal(balanced$summary$n_assumed, 23)
  expect_equal(balanced$summary$n_samples, 25)
  expect_equal(balanced$summary$n_valid, 25)
  expect_equal(balanced$settings$sampling_counts, data.frame(
    track_id = "1", n_observed = 2L, n_assumed = 23L, n_samples = 25L
  ))
  merged <- do.call(extract_to_track, c(arguments, list(
    sampling = "merged_steps"
  )))
  expect_equal(nrow(merged$samples), 37)
  expect_equal(merged$summary$n_assumed, 35)
  expect_true(any(abs(merged$samples$depth - 100) < 1e-10))
  expect_false(any(abs(balanced$samples$depth - 100) < 1e-10))
  expect_lte(max(diff(merged$samples$depth)), 10)
  expect_lte(max(as.numeric(diff(merged$samples$date), units = "hours")), 1)
  expect_error(do.call(extract_to_track, c(arguments, list(
    sampling = "merged_steps", max_samples = 30
  ))), "max_samples")
  expect_equal(nrow(do.call(extract_to_track, c(arguments, list(
    max_samples = 30
  )))$samples), 25)
  expect_error(extract_to_track(route, path, sampling = "regular"), "arg")
})


test_that("balanced distance or depth can set the common subdivision", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path)
  route <- data.frame(
    lon = c(0, 2), lat = 0, depth = c(0, 100),
    date = .track_clock(c(0, 2))
  )
  horizontal <- extract_to_track(
    route, path, sampling = "balanced", spacing_km = 30, depth_step_m = 50,
    time_step = 1, interpolate_depth = TRUE, interpolate_time = TRUE,
    quiet = TRUE
  )
  expect_equal(nrow(horizontal$samples), 9)
  expect_equal(horizontal$samples$depth, seq(0, 100, length.out = 9))
  expect_equal(horizontal$samples$date, .track_clock(seq(0, 2, by = 0.25)))
  vertical <- extract_to_track(
    route, path, spacing_km = 200, depth_step_m = 10, time_step = 1,
    interpolate_depth = TRUE, interpolate_time = TRUE, quiet = TRUE
  )
  expect_equal(nrow(vertical$samples), 11)
  expect_equal(vertical$samples$depth, seq(0, 100, by = 10))
  expect_equal(vertical$samples$date, .track_clock(seq(0, 2, by = 0.2)))
})


test_that("balanced calendar spans retain anchors and share inner cuts", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path)
  route <- data.frame(
    lon = c(0, 2), lat = 0,
    date = as.Date(c("2020-01-31", "2020-03-29"))
  )
  out <- extract_to_track(
    route, path, spacing_km = 60, time_step = 1, time_unit = "months",
    interpolate_time = TRUE, quiet = TRUE
  )
  expect_equal(nrow(out$samples), 5)
  expect_equal(out$samples$date, as.POSIXct(c(
    "2020-01-31 00:00:00", "2020-02-14 12:00:00", "2020-02-29 00:00:00",
    "2020-03-14 12:00:00", "2020-03-29 00:00:00"
  ), tz = "UTC"))
  expect_lte(max(diff(out$samples$distance_km)), 60)
  expect_equal(out$settings$sampling_counts$n_assumed, 3)
})


test_that("sampling counts retain tracks and missing results", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path, missing = TRUE)
  route <- data.frame(
    trip = c("b", "b", "a", "a"), lon = c(0, 2, 0, 1), lat = 0
  )
  observed <- extract_to_track(route, path, track_col = "trip", quiet = TRUE)
  expect_equal(observed$settings$sampling_counts, data.frame(
    track_id = c("b", "a"), n_observed = c(2L, 2L),
    n_assumed = c(0L, 0L), n_samples = c(2L, 2L)
  ))
  out <- extract_to_track(
    route, path, track_col = "trip", spacing_km = 112,
    statistics = c("mean", "max"), quiet = TRUE
  )
  expect_equal(out$settings$sampling_counts, data.frame(
    track_id = c("b", "a"), n_observed = c(2L, 2L),
    n_assumed = c(1L, 0L), n_samples = c(3L, 2L)
  ))
  expect_equal(out$summary$n_samples, c(3, 2))
  expect_equal(out$summary$n_observed, c(2, 2))
  expect_equal(out$summary$n_assumed, c(1, 0))
  expect_equal(out$summary$n_valid, c(2, 1))
  only <- extract_to_track(
    route, path, track_col = "trip", spacing_km = 112,
    statistics = NULL, quiet = TRUE
  )
  expect_null(only$summary)
  expect_equal(only$settings$sampling_counts, out$settings$sampling_counts)
  cells <- extract_to_track(
    route[1:2, ], path, sampling = "cells", quiet = TRUE
  )
  expect_equal(cells$settings$sampling_counts, data.frame(
    track_id = "1", n_observed = 2L, n_assumed = 3L, n_samples = 5L
  ))
})


test_that("buffered tracks share sampling choices and counts", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path)
  route <- data.frame(lon = c(0, 2), lat = 0, date = .track_clock(c(0, 3)))
  out <- extract_to_buffered_track(
    route, path, spacing_km = 112, time_step = 1, interpolate_time = TRUE,
    horizontal_radius_km = 1, infer_bounds = TRUE, quiet = TRUE
  )
  expect_identical(out$settings$sampling, "balanced")
  expect_equal(out$samples$date, .track_clock(0:3))
  expect_equal(out$settings$sampling_counts, data.frame(
    track_id = "1", n_observed = 2L, n_assumed = 2L, n_samples = 4L
  ))
  expect_equal(out$summary$n_assumed, 2)
  merged <- extract_to_buffered_track(
    route, path, sampling = "merged_steps", spacing_km = 112,
    time_step = 1, interpolate_time = TRUE, horizontal_radius_km = 1,
    infer_bounds = TRUE, statistics = NULL, quiet = TRUE
  )
  expect_null(merged$summary)
  expect_equal(merged$settings$sampling_counts$n_samples, 5)
  expect_equal(merged$settings$sampling_counts$n_assumed, 3)
  expect_error(extract_to_buffered_track(
    route, path, sampling = "regular", horizontal_radius_km = 1
  ), "arg")
})


test_that(
  "vertices reuse known point values a ... original row order",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    vertices <- data.frame(
      id = c("b", "a", "b", "a"),
      lon = c(
        0,
        2,
        1,
        1
      ),
      lat = 0,
      note = letters[1:4]
    )
    out <- extract_to_track(
      vertices,
      path,
      track_col = "id",
      sampling = "vertices",
      quiet = TRUE
    )
    expect_named(out, c("samples", "summary", "settings"))
    expect_equal(out$samples$temp, c(0, 20, 10, 10))
    expect_identical(out$samples$note, letters[1:4])
    expect_identical(out$samples$source_row, 1:4)
    expect_identical(out$summary$track_id, c("b", "a"))
    expect_equal(out$samples$sample_id, c(1, 1, 2, 2))
    expect_false(any(out$samples$is_interpolated))
  }
)


test_that(
  "balanced sampling inserts the route, ... and bounds spacing",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    out <- extract_to_track(
      data.frame(
        lon = c(0, 2),
        lat = 0,
        note = c(
          "start",
          "end"
        )
      ),
      path,
      spacing_km = 112,
      quiet = TRUE
    )
    expect_equal(out$samples$lon, c(0, 1, 2), tolerance = 1e-12)
    expect_equal(out$samples$temp, c(0, 10, 20))
    expect_identical(out$samples$is_interpolated, c(FALSE, TRUE, FALSE))
    expect_equal(out$samples$source_row, c(1, NA, 2))
    expect_identical(out$samples$note, c("start", NA, "end"))
    expect_lte(max(diff(out$samples$distance_km)), 112)
    expect_equal(
      sum(out$samples$support_km),
      222.638981586547,
      tolerance = 1e-08
    )
    expect_equal(out$summary$mean, 10)
    expect_identical(out$summary$status, "complete")
  }
)


test_that(
  "known statistics use endpoint suppo ... qual sample counts",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    vertices <- data.frame(lon = 0:2, lat = 0)
    out <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      statistics = c(
        "mean",
        "sd",
        "se",
        "min",
        "max"
      ),
      quiet = TRUE
    )
    expect_equal(out$summary$mean, 10)
    expect_equal(out$summary$sd, sqrt(80))
    expect_equal(out$summary$se, sqrt(30))
    expect_equal(out$summary$min, 0)
    expect_equal(out$summary$max, 20)
    expect_equal(out$summary$effective_n, 8 / 3)
    expect_identical(out$summary$variation_status, "available")
    equal <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      weighting = "equal",
      statistics = c("sd", "se"),
      quiet = TRUE
    )
    expect_equal(equal$summary$sd, 10)
    expect_equal(equal$summary$se, 10 / sqrt(3))
  }
)


test_that(
  "time support and distance support a ... distinct questions",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    vertices <- data.frame(
      lon = 0:2,
      lat = 0,
      date = .track_clock(
        c(
          0,
          1,
          4
        )
      )
    )
    out <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      weighting = "duration",
      quiet = TRUE
    )
    expect_equal(out$samples$support_seconds, c(1800, 7200, 5400))
    expect_equal(out$summary$mean, 12.5)
    expect_equal(out$summary$total_duration_seconds, 14400)
    expect_equal(out$summary$valid_fraction, 1)
  }
)


test_that(
  "partial coverage retains the full requested denominator",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    vertices <- data.frame(lon = c(0, 1, 3), lat = 0)
    rejected <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      quiet = TRUE
    )
    expect_equal(rejected$samples$temp, c(0, 10, NA))
    expect_equal(rejected$summary$valid_fraction, 2 / 3)
    expect_true(is.na(rejected$summary$mean))
    expect_identical(rejected$summary$status, "insufficient_coverage")
    expect_identical(rejected$summary$variation_status, "not_requested")
    accepted <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      min_valid_fraction = 2 / 3,
      statistics = c("mean", "min", "max"),
      quiet = TRUE
    )
    expect_equal(accepted$summary$mean, 7.5)
    expect_equal(accepted$summary$min, 0)
    expect_equal(accepted$summary$max, 10)
    expect_identical(accepted$summary$status, "partial")
    expect_message(
      extract_to_track(
        vertices,
        path,
        sampling = "vertices",
        time_report = FALSE,
        spatial_report = FALSE
      ),
      "Track report"
    )
  }
)


test_that(
  "all-missing and stationary tracks h ... cit summary states",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, missing = TRUE)
    vertices <- data.frame(lon = c(0, 0.1), lat = 0)
    out <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      min_valid_fraction = 0,
      quiet = TRUE
    )
    expect_identical(out$summary$status, "all_missing")
    expect_true(is.na(out$summary$mean))
    stationary <- data.frame(
      lon = c(1, 1),
      lat = 0,
      date = .track_clock(
        c(
          0,
          1
        )
      )
    )
    out <- extract_to_track(
      stationary,
      path,
      sampling = "vertices",
      quiet = TRUE
    )
    expect_true(is.na(out$summary$mean))
    expect_identical(out$summary$status, "empty")
    timed <- extract_to_track(
      stationary,
      path,
      sampling = "vertices",
      weighting = "duration",
      quiet = TRUE
    )
    expect_equal(timed$summary$mean, 10)
    expect_equal(timed$summary$total_distance_km, 0)
  }
)


test_that(
  "depth interpolation requires opt-in ... mples known layers",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, depth = c(0, 10, 20))
    vertices <- data.frame(lon = c(0, 2), lat = 0, depth = c(0, 20))
    expect_error(
      extract_to_track(vertices, path, spacing_km = 112),
      "interpolate_depth"
    )
    out <- extract_to_track(
      vertices,
      path,
      spacing_km = 112,
      interpolate_depth = TRUE,
      quiet = TRUE
    )
    expect_equal(out$samples$depth, c(0, 10, 20))
    expect_equal(out$samples$temp, c(0, 20, 40))
    expect_equal(out$summary$mean, 20)
    vertical <- data.frame(lon = 0, lat = 0, depth = c(0, 20))
    sampled <- extract_to_track(
      vertical,
      path,
      spacing_km = 1,
      depth_step_m = 10,
      interpolate_depth = TRUE,
      weighting = "equal",
      quiet = TRUE
    )
    expect_equal(sampled$samples$temp, c(0, 10, 20))
    expect_equal(sampled$samples$distance_km, c(0, 0, 0))
    expect_equal(sampled$summary$mean, 10)
  }
)


test_that(
  "time interpolation preserves UTC in ... mporal step limits",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, time = 0:2)
    vertices <- data.frame(
      lon = c(0, 2),
      lat = 0,
      date = .track_clock(
        c(
          0,
          2
        )
      )
    )
    expect_error(
      extract_to_track(vertices, path, spacing_km = 112),
      "interpolate_time"
    )
    out <- extract_to_track(
      vertices,
      path,
      spacing_km = 112,
      interpolate_time = TRUE,
      time_match = "exact",
      quiet = TRUE
    )
    expect_equal(out$samples$date, .track_clock())
    expect_equal(out$samples$temp, c(0, 1010, 2020))
    expect_equal(out$summary$mean, 1010)
    stationary <- data.frame(
      lon = 0,
      lat = 0,
      date = .track_clock(
        c(
          0,
          2
        )
      )
    )
    out <- extract_to_track(
      stationary,
      path,
      spacing_km = 2,
      interpolate_time = TRUE,
      time_step_seconds = 3600,
      weighting = "duration",
      quiet = TRUE
    )
    expect_equal(out$samples$temp, c(0, 1000, 2000))
    expect_equal(out$summary$mean, 1000)
  }
)


test_that(
  "POSIX local times keep instants and ... ns even when quiet",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, time = 0:2)
    local <- as.POSIXct(
      c("2020-01-01 01:00:00", "2020-01-01 03:00:00"),
      tz = "Europe/Paris"
    )
    vertices <- data.frame(lon = c(0, 2), lat = 0, date = local)
    out <- extract_to_track(
      vertices,
      path,
      spacing_km = 112,
      interpolate_time = TRUE,
      quiet = TRUE
    )
    expect_equal(as.numeric(out$samples$date), as.numeric(.track_clock()))
    expect_equal(out$samples$temp, c(0, 1010, 2020))
    vertices$date <- c("2020-01-01 00:00", "2020-01-01 02:00")
    expect_warning(
      out <- extract_to_track(
        vertices,
        path,
        spacing_km = 112,
        interpolate_time = TRUE,
        quiet = TRUE
      ),
      "UTC"
    )
    expect_equal(out$samples$temp, c(0, 1010, 2020))
  }
)


test_that(
  "missing dates and depths are retain ...  missing endpoints",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, depth = c(0, 10, 20), time = 0:2)
    vertices <- data.frame(
      lon = c(0, 2),
      lat = 0,
      depth = c(0, NA),
      date = .track_clock(c(0, NA))
    )
    out <- extract_to_track(
      vertices,
      path,
      spacing_km = 112,
      interpolate_depth = TRUE,
      interpolate_time = TRUE,
      quiet = TRUE
    )
    expect_equal(out$samples$depth, c(0, NA, NA))
    expect_equal(
      as.numeric(out$samples$date),
      c(
        as.numeric(.track_clock(0)),
        NA,
        NA
      )
    )
    expect_equal(out$samples$temp, c(0, NA, NA))
    expect_equal(out$summary$valid_fraction, 0.25)
    timed <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      weighting = "duration",
      quiet = TRUE
    )
    expect_true(is.na(timed$summary$valid_fraction))
    expect_identical(timed$summary$status, "empty")
  }
)


test_that(
  "surface bottom all and 2d follow th ...  selection methods",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, depth = c(0, 10, 20))
    vertices <- data.frame(
      lon = c(0, 2),
      lat = 0,
      depth = c(
        10,
        20
      )
    )
    surface <- extract_to_track(
      vertices,
      path,
      method = "surface",
      spacing_km = 112,
      quiet = TRUE
    )
    bottom <- extract_to_track(
      vertices,
      path,
      method = "bottom",
      spacing_km = 112,
      quiet = TRUE
    )
    expect_equal(surface$samples$temp, c(0, 10, 20))
    expect_equal(bottom$samples$temp, c(20, 30, 40))
    all <- extract_to_track(
      vertices,
      path,
      method = "all",
      spacing_km = 112,
      interpolate_depth = TRUE,
      quiet = TRUE
    )
    expect_equal(all$samples$surface_temp, c(0, 10, 20))
    expect_equal(all$samples$seabottom_temp, c(20, 30, 40))
    expect_equal(nrow(all$summary), 3)
    expect_error(
      extract_to_track(
        vertices,
        path,
        method = "2d",
        sampling = "vertices"
      ),
      "depth levels"
    )
  }
)


test_that(
  "mixed static temporal variables and ... es are independent",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, time = 0:2, extra_static = TRUE)
    vertices <- data.frame(lon = c(0, 2), lat = 0)
    static <- extract_to_track(
      vertices,
      path,
      var = "static",
      spacing_km = 112,
      quiet = TRUE
    )
    expect_equal(static$samples$static, rep(7, 3))
    expect_error(
      extract_to_track(vertices, path, var = "temp", spacing_km = 112),
      "date"
    )
    vertices$date <- .track_clock(c(0, 2))
    out <- extract_to_track(
      vertices,
      path,
      var = c("temp", "static"),
      spacing_km = 112,
      interpolate_time = TRUE,
      output_col = c(
        temp = "temperature",
        static = "index"
      ),
      quiet = TRUE
    )
    expect_equal(out$samples$temperature, c(0, 1010, 2020))
    expect_equal(out$samples$index, rep(7, 3))
    expect_identical(out$summary$output_column, c("temperature", "index"))
    expect_error(
      extract_to_track(vertices, path, spacing_km = 112),
      class = "ocean3d_error_var_detection"
    )
  }
)


test_that(
  "complementary singleton files are m ... ly in either order",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    .make_track_nc(paths[1L], time = 0)
    .make_track_nc(paths[2L], time = 2)
    vertices <- data.frame(
      lon = c(0, 2),
      lat = 0,
      date = .track_clock(
        c(
          0,
          2
        )
      )
    )
    for (files in list(paths, rev(paths))) {
      out <- extract_to_track(
        vertices,
        files,
        sampling = "vertices",
        quiet = TRUE
      )
      expect_equal(out$samples$temp, c(0, 2020))
      expect_equal(out$summary$mean, 1010)
    }
    table <- data.frame(path = paths)
    out <- extract_to_track(
      vertices,
      table,
      sampling = "vertices",
      file_col = "path",
      quiet = TRUE
    )
    expect_equal(out$samples$temp, c(0, 2020))
    handle <- ncdf4::nc_open(paths[1L])
    on.exit(ncdf4::nc_close(handle), add = TRUE)
    out <- extract_to_track(
      transform(
        vertices,
        date = .track_clock(
          c(
            0,
            0
          )
        )
      ),
      handle,
      sampling = "vertices",
      quiet = TRUE
    )
    expect_equal(out$samples$temp, c(0, 20))
    expect_equal(as.numeric(ncdf4::ncvar_get(handle, "temp"))[1L], 0)
  }
)


test_that(
  "custom columns matrices and named l ... xtraction contract",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    vertices <- data.frame(x = c(0, 2), y = 0, trip = c("haul", "haul"))
    out <- extract_to_track(
      vertices,
      path,
      lon_col = "x",
      lat_col = "y",
      track_col = "trip",
      spacing_km = 112,
      output_col = "temperature",
      quiet = TRUE
    )
    expect_equal(out$samples$temperature, c(0, 10, 20))
    expect_equal(out$samples$x, c(0, 1, 2))
    for (
      input in list(
        as.matrix(vertices[c("x", "y")]),
        as.list(
          vertices[c(
            "x",
            "y"
          )]
        )
      )
    ) {
      out <- extract_to_track(
        input,
        path,
        lon_col = "x",
        lat_col = "y",
        sampling = "vertices",
        quiet = TRUE
      )
      expect_equal(out$samples$temp, c(0, 20))
    }
  }
)


test_that(
  "sf lines transform CRS and accept e ... -vertex attributes",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, depth = c(0, 10, 20), time = 0:2)
    line <- sf::st_sf(
      id = "haul",
      depth = I(list(c(0, 20))),
      date = I(
        list(
          .track_clock(
            c(
              0,
              2
            )
          )
        )
      ),
      geometry = sf::st_sfc(
        sf::st_linestring(
          rbind(
            c(0, 0),
            c(2, 0)
          )
        ),
        crs = 4326
      )
    )
    for (input in list(line, sf::st_transform(line, 3857))) {
      out <- extract_to_track(
        input,
        path,
        track_col = "id",
        spacing_km = 112,
        interpolate_depth = TRUE,
        interpolate_time = TRUE,
        quiet = TRUE
      )
      expect_equal(out$samples$temp, c(0, 1020, 2040))
      expect_identical(out$samples$track_id, rep("haul", 3))
    }
    bad_list <- line
    bad_list$depth <- I(list(c(0, 10, 20)))
    expect_error(
      extract_to_track(bad_list, path, interpolate_depth = TRUE, quiet = TRUE),
      "one value or one per vertex in feature 1"
    )
    measured <- sf::st_sf(
      geometry = sf::st_sfc(
        sf::st_linestring(rbind(c(0, 0, 1), c(2, 0, 2)), dim = "XYM"),
        crs = 4326
      )
    )
    expect_error(extract_to_track(measured, path), "Measured")
    line3d <- sf::st_sf(
      date = .track_clock(0),
      geometry = sf::st_sfc(
        sf::st_linestring(
          rbind(
            c(
              0,
              0,
              0
            ),
            c(2, 0, 20)
          )
        ),
        crs = 4326
      )
    )
    expect_error(extract_to_track(line3d, path, spacing_km = 112), "z_as_depth")
    out <- extract_to_track(
      line3d,
      path,
      spacing_km = 112,
      z_as_depth = TRUE,
      interpolate_depth = TRUE,
      quiet = TRUE
    )
    expect_equal(out$samples$temp, c(0, 20, 40))
  }
)


test_that(
  "antimeridian routes take the short  ... r arcs stay finite",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, lon = -180:179, lat = c(0, 80, 85))
    out <- extract_to_track(
      data.frame(lon = c(179, -179), lat = 0), path,
      spacing_km = 112, distance_model = "sphere", quiet = TRUE
    )
    expect_equal(abs(out$samples$lon[2L]), 180)
    expect_equal(out$samples$distance_km[3L], 222.390202355, tolerance = 1e-08)
    polar <- extract_to_track(
      data.frame(lon = c(-45, 45), lat = 80), path,
      spacing_km = 1000, distance_model = "sphere", quiet = TRUE
    )
    expect_true(all(is.finite(c(polar$samples$lon, polar$samples$lat))))
    expect_gt(polar$samples$lat[2L], 80)
    expect_error(
      extract_to_track(
        data.frame(lon = c(0, 180), lat = 0), path,
        distance_model = "sphere", quiet = TRUE
      ),
      "Antipodal"
    )
  }
)


test_that(
  "default sampling retains vertices a ... ces memory budgets",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    vertices <- data.frame(lon = c(0, 2), lat = 0)
    out <- extract_to_track(vertices, path, quiet = TRUE)
    expect_null(out$settings$spacing_km)
    expect_equal(out$samples$source_row, c(1L, 2L))
    expect_equal(out$samples$temp, c(0, 20))
    expect_error(
      extract_to_track(vertices, path, sampling = "balanced"),
      "explicit"
    )
    expect_error(
      extract_to_track(
        vertices,
        path,
        spacing_km = 1e-09,
        max_samples = 100
      ),
      "max_samples"
    )
    expect_error(
      extract_to_track(vertices, path, max_samples = 1),
      "max_samples"
    )
    vertices$id <- c("a", "a")
    two <- rbind(vertices, transform(vertices, id = "b"))
    expect_error(
      extract_to_track(
        two,
        path,
        track_col = "id",
        spacing_km = 112,
        max_samples = 5
      ),
      "max_samples"
    )
  }
)


test_that(
  "time ordering and gap limits preven ... tended connections",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, time = 0:2)
    vertices <- data.frame(
      lon = c(0, 2),
      lat = 0,
      date = .track_clock(
        c(
          2,
          0
        )
      )
    )
    expect_error(
      extract_to_track(vertices, path, sampling = "vertices"),
      "traversal order"
    )
    vertices$date <- .track_clock(c(0, 2))
    expect_error(
      extract_to_track(
        vertices,
        path,
        sampling = "vertices",
        max_gap_hours = 1
      ),
      "time gap"
    )
    expect_error(
      extract_to_track(
        vertices,
        path,
        sampling = "vertices",
        max_segment_km = 100
      ),
      "segment"
    )
    vertices$date <- c("2020-01-01T00:00Z", "not a date")
    expect_error(
      extract_to_track(vertices, path, sampling = "vertices"),
      "not a date"
    )
  }
)


test_that(
  "invalid controls names and unsuppor ... tions fail clearly",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    vertices <- data.frame(lon = c(0, 2), lat = 0)
    for (bad in list(NA, 0, -1, Inf, "one", c(1, 2))) {
      expect_error(
        extract_to_track(vertices, path, spacing_km = bad),
        "spacing_km"
      )
    }
    for (bad in list(NA, "yes", c(TRUE, FALSE))) {
      expect_error(
        extract_to_track(vertices, path, interpolate_depth = bad),
        "interpolate_depth"
      )
    }
    for (bad in list(-1, 2, NA, "half")) {
      expect_error(
        extract_to_track(vertices, path, min_valid_fraction = bad),
        "min_valid_fraction"
      )
    }
    for (
      bad in list(
        "median",
        c("mean", "mean"),
        NA_character_,
        character()
      )
    ) {
      expect_error(
        extract_to_track(vertices, path, statistics = bad),
        "statistics"
      )
    }
    expect_error(
      extract_to_track(vertices, path, max_samples = 2.5),
      "integer"
    )
    expect_error(
      extract_to_track(vertices, path, depth_step_m = 10),
      "interpolation switch"
    )
    expect_error(
      extract_to_track(
        vertices,
        path,
        sampling = "vertices",
        spacing_km = 1
      ),
      "vertices sampling"
    )
    expect_error(
      extract_to_track(vertices, path, bad_option = TRUE),
      "point-selection"
    )
    expect_error(
      extract_to_track(vertices, path, output_col = "lon"),
      "overwrite"
    )
    expect_error(
      extract_to_track(
        transform(vertices, source_row = 1),
        path
      ),
      "reserved"
    )
    expect_error(
      extract_to_track(vertices, path, lon_col = "lat"),
      "distinct"
    )
    expect_error(
      extract_to_track(
        transform(vertices, track_id = "a"),
        path
      ),
      "reserved"
    )
    expect_error(extract_to_track(vertices[1, ], path), "two vertices")
    expect_error(extract_to_track(vertices[FALSE, ], path), "at least one")
    expect_error(
      extract_to_track(
        transform(vertices, lon = NA_real_),
        path
      ),
      "Split tracks"
    )
    expect_error(
      extract_to_track(
        transform(vertices, lat = 95),
        path
      ),
      "Latitudes"
    )
    expect_error(
      extract_to_track(vertices, path, var = "absent"),
      "No requested"
    )
  }
)


test_that(
  "samples-only mode and source settin ... elected provenance",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    out <- extract_to_track(
      data.frame(lon = c(0, 2), lat = 0),
      path,
      statistics = NULL,
      sampling = "vertices",
      quiet = TRUE
    )
    expect_null(out$summary)
    expect_equal(out$samples$temp, c(0, 20))
    expect_equal(out$settings$candidate_sources, path)
    expect_match(out$settings$integration, "approximate")
  }
)

test_that(
  "ambiguous geometries and malformed  ... butes are rejected",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    geometry <- sf::st_sfc(
      sf::st_linestring(rbind(c(0, 0), c(2, 0))),
      crs = 4326
    )
    line <- sf::st_sf(id = "one", geometry = geometry)
    expect_error(extract_to_track(sf::st_set_crs(line, NA), path), "CRS")
    points <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(0, 0))))
    sf::st_crs(points) <- 4326
    expect_error(extract_to_track(points, path), "LINESTRING")
    empty <- sf::st_sf(
      geometry = sf::st_sfc(
        sf::st_linestring(
          matrix(
            numeric(),
            ncol = 2
          )
        ),
        crs = 4326
      )
    )
    expect_error(extract_to_track(empty, path), "non-empty")
    expect_error(extract_to_track(line[FALSE, ], path), "at least one")
    repeated <- rbind(line, line)
    expect_error(
      extract_to_track(repeated, path, track_col = "id"),
      "distinct track ID"
    )
    line$depth <- I(list(c(1, 2, 3)))
    expect_error(extract_to_track(line, path), "one per vertex")
    measured <- sf::st_sf(
      geometry = sf::st_sfc(
        sf::st_linestring(
          rbind(
            c(
              0,
              0,
              1
            ),
            c(2, 0, 2)
          ),
          dim = "XYM"
        ),
        crs = 4326
      )
    )
    expect_error(extract_to_track(measured, path), "Measured")
    three <- sf::st_sf(
      depth = 2,
      geometry = sf::st_sfc(
        sf::st_linestring(
          rbind(
            c(
              0,
              0,
              1
            ),
            c(2, 0, 2)
          )
        ),
        crs = 4326
      )
    )
    expect_error(extract_to_track(three, path, z_as_depth = TRUE), "not both")
  }
)

test_that(
  "invalid schemas and variables cannot silently overwrite",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    vertices <- data.frame(lon = c(0, 2), lat = 0)
    expect_error(
      extract_to_track(
        transform(vertices, id = NA_character_),
        path,
        track_col = "id"
      ),
      "Track IDs"
    )
    expect_error(
      extract_to_track(
        transform(vertices, id = ""),
        path,
        track_col = "id"
      ),
      "Track IDs"
    )
    duplicate <- vertices
    names(duplicate) <- c("lon", "lon")
    expect_error(extract_to_track(duplicate, path), "unique")
    reserved <- cbind(as.matrix(vertices), source_row = 1)
    expect_error(extract_to_track(reserved, path), "reserved")
    for (bad in list("", NA_character_, 1, c("lon", "lat"))) {
      expect_error(
        extract_to_track(vertices, path, lon_col = bad),
        "lon_col"
      )
    }
    for (bad in list("", c("temp", "temp"), NA_character_, 1)) {
      expect_error(
        extract_to_track(vertices, path, var = bad),
        "`var`"
      )
    }
    expect_error(extract_to_track(vertices, list()), "source")
    expect_error(
      extract_to_track(
        vertices,
        path,
        var = c(
          "temp",
          "missing"
        )
      ),
      "not found: missing"
    )
    no_depth <- tempfile(fileext = ".nc")
    on.exit(unlink(no_depth), add = TRUE)
    .make_track_nc(no_depth, depth = c(0, 20))
    expect_error(extract_to_track(vertices, c(path, no_depth)), "incompatible")
    expect_error(
      extract_to_track(
        transform(vertices, depth = -1),
        no_depth
      ),
      "non-negative"
    )
  }
)

test_that(
  "singleton grids permit observations without densification",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, lon = 0, lat = 0)
    vertices <- data.frame(lon = c(0, 0), lat = 0)
    expect_equal(
      extract_to_track(vertices, path, quiet = TRUE)$samples$temp,
      c(0, 0)
    )
    expect_no_warning(
      out <- extract_to_track(
        vertices,
        path,
        spacing_km = 1,
        weighting = "equal",
        quiet = TRUE
      )
    )
    expect_equal(out$samples$temp, c(0, 0))
    expect_equal(out$summary$mean, 0)
  }
)

test_that(
  "time matching follows the original observation precision",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, time = c(0, 24, 48))
    vertices <- data.frame(
      lon = c(0, 2),
      lat = 0,
      date = as.Date(
        c(
          "2020-01-01",
          "2020-01-03"
        )
      )
    )
    out <- extract_to_track(
      vertices,
      path,
      spacing_km = 112,
      interpolate_time = TRUE,
      time_match = "exact",
      quiet = TRUE
    )
    expect_identical(out$settings$time_resolution, "daily")
    expect_equal(out$samples$temp, c(0, 24010, 48020))
    expect_error(
      extract_to_track(
        vertices,
        path,
        sampling = "vertices",
        time_resolution = "custom",
        custom_time_interval = as.difftime(
          24,
          units = "hours"
        ),
        quiet = TRUE
      ),
      "requires observation hours"
    )
    vertices$date <- as.POSIXct(vertices$date, tz = "UTC")
    custom <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      time_resolution = "custom",
      custom_time_interval = as.difftime(
        24,
        units = "hours"
      ),
      quiet = TRUE
    )
    expect_identical(custom$settings$time_resolution, "custom")
    expect_equal(custom$samples$temp, c(0, 48020))
  }
)

test_that(
  "SD can be unavailable while a valid mean is still reported",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, missing = TRUE)
    vertices <- data.frame(lon = c(0, 1), lat = 0)
    out <- extract_to_track(
      vertices,
      path,
      sampling = "vertices",
      min_valid_fraction = 0.5,
      statistics = c("mean", "sd", "se"),
      quiet = TRUE
    )
    expect_equal(out$summary$mean, 10)
    expect_true(is.na(out$summary$sd))
    expect_true(is.na(out$summary$se))
    expect_identical(
      out$summary$variation_status, "insufficient_variation_samples"
    )
  }
)


test_that("variation states distinguish requests, coverage and sample counts", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path)
  for (buffered in c(FALSE, TRUE)) {
    extract <- if (buffered) extract_to_buffered_track else extract_to_track
    settings <- list(nc = path, var = "temp", weighting = "equal", quiet = TRUE)
    if (buffered) {
      settings$horizontal_radius_km <- 1
      settings$infer_bounds <- TRUE
    }
    vertices <- data.frame(lon = 0:2, lat = 0)
    mean_only <- do.call(extract, c(list(data = vertices), settings))
    expect_identical(mean_only$summary$variation_status, "not_requested")
    expect_equal(mean_only$summary$mean, 10)
    for (requested in list("sd", "se", c("mean", "sd", "se"))) {
      out <- do.call(extract, c(
        list(data = vertices, statistics = requested), settings
      ))
      expect_identical(out$summary$variation_status, "available")
      if ("sd" %in% requested) expect_equal(out$summary$sd, 10)
      if ("se" %in% requested) expect_equal(out$summary$se, 10 / sqrt(3))
    }
    requested <- list(statistics = c("mean", "sd", "se"))
    rejected <- do.call(extract, c(
      list(data = data.frame(lon = c(0, 1, 3), lat = 0)), requested, settings
    ))
    expect_equal(rejected$summary$n_valid, 2)
    expect_equal(rejected$summary$valid_fraction, 2 / 3)
    expect_identical(rejected$summary$status, "insufficient_coverage")
    expect_identical(rejected$summary$variation_status, "not_available")
    expect_true(is.na(rejected$summary$sd))
    expect_true(is.na(rejected$summary$se))
    threshold <- if (buffered) {
      "track_min_valid_fraction"
    } else {
      "min_valid_fraction"
    }
    settings[[threshold]] <- 0.5
    singleton <- do.call(extract, c(
      list(data = data.frame(lon = c(1, 3), lat = 0)), requested, settings
    ))
    expect_identical(singleton$summary$status, "partial")
    expect_identical(
      singleton$summary$variation_status, "insufficient_variation_samples"
    )
    expect_equal(singleton$summary$mean, 10)
    expect_true(is.na(singleton$summary$sd))
    expect_true(is.na(singleton$summary$se))
    missing <- do.call(extract, c(
      list(data = data.frame(lon = 3:4, lat = 0)), requested, settings
    ))
    expect_identical(missing$summary$status, "all_missing")
    expect_identical(missing$summary$variation_status, "not_available")
    settings$weighting <- "duration"
    vertices$date <- .track_clock(c(0, 0, 0))
    stationary_time <- do.call(extract, c(
      list(data = vertices), requested, settings
    ))
    expect_identical(stationary_time$summary$status, "empty")
    expect_identical(stationary_time$summary$variation_status, "not_available")
  }
})


test_that("non-finite variation results are not marked available", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path)
  nc <- ncdf4::nc_open(path, write = TRUE)
  tryCatch(
    ncdf4::ncvar_put(nc, "temp", matrix(rep(c(-1e200, 1e200, 0), 2), 3)),
    finally = ncdf4::nc_close(nc)
  )
  out <- extract_to_track(
    data.frame(lon = 0:1, lat = 0), path, var = "temp", weighting = "equal",
    statistics = c("mean", "sd", "se"), quiet = TRUE
  )
  expect_identical(out$summary$status, "complete")
  expect_identical(out$summary$variation_status, "not_available")
  expect_equal(out$summary$mean, 0)
  expect_true(is.infinite(out$summary$sd))
  expect_true(is.infinite(out$summary$se))
})


test_that(
  "WGS84 geodesics have independent kn ... uatorial distances",
  {
    arc <- .track_geodesic(c(0, 1), c(0, 0))
    expect_equal(arc$length_km, 111.319490793274, tolerance = 1e-10)
    expect_equal(
      unname(.track_geodesic_positions(arc, 1, 0.5)),
      matrix(c(0.5, 0), nrow = 1),
      tolerance = 1e-10
    )
    expect_equal(.track_geodesic(c(1, 1), c(0, 0))$length_km, 0)
  }
)


test_that(
  "cell intervals preserve crossed cells and return visits",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    route <- data.frame(lon = c(0, 2, 0), lat = 0)
    out <- extract_to_track(route, path, sampling = "cells", quiet = TRUE)
    intervals <- subset(out$samples, sample_kind == "cell_interval")
    expect_equal(intervals$temp, c(0, 10, 20, 20, 10, 0))
    expect_equal(intervals$segment_id, c(1, 1, 1, 2, 2, 2))
    expect_equal(
      out$samples$source_row[!out$samples$is_interpolated],
      1:3
    )
    expect_equal(out$summary$mean, 10)
    expect_equal(
      out$summary$total_distance_km,
      445.277963173095,
      tolerance = 1e-09
    )
    expect_equal(sum(intervals$support_km), out$summary$total_distance_km)
    expect_equal(out$settings$sampling, "cells")
  }
)


test_that(
  "cells split changing depth and time ... ntal interpolation",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, depth = c(0, 20), time = c(0, 24))
    route <- data.frame(
      lon = c(1, 1),
      lat = 0,
      depth = c(0, 20),
      date = .track_clock(c(0, 24))
    )
    expect_error(
      extract_to_track(route, path, sampling = "cells"),
      "interpolat"
    )
    expect_warning(
      out <- extract_to_track(
        route,
        path,
        sampling = "cells",
        interpolate_depth = TRUE,
        interpolate_time = TRUE,
        time_resolution = "datetime",
        weighting = "duration",
        quiet = TRUE
      ),
      "Temporal warning"
    )
    intervals <- subset(out$samples, sample_kind == "cell_interval")
    expect_equal(intervals$temp, c(10, 24030))
    expect_equal(intervals$support_seconds, c(43200, 43200))
    expect_equal(out$summary$mean, 12020)
  }
)


test_that(
  "calendar steps clamp month ends and ... it sampling intent",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    route <- data.frame(
      lon = c(0, 2),
      lat = 0,
      date = as.Date(
        c(
          "2020-01-31",
          "2020-04-01"
        )
      )
    )
    out <- extract_to_track(
      route,
      path,
      interpolate_time = TRUE,
      time_step = 1,
      time_unit = "months",
      quiet = TRUE
    )
    expect_equal(
      format(out$samples$date, "%Y-%m-%d", tz = "UTC"),
      c("2020-01-31", "2020-02-29", "2020-03-29", "2020-04-01")
    )
    expect_error(
      extract_to_track(
        route,
        path,
        interpolate_time = TRUE,
        time_step = 0.5,
        time_unit = "years"
      ),
      "whole"
    )
  }
)




test_that(
  "marine masks examine connections, r ... ct unknown support",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    route <- data.frame(lon = c(0, 2), lat = 0)
    island <- .track_test_polygon(0.9, 1.1, -0.1, 0.1)
    domain <- .track_test_polygon(-1, 3, -1, 2)
    out <- extract_to_track(
      route,
      path,
      spacing_km = 50,
      land_mask = island,
      mask_domain = domain,
      quiet = TRUE
    )
    expect_equal(
      out$samples$temp[!out$samples$is_interpolated],
      c(0, 20)
    )
    expect_true(all(is.na(out$samples$temp[out$samples$is_interpolated])))
    expect_equal(out$summary$valid_fraction, 0)
    expect_true(is.na(out$summary$mean))
    expect_error(
      extract_to_track(
        route,
        path,
        land_mask = island,
        mask_domain = domain,
        land_policy = "error"
      ),
      "segment"
    )
    far <- .track_test_polygon(4, 5, -0.1, 0.1)
    unknown <- extract_to_track(route, path, land_mask = far, quiet = TRUE)
    expect_true(all(unknown$samples$segment_status == "unknown_mask_coverage"))
    expect_equal(unknown$summary$valid_fraction, 0)
    marine <- extract_to_track(
      route,
      path,
      land_mask = domain,
      mask_type = "marine",
      quiet = TRUE
    )
    expect_equal(marine$summary$mean, 10)
    raster <- terra::rast(
      nrows = 2,
      ncols = 3,
      xmin = -0.5,
      xmax = 2.5,
      ymin = -0.5,
      ymax = 1.5,
      crs = "EPSG:4326"
    )
    terra::values(raster) <- 1
    good <- extract_to_track(
      route,
      path,
      land_mask = raster,
      marine_values = 1,
      quiet = TRUE
    )
    expect_equal(good$summary$valid_fraction, 1)
    terra::values(raster)[5] <- NA
    missing <- extract_to_track(
      route,
      path,
      land_mask = raster,
      marine_values = 1,
      quiet = TRUE
    )
    expect_equal(missing$summary$valid_fraction, 0)
  }
)


test_that(
  "buffered tracks separate target int ... er/track summaries",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    route <- data.frame(lon = c(0, 2), lat = 0, target = c(1, 19))
    mean <- extract_to_buffered_track(
      route,
      path,
      var = "temp",
      horizontal_radius_km = 1,
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_equal(mean$samples$temp, c(0, 20))
    expect_equal(mean$summary$mean, 10)
    expect_error(
      extract_to_buffered_track(
        route,
        path,
        spacing_km = 100,
        statistic = "most_similar",
        target_col = "target",
        quiet = TRUE
      ),
      "interpolate_target"
    )
    selected <- extract_to_buffered_track(
      route,
      path,
      spacing_km = 100,
      interpolate_target = TRUE,
      statistic = "most_similar",
      target_col = "target",
      horizontal_radius_km = 60,
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_equal(selected$samples$target, c(1, 7, 13, 19))
    expect_equal(selected$samples$temp, c(0, 10, 10, 20))
    expect_equal(selected$settings$buffer_min_valid_fraction, 0)
    expect_equal(selected$settings$track_min_valid_fraction, 0.8)
    expect_true(
      all(
        vapply(
          selected$samples$temp_match_info,
          function(x) {
            nrow(x$target_match$matches) >=
              1
          },
          logical(1)
        )
      )
    )
  }
)


test_that(
  "masks treat coastal contact conserv ... nvalid connections",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    island <- .track_test_polygon(0.9, 1.1, 0, 0.1)
    domain <- .track_test_polygon(-1, 3, -1, 2)
    route <- data.frame(lon = c(0, 2, 2), lat = c(0, 0, 1))
    out <- extract_to_track(
      route,
      path,
      sampling = "cells",
      land_mask = island,
      mask_domain = domain,
      quiet = TRUE
    )
    intervals <- subset(out$samples, sample_kind == "cell_interval")
    expect_true(all(is.na(intervals$temp[intervals$segment_id == 1])))
    expect_true(
      all(
        is.finite(
          intervals$temp[intervals$segment_id ==
                           2]
        )
      )
    )
    expect_equal(
      out$samples$temp[!out$samples$is_interpolated],
      c(
        0,
        20,
        120
      )
    )
    expect_equal(out$summary$status, "insufficient_coverage")
    raster <- terra::rast(
      nrows = 1,
      ncols = 1,
      xmin = -0.5,
      xmax = 0.5,
      ymin = -0.5,
      ymax = 0.5,
      crs = "EPSG:4326"
    )
    terra::values(raster) <- 1
    expect_error(
      extract_to_track(route, path, land_mask = raster),
      "marine_values"
    )
    outside <- extract_to_track(
      route,
      path,
      land_mask = raster,
      marine_values = 1,
      quiet = TRUE
    )
    expect_equal(outside$summary$valid_fraction, 0)
  }
)


test_that(
  "buffered cells retain target NA and ... t extra permission",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    route <- data.frame(lon = c(0, 2), lat = 0, target = 10)
    out <- extract_to_buffered_track(
      route,
      path,
      sampling = "cells",
      statistic = "most_similar",
      target_col = "target",
      horizontal_radius_km = 120,
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_true(all(out$samples$target == 10))
    expect_true(all(out$samples$temp == 10))
    missing <- extract_to_buffered_track(
      transform(route, target = NA_real_),
      path,
      sampling = "cells",
      statistic = "most_similar",
      target_col = "target",
      horizontal_radius_km = 1,
      infer_bounds = TRUE,
      quiet = TRUE
    )
    expect_true(all(is.na(missing$samples$temp)))
    expect_true(is.na(missing$summary$mean))
  }
)
# Overlapping grids must not silently approximate source-switch boundaries.
test_that(
  "cells accept disjoint sources and r ... compatible overlap",
  {
    paths <- c(tempfile(fileext = ".nc"), tempfile(fileext = ".nc"))
    on.exit(unlink(paths))
    .make_track_nc(paths[1L], lon = 0:2)
    .make_track_nc(paths[2L], lon = 3:5)
    route <- data.frame(lon = c(0, 5), lat = 0)
    out <- extract_to_track(
      route,
      paths,
      var = "temp",
      sampling = "cells",
      quiet = TRUE
    )
    expect_equal(
      subset(out$samples, sample_kind == "cell_interval")$temp,
      c(0, 10, 20, 30, 40, 50)
    )
    .make_track_nc(paths[2L], lon = c(0.2, 1.2, 2.2))
    expect_error(
      extract_to_track(route, paths, var = "temp", sampling = "cells"),
      "overlapping different grids"
    )
  }
)
test_that(
  "cells detect thin cells that coarse ... r samples can miss",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, lon = c(0, 0.999, 1, 1.001, 2))
    route <- data.frame(lon = c(0, 2), lat = 0)
    cells <- extract_to_track(route, path, sampling = "cells", quiet = TRUE)
    expect_equal(
      subset(cells$samples, sample_kind == "cell_interval")$temp,
      c(0, 9.99, 10, 10.01, 20)
    )
    balanced <- extract_to_track(route, path, spacing_km = 100, quiet = TRUE)
    expect_false(any(balanced$samples$temp == 10))
  }
)

test_that(
  "cells use valid depth transitions and merge redundant cuts",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, depth = c(0, 10, 20))
    handle <- ncdf4::nc_open(path, write = TRUE)
    values <- ncdf4::ncvar_get(handle, "temp")
    values[, , 2L] <- NA_real_
    ncdf4::ncvar_put(handle, "temp", values)
    ncdf4::nc_close(handle)
    route <- data.frame(lon = c(1, 1), lat = 0, depth = c(0, 20))
    out <- extract_to_track(
      route,
      path,
      sampling = "cells",
      interpolate_depth = TRUE,
      weighting = "equal",
      quiet = TRUE
    )
    intervals <- subset(out$samples, sample_kind == "cell_interval")
    expect_equal(intervals$temp, c(10, 30))
    expect_equal(intervals$depth, c(5, 15))
    expect_equal(
      vapply(
        intervals$temp_match_info,
        function(x) x$depth_m,
        numeric(1)
      ),
      c(0, 20)
    )
    expect_equal(out$summary$mean, 20)
  }
)

test_that(
  "temporal limits leave unavailable d ... in the denominator",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, time = c(0, 24))
    route <- data.frame(
      lon = c(1, 1),
      lat = 0,
      date = .track_clock(c(0, 24))
    )
    expect_warning(
      out <- extract_to_track(
        route,
        path,
        sampling = "cells",
        interpolate_time = TRUE,
        weighting = "duration",
        max_time_diff = 0.25,
        quiet = TRUE
      ),
      "Temporal warning"
    )
    intervals <- subset(out$samples, sample_kind == "cell_interval")
    expect_equal(intervals$temp, c(10, NA, 24010))
    expect_equal(intervals$support_seconds, c(21600, 43200, 21600))
    expect_equal(out$summary$valid_fraction, 0.5)
    expect_true(is.na(out$summary$mean))
  }
)

test_that(
  "time steps preserve UTC instants across daylight saving",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    route <- data.frame(
      lon = c(0, 2),
      lat = 0,
      date = as.POSIXct(
        c("2020-03-29 01:30:00", "2020-03-29 04:30:00"),
        tz = "Europe/Paris"
      )
    )
    hours <- extract_to_track(
      route,
      path,
      interpolate_time = TRUE,
      time_step = 1,
      time_unit = "hours",
      quiet = TRUE
    )
    seconds <- extract_to_track(
      route,
      path,
      interpolate_time = TRUE,
      time_step = 3600,
      time_unit = "seconds",
      quiet = TRUE
    )
    expect_equal(
      format(hours$samples$date, "%H:%M", tz = "UTC"),
      c("00:30", "01:30", "02:30")
    )
    expect_equal(seconds$samples$date, hours$samples$date)
    route$date <- as.Date(c("2020-02-29", "2023-03-01"))
    years <- extract_to_track(
      route,
      path,
      interpolate_time = TRUE,
      time_step = 1,
      time_unit = "years",
      quiet = TRUE
    )
    expect_equal(
      format(years$samples$date, "%Y-%m-%d", tz = "UTC"),
      c("2020-02-29", "2021-02-28", "2022-02-28", "2023-02-28", "2023-03-01")
    )
  }
)

test_that(
  "invalid routes are distinct from va ... nmental candidates",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    route <- data.frame(lon = c(0, 2), lat = 0)
    island <- .track_test_polygon(0.9, 1.1, -0.1, 0.1)
    domain <- .track_test_polygon(-1, 3, -1, 2)
    out <- extract_to_track(
      route,
      path,
      spacing_km = 100,
      land_mask = island,
      mask_domain = domain,
      quiet = TRUE
    )
    i <- which(out$samples$is_interpolated)[1L]
    expect_equal(out$samples$temp_match_info[[i]]$status, "invalid_connection")
    expect_true(is.na(out$samples$temp[i]))
    expect_equal(
      out$samples$temp_match_info[[i]]$environmental_status,
      "selected"
    )
    expect_equal(
      out$samples$temp_match_info[[i]]$track_status,
      out$samples$segment_status[i]
    )
    buffered <- extract_to_buffered_track(
      route,
      path,
      spacing_km = 100,
      horizontal_radius_km = 10,
      infer_bounds = TRUE,
      land_mask = island,
      mask_domain = domain,
      quiet = TRUE
    )
    expect_equal(
      buffered$samples$temp_match_info[[i]]$status,
      "invalid_connection"
    )
    expect_equal(
      buffered$samples$temp_match_info[[i]]$environmental_status,
      "complete"
    )
    expect_equal(
      buffered$samples$temp_match_info[[i]]$track_status,
      buffered$samples$segment_status[i]
    )
    expect_equal(buffered$settings$buffer$horizontal_radius_km, 10)
    expect_identical(buffered$settings[["statistic"]], "mean")
    expect_identical(buffered$settings[["statistics"]], "mean")
    expect_identical(out$settings[["statistics"]], "mean")
    expect_true(all(nzchar(names(buffered$settings))))
    expect_true(buffered$settings$mask_checked)
  }
)
test_that(
  "cells normalize selection options before partitioning",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, time = 24)
    route <- data.frame(
      lon = c(0, 2),
      lat = 0,
      date = .track_clock(c(0, 48))
    )
    out <- extract_to_track(
      route,
      path,
      sampling = "cells",
      time_match = "ex",
      interpolate_time = TRUE,
      quiet = TRUE
    )
    expect_true(all(is.na(out$samples$temp)))
    expect_error(
      extract_to_track(
        route,
        path,
        sampling = "cells",
        spatial_match = "nea",
        max_distance_km = 10
      ),
      "Cells mode"
    )
    expect_error(
      extract_to_track(
        route,
        path,
        sampling = "cells",
        time_resolution = "cust"
      ),
      "Cells mode"
    )
    expect_error(extract_to_buffered_track(route, path, method = "all"), "all")
  }
)
test_that(
  "buffered tracks retain the original timezone warning",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, time = c(0, 24))
    route <- data.frame(
      lon = c(0, 2),
      lat = 0,
      date = "2020-01-01 00:00"
    )
    expect_warning(
      out <- extract_to_buffered_track(
        route,
        path,
        horizontal_radius_km = 10,
        infer_bounds = TRUE,
        quiet = TRUE
      ),
      "UTC"
    )
    expect_equal(as.numeric(out$samples$date), rep(1577836800, 2))
  }
)
test_that(
  "raster masks follow short antimeridian connections",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path, lon = 178:182)
    mask <- terra::rast(
      nrows = 2,
      ncols = 360,
      xmin = -180,
      xmax = 180,
      ymin = -1,
      ymax = 1,
      crs = "EPSG:4326"
    )
    terra::values(mask) <- 1
    xy <- terra::xyFromCell(mask, seq_len(terra::ncell(mask)))
    terra::values(mask)[abs(xy[, 1L]) < 1] <- 0
    route <- data.frame(lon = c(179, -179), lat = 0)
    out <- extract_to_track(
      route,
      path,
      land_mask = mask,
      marine_values = 1,
      quiet = TRUE
    )
    expect_true(all(out$samples$segment_status == "marine"))
    expect_equal(out$summary$mean, 1800)
    stationary <- extract_to_track(
      data.frame(lon = c(179, 179), lat = 0),
      path,
      land_mask = mask,
      marine_values = 1,
      weighting = "equal",
      quiet = TRUE
    )
    expect_true(all(stationary$samples$segment_status == "marine"))
    names(mask) <- "ID"
    land <- extract_to_track(
      data.frame(lon = c(-0.25, 0.25), lat = 0),
      path,
      land_mask = mask,
      marine_values = 1,
      quiet = TRUE
    )
    expect_true(all(land$samples$segment_status == "land_crossing"))
    levels(mask) <- data.frame(code = 0:1, label = c("land", "sea"))
    names(mask) <- "ID"
    categories <- extract_to_track(
      data.frame(lon = c(179, 179), lat = 0),
      path,
      land_mask = mask,
      marine_values = 1,
      weighting = "equal",
      quiet = TRUE
    )
    expect_true(all(categories$samples$segment_status == "marine"))
  }
)
test_that(
  "mask domains and longitude conventions are explicit",
  {
    path <- tempfile(fileext = ".nc")
    on.exit(unlink(path))
    .make_track_nc(path)
    route <- data.frame(lon = c(0, 2), lat = 0)
    raster <- terra::rast(
      nrows = 2,
      ncols = 3,
      xmin = -0.5,
      xmax = 2.5,
      ymin = -0.5,
      ymax = 1.5,
      crs = "EPSG:4326"
    )
    terra::values(raster) <- 1
    domain <- .track_test_polygon(-1, 3, -1, 2)
    expect_error(
      extract_to_track(
        route,
        path,
        land_mask = raster,
        marine_values = 1,
        mask_domain = domain
      ),
      "mask_domain"
    )
    invalid_domain <- sf::st_sfc(
      sf::st_linestring(
        matrix(c(0, 0, 2, 0), ncol = 2, byrow = TRUE)
      ),
      crs = 4326
    )
    expect_error(
      extract_to_track(
        route,
        path,
        land_mask = .track_test_polygon(4, 5, -1, 1),
        mask_domain = invalid_domain
      ),
      "mask_domain"
    )
    rotated <- terra::rast(
      nrows = 2,
      ncols = 3,
      xmin = 0,
      xmax = 360,
      ymin = -1,
      ymax = 1,
      crs = "EPSG:4326"
    )
    terra::values(rotated) <- 1
    expect_error(
      extract_to_track(route, path, land_mask = rotated, marine_values = 1),
      "Rotate"
    )
  }
)


test_that("tracks accept climatological time axes like point extraction", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path, time = c(0, 24))
  handle <- ncdf4::nc_open(path, write = TRUE)
  ncdf4::ncatt_put(handle, "time", "climatology", "climatology_bounds")
  ncdf4::nc_close(handle)
  route <- data.frame(
    lon = c(0, 2), lat = 0,
    date = as.POSIXct(c("2020-01-01", "2020-01-01"), tz = "UTC")
  )
  out <- extract_to_track(route, path, quiet = TRUE)
  expect_equal(out$samples$temp, c(0, 20))
})


test_that("a missing date column for a temporal variable matches points", {
  path <- tempfile(fileext = ".nc")
  on.exit(unlink(path))
  .make_track_nc(path, time = c(0, 24))
  route <- data.frame(lon = c(0, 2), lat = 0)
  err <- expect_error(extract_to_track(route, path, quiet = TRUE))
  expect_match(conditionMessage(err), "`temp` has a time dimension")
  expect_match(conditionMessage(err), "Column `date` is missing")
})
