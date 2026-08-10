# ------------------------------------------------------------------------------
# Unit tests for Copernicus helpers
#
# These tests are intentionally offline. They validate deterministic helper
# behaviour and input validation without contacting Copernicus services.
# ------------------------------------------------------------------------------

test_that("bbox validation accepts valid bounds", {
  expect_invisible(
    sharkabc3d:::.copernicus_validate_bbox(
      xmin = -5, xmax = 5,
      ymin = 35, ymax = 45
    )
  )

  expect_invisible(
    sharkabc3d:::.copernicus_validate_bbox(
      xmin = NULL, xmax = NULL,
      ymin = NULL, ymax = NULL
    )
  )
})

test_that("bbox validation rejects incomplete or invalid bounds", {
  expect_error(
    sharkabc3d:::.copernicus_validate_bbox(
      xmin = -5, xmax = 5,
      ymin = 35, ymax = NULL
    ),
    "must be supplied together"
  )

  expect_error(
    sharkabc3d:::.copernicus_validate_bbox(
      xmin = 5, xmax = -5,
      ymin = 35, ymax = 45
    ),
    "xmin < xmax"
  )

  expect_error(
    sharkabc3d:::.copernicus_validate_bbox(
      xmin = -5, xmax = 5,
      ymin = 45, ymax = 35
    ),
    "ymin < ymax"
  )

  expect_error(
    sharkabc3d:::.copernicus_validate_bbox(
      xmin = -5, xmax = 5,
      ymin = -95, ymax = 45
    ),
    "-90 <= ymin"
  )
})

test_that("depth validation accepts valid depth ranges", {
  expect_invisible(
    sharkabc3d:::.copernicus_validate_depth(
      depth_min = 0,
      depth_max = 100
    )
  )

  expect_invisible(
    sharkabc3d:::.copernicus_validate_depth(
      depth_min = NULL,
      depth_max = NULL
    )
  )
})

test_that("depth validation rejects invalid depth ranges", {
  expect_error(
    sharkabc3d:::.copernicus_validate_depth(
      depth_min = -1,
      depth_max = 100
    ),
    "must be >= 0"
  )

  expect_error(
    sharkabc3d:::.copernicus_validate_depth(
      depth_min = 100,
      depth_max = 50
    ),
    "depth_min.*depth_max"
  )
})

test_that("datetime parsing returns UTC ISO datetimes", {
  expect_identical(
    sharkabc3d:::.copernicus_datetime(
      as.Date("2020-01-01")
    ),
    "2020-01-01T00:00:00"
  )

  expect_identical(
    sharkabc3d:::.copernicus_datetime(
      "2020-01-01 12:30:00"
    ),
    "2020-01-01T12:30:00"
  )

  expect_identical(
    sharkabc3d:::.copernicus_datetime(
      "2020-01-01T12:30:00Z"
    ),
    "2020-01-01T12:30:00"
  )
})

test_that("datetime parsing rejects invalid input", {
  expect_error(
    sharkabc3d:::.copernicus_datetime(
      "not-a-date"
    ),
    "Could not parse"
  )

  expect_error(
    sharkabc3d:::.copernicus_datetime(
      c("2020-01-01", "2020-01-02")
    ),
    "must be a Date, POSIXt, or single character datetime"
  )
})

test_that("time ranges are normalized correctly", {
  out <- sharkabc3d:::.copernicus_time_range(
    start_datetime = "2020-01-01",
    end_datetime = NULL
  )

  expect_identical(out$start, "2020-01-01T00:00:00")
  expect_identical(out$end, "2020-01-01T00:00:00")

  out <- sharkabc3d:::.copernicus_time_range(
    start_datetime = "2020-01-01T12:00:00",
    end_datetime = "2020-01-02T12:00:00"
  )

  expect_identical(out$start, "2020-01-01T12:00:00")
  expect_identical(out$end, "2020-01-02T12:00:00")
})

test_that("time ranges reject invalid combinations", {
  expect_error(
    sharkabc3d:::.copernicus_time_range(
      start_datetime = NULL,
      end_datetime = "2020-01-02"
    ),
    "start_datetime"
  )

  expect_error(
    sharkabc3d:::.copernicus_time_range(
      start_datetime = "2020-01-02",
      end_datetime = "2020-01-01"
    ),
    "equal to or later"
  )
})

test_that("CERRA datasets are detected correctly", {
  expect_true(
    sharkabc3d:::.copernicus_is_cerra(
      "reanalysis-cerra-single-levels"
    )
  )

  expect_true(
    sharkabc3d:::.copernicus_is_cerra(
      "reanalysis-cerra-pressure-levels"
    )
  )

  expect_false(
    sharkabc3d:::.copernicus_is_cerra(
      "reanalysis-era5-single-levels"
    )
  )
})

test_that("copernicus_load validates source and dataset_id", {
  expect_error(
    copernicus_load(
      source = "invalid",
      dataset_id = "example"
    )
  )

  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = ""
    ),
    "dataset_id"
  )
})

test_that("copernicus_load validates variables and request", {
  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      variables = ""
    ),
    "variables"
  )

  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      request = list("reanalysis")
    ),
    "named list"
  )
})

test_that("copernicus_load validates force and quiet", {
  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      force = NA
    ),
    "force"
  )

  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      quiet = NA
    ),
    "quiet"
  )
})
