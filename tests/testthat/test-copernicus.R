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

test_that("copernicus_summarise validates inputs", {
  skip_if_not_installed("ncdf4")
  
  expect_error(
    copernicus_summarise(character()),
    "non-empty character vector"
  )
  
  expect_error(
    copernicus_summarise("does_not_exist.nc"),
    "not found"
  )
  
  tmp <- tempfile(fileext = ".txt")
  writeLines("x", tmp)
  on.exit(unlink(tmp), add = TRUE)
  
  expect_error(
    copernicus_summarise(tmp),
    "netCDF"
  )
})

test_that("copernicus_summarise rejects unsupported functions", {
  skip_if_not_installed("ncdf4")
  
  expect_error(
    copernicus_summarise(
      "does_not_matter.nc",
      fun = "median"
    )
  )
})

test_that("copernicus_summarise preserves lon lat depth and summarises time", {
  skip_if_not_installed("ncdf4")
  
  tmp <- tempfile("copernicus_summary_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  
  input <- file.path(tmp, "synthetic.nc")
  
  lon <- ncdf4::ncdim_def("longitude", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("latitude", "degrees_north", c(40, 41))
  depth <- ncdf4::ncdim_def("depth", "m", c(0, 100))
  time <- ncdf4::ncdim_def(
    "time",
    "days since 2000-01-01 00:00:00",
    0:2,
    unlim = TRUE
  )
  
  thetao <- ncdf4::ncvar_def(
    "thetao",
    "degrees_C",
    list(lon, lat, depth, time),
    missval = -9999,
    prec = "double"
  )
  
  nc <- ncdf4::nc_create(input, thetao)
  
  values <- array(
    c(
      # time 1
      1, 2, 3, 4, 5, 6, 7, 8,
      # time 2
      2, 3, 4, 5, 6, 7, 8, 9,
      # time 3
      3, 4, 5, 6, 7, 8, 9, 10
    ),
    dim = c(2, 2, 2, 3)
  )
  
  ncdf4::ncvar_put(nc, "thetao", values)
  ncdf4::nc_close(nc)
  
  out <- copernicus_summarise(
    input,
    fun = c("mean", "min", "max", "sd"),
    output_dir = tmp,
    filename = "summary.nc",
    force = TRUE,
    quiet = TRUE
  )
  
  expect_named(out, c("mean", "min", "max", "sd"))
  expect_true(all(file.exists(out)))
  
  # Mean
  nc_mean <- ncdf4::nc_open(out["mean"])
  on.exit(ncdf4::nc_close(nc_mean), add = TRUE)
  
  expect_true(all(c("longitude", "latitude", "depth") %in% names(nc_mean$dim)))
  expect_false("time" %in% names(nc_mean$dim))
  
  mean_values <- ncdf4::ncvar_get(nc_mean, "thetao")
  
  expected_mean <- apply(values, c(1, 2, 3), mean)
  
  expect_equal(mean_values, expected_mean)
  
  # Min
  nc_min <- ncdf4::nc_open(out["min"])
  min_values <- ncdf4::ncvar_get(nc_min, "thetao")
  ncdf4::nc_close(nc_min)
  
  expect_equal(min_values, apply(values, c(1, 2, 3), min))
  
  # Max
  nc_max <- ncdf4::nc_open(out["max"])
  max_values <- ncdf4::ncvar_get(nc_max, "thetao")
  ncdf4::nc_close(nc_max)
  
  expect_equal(max_values, apply(values, c(1, 2, 3), max))
  
  # SD
  nc_sd <- ncdf4::nc_open(out["sd"])
  sd_values <- ncdf4::ncvar_get(nc_sd, "thetao")
  ncdf4::nc_close(nc_sd)
  
  expect_equal(
    sd_values,
    apply(values, c(1, 2, 3), sd),
    tolerance = 1e-10
  )
}) 