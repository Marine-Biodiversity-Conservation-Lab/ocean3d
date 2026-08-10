# ------------------------------------------------------------------------------
# Live integration tests for Copernicus services
#
# These tests perform small real downloads and therefore require:
# - internet access
# - valid Copernicus/ECMWF credentials
# - accepted provider licences/policies
#
# They are skipped by default. To run them:
#
# Sys.setenv(RUN_COPERNICUS_LIVE_TESTS = "true")
# testthat::test_file("tests/testthat/test-copernicus-live.R")
#
# Disable again with:
# Sys.unsetenv("RUN_COPERNICUS_LIVE_TESTS")
# ------------------------------------------------------------------------------

skip_if(
  Sys.getenv("RUN_COPERNICUS_LIVE_TESTS") != "true",
  "Live Copernicus tests are disabled"
)

skip_if_not_installed("terra")

XMIN <- 0
XMAX <- 1
YMIN <- 38
YMAX <- 39

test_root <- file.path(tempdir(), "sharkabc3d_copernicus_live_tests")
dir.create(test_root, recursive = TRUE, showWarnings = FALSE)

check_raster_download <- function(path) {
  expect_true(is.character(path))
  expect_length(path, 1)
  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 0)

  r <- terra::rast(path)

  expect_s4_class(r, "SpatRaster")
  expect_gt(terra::ncell(r), 0)
  expect_gt(terra::nlyr(r), 0)

  invisible(r)
}

test_that("Copernicus Marine download works", {
  skip_if_not_installed("reticulate")

  output_dir <- file.path(test_root, "marine")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  path <- copernicus_load(
    source = "marine",
    dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
    variables = "thetao",
    start_datetime = "2020-01-01T12:00:00",
    end_datetime = "2020-01-01T12:00:00",
    xmin = XMIN, xmax = XMAX,
    ymin = YMIN, ymax = YMAX,
    depth_min = 0, depth_max = 10,
    output_dir = output_dir,
    filename = "marine_thetao_test.nc",
    force = TRUE,
    quiet = TRUE
  )

  r <- check_raster_download(path)

  expect_true(any(grepl("thetao", names(r), ignore.case = TRUE)))
})

test_that("CDS ERA5 download works", {
  skip_if_not_installed("ecmwfr")

  output_dir <- file.path(test_root, "cds_era5")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  path <- copernicus_load(
    source = "cds",
    dataset_id = "reanalysis-era5-single-levels",
    variables = "2m_temperature",
    xmin = XMIN, xmax = XMAX,
    ymin = YMIN, ymax = YMAX,
    request = list(
      product_type = "reanalysis",
      year = "2020",
      month = "01",
      day = "01",
      time = "12:00",
      data_format = "netcdf",
      download_format = "unarchived"
    ),
    output_dir = output_dir,
    filename = "era5_2m_temperature_test.nc",
    force = TRUE,
    quiet = TRUE
  )

  check_raster_download(path)
})

test_that("CDS CERRA download and local crop work", {
  skip_if_not_installed("ecmwfr")

  output_dir <- file.path(test_root, "cds_cerra")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  path <- copernicus_load(
    source = "cds",
    dataset_id = "reanalysis-cerra-single-levels",
    variables = "2m_temperature",
    xmin = XMIN, xmax = XMAX,
    ymin = YMIN, ymax = YMAX,
    request = list(
      level_type = "surface_or_atmosphere",
      data_type = "reanalysis",
      product_type = "analysis",
      year = "2020",
      month = "01",
      day = "01",
      time = "12:00",
      data_format = "grib"
    ),
    output_dir = output_dir,
    filename = "cerra_2m_temperature_cropped.grib",
    force = TRUE,
    quiet = TRUE
  )

  r <- check_raster_download(path)

  expect_lt(terra::ncol(r), 1069)
  expect_lt(terra::nrow(r), 1069)
})

test_that("ADS CAMS download works", {
  skip_if_not_installed("ecmwfr")

  output_dir <- file.path(test_root, "ads_cams")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  path <- copernicus_load(
    source = "ads",
    dataset_id = "cams-global-reanalysis-eac4",
    variables = "particulate_matter_2.5um",
    xmin = XMIN, xmax = XMAX,
    ymin = YMIN, ymax = YMAX,
    request = list(
      date = "2020-01-01/2020-01-01",
      time = "12:00",
      data_format = "netcdf",
      download_format = "unarchived"
    ),
    output_dir = output_dir,
    filename = "cams_pm25_test.nc",
    force = TRUE,
    quiet = TRUE
  )

  check_raster_download(path)
})