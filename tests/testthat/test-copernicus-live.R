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


test_that("Copernicus Marine multi-file temporal summaries work", {
  skip_if_not_installed("reticulate")
  skip_if_not_installed("ncdf4")

  output_dir <- file.path(test_root, "marine_summary")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Download three independent daily NetCDF files so that the live test checks
  # both real Copernicus files and aggregation across multiple file_paths.
  dates <- as.Date("2020-01-01") + 0:2

  input_paths <- vapply(dates, function(date) {
    date_chr <- format(date, "%Y-%m-%d")
    filename <- paste0("marine_thetao_", date_chr, ".nc")

    copernicus_load(
      source = "marine",
      dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
      variables = "thetao",
      start_datetime = paste0(date_chr, "T12:00:00"),
      end_datetime = paste0(date_chr, "T12:00:00"),
      xmin = XMIN, xmax = XMAX,
      ymin = YMIN, ymax = YMAX,
      depth_min = 0, depth_max = 10,
      output_dir = output_dir,
      filename = filename,
      force = TRUE,
      quiet = TRUE
    )
  }, character(1))

  expect_length(input_paths, 3)
  expect_true(all(file.exists(input_paths)))

  summary_paths <- copernicus_summarise(
    file_paths = input_paths,
    fun = c("mean", "min", "max", "sd"),
    output_dir = output_dir,
    filename = "marine_thetao_summary.nc",
    force = TRUE,
    quiet = TRUE
  )

  expect_named(summary_paths, c("mean", "min", "max", "sd"))
  expect_true(all(file.exists(summary_paths)))
  expect_true(all(file.info(summary_paths)$size > 0))

  # Independently read every real time slice with ncdf4. The expected summaries
  # are calculated here rather than reusing any sharkabc3d summarising helper.
  read_thetao_slices <- function(path) {
    nc <- ncdf4::nc_open(path)
    on.exit(ncdf4::nc_close(nc))

    var <- nc$var[["thetao"]]
    dim_names <- vapply(var$dim, `[[`, character(1), "name")
    time_idx <- match("time", tolower(dim_names))

    if (is.na(time_idx)) {
      stop("Live Marine file does not contain a time dimension.")
    }

    out_dim <- vapply(var$dim[-time_idx], `[[`, numeric(1), "len")
    n_time <- var$dim[[time_idx]]$len
    start <- rep(1, length(dim_names))
    count <- vapply(var$dim, `[[`, numeric(1), "len")
    count[time_idx] <- 1

    lapply(seq_len(n_time), function(tt) {
      start[time_idx] <- tt
      x <- ncdf4::ncvar_get(
        nc, "thetao",
        start = start,
        count = count,
        collapse_degen = FALSE
      )
      array(as.numeric(x), dim = out_dim)
    })
  }

  slices <- unlist(lapply(input_paths, read_thetao_slices), recursive = FALSE)
  expect_gte(length(slices), 3)

  slice_dim <- dim(slices[[1]])
  values <- vapply(slices, as.numeric, numeric(prod(slice_dim)))

  # Each row represents one fixed longitude × latitude × depth combination
  # across time. Summaries therefore collapse only the temporal dimension.
  expect_equal(nrow(values), prod(slice_dim))
  expect_equal(ncol(values), length(slices))

  # Independent reference summaries. Fully missing cells remain NA rather than
  # becoming Inf/-Inf for min/max.
  safe_mean <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  }

  safe_min <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    min(x, na.rm = TRUE)
  }

  safe_max <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  }

  safe_sd <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    stats::sd(x)
  }

  expected <- list(
    mean = array(apply(values, 1, safe_mean), dim = slice_dim),
    min = array(apply(values, 1, safe_min), dim = slice_dim),
    max = array(apply(values, 1, safe_max), dim = slice_dim),
    sd = array(apply(values, 1, safe_sd), dim = slice_dim)
  )

  # Check every statistic numerically and ensure the temporal dimension has
  # disappeared while longitude/latitude/depth remain.
  for (stat in names(summary_paths)) {
    nc <- ncdf4::nc_open(summary_paths[[stat]])

    expect_false(any(tolower(names(nc$dim)) %in% c("time", "valid_time")))
    expect_true(all(c("longitude", "latitude", "depth") %in% names(nc$dim)))
    expect_true("thetao" %in% names(nc$var))

    observed <- ncdf4::ncvar_get(nc, "thetao")
    ncdf4::nc_close(nc)

    expect_equal(
      observed,
      expected[[stat]],
      tolerance = 1e-5,
      info = paste("Incorrect live Marine", stat, "summary")
    )
  }
})


test_that("Copernicus Marine temporal summaries respect datetime filters", {
  skip_if_not_installed("reticulate")
  skip_if_not_installed("ncdf4")

  output_dir <- file.path(test_root, "marine_summary_filtered")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  dates <- as.Date("2020-01-01") + 0:2
  input_paths <- vapply(dates, function(date) {
    date_chr <- format(date, "%Y-%m-%d")

    copernicus_load(
      source = "marine",
      dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
      variables = "thetao",
      start_datetime = paste0(date_chr, "T12:00:00"),
      end_datetime = paste0(date_chr, "T12:00:00"),
      xmin = XMIN, xmax = XMAX,
      ymin = YMIN, ymax = YMAX,
      depth_min = 0, depth_max = 10,
      output_dir = output_dir,
      filename = paste0("marine_thetao_", date_chr, ".nc"),
      force = TRUE,
      quiet = TRUE
    )
  }, character(1))

  summary_path <- copernicus_summarise(
    file_paths = input_paths,
    fun = "mean",
    start_datetime = "2020-01-02T00:00:00",
    end_datetime = "2020-01-03T23:59:59",
    output_dir = output_dir,
    filename = "marine_thetao_filtered_mean.nc",
    force = TRUE,
    quiet = TRUE
  )

  expect_named(summary_path, "mean")
  expect_true(file.exists(summary_path))

  nc_out <- ncdf4::nc_open(summary_path)
  on.exit(ncdf4::nc_close(nc_out), add = TRUE)

  expect_equal(ncdf4::ncatt_get(nc_out, 0, "summary_time_steps")$value, 2)
  expect_match(
    ncdf4::ncatt_get(nc_out, 0, "summary_start_datetime")$value,
    "^2020-01-02"
  )
  expect_match(
    ncdf4::ncatt_get(nc_out, 0, "summary_end_datetime")$value,
    "^2020-01-03"
  )
})
