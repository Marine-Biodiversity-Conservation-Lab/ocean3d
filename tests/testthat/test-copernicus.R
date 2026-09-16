# ------------------------------------------------------------------------------
# Offline unit tests for Copernicus helpers
#
# No network access, authentication, or hard-coded package name should be
# required by this file. Real service tests belong in test-copernicus-live.R.
# ------------------------------------------------------------------------------

test_that("bbox validation accepts valid bounds", {
  expect_invisible(.copernicus_validate_bbox(-5, 5, 35, 45))
  expect_invisible(.copernicus_validate_bbox(NULL, NULL, NULL, NULL))
})

test_that("bbox validation rejects invalid bounds", {
  expect_error(.copernicus_validate_bbox(-5, 5, 35, NULL), "must be supplied together")
  expect_error(.copernicus_validate_bbox(5, -5, 35, 45), "xmin < xmax")
  expect_error(.copernicus_validate_bbox(-5, 5, 45, 35), "ymin < ymax")
  expect_error(.copernicus_validate_bbox(-5, 5, -95, 45), "-90 <= ymin")
})

test_that("depth validation works", {
  expect_invisible(.copernicus_validate_depth(0, 100))
  expect_invisible(.copernicus_validate_depth(NULL, NULL))
  expect_error(.copernicus_validate_depth(-1, 100), "must be >= 0")
  expect_error(.copernicus_validate_depth(100, 50), "depth_min.*depth_max")
})

test_that("datetime parsing returns normalized ISO datetimes", {
  expect_identical(.copernicus_datetime(as.Date("2020-01-01")), "2020-01-01T00:00:00")
  expect_identical(.copernicus_datetime("2020-01-01 12:30:00"), "2020-01-01T12:30:00")
  expect_identical(.copernicus_datetime("2020-01-01T12:30:00Z"), "2020-01-01T12:30:00")
})

test_that("datetime parsing rejects invalid input", {
  expect_error(.copernicus_datetime("not-a-date"), "Could not parse")
  expect_error(
    .copernicus_datetime(c("2020-01-01", "2020-01-02")),
    "must be a Date, POSIXt, or single character datetime"
  )
})

test_that("time ranges are normalized and validated", {
  out <- .copernicus_time_range("2020-01-01", NULL)
  expect_identical(out$start, "2020-01-01T00:00:00")
  expect_identical(out$end, "2020-01-01T00:00:00")
  
  out <- .copernicus_time_range(
    "2020-01-01T12:00:00",
    "2020-01-02T12:00:00"
  )
  expect_identical(out$start, "2020-01-01T12:00:00")
  expect_identical(out$end, "2020-01-02T12:00:00")
  
  expect_error(
    .copernicus_time_range(NULL, "2020-01-02"),
    "start_datetime"
  )
  expect_error(
    .copernicus_time_range("2020-01-02", "2020-01-01"),
    "equal to or later"
  )
})

test_that("CERRA datasets are detected correctly", {
  expect_true(.copernicus_is_cerra("reanalysis-cerra-single-levels"))
  expect_true(.copernicus_is_cerra("reanalysis-cerra-pressure-levels"))
  expect_false(.copernicus_is_cerra("reanalysis-era5-single-levels"))
})

# ------------------------------------------------------------------------------
# copernicus_load() validation
# ------------------------------------------------------------------------------

test_that("copernicus_load validates source and dataset_id", {
  expect_error(copernicus_load(source = "invalid", dataset_id = "example"))
  expect_error(
    copernicus_load(source = "marine", dataset_id = ""),
    "dataset_id"
  )
})

test_that("Marine requests require explicit variables", {
  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      variables = NULL
    ),
    "variables"
  )
  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      variables = ""
    ),
    "variables"
  )
})

test_that("copernicus_load validates request, force, quiet, and compression", {
  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      variables = "thetao",
      request = list("reanalysis")
    ),
    "named list"
  )
  
  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      variables = "thetao",
      force = NA
    ),
    "force"
  )
  
  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      variables = "thetao",
      quiet = NA
    ),
    "quiet"
  )
  
  expect_error(
    copernicus_load(
      source = "marine",
      dataset_id = "example",
      variables = "thetao",
      compression = 10
    ),
    "compression"
  )
})

# ------------------------------------------------------------------------------
# File splitting and organization
# ------------------------------------------------------------------------------

test_that("valid file structures are accepted", {
  expect_invisible(
    .copernicus_validate_file_structure(
      split_by = "day",
      organize_by = c("service", "dataset", "variable", "year", "month", "day"),
      variables = c("uo", "vo")
    )
  )
  
  expect_invisible(
    .copernicus_validate_file_structure(
      split_by = "week",
      organize_by = c("variable", "year", "week"),
      variables = "thetao"
    )
  )
  
  expect_invisible(
    .copernicus_validate_file_structure(
      split_by = "season",
      organize_by = c("variable", "year", "season"),
      variables = "thetao"
    )
  )
  
  expect_invisible(
    .copernicus_validate_file_structure(
      split_by = NULL,
      organize_by = c("service", "dataset", "variable"),
      variables = c("uo", "vo")
    )
  )
})

test_that("split_by controls time only", {
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "variable",
      variables = c("uo", "vo")
    ),
    "Variables do not need to be included"
  )
})

test_that("invalid organization structures are rejected", {
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "day",
      organize_by = c("year", "banana", "day"),
      variables = "thetao"
    ),
    "Unsupported `organize_by`"
  )
  
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "day",
      organize_by = c("year", "day", "day"),
      variables = "thetao"
    ),
    "duplicated"
  )
  
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "day",
      organize_by = c("year", "month", "week", "day"),
      variables = "thetao"
    ),
    "cannot combine"
  )
  
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "day",
      organize_by = c("day", "year"),
      variables = "thetao"
    ),
    "hierarchical order"
  )
  
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "month",
      organize_by = c("year", "month", "day"),
      variables = "thetao"
    ),
    "finer than or incompatible"
  )
  
  expect_error(
    .copernicus_validate_file_structure(
      split_by = NULL,
      organize_by = c("year", "day"),
      variables = "thetao"
    ),
    "split_by` is NULL"
  )
})

test_that("variable folders require requested variables", {
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "day",
      organize_by = c("variable", "year", "day"),
      variables = NULL
    ),
    "variables` was not supplied"
  )
})

test_that("hemisphere and concurrent_processes are validated", {
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "season",
      variables = "thetao",
      hemisphere = "east"
    ),
    "north.*south"
  )
  
  expect_error(
    .copernicus_validate_file_structure(
      concurrent_processes = 0,
      variables = c("uo", "vo")
    ),
    "integer >= 1"
  )
  
  expect_error(
    .copernicus_validate_file_structure(
      concurrent_processes = 2,
      variables = "thetao"
    ),
    "only one output file"
  )
  
  expect_error(
    .copernicus_validate_file_structure(
      split_by = "week",
      concurrent_processes = 2,
      variables = "thetao"
    ),
    "cannot currently be used"
  )
  
  expect_invisible(
    .copernicus_validate_file_structure(
      split_by = "day",
      concurrent_processes = 2,
      variables = "thetao"
    )
  )
})

# ------------------------------------------------------------------------------
# ISO weeks and meteorological seasons
# ------------------------------------------------------------------------------

test_that("ISO week helper handles ordinary and boundary dates", {
  info <- .copernicus_iso_week(as.Date("2026-07-15"))
  expect_identical(info$year, 2026L)
  expect_identical(info$week, 29L)
  expect_identical(info$label, "week_29")
  expect_identical(info$start, as.Date("2026-07-13"))
  
  boundary <- .copernicus_iso_week(as.Date("2021-01-01"))
  expect_identical(boundary$year, 2020L)
  expect_identical(boundary$week, 53L)
  expect_identical(boundary$label, "week_53")
})

test_that("meteorological seasons use hemisphere labels and season-year", {
  north <- .copernicus_season(as.Date("2025-12-15"), "north")
  south <- .copernicus_season(as.Date("2025-12-15"), "south")
  
  expect_identical(north$year, 2026L)
  expect_identical(north$block, "DJF")
  expect_identical(north$label, "winter_DJF")
  
  expect_identical(south$year, 2026L)
  expect_identical(south$block, "DJF")
  expect_identical(south$label, "summer_DJF")
  
  expect_identical(
    .copernicus_season(as.Date("2026-04-15"), "north")$label,
    "spring_MAM"
  )
  expect_identical(
    .copernicus_season(as.Date("2026-07-15"), "north")$label,
    "summer_JJA"
  )
  expect_identical(
    .copernicus_season(as.Date("2026-10-15"), "north")$label,
    "autumn_SON"
  )
})

test_that("custom week periods are clipped to the request", {
  periods <- .copernicus_custom_split_periods(
    start_datetime = "2026-07-15T00:00:00",
    end_datetime = "2026-07-21T23:59:59",
    split_by = "week",
    hemisphere = "north"
  )
  
  expect_length(periods, 2)
  expect_identical(periods[[1]]$label, "2026_week_29")
  expect_identical(periods[[1]]$start, "2026-07-15T00:00:00")
  expect_identical(periods[[1]]$end, "2026-07-19T23:59:59")
  expect_identical(periods[[2]]$label, "2026_week_30")
  expect_identical(periods[[2]]$start, "2026-07-20T00:00:00")
  expect_identical(periods[[2]]$end, "2026-07-21T23:59:59")
})

test_that("custom week periods use ISO year across New Year", {
  periods <- .copernicus_custom_split_periods(
    start_datetime = "2020-12-31T00:00:00",
    end_datetime = "2021-01-05T23:59:59",
    split_by = "week",
    hemisphere = "north"
  )
  
  expect_length(periods, 2)
  expect_identical(periods[[1]]$label, "2020_week_53")
  expect_identical(periods[[2]]$label, "2021_week_01")
})

test_that("custom season periods handle DJF across calendar years", {
  periods <- .copernicus_custom_split_periods(
    start_datetime = "2025-12-15T00:00:00",
    end_datetime = "2026-03-10T12:00:00",
    split_by = "season",
    hemisphere = "north"
  )
  
  expect_length(periods, 2)
  expect_identical(periods[[1]]$label, "2026_winter_DJF")
  expect_identical(periods[[1]]$start, "2025-12-15T00:00:00")
  expect_identical(periods[[1]]$end, "2026-02-28T23:59:59")
  expect_identical(periods[[2]]$label, "2026_spring_MAM")
  expect_identical(periods[[2]]$start, "2026-03-01T00:00:00")
  expect_identical(periods[[2]]$end, "2026-03-10T12:00:00")
})

test_that("custom splitting rejects unsupported modes or missing dates", {
  expect_error(
    .copernicus_custom_split_periods(
      "2026-01-01T00:00:00",
      "2026-01-02T00:00:00",
      "day",
      "north"
    ),
    "only used"
  )
  
  expect_error(
    .copernicus_custom_split_periods(
      NULL,
      NULL,
      "week",
      "north"
    ),
    "required"
  )
})

# ------------------------------------------------------------------------------
# Filename helpers and Toolbox JSON parsers
# ------------------------------------------------------------------------------

test_that("safe file components and period filenames are deterministic", {
  expect_identical(
    .copernicus_safe_component("thetao / daily:* test"),
    "thetao_daily_test"
  )
  
  expect_error(
    .copernicus_safe_component("   "),
    "valid file or directory name"
  )
  
  expect_identical(
    .copernicus_period_filename(
      NULL,
      "example-dataset",
      "thetao",
      "2026_week_29"
    ),
    "example-dataset_thetao_2026_week_29.nc"
  )
  
  expect_identical(
    .copernicus_period_filename(
      "my output.nc",
      "example-dataset",
      "thetao",
      "summer/JJA"
    ),
    "my_output_summer_JJA.nc"
  )
})

test_that("dry-run parser extracts valid planning metadata", {
  output <- c(
    "INFO - Selected dataset version",
    "[",
    "  {",
    "    \"filename\": \"thetao.nc\",",
    "    \"variables\": [\"thetao\"],",
    "    \"coordinates_extent\": [",
    "      {\"coordinate_id\": \"longitude\", \"minimum\": -1, \"maximum\": 0},",
    "      {\"coordinate_id\": \"latitude\", \"minimum\": 38, \"maximum\": 39}",
    "    ],",
    "    \"data_transfer_size\": \"20.95 MB\"",
    "  }",
    "]"
  )
  
  plan <- .copernicus_parse_dry_run(output)
  
  expect_length(plan, 1)
  expect_identical(plan[[1]]$filename, "thetao.nc")
  expect_identical(plan[[1]]$variables, list("thetao"))
  expect_identical(plan[[1]]$data_transfer_size, "20.95 MB")
})

test_that("dry-run parser rejects malformed plans", {
  expect_error(
    .copernicus_parse_dry_run(character()),
    "returned no output"
  )
  
  expect_error(
    .copernicus_parse_dry_run("INFO only"),
    "Could not find the JSON"
  )
  
  multi_variable <- c(
    "[",
    "  {",
    "    \"filename\": \"bad.nc\",",
    "    \"variables\": [\"uo\", \"vo\"],",
    "    \"coordinates_extent\": []",
    "  }",
    "]"
  )
  
  expect_error(
    .copernicus_parse_dry_run(multi_variable),
    "exactly one environmental variable"
  )
})


# ------------------------------------------------------------------------------
# Marine backend selection
# ------------------------------------------------------------------------------

test_that("auto prefers standalone when executable is available", {
  
  fake_executable <- tempfile(
    "copernicusmarine_"
  )
  
  file.create(
    fake_executable
  )
  
  on.exit(
    unlink(
      fake_executable
    ),
    add = TRUE
  )
  
  with_mocked_bindings(
    expect_identical(
      .copernicus_select_marine_backend(
        "auto"
      ),
      "standalone"
    ),
    .copernicus_find_executable = function(path = NULL) {
      fake_executable
    }
  )
})

test_that("explicit standalone works when executable is available", {
  
  fake_executable <- tempfile(
    "copernicusmarine_"
  )
  
  file.create(
    fake_executable
  )
  
  on.exit(
    unlink(
      fake_executable
    ),
    add = TRUE
  )
  
  with_mocked_bindings(
    expect_identical(
      .copernicus_select_marine_backend(
        "standalone"
      ),
      "standalone"
    ),
    .copernicus_find_executable = function(path = NULL) {
      fake_executable
    }
  )
})

test_that("explicit standalone fails clearly when executable is unavailable", {
  
  with_mocked_bindings(
    expect_error(
      .copernicus_select_marine_backend(
        "standalone"
      ),
      "copernicus_setup|copernicus_status"
    ),
    .copernicus_find_executable = function(path = NULL) NULL
  )
})

test_that("Python backend can be selected when reticulate is installed", {
  skip_if_not_installed("reticulate")
  
  expect_identical(
    .copernicus_select_marine_backend("python"),
    "python"
  )
})

test_that("auto falls back to Python when standalone is unavailable", {
  skip_if_not_installed("reticulate")
  
  with_mocked_bindings(
    expect_identical(
      .copernicus_select_marine_backend(
        "auto"
      ),
      "python"
    ),
    .copernicus_find_executable = function(path = NULL) NULL
  )
})

# ------------------------------------------------------------------------------
# copernicus_summarise()
# ------------------------------------------------------------------------------

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

test_that("copernicus_summarise preserves space and summarises time", {
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
      1, 2, 3, 4, 5, 6, 7, 8,
      2, 3, 4, 5, 6, 7, 8, 9,
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
  
  nc_mean <- ncdf4::nc_open(out[["mean"]])
  expect_true(all(c("longitude", "latitude", "depth") %in% names(nc_mean$dim)))
  expect_false("time" %in% names(nc_mean$dim))
  mean_values <- ncdf4::ncvar_get(nc_mean, "thetao")
  ncdf4::nc_close(nc_mean)
  expect_equal(mean_values, apply(values, c(1, 2, 3), mean))
  
  nc_min <- ncdf4::nc_open(out[["min"]])
  min_values <- ncdf4::ncvar_get(nc_min, "thetao")
  ncdf4::nc_close(nc_min)
  expect_equal(min_values, apply(values, c(1, 2, 3), min))
  
  nc_max <- ncdf4::nc_open(out[["max"]])
  max_values <- ncdf4::ncvar_get(nc_max, "thetao")
  ncdf4::nc_close(nc_max)
  expect_equal(max_values, apply(values, c(1, 2, 3), max))
  
  nc_sd <- ncdf4::nc_open(out[["sd"]])
  sd_values <- ncdf4::ncvar_get(nc_sd, "thetao")
  ncdf4::nc_close(nc_sd)
  expect_equal(
    sd_values,
    apply(values, c(1, 2, 3), sd),
    tolerance = 1e-10
  )
})
