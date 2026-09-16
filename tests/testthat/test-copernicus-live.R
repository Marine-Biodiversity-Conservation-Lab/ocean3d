# ------------------------------------------------------------------------------
# Live integration tests for Copernicus Marine
#
# These tests make REAL requests to Copernicus Marine. They are skipped by
# default and must never run during ordinary package checks.
#
# BASIC LIVE TESTS
# ----------------
# Small standalone-backend downloads covering:
#   1. one-variable download
#   2. multiple variables -> one NetCDF per variable
#   3. temporal splitting + directory organization
#   4. complete reuse when valid files already exist
#   5. partial resume when one expected file is missing
#
# Run with:
#
#   Sys.setenv(RUN_COPERNICUS_LIVE_TESTS = "true")
#   devtools::load_all()
#   testthat::test_file("tests/testthat/test-copernicus-live.R")
#
# Disable again with:
#
#   Sys.unsetenv("RUN_COPERNICUS_LIVE_TESTS")
#
#
# EXTENDED LIVE TESTS
# -------------------
# Week/season splitting and force/replacement checks make additional requests
# and are therefore opt-in:
#
#   Sys.setenv(RUN_COPERNICUS_EXTENDED_LIVE_TESTS = "true")
#
#
# PYTHON BACKEND TEST
# -------------------
# The optional Python backend is tested only when explicitly requested:
#
#   Sys.setenv(RUN_COPERNICUS_PYTHON_LIVE_TESTS = "true")
#
# This may provision a Python environment through reticulate.
#
# ------------------------------------------------------------------------------

skip_if(
  Sys.getenv("RUN_COPERNICUS_LIVE_TESTS") != "true",
  "Live Copernicus tests are disabled"
)

skip_if_not_installed("ncdf4")

# ------------------------------------------------------------------------------
# Shared configuration
# ------------------------------------------------------------------------------

DATASET_PHY <- "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m"
DATASET_CUR <- "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m"

XMIN <- -1
XMAX <- 0
YMIN <- 38
YMAX <- 39

START_DAY <- "2026-07-15T00:00:00"
END_DAY <- "2026-07-15T23:59:59"

test_root <- file.path(
  tempdir(),
  "copernicus_live_tests"
)

dir.create(
  test_root,
  recursive = TRUE,
  showWarnings = FALSE
)

# Use the standalone backend explicitly in the core live tests. This prevents
# an unavailable executable from silently turning these tests into Python tests.
standalone_executable <- .copernicus_find_executable()

skip_if(
  is.null(standalone_executable) ||
    !file.exists(standalone_executable),
  paste(
    "Standalone Copernicus Marine Toolbox not available.",
    "Run copernicus_setup() before the live tests."
  )
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

live_reset_dir <- function(name) {
  path <- file.path(
    test_root,
    name
  )
  
  if (dir.exists(path)) {
    unlink(
      path,
      recursive = TRUE,
      force = TRUE
    )
  }
  
  dir.create(
    path,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  path
}


live_check_netcdf <- function(
    path,
    expected_variable = NULL
) {
  
  expect_true(
    is.character(path)
  )
  
  expect_length(
    path,
    1
  )
  
  expect_true(
    file.exists(path)
  )
  
  expect_gt(
    file.info(path)$size,
    0
  )
  
  nc <- ncdf4::nc_open(
    path
  )
  
  on.exit(
    ncdf4::nc_close(nc),
    add = TRUE
  )
  
  expect_gt(
    length(nc$dim),
    0
  )
  
  expect_gt(
    length(nc$var),
    0
  )
  
  if (!is.null(expected_variable)) {
    expect_true(
      expected_variable %in% names(nc$var),
      info = paste(
        "Expected variable",
        expected_variable,
        "was not found in",
        path
      )
    )
  }
  
  invisible(TRUE)
}


live_file_variable <- function(path) {
  
  nc <- ncdf4::nc_open(
    path
  )
  
  on.exit(
    ncdf4::nc_close(nc),
    add = TRUE
  )
  
  candidates <- intersect(
    c(
      "thetao",
      "uo",
      "vo"
    ),
    names(nc$var)
  )
  
  if (length(candidates) != 1L) {
    testthat::fail(
      paste(
        "Live output should contain exactly one requested environmental variable:",
        path
      )
    )
  }
  
  candidates[[1]]
}


live_common_cur_request <- function(
    output_dir,
    force = FALSE,
    quiet = TRUE
) {
  
  copernicus_load(
    source = "marine",
    backend = "standalone",
    dataset_id = DATASET_CUR,
    variables = c(
      "uo",
      "vo"
    ),
    start_datetime = START_DAY,
    end_datetime = END_DAY,
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = output_dir,
    split_by = "day",
    organize_by = c(
      "service",
      "dataset",
      "variable",
      "year",
      "month",
      "day"
    ),
    force = force,
    quiet = quiet
  )
}


# ------------------------------------------------------------------------------
# BASIC LIVE TESTS
# ------------------------------------------------------------------------------

test_that("standalone Marine download returns a valid one-variable NetCDF", {
  
  output_dir <- live_reset_dir(
    "single_variable"
  )
  
  paths <- copernicus_load(
    source = "marine",
    backend = "standalone",
    dataset_id = DATASET_PHY,
    variables = "thetao",
    start_datetime = START_DAY,
    end_datetime = END_DAY,
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = output_dir,
    filename = "thetao_live.nc",
    force = TRUE,
    quiet = TRUE
  )
  
  expect_length(
    paths,
    1
  )
  
  live_check_netcdf(
    paths[[1]],
    expected_variable = "thetao"
  )
  
  expect_identical(
    live_file_variable(
      paths[[1]]
    ),
    "thetao"
  )
})


test_that("multiple Marine variables are split and organized independently", {
  
  output_dir <- live_reset_dir(
    "multiple_variables"
  )
  
  paths <- live_common_cur_request(
    output_dir = output_dir,
    force = TRUE
  )
  
  expect_length(
    paths,
    2
  )
  
  expect_true(
    all(
      file.exists(paths)
    )
  )
  
  variables <- vapply(
    paths,
    live_file_variable,
    character(1)
  )
  
  expect_setequal(
    variables,
    c(
      "uo",
      "vo"
    )
  )
  
  for (i in seq_along(paths)) {
    
    variable <- variables[[i]]
    
    live_check_netcdf(
      paths[[i]],
      expected_variable = variable
    )
    
    normalized <- normalizePath(
      paths[[i]],
      winslash = "/",
      mustWork = TRUE
    )
    
    expected_fragment <- paste0(
      "/marine/",
      DATASET_CUR,
      "/",
      variable,
      "/2026/07/15/"
    )
    
    expect_match(
      normalized,
      expected_fragment,
      fixed = TRUE
    )
  }
})


test_that("valid Marine outputs are reused when force is FALSE", {
  
  output_dir <- live_reset_dir(
    "reuse"
  )
  
  first <- live_common_cur_request(
    output_dir = output_dir,
    force = TRUE
  )
  
  expect_length(
    first,
    2
  )
  
  first <- sort(
    normalizePath(
      first,
      winslash = "/",
      mustWork = TRUE
    )
  )
  
  first_info <- file.info(
    first
  )
  
  second <- live_common_cur_request(
    output_dir = output_dir,
    force = FALSE
  )
  
  second <- sort(
    normalizePath(
      second,
      winslash = "/",
      mustWork = TRUE
    )
  )
  
  expect_identical(
    second,
    first
  )
  
  second_info <- file.info(
    second
  )
  
  expect_equal(
    as.numeric(second_info$size),
    as.numeric(first_info$size)
  )
  
  expect_equal(
    as.numeric(second_info$mtime),
    as.numeric(first_info$mtime)
  )
})


test_that("Marine partial resume downloads only the missing output", {
  
  output_dir <- live_reset_dir(
    "partial_resume"
  )
  
  first <- live_common_cur_request(
    output_dir = output_dir,
    force = TRUE
  )
  
  expect_length(
    first,
    2
  )
  
  variables <- vapply(
    first,
    live_file_variable,
    character(1)
  )
  
  keep_path <- first[
    variables == "uo"
  ][[1]]
  
  remove_path <- first[
    variables == "vo"
  ][[1]]
  
  keep_path <- normalizePath(
    keep_path,
    winslash = "/",
    mustWork = TRUE
  )
  
  keep_info_before <- file.info(
    keep_path
  )
  
  unlink(
    remove_path,
    force = TRUE
  )
  
  expect_false(
    file.exists(
      remove_path
    )
  )
  
  resumed <- live_common_cur_request(
    output_dir = output_dir,
    force = FALSE
  )
  
  expect_length(
    resumed,
    2
  )
  
  expect_true(
    all(
      file.exists(resumed)
    )
  )
  
  resumed_variables <- vapply(
    resumed,
    live_file_variable,
    character(1)
  )
  
  expect_setequal(
    resumed_variables,
    c(
      "uo",
      "vo"
    )
  )
  
  keep_path_after <- resumed[
    resumed_variables == "uo"
  ][[1]]
  
  keep_path_after <- normalizePath(
    keep_path_after,
    winslash = "/",
    mustWork = TRUE
  )
  
  expect_identical(
    keep_path_after,
    keep_path
  )
  
  keep_info_after <- file.info(
    keep_path_after
  )
  
  expect_equal(
    as.numeric(keep_info_after$size),
    as.numeric(keep_info_before$size)
  )
  
  expect_equal(
    as.numeric(keep_info_after$mtime),
    as.numeric(keep_info_before$mtime)
  )
  
  vo_path <- resumed[
    resumed_variables == "vo"
  ][[1]]
  
  live_check_netcdf(
    vo_path,
    expected_variable = "vo"
  )
})


test_that("Marine output can remain flat when organize_by is NULL", {
  
  output_dir <- live_reset_dir(
    "flat_output"
  )
  
  paths <- copernicus_load(
    source = "marine",
    backend = "standalone",
    dataset_id = DATASET_PHY,
    variables = "thetao",
    start_datetime = START_DAY,
    end_datetime = END_DAY,
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = output_dir,
    filename = "flat_thetao.nc",
    force = TRUE,
    quiet = TRUE
  )
  
  expect_length(
    paths,
    1
  )
  
  live_check_netcdf(
    paths[[1]],
    expected_variable = "thetao"
  )
  
  expect_identical(
    dirname(
      normalizePath(
        paths[[1]],
        winslash = "/",
        mustWork = TRUE
      )
    ),
    normalizePath(
      output_dir,
      winslash = "/",
      mustWork = TRUE
    )
  )
})


# ------------------------------------------------------------------------------
# EXTENDED LIVE TESTS
# ------------------------------------------------------------------------------

test_that("Marine ISO-week splitting works across a week boundary", {
  
  skip_if(
    Sys.getenv("RUN_COPERNICUS_EXTENDED_LIVE_TESTS") != "true",
    "Extended Copernicus live tests are disabled"
  )
  
  output_dir <- live_reset_dir(
    "week_split"
  )
  
  paths <- copernicus_load(
    source = "marine",
    backend = "standalone",
    dataset_id = DATASET_PHY,
    variables = "thetao",
    start_datetime = "2026-07-19T00:00:00",
    end_datetime = "2026-07-20T23:59:59",
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = output_dir,
    split_by = "week",
    organize_by = c(
      "variable",
      "year",
      "week"
    ),
    force = TRUE,
    quiet = TRUE
  )
  
  expect_length(
    paths,
    2
  )
  
  expect_true(
    all(
      file.exists(paths)
    )
  )
  
  normalized <- normalizePath(
    paths,
    winslash = "/",
    mustWork = TRUE
  )
  
  expect_true(
    any(
      grepl(
        "/2026/week_29/",
        normalized,
        fixed = TRUE
      )
    )
  )
  
  expect_true(
    any(
      grepl(
        "/2026/week_30/",
        normalized,
        fixed = TRUE
      )
    )
  )
  
  invisible(
    lapply(
      paths,
      live_check_netcdf,
      expected_variable = "thetao"
    )
  )
})


test_that("Marine meteorological-season splitting works across a boundary", {
  
  skip_if(
    Sys.getenv("RUN_COPERNICUS_EXTENDED_LIVE_TESTS") != "true",
    "Extended Copernicus live tests are disabled"
  )
  
  output_dir <- live_reset_dir(
    "season_split"
  )
  
  paths <- copernicus_load(
    source = "marine",
    backend = "standalone",
    dataset_id = DATASET_PHY,
    variables = "thetao",
    start_datetime = "2026-05-31T00:00:00",
    end_datetime = "2026-06-01T23:59:59",
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = output_dir,
    split_by = "season",
    organize_by = c(
      "variable",
      "year",
      "season"
    ),
    hemisphere = "north",
    force = TRUE,
    quiet = TRUE
  )
  
  expect_length(
    paths,
    2
  )
  
  normalized <- normalizePath(
    paths,
    winslash = "/",
    mustWork = TRUE
  )
  
  expect_true(
    any(
      grepl(
        "/2026/spring_MAM/",
        normalized,
        fixed = TRUE
      )
    )
  )
  
  expect_true(
    any(
      grepl(
        "/2026/summer_JJA/",
        normalized,
        fixed = TRUE
      )
    )
  )
  
  invisible(
    lapply(
      paths,
      live_check_netcdf,
      expected_variable = "thetao"
    )
  )
})


test_that("force TRUE replaces an existing standalone Marine output", {
  
  skip_if(
    Sys.getenv("RUN_COPERNICUS_EXTENDED_LIVE_TESTS") != "true",
    "Extended Copernicus live tests are disabled"
  )
  
  output_dir <- live_reset_dir(
    "force_replace"
  )
  
  first <- copernicus_load(
    source = "marine",
    backend = "standalone",
    dataset_id = DATASET_PHY,
    variables = "thetao",
    start_datetime = START_DAY,
    end_datetime = END_DAY,
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = output_dir,
    filename = "force_thetao.nc",
    force = TRUE,
    quiet = TRUE
  )
  
  first_path <- normalizePath(
    first[[1]],
    winslash = "/",
    mustWork = TRUE
  )
  
  first_mtime <- file.info(
    first_path
  )$mtime
  
  # Ensure common file systems can record a different modification time.
  Sys.sleep(
    1.2
  )
  
  second <- copernicus_load(
    source = "marine",
    backend = "standalone",
    dataset_id = DATASET_PHY,
    variables = "thetao",
    start_datetime = START_DAY,
    end_datetime = END_DAY,
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = output_dir,
    filename = "force_thetao.nc",
    force = TRUE,
    quiet = TRUE
  )
  
  second_path <- normalizePath(
    second[[1]],
    winslash = "/",
    mustWork = TRUE
  )
  
  expect_identical(
    second_path,
    first_path
  )
  
  expect_gt(
    as.numeric(
      file.info(
        second_path
      )$mtime
    ),
    as.numeric(
      first_mtime
    )
  )
})


# ------------------------------------------------------------------------------
# OPTIONAL PYTHON BACKEND TEST
# ------------------------------------------------------------------------------

test_that("Python and standalone backends return compatible tiny Marine data", {
  
  skip_if(
    Sys.getenv("RUN_COPERNICUS_PYTHON_LIVE_TESTS") != "true",
    "Python Copernicus live tests are disabled"
  )
  
  skip_if_not_installed(
    "reticulate"
  )
  
  standalone_dir <- live_reset_dir(
    "backend_compare_standalone"
  )
  
  python_dir <- live_reset_dir(
    "backend_compare_python"
  )
  
  standalone_path <- copernicus_load(
    source = "marine",
    backend = "standalone",
    dataset_id = DATASET_PHY,
    variables = "thetao",
    start_datetime = START_DAY,
    end_datetime = END_DAY,
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = standalone_dir,
    filename = "thetao_standalone.nc",
    force = TRUE,
    quiet = TRUE
  )
  
  python_path <- copernicus_load(
    source = "marine",
    backend = "python",
    dataset_id = DATASET_PHY,
    variables = "thetao",
    start_datetime = START_DAY,
    end_datetime = END_DAY,
    xmin = XMIN,
    xmax = XMAX,
    ymin = YMIN,
    ymax = YMAX,
    depth_min = 0,
    depth_max = 10,
    output_dir = python_dir,
    filename = "thetao_python.nc",
    force = TRUE,
    quiet = TRUE
  )
  
  live_check_netcdf(
    standalone_path[[1]],
    expected_variable = "thetao"
  )
  
  live_check_netcdf(
    python_path[[1]],
    expected_variable = "thetao"
  )
  
  nc_standalone <- ncdf4::nc_open(
    standalone_path[[1]]
  )
  
  on.exit(
    ncdf4::nc_close(
      nc_standalone
    ),
    add = TRUE
  )
  
  nc_python <- ncdf4::nc_open(
    python_path[[1]]
  )
  
  on.exit(
    ncdf4::nc_close(
      nc_python
    ),
    add = TRUE
  )
  
  standalone_values <- ncdf4::ncvar_get(
    nc_standalone,
    "thetao"
  )
  
  python_values <- ncdf4::ncvar_get(
    nc_python,
    "thetao"
  )
  
  expect_equal(
    dim(
      python_values
    ),
    dim(
      standalone_values
    )
  )
  
  expect_equal(
    python_values,
    standalone_values,
    tolerance = 1e-6
  )
})
