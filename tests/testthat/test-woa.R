test_that("woa_cache_dir returns a writable path", {
  path <- woa_cache_dir()
  expect_true(dir.exists(path))
  test_file <- file.path(path, ".write_test")
  writeLines("ok", test_file)
  expect_true(file.exists(test_file))
  unlink(test_file)
})

test_that("woa_cache_dir is stable across calls", {
  expect_identical(woa_cache_dir(), woa_cache_dir())
})

test_that("woa_download errors on an unknown variable, resolution or period", {
  # `args` go to woa_download().
  cases <- list(
    list(label = "unknown variable",
         args = list("nonsense_variable", period = "annual"),
         pattern = "Unknown variable"),
    list(label = "unknown resolution",
         args = list("temperature", resolution = "7"),
         pattern = "Unknown resolution"),
    list(label = "unknown period name",
         args = list("temperature", period = "yearly"),
         pattern = "period must be"),
    list(label = "unknown period number",
         args = list("temperature", period = 99),
         pattern = "Numeric period")
  )
  for (case in cases) {
    expect_error(do.call(woa_download, case$args), case$pattern,
                 label = case$label)
  }
})

test_that("woa_download caches files and skips re-download", {
  # curl is mocked, so no request reaches NCEI.
  downloads <- 0L
  testthat::local_mocked_bindings(
    curl_download = function(url, destfile, ...) {
      downloads <<- downloads + 1L
      writeBin(raw(1024), destfile)
      invisible(destfile)
    },
    curl_fetch_memory = function(url, handle, ...) {
      list(status_code = 200L, headers = charToRaw(
        "HTTP/1.1 200 OK\r\nContent-Length: 1024\r\n\r\n"
      ))
    },
    .package = "curl"
  )
  tmp <- tempfile("woa_cache_")
  on.exit(unlink(tmp, recursive = TRUE))

  expect_message(
    path <- woa_download("temperature", period = "annual", resolution = "1",
                         output_dir = tmp),
    "Downloading (0.0 MB): woa23_decav_t00_01.nc", fixed = TRUE
  )
  expect_identical(basename(path), "woa23_decav_t00_01.nc")
  expect_true(file.exists(path))
  expect_identical(downloads, 1L)

  # The second call finds the file cached at the remote size.
  expect_message(
    path2 <- woa_download("temperature", period = "annual", resolution = "1",
                          output_dir = tmp),
    "Cached: woa23_decav_t00_01.nc", fixed = TRUE
  )
  expect_identical(path2, path)
  expect_identical(downloads, 1L)
})

test_that("woa_download builds the NCEI URL for each resolution", {
  # Checked against the NCEI server with HEAD requests in #46.
  base <- paste0("https://www.ncei.noaa.gov/thredds-ocean/fileServer/",
                 "woa23/DATA/temperature/netcdf/decav/")
  cases <- list(
    list(resolution = "0.25", url = "0.25/woa23_decav_t00_04.nc"),
    list(resolution = 0.25, url = "0.25/woa23_decav_t00_04.nc"),
    list(resolution = "1", url = "1.00/woa23_decav_t00_01.nc"),
    list(resolution = "1.00", url = "1.00/woa23_decav_t00_01.nc"),
    list(resolution = 1, url = "1.00/woa23_decav_t00_01.nc"),
    list(resolution = "5", url = "5.00/woa23_decav_t00_5d.nc"),
    list(resolution = 5, url = "5.00/woa23_decav_t00_5d.nc")
  )
  url <- NULL
  testthat::local_mocked_bindings(
    curl_download = function(url, destfile, ...) {
      url <<- url
      writeBin(raw(1024), destfile)
      invisible(destfile)
    },
    .package = "curl"
  )
  for (case in cases) {
    tmp <- tempfile("woa_url_")
    woa_download("temperature", resolution = case$resolution,
                 output_dir = tmp, quiet = TRUE)
    unlink(tmp, recursive = TRUE)
    expect_identical(url, paste0(base, case$url),
                     label = paste("resolution", case$resolution))
  }
})

test_that("woa_download reports an unknown size when the server says 404", {
  # A 404 page has its own Content-Length, which isn't the file's size.
  testthat::local_mocked_bindings(
    curl_download = function(url, destfile, ...) {
      writeBin(raw(1024), destfile)
      invisible(destfile)
    },
    curl_fetch_memory = function(url, handle, ...) {
      list(status_code = 404L, headers = charToRaw(
        "HTTP/1.1 404 Not Found\r\nContent-Length: 431\r\n\r\n"
      ))
    },
    .package = "curl"
  )
  tmp <- tempfile("woa_404_")
  on.exit(unlink(tmp, recursive = TRUE))
  expect_message(
    woa_download("temperature", resolution = "1", output_dir = tmp),
    "Downloading (unknown size): woa23_decav_t00_01.nc", fixed = TRUE
  )
})
