#' test utility for creating small test netCDFs 
make_file <- function(path, values_by_depth) {
  layers <- lapply(names(values_by_depth), function(dn) {
    r <- terra::rast(nrows = 2, ncols = 2,
                      xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                      crs = "EPSG:4326")
    terra::values(r) <- values_by_depth[[dn]]
    names(r) <- dn
    r
  })
  stk <- terra::rast(layers)
  terra::writeCDF(stk, path, overwrite = TRUE, split = TRUE)
}

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

test_that("woa_download errors on unknown variable", {
  expect_error(
    woa_download("nonsense_variable", period = "annual"),
    "Unknown variable"
  )
})

test_that("woa_download errors on unknown resolution", {
  expect_error(
    woa_download("temperature", resolution = "7"),
    "Unknown resolution"
  )
})

test_that("woa_download errors on invalid period", {
  expect_error(
    woa_download("temperature", period = "yearly"),
    "period must be"
  )
  expect_error(
    woa_download("temperature", period = 99),
    "Numeric period"
  )
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
      list(headers = charToRaw(
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

test_that("woa_load_nc errors when file missing", {
  expect_error(woa_load_nc("/no/such/file.nc"), "File not found")
})

test_that(paste("woa_load_nc catches layer names that do not match WOA",
                "one-letter code for oceanographic variables"), {
  # Build two tiny synthetic NetCDFs matching the WOA layer-naming convention.
  tmp <- tempfile("woa_mon_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  f1 <- file.path(tmp, "m01.nc")
  make_file(f1, list(`do_sea_depth=0` = c(1, 2, 3, 4),
                     `do_sea_depth=100` = c(5, 6, 7, 8)))

  expect_error(
    woa_load_nc(file.path(tmp, "m01.nc"), field = "sea"),
    paste("check raster layer names. For WOA data, oceanographic variable",
          "one-letter code must be one of: t, s, o, O, A, n, p, i, I.")
  )
})

test_that(paste("woa_load_nc catches layer names that do not match WOA code",
                "for abbreviated fields"), {
  # Build two tiny synthetic NetCDFs matching the WOA layer-naming convention.
  tmp <- tempfile("woa_mon_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  f1 <- file.path(tmp, "m01.nc")
  make_file(f1, list(`o_inc_depth=0` = c(1, 2, 3, 4),
                     `o_inc_depth=100` = c(5, 6, 7, 8)))

  expect_error(
    woa_load_nc(file.path(tmp, "m01.nc")),
    paste("check raster layer names. For WOA data, field code must be one",
          "of: an, mn, dd, sd, se, oa, gp, sdo, sea.")
  )
})

test_that(paste("woa_load_nc provides descriptive error if standard field",
                "is missing from input netCDF"), {
  # Build a tiny synthetic NetCDF carrying only "an" layers.
  tmp <- tempfile("woa_mon_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  f1 <- file.path(tmp, "m01.nc")
  make_file(f1, list(`t_an_depth=0` = c(1, 2, 3, 4),
                     `t_an_depth=100` = c(5, 6, 7, 8)))

  # "sea" is a valid WOA field code, so it clears the `field` argument check,
  # but this file has no sea layers. Without the guard the empty subset fails
  # inside terra with "is.numeric(i) is not TRUE".
  expect_error(
    woa_load_nc(f1, field = "sea"),
    "No layers for field 'sea' in m01.nc. Fields present in this file: an.",
    fixed = TRUE
  )

  # A field that is present still loads normally
  expect_named(
    woa_load_nc(f1, field = "an"),
    c("t_an_depth=0", "t_an_depth=100")
  )
})