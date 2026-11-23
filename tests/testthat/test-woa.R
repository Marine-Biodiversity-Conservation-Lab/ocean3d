# A small WOA-style netCDF: each of `vars` on a 2 x 2 grid and two standard
# depths, with WOA's depth bounds and its one-step climatological time axis.
# Values are 1 to 8 for the first variable, 101 to 108 for the second, and so
# on, in storage order (lon fastest).
make_woa_file <- function(path, vars = "t_an") {
  dims <- list(
    lon = ncdf4::ncdim_def("lon", "degrees_east", c(0.5, 1.5)),
    lat = ncdf4::ncdim_def("lat", "degrees_north", c(0.5, 1.5)),
    depth = ncdf4::ncdim_def("depth", "meters", c(0, 100)),
    time = ncdf4::ncdim_def("time", "months since 0000-01-01 00:00:00", 6)
  )
  nv <- ncdf4::ncdim_def("nbounds", "", 1:2, create_dimvar = FALSE)
  defs <- c(
    lapply(vars, function(v) {
      ncdf4::ncvar_def(v, "degrees_celsius", dims, -9999)
    }),
    list(ncdf4::ncvar_def("depth_bnds", "meters", list(nv, dims$depth)))
  )
  nc <- ncdf4::nc_create(path, defs)
  for (i in seq_along(vars)) {
    ncdf4::ncvar_put(nc, vars[i], 1:8 + 100 * (i - 1))
  }
  ncdf4::ncvar_put(nc, "depth_bnds", c(0, 50, 50, 150))
  ncdf4::ncatt_put(nc, "depth", "bounds", "depth_bnds")
  ncdf4::ncatt_put(nc, "depth", "positive", "down")
  ncdf4::ncatt_put(nc, "time", "climatology", "climatology_bounds")
  ncdf4::nc_close(nc)
  path
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


test_that("woa_load_nc reads a field as a voxel with WOA's depth bounds", {
  f <- make_woa_file(tempfile(fileext = ".nc"), vars = c("t_an", "t_sd"))
  on.exit(unlink(f), add = TRUE)

  v <- woa_load_nc(f)
  expect_s4_class(v, "SpatVoxel")
  expect_identical(names(v), c("t_an_depth=0", "t_an_depth=100"))
  expect_identical(depths(v), c(0, 100))
  expect_equal(depth_bounds(v)$lower, c(0, 50))
  expect_equal(depth_bounds(v)$upper, c(50, 150))
  expect_true(terra::same.crs(v, "EPSG:4326"))
  # The variable's name is kept in terra::varnames(), not read from labels.
  expect_identical(unique(terra::varnames(v)), "t_an")
  # A WOA period is a climatology without a date.
  expect_false(terra::timeInfo(v)$time)
  # The top-left cell is lat 1.5, lon 0.5: stored values 3 (0 m) and 7.
  expect_equal(unname(terra::values(v)[1, ]), c(3, 7))

  # Another field from the same file.
  sd <- woa_load_nc(f, field = "sd")
  expect_identical(names(sd), c("t_sd_depth=0", "t_sd_depth=100"))
  expect_equal(unname(terra::values(sd)[1, ]), c(103, 107))
})

test_that("woa_load_nc rejects files and fields it can't read", {
  tmp <- tempfile("woa_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # Variable names that aren't a WOA variable code and field code.
  for (name in c("do_sea", "o_inc")) {
    f <- make_woa_file(file.path(tmp, paste0(name, ".nc")), vars = name)
    expect_error(woa_load_nc(f, field = "sea"),
                 "No World Ocean Atlas variables in")
  }

  # "sea" is a valid WOA field code, but this file has no sea layers.
  f <- make_woa_file(file.path(tmp, "m01.nc"), vars = "t_an")
  expect_error(
    woa_load_nc(f, field = "sea"),
    "No layers for field 'sea' in m01.nc. Fields present in this file: an.",
    fixed = TRUE
  )
  expect_error(woa_load_nc(f, field = "xx"), "field must be one of")
})
