# ---------------------------------------------------------------------------
# load_gebco_bathymetry: exercise the validation paths. The happy-path test
# writes a minimal global NetCDF with terra::writeCDF, which may not be
# available in all environments.
# ---------------------------------------------------------------------------

test_that("load_gebco_bathymetry errors when file is missing", {
  expect_error(load_gebco_bathymetry("/no/such/file.nc"), "File not found")
})

test_that("load_gebco_bathymetry errors when file is not .nc", {
  tmp <- tempfile(fileext = ".tif")
  file.create(tmp)
  on.exit(unlink(tmp))
  expect_error(load_gebco_bathymetry(tmp), "NetCDF")
})

test_that("load_gebco_bathymetry happy path and extent validation", {
  # Setup test raster
  tmp_good <- tempfile(fileext = ".nc")
  tmp_bad  <- tempfile(fileext = ".nc")
  on.exit(unlink(c(tmp_good, tmp_bad)))

  global <- terra::rast(nrows = 6, ncols = 12,
                        xmin = -180, xmax = 180,
                        ymin = -90, ymax = 90,
                        crs = "EPSG:4326")
  terra::values(global) <- seq_len(terra::ncell(global)) - 1000
  terra::varnames(global) <- "elevation"

  regional <- terra::rast(nrows = 4, ncols = 4,
                          xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                          crs = "EPSG:4326")
  terra::values(regional) <- seq_len(terra::ncell(regional))
  terra::varnames(regional) <- "elevation"

  ok <- tryCatch({
    terra::writeCDF(global, tmp_good, overwrite = TRUE, varname = "elevation")
    terra::writeCDF(regional, tmp_bad, overwrite = TRUE, varname = "elevation")
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "terra::writeCDF unavailable")

  r <- load_gebco_bathymetry(tmp_good)
  expect_s4_class(r, "SpatRaster")

  expect_error(load_gebco_bathymetry(tmp_bad), "global extent")
})

test_that("load_gebco_bathymetry checks variable name", {
  # Setup test raster
  tmp_good <- tempfile(fileext = ".nc")
  tmp_bad  <- tempfile(fileext = ".nc")
  on.exit(unlink(c(tmp_good, tmp_bad)))

  global <- terra::rast(nrows = 6, ncols = 12,
                        xmin = -180, xmax = 180,
                        ymin = -90, ymax = 90,
                        crs = "EPSG:4326")
  terra::values(global) <- seq_len(terra::ncell(global)) - 1000
  terra::varnames(global) <- "elevation"

  ok <- tryCatch({
    terra::writeCDF(global, tmp_good, overwrite = TRUE, varname = "elevation")
    terra::writeCDF(global, tmp_bad, overwrite = TRUE,
                    varname = "not_elevation")
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "terra::writeCDF unavailable")

  r <- load_gebco_bathymetry(tmp_good)
  expect_named(r, "elevation")

  expect_error(
    load_gebco_bathymetry(tmp_bad),
    paste0("Expected variable 'elevation' in NetCDF, ",
           "found: not_elevation")
  )
})