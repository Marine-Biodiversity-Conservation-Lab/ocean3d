# ---------------------------------------------------------------------------
# load_gebco_bathymetry: exercise the validation paths. The happy-path test
# writes a minimal global NetCDF with terra::writeCDF.
# ---------------------------------------------------------------------------

# A global elevation raster, and a regional one that doesn't cover the globe.
gebco_rasters <- function() {
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

  list(global = global, regional = regional)
}

test_that("load_gebco_bathymetry reads a global elevation file", {
  tmp_good <- tempfile(fileext = ".nc")
  on.exit(unlink(tmp_good))
  terra::writeCDF(gebco_rasters()$global, tmp_good, overwrite = TRUE,
                  varname = "elevation")

  r <- load_gebco_bathymetry(tmp_good)
  expect_s4_class(r, "SpatRaster")
  expect_named(r, "elevation")
})

test_that("load_gebco_bathymetry rejects files it can't read", {
  rasters <- gebco_rasters()
  not_nc <- tempfile(fileext = ".tif")
  file.create(not_nc)
  regional <- tempfile(fileext = ".nc")
  terra::writeCDF(rasters$regional, regional, overwrite = TRUE,
                  varname = "elevation")
  wrong_var <- tempfile(fileext = ".nc")
  terra::writeCDF(rasters$global, wrong_var, overwrite = TRUE,
                  varname = "not_elevation")
  on.exit(unlink(c(not_nc, regional, wrong_var)))

  cases <- list(
    list(label = "missing file", path = "/no/such/file.nc",
         pattern = "File not found"),
    list(label = "not a .nc file", path = not_nc, pattern = "NetCDF"),
    list(label = "regional extent", path = regional,
         pattern = "global extent"),
    list(label = "wrong variable name", path = wrong_var,
         pattern = paste0("Expected variable 'elevation' in NetCDF, ",
                          "found: not_elevation"))
  )
  for (case in cases) {
    expect_error(load_gebco_bathymetry(case$path), case$pattern,
                 label = case$label)
  }
})
