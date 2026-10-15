#' Load GEBCO bathymetry raster
#'
#' Load a GEBCO bathymetry raster from a NetCDF file using [terra::rast()].
#' Validates that the file is NetCDF format, contains an `elevation` variable,
#' and has global extent (-180 to 180, -90 to 90). Values are returned as-is
#' (GEBCO uses negative values for below sea level).
#'
#' @param file_path Character. Path to GEBCO bathymetry NetCDF file (e.g.,
#'   `"gebco_2025_sub_ice_topo/GEBCO_2025_sub_ice.nc"`).
#'
#' @returns SpatRaster with elevation values in metres.
#'
#' @examples
#' # A real GEBCO grid is a multi-gigabyte download, so this example builds a
#' # tiny stand-in with the same structure (global extent, `elevation`
#' # variable, negative values below sea level) and loads it.
#' nc_path <- file.path(tempdir(), "gebco_example.nc")
#'
#' template <- terra::rast(
#'   nrows = 36, ncols = 72,
#'   xmin = -180, xmax = 180, ymin = -90, ymax = 90,
#'   crs = "EPSG:4326"
#' )
#' terra::values(template) <- seq(-6000, 0, length.out = terra::ncell(template))
#' terra::varnames(template) <- "elevation"
#' terra::writeCDF(template, nc_path, varname = "elevation", overwrite = TRUE)
#'
#' bathy <- load_gebco_bathymetry(nc_path)
#' bathy
#'
#' # Depths are negative below sea level; flip the sign for use with
#' # vect_to_envelope(), which expects positive depths.
#' depth <- -bathy
#' terra::global(depth, "max", na.rm = TRUE)
#'
#' # In practice, point at a downloaded GEBCO NetCDF instead:
#' # bathy <- load_gebco_bathymetry(
#' #   "gebco_2025_sub_ice_topo/GEBCO_2025_sub_ice.nc"
#' # )
#'
#' unlink(nc_path)
#' @export
load_gebco_bathymetry <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Check file extension is .nc
  if (!grepl("\\.nc$", file_path, ignore.case = TRUE)) {
    stop("Expected a NetCDF (.nc) file, got: ", basename(file_path))
  }

  bathy <- terra::rast(file_path)

  # Check variable name is "elevation"
  if (!("elevation" %in% terra::varnames(bathy))) {
    stop(
      "Expected variable 'elevation' in NetCDF, found: ",
      toString(terra::varnames(bathy))
    )
  }

  # Check global extent (-180 to 180, -90 to 90)
  e <- as.vector(terra::ext(bathy))
  if (e[1] != -180 || e[2] != 180 || e[3] != -90 || e[4] != 90) {
    stop(
      "Expected global extent (-180, 180, -90, 90), got: (",
      toString(e), ")"
    )
  }

  bathy
}
