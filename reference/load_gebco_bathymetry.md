# Load GEBCO bathymetry raster

Load a GEBCO bathymetry raster from a NetCDF file using
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).
Validates that the file is NetCDF format, contains an `elevation`
variable, and has global extent (-180 to 180, -90 to 90). Values are
returned as-is (GEBCO uses negative values for below sea level).

## Usage

``` r
load_gebco_bathymetry(file_path)
```

## Arguments

- file_path:

  Character. Path to GEBCO bathymetry NetCDF file (e.g.,
  `"gebco_2025_sub_ice_topo/GEBCO_2025_sub_ice.nc"`).

## Value

SpatRaster with elevation values in metres.

## Examples

``` r
# A real GEBCO grid is a multi-gigabyte download, so this example builds a
# tiny stand-in with the same structure (global extent, `elevation`
# variable, negative values below sea level) and loads it.
nc_path <- file.path(tempdir(), "gebco_example.nc")

template <- terra::rast(
  nrows = 36, ncols = 72,
  xmin = -180, xmax = 180, ymin = -90, ymax = 90,
  crs = "EPSG:4326"
)
terra::values(template) <- seq(-6000, 0, length.out = terra::ncell(template))
terra::varnames(template) <- "elevation"
terra::writeCDF(template, nc_path, varname = "elevation", overwrite = TRUE)

bathy <- load_gebco_bathymetry(nc_path)
bathy
#> class       : SpatRaster
#> size        : 36, 72, 1  (nrow, ncol, nlyr)
#> dimensions  : longitude, latitude (72, 36}
#> resolution  : 5, 5  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> source      : gebco_example.nc
#> varname     : elevation
#> name        : elevation

# Depths are negative below sea level; flip the sign for use with
# vect_to_envelope(), which expects positive depths.
depth <- -bathy
terra::global(depth, "max", na.rm = TRUE)
#>            max
#> elevation 6000

# In practice, point at a downloaded GEBCO NetCDF instead:
# bathy <- load_gebco_bathymetry(
#   "gebco_2025_sub_ice_topo/GEBCO_2025_sub_ice.nc"
# )

unlink(nc_path)
```
