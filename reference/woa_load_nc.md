# Load a WOA NetCDF file

Load a WOA .nc file and select layers for a given statistical field.
Returns a SpatRaster with layer names already following the
`{variable}_depth={value}` convention used throughout this package
(native to WOA NetCDFs).

## Usage

``` r
woa_load_nc(file_path, field = "an")
```

## Arguments

- file_path:

  Character. Path to a WOA .nc file.

- field:

  Character. Statistical field to select. One of: `"an"` (objectively
  analyzed climatology), `"mn"` (statistical mean), `"dd"` (number of
  observations), `"sd"` (standard deviation), `"se"` (standard error),
  `"oa"` (mean minus climatology), `"gp"` (number of mean values within
  radius of influence), `"sdo"` (objectively analyzed standard
  deviation), `"sea"` (standard error of the analysis). See the WOA 2023
  Product Documentation:
  <https://repository.library.noaa.gov/view/noaa/70581>

## Value

SpatRaster with standardized depth layer names. Wrap it with
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
to use it with the `SpatVoxel` functions
([`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
[`mask()`](https://rspatial.github.io/terra/reference/mask.html),
[`extract_to_area()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_area.md),
...). The point, buffer and track extraction functions take the netCDF
file itself, not this raster.

## Examples

``` r
# Build a small WOA-style NetCDF to demonstrate field selection.
# Real files come from woa_download().
nc <- file.path(tempdir(), "woa_demo.nc")
demo <- terra::rast(nrows = 2, ncols = 2, nlyrs = 4,
                    xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                    crs = "EPSG:4326")
names(demo) <- c("t_an_depth=0", "t_an_depth=100",
                 "t_sd_depth=0", "t_sd_depth=100")
terra::values(demo) <- matrix(1:16, nrow = 4)
terra::writeCDF(demo, nc, overwrite = TRUE, split = TRUE)

# Objectively analyzed climatology (the default field)
temp <- woa_load_nc(nc)
names(temp)
#> [1] "t_an_depth=0"   "t_an_depth=100"

# Standard deviation layers from the same file
names(woa_load_nc(nc, field = "sd"))
#> [1] "t_sd_depth=0"   "t_sd_depth=100"

unlink(nc)

if (FALSE) { # \dontrun{
# Typical use on a downloaded file
f <- woa_download("temperature", period = "annual", resolution = "1")
temp <- woa_load_nc(f, field = "an")
terra::plot(temp[["t_an_depth=0"]])
} # }
```
