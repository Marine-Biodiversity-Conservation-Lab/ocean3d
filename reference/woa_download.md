# Download a WOA NetCDF file (with caching)

Download World Ocean Atlas 2023 NetCDF files from the NCEI THREDDS
server. Files are cached in
[`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md)
(or a user-supplied `output_dir`) and skipped on subsequent calls unless
`force = TRUE`.

## Usage

``` r
woa_download(
  variable,
  period = "annual",
  resolution = "0.25",
  decade = NULL,
  output_dir = NULL,
  force = FALSE,
  quiet = FALSE
)
```

## Arguments

- variable:

  Character. One of `"temperature"`, `"salinity"`, `"dissolved_oxygen"`,
  `"oxygen_saturation"`, `"AOU"`, `"nitrate"`, `"phosphate"`,
  `"silicate"`, `"density"`.

- period:

  Character or numeric. `"annual"` (default), `"monthly"` (all 12
  months), `"seasonal"` (4 seasons), or a numeric vector where 0 =
  annual, 1:12 = monthly, 13:16 = seasonal.

- resolution:

  Character or numeric. `"0.25"` (default), `"1"`, or `"5"` degrees.

- decade:

  Character. Decade code (e.g. `"decav"`, `"all"`). Defaults to the
  canonical decade for the variable.

- output_dir:

  Character. Destination directory. Defaults to
  [`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md).

- force:

  Logical. If `TRUE`, re-download even if the file exists. Default
  `FALSE`.

- quiet:

  Logical. Suppress download progress. Default `FALSE`.

## Value

Character vector of paths to downloaded .nc files.

## Details

On first use of the persistent cache in an interactive session, the
function prompts for consent to write to user filespace (per CRAN
Repository Policy). When the `curl` package is available, the estimated
size of each remote file is reported before downloading. A full set of
WOA files can be many gigabytes — supply `output_dir` to direct the
cache elsewhere (HPC scratch, external drive, etc.).

## Reading the files

Describe a downloaded file with
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md),
naming the field's variable, then read it with
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md).
Crop the proxy to your area first with
[`terra::crop()`](https://rspatial.github.io/terra/reference/crop.html),
so that only those cells are read from the file. Each file holds one
variable per statistical field, named `{code}_{field}`, such as `t_an`
for the objectively analyzed temperature. Leave out `var` and
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
lists the file's variables.

The variable codes are `t` (temperature), `s` (salinity), `o` (dissolved
oxygen), `O` (oxygen saturation), `A` (AOU), `n` (nitrate), `p`
(phosphate), `i` (silicate) and `I` (density). The field codes are `an`
(objectively analyzed climatology), `mn` (statistical mean), `dd`
(number of observations), `sd` (standard deviation), `se` (standard
error), `oa` (statistical mean minus the objectively analyzed mean),
`ma` (the period's objectively analyzed mean minus the annual one), `gp`
(number of grid squares with a statistical mean within the radius of
influence), `sdo` (objectively analyzed standard deviation) and `sea`
(standard error of the analysis). The article "Working with World Ocean
Atlas data", on the package website, shows each step. See the WOA 2023
Product Documentation:
<https://repository.library.noaa.gov/view/noaa/70581>.

The voxel carries WOA's own depth bounds
([`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md)),
its CRS (EPSG:4326) and no date, since a WOA file is a climatology for
one period. For the same reason, the 12 monthly (or 4 seasonal) files
can't be joined as one series with
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md):
open each file on its own. Because the voxel has depth bounds, the
`bounds` argument of
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) and the
other 3D functions doesn't apply to it. WOA's bounds are halfway between
its standard levels, which is what their default, `bounds = "midpoint"`,
gives.

## Examples

``` r
if (FALSE) { # \dontrun{
# Annual temperature climatology at 1-degree resolution, cached in a
# temporary directory so nothing is written to user filespace
tmp <- file.path(tempdir(), "woa")
f <- woa_download("temperature", period = "annual", resolution = "1",
                  output_dir = tmp)

# Describe the objectively analyzed mean without reading it, crop it to a
# study area, and read only that area as a SpatVoxel
temp <- as_proxy(f, var = "t_an")
v <- as_voxel(terra::crop(temp, terra::ext(-20, -10, 30, 40)))
depth_bounds(v)

# All 12 monthly files, downloaded in parallel, into the default cache
monthly <- woa_download("temperature", period = "monthly",
                        resolution = "1")
length(monthly)

# Specific periods: 0 = annual, 1:12 = monthly, 13:16 = seasonal
winter <- woa_download("dissolved_oxygen", period = 13, resolution = "1")

# Other variables use their own canonical decade codes. Note that
# 0.25-degree files are large (~1 GB each for salinity).
sal <- woa_download("salinity", period = "annual", resolution = "0.25")

# Re-download a file that is already cached
f <- woa_download("temperature", period = "annual", resolution = "1",
                  output_dir = tmp, force = TRUE, quiet = TRUE)
} # }
```
