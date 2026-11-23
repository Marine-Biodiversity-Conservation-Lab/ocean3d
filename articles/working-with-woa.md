# Working with World Ocean Atlas data

The [World Ocean
Atlas](https://www.ncei.noaa.gov/products/world-ocean-atlas) (WOA 2023)
from NOAA’s National Centers for Environmental Information gives
long-term averages (climatologies) of ocean variables such as
temperature, salinity and dissolved oxygen, on a global grid at standard
depth levels. This article shows how to:

- download a WOA file with
  [`woa_download()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_download.md);
- find the variables and statistical fields in it;
- read one field for a study area as a `SpatVoxel`, so that only that
  area is read from the file;
- compare fields, read monthly climatologies and other variables, and
  measure a volume of water.

``` r

library(ocean3d)
```

## Download a file

[`woa_download()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_download.md)
downloads WOA 2023 files from NCEI and returns their paths. Files are
kept in a cache
([`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md),
or `output_dir`), so a second call doesn’t download them again. The
first time a file is written to the cache in an interactive session,
you’re asked for permission.

``` r

f <- woa_download("temperature", period = "annual", resolution = "0.25")
basename(f)
#> [1] "woa23_decav_t00_04.nc"
```

You choose the file with three arguments:

- **`variable`**: `"temperature"`, `"salinity"`, `"dissolved_oxygen"`,
  `"oxygen_saturation"`, `"AOU"`, `"nitrate"`, `"phosphate"`,
  `"silicate"` or `"density"`.
- **`period`**: `"annual"` (the default), `"monthly"` (12 files),
  `"seasonal"` (4 files), or numbers: 0 for annual, 1 to 12 for the
  months, and 13 to 16 for the seasons January–March, April–June,
  July–September and October–December.
- **`resolution`**: `"0.25"` (the default), `"1"` or `"5"` degrees. A
  0.25° file is large (about 1 GB for salinity), so
  [`woa_download()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_download.md)
  reports each file’s size before downloading it.

Together, these make up the file name,
`woa23_{decade}_{code}{period}_{resolution}.nc`: here `decav` (the
average of all decades, the default for temperature), `t` (temperature),
`00` (annual) and `04` (0.25°).

## What a file holds

Each file holds one variable, such as temperature, as several
**statistical fields**. Each field is a netCDF variable named
`{code}_{field}`, such as `t_an`. Without `var`,
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
lists them:

``` r

as_proxy(f)
#> Error:
#> ! Automatic variable detection requires exactly one data variable in each
#> netCDF file (coordinate and bounds variables are ignored).
#> File: woa23_decav_t00_04.nc
#> Data variables found: t_an, t_mn, t_dd, t_sd, t_se, t_oa, t_gp, t_sdo, t_sea.
#> Please specify `var` manually.
```

The variable codes are:

| Variable                          | Code | Default decade |
|-----------------------------------|------|----------------|
| Temperature                       | `t`  | `decav`        |
| Salinity                          | `s`  | `decav`        |
| Dissolved oxygen                  | `o`  | `all`          |
| Oxygen saturation                 | `O`  | `all`          |
| Apparent oxygen utilization (AOU) | `A`  | `all`          |
| Nitrate                           | `n`  | `all`          |
| Phosphate                         | `p`  | `all`          |
| Silicate                          | `i`  | `all`          |
| Density                           | `I`  | `decav`        |

The codes are case-sensitive: `o` and `O`, and `i` and `I`, are
different variables. The fields are:

| Field | What it holds |
|----|----|
| `an` | Objectively analyzed mean: the climatology, with gaps between observations filled by objective analysis |
| `mn` | Statistical mean of the observations, only in grid squares with at least one |
| `dd` | Number of observations in each grid square and depth |
| `sd` | Standard deviation about the statistical mean |
| `se` | Standard error of the statistical mean |
| `oa` | Statistical mean minus the objectively analyzed mean |
| `ma` | The period’s objectively analyzed mean minus the annual one (monthly and seasonal files only) |
| `gp` | Number of grid squares with a statistical mean within the radius of influence |
| `sdo` | Objectively analyzed standard deviation |
| `sea` | Standard error of the objectively analyzed mean |

For most analyses, `an` is the field to use, since it has a value in
every ocean cell. `dd`, `sd` and `se` show how well a cell was sampled
and how much the observations vary. The [WOA 2023 Product
Documentation](https://repository.library.noaa.gov/view/noaa/70581)
describes each field.

## Describe a field without reading it

[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
describes one field without reading its values: the grid, the depth
levels with their bounds, the time step and the coordinate reference
system (CRS).

``` r

temp <- as_proxy(f, var = "t_an")
temp
#> class      : SpatProxy (experimental)
#> source     : woa23_decav_t00_04.nc
#> variable   : t_an (degrees_celsius), chunks 360 x 180 x 26 x 1
#> dimensions : 1440 x 720 x 102 x 1 (lon x lat x depth x time)
#> lon        : -179.875 to 179.875, 1440 cells, edges from bounds
#> lat        : -89.875 to 89.875, 720 cells, edges from bounds
#> depth      : 0 to 5500 m, 102 levels, shallow first, edges from bounds
#> time       : climatology, 1 step (months since 1955-01-01 00:00:00)
#> crs        : EPSG:4326
```

An annual file has 102 depth levels, down to 5500 m. Monthly and
seasonal files have 57, down to 1500 m. A WOA file is a climatology for
one period, so its single time step has no date.

## Read a study area

A 0.25° field covers 1440 x 720 cells at each depth. Crop the proxy to
your study area first, with
[`terra::crop()`](https://rspatial.github.io/terra/reference/crop.html):
this still reads nothing, but it narrows what
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
will read.

``` r

area <- terra::ext(-30, -10, 20, 40)
small <- terra::crop(temp, area)
small
#> class      : SpatProxy (experimental)
#> source     : woa23_decav_t00_04.nc
#> variable   : t_an (degrees_celsius), chunks 360 x 180 x 26 x 1
#> dimensions : 80 x 80 x 102 x 1 (lon x lat x depth x time)
#> lon        : -29.875 to -10.125, 80 cells, edges from bounds
#> lat        : 20.125 to 39.875, 80 cells, edges from bounds
#> depth      : 0 to 5500 m, 102 levels, shallow first, edges from bounds
#> time       : climatology, 1 step (months since 1955-01-01 00:00:00)
#> crs        : EPSG:4326
```

[`crop()`](https://rspatial.github.io/terra/reference/crop.html) takes
an extent, a raster, a `SpatVector` or an sf object, such as a species
range. It keeps the cells that overlap the area’s extent; to keep only
the cells inside a polygon, mask the voxel with it after reading. An
extent in the 0 to 360 longitude convention works on this -180 to 180
grid too.

[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
then reads the cropped cells from the file as a `SpatVoxel`, one layer
per depth:

``` r

v <- as_voxel(small)
dim(v)
#> [1]  80  80 102
head(depths(v))
#> [1]  0  5 10 15 20 25
head(depth_bounds(v), 4)
#>   depth lower upper
#> 1     0   0.0   2.5
#> 2     5   2.5   7.5
#> 3    10   7.5  12.5
#> 4    15  12.5  17.5
```

The voxel keeps WOA’s own depth bounds, the slab of water each level
stands for, and its CRS (EPSG:4326). Because it has depth bounds, the
`bounds` argument of
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) and the
other 3D functions doesn’t apply to it. WOA’s bounds are halfway between
its levels, which is what their default, `bounds = "midpoint"`, gives
anyway (see [Depth levels and the slabs they stand
for](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/depth-levels.md)).

## Compare statistical fields

Each field is read the same way. A short function reads one field of the
same file for the same area:

``` r

read_field <- function(field) {
  as_voxel(terra::crop(as_proxy(f, var = paste0("t_", field)), area))
}
an <- read_field("an")
mn <- read_field("mn")
dd <- read_field("dd")
sd <- read_field("sd")

pick <- which(depths(an) %in% c(0, 100, 500, 1000, 1500))
area_mean <- function(x) {
  round(terra::global(x[[pick]], "mean", na.rm = TRUE)[, 1], 3)
}
data.frame(
  depth = depths(an)[pick],
  analyzed_mean = area_mean(an),
  statistical_mean = area_mean(mn),
  standard_deviation = area_mean(sd),
  cells_with_data = area_mean(dd > 0)
)
#>   depth analyzed_mean statistical_mean standard_deviation cells_with_data
#> 1     0        20.652           20.585              1.954           0.998
#> 2   100        17.807           17.840              0.732           1.000
#> 3   500        11.940           11.992              0.291           0.997
#> 4  1000         8.495            8.548              0.273           0.983
#> 5  1500         5.955            5.974              0.218           0.975
```

In this area almost every ocean cell has observations at these depths,
so the analyzed and statistical means agree closely; in less sampled
regions they can differ more. To keep only cells with observations,
combine fields with terra, which keeps the voxel’s depths:

``` r

sampled <- terra::ifel(dd > 0, an, NA)
class(sampled)
#> [1] "SpatVoxel"
#> attr(,"package")
#> [1] "ocean3d"
```

## Monthly and seasonal climatologies

`period = "monthly"` downloads 12 files, one per month. Each is a
climatology for its month without a date, so they can’t be joined into
one series with
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md):

``` r

months <- woa_download("temperature", period = "monthly", resolution = "0.25")
as_proxy(months, var = "t_an")
#> Error:
#> ! woa23_decav_t01_04.nc has no dated time axis. Every file in a series
#> needs one, so its steps can be placed in time order.
```

Read each month on its own instead, here the mean surface temperature of
the study area:

``` r

surface <- vapply(months, function(p) {
  m <- as_voxel(terra::crop(as_proxy(p, var = "t_an"), area))
  terra::global(m[[1]], "mean", na.rm = TRUE)[, 1]
}, numeric(1))
data.frame(month = month.abb, surface_mean = round(unname(surface), 2))
#>    month surface_mean
#> 1    Jan        18.94
#> 2    Feb        18.36
#> 3    Mar        18.21
#> 4    Apr        18.51
#> 5    May        19.40
#> 6    Jun        20.74
#> 7    Jul        22.08
#> 8    Aug        23.21
#> 9    Sep        23.60
#> 10   Oct        22.95
#> 11   Nov        21.60
#> 12   Dec        20.24
```

Monthly and seasonal files also have the `ma` field, the month’s
analyzed mean minus the annual one. In August, the surface of this area
is on average 2.55 °C warmer than its annual mean:

``` r

august <- as_voxel(terra::crop(as_proxy(months[8], var = "t_ma"), area))
terra::global(august[[1]], "mean", na.rm = TRUE)
#>                 mean
#> t_ma_depth=0 2.55339
```

## Other variables

Every variable is read the same way; only the codes change. The code for
dissolved oxygen is `o`, and its default decade is `all`:

``` r

o <- woa_download("dissolved_oxygen", period = "annual", resolution = "1")
basename(o)
#> [1] "woa23_all_o00_01.nc"
oxygen <- as_proxy(o, var = "o_an")
oxygen
#> class      : SpatProxy (experimental)
#> source     : woa23_all_o00_01.nc
#> variable   : o_an (micromoles_per_kilogram), chunks 180 x 90 x 51 x 1
#> dimensions : 360 x 180 x 102 x 1 (lon x lat x depth x time)
#> lon        : -179.5 to 179.5, 360 cells, edges from bounds
#> lat        : -89.5 to 89.5, 180 cells, edges from bounds
#> depth      : 0 to 5500 m, 102 levels, shallow first, edges from bounds
#> time       : climatology, 1 step (months since 1965-01-01 00:00:00)
#> crs        : EPSG:4326
```

Check each variable’s units on its proxy: dissolved oxygen here is in
micromoles per kilogram.

## Measure a volume of water

A voxel works with the rest of the package. For example, the volume of
water between 10 and 15 °C in the study area, using WOA’s own depth
bounds:

``` r

band <- occupied(an, fun = function(x) x >= 10 & x <= 15)
volume(band)
#> [1] 2004705
```

The volume is in km³. With the voxel’s own bounds,
`volume(band, bounds = "top")` gives the same result.

## Where next

- [Extracting 3D environmental data for a single
  species](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/woa-environmental-extraction-single-species.md)
  uses WOA temperature and oxygen within a species’ 3D range.
- [`woa_download()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_download.md)
  lists the variable and field codes, and
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
  and
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  describe the proxy and how it is read.
- [`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md)
  shows where the files are kept, and
  [`woa_cache_clear()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_clear.md)
  removes them.
