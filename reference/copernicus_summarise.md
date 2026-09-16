# Summarise Copernicus netCDF data across time

Summarise one or more Copernicus netCDF files across their temporal
dimension while preserving all non-temporal dimensions (for example
longitude, latitude, depth, or projected x/y dimensions). Processing is
performed directly on the native netCDF data using `ncdf4`; conversion
to a `SpatRaster` is not required.

## Usage

``` r
copernicus_summarise(
  file_paths,
  fun = c("mean", "min", "max", "sd"),
  start_datetime = NULL,
  end_datetime = NULL,
  output_dir = NULL,
  filename = NULL,
  na.rm = TRUE,
  force = FALSE,
  quiet = FALSE
)
```

## Arguments

- file_paths:

  Character vector. Paths to one or more Copernicus netCDF (`.nc` or
  `.nc4`) files. This can be the output from
  [`copernicus_load()`](https://marine-biodiversity-conservation-lab.github.io/sharkabc3d/reference/copernicus_load.md).

- fun:

  Character vector of temporal summary statistics. Supported values are
  `"mean"`, `"min"`, `"max"`, and `"sd"`.

- start_datetime:

  Optional start datetime used to filter the available netCDF time steps
  before summarising. A `Date`, `POSIXt`, or character value interpreted
  in UTC.

- end_datetime:

  Optional end datetime. Defaults to `start_datetime` when only a start
  is supplied. If both are `NULL`, all available time steps are
  summarised.

- output_dir:

  Optional character. Directory in which summarised netCDF files are
  written. Defaults to the directory containing the first input file.

- filename:

  Optional character. Output filename. When multiple statistics are
  requested, the statistic is appended before the extension (for example
  `summary_mean.nc`, `summary_max.nc`).

- na.rm:

  Logical. If `TRUE` (default), missing values are ignored within each
  grid cell/depth combination. If `FALSE`, any missing temporal value
  produces a missing summary value at that location.

- force:

  Logical. Overwrite existing summary files. Default `FALSE`.

- quiet:

  Logical. Suppress progress messages. Default `FALSE`.

## Value

Named character vector containing paths to the summarised netCDF files,
with names corresponding to `fun`.

## Details

All variables containing the detected time dimension are summarised.
Variables without a time dimension are copied unchanged from the first
input file so that auxiliary coordinates, grid mappings, and other
static metadata remain available in the output.

The calculation is performed one time step at a time. This avoids
loading the complete multidimensional time series into memory and is
suitable for large Copernicus files. When multiple input files are
supplied they must have compatible non-temporal dimensions and
variables.

One netCDF file is written per requested summary statistic. The time
dimension is removed from the summarised variables; spatial and depth
dimensions are retained.

## See also

[`copernicus_load()`](https://marine-biodiversity-conservation-lab.github.io/sharkabc3d/reference/copernicus_load.md)

## Examples

``` r
if (FALSE) { # \dontrun{
summaries <- copernicus_summarise(
  file_paths = c("thetao_2020_01.nc", "thetao_2020_02.nc"),
  fun = c("mean", "min", "max", "sd"),
  output_dir = "summary"
)
} # }
```
