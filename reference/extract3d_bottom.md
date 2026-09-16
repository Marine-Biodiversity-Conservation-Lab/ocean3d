# Extract the bottom available layer from a 3D netCDF variable

For each observation, this function extracts the full vertical profile
at the nearest longitude, latitude and time cell, then returns the
deepest non-missing value. This is useful when the bathymetry of the
netCDF grid means that deeper layers are missing over shallow areas.

## Usage

``` r
extract3d_bottom(data, nc, var = NULL, ...)
```

## Arguments

- data:

  Object containing observation points. Can be a data frame, tibble,
  `sf` object with POINT geometries, matrix with column names, or named
  list.

- nc:

  netCDF source. Can be a file path, vector of file paths, list of file
  paths, opened `ncdf4` object, list of opened objects, or a data frame
  with a file column.

- var:

  Name of the variable to extract from the netCDF file. If `NULL`, the
  function tries to detect the variable automatically from each file.
  This only works when each netCDF file contains one single variable.

- ...:

  Additional arguments passed to
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md).

## Value

The input observations with extracted values appended. Data frames,
tibbles and `sf` objects retain their structure; matrices and named
lists are returned as data frames.

## Examples

``` r
nc_file <- system.file(
  "extdata", "example_3d.nc",
  package = "ocean3d"
)

observations <- data.frame(
  lon = c(0, 1),
  lat = c(40, 41),
  date = as.Date(c("2020-01-01", "2020-01-03"))
)

extract3d_bottom(
  data = observations,
  nc = nc_file,
  var = "temp"
)
#>   lon lat       date temp
#> 1   0  40 2020-01-01 10.0
#> 2   1  41 2020-01-03 13.5
```
