# Extract values from a 2D netCDF variable

Use this function when the variable has longitude and latitude
dimensions, and optionally a time dimension, but no depth dimension.

## Usage

``` r
extract2d(data, nc, var = NULL, ...)
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
  "extdata", "example_2d.nc",
  package = "ocean3d"
)

observations <- data.frame(
  lon = c(0, 1),
  lat = c(40, 41),
  date = as.Date(c("2020-01-01", "2020-01-03"))
)

extract2d(
  data = observations,
  nc = nc_file,
  var = "temp"
)
#>   lon lat       date temp
#> 1   0  40 2020-01-01 15.0
#> 2   1  41 2020-01-03 18.5
```
