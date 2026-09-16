# Extract nearest, surface and bottom values from a 3D netCDF variable

For each observation, this function extracts the full vertical profile
at the nearest longitude, latitude and time cell. It then returns three
summary columns: nearest valid depth, surface and bottom.

## Usage

``` r
extract3d_all(data, nc, var = NULL, ...)
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

The input observations with three extracted columns: `nearest_*`,
`surface_*` and `seabottom_*`. Data frames, tibbles and `sf` objects
retain their structure; matrices and named lists are returned as data
frames.

## Examples

``` r
nc_file <- system.file(
  "extdata", "example_3d.nc",
  package = "ocean3d"
)

observations <- data.frame(
  lon = c(0, 1),
  lat = c(40, 41),
  depth = c(60, 100),
  date = as.Date(c("2020-01-01", "2020-01-03"))
)

extract3d_all(
  data = observations,
  nc = nc_file,
  var = "temp"
)
#>   lon lat depth       date nearest_temp surface_temp seabottom_temp
#> 1   0  40    60 2020-01-01         15.0         20.0           10.0
#> 2   1  41   100 2020-01-03         13.5         23.5           13.5
```
