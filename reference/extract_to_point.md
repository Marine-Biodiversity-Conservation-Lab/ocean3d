# Extract values from netCDF files at point observations

A convenience wrapper over the build and intersect steps: it turns the
observations into
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
points with
[`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md),
describes the netCDF data with
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md),
and reads the value at every point with
[extract()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract-SpatProxy.md),
which reads each time step once for all the points in it. The same can
be done step by step; see the workflow article,
[`vignette("workflow-values-in-shapes")`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/workflow-values-in-shapes.md).

## Usage

``` r
extract_to_point(
  data = NULL,
  nc,
  lon = NULL,
  lat = NULL,
  depth = NULL,
  date = NULL,
  var = NULL,
  method = c("2d", "surface", "bottom", "nearest"),
  lon_col = "lon",
  lat_col = "lat",
  date_col = "date",
  depth_col = "depth",
  id_col = NULL,
  lon_dim = NULL,
  lat_dim = NULL,
  depth_dim = NULL,
  time_dim = NULL,
  time_match = c("nearest", "exact"),
  max_time_diff = NULL,
  output_prefix = NULL,
  output_col = NULL
)
```

## Arguments

- data:

  Optional observations: a data frame, tibble, `sf` object with POINT
  geometries, numeric matrix with column names, or named list. If
  `NULL`, give the coordinates directly with `lon`, `lat`, and `depth`
  and `date` as needed. `sf` points need a CRS and are transformed to
  longitude and latitude.

- nc:

  The netCDF data: a file path, several paths for a series of files on
  one grid (each with its own time steps, such as one file per day), a
  connection opened with
  [`ncdf4::nc_open()`](https://rdrr.io/pkg/ncdf4/man/nc_open.html), or a
  `SpatProxy` from
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md).
  Paths and connections are read through
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md);
  a `SpatProxy` is used as it is, so `var` and the dimension arguments
  don't apply to it.

- lon, lat, depth, date:

  Direct input when `data = NULL`: numeric longitudes, latitudes and
  depths (metres, positive down), and dates or date-times. Each has
  length 1 or the same length as the others. `depth` is needed for
  `method = "nearest"`.

- var:

  Name of the variable. `NULL` detects it when the file holds one data
  variable.

- method:

  `"2d"` (default) for a variable without depth or with a single depth
  level, or `"surface"`, `"bottom"` or `"nearest"` for one with depth
  levels (see Details).

- lon_col, lat_col, date_col, depth_col:

  Columns of `data` with the longitude, latitude, date and depth.
  `date_col` is needed when the variable has dated time steps (not a
  climatology); `depth_col` for `method = "nearest"`. The coordinate
  columns must be numeric.

- id_col:

  Optional column of observation identifiers. It is checked to exist;
  extraction doesn't use it.

- lon_dim, lat_dim, depth_dim, time_dim:

  Optional dimension names, passed to
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md).

- time_match:

  `"nearest"` (default) or `"exact"`: how observation times are matched
  to time steps, as in
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md).

- max_time_diff:

  Optional largest gap, in days, between an observation and its nearest
  time step with `time_match = "nearest"`; `NULL` (default) allows any
  gap.

- output_col, output_prefix:

  Name of the output column. If both are given, `output_col` is used. By
  default, the variable's name, or with `var` detected, `surface_<var>`,
  `seabottom_<var>` or `nearest_<var>` for those methods.

## Value

With `data`, the observations with the value column added. Data frames,
tibbles and `sf` objects keep their class; matrices and named lists
become data frames. With direct input, a numeric vector, one value per
point. A value is `NA` when the observation lies off the grid, has no
coordinates, has no time step within `max_time_diff`, or has no level
its `method` can take.

## Details

You choose the depth level with `method`, which follows the rules of the
same name in `extract()`: `"nearest"` takes the level nearest the
observation's depth, skipping levels without a value; `"surface"` takes
the shallowest level, even if it has no value; `"bottom"` takes the
deepest level with a value. `"2d"` is for a variable without depth. Each
observation takes the cell that contains it. Times are matched to the
nearest time step by default (see `time_match` and `max_time_diff`): a
date by its UTC day, a date-time as an instant.

## See also

[extract()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract-SpatProxy.md)
for the step this wraps;
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
and
[`extract_to_track()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_track.md).

## Examples

``` r
nc_file <- system.file("extdata", "example_3d.nc", package = "ocean3d")

# One observation: the nearest level with a value, on 3 January.
extract_to_point(nc = nc_file, lon = 0, lat = 40, depth = 60,
                 date = "2020-01-03", var = "temp", method = "nearest")
#> [1] 15.5

# A table of observations keeps its columns and gains one.
obs <- data.frame(id = c("a", "b"), lon = c(0, 1), lat = c(40, 41),
                  depth = c(60, 100),
                  date = as.Date(c("2020-01-03", "2020-01-03")))
extract_to_point(obs, nc_file, var = "temp", method = "bottom")
#>   id lon lat depth       date temp
#> 1  a   0  40    60 2020-01-03 10.5
#> 2  b   1  41   100 2020-01-03 13.5
```
