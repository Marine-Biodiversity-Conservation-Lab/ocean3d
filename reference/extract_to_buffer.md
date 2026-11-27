# Summarise environmental values around observations

Calculates one statistic per observation from the values in a buffer
around it: a circle of `horizontal_radius_km`, a depth window of
`depth_radius_m` above and below the observation, and a time window from
`time_before` to `time_after`. It is a convenience wrapper over the
three steps of an analysis, and the [workflow
article](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/workflow-values-in-shapes.md)
shows the same steps written out:

## Usage

``` r
extract_to_buffer(
  data,
  nc,
  var = NULL,
  horizontal_radius_km = 0,
  depth_radius_m = 0,
  time_before = 0,
  time_after = 0,
  time_unit = c("days", "hours", "minutes", "seconds"),
  statistic = c("mean", "sd", "min", "max"),
  weighting = c("weighted", "equal"),
  method = NULL,
  lon_col = "lon",
  lat_col = "lat",
  depth_col = "depth",
  date_col = "date",
  lon_dim = NULL,
  lat_dim = NULL,
  depth_dim = NULL,
  time_dim = NULL,
  time_match = c("nearest", "exact"),
  max_time_diff = NULL,
  output_col = NULL
)
```

## Arguments

- data:

  Observation table, matrix, named list or sf POINT object.

- nc:

  The data: a netCDF file path, several paths for a series of files on
  one grid, an open
  [`ncdf4::nc_open()`](https://rdrr.io/pkg/ncdf4/man/nc_open.html)
  connection, a
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)
  from
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md),
  or a
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  in memory.

- var:

  Name of the variable, for a netCDF source. `NULL` detects the only
  data variable of the file.

- horizontal_radius_km:

  Radius of the circle around each observation, in kilometres, measured
  on the WGS84 ellipsoid, below 10000 as in
  [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md).
  `0` (default) uses the observation's own cell.

- depth_radius_m:

  Distance in metres above and below the observation's depth. `0`
  (default) picks one level per cell with `method`.

- time_before, time_after:

  Time before and after the observation, in `time_unit`. Both `0`
  (default) read the one time step `time_match` matches. For
  observations with dates and no hours, the window covers whole UTC
  days, including the last.

- time_unit:

  `"days"` (default), `"hours"`, `"minutes"` or `"seconds"`.

- statistic:

  One of `"mean"` (default), `"sd"`, `"min"` and `"max"`.

- weighting:

  `"weighted"` (default) weights each value by the volume of its cell
  and level with a depth window, and by its cell's area without one (one
  level per cell) or for a variable without depth; `"equal"` counts each
  value once.

- method:

  How the level in each cell is chosen when there is no depth window:
  one of the rules of
  [`mask()`](https://rspatial.github.io/terra/reference/mask.html) for
  points, `"nearest"` (default; the level nearest the observation's
  depth that has a value), `"nearest_level"` (the nearest level, even if
  it has no value), `"bottom"` or `"surface"`. Leave `NULL` with a depth
  window. For a variable without depth, leave `NULL` or use `"2d"`.

- lon_col, lat_col, depth_col, date_col:

  Observation columns. The depth is needed with a depth window or a
  `"nearest"` rule, and the date for data with several time steps. Dates
  and date-times are read as in
  [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md).

- lon_dim, lat_dim, depth_dim, time_dim:

  Dimension names for a netCDF source, as in
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md).

- time_match, max_time_diff:

  Without a time window, how the time step is matched, as in
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md):
  the nearest step (default) or an exact one, and the largest gap in
  days for a nearest match.

- output_col:

  Name of the output column. `NULL` uses the variable's name.

## Value

`data` with two columns added: the statistic, named `output_col` or
after the variable, and `<output>_n`, the number of values it is built
on. sf input keeps its geometry. If `data` already has either column,
the function stops with an error; choose another `output_col`.

## Details

1.  **Build:** each observation becomes a point with
    [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md),
    and a buffer with
    [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)
    when `horizontal_radius_km` is positive.

2.  **Intersect:** the data for its time steps are read with
    [`crop()`](https://rspatial.github.io/terra/reference/crop.html) and
    [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
    (or
    [rast()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/rast-SpatProxy.md)
    for a variable without depth), and the values inside the buffer are
    kept with
    [`mask()`](https://rspatial.github.io/terra/reference/mask.html).

3.  **Summarise:**
    [`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md)
    gives the statistic of those values.

A cell or depth level counts in full when the buffer reaches it, as in
[`mask()`](https://rspatial.github.io/terra/reference/mask.html): a cell
the circle touches, and a level whose slab of water overlaps the depth
window. With `weighting = "weighted"` each value is weighted by the
volume of its cell and level (by the cell's area for a variable without
depth, or when `method` picks one level per cell), as
[`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md)
does.

Without a depth window, `method` picks the level in each cell: by
default the level nearest the observation's depth that has a value, as
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
does. Without a radius, the buffer is the observation's own cell.
Without a time window, the time step is the one `time_match` matches. A
time window pools the values of every time step inside it into one
statistic, each value keeping its own weight.

An observation outside the grid, without a date on data with several
time steps, without a depth where one is needed (or with a negative
one), or with no time step in its window or within `max_time_diff` gets
`NA`, and `<output>_n` is 0. One warning counts these observations, by
reason. On data with several time steps on one day, observations need
date-times or a time window: a date alone can't pick one of the steps,
so the function stops.

## See also

[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
[`extract_to_area()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_area.md),
and the [workflow
article](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/workflow-values-in-shapes.md).

## Examples

``` r
# Daily temperature on a 10 x 10 grid of 1 degree cells at four depths.
path <- tempfile(fileext = ".nc")
lon <- ncdf4::ncdim_def("lon", "degrees_east", seq(0.5, 9.5))
lat <- ncdf4::ncdim_def("lat", "degrees_north", seq(0.5, 9.5))
depth <- ncdf4::ncdim_def("depth", "m", c(5, 15, 30, 60))
time <- ncdf4::ncdim_def("time", "days since 2020-01-01", 0:2)
temp <- ncdf4::ncvar_def("temp", "degree_Celsius",
                         list(lon, lat, depth, time), -9999)
nc <- ncdf4::nc_create(path, temp)
ncdf4::ncvar_put(nc, temp, seq_len(10 * 10 * 4 * 3) / 10)
ncdf4::nc_close(nc)

obs <- data.frame(lon = 5.2, lat = 5.3, depth = 20,
                  date = as.Date("2020-01-02"))

# The mean within 150 km and 10 m above and below, on the observation's day.
extract_to_buffer(obs, path, horizontal_radius_km = 150,
                  depth_radius_m = 10)
#>   lon lat depth       date     temp temp_n
#> 1 5.2 5.3    20 2020-01-02 61.60395     26

# Pooled over a day either side.
extract_to_buffer(obs, path, horizontal_radius_km = 150,
                  depth_radius_m = 10, time_before = 1, time_after = 1)
#>   lon lat depth       date     temp temp_n
#> 1 5.2 5.3    20 2020-01-02 61.60395     78
```
