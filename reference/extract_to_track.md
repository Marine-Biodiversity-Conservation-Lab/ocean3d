# Summarise environmental values along tracks

Samples points along each track, such as a dive or a vessel route, reads
the value at each sample, and summarises the samples of each track. It
is a convenience wrapper over the three steps of an analysis, and the
[workflow
article](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/workflow-values-in-shapes.md)
shows the same steps written out:

## Usage

``` r
extract_to_track(
  data,
  nc,
  var = NULL,
  track_col = NULL,
  sampling = NULL,
  spacing_km = NULL,
  depth_step_m = NULL,
  time_step = NULL,
  time_unit = c("hours", "minutes", "seconds", "days", "months", "years"),
  max_samples = 100000L,
  horizontal_radius_km = 0,
  depth_radius_m = 0,
  method = NULL,
  statistic = "mean",
  weighting = c("distance", "duration", "equal"),
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

  Vertices of the tracks, in the order they were travelled: a data
  frame, matrix, named list or sf POINT object. sf points are
  transformed to longitude and latitude. To use sf lines, convert them
  first, e.g. with
  `as.data.frame(sf::st_coordinates(sf::st_transform(lines, 4326)))` and
  `lon_col = "X", lat_col = "Y", track_col = "L1"`.

- nc:

  The data: a netCDF file path, several paths for a series of files on
  one grid, an open
  [`ncdf4::nc_open()`](https://rdrr.io/pkg/ncdf4/man/nc_open.html)
  connection, or a
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)
  from
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md).

- var:

  Name of the variable, for a netCDF source. `NULL` detects the only
  data variable of the file.

- track_col:

  Column with track IDs. `NULL` (default) reads all of `data` as one
  track, with ID `"1"`. Each track needs at least two vertices.

- sampling, spacing_km, depth_step_m, time_step, time_unit, max_samples:

  Where the samples go along each track, as in
  [`sample_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/sample_3d.md):
  only the vertices (default without a step), or with a spacing in
  kilometres or a depth or time step, points in between as well.

- horizontal_radius_km:

  Radius of a circle around each sample, in kilometres, below 10000 as
  in
  [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md).
  `0` (default) reads the value at the sample itself. A buffer needs
  cells with a width, so data with a single longitude or latitude and no
  cell bounds (a station or a transect) takes no radius.

- depth_radius_m:

  Distance in metres above and below each sample's depth, for the
  buffers. Needed with a radius for a variable with depth, and used only
  with a radius.

- method:

  For a variable with depth and no radius, the depth level each sample
  takes, as in
  [extract()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract-SpatProxy.md):
  `"nearest"` (default; the level nearest the sample's depth that has a
  value), `"nearest_level"` (the nearest level, even if it has no
  value), `"bottom"` (the deepest level with a value) or `"surface"`.
  Leave `NULL` with a radius; for a variable without depth, leave `NULL`
  or use `"2d"`.

- statistic:

  One or more of `"mean"` (default), `"sd"`, `"min"` and `"max"`, of the
  samples of each track.

- weighting:

  How each sample is weighted in the track summaries: `"distance"`
  (default) by `support_km`, `"duration"` by `support_seconds` (needs
  dates), or `"equal"`.

- lon_col, lat_col, depth_col, date_col:

  Vertex columns. The depth is needed for a variable with depth, unless
  `method` is `"bottom"` or `"surface"`, and for `depth_step_m`. The
  date is needed for data with dated time steps, for `time_step` and for
  duration weights. Dates and date-times are read as in
  [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md).

- lon_dim, lat_dim, depth_dim, time_dim:

  Dimension names for a netCDF source, as in
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md).

- time_match, max_time_diff:

  How a sample's time is matched to a time step, as in
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md):
  the nearest step (default) or an exact one, and the largest gap in
  days for a nearest match.

- output_col:

  Name of the value column in `samples`. `NULL` uses the variable's
  name.

## Value

A list of two data frames:

- `samples`: one row per sample, in order along each track, with
  `track_id`, `sample_id`, the position (in `lon_col` and `lat_col`,
  with longitudes from -180 to 180), the depth and date when used (in
  `depth_col` and `date_col`; dates stay dates), `distance_km` from the
  track's start, `support_km`, `support_seconds`, `is_interpolated`
  (`FALSE` for the track's own vertices), and the value in `output_col`.
  With a radius, also `<output>_n`, the number of values in the sample's
  buffer.

- `summary`: one row per track, with `track_id`, a column per statistic,
  and `n`, the number of samples with a value.

## Details

1.  **Build:** the vertices of each track become a line with
    [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md),
    and
    [`sample_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/sample_3d.md)
    places points along it, with depth and time interpolated between
    vertices. With `horizontal_radius_km`,
    [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)
    turns each sample into a buffer.

2.  **Intersect:** without a radius,
    [extract()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract-SpatProxy.md)
    reads the value at each sample, in the cell that holds it and at the
    depth level `method` picks. With a radius, the data for the sample's
    time step are read with
    [`crop()`](https://rspatial.github.io/terra/reference/crop.html) and
    [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
    (or
    [rast()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/rast-SpatProxy.md)
    for a variable without depth), and
    [`mask()`](https://rspatial.github.io/terra/reference/mask.html)
    keeps the values inside its buffer: the cells the circle touches,
    and the levels whose slab of water overlaps the depth window.
    [`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md)
    gives their mean, weighted by cell volume.

3.  **Summarise:** the statistics of each track's sample values, each
    sample weighted by `weighting`.

Each sample is matched to the time step nearest its time, or equal to it
with `time_match = "exact"`, as
[extract()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract-SpatProxy.md)
matches points: a date by its UTC day, a date-time as an instant. A
sample gets `NA` when it is outside the grid, on a cell or level without
a value, has no date on data with dated time steps, or has no time step
that matches. One warning counts these samples, by reason. The track
summaries leave them out, and `n` is the number of samples with a value.

With `weighting = "distance"` each sample stands for half the distance
to the samples either side (`support_km`), so a stretch sampled more
densely doesn't count for more; with `"duration"` for half the time
(`support_seconds`). The weighted standard deviation is the one
[`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md)
uses. A track whose samples all have zero weight, such as a track that
doesn't move with distance weights, gets equal weights. Neighbouring
samples are usually correlated, so the standard deviation describes
their spread, not the uncertainty of the mean.

## See also

[`sample_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/sample_3d.md),
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md),
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

# A dive from the surface to 50 m over two days.
dive <- data.frame(lon = c(2.5, 6.5), lat = c(3.5, 5.5), depth = c(0, 50),
                   date = as.Date(c("2020-01-01", "2020-01-03")))

# The value at a sample at least every 50 km, and the mean along the dive.
track <- extract_to_track(dive, path, spacing_km = 50)
track$summary
#>   track_id mean  n
#> 1        1   42 11

# The mean within 30 km and 10 m above and below each sample instead.
extract_to_track(dive, path, spacing_km = 50, horizontal_radius_km = 30,
                 depth_radius_m = 10)$summary
#>   track_id    mean  n
#> 1        1 44.2459 11
```
