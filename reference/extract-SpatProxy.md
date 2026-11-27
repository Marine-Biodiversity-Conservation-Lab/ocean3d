# Values at many points from a SpatProxy

`extract()` reads the value of a gridded variable at each point of a
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md),
for many points at once. It gives the same values as calling
`mask(as_voxel(crop(proxy, p), time = p), p)` for each point `p`, much
faster: the points are matched to time steps together, and each time
step is read once, cropped to its points. The article
[`vignette("workflow-values-in-shapes")`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/workflow-values-in-shapes.md)
shows it in "Many points at once".

## Usage

``` r
# S4 method for class 'SpatProxy,SpatVector3D'
extract(
  x,
  y,
  method = c("nearest", "nearest_level", "bottom", "surface"),
  time_match = c("nearest", "exact"),
  max_time_diff = NULL,
  bind = FALSE,
  ...
)
```

## Arguments

- x:

  A
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md),
  from
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md):
  a netCDF file, a series of files, or a voxel in memory.

- y:

  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  points, from
  [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md).

- method:

  For a variable with depth, the depth level each point takes:
  `"nearest"` (default), `"nearest_level"`, `"bottom"` or `"surface"`,
  as in
  [mask()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md).
  `method` sets the depth level; `time_match` sets the time step.

- time_match, max_time_diff:

  How a point's time is matched to a time step, as in
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md):
  `"nearest"` (default), the step nearest in time, or `"exact"`; and for
  nearest matching, the largest gap in days. `NULL` (default) sets no
  limit. These set the time step; `method` sets the depth level.

- bind:

  `FALSE` (default) for a data frame of the values; `TRUE` for the
  points themselves with the values added, as terra's
  `extract(raster, vector, bind = TRUE)` returns the vector.

- ...:

  Not used.

## Value

With `bind = FALSE`, a data frame with one row per point, in the order
of `y`:

- `feature_id`, the point's ID;

- `value`, `NA` for a point off the grid, without a matching time step,
  or with no level its `method` can take;

- `time_step`, the date-time of the time step used, `NA` for a source
  without dated steps or a point with no matching step. A point off the
  grid keeps its step's time;

- `depth_level`, the depth of the level used in metres, `NA` for a
  variable without depth or when `value` is `NA`.

With `bind = TRUE`, `y` itself, a
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
with the same features in the same order and all their columns, plus
`value`, `time_step` and `depth_level` as above. They aren't named
`time` and `depth` in either form, because a point already has its own
`time` and `depth`. Points with a column named `value`, `time_step` or
`depth_level` stop with an error; rename the column first.

## Details

Each point takes the time step nearest its time, or the step equal to it
with `time_match = "exact"`, as
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
matches a point. With dates, a point's date matches by UTC day. A
climatology with one time step, such as a World Ocean Atlas annual
field, gives that step to every point, whatever its time. On a source
with dated time steps, even one, a point without a time gets `NA`, with
a warning.

Each point takes the cell that holds it and, for a variable with depth,
one depth level, chosen with `method` as
[mask()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md)
does for points:

- `"nearest"` (default): the level nearest the point's depth, skipping
  levels without a value, as `extract_to_point(method = "nearest")`
  does.

- `"nearest_level"`: the level nearest the point's depth, even if it has
  no value.

- `"bottom"`: the deepest level with a value.

- `"surface"`: the shallowest level.

For a variable without depth, each point takes the value in its cell,
and `method` doesn't apply.

A time step's points are read together while they are close enough. A
group of points is read in smaller groups of nearby points when one read
would exceed the read budget (the option `ocean3d.max_read_mb`, 1024 MB
by default), or would hold more than 10,000 cells and depth levels for
each point, so that points far apart aren't read with the water between
them. The values are the same.

## Experimental

`SpatProxy` and
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
are experimental. Their slots and behaviour may change without notice
while the core object is built.

## See also

[mask()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md)
and
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
for one shape at a time;
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md).

## Examples

``` r
p <- as_proxy(system.file("extdata", "example_3d.nc", package = "ocean3d"))
pts <- as_vector3d(data.frame(
  lon = c(0.1, 0.9), lat = c(40.2, 40.8), depth = c(5, 60),
  date = as.Date(c("2020-01-01", "2020-01-05"))
), time = "date")
terra::extract(p, pts)
#>   feature_id value  time_step depth_level
#> 1          1    20 2020-01-01           0
#> 2          2    19 2020-01-06          50

# The deepest level with a value at each point instead.
terra::extract(p, pts, method = "bottom")
#>   feature_id value  time_step depth_level
#> 1          1    10 2020-01-01         100
#> 2          2    14 2020-01-06         100

# The points themselves, with the values added to their columns.
as.data.frame(terra::extract(p, pts, bind = TRUE))
#>   feature_id vertex lon  lat depth       time value  time_step depth_level
#> 1          1      1 0.1 40.2     5 2020-01-01    20 2020-01-01           0
#> 2          2      1 0.9 40.8    60 2020-01-05    19 2020-01-06          50
```
