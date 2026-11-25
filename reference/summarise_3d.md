# Summarise a voxel's values across cells and depth levels

`summarise_3d()` computes statistics of every non-`NA` value of a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
across all its cells and depth levels at once. It is the summarise step
of an analysis: the values inside a shape are kept with
[`mask()`](https://rspatial.github.io/terra/reference/mask.html), and
`summarise_3d()` turns them into a mean, a standard deviation, a range
or a median. Compare
[`terra::global()`](https://rspatial.github.io/terra/reference/global.html),
which summarises each layer, here each depth level, on its own.

## Usage

``` r
summarise_3d(
  x,
  statistic = "mean",
  weighting = c("volume", "equal"),
  by = NULL,
  bounds = c("midpoint", "top")
)
```

## Arguments

- x:

  A
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
  typically `mask(v, shape)` for the values inside a shape, or a list of
  voxels, one per time step, as
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  gives for a time window. `x` can also be a one-layer `SpatRaster` of a
  variable without depth, as
  [rast()](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/rast-SpatProxy.md)
  reads from a
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md),
  or a list of them, one per time step. Envelopes have no values, so
  they aren't accepted; a netCDF file isn't either: describe it with
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
  and read one time step with
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  or `rast()`.

- statistic:

  One or more of `"mean"`, `"sd"`, `"min"`, `"max"` and `"median"`. Each
  becomes a column of the result, in this order.

- weighting:

  `"volume"` (default) weights each value by the volume of its cell and
  depth level, and each value of a one-layer `SpatRaster` by its cell's
  area; `"equal"` counts each value once. A row whose values all have
  zero volume, such as a voxel with a single level and no depth edges,
  or the deepest level under `bounds = "top"`, uses equal weights.

- by:

  `NULL` (default) gives one row per voxel, so one per time step for a
  list; `"depth"` one row per depth level of each voxel, with a `depth`
  column. Values from different time steps aren't combined, so a list
  shows the change from step to step. `"time"` means the same as `NULL`,
  and says in a call on a list that the rows are per time step. A
  one-layer `SpatRaster` has a single level, so `"depth"` is an error
  for it.

- bounds:

  As in
  [`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md):
  where the slab edges fall for a voxel that carries no depth edges;
  used with `weighting = "volume"` only. Voxels read with
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  from a
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)
  carry their edges, so `bounds` isn't used for them.

## Value

A data frame with one row per voxel, or per voxel and depth level with
`by = "depth"`, in the order of `x` and then of depth. Columns: `time`
(the voxel's
[`terra::time()`](https://rspatial.github.io/terra/reference/time.html),
UTC, or `NA` when it has none), `depth` (with `by = "depth"` only, in
metres), one column per `statistic`, and `n`, the number of values the
row is built on, so rows built on too few values can be dropped. A row
with no values has `NA` statistics and `n = 0`.

## Details

With `weighting = "volume"` (the default), each value is weighted by the
volume of its cell *and* depth level, from
`volume(x, per_cell = TRUE, per_level = TRUE)`: cells shrink towards the
poles, and depth levels range from a few metres thick to 100 m or more.
In a water column with thin warm levels near the surface and thick cold
ones below, the volume-weighted mean is colder than the equal-weighted
one. `weighting = "equal"` counts each value once, so `sd` and `median`
are then [`stats::sd()`](https://rdrr.io/r/stats/sd.html) and
[`stats::median()`](https://rdrr.io/r/stats/median.html).

With normalised weights \\w_i\\, the weighted standard deviation is the
square root of \\\sum w_i (v_i - \bar v)^2 / \sum w_i (1 - w_i)\\, the
sample standard deviation under equal weights and the same formula as in
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md).
The weighted median is the value at which the cumulative weight, over
the values sorted from smallest to largest, reaches half the total; when
it reaches half exactly, the two values either side are averaged, as
[`stats::median()`](https://rdrr.io/r/stats/median.html) does. `min` and
`max` ignore the weights.

## See also

[`mask()`](https://rspatial.github.io/terra/reference/mask.html) to keep
the values inside a shape;
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
for the weights;
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
to read one time step of a
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md).

## Examples

``` r
# One water column sampled at seven depths: thin warm levels near the
# surface and thick cold ones below.
r <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0,
                 ymax = 1, nlyrs = 7, crs = "EPSG:4326")
terra::values(r) <- t(c(20, 19, 18, 17, 16, 8, 4))
v <- as_voxel(r, depths = c(0, 10, 20, 30, 40, 100, 200), varname = "temp")

# The volume-weighted mean is colder: the deep levels stand for more water.
summarise_3d(v, statistic = c("mean", "sd", "min", "max"))
#>   time mean       sd min max n
#> 1 <NA> 10.2 6.271899   4  20 7
summarise_3d(v, statistic = "mean", weighting = "equal")
#>   time     mean n
#> 1 <NA> 14.57143 7

# One row per depth level.
summarise_3d(v, by = "depth")
#>   time depth mean n
#> 1 <NA>     0   20 1
#> 2 <NA>    10   19 1
#> 3 <NA>    20   18 1
#> 4 <NA>    30   17 1
#> 5 <NA>    40   16 1
#> 6 <NA>   100    8 1
#> 7 <NA>   200    4 1

# A list of voxels, one per time step, gives one row per step.
later <- v + 1
terra::time(v) <- rep(as.POSIXct("2020-01-01", tz = "UTC"), 7)
terra::time(later) <- rep(as.POSIXct("2020-01-02", tz = "UTC"), 7)
summarise_3d(list(v, later), statistic = c("mean", "median"))
#>         time mean median n
#> 1 2020-01-01 10.2      8 7
#> 2 2020-01-02 11.2      9 7
```
