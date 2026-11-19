# Depth levels and the slabs they stand for

A `SpatVoxel` samples the water column at a set of depth levels, such as
the World Ocean Atlas standard depths. A `SpatEnvelope` stores a
continuous `[depth_min, depth_max]` interval per cell. Converting
between the two with
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
or combining them with
[`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html),
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
or
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
needs a rule for which levels an interval occupies. This article sets
out that rule.

The examples use a two-cell footprint, so every result can be printed in
full.

``` r

library(ocean3d)

fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
terra::values(fp) <- c(1, 1)
```

## A depth level stands for a slab, not a knife edge

A cell occupies a level when its envelope overlaps the slab that level
stands for by more than a shared edge. This matters for thin intervals.
An envelope of `[10, 20]` contains none of `c(0, 100, 200)`, but it is
clearly within the water the 0 m level represents:

``` r

thin <- as_envelope(fp, depth_min = 10, depth_max = 20)
terra::values(envelope_to_voxel(thin, depths = c(0, 100, 200)))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]                1                 NA                 NA
#> [2,]                1                 NA                 NA
```

## Touching is not overlapping

By default the slab edges fall halfway between neighbouring levels, so
with levels at 0, 50, 150 and 200 m one edge falls at 100 m. An envelope
of `[0, 100]` ends exactly where the 150 m slab begins, so it occupies
the 0 m and 50 m levels and not the 150 m one:

``` r

levels <- c(0, 50, 150, 200)
upper <- as_envelope(fp, depth_min = 0, depth_max = 100)
lower <- as_envelope(fp, depth_min = 100, depth_max = 200)

terra::values(envelope_to_voxel(upper, depths = levels))
#>      presence_depth=0 presence_depth=50 presence_depth=150 presence_depth=200
#> [1,]                1                 1                 NA                 NA
#> [2,]                1                 1                 NA                 NA
```

Two envelopes that meet at 100 m therefore do not intersect.
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
gives the same answer for the envelopes themselves and for their voxels,
which land on disjoint levels:

``` r

terra::values(intersects_3d(upper, lower))
#>      intersects
#> [1,]      FALSE
#> [2,]      FALSE

terra::values(intersects_3d(
  envelope_to_voxel(upper, depths = levels),
  envelope_to_voxel(lower, depths = levels)
))
#>      intersects
#> [1,]      FALSE
#> [2,]      FALSE
```

## Where the slab edges fall: `bounds`

Where the slab edges fall is a property of the dataset the levels came
from, so `bounds` selects the convention:

- `"midpoint"` (the default) puts the edges halfway between neighbouring
  levels. This matches World Ocean Atlas levels. It is close for model
  output such as Copernicus Marine, which stores values at the middle of
  each layer, but it leaves out the water above the shallowest level and
  below the deepest.
- `"top"` runs each slab from its level down to the next, so each level
  is the top of its slab.

An interval of `[60, 70]` lands in a different level under each:

``` r

shallow <- as_envelope(fp, depth_min = 60, depth_max = 70)

terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200)))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]               NA                  1                 NA
#> [2,]               NA                  1                 NA
terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200),
                                bounds = "top"))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]                1                 NA                 NA
#> [2,]                1                 NA                 NA
```

The same `bounds` argument is accepted by
[`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html),
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md),
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
and
[`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md)
wherever they place an envelope on a voxel’s depth levels. See
[`?envelope_to_voxel`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
for the full rule, including how the shallowest and deepest levels are
clamped.

## A dataset’s own edges: `depth_bounds`

Both conventions work the edges out from the levels. When the dataset
gives its own edges, for example a netCDF file’s depth bounds, a
`SpatVoxel` can carry them: pass them to
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
or
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
as `depth_bounds`, a two-column table of top and bottom with one row per
level.
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
and the functions listed above then use these edges and ignore `bounds`,
and
[`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md)
reads them back:

``` r

edges <- cbind(c(0, 65, 150), c(65, 150, 250))
v <- envelope_to_voxel(shallow, depths = c(0, 100, 200), depth_bounds = edges)
terra::values(v)
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]                1                  1                 NA
#> [2,]                1                  1                 NA
depth_bounds(v)
#>   depth lower upper
#> 1     0     0    65
#> 2   100    65   150
#> 3   200   150   250
```

Here the 0 m level’s slab runs to 65 m, so `[60, 70]` occupies both the
0 m and the 100 m levels, where `"midpoint"` gave only the 100 m level.
