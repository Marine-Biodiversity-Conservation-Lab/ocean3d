# Convert Envelope 2.5D -\> Voxel 3D

Expand a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
to the
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
with an input of depth levels. A depth level belongs to a cell when the
slab of water it stands for shares some thickness with that cell's
`[depth_min, depth_max]` interval; depths outside the interval, and
cells that are `NA` in the envelope, are `NA` in every layer.

## Usage

``` r
envelope_to_voxel(
  x,
  depths,
  values = NULL,
  profile = NULL,
  bounds = c("top", "midpoint"),
  varname = "presence"
)
```

## Arguments

- x:

  SpatEnvelope, e.g. from
  [`as_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_envelope.md)
  or
  [`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md).

- depths:

  Array of values that can be coerced into numeric type, represents
  metres depth below sea level. Sorted shallow to deep, and
  deduplicated, before use.

- values:

  Optional numeric or SpatRaster of same CRS, resolution, extent as `x`.
  Defines values to write for voxel cells that are within depth
  intervals. `NULL` (the default) writes `1` at every occupied level,
  giving simple presence/absence.

- profile:

  Optional function distributing `values` across the depth dimension,
  taking `(ind, depths, n_depths)` and returning the weight to apply.
  `NULL` (the default) writes the same value at every depth, as
  [`profile_flat()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_profiles.md)
  does;
  [`profile_equal()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_profiles.md)
  divides the value equally across the depths a cell occupies. Any
  function meeting the contract in
  [voxel_profiles](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_profiles.md)
  works, including your own.

- bounds:

  How to place the edges of the slab each depth level stands for.
  `"top"` (the default) runs the slab from the level down to the next,
  so the level is its shallow edge; `"midpoint"` puts the edges halfway
  between neighbouring levels. Under both, the outer edges are clamped
  to `range(depths)` rather than extrapolated.

- varname:

  Name to use for depth layer, taking on form of
  `{varname}_depth={depths[i]}`

## Value

A
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
with one layer per depth in `depths`, on the grid of `x`, layers ordered
shallow to deep.

## Details

- `values` is the magnitude, which can be a single number for the whole
  grid or a SpatRaster carrying a value per cell. Default value is `1`,
  which makes the output a presence mask. Ex. `values` can be 2D
  fisheries effort rasters.

- `profile` decides how that magnitude is spread down the water column.
  Default `NULL` means the full value is written at every occupied
  level. Providing a `profile_*` function determines how `values` is
  transformed for each depth level. Ex.
  [`profile_equal()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_profiles.md)
  divides it evenly over the levels the cell occupies, so sum of all
  depth levels in the voxel equals `values`.

A depth level stands for a slab of water rather than a knife-edge, and a
cell occupies that level when its envelope overlaps the slab by more
than a shared edge. An interval of `[10, 20]` against levels `c(0, 100)`
therefore occupies the 0 m level, rather than falling between the levels
and coming back empty. Touching is not overlapping: an interval that
ends exactly where a slab begins, `[0, 100]` against the slab `100-200`,
does not occupy that level, so adjacent depth ranges are placed on
disjoint levels and do not intersect as voxels. This is the same rule
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
applies between two envelopes.

Where a slab's edges fall is a property of the dataset the levels came
from, not of the levels themselves, so `bounds` selects the convention:

- `"top"` (the default) runs each slab from its level down to the next,
  so a level names the top of what it stands for. Levels
  `c(0, 100, 200)` give slabs `0-100`, `100-200`.

- `"midpoint"` puts the edges halfway between neighbouring levels, so a
  level sits in the middle of what it stands for. This is the World
  Ocean Atlas convention: `0, 5, 10, ...` the 0 m layer covers `0-2.5`
  m, the 5 m layer `2.5-7.5` m, the 10 m layer `7.5-12.5` m, and so on.

Under either convention the slabs cover the span from the shallowest
level to the deepest and no further: outer edges are clamped to
`range(depths)` rather than extrapolated, which is why WOA's 0 m layer
starts at the surface rather than half a gap above it. The expansion is
therefore limited only at the two ends when a cell whose envelope lies
wholly above `min(depths)` or wholly below `max(depths)`. Under `"top"`
the deepest level has no next level to run to, so it stands for a
zero-thickness slab: an envelope is recorded there only when it reaches
below `max(depths)`, and one that ends exactly at `max(depths)` stops at
the level above. For the same reason a
[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
round trip drops each cell's deepest occupied level: the envelope's
`depth_max` names that level, which is the top of its slab.

## See also

[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md),
the reverse (and lossy) collapse.

## Examples

``` r
fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
terra::values(fp) <- c(1, NA)
e <- as_envelope(fp, depth_min = 50, depth_max = 200)

# presence at every standard depth the envelope covers
terra::values(envelope_to_voxel(e, depths = c(0, 50, 100, 200, 300)))
#>      presence_depth=0 presence_depth=50 presence_depth=100 presence_depth=200
#> [1,]               NA                 1                  1                 NA
#> [2,]               NA                NA                 NA                 NA
#>      presence_depth=300
#> [1,]                 NA
#> [2,]                 NA

# the same pattern, carrying a constant of your choosing instead of 1
terra::values(
  envelope_to_voxel(e, depths = c(0, 50, 100, 200, 300),
                    values = 5)
)
#>      presence_depth=0 presence_depth=50 presence_depth=100 presence_depth=200
#> [1,]               NA                 5                  5                 NA
#> [2,]               NA                NA                 NA                 NA
#>      presence_depth=300
#> [1,]                 NA
#> [2,]                 NA

# spread that value down the column instead: an even share at each of the
# depths the cell occupies, summing back to 1
terra::values(
  envelope_to_voxel(e, depths = c(0, 50, 100, 200, 300),
                    profile = profile_equal,
                    varname = "time")
)
#>      time_depth=0 time_depth=50 time_depth=100 time_depth=200 time_depth=300
#> [1,]           NA           0.5            0.5             NA             NA
#> [2,]           NA            NA             NA             NA             NA

# the same envelope read against WOA-style layer bounds: [60, 70] sits in the
# 50-150 m slab of the 100 m level, where the default puts it at 0-100 m
shallow <- as_envelope(fp, depth_min = 60, depth_max = 70)
terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200)))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]                1                 NA                 NA
#> [2,]               NA                 NA                 NA
terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200),
                                bounds = "midpoint"))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]               NA                  1                 NA
#> [2,]               NA                 NA                 NA

# a per-cell magnitude comes in as a raster: split each cell's effort evenly
# over the depths it occupies, conserving the cell total
effort <- terra::setValues(terra::rast(fp), c(100, NA))
terra::values(
  envelope_to_voxel(
    e, depths = c(0, 50, 100, 200, 300), values = effort,
    profile = profile_equal,
    varname = "effort"
  )
)
#>      effort_depth=0 effort_depth=50 effort_depth=100 effort_depth=200
#> [1,]             NA              50               50               NA
#> [2,]             NA              NA               NA               NA
#>      effort_depth=300
#> [1,]               NA
#> [2,]               NA
```
