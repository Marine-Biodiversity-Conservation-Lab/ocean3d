# Depth edges of a voxel's layers

Read the slab of water each layer of a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
stands for, when the voxel carries its own depth edges. Set them with
the `depth_bounds` argument of
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
or
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
for example from a netCDF file's depth bounds. A voxel read from a
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)
carries the proxy's edges, which are halfway between levels when the
file has no depth bounds (see
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)).

## Usage

``` r
depth_bounds(x)
```

## Arguments

- x:

  A
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).

## Value

A data frame with one row per layer of `x`, shallow to deep: `depth`,
and `lower` and `upper`, the top and bottom of the layer's slab, in
metres positive down. `NULL` when `x` doesn't carry edges for every
layer. The data frame can be passed back as `depth_bounds` to
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
or
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).

## Details

[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md),
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
and [`mask()`](https://rspatial.github.io/terra/reference/mask.html) use
these edges before their `bounds` argument. A voxel without edges, or
with edges for only some of its layers (for example after
[`c()`](https://rdrr.io/r/base/c.html) with a voxel that has none),
gives `NULL`, and those functions work the edges out from the depths
with `bounds` instead.

## See also

[`depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depths.md)
for the layer depths;
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
to set the edges.

## Examples

``` r
r <- terra::rast(nrows = 1, ncols = 1, nlyrs = 3, vals = 1)
v <- as_voxel(r, depths = c(0.5, 1.5, 3),
              depth_bounds = cbind(c(0, 1, 2), c(1, 2, 4)))
depth_bounds(v)
#>   depth lower upper
#> 1   0.5     0     1
#> 2   1.5     1     2
#> 3   3.0     2     4

# A layer subset gives the edges of the layers it kept.
depth_bounds(v[[2:3]])
#>   depth lower upper
#> 1   1.5     1     2
#> 2   3.0     2     4

# A voxel without edges gives NULL, so `bounds` applies.
depth_bounds(as_voxel(r, depths = c(0.5, 1.5, 3)))
#> NULL
```
