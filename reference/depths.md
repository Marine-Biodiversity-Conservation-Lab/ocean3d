# Depths of a voxel's layers

Read the depth of each layer of a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).
Depth is the layer index in a voxel, and each layer's depth is stored in
terra's own per-layer field,
[`terra::depth()`](https://rspatial.github.io/terra/reference/depth.html).
Layer names are labels for people to read: the depths are never read
from them, so renaming the layers doesn't change them.

## Usage

``` r
depths(x)
```

## Arguments

- x:

  A
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
  or a
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with its depths set.

## Value

Numeric vector, one depth per layer of `x`. An error when the depths
aren't set, for example after a terra operation that drops them (see
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)).

## Details

`x` can also be a plain
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
whose depths are set in
[`terra::depth()`](https://rspatial.github.io/terra/reference/depth.html),
such as a raster terra read from a netCDF file.

Depths are positive metres increasing downward.

## See also

[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md),
which sets the depths.

## Examples

``` r
r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
terra::values(r) <- runif(terra::ncell(r) * 3)
v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")

names(v)
#> [1] "temp_depth=0"   "temp_depth=100" "temp_depth=200"
depths(v)
#> [1]   0 100 200

# A per-layer summary has one row per layer, in layer order, so the voxel's
# depths label it.
per_depth <- terra::global(v, "mean", na.rm = TRUE)
per_depth$depth <- depths(v)
per_depth
#>                     mean depth
#> temp_depth=0   0.2381563     0
#> temp_depth=100 0.4578622   100
#> temp_depth=200 0.4748873   200
```
