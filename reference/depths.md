# Depths of a voxel's layers

Read the depth axis of a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).
Depth is the layer index in a voxel, and the depth itself is carried in
the layer name following the `{variable}_depth={value}` convention used
throughout the package; this parses it back out.

## Usage

``` r
depths(x)
```

## Arguments

- x:

  A
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
  a
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  or a character vector of layer names.

## Value

Numeric vector, one depth per layer of `x` (or per element, when `x` is
character). `NA` for any name that does not follow the convention; an
error when no name does.

## Details

`x` is normally a `SpatVoxel`, but any
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
whose layer names follow the convention works, as does a bare character
vector of layer names — useful for the row names of a
[`terra::global()`](https://rspatial.github.io/terra/reference/global.html)
result, which carry the layer names but not the raster.

Depths are positive metres increasing downward.

## See also

[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md),
which builds those layer names.

## Examples

``` r
r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
terra::values(r) <- runif(terra::ncell(r) * 3)
v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")

names(v)
#> [1] "temp_depth=0"   "temp_depth=100" "temp_depth=200"
depths(v)
#> [1]   0 100 200

# A character vector works too, so the depths of a per-layer summary can be
# recovered from its row names.
per_depth <- terra::global(v, "mean", na.rm = TRUE)
depths(rownames(per_depth))
#> [1]   0 100 200

# Names that do not follow the convention come back NA, as long as at least
# one name does.
depths(c("temp_depth=0", "not_a_depth_layer"))
#> [1]  0 NA
```
