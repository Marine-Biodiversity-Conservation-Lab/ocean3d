# 3D voxel model: one layer per standard depth level

A `SpatVoxel` is a
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
in which **depth is the layer index** and the cell values are the
variable (temperature, oxygen, presence, ...). Layer names follow the
`{variable}_depth={value}` convention and must be distinct and ordered
shallow to deep. The depth axis is grid-wide: every cell is sampled at
the same set of standard depths.

## Details

A voxel may have interior gaps — a cell can be NA at one depth and
non-NA at the depths above and below it.

## See also

[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
[SpatVolume](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVolume-class.md)

## Examples

``` r
# Two cells sampled at three standard depths. Build one with as_voxel(),
# which names the layers for you and sorts them shallow to deep.
r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1,
                 nlyrs = 3)
terra::values(r) <- cbind(c(12, 11), c(9, NA), c(6, 5))
v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")

names(v)          # depth lives in the layer name
#> [1] "temp_depth=0"   "temp_depth=100" "temp_depth=200"
terra::values(v)  # ...and the variable in the cell values
#>      temp_depth=0 temp_depth=100 temp_depth=200
#> [1,]           12              9              6
#> [2,]           11             NA              5

# Note cell 2: NA at 100 m, but with values above and below it. A voxel can
# hold that interior gap; a SpatEnvelope cannot.

# It is an ordinary SpatRaster underneath, so terra operations work on it
# unchanged. Re-wrapping with as_voxel() is cheap either way, since the
# constructor is idempotent.
methods::is(v, "SpatRaster")
#> [1] TRUE
terra::nlyr(terra::crop(v, terra::ext(0, 1, 0, 1)))
#> [1] 3
identical(names(as_voxel(v)), names(v))
#> [1] TRUE

# Validity is enforced on construction: every layer name must carry its
# depth, and the layers must run shallow to deep.
bad <- r
names(bad) <- c("a", "b", "c")
try(methods::new("SpatVoxel", bad))
#> Error in validObject(.Object) : 
#>   invalid class “SpatVoxel” object: every layer name must follow {variable}_depth={value}
```
