# 2.5D min-max envelope: exactly depth_min, depth_max

A `SpatEnvelope` is a
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with exactly two layers, `depth_min` and `depth_max`, in which **depth
is the cell value**. The variable is the object itself: the envelope
delimits the vertical extent of a spatial phenomenon (typically a
species range) rather than sampling a variable within it.

## Details

The vertical interval is per-cell and continuous, and is solid by
construction — an envelope cannot represent a gap in the vertical
distribution. Depths are positive metres increasing downward, and
`depth_max` is at least `depth_min` in every cell.

## See also

[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
[SpatVolume](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVolume-class.md)

## Examples

``` r
# A two-cell footprint, one cell present and one absent. Build an envelope
# with as_envelope() (from a raster footprint) or vect_to_envelope() (from
# polygons).
fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
terra::values(fp) <- c(1, NA)
e <- as_envelope(fp, depth_min = 50, depth_max = 200)

names(e)          # always exactly these two layers
#> [1] "depth_min" "depth_max"
terra::values(e)  # ...and depth is the cell value, in metres
#>      depth_min depth_max
#> [1,]        50       200
#> [2,]        NA        NA

# The interval is per-cell, so limits can vary across the grid.
terra::values(
  as_envelope(fp, depth_min = 50,
              depth_max = terra::setValues(terra::rast(fp), c(120, 300)))
)
#>      depth_min depth_max
#> [1,]        50       120
#> [2,]        NA        NA

# Any other set of layers is rejected: a bare footprint is not an envelope.
try(methods::new("SpatEnvelope", fp))
#> Error in validObject(.Object) : 
#>   invalid class “SpatEnvelope” object: layers must be exactly: depth_min, depth_max
```
