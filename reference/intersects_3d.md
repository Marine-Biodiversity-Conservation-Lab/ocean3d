# Do two objects share any 3D space?

`intersects_3d()` tests, cell by cell, whether `x` and `y` overlap in
3D. It is the yes/no form of
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
and takes the same inputs. It computes only the presence pattern, with
no cell areas and no volumes, so it is the cheap choice for richness and
tally maps.

## Usage

``` r
intersects_3d(x, y, ...)

# S4 method for class 'SpatEnvelope,SpatEnvelope'
intersects_3d(x, y, ...)

# S4 method for class 'SpatVoxel,SpatVoxel'
intersects_3d(x, y)

# S4 method for class 'SpatEnvelope,SpatVoxel'
intersects_3d(x, y, bounds = c("top", "midpoint"), ...)

# S4 method for class 'SpatVoxel,SpatEnvelope'
intersects_3d(x, y, bounds = c("top", "midpoint"), ...)

# S4 method for class 'SpatEnvelope,ANY'
intersects_3d(x, y, ...)

# S4 method for class 'SpatVoxel,ANY'
intersects_3d(x, y)

# S4 method for class 'ANY,SpatEnvelope'
intersects_3d(x, y, ...)

# S4 method for class 'ANY,SpatVoxel'
intersects_3d(x, y, ...)

# S4 method for class 'ANY,ANY'
intersects_3d(x, y, ...)
```

## Arguments

- x, y:

  The two objects. At least one must be a
  [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
  or a
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).
  The other can be either of those, a single-layer `SpatRaster`
  footprint, or polygons (`SpatVector`, `sf`, `sfc`). Rasters must share
  one grid (CRS, extent, resolution).

- ...:

  Arguments for the voxel methods. They are an error when both inputs
  are envelopes.

- bounds:

  Only when one input is an envelope and the other a voxel. How the
  envelope is placed on the voxel's depth levels: `"top"` (default) or
  `"midpoint"`. See
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).

## Value

A single-layer boolean `SpatRaster` named `intersects`, on the grid of
the 3D input: `TRUE`, `FALSE` or `NA` per cell as above.

## Details

Each cell gets one of three answers:

- `TRUE`: both objects are present and their depths overlap.

- `FALSE`: both are present but their depths do not overlap, or only one
  of them is present.

- `NA`: neither is present. There is nothing to compare.

This is how `terra` answers the 2D question for two rasters, so the
result sums and plots like any other boolean layer. Depth intervals that
only touch do not overlap. When `y` is a 2D footprint or polygons there
is no depth axis to compare, so the test is only whether both are
present.

The order of `x` and `y` does not matter.

## See also

[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
for the shared space itself;
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) to keep
a voxel's values inside a domain.

## Examples

``` r
fp <- terra::rast(nrows = 1, ncols = 4, xmin = 0, xmax = 4000,
                  ymin = 0, ymax = 1000,
                  crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")

# Four cells. A is present in cells 1-3 at 0-100 m. B is present in cells
# 1, 2 and 4, at 50-200 m in cell 1 and 300-400 m elsewhere.
a <- as_envelope(terra::setValues(fp, c(1, 1, 1, NA)),
                 depth_min = 0, depth_max = 100)
b <- as_envelope(terra::setValues(fp, c(1, 1, NA, 1)),
                 depth_min = terra::setValues(fp, c(50, 300, NA, 300)),
                 depth_max = terra::setValues(fp, c(200, 400, NA, 400)))

# Cell 1: overlap. Cell 2: both present, depths disjoint. Cell 3: A only.
# Cell 4: B only.
terra::values(intersects_3d(a, b))
#>      intersects
#> [1,]       TRUE
#> [2,]      FALSE
#> [3,]      FALSE
#> [4,]      FALSE

# Summing a stack of these gives a richness map. na.rm = TRUE counts each
# TRUE as 1 and each FALSE as 0.
richness <- sum(c(intersects_3d(a, b), intersects_3d(a, a)), na.rm = TRUE)
terra::values(richness)
#>      sum
#> [1,]   2
#> [2,]   1
#> [3,]   1
#> [4,]   0
```
