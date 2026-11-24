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

# S4 method for class 'SpatEnvelope,SpatVoxel'
intersects_3d(x, y, bounds = c("midpoint", "top"), ...)

# S4 method for class 'SpatVoxel,SpatEnvelope'
intersects_3d(x, y, bounds = c("midpoint", "top"), ...)

# S4 method for class 'ANY,ANY'
intersects_3d(x, y, ...)

# S4 method for class 'SpatVector3D,SpatVector3D'
intersects_3d(x, y, ...)
```

## Arguments

- x, y:

  The two objects: two of
  [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  and
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons or points, in any pairing except two sets of points, or an
  envelope or voxel with a 2D object (a single-layer `SpatRaster`
  footprint or 2D polygons: `SpatVector`, `sf`, `sfc`). Rasters must
  share one grid (CRS, extent, resolution). netCDF files aren't
  accepted; load them with terra and wrap them with
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  first.

- ...:

  Arguments for the voxel methods. They are an error when both inputs
  are envelopes, or when one is points and the other an envelope or a
  voxel.

- bounds:

  Only when one input is a voxel and the other an envelope or
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons. How the envelope or polygons are placed on the voxel's depth
  levels: `"midpoint"` (default) or `"top"`. See
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).
  Not used when the voxel carries its own depth edges (see
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)).
  To measure the result, pass the same `bounds` to
  [`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md).

## Value

A single-layer boolean `SpatRaster` named `intersects`, on the grid of
the 3D input: `TRUE`, `FALSE` or `NA` per cell as above. For two
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
sets, a list with one element per feature of `x`, named by its feature
ID: the feature IDs of the features of `y` it shares space with, in the
order of `y`, or `character(0)` for none.

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
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
polygons and points have depths, so they are compared as
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
compares them, after being put on the 3D object's grid. A cell is `TRUE`
where the voxel has a value at the level nearest a point's depth, or
where a point's depth lies within the envelope's interval.

Two
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
sets (polygons with polygons, or points with polygons) have no grid to
answer on, so the answer is a list instead, as
[`sf::st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
gives in 2D: for each feature of `x`, the features of `y` it shares 3D
space with. A pair is listed exactly when
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
gives it a polygon or point, so outlines that only share an edge, depth
ranges that only touch and times that don't overlap aren't listed.
Swapping `x` and `y` gives the same pairs listed the other way round.

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
