# Total 3D volume of a rasterized domain

Sums the occupied volume of a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
or a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
over the whole grid. Both representations describe a 3D domain over a 2D
grid, but they store depth in dual roles, so each computes the vertical
extent differently:

## Usage

``` r
volume(x, ...)

# S4 method for class 'SpatEnvelope'
volume(x, ...)

# S4 method for class 'SpatVoxel'
volume(x, bounds = c("top", "midpoint"))

# S4 method for class 'SpatRaster'
volume(x, ...)
```

## Arguments

- x:

  A
  [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
  or
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).
  A bare SpatRaster is rejected; build one with
  [`as_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_envelope.md),
  [`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
  or
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  first.

- ...:

  Arguments for the
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  method, which are an error for an envelope.

- bounds:

  SpatVoxel only. What each depth level stands for vertically: `"top"`
  (default) treats the level as the top of its slab, so the slab runs
  down to the next level; `"midpoint"` puts the level in the middle,
  with edges halfway to each neighbour (the World Ocean Atlas
  convention). Under `"top"` the deepest level has no next level and
  contributes no volume. Outer edges are clamped to `range(depths)`
  under both. Same argument, same meaning, as in
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).

## Value

Numeric of length 1. Total volume in km³.

## Details

- [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md):

  Depth is the cell value, and the interval is solid by construction.
  Volume is `sum(cell_area * (depth_max - depth_min))` over present
  cells.

- [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md):

  Depth is the layer index. Each occupied voxel contributes its cell
  area times the thickness of the slab its depth level stands for, so
  interior gaps are excluded rather than filled.

A voxel volume is quantized to the levels the voxel was sampled at, and
a voxel built by
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
counts a level as occupied on *any* overlap with the envelope. So a
range of 0-100 m discretized onto 0, 50, 100, 150, 200 m occupies three
levels spanning 0-150 m and reports more water than the envelope it came
from. Coarse levels overstate; the error shrinks as the levels tighten
around the range. Where an envelope's limits fall exactly on levels
under `bounds = "top"`, the two agree.

## See also

[`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md)
for the volume two domains share, and
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
for the shared domain itself;
[SpatVolume](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVolume-class.md)
for why the two representations dispatch separately.

## Examples

``` r
fp <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2000,
                  ymin = 0, ymax = 2000,
                  crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
terra::values(fp) <- c(1, 1, 1, NA)

# 2.5D: three 1 km² cells, each 0-200 m => 3 * 1 * 0.2 = 0.6 km³.
e <- as_envelope(fp, depth_min = 0, depth_max = 200)
volume(e)
#> [1] 0.6

# 3D: the same domain on three levels. Under "top" the slabs are 100 m,
# 100 m and 0 m, so the discretized volume matches exactly here.
v <- envelope_to_voxel(e, depths = c(0, 100, 200))
volume(v)
#> [1] 0.6

# An interior gap is volume a voxel can drop and an envelope cannot.
gappy <- v
gappy[[2]] <- terra::setValues(terra::rast(v[[2]]), NA)
gappy <- as_voxel(gappy)
volume(gappy)
#> [1] 0.3
volume(voxel_to_envelope(gappy))
#> [1] 0
```
