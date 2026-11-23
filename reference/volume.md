# 3D volume of a domain

Measures the occupied volume of a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
or
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
polygons, as a total, or cell by cell for envelopes and voxels.
Envelopes and voxels describe a 3D domain over a 2D grid, but they store
depth in dual roles, so each computes the vertical extent differently:

## Usage

``` r
volume(x, per_cell = FALSE, ...)

# S4 method for class 'SpatEnvelope'
volume(x, per_cell = FALSE, ...)

# S4 method for class 'SpatVoxel'
volume(x, per_cell = FALSE, bounds = c("midpoint", "top"))

# S4 method for class 'SpatVector3D'
volume(x, per_cell = FALSE, ...)

# S4 method for class 'SpatRaster'
volume(x, per_cell = FALSE, ...)
```

## Arguments

- x:

  A
  [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
  a
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
  or
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons. A bare SpatRaster is rejected; build one with
  [`as_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_envelope.md),
  [`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
  or
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  first. A netCDF file isn't accepted: describe it with
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
  and read one time step with
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md).

- per_cell:

  `FALSE` (default) for the total volume, `TRUE` for the volume of each
  cell. Not for
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons.

- ...:

  Arguments for the
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  method, which are an error for an envelope or polygons.

- bounds:

  SpatVoxel only. What each depth level stands for vertically:
  `"midpoint"` (default) puts the slab edges halfway to each
  neighbouring level (the World Ocean Atlas convention); `"top"` treats
  the level as the top of its slab, so the slab runs down to the next
  level. Under `"top"` the deepest level has no next level and
  contributes no volume. Outer edges are clamped to `range(depths)`
  under both. Same argument, same meaning, as in
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).
  Not used for a voxel that carries its own depth edges (see
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)).

## Value

With `per_cell = FALSE`, a number: the total volume in km³. With
`per_cell = TRUE`, a single-layer `SpatRaster` named `volume` on the
grid of `x`: the volume of each cell in km³, `NA` where `x` is absent
(no depth interval, or no occupied level). Its sum is the total.

## Details

- [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md):

  Depth is the cell value, and the interval is solid by construction.
  Volume is `sum(cell_area * (depth_max - depth_min))` over present
  cells.

- [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md):

  Depth is the layer index. Each occupied voxel contributes its cell
  area times the thickness of the slab its depth level stands for, so
  interior gaps are excluded rather than filled.

- [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons:

  Measured from the outlines, not a grid: each polygon's area on the
  WGS84 ellipsoid times `depth_max - depth_min`. Polygon edges are
  geodesics, so a polygon drawn along a grid's cell edges measures
  slightly differently from its envelope on that grid, whose cells are
  bounded by parallels. The time interval doesn't change the volume. The
  total adds up the polygons, so where they overlap the shared space is
  counted once for each; the overlap between two sets is
  `volume(intersect_3d(x, y))`. `per_cell = TRUE` is an error, since
  polygons have no cells; put them on a grid with
  [`vector3d_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vector3d_to_envelope.md)
  first. Points and lines have no volume, so they are an error;
  [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)
  turns them into polygons.

A voxel volume is quantized to the levels the voxel was sampled at, and
a voxel built by
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
counts a level as occupied on *any* overlap with the envelope. So a
range of 0-100 m discretized onto 0, 50, 100, 150, 200 m occupies three
levels spanning 0-150 m and reports more water than the envelope it came
from. Coarse levels overstate; the error shrinks as the levels tighten
around the range. Where an envelope's limits fall exactly on levels
under `bounds = "top"`, the two agree.

To measure the space two domains share, measure their intersection:
`volume(intersect_3d(a, b))`.

With one envelope and one voxel,
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
places the envelope on the voxel's depth levels, so the overlap is
measured on those levels while `volume()` of the envelope uses its own
depth limits. To compare the two on the same levels, measure the
envelope on them too:
`volume(envelope_to_voxel(a, depths(b), bounds = bounds,`
`depth_bounds = depth_bounds(b)), bounds = bounds)`. Pass the same
`bounds` to
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
and `volume()`; otherwise the overlap and the input volumes are measured
on different slabs and don't add up.

## See also

[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
for the space two domains share, which `volume()` then measures;
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
for whether two domains share any depth in each cell;
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

# 3D: the same domain on three levels. The slabs are 50 m, 100 m and 50 m,
# so the discretized volume matches exactly here.
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
#> [1] 0.6

# Per cell: the fourth cell is outside the domain, so it is NA.
terra::values(volume(e, per_cell = TRUE))
#>      volume
#> [1,]    0.2
#> [2,]    0.2
#> [3,]    0.2
#> [4,]     NA

# Overlap with a second domain at 100-300 m: intersect, then measure.
b <- as_envelope(fp, depth_min = 100, depth_max = 300)
ov <- intersect_3d(e, b)
volume(ov)
#> [1] 0.3
terra::values(volume(ov, per_cell = TRUE))
#>      volume
#> [1,]    0.1
#> [2,]    0.1
#> [3,]    0.1
#> [4,]     NA

# Polygons: a 10 km circle from 10 to 60 m deep, about
# pi * 10^2 * 0.05 = 15.7 km³.
detection <- as_vector3d(data.frame(lon = 150, lat = -20, depth = 35))
volume(buffer_3d(detection, radius_km = 10, depth_window = 25))
#> [1] 15.70776
```
