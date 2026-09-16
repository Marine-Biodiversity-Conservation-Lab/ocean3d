# The 3D space two objects share

`intersect_3d()` returns the part of space that is inside both `x` and
`y`. It answers the question "what do these two share?". To ask only "do
they share anything?", use
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md).
To keep one object's values where the other is present, use
[`mask()`](https://rspatial.github.io/terra/reference/mask.html).

## Usage

``` r
intersect_3d(x, y, ...)

# S4 method for class 'SpatEnvelope,SpatEnvelope'
intersect_3d(x, y, ...)

# S4 method for class 'SpatVoxel,SpatVoxel'
intersect_3d(x, y)

# S4 method for class 'SpatEnvelope,SpatVoxel'
intersect_3d(x, y, bounds = c("top", "midpoint"), ...)

# S4 method for class 'SpatVoxel,SpatEnvelope'
intersect_3d(x, y, bounds = c("top", "midpoint"), ...)

# S4 method for class 'SpatEnvelope,ANY'
intersect_3d(x, y, ...)

# S4 method for class 'SpatVoxel,ANY'
intersect_3d(x, y)

# S4 method for class 'ANY,SpatEnvelope'
intersect_3d(x, y, ...)

# S4 method for class 'ANY,SpatVoxel'
intersect_3d(x, y, ...)

# S4 method for class 'ANY,ANY'
intersect_3d(x, y, ...)
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

A
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
when both inputs are envelopes, or when one is an envelope and the other
is 2D. Otherwise a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
of presence. Either is on the grid of the 3D input.

## Details

At least one of `x` and `y` must be a 3D object: a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
or a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).
The other can be a 3D object too, or a 2D object. A 2D object restricts
the result horizontally and leaves its depths alone.

- Two envelopes:

  Cell by cell, the result runs from the deeper `depth_min` to the
  shallower `depth_max`. Intervals that only touch share no water, so
  that cell is empty (`NA`).

- Two voxels:

  Both must be sampled at the same depth levels. A level is in the
  result where both voxels occupy it.

- An envelope and a voxel:

  The envelope is first placed on the voxel's depth levels with
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).
  The result is a voxel.

- A 3D object and a `SpatRaster`:

  The raster is a footprint. Its non-`NA` cells are inside; its values
  are not read. It must have one layer and be on the same grid.

- A 3D object and polygons:

  A `SpatVector`, `sf` or `sfc` object. The polygons are rasterised onto
  the 3D object's grid. A cell is inside when its centre is. Polygons in
  another CRS are projected first.

The order of `x` and `y` does not matter.

The result is a domain, not a field. An envelope result carries the
shared depth interval. A voxel result carries presence: `1` where a
level is shared, `NA` elsewhere, in layers named
`presence_depth=<value>`. The cell values of a voxel input are not
carried over. To keep them, use
[`mask()`](https://rspatial.github.io/terra/reference/mask.html).

## terra's [`intersect()`](https://rdrr.io/r/base/sets.html)

[`terra::intersect()`](https://rspatial.github.io/terra/reference/intersect.html)
on two rasters is a 2D test that ignores depth. On a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
or a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
in either position, it is an error that points to
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
and `intersect_3d()`. Plain rasters, vectors and extents keep terra's
behaviour.

## See also

[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
for the yes/no form;
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) to keep
a voxel's values inside a domain;
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
for the volume of the result;
[`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md),
which is built on this function.

## Examples

``` r
fp <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3000,
                  ymin = 0, ymax = 1000,
                  crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
terra::values(fp) <- c(1, 1, 1)

# Two species: one at 0-100 m, one at 50-200 m everywhere.
a <- as_envelope(fp, depth_min = 0, depth_max = 100)
b <- as_envelope(fp, depth_min = 50, depth_max = 200)

# They share 50-100 m in every cell.
shared <- intersect_3d(a, b)
terra::values(shared)
#>      depth_min depth_max
#> [1,]        50       100
#> [2,]        50       100
#> [3,]        50       100
volume(shared)
#> [1] 0.15

# With a voxel, the result is a presence voxel on the voxel's levels.
bv <- envelope_to_voxel(b, depths = c(0, 50, 100, 150, 200))
names(intersect_3d(a, bv))
#> [1] "presence_depth=0"   "presence_depth=50"  "presence_depth=100"
#> [4] "presence_depth=150" "presence_depth=200"

# With polygons, the envelope is kept where the polygon is.
poly <- terra::vect("POLYGON ((0 0, 2000 0, 2000 1000, 0 1000, 0 0))",
                    crs = terra::crs(fp))
terra::values(intersect_3d(a, poly))
#>      depth_min depth_max
#> [1,]         0       100
#> [2,]         0       100
#> [3,]        NA        NA
```
