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
intersect_3d(x, y, bounds = c("midpoint", "top"), ...)

# S4 method for class 'SpatVoxel,SpatEnvelope'
intersect_3d(x, y, bounds = c("midpoint", "top"), ...)

# S4 method for class 'SpatEnvelope,ANY'
intersect_3d(x, y, ...)

# S4 method for class 'SpatVoxel,ANY'
intersect_3d(x, y)

# S4 method for class 'ANY,ANY'
intersect_3d(x, y, ...)

# S4 method for class 'SpatVector3D,SpatEnvelope'
intersect_3d(x, y, ...)

# S4 method for class 'SpatVector3D,SpatVoxel'
intersect_3d(x, y, ...)

# S4 method for class 'SpatEnvelope,SpatVector3D'
intersect_3d(x, y, ...)

# S4 method for class 'SpatVoxel,SpatVector3D'
intersect_3d(x, y, ...)

# S4 method for class 'SpatVector3D,SpatVector3D'
intersect_3d(x, y, ...)
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

[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
polygons when both inputs are polygons, and
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
points when one is points and the other polygons. A
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
when both inputs are envelopes, or when one is an envelope and the other
is 2D or a
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md).
Otherwise a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
of presence. Either is on the grid of the 3D input.

## Details

At least one of `x` and `y` must be a 3D object: a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
or
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
polygons or points. The other can be a 3D object too, or, with an
envelope or voxel, a 2D object. A 2D object restricts the result
horizontally and leaves its depths alone.

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

- An envelope or voxel and a `SpatRaster`:

  The raster is a footprint. Its non-`NA` cells are inside; its values
  are not read. It must have one layer and be on the same grid.

- An envelope or voxel and 2D polygons:

  A `SpatVector`, `sf` or `sfc` object. The polygons are rasterised onto
  the 3D object's grid. A cell is inside when a polygon touches it.
  Polygons in another CRS are projected first.

- An envelope and
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons:

  The polygons are first put on the envelope's grid with
  [`vector3d_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vector3d_to_envelope.md),
  then compared as an envelope. Each cell a polygon touches gets its
  depth range; where polygons overlap, a cell gets the shallowest top
  and the deepest bottom, so a depth gap between two polygons over one
  cell is filled.

- A voxel and
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons:

  The polygons are first put on the voxel's grid and depth levels with
  [`vector3d_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vector3d_to_voxel.md),
  one level at a time, then compared as a voxel. A depth gap between two
  polygons over one cell stays empty. When the voxel has no depth edges,
  you set each level's slab with `bounds`, as in
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).

- An envelope or voxel and
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons, either way:

  The polygons' time intervals are dropped, with a message. Polygons
  that miss the grid, or a set with no polygons, give an empty result.

- A voxel and
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  points:

  A point covers one cell and one depth level: the cell that holds it
  and the level whose depth is nearest the point's depth (a tie goes to
  the shallower level, as in
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)).
  The level is chosen by depth alone, so a point whose nearest level has
  no value in the voxel shares nothing there. When the voxel's depth
  edges give uneven slabs, the slab that contains a point can belong to
  another level; the nearest level is used. The result is a presence
  voxel on the voxel's levels, and `bounds` doesn't apply. Points off
  the grid share nothing. Points' times aren't compared, since a voxel
  shows one time step.

- An envelope and
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  points:

  In each cell where a point's depth lies within the envelope's interval
  (both ends included), the result holds an interval from the shallowest
  such point to the deepest. A lone point gives a zero-thickness
  interval, so the result's
  [`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
  is 0: points have no volume.

- [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  points and polygons:

  A point is inside a polygon when
  [`sf::st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
  finds its position in the outline (on WGS84 coordinates sf uses the
  sphere by default, where a point exactly on an edge may not count),
  its depth is within `depth_min` to `depth_max` and its time within the
  polygon's time interval. Depth and time include both ends (with dates,
  the interval covers all of its end day). A point or polygon without a
  time matches any time; one side with dates and the other with hours is
  an error. Each matching pair gives one point, numbered, with columns
  `feature_id_x` and `feature_id_y` naming the pair and the inputs'
  other columns dropped. Two sets of points share no space, so they are
  an error. Lines aren't supported;
  [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)
  turns them into polygons.

- Two sets of
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  polygons:

  Each pair of polygons, one from `x` and one from `y`, whose outlines,
  depth ranges and times overlap gives one polygon: the shared outline
  (as
  [`sf::st_intersection()`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html)),
  the shared depth range, and the shared time interval. Depth ranges
  that only touch share no water. Time intervals include both ends, so
  with dates, two that share a day overlap on it, and with hours, two
  that meet at one instant overlap at that instant (so a polygon at a
  single time can meet another). A polygon without a time overlaps any
  time and passes on the other's interval. One set can't have dates
  while the other has hours. The result's feature IDs are numbers;
  columns `feature_id_x` and `feature_id_y` say which polygons each came
  from, and the inputs' other columns are dropped. Pairs that share
  nothing give no polygon, so the result can have none. To measure the
  overlap, use `volume(intersect_3d(x, y))`.

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
for the volume of the result, in total or per cell. When one input is a
voxel, the envelope is placed on its levels, so compare the overlap with
`volume(envelope_to_voxel(...))` rather than with the envelope's own
volume.

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

# With SpatVector3D polygons, the depth ranges are compared too: a polygon
# from 80 to 300 m over the first two cells shares 80-100 m with `a`.
deep <- as_vector3d(
  sf::st_sf(geometry = sf::st_as_sfc(
    "POLYGON ((0 0, 2000 0, 2000 1000, 0 1000, 0 0))",
    crs = terra::crs(fp)
  )),
  depth_min = 80, depth_max = 300
)
terra::values(intersect_3d(a, deep))
#>      depth_min depth_max
#> [1,]        80       100
#> [2,]        80       100
#> [3,]        NA        NA

# A SpatVector3D point covers the cell that holds it and, on a voxel, the
# level nearest its depth: here the 50 m level of `bv` in the first cell
# (the point is projected onto the grid's CRS).
detection <- as_vector3d(data.frame(lon = 0.0045, lat = 0.0045, depth = 60))
terra::values(intersect_3d(detection, bv))
#>      presence_depth=0 presence_depth=50 presence_depth=100 presence_depth=150
#> [1,]               NA                 1                 NA                 NA
#> [2,]               NA                NA                 NA                 NA
#> [3,]               NA                NA                 NA                 NA
#>      presence_depth=200
#> [1,]                 NA
#> [2,]                 NA
#> [3,]                 NA

# Two sets of SpatVector3D polygons give the polygons they share.
square <- function(x0) {
  sf::st_polygon(list(rbind(c(x0, 0), c(x0 + 2, 0), c(x0 + 2, 2),
                            c(x0, 2), c(x0, 0))))
}
mpa <- as_vector3d(
  sf::st_sf(geometry = sf::st_sfc(square(150), crs = 4326)),
  depth_min = 0, depth_max = 200
)
habitat <- as_vector3d(
  sf::st_sf(geometry = sf::st_sfc(square(151), crs = 4326)),
  depth_min = 100, depth_max = 500
)
shared <- intersect_3d(mpa, habitat)
as.data.frame(shared)[c("feature_id_x", "feature_id_y", "depth_min",
                        "depth_max")]
#>   feature_id_x feature_id_y depth_min depth_max
#> 1            1            1       100       200
volume(shared)
#> [1] 2461.484
```
