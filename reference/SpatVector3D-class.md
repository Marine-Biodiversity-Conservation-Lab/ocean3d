# 3D points, lines and polygons with depth and time

A `SpatVector3D` holds features of one geometry type, reported by
[`geomtype()`](https://rspatial.github.io/terra/reference/geomtype.html):
`"points"` (a detection at 35 m), `"lines"` (a dive track or a vessel
route) or `"polygons"` (a marine protected area from 0 to 200 m).
Coordinates are longitude and latitude (WGS84), depths are positive
metres downward, and times are UTC. Build one with
[`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md)
from a table, sf object or
[`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html).

## Usage

``` r
# S4 method for class 'SpatVector3D'
geomtype(x, ...)

# S4 method for class 'SpatVector3D'
length(x)

# S4 method for class 'SpatVector3D'
show(object)

# S3 method for class 'SpatVector3D'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `SpatVector3D`.

- ...:

  Not used.

- object:

  A `SpatVector3D`.

- row.names, optional:

  Not used.

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): the
vertex table for points and lines, or the polygon table, with its sf
geometry column, for polygons.

## Details

Points and lines have a depth and optional time at each vertex. Polygons
have one depth range (`depth_min` to `depth_max`) and an optional time
interval for the whole shape.

It is the vector counterpart of
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
and
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md).
Points and lines are not stored as
[`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
or sf objects: those keep attributes per feature, so a line's depth and
time at each vertex would not survive their operations. A polygon's
depth range and time are per feature, so polygons keep their outlines in
an sf geometry column.

## Slots

- `type`:

  `"points"`, `"lines"` or `"polygons"`.

- `vertices`:

  Points and lines: a data frame with one row per vertex: `feature_id`,
  `vertex` (the order within the feature), `lon`, `lat`, `depth`, `time`
  (POSIXct, UTC, or `NA` when the features have no time), and any other
  columns from the input. Empty for polygons.

- `features`:

  Polygons: a data frame with one row per polygon: `feature_id`,
  `depth_min`, `depth_max`, `time_start` and `time_end` (POSIXct, UTC,
  or `NA` when the polygons have no time; with `time_resolution`
  `"daily"`, the interval includes all of the end day), `geometry` (an
  sf geometry column of POLYGON or MULTIPOLYGON, EPSG:4326), and any
  other columns from the input. Empty for points and lines.

- `time_resolution`:

  `"daily"` or `"datetime"` for the precision of the times, as for
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
  or `NA` when there is no time.

## See also

[`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md),
[`geomtype()`](https://rspatial.github.io/terra/reference/geomtype.html)

## Examples

``` r
obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35))
pts <- as_vector3d(obs)
geomtype(pts)
#> [1] "points"
pts
#> class      : SpatVector3D
#> geometry   : points 
#> features   : 2 
#> vertices   : 2 
#> depth (m)  : 10 to 35 
#> time       : none
```
