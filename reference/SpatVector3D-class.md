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

# S4 method for class 'SpatVector3D,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'SpatVector3D'
ext(x, ...)

# S4 method for class 'SpatVector3D'
crs(x, proj = FALSE, describe = FALSE, parse = FALSE)

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

- i:

  For `[`: the features to keep, by position (negative numbers drop), by
  a logical vector, or by feature ID (character). The result keeps the
  features in the order of `i`, with their IDs and other columns. For
  points and lines, every vertex of each chosen feature is kept.

- j, drop:

  Not used; `j` must be left out.

- proj, describe, parse:

  For `crs()`: as in
  [`terra::crs()`](https://rspatial.github.io/terra/reference/crs.html).

- row.names, optional:

  Not used.

## Value

`[`: a `SpatVector3D` of the same type with the chosen features, which
can be none.

`ext()`: the features' horizontal extent as a
[terra::SpatExtent](https://rspatial.github.io/terra/reference/SpatExtent-class.html)
in longitude and latitude: the shortest range of longitudes that covers
every feature, going either way round the globe. A shape that crosses
the antimeridian keeps its extent across it, with `xmax` past 180 (for
example 179 to 181), so a
[`terra::crop()`](https://rspatial.github.io/terra/reference/crop.html)
with it reads only the cells near the antimeridian. A polygon around a
pole covers every longitude, from -180 to 180, and reaches the pole.
Lines go the short way between vertices, as in
[`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md).
A `SpatVector3D` with no features has no extent, and `ext()` stops.

`crs()`: WGS84 (EPSG:4326), the CRS of every `SpatVector3D`, in the form
chosen with `proj`, `describe` or `parse`.

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
[`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md),
[`vector3d_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vector3d_to_envelope.md),
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

# One feature at a time, by position or ID, as in a loop over shapes.
pts[2]
#> class      : SpatVector3D
#> geometry   : points 
#> features   : 1 
#> vertices   : 1 
#> depth (m)  : 35 to 35 
#> time       : none
as.data.frame(pts["1"])
#>   feature_id vertex lon lat depth time
#> 1          1      1 150 -20    10 <NA>

# The horizontal extent and CRS let terra functions take the shape as an
# area, for example crop() of a SpatProxy.
terra::ext(pts)
#> SpatExtent : 150, 151, -21, -20 (xmin, xmax, ymin, ymax)
terra::crs(pts, describe = TRUE)$name
#> [1] "WGS 84"
```
