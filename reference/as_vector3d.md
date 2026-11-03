# Build 3D points or lines

Builds a
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
from a table with one row per vertex, or from an sf object with POINT or
LINESTRING geometries.

## Usage

``` r
as_vector3d(
  x,
  type = c("points", "lines"),
  depth = "depth",
  time = NULL,
  id = NULL,
  lon = "lon",
  lat = "lat"
)
```

## Arguments

- x:

  A data frame (or tibble, matrix, named list) with one row per vertex,
  or an sf object. sf points take their depth and time from columns; sf
  lines take them from the Z and M coordinates.

- type:

  `"points"` (default) or `"lines"`, for tables and sf points. With
  `"lines"`, the rows that share an `id` form one line, in row order. sf
  lines are always lines.

- depth:

  Column with the depth of each vertex, in positive metres downward. Not
  used for sf lines, whose depth is minus the Z coordinate (Z is height,
  as in other GIS tools).

- time:

  Optional column with the time of each vertex: `Date`, `POSIXct` or
  text, read as in
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md).
  For sf lines, the M coordinate in seconds since 1970-01-01 UTC, if
  present.

- id:

  Column with feature IDs. Required for lines from a table; for points,
  each row is a feature and IDs default to the row numbers.

- lon, lat:

  Longitude and latitude columns of a table.

## Value

A
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md).
Columns other than the coordinates, depth, time and ID are kept with
each vertex.

## See also

[`geomtype()`](https://rspatial.github.io/terra/reference/geomtype.html)

## Examples

``` r
obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35),
                  date = as.Date(c("2020-01-01", "2020-01-02")))
as_vector3d(obs, time = "date")
#> class      : SpatVector3D
#> geometry   : points 
#> features   : 2 
#> vertices   : 2 
#> depth (m)  : 10 to 35 
#> time       : 2020-01-01 to 2020-01-02 

dive <- data.frame(tag = "A", lon = c(150, 150.1, 150.2),
                   lat = -20, depth = c(0, 120, 15))
track <- as_vector3d(dive, type = "lines", id = "tag")
geomtype(track)
#> [1] "lines"
```
