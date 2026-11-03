# 3D points and lines with depth and time

A `SpatVector3D` holds features of one geometry type, reported by
[`geomtype()`](https://rspatial.github.io/terra/reference/geomtype.html):
`"points"` (a detection at 35 m) or `"lines"` (a dive track or a vessel
route). Each vertex has a longitude and latitude (WGS84), a depth in
positive metres downward, and optionally a time (UTC). Build one with
[`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md)
from a table or sf object.

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

## Details

It is the vector counterpart of
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
and
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md).
It is not built on
[`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
or sf: those keep attributes per feature, so a line's depth and time at
each vertex would not survive their operations.

## Slots

- `type`:

  `"points"` or `"lines"`.

- `vertices`:

  A data frame with one row per vertex: `feature_id`, `vertex` (the
  order within the feature), `lon`, `lat`, `depth`, `time` (POSIXct,
  UTC, or `NA` when the features have no time), and any other columns
  from the input.

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
