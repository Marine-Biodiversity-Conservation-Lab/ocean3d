# Buffer 3D points or lines into polygons

Draws a polygon around each point or line of a
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md),
with a depth range and an optional time interval, as
[`sf::st_buffer()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
does in 2D. A detection at 35 m with `radius_km = 10` and
`depth_window = 25` becomes a circle of 10 km radius from 10 to 60 m
deep.

## Usage

``` r
buffer_3d(x, radius_km, depth_window, time_window = NULL)
```

## Arguments

- x:

  A
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  of points or lines, from
  [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md).

- radius_km:

  Horizontal radius in kilometres: one positive number, below 10000.

- depth_window:

  Depth in metres added above and below each feature's depths: one
  number, 0 or more. The top stops at the surface (0 m). With `0`, the
  polygon spans just the feature's own depth range.

- time_window:

  Optional time added before and after each feature's times: a number of
  days, or a `difftime` (e.g. `as.difftime(6, units = "hours")`). Every
  vertex needs a time. If `x` has dates without hours, the window must
  be whole days, and the interval covers whole UTC days as in
  [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md).
  With `NULL` (default), each polygon keeps its feature's own times: a
  point's time, a line's first to last time, or no time. A line with
  times at some vertices but not others is an error either way.

## Value

A
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
of polygons, one per feature of `x`, with the same IDs. For points, the
other columns of `x` are kept. For lines, columns with one value along
each line (such as a tag name) are kept, and columns that change from
vertex to vertex (such as speed) are dropped.

## Details

The outline is a circle (for a line, the area within `radius_km` of it)
measured on the WGS84 ellipsoid, made with
[`terra::buffer()`](https://rspatial.github.io/terra/reference/buffer.html).
Outlines that cross the 180th meridian (the antimeridian) are split
there into two parts.

Each polygon has one depth range and one time interval, so a line's
polygon covers the line's whole depth and time range. Depth runs from
`depth_window` above the shallowest vertex to `depth_window` below the
deepest, and time from `time_window` before the first time to
`time_window` after the last. For a dive from the surface to 300 m, the
polygon spans 0 m to 300 m plus `depth_window` along the whole track. To
follow the dive more closely, split the line into shorter lines or
points and buffer those.

## See also

[`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md),
[`geomtype()`](https://rspatial.github.io/terra/reference/geomtype.html)

## Examples

``` r
obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(35, 120),
                  date = as.Date(c("2020-01-01", "2020-01-02")))
detections <- as_vector3d(obs, time = "date")

# 10 km circles, 25 m above and below each detection, a day either side.
buffers <- buffer_3d(detections, radius_km = 10, depth_window = 25,
                     time_window = 1)
as.data.frame(buffers)[c("feature_id", "depth_min", "depth_max",
                         "time_start", "time_end")]
#>   feature_id depth_min depth_max time_start   time_end
#> 1          1        10        60 2019-12-31 2020-01-02
#> 2          2        95       145 2020-01-01 2020-01-03

# A corridor 5 km either side of a dive track.
dive <- data.frame(tag = "A", lon = c(150, 150.1, 150.2), lat = -20,
                   depth = c(0, 120, 15))
track <- as_vector3d(dive, type = "lines", id = "tag")
buffer_3d(track, radius_km = 5, depth_window = 10)
#> class      : SpatVector3D
#> geometry   : polygons 
#> features   : 1 
#> depth (m)  : 0 to 130 
#> time       : none
```
