# Sample 3D lines into points

Turns each line of a
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md),
such as a dive track or a vessel route, into points along it: its
vertices, and with a spacing or step, points in between. Each point has
the depth and time interpolated along the segment it lies on. The points
can then go into the intersect step
([`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
[`mask()`](https://rspatial.github.io/terra/reference/mask.html)) or be
buffered with
[`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md).

## Usage

``` r
sample_3d(
  x,
  sampling = NULL,
  spacing_km = NULL,
  depth_step_m = NULL,
  time_step = NULL,
  time_unit = c("hours", "minutes", "seconds", "days", "months", "years"),
  max_samples = 100000L
)
```

## Arguments

- x:

  A
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  of lines, from
  [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md).

- sampling:

  How to place the points:

  - `"vertices"` keeps only the vertices.

  - `"balanced"` spaces points evenly along each segment, as many as the
    strictest step limit needs. With calendar steps, it keeps the step
    dates and spaces points evenly between them.

  - `"merged_steps"` places the evenly spaced points that the distance
    and depth limits need, and separately a point at each time step, and
    keeps both sets, so the points can be unevenly spaced.

  `NULL` (default) gives `"balanced"` when a step is given and
  `"vertices"` otherwise.

- spacing_km:

  Largest horizontal distance between points, in kilometres.

- depth_step_m:

  Largest change in depth between points, in metres.

- time_step, time_unit:

  Largest change in time between points, and its unit: seconds, minutes,
  hours, days, or whole calendar months or years. Needs a time at every
  vertex. If the lines have dates without a time of day, the step must
  be whole days, months or years. Calendar steps count from the first
  vertex of each segment, so a segment that starts on 31 January with
  monthly steps gets 29 February, 31 March, 30 April.

- max_samples:

  Largest total number of points, across all lines. `sample_3d()` stops
  with an error if the steps need more.

## Value

A
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
of points, one per sample, in order along each line. Each point's
feature ID is the line's ID and the sample number (`"A_3"`). Columns:
`track_id` (the line's ID), `sample_id` (1, 2, ... along the line),
`distance_km` (distance from the line's start), `support_km` and
`support_seconds` (half the length and duration of the segments either
side, as weights for summaries along the track; `NA` without times),
`is_interpolated` (`FALSE` for the line's own vertices), and the line's
columns that have one value along it. If the lines have dates without a
time of day, each point's time is the UTC date it falls on.

## Details

Points between vertices follow the WGS84 geodesic (the shortest path on
the ellipsoid), so a segment across the antimeridian goes the short way.
Depth and time change linearly with distance along each segment.
Distances are horizontal, not three-dimensional.

A line with a time at some vertices but not others is an error: its time
between those vertices isn't known. Give every vertex a time, or none.

## See also

[`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)
to turn the points into polygons;
[`extract_to_track()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_track.md),
which samples tracks with `sample_3d()` and summarises the values along
them.

## Examples

``` r
dive <- data.frame(tag = "A", lon = c(150, 150.1, 150.2), lat = -20,
                   depth = c(0, 120, 15),
                   time = as.POSIXct(c("2020-01-01 00:00",
                                       "2020-01-01 02:00",
                                       "2020-01-01 03:00"), tz = "UTC"))
track <- as_vector3d(dive, type = "lines", id = "tag", time = "time")

# A point at least every 2 km, with depth and time along each segment.
samples <- sample_3d(track, spacing_km = 2)
as.data.frame(samples)[c("track_id", "sample_id", "distance_km", "depth",
                         "time")]
#>    track_id sample_id distance_km depth                time
#> 1         A         1    0.000000   0.0 2020-01-01 00:00:00
#> 2         A         2    1.744118  20.0 2020-01-01 00:20:00
#> 3         A         3    3.488236  40.0 2020-01-01 00:40:00
#> 4         A         4    5.232354  60.0 2020-01-01 01:00:00
#> 5         A         5    6.976472  80.0 2020-01-01 01:20:00
#> 6         A         6    8.720590 100.0 2020-01-01 01:40:00
#> 7         A         7   10.464708 120.0 2020-01-01 02:00:00
#> 8         A         8   12.208827 102.5 2020-01-01 02:10:00
#> 9         A         9   13.952945  85.0 2020-01-01 02:20:00
#> 10        A        10   15.697063  67.5 2020-01-01 02:30:00
#> 11        A        11   17.441181  50.0 2020-01-01 02:40:00
#> 12        A        12   19.185299  32.5 2020-01-01 02:50:00
#> 13        A        13   20.929417  15.0 2020-01-01 03:00:00

# Buffer the samples to follow the dive more closely than one buffer
# around the whole line.
buffer_3d(samples, radius_km = 1, depth_window = 10)
#> class      : SpatVector3D
#> geometry   : polygons 
#> features   : 13 
#> depth (m)  : 0 to 130 
#> time       : 2020-01-01 00:00:00 UTC to 2020-01-01 03:00:00 UTC 
```
