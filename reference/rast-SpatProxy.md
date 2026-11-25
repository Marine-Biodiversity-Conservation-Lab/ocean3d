# Read a variable without depth from a SpatProxy

`rast()` reads time steps of a
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)
into memory as plain one-layer
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)s.
Use it for a variable without depth, such as sea-surface temperature or
chlorophyll: each step is one layer named after the variable. A variable
with depth stops with an error: read it with
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md),
which keeps its depth levels and their edges. So each kind of variable
has one reader.

## Usage

``` r
# S4 method for class 'SpatProxy'
rast(
  x,
  time = NULL,
  time_match = c("nearest", "exact"),
  max_time_diff = NULL,
  ...
)
```

## Arguments

- x:

  A
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md).

- time, time_match, max_time_diff:

  Which time steps to read, as in
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md):
  a step number, a date or date-time (matched with `time_match`, by
  default to the nearest step), a window of two dates or date-times, or
  a
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  with one feature.

- ...:

  Not used.

## Value

A one-layer
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
or a list of them, one per step, for a time window or a shape.

## Details

The proxy is read as
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
reads it. The time axis is decoded with CFtime rather than by terra,
which misreads the time of some files, and the date of each step is set
in
[`terra::time()`](https://rspatial.github.io/terra/reference/time.html).
A proxy cut down with `crop()` reads only its cells.

[`mask()`](https://rspatial.github.io/terra/reference/mask.html) by
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
polygons or points works on the result. It uses only their outlines or
positions, not their depths.
[`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md)
also takes the one-layer raster and weights each cell by its area.

## Experimental

`SpatProxy` and
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
are experimental. Their slots and behaviour may change without notice
while the core object is built.

## See also

[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
for a variable with depth.

## Examples

``` r
# A 2D temperature field on 1, 3 and 6 January 2020.
p <- as_proxy(system.file("extdata", "example_2d.nc", package = "ocean3d"))
r <- terra::rast(p, time = "2020-01-03")
r
#> class       : SpatRaster
#> size        : 2, 2, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : -0.5, 1.5, 39.5, 41.5  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> source(s)   : memory
#> varname     : temp
#> name        : temp
#> min value   : 15.5
#> max value   : 18.5
#> unit        : degrees_C
#> time        : 2020-01-03 00:00:00zUTC
terra::time(r)
#> [1] "2020-01-03 UTC"

# A window gives one raster per step inside it.
length(terra::rast(p, time = c("2020-01-01", "2020-01-03")))
#> [1] 2

# Summarised with each cell weighted by its area.
summarise_3d(r, statistic = c("mean", "sd"))
#>         time     mean       sd n
#> 1 2020-01-03 16.99266 1.290978 4
```
