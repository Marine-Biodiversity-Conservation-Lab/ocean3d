# Convert 3D polygons to an envelope on a grid

Puts
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
polygons onto the grid of `template` as a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
as
[`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
does for 2D polygons with one depth range. Each cell a polygon touches
gets that polygon's depth range. Cells that no polygon touches are `NA`,
and so are cells where `depth_min` equals `depth_max`, which have no
thickness.

## Usage

``` r
vector3d_to_envelope(x, template)
```

## Arguments

- x:

  A
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  of polygons, from
  [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md)
  or
  [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md).

- template:

  A `SpatRaster` whose grid (extent, resolution and CRS) the envelope
  uses; its values are ignored. For example, a
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  or a raster from
  [`create_study_raster()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/create_study_raster.md).
  Polygons are projected to its CRS if it isn't WGS84. Longitudes from 0
  to 360 work too.

## Value

A
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
on the grid of `template`. Stops with an error if it would have no
cells: when the polygons don't overlap `template`, or none of them has
`depth_max` deeper than `depth_min`.

## Details

Where polygons overlap, a cell gets the shallowest `depth_min` and the
deepest `depth_max` of the polygons that cover it. An envelope has one
depth range per cell, so polygons from 0 to 50 m and from 200 to 300 m
over the same cell give 0 to 300 m there.

An envelope has no time, so the polygons' time intervals are dropped,
with a message when any polygon has one. To make an envelope for one
period, keep only the polygons for that period first.

## See also

[`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
for 2D polygons,
[`as_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_envelope.md)

## Examples

``` r
template <- terra::rast(nrows = 4, ncols = 4, xmin = 150, xmax = 154,
                        ymin = -24, ymax = -20, crs = "EPSG:4326")
square <- function(x0, y0) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 2, y0), c(x0 + 2, y0 + 2),
                            c(x0, y0 + 2), c(x0, y0))))
}
areas <- sf::st_sf(
  top = c(0, 200), bottom = c(50, 300),
  geometry = sf::st_sfc(square(151, -23), square(152, -24), crs = 4326)
)
polygons <- as_vector3d(areas, depth_min = "top", depth_max = "bottom")

# The overlapping cell gets 0 to 300 m.
e <- vector3d_to_envelope(polygons, template)
terra::values(e)
#>       depth_min depth_max
#>  [1,]        NA        NA
#>  [2,]        NA        NA
#>  [3,]        NA        NA
#>  [4,]        NA        NA
#>  [5,]        NA        NA
#>  [6,]         0        50
#>  [7,]         0        50
#>  [8,]        NA        NA
#>  [9,]        NA        NA
#> [10,]         0        50
#> [11,]         0       300
#> [12,]       200       300
#> [13,]        NA        NA
#> [14,]        NA        NA
#> [15,]       200       300
#> [16,]       200       300
```
