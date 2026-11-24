# Convert 3D polygons to a voxel on a template's grid and depth levels

Puts
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
polygons onto the grid and depth levels of `template` as a presence
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md):
`1` at every cell and level a polygon covers, `NA` elsewhere. A polygon
covers a cell when it touches it, as in
[`vector3d_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vector3d_to_envelope.md),
and covers a level when its depth range shares some thickness with the
slab of water that level stands for, as in
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).
The slabs come from the template's depth edges (see
[`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md)),
or from `bounds` when it has none.

## Usage

``` r
vector3d_to_voxel(x, template, bounds = c("midpoint", "top"))
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

  A
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  whose grid, depth levels and depth edges the result uses; its values
  are ignored. Polygons are projected to its CRS if it isn't WGS84.
  Longitudes from 0 to 360 work too.

- bounds:

  How to place the edges of the slab each depth level stands for when
  `template` carries no depth edges: `"midpoint"` (default) or `"top"`,
  as in
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).

## Value

A presence
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
on the grid, depths and depth edges of `template`, with layers named
`presence_depth={value}`. Stops with an error if no cell and level would
be covered: when the polygons miss the grid, have no thickness, or lie
outside the template's levels.

## Details

Each level is filled on its own, so where two polygons cover one cell at
different depths, the levels between them stay empty: polygons from 0 to
50 m and from 200 to 300 m over the same cell cover the levels in those
two ranges and not the ones between.
[`vector3d_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vector3d_to_envelope.md)
gives that cell one range, 0 to 300 m, because an envelope holds one
interval per cell. Use this function when the next step works on a
voxel's levels, for example to compare the polygons with a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
of data;
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
and [`mask()`](https://rspatial.github.io/terra/reference/mask.html) use
it for polygons and a voxel.

Polygons with `depth_min` equal to `depth_max` have no thickness and
cover nothing. A voxel has no time, so the polygons' time intervals are
dropped, with a message when any polygon has one.

## See also

[`vector3d_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vector3d_to_envelope.md)
for one depth range per cell,
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)

## Examples

``` r
grid <- terra::rast(nrows = 4, ncols = 4, xmin = 150, xmax = 154,
                    ymin = -24, ymax = -20, crs = "EPSG:4326", nlyrs = 4,
                    vals = 1)
temp <- as_voxel(grid, depths = c(0, 100, 200, 300), varname = "temp")
square <- sf::st_polygon(list(rbind(c(151, -23), c(153, -23), c(153, -21),
                                    c(151, -21), c(151, -23))))
bands <- sf::st_sf(
  top = c(0, 200), bottom = c(50, 300),
  geometry = sf::st_sfc(square, square, crs = 4326)
)
polygons <- as_vector3d(bands, depth_min = "top", depth_max = "bottom")

# The 100 m level (50 to 150 m) lies between the two bands, so it stays
# empty; an envelope would fill it.
terra::values(vector3d_to_voxel(polygons, temp))[6, ]
#>   presence_depth=0 presence_depth=100 presence_depth=200 presence_depth=300 
#>                  1                 NA                  1                  1 
terra::values(envelope_to_voxel(vector3d_to_envelope(polygons, temp),
                                depths(temp)))[6, ]
#>   presence_depth=0 presence_depth=100 presence_depth=200 presence_depth=300 
#>                  1                  1                  1                  1 
```
