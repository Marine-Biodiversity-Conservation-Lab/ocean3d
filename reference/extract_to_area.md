# Extract a 3D raster to an area and a depth band

Keeps the values of a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
or
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)
inside an area polygon, and optionally a depth band, and crops the
result to the area. It is a shortcut for the build and intersect steps:
the data are described with
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md),
cropped to the area with
[`terra::crop()`](https://rspatial.github.io/terra/reference/crop.html)
and read with
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md),
and [`mask()`](https://rspatial.github.io/terra/reference/mask.html)
keeps the values inside the area, or inside
[SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
polygons from `min_depth` to `max_depth` built with
[`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md).
To summarise the result, pass it to
[`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md).

## Usage

``` r
extract_to_area(
  area,
  rast_3d,
  min_depth = NULL,
  max_depth = NULL,
  bounds = c("midpoint", "top"),
  time = NULL
)
```

## Arguments

- area:

  sf, sfc, or SpatVector polygons, with a CRS. Projected onto the data's
  grid, which may use longitudes from -180 to 180 or 0 to 360.

- rast_3d:

  A
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
  e.g. from
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  or
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
  or a
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)
  from
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md),
  of which only the cells under the area are read.

- min_depth, max_depth:

  The top and bottom of the depth band, in metres positive down. Give
  either or both; `NULL` for both (the default) keeps every level. With
  only one given, the band runs from the surface, or to the bottom of
  the deepest level. Under `bounds = "top"`, a voxel without depth edges
  gives its deepest level no thickness, so a band never reaches it.

- bounds:

  Where each level's depth range sits when the voxel has no depth edges:
  `"midpoint"` (default) or `"top"`, as in
  [`mask()`](https://rspatial.github.io/terra/reference/mask.html).
  Ignored when the voxel has depth edges, and without a band.

- time:

  The time step to read, as in
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md):
  needed only when the data have several time steps. A time window gives
  a list of voxels, one per step.

## Value

A
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
cropped to `area`, with the voxel's layers, `NA` outside the area and
the depth band. A list of them for a time window. A band that reaches
levels with no data in the area, such as water below the seafloor, gives
a voxel with no values and no warning.

## Details

A cell is inside the area when the polygon touches it. Without
`min_depth` and `max_depth`, every depth level of those cells has
values. With a band, a level has values where the band overlaps the
depth range of water the level represents by more than its edge. Those
depth ranges come from the voxel's depth edges when it has them (see
[`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md)),
and otherwise from `bounds`. Every layer of the voxel is returned;
levels outside the band, and cells outside the area, are `NA`. A band
that overlaps none of the levels gives a warning and a voxel with no
values.

Only the area's outline is used; its other columns are not returned.

To keep a voxel inside a depth window that changes from cell to cell,
such as a species range bounded by the seafloor, build a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
and mask with it instead: `mask(rast_3d, range_env)`. See
[mask-3d](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md).

## See also

[mask-3d](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md),
[`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md),
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md);
[`vignette("workflow-values-in-shapes")`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/workflow-values-in-shapes.md)
for the steps written out.

## Examples

``` r
# A 4x4 voxel over five standard depths.
r <- terra::rast(nrows = 4, ncols = 4, xmin = -10, xmax = 10,
                 ymin = -10, ymax = 10, nlyrs = 5, crs = "EPSG:4326")
terra::values(r) <- seq_len(terra::ncell(r) * 5)
v <- as_voxel(r, depths = c(0, 50, 100, 500, 1000), varname = "t_an")

area <- terra::vect("POLYGON ((-5 -5, 5 -5, 5 5, -5 5, -5 -5))",
                    crs = "EPSG:4326")

# 40-600 m overlaps the depth ranges of the 50, 100 and 500 m levels; the
# 0 and 1000 m layers are NA.
out <- extract_to_area(area, v, min_depth = 40, max_depth = 600)
terra::global(out, "notNA")
#>                 notNA
#> t_an_depth=0        0
#> t_an_depth=50       4
#> t_an_depth=100      4
#> t_an_depth=500      4
#> t_an_depth=1000     0
terra::ext(out)   # cropped to the polygon
#> SpatExtent : -5, 5, -5, 5 (xmin, xmax, ymin, ymax)

# The mean temperature in the area and band, weighted by volume.
summarise_3d(out)
#>   time     mean  n
#> 1 <NA> 49.32759 12
```
