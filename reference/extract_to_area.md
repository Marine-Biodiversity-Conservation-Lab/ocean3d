# Extract a 3D raster to an area and a depth band

Crop a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
to an area polygon and select the depth layers within a given depth
range. The nearest available depth layers to `min_depth` and `max_depth`
are used as the inclusive bounds, so the result always has at least one
layer; omit either to run to that end of the voxel.

## Usage

``` r
extract_to_area(area, rast_3d, min_depth = NULL, max_depth = NULL)
```

## Arguments

- area:

  sf, sfc, or SpatVector. Area polygon to crop the voxel to. Reprojected
  to the voxel's CRS when the two differ.

- rast_3d:

  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).
  Multi-depth raster, e.g. from
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  or
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).

- min_depth:

  Numeric. Shallowest depth in metres. `NULL` (the default) starts at
  the voxel's shallowest layer.

- max_depth:

  Numeric. Deepest depth in metres. `NULL` (the default) runs to the
  voxel's deepest layer.

## Value

A
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
cropped to `area` and filtered to the depth range.

## Details

This is the area counterpart of
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md):
the same "restrict a 3D source to a target geometry" operation, with a
polygon as the target rather than observation points.

To restrict a voxel to a species' *per-cell* depth window rather than
one depth band across the whole area, build a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
and mask with it instead: `mask(rast_3d, range_rast)`. See
[mask-3d](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md).

## See also

[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
for the point counterpart;
[`depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depths.md)
for the layer depths the bounds snap to.

## Examples

``` r
# A 4x4 voxel over five standard depths.
r <- terra::rast(nrows = 4, ncols = 4, xmin = -10, xmax = 10,
                 ymin = -10, ymax = 10, nlyrs = 5, crs = "EPSG:4326")
terra::values(r) <- seq_len(terra::ncell(r) * 5)
v <- as_voxel(r, depths = c(0, 50, 100, 500, 1000), varname = "t_an")

area <- terra::vect("POLYGON ((-5 -5, 5 -5, 5 5, -5 5, -5 -5))",
                    crs = "EPSG:4326")

# 40-600 m snaps to the 50, 100 and 500 m layers.
out <- extract_to_area(area, v, min_depth = 40, max_depth = 600)
names(out)
#> [1] "t_an_depth=50"  "t_an_depth=100" "t_an_depth=500"
terra::ext(out)   # cropped to the polygon
#> SpatExtent : -5, 5, -5, 5 (xmin, xmax, ymin, ymax)

# Omit the bounds to keep every layer.
names(extract_to_area(area, v))
#> [1] "t_an_depth=0"    "t_an_depth=50"   "t_an_depth=100"  "t_an_depth=500" 
#> [5] "t_an_depth=1000"
```
