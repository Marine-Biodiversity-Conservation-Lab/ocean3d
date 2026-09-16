# Keep a raster's values inside a 3D domain

These methods make
[`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html)
depth-aware. `mask(x, mask)` keeps the values of `x` where `mask` is
present and sets the rest to `NA`. When `mask` is a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
or a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
"present" is decided cell by cell **and** depth by depth. A value of `x`
stays only where the mask domain reaches that cell at that depth.

## Usage

``` r
# S4 method for class 'SpatVoxel,SpatEnvelope'
mask(x, mask, bounds = c("top", "midpoint"), ...)

# S4 method for class 'SpatVoxel,SpatVoxel'
mask(x, mask, ...)

# S4 method for class 'SpatEnvelope,SpatEnvelope'
mask(x, mask, ...)

# S4 method for class 'SpatEnvelope,SpatVoxel'
mask(x, mask, bounds = c("top", "midpoint"), ...)

# S4 method for class 'SpatRaster,SpatEnvelope'
mask(x, mask, ...)

# S4 method for class 'SpatRaster,SpatVoxel'
mask(x, mask, ...)
```

## Arguments

- x:

  The raster whose values are kept: a
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
  a
  [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
  or a plain `SpatRaster`.

- mask:

  The 3D domain to keep `x` inside: a
  [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
  or a
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  on the same grid as `x`.

- bounds:

  Only when one side is an envelope and the other a voxel. How the
  envelope is placed on the voxel's depth levels: `"top"` (default) or
  `"midpoint"`. See
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).

- ...:

  Passed on to
  [`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html).
  For example, `inverse = TRUE` keeps the values *outside* the domain
  instead.

## Value

`x`, with the same class and layers, with every value outside the domain
set to `NA`.

## Details

What each pairing does:

- voxel masked by an envelope:

  The envelope is first placed on the voxel's own depth levels with
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md).
  Each depth layer of the voxel is then masked by the matching level.
  This replaces the hand-written
  `terra::mask(v, envelope_to_voxel(e, depths(v)))` and cannot be
  pointed at the wrong depths.

- voxel masked by a voxel:

  Both must be sampled at the same depth levels. Each depth layer of `x`
  is masked by the matching layer of the mask, wherever the mask is
  occupied.

- envelope masked by an envelope or a voxel:

  A cell of `x` is kept where
  [`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
  finds the two domains overlap. Its depth interval is kept whole. To
  narrow the interval instead, use
  [`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md).

- plain `SpatRaster` masked by an envelope or a voxel:

  A cell is kept where the mask domain is present at any depth. Without
  this method terra would read the mask's first layer only.

A voxel or envelope masked by a plain `SpatRaster` or by polygons is
terra's own `mask()`. Every layer is masked by the same 2D pattern and
the class is kept. No method is added for those cases.

## See also

[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
for the shared domain itself;
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
for the yes/no test;
[`extract_to_area()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_area.md)
to restrict a voxel to a polygon and one depth band instead of a
per-cell depth window.

## Examples

``` r
fp <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3000,
                  ymin = 0, ymax = 1000,
                  crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
terra::values(fp) <- c(1, 1, 1)

# A temperature voxel on three levels, and a species found at 0-100 m in
# cells 1-2 only.
temp <- as_voxel(terra::rast(list(terra::setValues(fp, c(20, 21, 22)),
                                  terra::setValues(fp, c(15, 16, 17)),
                                  terra::setValues(fp, c(10, 11, 12)))),
                 depths = c(0, 100, 200), varname = "temp")
range_env <- as_envelope(terra::setValues(fp, c(1, 1, NA)),
                         depth_min = 0, depth_max = 100)

# The temperatures inside the species' range: cell 3 and the 200 m level
# are gone.
terra::values(mask(temp, range_env))
#>      temp_depth=0 temp_depth=100 temp_depth=200
#> [1,]           20             NA             NA
#> [2,]           21             NA             NA
#> [3,]           NA             NA             NA

# A 2D raster masked by the range keeps cells where the range is present
# at any depth.
terra::values(mask(fp, range_env))
#>      lyr.1
#> [1,]     1
#> [2,]     1
#> [3,]    NA
```
