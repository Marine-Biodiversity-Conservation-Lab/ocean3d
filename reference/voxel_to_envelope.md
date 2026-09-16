# Collapse Voxel 3D -\> Envelope 2.5D

Reduce a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
to the
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
that bounds it. A predicate `fun` is applied to the cell values at each
depth; for every cell, the shallowest depth at which the predicate is
`TRUE` becomes `depth_min` and the deepest becomes `depth_max`. Cells
where the predicate is never `TRUE` are `NA` in both layers.

## Usage

``` r
voxel_to_envelope(v, fun = function(x) !is.na(x))
```

## Arguments

- v:

  SpatVoxel (or a multi-depth SpatRaster with `{variable}_depth={value}`
  layer names).

- fun:

  Function taking a vector of cell values and returning a logical vector
  of the same length. `NA` results are treated as `FALSE`. Defaults to
  `\(x) !is.na(x)`.

## Value

A
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
with layers `depth_min` and `depth_max`, on the same grid as `v`.

## Details

The default predicate, `\(x) !is.na(x)`, gives the plain vertical extent
of the data: the shallowest and deepest depths at which the cell has any
value. Pass a different predicate to bound a subset of the values
instead, e.g. `\(x) x > 15` for the depths over which a cell exceeds 15
degrees.

**This conversion is lossy.** An envelope stores a single continuous
interval per cell, so any interior gap in the voxel is filled in: a cell
that satisfies `fun` at 0 m and 200 m but not at 100 m still yields the
envelope `[0, 200]`. Use
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
directly where interior gaps matter.

## See also

[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
the reverse expansion.

## Examples

``` r
# Two cells sampled at four standard depths.
r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1,
                 nlyrs = 4)
terra::values(r) <- cbind(c(12, NA), c(11, NA), c(NA, 8), c(6, NA))
v <- as_voxel(r, depths = c(0, 50, 100, 200), varname = "temp")

# Default predicate: the vertical extent of the data. Cell 1 has values at
# 0, 50 and 200 m, so it comes back as [0, 200] — the gap at 100 m is filled,
# because an envelope stores one continuous interval per cell.
terra::values(voxel_to_envelope(v))
#>      depth_min depth_max
#> [1,]         0       200
#> [2,]       100       100

# Another predicate bounds a subset of the values instead: here the depths
# over which a cell is warmer than 10 degrees. Cell 2 never qualifies, so it
# is NA in both layers.
terra::values(voxel_to_envelope(v, fun = function(x) x > 10))
#>      depth_min depth_max
#> [1,]         0        50
#> [2,]        NA        NA
```
