# 3D voxel model: one layer per standard depth level

A `SpatVoxel` is a
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
in which **depth is the layer index** and the cell values are the
variable (temperature, oxygen, presence, ...). Layer names follow the
`{variable}_depth={value}` convention and must be distinct and ordered
shallow to deep. The depth axis is grid-wide: every cell is sampled at
the same set of standard depths.

## Usage

``` r
# S4 method for class 'SpatVoxel'
c(x, ...)
```

## Arguments

- x, ...:

  For [`c()`](https://rdrr.io/r/base/c.html): a `SpatVoxel` and further
  rasters to combine, as with terra's
  [`c()`](https://rdrr.io/r/base/c.html). The result keeps the depth
  edges of every `SpatVoxel` among them.

## Details

A voxel may have interior gaps — a cell can be NA at one depth and
non-NA at the depths above and below it.

Each layer stands for a slab of water around its depth. A voxel can
carry the top and bottom of each slab, for example a netCDF file's depth
bounds, set with the `depth_bounds` argument of
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
or
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
and read with
[`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md).
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
[`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md),
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md),
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
and [`mask()`](https://rspatial.github.io/terra/reference/mask.html) use
these edges when the voxel has them for every layer, and work the edges
out from the depths with their `bounds` argument when it doesn't. Each
slab must contain its depth and have a positive height, and slabs must
not overlap; gaps between them are allowed.

terra operations that keep the `SpatVoxel` class (`[[`, `crop()`,
`mask()`, arithmetic, [`c()`](https://rdrr.io/r/base/c.html)) keep the
edges; an operation that returns a plain `SpatRaster` (`app()`, writing
to a file) drops them. Combining a voxel that has edges with one that
doesn't, with [`c()`](https://rdrr.io/r/base/c.html), leaves some layers
without edges, so `bounds` applies to the result.

## Slots

- `depth_edges`:

  A data frame with columns `depth`, `lower` and `upper`, in metres
  positive down: the slab each depth stands for, from `lower` (its top)
  to `upper` (its bottom), as in
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md).
  It has no rows when the voxel carries no edges. Each layer finds its
  slab by its depth, so the table may hold rows for depths the voxel no
  longer has.

## See also

[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
[SpatVolume](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVolume-class.md)

## Examples

``` r
# Two cells sampled at three standard depths. Build one with as_voxel(),
# which names the layers for you and sorts them shallow to deep.
r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
                 ymax = 1, nlyrs = 3)
terra::values(r) <- cbind(c(12, 11), c(9, NA), c(6, 5))
v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")

names(v)          # depth lives in the layer name
#> [1] "temp_depth=0"   "temp_depth=100" "temp_depth=200"
terra::values(v)  # ...and the variable in the cell values
#>      temp_depth=0 temp_depth=100 temp_depth=200
#> [1,]           12              9              6
#> [2,]           11             NA              5

# Note cell 2: NA at 100 m, but with values above and below it. A voxel can
# hold that interior gap; a SpatEnvelope cannot.

# It is an ordinary SpatRaster underneath, so terra operations work on it
# unchanged. Re-wrapping with as_voxel() is cheap either way, since the
# constructor is idempotent.
methods::is(v, "SpatRaster")
#> [1] TRUE
terra::nlyr(terra::crop(v, terra::ext(0, 1, 0, 1)))
#> [1] 3
identical(names(as_voxel(v)), names(v))
#> [1] TRUE

# Validity is enforced on construction: every layer name must carry its
# depth, and the layers must run shallow to deep.
bad <- r
names(bad) <- c("a", "b", "c")
try(methods::new("SpatVoxel", bad))
#> Error in validObject(.Object) : 
#>   invalid class “SpatVoxel” object: every layer name must follow {variable}_depth={value}
```
