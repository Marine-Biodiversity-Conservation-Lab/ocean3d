# 3D voxel model: one layer per standard depth level

A `SpatVoxel` is a
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
in which **depth is the layer index** and the cell values are the
variable (temperature, oxygen, presence, ...). Each layer's depth is
stored in terra's own per-layer field,
[`terra::depth()`](https://rspatial.github.io/terra/reference/depth.html),
in metres positive down; depths must be set, distinct and ordered
shallow to deep. Read them with
[`depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depths.md).
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
names the layers `{variable}_depth={value}` when it sets the depths, as
labels for people to read. The package never reads depths from layer
names, so renaming the layers doesn't change their depths. The depth
axis is grid-wide: every cell is sampled at the same set of standard
depths.

## Usage

``` r
# S4 method for class 'SpatVoxel'
c(x, ...)

# S4 method for class 'SpatVoxel'
subset(x, subset, negate = FALSE, NSE = FALSE, ...)

# S4 method for class 'SpatVoxel,numeric'
x[[i, j]] <- value

# S4 method for class 'SpatVoxel,character'
x[[i, j]] <- value

# S4 method for class 'SpatVoxel'
show(object)
```

## Arguments

- x:

  A `SpatVoxel`.

- ...:

  For [`c()`](https://rdrr.io/r/base/c.html): further rasters to combine
  with `x`, as with terra's [`c()`](https://rdrr.io/r/base/c.html). The
  result is a `SpatVoxel` with each layer's depth and the depth edges of
  every `SpatVoxel` among them, or a plain `SpatRaster` when the layers
  don't make a valid voxel (for example two layers at the same depth).
  For [`subset()`](https://rdrr.io/r/base/subset.html): passed to
  terra's [`subset()`](https://rdrr.io/r/base/subset.html).

- subset, negate, NSE:

  For [`subset()`](https://rdrr.io/r/base/subset.html), which `[[` and
  `$` also call: the layers to keep, as with terra's
  [`subset()`](https://rdrr.io/r/base/subset.html). The result is a
  `SpatVoxel` that keeps each layer's depth through later operations.

- i, j, value:

  For `[[<-`: the layer to replace, by number or name, and its new
  values, as with terra. A replaced layer keeps its depth.

- object:

  A `SpatVoxel`. `show()` prints the class, the depth levels, whether
  the voxel carries depth edges (see
  [`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md))
  and the time, then terra's description of the grid.

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

The depths stay through selecting layers (`[[`, `$`,
[`subset()`](https://rdrr.io/r/base/subset.html)), replacing a layer
with `[[<-`, [`c()`](https://rdrr.io/r/base/c.html), arithmetic,
`ifel()`, `crop()` and `mask()`. The `SpatVoxel` methods for
[`subset()`](https://rdrr.io/r/base/subset.html), `[[<-` and
[`c()`](https://rdrr.io/r/base/c.html) keep them where terra alone would
reset them to 0. A few terra operations drop them, such as `app()`,
which also returns a plain `SpatRaster`. The package's functions then
stop with an error; set the depths again with
`as_voxel(x, depths = ...)` or
[`terra::depth()`](https://rspatial.github.io/terra/reference/depth.html).
A plain `SpatRaster` (for example after `methods::as(v, "SpatRaster")`)
gets terra's own behaviour: after a layer subset, its next operation
resets the depths.

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
# which sets the depths and sorts the layers shallow to deep.
r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
                 ymax = 1, nlyrs = 3)
terra::values(r) <- cbind(c(12, 11), c(9, NA), c(6, 5))
v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")

depths(v)         # depth is stored per layer, in terra::depth()
#> [1]   0 100 200
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
identical(depths(as_voxel(v)), depths(v))
#> [1] TRUE

# Layer names are labels: renaming keeps the depths.
names(v) <- c("surface", "mid", "deep")
depths(v)
#> [1]   0 100 200

# Selecting layers keeps their depths through later operations.
depths(terra::crop(v[[2:3]], terra::ext(0, 1, 0, 1)) * 2)
#> [1] 100 200

# The class checks run on construction: the depths must be set, distinct
# and run shallow to deep.
try(methods::new("SpatVoxel", r))
#> Error in validObject(.Object) : 
#>   invalid class “SpatVoxel” object: depths must be set in terra::depth(); set them with as_voxel(x, depths = ...) or terra::depth()
```
