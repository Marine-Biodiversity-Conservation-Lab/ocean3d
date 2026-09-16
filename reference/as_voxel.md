# Create a voxel object

Wrap a multi-depth raster as a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md):
the validated 3D form used throughout the package, with one layer per
standard depth and layer names following the `{variable}_depth={value}`
convention.

## Usage

``` r
as_voxel(x, depths = NULL, varname = "value")
```

## Arguments

- x:

  SpatRaster with `{variable}_depth={value}` layer names, a list of
  single-depth SpatRasters, or an existing
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).

- depths:

  Optional numeric vector, one depth per layer, in metres. When
  supplied, layer names are (re)built from `depths` and `varname`,
  replacing any existing names. Required if `x` has no conforming layer
  names.

- varname:

  Character. Variable name used when building layer names from `depths`.
  Ignored when `depths` is `NULL`.

## Value

A
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
whose layers are ordered shallow to deep.

## Details

Most terra operations (`crop()`, `mask()`, `[[`, arithmetic) propagate
the class, so a `SpatVoxel` normally survives them; a few, such as
[`mean()`](https://rdrr.io/r/base/mean.html), return a plain
`SpatRaster` instead. Re-wrapping is cheap either way, since
`as_voxel()` is idempotent — given a valid `SpatVoxel` and no `depths`
it returns it untouched.

Propagating the class does not re-run the validity rules, so an
operation that changes the raster layer set can leave an object still
labelled `SpatVoxel` even if invalid. Passing an invalid `SpatVoxel` to
`as_voxel()` rebuilds from layer names, passing through the function as
if a plain multi-layer raster. Pass `depths` to rebuild the layer names.

Depths are positive metres increasing downward, matching the World Ocean
Atlas convention. Negative depths are an error rather than being
silently negated, since flipping the sign would change what the data
mean.

## Examples

``` r
r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
terra::values(r) <- runif(terra::ncell(r) * 3)

# build the layer names from a depth vector
v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")
names(v)
#> [1] "temp_depth=0"   "temp_depth=100" "temp_depth=200"

# already-conforming names are used as they stand, and sorted if needed
names(r) <- c("temp_depth=200", "temp_depth=0", "temp_depth=100")
names(as_voxel(r))
#> [1] "temp_depth=0"   "temp_depth=100" "temp_depth=200"
```
