# Create a voxel object

Wrap a multi-depth raster as a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md):
the validated 3D form used throughout the package, with one layer per
standard depth. Each layer's depth is stored in
[`terra::depth()`](https://rspatial.github.io/terra/reference/depth.html).
Given a
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md),
`as_voxel()` reads one time step of it into memory instead. Reading a
`SpatProxy` is experimental, like the class itself (see
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)),
and may change without notice while the core object is built.

## Usage

``` r
as_voxel(x, depths = NULL, varname = "value", depth_bounds = NULL, time = NULL)
```

## Arguments

- x:

  A
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md),
  a SpatRaster, a list of SpatRasters (one per depth, combined with
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html),
  which keeps their depths only when each has them set; otherwise pass
  `depths`), or an existing
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).
  This is how gridded data enter the `SpatVoxel` model: describe a
  netCDF file with
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
  and pass the proxy here, which reads one time step. Crop the proxy
  first with
  [`terra::crop()`](https://rspatial.github.io/terra/reference/crop.html)
  to read only part of the grid. netCDF paths aren't accepted directly.

- depths:

  Optional numeric vector, one depth per layer, in metres. When
  supplied, it sets the depths, and the layers are named
  `{varname}_depth={value}`, replacing any existing names. Required if
  `x` has no depths in
  [`terra::depth()`](https://rspatial.github.io/terra/reference/depth.html).
  Leave it out, and `depth_bounds` too, for a
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md),
  which gives its own.

- varname:

  Character. Variable name used in the layer names built from `depths`,
  or from a proxy's depths, where it defaults to the proxy's variable
  name; it is also stored in
  [`terra::varnames()`](https://rspatial.github.io/terra/reference/varnames.html),
  which renaming the layers keeps. Ignored otherwise.

- depth_bounds:

  Optional two-column matrix or data frame: the top and bottom of the
  slab of water each layer stands for, in metres positive down, one row
  per layer of `x` in its current order (for example a netCDF file's
  depth bounds). A data frame with columns `depth`, `lower` and `upper`,
  such as
  [`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md)
  of another voxel, is matched to the layers by depth instead. Each slab
  must contain its depth and have a positive height, and slabs must not
  overlap. `NULL` (the default) keeps the edges an existing `SpatVoxel`
  carries, unless `depths` renames its layers; a voxel without edges
  uses the `bounds` argument of the functions it is passed to.

- time:

  For a
  [SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)
  only: the time step to read, as a step number, or a date or date-time
  that matches one step exactly. A date such as `"2020-01-02"` matches
  the step on that day; with several steps a day, give a date-time.
  `NULL` (the default) reads the only step of a proxy that has one; a
  proxy with more needs `time`.

## Value

A
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
whose layers are ordered shallow to deep.

## Details

A voxel read from a proxy carries the proxy's depth edges (see
[`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md)):
the file's depth bounds, or, for a file without them (such as Copernicus
Marine), edges halfway between levels, from the surface to half a gap
below the deepest level. Because the voxel has edges, the `bounds`
argument of
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
and the other 3D functions doesn't apply to it, and each layer stands
for the same slab of water whichever layers are selected. To use other
edges, pass them with `as_voxel(v, depth_bounds = ...)`. A single level
without bounds has no edges, so `bounds` applies.

The depths come from `depths` when it is given, and otherwise from the
raster's own
[`terra::depth()`](https://rspatial.github.io/terra/reference/depth.html),
which terra also sets when it reads a netCDF file. A raster with neither
is an error: layer names are labels, never read as depths.

Most terra operations (`crop()`, `mask()`, `[[`, arithmetic,
[`c()`](https://rdrr.io/r/base/c.html)) keep the class and the depths; a
few, such as [`mean()`](https://rdrr.io/r/base/mean.html), return a
plain `SpatRaster` instead. Re-wrapping is cheap either way: given a
valid `SpatVoxel` and no `depths`, `as_voxel()` returns it untouched.

Keeping the class doesn't re-run the class checks, so an operation that
changes the layers can leave an object labelled `SpatVoxel` that isn't a
valid one (for example `v[[c(1, 1)]]`, two layers at the same depth).
`as_voxel()` checks such an object as it would a plain raster: unsorted
layers are sorted, and the rest is an error. A voxel whose depths terra
dropped, for example in `app()`, needs `depths` again (see
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)).

Depths are positive metres increasing downward, matching the World Ocean
Atlas convention. Negative depths are an error rather than being
silently negated, since flipping the sign would change what the data
mean.

## Examples

``` r
r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
terra::values(r) <- runif(terra::ncell(r) * 3)

# set the depths from a vector; the layers get readable names
v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")
depths(v)
#> [1]   0 100 200
names(v)
#> [1] "temp_depth=0"   "temp_depth=100" "temp_depth=200"

# depths already in terra::depth() are used as they stand, and sorted if
# needed
terra::depth(r) <- c(200, 0, 100)
depths(as_voxel(r))
#> [1]   0 100 200

# layer names aren't read as depths: a raster without terra::depth() needs
# `depths`
named <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2, vals = 1)
names(named) <- c("temp_depth=0", "temp_depth=50")
try(as_voxel(named))
#> Error : The raster has no depths in terra::depth(). Pass `depths` to as_voxel(), or set them with terra::depth(). Layer names aren't read as depths. terra drops the depths in some operations, such as app(); see ?`SpatVoxel-class`.

# carry each layer's slab of water, as a netCDF file's depth bounds give it
v <- as_voxel(r, depth_bounds = cbind(c(150, 0, 50), c(250, 50, 150)))
depth_bounds(v)
#>   depth lower upper
#> 1     0     0    50
#> 2   100    50   150
#> 3   200   150   250
```
