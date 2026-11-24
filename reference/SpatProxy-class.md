# A gridded dataset described but not loaded (experimental)

A `SpatProxy` describes one variable of a netCDF file without reading
its values: where the data are, how they are stored, and each axis in a
fixed form. In terra, "proxy" means the same thing: a `SpatVectorProxy`
describes vector data on disk. Build one with
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md).

## Usage

``` r
# S4 method for class 'SpatProxy'
crop(x, y, ...)

# S4 method for class 'SpatProxy'
show(object)
```

## Arguments

- x:

  A `SpatProxy`.

- y:

  For `crop()`: the area to keep, as a
  [terra::SpatExtent](https://rspatial.github.io/terra/reference/SpatExtent-class.html)
  or anything with an extent
  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  [terra::SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html),
  an sf object, or a
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  such as a buffer from
  [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)).
  The result keeps the cells whose longitude and latitude ranges overlap
  the area's extent, still without reading values, so
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  then reads only those cells from the file. A polygon is cut to its
  extent; to keep only the cells inside it, mask the voxel with it
  afterwards. An area in another CRS is projected first, and on a
  longitude-latitude grid an extent in the other longitude convention (0
  to 360 or -180 to 180) finds its cells too. An area that crosses the
  edge of the grid's longitudes (the antimeridian on a -180 to 180 grid,
  or 0 degrees on a 0 to 360 grid) is read from both ends of the grid:
  the cropped grid runs on across the edge, for example from 179 to 181
  degrees, so
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  gives one voxel and
  [`mask()`](https://rspatial.github.io/terra/reference/mask.html) with
  the same shape keeps the right cells. A
  [SpatVector3D](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVector3D-class.md)
  across the antimeridian has such an extent (see its `ext()`). A
  `SpatVector3D` with no features has no extent, so `crop()` stops.

- ...:

  Not used.

- object:

  A `SpatProxy`.

## Details

A `SpatProxy` is the planned core data object of ocean3d (SPEC.md, "Data
model").
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
reads one time step of it into memory as a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md),
with the proxy's depths, depth edges, CRS and date.

## Slots

- `source`:

  A list with `type` and where the data are: for `"netcdf"`, `path`, the
  full path of the file (a path rather than an open connection, so the
  object can be saved with
  [`base::saveRDS()`](https://rdrr.io/r/base/readRDS.html)); for
  `"series"`, `path`, the files in time order, and `steps`, the number
  of time steps in each; for `"voxel"`, `voxel`, the
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  in memory.

- `variable`:

  A list: `name`, `units`, `long_name`, `dims` (the variable's dimension
  names in ncdf4 order), `dim_lengths`, `axis_dims` (the dimension of
  each axis), `singleton_dims` (other dimensions of length 1) and
  `chunksizes` (`NULL` for unchunked storage).

- `axes`:

  A list with one element per axis, `NULL` for an axis the variable
  doesn't have. `lon`, `lat` and `depth` hold `dim`, `centres` and cell
  edges `lower` and `upper`, in storage order, and `edges`: where the
  edges come from (`"bounds"`, `"inferred"`, or `"none"` for a single
  cell without bounds). `depth` is in metres positive down and also has
  `storage` (`"shallow_first"` or `"deep_first"`). `time` holds `dim`,
  `cftime` (a CFtime object, with the file's bounds attached when it has
  them) and `edges` (`"bounds"` or `"none"`). A one-step climatology,
  such as a World Ocean Atlas period, has `climatology = TRUE`, its
  `units` and no `cftime`.

- `crs`:

  The coordinate reference system, as text.

## Experimental

`SpatProxy` and
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
are experimental. Their slots and behaviour may change without notice
while the core object is built.

## See also

[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
