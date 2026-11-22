# A gridded dataset described but not loaded (experimental)

A `SpatProxy` describes one variable of a netCDF file without reading
its values: where the data are, how they are stored, and each axis in a
fixed form. In terra, "proxy" means the same thing: a `SpatVectorProxy`
describes vector data on disk. Build one with
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md).

## Usage

``` r
# S4 method for class 'SpatProxy'
show(object)
```

## Arguments

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
