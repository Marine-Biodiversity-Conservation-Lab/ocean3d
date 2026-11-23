# Describe gridded data without reading its values (experimental)

`as_proxy()` describes one variable of gridded data and returns a
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md).
The data can be one netCDF file, a series of netCDF files, or a
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
already in memory. For a netCDF file it reads the metadata of the
variable and its coordinates and applies the CF metadata once, so later
steps can rely on each axis being in a fixed form. No values of the
variable are read.

## Usage

``` r
as_proxy(
  x,
  var = NULL,
  lon_dim = NULL,
  lat_dim = NULL,
  depth_dim = NULL,
  time_dim = NULL
)
```

## Arguments

- x:

  A netCDF file path or URL (e.g. OPeNDAP or NCZarr), several paths for
  a series of files, a connection opened with
  [`ncdf4::nc_open()`](https://rdrr.io/pkg/ncdf4/man/nc_open.html) (left
  open), or a
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).
  `var` and the dimension arguments apply to files.

- var:

  Name of the variable. `NULL` (the default) detects it when the file
  holds exactly one data variable.

- lon_dim, lat_dim, depth_dim, time_dim:

  Names of the dimensions for each axis. `NULL` (the default) detects
  them from CF metadata (`axis`, `standard_name`, `units`) or common
  names.

## Value

A
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md).

## Details

Besides one file, two sources give the same kind of proxy:

- **A series of files** (several paths): the same variable on the same
  grid and depths in every file, each with its own dated time axis, such
  as one file per day. The files are put in time order and their time
  axes joined into one, in the first file's units. A series stops with
  an error when files repeat or overlap each other's time steps, when a
  file has no dated time axis, or when the files differ in grid (spatial
  tiles) or in the variable's units. The time axis keeps its bounds only
  when every file has them.

- **A
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  already in memory**: one time step (the voxel's
  [`terra::time()`](https://rspatial.github.io/terra/reference/time.html)
  when set), with the voxel's cells, depths and depth edges; a voxel
  without edges gets slabs halfway between its depths, as
  `bounds = "midpoint"` gives. The time must be a date or a date-time; a
  date is read as midnight UTC. A voxel whose time is a year, a month or
  a plain number stops with an error: give it a date-time, or remove its
  time with `terra::time(v) <- NULL`.
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  returns the voxel itself.

For netCDF files, each axis is put in a fixed form:

- **Longitude and latitude:** cell edges come from the coordinate's CF
  `bounds` variable, or halfway between neighbouring centres.

- **Depth:** converted to metres positive down from the coordinate's
  `units` and `positive` attributes. Edges come from the CF bounds, or
  halfway between centres, cut at the surface and running half a gap
  below the deepest centre. A voxel read with
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  keeps these edges. A depth axis stored deep to shallow is kept in that
  order and recorded as such.

- **Time:** decoded with CFtime. Bounds are attached when the file has
  them; otherwise the axis has none (they aren't inferred). The time
  axis must be strictly increasing; a backwards or out-of-order axis is
  an error that suggests sorting the file, for example with CDO's
  `sorttimestamp`. Only the `standard`, `gregorian` and
  `proleptic_gregorian` calendars are supported. A climatological time
  axis with one step, such as a World Ocean Atlas period, is recorded as
  a climatology without a date; climatologies with more steps aren't
  supported yet.

- **CRS:** EPSG:4326, unless the variable's `grid_mapping` gives another
  geographic CRS as WKT (`crs_wkt` or `spatial_ref`).

Regular longitude-latitude grids with coordinate values are supported;
curvilinear, projected and rotated-pole grids are rejected. Bad cell
bounds (cells that don't contain their centre, have no width or overlap)
stop with an error, and so does a variable with no values along a
dimension (for example an unlimited time dimension with no records).

## Experimental

`SpatProxy` and `as_proxy()` are experimental. Their slots and behaviour
may change without notice while the core object is built.

## See also

[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md)

## Examples

``` r
x <- as_proxy(system.file("extdata", "example_3d.nc", package = "ocean3d"))
x
#> class      : SpatProxy (experimental)
#> source     : example_3d.nc
#> variable   : temp (degrees_C)
#> dimensions : 2 x 2 x 3 x 3 (lon x lat x depth x time)
#> lon        : 0 to 1, 2 cells, edges inferred
#> lat        : 40 to 41, 2 cells, edges inferred
#> depth      : 0 to 100 m, 3 levels, shallow first, edges inferred
#> time       : 2020-01-01 to 2020-01-06, 3 steps, standard calendar, no bounds
#> crs        : EPSG:4326

# Each axis in its fixed form: depth in metres, time as a CFtime object.
x@axes$depth$centres
#> [1]   0  50 100
x@axes$time$cftime
#> CF calendar:
#>   Origin  : 2020-01-01T00:00:00
#>   Units   : days
#>   Type    : standard
#> 
#> Time series:
#>   Elements: [2020-01-01 .. 2020-01-06] (average of 2.500000 days between 3 elements)
#>   Bounds  : not set

# Read one time step into memory as a SpatVoxel.
v <- as_voxel(x, time = 1)
depths(v)
#> [1]   0  50 100
```
