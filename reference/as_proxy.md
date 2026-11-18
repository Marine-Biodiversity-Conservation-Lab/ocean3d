# Describe a gridded netCDF variable without reading its values (experimental)

`as_proxy()` opens one netCDF file, reads the metadata of one variable
and its coordinates, and returns a
[SpatProxy](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatProxy-class.md).
The CF metadata is applied here, once, so later steps can rely on each
axis being in a fixed form. No values of the variable are read.

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

  A netCDF file path, or a connection opened with
  [`ncdf4::nc_open()`](https://rdrr.io/pkg/ncdf4/man/nc_open.html)
  (netCDF model). A connection is left open.

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

- **Longitude and latitude:** cell edges come from the coordinate's CF
  `bounds` variable, or halfway between neighbouring centres.

- **Depth:** converted to metres positive down from the coordinate's
  `units` and `positive` attributes. Edges come from the CF bounds, or
  halfway between centres, cut at the surface. A depth axis stored deep
  to shallow is kept in that order and recorded as such.

- **Time:** decoded with CFtime. Bounds are attached when the file has
  them; otherwise the axis has none (they aren't inferred). The time
  axis must be strictly increasing; a backwards or out-of-order axis is
  an error that suggests sorting the file, for example with CDO's
  `sorttimestamp`. Only the `standard`, `gregorian` and
  `proleptic_gregorian` calendars are supported, and climatological time
  axes aren't supported yet.

- **CRS:** EPSG:4326, unless the variable's `grid_mapping` gives another
  geographic CRS as WKT (`crs_wkt` or `spatial_ref`).

Regular longitude-latitude grids with coordinate values are supported;
curvilinear, projected and rotated-pole grids are rejected. Bad cell
bounds (cells that don't contain their centre, have no width or overlap)
stop with an error.

## Experimental

`SpatProxy` and `as_proxy()` are experimental. Their slots and behaviour
may change without notice while the core object is built, and no other
function accepts a `SpatProxy` yet.

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
```
