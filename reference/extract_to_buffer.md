# Summarise environmental values around observations

Calculate one summary per observation and variable within horizontal,
depth and/or temporal windows. Combine complementary netCDF sources
before assessing coverage. Observations and input order are retained.

## Usage

``` r
extract_to_buffer(
  data,
  nc,
  var = NULL,
  horizontal_radius_km = 0,
  lon_cell_radius = NULL,
  lat_cell_radius = NULL,
  depth_radius_m = 0,
  depth_layer_radius = NULL,
  time_before = 0,
  time_after = 0,
  time_unit = c("days", "hours", "minutes", "seconds", "months", "years"),
  statistic = "mean",
  weighting = c("weighted", "equal"),
  min_valid_fraction = 0.8,
  infer_bounds = FALSE,
  method = NULL,
  lon_col = "lon",
  lat_col = "lat",
  depth_col = "depth",
  date_col = "date",
  id_col = NULL,
  file_col = "file",
  lon_dim = NULL,
  lat_dim = NULL,
  depth_dim = NULL,
  time_dim = NULL,
  time_match = c("nearest", "exact"),
  time_resolution = c("auto", "datetime", "daily", "monthly", "yearly"),
  max_time_diff = NULL,
  spatial_match = c("within", "nearest"),
  max_distance_km = NULL,
  output_col = NULL,
  return_match_info = TRUE,
  quiet = FALSE
)
```

## Arguments

- data:

  Observation table, matrix, named list or sf POINT object.

- nc:

  netCDF path, paths, opened ncdf4 connection, list, or file table.

- var:

  Variable name or vector of names. NULL detects one data variable per
  file, excluding referenced bounds and coordinate variables.

- horizontal_radius_km:

  Horizontal geodesic radius in kilometres. Positive radii below 10000
  km are supported.

- lon_cell_radius, lat_cell_radius:

  Optional integer radii in grid cells. A radius of one on both axes
  requests a 3 by 3 neighbourhood. Mutually exclusive with a positive
  horizontal_radius_km.

- depth_radius_m:

  Symmetric depth radius in metres, positive downward.

- depth_layer_radius:

  Optional integer radius in physical depth order. Mutually exclusive
  with a positive depth_radius_m.

- time_before, time_after:

  Non-negative temporal distances before and after the observation. Both
  use time_unit.

- time_unit:

  One of days, hours, minutes, seconds, months or years. Days are 24
  elapsed hours. Months and years require whole numbers and shift the
  UTC calendar, clamping to the destination month's last day.

- statistic:

  One of "mean", "sd", "se", "min" or "max". SD uses the sample variance
  correction. SE is SD divided by the square root of the effective
  sample size. Min and max ignore weight magnitudes.

- weighting:

  "weighted" uses intersecting area, thickness and duration; "equal"
  gives equal weight to valid samples. Coverage uses physical support in
  both modes. A point-selected axis contributes a factor of one.

- min_valid_fraction:

  Minimum valid fraction of the entire requested support, between zero
  and one. Empty or entirely missing buffers return NA.

- infer_bounds:

  Allow missing bounds to be estimated using neighbouring coordinate
  midpoints, with half-gap extrapolation at outer edges. Inferred depth
  support is clipped at zero metres. A singleton needs explicit bounds.

- method:

  Depth selection on an unbuffered vertical axis: nearest, surface,
  bottom or 2d. NULL uses 2d without a depth axis and nearest with one.
  Leave NULL when requesting a depth interval or layer neighbourhood.

- lon_col, lat_col, depth_col, date_col:

  Observation coordinate columns. Dates are required only for temporal
  variables. Date and ISO date-only text centre time buffers at midnight
  UTC; POSIXt preserves instants.

- id_col:

  Optional observation identifier column used in diagnostics.

- file_col:

  File-path column when nc is a table.

- lon_dim, lat_dim, depth_dim, time_dim:

  Optional netCDF dimension overrides.

- time_match, time_resolution, max_time_diff:

  Point-time selection options on an unbuffered time axis; see
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md).
  max_time_diff is in days and cannot be combined with a temporal
  window.

- spatial_match, max_distance_km:

  Point-cell or index-seed selection options; see
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md).
  They do not expand a kilometre buffer.

- output_col:

  Output name for one variable, or a named character vector mapping
  multiple variable names to distinct output columns.

- return_match_info:

  Add a list column named

  \_match_info per variable, containing coverage, source, bounds and
  selection diagnostics.

- quiet:

  Hide explanatory summaries. Time-zone warnings remain active.

## Value

The normalised observation table with one numeric summary column per
variable and, optionally, a diagnostic list column per output. sf inputs
retain their geometry. Existing columns are never overwritten.

## Support and limitations

For normalised weights p, corrected variance is sum(p \* (x - mean)^2)
divided by sum(p \* (1 - p)); effective sample size is 1 / sum(p^2).
Equal weighting reproduces ordinary sample SD and SE. SD and SE require
at least two valid samples; otherwise they return NA with a diagnostic.
SE assumes independent samples and does not correct spatial or temporal
autocorrelation. It is not an automatic estimate of environmental
uncertainty. All statistics, including min and max, require the
requested valid coverage.

At least one window must be positive. Bounds describe represented
support; inferred intervals do not establish that instantaneous values
are interval means. The calculation treats each sample as constant over
its assigned support. No interpolation or nearest fallback expands an
empty buffer.

Geographic grids must have separable one-dimensional longitude and
latitude axes. Circles are geodesic on WGS84 with 720 segments. Each
cell is weighted by the area of the circle it covers, measured on the
WGS84 ellipsoid with terra. This is a polygonal approximation,
independent of sf_use_s2(). Unsupported projections, calendars,
climatologies or transformed vertical coordinates require explicit
preprocessing.

Complementary sources are combined before applying min_valid_fraction.
Overlapping selected support is an error, including duplicate values or
a missing copy. Sources for one variable must use matching value units.
Point-time selection is resolved across eligible sources before
aggregation. Index neighbourhoods require matching grid/layer
coordinates and bounds across sources. If requested indices extend
beyond available coverage, the result is NA because the full requested
physical support is unknown.

## See also

[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
[`extract_to_area()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_area.md)

## Examples

``` r
if (FALSE) { # \dontrun{
observations <- data.frame(lon = 0, lat = 38, depth = 60,
                           date = as.Date("2020-01-15"))
extract_to_buffer(observations, "temperature.nc", var = "thetao",
                  horizontal_radius_km = 5, method = "surface")
extract_to_buffer(observations, "temperature.nc", var = "thetao",
                  depth_radius_m = 10, time_before = 3)
} # }
```
