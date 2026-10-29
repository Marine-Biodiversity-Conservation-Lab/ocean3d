# Extract values from netCDFs to one or more point observations

This is the general point-extraction interface. It can extract
environmental values for one or multiple observations. The user chooses
the type of extraction with the `method` argument. The more specific
functions
[`extract2d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract2d.md),
[`extract3d_surface()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_surface.md),
[`extract3d_bottom()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_bottom.md),
[`extract3d_nearest()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_nearest.md)
and
[`extract3d_all()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_all.md)
are wrappers around this function.

## Usage

``` r
extract_to_point(
  data = NULL,
  nc,
  lon = NULL,
  lat = NULL,
  depth = NULL,
  date = NULL,
  var = NULL,
  method = c("2d", "surface", "bottom", "nearest", "all"),
  lon_col = "lon",
  lat_col = "lat",
  date_col = "date",
  depth_col = "depth",
  id_col = NULL,
  file_col = "file",
  lon_dim = NULL,
  lat_dim = NULL,
  depth_dim = NULL,
  time_dim = NULL,
  time_match = c("nearest", "exact"),
  max_time_diff = NULL,
  output_prefix = NULL,
  output_col = NULL,
  verbose = FALSE,
  time_resolution = c("auto", "datetime", "daily", "monthly", "yearly", "custom"),
  custom_time_interval = NULL,
  time_report = TRUE
)
```

## Arguments

- data:

  Optional object containing observation points. Can be a data frame,
  tibble, `sf` object with POINT geometries, matrix with column names,
  or named list. If `NULL`, coordinates can be supplied directly through
  `lon`, `lat`, and optionally `depth` and `date`. For `sf` objects,
  longitude and latitude are obtained from the POINT geometries;
  projected geometries are transformed internally to longitude and
  latitude when a CRS is available.

- nc:

  netCDF source. Can be a file path, vector of file paths, list of file
  paths, opened `ncdf4` object, list of opened objects, or a data frame
  with a file column.

- lon:

  Optional numeric longitude value or vector for direct point extraction
  when `data = NULL`.

- lat:

  Optional numeric latitude value or vector for direct point extraction
  when `data = NULL`.

- depth:

  Optional numeric observation depth value or vector for direct point
  extraction. Required when `data = NULL` and `method = "nearest"` or
  `method = "all"`.

- date:

  Optional observation date or vector of dates for direct point
  extraction. Used when `data = NULL` and the netCDF variable has a time
  dimension. Accepts Date and POSIXt objects, or character values in
  "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS" format (UTC). Direct input
  vectors may have length 1 or a common compatible length.

- var:

  Name of the variable to extract from the netCDF file. If `NULL`, the
  function tries to detect the variable automatically from each file.
  This only works when each netCDF file contains one single variable.

- method:

  Extraction method. Options are `"2d"`, `"surface"`, `"bottom"`,
  `"nearest"` and `"all"`. With `method = "all"`, the function returns
  nearest, surface and bottom outputs together.

- lon_col:

  Name of the longitude column in `data`.

- lat_col:

  Name of the latitude column in `data`.

- date_col:

  Name of the date column in `data`. Required only when the netCDF
  variable has a time dimension.

- depth_col:

  Name of the observation-depth column in `data`. Required only when
  `method = "nearest"` or `method = "all"`.

- id_col:

  Optional observation identifier column. It is not required for
  extraction, but it is checked when provided to help users detect
  mistakes.

- file_col:

  Column containing file paths when `nc` is a data frame.

- lon_dim:

  Optional name of the longitude dimension in the netCDF file.

- lat_dim:

  Optional name of the latitude dimension in the netCDF file.

- depth_dim:

  Optional name of the depth dimension in the netCDF file.

- time_dim:

  Optional name of the time dimension in the netCDF file.

- time_match:

  Character. How observation dates are matched to netCDF time values.
  `"nearest"` (default) selects the closest available date; `"exact"`
  requires an exact date match and returns `NA` when no match exists.

- max_time_diff:

  Optional non-negative numeric value giving the maximum allowed
  difference, in days, between an observation date and the nearest
  netCDF date when `time_match = "nearest"`. If `NULL` (default), no
  maximum temporal distance is imposed. Ignored when
  `time_match = "exact"`.

- output_prefix:

  Optional name used for the extracted output column. If omitted, `var`
  is used.

- output_col:

  Optional alias for `output_prefix`. Use this when you want to define
  the exact name of the output column directly. If both `output_col` and
  `output_prefix` are provided, `output_col` is used. With
  `method = "all"`, this name is used as the suffix in `nearest_*`,
  `surface_*` and `seabottom_*` columns.

- verbose:

  Logical. If `TRUE`, progress messages are printed every 100 rows.

- time_resolution:

  Precision used to compare observation times with dataset timestamps.
  "auto" uses daily comparison for Date objects and date-only text, and
  datetime comparison for POSIXt objects and date-time text. Explicit
  alternatives are "datetime", "daily", "monthly", "yearly" and
  "custom". Calendar comparisons use UTC. These options compare
  timestamp labels; they do not infer averaging intervals or calculate
  temporal summaries. Multiple dataset records mapping to the same
  comparison value cause an error. With nearest matching, monthly and
  yearly comparisons use the start of the labelled month or year.
  max_time_diff remains in days and applies to these comparison values.
  Use exact matching to require the same labelled day, month or year.

- custom_time_interval:

  Positive difftime value used as the warning interval when
  time_resolution = "custom". Matching uses full timestamps. This
  interval does not reject candidates; use max_time_diff for that.

- time_report:

  Logical. If TRUE, the default, display a temporal summary when
  returned values are displaced in time or outputs are missing. The
  return format is unchanged. FALSE hides the summary but does not
  suppress temporal warnings.

## Value

If `data` is provided, the input observations with extracted values
appended. Data frames, tibbles and `sf` objects retain their structure;
matrices and named lists are returned as data frames. If `data = NULL`,
direct extraction returns a numeric value for one point and one output,
a numeric vector for multiple points and one output, a named numeric
vector for one point and multiple outputs (for example,
`method = "all"`), or a data frame for multiple points and multiple
outputs.

## Details

Numeric netCDF time coordinates must provide units expressed as
milliseconds, seconds, minutes, hours or days since a reference
date-time. Singular and plural unit names are accepted.

Reference date-times are interpreted in UTC. A trailing "Z" or "UTC" is
accepted; numeric timezone offsets such as "+02:00" are not currently
supported.

Supported calendars are "standard", "gregorian" and
"proleptic_gregorian". The first two are supported only for reference
dates and timestamps on or after 1582-10-15. Unsupported units,
calendars or missing time units produce an error.

Temporal resolution controls comparison of timestamp labels. It does not
identify averaging intervals, calculate temporal means, or interpret
climatological periods.

## Multiple sources

With multiple sources, only cells containing the observation are
eligible. Explicit spatial cell bounds are preferred. Otherwise, edges
are inferred from adjacent centres; singleton axes without bounds
require coincidence. Longitude comparisons support equivalent
coordinates across the antimeridian. Among eligible sources the closest
temporal candidate with a valid value wins, with input order breaking
ties. Existing output columns are replaced, including by NA. Output
names cannot overwrite extraction inputs or combine different variables
in the same column.

## Custom temporal intervals

With time_resolution = "custom", observations require hours as for
datetime. For example, custom_time_interval = as.difftime(3, units =
"hours") warns above three hours of displacement. Exact matching
requires the same timestamp. Auto matching remains unchanged and uses
the observed cadence for diagnostics.

## Temporal diagnostics

Diagnostics are emitted after source selection, for values actually
returned. Daily, monthly and yearly comparisons warn on a change of
calendar period in UTC. Datetime comparisons warn when displacement
exceeds one observed cadence step, or when a nonzero displacement has no
reliable cadence. With a detected daily cadence, datetime comparison
also warns when the selected value falls on a different UTC day, even
below 24 hours. Explanatory warnings include the requested and selected
timestamps, their difference, the comparison resolution, cadence, row,
variable and file name. Reports explain non-warning displacements and
missing outputs without repeating examples already explained by a
temporal warning. Cadence inference does not change matching or impose a
maximum time gap. Singleton sources may use the observed sequence of
compatible files for the same variable. Time bounds describe intervals,
not sampling frequency.

## Examples

``` r
nc_file <- system.file(
  "extdata", "example_3d.nc",
  package = "ocean3d"
)

# Extract the value closest to a single observation
extract_to_point(
  nc = nc_file,
  lon = 0,
  lat = 40,
  depth = 60,
  date = "2020-01-03",
  var = "temp",
  method = "nearest"
)
#> [1] 15.5

# Extract values for multiple observations
extract_to_point(
  nc = nc_file,
  lon = c(0, 1),
  lat = c(40, 41),
  depth = c(60, 100),
  date = "2020-01-03",
  var = "temp",
  method = "nearest"
)
#> [1] 15.5 13.5

# Match a complete observation timestamp.
# The example file has a reading at midnight on 2020-01-03.
extract_to_point(
  nc = nc_file,
  lon = 0,
  lat = 40,
  depth = 60,
  date = as.POSIXct("2020-01-03 00:15:00", tz = "UTC"),
  var = "temp",
  method = "nearest",
  time_resolution = "datetime",
  max_time_diff = 30 / 1440
)
#> Warning: Temporal warning: 1 observation(s) need attention.
#> 
#> A different timestamp was selected, but its cadence is unknown.
#> Row 1; variable: temp; output: temp.
#> Observation: 2020-01-03 00:15:00 UTC
#> Selected data: 2020-01-03 00:00:00 UTC
#> Difference: 15 minutes earlier.
#> Comparison resolution: datetime.
#> Detected cadence: irregular; no reliable sampling interval.
#> The available timestamps do not establish a reliable sampling interval. The time difference cannot be assessed against a detected cadence. Use time_resolution = 'custom' with custom_time_interval if you know the appropriate warning interval.
#> File: example_3d.nc
#> [1] 15.5

# Compare monthly labels when the file has one record per month.
# Here the timestamp labels are illustrative.
monthly_labels <- as.Date(c("2020-01-15", "2020-02-15"))
observation_date <- as.Date("2020-02-28")
format(monthly_labels, "%Y-%m") ==
  format(observation_date, "%Y-%m")
#> [1] FALSE  TRUE

# For extraction from such a file, use:
# time_resolution = "monthly", time_match = "exact".
# Multiple records in the same month cause an error;
# they are not averaged automatically.
```
