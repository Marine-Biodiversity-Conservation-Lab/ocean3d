# Extract environmental values along tracks and lines

Sample ordered trajectories and reuse
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
for environmental selection. Return sampled positions and optional
summaries per track.

## Usage

``` r
extract_to_track(
  data,
  nc,
  var = NULL,
  track_col = NULL,
  sampling = NULL,
  spacing_km = NULL,
  interpolate_depth = FALSE,
  interpolate_time = FALSE,
  depth_step_m = NULL,
  time_step_seconds = NULL,
  method = NULL,
  statistics = "mean",
  weighting = c("distance", "duration", "equal"),
  min_valid_fraction = 0.8,
  lon_col = "lon",
  lat_col = "lat",
  depth_col = "depth",
  date_col = "date",
  output_col = NULL,
  time_resolution = c("auto", "datetime", "daily", "monthly", "yearly", "custom"),
  max_samples = 100000L,
  max_segment_km = NULL,
  max_gap_hours = NULL,
  quiet = FALSE,
  time_step = NULL,
  time_unit = c("hours", "minutes", "seconds", "days", "months", "years"),
  land_mask = NULL,
  marine_values = NULL,
  mask_type = c("land", "marine"),
  mask_domain = NULL,
  land_policy = c("skip", "error"),
  mask_step_km = 0.5,
  ...
)
```

## Arguments

- data:

  Ordered vertex table: a data frame, matrix, named list or sf POINT
  features, as for
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md).
  Rows are traversed in their supplied order within each track; they are
  never sorted by date. Convert sf LINESTRING features first, e.g. with
  `as.data.frame(sf::st_coordinates(sf::st_transform(lines, 4326)))` and
  `lon_col = "X", lat_col = "Y", track_col = "L1"`.

- nc, var:

  netCDF sources and selected variable names. See
  [`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md).

- track_col:

  Optional track identifier column. Without it, the input is one track.
  IDs are returned as text.

- sampling:

  NULL (default) uses 'balanced' when spacing_km, depth_step_m,
  time_step_seconds or time_step is supplied, and 'vertices' otherwise.
  'vertices' uses observations only; 'balanced' uses a shared
  subdivision satisfying all explicit step limits; 'merged_steps'
  combines distance/depth cuts with temporal cuts; 'cells' splits each
  connection at grid boundaries and samples interval midpoints,
  retaining original vertices. Repeated visits are retained. Cells mode
  currently uses spatial within matching.

- spacing_km:

  Maximum horizontal spacing in kilometres when inserting samples.
  Supply spacing or another step explicitly for balanced or
  merged_steps.

- interpolate_depth, interpolate_time:

  Allow linear interpolation of depth and elapsed UTC time between
  vertices. FALSE by default. Constant endpoint values are retained;
  missing endpoints produce missing intermediate values.

- depth_step_m, time_step_seconds:

  Optional maximum changes between samples. Require the corresponding
  interpolation switch. They supplement horizontal spacing, including
  for stationary or vertically moving tracks.

- method:

  Depth method passed to
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
  including 'all'. NULL uses nearest for variables with depth and 2d for
  variables without it.

- statistics:

  NULL for samples only, or distinct names from mean, sd, se, min and
  max. Summaries describe sampled environmental selections.

- weighting:

  'distance', 'duration' or 'equal'. Distance and duration assign half
  each adjacent segment's horizontal length or elapsed seconds to its
  endpoint samples. Equal gives each sample the same weight.

- min_valid_fraction:

  Minimum valid fraction of total requested track support for a summary,
  from zero to one. Missing values stay in its denominator.

- lon_col, lat_col, depth_col, date_col:

  Vertex coordinate columns. Dates are required for selected temporal
  variables, duration weights or time-gap limits.

- output_col:

  Output name or named variable-to-column mapping.

- time_resolution:

  Point-time precision; see
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md).
  With auto, the original observations determine precision, before
  interpolation.

- max_samples:

  Maximum total sampled positions, a positive integer.

- max_segment_km, max_gap_hours:

  Optional limits that reject long segments or known elapsed time gaps.
  Split disconnected trajectories explicitly. FALSE by default because
  geometric elevation may use a different convention.

- quiet:

  Hide explanatory reports; scientific and time-zone warnings remain.

- time_step, time_unit:

  Optional temporal step and units: seconds, minutes, hours, elapsed
  days, whole calendar months or years. Balanced sampling uses one equal
  subdivision per original segment for elapsed-time limits. Calendar
  steps advance from each preceding step in UTC, clamping to the
  destination month's end. Balanced sampling preserves these calendar
  anchors and subdivides each calendar span to satisfy distance/depth
  limits.

- land_mask:

  Optional sf polygon or single-layer terra SpatRaster mask. With no
  mask, marine connectivity is not checked. Geographic rasters must use
  longitudes from -180 to 180; rotate or reproject other conventions
  explicitly. Antimeridian connections are split.

- marine_values:

  Numeric raster codes representing marine cells, including the
  underlying codes of categorical rasters.

- mask_type:

  Meaning of polygon interiors: land or marine.

- mask_domain:

  Optional polygon domain for polygon masks. A land polygon alone cannot
  establish that its exterior is known sea. For rasters, coverage is
  defined by their extent and missing cells.

- land_policy:

  Skip invalid connections (default) or raise an error. Original
  observations remain; invalid support stays in the denominator.

- mask_step_km:

  Maximum segment length of the polyline used for mask predicates. This
  is an approximation to the geodesic; coastline accuracy also depends
  on the supplied mask. Touching land counts as an intersection.

- ...:

  Point-selection options: file_col, dimension overrides, time_match,
  max_time_diff, custom_time_interval, spatial_match, max_distance_km,
  and time_report, spatial_report, depth_report, verbose.

## Value

A list with samples (coordinates, original vertex attributes, source-row
references, sample support and extracted values), summary (one row per
track and output, or NULL), and settings. Tracks follow first appearance
order. Summary counts n_observed and n_assumed distinguish original
vertices from inserted positions; n_samples is their total, including
missing results and invalid connections. The same counts are available
once per track in settings\$sampling_counts, even when statistics =
NULL. Cell-interval midpoints count as assumed positions.
variation_status is not_requested when neither SD nor SE is requested,
available when the requested variation statistics are finite, or
insufficient_variation_samples when accepted support has fewer than two
valid samples. not_available covers rejected support or non-finite
statistics; consult status for coverage and support problems. status
follows extract_to_buffer(): empty (no known support), all_missing,
insufficient_coverage, partial or complete. Samples follow traversal
order. For vertices-only sampling the global input row order is
retained. Non-coordinate attributes of inserted samples are NA.
source_row indexes the input rows.

## Interpretation and limits

Positions follow WGS84 shortest geodesics. Distances are horizontal, not
three-dimensional travel distances. Antipodal endpoints need a waypoint.
Interpolated coordinates describe assumptions, not observations.

Balanced sampling approximates line integration. Cells mode finds
numerical intersections with rectilinear spatial cell boundaries and
selection changes in depth and time. It retains separate segment visits,
including returns to the same cell. Zero-support original vertices
remain available for inspection. Cells mode is not a volume-weighted
integration of physical averaging periods; the temporal rules still
compare timestamp labels. Custom temporal precision and spatial nearest
matching are currently unsupported in cells mode. Overlapping sources
with different horizontal grids require explicit alignment or another
sampling mode. Environmental values are selected using
extract_to_point(), never interpolated. Missing coordinates are rejected
and multipart lines must be split explicitly.

SD and SE use the same reliability-weight correction and effective
sample size as extract_to_buffer(). Adjacent samples are usually
correlated: SE is descriptive under an independence assumption, not an
uncertainty estimate. Densifying a track does not create independent
environmental observations. Distance-weighted summaries of stationary
tracks are NA; duration weighting requires known time support and
returns NA if any adjacent duration is unknown. settings records options
and candidate files. Actual selected records are stored in the
per-output match_info columns in samples.

## See also

[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)

## Examples

``` r
nc_file <- system.file("extdata", "example_2d.nc", package = "ocean3d")
track <- data.frame(
  lon = c(0, 1), lat = c(40, 40),
  date = as.Date(c("2020-01-01", "2020-01-01"))
)
result <- extract_to_track(track, nc_file,
  var = "temp", method = "2d",
  spacing_km = 20, quiet = TRUE
)
result$summary
#>   track_id output_column weighting n_samples n_observed n_assumed n_valid
#> 1        1          temp  distance         6          2         4       6
#>   valid_fraction total_distance_km total_duration_seconds effective_n   status
#> 1              1          85.39341                      0    5.555556 complete
#>   mean variation_status
#> 1 15.5    not_requested
```
