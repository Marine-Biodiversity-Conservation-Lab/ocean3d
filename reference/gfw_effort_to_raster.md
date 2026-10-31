# Rasterise a GFW effort tibble onto a target grid

Turn the long-format apparent-fishing-hours tibble returned by
[`gfwr::gfw_ais_fishing_hours()`](https://rdrr.io/pkg/gfwr/man/gfw_ais_fishing_hours.html)
(formerly
[`gfwr::get_raster()`](https://rdrr.io/pkg/gfwr/man/gfw_renamed.html))
into a multi-layer SpatRaster on the package's canonical study grid.
Each level of `layer_by` becomes its own layer, named `effort_<level>`
(e.g. `effort_drifting_longlines`).

## Usage

``` r
gfw_effort_to_raster(
  effort,
  grid = NULL,
  layer_by = "geartype",
  value = "Apparent Fishing Hours",
  fun = "sum"
)
```

## Arguments

- effort:

  Data frame. Output of
  [`gfwr::gfw_ais_fishing_hours()`](https://rdrr.io/pkg/gfwr/man/gfw_ais_fishing_hours.html)
  (a long-format tibble with at minimum `Lat`, `Lon`, a value, and a
  grouping column).

- grid:

  SpatRaster. Target grid (extent, resolution, CRS) — typically the same
  grid used for species ranges and WOA extraction. If `grid = NULL`,
  assume target grid from input effort data frame.

- layer_by:

  Character. Column in `effort` whose levels become layers. `NULL`
  produces a single-layer total-effort raster. Default `"geartype"`.

- value:

  Character. Column in `effort` to aggregate. Default
  `"Apparent Fishing Hours"`.

- fun:

  Character or function. Aggregation applied to records that fall into
  the same cell × layer level. Default `"sum"`.

## Value

A SpatRaster with one layer per `layer_by` level (or one layer if
`layer_by = NULL`). Layer names follow `effort_<level>`.

## Details

The input is expected to carry a cell centroid (`Lat`, `Lon`), a value
column (default `"Apparent Fishing Hours"`), and one categorical column
matching the API's `group_by` — for example `geartype` or `flag` (note
lower-case; this is what `gfwr` actually returns). Records that fall
into the same target cell × layer level are aggregated with `fun`.

## Examples

``` r
if (FALSE) { # \dontrun{
jpn_eez_id <- gfwr::gfw_region_id(region = "JPN", region_source = "EEZ")$id
effort <- gfwr::gfw_ais_fishing_hours(
  spatial_resolution  = "LOW",        # 0.01 deg
  temporal_resolution = "YEARLY",      # one bucket per year in [start, end]
  start_date          = "2022-01-01",
  end_date            = "2022-12-31",
  region              = jpn_eez_id,
  region_source       = "EEZ",
  group_by            = "GEARTYPE"
)

# Returns data.frame 
head(effort)

# Create multi-layer SpatRaster from data.frame, one for each geartype
effort_by_gear <- gfw_effort_to_raster(
  effort   = effort,
  layer_by = "geartype",
  value    = "Apparent Fishing Hours",
  fun      = "sum"
)

# Plot drifitng longlines effort
terra::plot(effort_by_gear$effort_drifting_longlines)
} # }
```
