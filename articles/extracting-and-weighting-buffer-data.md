# Extracting environmental data around observation points

Environmental observations are often associated with a single location,
depth and time. However, the conditions experienced by a marine animal
may extend beyond that exact point. An animal may move between
neighbouring cells, occupy several depths or encounter changing
temperatures over time. Similarly, a longline may span several
kilometres and remain deployed for hours. In these situations, we may
want to describe environmental conditions around the recorded location,
across relevant depths and over the duration of the event.

**How can we summarise these conditions, and why might two averages of
the same data give different results?**

In this tutorial, we use
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
with a small example of Mediterranean Sea potential temperature. We
build a buffer around one observation, compare two ways of calculating
its mean, and look at the steps
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
runs for us.

## Our research question

Imagine we have recorded a marine animal at:

- Longitude 0 degrees and latitude 38 degrees.
- A depth of 60 metres.
- Midnight UTC on 20 September 2026.

We would like to describe the temperature within 10 kilometres of this
location, between 50 and 70 metres deep, and from one day before to one
day after the observation. Our question concerns a **surrounding volume
over time**, not the single closest temperature value.

The data are a small netCDF subset of real modelled temperatures,
included with this article, so no download or API credentials are
needed. Their source is described at the end.

## Getting started

The article files are in `vignettes/articles/` in the source repository.
When running these examples by hand, use that directory as your working
directory.

``` r

library(ocean3d)

nc_file <- file.path("buffer-data", "temperature-subset.nc")
source(file.path("buffer-data", "figures.R")) # article_signature() below
```

We describe the file once with
[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md),
which reads its grid, depths and time steps but no values. The file has
no cell or depth bounds, so the edges of each cell and depth level are
placed halfway between neighbouring centres.

``` r

temperature <- as_proxy(nc_file, var = "thetao")
temperature
#> class      : SpatProxy (experimental)
#> source     : temperature-subset.nc
#> variable   : thetao (degrees_C)
#> dimensions : 9 x 8 x 6 x 5 (lon x lat x depth x time)
#> lon        : -0.166666 to 0.166667, 9 cells, edges inferred
#> lat        : 37.8542 to 38.1458, 8 cells, edges inferred
#> depth      : 46.6522 to 72.6237 m, 6 levels, shallow first, edges inferred
#> time       : 2026-09-18 to 2026-09-22, 5 steps, standard calendar, no bounds
#> crs        : EPSG:4326
```

Our observation goes in a table. Longitude and latitude are decimal
degrees, depth is in metres, positive downwards, and the time is in UTC.

``` r

observation <- data.frame(
  id = "observation_1",
  lon = 0,
  lat = 38,
  depth = 60,
  date = as.POSIXct("2026-09-20 00:00:00", tz = "UTC")
)
observation
#>              id lon lat depth       date
#> 1 observation_1   0  38    60 2026-09-20
```

## Step 1: a weighted mean

A horizontal radius of 10 km gives a circle around the location. A depth
radius of 10 m reaches from 50 to 70 m, and the time window runs one day
either side, so the buffer is a cylinder from 50 to 70 m over 48 hours.

``` r

weighted <- extract_to_buffer(
  observation, temperature,
  horizontal_radius_km = 10,
  depth_radius_m = 10,
  time_before = 1,
  time_after = 1,
  output_col = "temperature"
)
weighted[c("id", "temperature", "temperature_n")]
#>              id temperature temperature_n
#> 1 observation_1     15.9614           450
```

`temperature_n` is the number of values the mean is built on: every grid
cell the circle touches, at every depth level whose slab of water
overlaps 50 to 70 m, on each of the three daily time steps in the
window. A cell, level or time step counts in full once the buffer
reaches it.

Each value is weighted by the volume of water its cell and depth level
stand for. The depth levels in this part of the file thicken with depth,
from about 4.7 to 5.7 m, so a value from a deeper level counts for a
little more than one from a shallower level.

## Step 2: an equal mean

Let’s repeat the extraction, changing only the weighting.

``` r

equal <- extract_to_buffer(
  observation, temperature,
  horizontal_radius_km = 10,
  depth_radius_m = 10,
  time_before = 1,
  time_after = 1,
  weighting = "equal",
  output_col = "temperature"
)

data.frame(
  weighting = c("weighted", "equal"),
  mean = c(weighted$temperature, equal$temperature),
  n = c(weighted$temperature_n, equal$temperature_n)
)
#>   weighting     mean   n
#> 1  weighted 15.96140 450
#> 2     equal 15.98204 450
```

The same values give two slightly different means. With equal weighting,
every value counts once, so the thinner levels near the top of the
window count as much as the thicker ones below. Here the levels differ
little in thickness, so the means differ by only a few hundredths of a
degree. Use the weighted mean for the average temperature of the water
in the buffer, and the equal mean for the average of the grid values in
it.

## Step 3: the same steps written out

[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
is a shortcut for the build, intersect and summarise steps of an
analysis. Written out, the weighted mean above is:

``` r

# Build: the observation as a 3D point, then a 10 km circle from 50 to 70 m.
point <- as_vector3d(observation, time = "date")
buffer <- buffer_3d(point, radius_km = 10, depth_window = 10)

# Intersect: read the cells under the buffer for each day in the window, and
# keep the values inside it.
days <- as_voxel(crop(temperature, buffer),
                 time = as.POSIXct(c("2026-09-19", "2026-09-21"), tz = "UTC"))
inside <- lapply(days, function(v) suppressMessages(mask(v, buffer)))

# Summarise: one row per day, or all three days together.
summarise_3d(inside)
#>         time     mean   n
#> 1 2026-09-19 15.97945 150
#> 2 2026-09-20 15.98180 150
#> 3 2026-09-21 15.92295 150
summarise_3d(inside, by = "all")
#>   time    mean   n
#> 1 <NA> 15.9614 450
```

The pooled row is the mean
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
returned. Writing the steps out is useful when the same buffers are
summarised in several ways, or to see how the conditions change from one
day to the next.

## Step 4: other statistics

`statistic` also gives the standard deviation, the minimum or the
maximum. Each call gives one statistic per observation.

``` r

statistics <- c("mean", "sd", "min", "max")
data.frame(
  statistic = statistics,
  temperature = vapply(statistics, function(s) {
    extract_to_buffer(observation, temperature, horizontal_radius_km = 10,
                      depth_radius_m = 10, time_before = 1, time_after = 1,
                      statistic = s)$thetao
  }, numeric(1))
)
#>      statistic temperature
#> mean      mean  15.9613974
#> sd          sd   0.3928506
#> min        min  15.4408598
#> max        max  17.1730785
```

The standard deviation describes the spread of the values inside the
buffer. Values in neighbouring cells, levels and days are usually
correlated, so it isn’t a measure of the mean’s uncertainty.

## Step 5: buffers in fewer directions

A buffer doesn’t need to extend in all three directions. Without a depth
window, `method` picks one depth level in each cell the circle touches.
By default this is the level nearest the observation’s depth that has a
value, as
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
does; `"surface"` and `"bottom"` take the shallowest or deepest level
instead. Each cell then counts by its area. Without a time window, the
time step nearest the observation is used.

``` r

data.frame(
  method = c("nearest", "surface"),
  temperature = c(
    extract_to_buffer(observation, temperature,
                      horizontal_radius_km = 10)$thetao,
    extract_to_buffer(observation, temperature,
                      horizontal_radius_km = 10, method = "surface")$thetao
  )
)
#>    method temperature
#> 1 nearest    15.91890
#> 2 surface    17.28152
```

## When there is no value

An observation without a date, outside the grid, or with no time step in
its window gets `NA`, and its `_n` column is 0. One warning counts these
observations, by reason.

``` r

undated <- observation
undated$date <- as.POSIXct(NA, tz = "UTC")
extract_to_buffer(undated, temperature, horizontal_radius_km = 10)[
  c("id", "thetao", "thetao_n")
]
#> Warning: 1 of 1 observations get NA: 1 without a usable position, depth or
#> date.
#>              id thetao thetao_n
#> 1 observation_1     NA        0
```

To drop summaries built on too few values, filter on the `_n` column.

## Summary

We have described the conditions around an observation with a
horizontal, vertical and temporal buffer:

1.  We built a buffer around one observation.
2.  We compared a volume-weighted and an equal mean of the same values.
    Here they differ by a few hundredths of a degree, because the depth
    levels are similar in thickness.
3.  We wrote out the build, intersect and summarise steps that
    [`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
    runs.
4.  We calculated other statistics, and buffers in fewer directions.

### Next steps

- The [workflow
  article](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/workflow-values-in-shapes.md)
  shows the same steps for points, buffers, tracks and areas.
- Use
  [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
  or
  [`terra::extract()`](https://rspatial.github.io/terra/reference/extract.html)
  on a `SpatProxy`, when the question concerns a single cell and level
  rather than a surrounding buffer.

## Data sources and reproducibility

The example is a subset of **Mediterranean Sea Physics Reanalysis**,
product `MEDSEA_MULTIYEAR_PHY_006_004`, dataset
`cmems_mod_med_phy-temp_my_4.2km_P1D-m_202511`. It contains the original
`thetao` values and coordinate centres from a download made on 19
October 2026. The original metadata identifies CMCC as the producer and
describes the data as daily mean potential temperature.

E.U. Copernicus Marine Service Information;
<https://doi.org/10.48670/mds-00375>. The subset retains this
attribution. Source data remain subject to the [Copernicus Marine
Service
licence](https://marine.copernicus.eu/user-corner/service-commitments-and-licence).
The article does not imply endorsement by Copernicus Marine or CMCC.

The article’s
[`buffer-data/README.md`](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/blob/main/vignettes/articles/buffer-data/README.md)
records the subset and its preparation. The repository script
[`data-raw/buffer-article.R`](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/blob/main/data-raw/buffer-article.R)
regenerates the subset from the downloaded original. It keeps
neighbouring coordinates, so the inferred edges match the larger file.
It does not add bounds or alter the original temperature values.

## All arguments and defaults

Expand or collapse the complete extract_to_buffer() signature

This signature is generated from the current function. Arguments without
a default require input. Where a vector lists alternatives, its first
entry is the default; supply one alternative in a call.

``` r

extract_to_buffer(
  # Data and sources
  data, nc, var = NULL,

  # Horizontal buffer
  horizontal_radius_km = 0,

  # Depth buffer
  depth_radius_m = 0,

  # Time buffer
  time_before = 0, time_after = 0,
  time_unit = c("days", "hours", "minutes", "seconds"),

  # Statistics inside each buffer
  statistic = c("mean", "sd", "min", "max"), weighting = c("weighted", "equal"),

  # Depth selection
  method = NULL,

  # Observation and output columns
  lon_col = "lon", lat_col = "lat", depth_col = "depth", date_col = "date",
  output_col = NULL,

  # Source columns and dimensions
  lon_dim = NULL, lat_dim = NULL, depth_dim = NULL, time_dim = NULL,

  # Temporal selection
  time_match = c("nearest", "exact"), max_time_diff = NULL
)
```

See
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
for units, accepted inputs and the complete return description.
