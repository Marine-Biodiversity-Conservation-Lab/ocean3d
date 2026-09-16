# Extracting 3D Environmental Data for Sharks and Rays

## Extracting 3D Environmental Data for Sharks and Rays

This vignette demonstrates the World Ocean Atlas (WOA) environmental
extraction analysis, which summarises temperature and dissolved oxygen
conditions across the 3D habitat of every shark and ray species. Here we
implement this workflow using `ocean3d` package functions.

### Overview

The goal is to produce a single summary table with one row per species,
containing statistics (min, max, mean, range) for temperature and
dissolved oxygen across each species’ geographic range and depth limits.
The workflow has five steps:

1.  **Fetch species depth limits** from the IUCN Red List API
2.  **Fill missing depths** using genus-level means
3.  **Load species range polygons**
4.  **Prepare WOA rasters and the study grid** (download, annual means,
    align to grid)
5.  **Extract and summarise** environmental conditions per species

### Setup

``` r

library(ocean3d)
library(dplyr)
library(sf)
library(stringr)
library(terra)
library(here)
```

### Step 1: Fetch Species Depth Limits from IUCN Red List

``` r

# Your IUCN Red List API key
api_key <- Sys.getenv("IUCN_REDLIST_KEY")

# Fetch assessments (including depth limits) for all sharks and rays
depth_table <- fetch_species_assessments(
  api_key = api_key,
  group_code = "sharks_and_rays" # comprehensive group code
)

head(depth_table)
#> # A tibble: 6 x 4
#>   sis_id scientific_name          upper_depth_limit lower_depth_limit
#>    <dbl> <chr>                                <dbl>             <dbl>
#> 1  12345 Carcharodon carcharias                   0              1200
#> 2  12346 Sphyrna lewini                           0               512
#> ...

saveRDS(depth_table, "tmp/depth_table.RDS")
# fetch_species_assessments for the entire sharks_and_rays group can take 30 min
# readRDS("tmp/depth_table.RDS")
```

### Step 2: Fill Missing Depth Values

Some species lack depth data on the IUCN Red List. This is handled by
[`fill_missing_depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/fill_missing_depths.md).

[`fill_missing_depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/fill_missing_depths.md)
is designed for use inside
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) —
it takes upper/lower/genus vectors and returns a two-column tibble of
filled values that
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) unpacks
into new columns.

``` r

depth_table_complete <- depth_table |>
  mutate(fill_missing_depths(upper_depth_limit, lower_depth_limit, genus_name))

depth_table_missing <- depth_table_complete |>
  filter(is.na(upper_depth), is.na(lower_depth))
# Species whose entire genus lacks depth info will still have NAs; drop them.
depth_table_complete <- depth_table_complete |>
  filter(!is.na(upper_depth), !is.na(lower_depth))
```

### Step 3: Load Species Range Polygons

The IUCN SHARKS_RAYS_CHIMAERAS shapefile is licensed data that must be
downloaded manually from the [IUCN Red List spatial data
portal](https://www.iucnredlist.org/resources/spatial-data-download).
Point `iucn_shp` at your local copy. SQL filtering at read time avoids
loading the full multi-GB layer into memory.

``` r

iucn_shp <- here("/home/jay/Programming_Projects/Big_Data/SHARKS_RAYS_CHIMAERAS/SHARKS_RAYS_CHIMAERAS.shp")
layer_name <- sf::st_layers(iucn_shp)$name[1]
sis_ids_sql <- str_c(depth_table_complete$sis_id, collapse = ", ")
species_ranges <- sf::st_read(
  iucn_shp,
  query = str_glue(
    "SELECT * FROM \"{layer_name}\" ",
    "WHERE id_no IN ({sis_ids_sql}) ",
    "AND presence IN (1) AND origin IN (1, 2) AND seasonal IN (1, 2, 3)"
  ),
  quiet = TRUE
)
```

### Step 4: Prepare WOA Rasters

[`woa_download()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_download.md)
fetches (and caches) WOA files and
[`woa_load_nc()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_load_nc.md)
loads individual files. The monthly summaries are produced by a
`woa_summarise_monthly()` helper defined in this article (Step 5b)
rather than by the package.

#### 4a: Download WOA NetCDFs

[`woa_download()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_download.md)
caches files in
\[[`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md)\]
(defaults to the standard R user cache directory), so subsequent runs
skip the download. Pass `output_dir` if you prefer to keep files in the
project.

``` r

t_annual_path <- woa_download("temperature", period = "annual", resolution = "0.25")
t_monthly_paths <- woa_download("temperature", period = "monthly", resolution = "0.25")

do_annual_path <- woa_download("dissolved_oxygen", period = "annual", resolution = "1")
do_monthly_paths <- woa_download("dissolved_oxygen", period = "monthly", resolution = "1")
```

#### 4b: Load annual mean rasters

``` r

t_annual <- woa_load_nc(t_annual_path, field = "an")
do_annual <- woa_load_nc(do_annual_path, field = "an")
```

#### 4c: Prepare the study grid and bathymetry

Per-cell depth windows are computed by
[`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md),
which takes the seafloor as one more `depth_max` constraint and so
clamps the species’ nominal lower limit to the bed wherever the bed is
shallower. We pick one environmental raster as the canonical study grid;
every other raster (bathymetry and all WOA variables at other
resolutions) gets projected onto it so the entire downstream pipeline
operates on a single, consistent grid.

``` r

# Bathymetry from GEBCO (download once from
# https://www.gebco.net/data_and_products/gridded_bathymetry_data/).
bathy <- load_gebco_bathymetry("/home/jay/Programming_Projects/Big_Data/gebco_2025_sub_ice_topo/GEBCO_2025_sub_ice.nc")

# Use the temperature 0.25 deg grid as the canonical study grid.
study_grid <- t_annual[[1]]

# Bathymetry → positive depth, aligned to the study grid.
# Project before negating: negating first forces an eager pass over the full
# global 15 arc-second grid, which terra otherwise avoids by reading only the
# window the target grid needs.
seafloor <- (terra::project(bathy, study_grid) * -1) |>
  terra::clamp(lower = 0)
```

#### 4d: Align the annual rasters onto the study grid

Masking a range envelope onto an environmental raster requires the two
to share extent, resolution, and CRS. Project each WOA raster onto
`study_grid` up front (nearest-neighbour preserves the
`{variable}_depth={value}` layer naming and the standard depth bins
without interpolation across depths). We align the annual means here;
the monthly summaries are aligned in Step 5 once they have been
computed.

``` r

# Any raster not already on the study grid gets reprojected, then re-wrapped as
# a SpatVoxel: terra operations propagate the class but do not re-run its
# validity rules, so as_voxel() is how the depth axis is re-checked after a
# projection. It is idempotent, so wrapping defensively costs nothing.
align_to_grid <- function(r, grid) {
  if (terra::compareGeom(r[[1]], grid, stopOnError = FALSE, messages = FALSE)) {
    as_voxel(r)
  } else {
    as_voxel(terra::project(r, grid, method = "near"))
  }
}

t_annual  <- align_to_grid(t_annual, study_grid)
do_annual <- align_to_grid(do_annual, study_grid)
```

### Step 5: Extract and Summarise Environmental Conditions per Species

#### 5a: Quick look — sea surface temperature in a species’ range

Before the full pipeline, here is the simplest possible extraction: the
sea surface temperature (depth = 0) values where a species’ range
overlaps the WOA temperature raster.
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
puts the range envelope on `t_annual`’s depth levels so
[`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html)
can restrict it to the cells and depths the species occupies; selecting
the `t_an_depth=0` layer leaves just the surface.

We loop over every species and save its surface temperature values,
keyed by SIS id. Each element of `sst_values` is the vector of non-`NA`
SST cells falling inside that species’ range.

``` r

sst_values <- lapply(seq_len(nrow(depth_table_complete)), function(i) {
  sp_id <- depth_table_complete$sis_id[i]
  sp_range <- species_ranges[species_ranges$id_no == sp_id, ]
  sp_depths <- depth_table_complete[i, ]

  # Per-cell depth window for this species, clamped to the seafloor.
  range_rast <- vect_to_envelope(
    polygon = sp_range,
    template = study_grid,
    depth_min = sp_depths$upper_depth,
    depth_max = list(sp_depths$lower_depth, seafloor)
  )

  # Surface temperature across the species/WOA overlap. mask() is
  # depth-aware for an envelope: each depth layer keeps only the cells whose
  # window reaches that depth.
  in_range <- mask(t_annual, range_rast)
  terra::values(in_range[["t_an_depth=0"]], na.rm = TRUE)[, 1]
})
names(sst_values) <- depth_table_complete$sis_id

# Persist the per-species surface temperatures for use outside the vignette.
saveRDS(sst_values, "sst_values.rds")

# Per-species summary of the saved surface temperatures.
lapply(sst_values, summary)
```

#### 5b: Summarise monthly extremes and assemble the raster list

The helper below reduces the twelve monthly files to min/max/diff
rasters per depth layer, replacing the entire `data-raw/WOA.R` script
(~170 lines) with two function calls.

It lives in this article rather than in `ocean3d` because the package is
moving to a single generic `temporal_summarise()` that will cover both
WOA and Copernicus Marine data. Define it in your session before running
the next chunk.

``` r

woa_summarise_monthly <- function(monthly_dir = NULL, field = "an",
                                  files = NULL) {
  if (is.null(files)) {
    if (is.null(monthly_dir) || !dir.exists(monthly_dir)) {
      stop("Directory not found: ", monthly_dir)
    }
    files <- list.files(monthly_dir, pattern = "\\.nc$",
                        full.names = TRUE, ignore.case = TRUE)
  }
  if (length(files) == 0) {
    stop("No .nc files found in: ", monthly_dir)
  }

  monthly <- lapply(files, woa_load_nc, field = field)
  depth_names <- names(monthly[[1]])

  # For each depth layer, stack that layer across all months, then reduce.
  per_depth <- lapply(depth_names, function(dn) {
    stk <- terra::rast(lapply(monthly, function(x) x[[dn]]))
    mx <- terra::app(stk, fun = "max", na.rm = TRUE)
    mn <- terra::app(stk, fun = "min", na.rm = TRUE)
    df <- mx - mn
    names(mx) <- dn
    names(mn) <- dn
    names(df) <- dn
    list(max = mx, min = mn, diff = df)
  })

  list(
    min  = terra::rast(lapply(per_depth, `[[`, "min")),
    max  = terra::rast(lapply(per_depth, `[[`, "max")),
    diff = terra::rast(lapply(per_depth, `[[`, "diff"))
  )
}
```

``` r

t_monthly  <- woa_summarise_monthly(files = t_monthly_paths,  field = "an")
do_monthly <- woa_summarise_monthly(files = do_monthly_paths, field = "an")
# Each result is a list: $min, $max, $diff
```

Assemble every environmental raster into a single named list. The annual
means were already aligned in Step 4; the monthly summaries are aligned
here so the whole list shares the study grid.

``` r

env_rasters <- list(
  temperature_annual_mean       = t_annual,
  temperature_monthly_max       = t_monthly$max,
  temperature_monthly_min       = t_monthly$min,
  temperature_monthly_diff      = t_monthly$diff,
  dissolved_oxygen_annual_mean  = do_annual,
  dissolved_oxygen_monthly_max  = do_monthly$max,
  dissolved_oxygen_monthly_min  = do_monthly$min,
  dissolved_oxygen_monthly_diff = do_monthly$diff
) |>
  lapply(align_to_grid, grid = study_grid)
```

The original analysis used a ~90-line `foreach` loop that manually
called `woa_volume_extract()` eight times per species (4 rasters x 2
variables), computed summary statistics, and assembled output rows. That
collapses to the two helpers below, applied across species with
[`lapply()`](https://rdrr.io/r/base/lapply.html).

#### 5c: Two helpers — mask, then summarise

`summarise_range()` is the whole per-species step: for each
environmental raster, restrict it to the species’ 3D range and reduce
the result to one row of statistics.

Each raster gets its own occupancy voxel, because each may sit on a
different set of standard depths — `depths(v)` is read per raster rather
than assumed shared.
[`terra::values()`](https://rspatial.github.io/terra/reference/values.html)
returns a cells × depths matrix, so `n_cells` counts occupied **voxels**
(cell × depth) while `n_surface_cells` counts the map cells with data at
any depth.

``` r

# One row of statistics per named raster, prefixed with that name.
summarise_range <- function(range_rast, raster_list) {
  stopifnot(!is.null(names(raster_list)), all(names(raster_list) != ""))

  cols <- Map(function(v, nm) {
    vals <- terra::values(mask(v, range_rast))
    out <- c(
      min             = suppressWarnings(min(vals, na.rm = TRUE)),
      max             = suppressWarnings(max(vals, na.rm = TRUE)),
      mean            = suppressWarnings(mean(vals, na.rm = TRUE)),
      n_surface_cells = sum(rowSums(!is.na(vals)) > 0),
      n_cells         = sum(!is.na(vals)),
      n_depths        = ncol(vals)
    )
    # A species whose range misses the raster entirely gives -Inf / Inf.
    out[is.infinite(out)] <- NA_real_
    names(out) <- paste(nm, names(out), sep = "_")
    as.list(out)
  }, raster_list, names(raster_list))

  as.data.frame(unlist(cols, recursive = FALSE), stringsAsFactors = FALSE)
}
```

#### 5d: Summarise for a single species

``` r

i <- 1
sp_id <- depth_table_complete$sis_id[i]
sp_range <- species_ranges[species_ranges$id_no == sp_id, ]
sp_depths <- depth_table_complete[i, ]

range_rast <- vect_to_envelope(
  polygon = sp_range,
  template = study_grid,
  depth_min = sp_depths$upper_depth,
  depth_max = list(sp_depths$lower_depth, seafloor)
)

result <- summarise_range(
  range_rast = range_rast,
  raster_list = env_rasters
)

result
#> # A tibble: 1 x 32
#>   temperature_annual_mean_min temperature_annual_mean_max ...
```

#### 5e: Apply across all species

``` r

sf_use_s2(FALSE)

results <- lapply(seq_len(nrow(depth_table_complete)), function(i) {
  sp_id <- depth_table_complete$sis_id[i]
  sp_range <- species_ranges[species_ranges$id_no == sp_id, ]
  sp_depths <- depth_table_complete[i, ]

  row <- tryCatch({
    range_rast <- vect_to_envelope(
      polygon = sp_range,
      template = study_grid,
      depth_min = sp_depths$upper_depth,
      depth_max = list(sp_depths$lower_depth, seafloor)
    )
    summarise_range(range_rast, env_rasters)
  },
  error = function(e) {
    message("Error for ", sp_depths$scientific_name, ": ", e$message)
    NULL
  })
  if (!is.null(row)) row$sis_id <- sp_id
  row
})

summary_table <- results |>
  Filter(f = Negate(is.null)) |>
  bind_rows() |>
  merge(depth_table_complete, by = "sis_id")
```

#### 5f: Save output

``` r

write.csv(summary_table, here("outputs/species_environment_summary.csv"),
          row.names = FALSE)
```

### Visualisation

The two helpers below are defined in this article rather than exported
by the package. They are thin `ggplot2` wrappers used only for these
examples.

``` r

library(ggplot2)
library(tidyterra)

# Vertical depth profile of an environmental variable within a species range.
# At each depth layer the spatial mean, min, and max are computed across the
# cells where that depth falls inside the species' per-cell depth window
# (bathymetry-clamped). Layers that are all-NA are dropped.
plot_depth_profile <- function(species_name, range_rast, rast_3d) {
  masked <- mask(rast_3d, range_rast)
  layer_depths <- depths(masked)

  summary_per_depth <- lapply(seq_len(terra::nlyr(masked)), function(i) {
    v <- terra::values(masked[[i]])
    if (all(is.na(v))) return(NULL)
    data.frame(
      depth = layer_depths[i],
      mean  = mean(v, na.rm = TRUE),
      min   = min(v, na.rm = TRUE),
      max   = max(v, na.rm = TRUE)
    )
  })
  df <- do.call(rbind, Filter(Negate(is.null), summary_per_depth))
  if (is.null(df) || nrow(df) == 0) {
    stop("No cells overlap the species range and the raster's depth layers.",
         call. = FALSE)
  }

  ggplot(df, aes(x = mean, y = depth)) +
    geom_ribbon(aes(xmin = min, xmax = max), alpha = 0.2) +
    geom_line() +
    geom_point() +
    scale_y_reverse() +
    labs(title = species_name, x = "Value", y = "Depth (m)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "italic"))
}

# Map view of a species range with environmental values at the depth layer
# nearest `depth`.
plot_range_at_depth <- function(species_range, depth, rast_3d) {
  layer_depths <- depths(rast_3d)
  idx <- which.min(abs(layer_depths - depth))
  layer <- rast_3d[[idx]]
  actual_depth <- layer_depths[idx]

  range_vect <- if (inherits(species_range, "sf")) {
    terra::vect(species_range)
  } else {
    species_range
  }
  if (terra::crs(range_vect) != terra::crs(layer) &&
      terra::crs(range_vect) != "") {
    range_vect <- terra::project(range_vect, terra::crs(layer))
  }

  cropped <- terra::crop(layer, range_vect, mask = TRUE, touches = TRUE)
  range_sf <- if (inherits(species_range, "sf")) {
    species_range
  } else {
    sf::st_as_sf(range_vect)
  }

  ggplot() +
    geom_spatraster(data = cropped) +
    geom_sf(data = range_sf, fill = NA, colour = "black", linewidth = 0.3) +
    scale_fill_viridis_c(name = NULL, na.value = "transparent") +
    coord_sf() +
    labs(subtitle = paste0("Depth: ", actual_depth, " m")) +
    theme_minimal() +
    theme(axis.title = element_blank())
}
```

#### Depth profile

``` r

# Temperature depth profile within the white shark's range.
# Reuse the per-species vect_to_envelope() output from step 5e.
sp_range <- species_ranges[species_ranges$sci_name == "Carcharhinus longimanus", ]
sp_range_rast <- vect_to_envelope(
  polygon = sp_range, template = study_grid,
  depth_min = 0, depth_max = list(1200, seafloor)
)

# Is this really necessary 
plot_depth_profile(
  species_name = "Carcharhinus longimanus",
  range_rast = sp_range_rast,
  rast_3d = t_annual
)
```

#### Map at depth

``` r

# Map of temperature at 100m depth within the white shark's range
sp_range <- species_ranges[species_ranges$sci_name == "Carcharhinus longimanus", ]

plot_range_at_depth(
  species_range = sp_range,
  depth = 1000,
  rast_3d = t_annual
)
```

``` r

library(ggplot2)
library(gganimate)
library(tidyr)

sp_name <- "Carcharodon carcharias"
sp_range <- species_ranges[species_ranges$sci_name == sp_name, ]
sp_depths <- depth_table_complete[
  depth_table_complete$scientific_name == sp_name, ]

# Per-cell depth window, clamped to bathymetry.
range_rast <- vect_to_envelope(
  polygon = sp_range,
  template = study_grid,
  depth_min = sp_depths$upper_depth,
  depth_max = list(sp_depths$lower_depth, seafloor)
)

# Mask t_annual to the cells + depths the species actually occupies, then
# bounding-box crop so as.data.frame doesn't iterate the entire globe.
masked <- mask(t_annual, range_rast) |>
  terra::crop(terra::vect(sp_range))

long_df <- terra::as.data.frame(masked, xy = TRUE, na.rm = FALSE) |>
  pivot_longer(cols = -c(x, y), names_to = "layer", values_to = "value") |>
  mutate(depth = depths(layer)) |>
  filter(!is.na(value))

# Lock the legend to the global min/max so colour mapping is identical per frame.
v_limits <- range(long_df$value, na.rm = TRUE)

# Factor with numerically-sorted levels so frames play 0 m -> deepest.
long_df$depth <- factor(long_df$depth, levels = sort(unique(long_df$depth)))

p <- ggplot(long_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = sp_range, fill = NA, colour = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  scale_fill_viridis_c(name = NULL, limits = v_limits, na.value = "transparent") +
  coord_sf() +
  labs(subtitle = "Depth: {current_frame} m") +
  theme_minimal() +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
  ) +
  transition_manual(depth)

anim <- animate(
  p,
  nframes = nlevels(long_df$depth),
  fps = 4, width = 700, height = 400,
  renderer = gifski_renderer()
)

anim_save(here("output/range_at_depth_animation.gif"), anim)
anim
```

### Comparison with Original Analysis

| Aspect | Original (`explore_woa_take_2.qmd`) | `ocean3d` |
|----|----|----|
| Fetch depths from IUCN API | ~50 lines of custom API parsing | [`fetch_species_assessments()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/fetch_species_assessments.md) |
| Fill missing depths | ~30 lines of dplyr manipulation | [`fill_missing_depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/fill_missing_depths.md) |
| Load species ranges | ~10 lines with manual filtering | inline [`sf::st_read()`](https://r-spatial.github.io/sf/reference/st_read.html) with SQL filter |
| Download WOA files | manual URL lookup per variable | [`woa_download()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_download.md) (cached) |
| Create monthly summaries | `data-raw/WOA.R` (~170 lines) | `woa_summarise_monthly()` helper (this article) |
| Per-species extraction loop | ~90 lines foreach loop | [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md) + [`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html), wrapped in a `summarise_range()` helper (this article) |
| **Total** | **~350+ lines across multiple files** | **~30 lines of package calls** |
