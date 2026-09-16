
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ocean3d

<!-- badges: start -->
[![pkgcheck](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/workflows/pkgcheck/badge.svg)](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/actions?query=workflow%3Apkgcheck)
[![R-CMD-check](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/Marine-Biodiversity-Conservation-Lab/ocean3d/graph/badge.svg)](https://app.codecov.io/gh/Marine-Biodiversity-Conservation-Lab/ocean3d)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

**ocean3d** is an R package for three-dimensional marine spatial
analysis of species ranges, fisheries, and depth-stratified
oceanographic data. It grew out of work on sharks, rays, and chimaeras,
with the functions generalising to any marine taxon or spatial layer.

The three-dimensional marine environment poses a unique challenge for
spatial analyses. Most conventional GIS workflows represent space as a
two-dimensional plane. This abstraction fails to capture the range of
depths and vertical overlap of marine habitat. Fisheries also operate at
different depths depending on gear type and target species, creating
depth-specific patterns of threat exposure that require analysis in 3D
volume rather than 2D planes to quantify. `ocean3d` provides
documented workflows and reusable functions for:

- extracting values from depth-stratified oceanographic rasters within
  user-defined geographic areas and depth windows;
- characterising the abiotic habitat (e.g., temperature, dissolved
  oxygen) of marine species using IUCN Red List ranges and World Ocean
  Atlas 2023 climatologies;
- calculating three-dimensional overlap between species distributions,
  fisheries effort, and other spatial layers;
- producing depth-aware maps and depth-profile plots of species ranges,
  environmental conditions, and cumulative fishing pressure.

## How it works

All 3D analyses use a stacked-raster approach built on `terra` and `sf`.
Polygons (species ranges, fishery footprints) are rasterized onto a
common bathymetry-aware grid where each cell stores presence plus the
shallowest and deepest depths the feature occupies (clamped to the
seafloor). Volume overlap between two rasterized ranges is then computed
per-cell via raster algebra. Multi-depth environmental rasters (e.g.,
WOA temperature at 57 standard depths) follow a
`{variable}_depth={value}` layer-naming convention so that downstream
functions can select the correct layers for a given depth window.

## Installation

You can install the development version of ocean3d from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Marine-Biodiversity-Conservation-Lab/ocean3d")
```

Add an IUCN Red List API key to access species-assessment functions.
Sign up and acquire your API key at <https://api.iucnredlist.org/>:

``` r
usethis::edit_r_environ()
# Add the following line to the .Renviron file
# IUCN_REDLIST_KEY="your_iucn_api_key_here"
```

## Function overview

### Loading and preparing input data

- `load_gebco_bathymetry()` — load a GEBCO bathymetry NetCDF as a
  `SpatRaster`.
- `fetch_species_assessments()` — query the IUCN Red List API for
  taxonomy, Red List category, and depth limits, by SIS ID, scientific
  name, or comprehensive group code.
- `fill_missing_depths()` — fix swapped upper/lower depth values and
  fill missing values from genus-level means.

#### World Ocean Atlas 2023 utilities

- `woa_download()` — download WOA 2023 NetCDF files (temperature,
  salinity, dissolved oxygen, oxygen saturation, AOU, nitrate,
  phosphate, silicate, density) at 0.25°, 1°, or 5° resolution, with
  caching.
- `woa_load_nc()` — load a WOA NetCDF and select a statistical field
  (e.g., objectively analyzed climatology), returning a `SpatRaster`
  with the package’s standard `{variable}_depth={value}` layer names.
- `woa_cache_dir()` / `woa_cache_clear()` — manage the persistent WOA
  download cache.

#### Copernicus (Marine, CDS, ADS)

- `copernicus_setup()` — download and configure the standalone Toolbox
  executable (no separate Python installation needed).
- `copernicus_login()` — configure or verify authentication, once.
- `copernicus_status()` — report whether the Toolbox is available, and
  from which path and version.
- `copernicus_load()` — download a dataset from Copernicus Marine, the
  Climate Data Store, or the Atmosphere Data Store, subset by variable,
  time window, bounding box and depth range.
- `copernicus_summarise()` — reduce downloaded NetCDFs across time
  (mean, min, max, sd), preserving depth and the other non-temporal
  dimensions, working on the NetCDF directly via `ncdf4`.
- `copernicus_cache_dir()` / `copernicus_cache_clear()` — manage the
  persistent Copernicus download cache.

#### Global Fishing Watch (fisheries effort)

- `gfw_effort_to_raster()` — turn the long-format apparent-fishing-hours
  tibble from `gfwr::gfw_ais_fishing_hours()` into a multi-layer
  `SpatRaster`, one layer per gear (or other grouping).
  Extending that into 3D uses the general converters — `as_envelope()` and
  `envelope_to_voxel()` — driven by a gear-to-depth-band lookup. Those
  operating-depth priors are analysis assumptions rather than package data,
  so they live in the `gfw-fishing-effort-3d` article.

### Building the 3D study grid and rasterized ranges

- `create_study_raster()` — build an empty study-area `SpatRaster`
  covering the combined extent of one or more spatial inputs.

### Two 3D representations, and moving between them

A `SpatEnvelope` stores one continuous `[depth_min, depth_max]` interval
per cell; a `SpatVoxel` stores a value at each of a set of depth levels,
so it can hold a variable and interior gaps.

- `as_voxel()` — wrap a multi-depth `SpatRaster` (or a list of
  single-depth ones) as a validated `SpatVoxel`, ordered shallow to deep
  with `{variable}_depth={value}` layer names. 
- `as_envelope()` — define a 3D envelope with two rasters, with `depth_min` 
  and `depth_max` corresponding to the minimum and maximum depths present in 
  given cell. 
- `vect_to_envelope()` — rasterize a single species range or fishery
  footprint onto the study grid and attach per-cell `depth_min` /
  `depth_max`, returning a `SpatEnvelope`. 
- `envelope_to_voxel()` — expand an envelope onto a set of depth levels.
- `voxel_to_envelope()` — collapse a voxel to the envelope bounding it,
  using a predicate on the cell values. Lossy: interior gaps are filled.

### 3D spatial query

Three verbs ask how a `SpatEnvelope` or `SpatVoxel` relates to another
object in 3D. The other object can be a 3D object too, or a 2D one: a
`SpatRaster` footprint or polygons (`SpatVector`, `sf`). A 2D object
restricts the domain horizontally and leaves its depths alone.

- `intersect_3d(x, y)` — the 3D space the two share. Returns an
  envelope with the shared depth interval, or a presence voxel.
- `intersects_3d(x, y)` — whether they share any 3D space, cell by cell:
  `TRUE`, `FALSE` (both present but disjoint, or only one present), or
  `NA` (neither present). Sum a stack of these for a richness map.
- `mask(x, mask)` — `terra::mask()`, made depth-aware. Keeps the values
  of `x` where the mask domain reaches that cell *at that depth*.

`terra::intersect()` on a 3D object is an error that points to these,
because terra’s version ignores depth.

### 3D volume

Both dispatch on the 3D representation, so they take either a
`SpatEnvelope` or a `SpatVoxel`. A voxel’s volume sums the slab each
occupied depth level stands for, so interior gaps cost volume instead of
being filled in. Given one of each, the envelope is discretized onto the
voxel’s depth levels.

- `volume()` — total 3D volume (km³) of a rasterized domain.
- `calc_volume_overlap()` — per-cell depth intervals and volumes for two
  rasterized domains and their intersection (returns a 9-layer stack).
  Built on `intersect_3d()`.

### Environmental extraction to areas

- `depths()` — the depths a `SpatVoxel`’s layers stand for, parsed from
  the `{variable}_depth={value}` layer names. Also accepts a bare
  character vector of layer names.
- `extract_to_area()` — crop a `SpatVoxel` to an area polygon and select
  layers within a depth range.

To restrict a `SpatVoxel` to a species’ *per-cell* depth window —
preserving each cell’s vertical refuge — mask it with the range envelope:

``` r
mask(rast_3d, range_env)
```

See `vignette("woa-species-range-voxels")` for the full workflow, from a
range polygon through to summary statistics.

### Environmental extraction to point observations

- `extract_to_point()` — the general engine: match observation points
  (data frame, tibble, `sf` POINT, matrix, or named list) against one or
  more NetCDF files, snapping to the nearest longitude, latitude, time
  and depth cell. Returns the input with the extracted values appended,
  structure preserved.
- `extract2d()` — for variables with no depth dimension.
- `extract3d_nearest()` — the valid, non-missing depth layer closest to
  each observation’s depth, so a cell below the seabed does not silently
  return `NA`.
- `extract3d_surface()` / `extract3d_bottom()` — the first available
  depth layer, or the deepest non-missing one at that cell.
- `extract3d_all()` — nearest, surface and bottom in one pass, as three
  columns.



## Contributing

Contributions are welcome. `ocean3d` is being actively built out from a set
of lab analyses into a general, tested package, and there is plenty of
well-scoped work available. Start with:

- **[CONTRIBUTING.md](CONTRIBUTING.md)** — development setup, data requirements,
  the branching/coordination workflow, and coding conventions.
- **[SPEC.md](SPEC.md)** — the roadmap; every planned function with its priority
  and source. Unchecked items are up for grabs.
- **[Code of Conduct](CODE_OF_CONDUCT.md)**.

New contributors: the "Known rough edges / good first issues" section of
`CONTRIBUTING.md` lists concrete starting points.

## Acknowledgements

Thank you to the people that have inspired and collaborated on this work! 

David Ruiz-García, Rachel Aitchison, Wade VanderWright, Amanda Arnold, Dr. Samm Sherman, Dr. Alifa Haque.

## Citation

Jay H. Matsushiba, David Ruiz-García, Rachel M. Aitchison, Wade J. VanderWright, C. Samantha Sherman, Alifa Bintha Haque, Nicholas K. Dulvy. *ocean3d: An R Package for Three-Dimensional Spatial Analyses of Marine Habitat.*
