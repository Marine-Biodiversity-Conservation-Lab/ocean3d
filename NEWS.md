# ocean3d 0.1.1.9009

Housekeeping release addressing `goodpractice` checks (#58). Apart from the
changes below, behaviour is unchanged.

## Breaking changes

* `terra::intersect()` is no longer intercepted for `SpatEnvelope` and
  `SpatVoxel`. It previously raised an error pointing to `intersects_3d()` and
  `intersect_3d()`; it now falls through to terra's 2D method, which ignores
  depth. Use `intersects_3d()` / `intersect_3d()` for 3D questions.

## Bug fixes

* Point extraction functions (`extract_to_point()`, `extract2d()`,
  `extract3d_*()`) no longer drop netCDF sources when `nc` is a list that
  contains a data frame of file paths; previously the data frame replaced
  every source listed before it. A vector of paths inside a list is now
  accepted too.

## Minor improvements

* `dplyr` and `gfwr` move from `Imports` to `Suggests`; no package code uses
  them, so installing ocean3d no longer pulls them in.

# ocean3d 0.1.1.9008

The `copernicus_*` family is withdrawn for now. It is not abandoned: the code
is preserved on the `feat/7-copernicus` branch and in the history of this one,
and it will return once it is properly covered by tests (#50).

## Breaking changes

* `copernicus_setup()`, `copernicus_login()`, `copernicus_status()`,
  `copernicus_load()`, `copernicus_summarise()`, `copernicus_cache_dir()` and
  `copernicus_cache_clear()` are removed. Code calling them should pin to
  0.1.1.9007 or use the `feat/7-copernicus` branch.
* `ecmwfr` and `jsonlite` are dropped from `Imports`, and `reticulate` and
  `rstudioapi` from `Suggests`; they were only ever used by the Copernicus
  backends. `here` moves out of `Imports` but stays in `Suggests`, where the
  vignette articles still use it.

# ocean3d 0.1.1.9007

The package is renamed from `sharkabc3d` to `ocean3d`, to reflect that the 3D
workflows apply to any marine taxon or spatial layer, not just sharks and rays.

## Breaking changes

* `library(sharkabc3d)` becomes `library(ocean3d)`. No function names,
  arguments or class names changed — only the package name.
* The WOA download cache moves with the package name, from
  `tools::R_user_dir("sharkabc3d", "cache")/woa` to
  `tools::R_user_dir("ocean3d", "cache")/woa`. Previously downloaded WOA files
  are not found at the new location; either re-download them or move the old
  `woa` directory across. `woa_cache_dir()` reports the new path.
* Package `Title` and `Description` are broadened from "Shark and Ray Abiotic
  Covariates in 3 Dimensions" to "Three-Dimensional Marine Spatial Analysis".
* The repository and pkgdown site move to
  `Marine-Biodiversity-Conservation-Lab/ocean3d`; GitHub redirects the old
  URLs, but `devtools::install_github()` calls should be updated.

# sharkabc3d 0.1.1.9006

Copernicus environmental data acquisition. A new `copernicus_*` family
downloads, organises and summarises data from Copernicus Marine and the
ECMWF-operated Copernicus Data Stores.

## New functions

* `copernicus_load()` acquires data from Copernicus Marine (`source =
  "marine"`) or the CDS/ADS stores (`"cds"`, `"ads"`), filtered by variable,
  time range, bounding box and depth range, and returns the paths it wrote.
  Marine requests go through either the standalone Copernicus Marine Toolbox
  or the Python `copernicusmarine` package; `backend = "auto"` prefers the
  standalone executable and falls back to Python. Output is always one
  variable per file, optionally subdivided with `split_by` (year, season,
  month, ISO week, day, hour) and filed into a directory hierarchy with
  `organize_by`. Seasons are meteorological, with `hemisphere` choosing only
  the human-readable names. CDS and ADS requests are served by `ecmwfr`.
* `copernicus_summarise()` collapses one or more Copernicus netCDF files
  across time into `"mean"`, `"min"`, `"max"` and `"sd"` files, one per
  statistic, working a time step at a time so large downloads do not have to
  fit in memory. Spatial and depth dimensions are kept.
* `copernicus_setup()` installs the versioned standalone Toolbox into the
  user's application data directory — no Python, conda or `reticulate`
  needed. `copernicus_login()` configures credentials once (in the RStudio
  Terminal when run from RStudio, since external password prompts are
  unreliable in the console), and `copernicus_status()` reports what is
  installed and whether credentials are valid.
* `copernicus_cache_dir()` and `copernicus_cache_clear()` manage the
  persistent download cache under `tools::R_user_dir()`.

## Dependencies

* New `Imports`: `ecmwfr` and `jsonlite`, plus `dplyr`, `gfwr (>= 3.0)` and
  `here`, which were already used by the GFW and IUCN code but undeclared.
  `reticulate` and `rstudioapi` are in `Suggests`, needed only for the Python
  backend and the RStudio login path respectively.

# sharkabc3d 0.1.1.9005

The 3D object model. `SpatVoxel` / `SpatEnvelope` replace the ad-hoc "3D
raster" conventions, and the volumetric functions now dispatch on them.

## New classes and constructors

* `SpatVoxel` (one layer per standard depth, cell value is the variable) and
  `SpatEnvelope` (exactly `depth_min`/`depth_max`, cell value is depth): S4
  classes extending `terra::SpatRaster` with their own validity rules.
  `SpatVolume` is a class union over the two and the dispatch target for
  volumetric operations, replacing the virtual `SpatDepthRaster`.
* `as_voxel()` / `as_envelope()` construct them. Prefer these over
  `methods::new()`: they normalise as well as validate, and are idempotent, so
  re-wrapping after a terra operation that drops the class is cheap.
  `as_voxel()` accepts a conforming multi-depth `SpatRaster`, a list of
  single-depth rasters, or a `SpatVoxel`, and sorts layers shallow to deep.
  `as_envelope()` attaches depth limits to a 2D footprint: each limit is a
  single number or a per-cell single-layer `SpatRaster`, or `x` already
  carries the two depth layers.
* `voxel_to_envelope()` / `envelope_to_voxel()` convert between the two, each
  taking a `fun` — an occupancy predicate one way (default `\(x) !is.na(x)`),
  a per-voxel profile the other. `voxel_to_envelope()` is lossy: an envelope
  is one continuous interval per cell, so interior gaps are filled.
* `vect_to_envelope()` builds an envelope from polygons. `depth_min` /
  `depth_max` are lists mixing numerics and rasters; per cell the deepest
  `depth_min` and shallowest `depth_max` win, so bathymetry is just another
  constraint. Constraint rasters must match the template's CRS, resolution and
  extent, and their NAs propagate rather than falling back to the remaining
  constraints.
* `depths()`, exported: the depths a voxel's layers stand for, parsed from
  `{variable}_depth={value}` layer names. It accepts a `SpatVoxel`, any
  conforming `SpatRaster`, or a bare character vector, so
  `depths(rownames(terra::global(v, "sum")))` works.
* Depths are positive metres increasing downward throughout the 3D model.
  Negative depths are an error, not silently flipped.

## New spatial query verbs

Each takes a `SpatEnvelope` or `SpatVoxel` paired with another 3D object, a 2D
`SpatRaster` footprint, or polygons (`SpatVector`, `sf`, `sfc`). A 2D object
restricts the domain horizontally and leaves its depths alone.

* `intersect_3d(x, y)` returns the 3D space the two share: an envelope with the
  shared depth interval, or a presence voxel. It is the one place that decides
  where two domains overlap; `calc_volume_overlap()` is now built on it, with
  unchanged numbers.
* `intersects_3d(x, y)` is the yes/no form: `TRUE` where both are present and
  their depths overlap, `FALSE` where both are present but disjoint or only one
  is, `NA` where neither is. terra's own three-way answer, made depth-aware.
* `mask()` gains depth-aware methods on terra's generic, so the article idiom
  `terra::mask(v, envelope_to_voxel(range, depths(v)))` becomes
  `mask(v, range)`. A voxel masked by an envelope puts the envelope on the
  voxel's depth levels first; a voxel masked by a voxel requires the same
  levels. Mismatched levels used to mask by layer position, silently; that is
  now an error.
* `terra::intersect()` on a `SpatEnvelope` or `SpatVoxel` is now an error
  naming `intersects_3d()` / `intersect_3d()`. terra's version is a 2D "both
  have data" test, and because terra keeps the subclass on its result, two
  envelopes 100 m apart came back as a `SpatEnvelope` of `TRUE` that passed
  validity.
* The verbs re-check tagged inputs and point to `as_voxel()` / `as_envelope()`,
  since terra keeps the class tag through operations that change the layer set
  (`[[`, `app()`).

## Breaking changes

| before | after |
| --- | --- |
| `calc_volume(r)` | `volume(as_envelope(r))` |
| `calc_volume_overlap(range_rast_a = a, range_rast_b = b)` | `calc_volume_overlap(x = a, y = b)` |
| `count_3d_overlap(a, b)` | `intersects_3d(a, b)` |
| `extract_rast_range(range, v)` | `mask(v, range)` |
| `extract_rast_volume(area, min_depth = 40, max_depth = 600, rast_3d = v)` | `extract_to_area(area, v, min_depth = 40, max_depth = 600)` |
| `voxelize_range(poly, voxel, depth_min = 0, depth_max = 500)` | `vect_to_envelope(poly, template, depth_min = 0, depth_max = list(500, seafloor))` |
| `create_study_voxel(template, bathymetry, depths)` | no replacement object — pass the three separately (the `study_voxel` class and its `print()` method go too) |
| `summarise_species_environment()` | `summarise_range()` in `vignette("woa-environmental-extraction")` |
| `voxel_to_envelope(fun = "extent" \| "threshold")` | `voxel_to_envelope(fun = <predicate>)` |

* `volume()` (was `calc_volume()`, matching how terra names `area()`) and
  `calc_volume_overlap()` are S4 generics on the 3D classes, so they accept a
  `SpatVoxel` as well as a `SpatEnvelope`, with formals `x` / `y`. The class is
  now the contract: a bare `SpatRaster` carrying `depth_min`/`depth_max` layers
  is rejected rather than duck-typed. A voxel's volume sums the slab each
  occupied level stands for, so interior gaps cost volume rather than being
  filled in; given one envelope and one voxel, the envelope is discretized onto
  the voxel's levels. Both now check that their inputs share a grid instead of
  reusing the first argument's cell areas for the second. Envelope volumes are
  numerically unchanged.
* `count_3d_overlap()` is retired for `intersects_3d()`, which returns
  `TRUE`/`FALSE`/`NA` instead of `1`/`NA`. `sum(stack, na.rm = TRUE)` over a
  stack of results still gives a richness map, now with a real `0` where a
  domain is present but nothing overlaps it.
* `extract_to_area()` (was `extract_rast_volume()`) takes the target first and
  the source second, as `extract_to_point()` does, and its depth bounds are
  optional and named. It requires and returns a `SpatVoxel`, accepts `sfc`
  areas alongside `sf` and `SpatVector`, and compares CRS with
  `terra::same.crs()` rather than string equality on WKT.
* `vect_to_envelope()` requires `polygon` and `template` to already share a
  CRS rather than projecting silently.
* `summarise_species_environment()` was analysis rather than package
  machinery. Its column contract (`{name}_min`, `_max`, `_mean`,
  `_n_surface_cells`, `_n_cells`, `_n_depths`) is unchanged in the vignette
  helper.
* `R/extract.R` is retired — it re-derived by hand what the classes now
  guarantee. `methods` moves into Imports.

## Bug fixes

* `mask(v, range)` also fixes a depth bug in `extract_rast_range()`, which
  tested exact point containment, so an envelope falling between two depth
  levels occupied none of them: `[10, 20]` against depths `c(0, 100, 200)` came
  back empty. `envelope_to_voxel()` tests slab overlap instead, so the same
  envelope now occupies the 0 m level, and its `bounds` argument chooses where
  the slab edges fall (`"top"` or the World Ocean Atlas `"midpoint"`).
* Touching is not overlapping, for voxels as well as envelopes.
  `envelope_to_voxel()` used to record an envelope at a level whose slab it
  merely ended on, so `[0, 100]` and `[100, 200]` both occupied the 100 m
  level and their voxels intersected, while `intersect_3d()` on the same two
  envelopes said they did not. The occupancy test is now strict on the
  envelope's floor: a level is occupied only where the envelope shares some
  thickness of water with its slab. Adjacent depth ranges therefore land on
  disjoint levels, and `intersect_3d()`, `intersects_3d()`, `mask()` and
  `calc_volume_overlap()` give the envelope answer whichever types they are
  handed. Two consequences: an envelope ending exactly at a level no longer
  occupies that level, and under `"top"` an envelope ending at `max(depths)`
  is not recorded at the deepest level, which stands for a zero-thickness
  slab.
* `calc_volume_overlap()` returned its nine-layer stack tagged `SpatEnvelope`.
  It is now the plain `SpatRaster` its documentation describes.

## Documentation and internals

* New article `vignette("woa-species-range-voxels")`: a species range polygon
  to a `SpatEnvelope`, applied to a WOA-shaped `SpatVoxel`, through to the
  overlapping voxel values and their summaries. Unlike the other articles it is
  fully synthetic and runs, so it doubles as a worked check of the idioms.
* Helpers shared by `R/volume.R` and the new `R/intersect.R` move to
  `R/utils-3d.R`. The constructors work on plain copies of their input, so a
  voxel layer can serve as a footprint without its class tag reaching the
  depth-aware `mask()` methods.

# sharkabc3d 0.1.1.9004

* `gfw_effort_to_raster()` can infer grid resolution and extent from the input
  data frame.

# sharkabc3d 0.1.1.9003

* Breaking: `load_bathymetry()` is renamed `load_gebco_bathymetry()` and moves
  from `R/load_data.R` to `R/gebco_bathymetry.R`. The new name makes the
  expected data source explicit; the function is otherwise unchanged.
* Runnable examples for `load_gebco_bathymetry()` — building a small stand-in
  NetCDF rather than requiring a multi-gigabyte GEBCO download (#37) — and for
  `fill_missing_depths()`, which moves to `R/iucn_utils.R` alongside the
  `fetch_species_assessments()` depth limits it repairs.

# sharkabc3d 0.1.1.9002

* `extract_to_point()` generalises point extraction from netCDF files, taking
  data frames, tibbles, `sf` POINT objects, matrices, named lists, and bare
  longitude/latitude/depth/date values or vectors. It preserves input structure
  where appropriate and transforms projected `sf` coordinates to
  longitude/latitude when possible.
* Runnable synthetic 2D and 3D netCDF examples for `extract_to_point()`,
  `extract2d()`, `extract3d_surface()`, `extract3d_bottom()`,
  `extract3d_nearest()` and `extract3d_all()`.

# sharkabc3d 0.1.1.9001

* Retire `R/plot.R`. The `plot_*()` functions were used by a single article and
  are no longer exported; `plot_depth_profile()` and `plot_range_at_depth()`
  now live inline in the WOA environmental extraction articles. `ggplot2` and
  `tidyterra` move from Imports to Suggests.

# sharkabc3d 0.1.1.9000

* Deprecate `woa_nc_extract()`, folding it into `woa_load_nc()`, its only
  caller. `woa_load_nc()` now checks for layer names that don't follow the
  package convention.

# sharkabc3d 0.1.1

* Add NetCDF extraction to observation points (#6, @davidruizgarci)
