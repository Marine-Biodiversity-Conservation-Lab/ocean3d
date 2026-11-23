# ocean3d (development version)

## Breaking changes

`calc_volume_overlap()` is removed (#150). It returned a nine-layer
`SpatRaster` that the other 3D functions could only use as a flat 2D layer,
not as the 3D space it described.
Intersect with `intersect_3d()`, then measure with `volume()`, which has a
new `per_cell` argument: `volume(x, per_cell = TRUE)` returns a one-layer
`SpatRaster` named `volume` with the volume of each cell in km³, `NA` where
`x` is absent. Each layer of the old stack maps to:

| `calc_volume_overlap(a, b)` layer | Now |
|---|---|
| `depth_min_a`, `depth_max_a` | `a` itself (envelope), or `voxel_to_envelope(a)` (voxel) |
| `depth_min_overlap`, `depth_max_overlap` | `ov <- intersect_3d(a, b)` |
| `volume_a`, `volume_b` | `volume(a, per_cell = TRUE)`, `volume(b, per_cell = TRUE)` |
| `volume_overlap` | `volume(ov, per_cell = TRUE)` |
| `0` where both are present but share no depth | `both <- !is.na(volume(a, per_cell = TRUE)) & !is.na(volume(b, per_cell = TRUE))`, then `terra::ifel(both & is.na(vol_ov), 0, vol_ov)` with `vol_ov <- volume(ov, per_cell = TRUE)` |
| totals | `volume(a)`, `volume(b)`, `volume(ov)` |

Other changes to callers:

* With an envelope and a voxel, `calc_volume_overlap()` reported the volume of
  the envelope after placing it on the voxel's depth levels. `volume()` of the
  envelope uses its own depth limits. To compare it with the overlap on the
  same levels, use
  `volume(envelope_to_voxel(a, depths(b), bounds = bounds, depth_bounds = depth_bounds(b)), bounds = bounds)`.
* Cells with nothing in them are now `NA` everywhere. The old
  `volume_overlap` layer was the one place that used `0`, for cells where
  both domains are present but share no depth; build that `0` from the
  inputs as in the table.
* `volume()` of a domain with no cells now totals `0`, not `NA`.
* `per_cell` is the second argument of `volume()`, so a positional
  `volume(v, "top")` is now an error. Name it: `volume(v, bounds = "top")`.
* Pass the same `bounds` to `intersect_3d()` and `volume()`; otherwise the
  overlap and the input volumes are measured on different slabs and don't add
  up.

A `SpatVoxel` now stores each layer's depth in terra's own per-layer field,
`terra::depth()`, instead of in its layer names (#113). Before, the depth was
read from names like `temp_depth=100`, so renaming a layer broke the voxel.
Now the package never reads depths from layer names:

* Renaming a layer doesn't change its depth. `as_voxel()` still names layers
  `{varname}_depth={value}` when it sets the depths, as labels.
* `as_voxel()` takes the depths from `depths`, or else from the raster's
  `terra::depth()`, which terra also sets when it reads a netCDF file. A
  raster with neither is an error, even if its layers are named
  `{variable}_depth={value}`. Before, such names were read as depths.
* `depths()` reads `terra::depth()`, and is an error for a raster without
  depths. It no longer takes a character vector of layer names. To label a
  per-layer summary such as `terra::global(v, ...)`, use `depths(v)`: the
  summary has one row per layer, in layer order.
* Selecting layers (`v[[2:3]]`, `v$name`, `subset()`) and replacing one with
  `[[<-` keep each layer's depth through later operations. Without the
  `SpatVoxel` methods for these, terra would reset the depths to 0, for
  example in `v[[2:3]] * 2` or a crop of `v[[3]]`.
* `c()` of voxels keeps each layer's depth. When the layers don't make a
  valid voxel, for example two layers at the same depth, it returns a plain
  `SpatRaster`; before, it returned an invalid `SpatVoxel`.
* A few terra operations, such as `app()`, drop the depths and return a
  plain `SpatRaster`; the package's functions then stop with an error. Set
  the depths again with `as_voxel(x, depths = ...)`.

Code that set a voxel's depth through its layer names, or read depths back
from names, should pass `depths` to `as_voxel()` or set `terra::depth()`,
and read them with `depths()`. `terra::writeCDF()` and
`saveRDS(terra::wrap(v))` keep the depths. GeoTIFF doesn't, so a voxel read
back from GeoTIFF needs `depths`. All three return a plain `SpatRaster`, so
re-wrap with `as_voxel()` after reading.

`woa_load_nc()` is removed (#113). It loaded one field of a World Ocean Atlas
file as a plain `SpatRaster`. Read a WOA file like any other netCDF file:
describe the field with `as_proxy()`, naming its variable (`{code}_{field}`,
such as `t_an`; leave out `var` to list the file's variables), crop it to
your area with `terra::crop()`, and read it with `as_voxel()`:

```r
temp <- as_proxy(woa_download("temperature"), var = "t_an")
v <- as_voxel(terra::crop(temp, study_area))
```

Only the cropped cells are read from the file. The voxel carries WOA's own
depth bounds, its CRS (EPSG:4326) and no date, since a WOA file is a
climatology for one period. Because it has depth bounds, the `bounds`
argument of `volume()`, `mask()` and the other 3D functions doesn't apply
to it. WOA's bounds are halfway between its standard levels, so results with
the default `bounds = "midpoint"` don't change; `bounds = "top"` no longer
does anything. `woa_download()`'s help page lists the variable and field
codes, and the new article "Working with World Ocean Atlas data" shows how
to read each variable and statistical field, monthly and seasonal files
included.

`bounds` now defaults to `"midpoint"` instead of `"top"` in `volume()`,
`envelope_to_voxel()`, `intersect_3d()`,
`intersects_3d()` and `mask()` (#138). The slab of water each depth level
stands for now runs halfway to the neighbouring levels. This matches World
Ocean Atlas levels and is closer than `"top"` to Copernicus Marine layers.
Under the old default each level was the top of its slab and the deepest
level had no volume. Results change for voxels that don't carry their own
depth edges (see `depth_bounds` under New features):

* `volume()` gives different volumes. For example,
  levels at 0, 100 and 300 m now stand for 0-50, 50-200 and 200-300 m rather
  than 0-100, 100-300 and none.
* `envelope_to_voxel()` can place an envelope on different levels, and so
  can `mask()`, `intersect_3d()` and `intersects_3d()` when they combine an
  envelope with a voxel. An envelope from 60 to 70 m on levels at 0, 100 and
  200 m now occupies the 100 m level rather than the 0 m one.

Pass `bounds = "top"` to keep the old results.

Point extraction (`extract_to_point()`, `extract2d()`, `extract3d_*()`) now
checks the grid's coordinates instead of assuming one layout (#74). Some
calls that used to return a value now return `NA` or stop with an error:

* By default (`spatial_match = "within"`), a point must lie inside a grid
  cell. Points outside the grid, which used to get the value of the nearest
  edge cell, now return `NA`. This includes observations near a
  single-station file that are not at the station itself. Use
  `spatial_match = "nearest"` with `max_distance_km` to accept nearby cells.
* Observation depths must be metres, positive downward. Negative depths
  (heights) now stop with an error; convert them, and km or cm depths,
  before extraction.
* Depth coordinates are converted to metres from their units, and the
  `positive` attribute sets their direction. Surface, bottom and nearest are
  chosen by physical depth, whether the coordinate is stored shallow to
  deep or deep to shallow. The depth coordinate must have length units,
  including for `extract3d_surface()` and `extract3d_bottom()`; files
  without them stop with an error. A file with a single depth level can be
  read with `extract2d()`.
* With several files, a temporal tie between sources keeps input order under
  `"within"` matching, and goes to the closer cell centre under `"nearest"`
  matching.

Point, track and buffer extraction now check cell bounds the same way, and
stop with an error when they are bad (#113):

* Each cell must have a positive width, contain its centre and not overlap
  another cell, and a `bounds` attribute must name a usable bounds
  variable. This applies to longitude and latitude bounds everywhere, and
  to depth and time bounds in buffer extraction. Horizontal buffers also
  reject longitude cells that overlap after wrapping around 360°. The error
  names the file.
* Point and track extraction used to warn and infer longitude and latitude
  edges from the centres. They don't read time bounds.
* To read a file that now fails, fix its bounds, or remove the coordinate's
  `bounds` attribute so edges are inferred from the centres (for example
  with NCO: `ncatted -a bounds,lon,d,, file.nc`).

Point and track extraction, and buffer extraction without a depth window
(`depth_radius_m = 0`), used to accept a netCDF depth coordinate stored out
of order (for example 50, 0, 100 m), and rejected only duplicate depths.
They now stop with an error naming the coordinate and the file unless the
depths run shallow to deep or deep to shallow (#147). Sort such files along
the depth dimension before extraction. The error for longitude or latitude
stored out of order, which already stopped extraction, now names the
coordinate and the file too.

With `return_match_info = TRUE`, the `bounds` entry of the
`<output>_match_info` column from `extract_to_buffer()` used to label
longitude and latitude edges inferred from the centres
`"inferred_regular"` or `"inferred_irregular"`. Both now read
`"inferred_midpoints"`, the label depth and time already use, so the entry
no longer says whether the spacing is regular (#147). Code that checks for
the old labels should check for `"inferred_midpoints"`.

`vect_to_envelope()` now stops with an error when the envelope would have no
cells: no cell the polygons cover has `depth_max` deeper than `depth_min`.
It used to warn and return an empty envelope, so code that loops over
polygons and drops empty results has to handle the error instead. The
message says how many cells have the limits swapped (or a seafloor above
`depth_min`), equal, or missing (#114).

`envelope_to_voxel()` and `depth_bounds()` now stop with an error when given
a `SpatEnvelope` or `SpatVoxel` whose layers no longer match its class, for
example after a terra operation changed them (#97). They used to use such an
object without checking it. `envelope_to_voxel()` could then return a wrong
voxel, such as an empty one for an envelope whose two layers had been
swapped. `occupied()`, `voxel_to_envelope()`, `volume()`, `intersect_3d()`,
`intersects_3d()` and `mask()` already rejected such objects. Rebuild the
object with `as_envelope()` or `as_voxel()` first.

`envelope_to_voxel()` now stops with an error when `depths` repeats a
depth, as `as_voxel()` does (#97). It used to drop the repeats silently.
Pass each depth once, for example with `unique(depths)`.

`as_voxel()` on a `SpatProxy` now reads the step nearest to a date or
date-time in `time`, as `extract_to_point()` does (#115). It used to require
an exact match and stop otherwise. For the old behaviour, pass
`time_match = "exact"`; `max_time_diff` limits how far a nearest match may
be, in days.

## New features

* `as_voxel()` on a `SpatProxy` reads several time steps at once (#115):
  `time` can be a window of two dates or date-times, giving a list of
  voxels, one per step inside it, or a `SpatVector3D` with one feature,
  whose time (a point) or time interval (a polygon) selects the steps. A
  shape without a time reads the only step of a one-step source, such as a
  World Ocean Atlas climatology.

* New arguments `spatial_match`, `max_distance_km`, `spatial_report` and
  `depth_report`. Longitudes match across the 0-360 and -180-180
  conventions and across the antimeridian.
* Missing outputs are explained once, under the report for their cause
  (spatial, depth or time), with repeated reasons grouped by file. The time
  report now also explains outputs that are missing because no observation
  date was supplied for a file with a time dimension.
* New experimental class `SpatProxy`, built with `as_proxy()`, describes one
  variable of gridded data without reading its values (#113): its axes,
  chunk sizes and coordinate reference system (CRS). Depth is in metres
  positive down with each level's edges, and time is a CFtime object. The
  data can be:
  * one netCDF file;
  * a series of netCDF files, such as one per day, put in time order and
    joined into one time axis. The files must hold the same variable on the
    same grid and depths, and their time steps must not repeat or overlap;
  * a `SpatVoxel` already in memory, as one time step. Its time must be a
    date or a date-time.

  A climatology with one time step, such as a World Ocean Atlas file, is
  read without a date. `as_proxy()` stops with an error if a time axis
  isn't strictly increasing or is a climatology with more than one step.
* `as_voxel(proxy, time = ...)` reads one time step of a `SpatProxy` into a
  `SpatVoxel` (#113), chosen by step number or date. The voxel's extent
  comes from the cell edges, its CRS is always set, its depths and depth
  edges come from the proxy, and its date is in `terra::time()`. For a file
  without depth bounds, the depth edges are halfway between levels, so the
  voxel's `bounds` argument doesn't apply. Grids that are regular apart from
  rounding in the stored coordinates are treated as regular; other irregular
  grids stop with an error.
* `terra::crop()` of a `SpatProxy` keeps the cells that overlap an area (an
  extent, raster, vector or sf object), still without reading values, so
  `as_voxel()` then reads only those cells from the file (#113). An area in
  another CRS is projected first, and an extent in the other longitude
  convention (0 to 360 or -180 to 180) finds its cells. An area that crosses
  the edge of the grid's longitudes stops; crop each side separately.

  `SpatProxy`, `as_proxy()`, cropping a `SpatProxy` and reading one with
  `as_voxel()` are experimental and may change without notice while the core
  object is built.
* New class `SpatVector3D` for 3D points, lines and polygons, built with
  `as_vector3d()` from a table, sf object or `terra::SpatVector`. Points
  and lines have a depth (positive metres down) and an optional time at
  each vertex (for sf lines, depth is minus Z and time is M). Polygons have
  a depth range, `depth_min` to `depth_max`, and an optional time interval
  per feature, such as a protected area from 0 to 200 m. `geomtype()`
  reports the type (#114).
* New `buffer_3d()` turns each `SpatVector3D` point or line into a polygon
  with a depth range and an optional time interval, such as a circle of
  10 km radius from 10 to 60 m deep around a detection at 35 m. `volume()`
  on points or lines now suggests it (#114).
* New `vector3d_to_envelope()` puts `SpatVector3D` polygons onto a grid as
  a `SpatEnvelope`, as `vect_to_envelope()` does for 2D polygons. Where
  polygons overlap, a cell gets the shallowest top and deepest bottom of
  the polygons covering it, so polygons from 0 to 50 m and from 200 to
  300 m give 0 to 300 m. Cells where top equals bottom are dropped, and so
  are polygon times, with a message (#114).
* `volume()` measures `SpatVector3D` polygons: the sum of each polygon's
  area times its depth range, in km³. Space where polygons overlap is
  counted once for each polygon, so the total isn't the volume of their
  union (#115).
* `intersect_3d()` and `intersects_3d()` accept `SpatVector3D` polygons
  with a `SpatEnvelope` or `SpatVoxel`, in either order, and `mask()`
  accepts them as the mask. The polygons are first put on the raster's
  grid, so a depth gap between overlapping polygons over one cell is filled
  in. `mask(voxel, polygons)` keeps the voxel's values in the cells and
  depth levels the polygons reach. Polygons that miss the grid give an
  empty result, not an error (#115).
* `intersect_3d()` of two sets of `SpatVector3D` polygons gives a polygon
  for each pair whose outlines, depth ranges and time intervals overlap,
  with columns `feature_id_x` and `feature_id_y` naming the pair. Measure
  the overlap with `volume(intersect_3d(x, y))` (#115).
* `intersects_3d()` of two sets of `SpatVector3D` polygons gives, for each
  polygon of `x`, the feature IDs of the polygons of `y` it shares 3D space
  with, as a list named by the IDs of `x` (#115).
* A `SpatVoxel` can carry the top and bottom of the slab of water each layer
  stands for, such as a netCDF file's depth bounds: pass `depth_bounds` to
  `as_voxel()` or `envelope_to_voxel()`, and read them with the new
  `depth_bounds()`. `volume()`, `calc_volume_overlap()`,
  `envelope_to_voxel()`, `intersect_3d()`, `intersects_3d()` and `mask()`
  then use these edges instead of `bounds`. Results that stay a `SpatVoxel`
  keep the edges, including after `c()`; `app()` and writing to a file drop
  them. For such a voxel, `voxel_to_envelope()` gives each cell an envelope
  from the top of its shallowest selected slab to the bottom of its deepest,
  and the depth limits from `calc_volume_overlap()` follow the same rule
  (#138).

## Bug fixes

* `as_proxy()`, `extract_to_point()`, `extract_to_buffer()` and
  `extract_to_track()` now open remote sources given as a URL, such as
  OPeNDAP or NCZarr addresses. They used to stop with "netCDF file does
  not exist" for any path that was not a local file (#160).
* netCDF time axes are decoded with the CFtime package. Time units that
  used to stop extraction are now read: a trailing `UTC` or `Z`
  (`days since 2020-01-01 UTC`), origins written `2020-1-1`,
  `00:00:0.0` or without seconds, years before 1000 with the
  `proleptic_gregorian` calendar, and a single time step in
  `months since` units, as in World Ocean Atlas files (#88). A time axis
  of several months or years still stops with an error, because their
  length is ambiguous.
* Reference dates before 1582-10-15 on the `standard` and `gregorian`
  calendars are accepted, such as `hours since 1-1-1 00:00:0.0` in
  older NCEP reanalysis files; CFtime counts the Julian days before the
  calendar change. Time steps before 1582-10-15 on these calendars still
  stop with an error.
* Time steps stored as fractional days (for example 10-minute data in
  `days since 1970-01-01`) decode exactly, so exact time matching no
  longer returns `NA` for observations on a file timestamp (#90).
* A file passed as a connection opened with
  `ncdf4::nc_open(readunlim = FALSE)` now gives the same result as its
  path when the time dimension is unlimited (#92).
* Files on calendars other than `standard`, `gregorian` and
  `proleptic_gregorian` (such as `noleap` or `360_day`) still stop, now with
  an error that names the calendar.
* Matching observation dates to long time axes is fast again: 500 dates
  against a 10-year daily file take about 0.35 s, down from about 11 s
  (#89).
* A time step that appears more than once (for example after joining two
  files that overlap) no longer stops the whole extraction. It is an error
  only when an observation selects it, and the error names the file and
  the timestamp (#91).
* An observation exactly halfway between two time steps takes the one on
  its own UTC day. Midnight date-times, as `readxl::read_excel()` returns
  for Excel dates, used to get the previous day on files stamped at noon
  (#71).
* Files whose time bounds are back to back (one interval ends where the
  next starts, as in monthly means) no longer print CFtime's "Offsets not
  monotonically increasing" warning in point and buffer extraction. Time
  bounds that don't contain their time step still stop buffer extraction
  with an error (#113).
* The `axis` and `standard_name` attributes of ordinary dimension
  coordinates are now read when longitude, latitude, depth and time are
  detected in point, buffer and track extraction. Files whose dimension
  names don't say which axis is which (such as `x` labelled latitude) used
  to be read along the wrong axes without a warning, and an `altitude`
  dimension with `axis = "Z"` stopped with an "Unrecognised non-singleton
  dimension" error unless `depth_dim` was given (#93).
* Curvilinear grids, such as native NEMO output with two-dimensional
  `nav_lon` / `nav_lat` coordinates, now stop with an error when those
  coordinates are identified only by their units (`degrees_east`,
  `degrees_north`). They used to be read with the grid indices as longitude
  and latitude. Regrid such data to a regular longitude/latitude grid
  first (#93).

## Minor improvements

* `vect_to_envelope()` now projects polygons in another CRS onto the CRS of
  `template`, as the 3D query functions do with polygons (#97). It used to
  stop with an error asking for the polygons to be projected first. Depth
  rasters in `depth_min` and `depth_max` must still be on the template's
  grid, CRS included.

* Functions that take a `SpatEnvelope` check it faster when it is held in
  memory. Envelopes read from a file, or with more than a million cells, are
  checked as before. On 9-cell test grids, `as_envelope()` is about 5 times
  faster, `intersect_3d()` about 2 times, `volume()` about 1.6 times and
  `calc_volume_overlap()` about 1.2 times (#154).

* `as_envelope()`'s depth errors now come from the `SpatEnvelope` validity
  check, so they start with `invalid class "SpatEnvelope" object:` and use
  its wording: "depth_max must be at least depth_min in every cell", and for
  negative depths "depths are positive metres increasing downward, but
  `depth_min` has negative values; if these are elevations, flip their sign
  first" (#154).

## Dependencies

* CFtime (>= 1.7.3) is now imported, for decoding CF time units.

# ocean3d 0.1.1.9009

Housekeeping release addressing `goodpractice` checks (#58). Apart from the
changes below, behaviour is unchanged.

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
