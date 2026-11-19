# Changelog

## ocean3d (development version)

### Breaking changes

`bounds` now defaults to `"midpoint"` instead of `"top"` in
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
[`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md),
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
and [`mask()`](https://rspatial.github.io/terra/reference/mask.html)
([\#138](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/138)).
The slab of water each depth level stands for now runs halfway to the
neighbouring levels. This matches World Ocean Atlas levels and is closer
than `"top"` to Copernicus Marine layers. Under the old default each
level was the top of its slab and the deepest level had no volume.
Results change for voxels that don’t carry their own depth edges (see
`depth_bounds` under New features):

- [`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
  and
  [`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md)
  give different volumes. For example, levels at 0, 100 and 300 m now
  stand for 0-50, 50-200 and 200-300 m rather than 0-100, 100-300 and
  none.
- [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
  can place an envelope on different levels, and so can
  [`mask()`](https://rspatial.github.io/terra/reference/mask.html),
  [`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
  and
  [`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
  when they combine an envelope with a voxel. An envelope from 60 to 70
  m on levels at 0, 100 and 200 m now occupies the 100 m level rather
  than the 0 m one.

Pass `bounds = "top"` to keep the old results.

Point extraction
([`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
[`extract2d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract2d.md),
`extract3d_*()`) now checks the grid’s coordinates instead of assuming
one layout
([\#74](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/74)).
Some calls that used to return a value now return `NA` or stop with an
error:

- By default (`spatial_match = "within"`), a point must lie inside a
  grid cell. Points outside the grid, which used to get the value of the
  nearest edge cell, now return `NA`. This includes observations near a
  single-station file that are not at the station itself. Use
  `spatial_match = "nearest"` with `max_distance_km` to accept nearby
  cells.
- Observation depths must be metres, positive downward. Negative depths
  (heights) now stop with an error; convert them, and km or cm depths,
  before extraction.
- Depth coordinates are converted to metres from their units, and the
  `positive` attribute sets their direction. Surface, bottom and nearest
  are chosen by physical depth, whether the coordinate is stored shallow
  to deep or deep to shallow. The depth coordinate must have length
  units, including for
  [`extract3d_surface()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_surface.md)
  and
  [`extract3d_bottom()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_bottom.md);
  files without them stop with an error. A file with a single depth
  level can be read with
  [`extract2d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract2d.md).
- With several files, a temporal tie between sources keeps input order
  under `"within"` matching, and goes to the closer cell centre under
  `"nearest"` matching.

Point, track and buffer extraction now check cell bounds the same way,
and stop with an error when they are bad
([\#113](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/113)):

- Each cell must have a positive width, contain its centre and not
  overlap another cell, and a `bounds` attribute must name a usable
  bounds variable. This applies to longitude and latitude bounds
  everywhere, and to depth and time bounds in buffer extraction.
  Horizontal buffers also reject longitude cells that overlap after
  wrapping around 360°. The error names the file.
- Point and track extraction used to warn and infer longitude and
  latitude edges from the centres. They don’t read time bounds.
- To read a file that now fails, fix its bounds, or remove the
  coordinate’s `bounds` attribute so edges are inferred from the centres
  (for example with NCO: `ncatted -a bounds,lon,d,, file.nc`).

Point and track extraction, and buffer extraction without a depth window
(`depth_radius_m = 0`), used to accept a netCDF depth coordinate stored
out of order (for example 50, 0, 100 m), and rejected only duplicate
depths. They now stop with an error naming the coordinate and the file
unless the depths run shallow to deep or deep to shallow
([\#147](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/147)).
Sort such files along the depth dimension before extraction. The error
for longitude or latitude stored out of order, which already stopped
extraction, now names the coordinate and the file too.

With `return_match_info = TRUE`, the `bounds` entry of the
`<output>_match_info` column from
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
used to label longitude and latitude edges inferred from the centres
`"inferred_regular"` or `"inferred_irregular"`. Both now read
`"inferred_midpoints"`, the label depth and time already use, so the
entry no longer says whether the spacing is regular
([\#147](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/147)).
Code that checks for the old labels should check for
`"inferred_midpoints"`.

[`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
now stops with an error when the envelope would have no cells: no cell
the polygons cover has `depth_max` deeper than `depth_min`. It used to
warn and return an empty envelope, so code that loops over polygons and
drops empty results has to handle the error instead. The message says
how many cells have the limits swapped (or a seafloor above
`depth_min`), equal, or missing
([\#114](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/114)).

### New features

- New arguments `spatial_match`, `max_distance_km`, `spatial_report` and
  `depth_report`. Longitudes match across the 0-360 and -180-180
  conventions and across the antimeridian.
- Missing outputs are explained once, under the report for their cause
  (spatial, depth or time), with repeated reasons grouped by file. The
  time report now also explains outputs that are missing because no
  observation date was supplied for a file with a time dimension.
- New experimental class `SpatProxy`, built with
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md),
  describes one variable of a netCDF file without reading its values
  ([\#113](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/113)):
  its axes, chunk sizes and coordinate reference system (CRS). Depth is
  in metres positive down and time is a CFtime object.
  [`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
  stops with an error if the time axis isn’t strictly increasing or is a
  climatology. No other function accepts a `SpatProxy` yet, and it may
  change without notice.
- New class `SpatVector3D` for 3D points, lines and polygons, built with
  [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md)
  from a table, sf object or
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html).
  Points and lines have a depth (positive metres down) and an optional
  time at each vertex (for sf lines, depth is minus Z and time is M).
  Polygons have a depth range, `depth_min` to `depth_max`, and an
  optional time interval per feature, such as a protected area from 0 to
  200 m.
  [`geomtype()`](https://rspatial.github.io/terra/reference/geomtype.html)
  reports the type
  ([\#114](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/114)).
- New
  [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)
  turns each `SpatVector3D` point or line into a polygon with a depth
  range and an optional time interval, such as a circle of 10 km radius
  from 10 to 60 m deep around a detection at 35 m.
  [`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
  on points or lines now suggests it
  ([\#114](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/114)).
- New
  [`vector3d_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vector3d_to_envelope.md)
  puts `SpatVector3D` polygons onto a grid as a `SpatEnvelope`, as
  [`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
  does for 2D polygons. Where polygons overlap, a cell gets the
  shallowest top and deepest bottom of the polygons covering it, so
  polygons from 0 to 50 m and from 200 to 300 m give 0 to 300 m. Cells
  where top equals bottom are dropped, and so are polygon times, with a
  message
  ([\#114](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/114)).
- A `SpatVoxel` can carry the top and bottom of the slab of water each
  layer stands for, such as a netCDF file’s depth bounds: pass
  `depth_bounds` to
  [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  or
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
  and read them with the new
  [`depth_bounds()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depth_bounds.md).
  [`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
  [`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md),
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md),
  [`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
  [`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
  and [`mask()`](https://rspatial.github.io/terra/reference/mask.html)
  then use these edges instead of `bounds`. Results that stay a
  `SpatVoxel` keep the edges, including after
  [`c()`](https://rdrr.io/r/base/c.html);
  [`app()`](https://rspatial.github.io/terra/reference/app.html) and
  writing to a file drop them. For such a voxel,
  [`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
  gives each cell an envelope from the top of its shallowest selected
  slab to the bottom of its deepest, and the depth limits from
  [`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md)
  follow the same rule
  ([\#138](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/138)).

### Bug fixes

- netCDF time axes are decoded with the CFtime package. Time units that
  used to stop extraction are now read: a trailing `UTC` or `Z`
  (`days since 2020-01-01 UTC`), origins written `2020-1-1`, `00:00:0.0`
  or without seconds, years before 1000 with the `proleptic_gregorian`
  calendar, and a single time step in `months since` units, as in World
  Ocean Atlas files
  ([\#88](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/88)).
  A time axis of several months or years still stops with an error,
  because their length is ambiguous.
- Reference dates before 1582-10-15 on the `standard` and `gregorian`
  calendars are accepted, such as `hours since 1-1-1 00:00:0.0` in older
  NCEP reanalysis files; CFtime counts the Julian days before the
  calendar change. Time steps before 1582-10-15 on these calendars still
  stop with an error.
- Time steps stored as fractional days (for example 10-minute data in
  `days since 1970-01-01`) decode exactly, so exact time matching no
  longer returns `NA` for observations on a file timestamp
  ([\#90](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/90)).
- A file passed as a connection opened with
  `ncdf4::nc_open(readunlim = FALSE)` now gives the same result as its
  path when the time dimension is unlimited
  ([\#92](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/92)).
- Files on calendars other than `standard`, `gregorian` and
  `proleptic_gregorian` (such as `noleap` or `360_day`) still stop, now
  with an error that names the calendar.
- Matching observation dates to long time axes is fast again: 500 dates
  against a 10-year daily file take about 0.35 s, down from about 11 s
  ([\#89](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/89)).
- A time step that appears more than once (for example after joining two
  files that overlap) no longer stops the whole extraction. It is an
  error only when an observation selects it, and the error names the
  file and the timestamp
  ([\#91](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/91)).
- An observation exactly halfway between two time steps takes the one on
  its own UTC day. Midnight date-times, as `readxl::read_excel()`
  returns for Excel dates, used to get the previous day on files stamped
  at noon
  ([\#71](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/71)).
- Files whose time bounds are back to back (one interval ends where the
  next starts, as in monthly means) no longer print CFtime’s “Offsets
  not monotonically increasing” warning in point and buffer extraction.
  Time bounds that don’t contain their time step still stop buffer
  extraction with an error
  ([\#113](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/113)).
- The `axis` and `standard_name` attributes of ordinary dimension
  coordinates are now read when longitude, latitude, depth and time are
  detected in point, buffer and track extraction. Files whose dimension
  names don’t say which axis is which (such as `x` labelled latitude)
  used to be read along the wrong axes without a warning, and an
  `altitude` dimension with `axis = "Z"` stopped with an “Unrecognised
  non-singleton dimension” error unless `depth_dim` was given
  ([\#93](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/93)).
- Curvilinear grids, such as native NEMO output with two-dimensional
  `nav_lon` / `nav_lat` coordinates, now stop with an error when those
  coordinates are identified only by their units (`degrees_east`,
  `degrees_north`). They used to be read with the grid indices as
  longitude and latitude. Regrid such data to a regular
  longitude/latitude grid first
  ([\#93](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/93)).

### Dependencies

- CFtime (\>= 1.7.3) is now imported, for decoding CF time units.

## ocean3d 0.1.1.9009

Housekeeping release addressing `goodpractice` checks
([\#58](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/58)).
Apart from the changes below, behaviour is unchanged.

### Bug fixes

- Point extraction functions
  ([`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
  [`extract2d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract2d.md),
  `extract3d_*()`) no longer drop netCDF sources when `nc` is a list
  that contains a data frame of file paths; previously the data frame
  replaced every source listed before it. A vector of paths inside a
  list is now accepted too.

### Minor improvements

- `dplyr` and `gfwr` move from `Imports` to `Suggests`; no package code
  uses them, so installing ocean3d no longer pulls them in.

## ocean3d 0.1.1.9008

The `copernicus_*` family is withdrawn for now. It is not abandoned: the
code is preserved on the `feat/7-copernicus` branch and in the history
of this one, and it will return once it is properly covered by tests
([\#50](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues/50)).

### Breaking changes

- `copernicus_setup()`, `copernicus_login()`, `copernicus_status()`,
  `copernicus_load()`, `copernicus_summarise()`,
  `copernicus_cache_dir()` and `copernicus_cache_clear()` are removed.
  Code calling them should pin to 0.1.1.9007 or use the
  `feat/7-copernicus` branch.
- `ecmwfr` and `jsonlite` are dropped from `Imports`, and `reticulate`
  and `rstudioapi` from `Suggests`; they were only ever used by the
  Copernicus backends. `here` moves out of `Imports` but stays in
  `Suggests`, where the vignette articles still use it.

## ocean3d 0.1.1.9007

The package is renamed from `sharkabc3d` to `ocean3d`, to reflect that
the 3D workflows apply to any marine taxon or spatial layer, not just
sharks and rays.

### Breaking changes

- [`library(sharkabc3d)`](https://rdrr.io/r/base/library.html) becomes
  [`library(ocean3d)`](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d).
  No function names, arguments or class names changed — only the package
  name.
- The WOA download cache moves with the package name, from
  `tools::R_user_dir("sharkabc3d", "cache")/woa` to
  `tools::R_user_dir("ocean3d", "cache")/woa`. Previously downloaded WOA
  files are not found at the new location; either re-download them or
  move the old `woa` directory across.
  [`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md)
  reports the new path.
- Package `Title` and `Description` are broadened from “Shark and Ray
  Abiotic Covariates in 3 Dimensions” to “Three-Dimensional Marine
  Spatial Analysis”.
- The repository and pkgdown site move to
  `Marine-Biodiversity-Conservation-Lab/ocean3d`; GitHub redirects the
  old URLs, but `devtools::install_github()` calls should be updated.
