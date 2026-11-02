# Data model trial: stars, ncdfCF and the current ncdf4 code (#112)

Trial run on 2 November 2026 on a machine with 31 GB of memory. Package
versions, scripts and raw results are in
[`data-model-trial/`](data-model-trial/).

## Decision

**Keep the package's own core object (its description of a gridded dataset,
without the values), read it lazily (only the values each step needs) with
ncdf4, and use CFtime for time, as SPEC.md already plans. Don't make a stars
proxy the core object, and don't depend on ncdfCF yet. Read remote Zarr
stores (a cloud storage format) with `stars::read_mdim()` as an optional
dependency, with the package correcting its coordinates.**

Agreed by Jay Matsushiba on 2 November 2026 and recorded in SPEC.md
("Planned: one model"). "CF" below means the Climate and Forecast
conventions for netCDF metadata; "as doubles" is a file's size once loaded
into R.

- **stars objects that point to a file without reading it (proxies) read
  far more than needed.** Extracting 2 points took 198 s, or failed asking
  for 51 GB of memory. Reading cell by cell took 266–333 s for 500 points
  (package: 11–49 s). stars couldn't hand slices of the real Copernicus file
  to terra for volume. Its second reader, `read_mdim()`, got 17 of 29 cases
  wrong, and its 500-point results differ from every other reader on all
  three datasets.
- **ncdfCF reads only what is asked for and gets 25 of 29 cases right (as
  does stars `read_ncdf()`, which ignores explicit bounds), but its interface
  still changes between releases.** Its values and volumes matched the
  package's, and its speed was similar. It is marked experimental, changed
  its interface in 4 of its last 8 releases, has one maintainer, and has
  quirks the package would have to work around (section 4). It is a useful
  reference for tests; revisit at its 1.0 release or if native model grids
  come into scope.
- **CFtime fixes the time cases.** Every time failure in the package's
  current reader (T03–T12) is decoded correctly by both readers that use
  CFtime (stars `read_ncdf()` and ncdfCF), and `CFtime::indexOf()` matched
  500 dates to a 10-year daily axis in 0.08 s. This supports #111 as planned.
- **Some CF handling stays in the package whichever reader is used.** None
  of the readers tested converts depth units to metres (D02) or turns
  upward elevation into depth (D03): stars drops the information and
  ncdfCF only exposes it. Choosing between files in a series (#94, #96) and
  the nearest-layer and nearest-time rules stay package code with each of
  them.
- **terra's own netCDF reading is wrong without warning on some files** (19
  of 29 cases right): `360_day` dates move by a day, explicit cell bounds
  are ignored, and irregular grids become index positions. So ncdf4 also
  reads the slices that become voxels.
- **The package's buffer extraction is slow:** 76–434 s for one 200 km
  buffer that a single block read handles in 4–7 s. This is the current
  code, not a reader limit.
- **Native model grids are out of scope for now** (decided by Jay
  Matsushiba, 2 November 2026). Inputs are regular lon/lat grids with fixed
  depth levels. NEMO-style curvilinear grids and terrain-following (sigma)
  depths, the main things ncdfCF handles that the package doesn't, are
  regridded before use.

The rest of this note gives the results each point rests on.

## What was compared

| Option | How it reads |
|---|---|
| ocean3d | The package's current code: `extract_to_point()`, `extract_to_buffer()` and their ncdf4 helpers; terra for volume |
| stars `read_ncdf()` | ncmeta and RNetCDF, `proxy = TRUE` |
| stars `read_mdim()` | GDAL's multidimensional API, `proxy = TRUE` |
| ncdfCF | `open_ncdf()`, then `subset()` (lazy) and `raw()` / `terra()` |
| terra (correctness only) | `terra::rast(file)` through GDAL's netCDF driver, with `depth()` and `time()`; how `woa_load_nc()` reads now |

For stars and ncdfCF, the code the package would still write on top (picking
the file for a date, the nearest cell on each axis, the cells inside a
buffer) is in `bench.R`. It uses only what each reader returns.

## 1. Correctness

29 small netCDF files (`cases.R`), each testing one CF feature, with the
expected value at a query cell. For #88–#96 the expected value is the one
written in the issue; otherwise it is the package's current, tested result.
The package (last column) fails in the same way as described in each issue
from #88 to #96 that a case covers, which confirms the cases match the
issues.

"ok": right cell. "wrong": wrong cell or no value. "error": the reader
stopped. stars and ncdfCF are scored on their coordinates as returned, with
no unit conversion or direction flip added.

| Case | From | stars `read_ncdf` | stars `read_mdim` | ncdfCF | terra | ocean3d |
|---|---|---|---|---|---|---|
| T01 days since, steps 0, 2, 5 | fixture | ok | wrong¹ | ok | ok | ok |
| T02 origin with an hour (06:00) | test | ok | wrong¹ | ok | ok | ok |
| T03 `noleap` calendar | test | ok | ok | ok | ok | error |
| T04 `360_day` calendar | test | ok | ok | ok | wrong⁹ | error |
| T05 `months since`, one step (WOA23) | #88 | ok | ok | ok | wrong⁹ | error |
| T06 `... UTC` suffix (GHRSST) | #88 | ok | wrong¹ | ok | ok | error |
| T07 `00:00:0.0` origin (NCEP, ERA5) | #88 | ok | wrong¹ | ok | ok | error |
| T08 unpadded `2020-1-1` (PSL ERSST) | #88 | ok | wrong¹ | ok | ok | error |
| T09 no seconds `00:00` (NCEI ERSST) | #88 | ok | wrong¹ | ok | ok | error |
| T10 year 0001, `proleptic_gregorian` | #88 | ok | ok | ok | ok | error |
| T11 10-minute steps as fractional days | #90 | ok | wrong¹ | ok | wrong⁹ | wrong (32 NA) |
| T12 repeated timestamp | #91 | ok | wrong¹ | ok | ok | error |
| T13 UNLIMITED time dimension | #92 | ok | wrong¹ | ok | ok | ok |
| T14 monthly means with `time_bnds` | test | ok² | wrong¹ | ok | ok | ok |
| D01 depth (m), `positive = "down"` | fixture | ok | ok | wrong³ | ok | ok |
| D02 depth in cm | test | wrong⁴ | wrong⁴ | wrong⁴ | wrong | ok |
| D03 elevation, `positive = "up"` (Copernicus Zarr) | test | wrong⁵ | wrong⁵ | wrong⁵ | wrong | error |
| D04 depth stored deepest first | test | ok | wrong⁶ | wrong³ | ok | ok |
| D05 `altitude` with `axis = "Z"` | #93 | ok | ok | ok | ok | error |
| D06 explicit depth bounds | test | wrong⁷ | wrong⁶ | ok | wrong⁷ | differs⁸ |
| D07 depth with `standard_name`, `axis = "Z"` (Copernicus) | layout | ok | ok | ok | ok | ok |
| H01 0–360 file, query −90 | test | ok | ok | ok | ok | ok |
| H02 −180–180 file, query 270 | test | ok | ok | ok | ok | ok |
| H03 explicit longitude bounds | test | wrong⁷ | wrong⁷ | ok | wrong⁷ | ok |
| H04 axes named `x`/`y`, CF labels swapped | #93 | ok | wrong | ok | wrong | wrong |
| H05 units-only 2D `nav_lon`/`nav_lat` | #93 | ok (error) | wrong | ok (error) | wrong | wrong |
| H06 coarse longitudes 0, 200, 300 | #95 | ok | ok | ok | wrong⁹ | error |
| H07 latitude stored north to south | test | ok | ok | ok | ok | ok |
| V01 `short` with `scale_factor`, `add_offset`, `_FillValue` | layout | ok | ok | ok | ok | ok |
| **Right** | | **25** | **12** | **25** | **19** | **13** |

1. `read_mdim()` decodes time only when the axis has a `calendar` attribute.
   Without one (T01, T02, T06–T09, T11–T14), time is a plain number. CF
   makes `calendar` optional (default `standard`), so files may leave it
   out.
2. Right cell, but `read_ncdf()` ignores `time_bnds`: it builds intervals
   from neighbouring time values ([16 Jan, 15 Feb) instead of [1 Jan, 1 Feb)).
   ncdfCF reads the bounds.
3. ncdfCF classes a `depth` axis that has only units and `positive` as a
   plain numeric axis, not a vertical one, although CF says `positive`
   identifies a vertical coordinate. With `standard_name` or `axis = "Z"`
   (D07, as in Copernicus files) it is recognised.
4. No reader converts the units. stars keeps 5000 [cm], ncdfCF exposes the
   `units` attribute, and neither gives metres.
5. No reader turns elevation into depth. ncdfCF exposes `positive = "up"`;
   stars drops it.
6. `read_mdim()` treats each coordinate as the start of a cell, so centres
   move half a step (5000 cm becomes 7500, 100 m becomes 75 m). With explicit
   bounds it uses their midpoints as cell starts.
7. Explicit bounds variables (`lon_bnds`, `depth_bnds`) are ignored; edges
   are inferred from the centres.
8. `method = "nearest"` picks the nearest layer centre by design and doesn't
   use depth bounds; the case expects the layer whose bounds contain 35 m.
9. terra (GDAL's netCDF driver) converts `360_day` dates to R dates (30
   February becomes 1 March), returns `months since` as a bare number,
   rounds fractional-day steps off (16 of 144 exact matches fail, as in
   #90), and reads irregular coordinates (0, 200, 300) as index positions
   because its grids must be regular. It also can't place a file with a
   single latitude row, so the longitude cases use two rows.

**Multi-file cases (#94, #96).** Of the routes tested, none chooses between
files in a series: `read_ncdf()`, `read_mdim()` and ncdfCF each read one file
(ncdfCF can join files along time afterwards with `append()`). stars'
`read_stars()` can stack files along a dimension, but wasn't tested: it
reads netCDF through the same GDAL driver as terra, and stacking files
doesn't apply the rules #94 and #96 need (choosing between a static file and
a time series, not borrowing a time step from an unrelated tile). Those
rules stay in package code.

**#92** (a caller-opened ncdf4 connection with `readunlim = FALSE`) is
specific to accepting ncdf4 connections; stars and ncdfCF open the file
themselves and read the UNLIMITED axis correctly (T13).

**#89 (speed of date matching).** 500 dates against a 3,650-step daily axis:
`extract2d()` 10.9 s end to end, `CFtime::indexOf()` 0.08 s, both 500 of 500
right (`time_matching.R`).

**terra as a reader.** terra is in the table because it already reads
netCDF for the voxel tools (`woa_load_nc()` calls `terra::rast(file)`). It
fills `terra::depth()` and `terra::time()` from the file and reads regular
Copernicus-style grids correctly (its volumes matched on the real file), but
its errors on the other cases are silent: a `360_day` date moves to the next
day, bounds that change cell edges (and so volumes) are ignored, irregular
grids become index positions. Using it for voxel slices would need its grid
checked against the core object's axes on every read, with ncdf4 as the
fallback: two read paths with two sets of CF rules. So ncdf4 reads all
netCDF files that have depth or time, including voxel slices, and terra works
on the slices once they are in memory.

## 2. Memory and speed

Peak memory is for the whole R process, including about 200–240 MB for R
and the packages (the "load" row). How it was measured is in the
[README](data-model-trial/README.md).

Three datasets, all larger than memory or split across files:

- **Copernicus Marine download (main test):** global daily temperature at
  1/12°, cut to 150°E–150°W, 30°S–30°N, the top 47 m and two years (June
  2022 to May 2024). 16.9 GB on disk, **54.7 GB once loaded**, stored in
  blocks of 23 days.
- **Synthetic:** a 0.25° global file shaped like a Copernicus download,
  stored one day per block, the slowest layout for time series. 8.7 GB on
  disk, **60 GB once loaded**.
- **NOAA OISST v2.1**, daily global 0.25° sea surface temperature, one file
  per year (2022, 2023; 470 MB each, 2D + time), as a multi-file series.

Seconds / peak MB. **Copernicus Marine download:**

| Step | ocean3d | stars `read_ncdf` | stars `read_mdim` | ncdfCF |
|---|---|---|---|---|
| Load packages only | 4.5 / 234 | 3.2 / 215 | 3.9 / 214 | 3.5 / 204 |
| Open without reading values | 4.3 / 243 | 4.1 / 252 | 3.9 / 226 | 3.9 / 233 |
| 500 points (place, depth, day) | 11 / 277 | 333 / 252 | 30 / 256 ⚠ | 47 / 329 |
| same, stars' `st_extract()` on the proxy, 2 points | | 198 / 427 | fails (51.0 GB) | |
| One point, 731 days | 7 / 280 | 7 / 253 | 4 / 236 | 4 / 244 |
| Volume and overlap, one day, via terra | 37 / 1,415 | fails² | fails² | 17 / 2,429 |
| Buffer: 200 km × ±20 m × ±15 days | 76 / 806 | 6 / 288 | 4 / 273 | 4 / 279 |
| `subset()` of 380 million cells (3 GB) | | | | 3.9 / 233 |
| then read 2 days of it | | | | 4.1 / 292 |

2. The slice reads, but stars won't convert it to terra: the file's float32
   coordinates make stars treat the grid as irregular, and
   `as(x, "SpatRaster")` stops with "only regular rasters can be
   converted". The same happens with the Zarr store (section 3), where the
   SpatRaster had to be built from the array by hand.

**Synthetic file:**

| Step | ocean3d | stars `read_ncdf` | stars `read_mdim` | ncdfCF |
|---|---|---|---|---|
| Load packages only | 3.3 / 236 | 2.7 / 214 | 2.7 / 214 | 2.5 / 204 |
| Open without reading values | 3.4 / 244 | 3.6 / 250 | 2.9 / 227 | 2.9 / 230 |
| 500 points (place, depth, day) | 49 / 302 | 277 / 250 | 36 / 258 ⚠ | 33 / 385 |
| same, stars' `st_extract()` on the proxy, 2 points | | 243 / 1,020 | fails (56.4 GB) | |
| One point, 730 days | 11 / 303 | 12 / 296 | 10 / 286 | 11 / 310 |
| Volume and overlap, one day, via terra | 12 / 1,866 | 13 / 2,248 | 11 / 2,094 | 14 / 2,774 |
| Buffer: 200 km × ±50 m × ±15 days | 434 / 419 | 7 / 273 | 5 / 252 | 6 / 276 |
| `subset()` of 757 million cells (6 GB) | | | | 3.4 / 233 |
| then read 2 days of it | | | | 3.3 / 337 |

**OISST, two files.** Points and the series span both files; the
buffer's ±15-day window crosses the year boundary.

| Step | ocean3d | stars `read_ncdf` | stars `read_mdim` | ncdfCF |
|---|---|---|---|---|
| Load packages only | 3.3 / 236 | 3.1 / 215 | 3.2 / 215 | 2.8 / 205 |
| Open both files | 3.5 / 247 | 4.6 / 251 | 3.2 / 227 | 3.2 / 230 |
| 500 points | 11 / 288 | 266 / 255 | 42 / 239 ⚠ | 27 / 330 |
| One point, 730 days | 15 / 339 | 15 / 277 | 12 / 259 | 12 / 248 |
| Buffer across the two files | 92 / 384 | 6 / 272 | 4 / 253 | 4 / 263 |

Results agree between readers except where marked:

- **Points:** `read_ncdf()` and ncdfCF returned the same 500 values on all
  three datasets. The package matched them, except for one extra coastal
  value on the Copernicus file (its nearest-depth rule skips missing
  layers). ⚠ `read_mdim()` differed on every dataset: its half-cell
  coordinate shift (note 6) picks a neighbouring layer or cell for some
  points. The values are in `results/benchmarks.tsv`.
- **Volume:** every reader that completed gives the same volumes
  (Copernicus: 6,438 km³ at 10–20 °C, 360,362 km³ at 15–25 °C; synthetic:
  116,835,900 km³, 111,867,900 km³ and 69,352,360 km³ overlap). Peak memory
  is one day's slice plus the package's terra steps, so it follows the
  slice size, not the file.
- **Series and buffers:** identical series means. The package's buffer
  means differ from a plain average of the cells inside the radius by under
  0.01 °C, because it weights cells by how much of each lies in the buffer.
  With a depth window reaching below the file's deepest level it
  returns `NA` with `insufficient_coverage` (less than 80% of the buffer has
  data), where the plain reads silently average what is there.

What this shows:

- **Every option opens without reading values** and keeps peak memory near
  the size of what each step reads (under 400 MB for points and series on
  files of 55–60 GB as doubles).
- **stars' own point extraction isn't lazy enough to use** (198–243 s for 2
  points; allocations of 51–56 GB on `read_mdim()`). Reading cell by cell
  through `read_ncdf(ncsub = )` works but costs about 0.5 s per call.
- **stars slices of real Copernicus data can't go to terra** without
  rebuilding the raster by hand (note 2). Volume and overlap, which run on
  terra, are where this matters.
- **The package's `extract_to_buffer()` is slow:** 76 s (Copernicus), 92 s
  (OISST) and 434 s (synthetic) for one buffer that plain block reads handle
  in 4–7 s. This is the current code, not a reader limit; the block reads
  planned for #115 address it.

## 3. Remote Zarr (optional)

The real product, `cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m` (global
1/12°, 50 levels, daily since June 2022: 4320 × 2041 × 50 × 1625 cells, 2.9
TB as float32), read with `read_mdim()` from Copernicus Marine's public ARCO
stores, without a login (`zarr.R`).

| Step | Seconds / peak MB | Result |
|---|---|---|
| Open (`timeChunked` store) | 14 / 216 | dimensions and decoded time read; no values |
| Surface series, one point, 1,625 days (`geoChunked`) | 89 / 220 | 1,625 values |
| One day's global surface map (`timeChunked`) | 79 / 532 | 8.8 million values |
| Volume and overlap, one day, 50 levels, 10° × 10° | 310 / 480 | 199,227 km³ at 10–20 °C, 242,608 km³ at 15–25 °C, 85,121 km³ overlap¹ |
| ncdf4 (netCDF-C NCZarr) on the same store | fails | `nc_open` can't open it, so ncdfCF and the package's ncdf4 code have no Zarr route here |

1. The process exited with status 1 after printing these volumes; the
   cause wasn't investigated. The SpatRaster was built from the array by
   hand: the store's float32 coordinates make stars treat the grid as
   irregular, and `as(x, "SpatRaster")` then refuses it ("only regular
   rasters can be converted").

- Reads stay lazy: memory follows the slice read, not the 2.9 TB store.
- The vertical axis is `elevation` with `positive = "up"` (−5728 to −0.494
  m), returned as is, and time comes back as `Date`. The package would flip
  elevation to depth (D03) and keep sub-daily time itself.
- Speed depends on the store layout: a time series is read from
  `geoChunked` (134 days per chunk) and maps from `timeChunked` (one day per
  chunk). The core object would need to know both stores.

## 4. ncdfCF: lazy `subset()` and interface stability

- **Lazy:** yes, for regular lon/lat grids. `subset()` returns a
  file-backed `CFVariable` with start and count and reads nothing until
  `raw()`, `terra()` or similar; peak memory was 233 MB for a 6 GB subset.
  On grids with 2D (auxiliary) lon/lat it reads straight away, because it
  interpolates.
- **Quirks the package would work around:**
  - `subset()` selects the cells whose *coordinates* fall in a range. A
    single value, or a range between two cell centres, returns `NULL`, so a
    point needs a range around its cell's coordinate.
  - Time ranges are half-open even with `rightmost.closed = TRUE`:
    `time = c("2022-01-01", "2023-12-31")` leaves out 31 December.
  - `terra()` returns a `SpatRasterDataset` without a CRS, so terra
    computes cell areas (and volumes) on a plane until the CRS is set.
  - `profile()` reads the whole depth × time column at a point, which on
    files chunked by time step means decompressing every chunk.
- **Stability:** 11 releases from June 2024 to March 2026, with interface changes
  in 0.4 (`CFData` renamed `CFArray`), 0.6 (`subset()` signature), 0.7
  (`CFArray` merged into `CFVariable`) and 0.8 (indexing must include
  length-1 axes, `CFResource` renamed `NCResource`). One maintainer (who
  also maintains CFtime); built on RNetCDF, not ncdf4. On CRAN on 2
  November 2026, only stars listed it, as an optional dependency.
