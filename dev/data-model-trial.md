# Data model trial: stars, ncdfCF and the current ncdf4 code (#112)

Trial run on 2 November 2026 with stars 0.7.3, ncdfCF 0.8.2, CFtime 1.7.3,
ncdf4 1.24, terra 1.9.50 and GDAL 3.8.4, on R 4.5.3 (Linux, 8 cores, 31 GB
RAM). Scripts and raw results are in [`data-model-trial/`](data-model-trial/).

## Recommendation

**Keep the package's own core object with lazy ncdf4 reads, and use CFtime
for time, as SPEC.md already plans. Don't make a `stars_proxy` the core
object, and don't depend on ncdfCF yet. Read remote Zarr stores with
`stars::read_mdim()` as an optional source (Suggests), with the package
correcting its coordinates.**

- **stars proxies don't stay lazy for the operations the package needs.**
  `st_extract()` on a proxy took 243 s and 1 GB for 2 points on the large
  file (`read_ncdf()`), and tried to allocate 56.4 GB on a `read_mdim()`
  proxy. Reads through `read_ncdf(ncsub = )` return the right values but
  keep the full time axis in the dimensions, and cost about 0.5 s per call
  (266–277 s for 500 points against 11–49 s for the package). `read_mdim()`
  gets 17 of the 29 CF cases wrong, and its errors show up on real data: the
  500-point means differ from every other reader on both large files.
- **ncdfCF reads lazily and gets 25 of the 29 CF cases right (as does stars
  `read_ncdf()`, which ignores explicit bounds), but its API is still
  moving.** `subset()` doesn't read values: subsetting 757 million cells
  (6 GB) used no more memory than opening the file. For 500 points it took
  33 s on the large file (package 49 s) and 27 s on OISST (package 11 s).
  But it is marked experimental,
  has broken its API in 4 of its last 8 releases, has one maintainer, and has
  quirks the package would have to work around (below). It is a good
  reference to test against, and the dependency question can be revisited at
  its 1.0 release.
- **CFtime fixes the time cases.** Every time failure in the package's
  current reader (T03–T12) is decoded correctly by both readers that use
  CFtime (stars `read_ncdf()` and ncdfCF), and `CFtime::indexOf()` matched
  500 dates to a 10-year daily axis in 0.08 s. This supports #111 as planned.
- **Some CF handling stays in the package whichever reader is used.** No
  reader converts depth units to metres (D02) or flips `positive = "up"`
  axes to depth (D03): stars drops the information and ncdfCF only exposes
  it. Choosing between files in a series (#94, #96) and the nearest-layer
  and nearest-time rules are package logic in every option.

The rest of this note gives the results each point rests on.

## What was compared

| Option | How it reads |
|---|---|
| ocean3d | The package's current code: `extract_to_point()`, `extract_to_buffer()` and their ncdf4 helpers; terra for volume |
| stars `read_ncdf()` | ncmeta and RNetCDF, `proxy = TRUE` |
| stars `read_mdim()` | GDAL's multidimensional API, `proxy = TRUE` |
| ncdfCF | `open_ncdf()`, then `subset()` (lazy) and `raw()` / `terra()` |

For stars and ncdfCF, the code the package would still write on top (picking
the file for a date, the nearest cell on each axis, the cells inside a
buffer) is in `bench.R`. It uses only what each reader returns.

## 1. Correctness

29 small netCDF files (`cases.R`), each testing one CF feature, with the
expected value at a query cell. For #88–#96 the expected value is the one
written in the issue; otherwise it is the package's current, tested result.
The package (last column) reproduces every failure described in #88–#96,
which confirms the cases match the issues.

"ok": right cell. "wrong": wrong cell or no value. "error": the reader
stopped. stars and ncdfCF are scored on their coordinates as returned, with
no unit conversion or direction flip added.

| Case | From | stars `read_ncdf` | stars `read_mdim` | ncdfCF | ocean3d |
|---|---|---|---|---|---|
| T01 days since, steps 0, 2, 5 | fixture | ok | wrong¹ | ok | ok |
| T02 origin with an hour (06:00) | test | ok | wrong¹ | ok | ok |
| T03 `noleap` calendar | test | ok | ok | ok | error |
| T04 `360_day` calendar | test | ok | ok | ok | error |
| T05 `months since`, one step (WOA23) | #88 | ok | ok | ok | error |
| T06 `... UTC` suffix (GHRSST) | #88 | ok | wrong¹ | ok | error |
| T07 `00:00:0.0` origin (NCEP, ERA5) | #88 | ok | wrong¹ | ok | error |
| T08 unpadded `2020-1-1` (PSL ERSST) | #88 | ok | wrong¹ | ok | error |
| T09 no seconds `00:00` (NCEI ERSST) | #88 | ok | wrong¹ | ok | error |
| T10 year 0001, `proleptic_gregorian` | #88 | ok | ok | ok | error |
| T11 10-minute steps as fractional days | #90 | ok | wrong¹ | ok | wrong (32 NA) |
| T12 repeated timestamp | #91 | ok | wrong¹ | ok | error |
| T13 UNLIMITED time dimension | #92 | ok | wrong¹ | ok | ok |
| T14 monthly means with `time_bnds` | test | ok² | wrong¹ | ok | ok |
| D01 depth (m), `positive = "down"` | fixture | ok | ok | wrong³ | ok |
| D02 depth in cm | test | wrong⁴ | wrong⁴ | wrong⁴ | ok |
| D03 elevation, `positive = "up"` (Copernicus Zarr) | test | wrong⁵ | wrong⁵ | wrong⁵ | error |
| D04 depth stored deepest first | test | ok | wrong⁶ | wrong³ | ok |
| D05 `altitude` with `axis = "Z"` | #93 | ok | ok | ok | error |
| D06 explicit depth bounds | test | wrong⁷ | wrong⁶ | ok | differs⁸ |
| D07 depth with `standard_name`, `axis = "Z"` (Copernicus) | layout | ok | ok | ok | ok |
| H01 0–360 file, query −90 | test | ok | ok | ok | ok |
| H02 −180–180 file, query 270 | test | ok | ok | ok | ok |
| H03 explicit longitude bounds | test | wrong⁷ | wrong⁷ | ok | ok |
| H04 axes named `x`/`y`, CF labels swapped | #93 | ok | wrong | ok | wrong |
| H05 units-only 2D `nav_lon`/`nav_lat` | #93 | ok (error) | wrong | ok (error) | wrong |
| H06 coarse longitudes 0, 200, 300 | #95 | ok | ok | ok | error |
| H07 latitude stored north to south | test | ok | ok | ok | ok |
| V01 `short` with `scale_factor`, `add_offset`, `_FillValue` | layout | ok | ok | ok | ok |
| **Right** | | **25** | **12** | **25** | **13** |

1. `read_mdim()` decodes time only when the axis has a `calendar` attribute.
   Without one (T01, T02, T06–T09, T11–T14), time is a plain number. CF
   makes `calendar` optional (default `standard`), and many files leave it
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

**Multi-file cases (#94, #96).** Neither stars nor ncdfCF chooses between
files in a series: `read_mdim()` reads one file, and `read_ncdf()` and ncdfCF
read one file at a time (ncdfCF can `append()` along time after reading). The
static-file and cadence rules from #94 and #96 stay in package code with any
reader.

**#92** (a caller-opened ncdf4 connection with `readunlim = FALSE`) is
specific to accepting ncdf4 connections; stars and ncdfCF open the file
themselves and read the UNLIMITED axis correctly (T13).

**#89 (speed of date matching).** 500 dates against a 3,650-step daily axis:
`extract2d()` 10.9 s end to end, `CFtime::indexOf()` 0.08 s, both 500 of 500
right (`time_matching.R`).

## 2. Memory and speed

Each step ran in a new R process, so "peak" is the whole process's peak
resident memory (GNU `time`), including R and the loaded packages (the
"load" row). Each process had a 12 GB address-space limit.

**Large file.** Synthetic, shaped like a Copernicus Marine download
(`make_large_file.R`): 1440 × 720 at 0.25°, 10 depths, 730 daily steps,
`short` with `scale_factor`, `_FillValue` on land, deflate, one chunk per
depth and day. 8.7 GB on disk, 15 GB packed and **60 GB as doubles**, about
twice the machine's RAM. No real multi-year 3D daily file could be
downloaded without Copernicus credentials, so the real product was read as
Zarr instead (section 3).

Seconds / peak MB:

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

**Real files, multi-file series.** NOAA OISST v2.1, daily global 0.25° sea
surface temperature, one file per year (2022, 2023; 470 MB each, 2D + time).
Points and the series span both files; the buffer's ±15-day window crosses
the year boundary.

| Step | ocean3d | stars `read_ncdf` | stars `read_mdim` | ncdfCF |
|---|---|---|---|---|
| Load packages only | 3.3 / 236 | 3.1 / 215 | 3.2 / 215 | 2.8 / 205 |
| Open both files | 3.5 / 247 | 4.6 / 251 | 3.2 / 227 | 3.2 / 230 |
| 500 points | 11 / 288 | 266 / 255 | 42 / 239 ⚠ | 27 / 330 |
| One point, 730 days | 15 / 339 | 15 / 277 | 12 / 259 | 12 / 248 |
| Buffer across the two files | 92 / 384 | 6 / 272 | 4 / 253 | 4 / 263 |

Results agree between readers except where marked:

- **Points:** ocean3d, `read_ncdf()` and ncdfCF return the same 500 values
  (mean 15.8731 on the synthetic file; 18.4863 with 115 land `NA`s on
  OISST). ⚠ `read_mdim()` gives 16.2477 and 18.4551: its half-cell
  coordinate shift (note 6) picks a neighbouring depth layer or cell for
  some points.
- **Volume:** all four give the same volumes (116,835,900 km³ at 10–20 °C,
  111,867,900 km³ at 15–25 °C, 69,352,360 km³ overlap). Peak memory is one
  day's 10.4 million cells plus the package's terra steps, so it follows
  the slice size, not the file.
- **Series and buffers:** identical series means. The package's buffer
  (25.50317 synthetic, 24.902 OISST) differs from the plain equal-weight
  mean of the cells inside the radius (25.50388, 24.89858) because it
  weights cells by their overlap with the buffer.

What this shows:

- **Every option opens without reading values** and keeps peak memory near
  the size of what each step reads (under 400 MB for points, series and
  buffers on a 60 GB file).
- **stars' own point extraction isn't lazy enough to use** (243 s for 2
  points; a 56.4 GB allocation on `read_mdim()`). Reading cell by cell
  through `read_ncdf(ncsub = )` works but costs about 0.5 s per call.
- **The package's `extract_to_buffer()` is slow:** 434 s for one buffer
  that the other readers' plain reads handle in 5–7 s. This is the current
  code, not a reader limit, and is worth its own issue.
- **Converting a 3D stars object to terra** needs `as(x, "SpatRaster")`;
  `terra::rast()` fails on a depth dimension ("1:dd[4]: NA/NaN argument").

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
- The path must be written `ZARR:"/vsicurl/<url>"` with
  `normalize_path = FALSE`; a plain `/vsicurl/` path fails (HTTP 403 on the
  bucket listing).
- The vertical axis is `elevation` with `positive = "up"` (−5728 to −0.494
  m), returned as is, and time comes back as `Date`. The package would flip
  elevation to depth (D03) and keep sub-daily time itself.
- Speed depends on the store layout: a time series is read from
  `geoChunked` (134 days per chunk) and maps from `timeChunked` (one day per
  chunk). The core object would need to know both stores.

## 4. ncdfCF: lazy `subset()` and API stability

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
- **Stability:** 11 releases from June 2024 to March 2026, with API breaks
  in 0.4 (`CFData` renamed `CFArray`), 0.6 (`subset()` signature), 0.7
  (`CFArray` merged into `CFVariable`) and 0.8 (indexing must include
  length-1 axes, `CFResource` renamed `NCResource`). One maintainer (who
  also maintains CFtime); built on RNetCDF, not ncdf4. Its only CRAN
  reverse dependency is stars (Suggests).

## What changes in SPEC.md

"Planned: one model → Direction: a hybrid" and "First step" are updated with
this decision once it is agreed (#112, step 5 of the order of work).
