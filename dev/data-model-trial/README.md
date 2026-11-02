# Data model trial scripts (#112)

Scripts behind [`dev/data-model-trial.md`](../data-model-trial.md). Run them
from the package root.

The trial ran on 2 November 2026 with stars 0.7.3, ncdfCF 0.8.2, CFtime
1.7.3, ncdf4 1.24, terra 1.9.50 and GDAL 3.8.4, on R 4.5.3 (Linux, 8 cores,
31 GB RAM). Each benchmark step ran in a new R process, so peak memory is the
whole process's peak resident memory as reported by GNU `time`, including R
and the loaded packages (the "load" rows). Each process had a 12 GB
address-space limit.

The Copernicus Marine download is
`cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m`: 721 × 721 × 18 × 731 cells
(0.5–47 m, 1 June 2022 to 31 May 2024), `float`, compressed, chunks of 46 ×
46 cells × 2 levels × 23 days, longitudes 150–210. The synthetic file is
1440 × 720 × 10 × 730 cells, `short` with `scale_factor`, `_FillValue` on
land, one chunk per depth and day. stars, ncdfCF and CFtime aren't package dependencies,
so the scripts take a library path where they are installed:

```r
install.packages(c("stars", "ncdfCF", "CFtime"), lib = "<lib>")
```

| Script | What it does |
|---|---|
| `cases.R` | The CF cases (test fixtures and #88–#96), each a small netCDF file with the expected result |
| `correctness.R` | Reads every case with stars (`read_ncdf()`, `read_mdim()`), ncdfCF and the package; writes `results/correctness.csv` |
| `make_large_file.R` | Writes the synthetic 4D daily file used for the timings (8.7 GB on disk, 60 GB as doubles) |
| `bench.R` | One reader, one step, one dataset, in its own process |
| `run_bench.sh` | Runs `bench.R` for each reader and step; appends to `results/benchmarks.tsv` |
| `measure.sh` | Wall time and peak resident memory of one process (GNU `time`), with a 12 GB address-space limit |
| `zarr.R`, `run_zarr.sh` | The same steps on the real Copernicus Marine product as a remote Zarr store |
| `time_matching.R` | #89: matching 500 dates to a 10-year daily axis, package against CFtime |

```bash
Rscript dev/data-model-trial/correctness.R <lib>
Rscript dev/data-model-trial/make_large_file.R <data>/synthetic_4d.nc
# OISST: https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.mean.{2022,2023}.nc into <data>/oisst/
# Copernicus Marine (needs a login; 16.9 GB):
copernicusmarine subset --dataset-id cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m \
  --variable thetao --minimum-longitude 150 --maximum-longitude -150 \
  --minimum-latitude -30 --maximum-latitude 30 --minimum-depth 0 --maximum-depth 50 \
  --start-datetime 2022-06-01 --end-datetime 2024-05-31 --netcdf-compression-level 1 \
  --output-directory <data>/cmems --output-filename cmems_thetao_pacific_2022-2024.nc
dev/data-model-trial/run_bench.sh <lib> <data> cmems "ocean3d stars_ncdf stars_mdim ncdfCF" "load open points series volume buffer"
dev/data-model-trial/run_bench.sh <lib> <data> synthetic "ocean3d stars_ncdf stars_mdim ncdfCF" "load open points series volume buffer"
dev/data-model-trial/run_bench.sh <lib> <data> oisst "ocean3d stars_ncdf stars_mdim ncdfCF" "load open points series buffer"
dev/data-model-trial/run_zarr.sh <lib> "ncdf4 open series map volume"
Rscript dev/data-model-trial/time_matching.R <lib>
```
