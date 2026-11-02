# Data model trial scripts (#112)

Scripts behind [`dev/data-model-trial.md`](../data-model-trial.md). Run them
from the package root. stars, ncdfCF and CFtime aren't package dependencies,
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
dev/data-model-trial/run_bench.sh <lib> <data> synthetic "ocean3d stars_ncdf stars_mdim ncdfCF" "load open points series volume buffer"
dev/data-model-trial/run_bench.sh <lib> <data> oisst "ocean3d stars_ncdf stars_mdim ncdfCF" "load open points series buffer"
dev/data-model-trial/run_zarr.sh <lib> "ncdf4 open series map volume"
Rscript dev/data-model-trial/time_matching.R <lib>
```
