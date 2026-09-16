# Download a WOA NetCDF file (with caching)

Download World Ocean Atlas 2023 NetCDF files from the NCEI THREDDS
server. Files are cached in
[`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md)
(or a user-supplied `output_dir`) and skipped on subsequent calls unless
`force = TRUE`.

## Usage

``` r
woa_download(
  variable,
  period = "annual",
  resolution = "0.25",
  decade = NULL,
  output_dir = NULL,
  force = FALSE,
  quiet = FALSE
)
```

## Arguments

- variable:

  Character. One of `"temperature"`, `"salinity"`, `"dissolved_oxygen"`,
  `"oxygen_saturation"`, `"AOU"`, `"nitrate"`, `"phosphate"`,
  `"silicate"`, `"density"`.

- period:

  Character or numeric. `"annual"` (default), `"monthly"` (all 12
  months), `"seasonal"` (4 seasons), or a numeric vector where 0 =
  annual, 1:12 = monthly, 13:16 = seasonal.

- resolution:

  Character or numeric. `"0.25"` (default), `"1"`, or `"5"` degrees.

- decade:

  Character. Decade code (e.g. `"decav"`, `"all"`). Defaults to the
  canonical decade for the variable.

- output_dir:

  Character. Destination directory. Defaults to
  [`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md).

- force:

  Logical. If `TRUE`, re-download even if the file exists. Default
  `FALSE`.

- quiet:

  Logical. Suppress download progress. Default `FALSE`.

## Value

Character vector of paths to downloaded .nc files.

## Details

On first use of the persistent cache in an interactive session, the
function prompts for consent to write to user filespace (per CRAN
Repository Policy). When the `curl` package is available, the estimated
size of each remote file is reported before downloading. A full set of
WOA files can be many gigabytes — supply `output_dir` to direct the
cache elsewhere (HPC scratch, external drive, etc.).

## Examples

``` r
if (FALSE) { # \dontrun{
# Annual temperature climatology at 1-degree resolution, cached in a
# temporary directory so nothing is written to user filespace
tmp <- file.path(tempdir(), "woa")
f <- woa_download("temperature", period = "annual", resolution = "1",
                  output_dir = tmp)
woa_load_nc(f)

# All 12 monthly files, downloaded in parallel, into the default cache
monthly <- woa_download("temperature", period = "monthly",
                        resolution = "1")
length(monthly)

# Specific periods: 0 = annual, 1:12 = monthly, 13:16 = seasonal
winter <- woa_download("dissolved_oxygen", period = 13, resolution = "1")

# Other variables use their own canonical decade codes. Note that
# 0.25-degree files are large (~1 GB each for salinity).
sal <- woa_download("salinity", period = "annual", resolution = "0.25")

# Re-download a file that is already cached
f <- woa_download("temperature", period = "annual", resolution = "1",
                  output_dir = tmp, force = TRUE, quiet = TRUE)
} # }
```
