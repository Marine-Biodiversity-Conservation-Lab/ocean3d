# Load environmental data from Copernicus services

Acquire environmental data from Copernicus Marine or ECMWF-operated
Copernicus Data Stores and save the resulting files locally. Supported
sources are Copernicus Marine (`"marine"`), Climate Data Store
(`"cds"`), and Atmosphere Data Store (`"ads"`).

## Usage

``` r
copernicus_load(
  source = c("marine", "cds", "ads"),
  backend = c("auto", "standalone", "python"),
  dataset_id,
  variables = NULL,
  start_datetime = NULL,
  end_datetime = start_datetime,
  xmin = NULL,
  xmax = NULL,
  ymin = NULL,
  ymax = NULL,
  depth_min = NULL,
  depth_max = NULL,
  request = list(),
  output_dir = NULL,
  filename = NULL,
  split_by = NULL,
  organize_by = NULL,
  hemisphere = c("north", "south"),
  concurrent_processes = NULL,
  compression = NULL,
  dataset_version = NULL,
  dataset_part = NULL,
  force = FALSE,
  quiet = FALSE
)
```

## Arguments

- source:

  Character. Copernicus service: `"marine"`, `"cds"`, or `"ads"`.

- backend:

  Character. Backend used for Copernicus Marine requests: `"auto"`,
  `"standalone"`, or `"python"`.

- dataset_id:

  Character. Provider dataset identifier.

- variables:

  Character vector containing one or more environmental variables to
  download. For Copernicus Marine this must be supplied explicitly. Each
  requested variable is stored in separate NetCDF files.

- start_datetime:

  Optional start datetime.

- end_datetime:

  Optional end datetime.

- xmin, xmax, ymin, ymax:

  Optional geographic bounding box.

- depth_min, depth_max:

  Optional depth range in metres, positive downward.

- request:

  Named list of provider- or dataset-specific request arguments.

- output_dir:

  Character. Root destination directory.

- filename:

  Optional base output filename.

- split_by:

  Optional character defining temporal file subdivision for Copernicus
  Marine. Supported values are `"year"`, `"season"`, `"month"`,
  `"week"`, `"day"`, and `"hour"`. Variables are always split
  automatically.

- organize_by:

  Optional character vector defining the directory hierarchy. Supported
  components are `"service"`, `"dataset"`, `"variable"`, `"year"`,
  `"season"`, `"month"`, `"week"`, `"day"`, and `"hour"`.

- hemisphere:

  Character. `"north"` or `"south"`. Used only to assign human-readable
  names to meteorological seasons.

- concurrent_processes:

  Optional integer \>= 1.

- compression:

  Optional integer from 0 to 9.

- dataset_version:

  Optional Copernicus Marine dataset version.

- dataset_part:

  Optional Copernicus Marine dataset part.

- force:

  Logical. Replace existing output files.

- quiet:

  Logical. Suppress progress messages where possible.

## Value

Character vector containing paths to downloaded files.

## Details

Copernicus Marine requests can use either the official standalone
Copernicus Marine Toolbox or the Python `copernicusmarine` package
through `reticulate`. The standalone Toolbox is the recommended backend
because it does not require users to install or configure Python.

With `backend = "auto"`, the standalone Toolbox is used when available.
The Python backend is used only when the standalone executable is not
available. Backend selection occurs before the request starts; a failed
request is never automatically retried using another backend.

For Copernicus Marine, each output NetCDF file contains exactly one
requested environmental variable. When several variables are supplied,
they are automatically separated into different files. `split_by`
therefore controls temporal subdivision only.

Downloaded Copernicus Marine files can optionally be arranged into a
directory hierarchy with `organize_by`. For example,

`c("service", "dataset", "variable", "year", "month", "day")`

stores each variable independently below the service and dataset,
followed by calendar year, month, and day.

Temporal organization cannot be finer than `split_by`. For example,
`split_by = "month"` cannot be combined with
`organize_by = c("year", "month", "day")`, because one monthly output
file can contain multiple days. Such incompatible combinations are
detected before a download starts.

Week splitting and organization use ISO weeks. ISO weeks start on Monday
and week 1 is the week containing 4 January. The associated year is the
ISO week-year rather than necessarily the ordinary calendar year.

Season splitting uses meteorological rather than astronomical seasons.
Meteorological seasons are fixed three-month periods beginning on the
first day of a month:

- DJF: December-February

- MAM: March-May

- JJA: June-August

- SON: September-November

These temporal boundaries are identical in both hemispheres.
`hemisphere` determines only the human-readable season name. Thus DJF is
`"winter_DJF"` in the Northern Hemisphere and `"summer_DJF"` in the
Southern Hemisphere. DJF is assigned to the year containing January and
February; for example, December 2025 belongs to DJF 2026.

The official Copernicus Marine Toolbox natively supports splitting by
year, month, day, and hour. Week and meteorological-season splitting are
managed by this package by issuing one appropriately bounded Toolbox
request per temporal period.

CDS and ADS requests use `ecmwfr`. File organization through
`organize_by` is currently implemented for Copernicus Marine only.

Authentication is managed through the official provider clients. For
Copernicus Marine, run
[`copernicus_login()`](https://marine-biodiversity-conservation-lab.github.io/sharkabc3d/reference/copernicus_login.md)
once before the first request. The credentials created by the official
Copernicus Marine Toolbox are shared by both the standalone and Python
backends, so the same login can be used with `backend = "standalone"` or
`backend = "python"`.

If the Python backend does not find valid credentials, it asks the user
to configure authentication with
[`copernicus_login()`](https://marine-biodiversity-conservation-lab.github.io/sharkabc3d/reference/copernicus_login.md)
rather than attempting an interactive password prompt through
`reticulate`, because such prompts are not reliably supported in R
sessions.

## Examples

``` r
if (FALSE) { # \dontrun{
files <- copernicus_load(
  source = "marine",
  dataset_id = "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m",
  variables = "thetao",
  start_datetime = "2026-07-01",
  end_datetime = "2026-07-31",
  xmin = -6,
  xmax = 10,
  ymin = 35,
  ymax = 45,
  depth_min = 0,
  depth_max = 500,
  split_by = "day",
  organize_by = c(
    "service",
    "dataset",
    "variable",
    "year",
    "month",
    "day"
  )
)
} # }
```
