# Set up the Copernicus Marine Toolbox

Download and configure a standalone version of the official Copernicus
Marine Toolbox for use by this package. The standalone executable does
not require a separate Python installation.

## Usage

``` r
copernicus_setup(
  version = .COPERNICUS_TOOLBOX_VERSION,
  force = FALSE,
  quiet = FALSE
)
```

## Arguments

- version:

  Character. Copernicus Marine Toolbox version to install. Defaults to
  the version tested by this package.

- force:

  Logical. If `TRUE`, reinstall the executable even when the requested
  version is already available.

- quiet:

  Logical. Suppress setup messages where possible.

## Value

Invisibly, the result of
[`copernicus_status()`](https://marine-biodiversity-conservation-lab.github.io/sharkabc3d/reference/copernicus_status.md)
for the installed executable.

## Details

The executable is downloaded from the official versioned Copernicus
Marine Toolbox release and stored in the user's platform-specific
application data directory. Nothing is downloaded when this package
itself is installed or loaded.

## Examples

``` r
if (FALSE) { # \dontrun{
copernicus_setup()
} # }
```
