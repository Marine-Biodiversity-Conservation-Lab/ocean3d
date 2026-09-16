# Check Copernicus Marine Toolbox availability

Check whether the official Copernicus Marine Toolbox standalone
executable is available and report the executable path and version.

## Usage

``` r
copernicus_status(path = NULL)
```

## Arguments

- path:

  Optional character. Explicit path to a Copernicus Marine Toolbox
  executable.

## Value

A named list containing `available`, `path`, and `version`.

## Details

The executable is searched for in the package-managed Copernicus
directory and on the system `PATH`. An explicit executable path can also
be supplied.

## Examples

``` r
copernicus_status()
#> Copernicus Marine Toolbox standalone executable was not found.
#> Run `copernicus_setup()` to install the supported standalone Toolbox.
```
