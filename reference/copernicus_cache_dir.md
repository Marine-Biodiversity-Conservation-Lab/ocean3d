# Copernicus cache directory

Returns the path to the package's persistent cache directory for
downloaded Copernicus environmental data. Uses
[`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) so the
location survives across sessions and follows platform conventions.

## Usage

``` r
copernicus_cache_dir()
```

## Value

Character. Path to cache directory (created if missing).

## Examples

``` r
copernicus_cache_dir()
#> [1] "/home/runner/.cache/R/ocean3d/copernicus"
```
