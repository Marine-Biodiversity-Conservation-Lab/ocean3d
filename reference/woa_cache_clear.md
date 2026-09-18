# Clear the WOA cache

Remove all cached WOA NetCDF files.

## Usage

``` r
woa_cache_clear(confirm = TRUE)
```

## Arguments

- confirm:

  Logical. Require interactive confirmation. Default `TRUE`.

## Value

Invisibly, `TRUE` on success.

## Examples

``` r
# Redirect the cache into a temporary directory so this example does not
# delete real downloads. In normal use you would skip this step.
old <- Sys.getenv("R_USER_CACHE_DIR", unset = NA)
Sys.setenv(R_USER_CACHE_DIR = tempdir())

cache <- woa_cache_dir()
file.create(file.path(cache, "woa23_decav_t00_01.nc"))
#> [1] TRUE
list.files(cache)
#> [1] "woa23_decav_t00_01.nc"

# confirm = FALSE deletes without prompting (e.g. from a script)
woa_cache_clear(confirm = FALSE)
dir.exists(cache)
#> [1] FALSE

# The cache is recreated empty by the next call that needs it
woa_cache_dir()
#> [1] "/tmp/RtmpZgGnIj/R/ocean3d/woa"

if (is.na(old)) Sys.unsetenv("R_USER_CACHE_DIR") else
  Sys.setenv(R_USER_CACHE_DIR = old)

if (FALSE) { # \dontrun{
# Called with no arguments it asks for confirmation when interactive
woa_cache_clear()
} # }
```
