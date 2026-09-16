# Fill missing depth values

Fix swapped upper/lower depth values and fill NAs using genus-level
means. IUCN Red List assessments (see
[`fetch_species_assessments()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/fetch_species_assessments.md))
sometimes record depth limits the wrong way around or omit them
entirely. Designed for use inside
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) —
returns a two-column tibble (`upper_depth` and `lower_depth`) that can
be unpacked with `mutate()`.

## Usage

``` r
fill_missing_depths(upper, lower, genus, method = "genus_mean")
```

## Arguments

- upper:

  Numeric vector. Upper (shallower) depth limit values, possibly with
  NAs or swapped values.

- lower:

  Numeric vector. Lower (deeper) depth limit values, possibly with NAs
  or swapped values.

- genus:

  Character vector. Genus names, used to compute genus-level mean depths
  for filling NAs.

- method:

  Character. Method for filling missing values. Currently only
  `"genus_mean"` is supported. Default `"genus_mean"`.

## Value

A tibble with columns `upper_depth` and `lower_depth`, suitable for use
with
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).

## Examples

``` r
# Depth limits as they can arrive from an IUCN assessment: the first row is
# swapped (upper deeper than lower) and the third is missing entirely.
assessments <- data.frame(
  genus_name = c("Carcharhinus", "Carcharhinus", "Carcharhinus", "Sphyrna"),
  upper_depth_limit = c(280, 0, NA, 0),
  lower_depth_limit = c(0, 100, NA, 512)
)

fill_missing_depths(
  assessments$upper_depth_limit,
  assessments$lower_depth_limit,
  assessments$genus_name
)
#>   upper_depth lower_depth
#> 1           0         280
#> 2           0         100
#> 3           0         190
#> 4           0         512

# A genus with no non-missing values anywhere has nothing to fill from:
fill_missing_depths(c(NA, NA), c(NA, NA), c("Raja", "Raja"))
#>   upper_depth lower_depth
#> 1         NaN         NaN
#> 2         NaN         NaN

# Typical use, unpacking both columns inside dplyr::mutate():
if (FALSE) { # \dontrun{
assessments <- assessments %>%
  mutate(fill_missing_depths(upper_depth_limit, lower_depth_limit, genus_name))
} # }
```
