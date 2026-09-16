# Reduce a variable voxel to a presence voxel

A
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
holds a variable, such as temperature, oxygen or fishing effort. It does
not say which cells count as occupied. `occupied()` applies a predicate
and returns presence: 1 where `fun` is `TRUE`, `NA` elsewhere. The grid
and the depth levels stay the same.

## Usage

``` r
occupied(x, fun = function(v) !is.na(v))
```

## Arguments

- x:

  A
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).

- fun:

  Function taking a vector of cell values and returning a logical vector
  of the same length. `NA` counts as `FALSE`. Defaults to
  `function(v) !is.na(v)`, which treats any recorded value as occupied.

## Value

A
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
on the depth levels of `x`. Values are 1 (occupied) or `NA`. Layers are
named `presence_depth={value}`.

## Details

The 3D query verbs read presence directly and take no predicate.
Threshold each voxel first, so both sides of a query use their own
cutoff:

    intersect_3d(occupied(temp, function(v) v > 15),
                 occupied(effort, function(v) v > 0))

`fun` receives one depth layer's cell values as a numeric vector. It
must return a logical vector of the same length. A predicate that
summarises the layer, such as `function(v) v > mean(v, na.rm = TRUE)`,
therefore compares each cell against its own depth.
[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
uses the same contract.

Running `occupied()` on a presence voxel returns it unchanged.

## See also

[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
collapses the result to a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md).
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md),
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md),
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
and [`mask()`](https://rspatial.github.io/terra/reference/mask.html)
take a presence voxel.

## Examples

``` r
r <- terra::rast(nrows = 1, ncols = 3, nlyrs = 2)
terra::values(r) <- cbind(c(12, 18, 4), c(9, 16, 2))
v <- as_voxel(r, depths = c(0, 100), varname = "temp")

# Occupied where the water is warmer than 10 degrees.
warm <- occupied(v, function(x) x > 10)
names(warm)
#> [1] "presence_depth=0"   "presence_depth=100"
terra::values(warm)
#>      presence_depth=0 presence_depth=100
#> [1,]                1                 NA
#> [2,]                1                  1
#> [3,]               NA                 NA

# The verbs need no predicate.
volume(warm)
#> [1] 34004375

# A presence voxel is returned unchanged.
identical(terra::values(occupied(warm)), terra::values(warm))
#> [1] TRUE

# Each side of a query uses its own cutoff.
effort <- as_voxel(terra::setValues(r, cbind(c(0, 3, 0), c(0, 2, 0))),
                   depths = c(0, 100), varname = "effort")
terra::values(intersect_3d(warm, occupied(effort, function(x) x > 0)))
#>      presence_depth=0 presence_depth=100
#> [1,]               NA                 NA
#> [2,]                1                  1
#> [3,]               NA                 NA
```
