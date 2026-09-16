# Coerce a 2D footprint to a 2.5D min-max envelope

Build a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
— one continuous `[depth_min, depth_max]` interval per grid cell — from
a 2D horizontal footprint and the depth limits that apply to it. This is
the general form of the conversion
[`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
performs for polygons: any 2D raster whose non-`NA` cells mark presence
becomes a 3D domain once depth limits are attached to it.

## Usage

``` r
as_envelope(x, depth_min, depth_max)
```

## Arguments

- x:

  SpatRaster. Either a single-layer footprint whose non-`NA` cells are
  present, or a two-layer raster already named `depth_min`, `depth_max`.

- depth_min, depth_max:

  Shallowest and deepest depth in metres, positive down. Each is either
  a single number applying to the whole footprint, or a single-layer
  SpatRaster on the grid of `x` giving the limit per cell. Omit both
  when `x` already carries the two depth layers.

## Value

A
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
with layers `depth_min` and `depth_max`, on the grid of `x`.

## Details

`x` is the footprint. Its values are not read, only their non-`NA`
pattern: a rasterized species range, a Global Fishing Watch effort
layer, or a plain presence mask all work. A SpatRaster that *already*
has exactly the two layers `depth_min` and `depth_max` (the output of
[`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md),
or a
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
itself) is promoted directly instead, in which case `depth_min` and
`depth_max` must be omitted.

Depths are positive metres increasing downward, following the package's
depth sign convention. Flip GEBCO-style elevation first, e.g.
`terra::clamp(-terra::project(bathy, x), lower = 0)`, which gives a
positive seafloor depth with land clamped to 0.

## Examples

``` r
# A 2x2 footprint: three cells present (any non-NA value), one absent.
fp <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
terra::values(fp) <- c(1, 1, 1, NA)

# A species recorded between 0 and 200 m.
e <- as_envelope(fp, depth_min = 0, depth_max = 200)
terra::values(e)
#>      depth_min depth_max
#> [1,]         0       200
#> [2,]         0       200
#> [3,]         0       200
#> [4,]        NA        NA

# Per-cell limits are allowed too, as single-layer rasters.
dmax <- terra::setValues(terra::rast(fp), c(100, 200, 300, 400))
terra::values(as_envelope(fp, depth_min = 10, depth_max = dmax))
#>      depth_min depth_max
#> [1,]        10       100
#> [2,]        10       200
#> [3,]        10       300
#> [4,]        NA        NA
```
