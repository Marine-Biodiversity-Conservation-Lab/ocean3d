# Per-cell 3D volume overlap between two rasterized domains

Computes the depth interval and volume of each domain and of their
intersection, cell by cell. The intersection is
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md);
this function measures it. Accepts any combination of
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
and
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md);
when the two differ, the envelope is discretized onto the voxel's depth
levels with
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
rather than the voxel being collapsed, because
[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
fills interior gaps and would overstate the overlap.

## Usage

``` r
calc_volume_overlap(x, y, ...)

# S4 method for class 'SpatEnvelope,SpatEnvelope'
calc_volume_overlap(x, y, ...)

# S4 method for class 'SpatVoxel,SpatVoxel'
calc_volume_overlap(x, y, bounds = c("top", "midpoint"))

# S4 method for class 'SpatEnvelope,SpatVoxel'
calc_volume_overlap(x, y, bounds = c("top", "midpoint"), ...)

# S4 method for class 'SpatVoxel,SpatEnvelope'
calc_volume_overlap(x, y, bounds = c("top", "midpoint"), ...)

# S4 method for class 'SpatRaster,SpatRaster'
calc_volume_overlap(x, y, ...)
```

## Arguments

- x, y:

  The two domains, each a
  [SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
  or
  [SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md).

- ...:

  Arguments for the voxel methods, which are an error for a pair of
  envelopes.

- bounds:

  Voxel methods only. See
  [`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md).
  Also governs how an envelope is discretized when the two inputs
  differ.

## Value

Multi-layer SpatRaster with 9 layers:

- depth_min_a, depth_max_a:

  Depth limits of `x` (m)

- depth_min_b, depth_max_b:

  Depth limits of `y` (m)

- depth_min_overlap, depth_max_overlap:

  Depth limits of the intersection (m). NA where the domains do not
  overlap in depth.

- volume_a, volume_b:

  Per-cell volume of each domain (km³)

- volume_overlap:

  Per-cell overlap volume (km³). `0` where both are present but their
  depths do not intersect, NA where either is absent.

For voxel inputs the depth layers are *summary bounds* — the shallowest
and deepest occupied level. A voxel with an interior gap is not solid
between them, so read the volume layers, not the depth layers, for
occupied volume. The result is a plain SpatRaster: it is neither an
envelope nor a voxel.

## Details

Two voxels must be sampled at the same depth levels, and all inputs must
be on the same grid. To ask only whether two domains overlap, use
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md);
it computes no volumes.

## See also

[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
for the shared domain itself, and
[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
for its total volume;
[`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
when only the presence of overlap is needed.

## Examples

``` r
fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2000,
                  ymin = 0, ymax = 1000,
                  crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
terra::values(fp) <- c(1, 1)

a <- as_envelope(fp, depth_min = 0, depth_max = 100)
b <- as_envelope(fp, depth_min = 50, depth_max = 200)

ov <- calc_volume_overlap(a, b)
names(ov)
#> [1] "depth_min_a"       "depth_max_a"       "depth_min_b"      
#> [4] "depth_max_b"       "depth_min_overlap" "depth_max_overlap"
#> [7] "volume_a"          "volume_b"          "volume_overlap"   
terra::global(ov[[c("volume_a", "volume_b", "volume_overlap")]], "sum",
              na.rm = TRUE)
#>                sum
#> volume_a       0.2
#> volume_b       0.3
#> volume_overlap 0.1

# Mixed input: the envelope is discretized onto the voxel's levels.
bv <- envelope_to_voxel(b, depths = c(0, 50, 100, 150, 200))
terra::global(calc_volume_overlap(a, bv)[["volume_overlap"]], "sum",
              na.rm = TRUE)
#>                sum
#> volume_overlap 0.1
```
