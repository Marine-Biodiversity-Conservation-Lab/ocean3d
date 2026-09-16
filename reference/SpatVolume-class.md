# Union of the package's 3D domain representations

`SpatVolume` is a class union over
[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
and
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md).
It exists as a dispatch target for operations that are meaningful on
either representation because both determine a **3D domain over a 2D
grid**: volume, volumetric overlap, vertical extent, printing.

## Details

It deliberately asserts no shared structure. The two members store depth
in dual roles (layer index vs. cell value), so any operation that
touches the depth axis directly should dispatch on the concrete class
instead.

## See also

[SpatVoxel](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
and
[SpatEnvelope](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md),
the two members;
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
and
[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
convert between them.

## Examples

``` r
# The same 3D domain in both representations.
fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
terra::values(fp) <- c(1, 1)
e <- as_envelope(fp, depth_min = 0, depth_max = 200)
v <- envelope_to_voxel(e, depths = c(0, 100, 200))

# Neither inherits from the other; both belong to the union.
methods::is(e, "SpatVolume")
#> [1] TRUE
methods::is(v, "SpatVolume")
#> [1] TRUE
methods::is(v, "SpatEnvelope")
#> [1] FALSE

# `SpatVolume` is a dispatch target, so a generic written against it accepts
# either representation.
setGeneric("n_depth_layers", function(x) standardGeneric("n_depth_layers"))
#> [1] "n_depth_layers"
setMethod("n_depth_layers", "SpatVolume", function(x) terra::nlyr(x))
n_depth_layers(e)
#> [1] 2
n_depth_layers(v)
#> [1] 3

# But the depth axis itself is stored differently in the two, so anything
# reading depths must dispatch on the concrete class rather than the union.
names(e)
#> [1] "depth_min" "depth_max"
names(v)
#> [1] "presence_depth=0"   "presence_depth=100" "presence_depth=200"
```
