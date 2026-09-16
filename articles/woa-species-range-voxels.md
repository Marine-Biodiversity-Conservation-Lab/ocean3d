# Extracting Voxel Values Within a Species' 3D Range

## Extracting Voxel Values Within a Species’ 3D Range

This article walks the core operation the package exists for: take a
species range **polygon**, turn it into a 3D domain, and read the values
a depth-stratified environmental raster carries inside that domain.

Three objects and two conversions do the whole job:

| Object | What it is |
|----|----|
| `SpatEnvelope` | The species’ 3D range. Two layers, `depth_min` and `depth_max`; **depth is the cell value**, so the vertical interval varies per cell. |
| `SpatVoxel` | The environmental raster. One layer per standard depth, named `{variable}_depth={value}`; **depth is the layer index**. |
| occupancy voxel | The envelope re-expressed on the raster’s depth axis, so the two can be combined layer by layer. |

Unlike the other articles, everything here is synthetic and runs — no
downloads, no local data files. The grids are small enough to print, so
you can check every step against the numbers.

``` r

library(ocean3d)
```

### Step 1: A study grid and a seafloor

The study grid defines the geometry everything else is aligned to. The
seafloor is depth in **positive metres increasing downward**, the
package’s convention throughout — GEBCO-style elevation has to be
flipped before it gets here.

``` r

grid <- terra::rast(
  nrows = 20, ncols = 20,
  xmin = 0, xmax = 20, ymin = 0, ymax = 20,
  crs = "EPSG:4326"
)

# A shelf that deepens from north to south, 20 m down to 900 m.
seafloor <- terra::rast(grid)
terra::values(seafloor) <- rep(seq(20, 900, length.out = 20), each = 20)

seafloor
#> class       : SpatRaster
#> size        : 20, 20, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 20, 0, 20  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> source(s)   : memory
#> name        : lyr.1
#> min value   :    20
#> max value   :   900
```

### Step 2: A WOA-shaped environmental raster, as a `SpatVoxel`

World Ocean Atlas products arrive as one layer per standard depth. Here
that shape is built by hand: temperature falling with depth, on a
latitudinal gradient.
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
attaches the layer names and validates the depth axis.

``` r

standard_depths <- c(0, 50, 100, 200, 500)

lat <- terra::init(terra::rast(grid), "y")
t_layers <- lapply(standard_depths, function(d) 26 - 0.03 * d - 0.4 * lat)

t_annual <- as_voxel(
  terra::rast(t_layers),
  depths = standard_depths,
  varname = "t_an"
)

names(t_annual)
#> [1] "t_an_depth=0"   "t_an_depth=50"  "t_an_depth=100" "t_an_depth=200"
#> [5] "t_an_depth=500"
```

[`depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depths.md)
reads that axis back out. It is the accessor to reach for whenever you
need the numeric depths rather than the layer names:

``` r

depths(t_annual)
#> [1]   0  50 100 200 500
```

### Step 3: The species range polygon, as a `SpatEnvelope`

[`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
rasterises the polygon onto the study grid and attaches the depth
limits, giving the species’ 3D domain.

The seafloor is not a special argument — it is simply one more
`depth_max` constraint. Per cell the **shallowest** `depth_max` wins, so
a species with a nominal 600 m limit is clamped to the bed wherever the
bed is shallower. Cells where the bed sits above `depth_min` have no
water column left and drop out entirely.

``` r

poly <- terra::vect(
  "POLYGON ((4 2, 16 2, 16 18, 4 18, 4 2))",
  crs = "EPSG:4326"
)

range_env <- vect_to_envelope(
  polygon   = poly,
  template  = grid,
  depth_min = 0,
  depth_max = list(600, seafloor)
)

range_env
#> class       : SpatRaster
#> size        : 20, 20, 2  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 20, 0, 20  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> source(s)   : memory
#> names       : depth_min,  depth_max
#> min values  :         0, 112.631579
#> max values  :         0,        600
```

Depth is the cell value here, and the interval is per cell — the
northern rows are clamped to a shallow bed, the southern rows reach the
full 600 m:

``` r

range(terra::values(range_env[["depth_max"]]), na.rm = TRUE)
#> [1] 112.6316 600.0000
```

### Step 4: Put the envelope on the voxel’s depth axis

The envelope and the voxel store depth in different roles, so they
cannot be combined directly.
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
re-expresses the envelope on the voxel’s depth levels, giving a `1`/`NA`
occupancy voxel with one layer per depth:

``` r

occ <- envelope_to_voxel(range_env, depths = depths(t_annual))

# Cells occupied at each standard depth.
terra::global(occ, "sum", na.rm = TRUE)
#>                    sum
#> presence_depth=0   192
#> presence_depth=50  192
#> presence_depth=100 192
#> presence_depth=200 168
#> presence_depth=500  84
```

#### A depth level stands for a slab, not a knife edge

A cell occupies a level when its envelope overlaps the slab that level
stands for by more than a shared edge. This matters for thin intervals.
An envelope of `[10, 20]` contains none of `c(0, 100, 200)`, but it is
clearly within the water the 0 m level represents:

``` r

fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
terra::values(fp) <- c(1, 1)

thin <- as_envelope(fp, depth_min = 10, depth_max = 20)
terra::values(envelope_to_voxel(thin, depths = c(0, 100, 200)))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]                1                 NA                 NA
#> [2,]                1                 NA                 NA
```

Touching is not overlapping. An envelope of `[0, 100]` ends exactly
where the 100 m slab begins, so it occupies the 0 m level alone, and two
envelopes that meet at 100 m land on disjoint levels and do not
intersect as voxels, the same answer
[`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
gives for the envelopes themselves:

``` r

upper <- as_envelope(fp, depth_min = 0, depth_max = 100)
lower <- as_envelope(fp, depth_min = 100, depth_max = 200)
terra::values(envelope_to_voxel(upper, depths = c(0, 100, 200)))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]                1                 NA                 NA
#> [2,]                1                 NA                 NA
terra::values(intersects_3d(envelope_to_voxel(upper, depths = c(0, 100, 200)),
                            envelope_to_voxel(lower, depths = c(0, 100, 200))))
#>      intersects
#> [1,]      FALSE
#> [2,]      FALSE
```

Where the slab edges fall is a property of the dataset the levels came
from, so `bounds` selects the convention. `"top"` (the default) runs
each slab from its level down to the next; `"midpoint"` puts the edges
halfway between neighbours, which is what WOA means. An interval of
`[60, 70]` lands in a different level under each:

``` r

shallow <- as_envelope(fp, depth_min = 60, depth_max = 70)

terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200)))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]                1                 NA                 NA
#> [2,]                1                 NA                 NA
terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200),
                                bounds = "midpoint"))
#>      presence_depth=0 presence_depth=100 presence_depth=200
#> [1,]               NA                  1                 NA
#> [2,]               NA                  1                 NA
```

### Step 5: Read the values inside the range

With both objects on the same depth axis, masking is ordinary `terra`:

``` r

in_range <- terra::mask(t_annual, occ)

names(in_range)
#> [1] "t_an_depth=0"   "t_an_depth=50"  "t_an_depth=100" "t_an_depth=200"
#> [5] "t_an_depth=500"
```

The layer names survive, so a single depth can still be pulled out by
name — sea surface temperature within the range, for instance:

``` r

sst <- in_range[["t_an_depth=0"]]
summary(terra::values(sst, na.rm = TRUE)[, 1])
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    19.0    20.5    22.0    22.0    23.5    25.0
```

#### Vertical profile

[`depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depths.md)
also accepts a character vector, which is what makes a per-layer
[`terra::global()`](https://rspatial.github.io/terra/reference/global.html)
result easy to turn into a profile — the layer names are on the rows,
not on a raster:

``` r

per_depth <- terra::global(in_range, c("mean", "min", "max"), na.rm = TRUE)
per_depth$depth <- depths(rownames(per_depth))
per_depth$n_cells <- terra::global(!is.na(in_range), "sum", na.rm = TRUE)[, 1]

per_depth[, c("depth", "n_cells", "mean", "min", "max")]
#>                depth n_cells mean  min  max
#> t_an_depth=0       0     192 22.0 19.0 25.0
#> t_an_depth=50     50     192 20.5 17.5 23.5
#> t_an_depth=100   100     192 19.0 16.0 22.0
#> t_an_depth=200   200     168 16.4 13.8 19.0
#> t_an_depth=500   500      84  8.8  7.6 10.0
```

The cell count falls with depth: the range is clamped to the bed, so the
northern rows drop out of the deeper levels.

#### Summary statistics across the whole 3D range

Summaries over the masked voxel are a few lines of base R over
[`terra::values()`](https://rspatial.github.io/terra/reference/values.html),
which returns a cells × depths matrix. `n_cells` counts occupied
**voxels** (cell × depth), while `n_surface_cells` counts the map cells
with data at any depth:

``` r

summarise_voxel <- function(v, name) {
  vals <- terra::values(v)
  out <- c(
    min             = suppressWarnings(min(vals, na.rm = TRUE)),
    max             = suppressWarnings(max(vals, na.rm = TRUE)),
    mean            = suppressWarnings(mean(vals, na.rm = TRUE)),
    n_surface_cells = sum(rowSums(!is.na(vals)) > 0),
    n_cells         = sum(!is.na(vals)),
    n_depths        = terra::nlyr(v)
  )
  # An all-NA range gives -Inf / Inf from min() and max().
  out[is.infinite(out)] <- NA_real_
  names(out) <- paste(name, names(out), sep = "_")
  as.data.frame(as.list(out))
}

summarise_voxel(in_range, "temperature")
#>   temperature_min temperature_max temperature_mean temperature_n_surface_cells
#> 1             7.6              25         18.48116                         192
#>   temperature_n_cells temperature_n_depths
#> 1                 828                    5
```

Across several variables, and several species, this is an
[`lapply()`](https://rdrr.io/r/base/lapply.html) over the same three
lines:

``` r

# A second variable, on its own set of standard depths.
do_annual <- as_voxel(
  terra::rast(lapply(c(0, 100, 300), function(d) 300 - 0.2 * d + 2 * lat)),
  depths = c(0, 100, 300),
  varname = "o_an"
)

env_rasters <- list(temperature = t_annual, oxygen = do_annual)

do.call(cbind, unname(Map(function(v, nm) {
  summarise_voxel(mask(v, range_env), nm)
}, env_rasters, names(env_rasters))))
#>   temperature_min temperature_max temperature_mean temperature_n_surface_cells
#> 1             7.6              25         18.48116                         192
#>   temperature_n_cells temperature_n_depths oxygen_min oxygen_max oxygen_mean
#> 1                 828                    5        245        335    295.9302
#>   oxygen_n_surface_cells oxygen_n_cells oxygen_n_depths
#> 1                    192            516               3
```

[`mask()`](https://rspatial.github.io/terra/reference/mask.html) puts
the envelope on each raster’s own depth levels before it masks, so each
environmental raster gets its own occupancy voxel. Note `oxygen` reports
three depths against temperature’s five: the levels are read from the
raster being masked, never assumed shared.

### Step 6: Collapsing back to an envelope

[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
is the return trip: it reduces a voxel to the envelope bounding it. A
predicate decides what counts as present, so it can bound a subset of
the values rather than just the data extent — here, the depths over
which the range is warmer than 20 °C:

``` r

warm <- voxel_to_envelope(in_range, fun = function(x) x > 20)
range(terra::values(warm[["depth_max"]]), na.rm = TRUE)
#> [1]   0 100
```

This conversion is **lossy**. An envelope holds one continuous interval
per cell, so any interior gap in the voxel is filled in — a cell
satisfying the predicate at 0 m and 200 m but not at 100 m still comes
back as `[0, 200]`. Where interior gaps matter, stay in the voxel.

### Summary

``` r

# 1. the species' 3D range
range_env <- vect_to_envelope(poly, grid, depth_min = 0,
                              depth_max = list(600, seafloor))

# 2. the environmental raster
t_annual <- as_voxel(stack, depths = standard_depths, varname = "t_an")

# 3. the values inside the range
in_range <- mask(t_annual, range_env)
```

Step 3 is the one line to reach for whenever a voxel needs restricting
to a species’ per-cell depth window.
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) is
depth-aware for the package’s 3D classes: it places the envelope on the
voxel’s own depth levels, as
[`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
did by hand above, and masks layer by layer.

``` r

mask(rast_3d, range_env)
```

To restrict a voxel to a polygon and a single depth band across the
whole area instead — no per-cell window — use
[`extract_to_area()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_area.md).
