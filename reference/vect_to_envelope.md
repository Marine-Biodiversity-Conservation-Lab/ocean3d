# Convert SpatVector or sf to SpatEnvelope

Rasterize polygons onto a template raster grid and assign per-cell depth
limits, outputting an SpatEnvelope object. The depths are clamped by
depth_min (shallowest) and depth_max (deepest), and/or by rasters (ex.
bathymetry) that are the same coordinate and resolution as the template.
Where a constraint raster is NA at a cell its limit there is unknown, so
the cell is dropped rather than falling back to the remaining
constraints.

## Usage

``` r
vect_to_envelope(polygon, template, depth_min, depth_max)
```

## Arguments

- polygon:

  sf or SpatVector, for example species ranges or fishery footprints.

- template:

  SpatRaster that defines the horizontal grid of the SpatEnvelope
  output.

- depth_min:

  List. Can contain both numeric and SpatRasters that match CRS,
  resolution, extent of
  ``` template`` and contain numeric values. For each cell, the maximum across the  ```depth_min\`
  list parameters is used as output SpatEnvelope depth_min layer cell
  value.

- depth_max:

  List. Can contain both numeric and SpatRasters that match CRS,
  resolution, extent of
  ``` template`` and contain numeric values. For each cell, the minimum across the  ```depth_max`list parameters is used as output SpatEnvelope depth_max layer cell value. Note that when`depth_min`is exactly`depth_max\`,
  these cell values are dropped and replaced with NA.

## Value

SpatEnvelope, with depth_min and depth_max layers.

## See also

[`as_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_envelope.md),
which attaches depth limits to a footprint that is already a raster.

## Examples

``` r
# A 4x4 study grid, and a range polygon covering its middle. `polygon` and
# `template` must already share a CRS: this function does not reproject.
template <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
                        ymin = 0, ymax = 4, crs = "EPSG:4326")
poly <- terra::vect("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))", crs = "EPSG:4326")

# A species recorded between 50 m and 200 m throughout its range. Cells
# outside the polygon are NA.
e <- vect_to_envelope(poly, template, depth_min = 50, depth_max = 200)
terra::values(e)
#>       depth_min depth_max
#>  [1,]        NA        NA
#>  [2,]        NA        NA
#>  [3,]        NA        NA
#>  [4,]        NA        NA
#>  [5,]        NA        NA
#>  [6,]        50       200
#>  [7,]        50       200
#>  [8,]        NA        NA
#>  [9,]        NA        NA
#> [10,]        50       200
#> [11,]        50       200
#> [12,]        NA        NA
#> [13,]        NA        NA
#> [14,]        NA        NA
#> [15,]        NA        NA
#> [16,]        NA        NA

# Bathymetry is not a special argument: pass the seafloor as a second
# `depth_max` and, per cell, the shallower of the two wins. This seafloor
# deepens from north to south, so the two rows the polygon covers meet the
# bed very differently.
seafloor <- terra::setValues(terra::rast(template),
                             rep(c(20, 45, 250, 400), each = 4))
e2 <- vect_to_envelope(poly, template,
                       depth_min = 50, depth_max = list(300, seafloor))
terra::values(e2)
#>       depth_min depth_max
#>  [1,]        NA        NA
#>  [2,]        NA        NA
#>  [3,]        NA        NA
#>  [4,]        NA        NA
#>  [5,]        NA        NA
#>  [6,]        NA        NA
#>  [7,]        NA        NA
#>  [8,]        NA        NA
#>  [9,]        NA        NA
#> [10,]        50       250
#> [11,]        50       250
#> [12,]        NA        NA
#> [13,]        NA        NA
#> [14,]        NA        NA
#> [15,]        NA        NA
#> [16,]        NA        NA

# Two things happened. Where the bed is at 250 m, the nominal 300 m limit is
# truncated to it. Where the bed is at 45 m it sits above `depth_min`, so
# there is no water column left to occupy and the cell drops out entirely —
# leaving 2 of the polygon's 4 cells.
sum(!is.na(terra::values(e2[["depth_min"]])))
#> [1] 2
```
