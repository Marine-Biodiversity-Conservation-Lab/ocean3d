# Environmental values inside points, buffers, tracks and areas

This article shows how to get environmental values inside 3D shapes,
such as tag detections, buffers around them, a tag’s track, or a
protected area down to a set depth. Each example uses the same three
steps:

1.  **Build** the shapes with
    [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md),
    [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)
    or
    [`sample_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/sample_3d.md).
2.  **Intersect** each shape with the data:
    [`crop()`](https://rspatial.github.io/terra/reference/crop.html) the
    dataset to the shape, read its time steps with
    [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md),
    and keep the values inside the shape with
    [`mask()`](https://rspatial.github.io/terra/reference/mask.html).
3.  **Summarise** the values with
    [`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md).

You loop over the shapes yourself, so you can stop after any step and
look at the result. Each example ends with a table of values per shape,
which you can join back to your observations. Volumes are in km³
throughout.

## Example data

The examples use a small synthetic temperature dataset: a 0.05° grid off
eastern Australia, six depth levels from 0 to 200 m, and five daily time
steps. Temperature falls with depth, rises to the north and warms a
little each day. The code below writes it to a netCDF file, as you might
download from Copernicus Marine.

``` r

library(ocean3d)

lon <- seq(150.025, 150.975, by = 0.05)
lat <- seq(-20.975, -20.025, by = 0.05)
depth <- c(0, 10, 25, 50, 100, 200)
time <- 0:4
grid <- expand.grid(lon = lon, lat = lat, depth = depth, time = time)
temp <- 28 - 0.06 * grid$depth + 2 * (grid$lat + 20.5) + 0.2 * grid$time

dims <- list(
  ncdf4::ncdim_def("lon", "degrees_east", lon),
  ncdf4::ncdim_def("lat", "degrees_north", lat),
  ncdf4::ncdim_def("depth", "m", depth),
  ncdf4::ncdim_def("time", "days since 2020-03-01", time)
)
var <- ncdf4::ncvar_def("temp", "degC", dims, missval = -9999,
                        prec = "double")
path <- tempfile(fileext = ".nc")
nc <- ncdf4::nc_create(path, var)
ncdf4::ncvar_put(nc, var, array(temp, lengths(lapply(dims, `[[`, "vals"))))
ncdf4::nc_close(nc)
```

[`as_proxy()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_proxy.md)
describes the file without reading its values. Each later step reads
only the cells and time steps it needs, at every depth level.

``` r

proxy <- as_proxy(path)
proxy
#> class      : SpatProxy (experimental)
#> source     : file2306e43bbc1.nc
#> variable   : temp (degC)
#> dimensions : 20 x 20 x 6 x 5 (lon x lat x depth x time)
#> lon        : 150.025 to 150.975, 20 cells, edges inferred
#> lat        : -20.975 to -20.025, 20 cells, edges inferred
#> depth      : 0 to 200 m, 6 levels, shallow first, edges inferred
#> time       : 2020-03-01 to 2020-03-05, 5 steps, standard calendar, no bounds
#> crs        : EPSG:4326
```

## Points: tag detections

Three detections, each with a position, a depth and a date, become
`SpatVector3D` points:

``` r

obs <- data.frame(
  tag = c("A", "B", "C"),
  lon = c(150.21, 150.52, 150.80),
  lat = c(-20.31, -20.56, -20.82),
  depth = c(12, 48, 140),
  date = as.Date(c("2020-03-01", "2020-03-02", "2020-03-04"))
)
detections <- as_vector3d(obs, time = "date", id = "tag")
detections
#> class      : SpatVector3D
#> geometry   : points 
#> features   : 3 
#> vertices   : 3 
#> depth (m)  : 12 to 140 
#> time       : 2020-03-01 to 2020-03-04
```

For one detection,
[`crop()`](https://rspatial.github.io/terra/reference/crop.html) keeps
the grid cell it falls in, and
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
reads the time step nearest its date. With a shape as `time`,
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
always returns a list, here with one voxel.
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) then
keeps the value at the depth level nearest the detection’s depth:

``` r

a <- detections[1]
v <- as_voxel(crop(proxy, a), time = a)[[1]]
summarise_3d(mask(v, a), statistic = "mean")
#>         time  mean n
#> 1 2020-03-01 27.75 1
```

Detection A at 12 m takes the 10 m level.
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
gives the same value from the file:

``` r

extract_to_point(obs[1, ], path, var = "temp", date_col = "date",
                 method = "nearest")$temp
#> [1] 27.75
```

Both pick the nearest time step: on this data they give the same values
for all three detections, for a time halfway between two steps, and for
a date after the last step.
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
warns when the step it picks is on another day;
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
doesn’t, so set `max_time_diff` in
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
to limit how far away the step may be.

The two differ when the nearest depth level has no value.
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
skips such levels and takes the nearest level with a value. By default,
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) takes
the nearest level whatever its value, so such a detection gets `NA`.
Pass `method = "nearest_value"` to skip empty levels as
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
does, or `method = "bottom"` for the deepest level with a value, as
`extract_to_point(method = "bottom")` gives.

This file has a value at every level, so the next chunk removes the 10 m
level from a copy of the voxel. Detection A at 12 m then gets `NA` by
default. With `"nearest_value"` it takes the 0 m level, 12 m away, which
is nearer than the 25 m level. `"bottom"` takes the 200 m level, the
same value
[`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
gives:

``` r

gappy <- v
gappy[[2]] <- NA_real_   # the 10 m level
summarise_3d(mask(gappy, a), statistic = "mean")
#>         time mean n
#> 1 2020-03-01   NA 0
summarise_3d(mask(gappy, a, method = "nearest_value"), statistic = "mean")
#>         time  mean n
#> 1 2020-03-01 28.35 1
summarise_3d(mask(v, a, method = "bottom"), statistic = "mean")
#>         time  mean n
#> 1 2020-03-01 16.35 1
extract_to_point(obs[1, ], path, var = "temp", date_col = "date",
                 method = "bottom")$temp
#> [1] 16.35
```

## Surface data: a variable without depth

Sea-surface temperature, chlorophyll and many satellite products have no
depth. Read them with
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)
on the proxy, which gives a plain one-layer raster with the step’s date.
The other steps are the same as for a voxel, except that
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) uses
only the detection’s position, not its depth, and
[`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md)
weights each cell by its area. Here the surface file holds the top depth
level of the example temperature.

``` r

sst <- ncdf4::ncvar_def("sst", "degC", dims[c(1, 2, 4)], missval = -9999,
                        prec = "double")
sst_path <- tempfile(fileext = ".nc")
nc <- ncdf4::nc_create(sst_path, sst)
surface <- array(temp, lengths(lapply(dims, `[[`, "vals")))[, , 1, ]
ncdf4::ncvar_put(nc, sst, surface)
ncdf4::nc_close(nc)

sst_proxy <- as_proxy(sst_path)
r <- terra::rast(crop(sst_proxy, a), time = a)[[1]]
summarise_3d(mask(r, a), statistic = "mean")
#>         time  mean n
#> 1 2020-03-01 28.35 1
```

The result is the surface temperature in the cell of detection A on its
day, from one value (`n = 1`). It is warmer than the 10 m value above,
because the example water cools with depth.

## Buffers: the water around each detection

[`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md)
turns each detection into a polygon: here a circle of 5 km radius, from
20 m above to 20 m below the detection, and one day either side.

``` r

buffers <- buffer_3d(detections, radius_km = 5, depth_window = 20,
                     time_window = 1)
as.data.frame(buffers)[c("feature_id", "depth_min", "depth_max",
                         "time_start", "time_end")]
#>   feature_id depth_min depth_max time_start   time_end
#> 1          A         0        32 2020-02-29 2020-03-02
#> 2          B        28        68 2020-03-01 2020-03-03
#> 3          C       120       160 2020-03-03 2020-03-05
```

The steps are the same as for one detection, run for each buffer in
turn.
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
returns one voxel for each time step inside a buffer’s window: three for
B and C, and two for A, whose window starts before the first day of
data. Given that list of voxels,
[`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md)
gives one row per time step. By default the mean is weighted by the
volume of water each value stands for, and `n` counts the values used.

[`mask()`](https://rspatial.github.io/terra/reference/mask.html) uses
only a buffer’s outline and depth range, because
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
has already used its time. It prints a message saying so each time,
which [`suppressMessages()`](https://rdrr.io/r/base/message.html) hides
inside the loop.

``` r

buffer_rows <- lapply(seq_along(buffers), function(i) {
  b <- buffers[i]
  voxels <- as_voxel(crop(proxy, b), time = b)
  inside <- lapply(voxels, function(v) suppressMessages(mask(v, b)))
  rows <- summarise_3d(inside, statistic = c("mean", "sd"))
  cbind(feature_id = as.data.frame(b)$feature_id, rows)
})
do.call(rbind, buffer_rows)
#>   feature_id       time     mean        sd  n
#> 1          A 2020-03-01 27.36252 0.5863953 24
#> 2          A 2020-03-02 27.56252 0.5863953 24
#> 3          B 2020-03-01 25.38426 0.7445443 16
#> 4          B 2020-03-02 25.58426 0.7445443 16
#> 5          B 2020-03-03 25.78426 0.7445443 16
#> 6          C 2020-03-03 18.32145 3.1053138 12
#> 7          C 2020-03-04 18.52145 3.1053138 12
#> 8          C 2020-03-05 18.72145 3.1053138 12
```

Each row summarises one time step;
[`summarise_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/summarise_3d.md)
doesn’t pool values across steps.

## Tracks: buffers along a route

A track is a line with a depth and a time at each vertex. This tag
records three positions over 30 hours:

``` r

route <- data.frame(
  tag = "D",
  lon = c(150.15, 150.45, 150.85),
  lat = c(-20.85, -20.50, -20.40),
  depth = c(5, 60, 20),
  time = as.POSIXct(c("2020-03-02 00:00", "2020-03-02 12:00",
                      "2020-03-03 06:00"), tz = "UTC")
)
track <- as_vector3d(route, type = "lines", id = "tag", time = "time")
```

A single buffer around the line would span its full depth and time
range, from 5 m to 60 m over two days. To follow the route more closely,
[`sample_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/sample_3d.md)
places points along it, here at most 10 km apart. Each point gets the
depth and time interpolated between the vertices either side, and
`support_km`, the length of route it stands for: half the distance to
the points either side.

``` r

samples <- sample_3d(track, spacing_km = 10)
as.data.frame(samples)[c("feature_id", "distance_km", "depth", "time",
                         "support_km")]
#>    feature_id distance_km depth                time support_km
#> 1         D_1    0.000000     5 2020-03-02 00:00:00   4.978544
#> 2         D_2    9.957087    16 2020-03-02 02:24:00   9.957087
#> 3         D_3   19.914175    27 2020-03-02 04:48:00   9.957087
#> 4         D_4   29.871262    38 2020-03-02 07:12:00   9.957087
#> 5         D_5   39.828350    49 2020-03-02 09:36:00   9.957087
#> 6         D_6   49.785437    60 2020-03-02 12:00:00   9.296733
#> 7         D_7   58.421816    52 2020-03-02 15:36:00   8.636379
#> 8         D_8   67.058195    44 2020-03-02 19:12:00   8.636379
#> 9         D_9   75.694574    36 2020-03-02 22:48:00   8.636379
#> 10       D_10   84.330953    28 2020-03-03 02:24:00   8.636379
#> 11       D_11   92.967333    20 2020-03-03 06:00:00   4.318190
```

You then buffer each point and loop over the buffers as above. Here each
buffer is a circle of 3 km radius, from 10 m above to 10 m below its
point. The loop passes the point, not the buffer, to
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md),
which then reads the time step nearest the point’s time. A buffer would
need a time window wide enough to contain a step; with daily data and a
window of a few hours,
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
stops because the window holds no step. Each row also keeps its point’s
`support_km`, so the weights stay with their rows:

``` r

sample_buffers <- buffer_3d(samples, radius_km = 3, depth_window = 10)
local <- do.call(rbind, lapply(seq_along(sample_buffers), function(i) {
  b <- sample_buffers[i]
  v <- as_voxel(crop(proxy, b), time = samples[i])[[1]]
  cbind(feature_id = as.data.frame(b)$feature_id,
        support_km = as.data.frame(b)$support_km,
        summarise_3d(suppressMessages(mask(v, b)), statistic = "mean"))
}))
local
#>    feature_id support_km       time     mean n
#> 1         D_1   4.978544 2020-03-02 27.07144 8
#> 2         D_2   9.957087 2020-03-02 26.44616 8
#> 3         D_3   9.957087 2020-03-02 26.64616 8
#> 4         D_4   9.957087 2020-03-02 25.42175 8
#> 5         D_5   9.957087 2020-03-02 25.10001 4
#> 6         D_6   9.296733 2020-03-02 25.20001 4
#> 7         D_7   8.636379 2020-03-03 25.41667 3
#> 8         D_8   8.636379 2020-03-03 26.02175 8
#> 9         D_9   8.636379 2020-03-03 26.02175 8
#> 10       D_10   8.636379 2020-03-03 26.10508 6
#> 11       D_11   4.318190 2020-03-03 27.44616 8
```

Each row is a local summary: the conditions inside one buffer. `n`
varies from buffer to buffer, here from 3 to 8 values, because the
buffers are smaller than the grid cells and the depth window catches one
or two levels. A buffer with no water, for example on land near a coast,
has a `NA` mean and `n = 0`.

To summarise the whole track, average the local means weighted by
`support_km`. This estimates the mean along the distance travelled, from
the water the buffers cover, so a stretch of route sampled more densely
doesn’t count for more. To weight by time spent instead, use the
`support_seconds` column. `na.rm = TRUE` leaves out buffers with no
water:

``` r

stats::weighted.mean(local$mean, local$support_km, na.rm = TRUE)
#> [1] 25.96283
```

The weights ignore `n`, so a buffer with 3 values counts as much as one
with 8. This is a mean of the local means, not a mean over all the water
the buffers cover. These buffers don’t overlap, but neighbouring ones
can use the same grid cells, so the local summaries aren’t independent.

[`sample_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/sample_3d.md)
doesn’t yet check that each segment stays in water, so check that your
samples aren’t on land before summarising them (#105).

## Areas: a protected area down to 50 m

An area is a polygon with a depth range. This one covers the middle of
the grid, from the surface to 50 m:

``` r

mpa <- sf::st_sf(
  name = "reef",
  geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(150.3, -20.7), c(150.7, -20.7), c(150.7, -20.3), c(150.3, -20.3),
    c(150.3, -20.7)
  ))), crs = 4326)
)
reef <- as_vector3d(mpa, depth_min = 0, depth_max = 50)
volume(reef) # km³
#> [1] 92.38822
```

The area has no time, so pick the time step with a date. `by = "depth"`
gives one row per depth level:

``` r

v <- as_voxel(crop(proxy, reef), time = "2020-03-03")
summarise_3d(mask(v, reef), statistic = c("mean", "min", "max"),
             by = "depth")
#>         time depth     mean   min   max  n
#> 1 2020-03-03     0 28.40017 28.05 28.75 64
#> 2 2020-03-03    10 27.80017 27.45 28.15 64
#> 3 2020-03-03    25 26.90017 26.55 27.25 64
#> 4 2020-03-03    50 25.40017 25.05 25.75 64
#> 5 2020-03-03   100       NA    NA    NA  0
#> 6 2020-03-03   200       NA    NA    NA  0
```

The 50 m level counts because the water it stands for, from 37.5 m to 75
m, overlaps the area’s depth range;
[`vignette("depth-levels")`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/articles/depth-levels.md)
describes this rule. The 100 m and 200 m levels are outside the area, so
their rows have `NA` and `n = 0`.

`extract_to_area(mpa, proxy, min_depth = 0, max_depth = 50, time = "2020-03-03")`
runs the same
[`crop()`](https://rspatial.github.io/terra/reference/crop.html),
[`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
and [`mask()`](https://rspatial.github.io/terra/reference/mask.html)
steps in one call, and returns the same voxel.

### Down to the seafloor

A depth range that follows the seafloor, such as a species found from
the surface down to 200 m or the bed, whichever is shallower, or a layer
within 20 m of the bed, varies from cell to cell.
[`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
builds it on the data’s grid, taking each depth limit as a number or a
raster. Here a synthetic bathymetry deepens from 20 m in the west to 320
m in the east. Bathymetry such as GEBCO is elevation, negative below sea
level, so flip its sign and set land to 0 before use, then resample it
onto the voxel’s grid:

``` r

elevation <- terra::rast(nrows = 10, ncols = 10, xmin = 150, xmax = 151,
                         ymin = -21, ymax = -20, crs = "EPSG:4326")
elevation <- -(20 + 300 * (terra::init(elevation, "x") - 150))
seafloor <- terra::resample(terra::clamp(-elevation, lower = 0), v)

to_bed <- vect_to_envelope(mpa, v, depth_min = 0,
                           depth_max = list(200, seafloor))
volume(to_bed) # km³
#> [1] 307.1901
summarise_3d(mask(v, to_bed), statistic = c("mean", "min", "max"))
#>         time     mean   min   max   n
#> 1 2020-03-03 22.04723 16.05 28.75 360

near_bed <- vect_to_envelope(mpa, v, depth_min = list(seafloor - 20),
                             depth_max = seafloor)
summarise_3d(mask(v, near_bed), statistic = c("mean", "min", "max"))
#>         time     mean   min   max  n
#> 1 2020-03-03 18.65017 16.05 22.75 72
```

In each cell, the deepest `depth_min` and the shallowest `depth_max`
given apply. Use this instead of
[`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md)
when the range’s top or bottom depends on the seafloor or another
surface.

## Through time: habitat volume per day

[`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
on the proxy measures each time step in turn, in km³. Pass `fun` to
count only some values as habitat, here water between 24 and 26 °C:

``` r

volume(proxy, fun = function(x) x >= 24 & x <= 26)
#>         time volume_km3
#> 1 2020-03-01   490.6675
#> 2 2020-03-02   424.1851
#> 3 2020-03-03   357.7448
#> 4 2020-03-04   302.8561
#> 5 2020-03-05   259.5079
```

The habitat shrinks from day to day because the synthetic water warms
each day.

For the habitat inside the area only, keep the habitat in one day’s
voxel with
[`occupied()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/occupied.md).
Masking it by the area would keep whole depth levels, so the 50 m level
would count water down to 75 m, below the area’s floor. Instead, turn
the habitat into a depth range per cell with
[`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
and intersect it with the area, which cuts it at 50 m:

``` r

habitat <- occupied(v, function(x) x >= 24 & x <= 26)
inside <- intersect_3d(voxel_to_envelope(habitat), reef)
volume(inside) # km³
#> [1] 23.097
```

## Putting the rows together

The buffer loop ends with a list of data frames, one per buffer, each
starting with the buffer’s `feature_id`. `do.call(rbind, ...)` joins
them into one table, which you can then join back to the observations by
ID. In the result, `date` is the detection date and `time` is the time
step the values come from:

``` r

results <- do.call(rbind, buffer_rows)
merge(obs, results, by.x = "tag", by.y = "feature_id")
#>   tag    lon    lat depth       date       time     mean        sd  n
#> 1   A 150.21 -20.31    12 2020-03-01 2020-03-01 27.36252 0.5863953 24
#> 2   A 150.21 -20.31    12 2020-03-01 2020-03-02 27.56252 0.5863953 24
#> 3   B 150.52 -20.56    48 2020-03-02 2020-03-01 25.38426 0.7445443 16
#> 4   B 150.52 -20.56    48 2020-03-02 2020-03-02 25.58426 0.7445443 16
#> 5   B 150.52 -20.56    48 2020-03-02 2020-03-03 25.78426 0.7445443 16
#> 6   C 150.80 -20.82   140 2020-03-04 2020-03-03 18.32145 3.1053138 12
#> 7   C 150.80 -20.82   140 2020-03-04 2020-03-04 18.52145 3.1053138 12
#> 8   C 150.80 -20.82   140 2020-03-04 2020-03-05 18.72145 3.1053138 12
```
