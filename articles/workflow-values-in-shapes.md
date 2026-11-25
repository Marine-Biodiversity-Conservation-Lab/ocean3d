# Environmental values inside points, buffers and areas

This article shows how to get environmental values inside 3D shapes,
such as tag detections, buffers around them, or a protected area down to
a set depth. Each example uses the same three steps:

1.  **Build** the shapes with
    [`as_vector3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_vector3d.md)
    or
    [`buffer_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/buffer_3d.md).
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
#> source     : file22a9108c62f.nc
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
skips depth levels with no value when it looks for the nearest one.
[`mask()`](https://rspatial.github.io/terra/reference/mask.html) uses
the nearest level whatever its value, so such a detection gets `NA`.

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
