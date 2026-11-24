# SpatVector3D points, lines and polygons (#114).

dive <- function() {
  data.frame(
    tag = rep(c("A", "B"), each = 3),
    lon = c(150, 150.1, 150.2, 160, 160.1, 160.2), lat = -20,
    depth = c(0, 120, 15, 5, 50, 10),
    when = as.POSIXct("2020-01-01", tz = "UTC") + 3600 * c(0:2, 0:2),
    speed = 1:6
  )
}

# The dive() tracks as sf lines: Z is height (minus depth), M is time.
dive_sf <- function() {
  d <- dive()
  lines <- lapply(split(d, d$tag), function(rows) {
    sf::st_linestring(
      cbind(rows$lon, rows$lat, -rows$depth, as.numeric(rows$when)),
      dim = "XYZM"
    )
  })
  sf::st_sf(tag = names(lines), geometry = sf::st_sfc(lines, crs = 4326))
}

vertex_cols <- c("feature_id", "vertex", "lon", "lat", "depth", "time")

test_that("points keep one depth, time and attributes per point", {
  obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35),
                    date = as.Date(c("2020-01-01", "2020-01-02")),
                    species = "shark")
  p <- as_vector3d(obs, time = "date")

  expect_s4_class(p, "SpatVector3D")
  expect_identical(geomtype(p), "points")
  expect_identical(length(p), 2L)
  v <- as.data.frame(p)
  expect_equal(v$depth, c(10, 35))
  expect_equal(v$time, as.POSIXct(c("2020-01-01", "2020-01-02"), tz = "UTC"))
  expect_identical(v$species, c("shark", "shark"))
  expect_identical(p@time_resolution, "daily")
})

test_that("lines join each id's rows in row order", {
  d <- dive()
  shuffled <- d[c(1, 4, 2, 5, 3, 6), ]
  tr <- as_vector3d(shuffled, type = "lines", id = "tag", time = "when")

  expect_identical(geomtype(tr), "lines")
  expect_identical(length(tr), 2L)
  v <- as.data.frame(tr)
  expect_identical(v$feature_id, rep(c("A", "B"), each = 3))
  expect_identical(v$vertex, rep(1:3, 2))
  expect_equal(v$depth, d$depth)
  expect_identical(v$speed, 1:6)
  expect_identical(tr@time_resolution, "datetime")
})

test_that("invalid points and lines are rejected with a clear error", {
  # `args` go to as_vector3d(). Errors are ocean3d_error_vector3d unless
  # `unclassed` is TRUE; `fixed` matches `pattern` as a fixed string.
  obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35))
  cases <- list(
    list(label = "missing depth",
         args = list(transform(obs, depth = c(10, NA))),
         pattern = "every vertex needs a depth", fixed = TRUE),
    list(label = "negative depth",
         args = list(transform(obs, depth = c(-5, 10))),
         pattern = "positive metres", fixed = TRUE),
    list(label = "latitude past the pole",
         args = list(transform(obs, lat = c(-20, -95))),
         pattern = "latitudes within", fixed = TRUE),
    list(label = "factor depth",
         args = list(transform(obs, depth = factor(c("10", "35")))),
         pattern = "`depth` must be numeric, but it is factor"),
    list(label = "text longitude",
         args = list(transform(obs, lon = c("150", "151"))),
         pattern = "`lon` must be numeric, but it is character"),
    list(label = "line with one vertex",
         args = list(dive()[1, ], type = "lines", id = "tag", time = "when"),
         pattern = "at least two vertices"),
    list(label = "lines without `id`", args = list(dive(), type = "lines"),
         pattern = "Lines need `id`"),
    list(label = "missing line ID",
         args = list(local({
           d <- dive()
           d$tag[2] <- NA
           d
         }), type = "lines", id = "tag"),
         pattern = "feature IDs must be non-missing"),
    list(label = "repeated point IDs",
         args = list(transform(obs, id = c(1, 1)), id = "id"),
         pattern = "each point needs its own ID"),
    list(label = "time column without `time`",
         args = list(transform(obs, time = Sys.time())),
         pattern = "time = \"time\"", fixed = TRUE, unclassed = TRUE),
    list(label = "depth range for points",
         args = list(dive(), depth_min = 0, depth_max = 10),
         pattern = "are for polygons"),
    # IDs must be unique for sf lines.
    list(label = "sf lines sharing an ID",
         args = list(local({
           s <- dive_sf()
           s$tag <- c("A", "A")
           s
         }), id = "tag"),
         pattern = "separate lines would be joined"),
    list(label = "sf line without an ID",
         args = list(local({
           s <- dive_sf()
           s$tag <- c("A", NA)
           s
         }), id = "tag"),
         pattern = "its own non-missing ID"),
    list(label = "sf lines without Z coordinates",
         args = list(sf::st_sf(geometry = sf::st_sfc(
           sf::st_linestring(matrix(c(150, -20, 151, -21), ncol = 2,
                                    byrow = TRUE)),
           crs = 4326
         ))),
         pattern = "Z coordinates"),
    # Rejected, not dropped.
    list(label = "sf line with an empty geometry",
         args = list(local({
           s <- dive_sf()
           sf::st_geometry(s)[[2]] <- sf::st_linestring(
             matrix(numeric(0), ncol = 4), dim = "XYZM"
           )
           s
         }), id = "tag"),
         pattern = "Row 2 of `x` has an empty geometry", fixed = TRUE),
    # An attribute named like a vertex column stops, not vanishes.
    list(label = "sf line attribute named `time`",
         args = list(local({
           s <- dive_sf()
           s$time <- as.Date(c("2020-01-01", "2020-01-02"))
           s
         }), id = "tag"),
         pattern = "Column names time are used by SpatVector3D", fixed = TRUE)
  )
  for (case in cases) {
    class <- if (isTRUE(case$unclassed)) NULL else "ocean3d_error_vector3d"
    expect_error(do.call(as_vector3d, case$args), case$pattern, class = class,
                 fixed = isTRUE(case$fixed), label = case$label)
  }
})

test_that("sf lines take depth from minus Z and time from M", {
  from_sf <- as_vector3d(dive_sf(), id = "tag")
  from_table <- as_vector3d(dive(), type = "lines", id = "tag", time = "when")
  expect_identical(geomtype(from_sf), "lines")
  expect_equal(as.data.frame(from_sf)[vertex_cols],
               as.data.frame(from_table)[vertex_cols],
               ignore_attr = TRUE)
})

test_that("sf points take depth from a column", {
  obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35))
  p <- as_vector3d(sf::st_as_sf(obs, coords = c("lon", "lat"), crs = 4326))
  expect_identical(geomtype(p), "points")
  expect_equal(as.data.frame(p)[c("lon", "lat", "depth")], obs,
               ignore_attr = TRUE)
})

test_that("sf input in a geographic CRS other than WGS84 is transformed", {
  # EPSG:4807 (NTF Paris) is longitude/latitude, but in grads from the Paris
  # meridian, so its coordinates differ from WGS84's.
  paris <- sf::st_transform(dive_sf(), 4807)
  expect_true(sf::st_is_longlat(paris))
  expect_gt(max(abs(sf::st_coordinates(paris)[, "X"] - dive()$lon)), 1)
  lines <- as.data.frame(as_vector3d(paris, id = "tag"))
  expect_equal(lines$lon, dive()$lon, tolerance = 1e-6)
  expect_equal(lines$lat, dive()$lat, tolerance = 1e-6)

  obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35))
  points_paris <- sf::st_transform(
    sf::st_as_sf(obs, coords = c("lon", "lat"), crs = 4326), 4807
  )
  points <- as.data.frame(as_vector3d(points_paris))
  expect_equal(points$lon, obs$lon, tolerance = 1e-6)
  expect_equal(points$lat, obs$lat, tolerance = 1e-6)
})

test_that("points and lines have no volume but are 3D shapes", {
  tr <- as_vector3d(dive(), type = "lines", id = "tag")
  expect_true(methods::is(tr, "SpatVolume"))
  expect_error(volume(tr), "lines, which have no volume",
               class = "ocean3d_error_volume_type")
  expect_error(volume(tr), "buffer_3d() turns lines into polygons",
               class = "ocean3d_error_volume_type", fixed = TRUE)
})

# Polygons ------------------------------------------------------------------

# Two one-degree squares off Queensland, with depth ranges and time intervals.
areas <- function() {
  square <- function(x0) {
    sf::st_polygon(list(rbind(
      c(x0, -20), c(x0 + 1, -20), c(x0 + 1, -21), c(x0, -21), c(x0, -20)
    )))
  }
  sf::st_sf(
    name = c("reef", "bank"), top = c(0, 10), bottom = c(200, 50),
    from = as.Date(c("2020-01-01", "2020-02-01")),
    to = as.Date(c("2020-03-01", "2020-02-15")),
    status = c("no-take", "partial"),
    geometry = sf::st_sfc(square(150), square(152), crs = 4326)
  )
}

test_that("polygons take a depth range from numbers and keep attributes", {
  p <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  expect_s4_class(p, "SpatVector3D")
  expect_true(methods::is(p, "SpatVolume"))
  expect_identical(geomtype(p), "polygons")
  expect_identical(length(p), 2L)
  f <- as.data.frame(p)
  expect_identical(f$feature_id, c("1", "2"))
  expect_equal(f$depth_min, c(0, 0))
  expect_equal(f$depth_max, c(200, 200))
  expect_true(all(is.na(f$time_start)) && all(is.na(f$time_end)))
  expect_identical(f$status, c("no-take", "partial"))
  expect_s3_class(f$geometry, "sfc")
  expect_equal(sf::st_crs(f$geometry), sf::st_crs(4326))
  expect_identical(p@time_resolution, NA_character_)
})

test_that("polygons take depths, IDs and a time interval from columns", {
  p <- as_vector3d(areas(), depth_min = "top", depth_max = "bottom",
                   time = c("from", "to"), id = "name")
  f <- as.data.frame(p)
  expect_identical(f$feature_id, c("reef", "bank"))
  expect_equal(f$depth_min, c(0, 10))
  expect_equal(f$depth_max, c(200, 50))
  expect_equal(f$time_start,
               as.POSIXct(c("2020-01-01", "2020-02-01"), tz = "UTC"))
  expect_equal(f$time_end,
               as.POSIXct(c("2020-03-01", "2020-02-15"), tz = "UTC"))
  expect_identical(p@time_resolution, "daily")
  # Columns used as arguments aren't repeated as attributes.
  expect_false(any(c("name", "top", "bottom", "from", "to") %in% names(f)))
  expect_output(show(p), "polygons")
  expect_output(show(p), "0 to 200")
  expect_output(show(p), "2020-01-01 to 2020-03-01")
})

test_that("one time column gives each polygon a single time", {
  f <- as.data.frame(
    as_vector3d(areas(), depth_min = 0, depth_max = 50, time = "from")
  )
  expect_equal(f$time_start, f$time_end)
  expect_equal(f$time_start,
               as.POSIXct(c("2020-01-01", "2020-02-01"), tz = "UTC"))
})

test_that("polygon start and end times keep their hours and one resolution", {
  a <- areas()
  a$from <- as.POSIXct(c("2020-01-01 06:00", "2020-02-01 18:30"), tz = "UTC")
  a$to <- as.POSIXct(c("2020-01-01 09:15", "2020-02-02 00:00"), tz = "UTC")
  p <- as_vector3d(a, depth_min = 0, depth_max = 50, time = c("from", "to"))
  f <- as.data.frame(p)
  expect_equal(f$time_start, a$from)
  expect_equal(f$time_end, a$to)
  expect_identical(p@time_resolution, "datetime")
})

test_that("polygons come from SpatVector, MULTIPOLYGON and other CRSs", {
  from_sf <- as.data.frame(as_vector3d(areas(), depth_min = 0,
                                       depth_max = 50))
  from_terra <- as.data.frame(as_vector3d(terra::vect(areas()),
                                          depth_min = 0, depth_max = 50))
  expect_equal(from_terra$depth_max, from_sf$depth_max)
  expect_equal(as.numeric(sf::st_bbox(from_terra$geometry)),
               as.numeric(sf::st_bbox(from_sf$geometry)))

  multi <- sf::st_cast(areas(), "MULTIPOLYGON")
  expect_identical(geomtype(as_vector3d(multi, depth_min = 0,
                                        depth_max = 50)), "polygons")

  mercator <- sf::st_transform(areas(), 3857)
  back <- as.data.frame(as_vector3d(mercator, depth_min = 0, depth_max = 50))
  expect_equal(as.numeric(sf::st_bbox(back$geometry)),
               as.numeric(sf::st_bbox(from_sf$geometry)), tolerance = 1e-6)
})

test_that("invalid polygon input is rejected with a clear error", {
  # Each row's polygons are areas(), changed by `change` if it has one, and
  # `args` are the other arguments to as_vector3d(). `fixed` matches
  # `pattern` as a fixed string.
  depths <- list(depth_min = 0, depth_max = 50)
  cases <- list(
    list(label = "no depths", args = list(),
         pattern = "Polygons need `depth_min` and `depth_max`"),
    list(label = "two depth_max values",
         args = list(depth_min = 0, depth_max = c(50, 60)),
         pattern = "must each be one number, the name of a column"),
    list(label = "depth_min deeper than depth_max",
         args = list(depth_min = 300, depth_max = 200),
         pattern = "depth_min must not be deeper than depth_max"),
    list(label = "negative depth_min",
         args = list(depth_min = -5, depth_max = 200),
         pattern = "positive metres"),
    list(label = "text depth column",
         change = function(a) {
           a$top_text <- c("0", "10")
           a
         },
         args = list(depth_min = "top_text", depth_max = 200),
         pattern = "`top_text` must be numeric, but it is character"),
    list(label = "repeated polygon IDs",
         change = function(a) {
           a$name <- "reef"
           a
         },
         args = c(depths, list(id = "name")),
         pattern = "each polygon needs its own ID"),
    list(label = "start after end",
         args = c(depths, list(time = c("to", "from"))),
         pattern = "time_start must not be after time_end"),
    list(label = "three time columns",
         args = c(depths, list(time = c("from", "to", "from"))),
         pattern = "one column (a single time) or two", fixed = TRUE),
    # A Date start with a POSIXct end would otherwise lose the end's hours.
    list(label = "Date start with a POSIXct end",
         change = function(a) {
           a$to <- as.POSIXct(c("2020-01-01 09:15", "2020-02-02 00:00"),
                              tz = "UTC")
           a
         },
         args = c(depths, list(time = c("from", "to"))),
         pattern = "`from` and `to` must both be dates or both have hours",
         fixed = TRUE),
    list(label = "attribute named depth_min",
         change = function(a) {
           a$depth_min <- 1
           a
         },
         args = depths,
         pattern = "Column names depth_min are used by SpatVector3D",
         fixed = TRUE),
    list(label = "empty geometry",
         change = function(a) {
           sf::st_geometry(a)[[2]] <- sf::st_polygon()
           a
         },
         args = depths, pattern = "Row 2 of `x` has an empty geometry",
         fixed = TRUE)
  )
  for (case in cases) {
    a <- if (is.null(case$change)) areas() else case$change(areas())
    expect_error(do.call(as_vector3d, c(list(a), case$args)), case$pattern,
                 class = "ocean3d_error_vector3d", fixed = isTRUE(case$fixed),
                 label = case$label)
  }
})

# Subsetting, extent and CRS -----------------------------------------------

test_that("`[` picks features by position, logical or ID, in that order", {
  obs <- data.frame(lon = c(150, 151, 152), lat = c(-20, -21, -22),
                    depth = c(10, 35, 50), species = c("a", "b", "c"))
  pts <- as_vector3d(obs)
  lines <- as_vector3d(dive(), type = "lines", id = "tag")
  polygons <- as_vector3d(areas(), id = "name", depth_min = "top",
                          depth_max = "bottom", time = c("from", "to"))

  # The same features whichever way they are picked, in the order asked,
  # with their IDs and other columns.
  picks <- list(c(3, 1), c(TRUE, FALSE, TRUE), c("3", "1"), -2)
  for (i in picks) {
    sub <- pts[i]
    expect_s4_class(sub, "SpatVector3D")
    expect_identical(geomtype(sub), "points")
    want <- if (is.logical(i) || is.numeric(i) && any(i < 0)) {
      c("1", "3")
    } else {
      c("3", "1")
    }
    expect_identical(as.data.frame(sub)$feature_id, want)
    expect_identical(as.data.frame(sub)$species, letters[as.integer(want)])
  }

  # Lines keep every vertex of each chosen line, in line then vertex order.
  b_then_a <- as.data.frame(lines[c("B", "A")])
  expect_identical(b_then_a$feature_id, rep(c("B", "A"), each = 3))
  expect_identical(b_then_a$vertex, rep(1:3, 2))
  expect_identical(b_then_a$speed, c(4:6, 1:3))
  expect_identical(length(lines[1]), 1L)

  # Polygons keep their outline, depth range, time and columns.
  bank <- polygons["bank"]
  expect_identical(as.data.frame(bank)$feature_id, "bank")
  expect_identical(as.data.frame(bank)$status, "partial")
  expect_identical(bank@time_resolution, "daily")
  expect_equal(volume(bank), volume(polygons) - volume(polygons[1]))

  # No feature is a valid, empty result.
  none <- pts[0]
  expect_identical(length(none), 0L)
  expect_s4_class(none, "SpatVector3D")
  expect_identical(length(polygons[c(FALSE, FALSE)]), 0L)

  # A feature can be picked once, whatever the geometry type.
  bad <- list(
    list(call = quote(pts[1, 2]), pattern = "by feature only"),
    list(call = quote(pts["z"]), pattern = "No feature has ID z"),
    list(call = quote(pts[4]), pattern = "whole numbers from 1 to 3"),
    list(call = quote(pts[c(-1, 2)]), pattern = "all positive or all negative"),
    list(call = quote(pts[c(TRUE, FALSE)]), pattern = "one TRUE or FALSE"),
    list(call = quote(pts[c(TRUE, NA, FALSE)]), pattern = "1 of them NA"),
    list(call = quote(pts[c(1, 1, 2)]), pattern = "repeats feature 1"),
    list(call = quote(lines[c("A", "A")]), pattern = "repeats feature A"),
    list(call = quote(polygons[c(2, 2)]), pattern = "repeats feature bank"),
    list(call = quote(pts[list(1)]), pattern = "got list"),
    list(call = quote(none[1]), pattern = "no features to pick from")
  )
  for (case in bad) {
    expect_error(eval(case$call), case$pattern,
                 class = "ocean3d_error_vector3d")
  }
})

test_that("ext() and crs() give the horizontal extent in WGS84", {
  pts <- as_vector3d(data.frame(lon = c(150, 152), lat = c(-22, -20),
                                depth = c(10, 35)))
  expect_equal(as.vector(terra::ext(pts)),
               c(xmin = 150, xmax = 152, ymin = -22, ymax = -20))
  lines <- as_vector3d(dive(), type = "lines", id = "tag")
  expect_equal(as.vector(terra::ext(lines)),
               c(xmin = 150, xmax = 160.2, ymin = -20, ymax = -20))
  polygons <- as_vector3d(areas(), depth_min = "top", depth_max = "bottom")
  expect_equal(as.vector(terra::ext(polygons)),
               c(xmin = 150, xmax = 153, ymin = -21, ymax = -20))
  expect_true(terra::same.crs(terra::crs(pts), "EPSG:4326"))
  expect_identical(terra::crs(pts, describe = TRUE)$code, "4326")
  expect_error(terra::ext(pts[0]), "no features, so it has no extent",
               class = "ocean3d_error_vector3d")
})

test_that("ext() stays continuous across the antimeridian and round a pole", {
  lon_range <- function(x) as.vector(terra::ext(x))[c("xmin", "xmax")]
  # A line from 179 to -179 goes the short way, through 180.
  seam <- as_vector3d(data.frame(id = "A", lon = c(179, -179), lat = 0,
                                 depth = 10), type = "lines", id = "id")
  expect_equal(lon_range(seam), c(xmin = 179, xmax = 181))
  # A buffer split at the antimeridian, and points on either side of it.
  buffer <- buffer_3d(seam, radius_km = 20, depth_window = 5)
  got <- lon_range(buffer)
  expect_true(got[["xmin"]] > 178.7 && got[["xmin"]] < 179 &&
                got[["xmax"]] > 181 && got[["xmax"]] < 181.3)
  points <- as_vector3d(data.frame(lon = c(170, -170, 175), lat = 0,
                                   depth = 10))
  expect_equal(lon_range(points), c(xmin = 170, xmax = 190))
  # Features that leave no gap round the globe cover every longitude.
  spread <- as_vector3d(data.frame(id = "A", lon = c(-170, -50, 70, 190),
                                   lat = 0, depth = 10),
                        type = "lines", id = "id")
  expect_equal(lon_range(spread), c(xmin = -180, xmax = 180))
  # A buffer around the North Pole covers every longitude up to the pole.
  pole <- buffer_3d(as_vector3d(data.frame(lon = 0, lat = 89.9, depth = 10)),
                    radius_km = 50, depth_window = 5)
  e <- as.vector(terra::ext(pole))
  expect_equal(e[c("xmin", "xmax", "ymax")],
               c(xmin = -180, xmax = 180, ymax = 90))
})

# Buffers -------------------------------------------------------------------

area_km2 <- function(polygons) {
  as.numeric(sf::st_area(as.data.frame(polygons)$geometry)) / 1e6
}

test_that("a point becomes a circle with the depth and time window", {
  obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(35, 10),
                    date = as.Date(c("2020-01-01", "2020-01-02")),
                    species = "shark")
  b <- buffer_3d(as_vector3d(obs, time = "date"), radius_km = 10,
                 depth_window = 25, time_window = 1)

  expect_identical(geomtype(b), "polygons")
  f <- as.data.frame(b)
  expect_identical(f$feature_id, c("1", "2"))
  expect_equal(area_km2(b), rep(pi * 10^2, 2), tolerance = 0.01)
  centres <- sf::st_sfc(sf::st_point(c(150, -20)), sf::st_point(c(151, -21)),
                        crs = 4326)
  expect_true(all(sf::st_within(centres, f$geometry, sparse = FALSE)[
    cbind(1:2, 1:2)
  ]))
  # The second top would be -15 m, so it stops at the surface.
  expect_equal(f$depth_min, c(10, 0))
  expect_equal(f$depth_max, c(60, 35))
  expect_equal(f$time_start,
               as.POSIXct(c("2019-12-31", "2020-01-01"), tz = "UTC"))
  expect_equal(f$time_end,
               as.POSIXct(c("2020-01-02", "2020-01-03"), tz = "UTC"))
  expect_identical(b@time_resolution, "daily")
  expect_identical(f$species, c("shark", "shark"))
})

test_that("without a time window, polygons keep the features' own times", {
  obs <- data.frame(lon = 150, lat = -20, depth = 35,
                    when = as.POSIXct("2020-01-01 06:00", tz = "UTC"))
  f <- as.data.frame(buffer_3d(as_vector3d(obs, time = "when"), 1, 0))
  expect_equal(f$time_start, obs$when)
  expect_equal(f$time_end, obs$when)
  expect_equal(c(f$depth_min, f$depth_max), c(35, 35))

  untimed <- buffer_3d(as_vector3d(obs[1:3]), 1, 5)
  expect_true(is.na(as.data.frame(untimed)$time_start))
  expect_true(is.na(untimed@time_resolution))
})

test_that("a line's polygon covers its whole depth and time range", {
  tracks <- as_vector3d(dive(), type = "lines", id = "tag", time = "when")
  b <- buffer_3d(tracks, radius_km = 5, depth_window = 10,
                 time_window = as.difftime(30, units = "mins"))
  f <- as.data.frame(b)

  expect_identical(f$feature_id, c("A", "B"))
  expect_equal(f$depth_min, c(0, 0))
  expect_equal(f$depth_max, c(130, 60))
  t0 <- as.POSIXct("2020-01-01", tz = "UTC")
  expect_equal(f$time_start, rep(t0 - 1800, 2))
  expect_equal(f$time_end, rep(t0 + 2 * 3600 + 1800, 2))
  # Area within 5 km of a line: about 2 r L + pi r^2.
  length_km <- as.numeric(sf::st_length(sf::st_sfc(
    sf::st_linestring(cbind(c(150, 150.1, 150.2), -20)), crs = 4326
  ))) / 1000
  expect_equal(area_km2(b)[1], 2 * 5 * length_km + pi * 25, tolerance = 0.01)
  # speed changes along each line, so it has no one value per polygon.
  expect_false("speed" %in% names(f))
})

test_that("line attributes with one value per line are kept", {
  d <- dive()
  d$ship <- rep(c("X", "Y"), each = 3)
  f <- as.data.frame(buffer_3d(as_vector3d(d, type = "lines", id = "tag"),
                               5, 10))
  expect_identical(f$ship, c("X", "Y"))
})

test_that("buffers across the antimeridian go the short way", {
  short <- data.frame(id = "a", lon = c(179.9, -179.9), lat = 0, depth = 0)
  b <- buffer_3d(as_vector3d(short, type = "lines", id = "id"), 5, 0)
  length_km <- 0.2 * 111.32
  expect_equal(area_km2(b), 2 * 5 * length_km + pi * 25, tolerance = 0.02)

  point <- data.frame(lon = 180, lat = 0, depth = 0)
  circle <- buffer_3d(as_vector3d(point), 10, 0)
  expect_equal(area_km2(circle), pi * 100, tolerance = 0.01)
  expect_true(all(sf::st_is_valid(as.data.frame(circle)$geometry)))
})

test_that("buffer_3d() rejects inputs it can't buffer", {
  # `args` go to buffer_3d(); `fixed` matches `pattern` as a fixed string.
  pts <- as_vector3d(data.frame(lon = 150, lat = -20, depth = 5))
  dated <- as_vector3d(data.frame(lon = 150, lat = -20, depth = 5,
                                  date = as.Date("2020-01-01")),
                       time = "date")
  partly_timed <- local({
    d <- dive()
    d$when[2] <- NA
    as_vector3d(d, type = "lines", id = "tag", time = "when")
  })
  partly_timed_msg <- "Line A of `x` has a time at some vertices but not others"
  cases <- c(
    list(
      list(label = "a data frame",
           args = list(data.frame(lon = 1, lat = 1, depth = 1), 1, 1),
           pattern = "must be SpatVector3D points or lines"),
      list(label = "polygons", args = list(buffer_3d(pts, 1, 1), 1, 1),
           pattern = "already have a depth range")
    ),
    lapply(list(0, -1, 10000, c(1, 2), "5", NA_real_), function(radius) {
      list(label = paste("radius_km", deparse(radius)),
           args = list(pts, radius, 1), pattern = "`radius_km` must be")
    }),
    lapply(list(-1, c(1, 2), NA_real_, "5"), function(window) {
      list(label = paste("depth_window", deparse(window)),
           args = list(pts, 1, window), pattern = "`depth_window` must be")
    }),
    list(
      list(label = "time_window for untimed points",
           args = list(pts, 1, 1, time_window = 1),
           pattern = paste("needs a time at every vertex.",
                           "Feature 1 of `x` has no time")),
      list(label = "negative time_window",
           args = list(pts, 1, 1, time_window = -1),
           pattern = "`time_window` must be"),
      list(label = "part-day time_window for dates",
           args = list(dated, 1, 1, time_window = 0.5),
           pattern = "must be whole days"),
      list(label = "no features",
           args = list(as_vector3d(data.frame(lon = numeric(), lat = numeric(),
                                              depth = numeric())), 1, 1),
           pattern = "has no points or lines"),
      # Columns named like polygon columns stop, not duplicate.
      list(label = "attribute named geometry",
           args = list(as_vector3d(data.frame(lon = 150, lat = -20, depth = 5,
                                              geometry = "x")), 5, 0),
           pattern = "Column names geometry are used by SpatVector3D"),
      list(label = "line timed at some vertices, no time_window",
           args = list(partly_timed, 5, 0, time_window = NULL),
           pattern = partly_timed_msg, fixed = TRUE),
      list(label = "line timed at some vertices, with a time_window",
           args = list(partly_timed, 5, 0, time_window = 1),
           pattern = partly_timed_msg, fixed = TRUE)
    )
  )
  for (case in cases) {
    expect_error(do.call(buffer_3d, case$args), case$pattern,
                 class = "ocean3d_error_vector3d", fixed = isTRUE(case$fixed),
                 label = case$label)
  }

  # A window computed in floating point that is a whole day to within
  # rounding counts as one.
  near_day <- buffer_3d(dated, 1, 1, time_window = 1 + 1e-12)
  expect_equal(as.data.frame(near_day)$time_end,
               as.POSIXct("2020-01-02", tz = "UTC"))
})

test_that("vertices stored out of order are buffered in line order", {
  # Line A is a V, so joining its vertices in another order changes the shape.
  d <- dive()
  d$lat[2] <- -20.5
  tracks <- as_vector3d(d, type = "lines", id = "tag")
  shuffled <- tracks
  shuffled@vertices <- tracks@vertices[c(3, 1, 2, 6, 4, 5), ]
  expect_equal(as.data.frame(buffer_3d(shuffled, 5, 10)),
               as.data.frame(buffer_3d(tracks, 5, 10)))
})

# Envelopes ----------------------------------------------------------------

# A 0.25-degree grid over both areas() squares and the gap between them.
envelope_grid <- function() {
  terra::rast(nrows = 4, ncols = 12, xmin = 150, xmax = 153, ymin = -21,
              ymax = -20, crs = "EPSG:4326")
}

test_that("polygons with one depth range match vect_to_envelope()", {
  p <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  e <- vector3d_to_envelope(p, envelope_grid())
  expected <- vect_to_envelope(terra::vect(areas()), envelope_grid(),
                               depth_min = 0, depth_max = 200)
  expect_s4_class(e, "SpatEnvelope")
  expect_equal(terra::values(e), terra::values(expected),
               ignore_attr = TRUE)
})

test_that("each cell gets the polygons' own depths, outer range on overlap", {
  p <- as_vector3d(areas(), depth_min = "top", depth_max = "bottom")
  v <- terra::values(vector3d_to_envelope(p, envelope_grid()))
  expect_setequal(stats::na.omit(v[, "depth_min"]), c(0, 10))
  expect_setequal(stats::na.omit(v[, "depth_max"]), c(200, 50))

  # Two squares over the same cells, 0 to 50 m and 200 to 300 m.
  stacked <- areas()
  stacked$geometry[2] <- stacked$geometry[1]
  stacked$top <- c(0, 200)
  stacked$bottom <- c(50, 300)
  s <- as_vector3d(stacked, depth_min = "top", depth_max = "bottom")
  v <- stats::na.omit(terra::values(vector3d_to_envelope(s, envelope_grid())))
  expect_true(nrow(v) > 0)
  expect_true(all(v[, "depth_min"] == 0 & v[, "depth_max"] == 300))
})

test_that("polygon times are dropped with a message", {
  timed <- as_vector3d(areas(), depth_min = 0, depth_max = 200,
                       time = c("from", "to"))
  expect_message(vector3d_to_envelope(timed, envelope_grid()),
                 "times are dropped", class = "ocean3d_message_time_dropped")
  untimed <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  expect_no_message(vector3d_to_envelope(untimed, envelope_grid()))
})

test_that("polygons are projected onto a template in another CRS", {
  p <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  grid <- terra::project(envelope_grid(), "EPSG:3857")
  e <- vector3d_to_envelope(p, grid)
  expected <- vect_to_envelope(terra::project(terra::vect(areas()), grid),
                               grid, depth_min = 0, depth_max = 200)
  expect_equal(terra::values(e), terra::values(expected),
               ignore_attr = TRUE)
})

test_that("a buffered detection becomes an envelope with its depth window", {
  obs <- data.frame(lon = 151.5, lat = -20.5, depth = 35)
  b <- buffer_3d(as_vector3d(obs), radius_km = 20, depth_window = 25)
  v <- stats::na.omit(terra::values(vector3d_to_envelope(b, envelope_grid())))
  expect_true(nrow(v) > 0)
  expect_true(all(v[, "depth_min"] == 10 & v[, "depth_max"] == 60))
})

test_that("vector3d_to_envelope() rejects inputs it can't convert", {
  # `x` and `template` go to vector3d_to_envelope(); `fixed` matches
  # `pattern` as a fixed string.
  input <- "ocean3d_error_envelope_input"
  empty <- "ocean3d_error_envelope_empty"
  pts <- as_vector3d(data.frame(lon = 151, lat = -20.5, depth = 5))
  p <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  cases <- list(
    list(label = "points", x = pts, template = envelope_grid(),
         class = "ocean3d_error_vector3d_type",
         pattern = "these are points. buffer_3d() turns points into polygons",
         fixed = TRUE),
    list(label = "sf polygons", x = areas(), template = envelope_grid(),
         class = input, pattern = "must be SpatVector3D polygons"),
    list(label = "matrix template", x = p, template = matrix(1),
         class = input, pattern = "`template` must be a SpatRaster"),
    list(label = "template away from the polygons", x = p,
         template = terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 1,
                                ymin = 0, ymax = 1, crs = "EPSG:4326"),
         class = empty, pattern = "don't overlap any cell"),
    list(label = "no polygons",
         x = local({
           none <- p
           none@features <- p@features[0, ]
           none
         }),
         template = envelope_grid(), class = input,
         pattern = "`x` has no polygons"),
    list(label = "template without a CRS", x = p,
         template = terra::rast(nrows = 2, ncols = 2, xmin = 150, xmax = 153,
                                ymin = -21, ymax = -20, crs = ""),
         class = input, pattern = "`template` has no CRS"),
    list(label = "polygons with no thickness",
         x = buffer_3d(pts, radius_km = 20, depth_window = 0),
         template = envelope_grid(), class = empty,
         pattern = paste0("Every polygon has `depth_min` equal to ",
                          "`depth_max`.*buffer_3d"))
  )
  for (case in cases) {
    expect_error(vector3d_to_envelope(case$x, case$template), case$pattern,
                 class = case$class, fixed = isTRUE(case$fixed),
                 label = case$label)
  }
})

test_that("cells with no thickness are dropped", {
  flat <- areas()
  flat$top <- c(30, 10)
  flat$bottom <- c(30, 50)
  p <- as_vector3d(flat, depth_min = "top", depth_max = "bottom")
  v <- stats::na.omit(terra::values(vector3d_to_envelope(p, envelope_grid())))
  expect_true(all(v[, "depth_min"] == 10 & v[, "depth_max"] == 50))

  # The flat square under the thick one keeps the thick one's range.
  flat$geometry[1] <- flat$geometry[2]
  p <- as_vector3d(flat, depth_min = "top", depth_max = "bottom")
  v <- stats::na.omit(terra::values(vector3d_to_envelope(p, envelope_grid())))
  expect_true(all(v[, "depth_min"] == 10 & v[, "depth_max"] == 50))
})

test_that("envelope depths keep their decimals", {
  p <- as_vector3d(areas(), depth_min = 10.1, depth_max = 35.1)
  v <- stats::na.omit(terra::values(vector3d_to_envelope(p, envelope_grid())))
  expect_identical(unique(v[, "depth_min"]), 10.1)
  expect_identical(unique(v[, "depth_max"]), 35.1)
})

test_that("polygons land on grids with longitudes past 180", {
  # A circle across the antimeridian is split in two; both halves count.
  b <- buffer_3d(as_vector3d(data.frame(lon = 180, lat = 0, depth = 10)),
                 radius_km = 300, depth_window = 5)
  b0 <- buffer_3d(as_vector3d(data.frame(lon = 0, lat = 0, depth = 10)),
                  radius_km = 300, depth_window = 5)
  cells <- function(polygons, xmin) {
    grid <- terra::rast(nrows = 10, ncols = 40, xmin = xmin,
                        xmax = xmin + 20, ymin = -5, ymax = 5,
                        crs = "EPSG:4326")
    e <- vector3d_to_envelope(polygons, grid)
    sum(!is.na(terra::values(e)[, "depth_min"]))
  }
  expected <- cells(b0, -10)
  expect_gt(expected, 0)
  expect_identical(cells(b, 170), expected)
  expect_identical(cells(b, -190), expected)

  # An outline stored past 180 lands on a -180 to 180 grid too.
  past <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(178, -1), c(182, -1), c(182, 1), c(178, 1), c(178, -1)
  ))), crs = 4326))
  p <- as_vector3d(past, depth_min = 0, depth_max = 10)
  west <- terra::rast(nrows = 2, ncols = 10, xmin = -180, xmax = -170,
                      ymin = -1, ymax = 1, crs = "EPSG:4326")
  # It covers -180 to -178 on this 1-degree grid: 2 columns by 2 rows.
  v <- terra::values(vector3d_to_envelope(p, west))[, "depth_min"]
  expect_identical(sum(!is.na(v)), 4L)
})

# ---- vector3d_to_voxel ------------------------------------------------------

# A voxel on envelope_grid() with levels 0, 100, 200 and 300 m.
voxel_grid <- function(depth_bounds = NULL) {
  r <- terra::rast(envelope_grid(), nlyrs = 4, vals = 1)
  as_voxel(r, depths = c(0, 100, 200, 300), varname = "temp",
           depth_bounds = depth_bounds)
}

# Two squares over the same cells, 0 to 50 m and 200 to 300 m.
stacked_bands <- function() {
  stacked <- areas()
  stacked$geometry[2] <- stacked$geometry[1]
  stacked$top <- c(0, 200)
  stacked$bottom <- c(50, 300)
  as_vector3d(stacked, depth_min = "top", depth_max = "bottom")
}

test_that("a depth gap between polygons over one cell stays empty", {
  v <- vector3d_to_voxel(stacked_bands(), voxel_grid())
  expect_s4_class(v, "SpatVoxel")
  expect_identical(depths(v), c(0, 100, 200, 300))
  expect_identical(names(v), paste0("presence_depth=", c(0, 100, 200, 300)))
  vals <- terra::values(v)
  covered <- !is.na(vals[, 1])
  expect_true(any(covered))
  # Midpoint slabs: 0-50, 50-150, 150-250 and 250-300 m. The 100 m level
  # lies between the bands.
  expect_true(all(vals[covered, 1] == 1))
  expect_true(all(is.na(vals[covered, 2])))
  expect_true(all(vals[covered, 3:4] == 1))
  expect_true(all(is.na(vals[!covered, ])))

  # An envelope can't hold the gap, so the same polygons fill the 100 m level.
  filled <- envelope_to_voxel(vector3d_to_envelope(stacked_bands(),
                                                   voxel_grid()),
                              depths(voxel_grid()))
  expect_true(all(terra::values(filled)[covered, 2] == 1))
})

test_that("one depth range lands on the same levels as an envelope", {
  p <- as_vector3d(areas(), depth_min = "top", depth_max = "bottom")
  for (bounds in c("midpoint", "top")) {
    v <- vector3d_to_voxel(p, voxel_grid(), bounds = bounds)
    expected <- envelope_to_voxel(vector3d_to_envelope(p, voxel_grid()),
                                  depths(voxel_grid()), bounds = bounds)
    expect_equal(terra::values(v), terra::values(expected),
                 ignore_attr = TRUE)
  }
})

test_that("the template's depth edges are used before `bounds`", {
  # Slabs 0-20, 20-180, 180-220 and 220-400 m: the 100 m level now reaches
  # the 0-50 m band, and the 200 m level only the 200-300 m band.
  edged <- voxel_grid(depth_bounds = cbind(c(0, 20, 180, 220),
                                           c(20, 180, 220, 400)))
  v <- vector3d_to_voxel(stacked_bands(), edged, bounds = "top")
  vals <- terra::values(v)
  covered <- !is.na(vals[, 1])
  expect_true(all(vals[covered, 1:4] == 1))
  expect_identical(depth_bounds(v), depth_bounds(edged))

  shallow <- as_vector3d(areas(), depth_min = 0, depth_max = 10)
  vals <- terra::values(vector3d_to_voxel(shallow, edged))
  expect_true(all(is.na(vals[, 2:4])))
})

test_that("polygons land on a voxel in another CRS or past 180 degrees", {
  p <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  # Projecting drops terra's depths, so they are set again.
  projected <- as_voxel(terra::project(voxel_grid(), "EPSG:3857"),
                        depths = depths(voxel_grid()))
  v <- vector3d_to_voxel(p, projected)
  expected <- vect_to_envelope(terra::project(terra::vect(areas()), projected),
                               projected, depth_min = 0, depth_max = 200)
  expect_identical(is.na(terra::values(v)[, 1]),
                   is.na(terra::values(expected)[, 1]))

  past <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(178, -1), c(182, -1), c(182, 1), c(178, 1), c(178, -1)
  ))), crs = 4326))
  p <- as_vector3d(past, depth_min = 0, depth_max = 10)
  west <- terra::rast(nrows = 2, ncols = 10, xmin = -180, xmax = -170,
                      ymin = -1, ymax = 1, crs = "EPSG:4326", nlyrs = 2,
                      vals = 1)
  west <- as_voxel(west, depths = c(0, 50))
  # It covers -180 to -178 on this 1-degree grid: 2 columns by 2 rows.
  v <- terra::values(vector3d_to_voxel(p, west))
  expect_identical(sum(!is.na(v[, 1])), 4L)
})

test_that("polygon times are dropped with a message for a voxel too", {
  timed <- as_vector3d(areas(), depth_min = 0, depth_max = 200,
                       time = c("from", "to"))
  expect_message(vector3d_to_voxel(timed, voxel_grid()),
                 "times are dropped", class = "ocean3d_message_time_dropped")
  untimed <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  expect_no_message(vector3d_to_voxel(untimed, voxel_grid()))
})

test_that("vector3d_to_voxel() rejects inputs it can't convert", {
  input <- "ocean3d_error_voxel_input"
  empty <- "ocean3d_error_envelope_empty"
  pts <- as_vector3d(data.frame(lon = 151, lat = -20.5, depth = 5))
  p <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  none <- p
  none@features <- p@features[0, ]
  far <- as_voxel(terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 1,
                              ymin = 0, ymax = 1, crs = "EPSG:4326",
                              nlyrs = 2, vals = 1),
                  depths = c(0, 50))
  cases <- list(
    list(label = "points", x = pts, template = voxel_grid(),
         class = "ocean3d_error_vector3d_type",
         pattern = "these are points. buffer_3d() turns points into polygons",
         fixed = TRUE),
    list(label = "sf polygons", x = areas(), template = voxel_grid(),
         class = input, pattern = "must be SpatVector3D polygons"),
    list(label = "envelope template", x = p, template = envelope_50_200(),
         class = "ocean3d_error_not_3d",
         pattern = "`template` must be a SpatVoxel"),
    list(label = "plain raster template", x = p, template = envelope_grid(),
         class = "ocean3d_error_not_3d",
         pattern = "`template` must be a SpatVoxel"),
    list(label = "no polygons", x = none, template = voxel_grid(),
         class = input, pattern = "`x` has no polygons"),
    list(label = "template without a CRS", x = p,
         template = as_voxel(terra::rast(nrows = 2, ncols = 2, xmin = 150,
                                         xmax = 153, ymin = -21, ymax = -20,
                                         crs = "", nlyrs = 2, vals = 1),
                             depths = c(0, 50)),
         class = "ocean3d_error_envelope_input",
         pattern = "`template` has no CRS"),
    list(label = "template away from the polygons", x = p, template = far,
         class = empty, pattern = "cover no cell and depth level"),
    list(label = "polygons with no thickness",
         x = buffer_3d(pts, radius_km = 20, depth_window = 0),
         template = voxel_grid(), class = empty,
         pattern = "cover no cell and depth level"),
    list(label = "polygons below every level",
         x = as_vector3d(areas(), depth_min = 500, depth_max = 600),
         template = voxel_grid(), class = empty,
         pattern = "cover no cell and depth level")
  )
  for (case in cases) {
    expect_error(vector3d_to_voxel(case$x, case$template), case$pattern,
                 class = case$class, fixed = isTRUE(case$fixed),
                 label = case$label)
  }
})

# ---- Raster depth surfaces (#117) --------------------------------------------

# A seafloor deepening across make_template()'s 5 x 5 grid, from 100 to 340 m.
seafloor_5x5 <- function() {
  terra::setValues(make_template(), seq(100, 340, by = 10))
}

test_that("a seafloor surface gives each cell its own depth on a grid", {
  seafloor <- seafloor_5x5()
  p <- as_vector3d(make_polygon(), depth_min = 0, depth_max = seafloor)
  expect_equal(terra::values(vector3d_to_envelope(p, make_template())),
               terra::values(vect_to_envelope(make_polygon(), make_template(),
                                              depth_min = 0,
                                              depth_max = seafloor)))
  # The column holds the deepest value over the polygon's cells.
  touched <- terra::values(vect_to_envelope(make_polygon(), make_template(),
                                            0, seafloor))[, 2]
  expect_equal(as.data.frame(p)$depth_max, max(touched, na.rm = TRUE))
  expect_length(p[1]@surfaces, 1L)
})

test_that("a constant surface gives the same results as the number", {
  constant <- terra::setValues(make_template(), 200)
  surface <- as_vector3d(make_polygon(), depth_min = 0, depth_max = constant)
  number <- as_vector3d(make_polygon(), depth_min = 0, depth_max = 200)
  expect_equal(terra::values(vector3d_to_envelope(surface, make_template())),
               terra::values(vector3d_to_envelope(number, make_template())))
  levels <- as_voxel(terra::rast(make_template(), nlyrs = 4, vals = 1),
                     depths = c(0, 100, 200, 300))
  expect_equal(terra::values(vector3d_to_voxel(surface, levels)),
               terra::values(vector3d_to_voxel(number, levels)))
  # Sampled on the surface's cells against the geodesic outline.
  expect_equal(volume(surface), volume(number), tolerance = 1e-3)
  # Saved and read back, the surface still works.
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(surface, file)
  expect_equal(volume(readRDS(file)), volume(surface))
})

test_that("a layer above the seafloor follows the bed", {
  seafloor <- seafloor_5x5()
  near_bed <- as_vector3d(make_polygon(), depth_min = seafloor - 20,
                          depth_max = seafloor)
  e <- vector3d_to_envelope(near_bed, make_template())
  thickness <- terra::values(e$depth_max - e$depth_min)[, 1]
  expect_equal(unique(stats::na.omit(thickness)), 20)
  expect_equal(volume(near_bed),
               volume(as_vector3d(make_polygon(), depth_min = 0,
                                  depth_max = 20)),
               tolerance = 1e-3)
  levels <- as_voxel(terra::rast(make_template(), nlyrs = 4, vals = 1),
                     depths = c(0, 100, 200, 300))
  # The deepest cells reach the deepest level, as the bed deepens.
  cover <- terra::values(mask(levels, near_bed))
  expect_gt(sum(!is.na(cover[, 4])), 0)
  expect_true(all(is.na(cover[, 1])))
})

test_that("cells without water are left out, as in vect_to_envelope()", {
  # From 200 m to the bed: the shelf cells shallower than 200 m drop out.
  seafloor <- seafloor_5x5()
  shelf <- as_vector3d(make_polygon(), depth_min = 200, depth_max = seafloor)
  expect_equal(terra::values(vector3d_to_envelope(shelf, make_template())),
               terra::values(vect_to_envelope(make_polygon(), make_template(),
                                              200, seafloor)))
  expect_equal(as.data.frame(shelf)$depth_min, 200)
  # Land (a negative depth) adds no volume.
  land <- terra::setValues(make_template(), c(rep(-5, 12), rep(100, 13)))
  coast <- as_vector3d(make_polygon(), depth_min = 0, depth_max = land)
  ocean <- terra::setValues(make_template(), c(rep(NA, 12), rep(100, 13)))
  expect_equal(volume(coast),
               volume(as_vector3d(make_polygon(), depth_min = 0,
                                  depth_max = ocean)))
})

test_that("a surface is rejected where it can't give a depth range", {
  seafloor <- seafloor_5x5()
  far <- terra::rast(nrows = 2, ncols = 2, xmin = 20, xmax = 22, ymin = 20,
                     ymax = 22, crs = "EPSG:4326", vals = 100)
  bare <- terra::rast(matrix(100, 5, 5))
  terra::crs(bare) <- ""
  cases <- list(
    list(depth_min = 400, depth_max = seafloor,
         pattern = "Polygon 1 has no water between"),
    list(depth_min = 0, depth_max = far,
         pattern = "Polygon 1 has no water between"),
    list(depth_min = 0, depth_max = bare,
         pattern = "`depth_max` has no CRS")
  )
  for (case in cases) {
    expect_error(as_vector3d(make_polygon(), depth_min = case$depth_min,
                             depth_max = case$depth_max),
                 case$pattern, class = "ocean3d_error_vector3d")
  }
  expect_error(as_vector3d(make_polygon(), depth_min = 0,
                           depth_max = c(seafloor, seafloor)),
               "single-layer SpatRaster", class = "ocean3d_error_grid_value")
  coarse <- terra::aggregate(seafloor, 5)
  expect_error(as_vector3d(make_polygon(), depth_min = coarse,
                           depth_max = seafloor),
               "does not align with", class = "ocean3d_error_grid_mismatch")
})

test_that("only the surface under the polygons is kept", {
  big <- terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 50, ymin = 0,
                     ymax = 50, crs = "EPSG:4326", vals = 300)
  small <- as_vector3d(make_polygon(), depth_min = 0, depth_max = big)
  kept <- terra::unwrap(small@surfaces$depth_max)
  expect_lt(terra::ncell(kept), terra::ncell(big) / 10)
  expect_equal(volume(small),
               volume(as_vector3d(make_polygon(), depth_min = 0,
                                  depth_max = 300)),
               tolerance = 1e-3)
})

test_that("a surface on another grid keeps its decimals", {
  fine <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 5, ymin = 0,
                      ymax = 5, crs = "EPSG:4326",
                      vals = seq(100.1, 595.1, by = 5))
  p <- as_vector3d(make_polygon(), depth_min = 0, depth_max = fine)
  finer <- terra::disagg(fine, 2)
  bottom <- stats::na.omit(terra::values(vector3d_to_envelope(p, finer))[, 2])
  # Bilinear on a finer grid: every value is a multiple of 0.05 m.
  expect_true(all(abs(bottom * 20 - round(bottom * 20)) < 1e-9))
  # A coarser grid averages the cells under each of its cells.
  coarse_bottom <- terra::values(vector3d_to_envelope(p, make_template()))[, 2]
  averaged <- terra::values(terra::aggregate(fine, 2, "mean"))[, 1]
  expect_equal(coarse_bottom[!is.na(coarse_bottom)],
               averaged[!is.na(coarse_bottom)])
})

test_that("a surface reaches a voxel cropped across the antimeridian", {
  # A one-degree band around the globe, cropped across 180 degrees: the
  # voxel runs from 179 to 181, past the surface's longitudes.
  band <- terra::rast(nrows = 4, ncols = 360, xmin = -180, xmax = 180,
                      ymin = -2, ymax = 2, crs = "EPSG:4326", nlyrs = 2,
                      vals = 1)
  levels <- as_voxel(band, depths = c(0, 100))
  seafloor <- terra::setValues(terra::rast(band[[1]]), 150)
  buffer <- buffer_3d(as_vector3d(data.frame(lon = 180, lat = 0, depth = 10)),
                      radius_km = 50, depth_window = 5)
  outline <- sf::st_sf(geometry = as.data.frame(buffer)$geometry)
  shape <- as_vector3d(outline, depth_min = 0, depth_max = seafloor)
  cropped <- as_voxel(crop(as_proxy(levels), shape))
  expect_equal(as.vector(terra::ext(cropped))[c("xmin", "xmax")],
               c(xmin = 179, xmax = 181))
  plain <- as_vector3d(outline, depth_min = 0, depth_max = 150)
  expect_equal(terra::values(mask(cropped, shape)),
               terra::values(mask(cropped, plain)))
})

test_that("a surface outside a voxel's levels covers nothing quietly", {
  levels <- as_voxel(terra::rast(make_template(), nlyrs = 4, vals = 1),
                     depths = c(0, 100, 200, 300))
  deep <- as_vector3d(make_polygon(), depth_min = 500,
                      depth_max = terra::setValues(make_template(), 600))
  expect_error(vector3d_to_voxel(deep, levels),
               "cover no cell and depth level",
               class = "ocean3d_error_envelope_empty")
  expect_no_warning(masked <- mask(levels, deep))
  expect_true(all(is.na(terra::values(masked))))
})

test_that("polygons with a surface aren't compared with other shapes", {
  surface <- as_vector3d(make_polygon(), depth_min = 0,
                         depth_max = seafloor_5x5())
  points <- as_vector3d(data.frame(lon = 2.5, lat = 2.5, depth = 10))
  for (y in list(polygons_0_100(), points)) {
    expect_error(intersect_3d(surface, y), "work with grids only",
                 class = "ocean3d_error_vector3d")
    expect_error(intersects_3d(y, surface), "work with grids only",
                 class = "ocean3d_error_vector3d")
  }
})
