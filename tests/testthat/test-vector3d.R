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

test_that("invalid vertices are rejected with a clear error", {
  obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35))
  cases <- list(
    list(x = transform(obs, depth = c(10, NA)),
         why = "every vertex needs a depth"),
    list(x = transform(obs, depth = c(-5, 10)), why = "positive metres"),
    list(x = transform(obs, lat = c(-20, -95)), why = "latitudes within")
  )
  for (case in cases) {
    err <- expect_error(as_vector3d(case$x), class = "ocean3d_error_vector3d")
    expect_match(conditionMessage(err), case$why, fixed = TRUE)
  }

  one_vertex <- dive()[1, ]
  expect_error(
    as_vector3d(one_vertex, type = "lines", id = "tag", time = "when"),
    "at least two vertices", class = "ocean3d_error_vector3d"
  )
  expect_error(as_vector3d(dive(), type = "lines"),
               "Lines need `id`", class = "ocean3d_error_vector3d")
  expect_error(
    as_vector3d(transform(obs, id = c(1, 1)), id = "id"),
    "each point needs its own ID", class = "ocean3d_error_vector3d"
  )
  expect_error(as_vector3d(transform(obs, time = Sys.time())),
               "time = \"time\"", fixed = TRUE)
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

test_that("coordinates and depths must be numeric, not factors or text", {
  obs <- data.frame(lon = c(150, 151), lat = c(-20, -21))
  obs$depth <- factor(c("10", "35"))
  expect_error(as_vector3d(obs), "`depth` must be numeric, but it is factor",
               class = "ocean3d_error_vector3d")
  obs$depth <- c(10, 35)
  obs$lon <- c("150", "151")
  expect_error(as_vector3d(obs), "`lon` must be numeric, but it is character",
               class = "ocean3d_error_vector3d")
})

test_that("feature IDs must be present, and unique for sf lines", {
  d <- dive()
  d$tag[2] <- NA
  expect_error(as_vector3d(d, type = "lines", id = "tag"),
               "feature IDs must be non-missing",
               class = "ocean3d_error_vector3d")

  s <- dive_sf()
  s$tag <- c("A", "A")
  expect_error(as_vector3d(s, id = "tag"),
               "separate lines would be joined",
               class = "ocean3d_error_vector3d")
  s$tag <- c("A", NA)
  expect_error(as_vector3d(s, id = "tag"),
               "its own non-missing ID",
               class = "ocean3d_error_vector3d")
})

test_that("sf lines without Z coordinates need a depth", {
  flat <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(matrix(c(150, -20, 151, -21), ncol = 2, byrow = TRUE)),
    crs = 4326
  ))
  expect_error(as_vector3d(flat), "Z coordinates",
               class = "ocean3d_error_vector3d")
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

test_that("sf lines with empty geometries are rejected, not dropped", {
  s <- dive_sf()
  sf::st_geometry(s)[[2]] <- sf::st_linestring(matrix(numeric(0), ncol = 4),
                                               dim = "XYZM")
  expect_error(as_vector3d(s, id = "tag"),
               "Row 2 of `x` has an empty geometry",
               class = "ocean3d_error_vector3d", fixed = TRUE)
})

test_that("sf line attributes named like vertex columns stop, not vanish", {
  s <- dive_sf()
  s$time <- as.Date(c("2020-01-01", "2020-01-02"))
  expect_error(as_vector3d(s, id = "tag"),
               "Column names time are used by SpatVector3D",
               class = "ocean3d_error_vector3d", fixed = TRUE)
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

  # A Date start with a POSIXct end would otherwise lose the end's hours.
  a$from <- as.Date(c("2020-01-01", "2020-02-01"))
  expect_error(
    as_vector3d(a, depth_min = 0, depth_max = 50, time = c("from", "to")),
    "`from` and `to` must both be dates or both have hours",
    class = "ocean3d_error_vector3d", fixed = TRUE
  )
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
  a <- areas()
  expect_error(as_vector3d(a), "Polygons need `depth_min` and `depth_max`",
               class = "ocean3d_error_vector3d")
  expect_error(as_vector3d(a, depth_min = 0, depth_max = c(50, 60)),
               "must each be one number or the name of a column",
               class = "ocean3d_error_vector3d")
  expect_error(as_vector3d(a, depth_min = 300, depth_max = 200),
               "depth_min must not be deeper than depth_max",
               class = "ocean3d_error_vector3d")
  expect_error(as_vector3d(a, depth_min = -5, depth_max = 200),
               "positive metres", class = "ocean3d_error_vector3d")
  a$top_text <- c("0", "10")
  expect_error(as_vector3d(a, depth_min = "top_text", depth_max = 200),
               "`top_text` must be numeric, but it is character",
               class = "ocean3d_error_vector3d")

  a <- areas()
  a$name <- "reef"
  expect_error(as_vector3d(a, depth_min = 0, depth_max = 50, id = "name"),
               "each polygon needs its own ID",
               class = "ocean3d_error_vector3d")

  a <- areas()
  expect_error(as_vector3d(a, depth_min = 0, depth_max = 50,
                           time = c("to", "from")),
               "time_start must not be after time_end",
               class = "ocean3d_error_vector3d")
  expect_error(as_vector3d(a, depth_min = 0, depth_max = 50,
                           time = c("from", "to", "from")),
               "one column (a single time) or two", fixed = TRUE,
               class = "ocean3d_error_vector3d")

  a$depth_min <- 1
  expect_error(as_vector3d(a, depth_min = 0, depth_max = 50),
               "Column names depth_min are used by SpatVector3D",
               class = "ocean3d_error_vector3d", fixed = TRUE)

  a <- areas()
  sf::st_geometry(a)[[2]] <- sf::st_polygon()
  expect_error(as_vector3d(a, depth_min = 0, depth_max = 50),
               "Row 2 of `x` has an empty geometry", fixed = TRUE,
               class = "ocean3d_error_vector3d")
})

test_that("depth_min and depth_max are only for polygons", {
  expect_error(as_vector3d(dive(), depth_min = 0, depth_max = 10),
               "are for polygons", class = "ocean3d_error_vector3d")
})

test_that("volume() on polygons says it isn't supported yet", {
  p <- as_vector3d(areas(), depth_min = 0, depth_max = 50)
  expect_error(volume(p), "doesn't support SpatVector3D polygons yet",
               class = "ocean3d_error_volume_type")
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
  pts <- as_vector3d(data.frame(lon = 150, lat = -20, depth = 5))
  err <- "ocean3d_error_vector3d"
  expect_error(buffer_3d(data.frame(lon = 1, lat = 1, depth = 1), 1, 1),
               "must be SpatVector3D points or lines", class = err)
  polygons <- buffer_3d(pts, 1, 1)
  expect_error(buffer_3d(polygons, 1, 1), "already have a depth range",
               class = err)
  for (radius in list(0, -1, 10000, c(1, 2), "5", NA_real_)) {
    expect_error(buffer_3d(pts, radius, 1), "`radius_km` must be",
                 class = err)
  }
  for (window in list(-1, c(1, 2), NA_real_, "5")) {
    expect_error(buffer_3d(pts, 1, window), "`depth_window` must be",
                 class = err)
  }
  expect_error(buffer_3d(pts, 1, 1, time_window = 1),
               "needs a time at every vertex. Feature 1 of `x` has no time",
               class = err)
  expect_error(buffer_3d(pts, 1, 1, time_window = -1),
               "`time_window` must be", class = err)
  dated <- as_vector3d(data.frame(lon = 150, lat = -20, depth = 5,
                                  date = as.Date("2020-01-01")),
                       time = "date")
  expect_error(buffer_3d(dated, 1, 1, time_window = 0.5),
               "must be whole days", class = err)
  # A window computed in floating point that is a whole day to within
  # rounding counts as one.
  near_day <- buffer_3d(dated, 1, 1, time_window = 1 + 1e-12)
  expect_equal(as.data.frame(near_day)$time_end,
               as.POSIXct("2020-01-02", tz = "UTC"))

  empty <- as_vector3d(data.frame(lon = numeric(), lat = numeric(),
                                  depth = numeric()))
  expect_error(buffer_3d(empty, 1, 1), "has no points or lines", class = err)
})

test_that("columns named like polygon columns stop, not duplicate", {
  obs <- data.frame(lon = 150, lat = -20, depth = 5, geometry = "x")
  expect_error(buffer_3d(as_vector3d(obs), 5, 0),
               "Column names geometry are used by SpatVector3D",
               class = "ocean3d_error_vector3d")
})

test_that("a line timed at only some vertices stops", {
  d <- dive()
  d$when[2] <- NA
  tracks <- as_vector3d(d, type = "lines", id = "tag", time = "when")
  for (window in list(NULL, 1)) {
    expect_error(buffer_3d(tracks, 5, 0, time_window = window),
                 "Line A of `x` has a time at some vertices but not others",
                 class = "ocean3d_error_vector3d", fixed = TRUE)
  }
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
  err <- "ocean3d_error_envelope_input"
  empty <- "ocean3d_error_envelope_empty"
  pts <- as_vector3d(data.frame(lon = 151, lat = -20.5, depth = 5))
  expect_error(vector3d_to_envelope(pts, envelope_grid()),
               "These are points; buffer_3d() turns them into polygons",
               class = err, fixed = TRUE)
  expect_error(vector3d_to_envelope(areas(), envelope_grid()),
               "must be SpatVector3D polygons", class = err)
  p <- as_vector3d(areas(), depth_min = 0, depth_max = 200)
  expect_error(vector3d_to_envelope(p, matrix(1)),
               "`template` must be a SpatRaster", class = err)
  far <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 1, ymin = 0,
                     ymax = 1, crs = "EPSG:4326")
  expect_error(vector3d_to_envelope(p, far), "don't overlap any cell",
               class = empty)
  none <- p
  none@features <- p@features[0, ]
  expect_error(vector3d_to_envelope(none, envelope_grid()),
               "`x` has no polygons", class = err)
  no_crs <- terra::rast(nrows = 2, ncols = 2, xmin = 150, xmax = 153,
                        ymin = -21, ymax = -20, crs = "")
  expect_error(vector3d_to_envelope(p, no_crs), "`template` has no CRS",
               class = err)
  flat <- buffer_3d(pts, radius_km = 20, depth_window = 0)
  expect_error(vector3d_to_envelope(flat, envelope_grid()),
               "Every polygon has `depth_min` equal to `depth_max`.*buffer_3d",
               class = empty)
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
