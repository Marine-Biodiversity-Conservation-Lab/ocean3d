# SpatVector3D points and lines (#114).

dive <- function() {
  data.frame(
    tag = rep(c("A", "B"), each = 3),
    lon = c(150, 150.1, 150.2, 160, 160.1, 160.2), lat = -20,
    depth = c(0, 120, 15, 5, 50, 10),
    when = as.POSIXct("2020-01-01", tz = "UTC") + 3600 * c(0:2, 0:2),
    speed = 1:6
  )
}

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
    expect_error(as_vector3d(case$x), class = "ocean3d_error_vector3d")
    err <- tryCatch(as_vector3d(case$x),
                    error = function(e) conditionMessage(e))
    expect_match(err, case$why, fixed = TRUE)
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
    "Each point needs its own ID", class = "ocean3d_error_vector3d"
  )
  expect_error(as_vector3d(transform(obs, time = Sys.time())),
               "time = \"time\"", fixed = TRUE)
})

test_that("points and lines have no volume but are 3D shapes", {
  tr <- as_vector3d(dive(), type = "lines", id = "tag")
  expect_true(methods::is(tr, "SpatVolume"))
  expect_error(volume(tr), "lines, which have no volume",
               class = "ocean3d_error_volume_type")
})
