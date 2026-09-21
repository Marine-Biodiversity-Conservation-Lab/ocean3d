# Synthetic GFW-shaped effort tibble. Mirrors the columns returned by
# gfwr::gfw_ais_fishing_hours() with group_by = "GEARTYPE" so the test fixtures
# look like the real input the function receives.
make_effort <- function(rows) {
  out <- data.frame(
    Lat = rows$lat,
    Lon = rows$lon,
    `Time Range` = rows$year %||% 2022L,
    geartype = rows$geartype,
    `Vessel IDs` = rows$vessels %||% 1L,
    `Apparent Fishing Hours` = rows$hours,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# A 5x5 study grid covering 0..5 in both dimensions, EPSG:4326. Cells are
# 1 deg square so longitude/latitude maps to integer cell indices cleanly.
make_grid <- function(res = 1, xmin = 0, xmax = 5, ymin = 0, ymax = 5,
                      crs = "EPSG:4326") {
  terra::rast(
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    resolution = res, crs = crs
  )
}

test_that("per-gear totals round-trip from effort table to raster", {
  # Create test data
  effort <- make_effort(list(
    lat      = c(0.5, 1.5, 0.5, 2.5, 3.5),
    lon      = c(0.5, 1.5, 2.5, 3.5, 4.5),
    geartype = c("trawlers", "trawlers", "drifting_longlines",
                 "drifting_longlines", "set_gillnets"),
    hours    = c(10,  5,  3,  7,  2)
  ))
  grid <- make_grid()

  out <- gfw_effort_to_raster(effort, grid)

  totals <- terra::global(out, "sum", na.rm = TRUE)$sum
  names(totals) <- names(out)

  expect_identical(totals[["effort_trawlers"]], 15)
  expect_identical(totals[["effort_drifting_longlines"]], 10)
  expect_identical(totals[["effort_set_gillnets"]], 2)
})

test_that("output has one layer per gear class, named effort_<level>", {
  # Create test data
  effort <- make_effort(list(
    lat      = c(0.5, 1.5, 2.5),
    lon      = c(0.5, 1.5, 2.5),
    geartype = c("trawlers", "drifting_longlines", "set_gillnets"),
    hours    = c(1, 1, 1)
  ))
  grid <- make_grid()

  out <- gfw_effort_to_raster(effort, grid)

  expect_identical(terra::nlyr(out), 3)
  expect_setequal(
    names(out),
    c("effort_trawlers", "effort_drifting_longlines", "effort_set_gillnets")
  )
})

test_that("multiple records in the same cell are aggregated by `fun`", {
  # Create test data
  # Three trawler records all fall inside cell (0.5, 0.5)
  effort <- make_effort(list(
    lat      = c(0.5, 0.6, 0.7),
    lon      = c(0.5, 0.6, 0.7),
    geartype = c("trawlers", "trawlers", "trawlers"),
    hours    = c(2, 3, 5)
  ))
  grid <- make_grid()

  summed <- gfw_effort_to_raster(effort, grid, fun = "sum")
  meaned <- gfw_effort_to_raster(effort, grid, fun = "mean")

  expect_identical(terra::global(summed, "sum", na.rm = TRUE)$sum, 10)
  expect_identical(
    terra::global(meaned, "max", na.rm = TRUE)$max,
    mean(c(2, 3, 5))
  )
})

test_that("layer_by = NULL produces a single total-effort raster named 'effort'", {
  # Create test data
  effort <- make_effort(list(
    lat      = c(0.5, 1.5, 2.5),
    lon      = c(0.5, 1.5, 2.5),
    geartype = c("trawlers", "drifting_longlines", "set_gillnets"),
    hours    = c(2, 3, 5)
  ))
  grid <- make_grid()

  out <- gfw_effort_to_raster(effort, grid, layer_by = NULL)

  expect_identical(terra::nlyr(out), 1)
  expect_identical(names(out), "effort")
  expect_identical(terra::global(out, "sum", na.rm = TRUE)$sum, 10)
})

test_that("missing required columns produce an informative error", {
  # Create test data
  bad <- data.frame(
    Lat = 0.5, Lon = 0.5, geartype = "trawlers",
    stringsAsFactors = FALSE
  )
  grid <- make_grid()

  expect_error(
    gfw_effort_to_raster(bad, grid),
    "Apparent Fishing Hours"
  )
})

test_that("custom `value` and `layer_by` columns are honoured", {
  # Create test data
  effort <- data.frame(
    Lat = c(0.5, 1.5),
    Lon = c(0.5, 1.5),
    flag = c("BGD", "IND"),
    hours_total = c(4, 6),
    stringsAsFactors = FALSE
  )
  grid <- make_grid()

  out <- gfw_effort_to_raster(
    effort, grid,
    layer_by = "flag",
    value = "hours_total"
  )

  expect_setequal(names(out), c("effort_BGD", "effort_IND"))
  totals <- terra::global(out, "sum", na.rm = TRUE)$sum
  names(totals) <- names(out)
  expect_identical(totals[["effort_BGD"]], 4)
  expect_identical(totals[["effort_IND"]], 6)
})

test_that("points in EPSG:4326 are reprojected onto a non-4326 grid", {
  # Same three trawler records as round-trip test, but rasterized onto a
  # Mollweide grid. The total fishing hours should be conserved across the
  # CRS change.
  effort <- make_effort(list(
    lat      = c(0.5, 1.5, 2.5),
    lon      = c(0.5, 1.5, 2.5),
    geartype = c("trawlers", "trawlers", "trawlers"),
    hours    = c(4, 5, 6)
  ))

  moll_grid <- terra::project(make_grid(), "+proj=moll")

  out <- gfw_effort_to_raster(effort, moll_grid)

  expect_identical(terra::crs(out), terra::crs(moll_grid))
  expect_identical(
    terra::global(out, "sum", na.rm = TRUE)$sum,
    15
  )
})

test_that("gfw_effort_to_raster assumes extent and resolution from input gfwr effort object", {
  effort <- make_effort(list(
    lat      = c(0.5, 1.5, 2.5),
    lon      = c(10.5, 11.5, 12.5),
    geartype = c("trawlers", "trawlers", "trawlers"),
    hours    = c(4, 5, 6)
  ))

  out <- gfw_effort_to_raster(effort)

  # Check SpatRaster
  expect_s4_class(out, "SpatRaster")
  # Check correct resolution
  expect_identical(res(out), c(1,1))
  # Check correct extent
  expect_true(terra::identical(ext(out), ext(10, 13, 0, 3))) # nolint: expect_identical_linter. expect_identical() doesn't work for terra objects, since they are pointers in memory
})

test_that("gfw_effort_to_raster cannot assume extent and resolution from irregular interval between effort values", {
  bad_effort <- make_effort(list(
    lat      = c(0.5, 1.6, 2.5),
    lon      = c(10.5, 11.6, 12.5),
    geartype = c("trawlers", "trawlers", "trawlers"),
    hours    = c(4, 5, 6)
  ))

  expect_error(gfw_effort_to_raster(bad_effort), "Unable to assume grid resolution from effort data frame. Please provide grid.")
})
