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

# The study grid is make_template() (helper-example_objects.R): 5x5 over
# 0..5 in both dimensions, EPSG:4326. Cells are 1 deg square so
# longitude/latitude maps to integer cell indices cleanly.

test_that("layer totals round-trip from effort table to raster", {
  # `args` are the other arguments to gfw_effort_to_raster(); `totals` are the
  # summed effort in each layer, named by layer.
  cases <- list(
    list(label = "per gear",
         effort = make_effort(list(
           lat      = c(0.5, 1.5, 0.5, 2.5, 3.5),
           lon      = c(0.5, 1.5, 2.5, 3.5, 4.5),
           geartype = c("trawlers", "trawlers", "drifting_longlines",
                        "drifting_longlines", "set_gillnets"),
           hours    = c(10,  5,  3,  7,  2)
         )),
         args = list(),
         totals = c(effort_trawlers = 15, effort_drifting_longlines = 10,
                    effort_set_gillnets = 2)),
    # Custom `value` and `layer_by` columns are honoured.
    list(label = "custom value and layer_by columns",
         effort = data.frame(
           Lat = c(0.5, 1.5),
           Lon = c(0.5, 1.5),
           flag = c("BGD", "IND"),
           hours_total = c(4, 6),
           stringsAsFactors = FALSE
         ),
         args = list(layer_by = "flag", value = "hours_total"),
         totals = c(effort_BGD = 4, effort_IND = 6))
  )
  for (case in cases) {
    out <- do.call(gfw_effort_to_raster,
                   c(list(case$effort, make_template()), case$args))

    # expect_setequal() has no `label`; the same names in any order.
    expect_identical(sort(names(out)), sort(names(case$totals)),
                     label = case$label)
    totals <- terra::global(out, "sum", na.rm = TRUE)$sum
    names(totals) <- names(out)
    for (layer in names(case$totals)) {
      expect_identical(totals[[layer]], case$totals[[layer]],
                       label = paste0(case$label, ": ", layer))
    }
  }
})

test_that("output has one layer per gear class, named effort_<level>", {
  # Create test data
  effort <- make_effort(list(
    lat      = c(0.5, 1.5, 2.5),
    lon      = c(0.5, 1.5, 2.5),
    geartype = c("trawlers", "drifting_longlines", "set_gillnets"),
    hours    = c(1, 1, 1)
  ))
  grid <- make_template()

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
  grid <- make_template()

  summed <- gfw_effort_to_raster(effort, grid, fun = "sum")
  meaned <- gfw_effort_to_raster(effort, grid, fun = "mean")

  expect_identical(terra::global(summed, "sum", na.rm = TRUE)$sum, 10)
  expect_identical(
    terra::global(meaned, "max", na.rm = TRUE)$max,
    mean(c(2, 3, 5))
  )
})

test_that("layer_by = NULL gives one total-effort layer named 'effort'", {
  # Create test data
  effort <- make_effort(list(
    lat      = c(0.5, 1.5, 2.5),
    lon      = c(0.5, 1.5, 2.5),
    geartype = c("trawlers", "drifting_longlines", "set_gillnets"),
    hours    = c(2, 3, 5)
  ))
  grid <- make_template()

  out <- gfw_effort_to_raster(effort, grid, layer_by = NULL)

  expect_identical(terra::nlyr(out), 1)
  expect_named(out, "effort")
  expect_identical(terra::global(out, "sum", na.rm = TRUE)$sum, 10)
})

test_that("missing required columns produce an informative error", {
  # Create test data
  bad <- data.frame(
    Lat = 0.5, Lon = 0.5, geartype = "trawlers",
    stringsAsFactors = FALSE
  )
  grid <- make_template()

  expect_error(
    gfw_effort_to_raster(bad, grid),
    "Apparent Fishing Hours"
  )
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

  moll_grid <- terra::project(make_template(), "+proj=moll")

  out <- gfw_effort_to_raster(effort, moll_grid)

  expect_identical(terra::crs(out), terra::crs(moll_grid))
  expect_identical(
    terra::global(out, "sum", na.rm = TRUE)$sum,
    15
  )
})

test_that("gfw_effort_to_raster() infers extent and resolution if no grid", {
  # both are assumed from the input gfwr effort object itself
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
  # expect_identical() doesn't work for terra objects, since they are pointers
  # in memory
  # nolint start: expect_identical_linter.
  expect_true(terra::identical(ext(out), ext(10, 13, 0, 3)))
  # nolint end
})

test_that("gfw_effort_to_raster() errors when effort spacing is irregular", {
  # an irregular interval between effort values leaves no extent and
  # resolution to assume, so a grid must be provided
  bad_effort <- make_effort(list(
    lat      = c(0.5, 1.6, 2.5),
    lon      = c(10.5, 11.6, 12.5),
    geartype = c("trawlers", "trawlers", "trawlers"),
    hours    = c(4, 5, 6)
  ))

  expect_error(
    gfw_effort_to_raster(bad_effort),
    paste0("Unable to assume grid resolution from effort data frame. ",
           "Please provide grid.")
  )
})
