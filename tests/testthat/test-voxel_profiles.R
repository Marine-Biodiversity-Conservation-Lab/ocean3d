# Vertical profiles, tested on their own. Each is handed the arguments
# envelope_to_voxel() would hand it and checked on the weights it returns, so a
# shape can be written or changed without reasoning about voxel assembly.

# profile_flat() ----

test_that("profile_flat() weights every level the same", {
  ind <- make_occupancy()
  w <- profile_flat(ind, c(0, 100, 200, 300), sum(ind))

  # a bare 1 broadcasts across cells and levels alike, leaving the value whole
  expect_identical(w, 1)
})

# profile_equal() ----

test_that("profile_equal() splits a cell's value over the levels it occupies", {
  ind <- make_occupancy()
  w <- profile_equal(ind, c(0, 100, 200, 300), sum(ind))

  expect_s4_class(w, "SpatRaster")
  expect_identical(terra::nlyr(w), 1)
  # cells occupy 1, 2 and 4 levels, so carry 1, 1/2 and 1/4 at each
  expect_identical(unname(terra::values(w)[1:3, 1]), c(1, 0.5, 0.25))
})

test_that("profile_equal() weights sum to 1 across each cell's occupied levels", {
  ind <- make_occupancy()
  w <- profile_equal(ind, c(0, 100, 200, 300), sum(ind))

  # the weight lands only where the cell is occupied, so this is the total the
  # cell's value is spread over — 1 wherever there is anywhere to put it
  totals <- terra::values(sum(ind * w))
  expect_identical(unname(totals[1:3, 1]), c(1, 1, 1))
})

test_that("profile_equal() gives NA, not Inf, where no level is occupied", {
  ind <- make_occupancy()
  w <- profile_equal(ind, c(0, 100, 200, 300), sum(ind))

  # cell 4 occupies nothing: dividing by its zero count must not reach 1/0
  expect_true(is.na(terra::values(w)[4, 1]))
  expect_false(any(is.infinite(terra::values(w))))
})

test_that("profile_equal() carries NA cells of the envelope through", {
  # a cell absent from the envelope is NA in the occupancy stack, and stays NA
  ind <- make_occupancy(cells = list(c(TRUE, TRUE, FALSE, FALSE),
                                     c(NA, NA, NA, NA),
                                     c(TRUE, TRUE, TRUE, TRUE),
                                     c(FALSE, FALSE, FALSE, FALSE)))
  w <- profile_equal(ind, c(0, 100, 200, 300), sum(ind))

  expect_true(is.na(terra::values(w)[2, 1]))
  expect_identical(unname(terra::values(w)[c(1, 3), 1]), c(0.5, 0.25))
})

test_that("profile_equal() does not depend on the depth values themselves", {
  # only the count of occupied levels matters, not how deep or how spaced
  ind <- make_occupancy()
  n <- sum(ind)

  near <- profile_equal(ind, c(0, 1, 2, 3), n)
  far <- profile_equal(ind, c(0, 500, 1000, 4000), n)

  expect_identical(terra::values(near), terra::values(far))
})

# custom profiles ----

test_that("envelope_to_voxel() accepts a profile written by the caller", {
  fp <- make_footprint(c(1, 1, 1, 1))
  envel <- as_envelope(fp, depth_min = 0, depth_max = 350)
  depths <- c(0, 100, 200, 300)

  # weight the top two levels twice as heavily as the deeper ones, normalised
  # so each cell still sums back to its value
  surface_weighted <- function(ind, depths, n_depths) {
    w <- terra::rast(lapply(depths, function(d) {
      terra::setValues(terra::rast(ind[[1]]), ifelse(d <= 100, 2, 1))
    }))
    w <- w * ind
    w / sum(w)
  }

  vals <- terra::values(
    envelope_to_voxel(envel, depths = depths, values = 6,
                      profile = surface_weighted)
  )

  # shares of 2:2:1:1 out of 6, and the cell total is preserved
  expect_identical(unname(vals[1, ]), c(2, 2, 1, 1))
  expect_identical(unname(rowSums(vals, na.rm = TRUE)), rep(6, 4))
})

test_that("a custom profile may name its arguments however it likes", {
  # the profile is called positionally, so these formals are not a contract
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 350)

  vals <- terra::values(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300), values = 8,
                      profile = function(occ, z, n) 1 / n)
  )

  expect_identical(unname(vals[1, ]), rep(2, 4))
})

test_that("a custom profile sees the depths it is being asked about", {
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 350)
  seen <- NULL

  envelope_to_voxel(envel, depths = c(0, 100, 200, 300),
                    profile = function(ind, depths, n_depths) {
                      seen <<- depths
                      1
                    })

  expect_identical(seen, c(0, 100, 200, 300))
})

test_that("a custom profile may return one weight per depth", {
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 350)
  depths <- c(0, 100, 200, 300)

  # a full stack: half the value at the surface layer, whole value below
  vals <- terra::values(
    envelope_to_voxel(envel, depths = depths, values = 10,
                      profile = function(ind, depths, n_depths) {
                        terra::rast(lapply(c(0.5, 1, 1, 1), function(k) {
                          terra::setValues(terra::rast(ind[[1]]), k)
                        }))
                      })
  )

  expect_identical(unname(vals[1, ]), c(5, 10, 10, 10))
})

# profile return values ----

test_that("envelope_to_voxel() rejects a profile returning the wrong layer count", {
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 350)

  expect_error(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300),
                      profile = function(ind, depths, n_depths) ind[[1:2]]),
    "`profile` must return a SpatRaster with 1 layer or one per depth (4); got 2",
    fixed = TRUE
  )
})

test_that("envelope_to_voxel() rejects a profile returning an off-grid raster", {
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 350)
  other <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3,
                       ymin = 0, ymax = 3, vals = 1)

  expect_error(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300),
                      profile = function(ind, depths, n_depths) other),
    "`profile` must return a SpatRaster on the same grid",
    fixed = TRUE
  )
})

test_that("envelope_to_voxel() rejects a profile returning a non-weight", {
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 350)
  depths <- c(0, 100, 200, 300)

  expect_error(
    envelope_to_voxel(envel, depths = depths,
                      profile = function(ind, depths, n_depths) "half"),
    "`profile` must return a single non-missing number or a SpatRaster",
    fixed = TRUE
  )
  expect_error(
    envelope_to_voxel(envel, depths = depths,
                      profile = function(ind, depths, n_depths) c(1, 2, 3)),
    "of length 3",
    fixed = TRUE
  )
  expect_error(
    envelope_to_voxel(envel, depths = depths,
                      profile = function(ind, depths, n_depths) NULL),
    "`profile` must return a single non-missing number or a SpatRaster",
    fixed = TRUE
  )
})

test_that("a failing profile is reported as the profile's failure", {
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 350)

  expect_error(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300),
                      profile = function(ind, depths, n_depths) stop("nope")),
    "`profile` failed when called with (ind, depths, n_depths): nope",
    fixed = TRUE
  )

  # a profile that cannot take the three arguments fails the same way, rather
  # than surfacing R's bare "unused argument" from inside the call
  expect_error(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300),
                      profile = function(ind) 1),
    "`profile` failed when called with (ind, depths, n_depths)",
    fixed = TRUE
  )
})
