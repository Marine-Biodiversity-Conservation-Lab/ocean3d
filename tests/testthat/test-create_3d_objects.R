test_that("SpatVoxel and SpatEnvelope are SpatRasters and SpatVolumes", {
  v <- as_voxel(make_multidepth_rast())

  expect_s4_class(v, "SpatVoxel")
  expect_s4_class(v, "SpatRaster")
  expect_s4_class(v, "SpatVolume")
  # terra operations still work on the subclass
  expect_identical(terra::nlyr(v), 4)
})

test_that("new() refuses an object that breaks the class rules", {
  unsorted <- make_multidepth_rast()
  terra::depth(unsorted) <- c(300, 200, 100, 0)
  # Names in the old pattern, but no terra::depth().
  unset <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2, vals = 1)
  names(unset) <- c("t_depth=0", "t_depth=100")
  repeated <- make_multidepth_rast()
  terra::depth(repeated) <- c(0, 100, 100, 300)
  no_values <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2)
  names(no_values) <- c("depth_min", "depth_max")

  cases <- list(
    list(label = "voxel without depths", class = "SpatVoxel",
         object = unset, pattern = "must be set in terra::depth"),
    list(label = "voxel with a repeated depth", class = "SpatVoxel",
         object = repeated, pattern = "must be distinct"),
    list(label = "voxel with unsorted layers", class = "SpatVoxel",
         object = unsorted, pattern = "shallow to deep"),
    list(label = "envelope without the two depth layers",
         class = "SpatEnvelope", object = make_multidepth_rast(),
         pattern = "exactly: depth_min, depth_max"),
    list(label = "envelope without cell values", class = "SpatEnvelope",
         object = no_values, pattern = "must have values")
  )
  for (case in cases) {
    expect_error(methods::new(case$class, case$object), case$pattern,
                 label = case$label)
  }
})

test_that("a voxel's layer names are only labels", {
  # Only terra::depth() counts, whatever the names.
  free <- make_multidepth_rast()
  names(free) <- c("a", "b", "c", "d")
  expect_true(methods::validObject(methods::new("SpatVoxel", free)))
})

# as_envelope() ----

test_that("as_envelope() builds a SpatEnvelope on the footprint's grid", {
  # each depth limit is a constant or a single-layer raster on the footprint's
  # grid, and only the footprint's cells carry it
  fp <- make_footprint()

  cases <- list(
    list(label = "constant limits", fp = fp, dmin = 0, dmax = 200,
         min = c(0, 0, 0, NA), max = c(200, 200, 200, NA)),
    list(label = "per-cell limits", fp = fp,
         dmin = make_footprint(c(0, 50, 100, 150)),
         dmax = make_footprint(c(100, 200, 300, 400)),
         min = c(0, 50, 100, NA), max = c(100, 200, 300, NA)),
    list(label = "a constant and a per-cell limit", fp = fp, dmin = 10,
         dmax = make_footprint(c(100, 200, 300, 400)),
         min = c(10, 10, 10, NA), max = c(100, 200, 300, NA)),
    list(label = "empty footprint", fp = make_footprint(c(NA, NA, NA, NA)),
         dmin = 0, dmax = 200, min = rep(NA_real_, 4), max = rep(NA_real_, 4))
  )
  for (case in cases) {
    e <- as_envelope(case$fp, depth_min = case$dmin, depth_max = case$dmax)

    expect_true(methods::is(e, "SpatEnvelope"), label = case$label)
    expect_true(methods::is(e, "SpatVolume"), label = case$label)
    expect_named(e, c("depth_min", "depth_max"), label = case$label)
    expect_true(terra::compareGeom(e, case$fp, stopOnError = FALSE),
                label = case$label)
    vals <- terra::values(e)
    expect_identical(unname(vals[, "depth_min"]), case$min, label = case$label)
    expect_identical(unname(vals[, "depth_max"]), case$max, label = case$label)
  }
})

test_that("as_envelope() reads only the non-NA pattern of the footprint", {
  # Values differ wildly, including a zero and a negative; presence does not.
  a <- as_envelope(make_footprint(c(1, 5, 0, NA)), 10, 20)
  b <- as_envelope(make_footprint(c(-7, 0.5, 999, NA)), 10, 20)

  expect_identical(terra::values(a), terra::values(b))
})

test_that("as_envelope() promotes a raster that already has the depth layers", {
  fp <- make_footprint()
  plain <- c(
    terra::setValues(terra::rast(fp), c(0, 0, 0, NA)),
    terra::setValues(terra::rast(fp), c(200, 200, 200, NA))
  )
  names(plain) <- c("depth_min", "depth_max")

  e <- as_envelope(plain)

  expect_s4_class(e, "SpatEnvelope")
  expect_identical(terra::values(e), terra::values(plain))
})

test_that("as_envelope() is idempotent on its own output", {
  e <- as_envelope(make_footprint(), 0, 200)

  expect_identical(terra::values(as_envelope(e)), terra::values(e))
})

test_that("as_envelope() says how a per-cell limit is off the grid", {
  fp <- make_footprint()
  other_crs <- terra::rast(fp, vals = 100)
  terra::crs(other_crs) <- "EPSG:3857"
  cases <- list(
    list(limit = c(terra::rast(fp, vals = 100), terra::rast(fp, vals = 200)),
         pattern = "`depth_max` must be a single-layer SpatRaster; got 2"),
    list(limit = other_crs,
         pattern = "`depth_max` has a different CRS than `x`"),
    list(limit = terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 2,
                             ymin = 0, ymax = 2, vals = 100),
         pattern = "`depth_max` does not align with `x`")
  )
  for (case in cases) {
    err <- expect_error(as_envelope(fp, 0, case$limit),
                        class = "ocean3d_error_grid_value")
    expect_match(conditionMessage(err), case$pattern, fixed = TRUE)
  }
})

test_that("as_envelope() rejects unusable input and says what is wrong", {
  fp <- make_footprint()
  e <- as_envelope(fp, 0, 200)
  inverted_dmax <- terra::setValues(terra::rast(fp), c(300, 300, 5, NA))
  poly <- terra::vect("POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))")

  cases <- list(
    list(label = "depth args for an envelope",
         args = list(e, depth_min = 0, depth_max = 100),
         pattern = "already has depth_min and depth_max"),
    list(label = "bare footprint without depths", args = list(fp),
         pattern = "are required"),
    list(label = "multi-layer footprint",
         args = list(make_multidepth_rast(), 0, 200),
         pattern = "single-layer footprint"),
    # names the function that does take polygons
    list(label = "SpatVector", args = list(poly, 0, 200),
         pattern = "vect_to_envelope"),
    list(label = "non-raster", args = list(1:10, 0, 200),
         pattern = "must be a SpatRaster"),
    list(label = "inverted constant limits",
         args = list(fp, depth_min = 200, depth_max = 0),
         pattern = "depth_max must be at least depth_min"),
    list(label = "inverted per-cell limit",
         args = list(fp, depth_min = 10, depth_max = inverted_dmax),
         pattern = "depth_max must be at least depth_min"),
    # positive-down convention
    list(label = "depth_min of -50",
         args = list(fp, depth_min = -50, depth_max = 200),
         pattern = "positive metres increasing downward"),
    list(label = "depth_min of -10",
         args = list(fp, depth_min = -10, depth_max = 200),
         pattern = "positive metres"),
    # limits that aren't scalars or rasters
    list(label = "vector limit", args = list(fp, c(0, 100), 200),
         pattern = "single non-missing number"),
    list(label = "NA limit", args = list(fp, NA, 200),
         pattern = "single non-missing number"),
    list(label = "character limit", args = list(fp, "0", 200),
         pattern = "single non-missing number")
  )
  for (case in cases) {
    expect_error(do.call(as_envelope, case$args), case$pattern,
                 label = case$label)
  }
})

# as_voxel() ----

test_that("as_voxel() wraps a conforming multi-depth SpatRaster", {
  v <- as_voxel(make_multidepth_rast())

  expect_s4_class(v, "SpatVoxel")
  expect_s4_class(v, "SpatVolume")
  expect_named(v, paste0("temp_depth=", c(0, 100, 200, 300)))
})

test_that("as_voxel() is idempotent", {
  v <- as_voxel(make_multidepth_rast())

  expect_s4_class(as_voxel(v), "SpatVoxel")
  expect_named(as_voxel(v), names(v))
})

test_that("as_voxel() sorts layers shallow to deep, values and all", {
  # new() rejects unsorted input; the constructor repairs it, and the 0 m layer
  # must still hold the values it had before the sort
  cases <- list(
    list(label = "four levels",
         x = make_multidepth_rast(depths = c(300, 200, 100, 0)),
         depths = c(0, 100, 200, 300), at_0 = c(NA, NA, 8, NA)),
    list(label = "two levels",
         x = make_multidepth_rast(depths = c(300, 0),
                                  vals = list(c(1, 2, 3, 4), c(5, 6, 7, 8))),
         depths = c(0, 300), at_0 = c(5, 6, 7, 8))
  )
  for (case in cases) {
    v <- as_voxel(case$x)

    expect_true(methods::is(v, "SpatVoxel"), label = case$label)
    expect_named(v, paste0("temp_depth=", case$depths), label = case$label)
    expect_identical(unname(terra::values(v[["temp_depth=0"]])[, 1]),
                     case$at_0, label = case$label)
  }
})

test_that("as_voxel() builds layer names from `depths`", {
  r <- make_multidepth_rast()
  names(r) <- c("a", "b", "c", "d")

  v <- as_voxel(r, depths = c(0, 100, 200, 300), varname = "temp")
  expect_named(v, paste0("temp_depth=", c(0, 100, 200, 300)))
})

test_that("as_voxel() accepts a list of single-depth SpatRasters", {
  r <- make_multidepth_rast()
  lst <- lapply(seq_len(terra::nlyr(r)), function(i) r[[i]])

  # Each raster has its depth, which terra keeps when it combines them.
  v <- as_voxel(lst)
  expect_s4_class(v, "SpatVoxel")
  expect_identical(depths(v), c(0, 100, 200, 300))

  # Rasters without depths need `depths`.
  bare <- lapply(1:2, function(i) terra::rast(nrows = 2, ncols = 2, vals = i))
  expect_error(as_voxel(bare), class = "ocean3d_error_no_depths")
  expect_identical(depths(as_voxel(bare, depths = c(0, 50))), c(0, 50))
})

test_that("as_voxel() reads depths from terra::depth(), never from names", {
  # terra::depth(), whatever the names.
  r <- make_multidepth_rast()
  terra::depth(r) <- c(5, 15, 25, 35)
  v <- as_voxel(r)
  expect_identical(depths(v), c(5, 15, 25, 35))
  expect_identical(names(v), names(r))

  # Names in the old `{variable}_depth={value}` pattern without
  # terra::depth(): an error that says how to fix it.
  named <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2, vals = 1)
  names(named) <- c("t_depth=0", "t_depth=100")
  err <- expect_error(as_voxel(named), class = "ocean3d_error_no_depths")
  expect_match(conditionMessage(err), "Pass `depths` to as_voxel\\(\\)")
  expect_match(conditionMessage(err), "Layer names aren't read as depths")

  # A single layer set to 0 m is set; this also fails if terra stops
  # flagging the depths it holds.
  surface <- terra::rast(nrows = 2, ncols = 2, vals = 1)
  terra::depth(surface) <- 0
  expect_identical(depths(as_voxel(surface)), 0)
})

test_that("as_voxel() leaves the caller's raster as it was", {
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3, vals = 1)
  names(r) <- c("a", "b", "c")
  v <- as_voxel(r, depths = c(0, 50, 100), varname = "t")
  expect_identical(depths(v), c(0, 50, 100))
  expect_identical(terra::depth(r), c(0, 0, 0))
  expect_identical(names(r), c("a", "b", "c"))
  expect_error(as_voxel(r), class = "ocean3d_error_no_depths")

  unsorted <- make_multidepth_rast(depths = c(200, 0, 100, 300))
  as_voxel(unsorted)
  expect_identical(terra::depth(unsorted), c(200, 0, 100, 300))
})

test_that("as_voxel() rejects unusable input and says how to fix it", {
  r <- make_multidepth_rast()

  cases <- list(
    # rejected rather than negated
    list(label = "negative depths",
         args = list(r, depths = c(0, -100, -200, -300)),
         pattern = "positive metres increasing downward"),
    list(label = "duplicate depths",
         args = list(r, depths = c(0, 100, 100, 300)),
         pattern = "duplicate depth"),
    list(label = "mismatched `depths` length",
         args = list(r, depths = c(0, 100)), pattern = "one value per layer",
         class = "ocean3d_error_voxel_input"),
    list(label = "text `depths`", args = list(r, depths = letters[1:4]),
         pattern = "must be numeric metres",
         class = "ocean3d_error_voxel_input"),
    list(label = "NA in `depths`", args = list(r, depths = c(0, NA, 200, 300)),
         pattern = "some are NA", class = "ocean3d_error_voxel_input"),
    list(label = "non-raster", args = list("not a raster"),
         pattern = "must be a SpatRaster")
  )
  for (case in cases) {
    err <- expect_error(do.call(as_voxel, case$args), case$pattern,
                        label = case$label)
    if (!is.null(case$class)) {
      expect_s3_class(err, case$class)
    }
  }
})

test_that("as_voxel() repairs a tagged SpatVoxel with unsorted layers", {
  # package functions advise that as_voxel() can be used to repair SpatVoxel
  # so it should be able to repair when possible from the input invalid
  # SpatVoxel
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = 3)
  terra::values(r) <- cbind(c(1, 1), c(2, 2), c(3, 3))
  v <- as_voxel(r, depths = c(0, 100, 200), varname = "presence")

  reversed <- v[[c(3, 1, 2)]]
  expect_s4_class(reversed, "SpatVoxel")
  expect_false(isTRUE(methods::validObject(reversed, test = TRUE)))
  expect_error(volume(reversed), "not a valid one")

  fixed <- as_voxel(reversed)
  expect_true(methods::validObject(fixed, test = TRUE))
  expect_named(fixed, names(v))
  expect_identical(terra::values(fixed), terra::values(v))
  expect_identical(volume(fixed), volume(v))
})

# voxel_to_envelope() ----
test_that("voxel_to_envelope() rejects a non-SpatVoxel or a non-function fun", {
  # a plain multi-depth SpatRaster is refused until wrapped by as_voxel()
  not_voxel <- make_multidepth_rast()
  expect_error(voxel_to_envelope(not_voxel), "must be a SpatVoxel")
  expect_s4_class(as_voxel(not_voxel), "SpatVoxel")

  v <- as_voxel(make_multidepth_rast())
  expect_error(voxel_to_envelope(v, "extent"), "must be a function")
})

test_that("voxel_to_envelope() bounds the depths where `fun` holds", {
  v <- as_voxel(make_multidepth_rast())

  cases <- list(
    list(label = "default: the non-NA vertical extent", v = v, fun = NULL,
         min = c(0, 100, 0, NA), max = c(200, 200, 300, NA)),
    # only cell 2 exceeds 15, at 100 m and 200 m
    list(label = "x > 15", v = v, fun = function(x) x > 15,
         min = c(NA, 100, NA, NA), max = c(NA, 200, NA, NA)),
    list(label = "a predicate that never holds", v = v,
         fun = function(x) x > 1e6,
         min = rep(NA_real_, 4), max = rep(NA_real_, 4)),
    # a one-layer voxel collapses to a zero-thickness envelope at that depth
    list(label = "single-depth voxel",
         v = as_voxel(make_multidepth_rast(depths = 50,
                                           vals = list(c(10, NA, 5, NA)))),
         fun = NULL, min = c(50, NA, 50, NA), max = c(50, NA, 50, NA))
  )
  for (case in cases) {
    e <- if (is.null(case$fun)) {
      voxel_to_envelope(case$v)
    } else {
      voxel_to_envelope(case$v, case$fun)
    }

    expect_true(methods::is(e, "SpatEnvelope"), label = case$label)
    expect_named(e, c("depth_min", "depth_max"), label = case$label)
    vals <- terra::values(e)
    expect_identical(unname(vals[, "depth_min"]), case$min, label = case$label)
    expect_identical(unname(vals[, "depth_max"]), case$max, label = case$label)
  }
})

# envelope_to_voxel() ----

# Local helper: runs one row of a table. An error that isn't a failed
# expectation would end the test without saying which row raised it, so the
# row's label is added to its message.
in_case <- function(label, code) {
  withCallingHandlers(code, error = function(e) {
    if (!inherits(e, "expectation")) {
      stop("in row \"", label, "\": ", conditionMessage(e), call. = FALSE)
    }
  })
}

test_that("envelope_to_voxel() rejects invalid input", {
  # a multi-depth raster is a plausible-looking but wrong input: it is the
  # voxel form of a 3D domain, not the envelope form
  bad_envel <- make_multidepth_rast()
  good_envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)
  fp <- make_footprint()

  good_depths <- c(0, 50, 100, 150, 200, 300, 600)
  profile_msg <- paste0("`profile` needs to be NULL or a function of ",
                        "(ind, depths, n_depths); got ")

  cases <- list(
    list(label = "`x` not an envelope",
         args = list(x = bad_envel, depths = good_depths),
         pattern = "`x` needs to be of `SpatEnvelope` class"),
    list(label = "non-numeric depths",
         args = list(x = good_envel, depths = c("hello", 100, "depths")),
         pattern = "`depths` needs to be array coercible to numeric type"),
    list(label = "empty depths",
         args = list(x = good_envel, depths = numeric(0)),
         pattern = "`depths` needs to be array coercible to numeric type"),
    # positive-down convention
    list(label = "negative depths",
         args = list(x = good_envel, depths = c(-100, 0, 100)),
         pattern = "positive metres increasing downward"),
    # A shape is passed as the function itself, so a name that used to select
    # one is now just the wrong type.
    list(label = "profile as a name",
         args = list(x = good_envel, depths = good_depths, profile = "equal"),
         pattern = paste0(profile_msg, "character")),
    list(label = "profile as a number",
         args = list(x = good_envel, depths = good_depths, profile = 100),
         pattern = paste0(profile_msg, "numeric")),
    list(label = "profile NA",
         args = list(x = good_envel, depths = good_depths, profile = NA),
         pattern = "`profile` needs to be NULL or a function"),
    list(label = "unknown bounds",
         args = list(x = good_envel, depths = c(0, 100), bounds = "gaussian"),
         pattern = "should be one of"),
    list(label = "numeric varname",
         args = list(x = good_envel, depths = good_depths, varname = 100),
         pattern = "`varname` needs to be a string"),
    list(label = "two varnames",
         args = list(x = good_envel, depths = good_depths,
                     varname = c("presence", "absence")),
         pattern = "`varname` needs to be a string"),
    list(label = "values on another grid",
         args = list(x = good_envel, depths = c(0, 100),
                     values = terra::rast(nrows = 3, ncols = 3)),
         pattern = "`values` does not align with `x`"),
    list(label = "multi-layer values",
         args = list(x = good_envel, depths = c(0, 100),
                     values = c(terra::rast(fp), terra::rast(fp))),
         pattern = "must be a single-layer SpatRaster"),
    list(label = "vector values",
         args = list(x = good_envel, depths = c(0, 100), values = c(1, 2)),
         pattern = "single non-missing number")
  )
  for (case in cases) {
    expect_error(do.call(envelope_to_voxel, case$args), case$pattern,
                 fixed = TRUE, label = case$label)
  }

  # the supplied shapes and the default all have to survive the profile check
  expect_no_error(
    envelope_to_voxel(x = good_envel, depths = good_depths, profile = NULL)
  )
  expect_no_error(
    envelope_to_voxel(x = good_envel, depths = good_depths,
                      profile = profile_flat)
  )
  expect_no_error(
    envelope_to_voxel(x = good_envel, depths = good_depths,
                      profile = profile_equal)
  )
})

test_that("envelope_to_voxel() accepts depths that merely coerce to numeric", {
  good_envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)

  # character and integer depths are coerced rather than rejected as the wrong
  # type, and the coerced values are what the layer names are built from
  expect_named(
    envelope_to_voxel(good_envel, depths = c("0", "100")),
    paste0("presence_depth=", c(0, 100))
  )
  expect_named(
    envelope_to_voxel(good_envel, depths = 0:3),
    paste0("presence_depth=", 0:3)
  )
})

test_that("envelope_to_voxel() occupies every slab the envelope overlaps", {
  # Each level stands for a slab: under "top" from the level to the next one,
  # under "midpoint" between the midpoints to its neighbours, and the deepest
  # level is a point in both. A level is occupied when its slab shares water
  # with the envelope. Touching is not overlapping. `expected` holds the first
  # rows of the voxel's values, one row per cell; NULL bounds is the default,
  # "midpoint". Rows about the "top" convention set it. `env` holds the
  # as_envelope() arguments.
  full <- make_footprint(c(1, 1, 1, 1))
  cases <- list(
    # the case the point-in-interval rule missed: [10, 20] contains none of
    # the levels, but lies inside the slab the 0 m level stands for
    list(label = "between two levels", env = list(full, 10, 20),
         depths = c(0, 100), expected = rbind(c(1L, NA))),
    # starts inside the 50-100 slab and ends inside 200-300
    list(label = "across several slabs", env = list(full, 70, 210),
         depths = c(0, 50, 100, 200, 300),
         expected = rbind(c(NA, 1L, 1L, 1L, NA))),
    # wholly inside one slab: the level below owns its own top
    list(label = "inside one slab", env = list(full, 60, 99),
         depths = c(0, 50, 100, 200), bounds = "top",
         expected = rbind(c(NA, 1L, NA, NA))),
    # ending exactly at 100 only touches the 100-200 slab
    list(label = "ending on a level", env = list(full, 60, 100),
         depths = c(0, 50, 100, 200), bounds = "top",
         expected = rbind(c(NA, 1L, NA, NA))),
    # reaching past 100 by any amount brings that level in
    list(label = "reaching past a level", env = list(full, 60, 101),
         depths = c(0, 50, 100, 200), bounds = "top",
         expected = rbind(c(NA, 1L, 1L, NA))),
    # Adjacent ranges [0, 100] and [100, 200] share the boundary at 100 m and
    # nothing else, so they occupy no common level whenever a slab edge falls
    # at 100 m: under "top" the 100 level itself, under "midpoint" halfway
    # between 50 and 150. Under "top" the shared level goes to the deeper range
    # alone, and the deepest level is a point, so [100, 200] holds 100 only.
    list(label = "upper of two adjacent ranges, top",
         env = list(full, 0, 100), depths = c(0, 100, 200), bounds = "top",
         expected = rbind(c(1L, NA, NA))),
    list(label = "lower of two adjacent ranges, top",
         env = list(full, 100, 200), depths = c(0, 100, 200), bounds = "top",
         expected = rbind(c(NA, 1L, NA))),
    list(label = "upper of two adjacent ranges, midpoint edge at 100",
         env = list(full, 0, 100), depths = c(0, 50, 150, 200),
         bounds = "midpoint",
         expected = rbind(c(1L, 1L, NA, NA))),
    list(label = "lower of two adjacent ranges, midpoint edge at 100",
         env = list(full, 100, 200), depths = c(0, 50, 150, 200),
         bounds = "midpoint",
         expected = rbind(c(NA, NA, 1L, 1L))),
    # With evenly spaced levels under "midpoint" the 100 m boundary falls
    # inside the 100 m slab, so both ranges share water with it and both
    # occupy it. That is the discretisation, not a touch.
    list(label = "upper of two adjacent ranges, even midpoint",
         env = list(full, 0, 100), depths = c(0, 50, 100, 150, 200),
         bounds = "midpoint", expected = rbind(c(1L, 1L, 1L, NA, NA))),
    list(label = "lower of two adjacent ranges, even midpoint",
         env = list(full, 100, 200), depths = c(0, 50, 100, 150, 200),
         bounds = "midpoint", expected = rbind(c(NA, NA, 1L, 1L, 1L))),
    # overlap widens the old point rule and never narrows it: every level
    # whose own value lies inside the envelope (50, 100, 200) is occupied
    list(label = "levels inside the envelope", env = list(full, 50, 250),
         depths = c(0, 50, 100, 200, 300),
         expected = rbind(c(NA, 1L, 1L, 1L, NA))),
    # 50 heads the slab 50-100, which the envelope runs through; 200 heads
    # 200-300, which it only touches
    list(label = "a level at the envelope's top, not its bottom",
         env = list(full, 50, 200), depths = c(0, 50, 100, 200, 300),
         bounds = "top", expected = rbind(c(NA, 1L, 1L, NA, NA))),
    # reaching the deepest level is enough to be recorded at it
    list(label = "reaching the deepest level", env = list(full, 290, 400),
         depths = c(0, 100, 200, 300), bounds = "top",
         expected = rbind(c(NA, NA, 1L, 1L))),
    # the deepest slab is closed at max(depths), so an envelope sitting
    # exactly on the last level is caught rather than falling off the bottom
    list(label = "on the deepest level, midpoint", env = list(full, 300, 300),
         depths = c(0, 100, 200, 300), bounds = "midpoint",
         expected = rbind(c(NA, NA, NA, 1L))),
    # [60, 70] sits in the 0-100 slab under "top", but in the 50-150 slab of
    # the 100 m level once the edges move to the midpoints
    list(label = "[60, 70] under top", env = list(full, 60, 70),
         depths = c(0, 100, 200), bounds = "top",
         expected = rbind(c(1L, NA, NA))),
    list(label = "[60, 70] under midpoint", env = list(full, 60, 70),
         depths = c(0, 100, 200), bounds = "midpoint",
         expected = rbind(c(NA, 1L, NA))),
    # by default the slabs are 0-50, 50-150, 150-250 and 250-300 m, so
    # [60, 210] holds the 100 and 200 m levels; under "top" it would also
    # hold 0 m
    list(label = "default bounds", env = list(full, 60, 210),
         depths = c(0, 100, 200, 300), expected = rbind(c(NA, 1L, 1L, NA)),
         differs_from_top = TRUE),
    list(label = "[60, 210] under top", env = list(full, 60, 210),
         depths = c(0, 100, 200, 300), bounds = "top",
         expected = rbind(c(1L, 1L, 1L, NA))),
    # cell 4 of make_footprint() is NA, so it is NA at every depth. [0, 200]
    # shares water with all three default slabs (0-50, 50-150, 150-200 m).
    list(label = "a cell absent from the envelope",
         env = list(make_footprint(), 0, 200), depths = c(0, 100, 200),
         expected = rbind(c(1L, 1L, 1L), c(1L, 1L, 1L), c(1L, 1L, 1L),
                          c(NA, NA, NA))),
    # each range holds the slabs it runs through and not the one it ends on:
    # [0, 100], [100, 300], [200, 300] and an absent cell
    list(label = "per-cell depth limits",
         env = list(make_footprint(), make_footprint(c(0, 100, 200, 0)),
                    make_footprint(c(100, 300, 300, 300))),
         depths = c(0, 100, 200, 300), bounds = "top",
         expected = rbind(c(1L, NA, NA, NA), c(NA, 1L, 1L, NA),
                          c(NA, NA, 1L, NA), c(NA, NA, NA, NA)))
  )
  for (case in cases) in_case(case$label, {
    envel <- do.call(as_envelope, case$env)
    # a NULL bounds is left out, so the default applies
    args <- Filter(Negate(is.null),
                   list(envel, depths = case$depths, bounds = case$bounds))
    warnings <- capture_warnings(v <- do.call(envelope_to_voxel, args))
    expect_identical(warnings, character(), label = case$label)

    expect_true(methods::is(v, "SpatVoxel"), label = case$label)
    expect_named(v, paste0("presence_depth=", case$depths), label = case$label)
    vals <- terra::values(v)
    expect_identical(unname(vals[seq_len(nrow(case$expected)), , drop = FALSE]),
                     case$expected, label = case$label)
    # the default bounds are "midpoint"
    if (is.null(case$bounds)) {
      mid <- envelope_to_voxel(envel, depths = case$depths, bounds = "midpoint")
      expect_identical(vals, terra::values(mid), label = case$label)
    }
    if (isTRUE(case$differs_from_top)) {
      top <- envelope_to_voxel(envel, depths = case$depths, bounds = "top")
      expect_false(identical(vals, terra::values(top)), label = case$label)
    }
  })
})

test_that("envelope_to_voxel() warns about an envelope outside the levels", {
  # An open-ended deepest slab would swallow any range below the grid and
  # report it as resolved; keeping it a point leaves the gap visible. The
  # shallowest level is the top of the first slab, so a cell entirely above it
  # is out of reach in the same way.
  full <- make_footprint(c(1, 1, 1, 1))
  cases <- list(
    list(label = "below the deepest level", env = list(full, 310, 400),
         depths = c(0, 100, 200, 300), bounds = "top"),
    list(label = "far below, top", env = list(full, 500, 600),
         depths = c(0, 100, 200), bounds = "top"),
    list(label = "far below, midpoint", env = list(full, 500, 600),
         depths = c(0, 100, 200), bounds = "midpoint"),
    list(label = "above the shallowest level", env = list(full, 5, 40),
         depths = c(50, 100, 200))
  )
  for (case in cases) in_case(case$label, {
    envel <- do.call(as_envelope, case$env)
    # a NULL bounds is left out, so the default applies
    args <- Filter(Negate(is.null),
                   list(envel, depths = case$depths, bounds = case$bounds))
    expect_warning(v <- do.call(envelope_to_voxel, args),
                   class = "ocean3d_warning_envelope_outside_slabs",
                   label = case$label)
    expect_true(all(is.na(terra::values(v))), label = case$label)
  })
})

test_that(".depth_slabs() gives each level the slab it stands for", {
  cases <- list(
    # each slab runs to the next level, and the deepest level has no next one
    # and collapses to a point
    list(depths = c(0, 100, 200), bounds = "top",
         lower = c(0, 100, 200), upper = c(100, 200, 200)),
    # SPEC.md: WOA's 0 m layer has bounds 0-2.5, the 5 m layer 2.5-7.5, the
    # 10 m layer 7.5-12.5
    list(depths = c(0, 5, 10, 15), bounds = "midpoint",
         lower = c(0, 2.5, 7.5, 12.5), upper = c(2.5, 7.5, 12.5, 15)),
    # the outer edges are clamped rather than extrapolated: half a gap above
    # WOA's shallowest level would be a negative depth, and half a gap below
    # the deepest would invent reach the levels never had
    list(depths = c(0, 100, 200), bounds = "midpoint",
         lower = c(0, 50, 150), upper = c(50, 150, 200)),
    # a lone level collapses to a point either way
    list(depths = 100, bounds = "top", lower = 100, upper = 100),
    list(depths = 100, bounds = "midpoint", lower = 100, upper = 100)
  )
  for (case in cases) {
    label <- paste(case$bounds, toString(case$depths))
    slab <- .depth_slabs(case$depths, case$bounds)
    expect_identical(slab$lower, case$lower, label = label)
    expect_identical(slab$upper, case$upper, label = label)
  }
})

test_that(".depth_slabs() tiles the depth span without gaps or overlaps", {
  # each slab picks up exactly where the last left off, so the levels partition
  # [min(depths), max(depths)]
  for (bounds in c("top", "midpoint")) {
    slab <- .depth_slabs(c(0, 50, 100, 400), bounds)
    expect_identical(slab$lower[-1], slab$upper[-length(slab$upper)],
                     info = bounds)
  }
})

test_that("envelope_to_voxel() writes `values` and spreads them by `profile`", {
  # `values` replaces presence at occupied levels; a profile divides it across
  # them. `expected` holds rows `at` of the voxel's values, and `sums` each
  # cell's total over depth. `env` holds the as_envelope() arguments.
  depths <- c(0, 100, 200, 300)
  full <- make_footprint(c(1, 1, 1, 1))
  cases <- list(
    # same occupancy as presence, carrying 5 instead of 1
    list(label = "scalar values",
         env = list(make_footprint(), 0,
                    make_footprint(c(150, 350, 350, 350))),
         depths = depths, values = 5, varname = "score",
         at = 1:2, expected = rbind(c(5L, 5L, NA, NA), rep(5L, 4))),
    list(label = "per-cell values", env = list(full, 50, 250),
         depths = c(0, 50, 100, 200, 300), values = c(10, 20, 30, 40),
         varname = "effort", at = c(1, 4),
         expected = rbind(c(NA, 10, 10, 10, NA), c(NA, 40, 40, 40, NA))),
    # the two distinct intervals occupy 2 and 4 levels, so get different
    # shares
    list(label = "profile_equal splits presence",
         env = list(make_footprint(), 0,
                    make_footprint(c(150, 350, 350, 350))),
         depths = depths, profile = profile_equal, varname = "time",
         at = 1:2, expected = rbind(c(0.5, 0.5, NA, NA), rep(0.25, 4))),
    # cells occupy 2, 3, 4 and 4 levels; each still sums to its input
    list(label = "profile_equal conserves each cell's value",
         env = list(full, 0, make_footprint(c(150, 250, 350, 350))),
         depths = depths,
         values = c(10, 20, 30, 40), profile = profile_equal,
         varname = "effort", at = 1, expected = rbind(c(5, 5, NA, NA)),
         sums = c(10, 20, 30, 40)),
    # full value at all four levels: a footprint map, not a conserved total
    list(label = "no profile does not divide the value",
         env = list(full, 0, 350), depths = depths, values = c(10, 20, 30, 40),
         varname = "presence", at = 1, expected = rbind(rep(10, 4)),
         sums = c(40, 80, 120, 160))
  )
  for (case in cases) in_case(case$label, {
    envel <- do.call(as_envelope, case$env)
    values <- case$values
    if (length(values) > 1) {
      values <- make_footprint(values)
    }
    args <- list(envel, depths = case$depths, values = values,
                 profile = case$profile, varname = case$varname)
    v <- do.call(envelope_to_voxel, args)

    expect_named(v, paste0(case$varname, "_depth=", case$depths),
                 label = case$label)
    vals <- terra::values(v)
    expect_identical(unname(vals[case$at, , drop = FALSE]), case$expected,
                     label = case$label)
    if (!is.null(case$sums)) {
      expect_identical(unname(rowSums(vals, na.rm = TRUE)), case$sums,
                       label = case$label)
    }
  })
})

test_that("envelope_to_voxel() presence is the same as values = 1", {
  fp <- make_footprint()
  dmax <- terra::setValues(terra::rast(fp), c(100, 300, 300, 300))
  envel <- as_envelope(fp, depth_min = 0, depth_max = dmax)
  depths <- c(0, 100, 200, 300)

  expect_identical(
    terra::values(envelope_to_voxel(envel, depths)),
    terra::values(envelope_to_voxel(envel, depths, values = 1))
  )
})

test_that("envelope_to_voxel() sorts and deduplicates the requested depths", {
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)),
                       depth_min = 0, depth_max = 300)

  v <- envelope_to_voxel(envel, depths = c(200, 0, 100, 200))

  expect_named(v, paste0("presence_depth=", c(0, 100, 200)))
})

test_that("envelope_to_voxel() round-trips a voxel built on the same depths", {
  # a gap-free voxel is the case where the two conversions are inverses
  depths <- c(0, 100, 200, 300)
  solid <- make_multidepth_rast(depths = depths,
                                vals = list(c(10, NA, 5, NA),
                                            c(11, 20, 6, NA),
                                            c(12, 21, 7, NA),
                                            c(NA, 22, 8, NA)))
  v <- as_voxel(solid)

  # By default each level sits in the middle of its slab, so the envelope's
  # limits (the shallowest and deepest occupied levels) fall inside the end
  # slabs and every occupied level comes back. The values themselves are not
  # recovered, since an envelope carries only the depth limits.
  # By default each level sits in the middle of its slab, so the envelope's
  # limits (the shallowest and deepest occupied levels) fall inside the end
  # slabs and every occupied level comes back. The values themselves are not
  # recovered, since an envelope carries only the depth limits.
  back <- envelope_to_voxel(voxel_to_envelope(v), depths = depths,
                            varname = "temp")
  expect_named(back, names(v))
  expect_identical(!is.na(terra::values(back)), !is.na(terra::values(v)))

  # Under "top" a level names the top of its slab, so the envelope only
  # touches the deepest occupied level's slab, and touching is not
  # overlapping: the round trip drops that level in each cell.
  back <- envelope_to_voxel(voxel_to_envelope(v), depths = depths,
                            varname = "temp", bounds = "top")
  present <- !is.na(terra::values(back))
  # 0-200 -> 0, 100
  expect_identical(unname(present[1, ]), c(TRUE, TRUE, FALSE, FALSE))
  # 100-300 -> 100, 200
  expect_identical(unname(present[2, ]), c(FALSE, TRUE, TRUE, FALSE))
  # 0-300 -> 0, 100, 200
  expect_identical(unname(present[3, ]), c(TRUE, TRUE, TRUE, FALSE))
  # absent
  expect_identical(unname(present[4, ]), rep(FALSE, 4))
})

# vect_to_envelope() ----

test_that("vect_to_envelope() rejects unusable input", {
  input <- "ocean3d_error_envelope_input"
  grid <- "ocean3d_error_grid_value"
  # a raster that doesn't overlap make_polygon()
  elsewhere <- terra::rast(nrows = 5, ncols = 5, xmin = 5, xmax = 10,
                           ymin = 5, ymax = 10, crs = "EPSG:4326")
  # a template in another CRS
  mercator <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5,
                          ymin = 0, ymax = 5, crs = "EPSG:3857")
  # depth rasters whose resolution, extent or CRS differ from the template's;
  # with a different CRS but the same extent numbers, terra would otherwise
  # return a wrong answer rather than fail
  bad_res <- terra::rast(nrows = 9, ncols = 9, xmin = 0, xmax = 5,
                         ymin = 0, ymax = 5, vals = 50, crs = "EPSG:4326")
  bad_ext <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 9,
                         ymin = 0, ymax = 9, vals = 50, crs = "EPSG:4326")
  bad_crs <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5,
                         ymin = 0, ymax = 5, vals = 50, crs = "EPSG:3857")
  poly <- make_polygon()
  tmpl <- make_template()

  cases <- list(
    list(args = list("test", tmpl, 0, 10), class = input,
         pattern = "`polygon` must be a SpatVector, sf or sfc object"),
    list(args = list(poly, "test", 0, 10), class = input,
         pattern = "`template` must be a SpatRaster"),
    list(args = list(poly, elsewhere, 0, 10),
         class = "ocean3d_error_envelope_empty",
         pattern = "The polygons don't overlap any cell of `template`"),
    list(args = list(poly, mercator, 0, 10), class = input,
         pattern = "`polygon` and `template` have different CRS."),
    list(args = list(poly, tmpl, 0, list(10, bad_res)), class = grid,
         pattern = "does not align with `template`"),
    list(args = list(poly, tmpl, 0, list(10, bad_ext)), class = grid,
         pattern = "does not align with `template`"),
    list(args = list(poly, tmpl, 0, list(10, bad_crs)), class = grid,
         pattern = "different CRS than `template`"),
    list(args = list(poly, tmpl, list(), 10), class = input,
         pattern = "`depth_min` is empty"),
    list(args = list(poly, tmpl, 0, list()), class = input,
         pattern = "`depth_max` is empty"),
    list(args = list(poly, tmpl, 0, list(10, "deep")), class = grid,
         pattern = "element 2 must be a single non-missing number"),
    # `na.rm = TRUE` would drop an NA silently, widening rather than narrowing.
    list(args = list(poly, tmpl, 0, list(10, NA_real_)), class = grid,
         pattern = "element 2 must be a single non-missing number")
  )
  for (case in cases) {
    expect_error(do.call(vect_to_envelope, case$args), case$pattern,
                 class = case$class, label = case$pattern)
  }

  # TODO: implement check for the depth_min and depth_max inputs
})

test_that("vect_to_envelope() accepts scalar depth_min and depth_max", {
  r <- make_template()

  expected_result <- rast(c(
    depth_min = rast(r, vals = 10),
    depth_max = rast(r, vals = 20)
  ))
  expected_result$depth_min <- expected_result$depth_min %>%
    mask(make_polygon())
  expected_result$depth_max <- expected_result$depth_max %>%
    mask(make_polygon())
  expected_result <- as_envelope(expected_result)

  result <- vect_to_envelope(make_polygon(), make_template(),
                             depth_min = 10, depth_max = 20)

  # terra objects are pointers, not identical according to expect_identical()
  # nolint start: expect_identical_linter.
  expect_true(terra::identical(result, expected_result))
  # nolint end
})

test_that("vect_to_envelope() takes the narrowest limits of a depth list", {
  # depth_min = c(10, r): each cell takes the deeper of the two;
  # depth_max = c(10, r): each cell takes the shallower of the two
  r <- make_template()
  terra::values(r) <- seq_along(values(r))

  inside <- mask(r, make_polygon())
  deepest_min <- max(inside, 10)
  names(deepest_min) <- "depth_min"
  deepest_min$depth_max <- rast(r, vals = 25)

  shallowest_max <- rast()
  shallowest_max$depth_min <- rast(r, vals = 0)
  shallowest_max$depth_max <- min(inside, 10)

  cases <- list(
    list(layer = "depth_min", depth_min = c(10, r), depth_max = 25,
         expected = deepest_min),
    list(layer = "depth_max", depth_min = 0, depth_max = c(10, r),
         expected = shallowest_max)
  )
  for (case in cases) {
    result <- vect_to_envelope(make_polygon(), make_template(),
                               depth_min = case$depth_min,
                               depth_max = case$depth_max)
    expected_result <- as_envelope(case$expected)

    # terra objects are pointers, not identical according to expect_identical()
    # nolint start: expect_identical_linter.
    expect_true(terra::identical(result[case$layer],
                                 expected_result[case$layer]),
                label = case$layer)
    # nolint end
  }
})

# TODO: deal with case where depth_max is shallower than depth_min
test_that("vect_to_envelope() fills NA where depth_max < depth_min", {
  r <- make_template()
  terra::values(r) <- seq_along(values(r))

  # minimum depth of envelope is 10
  depth_min_rast <- rast(r, vals = 10)
  # maximum depth is either 16 or the raster r value
  depth_max_rast <- r %>%
    mask(make_polygon()) %>%
    min(16)

  # valid cells are where depth_max is greater than depth_min
  diffs <- depth_max_rast - depth_min_rast
  msk <- ifel(diffs > 0, 1, NA)

  depth_min_rast <- mask(depth_min_rast, msk)
  depth_max_rast <- mask(depth_max_rast, msk)

  expected_result <- rast(c(
    depth_min = depth_min_rast,
    depth_max = depth_max_rast
  )) %>% as_envelope()

  result <- vect_to_envelope(make_polygon(), make_template(),
                             depth_min = 10, depth_max = c(16, r))

  # terra objects are pointers, not identical according to expect_identical()
  # nolint start: expect_identical_linter.
  expect_true(terra::identical(result, expected_result))
  # nolint end
})

test_that("vect_to_envelope() stops when no cell would be left", {
  r <- make_template()
  terra::values(r) <- seq_along(values(r))
  empty <- "ocean3d_error_envelope_empty"

  expect_error(
    vect_to_envelope(make_polygon(), make_template(),
                     depth_min = 50, depth_max = c(10, r)),
    "have `depth_max` shallower than `depth_min` \\(the limits may be swapped",
    class = empty
  )
  expect_error(
    vect_to_envelope(make_polygon(), make_template(),
                     depth_min = 30, depth_max = 30),
    "have them equal, so no thickness", class = empty
  )
  unknown <- make_template()
  terra::values(unknown) <- NA_real_
  expect_error(
    vect_to_envelope(make_polygon(), make_template(),
                     depth_min = 0, depth_max = list(100, unknown)),
    "have no `depth_min` or `depth_max`", class = empty
  )

  # Cells that fail in different ways are counted separately.
  mixed <- terra::setValues(make_template(), 30)
  covered <- which(!is.na(terra::values(terra::mask(mixed, make_polygon()))))
  mixed[covered[1]] <- 5
  err <- expect_error(
    vect_to_envelope(make_polygon(), make_template(),
                     depth_min = 30, depth_max = list(100, mixed)),
    class = empty
  )
  expect_match(conditionMessage(err), "1 has `depth_max` shallower",
               fixed = TRUE)
  expect_match(conditionMessage(err), "have them equal", fixed = TRUE)
})

test_that("vect_to_envelope() keeps decimal depths when terra writes to disk", {
  old <- terra::terraOptions(print = FALSE)$todisk
  terra::terraOptions(todisk = TRUE)
  on.exit(terra::terraOptions(todisk = old))
  out <- vect_to_envelope(make_polygon(), make_template(),
                          depth_min = 10.1, depth_max = 35.1)
  v <- stats::na.omit(terra::values(out))
  expect_identical(unique(v[, "depth_min"]), 10.1)
  expect_identical(unique(v[, "depth_max"]), 35.1)
  # terra's own default is restored afterwards.
  expect_identical(terra::terraOptions(print = FALSE)$datatype, "FLT4S")
})

test_that("vect_to_envelope() returns a valid SpatEnvelope, visibly", {
  # sf, sfc and SpatVector polygons are all accepted. The constructor assigns
  # to `out` last; without a trailing bare `out` the result would come back
  # invisibly and never print at the console.
  polygons <- list(sf = make_polygon(),
                   sfc = sf::st_geometry(make_polygon()),
                   SpatVector = terra::vect(make_polygon()))
  for (form in names(polygons)) {
    out <- withVisible(
      vect_to_envelope(polygons[[form]], make_template(), 0, 10)
    )

    expect_true(out$visible, label = form)
    expect_true(methods::is(out$value, "SpatEnvelope"), label = form)
    expect_true(isTRUE(methods::validObject(out$value)), label = form)
    expect_named(out$value, c("depth_min", "depth_max"), label = form)
  }
})

# Ported from the retired voxelize_range(): the seafloor is now just another
# `depth_max` constraint rather than a dedicated `bathymetry` argument.
test_that("vect_to_envelope() clamps depth_max to a shallower seafloor", {
  seafloor <- make_template(seq(10, 250, length.out = 25))

  out <- vect_to_envelope(make_polygon(), make_template(),
                          depth_min = 50, depth_max = list(200, seafloor))

  got <- terra::values(out$depth_max)
  got <- got[!is.na(got)]
  expect_true(all(got <= 200))
  expect_true(any(got < 200))
})

# Ported from the retired voxelize_range(): cells whose bed sits above the
# shallowest depth the animal occupies carry no envelope at all.
test_that("vect_to_envelope() drops cells whose seafloor is above depth_min", {
  seafloor <- make_template(seq(10, 250, length.out = 25))

  out <- vect_to_envelope(make_polygon(), make_template(),
                          depth_min = 100, depth_max = list(200, seafloor))

  kept <- !is.na(terra::values(out$depth_min))
  expect_true(any(kept))
  expect_true(all(terra::values(seafloor)[kept] > 100))
})

test_that("vect_to_envelope() keeps NA where a depth_max raster is NA", {
  # NA in the depth_max parameter (ex. bathymetry) means there is no valid
  # depth at that cell. These should stay NA in the output as well.
  vals <- seq(10, 250, length.out = 25)
  gaps <- c(7, 13, 19)              # cells inside the polygon footprint
  vals[gaps] <- NA
  seafloor <- make_template(vals)
  new <- vect_to_envelope(make_polygon(), make_template(),
                          depth_min = 50, depth_max = list(200, seafloor))

  expect_true(all(is.na(terra::values(new$depth_max)[gaps])))
  expect_true(all(is.na(terra::values(new$depth_min)[gaps])))
})

test_that("constructors accept tagged input as a footprint or template", {
  field <- make_field()
  # A voxel layer used as a footprint is still tagged SpatVoxel.
  e <- as_envelope(field[[1]], depth_min = 0, depth_max = 100)
  expect_identical(class(e)[[1]], "SpatEnvelope")
  expect_identical(terra::values(e)[, "depth_min"], rep(0, 9))

  poly <- sf::st_as_sf(make_left_polygon())
  e2 <- vect_to_envelope(poly, field, depth_min = 0, depth_max = 100)
  expect_identical(class(e2)[[1]], "SpatEnvelope")
  expect_identical(terra::values(e2)[, "depth_max"],
                   c(100, 100, NA, 100, 100, NA, 100, 100, NA))
})

# ---- occupied() -------------------------------------------------------------

# Local helper: a 1x3 voxel over two standard depths, named by variable.
make_named_voxel <- function(vals, varname, depths = c(0, 100)) {
  r <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                   nlyrs = length(depths))
  terra::values(r) <- vals
  as_voxel(r, depths = depths, varname = varname)
}

test_that("occupied() turns a voxel into a 1-or-NA presence voxel", {
  # Presence is 1 where `fun` holds, applied one depth layer at a time, and NA
  # elsewhere. A presence voxel comes back unchanged.
  v <- as_voxel(make_multidepth_rast())
  vals <- terra::values(v)
  once <- occupied(v)

  cases <- list(
    list(label = "default: where the variable was recorded", x = v,
         fun = NULL, expected = ifelse(is.na(vals), NA, 1)),
    list(label = "x > 10", x = v, fun = function(x) x > 10,
         expected = ifelse(!is.na(vals) & vals > 10, 1, NA)),
    list(label = "a presence voxel", x = once, fun = NULL,
         expected = terra::values(once))
  )
  for (case in cases) {
    out <- if (is.null(case$fun)) {
      occupied(case$x)
    } else {
      occupied(case$x, case$fun)
    }

    expect_true(methods::is(out, "SpatVoxel"), label = case$label)
    expect_named(out, paste0("presence_depth=", c(0, 100, 200, 300)),
                 label = case$label)
    expect_identical(depths(out), depths(v), label = case$label)
    expect_identical(unname(terra::values(out)), unname(case$expected),
                     label = case$label)
  }
})

test_that("occupied() rejects input that is not a voxel", {
  expect_error(occupied(make_footprint()), "must be a SpatVoxel")
  expect_error(occupied(as_envelope(make_footprint(), 0, 100)),
               "must be a SpatVoxel")
})

test_that("occupied() gives a reducing predicate its per-depth meaning", {
  # `fun` gets each layer's values as a numeric vector, so the predicate
  # compares each cell against its own depth's mean. On a one-layer
  # SpatRaster, terra's mean() returns the layer itself, making
  # `x > mean(x)` all-FALSE.
  v <- make_named_voxel(cbind(c(1, 2, 3), c(10, 20, 30)), "temp")
  above_mean <- function(x) x > mean(x, na.rm = TRUE)

  out <- occupied(v, above_mean)

  # Above each layer's own mean (2 and 20): the third cell in both layers.
  expect_identical(unname(terra::values(out)),
               unname(cbind(c(NA, NA, 1), c(NA, NA, 1))))
  expect_gt(volume(out), 0)

  # ...and it agrees with the same predicate through voxel_to_envelope().
  expect_identical(terra::values(voxel_to_envelope(v, fun = above_mean)),
               terra::values(voxel_to_envelope(out)))
})

test_that("occupied() lets each side of a query use its own predicate", {
  temp <- make_named_voxel(cbind(c(20, 20, 4), c(20, 4, 4)), "temp")
  effort <- make_named_voxel(cbind(c(0, 5, 5), c(5, 5, 0)), "effort")

  # Warm water and fished. One shared predicate cannot express both cutoffs.
  out <- intersect_3d(occupied(temp, function(v) v > 15),
                      occupied(effort, function(v) v > 0))

  expect_identical(unname(terra::values(out)),
               unname(cbind(c(NA, 1, NA), c(1, NA, NA))))
})

test_that("tagged 3D objects that a terra operation broke are refused", {
  # terra keeps the class tag through `[[`, `names<-` and arithmetic without
  # re-running validity, so an object can be labelled SpatVoxel or
  # SpatEnvelope and not be one. The verbs, and the idempotent fast path of
  # as_voxel(), must not wave it through on the tag.
  e <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)
  # swapping the two layers yields an envelope whose interval runs backwards;
  # volume() would otherwise report a negative volume for it
  inverted <- e[[c("depth_max", "depth_min")]]
  names(inverted) <- c("depth_min", "depth_max")

  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = 3)
  terra::values(r) <- 1
  v <- as_voxel(r, depths = c(0, 100, 200), varname = "presence")
  expect_true(methods::validObject(v, test = TRUE))
  renamed <- v
  names(renamed) <- c("a", "b", "c")

  r4 <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
  terra::values(r4) <- 1
  terra::depth(r4) <- c(0, 100, 200)

  cases <- list(
    list(label = "inverted envelope", object = inverted,
         class = "SpatEnvelope", call = volume, pattern = "not a valid one"),
    # arithmetic can lift depths above the surface
    list(label = "envelope above the surface", object = e - 500,
         class = "SpatEnvelope"),
    # selecting the first layer twice leaves depths 0, 0, 100; as_voxel()
    # refuses it with the message a plain raster would get
    list(label = "voxel with a repeated layer", object = v[[c(1, 1, 2)]],
         class = "SpatVoxel", call = as_voxel, pattern = "duplicate depth"),
    # a reversed voxel would otherwise collapse to an inverted envelope
    list(label = "reversed voxel", object = as_voxel(r4)[[c(3, 2, 1)]],
         class = "SpatVoxel", call = voxel_to_envelope,
         pattern = "not a valid one")
  )
  for (case in cases) {
    expect_true(methods::is(case$object, case$class), label = case$label)
    expect_false(isTRUE(methods::validObject(case$object, test = TRUE)),
                 label = case$label)
    if (!is.null(case$call)) {
      expect_error(case$call(case$object), case$pattern, label = case$label)
    }
  }

  # Renaming doesn't break a voxel: the depths are in terra::depth(). Passing
  # `depths` gives readable names again.
  expect_true(methods::validObject(renamed, test = TRUE))
  expect_identical(depths(as_voxel(renamed)), c(0, 100, 200))
  rebuilt <- as_voxel(renamed, depths = c(0, 100, 200), varname = "presence")
  expect_true(methods::validObject(rebuilt, test = TRUE))
  expect_named(rebuilt, names(v))
})

test_that("envelope depth rules hold for file-backed rasters", {
  # The depth rules are read from memory when the raster is in memory and
  # summarised block by block when it is on disk; both must catch the same
  # breaches.
  paths <- character()
  on.exit(unlink(paths))
  on_disk <- function(dmin, dmax) {
    r <- terra::rast(make_footprint(), nlyrs = 2)
    terra::values(r) <- cbind(dmin, dmax)
    names(r) <- c("depth_min", "depth_max")
    path <- tempfile(fileext = ".tif")
    paths <<- c(paths, path)
    terra::writeRaster(r, path)
    out <- terra::rast(path)
    expect_false(any(terra::inMemory(out)))
    out
  }

  good <- on_disk(c(0, 10, 20, NA), c(100, 10, 50, NA))
  expect_s4_class(as_envelope(good), "SpatEnvelope")
  expect_true(methods::validObject(methods::new("SpatEnvelope", good),
                                   test = TRUE))

  cases <- list(
    list(label = "negative depth_min", object = on_disk(c(-5, 0, 0, NA),
                                                        c(10, 10, 10, NA)),
         constructor = "positive metres increasing downward",
         validity = "positive metres increasing downward"),
    list(label = "inverted cell", object = on_disk(c(0, 50, 0, NA),
                                                   c(10, 20, 10, NA)),
         constructor = "depth_max must be at least depth_min",
         validity = "depth_max must be at least depth_min")
  )
  for (case in cases) {
    expect_error(as_envelope(case$object), case$constructor,
                 label = case$label)
    expect_error(methods::new("SpatEnvelope", case$object), case$validity,
                 label = case$label)
  }
})


# ---- depth edges (#138) -----------------------------------------------------

# Copernicus-style levels at cell centres, with edges that aren't halfway.
edge_depths <- c(0.5, 1.5, 3)
edge_bounds <- cbind(c(0, 1, 2), c(1, 2, 4))

test_that("as_voxel() stores depth edges keyed by depth", {
  r <- terra::rast(make_grid(), nlyrs = 3)
  terra::values(r) <- 1
  # Given out of order: the edges follow their layers when they are sorted.
  names(r) <- paste0("t_depth=", edge_depths[c(3, 1, 2)])
  terra::depth(r) <- edge_depths[c(3, 1, 2)]
  v <- as_voxel(r, depth_bounds = edge_bounds[c(3, 1, 2), ])

  expect_identical(names(v), paste0("t_depth=", edge_depths))
  expect_identical(names(v@depth_edges), c("depth", "lower", "upper"))
  expect_equal(v@depth_edges$depth, edge_depths[c(3, 1, 2)])
  expect_equal(edges_of(v), edge_bounds)

  # Without depth_bounds a voxel carries none; a data frame works too.
  expect_identical(nrow(as_voxel(r)@depth_edges), 0L)
  df <- as_voxel(r, depth_bounds = as.data.frame(edge_bounds[c(3, 1, 2), ]))
  expect_equal(edges_of(df), edge_bounds)

  # A voxel's own table goes back in, matched by depth whatever its order.
  again <- as_voxel(make_voxel(2, edge_depths),
                    depth_bounds = v@depth_edges[c(2, 3, 1), ])
  expect_equal(edges_of(again), edge_bounds)

  # Depths that round when written into layer names still find their edges.
  tenths <- as_voxel(terra::rast(nrows = 2, ncols = 2, nlyrs = 4, vals = 1),
                     depths = (0:3) * 0.1,
                     depth_bounds = cbind((0:3) * 0.1, c((1:3) * 0.1, 0.4)))
  expect_equal(edges_of(tenths), cbind((0:3) * 0.1, c((1:3) * 0.1, 0.4)))
})

test_that("terra operations keep a voxel's depths", {
  v <- make_voxel(1, c(0, 100, 200))
  poly <- terra::vect("POLYGON ((0 0, 2000 0, 2000 2000, 0 2000, 0 0))",
                      crs = terra::crs(v))
  ops <- list(
    crop = function(x) terra::crop(x, terra::ext(0, 1000, 0, 1000)),
    polygon = function(x) terra::crop(x, poly, mask = TRUE),
    mask = function(x) terra::mask(x, poly),
    times = function(x) x * 2,
    ifel = function(x) terra::ifel(x > 0, 1, NA)
  )
  # On the whole voxel, and after selecting layers: terra alone clears its
  # depth flag on a layer subset, so the next operation would reset them.
  picks <- list(
    all = list(x = v, depths = c(0, 100, 200)),
    numbers = list(x = v[[2:3]], depths = c(100, 200)),
    one = list(x = v[[3]], depths = 200),
    name = list(x = v[["presence_depth=100"]], depths = 100),
    dollar = list(x = v$`presence_depth=100`, depths = 100),
    logical = list(x = v[[c(TRUE, FALSE, TRUE)]], depths = c(0, 200)),
    subset = list(x = subset(v, 2:3), depths = c(100, 200))
  )
  for (pick in names(picks)) {
    for (op in names(ops)) {
      out <- ops[[op]](picks[[pick]]$x)
      expect_s4_class(out, "SpatVoxel")
      expect_identical(depths(out), picks[[pick]]$depths,
                       label = paste(pick, op))
    }
  }
  # terra's non-standard subset() too.
  named <- v
  names(named) <- c("top", "mid", "deep")
  expect_identical(depths(subset(named, c(top, deep), NSE = TRUE) * 2),
                   c(0, 200))

  # Replacing a layer with `[[<-`, by number or name, keeps its depth, even
  # when the new values have none.
  blank <- terra::setValues(terra::rast(methods::as(v[[2]], "SpatRaster")), 0)
  replaced <- v
  replaced[[2]] <- blank
  replaced[["presence_depth=200"]] <- blank
  expect_s4_class(replaced, "SpatVoxel")
  expect_identical(depths(replaced), c(0, 100, 200))
  expect_identical(depths(v), c(0, 100, 200))

  # app() drops the depths and the class; as_voxel() then needs `depths`.
  applied <- terra::app(v, function(x) x * 2)
  expect_false(methods::is(applied, "SpatVoxel"))
  expect_error(as_voxel(applied), class = "ocean3d_error_no_depths")
  expect_identical(depths(as_voxel(applied, depths = c(0, 100, 200))),
                   c(0, 100, 200))

  # terra's own c() resets depths to 0 and edits its inputs' depths in place;
  # the SpatVoxel method puts them back and leaves the inputs alone.
  shallow <- v[[1]]
  deep <- v[[2:3]]
  joined <- c(shallow, deep)
  expect_s4_class(joined, "SpatVoxel")
  expect_identical(depths(joined), c(0, 100, 200))
  expect_identical(depths(shallow), 0)
  expect_identical(depths(deep), c(100, 200))
  # Two layers at one depth aren't a voxel: c() returns a plain raster.
  expect_false(methods::is(c(shallow, shallow), "SpatVoxel"))
})

test_that("depth edges survive terra operations that keep the class", {
  v <- as_voxel(make_voxel(1, edge_depths), depth_bounds = edge_bounds)

  # A subset keeps the rows for the layers it kept.
  expect_equal(edges_of(v[[2:3]]), edge_bounds[2:3, ])
  expect_equal(edges_of(as_voxel(v[[c(1, 3)]])), edge_bounds[c(1, 3), ])
  expect_equal(edges_of(terra::crop(v, terra::ext(0, 1, 0, 1))), edge_bounds)
  expect_equal(edges_of(v * 2), edge_bounds)

  # c() keeps the edges of every voxel it combines, each built with only its
  # own layers' edges.
  top <- as_voxel(make_voxel(1, edge_depths[1]),
                  depth_bounds = edge_bounds[1, , drop = FALSE])
  rest <- as_voxel(make_voxel(1, edge_depths[2:3]),
                   depth_bounds = edge_bounds[2:3, ])
  joined <- c(top, rest)
  expect_true(methods::validObject(joined))
  expect_equal(edges_of(joined), edge_bounds)
  expect_equal(volume(joined), volume(v))
  # Given deepest first, the layers are sorted and keep their edges.
  reversed <- c(rest, top)
  expect_s4_class(reversed, "SpatVoxel")
  expect_identical(depths(reversed), edge_depths)
  expect_equal(edges_of(reversed), edge_bounds)
  # Slabs that overlap between voxels aren't a voxel: a plain raster.
  clash <- as_voxel(make_voxel(1, 0.8), depth_bounds = cbind(0.6, 1.2))
  expect_false(methods::is(c(top, clash), "SpatVoxel"))
  # With a voxel that has none, some layers have no edges, so `bounds`
  # applies to all of them.
  mixed <- c(top, make_voxel(1, edge_depths[2:3]))
  expect_true(methods::validObject(mixed))
  expect_null(edges_of(mixed))
  expect_equal(volume(mixed), volume(make_voxel(1, edge_depths)))
  # A subset's rows for layers it dropped don't stand in for another voxel's
  # edges at the same depth: v[[1]] still holds the 1.5 m row (1-2 m).
  other <- as_voxel(make_voxel(1, edge_depths[2]),
                    depth_bounds = cbind(1, 2.5))
  expect_equal(edges_of(c(v[[1]], other)), rbind(c(0, 1), c(1, 2.5)))

  # Re-wrapping keeps them; renaming the layers to other depths drops them.
  expect_identical(as_voxel(v), v)
  expect_identical(nrow(as_voxel(v, depths = c(1, 2, 5))@depth_edges), 0L)
  # New edges replace the old ones.
  wider <- cbind(c(0, 1, 2), c(1, 2, 6))
  expect_equal(edges_of(as_voxel(v, depth_bounds = wider)), wider)
})

test_that("invalid depth edges stop", {
  r <- make_voxel(1, edge_depths)
  bad <- list(
    shape = edge_bounds[1:2, ],
    contains = cbind(c(0, 2, 2), c(1, 2.5, 4)),
    overlap = cbind(c(0, 0.8, 2), c(1, 2, 4)),
    zero_width = cbind(c(0, 1, 3), c(1, 2, 3)),
    infinite = cbind(c(0, 1, 2), c(1, 2, Inf)),
    negative = cbind(c(-1, 1, 2), c(1, 2, 4)),
    no_row = data.frame(depth = c(0.5, 1.5), lower = c(0, 1), upper = c(1, 2))
  )
  messages <- c(
    shape = "one row per depth: got 2 x 2 for 3 depths",
    contains = "positive widths and must contain their centres",
    overlap = "overlap",
    zero_width = "positive widths",
    infinite = "positive widths",
    negative = "positive metres",
    no_row = "no row for depth 3"
  )
  for (case in names(bad)) {
    err <- expect_error(as_voxel(r, depth_bounds = bad[[case]]),
                        class = "ocean3d_error_bounds")
    expect_match(conditionMessage(err), "Bounds for `depth_bounds`",
                 fixed = TRUE)
    expect_match(conditionMessage(err), messages[[case]], fixed = TRUE)
  }

  # Gaps between slabs are allowed.
  gap <- cbind(c(0, 1.2, 2), c(1, 2, 4))
  expect_equal(edges_of(as_voxel(r, depth_bounds = gap)), gap)

  # The validity method checks a voxel whose slot was set by hand.
  v <- as_voxel(r, depth_bounds = edge_bounds)
  v@depth_edges$upper[1] <- 1.2
  expect_error(methods::validObject(v), "depth_edges overlap")
  v@depth_edges <- data.frame(depth = 1)
  expect_error(methods::validObject(v),
               "numeric columns depth, lower and upper")
})

test_that("envelope_to_voxel() places levels with given depth edges", {
  # [1.1, 1.3] lies in the 1.5 m level's slab under "midpoint" (1-2.25), but
  # in the 0.5 m level's slab with these edges.
  fp <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0,
                    ymax = 1)
  terra::values(fp) <- 1
  e <- as_envelope(fp, depth_min = 1.1, depth_max = 1.3)
  shifted <- cbind(c(0, 1.4, 2), c(1.4, 2, 4))

  expect_equal(as.vector(terra::values(
    envelope_to_voxel(e, edge_depths, bounds = "midpoint")
  )), c(NA, 1, NA))
  v <- envelope_to_voxel(e, edge_depths, depth_bounds = shifted)
  expect_equal(as.vector(terra::values(v)), c(1, NA, NA))
  expect_equal(edges_of(v), shifted)

  expect_error(envelope_to_voxel(e, edge_depths, depth_bounds = shifted[1:2, ]),
               class = "ocean3d_error_bounds")
  # A repeated depth with edges is an error, not two rows for one level.
  expect_error(
    envelope_to_voxel(e, c(0.5, 1.5, 1.5),
                      depth_bounds = cbind(c(0, 1, 1), c(1, 2, 2))),
    class = "ocean3d_error_axis_order"
  )
})

test_that("envelope_to_voxel() leaves out a slab the envelope only touches", {
  fp <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0,
                    ymax = 1)
  terra::values(fp) <- 1

  # [4, 5] starts where the deepest stored slab (2-4 m) ends: no shared water.
  below <- as_envelope(fp, depth_min = 4, depth_max = 5)
  expect_warning(
    v <- envelope_to_voxel(below, edge_depths, depth_bounds = edge_bounds),
    class = "ocean3d_warning_envelope_outside_slabs"
  )
  expect_true(all(is.na(terra::values(v))))
  # So a voxel and an envelope that only meet there share no volume.
  edged <- as_voxel(make_voxel(1, edge_depths), depth_bounds = edge_bounds)
  touching <- make_range_rast(rep(4, 9), rep(5, 9))
  expect_warning(
    overlap <- intersect_3d(touching, edged),
    class = "ocean3d_warning_envelope_outside_slabs"
  )
  # The converted envelope is empty, so the overlap is absent (NA) in every
  # cell, never a positive volume.
  expect_true(all(is.na(terra::values(volume(overlap, per_cell = TRUE)))))

  # An envelope in a gap between stored slabs gets the same warning, which
  # says it can fall between them.
  gappy <- cbind(c(0, 1.2, 2), c(1, 2, 4))
  between <- as_envelope(fp, depth_min = 1.05, depth_max = 1.15)
  expect_warning(
    envelope_to_voxel(between, edge_depths, depth_bounds = gappy),
    "above, below or between them"
  )
})

test_that("voxel_to_envelope() runs from the top to the bottom slab edge", {
  # Cell 1 occupies all three levels, cell 2 only 1.5 m, cell 3 1.5 m and 3 m
  # (a gap is filled), the rest none.
  vals <- cbind(c(1, NA, NA, rep(NA, 6)), c(1, 1, 1, rep(NA, 6)),
                c(1, NA, 1, rep(NA, 6)))
  v <- as_voxel(make_voxel(vals, edge_depths), depth_bounds = edge_bounds)

  e <- terra::values(voxel_to_envelope(v))
  expect_equal(unname(e[1:3, ]), rbind(c(0, 4), c(1, 2), c(1, 4)))
  expect_true(all(is.na(e[4:9, ])))

  # Without edges, the layer depths, as before.
  plain <- terra::values(voxel_to_envelope(make_voxel(vals, edge_depths)))
  expect_equal(unname(plain[1:3, ]), rbind(c(0.5, 3), c(1.5, 1.5), c(1.5, 3)))

  # With contiguous edges the round trip recovers every occupied level of the
  # cells without a gap (cell 3's gap is filled, as always).
  back <- envelope_to_voxel(voxel_to_envelope(v), edge_depths,
                            depth_bounds = edge_bounds)
  expect_identical(!is.na(terra::values(back))[1:2, ],
                   !is.na(terra::values(v))[1:2, ])
})

test_that("occupied() keeps a voxel's depth edges", {
  v <- as_voxel(make_voxel(seq_len(27), edge_depths),
                depth_bounds = edge_bounds)
  kept <- occupied(v, function(x) x > 0)
  expect_equal(edges_of(kept), edge_bounds)
  expect_equal(volume(kept), volume(v))
})
