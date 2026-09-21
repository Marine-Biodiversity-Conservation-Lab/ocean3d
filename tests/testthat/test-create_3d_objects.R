test_that("SpatVoxel and SpatEnvelope are SpatRasters and SpatVolumes", {
  v <- methods::new("SpatVoxel", make_multidepth_rast())

  expect_s4_class(v, "SpatVoxel")
  expect_s4_class(v, "SpatRaster")
  expect_s4_class(v, "SpatVolume")
  # terra operations still work on the subclass
  expect_identical(terra::nlyr(v), 4)
})

test_that("SpatVoxel validity rejects non-conforming layer names and order", {
  bad_names <- make_multidepth_rast()
  names(bad_names) <- c("a", "b", "c", "d")
  expect_error(methods::new("SpatVoxel", bad_names), "\\{variable\\}_depth=")

  unsorted <- make_multidepth_rast(depths = c(300, 200, 100, 0))
  expect_error(methods::new("SpatVoxel", unsorted), "shallow to deep")
})

test_that("SpatEnvelope validity requires exactly depth_min, depth_max", {
  expect_error(methods::new("SpatEnvelope", make_multidepth_rast()),
               "exactly: depth_min, depth_max")
})

# as_envelope() ----

test_that("as_envelope() builds a SpatEnvelope from a footprint and constant depths", {
  e <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)

  expect_s4_class(e, "SpatEnvelope")
  expect_s4_class(e, "SpatVolume")
  expect_named(e, c("depth_min", "depth_max"))

  vals <- terra::values(e)
  expect_identical(unname(vals[, "depth_min"]), c(0, 0, 0, NA))
  expect_identical(unname(vals[, "depth_max"]), c(200, 200, 200, NA))
})

test_that("as_envelope() reads only the non-NA pattern of the footprint", {
  # Values differ wildly, including a zero and a negative; presence does not.
  a <- as_envelope(make_footprint(c(1, 5, 0, NA)), 10, 20)
  b <- as_envelope(make_footprint(c(-7, 0.5, 999, NA)), 10, 20)

  expect_identical(terra::values(a), terra::values(b))
})

test_that("as_envelope() keeps the grid of the footprint", {
  fp <- make_footprint()
  e <- as_envelope(fp, 0, 200)

  expect_true(terra::compareGeom(e, fp, stopOnError = FALSE))
})

test_that("as_envelope() accepts per-cell depth limits as single-layer rasters", {
  fp <- make_footprint()
  dmin <- terra::setValues(terra::rast(fp), c(0, 50, 100, 150))
  dmax <- terra::setValues(terra::rast(fp), c(100, 200, 300, 400))

  vals <- terra::values(as_envelope(fp, depth_min = dmin, depth_max = dmax))

  expect_identical(unname(vals[, "depth_min"]), c(0, 50, 100, NA))
  expect_identical(unname(vals[, "depth_max"]), c(100, 200, 300, NA))
})

test_that("as_envelope() mixes a constant limit with a per-cell limit", {
  fp <- make_footprint()
  dmax <- terra::setValues(terra::rast(fp), c(100, 200, 300, 400))

  vals <- terra::values(as_envelope(fp, depth_min = 10, depth_max = dmax))

  expect_identical(unname(vals[, "depth_min"]), c(10, 10, 10, NA))
  expect_identical(unname(vals[, "depth_max"]), c(100, 200, 300, NA))
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

test_that("as_envelope() promotes voxel_to_envelope() output unchanged", {
  v <- methods::new("SpatVoxel", make_multidepth_rast())
  e <- voxel_to_envelope(v)

  expect_identical(terra::values(as_envelope(e)), terra::values(e))
})

test_that("as_envelope() rejects depth arguments alongside existing depth layers", {
  e <- as_envelope(make_footprint(), 0, 200)

  expect_error(as_envelope(e, depth_min = 0, depth_max = 100),
               "already has depth_min and depth_max")
})

test_that("as_envelope() requires depth limits for a bare footprint", {
  expect_error(as_envelope(make_footprint()), "are required")
})

test_that("as_envelope() rejects a multi-layer footprint", {
  expect_error(as_envelope(make_multidepth_rast(), 0, 200), "single-layer footprint")
})

test_that("as_envelope() rejects vector geometry with a pointer to vect_to_envelope()", {
  poly <- terra::vect("POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))")

  expect_error(as_envelope(poly, 0, 200), "vect_to_envelope")
})

test_that("as_envelope() rejects a non-raster x", {
  expect_error(as_envelope(1:10, 0, 200), "must be a SpatRaster")
})

test_that("as_envelope() rejects inverted depth limits", {
  expect_error(as_envelope(make_footprint(), depth_min = 200, depth_max = 0),
               "at least `depth_min`")

  fp <- make_footprint()
  dmax <- terra::setValues(terra::rast(fp), c(300, 300, 5, NA))
  expect_error(as_envelope(fp, depth_min = 10, depth_max = dmax),
               "at least `depth_min`")
})

test_that("as_envelope() rejects negative depths (positive-down convention)", {
  expect_error(as_envelope(make_footprint(), depth_min = -50, depth_max = 200),
               "positive metres increasing downward")
})

test_that("as_envelope() rejects depth limits that are not scalars or rasters", {
  fp <- make_footprint()

  expect_error(as_envelope(fp, depth_min = c(0, 100), depth_max = 200),
               "single non-missing number")
  expect_error(as_envelope(fp, depth_min = NA, depth_max = 200),
               "single non-missing number")
  expect_error(as_envelope(fp, depth_min = "0", depth_max = 200),
               "single non-missing number")
})

test_that("as_envelope() handles an empty footprint", {
  empty <- make_footprint(c(NA, NA, NA, NA))
  e <- as_envelope(empty, 0, 200)

  expect_s4_class(e, "SpatEnvelope")
  expect_true(all(is.na(terra::values(e))))
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

test_that("as_voxel() sorts layers shallow to deep instead of rejecting them", {
  # new() rejects unsorted input; the constructor repairs it
  unsorted <- make_multidepth_rast(depths = c(300, 200, 100, 0))
  expect_error(methods::new("SpatVoxel", unsorted), "shallow to deep")

  v <- as_voxel(unsorted)
  expect_s4_class(v, "SpatVoxel")
  expect_named(v, paste0("temp_depth=", c(0, 100, 200, 300)))
})

test_that("as_voxel() sorting carries the values with the layers", {
  unsorted <- make_multidepth_rast(
    depths = c(300, 0),
    vals = list(c(1, 2, 3, 4), c(5, 6, 7, 8))
  )
  v <- as_voxel(unsorted)

  expect_named(v, paste0("temp_depth=", c(0, 300)))
  # the 0 m layer must still hold the values it had before the sort
  expect_identical(unname(terra::values(v[["temp_depth=0"]])[, 1]), c(5, 6, 7, 8))
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

  expect_s4_class(as_voxel(lst), "SpatVoxel")
})

test_that("as_voxel() rejects non-conforming names and says how to fix it", {
  r <- make_multidepth_rast()
  names(r) <- c("a", "b", "c", "d")

  expect_error(as_voxel(r), "do not follow")
  expect_error(as_voxel(r), "Pass `depths`")
})

test_that("as_voxel() rejects negative depths rather than negating them", {
  expect_error(
    as_voxel(make_multidepth_rast(), depths = c(0, -100, -200, -300)),
    "positive metres increasing downward"
  )
})

test_that("as_voxel() rejects duplicate depths", {
  expect_error(as_voxel(make_multidepth_rast(), depths = c(0, 100, 100, 300)),
               "duplicate depth")
})

test_that("as_voxel() rejects a mismatched `depths` length", {
  expect_error(as_voxel(make_multidepth_rast(), depths = c(0, 100)),
               "one value per layer")
})

test_that("as_voxel() rejects a non-raster input", {
  expect_error(as_voxel("not a raster"), "must be a SpatRaster")
})

test_that("as_voxel() catches a SpatVoxel that a terra operation broke", {
  # terra propagates the class tag through `[[` without re-running validity,
  # so an object can be labelled SpatVoxel and not be one. The idempotent
  # fast path must not wave it through.
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = 3)
  terra::values(r) <- 1
  v <- as_voxel(r, depths = c(0, 100, 200), varname = "presence")
  expect_true(methods::validObject(v, test = TRUE))

  # selecting the first layer twice leaves depths 0, 0, 100
  broken <- v[[c(1, 1, 2)]]
  expect_s4_class(broken, "SpatVoxel")
  expect_false(isTRUE(methods::validObject(broken, test = TRUE)))

  # as_voxel() refuses it, with the message a plain raster would get
  expect_error(as_voxel(broken), "duplicate depth")

  # and the same when the depth is lost from the layer names
  renamed <- v
  names(renamed) <- c("a", "b", "c")
  expect_s4_class(renamed, "SpatVoxel")
  expect_error(as_voxel(renamed), "_depth=")

  # a lost name is recoverable by supplying `depths`, which rebuilds the names
  rebuilt <- as_voxel(renamed, depths = c(0, 100, 200), varname = "presence")
  expect_true(methods::validObject(rebuilt, test = TRUE))
  expect_named(rebuilt, names(v))
})

test_that("as_voxel() repairs a tagged SpatVoxel whose layers are out of order", {
  # package functions advise that as_voxel() can be used to repair SpatVoxel
  # so it should be able to repair when possible from the input invalid SpatVoxel
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

test_that("as_voxel() output is accepted by voxel_to_envelope()", {
  e <- voxel_to_envelope(as_voxel(make_multidepth_rast()))
  expect_s4_class(e, "SpatEnvelope")
})

# voxel_to_envelope() ---- 
test_that("voxel_to_envelope() gracefully rejects non SpatVoxel object as input param", {
  not_voxel <- make_multidepth_rast()
  expect_error(voxel_to_envelope(not_voxel))

  voxel <- as_voxel(not_voxel)
  expect_s4_class(voxel, "SpatVoxel")
})

test_that("voxel_to_envelope() rejects a non-function fun", {
  v <- methods::new("SpatVoxel", make_multidepth_rast())
  expect_error(voxel_to_envelope(v, "extent"), "must be a function")
})


test_that("voxel_to_envelope() defaults to the non-NA vertical extent", {
  v <- methods::new("SpatVoxel", make_multidepth_rast())
  e <- voxel_to_envelope(v)

  expect_s4_class(e, "SpatEnvelope")
  expect_named(e, c("depth_min", "depth_max"))

  vals <- terra::values(e)
  expect_identical(unname(vals[, "depth_min"]), c(0, 100, 0, NA))
  expect_identical(unname(vals[, "depth_max"]), c(200, 200, 300, NA))
})

test_that("voxel_to_envelope() fills interior gaps (lossy collapse)", {
  # cell 1 has values at 0 m and 200 m but not at 100 m; the envelope is solid
  v <- methods::new("SpatVoxel", make_multidepth_rast())
  vals <- terra::values(voxel_to_envelope(v))

  expect_identical(unname(vals[1, "depth_min"]), 0)
  expect_identical(unname(vals[1, "depth_max"]), 200)
})

test_that("voxel_to_envelope() bounds the depths where the predicate holds", {
  v <- methods::new("SpatVoxel", make_multidepth_rast())
  vals <- terra::values(voxel_to_envelope(v, function(x) x > 15))

  # only cell 2 exceeds 15, at 100 m and 200 m
  expect_identical(unname(vals[, "depth_min"]), c(NA, 100, NA, NA))
  expect_identical(unname(vals[, "depth_max"]), c(NA, 200, NA, NA))
})

test_that("voxel_to_envelope() returns NA where the predicate never holds", {
  v <- methods::new("SpatVoxel", make_multidepth_rast())
  vals <- terra::values(voxel_to_envelope(v, function(x) x > 1e6))

  expect_true(all(is.na(vals)))
})

test_that("voxel_to_envelope() handles a single-depth voxel", {
  # a one-layer voxel collapses to a zero-thickness envelope at that depth
  v <- methods::new("SpatVoxel",
                    make_multidepth_rast(depths = 50,
                                         vals = list(c(10, NA, 5, NA))))
  vals <- terra::values(voxel_to_envelope(v))

  expect_identical(unname(vals[, "depth_min"]), c(50, NA, 50, NA))
  expect_identical(unname(vals[, "depth_max"]), c(50, NA, 50, NA))
})

# envelope_to_voxel() ----
test_that("envelope_to_voxel() rejects invalid input param types", {
  # a multi-depth raster is a plausible-looking but wrong input: it is the
  # voxel form of a 3D domain, not the envelope form
  bad_envel <- make_multidepth_rast()
  good_envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)

  bad_depths <- c("hello", 100, "depths")
  good_depths <- c(0, 50, 100, 150, 200, 300, 600)

  # check for param x valid SpatEnvelope
  expect_error(
    envelope_to_voxel(x = bad_envel, depths = good_depths),
    "`x` needs to be of `SpatEnvelope` class",
    fixed = TRUE
  )

  # check for param depths valid numeric
  expect_error(
    envelope_to_voxel(x = good_envel, depths = bad_depths),
    "`depths` needs to be array coercible to numeric type",
    fixed = TRUE
  )
  expect_error(
    envelope_to_voxel(x = good_envel, depths = numeric(0)),
    "`depths` needs to be array coercible to numeric type",
    fixed = TRUE
  )

  # check for param profile a function. A shape is passed as the function
  # itself, so a name that used to select one is now just the wrong type.
  expect_error(
    envelope_to_voxel(x = good_envel, depths = good_depths, profile = "equal"),
    "`profile` needs to be NULL or a function of (ind, depths, n_depths); got character",
    fixed = TRUE
  )
  expect_error(
    envelope_to_voxel(x = good_envel, depths = good_depths, profile = 100),
    "`profile` needs to be NULL or a function of (ind, depths, n_depths); got numeric",
    fixed = TRUE
  )
  expect_error(
    envelope_to_voxel(x = good_envel, depths = good_depths, profile = NA),
    "`profile` needs to be NULL or a function",
    fixed = TRUE
  )

  # the supplied shapes and the default all have to survive the same check
  expect_no_error(
    envelope_to_voxel(x = good_envel, depths = good_depths, profile = NULL)
  )
  expect_no_error(
    envelope_to_voxel(x = good_envel, depths = good_depths, profile = profile_flat)
  )
  expect_no_error(
    envelope_to_voxel(x = good_envel, depths = good_depths, profile = profile_equal)
  )

  # check for param varname valid string
  expect_error(
    envelope_to_voxel(x = good_envel, depths = good_depths, varname = 100),
    "`varname` needs to be a string",
    fixed = TRUE
  )
  expect_error(
    envelope_to_voxel(x = good_envel, depths = good_depths,
                      varname = c("presence", "absence")),
    "`varname` needs to be a string",
    fixed = TRUE
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

test_that("envelope_to_voxel() generates voxel from envelope", {
  temp_r <- make_footprint(c(1,1,1,1))
  # Create 2x2 simple envelope with depth_min and depth_max layers. 
  envel <- as_envelope(temp_r, depth_min = 70, depth_max = 210)

  # Defined depth layers 
  a_depths <- c(0, 50, 100, 150, 200, 300, 400)

  # Expected result. Each level stands for the slab running down to the next,
  # so [70, 210] occupies every slab it overlaps: it starts inside 50-100, so
  # 50 is occupied though 70 > 50, and ends inside 200-300, so 200 is occupied
  # and 300 is not.
  occupied <- c(50, 100, 150, 200)
  expected_vox <- lapply(a_depths, function(x) {
    if (x %in% occupied) {
      r <- terra::setValues(temp_r, 1L)
    } else {
      r <- terra::setValues(temp_r, NA_integer_)
    }
    names(r) <- paste0("test_depth=", x)
    r
  }) %>% as_voxel()

  got <- envelope_to_voxel(x = envel, depths = a_depths, varname = "test")

  expect_s4_class(got, "SpatVoxel")
  expect_named(got, names(expected_vox))
  expect_identical(terra::values(got), terra::values(expected_vox))
})

test_that("envelope_to_voxel() catches an envelope that falls between levels", {
  # the case the point-in-interval rule missed: [10, 20] contains none of the
  # levels, but lies inside the slab the 0 m level stands for
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 10, depth_max = 20)

  expect_no_warning(v <- envelope_to_voxel(envel, depths = c(0, 100)))
  expect_identical(unname(terra::values(v)[1, ]), c(1L, NA_integer_))
})

test_that("envelope_to_voxel() occupies every slab the envelope overlaps", {
  # [70, 210] starts inside the 50-100 slab and ends inside 200-300
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 70, depth_max = 210)
  vals <- terra::values(envelope_to_voxel(envel, depths = c(0, 50, 100, 200, 300)))

  expect_identical(unname(vals[1, ]), c(NA_integer_, 1L, 1L, 1L, NA_integer_))
})

test_that("envelope_to_voxel() gives a boundary depth to the slab it heads", {
  # an envelope sitting wholly inside one slab occupies that slab alone; the
  # level below owns its own top, so it is not dragged in
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 60, depth_max = 99)
  vals <- terra::values(envelope_to_voxel(envel, depths = c(0, 50, 100, 200)))

  expect_identical(unname(vals[1, ]), c(NA_integer_, 1L, NA_integer_, NA_integer_))

  # extending it to end exactly at 100 still leaves that level out: the
  # envelope only touches the 100-200 slab and shares no water with it
  envel2 <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 60, depth_max = 100)
  vals2 <- terra::values(envelope_to_voxel(envel2, depths = c(0, 50, 100, 200)))

  expect_identical(unname(vals2[1, ]), c(NA_integer_, 1L, NA_integer_, NA_integer_))

  # reaching past 100 by any amount brings that level in
  envel3 <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 60, depth_max = 101)
  vals3 <- terra::values(envelope_to_voxel(envel3, depths = c(0, 50, 100, 200)))

  expect_identical(unname(vals3[1, ]), c(NA_integer_, 1L, 1L, NA_integer_))
})

test_that("envelope_to_voxel() puts adjacent depth ranges on disjoint levels", {
  # touching is not overlapping: [0, 100] and [100, 200] share the boundary at
  # 100 m and nothing else, so no level is occupied by both. The levels must
  # put a slab edge at 100 m for this to be expressible at all: under "top"
  # that is the 100 level itself, under "midpoint" it is halfway between 50
  # and 150.
  fp <- make_footprint(c(1, 1, 1, 1))
  upper <- as_envelope(fp, depth_min = 0, depth_max = 100)
  lower <- as_envelope(fp, depth_min = 100, depth_max = 200)

  for (bounds in c("top", "midpoint")) {
    depths <- if (bounds == "top") c(0, 100, 200) else c(0, 50, 150, 200)
    u <- terra::values(envelope_to_voxel(upper, depths, bounds = bounds))
    l <- terra::values(envelope_to_voxel(lower, depths, bounds = bounds))
    both <- !is.na(u) & !is.na(l)
    expect_false(any(both), info = bounds)
  }

  # under "top" the shared level is claimed by the deeper range alone. The
  # deepest level stands for a zero-thickness slab, so a range that ends
  # exactly on it is not recorded there either: [100, 200] holds 100 only.
  u <- terra::values(envelope_to_voxel(upper, c(0, 100, 200)))
  l <- terra::values(envelope_to_voxel(lower, c(0, 100, 200)))
  expect_identical(unname(u[1, ]), c(1L, NA, NA))
  expect_identical(unname(l[1, ]), c(NA, 1L, NA))

  # under "midpoint" the 100 m boundary falls inside the 100 m slab when the
  # levels are evenly spaced, so both ranges genuinely share water with that
  # slab and both occupy it. That is the discretisation, not a touch.
  depths <- c(0, 50, 100, 150, 200)
  u <- terra::values(envelope_to_voxel(upper, depths, bounds = "midpoint"))
  l <- terra::values(envelope_to_voxel(lower, depths, bounds = "midpoint"))
  expect_identical(unname(u[1, ]), c(1L, 1L, 1L, NA, NA))
  expect_identical(unname(l[1, ]), c(NA, NA, 1L, 1L, 1L))
})

test_that("envelope_to_voxel() keeps every level the old point rule found", {
  # overlap is a widening of the rule, never a narrowing: a level whose own
  # value lies inside the envelope is still occupied
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 50, depth_max = 250)
  depths <- c(0, 50, 100, 200, 300)
  vals <- terra::values(envelope_to_voxel(envel, depths = depths))

  inside <- depths >= 50 & depths <= 250
  expect_false(any(is.na(unname(vals[1, ]))[inside]))
})

test_that("envelope_to_voxel() treats the deepest level as a point, not a floor", {
  # an open-ended deepest slab would swallow any range below the grid and
  # report it as resolved; keeping it a point leaves the gap visible
  below <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 310, depth_max = 400)
  expect_warning(v <- envelope_to_voxel(below, depths = c(0, 100, 200, 300)),
                 "lying entirely outside")
  expect_true(all(is.na(terra::values(v))))

  # reaching the deepest level is enough to be recorded at it
  touching <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 290, depth_max = 400)
  vals <- terra::values(envelope_to_voxel(touching, depths = c(0, 100, 200, 300)))
  expect_identical(unname(vals[1, ]), c(NA, NA, 1L, 1L))
})

test_that(".depth_slabs() runs each slab to the next level under 'top'", {
  slab <- .depth_slabs(c(0, 100, 200), "top")

  expect_identical(slab$lower, c(0, 100, 200))
  # the deepest level has no next one and collapses to a point
  expect_identical(slab$upper, c(100, 200, 200))
})

test_that(".depth_slabs() halves the gaps under 'midpoint', as SPEC.md's WOA does", {
  # SPEC.md: WOA's 0 m layer has bounds 0-2.5, the 5 m layer 2.5-7.5, the 10 m
  # layer 7.5-12.5
  slab <- .depth_slabs(c(0, 5, 10, 15), "midpoint")

  expect_identical(slab$lower, c(0, 2.5, 7.5, 12.5))
  expect_identical(slab$upper, c(2.5, 7.5, 12.5, 15))
})

test_that(".depth_slabs() clamps the outer edges rather than extrapolating", {
  # half a gap above WOA's shallowest level would be a negative depth, and half
  # a gap below the deepest would invent reach the levels never had
  for (bounds in c("top", "midpoint")) {
    slab <- .depth_slabs(c(0, 100, 200), bounds)
    expect_identical(min(slab$lower), 0, info = bounds)
    expect_identical(max(slab$upper), 200, info = bounds)
  }
})

test_that(".depth_slabs() collapses a lone level to a point either way", {
  for (bounds in c("top", "midpoint")) {
    slab <- .depth_slabs(100, bounds)
    expect_identical(slab$lower, 100, info = bounds)
    expect_identical(slab$upper, 100, info = bounds)
  }
})

test_that(".depth_slabs() tiles the depth span without gaps or overlaps", {
  # each slab picks up exactly where the last left off, so the levels partition
  # [min(depths), max(depths)]
  for (bounds in c("top", "midpoint")) {
    slab <- .depth_slabs(c(0, 50, 100, 400), bounds)
    expect_identical(slab$lower[-1], slab$upper[-length(slab$upper)], info = bounds)
  }
})

test_that("envelope_to_voxel() reads an envelope against the chosen bounds", {
  # [60, 70] sits in the 0-100 slab under 'top', but in the 50-150 slab of the
  # 100 m level once the edges move to the midpoints
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 60, depth_max = 70)
  depths <- c(0, 100, 200)

  top <- terra::values(envelope_to_voxel(envel, depths = depths))
  mid <- terra::values(envelope_to_voxel(envel, depths = depths,
                                         bounds = "midpoint"))

  expect_identical(unname(top[1, ]), c(1L, NA, NA))
  expect_identical(unname(mid[1, ]), c(NA, 1L, NA))
})

test_that("envelope_to_voxel() defaults to the 'top' bounds", {
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 60, depth_max = 210)
  depths <- c(0, 100, 200, 300)

  expect_identical(
    terra::values(envelope_to_voxel(envel, depths = depths)),
    terra::values(envelope_to_voxel(envel, depths = depths, bounds = "top"))
  )
})

test_that("envelope_to_voxel() still reaches the deepest level under 'midpoint'", {
  # the deepest slab is closed at max(depths), so an envelope sitting exactly
  # on the last level is caught rather than falling off the bottom
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 300, depth_max = 300)
  vals <- terra::values(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300), bounds = "midpoint")
  )

  expect_identical(unname(vals[1, ]), c(NA, NA, NA, 1L))
})

test_that("envelope_to_voxel() warns outside the depth span under either bounds", {
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 500, depth_max = 600)

  for (bounds in c("top", "midpoint")) {
    expect_warning(
      envelope_to_voxel(envel, depths = c(0, 100, 200), bounds = bounds),
      "lying entirely outside"
    )
  }
})

test_that("envelope_to_voxel() rejects an unknown bounds convention", {
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)

  expect_error(
    envelope_to_voxel(envel, depths = c(0, 100), bounds = "gaussian"),
    "should be one of"
  )
})

test_that("envelope_to_voxel() occupies the level at its top limit, not at its bottom", {
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 50, depth_max = 200)
  vals <- terra::values(envelope_to_voxel(envel, depths = c(0, 50, 100, 200, 300)))

  # 50 heads the slab 50-100, which the envelope runs through, so it is
  # occupied. 200 heads the slab 200-300, which the envelope only touches, so
  # it is not: touching is not overlapping.
  expect_identical(unname(vals[1, ]), c(NA, 1L, 1L, NA, NA))
})

test_that("envelope_to_voxel() leaves cells absent from the envelope empty", {
  # cell 4 of make_footprint() is NA, so it is NA at every depth
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)
  vals <- terra::values(envelope_to_voxel(envel, depths = c(0, 100, 200)))

  expect_identical(unname(vals[, 1]), c(1L, 1L, 1L, NA))
  expect_true(all(is.na(vals[4, ])))
})

test_that("envelope_to_voxel() honours per-cell depth limits", {
  fp <- make_footprint()
  dmin <- terra::setValues(terra::rast(fp), c(0, 100, 200, 0))
  dmax <- terra::setValues(terra::rast(fp), c(100, 300, 300, 300))
  envel <- as_envelope(fp, depth_min = dmin, depth_max = dmax)

  vals <- terra::values(envelope_to_voxel(envel, depths = c(0, 100, 200, 300)))

  # each range holds the slabs it runs through and not the one it ends on
  expect_identical(unname(vals[1, ]), c(1L, NA, NA, NA))   # [0, 100]
  expect_identical(unname(vals[2, ]), c(NA, 1L, 1L, NA))    # [100, 300]
  expect_identical(unname(vals[3, ]), c(NA, NA, 1L, NA))   # [200, 300]
  expect_true(all(is.na(vals[4, ])))                  # absent footprint
})

test_that("envelope_to_voxel() writes a constant `values` at every occupied level", {
  fp <- make_footprint()
  dmax <- terra::setValues(terra::rast(fp), c(150, 350, 350, 350))
  envel <- as_envelope(fp, depth_min = 0, depth_max = dmax)

  vals <- terra::values(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300), values = 5,
                      varname = "score")
  )

  # same occupancy as presence, carrying 5 instead of 1
  expect_identical(unname(vals[1, ]), c(5L, 5L, NA, NA))
  expect_identical(unname(vals[2, ]), rep(5L, 4))
})

test_that("envelope_to_voxel() writes a per-cell `values` raster", {
  fp <- make_footprint(c(1, 1, 1, 1))
  envel <- as_envelope(fp, depth_min = 50, depth_max = 250)
  effort <- terra::setValues(terra::rast(fp), c(10, 20, 30, 40))

  vals <- terra::values(
    envelope_to_voxel(envel, depths = c(0, 50, 100, 200, 300), values = effort,
                      varname = "effort")
  )

  expect_identical(unname(vals[1, ]), c(NA, 10, 10, 10, NA))
  expect_identical(unname(vals[4, ]), c(NA, 40, 40, 40, NA))
})

test_that("envelope_to_voxel() profile_equal divides across the occupied depths", {
  fp <- make_footprint()
  dmax <- terra::setValues(terra::rast(fp), c(150, 350, 350, 350))
  envel <- as_envelope(fp, depth_min = 0, depth_max = dmax)

  # the two distinct intervals occupy 2 and 4 levels, so get different shares
  v <- envelope_to_voxel(envel, depths = c(0, 100, 200, 300),
                         profile = profile_equal, varname = "time")
  vals <- terra::values(v)

  expect_named(v, paste0("time_depth=", c(0, 100, 200, 300)))
  expect_identical(unname(vals[1, ]), c(0.5, 0.5, NA, NA))
  expect_identical(unname(vals[2, ]), rep(0.25, 4))
})

test_that("envelope_to_voxel() profile_equal conserves each cell's value", {
  fp <- make_footprint(c(1, 1, 1, 1))
  dmax <- terra::setValues(terra::rast(fp), c(150, 250, 350, 350))
  envel <- as_envelope(fp, depth_min = 0, depth_max = dmax)
  effort <- terra::setValues(terra::rast(fp), c(10, 20, 30, 40))

  vals <- terra::values(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300), values = effort,
                      profile = profile_equal, varname = "effort")
  )

  # cells occupy 2, 3, 4 and 4 levels respectively; each still sums to its input
  expect_identical(unname(vals[1, ]), c(5, 5, NA, NA))
  expect_identical(unname(rowSums(vals, na.rm = TRUE)), c(10, 20, 30, 40))
})

test_that("envelope_to_voxel() without a profile does not divide the value", {
  fp <- make_footprint(c(1, 1, 1, 1))
  envel <- as_envelope(fp, depth_min = 0, depth_max = 350)
  effort <- terra::setValues(terra::rast(fp), c(10, 20, 30, 40))

  vals <- terra::values(
    envelope_to_voxel(envel, depths = c(0, 100, 200, 300), values = effort)
  )

  # full value at all four levels: a footprint map, not a conserved total
  expect_identical(unname(vals[1, ]), rep(10, 4))
  expect_identical(unname(rowSums(vals, na.rm = TRUE)), c(40, 80, 120, 160))
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

test_that("envelope_to_voxel() rejects an unusable `values`", {
  fp <- make_footprint(c(1, 1, 1, 1))
  envel <- as_envelope(fp, depth_min = 0, depth_max = 300)

  expect_error(
    envelope_to_voxel(envel, depths = c(0, 100),
                      values = terra::rast(nrows = 3, ncols = 3)),
    "same grid", fixed = TRUE
  )
  expect_error(
    envelope_to_voxel(envel, depths = c(0, 100),
                      values = c(terra::rast(fp), terra::rast(fp))),
    "must be a single-layer SpatRaster", fixed = TRUE
  )
  expect_error(
    envelope_to_voxel(envel, depths = c(0, 100), values = c(1, 2)),
    "single non-missing number", fixed = TRUE
  )
})

test_that("envelope_to_voxel() sorts and deduplicates the requested depths", {
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 0, depth_max = 300)

  v <- envelope_to_voxel(envel, depths = c(200, 0, 100, 200))

  expect_named(v, paste0("presence_depth=", c(0, 100, 200)))
})

test_that("envelope_to_voxel() rejects negative depths (positive-down convention)", {
  envel <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)

  expect_error(envelope_to_voxel(envel, depths = c(-100, 0, 100)),
               "positive metres increasing downward")
})

test_that("envelope_to_voxel() warns when an envelope lies below every depth level", {
  # the slabs reach only as deep as max(depths), so nothing can hold this cell
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 500, depth_max = 600)

  expect_warning(v <- envelope_to_voxel(envel, depths = c(0, 100, 200)),
                 "lying entirely outside")
  expect_true(all(is.na(terra::values(v))))
})

test_that("envelope_to_voxel() warns when an envelope lies above every depth level", {
  # the shallowest level is the top of the first slab, so a cell entirely
  # above it is out of reach in the same way
  envel <- as_envelope(make_footprint(c(1, 1, 1, 1)), depth_min = 5, depth_max = 40)

  expect_warning(v <- envelope_to_voxel(envel, depths = c(50, 100, 200)),
                 "lying entirely outside")
  expect_true(all(is.na(terra::values(v))))
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

  back <- envelope_to_voxel(voxel_to_envelope(v), depths = depths, varname = "temp")

  # voxel_to_envelope() reports each cell's deepest occupied level as its
  # depth_max, and under "top" a level names the top of its slab. The envelope
  # therefore only touches that last slab, and touching is not overlapping, so
  # the round trip recovers every occupied level except the deepest one per
  # cell. The values themselves are not recovered either, since an envelope
  # carries only the depth limits.
  expect_named(back, names(v))
  present <- !is.na(terra::values(back))
  expect_identical(unname(present[1, ]), c(TRUE, TRUE, FALSE, FALSE))   # 0-200 -> 0, 100
  expect_identical(unname(present[2, ]), c(FALSE, TRUE, TRUE, FALSE))   # 100-300 -> 100, 200
  expect_identical(unname(present[3, ]), c(TRUE, TRUE, TRUE, FALSE))    # 0-300 -> 0, 100, 200
  expect_identical(unname(present[4, ]), rep(FALSE, 4))                 # absent
})

# vect_to_envelope() ----

test_that("vect_to_envelope() checks for correct params", {
  expect_error(
    vect_to_envelope("test", make_template(), 0, 10),
    "`polygon` needs to be of class SpatVector, sf, or sfc"
  )
  expect_error(
    vect_to_envelope(make_polygon(), "test", 0, 10),
    "`template` needs to be of class SpatRaster"
  )
})

test_that("vect_to_envelope() provides meaningful message when polygon and template don't overlap", {
  # Create raster that doesn't overlap with make_polygon()
  template_r <- rast(
    nrows = 5, ncols = 5, 
    xmin = 5, xmax = 10, ymin = 5, ymax = 10, crs = "EPSG:4326"
  )
  expect_error(
    vect_to_envelope(make_polygon(), template_r, 0, 10),
    "All output cell values are NA. `polygon` and `template` may not spatially overlap."
  )
})

test_that("vect_to_envelope() checks for agreement of CRS between inputs", {
  template_r <- rast(
    nrows = 5, ncols = 5, 
    xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "EPSG:3857"
  )
  expect_error(
    vect_to_envelope(make_polygon(), template_r, 0, 10),
    "`polygon` and `template` have different CRS."
  )

  # TODO: implement check for the depth_min and depth_max inputs
})

test_that("vect_to_envelope() takes single numeric value for depth_min and depth_max", {
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

  result <- vect_to_envelope(make_polygon(), make_template(), depth_min = 10, depth_max = 20) 

  expect_true(terra::identical(result, expected_result)) # nolint: expect_identical_linter. terra objects are pointers, not identical according to expect_identical()
})

test_that("vect_to_envelope() correctly takes deeper minimum depth in the depth_min list params", {
  r <- make_template()
  terra::values(r) <- seq_along(values(r))

  expected_result <- r %>%
    mask(make_polygon()) %>%
    max(10) 
  names(expected_result) <- "depth_min"
  expected_result$depth_max <- rast(r, vals = 25)
  expected_result <- as_envelope(expected_result)

  result <- vect_to_envelope(make_polygon(), make_template(), depth_min = c(10, r), depth_max = 25) 

  expect_true(terra::identical(result["depth_min"], expected_result["depth_min"])) # nolint: expect_identical_linter. terra objects are pointers, not identical according to expect_identical()
})

test_that("vect_to_envelope() correctly takes shallower maximum depth in the depth_max list params", {
  r <- make_template()
  terra::values(r) <- seq_along(values(r))

  expected_result <- rast()
  expected_result$depth_min <- rast(r, vals = 0)
  expected_result$depth_max <- r %>%
    mask(make_polygon()) %>%
    min(10) 
  expected_result <- as_envelope(expected_result)

  result <- vect_to_envelope(make_polygon(), make_template(), depth_min = 0, depth_max = c(10, r)) 

  expect_true(terra::identical(result["depth_max"], expected_result["depth_max"])) # nolint: expect_identical_linter. terra objects are pointers, not identical according to expect_identical()
})

# TODO: deal with case where depth_max is shallower than depth_min
test_that("vect_to_envelope() fills NA where depth_max is shallower than depth_min", {
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

  result <- vect_to_envelope(make_polygon(), make_template(), depth_min = 10, depth_max = c(16, r)) 

  expect_true(terra::identical(result, expected_result)) # nolint: expect_identical_linter. terra objects are pointers, not identical according to expect_identical()
})

test_that("vect_to_envelope() catches case where all cells depth_min are deeper than depth_max", {
  r <- make_template()
  terra::values(r) <- seq_along(values(r))

  expected_result <- r %>%
    mask(make_polygon()) %>%
    min(10) 
  names(expected_result) <- "depth_max"

  expect_warning(
    vect_to_envelope(make_polygon(), make_template(), depth_min = 50, depth_max = c(10, r)), 
    "All depth_max values are shallower than depth_min. Check `depth_min` and `depth_max`, you may have swapped these two."
  )
})

test_that("vect_to_envelope() accepts SpatVector and sfc as well as sf", {
  expect_s4_class(
    vect_to_envelope(terra::vect(make_polygon()), make_template(), 0, 10),
    "SpatEnvelope"
  )
  expect_s4_class(
    vect_to_envelope(sf::st_geometry(make_polygon()), make_template(), 0, 10),
    "SpatEnvelope"
  )
})

test_that("vect_to_envelope() returns a valid SpatEnvelope, visibly", {
  out <- vect_to_envelope(make_polygon(), make_template(), 0, 10)

  expect_s4_class(out, "SpatEnvelope")
  expect_true(isTRUE(methods::validObject(out)))
  expect_named(out, c("depth_min", "depth_max"))
  # The constructor assigns to `out` last; without a trailing bare `out` the
  # result would come back invisibly and never print at the console.
  expect_true(withVisible(vect_to_envelope(make_polygon(), make_template(), 0, 10))$visible)
})

# Ported from the retired voxelize_range(): the seafloor is now just another
# `depth_max` constraint rather than a dedicated `bathymetry` argument.
test_that("vect_to_envelope() clamps depth_max to the seafloor where it is shallower", {
  seafloor <- make_template(seq(10, 250, length.out = 25))

  out <- vect_to_envelope(make_polygon(), make_template(),
                          depth_min = 50, depth_max = list(200, seafloor))

  got <- terra::values(out$depth_max)
  got <- got[!is.na(got)]
  expect_true(all(got <= 200))
  expect_true(any(got < 200))
})

# Ported from the retired voxelize_range(): cells whose bed sits above the shallowest depth
# the animal occupies carry no envelope at all.
test_that("vect_to_envelope() drops cells where the seafloor is shallower than depth_min", {
  seafloor <- make_template(seq(10, 250, length.out = 25))

  out <- vect_to_envelope(make_polygon(), make_template(),
                          depth_min = 100, depth_max = list(200, seafloor))

  kept <- !is.na(terra::values(out$depth_min))
  expect_true(any(kept))
  expect_true(all(terra::values(seafloor)[kept] > 100))
})

# Ported from the retired voxelize_range()'s bathymetry CRS / resolution checks, which now
# apply to every SpatRaster in a depth list.
test_that("vect_to_envelope() errors on a depth raster that does not align with template", {
  bad_res <- terra::rast(nrows = 9, ncols = 9, xmin = 0, xmax = 5,
                         ymin = 0, ymax = 5, vals = 50, crs = "EPSG:4326")
  expect_error(
    vect_to_envelope(make_polygon(), make_template(), 0, list(10, bad_res)),
    "does not align with `template`"
  )

  bad_ext <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 9,
                         ymin = 0, ymax = 9, vals = 50, crs = "EPSG:4326")
  expect_error(
    vect_to_envelope(make_polygon(), make_template(), 0, list(10, bad_ext)),
    "does not align with `template`"
  )
})

test_that("vect_to_envelope() errors on a depth raster with a different CRS", {
  # Same extent numbers, different CRS: terra otherwise returns a wrong answer
  # here rather than failing.
  bad_crs <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5,
                         ymin = 0, ymax = 5, vals = 50, crs = "EPSG:3857")
  expect_error(
    vect_to_envelope(make_polygon(), make_template(), 0, list(10, bad_crs)),
    "different CRS than `template`"
  )
})

test_that("vect_to_envelope() rejects empty and unusable depth inputs", {
  expect_error(
    vect_to_envelope(make_polygon(), make_template(), list(), 10),
    "`depth_min` is empty"
  )
  expect_error(
    vect_to_envelope(make_polygon(), make_template(), 0, list()),
    "`depth_max` is empty"
  )
  expect_error(
    vect_to_envelope(make_polygon(), make_template(), 0, list(10, "deep")),
    "must be numeric or a SpatRaster"
  )
  # `na.rm = TRUE` would drop an NA silently, widening rather than narrowing.
  expect_error(
    vect_to_envelope(make_polygon(), make_template(), 0, list(10, NA_real_)),
    "is NA"
  )
})

test_that("vect_to_envelope() returns NA values for cells where input raster depth_max is NA", {
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

# Masking a voxel by an envelope is the package's "restrict a 3D raster to a
# species' per-cell depth window" operation:
#   terra::mask(v, envelope_to_voxel(e, depths(v)))
# It replaced the retired extract_rast_range(), so its behaviour is pinned here.
test_that("masking a voxel by an envelope respects the per-cell depth window", {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                   nlyrs = 3)
  terra::values(r) <- cbind(rep(1, 4), rep(2, 4), rep(3, 4))
  v <- as_voxel(r, depths = c(0, 100, 500), varname = "tan")

  # Cells 1-2 end at 100 m; cells 3-4 end at 500 m.
  fp <- terra::rast(v[[1]])
  terra::values(fp) <- 1
  e <- as_envelope(fp, depth_min = 0,
                   depth_max = terra::setValues(terra::rast(fp),
                                                c(100, 100, 500, 500)))

  out <- terra::mask(v, envelope_to_voxel(e, depths(v)))

  # Layer names survive the mask, so callers can still index by depth.
  expect_named(out, names(v))
  # Every cell holds the 0 m slab. Only the last two run into the 100-500
  # slab; the first two end where it begins and touching is not overlapping.
  # The 500 m level heads a zero-thickness slab that no cell reaches below.
  expect_identical(sum(!is.na(terra::values(out[["tan_depth=0"]]))), 4L)
  expect_identical(sum(!is.na(terra::values(out[["tan_depth=100"]]))), 2L)
  expect_identical(sum(!is.na(terra::values(out[["tan_depth=500"]]))), 0L)
  # Values themselves are untouched where the cell is in window.
  expect_identical(unique(terra::values(out[["tan_depth=0"]])[, 1]), 1)
})

test_that("masking keeps a cell whose envelope falls between two depth levels", {
  # A [10, 20] envelope contains none of c(0, 100, 200). Slab overlap places it
  # in the 0 m level; the retired point-containment test dropped the cell
  # entirely, which is the bug this idiom fixes.
  r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1,
                   nlyrs = 3)
  terra::values(r) <- cbind(c(1, 1), c(2, 2), c(3, 3))
  v <- as_voxel(r, depths = c(0, 100, 200), varname = "tan")

  fp <- terra::rast(v[[1]])
  terra::values(fp) <- 1
  e <- as_envelope(fp, depth_min = 10, depth_max = 20)

  out <- terra::mask(v, envelope_to_voxel(e, depths(v)))

  expect_identical(sum(!is.na(terra::values(out[["tan_depth=0"]]))), 2L)
  expect_true(all(is.na(terra::values(out[["tan_depth=100"]]))))
  expect_true(all(is.na(terra::values(out[["tan_depth=200"]]))))
})

test_that("masking drops cells that are NA in the envelope", {
  r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1,
                   nlyrs = 2)
  terra::values(r) <- cbind(c(1, 1), c(2, 2))
  v <- as_voxel(r, depths = c(0, 100), varname = "tan")

  fp <- terra::rast(v[[1]])
  terra::values(fp) <- c(1, NA)     # second cell absent
  e <- as_envelope(fp, depth_min = 0, depth_max = 150)

  out <- terra::mask(v, envelope_to_voxel(e, depths(v)))

  expect_identical(sum(!is.na(terra::values(out))), 2L)  # one cell x two depths
})

# ---- occupied() -------------------------------------------------------------

# Local helper: a 1x3 voxel over two standard depths, named by variable.
make_named_voxel <- function(vals, varname, depths = c(0, 100)) {
  r <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                   nlyrs = length(depths))
  terra::values(r) <- vals
  as_voxel(r, depths = depths, varname = varname)
}

test_that("occupied() turns a variable voxel into a 1-or-NA presence voxel", {
  v <- as_voxel(make_multidepth_rast())

  out <- occupied(v)

  expect_s4_class(out, "SpatVoxel")
  expect_named(out, paste0("presence_depth=", c(0, 100, 200, 300)))
  expect_identical(depths(out), depths(v))
  # Presence is 1 where the variable was recorded, NA where it was not.
  expect_identical(unname(terra::values(out)),
               unname(ifelse(is.na(terra::values(v)), NA, 1)))
})

test_that("occupied() applies the predicate one depth layer at a time", {
  v <- as_voxel(make_multidepth_rast())

  out <- occupied(v, function(x) x > 10)

  vals <- terra::values(v)
  expected <- ifelse(!is.na(vals) & vals > 10, 1, NA)
  expect_identical(unname(terra::values(out)), unname(expected))
})

test_that("occupied() is idempotent on a presence voxel", {
  v <- as_voxel(make_multidepth_rast())
  once <- occupied(v)

  expect_identical(terra::values(occupied(once)), terra::values(once))
  expect_named(occupied(once), names(once))
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

test_that("SpatVoxel validity rejects duplicate depths, as as_voxel() does", {
  raw <- make_multidepth_rast()

  # as_voxel() has always refused duplicates.
  expect_error(as_voxel(raw[[c(1, 1, 2)]]), "duplicate depth")

  # terra keeps the class tag through subsetting, so an invalid voxel can
  # still reach the verbs. Validity must catch what the constructor refuses.
  dup <- as_voxel(raw)[[c(1, 1, 2)]]
  expect_s4_class(dup, "SpatVoxel")
  expect_false(isTRUE(methods::validObject(dup, test = TRUE)))

  # Otherwise volume() counts the repeated slab and returns a wrong number.
  expect_error(volume(dup), "not a valid one")
})

test_that("SpatEnvelope validity rejects inverted and negative depth intervals", {
  e <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)

  # terra keeps the class tag through subsetting, so swapping the two layers
  # yields a tagged envelope whose interval runs backwards.
  inv <- e[[c("depth_max", "depth_min")]]
  names(inv) <- c("depth_min", "depth_max")
  expect_s4_class(inv, "SpatEnvelope")
  expect_false(isTRUE(methods::validObject(inv, test = TRUE)))

  # Arithmetic keeps the tag too, and can lift depths above the surface.
  expect_false(isTRUE(methods::validObject(e - 500, test = TRUE)))

  # as_envelope() has always refused both.
  expect_error(as_envelope(make_footprint(), depth_min = 200, depth_max = 0),
               "at least")
  expect_error(as_envelope(make_footprint(), depth_min = -10, depth_max = 200),
               "positive metres")
})

test_that("volume() and voxel_to_envelope() reject invalid 3D input", {
  e <- as_envelope(make_footprint(), depth_min = 0, depth_max = 200)
  inv <- e[[c("depth_max", "depth_min")]]
  names(inv) <- c("depth_min", "depth_max")

  # Otherwise volume() reports a negative volume for an impossible domain.
  expect_error(volume(inv), "not a valid one")

  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
  terra::values(r) <- 1
  names(r) <- paste0("t_depth=", c(0, 100, 200))
  v <- as_voxel(r)

  # A reversed voxel would otherwise collapse to an inverted envelope.
  expect_error(voxel_to_envelope(v[[c(3, 2, 1)]]), "not a valid one")
})
