# R/utils-3d.R
#
# Internal helpers shared by the operations on the package's two 3D
# representations: R/volume.R (volume and volumetric overlap) and
# R/intersect.R (spatial query). Nothing here is exported.
#
# SpatEnvelope and SpatVoxel both describe a 3D domain over a 2D grid, but they
# store depth in dual roles -- the envelope as a per-cell value, the voxel as
# the layer index. The helpers below are the pieces both files need to read
# that domain: occupancy, presence, depth bounds, cell areas, and the checks
# that two inputs can be combined at all.


# Internal: cell areas in km^2. terra::cellSize() returns m^2 by default, and
# every volume here is km^3, so the unit change happens once, in one place.
.cell_area_km2 <- function(x) {
  terra::cellSize(x[[1]], unit = "km")
}

# Internal: two rasters can only be combined cell by cell if they are on the
# same grid. Checks `x` against the grid of `template`, naming them with the
# labels `arg` and `grid` (with backticks, e.g. "`mask`"). A CRS mismatch is
# told apart from an extent or resolution mismatch, since the fix differs
# (project, or resample). The one producer of `ocean3d_error_grid_mismatch`,
# for the query verbs, the constructors (through .grid_value()) and the 2D
# side of a query (.footprint()).
.check_same_grid <- function(x, template, arg = "`y`", grid = "`x`") {
  problem <- if (!terra::same.crs(x, template)) {
    "has a different CRS than {grid}. Project it onto that grid first."
  } else if (!isTRUE(terra::compareGeom(x, template, crs = FALSE,
                                        stopOnError = FALSE,
                                        messages = FALSE))) {
    paste("does not align with {grid}; its extent or resolution differs.",
          "Resample it onto that grid first.")
  }
  if (!is.null(problem)) {
    rlang::abort(glue::glue("{arg} ", problem),
                 class = "ocean3d_error_grid_mismatch", call = NULL)
  }
  invisible(TRUE)
}

# Internal: reject arguments a method cannot act on, rather than letting the
# generic's `...` swallow them. Passing `bounds` to a pair of envelopes is a
# misunderstanding worth naming, not a no-op.
.no_extra_args <- function(...) {
  extra <- list(...)
  if (length(extra) == 0L) return(invisible(TRUE))
  nms <- names(extra)
  nms <- if (is.null(nms)) rep("", length(extra)) else nms
  nms[nms == ""] <- "<unnamed>"
  rlang::abort(glue::glue(
    "unused argument(s): {toString(nms)}. `bounds` applies only when one ",
    "input is a SpatVoxel; envelopes and polygons carry their own depth ",
    "intervals."
  ), class = "ocean3d_error_unused_argument", call = NULL)
}

# Internal: a voxel's occupancy as a 0/1 stack, one layer per depth, ready to
# multiply by slab thicknesses. occupied() and voxel_to_envelope() pass the
# caller's `fun`; the query verbs read presence as it stands.
.voxel_occupancy <- function(x, fun = function(v) !is.na(v)) {
  if (!is.function(fun)) {
    stop("`fun` must take a vector of cell values and return a logical ",
         "vector of the same length.", call. = FALSE)
  }
  n_cells <- terra::ncell(x)

  occ <- terra::rast(lapply(seq_len(terra::nlyr(x)), function(i) {
    # `[[` keeps the SpatVoxel tag, and a tagged raster dispatches to the
    # depth-aware methods in R/intersect.R, so strip it first. `fun` reads
    # values rather than the SpatRaster, so ordinary R predicates work: on a
    # one-layer SpatRaster, `v > mean(v)` compares the layer against itself.
    layer <- methods::as(x[[i]], "SpatRaster")  # remove the 3D class
    hit <- as.logical(fun(terra::values(layer)[, 1]))
    if (length(hit) != n_cells) {
      stop("`fun` must return one TRUE/FALSE per cell value; got ",
           length(hit), " values for ", n_cells, " cells.", call. = FALSE)
    }
    # A rejected cell and an NA cell are both unoccupied. NA becomes 0 so
    # that it does not spoil the layer sum.
    hit[is.na(hit)] <- FALSE
    terra::setValues(terra::rast(layer), as.numeric(hit))
  }))
  names(occ) <- names(x)
  # terra tags the stack with the voxel's class; occupancy is not a voxel.
  methods::as(occ, "SpatRaster")  # remove the 3D class
}

# Internal: the shared depth levels of two voxels, and their shared depth
# edges: whichever voxel carries them, or NULL when neither does. Resampling
# one onto the other would be a guess about the vertical axis the caller is
# better placed to make, so different levels, or different edges for the same
# levels, are an error.
.shared_depths <- function(x, y) {
  dx <- depths(x)
  dy <- depths(y)
  if (!isTRUE(all.equal(dx, dy))) {
    stop("`x` and `y` must be sampled at the same depth levels; got ",
         toString(dx), " and ", toString(dy),
         ". Rebuild one on the other's depths.", call. = FALSE)
  }
  bx <- depth_bounds(x)
  by <- depth_bounds(y)
  if (!is.null(bx) && !is.null(by) &&
        !isTRUE(all.equal(bx[c("lower", "upper")], by[c("lower", "upper")]))) {
    rlang::abort(glue::glue(
      "`x` and `y` carry different depth edges for the same depth levels. ",
      "Rebuild one with the other's edges (`depth_bounds` in as_voxel())."
    ), class = "ocean3d_error_depth_edges_differ", call = NULL)
  }
  list(depths = dx, depth_bounds = if (is.null(bx)) by else bx)
}

# Internal: a 3D input must be one of the classes `class` allows, and satisfy
# that class's validity rules, before it is read. The class boundary is the
# contract: a raster that merely carries the right layer names isn't one of
# the package's 3D representations. terra also keeps the SpatVoxel /
# SpatEnvelope tag through operations that change the layer set, so an object
# can arrive tagged but not valid; the checks here say how to build or
# rebuild it, instead of failing deeper in on a missing layer name.
.check_3d <- function(x, arg, class = c("SpatVoxel", "SpatEnvelope")) {
  builders <- c(
    SpatVoxel = "as_voxel() for a multi-depth voxel",
    SpatEnvelope = "as_envelope() or vect_to_envelope() for a min-max envelope"
  )
  is_class <- vapply(class, function(cl) methods::is(x, cl), logical(1))
  if (!any(is_class)) {
    wanted <- glue::glue_collapse(class, sep = ", ", last = " or ")
    got <- paste(class(x), collapse = "/")
    hint <- glue::glue_collapse(builders[class], sep = ", or with ")
    rlang::abort(glue::glue(
      "`{arg}` must be a {wanted}; got {got}. Build one with {hint}."
    ), class = "ocean3d_error_not_3d", call = NULL)
  }
  ok <- methods::validObject(x, test = TRUE)
  if (!isTRUE(ok)) {
    tagged <- class[is_class][[1]]
    rebuild <- c(SpatVoxel = "as_voxel()", SpatEnvelope = "as_envelope()")
    rlang::abort(glue::glue(
      "`{arg}` is tagged {tagged} but is not a valid one: {ok}. A terra ",
      "operation probably changed its layers; rebuild it with ",
      "{rebuild[[tagged]]}."
    ), class = "ocean3d_error_invalid_3d", call = NULL)
  }
  invisible(TRUE)
}

# Internal: one depth layer of an envelope as a plain SpatRaster. `[[` keeps
# the SpatEnvelope tag on a single layer, which is both invalid and, since a
# tagged raster dispatches to the depth-aware methods in R/intersect.R,
# hazardous to compute with. Internals that compute with an envelope's
# layers read them through here; the validity method's in-memory check reads
# the cell values directly.
.envelope_layer <- function(x, which = c("depth_min", "depth_max")) {
  which <- match.arg(which)
  methods::as(x[[which]], "SpatRaster")  # remove the 3D class
}

# Internal: where `x` is present, as a plain single-layer logical SpatRaster
# named "footprint": TRUE where present, FALSE elsewhere, never NA. Callers
# that mask with it wrap it in terra::ifel(fp, 1, NA).
#
# An envelope cell is present when both depth layers are non-NA, a voxel cell
# when it is occupied at any level. The 2D side of a query is put on
# `template`'s grid: a plain SpatRaster is read by its non-NA pattern and must
# be single-layer and already on the grid. Polygons (SpatVector, sf, sfc) are
# projected onto the template's CRS when they differ, then rasterised with the
# rule vect_to_envelope() uses: terra::mask()'s default, under which every
# cell a polygon touches is inside. Anything else is rejected, naming `arg`,
# with the list of accepted types.
.footprint <- function(x, template = NULL, arg = "y") {
  if (inherits(x, c("sf", "sfc"))) x <- terra::vect(x)

  out <- if (methods::is(x, "SpatEnvelope")) {
    !is.na(.envelope_layer(x, "depth_min")) &
      !is.na(.envelope_layer(x, "depth_max"))
  } else if (methods::is(x, "SpatVoxel")) {
    # The same rule as .voxel_occupancy(): a cell is occupied at a level where
    # it is not NA. allNA() reads that in one pass, without pulling each
    # layer into R.
    !terra::allNA(methods::as(x, "SpatRaster"))  # remove the 3D class
  } else if (inherits(x, "SpatVector")) {
    if (!terra::same.crs(x, template)) {
      x <- terra::project(x, terra::crs(template))
    }
    # remove the 3D class
    ones <- terra::setValues(
      terra::rast(methods::as(template[[1]], "SpatRaster")), 1
    )
    !is.na(terra::mask(ones, x))
  } else if (inherits(x, "SpatRaster")) {
    if (terra::nlyr(x) != 1) {
      rlang::abort(glue::glue(
        "`{arg}` must be a single-layer footprint when it is a plain ",
        "SpatRaster; got {terra::nlyr(x)} layers. If it is a multi-depth ",
        "raster, wrap it with as_voxel() so its depth axis is used."
      ), class = "ocean3d_error_footprint_layers", call = NULL)
    }
    # The 3D side is `x` unless the inputs were swapped to put it first.
    grid <- if (arg == "x") "`y`" else "`x`"
    .check_same_grid(x, template, glue::glue("`{arg}`"), grid)
    !is.na(x)
  } else {
    got <- paste(class(x), collapse = "/")
    rlang::abort(glue::glue(
      "`{arg}` must be a SpatEnvelope, a SpatVoxel, a single-layer ",
      "SpatRaster footprint, or polygons (SpatVector, sf, sfc); got {got}."
    ), class = "ocean3d_error_footprint_type", call = NULL)
  }
  names(out) <- "footprint"
  methods::as(out, "SpatRaster")  # remove the 3D class
}
