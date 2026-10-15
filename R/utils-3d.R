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

# Internal: two domains can only be combined cell by cell if they are on the
# same grid. Nothing upstream guarantees this -- the old code took cell areas
# from the first argument and silently reused them for the second.
.check_same_grid <- function(x, y, names = c("x", "y")) {
  if (!terra::compareGeom(x, y, stopOnError = FALSE)) {
    stop("`", names[1], "` and `", names[2], "` must be on the same grid ",
         "(CRS, extent, resolution). Project or resample one onto the other ",
         "first.", call. = FALSE)
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
  stop("unused argument(s): ", toString(nms),
       ". `bounds` applies to SpatVoxel input only; a SpatEnvelope ",
       "carries its depth interval per cell.", call. = FALSE)
}

# Internal: a voxel's occupancy as a 0/1 stack, one layer per depth, ready to
# multiply by slab thicknesses. Only occupied() passes a non-default `fun`;
# the query verbs read presence as it stands.
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
    layer <- .as_plain_raster(x[[i]])
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
  .as_plain_raster(occ)
}

# Internal: vertical extent each depth level stands for, in metres.
# .depth_slabs() (create_3d_objects.R) owns the top/midpoint conventions; this
# is only their height, so the two stay in step by construction.
.slab_thickness <- function(depths, bounds) {
  slab <- .depth_slabs(depths, bounds)
  slab$upper - slab$lower
}

# Internal: per-cell occupied volume of a voxel, km^3. Summing thickness over
# occupied layers is what makes an interior gap cost volume, which is the whole
# reason the voxel path is not just voxel_to_envelope() plus the 2.5D formula.
.voxel_cell_volume <- function(occ, depths, bounds, area_km2) {
  thickness_m <- .slab_thickness(depths, bounds)
  # The raster stays on the left of the operator: terra mishandles
  # `numeric <op> SpatRaster`. The vector is recycled layer-wise.
  occupied_m <- sum(occ * thickness_m)
  area_km2 * (occupied_m / 1000)
}

# Internal: 1 where a cell is occupied at any depth, NA elsewhere -- the voxel
# equivalent of a non-NA envelope cell.
.voxel_present <- function(occ) {
  terra::ifel(sum(occ) > 0, 1, NA)
}

# Internal: shallowest and deepest occupied level of an occupancy stack.
# These are summary bounds only: a voxel with an interior gap is not solid
# between them, which is why volumes are computed from the occupancy itself.
.occupancy_bounds <- function(occ) {
  voxel_to_envelope(as_voxel(occ), fun = function(v) v > 0)
}

# Internal: the shared depth levels of two voxels. Resampling one onto the
# other would be a guess about the vertical axis the caller is better placed to
# make, so a mismatch is an error.
.shared_depths <- function(x, y) {
  dx <- .parse_depth_layers(x)
  dy <- .parse_depth_layers(y)
  if (!isTRUE(all.equal(dx, dy))) {
    stop("`x` and `y` must be sampled at the same depth levels; got ",
         toString(dx), " and ", toString(dy),
         ". Rebuild one on the other's depths.", call. = FALSE)
  }
  dx
}

# Internal: bring an envelope up to a voxel's depth levels. This direction, not
# the reverse: envelope_to_voxel() quantizes, but voxel_to_envelope() fills
# interior gaps, which would overstate an overlap.
.promote_to_voxel <- function(e, v, bounds) {
  envelope_to_voxel(e, depths = .parse_depth_layers(v), bounds = bounds)
}

# Internal: the class boundary is the contract. A raster that merely happens to
# carry the right layer names is not one of the package's 3D representations --
# duck-typing was how these functions used to accept a voxel silently and hand
# back nonsense for it.
.stop_not_volume <- function() {
  stop("input must be a SpatEnvelope or a SpatVoxel, not a bare SpatRaster. ",
       "Build one with as_envelope() or vect_to_envelope() for a 2.5D ",
       "min-max envelope, or with as_voxel() for a multi-depth voxel model.",
       call. = FALSE)
}


# Internal: a 3D input must satisfy its class's validity rules before it is
# read. terra keeps the SpatVoxel / SpatEnvelope tag through operations that
# change the layer set, so an object can arrive tagged but not valid; the
# query verbs check here and say how to rebuild it, instead of failing deeper
# in on a missing layer name.
.check_3d <- function(x, arg) {
  ok <- methods::validObject(x, test = TRUE)
  if (!isTRUE(ok)) {
    stop("`", arg, "` is tagged ", class(x)[[1]], " but is not a valid one: ",
         ok, ". A terra operation probably changed its layers; rebuild it ",
         "with as_voxel() or as_envelope().", call. = FALSE)
  }
  invisible(TRUE)
}

# Internal: one depth layer of an envelope as a plain SpatRaster. `[[` keeps
# the SpatEnvelope tag on a single layer, which is both invalid and, since a
# tagged raster dispatches to the depth-aware methods in R/intersect.R,
# hazardous to compute with. Every internal that reads an envelope's layers
# goes through here.
.envelope_layer <- function(x, which = c("depth_min", "depth_max")) {
  which <- match.arg(which)
  .as_plain_raster(x[[which]])
}

# Internal: where an envelope is present, as a plain logical layer without NA.
.envelope_present <- function(x) {
  !is.na(.envelope_layer(x, "depth_min")) &
    !is.na(.envelope_layer(x, "depth_max"))
}

# Internal: drop the SpatVoxel / SpatEnvelope class tag. terra propagates the
# subclass through nearly every operation (`mask()`, `[[`, arithmetic,
# `ifel()`, `c()`), so anything that is not a 3D domain -- a predicate layer, a
# footprint, an overlap stack -- must be stripped explicitly or it comes back
# wearing a class it does not satisfy. `as()` copies; it does not alias.
.as_plain_raster <- function(x) {
  methods::as(x, "SpatRaster")
}

# Internal: where a 3D object is present at any depth, as a plain single-layer
# SpatRaster: 1 where present, NA elsewhere. An envelope cell is present when
# both depth layers are non-NA; a voxel cell when it is occupied at any level.
.footprint <- function(x) {
  out <- if (methods::is(x, "SpatEnvelope")) {
    terra::ifel(.envelope_present(x), 1, NA)
  } else if (methods::is(x, "SpatVoxel")) {
    .voxel_present(.voxel_occupancy(x))
  } else {
    stop("internal: .footprint() expects a SpatEnvelope or SpatVoxel.",
         call. = FALSE)
  }
  names(out) <- "footprint"
  .as_plain_raster(out)
}

# Internal: the 2D side of a spatial query as a footprint on `template`'s grid,
# 1 where present and NA elsewhere. A SpatRaster is read by its non-NA pattern
# and must be single-layer and already on the grid. Polygons (SpatVector, sf,
# sfc) are projected onto the template's CRS when they differ, then rasterised
# with the same cell-centre rule vect_to_envelope() uses: a cell is inside when
# its centre is. Anything else is rejected with the list of accepted types.
.footprint_of <- function(y, template, arg = "y") {
  if (inherits(y, c("sf", "sfc"))) y <- terra::vect(y)

  if (inherits(y, "SpatVector")) {
    if (!terra::same.crs(y, template)) {
      y <- terra::project(y, terra::crs(template))
    }
    ones <- terra::setValues(terra::rast(.as_plain_raster(template[[1]])), 1)
    out <- terra::mask(ones, y)
  } else if (inherits(y, "SpatRaster")) {
    if (terra::nlyr(y) != 1) {
      stop("`", arg, "` must be a single-layer footprint when it is a plain ",
           "SpatRaster; got ", terra::nlyr(y), " layers. If it is a ",
           "multi-depth raster, wrap it with as_voxel() so its depth axis is ",
           "used.", call. = FALSE)
    }
    .check_same_grid(template, y)
    out <- terra::ifel(is.na(y), NA, 1)
  } else {
    stop("`", arg, "` must be a SpatEnvelope, a SpatVoxel, a single-layer ",
         "SpatRaster footprint, or polygons (SpatVector, sf, sfc); got ",
         paste(class(y), collapse = "/"), ".", call. = FALSE)
  }
  names(out) <- "footprint"
  .as_plain_raster(out)
}
