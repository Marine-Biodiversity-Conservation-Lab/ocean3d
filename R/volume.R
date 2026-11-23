# R/volume.R
#
# Volume for the package's 3D representations. The generic is declared in
# R/AllGenerics.R; this file holds the arithmetic.
#
# SpatEnvelope and SpatVoxel both describe a 3D domain over a 2D grid, but they
# store depth in dual roles -- the envelope as a per-cell value, the voxel as
# the layer index -- so the vertical extent is derived differently in each and
# each gets its own method rather than a shared SpatVolume one. See
# [SpatVolume-class].
#
# The internals both this file and R/intersect.R rely on (occupancy, presence,
# cell areas, grid and depth checks) live in R/utils-3d.R. Where two domains
# overlap is intersect_3d()'s job (R/intersect.R); volume() measures the result
# like any other domain.

# Internal: what both methods return from their per-cell volume raster (km³,
# NA outside the domain): the raster itself, or its total. The generic has
# already checked `per_cell`.
.volume_result <- function(vol_rast, per_cell) {
  if (!per_cell) {
    # terra sums an all-NA raster to NA; an empty domain holds no water.
    total <- terra::global(vol_rast, "sum", na.rm = TRUE)[[1]]
    return(if (is.na(total)) 0 else total)
  }
  names(vol_rast) <- "volume"
  methods::as(vol_rast, "SpatRaster")  # remove the 3D class
}

#' @rdname volume
#' @export
setMethod("volume", "SpatEnvelope", function(x, per_cell = FALSE, ...) {
  .no_extra_args(...)
  .check_3d(x, "x")

  depth_extent <- .envelope_layer(x, "depth_max") -
    .envelope_layer(x, "depth_min")

  # Volume per cell in km³ (depth in m, converted to km). NA depths, outside
  # the domain, give NA.
  .volume_result(.cell_area_km2(x) * (depth_extent / 1000), per_cell)
})

#' @rdname volume
#' @export
setMethod(
  "volume", "SpatVoxel",
  function(x, per_cell = FALSE, bounds = c("midpoint", "top")) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")

    occ <- .voxel_occupancy(x)
    slab <- .depth_slabs(depths(x), bounds, depth_bounds(x))
    # Summing thickness over occupied levels only is what makes an interior
    # gap cost volume, unlike voxel_to_envelope() plus the envelope formula.
    thickness_m <- slab$upper - slab$lower
    # The raster stays on the left of the operator: terra mishandles
    # `numeric <op> SpatRaster`. The vector is recycled layer-wise.
    occupied_m <- sum(occ * thickness_m)
    # Occupancy is 0 at unoccupied levels, so a cell with none would read 0;
    # masking makes it NA, as outside an envelope.
    vol_rast <- terra::ifel(.footprint(x),
                            .cell_area_km2(x) * (occupied_m / 1000), NA)

    .volume_result(vol_rast, per_cell)
  }
)

# SpatVector3D polygons are measured from their outlines, not a grid
# (SPEC.md, "Summarise: volume calculations"): geodesic area times the depth
# range, with terra::expanse(). It agrees with terra::cellSize() for a polygon
# drawn along a cell's edges to within about 0.02% at low latitudes: the
# polygon's edges are geodesics, while cells are bounded by parallels.
# Points and lines have no volume.
#' @rdname volume
#' @export
setMethod("volume", "SpatVector3D", function(x, per_cell = FALSE, ...) {
  .no_extra_args(...)
  if (x@type != "polygons") {
    rlang::abort(glue::glue(
      "These are {x@type}, which have no volume. Volume is defined for ",
      "shapes with a depth range: voxels, envelopes and polygons. ",
      "buffer_3d() turns {x@type} into polygons."
    ), class = "ocean3d_error_volume_type", call = NULL)
  }
  if (per_cell) {
    rlang::abort(
      "Polygons have no cells, so `per_cell = TRUE` doesn't apply.",
      class = "ocean3d_error_volume_per_cell", call = NULL
    )
  }
  f <- x@features
  if (nrow(f) == 0L) return(0)
  area_km2 <- terra::expanse(terra::vect(f$geometry), unit = "km")
  sum(area_km2 * (f$depth_max - f$depth_min) / 1000)
})

# Rejected here rather than by S4's default "unable to find an inherited
# method", which says nothing about how to build the right object.
#' @rdname volume
#' @export
setMethod("volume", "SpatRaster",
          function(x, per_cell = FALSE, ...) .check_3d(x, "x"))

# A SpatProxy is measured one time step at a time, reading each step through
# .proxy_to_voxel(); a step larger than the read budget is read in groups of
# depth layers. A voxel read from a proxy always carries its depth edges, so
# each layer stands for the same slab whichever layers are read with it, and
# the groups' volumes add up to the whole step's.
#' @rdname volume
#' @param fun SpatProxy only. Which values count as occupied, as in
#'   [occupied()]: `NULL` (default) counts every non-`NA` value; otherwise a
#'   function of a layer's values returning `TRUE` where the cell is
#'   occupied, for example `function(v) v > 10 & v < 15` for the water
#'   between 10 and 15 degrees.
#' @section Memory:
#' `volume()` on a [SpatProxy-class] holds one time step in memory at a time.
#' When one step is larger than the read budget, the option
#' `ocean3d.max_read_mb` (default 1024, in megabytes of doubles: cells times
#' layers times 8 bytes), the step is read in groups of depth layers that
#' each fit the budget, and the groups' volumes are added. Set the option
#' lower on a machine with little memory, or higher for faster reads.
#' @section Overlap through time:
#' `intersect_3d()` and the other query functions take voxels, not a
#' `SpatProxy`. To measure the overlap of a time-varying dataset with a
#' shape, read each step and intersect it:
#' `volume(intersect_3d(as_voxel(x, time = i), shape))` for each step `i`.
#' @export
setMethod("volume", "SpatProxy", function(x, per_cell = FALSE, fun = NULL,
                                          ...) {
  .no_extra_args(...)
  variable <- x@variable
  axis <- x@axes$time
  n_steps <- if (is.null(axis)) 1L else variable$dim_lengths[[axis$dim]]
  dated <- !is.null(axis) && !isTRUE(axis$climatology)
  times <- if (dated) axis$cftime$as_timestamp(asPOSIX = TRUE) else
    .POSIXct(rep(NA_real_, n_steps), tz = "UTC")

  # Groups of depth layers (first and last index, storage order) that fit
  # the read budget; a proxy without a depth axis gets one group of nothing,
  # and .proxy_to_voxel() stops for it.
  depth_dim <- variable$axis_dims$depth
  n_depth <- if (is.null(depth_dim)) 1L else variable$dim_lengths[[depth_dim]]
  n_cells <- prod(variable$dim_lengths[unlist(variable$axis_dims[c("lon",
                                                                   "lat")])])
  budget <- getOption("ocean3d.max_read_mb", 1024)
  per_group <- max(1L, min(n_depth, floor(budget * 1024^2 / (8 * n_cells))))
  groups <- if (is.null(depth_dim)) {
    list(NULL)
  } else {
    lapply(seq(1L, n_depth, by = per_group), function(first) {
      c(first, min(first + per_group - 1L, n_depth))
    })
  }

  layers <- lapply(seq_len(n_steps), function(step) {
    parts <- lapply(groups, function(span) {
      v <- .proxy_to_voxel(x, if (is.null(axis)) NULL else step,
                           variable$name, depth = span)
      if (!is.null(fun)) v <- occupied(v, fun)
      volume(v, per_cell = TRUE)
    })
    parts <- terra::rast(parts)
    # NA where no group has an occupied level, as volume() gives.
    terra::ifel(terra::allNA(parts), NA, sum(parts, na.rm = TRUE))
  })

  if (!per_cell) {
    totals <- vapply(layers, .volume_result, numeric(1), per_cell = FALSE)
    return(data.frame(time = times, volume_km3 = totals))
  }
  out <- terra::rast(layers)
  names(out) <- rep("volume", n_steps)
  if (dated) terra::time(out) <- times
  out
})
