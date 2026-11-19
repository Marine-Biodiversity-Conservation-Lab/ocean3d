# R/volume.R
#
# Volume for the package's two 3D representations. The generic is declared in
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
# NA outside the domain): the raster itself, or its total.
.volume_result <- function(vol_rast, per_cell) {
  if (!rlang::is_bool(per_cell)) {
    rlang::abort("`per_cell` must be TRUE or FALSE.",
                 class = "ocean3d_error_per_cell", call = NULL)
  }
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
    vol_rast <- terra::mask(.cell_area_km2(x) * (occupied_m / 1000),
                            .voxel_present(occ))

    .volume_result(vol_rast, per_cell)
  }
)

# Rejected here rather than by S4's default "unable to find an inherited
# method", which says nothing about how to build the right object.
#' @rdname volume
#' @export
setMethod("volume", "SpatRaster", function(x, ...) .stop_not_volume())
