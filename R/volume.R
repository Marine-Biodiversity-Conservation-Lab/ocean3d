# R/volume.R
#
# Volume and volumetric overlap for the package's two 3D representations. The
# generics are declared in R/AllGenerics.R; this file holds the arithmetic.
#
# SpatEnvelope and SpatVoxel both describe a 3D domain over a 2D grid, but they
# store depth in dual roles -- the envelope as a per-cell value, the voxel as
# the layer index -- so the vertical extent is derived differently in each and
# each gets its own method rather than a shared SpatVolume one. See
# [SpatVolume-class].
#
# The internals both this file and R/intersect.R rely on (occupancy, presence,
# cell areas, grid and depth checks) live in R/utils-3d.R. Where two domains
# overlap is intersect_3d()'s job (R/intersect.R); this file only measures the
# result.

# ---- volume ------------------------------------------------------------

#' @rdname volume
#' @export
setMethod("volume", "SpatEnvelope", function(x, ...) {
  .no_extra_args(...)
  .check_3d(x, "x")

  depth_extent <- .envelope_layer(x, "depth_max") -
    .envelope_layer(x, "depth_min")

  # Volume per cell in km³ (depth in m, converted to km)
  vol_rast <- .cell_area_km2(x) * (depth_extent / 1000)

  terra::global(vol_rast, "sum", na.rm = TRUE)[[1]]
})

#' @rdname volume
#' @export
setMethod(
  "volume", "SpatVoxel",
  function(x, bounds = c("top", "midpoint")) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")

    depths <- .parse_depth_layers(x)
    occ <- .voxel_occupancy(x)
    vol_rast <- .voxel_cell_volume(occ, depths, bounds, .cell_area_km2(x))

    terra::global(vol_rast, "sum", na.rm = TRUE)[[1]]
  }
)

# Rejected here rather than by S4's default "unable to find an inherited
# method", which says nothing about how to build the right object.
#' @rdname volume
#' @export
setMethod("volume", "SpatRaster", function(x, ...) .stop_not_volume())

# ---- calc_volume_overlap ----------------------------------------------------

# Internal: a depth layer of an envelope, renamed for the overlap stack.
.named_layer <- function(e, which, name) {
  out <- .envelope_layer(e, which)
  names(out) <- name
  out
}

#' @rdname calc_volume_overlap
#' @export
setMethod(
  "calc_volume_overlap", c("SpatEnvelope", "SpatEnvelope"),
  function(x, y, ...) {
    .no_extra_args(...)
    # Validates both inputs and checks they share a grid.
    shared <- intersect_3d(x, y)

    cell_area_km2 <- .cell_area_km2(x)

    dmin_a <- .named_layer(x, "depth_min", "depth_min_a")
    dmax_a <- .named_layer(x, "depth_max", "depth_max_a")
    dmin_b <- .named_layer(y, "depth_min", "depth_min_b")
    dmax_b <- .named_layer(y, "depth_max", "depth_max_b")
    overlap_min <- .named_layer(shared, "depth_min", "depth_min_overlap")
    overlap_max <- .named_layer(shared, "depth_max", "depth_max_overlap")

    # Per-cell volume for A and B
    vol_a <- cell_area_km2 * ((dmax_a - dmin_a) / 1000)
    names(vol_a) <- "volume_a"
    vol_b <- cell_area_km2 * ((dmax_b - dmin_b) / 1000)
    names(vol_b) <- "volume_b"

    # The shared interval is NA wherever the domains do not meet, so the fill
    # below is what distinguishes the two ways a cell can have no overlap: 0
    # where both are present but vertically disjoint, NA where one is absent.
    both_present <- .envelope_present(x) & .envelope_present(y)
    overlap_depth <- overlap_max - overlap_min
    overlap_depth <- terra::ifel(is.na(overlap_depth) & both_present, 0,
                                 overlap_depth)
    vol_overlap <- cell_area_km2 * (overlap_depth / 1000)
    names(vol_overlap) <- "volume_overlap"

    .as_plain_raster(c(dmin_a, dmax_a, dmin_b, dmax_b, overlap_min, overlap_max,
                       vol_a, vol_b, vol_overlap))
  }
)

#' @rdname calc_volume_overlap
#' @export
setMethod(
  "calc_volume_overlap", c("SpatVoxel", "SpatVoxel"),
  function(x, y, bounds = c("top", "midpoint")) {
    bounds <- match.arg(bounds)
    # Validates both inputs and checks they share a grid and depth levels.
    shared <- intersect_3d(x, y)
    depths <- .parse_depth_layers(x)

    occ_a <- .voxel_occupancy(x)
    occ_b <- .voxel_occupancy(y)
    # `shared` is presence (1/NA). Reading it back as occupancy gives the
    # 1/0 stack the volume arithmetic needs.
    occ_overlap <- .voxel_occupancy(shared)

    cell_area_km2 <- .cell_area_km2(x)

    present_a <- .voxel_present(occ_a)
    present_b <- .voxel_present(occ_b)
    both_mask <- terra::mask(present_a, present_b)

    # Summary depth bounds, not a claim that the domain is solid between them.
    bounds_a <- .occupancy_bounds(occ_a)
    bounds_b <- .occupancy_bounds(occ_b)
    bounds_overlap <- .occupancy_bounds(occ_overlap)

    dmin_a <- .named_layer(bounds_a, "depth_min", "depth_min_a")
    dmax_a <- .named_layer(bounds_a, "depth_max", "depth_max_a")
    dmin_b <- .named_layer(bounds_b, "depth_min", "depth_min_b")
    dmax_b <- .named_layer(bounds_b, "depth_max", "depth_max_b")
    overlap_min <- .named_layer(bounds_overlap, "depth_min",
                                "depth_min_overlap")
    overlap_max <- .named_layer(bounds_overlap, "depth_max",
                                "depth_max_overlap")

    # Volumes come from the occupancy stacks, so interior gaps are excluded
    # even though the depth bounds above span them.
    vol_a <- terra::mask(
      .voxel_cell_volume(occ_a, depths, bounds, cell_area_km2), present_a)
    names(vol_a) <- "volume_a"
    vol_b <- terra::mask(
      .voxel_cell_volume(occ_b, depths, bounds, cell_area_km2), present_b)
    names(vol_b) <- "volume_b"

    # Already 0 where both are present but share no occupied level; masking to
    # both_mask makes it NA where either is absent, matching the 2.5D method.
    vol_overlap <- terra::mask(
      .voxel_cell_volume(occ_overlap, depths, bounds, cell_area_km2),
      both_mask)
    names(vol_overlap) <- "volume_overlap"

    .as_plain_raster(c(dmin_a, dmax_a, dmin_b, dmax_b, overlap_min, overlap_max,
                       vol_a, vol_b, vol_overlap))
  }
)

#' @rdname calc_volume_overlap
#' @export
setMethod(
  "calc_volume_overlap", c("SpatEnvelope", "SpatVoxel"),
  function(x, y, bounds = c("top", "midpoint"), ...) {
    bounds <- match.arg(bounds)
    calc_volume_overlap(.promote_to_voxel(x, y, bounds), y,
                        bounds = bounds, ...)
  }
)

#' @rdname calc_volume_overlap
#' @export
setMethod(
  "calc_volume_overlap", c("SpatVoxel", "SpatEnvelope"),
  function(x, y, bounds = c("top", "midpoint"), ...) {
    bounds <- match.arg(bounds)
    calc_volume_overlap(x, .promote_to_voxel(y, x, bounds),
                        bounds = bounds, ...)
  }
)

# One guard catches every pairing with a bare SpatRaster: any signature
# naming the SpatVolume union would tie with this one and lose.
#' @rdname calc_volume_overlap
#' @export
setMethod("calc_volume_overlap", c("SpatRaster", "SpatRaster"),
          function(x, y, ...) .stop_not_volume())
