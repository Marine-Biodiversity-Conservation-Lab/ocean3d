# Accessors for the package's 3D classes.

#' Depths of a voxel's layers
#'
#' Read the depth of each layer of a [SpatVoxel-class]. Depth is the layer
#' index in a voxel, and each layer's depth is stored in terra's own per-layer
#' field, [terra::depth()]. Layer names are labels for people to read: the
#' depths are never read from them, so renaming the layers doesn't change
#' them.
#'
#' `x` can also be a plain [terra::SpatRaster] whose depths are set in
#' [terra::depth()], such as a raster terra read from a netCDF file.
#'
#' Depths are positive metres increasing downward.
#'
#' @param x A [SpatVoxel-class], or a [terra::SpatRaster] with its depths set.
#'
#' @returns Numeric vector, one depth per layer of `x`. An error when the
#'   depths aren't set, for example after a terra operation that drops them
#'   (see [SpatVoxel-class]).
#'
#' @seealso [as_voxel()], which sets the depths.
#'
#' @examples
#' r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
#' terra::values(r) <- runif(terra::ncell(r) * 3)
#' v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")
#'
#' names(v)
#' depths(v)
#'
#' # A per-layer summary has one row per layer, in layer order, so the voxel's
#' # depths label it.
#' per_depth <- terra::global(v, "mean", na.rm = TRUE)
#' per_depth$depth <- depths(v)
#' per_depth
#' @export
depths <- function(x) {
  if (!methods::is(x, "SpatRaster")) {
    rlang::abort(glue::glue(
      "`x` must be a SpatVoxel or a SpatRaster; got ",
      "{paste(class(x), collapse = '/')}."
    ), class = "ocean3d_error_depths_input", call = NULL)
  }
  if (!.depth_set(x)) .abort_no_depths()
  terra::depth(x)
}

# Internal: whether terra holds the raster's depths. terra sets a flag when
# depths are set with terra::depth() or read from a netCDF file, even when
# they are all 0, and clears it whenever an operation resets them to 0
# (app(), `[[<-`, c(), and anything after a layer subset of a raster in
# memory). Unset depths also read as 0, so the flag is the only way to tell a
# 0 m layer from a layer without a depth. A layer subset of a raster in
# memory keeps the values but clears the flag; the SpatVoxel method for
# subset() sets it again. terra has no exported accessor for the flag. If a
# terra release renames it, every raster reads as having no depths, so
# as_voxel() and depths() stop rather than return wrong depths.
.depth_set <- function(x) {
  isTRUE(tryCatch(x@pntr$hasDepth, error = function(e) FALSE))
}

# Internal: error, class `ocean3d_error_no_depths`, for a raster whose depths
# aren't set in terra::depth(). Raised by depths() and as_voxel().
.abort_no_depths <- function() {
  rlang::abort(paste(
    "The raster has no depths in terra::depth(). Pass `depths` to",
    "as_voxel(), or set them with terra::depth(). Layer names aren't read as",
    "depths. terra drops the depths in some operations, such as app(); see",
    "?`SpatVoxel-class`."
  ), class = "ocean3d_error_no_depths", call = NULL)
}

#' Depth edges of a voxel's layers
#'
#' Read the slab of water each layer of a [SpatVoxel-class] stands for, when
#' the voxel carries its own depth edges. Set them with the `depth_bounds`
#' argument of [as_voxel()] or [envelope_to_voxel()], for example from a
#' netCDF file's depth bounds. A voxel read from a [SpatProxy-class] carries
#' the proxy's edges, which are halfway between levels when the file has no
#' depth bounds (see [as_voxel()]).
#'
#' [volume()], [envelope_to_voxel()],
#' [voxel_to_envelope()], [intersect_3d()], [intersects_3d()] and [mask()]
#' use these edges before their `bounds` argument. A voxel without edges, or
#' with edges for only some of its layers (for example after `c()` with a
#' voxel that has none), gives `NULL`, and those functions work the edges out
#' from the depths with `bounds` instead.
#'
#' @param x A [SpatVoxel-class].
#'
#' @returns A data frame with one row per layer of `x`, shallow to deep:
#'   `depth`, and `lower` and `upper`, the top and bottom of the layer's slab,
#'   in metres positive down. `NULL` when `x` doesn't carry edges for every
#'   layer. The data frame can be passed back as `depth_bounds` to
#'   [as_voxel()] or [envelope_to_voxel()].
#'
#' @seealso [depths()] for the layer depths; [as_voxel()] to set the edges.
#'
#' @examples
#' r <- terra::rast(nrows = 1, ncols = 1, nlyrs = 3, vals = 1)
#' v <- as_voxel(r, depths = c(0.5, 1.5, 3),
#'               depth_bounds = cbind(c(0, 1, 2), c(1, 2, 4)))
#' depth_bounds(v)
#'
#' # A layer subset gives the edges of the layers it kept.
#' depth_bounds(v[[2:3]])
#'
#' # A voxel without edges gives NULL, so `bounds` applies.
#' depth_bounds(as_voxel(r, depths = c(0.5, 1.5, 3)))
#' @export
depth_bounds <- function(x) {
  if (!methods::is(x, "SpatVoxel")) {
    rlang::abort(glue::glue(
      "`x` must be a SpatVoxel; got {paste(class(x), collapse = '/')}. ",
      "Build one with as_voxel()."
    ), class = "ocean3d_error_not_voxel", call = NULL)
  }
  # The table is keyed by depth and may hold rows for dropped layers, in any
  # order, or lack rows for some layers: return the rows for the current
  # layers, or NULL if any is missing.
  edges <- x@depth_edges
  row <- .edge_rows(edges$depth, depths(x))
  if (!nrow(edges) || anyNA(row)) return(NULL)
  edges <- edges[row, , drop = FALSE]
  rownames(edges) <- NULL
  edges
}
