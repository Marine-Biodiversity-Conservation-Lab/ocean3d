# Accessors for the package's 3D classes.

#' Depths of a voxel's layers
#'
#' Read the depth axis of a [SpatVoxel-class]. Depth is the layer index in a
#' voxel, and each layer's depth is stored in terra's own per-layer field,
#' [terra::depth()]. Layer names are free.
#'
#' `x` is normally a `SpatVoxel`. A plain [terra::SpatRaster] works too: its
#' [terra::depth()] when set, otherwise the depths in layer names that follow
#' the `{variable}_depth={value}` pattern (the names [as_voxel()] gives, and
#' that terra gives layers read from netCDF). A bare character vector of layer
#' names is parsed the same way, which is useful for the row names of a
#' [terra::global()] result.
#'
#' Depths are positive metres increasing downward.
#'
#' @param x A [SpatVoxel-class], a [terra::SpatRaster], or a character vector of
#'   layer names.
#'
#' @returns Numeric vector, one depth per layer of `x` (or per element, when `x`
#'   is character). When depths come from names: `NA` for any name that doesn't
#'   follow the pattern, and an error when no name does.
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
#' # A character vector works too, so the depths of a per-layer summary can be
#' # recovered from its row names.
#' per_depth <- terra::global(v, "mean", na.rm = TRUE)
#' depths(rownames(per_depth))
#'
#' # Names that do not follow the convention come back NA, as long as at least
#' # one name does.
#' depths(c("temp_depth=0", "not_a_depth_layer"))
#' @export
depths <- function(x) {
  if (is.character(x)) {
    return(.parse_depth_layers_names(x))
  }
  if (!methods::is(x, "SpatRaster")) {
    stop("`x` must be a SpatVoxel, a SpatRaster, or a character vector of ",
         "layer names; got ", paste(class(x), collapse = "/"), ".",
         call. = FALSE)
  }
  if (methods::is(x, "SpatVoxel") || .depth_set(x)) return(terra::depth(x))
  .parse_depth_layers_names(names(x))
}

# Internal: whether a plain SpatRaster's depths are set with terra::depth().
# terra flags a raster whose depths were set, but a layer subset keeps the
# values and drops the flag, so any non-zero depth counts as set too. Unset
# depths read as 0. A SpatVoxel always has its depths set by as_voxel(), so
# callers don't ask this of one.
.depth_set <- function(x) {
  flag <- tryCatch(x@pntr$hasDepth, error = function(e) FALSE)
  isTRUE(flag) || any(terra::depth(x) != 0)
}

#' Depth edges of a voxel's layers
#'
#' Read the slab of water each layer of a [SpatVoxel-class] stands for, when
#' the voxel carries its own depth edges. Set them with the `depth_bounds`
#' argument of [as_voxel()] or [envelope_to_voxel()], for example from a
#' netCDF file's depth bounds.
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

# Internal: parse numeric depths from layer names following the
# `{variable}_depth={value}` pattern: the names terra gives layers read from
# netCDF, and that voxels had before depths moved to terra::depth(). Returns a
# numeric vector the same length as its input, NA for any name that does not
# match. With `error = TRUE` (the default) it stops when no name matches.
# Used by depths() for plain rasters without terra::depth(), and by
# as_voxel() with `error = FALSE`, to report bad names itself.
.parse_depth_layers_names <- function(layer_names, error = TRUE) {
  depths <- suppressWarnings(
    as.numeric(stringr::str_extract(layer_names, "(?<=_depth=)-?[0-9.]+"))
  )
  if (error && all(is.na(depths))) {
    stop(
      "No layer names match the '{variable}_depth={value}' convention. ",
      "Got: ", toString(utils::head(layer_names)),
      call. = FALSE
    )
  }
  depths
}
