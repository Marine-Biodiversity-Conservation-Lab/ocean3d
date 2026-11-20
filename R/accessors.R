# Accessors for the package's 3D classes.

#' Depths of a voxel's layers
#'
#' Read the depth axis of a [SpatVoxel-class]. Depth is the layer index in a
#' voxel, and the depth itself is carried in the layer name following the
#' `{variable}_depth={value}` convention used throughout the package; this
#' parses it back out.
#'
#' `x` is normally a `SpatVoxel`, but any [terra::SpatRaster] whose layer names
#' follow the convention works, as does a bare character vector of layer names —
#' useful for the row names of a [terra::global()] result, which carry the layer
#' names but not the raster.
#'
#' Depths are positive metres increasing downward.
#'
#' @param x A [SpatVoxel-class], a [terra::SpatRaster], or a character vector of
#'   layer names.
#'
#' @returns Numeric vector, one depth per layer of `x` (or per element, when `x`
#'   is character). `NA` for any name that does not follow the convention; an
#'   error when no name does.
#'
#' @seealso [as_voxel()], which builds those layer names.
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
  .parse_depth_layers_names(names(x))
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
# `{variable}_depth={value}` convention used throughout the package. Returns a
# numeric vector the same length as its input, NA for any name that does not
# match. With `error = TRUE` (the default) it stops when no name matches.
# Package code reads a raster's depths with depths(); only the SpatVoxel
# validity method and as_voxel() call this directly, with `error = FALSE`, so
# they can report bad names themselves.
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
