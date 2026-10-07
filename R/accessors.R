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
  .parse_depth_layers(x)
}

# Internal: parse numeric depths from the `{variable}_depth={value}` layer
# naming convention used throughout the package. Returns a numeric vector the
# same length as its input, NA for any name that does not match. With
# `error = TRUE` (the default) it stops when no name matches; validity methods
# pass `error = FALSE` so they can report the failure themselves.
#
# The `_names()` form takes the layer names directly, for callers that have
# names but no raster.
.parse_depth_layers <- function(rast, error = TRUE) {
  .parse_depth_layers_names(layer_names = names(rast), error = error)
}

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
