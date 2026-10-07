#' Extract a 3D raster to an area and a depth band
#'
#' Crop a [SpatVoxel-class] to an area polygon and select the depth layers
#' within a given depth range. The nearest available depth layers to
#' `min_depth` and `max_depth` are used as the inclusive bounds, so the result
#' always has at least one layer; omit either to run to that end of the voxel.
#'
#' This is the area counterpart of [extract_to_point()]: the same "restrict a
#' 3D source to a target geometry" operation, with a polygon as the target
#' rather than observation points.
#'
#' To restrict a voxel to a species' *per-cell* depth window rather than one
#' depth band across the whole area, build a [SpatEnvelope-class] and mask with
#' it instead: `mask(rast_3d, range_rast)`. See [mask-3d].
#'
#' @param area sf, sfc, or SpatVector. Area polygon to crop the voxel to.
#'   Reprojected to the voxel's CRS when the two differ.
#' @param rast_3d [SpatVoxel-class]. Multi-depth raster, e.g. from [as_voxel()]
#'   or [envelope_to_voxel()].
#' @param min_depth Numeric. Shallowest depth in metres. `NULL` (the default)
#'   starts at the voxel's shallowest layer.
#' @param max_depth Numeric. Deepest depth in metres. `NULL` (the default) runs
#'   to the voxel's deepest layer.
#'
#' @returns A [SpatVoxel-class] cropped to `area` and filtered to the depth
#'   range.
#'
#' @seealso [extract_to_point()] for the point counterpart; [depths()] for the
#'   layer depths the bounds snap to.
#'
#' @examples
#' # A 4x4 voxel over five standard depths.
#' r <- terra::rast(nrows = 4, ncols = 4, xmin = -10, xmax = 10,
#'                  ymin = -10, ymax = 10, nlyrs = 5, crs = "EPSG:4326")
#' terra::values(r) <- seq_len(terra::ncell(r) * 5)
#' v <- as_voxel(r, depths = c(0, 50, 100, 500, 1000), varname = "t_an")
#'
#' area <- terra::vect("POLYGON ((-5 -5, 5 -5, 5 5, -5 5, -5 -5))",
#'                     crs = "EPSG:4326")
#'
#' # 40-600 m snaps to the 50, 100 and 500 m layers.
#' out <- extract_to_area(area, v, min_depth = 40, max_depth = 600)
#' names(out)
#' terra::ext(out)   # cropped to the polygon
#'
#' # Omit the bounds to keep every layer.
#' names(extract_to_area(area, v))
#' @export
extract_to_area <- function(area, rast_3d, min_depth = NULL, max_depth = NULL) {
  # The voxel is checked before `area` is touched, so a malformed depth axis
  # is reported as such rather than as a downstream geometry error.
  if (!methods::is(rast_3d, "SpatVoxel")) {
    stop("Input error for extract_to_area(): `rast_3d` needs to be of ",
         "`SpatVoxel` class. Wrap a multi-depth SpatRaster with as_voxel().",
         call. = FALSE)
  }
  if (!inherits(area, c("sf", "sfc", "SpatVector"))) {
    stop("Input error for extract_to_area(): `area` needs to be of class ",
         "sf, sfc, or SpatVector; got ", paste(class(area), collapse = "/"),
         ".", call. = FALSE)
  }

  d <- depths(rast_3d)

  # NULL means "run to that end of the voxel", so each bound falls back to the
  # depth that already sits there and the nearest-layer snap is unchanged.
  idx_min <- if (is.null(min_depth)) 1L else which.min(abs(d - min_depth))
  idx_max <- if (is.null(max_depth)) {
    length(d)
  } else {
    which.min(abs(d - max_depth))
  }
  idx <- sort(c(idx_min, idx_max))
  selected <- rast_3d[[seq.int(idx[1], idx[2])]]

  if (inherits(area, "sf") || inherits(area, "sfc")) {
    area <- terra::vect(area)
  }
  if (terra::crs(area) != "" && !terra::same.crs(area, selected)) {
    area <- terra::project(area, terra::crs(selected))
  }

  as_voxel(terra::crop(selected, area, mask = TRUE, touches = TRUE))
}
