#' Extract a 3D raster to an area and a depth band
#'
#' Keeps the values of a [SpatVoxel-class] inside an area polygon and a depth
#' band, and crops the voxel to the area. It is a shortcut for the build and
#' intersect steps: the area becomes [SpatVector3D-class] polygons from
#' `min_depth` to `max_depth` with [as_vector3d()], the voxel is cropped to
#' them with [terra::crop()], and [mask()] keeps the values inside. To
#' summarise the result, pass it to [summarise_3d()].
#'
#' A depth level is kept where the band shares some thickness with the slab
#' of water the level stands for: the voxel's own depth edges when it has
#' them (see [depth_bounds()]), otherwise the slabs `bounds` gives. A band
#' that only touches a slab's edge doesn't keep that level. Every layer of
#' the voxel is returned; levels outside the band, and cells outside the
#' area, are `NA`. A cell is inside the area when the polygon touches it.
#'
#' `rast_3d` can also be a [SpatProxy-class]. Only the cells under the area
#' are then read, for the time step chosen with `time` (see [as_voxel()]).
#'
#' To keep a voxel inside a depth window that changes from cell to cell, such
#' as a species range bounded by the seafloor, build a [SpatEnvelope-class]
#' and mask with it instead: `mask(rast_3d, range_env)`. See [mask-3d].
#'
#' @param area sf, sfc, or SpatVector polygons. Projected to WGS84 by
#'   [as_vector3d()], then placed on the voxel's grid.
#' @param rast_3d A [SpatVoxel-class], e.g. from [as_voxel()] or
#'   [envelope_to_voxel()], or a [SpatProxy-class] from [as_proxy()].
#' @param min_depth,max_depth The top and bottom of the depth band, in metres
#'   positive down. `NULL` (the default) runs to that end of the voxel: the
#'   surface for `min_depth`, and the bottom of the deepest slab for
#'   `max_depth`.
#' @param bounds How a voxel without depth edges places each level's slab:
#'   `"midpoint"` (default) or `"top"`, as in [mask()]. Not used when the
#'   voxel has depth edges.
#' @param time For a [SpatProxy-class] only: the time step to read, as in
#'   [as_voxel()]. A time window gives a list of voxels, one per step.
#'
#' @returns A [SpatVoxel-class] cropped to `area`, with the voxel's layers,
#'   `NA` outside the area and the depth band. A list of them for a
#'   [SpatProxy-class] read over a time window.
#'
#' @seealso [mask-3d], [summarise_3d()], [volume()];
#'   `vignette("workflow-values-in-shapes")` for the steps written out.
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
#' # 40-600 m shares water with the slabs of the 50, 100 and 500 m levels.
#' out <- extract_to_area(area, v, min_depth = 40, max_depth = 600)
#' terra::global(out, "notNA")
#' terra::ext(out)   # cropped to the polygon
#'
#' # The mean temperature in the area and band, weighted by volume.
#' summarise_3d(out)
#' @export
extract_to_area <- function(area, rast_3d, min_depth = NULL, max_depth = NULL,
                            bounds = c("midpoint", "top"), time = NULL) {
  bounds <- match.arg(bounds)
  if (!inherits(area, c("sf", "sfc", "SpatVector"))) {
    rlang::abort(glue::glue(
      "`area` must be sf, sfc, or SpatVector polygons; got ",
      "{paste(class(area), collapse = '/')}."
    ), class = "ocean3d_error_area_input", call = NULL)
  }
  if (!methods::is(rast_3d, "SpatVoxel") &&
        !methods::is(rast_3d, "SpatProxy")) {
    rlang::abort(glue::glue(
      "`rast_3d` must be a SpatVoxel or a SpatProxy; got ",
      "{paste(class(rast_3d), collapse = '/')}. Wrap a multi-depth SpatRaster ",
      "with as_voxel()."
    ), class = "ocean3d_error_area_input", call = NULL)
  }
  # Only the cells under the area are kept, or read from a proxy, which
  # projects the area itself; a voxel in another CRS needs the area projected.
  target <- if (inherits(area, "SpatVector")) area else terra::vect(area)
  if (methods::is(rast_3d, "SpatVoxel") && !terra::same.crs(target, rast_3d)) {
    target <- terra::project(target, terra::crs(rast_3d))
  }
  voxels <- as_voxel(crop(rast_3d, target), time = time)

  # By default the band runs from the surface to the deepest water the voxel
  # stands for; every voxel of a time window has the same levels.
  first <- if (is.list(voxels)) voxels[[1L]] else voxels
  if (is.null(min_depth)) min_depth <- 0
  if (is.null(max_depth)) {
    max_depth <- max(depths(first), depth_bounds(first)$upper)
  }
  if (inherits(area, "sfc")) area <- sf::st_sf(geometry = area)
  shape <- as_vector3d(area, depth_min = min_depth, depth_max = max_depth)
  if (is.list(voxels)) return(lapply(voxels, mask, shape, bounds = bounds))
  mask(voxels, shape, bounds = bounds)
}
