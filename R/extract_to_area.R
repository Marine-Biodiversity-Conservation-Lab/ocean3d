#' Extract a 3D raster to an area and a depth band
#'
#' Keeps the values of a [SpatVoxel-class] or [SpatProxy-class] inside an
#' area polygon, and optionally a depth band, and crops the result to the
#' area. It is a shortcut for the build and intersect steps: the data are
#' described with [as_proxy()], cropped to the area with [terra::crop()] and
#' read with [as_voxel()], and [mask()] keeps the values inside the area, or
#' inside [SpatVector3D-class] polygons from `min_depth` to `max_depth`
#' built with [as_vector3d()]. To summarise the result, pass it to
#' [summarise_3d()].
#'
#' A cell is inside the area when the polygon touches it. Without
#' `min_depth` and `max_depth`, every depth level of those cells has values.
#' With a band, a level has values where the band overlaps the depth range of
#' water the level represents by more than its edge. Those depth ranges come
#' from the voxel's depth edges when it has them (see [depth_bounds()]), and
#' otherwise from `bounds`. Every layer of the voxel is returned; levels
#' outside the band, and cells outside the area, are `NA`. A band that
#' overlaps none of the levels gives a warning and a voxel with no values.
#'
#' Only the area's outline is used; its other columns are not returned.
#'
#' To keep a voxel inside a depth window that changes from cell to cell, such
#' as a species range bounded by the seafloor, build a [SpatEnvelope-class]
#' and mask with it instead: `mask(rast_3d, range_env)`. See [mask-3d].
#'
#' @param area sf, sfc, or SpatVector polygons, with a CRS. Projected onto the
#'   data's grid, which may use longitudes from -180 to 180 or 0 to 360.
#' @param rast_3d A [SpatVoxel-class], e.g. from [as_voxel()] or
#'   [envelope_to_voxel()], or a [SpatProxy-class] from [as_proxy()], of
#'   which only the cells under the area are read.
#' @param min_depth,max_depth The top and bottom of the depth band, in metres
#'   positive down. Give either or both; `NULL` for both (the default) keeps
#'   every level. With only one given, the band runs from the surface, or to
#'   the bottom of the deepest level. Under `bounds = "top"`, a voxel without
#'   depth edges gives its deepest level no thickness, so a band never
#'   reaches it.
#' @param bounds Where each level's depth range sits when the voxel has no
#'   depth edges: `"midpoint"` (default) or `"top"`, as in [mask()]. Ignored
#'   when the voxel has depth edges, and without a band.
#' @param time The time step to read, as in [as_voxel()]: needed only when
#'   the data have several time steps. A time window gives a list of voxels,
#'   one per step.
#'
#' @returns A [SpatVoxel-class] cropped to `area`, with the voxel's layers,
#'   `NA` outside the area and the depth band. A list of them for a time
#'   window. A band that reaches levels with no data in the area, such as
#'   water below the seafloor, gives a voxel with no values and no warning.
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
#' # 40-600 m overlaps the depth ranges of the 50, 100 and 500 m levels; the
#' # 0 and 1000 m layers are NA.
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
  got <- function(x) paste(class(x), collapse = "/")
  problem <- if (!inherits(area, c("sf", "sfc", "SpatVector"))) {
    glue::glue(
      "`area` must be sf, sfc, or SpatVector polygons; got {got(area)}."
    )
  } else if (is.na(sf::st_crs(sf::st_as_sf(area)))) {
    "`area` has no CRS, so it can't be placed on the data. Set one first."
  } else if (!methods::is(rast_3d, "SpatVoxel") &&
               !methods::is(rast_3d, "SpatProxy")) {
    glue::glue(
      "`rast_3d` must be a SpatVoxel or a SpatProxy; got {got(rast_3d)}. ",
      "Wrap a multi-depth SpatRaster with as_voxel()."
    )
  }
  if (!is.null(problem)) {
    rlang::abort(problem, class = "ocean3d_error_area_input", call = NULL)
  }
  # Only the outline is used, so the area's columns can't clash with the
  # ones SpatVector3D uses.
  outline <- sf::st_sf(geometry = sf::st_geometry(sf::st_as_sf(area)))
  # A voxel goes through the same path as a file: crop() of a proxy keeps
  # every cell the area's extent overlaps, projects the area and finds it on
  # a 0-360 grid.
  proxy <- rast_3d
  if (methods::is(rast_3d, "SpatVoxel")) proxy <- as_proxy(rast_3d)
  voxels <- as_voxel(crop(proxy, outline), time = time)
  window <- is.list(voxels)
  if (!window) voxels <- list(voxels)

  # The deepest water the levels stand for; every voxel of a time window has
  # the same levels. A shape reaching past it covers every level under
  # either `bounds`, including the deepest level, which "top" gives no
  # thickness.
  deepest <- max(depths(voxels[[1L]]), depth_bounds(voxels[[1L]])$upper)
  whole <- as_vector3d(outline, depth_min = 0, depth_max = deepest + 1)
  # A band runs from the surface, or to the deepest level, when one end is
  # left out. A min_depth below that gives a band of no thickness, which
  # overlaps no level.
  band <- !is.null(min_depth) || !is.null(max_depth)
  shape <- whole
  if (band) {
    if (is.null(min_depth)) min_depth <- 0
    if (is.null(max_depth)) max_depth <- max(deepest, min_depth)
    shape <- as_vector3d(outline, depth_min = min_depth, depth_max = max_depth)
    # Whether the band overlaps a level is a question about depths, not
    # values: a copy filled with 1s, on the same levels and depth edges,
    # keeps values only on the levels the band reaches. A band that reaches
    # a level with no data gives an empty result without a warning.
    first <- voxels[[1L]]
    filled <- terra::setValues(methods::as(first, "SpatRaster"), 1)
    filled <- as_voxel(filled, depths = depths(first),
                       depth_bounds = depth_bounds(first))
    reached <- terra::global(mask(filled, shape, bounds = bounds), "notNA")
    if (all(reached[, 1] == 0)) {
      rlang::warn(glue::glue(
        "The band from {min_depth} to {max_depth} m overlaps none of the ",
        "voxel's depth levels, so no values are kept. Widen the band."
      ), class = "ocean3d_warning_area_band")
    }
  }
  out <- lapply(voxels, mask, shape, bounds = bounds)
  if (window) out else out[[1L]]
}
