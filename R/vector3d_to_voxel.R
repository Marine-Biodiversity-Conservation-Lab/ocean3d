# SpatVector3D polygons to a presence SpatVoxel on a voxel's levels (#115).
# Each level is filled on its own, so a depth gap between two polygons over
# one cell stays a gap; vector3d_to_envelope() merges them into one range.

#' Convert 3D polygons to a voxel on another voxel's depth levels
#'
#' Puts [SpatVector3D-class] polygons onto the grid and depth levels of
#' `template` as a presence [SpatVoxel-class]: `1` at every cell and level a
#' polygon covers, `NA` elsewhere. A polygon covers a cell when it touches it,
#' as in [vector3d_to_envelope()], and covers a level when its depth range
#' shares some thickness with the slab of water that level stands for, as in
#' [envelope_to_voxel()]. The slabs come from the template's depth edges (see
#' [depth_bounds()]), or from `bounds` when it has none.
#'
#' Each level is filled on its own, so where polygons overlap over a cell,
#' the levels between them stay empty: polygons from 0 to 50 m and from 200
#' to 300 m over the same cell cover the levels in those two ranges and not
#' the ones between. [vector3d_to_envelope()] gives that cell one range, 0 to
#' 300 m, because an envelope holds one interval per cell. Use this function
#' when the next step works on a voxel's levels, for example to compare the
#' polygons with a [SpatVoxel-class] of data; [intersect_3d()],
#' [intersects_3d()] and [mask()] use it for polygons and a voxel.
#'
#' Polygons with `depth_min` equal to `depth_max` have no thickness and cover
#' nothing. A voxel has no time, so the polygons' time intervals are dropped,
#' with a message when any polygon has one.
#'
#' @param x A [SpatVector3D-class] of polygons, from [as_vector3d()] or
#'   [buffer_3d()].
#' @param template A [SpatVoxel-class] whose grid, depth levels and depth
#'   edges the result uses; its values are ignored. Polygons are projected to
#'   its CRS if it isn't WGS84. Longitudes from 0 to 360 work too.
#' @param bounds How to place the edges of the slab each depth level stands
#'   for when `template` carries no depth edges: `"midpoint"` (default) or
#'   `"top"`, as in [envelope_to_voxel()].
#'
#' @return A presence [SpatVoxel-class] on the grid, depths and depth edges of
#'   `template`, with layers named `presence_depth={value}`. Stops with an
#'   error if no cell and level would be covered: when the polygons miss the
#'   grid, have no thickness, or lie outside the template's levels.
#' @seealso [vector3d_to_envelope()] for one depth range per cell,
#'   [envelope_to_voxel()]
#' @examples
#' grid <- terra::rast(nrows = 4, ncols = 4, xmin = 150, xmax = 154,
#'                     ymin = -24, ymax = -20, crs = "EPSG:4326", nlyrs = 4,
#'                     vals = 1)
#' temp <- as_voxel(grid, depths = c(0, 100, 200, 300), varname = "temp")
#' square <- sf::st_polygon(list(rbind(c(151, -23), c(153, -23), c(153, -21),
#'                                     c(151, -21), c(151, -23))))
#' bands <- sf::st_sf(
#'   top = c(0, 200), bottom = c(50, 300),
#'   geometry = sf::st_sfc(square, square, crs = 4326)
#' )
#' polygons <- as_vector3d(bands, depth_min = "top", depth_max = "bottom")
#'
#' # The 100 m level (50 to 150 m) lies between the two bands, so it stays
#' # empty; an envelope would fill it.
#' terra::values(vector3d_to_voxel(polygons, temp))[6, ]
#' terra::values(envelope_to_voxel(vector3d_to_envelope(polygons, temp),
#'                                 depths(temp)))[6, ]
#' @export
vector3d_to_voxel <- function(x, template, bounds = c("midpoint", "top")) {
  bounds <- match.arg(bounds)
  if (!methods::is(x, "SpatVector3D")) {
    .abort_voxel_input("`x` must be SpatVector3D polygons.")
  }
  .check_polygons(x)
  .check_3d(template, "template", "SpatVoxel")
  f <- x@features
  if (nrow(f) == 0L) {
    .abort_voxel_input("`x` has no polygons to put on the voxel's levels.")
  }
  .inform_time_dropped(x)

  depths <- depths(template)
  edges <- depth_bounds(template)
  slab <- .depth_slabs(depths, bounds, edges)
  # Slabs worked out from `depths` are half-open with the deepest closed,
  # given edges all half-open: the rule of .promote_to_voxel(), so a polygon
  # and the same range as an envelope land on the same levels.
  closed <- is.null(edges)
  deepest <- length(depths)
  # A polygon with no thickness covers nothing, as in vector3d_to_envelope().
  thick <- f$depth_max > f$depth_min

  # Only the grid is used. Copies shifted by 360 degrees put outlines on a
  # grid in either longitude convention, as in vector3d_to_envelope().
  grid <- methods::as(template[[1]], "SpatRaster")  # remove the 3D class
  polygons <- terra::vect(sf::st_sf(geometry = f$geometry))
  if (!terra::same.crs(polygons, grid)) {
    polygons <- terra::project(polygons, terra::crs(grid))
  }
  copies <- 1L
  if (terra::is.lonlat(grid)) {
    polygons <- rbind(polygons, terra::shift(polygons, dx = 360),
                      terra::shift(polygons, dx = -360))
    copies <- 3L
  }
  occ <- terra::rast(lapply(seq_len(deepest), function(i) {
    covers <- thick & f$depth_max > slab$lower[i] &
      if (closed && i == deepest) f$depth_min <= slab$upper[i] else
        f$depth_min < slab$upper[i]
    if (!any(covers)) return(terra::setValues(terra::rast(grid), 0))
    .footprint(polygons[which(rep(covers, copies))], grid)
  }))
  if (!any(terra::global(occ, "max", na.rm = TRUE)[, 1] > 0)) {
    .abort_envelope_empty(paste0(
      "The polygons cover no cell and depth level of `template`: they miss ",
      "its grid, have no thickness, or lie outside its depth levels."
    ))
  }
  .presence_voxel(occ, depths, edges)
}
