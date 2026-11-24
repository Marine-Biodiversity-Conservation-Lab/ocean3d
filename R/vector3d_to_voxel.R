# SpatVector3D polygons to a presence SpatVoxel on a voxel's levels (#115).
# Each level is filled on its own, so a depth gap between two polygons over
# one cell stays a gap; vector3d_to_envelope() merges them into one range.

#' Convert 3D polygons to a voxel on a template's grid and depth levels
#'
#' Puts [SpatVector3D-class] polygons onto the grid and depth levels of
#' `template` as a presence [SpatVoxel-class]: `1` at every cell and level a
#' polygon covers, `NA` elsewhere. A polygon covers a cell when it touches it,
#' as in [vector3d_to_envelope()], and covers a level when its depth range
#' shares some thickness with the slab of water that level stands for, as in
#' [envelope_to_voxel()]. The slabs come from the template's depth edges (see
#' [depth_bounds()]), or from `bounds` when it has none.
#'
#' Each level is filled on its own, so where two polygons cover one cell at
#' different depths, the levels between them stay empty: polygons from 0 to
#' 50 m and from 200 to 300 m over the same cell cover the levels in those
#' two ranges and not the ones between. [vector3d_to_envelope()] gives that
#' cell one range, 0 to 300 m, because an envelope holds one interval per
#' cell. Use this function
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
  .check_geomtype(x)
  .check_3d(template, "template", "SpatVoxel")
  f <- x@features
  if (nrow(f) == 0L) {
    .abort_voxel_input("`x` has no polygons to put on the voxel's levels.")
  }
  .inform_time_dropped(x)

  depths <- depths(template)
  edges <- depth_bounds(template)
  slab <- .depth_slabs(depths, bounds, edges)
  # Only the grid is used; the class tag would reach the depth-aware methods.
  grid <- methods::as(template[[1]], "SpatRaster")  # remove the 3D class
  polygons <- .outlines_on_grid(x, grid)

  # The polygons covering each level, by the rule envelope_to_voxel() uses,
  # so one polygon lands on the same levels as its range as an envelope. A
  # polygon with no thickness covers nothing, as in vector3d_to_envelope().
  thick <- f$depth_max > f$depth_min
  covers <- lapply(seq_along(depths), function(i) {
    which(thick & .in_slab(f$depth_min, f$depth_max, slab, i,
                           closed = is.null(edges)))
  })
  # Each set of polygons is rasterized once: a polygon set with one depth
  # range covers a run of levels, often dozens on a World Ocean Atlas grid.
  sets <- unique(covers)
  footprints <- lapply(sets, function(rows) {
    if (!length(rows)) return(terra::setValues(terra::rast(grid), 0))
    .footprint(polygons[polygons$row %in% rows], grid)
  })
  occ <- terra::rast(footprints[match(covers, sets)])
  if (!any(terra::global(occ, "max", na.rm = TRUE)[, 1] > 0)) {
    .abort_envelope_empty(paste0(
      "The polygons cover no cell and depth level of `template`: they miss ",
      "its grid, have no thickness, or lie outside its depth levels."
    ))
  }
  .presence_voxel(occ, depths, edges)
}
