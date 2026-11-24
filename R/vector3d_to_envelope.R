# SpatVector3D polygons to a SpatEnvelope on a template grid (#114).

#' Convert 3D polygons to an envelope on a grid
#'
#' Puts [SpatVector3D-class] polygons onto the grid of `template` as a
#' [SpatEnvelope-class], as [vect_to_envelope()] does for 2D polygons with
#' one depth range. Each cell a polygon touches gets that polygon's depth
#' range. Cells that no polygon touches are `NA`, and so are cells where
#' `depth_min` equals `depth_max`, which have no thickness.
#'
#' Where polygons overlap, a cell gets the shallowest `depth_min` and the
#' deepest `depth_max` of the polygons that cover it. An envelope has one
#' depth range per cell, so polygons from 0 to 50 m and from 200 to 300 m
#' over the same cell give 0 to 300 m there.
#'
#' An envelope has no time, so the polygons' time intervals are dropped, with
#' a message when any polygon has one. To make an envelope for one period,
#' keep only the polygons for that period first.
#'
#' @param x A [SpatVector3D-class] of polygons, from [as_vector3d()] or
#'   [buffer_3d()].
#' @param template A `SpatRaster` whose grid (extent, resolution and CRS) the
#'   envelope uses; its values are ignored. For example, a [SpatVoxel-class]
#'   or a raster from [create_study_raster()]. Polygons are projected to its
#'   CRS if it isn't WGS84. Longitudes from 0 to 360 work too.
#'
#' @return A [SpatEnvelope-class] on the grid of `template`. Stops with an
#'   error if it would have no cells: when the polygons don't overlap
#'   `template`, or none of them has `depth_max` deeper than `depth_min`.
#' @seealso [vect_to_envelope()] for 2D polygons, [as_envelope()]
#' @examples
#' template <- terra::rast(nrows = 4, ncols = 4, xmin = 150, xmax = 154,
#'                         ymin = -24, ymax = -20, crs = "EPSG:4326")
#' square <- function(x0, y0) {
#'   sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 2, y0), c(x0 + 2, y0 + 2),
#'                             c(x0, y0 + 2), c(x0, y0))))
#' }
#' areas <- sf::st_sf(
#'   top = c(0, 200), bottom = c(50, 300),
#'   geometry = sf::st_sfc(square(151, -23), square(152, -24), crs = 4326)
#' )
#' polygons <- as_vector3d(areas, depth_min = "top", depth_max = "bottom")
#'
#' # The overlapping cell gets 0 to 300 m.
#' e <- vector3d_to_envelope(polygons, template)
#' terra::values(e)
#' @export
vector3d_to_envelope <- function(x, template) {
  if (!methods::is(x, "SpatVector3D")) {
    .abort_envelope_input("`x` must be SpatVector3D polygons.")
  }
  .check_geomtype(x)
  if (!inherits(template, "SpatRaster")) {
    .abort_envelope_input(
      "`template` must be a SpatRaster that defines the grid."
    )
  }
  f <- x@features
  if (nrow(f) == 0L) {
    .abort_envelope_input("`x` has no polygons to put on the grid.")
  }
  if (all(f$depth_min == f$depth_max)) {
    .abort_envelope_empty(paste0(
      "Every polygon has `depth_min` equal to `depth_max`, so the envelope ",
      "would be empty. Give them a depth range, for example a ",
      "`depth_window` above 0 in buffer_3d()."
    ))
  }
  .inform_time_dropped(x)

  # Only the grid is used. Dropping the class tag keeps a voxel template's
  # layers away from the depth-aware methods in R/intersect.R.
  grid <- methods::as(template[[1]], "SpatRaster")  # remove the 3D class
  polygons <- .outlines_on_grid(x, grid, c("depth_min", "depth_max"))
  # Each cell's shallowest top and deepest bottom over the polygons touching
  # it, picked as terra::mask() picks cells in vect_to_envelope(). terra
  # rasterizes polygon values as 4-byte floats (since 1.9-60 whatever the
  # datatype asked for), which would turn 10.1 m into 10.1000004 m. So each
  # polygon's rank among the depths is rasterized instead: ranks are whole
  # numbers, which floats hold exactly, and the lowest or highest rank picks
  # the same polygon as the lowest or highest depth. The depths are then
  # looked up as doubles.
  extreme <- function(field, fun) {
    given <- terra::values(polygons)[[field]]
    levels <- sort(unique(given))
    polygons$rank <- match(given, levels)
    rank <- terra::rasterize(polygons, grid, field = "rank", fun = fun,
                             touches = TRUE)
    terra::setValues(rank, levels[terra::values(rank)[, 1]])
  }
  top <- extreme("depth_min", "min")
  bottom <- extreme("depth_max", "max")
  # vect_to_envelope() keeps the cells the polygons touch, drops cells with
  # no thickness, and stops if none are left.
  vect_to_envelope(polygons, grid, depth_min = top, depth_max = bottom)
}
