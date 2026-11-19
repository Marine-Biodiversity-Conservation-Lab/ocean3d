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
#' @return A [SpatEnvelope-class] on the grid of `template`.
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
  if (!methods::is(x, "SpatVector3D") || x@type != "polygons") {
    hint <- if (methods::is(x, "SpatVector3D")) {
      glue::glue(" These are {x@type}; buffer_3d() turns them into polygons.")
    } else {
      ""
    }
    .abort_envelope_input(
      glue::glue("`x` must be SpatVector3D polygons.{hint}")
    )
  }
  if (!inherits(template, "SpatRaster")) {
    .abort_envelope_input(
      "`template` must be a SpatRaster that defines the grid."
    )
  }
  if (terra::crs(template) == "") {
    .abort_envelope_input(paste0(
      "`template` has no CRS, so the polygons can't be placed on it. Set one ",
      "with terra::crs()."
    ))
  }
  f <- x@features
  if (any(!is.na(f$time_start))) {
    rlang::inform(
      paste0("The polygons' times are dropped: an envelope keeps only their ",
             "area and depth range."),
      class = "ocean3d_message_time_dropped"
    )
  }

  # Only the grid is used. Dropping the class tag keeps a voxel template's
  # layers away from the depth-aware methods in R/intersect.R.
  grid <- .as_plain_raster(template[[1]])
  polygons <- terra::vect(sf::st_sf(depth_min = f$depth_min,
                                    depth_max = f$depth_max,
                                    geometry = f$geometry))
  if (!terra::same.crs(polygons, grid)) {
    polygons <- terra::project(polygons, terra::crs(grid))
  }
  # Outlines are stored in [-180, 180]. On a longitude grid that runs past
  # 180 (0 to 360) or below -180, copies shifted by 360 degrees put them, and
  # both halves of a shape split at the antimeridian, where the grid is.
  if (terra::is.lonlat(grid)) {
    grid_x <- terra::ext(grid)[1:2]
    shifts <- c(if (grid_x[2] > 180) 360, if (grid_x[1] < -180) -360)
    polygons <- do.call(rbind, c(list(polygons), lapply(shifts, function(dx) {
      terra::shift(polygons, dx = dx)
    })))
  }
  # Each cell's shallowest top and deepest bottom over the polygons touching
  # it, picked as terra::mask() picks cells in vect_to_envelope(). Depths stay
  # doubles: terra's default of 4-byte floats would turn 10.1 m into
  # 10.1000004 m.
  doubles <- list(datatype = "FLT8S")
  top <- terra::rasterize(polygons, grid, field = "depth_min", fun = "min",
                          touches = TRUE, wopt = doubles)
  bottom <- terra::rasterize(polygons, grid, field = "depth_max", fun = "max",
                             touches = TRUE, wopt = doubles)
  # vect_to_envelope() keeps the cells the polygons touch, drops cells with
  # no thickness, and stops if none are left.
  vect_to_envelope(polygons, grid, depth_min = top, depth_max = bottom)
}
