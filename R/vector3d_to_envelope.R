# SpatVector3D polygons to a SpatEnvelope on a template grid (#114).

#' Convert 3D polygons to an envelope on a grid
#'
#' Puts [SpatVector3D-class] polygons onto the grid of `template` as a
#' [SpatEnvelope-class]: each cell a polygon touches gets that polygon's
#' depth range, as [vect_to_envelope()] does for 2D polygons with one depth
#' range. Cells no polygon touches are `NA`, and so are cells where
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
#'   envelope uses, such as a [SpatVoxel-class] or a raster from
#'   [create_study_raster()]. Only its grid is used, not its values. Polygons
#'   are projected to its CRS if it isn't WGS84.
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
    rlang::abort(glue::glue("`x` must be SpatVector3D polygons.{hint}"),
                 class = "ocean3d_error_vector3d_to_envelope", call = NULL)
  }
  if (!inherits(template, "SpatRaster")) {
    rlang::abort("`template` must be a SpatRaster that defines the grid.",
                 class = "ocean3d_error_vector3d_to_envelope", call = NULL)
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
  # Every cell a polygon touches, as terra::mask() picks them in
  # vect_to_envelope().
  top <- terra::rasterize(polygons, grid, field = "depth_min", fun = "min",
                          touches = TRUE)
  bottom <- terra::rasterize(polygons, grid, field = "depth_max", fun = "max",
                             touches = TRUE)
  if (terra::global(top, "notNA")[1, 1] == 0) {
    rlang::abort(
      paste0("The polygons don't touch any cell of `template`. Check that ",
             "they overlap it."),
      class = "ocean3d_error_vector3d_to_envelope", call = NULL
    )
  }
  # A cell with no thickness has no volume, so it is dropped, as in
  # vect_to_envelope().
  thick <- terra::ifel(bottom > top, 1, NA)
  if (terra::global(thick, "notNA")[1, 1] == 0) {
    rlang::abort(
      paste0("Every cell the polygons touch has `depth_min` equal to ",
             "`depth_max`, so the envelope would be empty. Give the ",
             "polygons a depth range, e.g. a `depth_window` above 0 in ",
             "buffer_3d()."),
      class = "ocean3d_error_vector3d_to_envelope", call = NULL
    )
  }
  layers <- terra::mask(c(top, bottom), thick)
  names(layers) <- c("depth_min", "depth_max")
  as_envelope(layers)
}
