# R/AllClasses.R

#' @importFrom methods setClass setClassUnion setValidity new validObject is
#' @importClassesFrom terra SpatRaster
NULL

#' 3D voxel model: one layer per standard depth level
#'
#' A `SpatVoxel` is a [terra::SpatRaster] in which **depth is the layer index**
#' and the cell values are the variable (temperature, oxygen, presence, ...).
#' Layer names follow the `{variable}_depth={value}` convention and must be
#' distinct and ordered shallow to deep. The depth axis is grid-wide: every
#' cell is sampled at the same set of standard depths.
#'
#' A voxel may have interior gaps — a cell can be NA at one depth and non-NA
#' at the depths above and below it.
#'
#' @seealso [SpatEnvelope-class], [SpatVolume-class]
#'
#' @examples
#' # Two cells sampled at three standard depths. Build one with as_voxel(),
#' # which names the layers for you and sorts them shallow to deep.
#' r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
#'                  ymax = 1, nlyrs = 3)
#' terra::values(r) <- cbind(c(12, 11), c(9, NA), c(6, 5))
#' v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")
#'
#' names(v)          # depth lives in the layer name
#' terra::values(v)  # ...and the variable in the cell values
#'
#' # Note cell 2: NA at 100 m, but with values above and below it. A voxel can
#' # hold that interior gap; a SpatEnvelope cannot.
#'
#' # It is an ordinary SpatRaster underneath, so terra operations work on it
#' # unchanged. Re-wrapping with as_voxel() is cheap either way, since the
#' # constructor is idempotent.
#' methods::is(v, "SpatRaster")
#' terra::nlyr(terra::crop(v, terra::ext(0, 1, 0, 1)))
#' identical(names(as_voxel(v)), names(v))
#'
#' # Validity is enforced on construction: every layer name must carry its
#' # depth, and the layers must run shallow to deep.
#' bad <- r
#' names(bad) <- c("a", "b", "c")
#' try(methods::new("SpatVoxel", bad))
#' @export
setClass("SpatVoxel", contains = "SpatRaster")

#' 2.5D min-max envelope: exactly depth_min, depth_max
#'
#' A `SpatEnvelope` is a [terra::SpatRaster] with exactly two layers,
#' `depth_min` and `depth_max`, in which **depth is the cell value**. The
#' variable is the object itself: the envelope delimits the vertical extent of
#' a spatial phenomenon (typically a species range) rather than sampling a
#' variable within it.
#'
#' The vertical interval is per-cell and continuous, and is solid by
#' construction — an envelope cannot represent a gap in the vertical
#' distribution. Depths are positive metres increasing downward, and
#' `depth_max` is at least `depth_min` in every cell.
#'
#' @seealso [SpatVoxel-class], [SpatVolume-class]
#'
#' @examples
#' # A two-cell footprint, one cell present and one absent. Build an envelope
#' # with as_envelope() (from a raster footprint) or vect_to_envelope() (from
#' # polygons).
#' fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
#'                   ymax = 1)
#' terra::values(fp) <- c(1, NA)
#' e <- as_envelope(fp, depth_min = 50, depth_max = 200)
#'
#' names(e)          # always exactly these two layers
#' terra::values(e)  # ...and depth is the cell value, in metres
#'
#' # The interval is per-cell, so limits can vary across the grid.
#' terra::values(
#'   as_envelope(fp, depth_min = 50,
#'               depth_max = terra::setValues(terra::rast(fp), c(120, 300)))
#' )
#'
#' # Any other set of layers is rejected: a bare footprint is not an envelope.
#' try(methods::new("SpatEnvelope", fp))
#' @export
setClass("SpatEnvelope", contains = "SpatRaster")

#' 3D points and lines with depth and time
#'
#' A `SpatVector3D` holds features of one geometry type, reported by
#' [geomtype()]: `"points"` (a detection at 35 m) or `"lines"` (a dive track
#' or a vessel route). Each vertex has a longitude and latitude (WGS84), a
#' depth in positive metres downward, and optionally a time (UTC). Build one
#' with [as_vector3d()]; convert to sf with `sf::st_as_sf()`.
#'
#' It is the vector counterpart of [SpatVoxel-class] and [SpatEnvelope-class].
#' It is not built on `terra::SpatVector` or sf: those keep attributes per
#' feature, so a line's depth and time at each vertex would not survive their
#' operations.
#'
#' @slot type `"points"` or `"lines"`.
#' @slot vertices A data frame with one row per vertex: `feature_id`, `vertex`
#'   (the order within the feature), `lon`, `lat`, `depth`, `time` (POSIXct,
#'   UTC, or `NA` when the features have no time), and any other columns from
#'   the input.
#' @slot time_resolution `"daily"` or `"datetime"` for the precision of the
#'   times, as for [extract_to_point()], or `NA` when there is no time.
#'
#' @seealso [as_vector3d()], [geomtype()]
#' @examples
#' obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35))
#' pts <- as_vector3d(obs)
#' geomtype(pts)
#' pts
#' @export
setClass("SpatVector3D", slots = c(
  type = "character",
  vertices = "data.frame",
  time_resolution = "character"
))

setValidity("SpatVector3D", function(object) {
  v <- object@vertices
  if (length(object@type) != 1L || !object@type %in% c("points", "lines")) {
    return("type must be \"points\" or \"lines\"")
  }
  needed <- c("feature_id", "vertex", "lon", "lat", "depth", "time")
  if (!all(needed %in% names(v))) {
    return(paste("vertices must have columns", toString(needed)))
  }
  if (any(!is.finite(v$lon)) || any(!is.finite(v$lat)) ||
        any(abs(v$lat) > 90)) {
    return(paste("longitudes and latitudes must be finite,",
                 "latitudes within [-90, 90]"))
  }
  if (any(!is.finite(v$depth))) return("every vertex needs a depth")
  if (any(v$depth < 0)) return("depths are positive metres increasing downward")
  if (!inherits(v$time, "POSIXct")) return("time must be POSIXct (UTC)")
  if (object@type == "lines" && any(table(v$feature_id) < 2L)) {
    return("each line needs at least two vertices")
  }
  TRUE
})

#' Union of the package's 3D domain representations
#'
#' `SpatVolume` is a class union over [SpatVoxel-class],
#' [SpatEnvelope-class] and [SpatVector3D-class]. It exists as a dispatch
#' target for operations that are meaningful on any 3D representation, such
#' as volume, volumetric overlap and vertical extent. Only shapes with a
#' depth range have volume: voxels, envelopes and (later) polygons, not
#' points or lines.
#'
#' It deliberately asserts no shared structure. The raster members store depth
#' in dual roles (layer index vs. cell value), and `SpatVector3D` stores it per
#' vertex, so any operation that touches the depth axis directly should
#' dispatch on the concrete class instead.
#'
#' @seealso [SpatVoxel-class], [SpatEnvelope-class] and [SpatVector3D-class],
#'   the members; [envelope_to_voxel()] and [voxel_to_envelope()] convert
#'   between the raster members.
#'
#' @examples
#' # The same 3D domain in both representations.
#' fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
#'                   ymax = 1)
#' terra::values(fp) <- c(1, 1)
#' e <- as_envelope(fp, depth_min = 0, depth_max = 200)
#' v <- envelope_to_voxel(e, depths = c(0, 100, 200))
#'
#' # Neither inherits from the other; both belong to the union.
#' methods::is(e, "SpatVolume")
#' methods::is(v, "SpatVolume")
#' methods::is(v, "SpatEnvelope")
#'
#' # `SpatVolume` is a dispatch target, so a generic written against it accepts
#' # either representation.
#' setGeneric("n_depth_layers", function(x) standardGeneric("n_depth_layers"))
#' setMethod("n_depth_layers", "SpatVolume", function(x) terra::nlyr(x))
#' n_depth_layers(e)
#' n_depth_layers(v)
#'
#' # But the depth axis itself is stored differently in the two, so anything
#' # reading depths must dispatch on the concrete class rather than the union.
#' names(e)
#' names(v)
#' @export
setClassUnion("SpatVolume", c("SpatVoxel", "SpatEnvelope", "SpatVector3D"))

setValidity("SpatVoxel", function(object) {
  d <- .parse_depth_layers(object, error = FALSE)
  if (anyNA(d)) return("every layer name must follow {variable}_depth={value}")
  if (is.unsorted(d)) return("depth layers must be ordered shallow to deep")
  if (anyDuplicated(d)) return("depth layers must be distinct")
  if (any(d < 0)) return("depths are positive metres increasing downward")
  TRUE
})

setValidity("SpatEnvelope", function(object) {
  if (!identical(names(object), c("depth_min", "depth_max")))
    return("layers must be exactly: depth_min, depth_max")

  # Depth is the cell value here, not the layer name, so the interval rules
  # can only be checked by reading the values. That costs a pass over the
  # grid, which the sibling class avoids because its depths live in the layer
  # names. It is worth paying: without it a tagged raster that terra produced
  # passes as a valid envelope, and volume() returns a negative number.
  min_of <- function(r) {
    v <- terra::global(r, "min", na.rm = TRUE)[1, 1]
    if (is.null(v) || is.na(v)) NA_real_ else v
  }
  shallowest <- min_of(.envelope_layer(object, "depth_min"))
  if (!is.na(shallowest) && shallowest < 0)
    return("depths are positive metres increasing downward")
  thinnest <- min_of(.envelope_layer(object, "depth_max") -
                       .envelope_layer(object, "depth_min"))
  if (!is.na(thinnest) && thinnest < 0)
    return("depth_max must be at least depth_min in every cell")
  TRUE
})
