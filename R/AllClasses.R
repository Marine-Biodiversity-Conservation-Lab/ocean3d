# R/AllClasses.R

#' @importFrom methods setClass setClassUnion setValidity new validObject is
#' @importFrom methods show
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

#' 3D points, lines and polygons with depth and time
#'
#' A `SpatVector3D` holds features of one geometry type, reported by
#' [geomtype()]: `"points"` (a detection at 35 m), `"lines"` (a dive track or
#' a vessel route) or `"polygons"` (a marine protected area from 0 to 200 m).
#' Coordinates are longitude and latitude (WGS84), depths are positive metres
#' downward, and times are UTC. Build one with [as_vector3d()] from a table,
#' sf object or `terra::SpatVector`.
#'
#' Points and lines have a depth and optional time at each vertex. Polygons
#' have one depth range (`depth_min` to `depth_max`) and an optional time
#' interval for the whole shape.
#'
#' It is the vector counterpart of [SpatVoxel-class] and [SpatEnvelope-class].
#' Points and lines are not stored as `terra::SpatVector` or sf objects:
#' those keep attributes per feature, so a line's depth and time at each
#' vertex would not survive their operations. A polygon's depth range and
#' time are per feature, so polygons keep their outlines in an sf geometry
#' column.
#'
#' @slot type `"points"`, `"lines"` or `"polygons"`.
#' @slot vertices Points and lines: a data frame with one row per vertex:
#'   `feature_id`, `vertex` (the order within the feature), `lon`, `lat`,
#'   `depth`, `time` (POSIXct, UTC, or `NA` when the features have no time),
#'   and any other columns from the input. Empty for polygons.
#' @slot features Polygons: a data frame with one row per polygon:
#'   `feature_id`, `depth_min`, `depth_max`, `time_start` and `time_end`
#'   (POSIXct, UTC, or `NA` when the polygons have no time), `geometry` (an
#'   sf geometry column of POLYGON or MULTIPOLYGON, EPSG:4326), and any other
#'   columns from the input. Empty for points and lines.
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
  features = "data.frame",
  time_resolution = "character"
))

setValidity("SpatVector3D", function(object) {
  types <- c("points", "lines", "polygons")
  if (length(object@type) != 1L || !object@type %in% types) {
    return("type must be \"points\", \"lines\" or \"polygons\"")
  }
  if (object@type == "polygons") {
    f <- object@features
    needed <- c("feature_id", "depth_min", "depth_max", "time_start",
                "time_end", "geometry")
    if (!all(needed %in% names(f))) {
      return(paste("features must have columns", toString(needed)))
    }
    if (anyNA(f$feature_id) || any(!nzchar(f$feature_id))) {
      return("feature IDs must be non-missing and non-empty")
    }
    if (anyDuplicated(f$feature_id)) {
      return("each polygon needs its own ID; IDs repeat")
    }
    if (!inherits(f$geometry, "sfc") ||
          !all(sf::st_geometry_type(f$geometry) %in%
                 c("POLYGON", "MULTIPOLYGON"))) {
      return("geometry must be an sf column of POLYGON or MULTIPOLYGON")
    }
    if (sf::st_crs(f$geometry) != sf::st_crs(4326)) {
      return("geometry must be in WGS84 (EPSG:4326)")
    }
    if (any(sf::st_is_empty(f$geometry))) {
      return("every polygon needs an outline; some geometries are empty")
    }
    if (any(!is.finite(f$depth_min)) || any(!is.finite(f$depth_max))) {
      return("every polygon needs a depth_min and depth_max")
    }
    if (any(f$depth_min < 0)) {
      return("depths are positive metres increasing downward")
    }
    if (any(f$depth_min > f$depth_max)) {
      return("depth_min must not be deeper than depth_max")
    }
    if (!inherits(f$time_start, "POSIXct") ||
          !inherits(f$time_end, "POSIXct")) {
      return("time_start and time_end must be POSIXct (UTC)")
    }
    if (any(is.na(f$time_start) != is.na(f$time_end))) {
      return("a polygon's time needs both a start and an end, or neither")
    }
    if (any(f$time_start > f$time_end, na.rm = TRUE)) {
      return("time_start must not be after time_end")
    }
    return(TRUE)
  }
  v <- object@vertices
  needed <- c("feature_id", "vertex", "lon", "lat", "depth", "time")
  if (!all(needed %in% names(v))) {
    return(paste("vertices must have columns", toString(needed)))
  }
  if (anyNA(v$feature_id) || any(!nzchar(v$feature_id))) {
    return("feature IDs must be non-missing and non-empty")
  }
  if (object@type == "points" && anyDuplicated(v$feature_id)) {
    return("each point needs its own ID; IDs repeat")
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
#' depth range have volume: voxels, envelopes and polygons, not points or
#' lines.
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
