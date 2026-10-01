# R/intersect.R
#
# Spatial query for the package's 3D representations: intersect_3d() gives
# the domain two objects share, intersects_3d() says whether there is one,
# and depth-aware mask() methods keep one object's values inside another's
# domain. The first two generics are declared in R/AllGenerics.R; mask()
# attaches to terra's generic.
#
# Dispatch is on concrete classes only. A method on the SpatVolume class union
# loses to a SpatRaster method for the same generic (S4 resolves that tie in
# favour of the ordinary superclass), so every signature here names
# SpatEnvelope, SpatVoxel, SpatRaster or ANY. The ANY methods take the 2D
# side of a query -- a footprint raster or polygons -- and reject anything
# else with the list of accepted types.
#
# terra propagates the subclass through its operations, so every result that
# is not a 3D domain is stripped with .as_plain_raster() before it is returned,
# and every 3D result is rebuilt through as_envelope() or as_voxel() so the
# validity rules run.

# Internal: the error for a query in which neither side is a 3D object.
.stop_neither_3d <- function() {
  stop("neither `x` nor `y` is a SpatEnvelope or a SpatVoxel. A 3D query ",
       "needs at least one 3D object; build one with as_envelope(), ",
       "vect_to_envelope() or as_voxel().", call. = FALSE)
}

# Internal: a voxel's occupancy as a presence domain: 1 where occupied, NA
# elsewhere, on the depth levels of `x`.
.presence_voxel <- function(occ, depths) {
  out <- terra::ifel(occ > 0, 1, NA)
  as_voxel(.as_plain_raster(out), depths = depths, varname = "presence")
}

# ---- intersect_3d -----------------------------------------------------------

#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatEnvelope", "SpatEnvelope"),
  function(x, y, ...) {
    .no_extra_args(...)
    .check_3d(x, "x")
    .check_3d(y, "y")
    .check_same_grid(x, y)

    x_min <- .envelope_layer(x, "depth_min")
    x_max <- .envelope_layer(x, "depth_max")
    y_min <- .envelope_layer(y, "depth_min")
    y_max <- .envelope_layer(y, "depth_max")

    # The deeper of the two tops and the shallower of the two bottoms. Where
    # either envelope is absent the comparison is NA and the cell stays NA.
    dmin <- terra::ifel(x_min > y_min, x_min, y_min)
    dmax <- terra::ifel(x_max < y_max, x_max, y_max)

    # Strict `>`: intervals that only touch share no water.
    keep <- terra::ifel(dmax > dmin, 1, NA)
    out <- c(terra::mask(dmin, keep), terra::mask(dmax, keep))
    names(out) <- c("depth_min", "depth_max")
    as_envelope(out)
  }
)

#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatVoxel", "SpatVoxel"),
  function(x, y) {
    .check_3d(x, "x")
    .check_3d(y, "y")
    .check_same_grid(x, y)
    depths <- .shared_depths(x, y)

    occ <- .voxel_occupancy(x) * .voxel_occupancy(y)
    .presence_voxel(occ, depths)
  }
)

# Mixed input: the envelope goes onto the voxel's levels, never the reverse,
# because voxel_to_envelope() fills interior gaps and would overstate the
# shared domain.
#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatEnvelope", "SpatVoxel"),
  function(x, y, bounds = c("top", "midpoint"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(y, "y")
    intersect_3d(.promote_to_voxel(x, y, bounds), y, ...)
  }
)

#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatVoxel", "SpatEnvelope"),
  function(x, y, bounds = c("top", "midpoint"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(y, "y")
    intersect_3d(x, .promote_to_voxel(y, x, bounds), ...)
  }
)

# 2D input: a footprint raster or polygons restrict the domain horizontally
# and leave its depths alone.
#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatEnvelope", "ANY"),
  function(x, y, ...) {
    .no_extra_args(...)
    .check_3d(x, "x")
    fp <- .footprint_of(y, x)
    as_envelope(terra::mask(.as_plain_raster(x), fp))
  }
)

#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatVoxel", "ANY"),
  function(x, y) {
    .check_3d(x, "x")
    fp <- .footprint_of(y, x)
    occ <- terra::mask(.voxel_occupancy(x), fp, updatevalue = 0)
    .presence_voxel(occ, .parse_depth_layers(x))
  }
)

# Intersection is symmetric, so a 2D object in first position is the same
# query the other way round.
#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("ANY", "SpatEnvelope"),
          function(x, y, ...) intersect_3d(y, x, ...))

#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("ANY", "SpatVoxel"),
          function(x, y, ...) intersect_3d(y, x, ...))

#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("ANY", "ANY"),
          function(x, y, ...) .stop_neither_3d())

# ---- intersects_3d ----------------------------------------------------------

# Internal: the tri-state answer from a per-cell hit layer and the two
# presence layers. `hit` is NA wherever either side is absent, so it is
# flattened to FALSE before the presence logic runs; "neither present" is
# applied last because it is the only case that stays NA.
.intersects_answer <- function(hit, present_x, present_y) {
  hit <- terra::ifel(is.na(hit), FALSE, hit)
  out <- terra::ifel(!present_x & !present_y, NA, present_x & present_y & hit)
  # ifel() yields 1/0/NA; only as.bool() makes values() give TRUE/FALSE/NA.
  out <- terra::as.bool(out)
  names(out) <- "intersects"
  .as_plain_raster(out)
}

#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatEnvelope", "SpatEnvelope"),
  function(x, y, ...) {
    .no_extra_args(...)
    .check_3d(x, "x")
    .check_3d(y, "y")
    .check_same_grid(x, y)

    x_min <- .envelope_layer(x, "depth_min")
    x_max <- .envelope_layer(x, "depth_max")
    y_min <- .envelope_layer(y, "depth_min")
    y_max <- .envelope_layer(y, "depth_max")

    overlap_min <- terra::ifel(x_min > y_min, x_min, y_min)
    overlap_max <- terra::ifel(x_max < y_max, x_max, y_max)
    # Strict `>`: intervals that only touch do not intersect.
    .intersects_answer(overlap_max > overlap_min,
                       .envelope_present(x), .envelope_present(y))
  }
)

#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatVoxel", "SpatVoxel"),
  function(x, y) {
    .check_3d(x, "x")
    .check_3d(y, "y")
    .check_same_grid(x, y)
    .shared_depths(x, y)

    occ_x <- .voxel_occupancy(x)
    occ_y <- .voxel_occupancy(y)
    .intersects_answer(sum(occ_x * occ_y) > 0, sum(occ_x) > 0, sum(occ_y) > 0)
  }
)

#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatEnvelope", "SpatVoxel"),
  function(x, y, bounds = c("top", "midpoint"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(y, "y")
    intersects_3d(.promote_to_voxel(x, y, bounds), y, ...)
  }
)

#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatVoxel", "SpatEnvelope"),
  function(x, y, bounds = c("top", "midpoint"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(y, "y")
    intersects_3d(x, .promote_to_voxel(y, x, bounds), ...)
  }
)

# 2D input has no depth axis, so the test is whether both are present.
#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatEnvelope", "ANY"),
  function(x, y, ...) {
    .no_extra_args(...)
    .check_3d(x, "x")
    present_x <- .envelope_present(x)
    present_y <- !is.na(.footprint_of(y, x))
    .intersects_answer(present_x & present_y, present_x, present_y)
  }
)

#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatVoxel", "ANY"),
  function(x, y) {
    .check_3d(x, "x")
    present_x <- sum(.voxel_occupancy(x)) > 0
    present_y <- !is.na(.footprint_of(y, x))
    .intersects_answer(present_x & present_y, present_x, present_y)
  }
)

#' @rdname intersects_3d
#' @export
setMethod("intersects_3d", c("ANY", "SpatEnvelope"),
          function(x, y, ...) intersects_3d(y, x, ...))

#' @rdname intersects_3d
#' @export
setMethod("intersects_3d", c("ANY", "SpatVoxel"),
          function(x, y, ...) intersects_3d(y, x, ...))

#' @rdname intersects_3d
#' @export
setMethod("intersects_3d", c("ANY", "ANY"),
          function(x, y, ...) .stop_neither_3d())

# ---- mask -------------------------------------------------------------------

#' Keep a raster's values inside a 3D domain
#'
#' These methods make [terra::mask()] depth-aware. `mask(x, mask)` keeps the
#' values of `x` where `mask` is present and sets the rest to `NA`. When
#' `mask` is a [SpatEnvelope-class] or a [SpatVoxel-class], "present" is
#' decided cell by cell **and** depth by depth. A value of `x` stays only
#' where the mask domain reaches that cell at that depth.
#'
#' What each pairing does:
#' \describe{
#'   \item{voxel masked by an envelope}{The envelope is first placed on the
#'     voxel's own depth levels with [envelope_to_voxel()]. Each depth layer
#'     of the voxel is then masked by the matching level. This replaces the
#'     hand-written `terra::mask(v, envelope_to_voxel(e, depths(v)))` and
#'     cannot be pointed at the wrong depths.}
#'   \item{voxel masked by a voxel}{Both must be sampled at the same depth
#'     levels. Each depth layer of `x` is masked by the matching layer of the
#'     mask, wherever the mask is occupied.}
#'   \item{envelope masked by an envelope or a voxel}{A cell of `x` is kept
#'     where [intersects_3d()] finds the two domains overlap. Its depth
#'     interval is kept whole. To narrow the interval instead, use
#'     [intersect_3d()].}
#'   \item{plain `SpatRaster` masked by an envelope or a voxel}{A cell is kept
#'     where the mask domain is present at any depth. Without this method
#'     terra would read the mask's first layer only.}
#' }
#'
#' A voxel or envelope masked by a plain `SpatRaster` or by polygons is
#' terra's own `mask()`. Every layer is masked by the same 2D pattern and the
#' class is kept. No method is added for those cases.
#'
#' @param x The raster whose values are kept: a [SpatVoxel-class], a
#'   [SpatEnvelope-class], or a plain `SpatRaster`.
#' @param mask The 3D domain to keep `x` inside: a [SpatEnvelope-class] or a
#'   [SpatVoxel-class] on the same grid as `x`.
#' @param bounds Only when one side is an envelope and the other a voxel. How
#'   the envelope is placed on the voxel's depth levels: `"top"` (default) or
#'   `"midpoint"`. See [envelope_to_voxel()].
#' @param ... Passed on to [terra::mask()]. For example, `inverse = TRUE`
#'   keeps the values *outside* the domain instead.
#'
#' @returns `x`, with the same class and layers, with every value outside the
#'   domain set to `NA`.
#'
#' @seealso [intersect_3d()] for the shared domain itself; [intersects_3d()]
#'   for the yes/no test; [extract_to_area()] to restrict a voxel to a polygon
#'   and one depth band instead of a per-cell depth window.
#'
#' @examples
#' fp <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3000,
#'                   ymin = 0, ymax = 1000,
#'                   crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
#' terra::values(fp) <- c(1, 1, 1)
#'
#' # A temperature voxel on three levels, and a species found at 0-100 m in
#' # cells 1-2 only.
#' temp <- as_voxel(terra::rast(list(terra::setValues(fp, c(20, 21, 22)),
#'                                   terra::setValues(fp, c(15, 16, 17)),
#'                                   terra::setValues(fp, c(10, 11, 12)))),
#'                  depths = c(0, 100, 200), varname = "temp")
#' range_env <- as_envelope(terra::setValues(fp, c(1, 1, NA)),
#'                          depth_min = 0, depth_max = 100)
#'
#' # The temperatures inside the species' range: cell 3 and the 200 m level
#' # are gone.
#' terra::values(mask(temp, range_env))
#'
#' # A 2D raster masked by the range keeps cells where the range is present
#' # at any depth.
#' terra::values(mask(fp, range_env))
#' @name mask-3d
NULL

# Internal: terra's own mask(), reached on plain copies so a method here never
# dispatches back to itself.
.mask_plain <- function(x, keep, ...) {
  terra::mask(.as_plain_raster(x), .as_plain_raster(keep), ...)
}

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatVoxel", "SpatEnvelope"),
  function(x, mask, bounds = c("top", "midpoint"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(mask, "mask")
    .check_same_grid(x, mask, c("x", "mask"))
    keep <- .promote_to_voxel(mask, x, bounds)
    as_voxel(.mask_plain(x, keep, ...))
  }
)

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatVoxel", "SpatVoxel"),
  function(x, mask, ...) {
    .check_3d(x, "x")
    .check_3d(mask, "mask")
    .check_same_grid(x, mask, c("x", "mask"))
    .shared_depths(x, mask)
    keep <- .presence_voxel(.voxel_occupancy(mask),
                            .parse_depth_layers(mask))
    as_voxel(.mask_plain(x, keep, ...))
  }
)

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatEnvelope", "SpatEnvelope"),
  function(x, mask, ...) {
    .check_3d(x, "x")
    .check_3d(mask, "mask")
    keep <- terra::ifel(intersects_3d(x, mask), 1, NA)
    as_envelope(.mask_plain(x, keep, ...))
  }
)

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatEnvelope", "SpatVoxel"),
  function(x, mask, bounds = c("top", "midpoint"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(mask, "mask")
    keep <- terra::ifel(intersects_3d(x, mask, bounds = bounds), 1, NA)
    as_envelope(.mask_plain(x, keep, ...))
  }
)

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatRaster", "SpatEnvelope"),
  function(x, mask, ...) {
    .check_3d(mask, "mask")
    .check_same_grid(x, mask, c("x", "mask"))
    terra::mask(x, .footprint(mask), ...)
  }
)

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatRaster", "SpatVoxel"),
  function(x, mask, ...) {
    .check_3d(mask, "mask")
    .check_same_grid(x, mask, c("x", "mask"))
    terra::mask(x, .footprint(mask), ...)
  }
)