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
# is not a 3D domain has the class removed with methods::as(x, "SpatRaster")
# before it is returned (see R/AllClasses.R), and every 3D result is rebuilt
# through as_envelope() or as_voxel() so the validity rules run.

# Internal: the one producer of `ocean3d_error_query_input`, for a pair no
# method takes: neither side is a 3D object, SpatVector3D polygons or points
# meet something other than an envelope, a voxel or a SpatVector3D, or two
# sets of points meet (points have no extent, so they share no space).
.stop_neither_3d <- function(x, y) {
  vector <- methods::is(x, "SpatVector3D") || methods::is(y, "SpatVector3D")
  both <- methods::is(x, "SpatVector3D") && methods::is(y, "SpatVector3D")
  message <- if (both && x@type == "points" && y@type == "points") {
    glue::glue(
      "Two sets of points share no space; compare points with polygons, an ",
      "envelope or a voxel, or buffer one set with buffer_3d()."
    )
  } else if (vector) {
    other <- if (methods::is(x, "SpatVector3D")) y else x
    glue::glue(
      "SpatVector3D polygons or points can be compared with a SpatEnvelope, ",
      "a SpatVoxel or a SpatVector3D (polygons, or points with polygons), ",
      "not with {class(other)[1L]}. Build 3D polygons from 2D ones with ",
      "as_vector3d()."
    )
  } else {
    glue::glue(
      "neither `x` nor `y` is a SpatEnvelope, a SpatVoxel or SpatVector3D ",
      "polygons or points. A 3D query needs at least one 3D object; build ",
      "one with as_envelope(), vect_to_envelope(), as_voxel() or ",
      "as_vector3d()."
    )
  }
  rlang::abort(message, class = "ocean3d_error_query_input", call = NULL)
}

# Internal: a voxel's occupancy as a presence domain: 1 where occupied, NA
# elsewhere, on the depth levels of `x`, carrying `depth_bounds` if given.
.presence_voxel <- function(occ, depths, depth_bounds = NULL) {
  out <- terra::ifel(occ > 0, 1, NA)
  # remove the 3D class
  as_voxel(methods::as(out, "SpatRaster"), depths = depths,
           varname = "presence", depth_bounds = depth_bounds)
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
    .check_same_grid(y, x)

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
    .check_same_grid(y, x)
    shared <- .shared_depths(x, y)

    occ <- .voxel_occupancy(x) * .voxel_occupancy(y)
    .presence_voxel(occ, shared$depths, shared$depth_bounds)
  }
)

# Mixed input: the envelope goes onto the voxel's levels, never the reverse,
# because voxel_to_envelope() fills interior gaps and would overstate the
# shared domain.
#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatEnvelope", "SpatVoxel"),
  function(x, y, bounds = c("midpoint", "top"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(y, "y")
    promoted <- .promote_to_voxel(x, depths(y), bounds, depth_bounds(y))
    intersect_3d(promoted, y, ...)
  }
)

#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatVoxel", "SpatEnvelope"),
  function(x, y, bounds = c("midpoint", "top"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(y, "y")
    promoted <- .promote_to_voxel(y, depths(x), bounds, depth_bounds(x))
    intersect_3d(x, promoted, ...)
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
    # Assigned before ifel() takes it: an error raised while S4 evaluates an
    # argument for dispatch is reworded.
    fp <- .footprint(y, x)
    fp <- terra::ifel(fp, 1, NA)
    # remove the 3D class
    as_envelope(terra::mask(methods::as(x, "SpatRaster"), fp))
  }
)

#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatVoxel", "ANY"),
  function(x, y) {
    .check_3d(x, "x")
    fp <- .footprint(y, x)
    fp <- terra::ifel(fp, 1, NA)
    occ <- terra::mask(.voxel_occupancy(x), fp, updatevalue = 0)
    .presence_voxel(occ, depths(x), depth_bounds(x))
  }
)

# Internal: whether `x` is a gridded 3D object. Not the SpatVolume union,
# which also holds SpatVector3D polygons.
.is_grid_3d <- function(x) {
  methods::is(x, "SpatEnvelope") || methods::is(x, "SpatVoxel")
}

# Intersection is symmetric, so a 2D object in first position is the same
# query the other way round. Both inputs are checked here first, so an error
# names the argument as the caller passed it. Any other pair has no 3D object
# to answer on.
#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("ANY", "ANY"), function(x, y, ...) {
  if (.is_grid_3d(y)) {
    .check_3d(y, "y")
    .footprint(x, y, "x")
    return(intersect_3d(y, x, ...))
  }
  .stop_neither_3d(x, y)
})

# ---- intersects_3d ----------------------------------------------------------

# Two objects intersect in a cell where intersect_3d() keeps any of it. Only
# the envelope/voxel pairs need methods of their own: the envelope goes onto
# the voxel's levels first, as in intersect_3d(), so that its footprint is
# the promoted one and `bounds` reaches the promotion.
#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatEnvelope", "SpatVoxel"),
  function(x, y, bounds = c("midpoint", "top"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(y, "y")
    promoted <- .promote_to_voxel(x, depths(y), bounds, depth_bounds(y))
    intersects_3d(promoted, y, ...)
  }
)

#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatVoxel", "SpatEnvelope"),
  function(x, y, bounds = c("midpoint", "top"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(y, "y")
    promoted <- .promote_to_voxel(y, depths(x), bounds, depth_bounds(x))
    intersects_3d(x, promoted, ...)
  }
)

# Every other pair, in either order. intersect_3d() checks the inputs and
# raises `ocean3d_error_query_input` for a pair with no 3D object, so it runs
# first. It builds the shared domain, which can be empty; an empty envelope
# or voxel is valid, so that costs one more validity check and nothing else.
#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("ANY", "ANY"),
  function(x, y, ...) {
    # Polygons go on the other object's grid (on a voxel's levels with
    # `bounds`, or as an envelope), and points on a voxel's grid as a
    # presence voxel, then take that pair's methods. Points with an envelope
    # take the general path below: intersect_3d() gives the part of the
    # envelope they reach. Assigned first: an error raised while S4
    # evaluates an argument for dispatch loses its class.
    if (methods::is(x, "SpatVector3D") && .is_grid_3d(y) &&
          .placed_on_grid(x, y)) {
      shape <- .vector3d_on_grid(x, y, ...)
      return(intersects_3d(shape, y))
    }
    if (methods::is(y, "SpatVector3D") && .is_grid_3d(x) &&
          .placed_on_grid(y, x)) {
      shape <- .vector3d_on_grid(y, x, ...)
      return(intersects_3d(x, shape))
    }

    shared <- intersect_3d(x, y, ...)
    grid <- if (.is_grid_3d(x)) x else y
    present_x <- .footprint(x, grid)
    present_y <- .footprint(y, grid)
    # The shared domain lies where both are present, so its footprint is the
    # answer there; FALSE where only one is, NA where neither is.
    out <- terra::ifel(present_x | present_y, .footprint(shared), NA)
    # ifel() yields 1/0/NA; only as.bool() makes values() give TRUE/FALSE/NA.
    out <- terra::as.bool(out)
    names(out) <- "intersects"
    methods::as(out, "SpatRaster")  # remove the 3D class
  }
)

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
#'   \item{voxel masked by [SpatVector3D-class] polygons}{The polygons are
#'     put on the voxel's grid and levels with [vector3d_to_voxel()], then
#'     used as a voxel mask as above: the voxel keeps its values in the cells
#'     and levels the polygons reach, and a depth gap between polygons over
#'     one cell stays empty. When the voxel has no depth edges, you set each
#'     level's slab with `bounds`, as in [envelope_to_voxel()]. The time
#'     intervals are dropped, with a message.}
#'   \item{envelope or plain `SpatRaster` masked by [SpatVector3D-class]
#'     polygons}{The polygons are put on the grid of `x` with
#'     [vector3d_to_envelope()], then used as an envelope mask as above. Where
#'     polygons overlap, a cell gets the shallowest top and the deepest
#'     bottom, so a depth gap between them is filled. The time intervals are
#'     dropped, with a message.}
#'   \item{anything masked by [SpatVector3D-class] points}{A point covers
#'     one cell and one depth level, as in [intersect_3d()]. A voxel keeps
#'     its value in each point's cell at one level, which `method` sets;
#'     `bounds` is an error with points. An envelope
#'     keeps its interval whole in the cells where a point lies within it. A
#'     plain `SpatRaster` keeps the cells that hold a point. Lines aren't
#'     supported; [buffer_3d()] turns them into polygons.}
#' }
#'
#' A voxel or envelope masked by a plain `SpatRaster` or by polygons is
#' terra's own `mask()`. Every layer is masked by the same 2D pattern and the
#' class is kept. No method is added for those cases.
#'
#' @param x The raster whose values are kept: a [SpatVoxel-class], a
#'   [SpatEnvelope-class], or a plain `SpatRaster`. A netCDF file isn't
#'   accepted: describe it with [as_proxy()] and read one time step with
#'   [as_voxel()].
#' @param mask The 3D domain to keep `x` inside: a [SpatEnvelope-class] or a
#'   [SpatVoxel-class] on the same grid as `x`, or [SpatVector3D-class]
#'   polygons or points, which are put on that grid.
#' @param bounds Only when `x` is a voxel and `mask` an envelope or
#'   [SpatVector3D-class] polygons. How the envelope or polygons are placed
#'   on the voxel's depth levels: `"midpoint"` (default) or `"top"`. See
#'   [envelope_to_voxel()]. Not used when the voxel carries its own depth
#'   edges (see [SpatVoxel-class]).
#' @param method Only when `x` is a voxel and `mask` [SpatVector3D-class]
#'   points: which depth level each point keeps.
#'   * `"nearest"` (default): the level nearest the point's depth, whatever
#'     the voxel holds there, as [intersect_3d()] places points. A point
#'     whose nearest level has no value keeps nothing.
#'   * `"nearest_value"`: the nearest level that has a value in the point's
#'     cell, as `extract_to_point(method = "nearest")` picks it. Levels
#'     without a value are skipped.
#'   * `"bottom"`: the deepest level that has a value in the point's cell,
#'     whatever the point's depth, as `extract_to_point(method = "bottom")`.
#'   * `"surface"`: the shallowest level, whatever the point's depth and
#'     whether or not it has a value, as
#'     `extract_to_point(method = "surface")`.
#'
#'   For the two nearest rules, a point halfway between two levels takes the
#'   shallower one. A point whose cell has no value at any level keeps
#'   nothing under `"nearest_value"` and `"bottom"`. These rules match
#'   [extract_to_point()]'s for the depth level. The cell is the one that
#'   holds the point; for a point exactly on the edge between two cells it
#'   can differ from the cell [extract_to_point()] picks.
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
  # remove the 3D class
  terra::mask(methods::as(x, "SpatRaster"), methods::as(keep, "SpatRaster"),
              ...)
}

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatVoxel", "SpatEnvelope"),
  function(x, mask, bounds = c("midpoint", "top"), ...) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")
    .check_3d(mask, "mask")
    .check_same_grid(mask, x, "`mask`")
    keep <- .promote_to_voxel(mask, depths(x), bounds, depth_bounds(x))
    as_voxel(.mask_plain(x, keep, ...), depth_bounds = depth_bounds(x))
  }
)

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatVoxel", "SpatVoxel"),
  function(x, mask, ...) {
    .check_3d(x, "x")
    .check_3d(mask, "mask")
    .check_same_grid(mask, x, "`mask`")
    # The result keeps the pair's shared edges, whichever voxel carries them,
    # as intersect_3d() does.
    shared <- .shared_depths(x, mask)
    keep <- .presence_voxel(.voxel_occupancy(mask), shared$depths)
    as_voxel(.mask_plain(x, keep, ...), depth_bounds = shared$depth_bounds)
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
  function(x, mask, bounds = c("midpoint", "top"), ...) {
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
    .check_same_grid(mask, x, "`mask`")
    terra::mask(x, terra::ifel(.footprint(mask), 1, NA), ...)
  }
)

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatRaster", "SpatVoxel"),
  function(x, mask, ...) {
    .check_3d(mask, "mask")
    .check_same_grid(mask, x, "`mask`")
    terra::mask(x, terra::ifel(.footprint(mask), 1, NA), ...)
  }
)

# ---- SpatVector3D polygons and points ---------------------------------------

# Internal: a SpatVector3D on the grid of `template`, so the queries above
# compare it with an envelope or voxel cell by cell, as SPEC.md ("Intersect:
# spatial query") plans. Lines are an error. Polygons on a voxel become a
# presence voxel on its levels, as vector3d_to_voxel() makes it, so a depth
# gap between two polygons over one cell stays empty; `...` is that
# function's `bounds`. Elsewhere polygons become an envelope as
# vector3d_to_envelope() makes it, and `...` is an error, as for two
# envelopes. No polygons, polygons that miss the grid and polygons with no
# thickness give an empty result instead of the converters' errors: sharing
# no space is an answer here, as it is for two envelopes, and intersect_3d()
# of two polygon sets can return no polygons.
#
# A point covers one cell and one depth level (#115). On a voxel, the points
# become a presence voxel on its levels: 1 at each point's cell and one
# level, NA elsewhere. `.level_rule` picks the level, through
# .summarise_depth_profile(), the rule extract_to_point() uses: "nearest"
# (the default, and the only rule intersect_3d() and intersects_3d() use)
# takes the level whose depth is nearest the point's depth, a tie going to
# the shallower level, whatever the voxel holds there. mask() can also pass
# "nearest_value", "bottom" or "surface", which read the cell's values. The
# argument has a dotted name so the generics' `...`, which reaches here,
# can't set it: `method` given to intersect_3d() stays an unused argument.
# On an envelope, the result
# is the envelope's part the points reach: in each cell holding a point whose
# depth lies within the cell's interval (both ends included), the shallowest
# and deepest such depth, so a lone point is a zero-thickness interval; NA
# elsewhere. Points off the grid aren't present. Callers assign the result
# before passing it on: an error raised while S4 evaluates an argument for
# dispatch loses its class.
.vector3d_on_grid <- function(x, template, ...,
                              .level_rule = "nearest") {
  .check_geomtype(x, c("polygons", "points"))
  # remove the 3D class
  grid <- methods::as(template[[1]], "SpatRaster")
  voxel <- methods::is(template, "SpatVoxel")
  if (x@type == "polygons" && voxel) {
    empty <- function(...) {
      occ <- terra::rast(grid, nlyrs = terra::nlyr(template), vals = 0)
      .presence_voxel(occ, depths(template), depth_bounds(template))
    }
    if (length(x) == 0L) return(empty())
    return(tryCatch(vector3d_to_voxel(x, template, ...),
                    ocean3d_error_envelope_empty = empty))
  }
  .no_extra_args(...)
  # An envelope with no cells, filled below where points reach it.
  out <- terra::rast(grid, nlyrs = 2, vals = NA_real_)
  names(out) <- c("depth_min", "depth_max")
  if (x@type == "polygons") {
    if (length(x) == 0L) return(as_envelope(out))
    return(tryCatch(vector3d_to_envelope(x, template),
                    ocean3d_error_envelope_empty = function(e) {
                      as_envelope(out)
                    }))
  }
  v <- x@vertices
  cells <- .point_cells(x, grid)
  on <- !is.na(cells)
  if (voxel) {
    d <- depths(template)
    kind <- c(nearest = "nearest", nearest_value = "nearest",
              bottom = "seabottom", surface = "surface")[[.level_rule]]
    # "nearest" reads no values, so one profile in which every level has a
    # value serves every point. The other rules read one column of values
    # per distinct cell. Each distinct profile and depth goes through
    # .summarise_depth_profile() once, so memory and time follow the cells
    # and depths, not the number of points. With uneven slabs (stored depth
    # edges), the slab that contains a point can belong to another level;
    # the nearest level is still used.
    level <- rep(NA_integer_, length(cells))
    if (any(on)) {
      if (.level_rule == "nearest") {
        profiles <- matrix(1, 1L, length(d))
        row <- rep(1L, sum(on))
      } else {
        cell_ids <- unique(cells[on])
        # remove the 3D class
        profiles <- as.matrix(terra::extract(
          methods::as(template, "SpatRaster"), cell_ids
        ))
        row <- match(cells[on], cell_ids)
      }
      target <- if (kind == "nearest") v$depth[on] else rep(NA_real_, sum(on))
      key <- paste(row, target)
      first <- !duplicated(key)
      picked <- vapply(which(first), function(k) {
        profile <- profiles[row[k], ]
        .summarise_depth_profile(profile, d, target[k])$index[[kind]]
      }, integer(1))
      level[on] <- picked[match(key, key[first])]
    }
    on <- on & !is.na(level)
    occ <- matrix(0, terra::ncell(grid), length(d))
    occ[cbind(cells[on], level[on])] <- 1
    occ <- terra::setValues(terra::rast(grid, nlyrs = length(d)), occ)
    return(.presence_voxel(occ, d, depth_bounds(template)))
  }
  limits <- terra::values(template, mat = TRUE)[cells[on], , drop = FALSE]
  depth <- v$depth[on]
  inside <- !is.na(limits[, "depth_min"]) & limits[, "depth_min"] <= depth &
    depth <= limits[, "depth_max"]
  if (!any(inside)) return(as_envelope(out))
  cell <- factor(cells[on][inside])
  reached <- data.frame(
    cell = as.integer(levels(cell)),
    depth_min = tapply(depth[inside], cell, min),
    depth_max = tapply(depth[inside], cell, max)
  )
  out[reached$cell] <- as.matrix(reached[c("depth_min", "depth_max")])
  as_envelope(out)
}

# Internal: whether a SpatVector3D goes on a grid object's grid before the
# two are compared, as .vector3d_on_grid() puts it. Polygons always do, and
# points on a voxel; points with an envelope are compared directly, since
# intersect_3d() gives the part of the envelope they reach. Shared by
# intersects_3d() and mask().
.placed_on_grid <- function(x, grid) {
  x@type == "polygons" || methods::is(grid, "SpatVoxel")
}

# The shape goes on the grid, with `...` (`bounds` for polygons on a voxel's
# levels, an error otherwise), then takes the methods for that pair. On an
# envelope, the part points reach is the answer itself (a second envelope
# intersection would drop a lone point's zero-thickness interval).
#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("SpatVector3D", "SpatEnvelope"),
          function(x, y, ...) {
            .check_3d(y, "y")
            on_grid <- .vector3d_on_grid(x, y, ...)
            if (x@type == "points") return(on_grid)
            intersect_3d(on_grid, y)
          })

#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("SpatVector3D", "SpatVoxel"),
          function(x, y, ...) {
            # Points take the level nearest their depth, so `bounds` has
            # nothing to place; .vector3d_on_grid() rejects it.
            on_grid <- .vector3d_on_grid(x, y, ...)
            intersect_3d(on_grid, y)
          })

#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("SpatEnvelope", "SpatVector3D"),
          function(x, y, ...) intersect_3d(y, x, ...))

#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("SpatVoxel", "SpatVector3D"),
          function(x, y, ...) intersect_3d(y, x, ...))

# Two sets of polygons: every pair whose outlines, depth ranges and times
# overlap gives one polygon, as sf::st_intersection() does in 2D. Points and
# polygons: every point inside a polygon's outline, depth range and time
# interval gives one point. Two sets of points share no space, and lines
# aren't taken.
#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatVector3D", "SpatVector3D"),
  function(x, y, ...) {
    .no_extra_args(...)
    .check_geomtype(x, c("polygons", "points"))
    .check_geomtype(y, c("polygons", "points"))
    if (x@type == "points" && y@type == "points") .stop_neither_3d(x, y)
    resolution <- unique(stats::na.omit(c(x@time_resolution,
                                          y@time_resolution)))
    if (length(resolution) > 1L) {
      rlang::abort(paste0(
        "One set has dates and the other has hours, so their ",
        "times can't be compared. Rebuild one with as_vector3d() from ",
        "times of the same kind."
      ), class = "ocean3d_error_time_resolution_differ", call = NULL)
    }
    resolution <- if (length(resolution)) resolution else NA_character_
    if (x@type == "points" || y@type == "points") {
      return(.points_in_polygons(x, y, resolution))
    }
    as_pieces <- function(v, suffix) {
      f <- v@features
      out <- sf::st_sf(f[c("feature_id", "depth_min", "depth_max",
                           "time_start", "time_end")],
                       geometry = f$geometry)
      names(out)[1:5] <- paste0(names(out)[1:5], suffix)
      # Depths and times hold for the whole polygon, which tells sf not to
      # warn that they might not.
      sf::st_agr(out) <- "constant"
      out
    }
    pairs <- sf::st_intersection(as_pieces(x, "_x"), as_pieces(y, "_y"))

    depth_min <- pmax(pairs$depth_min_x, pairs$depth_min_y)
    depth_max <- pmin(pairs$depth_max_x, pairs$depth_max_y)
    # A polygon without a time is there at any time, so the pair takes the
    # other polygon's interval. Intervals are closed: with dates, both
    # include their end day, so two that share a day overlap on it.
    time_start <- pmax(pairs$time_start_x, pairs$time_start_y, na.rm = TRUE)
    time_end <- pmin(pairs$time_end_x, pairs$time_end_y, na.rm = TRUE)
    # Outlines that only share an edge or corner give lines or points,
    # which have no area; a collection keeps its polygon parts.
    geometry <- lapply(sf::st_geometry(pairs), function(g) {
      if (inherits(g, "GEOMETRYCOLLECTION")) {
        # On one geometry, st_collection_extract() and st_union() give a
        # geometry, not a column, so the parts go in a column first.
        parts <- sf::st_sfc(sf::st_collection_extract(g, "POLYGON"))
        g <- sf::st_union(parts)[[1]]
      }
      g
    })
    geometry <- sf::st_sfc(geometry, crs = 4326)
    polygonal <- sf::st_geometry_type(geometry) %in%
      c("POLYGON", "MULTIPOLYGON") & !sf::st_is_empty(geometry)
    # Strict `>`, as for envelopes: depth ranges that only touch share no
    # water.
    keep <- polygonal & depth_max > depth_min &
      (is.na(time_start) | time_start <= time_end)

    features <- data.frame(
      feature_id = as.character(seq_len(sum(keep))),
      depth_min = depth_min[keep],
      depth_max = depth_max[keep],
      time_start = time_start[keep],
      time_end = time_end[keep],
      feature_id_x = pairs$feature_id_x[keep],
      feature_id_y = pairs$feature_id_y[keep],
      stringsAsFactors = FALSE
    )
    features$geometry <- geometry[keep]
    .new_vector3d("polygons", resolution, features = features)
  }
)

# Internal: the points of one SpatVector3D inside the polygons of the other,
# for intersect_3d(). A point is inside a polygon when its position is within
# the outline, its depth within the depth range and its time within the time
# interval, all with both ends included. With dates, a polygon's interval
# covers all of its end day. A point or polygon without a time matches any
# time. One point per matching pair, numbered, with `feature_id_x` and
# `feature_id_y` naming the pair as the polygon-pair result does; the inputs'
# other columns are dropped. `resolution` is the pair's shared time
# resolution, already checked.
.points_in_polygons <- function(x, y, resolution) {
  points_first <- x@type == "points"
  points <- if (points_first) x else y
  polygons <- if (points_first) y else x
  v <- points@vertices
  f <- polygons@features
  # Outlines and positions are both WGS84, so sf tests them with its own
  # model for geographic coordinates (the sphere, by default).
  positions <- sf::st_geometry(
    sf::st_as_sf(v[c("lon", "lat")], coords = c("lon", "lat"), crs = 4326)
  )
  hits <- sf::st_intersects(positions, f$geometry)
  i <- rep(seq_len(nrow(v)), lengths(hits))
  j <- unlist(hits)
  # Depths within the range, both ends included: a point at a polygon's
  # bottom is inside it.
  inside <- v$depth[i] >= f$depth_min[j] & v$depth[i] <= f$depth_max[j]
  # Times within the interval, both ends included, the rule the polygon pairs
  # use. With dates, every time is stored at midnight UTC and the two sides
  # share one resolution (checked by the caller), so `time <= time_end` also
  # takes a point on the interval's end day. A missing time on either side
  # matches any time.
  time <- v$time[i]
  in_time <- is.na(time) | is.na(f$time_start[j]) |
    (time >= f$time_start[j] & time <= f$time_end[j])
  keep <- inside & in_time
  i <- i[keep]
  j <- j[keep]
  vertices <- data.frame(
    feature_id = as.character(seq_along(i)),
    vertex = rep(1L, length(i)),
    lon = v$lon[i], lat = v$lat[i], depth = v$depth[i], time = v$time[i],
    feature_id_x = if (points_first) v$feature_id[i] else f$feature_id[j],
    feature_id_y = if (points_first) f$feature_id[j] else v$feature_id[i],
    stringsAsFactors = FALSE
  )
  # The points keep their times, so the result has their resolution.
  .new_vector3d("points", if (all(is.na(v$time[i]))) NA_character_ else
                  resolution, vertices = vertices)
}

# Two SpatVector3D sets: for each feature of `x`, the features of `y` it
# shares 3D space with. The pairs come from intersect_3d(), so the two always
# agree. For two polygon sets that builds every shared outline, which a
# yes/no answer doesn't need; testing outlines with a predicate that ignores
# shared edges, then depths and times, would be faster for large sets.
#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatVector3D", "SpatVector3D"),
  function(x, y, ...) {
    # Assigned first: an error raised while S4 evaluates an argument for
    # dispatch loses its class.
    shared <- intersect_3d(x, y, ...)
    shared <- as.data.frame(shared)
    ids_x <- as.data.frame(x)$feature_id
    ids_y <- as.data.frame(y)$feature_id
    hits <- split(shared$feature_id_y, factor(shared$feature_id_x, ids_x))
    # In the order of `y`, as sf::st_intersects() lists them.
    lapply(hits, function(h) ids_y[ids_y %in% h])
  }
)

# One signature serves voxels, envelopes and plain rasters. Polygons go on
# a voxel's levels with `bounds` and on other grids as an envelope, and
# points on a voxel as a presence voxel, so mask() then dispatches to the
# methods above for that class of `x`. An envelope keeps the cells where a
# point lies within its interval, and a plain raster the cells that hold a
# point. `...` is terra's.
#' @rdname mask-3d
#' @export
setMethod("mask", c("SpatRaster", "SpatVector3D"),
          function(x, mask, bounds = c("midpoint", "top"),
                   method = c("nearest", "nearest_value", "bottom",
                              "surface"),
                   ...) {
            # Checked first: below, an error raised while S4 evaluates
            # ifel()'s argument would lose its class.
            .check_geomtype(mask, c("polygons", "points"))
            # `method` applies only to points on a voxel, so it goes on only
            # when given; elsewhere it is the unused-argument error.
            rules <- c("nearest", "nearest_value", "bottom", "surface")
            if (missing(method)) {
              method <- NULL
            } else if (!rlang::is_string(method) || !method %in% rules) {
              rlang::abort(glue::glue(
                "`method` must be one of {toString(dQuote(rules, FALSE))}."
              ), class = "ocean3d_error_mask_method", call = NULL)
            }
            # `bounds` goes on only when given, so wherever there are no
            # levels to place it is the unused-argument error, as elsewhere.
            if (!missing(bounds) && !.placed_on_grid(mask, x)) {
              .no_extra_args(bounds = bounds)
            }
            if (!is.null(method) &&
                  !(mask@type == "points" && methods::is(x, "SpatVoxel"))) {
              .no_extra_args(method = method)
            }
            if (.placed_on_grid(mask, x)) {
              shape <- if (missing(bounds)) {
                .vector3d_on_grid(mask, x,
                                  .level_rule = method %||% "nearest")
              } else {
                .vector3d_on_grid(mask, x, bounds = bounds)
              }
              return(mask(x, shape, ...))
            }
            if (methods::is(x, "SpatEnvelope")) {
              .check_3d(x, "x")
              keep <- terra::ifel(intersects_3d(x, mask), 1, NA)
              return(as_envelope(.mask_plain(x, keep, ...)))
            }
            terra::mask(x, terra::ifel(.footprint(mask, x, "mask"), 1, NA),
                        ...)
          })

# ---- terra::intersect() guard -----------------------------------------------

#' @rdname intersect_3d
#' @section terra's `intersect()`:
#' `terra::intersect()` on two rasters is a 2D test that ignores depth. On a
#' [SpatEnvelope-class] or a [SpatVoxel-class], in either position, it is an
#' error that points to [intersects_3d()] and `intersect_3d()`. Plain
#' rasters, vectors and extents keep terra's behaviour.
#' @name intersect-guard
#' @aliases intersect,SpatEnvelope,SpatEnvelope-method
#' @aliases intersect,SpatEnvelope,SpatVoxel-method
#' @aliases intersect,SpatVoxel,SpatEnvelope-method
#' @aliases intersect,SpatVoxel,SpatVoxel-method
#' @aliases intersect,SpatEnvelope,SpatRaster-method
#' @aliases intersect,SpatVoxel,SpatRaster-method
#' @aliases intersect,SpatRaster,SpatEnvelope-method
#' @aliases intersect,SpatRaster,SpatVoxel-method
#' @exportMethod intersect
NULL

# terra's intersect() for two rasters is a 2D "both have data" test, and terra
# propagates the subclass, so intersect(envelope_a, envelope_b) came back as a
# SpatEnvelope holding TRUE/FALSE that ignored depth and still passed
# validity. On a 3D object that answer is misleading under any class tag, so
# it is an error that names the two replacements. Plain rasters, vectors and
# extents keep terra's methods: none of the signatures below match them, and
# terra's generic has no `...`, so these methods take exactly (x, y).
for (sig in list(c("SpatEnvelope", "SpatEnvelope"),
                 c("SpatEnvelope", "SpatVoxel"),
                 c("SpatVoxel", "SpatEnvelope"),
                 c("SpatVoxel", "SpatVoxel"),
                 c("SpatEnvelope", "SpatRaster"),
                 c("SpatVoxel", "SpatRaster"),
                 c("SpatRaster", "SpatEnvelope"),
                 c("SpatRaster", "SpatVoxel"))) {
  setMethod("intersect", sig, function(x, y) {
    stop("intersect() is terra's 2D test and ignores the depth axis. For a ",
         "SpatEnvelope or SpatVoxel, use intersects_3d() to ask whether two ",
         "objects share 3D space, or intersect_3d() for the space they share.",
         call. = FALSE)
  })
}
rm(sig)
