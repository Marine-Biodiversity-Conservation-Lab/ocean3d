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
# method takes: neither side is a 3D object, or SpatVector3D polygons meet
# something other than an envelope, a voxel or other polygons.
.stop_neither_3d <- function(x, y) {
  vector <- methods::is(x, "SpatVector3D") || methods::is(y, "SpatVector3D")
  message <- if (vector) {
    other <- if (methods::is(x, "SpatVector3D")) y else x
    glue::glue(
      "SpatVector3D polygons can be compared with a SpatEnvelope, a ",
      "SpatVoxel or other SpatVector3D polygons, not with ",
      "{class(other)[1L]}. Build 3D polygons from 2D ones with as_vector3d()."
    )
  } else {
    paste0(
      "neither `x` nor `y` is a SpatEnvelope, a SpatVoxel or SpatVector3D ",
      "polygons. A 3D query needs at least one 3D object; build one with ",
      "as_envelope(), vect_to_envelope(), as_voxel() or as_vector3d()."
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
    intersect_3d(.promote_to_voxel(x, y, bounds), y, ...)
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

# Intersection is symmetric, so a 2D object in first position is the same
# query the other way round. Any other pair has no 3D object to answer on.
#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("ANY", "ANY"), function(x, y, ...) {
  if (methods::is(y, "SpatEnvelope") || methods::is(y, "SpatVoxel")) {
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
    intersects_3d(.promote_to_voxel(x, y, bounds), y, ...)
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
    intersects_3d(x, .promote_to_voxel(y, x, bounds), ...)
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
    is_grid <- function(v) {
      methods::is(v, "SpatEnvelope") || methods::is(v, "SpatVoxel")
    }
    # Polygons go on the other object's grid as an envelope, then take the
    # envelope methods, including the promotion onto a voxel's levels.
    # Assigned first: an error raised while S4 evaluates an argument for
    # dispatch loses its class.
    if (methods::is(x, "SpatVector3D") && is_grid(y)) {
      polygons <- .polygons_on_grid(x, y)
      return(intersects_3d(polygons, y, ...))
    }
    if (methods::is(y, "SpatVector3D") && is_grid(x)) {
      polygons <- .polygons_on_grid(y, x)
      return(intersects_3d(x, polygons, ...))
    }

    shared <- intersect_3d(x, y, ...)
    grid <- if (is_grid(x)) x else y
    present_x <- .footprint(x, grid, "x")
    present_y <- .footprint(y, grid, "y")
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
#'   \item{anything masked by [SpatVector3D-class] polygons}{The polygons
#'     are put on the grid of `x` with [vector3d_to_envelope()], then used as
#'     an envelope mask as above. A voxel keeps its values in the cells and
#'     levels the polygons reach. Where polygons overlap, a cell gets the
#'     shallowest top and the deepest bottom, so a depth gap between them is
#'     kept too. The time intervals are dropped, with a message.}
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
#'   polygons, which are put on that grid.
#' @param bounds Only when one side is an envelope and the other a voxel. How
#'   the envelope is placed on the voxel's depth levels: `"midpoint"`
#'   (default) or `"top"`. See [envelope_to_voxel()]. Not used when the voxel
#'   carries its own depth edges (see [SpatVoxel-class]).
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
    .check_same_grid(x, mask, c("x", "mask"))
    keep <- .promote_to_voxel(mask, x, bounds)
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
    .check_same_grid(x, mask, c("x", "mask"))
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
    .check_same_grid(x, mask, c("x", "mask"))
    terra::mask(x, terra::ifel(.footprint(mask), 1, NA), ...)
  }
)

#' @rdname mask-3d
#' @export
setMethod(
  "mask", c("SpatRaster", "SpatVoxel"),
  function(x, mask, ...) {
    .check_3d(mask, "mask")
    .check_same_grid(x, mask, c("x", "mask"))
    terra::mask(x, terra::ifel(.footprint(mask), 1, NA), ...)
  }
)

# ---- SpatVector3D polygons --------------------------------------------------

# Internal: SpatVector3D polygons on the grid of `template` as an envelope,
# so the queries above compare them with an envelope or voxel cell by cell,
# as SPEC.md ("Intersect: spatial query") plans. No polygons, polygons that
# miss the grid and polygons with no thickness give an empty envelope instead
# of vector3d_to_envelope()'s errors: sharing no space is an answer here, as
# it is for two envelopes, and intersect_3d() of two polygon sets can return
# no polygons. Callers assign its result before passing it on: an error
# raised while S4 evaluates an argument for dispatch loses its class.
.polygons_on_grid <- function(x, template) {
  # Checked first, so empty points or lines are still an error.
  .check_polygons(x)
  empty <- function(...) {
    # remove the 3D class
    grid <- methods::as(template[[1]], "SpatRaster")
    out <- terra::rast(grid, nlyrs = 2, vals = NA_real_)
    names(out) <- c("depth_min", "depth_max")
    as_envelope(out)
  }
  if (length(x) == 0L) return(empty())
  tryCatch(vector3d_to_envelope(x, template),
           ocean3d_error_envelope_empty = empty)
}

#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("SpatVector3D", "SpatEnvelope"),
          function(x, y, ...) {
            polygons <- .polygons_on_grid(x, y)
            intersect_3d(polygons, y, ...)
          })

#' @rdname intersect_3d
#' @export
setMethod("intersect_3d", c("SpatVector3D", "SpatVoxel"),
          function(x, y, ...) {
            polygons <- .polygons_on_grid(x, y)
            intersect_3d(polygons, y, ...)
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
# overlap gives one polygon, as sf::st_intersection() does in 2D.
#' @rdname intersect_3d
#' @export
setMethod(
  "intersect_3d", c("SpatVector3D", "SpatVector3D"),
  function(x, y, ...) {
    .no_extra_args(...)
    .check_polygons(x)
    .check_polygons(y)
    resolution <- unique(stats::na.omit(c(x@time_resolution,
                                          y@time_resolution)))
    if (length(resolution) > 1L) {
      rlang::abort(paste0(
        "One set of polygons has dates and the other has hours, so their ",
        "times can't be compared. Rebuild one with as_vector3d() from ",
        "times of the same kind."
      ), class = "ocean3d_error_time_resolution_differ", call = NULL)
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
    .new_vector3d("polygons",
                  if (length(resolution)) resolution else NA_character_,
                  features = features)
  }
)

# Two sets of polygons: for each polygon of `x`, the polygons of `y` it
# shares 3D space with. The pairs come from intersect_3d(), so the two always
# agree. That builds every shared outline, which a yes/no answer doesn't
# need; testing outlines with a predicate that ignores shared edges, then
# depths and times, would be faster for large sets.
#' @rdname intersects_3d
#' @export
setMethod(
  "intersects_3d", c("SpatVector3D", "SpatVector3D"),
  function(x, y, ...) {
    shared <- as.data.frame(intersect_3d(x, y, ...))
    ids_x <- as.data.frame(x)$feature_id
    ids_y <- as.data.frame(y)$feature_id
    hits <- split(shared$feature_id_y, factor(shared$feature_id_x, ids_x))
    # In the order of `y`, as sf::st_intersects() lists them.
    lapply(hits, function(h) ids_y[ids_y %in% h])
  }
)

# One signature serves voxels, envelopes and plain rasters: the polygons go
# on the grid of `x`, and mask() then dispatches to the envelope methods
# above for that class of `x`.
#' @rdname mask-3d
#' @export
setMethod("mask", c("SpatRaster", "SpatVector3D"),
          function(x, mask, ...) {
            polygons <- .polygons_on_grid(mask, x)
            mask(x, polygons, ...)
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
