# Constructors 

#' Coerce a 2D footprint to a 2.5D min-max envelope
#'
#' Build a [SpatEnvelope-class] — one continuous `[depth_min, depth_max]`
#' interval per grid cell — from a 2D horizontal footprint and the depth limits
#' that apply to it. This is the general form of the conversion
#' [vect_to_envelope()] performs for polygons: any 2D raster whose non-`NA`
#' cells mark presence becomes a 3D domain once depth limits are attached to
#' it.
#'
#' `x` is the footprint. Its values are not read, only their non-`NA` pattern:
#' a rasterized species range, a Global Fishing Watch effort layer, or a plain
#' presence mask all work. A SpatRaster that *already* has exactly the two
#' layers `depth_min` and `depth_max` (the output of [vect_to_envelope()], or a
#' [SpatEnvelope-class] itself) is promoted directly instead, in which case
#' `depth_min` and `depth_max` must be omitted.
#'
#' Depths are positive metres increasing downward, following the package's
#' depth sign convention. Flip GEBCO-style elevation first, e.g.
#' `terra::clamp(-terra::project(bathy, x), lower = 0)`, which gives a positive
#' seafloor depth with land clamped to 0.
#'
#' @param x SpatRaster. Either a single-layer footprint whose non-`NA` cells are
#'   present, or a two-layer raster already named `depth_min`, `depth_max`.
#' @param depth_min,depth_max Shallowest and deepest depth in metres, positive
#'   down. Each is either a single number applying to the whole footprint, or a
#'   single-layer SpatRaster on the grid of `x` giving the limit per cell. Omit
#'   both when `x` already carries the two depth layers.
#'
#' @returns A [SpatEnvelope-class] with layers `depth_min` and `depth_max`, on
#'   the grid of `x`.
#'
#' @examples
#' # A 2x2 footprint: three cells present (any non-NA value), one absent.
#' fp <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2,
#'                   ymin = 0, ymax = 2)
#' terra::values(fp) <- c(1, 1, 1, NA)
#'
#' # A species recorded between 0 and 200 m.
#' e <- as_envelope(fp, depth_min = 0, depth_max = 200)
#' terra::values(e)
#'
#' # Per-cell limits are allowed too, as single-layer rasters.
#' dmax <- terra::setValues(terra::rast(fp), c(100, 200, 300, 400))
#' terra::values(as_envelope(fp, depth_min = 10, depth_max = dmax))
#' @export
as_envelope <- function(x, depth_min, depth_max) {
  if (inherits(x, "sf") || inherits(x, "sfc") || inherits(x, "SpatVector")) {
    stop("`x` must be a SpatRaster footprint, not vector geometry. ",
         "Rasterize the polygons onto the study grid first, or use ",
         "vect_to_envelope(), which does both.", call. = FALSE)
  }
  if (!inherits(x, "SpatRaster")) {
    stop("`x` must be a SpatRaster.", call. = FALSE)
  }
  # Work on a plain raster whatever class tag `x` arrived with: a footprint
  # taken from a voxel is still tagged SpatVoxel, and a tagged intermediate
  # would dispatch to the depth-aware mask() methods in R/intersect.R.
  x <- methods::as(x, "SpatRaster")  # remove the 3D class

  # A raster already carrying the two depth layers is promoted, not rebuilt.
  already_envelope <- identical(names(x), c("depth_min", "depth_max"))
  supplied <- !missing(depth_min) || !missing(depth_max)

  if (already_envelope) {
    if (supplied) {
      stop("`x` already has depth_min and depth_max layers; omit the ",
           "`depth_min` and `depth_max` arguments.", call. = FALSE)
    }
    dmin <- x[["depth_min"]]
    dmax <- x[["depth_max"]]
  } else {
    if (missing(depth_min) || missing(depth_max)) {
      stop("`depth_min` and `depth_max` are required unless `x` already has ",
           "layers named depth_min and depth_max. Got layers: ",
           toString(names(x)), call. = FALSE)
    }
    if (terra::nlyr(x) != 1) {
      stop("`x` must be a single-layer footprint when depth limits are ",
           "supplied; got ", terra::nlyr(x), " layers.", call. = FALSE)
    }

    # Only the non-NA pattern of the footprint matters, never its values.
    present <- terra::ifel(is.na(x), NA, 1)
    dmin <- terra::mask(.as_depth_layer(depth_min, "depth_min", x), present)
    dmax <- terra::mask(.as_depth_layer(depth_max, "depth_max", x), present)
  }

  out <- c(dmin, dmax)
  names(out) <- c("depth_min", "depth_max")

  # new() runs the validity method, which checks the depth rules.
  methods::new("SpatEnvelope", out)
}

# Internal: resolve a depth limit given either as a single number or as a
# single-layer SpatRaster into a named layer on the grid of `template`.
.as_depth_layer <- function(value, name, template) {
  if (inherits(value, "SpatRaster")) {
    if (terra::nlyr(value) != 1) {
      stop("`", name, "` must be a single-layer SpatRaster; got ",
           terra::nlyr(value), " layers.", call. = FALSE)
    }
    if (!terra::compareGeom(value, template, stopOnError = FALSE)) {
      stop("`", name, "` must be on the same grid (CRS, extent, resolution) ",
           "as `x`.", call. = FALSE)
    }
    out <- value
  } else {
    if (!is.numeric(value) || length(value) != 1 || is.na(value)) {
      stop("`", name, "` must be a single non-missing number or a ",
           "single-layer SpatRaster.", call. = FALSE)
    }
    out <- terra::setValues(terra::rast(template), value)
  }
  names(out) <- name
  out
}

#' Create a voxel object
#'
#' Wrap a multi-depth raster as a [SpatVoxel-class]: the validated 3D form used
#' throughout the package, with one layer per standard depth. Each layer's
#' depth is stored in [terra::depth()]. Given a [SpatProxy-class],
#' `as_voxel()` reads one time step of it into memory instead.
#'
#' A voxel read from a proxy carries the proxy's depth edges (see
#' [depth_bounds()]): the file's depth bounds, or, for a file without them
#' (such as Copernicus Marine), edges halfway between levels, from the
#' surface to half a gap below the deepest level. Because the voxel has
#' edges, the `bounds` argument of [volume()] and the other 3D functions
#' doesn't apply to it, and each layer stands for the same slab of water
#' whichever layers are selected. To use other edges, pass them with
#' `as_voxel(v, depth_bounds = ...)`. A single level without bounds has no
#' edges, so `bounds` applies.
#'
#' The depths come from `depths` when it is given, and otherwise from the
#' raster's own [terra::depth()], which terra also sets when it reads a netCDF
#' file. A raster with neither is an error: layer names are labels, never
#' read as depths.
#'
#' Most terra operations (`crop()`, `mask()`, `[[`, arithmetic, `c()`) keep
#' the class and the depths; a few, such as `mean()`, return a plain
#' `SpatRaster` instead. Re-wrapping is cheap either way: given a valid
#' `SpatVoxel` and no `depths`, `as_voxel()` returns it untouched.
#'
#' Keeping the class doesn't re-run the class checks, so an operation that
#' changes the layers can leave an object labelled `SpatVoxel` that isn't a
#' valid one (for example `v[[c(1, 1)]]`, two layers at the same depth).
#' `as_voxel()` checks such an object as it would a plain raster: unsorted
#' layers are sorted, and the rest is an error. A voxel whose depths terra
#' dropped, for example in `app()`, needs `depths` again (see
#' [SpatVoxel-class]).
#'
#' Depths are positive metres increasing downward, matching the World Ocean
#' Atlas convention. Negative depths are an error rather than being silently
#' negated, since flipping the sign would change what the data mean.
#'
#' @param x A [SpatProxy-class], a SpatRaster, a list of SpatRasters (one per
#'   depth, combined with [terra::rast()], which keeps their depths only when
#'   each has them set; otherwise pass `depths`), or an existing
#'   [SpatVoxel-class]. This is how gridded data
#'   enter the `SpatVoxel` model: describe a netCDF file with [as_proxy()] and
#'   pass the proxy here, which reads one time step. netCDF paths aren't
#'   accepted directly.
#' @param depths Optional numeric vector, one depth per layer, in metres. When
#'   supplied, it sets the depths, and the layers are named
#'   `{varname}_depth={value}`, replacing any existing names. Required if `x`
#'   has no depths in [terra::depth()]. Leave it out, and `depth_bounds` too,
#'   for a [SpatProxy-class], which gives its own.
#' @param varname Character. Variable name used in the layer names built from
#'   `depths`, or from a proxy's depths, where it defaults to the proxy's
#'   variable name. Ignored otherwise.
#' @param time For a [SpatProxy-class] only: the time step to read, as a step
#'   number or a date or date-time matched exactly at its own resolution (so
#'   `"2020-01-02"` selects that day). `NULL` (the default) reads the only
#'   step of a proxy that has one; a proxy with more needs `time`.
#' @param depth_bounds Optional two-column matrix or data frame: the top and
#'   bottom of the slab of water each layer stands for, in metres positive
#'   down, one row per layer of `x` in its current order (for example a netCDF
#'   file's depth bounds). A data frame with columns `depth`, `lower` and
#'   `upper`, such as [depth_bounds()] of another voxel, is matched to the
#'   layers by depth instead. Each slab must contain its depth and have a
#'   positive height, and slabs must not overlap. `NULL` (the default) keeps
#'   the edges an existing `SpatVoxel` carries, unless `depths` renames its
#'   layers; a voxel without edges uses the `bounds` argument of the functions
#'   it is passed to.
#'
#' @returns A [SpatVoxel-class] whose layers are ordered shallow to deep.
#' @examples
#' r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
#' terra::values(r) <- runif(terra::ncell(r) * 3)
#'
#' # set the depths from a vector; the layers get readable names
#' v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")
#' depths(v)
#' names(v)
#'
#' # depths already in terra::depth() are used as they stand, and sorted if
#' # needed
#' terra::depth(r) <- c(200, 0, 100)
#' depths(as_voxel(r))
#'
#' # layer names aren't read as depths: a raster without terra::depth() needs
#' # `depths`
#' named <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2, vals = 1)
#' names(named) <- c("temp_depth=0", "temp_depth=50")
#' try(as_voxel(named))
#'
#' # carry each layer's slab of water, as a netCDF file's depth bounds give it
#' v <- as_voxel(r, depth_bounds = cbind(c(150, 0, 50), c(250, 50, 150)))
#' depth_bounds(v)
#' @export
as_voxel <- function(x, depths = NULL, varname = "value",
                     depth_bounds = NULL, time = NULL) {
  # A SpatProxy is read one time step at a time; it gives its own depths and
  # edges, and `time` applies only to it.
  proxy <- methods::is(x, "SpatProxy")
  misplaced <- if (proxy) {
    if (!is.null(depths) || !is.null(depth_bounds)) {
      paste("A SpatProxy gives its own depths and edges; omit `depths`",
            "and `depth_bounds`.")
    }
  } else if (!is.null(time)) {
    "`time` selects a time step of a SpatProxy; `x` isn't one."
  }
  if (!is.null(misplaced)) {
    rlang::abort(misplaced, class = "ocean3d_error_voxel_input", call = NULL)
  }
  if (proxy) {
    if (missing(varname)) varname <- x@variable$name
    return(.proxy_to_voxel(x, time, varname))
  }
  # Idempotent: a valid SpatVoxel comes back as it is, so re-wrapping is cheap
  # and callers can do it defensively. terra propagates the class tag through
  # `[[`, `names<-` and the like without re-running validity, so a tagged
  # object may no longer be a voxel. `as_voxel()` repairs what is unambiguous 
  # (layer order) and rejects the rest with the constructor's own messages. 
  edges <- NULL
  was_voxel <- methods::is(x, "SpatVoxel")
  if (was_voxel) {
    if (is.null(depths) && is.null(depth_bounds) &&
          isTRUE(methods::validObject(x, test = TRUE))) {
      return(x)
    }
    # Edges are keyed by depth, so they still apply unless the layers are
    # given other depths.
    if (is.null(depths)) edges <- x@depth_edges
    x <- methods::as(x, "SpatRaster")  # remove the 3D class
  }

  if (is.list(x)) x <- terra::rast(x)
  if (!methods::is(x, "SpatRaster")) {
    stop("`x` must be a SpatRaster or a list of SpatRasters; got ",
         paste(class(x), collapse = "/"), ".", call. = FALSE)
  }

  # The depths come from `depths`, else from terra::depth(); layer names are
  # labels and never read. terra::depth() is set only on a copy of `x`,
  # because terra sets it in place and `x` may still be the caller's raster.
  if (!is.null(depths)) {
    problem <- if (length(depths) != terra::nlyr(x)) {
      glue::glue("must have one value per layer: got {length(depths)} for ",
                 "{terra::nlyr(x)} layers")
    } else if (!is.numeric(depths)) {
      "must be numeric metres"
    } else if (anyNA(depths)) {
      "must give every layer a depth; some are NA"
    }
    if (!is.null(problem)) {
      rlang::abort(glue::glue("`depths` {problem}."),
                   class = "ocean3d_error_voxel_input", call = NULL)
    }
    d <- as.numeric(depths)
    # Readable names; `names<-` copies the raster first.
    names(x) <- paste0(varname, "_depth=", d)
    terra::depth(x) <- d
  } else if (.depth_set(x)) {
    d <- terra::depth(x)
  } else {
    .abort_no_depths()
  }
  if (any(d < 0)) {
    stop("depths are positive metres increasing downward; got ",
         toString(d[d < 0]),
         ". Negate the depths in `depths` or terra::depth().",
         call. = FALSE)
  }
  if (anyDuplicated(d)) {
    stop("duplicate depth(s): ", toString(unique(d[duplicated(d)])),
         ". Each layer must be a distinct depth.", call. = FALSE)
  }

  if (!is.null(depth_bounds)) {
    edges <- .depth_edge_table(d, depth_bounds)
  } else if (!is.null(edges) && nrow(edges)) {
    # Edges the voxel had, for the layers it still has.
    row <- .edge_rows(edges$depth, d)
    edges <- edges[row[!is.na(row)], , drop = FALSE]
  }

  # Repair what is unambiguously repairable rather than rejecting it. The
  # subset is a new raster, and a layer subset clears terra's depth flag, so
  # the depths are set again.
  if (is.unsorted(d)) {
    x <- x[[order(d)]]
    d <- sort(d)
    terra::depth(x) <- d
  }

  # Without edges, the class prototype gives the empty table.
  if (is.null(edges)) return(methods::new("SpatVoxel", x))
  rownames(edges) <- NULL
  methods::new("SpatVoxel", x, depth_edges = edges)
}

#' @rdname SpatVoxel-class
#' @param x A `SpatVoxel`.
#' @param ... For `c()`: further rasters to combine with `x`, as with terra's
#'   `c()`. The result is a `SpatVoxel` with each layer's depth and the depth
#'   edges of every `SpatVoxel` among them, or a plain `SpatRaster` when the
#'   layers don't make a valid voxel (for example two layers at the same
#'   depth). For `subset()`: passed to terra's `subset()`.
#' @export
setMethod("c", "SpatVoxel", function(x, ...) {
  args <- list(x, ...)
  is_raster <- vapply(args, methods::is, logical(1), "SpatRaster")
  # Each raster's depths, in order, read before combining.
  d <- unlist(lapply(args[is_raster], function(a) {
    tryCatch(depths(a), error = function(e) rep(NA_real_, terra::nlyr(a)))
  }))
  # terra edits its inputs' depths in place when it combines them, so it gets
  # plain copies, and the result is tagged as a voxel again afterwards.
  args[is_raster] <- lapply(args[is_raster], function(a) {
    terra::deepcopy(methods::as(a, "SpatRaster"))
  })
  out <- do.call(c, args)
  if (anyNA(d) || length(d) != terra::nlyr(out)) return(out)
  terra::depth(out) <- d
  # Only each voxel's rows for the layers it still has, so a row left over
  # from a dropped layer can't stand in for another voxel's edges.
  voxels <- Filter(function(a) methods::is(a, "SpatVoxel"), list(x, ...))
  edges <- do.call(rbind, lapply(voxels, function(v) {
    row <- .edge_rows(v@depth_edges$depth, depths(v))
    v@depth_edges[row[!is.na(row)], , drop = FALSE]
  }))
  edges <- edges[!duplicated(edges$depth), , drop = FALSE]
  rownames(edges) <- NULL
  # new() runs the validity checks; layers that aren't a valid voxel stay a
  # plain raster.
  tryCatch(methods::new("SpatVoxel", out, depth_edges = edges),
           error = function(e) out)
})

#' @rdname SpatVoxel-class
#' @param subset,negate,NSE For `subset()`, which `[[` and `$` also call: the
#'   layers to keep, as with terra's `subset()`. The result is a `SpatVoxel`
#'   that keeps each layer's depth through later operations.
#' @export
# `NSE` keeps terra's argument name.
setMethod("subset", "SpatVoxel", function(x, subset, negate = FALSE,
                                          NSE = FALSE, ...) { # nolint
  # terra's method reads a non-standard expression with substitute(), which
  # would see only this method's argument, so the layers are found here, as
  # terra does, and passed on as indices.
  if (isTRUE(NSE)) {
    layers <- stats::setNames(as.list(seq_len(terra::nlyr(x))), names(x))
    subset <- if (missing(subset)) unlist(layers) else
      eval(substitute(subset), layers, parent.frame())
  }
  out <- methods::callNextMethod(x, subset = subset, negate = negate,
                                 NSE = FALSE, ...)
  # A layer subset of a raster in memory keeps the depths but clears terra's
  # depth flag, so the next operation (arithmetic, crop(), mask(), ...) would
  # reset them to 0. The subset is a new raster, so setting its depths again
  # doesn't touch `x`.
  if (.depth_set(x) && !.depth_set(out)) terra::depth(out) <- terra::depth(out)
  out
})

#' @rdname SpatVoxel-class
#' @param i,j,value For `[[<-`: the layer to replace, by number or name, and
#'   its new values, as with terra. A replaced layer keeps its depth.
#' @export
setReplaceMethod("[[", c("SpatVoxel", "numeric"), function(x, i, j, value) {
  .replace_layer(x, methods::callNextMethod())
})

#' @rdname SpatVoxel-class
#' @export
setReplaceMethod("[[", c("SpatVoxel", "character"), function(x, i, j, value) {
  .replace_layer(x, methods::callNextMethod())
})

# Internal: the result `out` of terra's `[[<-` on voxel `x`, with x's depths
# and class back. terra gives a replaced layer the depth of the new values,
# which is 0 when they have none, so each replaced layer keeps its own depth
# instead. Replacing by name returns a plain SpatRaster, so the result goes
# back into `x`, whose only other slot is depth_edges; like terra's own
# operations, this doesn't re-run the class checks. A layer added by name has
# no depth to keep, so that result is left as terra returns it.
.replace_layer <- function(x, out) {
  d <- if (.depth_set(x)) terra::depth(x)
  if (length(d) != terra::nlyr(out)) return(out)
  terra::depth(out) <- d
  x@pntr <- out@pntr
  x
}

# Internal: a depth_edges table (depth, lower, upper) from `depth_bounds`: a
# two-column matrix or data frame of top and bottom with one row per depth, in
# the order of `depths`, or a table with columns depth, lower and upper (such
# as depth_bounds() of a voxel), matched to `depths` by depth. The slabs are
# checked as every depth axis is, by .axis_intervals(); problems stop with
# `ocean3d_error_bounds`.
.depth_edge_table <- function(depths, depth_bounds) {
  if (is.data.frame(depth_bounds) &&
        all(c("depth", "lower", "upper") %in% names(depth_bounds))) {
    row <- .edge_rows(depth_bounds$depth, depths)
    if (anyNA(row)) {
      .abort_bounds("depth_bounds", NULL, glue::glue(
        "have no row for depth {toString(depths[is.na(row)])}"
      ))
    }
    m <- cbind(depth_bounds$lower[row], depth_bounds$upper[row])
  } else {
    m <- if (is.data.frame(depth_bounds)) as.matrix(depth_bounds) else
      depth_bounds
  }
  if (!is.matrix(m) || !is.numeric(m) || ncol(m) != 2L ||
        nrow(m) != length(depths)) {
    shape <- if (is.matrix(m)) paste(dim(m), collapse = " x ") else class(m)[1]
    .abort_bounds("depth_bounds", NULL, glue::glue(
      "must be a numeric two-column matrix or data frame (top, bottom) with ",
      "one row per depth: got {shape} for {length(depths)} depths"
    ))
  }
  o <- order(depths)
  .axis_intervals(depths[o], t(m[o, , drop = FALSE]), "depth", "depth_bounds")
  lower <- pmin(m[, 1], m[, 2])
  if (any(lower < 0)) {
    .abort_bounds("depth_bounds", NULL,
                  "must be positive metres increasing downward")
  }
  data.frame(depth = depths, lower = lower, upper = pmax(m[, 1], m[, 2]))
}

#' Convert Envelope 2.5D -> Voxel 3D
#' 
#' Expand a [SpatEnvelope-class] to the [SpatVoxel-class], with an input of
#' depth levels. A depth level belongs to a cell when the slab of water it
#' stands for shares some thickness with that cell's `[depth_min, depth_max]`
#' interval; depths outside the interval, and cells that are `NA` in the
#' envelope, are `NA` in every layer.
#'
#' \itemize{
#'   \item `values` is the magnitude, which can be a single number for the
#'     whole grid or a SpatRaster carrying a value per cell. Default value is
#'     `1`, which makes the output a presence mask. Ex. `values` can be 2D
#'     fisheries effort rasters.
#'   \item `profile` decides how that magnitude is spread down the water
#'     column. Default `NULL` means the full value is written at every
#'     occupied level. Providing a `profile_*` function determines how
#'     `values` is transformed for each depth level. Ex. [profile_equal()]
#'     divides it evenly over the levels the cell occupies, so sum of all
#'     depth levels in the voxel equals `values`.
#' }
#'
#' A depth level stands for a slab of water rather than a knife-edge, and a cell
#' occupies that level when its envelope overlaps the slab by more than a
#' shared edge. An interval of `[10, 20]` against levels `c(0, 100)` therefore
#' occupies the 0 m level, rather than falling between the levels and coming
#' back empty. Touching is not overlapping: an interval that ends exactly where
#' a slab begins, `[0, 100]` against the slab `100-200`, does not occupy that
#' level, so adjacent depth ranges are placed on disjoint levels and do not
#' intersect as voxels. This is the same rule [intersect_3d()] applies between
#' two envelopes.
#'
#' Where a slab's edges fall is a property of the dataset the levels came from,
#' not of the levels themselves, so `bounds` selects the convention:
#'
#' \itemize{
#'   \item `"midpoint"` (the default) puts the edges halfway between
#'     neighbouring levels. This is the World Ocean Atlas convention:
#'     `0, 5, 10, ...` the 0 m layer covers `0-2.5` m, the 5 m layer
#'     `2.5-7.5` m, the 10 m layer `7.5-12.5` m, and so on. It is close for
#'     model output such as Copernicus Marine, but leaves out the water above
#'     the shallowest level and below the deepest.
#'   \item `"top"` runs each slab from its level down to the next, so a level
#'     names the top of what it stands for. Levels `c(0, 100, 200)` give slabs
#'     `0-100`, `100-200`.
#' }
#'
#' When the dataset gives the real edges (a netCDF file's depth bounds), pass
#' them as `depth_bounds` instead.
#'
#' Under either convention the slabs cover the span from the shallowest level to
#' the deepest and no further: outer edges are clamped to `range(depths)` rather
#' than extrapolated, which is why WOA's 0 m layer starts at the surface rather
#' than half a gap above it. The expansion is therefore limited only at the two
#' ends when a cell whose envelope lies wholly above `min(depths)` or wholly
#' below `max(depths)`. Under `"top"` the deepest level has no next level to
#' run to, so it stands for a zero-thickness slab: an envelope is recorded
#' there only when it reaches below `max(depths)`, and one that ends exactly at
#' `max(depths)` stops at the level above. For the same reason a
#' [voxel_to_envelope()] round trip drops each cell's deepest occupied level:
#' the envelope's `depth_max` names that level, which is the top of its slab.
#'
#' @param x SpatEnvelope, e.g. from [as_envelope()] or [voxel_to_envelope()].
#' @param depths Array of values that can be coerced into numeric type,
#'   represents metres depth below sea level. Sorted shallow to deep, and
#'   deduplicated, before use.
#' @param values Optional numeric or SpatRaster of same CRS, resolution, extent
#'   as `x`. Defines values to write for voxel cells that are within depth
#'   intervals.
#'   `NULL` (the default) writes `1` at every occupied level, giving simple 
#'   presence/absence.
#' @param profile Optional function distributing `values` across the depth
#'   dimension, taking `(ind, depths, n_depths)` and returning the weight to
#'   apply. `NULL` (the default) writes the same value at every depth, as
#'   [profile_flat()] does; [profile_equal()] divides the value equally across
#'   the depths a cell occupies. Any function meeting the contract in
#'   [voxel_profiles] works, including your own.
#' @param bounds How to place the edges of the slab each depth level stands for.
#'   `"midpoint"` (the default) puts the edges halfway between neighbouring
#'   levels; `"top"` runs the slab from the level down to the next, so the
#'   level is its shallow edge. Under both, the outer edges are clamped to
#'   `range(depths)` rather than extrapolated. Not used when `depth_bounds`
#'   is given.
#' @param varname Variable name for the layer names, which take the form
#'   `{varname}_depth={depths[i]}`. The depths themselves are stored in
#'   [terra::depth()].
#' @param depth_bounds Optional two-column matrix or data frame: the top and
#'   bottom of the slab of water each depth stands for, in metres, one row per
#'   element of `depths` (or a `depth`, `lower`, `upper` table, as in
#'   [as_voxel()]). When given, these slabs are used instead of `bounds`, the
#'   returned voxel carries them, and `depths` must not repeat a depth.
#'
#' @returns A [SpatVoxel-class] with one layer per depth in `depths`, on the
#'   grid of `x`, layers ordered shallow to deep.
#'
#' @seealso [voxel_to_envelope()], the reverse (and lossy) collapse.
#'
#' @examples
#' fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2,
#'                   ymin = 0, ymax = 1)
#' terra::values(fp) <- c(1, NA)
#' e <- as_envelope(fp, depth_min = 50, depth_max = 200)
#'
#' # presence at every standard depth the envelope covers
#' terra::values(envelope_to_voxel(e, depths = c(0, 50, 100, 200, 300)))
#'
#' # the same pattern, carrying a constant of your choosing instead of 1
#' terra::values(
#'   envelope_to_voxel(e, depths = c(0, 50, 100, 200, 300),
#'                     values = 5)
#' )
#'
#' # spread that value down the column instead: an even share at each of the
#' # depths the cell occupies, summing back to 1
#' terra::values(
#'   envelope_to_voxel(e, depths = c(0, 50, 100, 200, 300),
#'                     profile = profile_equal,
#'                     varname = "time")
#' )
#'
#' # where a slab's edges fall decides the level: [60, 70] sits in the 50-150 m
#' # slab of the 100 m level by default (WOA-style), but in the 0-100 m slab of
#' # the 0 m level under "top"
#' shallow <- as_envelope(fp, depth_min = 60, depth_max = 70)
#' terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200)))
#' terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200),
#'                                 bounds = "top"))
#'
#' # or give the dataset's own edges
#' terra::values(envelope_to_voxel(
#'   shallow, depths = c(0, 100, 200),
#'   depth_bounds = cbind(c(0, 65, 150), c(65, 150, 250))
#' ))
#'
#' # a per-cell magnitude comes in as a raster: split each cell's effort evenly
#' # over the depths it occupies, conserving the cell total
#' effort <- terra::setValues(terra::rast(fp), c(100, NA))
#' terra::values(
#'   envelope_to_voxel(
#'     e, depths = c(0, 50, 100, 200, 300), values = effort,
#'     profile = profile_equal,
#'     varname = "effort"
#'   )
#' )
#' @export
envelope_to_voxel <- function(x, depths, values = NULL, profile = NULL,
                              bounds = c("midpoint", "top"),
                              varname = "presence", depth_bounds = NULL) {
  if (!is(x, "SpatEnvelope")) {
    stop("Input error for envelope_to_voxel(): `x` needs to be of ",
         "`SpatEnvelope` class.", call. = FALSE)
  }
  # Coercion, not class: an integer vector, or a character vector of numbers,
  # is as good as a double here — "hello" is not a depth.
  depths <- suppressWarnings(as.numeric(depths))
  if (length(depths) == 0 || anyNA(depths)) {
    stop("Input error for envelope_to_voxel(): `depths` needs to be array ",
         "coercible to numeric type.", call. = FALSE)
  }
  # match.arg() is the single source of truth for the accepted conventions:
  # they are named once, in the signature, where .depth_slabs() also reads them.
  bounds <- match.arg(bounds)
  if (!is.character(varname) || length(varname) != 1 || is.na(varname)) {
    stop("Input error for envelope_to_voxel(): `varname` needs to be a string.",
         call. = FALSE)
  }
  if (any(depths < 0)) {
    stop("Input error for envelope_to_voxel(): depths are positive metres ",
         "increasing downward; got ", toString(depths[depths < 0]),
         ".", call. = FALSE)
  }
  # Checked up front rather than on first use, so a bad profile is an error
  # even on an envelope where no cell turns out to be occupied. What the
  # function returns is checked later, once there is an `ind` to size it
  # against.
  if (!is.null(profile) && !is.function(profile)) {
    stop("Input error for envelope_to_voxel(): `profile` needs to be NULL or ",
         "a function of (ind, depths, n_depths); got ",
         paste(class(profile), collapse = "/"), ". See ?voxel_profiles.",
         call. = FALSE)
  }

  # A voxel's layer axis is one layer per depth, shallow to deep, so the
  # requested levels are put in that form up front. Given edges come along;
  # with edges, a repeated depth is an error from the edge checks.
  if (!is.null(depth_bounds)) {
    depth_bounds <- .depth_edge_table(depths, depth_bounds)
    depth_bounds <- depth_bounds[order(depth_bounds$depth), , drop = FALSE]
  }
  depths <- sort(unique(depths))

  # Step 1: occupancy. Each level stands for a slab of water, `bounds` (or
  # the given edges) says where its edges fall, and a cell occupies the level
  # when its envelope shares some thickness of water with that slab. Strict
  # `>` on the envelope's floor: an envelope that merely ends where a slab
  # begins shares no water with it, the same rule intersect_3d() applies
  # between envelopes, so adjacent depth ranges never meet on a level. Slabs
  # are half-open so a shared edge is claimed by the deeper level alone. The
  # deepest slab worked out from `depths` is closed, since under "top" it has
  # no height; given slabs have a height, so they all stay half-open.
  # Raster stays left: terra mishandles `scalar <op> raster`.
  slab <- .depth_slabs(depths, bounds, depth_bounds)
  closed <- is.null(depth_bounds)
  # Plain layers, so the occupancy stack built from them is a plain SpatRaster.
  lower <- .envelope_layer(x, "depth_min")
  upper <- .envelope_layer(x, "depth_max")
  deepest <- length(depths)
  ind <- terra::rast(lapply(seq_along(depths), function(i) {
    if (closed && i == deepest) {
      (lower <= slab$upper[i]) & (upper > slab$lower[i])
    } else {
      (lower < slab$upper[i]) & (upper > slab$lower[i])
    }
  }))
  # number of depth levels present per vertical column, cell 
  n_depths <- sum(ind)

  # Slabs worked out from `depths` tile `[min(depths), max(depths)]`, so a
  # cell is left empty only when its envelope lies wholly above the shallowest
  # level or wholly below the deepest: a range the levels cannot reach. Given
  # slabs can leave gaps, so there an envelope can also fall between them.
  # `n_depths` is NA where the envelope is, so this counts only real cells.
  n_missed <- terra::global(n_depths == 0, "sum", na.rm = TRUE)[1, 1]
  if (!is.na(n_missed) && n_missed > 0) {
    cells <- if (n_missed == 1) "cell has" else "cells have"
    where <- if (closed) {
      paste("the slabs of `depths`: it lies wholly above or below them.",
            "These cells are empty in the voxel. Extend `depths` to span them.")
    } else {
      paste("the slabs in `depth_bounds`: it lies above, below or between",
            "them. These cells are empty in the voxel.")
    }
    rlang::warn(
      glue::glue("{n_missed} {cells} an envelope that shares no water with ",
                 "{where}"),
      class = "ocean3d_warning_envelope_outside_slabs"
    )
  }

  # Step 2: values. The magnitude each cell carries, and how `profile` spreads
  # it down the levels that cell occupies.
  v <- .cell_values(values, x)

  # NULL means the flat shape, which is a profile like any other. Called
  # positionally, so a caller's own profile is free to name its arguments
  # however it likes.
  if (is.null(profile)) profile <- profile_flat
  w <- .profile_weights(profile, ind, depths, n_depths)

  # Raster stays on the left throughout: terra broadcasts the single-layer `v`
  # and `w` across every depth layer of `ind`.
  out <- ind * v * w
  # Multiplication left out-of-band cells at 0 rather than NA.
  out <- terra::mask(out, ind, maskvalues = 0)
  as_voxel(out, depths = depths, varname = varname,
           depth_bounds = depth_bounds)
}

# Internal: the slab of water each depth level stands for, as parallel vectors
# of its shallow and deep edge. Slabs are half-open, `[lower, upper)`, except
# the deepest, which is closed so an envelope that starts exactly on the last
# level and reaches below it is still caught.
#
# Outer edges are clamped to `range(depths)` rather than extrapolated, under
# both conventions. Extrapolating would invent reach the levels were never
# given: at the top it can run above the surface (WOA's shallowest level is
# 0 m, and half a gap above that is a negative depth), and at the bottom it
# would quietly swallow every range below the grid, turning a range the levels
# cannot express into a plausible-looking answer instead of the warning the
# caller needs to see.
#
# `depth_bounds`, a voxel's stored edges (a table with `lower` and `upper` in
# the order of `depths`, as depth_bounds() gives it), wins over `bounds`
# when given.
.depth_slabs <- function(depths, bounds, depth_bounds = NULL) {
  if (!is.null(depth_bounds)) {
    return(list(lower = depth_bounds$lower, upper = depth_bounds$upper))
  }
  n <- length(depths)
  switch(
    bounds,
    # The level names the top of what it stands for, so the slab runs down to
    # the next level. The deepest has no next and collapses to a point.
    top = list(lower = depths, upper = c(depths[-1], depths[n])),
    # The level sits in the middle of what it stands for, so the edges fall
    # halfway to each neighbour. This is the World Ocean Atlas convention.
    midpoint = {
      mids <- (depths[-n] + depths[-1]) / 2
      list(lower = c(depths[1], mids), upper = c(mids, depths[n]))
    },
    stop("Unhandled `bounds`: \"", bounds, "\".", call. = FALSE)
  ) # nolint: unreachable_code_linter. switch() default, not dead code.
}

# Internal: resolve `values` into the magnitude each cell carries — a single
# number applying to the whole grid, or a single-layer raster on its grid.
.cell_values <- function(values, template) {
  if (is.null(values)) return(1)

  if (inherits(values, "SpatRaster")) {
    if (terra::nlyr(values) != 1) {
      stop("`values` must be a single-layer SpatRaster; got ",
           terra::nlyr(values), " layers.", call. = FALSE)
    }
    if (!terra::compareGeom(values, template, stopOnError = FALSE)) {
      stop("`values` must be on the same grid (CRS, extent, resolution) as ",
           "`x`.", call. = FALSE)
    }
    return(values)
  }

  if (!is.numeric(values) || length(values) != 1 || is.na(values)) {
    stop("`values` must be a single non-missing number or a single-layer ",
         "SpatRaster.", call. = FALSE)
  }
  values
}

#' Collapse Voxel 3D -> Envelope 2.5D
#'
#' Reduce a [SpatVoxel-class] to the [SpatEnvelope-class] that bounds it. A
#' predicate `fun` is applied to the cell values at each depth; for every cell,
#' the shallowest depth at which the predicate is `TRUE` becomes `depth_min`
#' and the deepest becomes `depth_max`. Cells where the predicate is never
#' `TRUE` are `NA` in both layers.
#'
#' The default predicate, `\(x) !is.na(x)`, gives the plain vertical extent of
#' the data: the shallowest and deepest depths at which the cell has any value.
#' Pass a different predicate to bound a subset of the values instead, e.g.
#' `\(x) x > 15` for the depths over which a cell exceeds 15 degrees.
#'
#' When the voxel carries depth edges (see [SpatVoxel-class]), each cell's
#' envelope runs from the top of the slab of its shallowest selected layer
#' (where `fun` is `TRUE`) to the bottom of its deepest one's, so it covers
#' the water those layers stand for. A voxel
#' without edges gives the depths of the layers themselves.
#'
#' **This conversion is lossy.** An envelope stores a single continuous
#' interval per cell, so any interior gap in the voxel is filled in: a cell that
#' satisfies `fun` at 0 m and 200 m but not at 100 m still yields the envelope
#' `[0, 200]`. Use [SpatVoxel-class] directly where interior gaps matter.
#'
#' @param v A [SpatVoxel-class]. A netCDF file isn't accepted: describe it
#'   with [as_proxy()] and read one time step with [as_voxel()].
#' @param fun Function taking a vector of cell values and returning a logical
#'   vector of the same length. `NA` results are treated as `FALSE`. Defaults to
#'   `\(x) !is.na(x)`.
#'
#' @returns A [SpatEnvelope-class] with layers `depth_min` and `depth_max`, on
#'   the same grid as `v`.
#'
#' @seealso [envelope_to_voxel()], the reverse expansion.
#'
#' @examples
#' # Two cells sampled at four standard depths.
#' r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2,
#'                  ymin = 0, ymax = 1, nlyrs = 4)
#' terra::values(r) <- cbind(c(12, NA), c(11, NA), c(NA, 8), c(6, NA))
#' v <- as_voxel(r, depths = c(0, 50, 100, 200), varname = "temp")
#'
#' # Default predicate: the vertical extent of the data. Cell 1 has values at
#' # 0, 50 and 200 m, so it comes back as [0, 200] — the gap at 100 m is
#' # filled, because an envelope stores one continuous interval per cell.
#' terra::values(voxel_to_envelope(v))
#'
#' # Another predicate bounds a subset of the values instead: here the depths
#' # over which a cell is warmer than 10 degrees. Cell 2 never qualifies, so it
#' # is NA in both layers.
#' terra::values(voxel_to_envelope(v, fun = function(x) x > 10))
#'
#' # With depth edges, the envelope covers the slabs: cell 1 runs from the top
#' # of the 0 m slab to the bottom of the 200 m slab.
#' edged <- as_voxel(v, depth_bounds = cbind(c(0, 25, 75, 150),
#'                                           c(25, 75, 150, 250)))
#' terra::values(voxel_to_envelope(edged))
#' @export
voxel_to_envelope <- function(v, fun = function(x) !is.na(x)) {
  if (!is.function(fun)) {
    stop("`fun` must be a function returning TRUE/FALSE for a cell value. ",
         "The former \"extent\" behaviour is the default, ",
         "function(x) !is.na(x).",
         call. = FALSE)
  }
  if (!is(v, "SpatVoxel")) {
    stop("`v` must be a SpatVoxel object.")
  }
  .check_3d(v, "v")
  depths <- depths(v)
  n_cells <- terra::ncell(v)
  n_depths <- length(depths)

  # Hit matrix (cells x depths): does `fun` hold for this cell at this depth?
  # The stack is read once and `fun` is still applied one depth layer at a
  # time, so a predicate that reduces over a layer (e.g. `\(x) x > mean(x)`)
  # keeps its per-depth meaning. `fun` is evaluated on the values rather than
  # on the SpatRaster so that any ordinary R predicate works, not only
  # terra-aware ones.
  vals <- terra::values(v)
  hit <- matrix(FALSE, nrow = n_cells, ncol = n_depths)
  for (i in seq_len(n_depths)) {
    h <- as.logical(fun(vals[, i]))
    if (length(h) != n_cells) {
      stop("`fun` must return one TRUE/FALSE per cell value; got ",
           length(h), " values for ", n_cells, " cells.",
           call. = FALSE)
    }
    h[is.na(h)] <- FALSE
    hit[, i] <- h
  }
  rm(vals)

  # `SpatVoxel` validity guarantees depths are non-NA and sorted shallow to
  # deep, so each row's first and last TRUE column are exactly the envelope
  # bounds. `max.col()` finds them at C level, which avoids building one
  # full-size intermediate raster per depth.
  first_hit <- max.col(hit, ties.method = "first")
  last_hit <- max.col(hit, ties.method = "last")

  # An all-FALSE row still yields a column index, so cells where the predicate
  # never holds must go back to NA rather than point at the shallowest depth.
  never <- !hit[cbind(seq_len(n_cells), first_hit)]
  # With stored edges, the top of the shallowest slab and the bottom of the
  # deepest; without, the layer depths.
  edges <- depth_bounds(v)
  depth_min <- if (is.null(edges)) depths[first_hit] else
    edges$lower[first_hit]
  depth_max <- if (is.null(edges)) depths[last_hit] else edges$upper[last_hit]
  depth_min[never] <- NA_real_
  depth_max[never] <- NA_real_

  # remove the 3D class
  out <- terra::setValues(
    terra::rast(methods::as(v[[1]], "SpatRaster"), nlyrs = 2),
    cbind(depth_min, depth_max)
  )
  names(out) <- c("depth_min", "depth_max")

  # new() runs the validity method.
  methods::new("SpatEnvelope", out)
}

#' Reduce a variable voxel to a presence voxel
#'
#' A [SpatVoxel-class] holds a variable, such as temperature, oxygen or
#' fishing effort. It does not say which cells count as occupied.
#' `occupied()` applies a predicate and returns presence: 1 where `fun` is
#' `TRUE`, `NA` elsewhere. The grid and the depth levels stay the same.
#'
#' The 3D query verbs read presence directly and take no predicate. Threshold
#' each voxel first, so both sides of a query use their own cutoff:
#'
#' ```
#' intersect_3d(occupied(temp, function(v) v > 15),
#'              occupied(effort, function(v) v > 0))
#' ```
#'
#' `fun` receives one depth layer's cell values as a numeric vector. It must
#' return a logical vector of the same length. A predicate that summarises the
#' layer, such as `function(v) v > mean(v, na.rm = TRUE)`, therefore compares
#' each cell against its own depth. [voxel_to_envelope()] uses the same
#' contract.
#'
#' Running `occupied()` on a presence voxel returns it unchanged.
#'
#' @param x A [SpatVoxel-class]. A netCDF file isn't accepted: describe it
#'   with [as_proxy()] and read one time step with [as_voxel()].
#' @param fun Function taking a vector of cell values and returning a logical
#'   vector of the same length. `NA` counts as `FALSE`. Defaults to
#'   `function(v) !is.na(v)`, which treats any recorded value as occupied.
#'
#' @returns A [SpatVoxel-class] on the depth levels of `x`. Values are 1
#'   (occupied) or `NA`. Layers are named `presence_depth={value}`.
#'
#' @seealso [voxel_to_envelope()] collapses the result to a
#'   [SpatEnvelope-class]. [volume()], [intersect_3d()], [intersects_3d()] and
#'   [mask()] take a presence voxel.
#'
#' @examples
#' r <- terra::rast(nrows = 1, ncols = 3, nlyrs = 2)
#' terra::values(r) <- cbind(c(12, 18, 4), c(9, 16, 2))
#' v <- as_voxel(r, depths = c(0, 100), varname = "temp")
#'
#' # Occupied where the water is warmer than 10 degrees.
#' warm <- occupied(v, function(x) x > 10)
#' names(warm)
#' terra::values(warm)
#'
#' # The verbs need no predicate.
#' volume(warm)
#'
#' # A presence voxel is returned unchanged.
#' identical(terra::values(occupied(warm)), terra::values(warm))
#'
#' # Each side of a query uses its own cutoff.
#' effort <- as_voxel(terra::setValues(r, cbind(c(0, 3, 0), c(0, 2, 0))),
#'                    depths = c(0, 100), varname = "effort")
#' terra::values(intersect_3d(warm, occupied(effort, function(x) x > 0)))
#' @export
occupied <- function(x, fun = function(v) !is.na(v)) {
  if (!methods::is(x, "SpatVoxel")) {
    stop("`x` must be a SpatVoxel; got ", paste(class(x), collapse = "/"),
         ". Build one with as_voxel() or envelope_to_voxel().", call. = FALSE)
  }
  .check_3d(x, "x")
  .presence_voxel(.voxel_occupancy(x, fun), depths(x),
                  depth_bounds(x))
}

# Internal: validate and normalise the `depth_min` / `depth_max` inputs of
# `vect_to_envelope()`. Accepts a list, a numeric vector, or a SpatRaster
# Returns a list whose elements can be used in `do.call("max", ...)` /
# `do.call("min", ...)`.
.normalize_depth_inputs <- function(x, template, arg) {
  if (inherits(x, "SpatRaster")) {
    x <- list(x)
  } else if (!is.list(x)) {
    x <- as.list(x)
  }

  if (length(x) == 0L) {
    .abort_envelope_input(glue::glue(
      "`{arg}` is empty. Supply at least one numeric value or SpatRaster."
    ))
  }

  for (i in seq_along(x)) {
    el <- x[[i]]
    label <- if (length(x) == 1L) {
      glue::glue("`{arg}`")
    } else {
      glue::glue("`{arg}` element {i}")
    }

    if (inherits(el, "SpatRaster")) {
      if (!terra::same.crs(el, template)) {
        .abort_envelope_input(glue::glue(
          "{label} has a different CRS than `template`. ",
          "Project it onto the template grid first."
        ))
      }
      aligned <- terra::compareGeom(el, template, crs = FALSE,
                                    stopOnError = FALSE, messages = FALSE)
      if (!isTRUE(aligned)) {
        .abort_envelope_input(glue::glue(
          "{label} does not align with `template`; its extent or ",
          "resolution differs. Resample it onto the template grid first."
        ))
      }
    } else if (!is.numeric(el)) {
      .abort_envelope_input(glue::glue(
        "{label} must be numeric or a SpatRaster, not {class(el)[1]}."
      ))
    } else if (anyNA(el)) {
      # A scalar NA carries no information but, under `na.rm = FALSE`, would
      # blank every cell of the output. Reject it here where we can say which
      # argument was at fault instead of reporting an empty envelope later.
      .abort_envelope_input(glue::glue(
        "{label} is NA. Drop it, or supply a real depth constraint."
      ))
    }
  }

  x
}

# Internal: collapse the validated `depth_min` / `depth_max` constraints of
# `vect_to_envelope()` into one SpatRaster on the `template` grid, applying
# `fun` across them cell by cell ("max" for depth_min, "min" for depth_max, so
# that every constraint narrows the envelope).
.combine_depths <- function(x, template, fun) {
  # identify which inputs are SpatRaster, which are not 
  is_rast <- vapply(x, inherits, logical(1), what = "SpatRaster")

  # skips the reordered do.call when there are no rasters in input list x
  # creates raster with single value across
  if (!any(is_rast)) {
    return(terra::setValues(terra::rast(template[[1]]), do.call(fun, x)))
  }
  # terra has min() and max() functions that can apply to a list that includes
  # both SpatRaster and numeric values. However, it only works when the
  # SpatRaster is first in the list; the base::min is called instead if the
  # first element is a numeric object. This do.call line takes the reordered
  # input list, with the SpatRaster objects first.

  # `na.rm = FALSE` so that cell where constraints / depth raster
  # has NA value keeps it in the output raster. Avoids filling the value from 
  # remaining constraints
  do.call(fun, c(x[is_rast], x[!is_rast], list(na.rm = FALSE)))
}

# The one producer of `ocean3d_error_envelope_input`: input that can't be
# put on a grid as an envelope, from vect_to_envelope() or
# vector3d_to_envelope().
.abort_envelope_input <- function(message) {
  rlang::abort(message, class = "ocean3d_error_envelope_input", call = NULL)
}

# The one producer of `ocean3d_error_envelope_empty`: an envelope with no
# cells left, because the shapes miss the grid or have no thickness.
.abort_envelope_empty <- function(message) {
  rlang::abort(message, class = "ocean3d_error_envelope_empty", call = NULL)
}

#' Convert SpatVector or sf to SpatEnvelope
#' 
#' Rasterize polygons onto a template raster grid and assign per-cell depth
#' limits, outputting an SpatEnvelope object. The depths are clamped by
#' depth_min (shallowest) and depth_max (deepest), and/or by rasters (ex.
#' bathymetry) that are the same coordinate and resolution as the template.
#' Where a constraint raster is NA at a cell its limit there is unknown, so the
#' cell is dropped rather than falling back to the remaining constraints.
#' 
#' @param polygon sf or SpatVector, for example species ranges or fishery
#'   footprints.
#' @param template SpatRaster that defines the horizontal grid of the
#'   SpatEnvelope output.
#' @param depth_min List. Can contain both numeric and SpatRasters that match
#'   CRS, resolution, extent of `template` and contain numeric values. For each
#'   cell, the maximum across the `depth_min` list parameters is used as output
#'   SpatEnvelope depth_min layer cell value.
#' @param depth_max List. Can contain both numeric and SpatRasters that match
#'   CRS, resolution, extent of `template` and contain numeric values. For each
#'   cell, the minimum across the `depth_max` list parameters is used as output
#'   SpatEnvelope depth_max layer cell value. Cells where it isn't deeper
#'   than `depth_min` are `NA`.
#' 
#' @returns SpatEnvelope, with depth_min and depth_max layers. Stops with an
#'   error if the envelope would have no cells: when the polygons don't
#'   overlap `template`, or no cell they cover has `depth_max` deeper than
#'   `depth_min`.
#'
#' @seealso [as_envelope()], which attaches depth limits to a footprint that is
#'   already a raster.
#'
#' @examples
#' # A 4x4 study grid, and a range polygon covering its middle. `polygon` and
#' # `template` must already share a CRS: this function does not reproject.
#' template <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
#'                         ymin = 0, ymax = 4, crs = "EPSG:4326")
#' poly <- terra::vect("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))", crs = "EPSG:4326")
#'
#' # A species recorded between 50 m and 200 m throughout its range. Cells
#' # outside the polygon are NA.
#' e <- vect_to_envelope(poly, template, depth_min = 50, depth_max = 200)
#' terra::values(e)
#'
#' # Bathymetry is not a special argument: pass the seafloor as a second
#' # `depth_max` and, per cell, the shallower of the two wins. This seafloor
#' # deepens from north to south, so the two rows the polygon covers meet the
#' # bed very differently.
#' seafloor <- terra::setValues(terra::rast(template),
#'                              rep(c(20, 45, 250, 400), each = 4))
#' e2 <- vect_to_envelope(poly, template,
#'                        depth_min = 50, depth_max = list(300, seafloor))
#' terra::values(e2)
#'
#' # Two things happened. Where the bed is at 250 m, the nominal 300 m limit is
#' # truncated to it. Where the bed is at 45 m it sits above `depth_min`, so
#' # there is no water column left to occupy and the cell drops out entirely —
#' # leaving 2 of the polygon's 4 cells.
#' sum(!is.na(terra::values(e2[["depth_min"]])))
#' @export
vect_to_envelope <- function(polygon, template, depth_min, depth_max) {
  # check params
  if (!inherits(polygon, c("sf", "sfc", "SpatVector"))) {
    .abort_envelope_input("`polygon` must be a SpatVector, sf or sfc object.")
  }
  if (inherits(polygon, "sf") || inherits(polygon, "sfc")) {
    polygon <- terra::vect(polygon)
  }
  if (!is(template, "SpatRaster")) {
    .abort_envelope_input("`template` must be a SpatRaster.")
  }
  # Only the grid is used. Drop any class tag so a voxel can serve as the
  # template without its layers reaching the depth-aware mask() methods.
  template <- methods::as(template[[1]], "SpatRaster")  # remove the 3D class
  if (!terra::same.crs(crs(polygon), crs(template))) {
    .abort_envelope_input(paste0(
      "`polygon` and `template` have different CRS. Project `polygon` onto ",
      "the CRS of `template` first."
    ))
  }

  # Depths stay 8-byte doubles if terra writes to disk; its default 4-byte
  # floats would turn 10.1 m into 10.1000004 m.
  datatype <- terra::terraOptions(print = FALSE)$datatype
  terra::terraOptions(datatype = "FLT8S")
  on.exit(terra::terraOptions(datatype = datatype), add = TRUE)

  # check the depth inputs for CRS, resolution, extent matching and reject NA
  # params
  depth_min <- .normalize_depth_inputs(depth_min, template, "depth_min")
  depth_max <- .normalize_depth_inputs(depth_max, template, "depth_max")

  # take first layer in case multi-layer raster
  all_one <- template[[1]]
  # make all cell values 1, since mask leaves NAs alone 
  terra::values(all_one) <- 1
  # intersect between polygon and template
  masked <- terra::mask(all_one, polygon) 
  if (terra::global(masked, "notNA")[1, 1] == 0) {
    .abort_envelope_empty(paste0(
      "The polygons don't overlap any cell of `template`, so the envelope ",
      "would be empty. Check that they cover the same area."
    ))
  }

  # Creating depth_min, depth_max rasters that form the SpatEnvelope output
  # - value that narrows the envelope the most is selected
  # - preserves NA values from the constraints, rather than fall back on other
  #   constraints
  # deepest of depth_min constraints
  depth_min_rast <- .combine_depths(depth_min, template, "max") %>%
    mask(masked)
  # shallowest of depth_max constraints
  depth_max_rast <- .combine_depths(depth_max, template, "min") %>%
    mask(masked)

  # Keep cells where depth_max is deeper than depth_min. A cell with no
  # thickness has no volume, and one with depth_max shallower is invalid.
  diffs <- depth_max_rast - depth_min_rast
  msk <- ifel(diffs > 0, 1, NA)
  if (terra::global(msk, "notNA")[1, 1] == 0) {
    # Say how the covered cells fail, since they can fail in different ways.
    count <- function(r) {
      n <- terra::global(r, "sum", na.rm = TRUE)[1, 1]
      if (is.na(n)) 0 else n
    }
    counts <- c(
      shallower = count(diffs < 0),
      equal = count(diffs == 0),
      missing = count(!is.na(masked) & is.na(diffs))
    )
    verb <- ifelse(counts == 1, "has", "have")
    parts <- c(
      if (counts[["shallower"]] > 0) glue::glue(
        "{counts[['shallower']]} {verb[['shallower']]} `depth_max` shallower ",
        "than `depth_min` (the limits may be swapped, or a seafloor raster is ",
        "above `depth_min`)"
      ),
      if (counts[["equal"]] > 0) {
        glue::glue(
          "{counts[['equal']]} {verb[['equal']]} them equal, so no thickness"
        )
      },
      if (counts[["missing"]] > 0) {
        glue::glue(
          "{counts[['missing']]} {verb[['missing']]} no `depth_min` or ",
          "`depth_max`"
        )
      }
    )
    parts <- glue::glue_collapse(parts, sep = ", ", last = " and ")
    .abort_envelope_empty(glue::glue(
      "No cell has `depth_max` deeper than `depth_min`, so the envelope ",
      "would be empty. Of the {sum(counts)} cells the polygons cover, ",
      "{parts}."
    ))
  }

  # mask cells that are not valid, depth_min is deeper than depth_max
  depth_min_rast <- mask(depth_min_rast, msk)
  depth_max_rast <- mask(depth_max_rast, msk)

  # assign SpatEnvelope class to output
  out <- rast(c(depth_min = depth_min_rast, depth_max = depth_max_rast)) %>%
    as_envelope()
  out
}

#' Create a study area raster grid
#'
#' Build an empty raster covering the combined extent of one or more spatial
#' objects. Useful for defining the common grid before rasterizing species
#' ranges and fishery footprints.
#'
#' @param layers List of sf, sfc, SpatVector, or SpatRaster objects. The output
#'   extent will cover all objects.
#' @param res Numeric vector of length 1 or 2. Cell resolution in units of
#'   `crs` (degrees for lon/lat). Default `0.01` (~1 km at equator).
#' @param crs Character. Coordinate reference system. Default `"EPSG:4326"`.
#'
#' @returns An empty SpatRaster with the computed extent, resolution, and CRS.
#'
#' @examples
#' # Two features in different places: the grid has to cover both.
#' a <- terra::vect("POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))", crs = "EPSG:4326")
#' b <- terra::vect("POLYGON ((3 3, 5 3, 5 5, 3 5, 3 3))", crs = "EPSG:4326")
#'
#' grid <- create_study_raster(list(a, b), res = 0.5)
#' terra::ext(grid)
#' dim(grid)              # rows, columns, layers
#' terra::hasValues(grid) # FALSE: the grid is empty, it only defines geometry
#'
#' # Inputs are projected to `crs` before their extents are combined, so a list
#' # mixing coordinate systems is fine.
#' create_study_raster(list(a, terra::project(b, "EPSG:3857")), res = 0.5)
#' @export
create_study_raster <- function(layers, res = 0.01, crs = "EPSG:4326") {
  extents <- lapply(layers, function(x) {
    if (inherits(x, "sf") || inherits(x, "sfc")) {
      x <- terra::vect(x)
    }
    if (inherits(x, "SpatVector")) {
      x <- terra::project(x, crs)
    } else if (inherits(x, "SpatRaster")) {
      x <- terra::project(x, crs)
    }
    terra::ext(x)
  })

  combined <- extents[[1]]
  for (i in seq_along(extents)[-1]) {
    combined <- terra::union(combined, extents[[i]])
  }

  terra::rast(combined, res = res, crs = crs)
}
