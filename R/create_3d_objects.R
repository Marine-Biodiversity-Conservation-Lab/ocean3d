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
#' fp <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
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
  x <- .as_plain_raster(x)

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

  # Check dmin and dmax layers 
  min_of <- function(r) {
    v <- terra::global(r, "min", na.rm = TRUE)[1, 1]
    if (is.null(v) || is.na(v)) NA_real_ else v
  }
  shallowest <- min_of(dmin)
  if (!is.na(shallowest) && shallowest < 0) {
    stop("depths are positive metres increasing downward; `depth_min` has ",
         "negative values. Flip the sign of elevation data first.",
         call. = FALSE)
  }
  thinnest <- min_of(dmax - dmin)
  if (!is.na(thinnest) && thinnest < 0) {
    stop("`depth_max` must be at least `depth_min` in every cell.",
         call. = FALSE)
  }

  out <- c(dmin, dmax)
  names(out) <- c("depth_min", "depth_max")

  out <- methods::new("SpatEnvelope", out)
  methods::validObject(out)
  out
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
#' throughout the package, with one layer per standard depth and layer names
#' following the `{variable}_depth={value}` convention.
#'
#' Most terra operations (`crop()`, `mask()`, `[[`, arithmetic) propagate the
#' class, so a `SpatVoxel` normally survives them; a few, such as `mean()`,
#' return a plain `SpatRaster` instead. Re-wrapping is cheap either way, since
#' `as_voxel()` is idempotent — given a valid `SpatVoxel` and no `depths` it
#' returns it untouched.
#'
#' Propagating the class does not re-run the validity rules, so an operation
#' that changes the raster layer set can leave an object still labelled `SpatVoxel`
#' even if invalid. Passing an invalid `SpatVoxel` to `as_voxel()` rebuilds
#' from layer names, passing through the function as if a plain multi-layer raster. 
#' Pass `depths` to rebuild the layer names.
#'
#' Depths are positive metres increasing downward, matching the World Ocean
#' Atlas convention. Negative depths are an error rather than being silently
#' negated, since flipping the sign would change what the data mean.
#'
#' @param x SpatRaster with `{variable}_depth={value}` layer names, a list of
#'   single-depth SpatRasters, or an existing [SpatVoxel-class].
#' @param depths Optional numeric vector, one depth per layer, in metres. When
#'   supplied, layer names are (re)built from `depths` and `varname`, replacing
#'   any existing names. Required if `x` has no conforming layer names.
#' @param varname Character. Variable name used when building layer names from
#'   `depths`. Ignored when `depths` is `NULL`.
#'
#' @returns A [SpatVoxel-class] whose layers are ordered shallow to deep.
#' @examples
#' r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
#' terra::values(r) <- runif(terra::ncell(r) * 3)
#'
#' # build the layer names from a depth vector
#' v <- as_voxel(r, depths = c(0, 100, 200), varname = "temp")
#' names(v)
#'
#' # already-conforming names are used as they stand, and sorted if needed
#' names(r) <- c("temp_depth=200", "temp_depth=0", "temp_depth=100")
#' names(as_voxel(r))
#' @export
as_voxel <- function(x, depths = NULL, varname = "value") {
  # Idempotent: a valid SpatVoxel comes back as it is, so re-wrapping is cheap
  # and callers can do it defensively. terra propagates the class tag through
  # `[[`, `names<-` and the like without re-running validity, so a tagged
  # object may no longer be a voxel. `as_voxel()` repairs what is unambiguous 
  # (layer order) and rejects the rest with the constructor's own messages. 
  if (methods::is(x, "SpatVoxel") && is.null(depths)) {
    if (isTRUE(methods::validObject(x, test = TRUE))) return(x)
    x <- .as_plain_raster(x)
  }

  if (is.list(x)) x <- terra::rast(x)
  if (!methods::is(x, "SpatRaster")) {
    stop("`x` must be a SpatRaster or a list of SpatRasters; got ",
         paste(class(x), collapse = "/"), ".", call. = FALSE)
  }

  if (!is.null(depths)) {
    if (length(depths) != terra::nlyr(x)) {
      stop("`depths` must have one value per layer: got ", length(depths),
           " for ", terra::nlyr(x), " layers.", call. = FALSE)
    }
    if (!is.numeric(depths)) {
      stop("`depths` must be numeric metres.", call. = FALSE)
    }
    names(x) <- paste0(varname, "_depth=", depths)
  }

  d <- .parse_depth_layers(x, error = FALSE)
  if (anyNA(d)) {
    stop("layer(s) ", toString(names(x)[is.na(d)]),
         " do not follow the '{variable}_depth={value}' convention. ",
         "Pass `depths` to build the layer names instead.", call. = FALSE)
  }
  if (any(d < 0)) {
    stop("depths are positive metres increasing downward; got ",
         toString(d[d < 0]),
         ". Negate the depths in the layer names or in `depths`.",
         call. = FALSE)
  }
  if (anyDuplicated(d)) {
    stop("duplicate depth(s): ", toString(unique(d[duplicated(d)])),
         ". Each layer must be a distinct depth.", call. = FALSE)
  }

  # Repair what is unambiguously repairable rather than rejecting it.
  if (is.unsorted(d)) x <- x[[order(d)]]

  methods::new("SpatVoxel", x)
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
#'   \item `values` is the magnitude, which can be a single number for the whole grid or a
#'     SpatRaster carrying a value per cell. Default value is `1`, which makes the
#'     output a presence mask. Ex. `values` can be 2D fisheries effort rasters. 
#'   \item `profile` decides how that magnitude is spread down the water
#'     column. Default `NULL` means the full value is written at every occupied level. 
#'     Providing a `profile_*` function determines how `values` is transformed for each 
#'     depth level. Ex. [profile_equal()] divides it evenly over the levels the cell occupies,
#'     so sum of all depth levels in the voxel equals `values`.
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
#'   \item `"top"` (the default) runs each slab from its level down to the next,
#'     so a level names the top of what it stands for. Levels `c(0, 100, 200)`
#'     give slabs `0-100`, `100-200`.
#'   \item `"midpoint"` puts the edges halfway between neighbouring levels, so a
#'     level sits in the middle of what it stands for. This is the World Ocean
#'     Atlas convention: `0, 5, 10, ...` the 0 m layer covers `0-2.5` m, the 5 m 
#'     layer `2.5-7.5` m, the 10 m layer `7.5-12.5` m, and so on.
#' }
#'
#' Under either convention the slabs cover the span from the shallowest level to
#' the deepest and no further: outer edges are clamped to `range(depths)` rather
#' than extrapolated, which is why WOA's 0 m layer starts at the surface rather
#' than half a gap above it. The expansion is therefore limited only at the two
#' ends when a cell whose envelope lies wholly above `min(depths)` or wholly below
#' `max(depths)`. Under `"top"` the deepest level has no next level to run to,
#' so it stands for a zero-thickness slab: an envelope is recorded there only
#' when it reaches below `max(depths)`, and one that ends exactly at
#' `max(depths)` stops at the level above. For the same reason a
#' [voxel_to_envelope()] round trip drops each cell's deepest occupied level:
#' the envelope's `depth_max` names that level, which is the top of its slab.
#'
#' @param x SpatEnvelope, e.g. from [as_envelope()] or [voxel_to_envelope()].
#' @param depths Array of values that can be coerced into numeric type,
#'   represents metres depth below sea level. Sorted shallow to deep, and
#'   deduplicated, before use.
#' @param values Optional numeric or SpatRaster of same CRS, resolution, extent as `x`. 
#'   Defines values to write for voxel cells that are within depth intervals. 
#'   `NULL` (the default) writes `1` at every occupied level, giving simple 
#'   presence/absence.
#' @param profile Optional function distributing `values` across the depth
#'   dimension, taking `(ind, depths, n_depths)` and returning the weight to
#'   apply. `NULL` (the default) writes the same value at every depth, as
#'   [profile_flat()] does; [profile_equal()] divides the value equally across
#'   the depths a cell occupies. Any function meeting the contract in
#'   [voxel_profiles] works, including your own.
#' @param bounds How to place the edges of the slab each depth level stands for.
#'   `"top"` (the default) runs the slab from the level down to the next, so the
#'   level is its shallow edge; `"midpoint"` puts the edges halfway between
#'   neighbouring levels. Under both, the outer edges are clamped to `range(depths)` 
#'   rather than extrapolated.
#' @param varname Name to use for depth layer, taking on form of `{varname}_depth={depths[i]}`
#'
#' @returns A [SpatVoxel-class] with one layer per depth in `depths`, on the
#'   grid of `x`, layers ordered shallow to deep.
#'
#' @seealso [voxel_to_envelope()], the reverse (and lossy) collapse.
#'
#' @examples
#' fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
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
#' # the same envelope read against WOA-style layer bounds: [60, 70] sits in the
#' # 50-150 m slab of the 100 m level, where the default puts it at 0-100 m
#' shallow <- as_envelope(fp, depth_min = 60, depth_max = 70)
#' terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200)))
#' terra::values(envelope_to_voxel(shallow, depths = c(0, 100, 200),
#'                                 bounds = "midpoint"))
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
                              bounds = c("top", "midpoint"),
                              varname = "presence") {
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
  # requested levels are put in that form up front.
  depths <- sort(unique(depths))

  # Step 1: occupancy. Each level stands for a slab of water, `bounds` says
  # where its edges fall, and a cell occupies the level when its envelope
  # shares some thickness of water with that slab. Strict `>` on the
  # envelope's floor: an envelope that merely ends where a slab begins shares
  # no water with it, the same rule intersect_3d() applies between envelopes,
  # so adjacent depth ranges never meet on a level. Slabs are half-open so a
  # shared edge is claimed by the deeper level alone, the deepest closed.
  # Raster stays left: terra mishandles `scalar <op> raster`.
  slab <- .depth_slabs(depths, bounds)
  # Plain layers, so the occupancy stack and everything derived from it is a
  # plain SpatRaster rather than an object still tagged SpatEnvelope. terra
  # propagates the class through `[[` and arithmetic, and a tagged
  # intermediate would reach the depth-aware mask() methods in R/intersect.R.
  lower <- .as_plain_raster(x[["depth_min"]])
  upper <- .as_plain_raster(x[["depth_max"]])
  deepest <- length(depths)
  ind <- terra::rast(lapply(seq_along(depths), function(i) {
    if (i == deepest) {
      (lower <= slab$upper[i]) & (upper > slab$lower[i])
    } else {
      (lower < slab$upper[i]) & (upper > slab$lower[i])
    }
  }))
  # number of depth levels present per vertical column, cell 
  n_depths <- sum(ind)

  # The slabs tile `[min(depths), max(depths)]`, so a cell is left empty only
  # when its envelope lies wholly above the shallowest level or wholly below
  # the deepest — a range the levels cannot reach, not one that fell between
  # them.
  # `n_depths` is NA where the envelope is, so this counts only real cells.
  n_missed <- terra::global(n_depths == 0, "sum", na.rm = TRUE)[1, 1]
  if (!is.na(n_missed) && n_missed > 0) {
    warning(n_missed, " cell(s) have an envelope lying entirely outside ",
            "`depths` and are empty in the voxel. Extend `depths` to span ",
            "them.", call. = FALSE)
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
  names(out) <- paste0(varname, "_depth=", depths)

  as_voxel(out)
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
.depth_slabs <- function(depths, bounds) {
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
#' **This conversion is lossy.** An envelope stores a single continuous
#' interval per cell, so any interior gap in the voxel is filled in: a cell that
#' satisfies `fun` at 0 m and 200 m but not at 100 m still yields the envelope
#' `[0, 200]`. Use [SpatVoxel-class] directly where interior gaps matter.
#'
#' @param v SpatVoxel (or a multi-depth SpatRaster with
#'   `{variable}_depth={value}` layer names).
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
#' r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1,
#'                  nlyrs = 4)
#' terra::values(r) <- cbind(c(12, NA), c(11, NA), c(NA, 8), c(6, NA))
#' v <- as_voxel(r, depths = c(0, 50, 100, 200), varname = "temp")
#'
#' # Default predicate: the vertical extent of the data. Cell 1 has values at
#' # 0, 50 and 200 m, so it comes back as [0, 200] — the gap at 100 m is filled,
#' # because an envelope stores one continuous interval per cell.
#' terra::values(voxel_to_envelope(v))
#'
#' # Another predicate bounds a subset of the values instead: here the depths
#' # over which a cell is warmer than 10 degrees. Cell 2 never qualifies, so it
#' # is NA in both layers.
#' terra::values(voxel_to_envelope(v, fun = function(x) x > 10))
#' @export
voxel_to_envelope <- function(v, fun = function(x) !is.na(x)) {
  if (!is.function(fun)) {
    stop("`fun` must be a function returning TRUE/FALSE for a cell value. ",
         "The former \"extent\" behaviour is the default, function(x) !is.na(x).",
         call. = FALSE)
  }
  if (!is(v, "SpatVoxel")) {
    stop("`v` must be a SpatVoxel object.")
  }
  .check_3d(v, "v")
  depths <- .parse_depth_layers(v)
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
  depth_min <- depths[first_hit]
  depth_max <- depths[last_hit]
  depth_min[never] <- NA_real_
  depth_max[never] <- NA_real_

  out <- terra::setValues(terra::rast(.as_plain_raster(v[[1]]), nlyrs = 2),
                          cbind(depth_min, depth_max))
  names(out) <- c("depth_min", "depth_max")

  out <- methods::new("SpatEnvelope", out)
  methods::validObject(out)
  out
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
#' @param x A [SpatVoxel-class].
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
  .presence_voxel(.voxel_occupancy(x, fun), .parse_depth_layers(x))
}

# Internal: validate and normalise the `depth_min` / `depth_max` inputs of
# `vect_to_envelope()`. Accepts a list, a numeric vector, or a SpatRaster
# Returns a list whose elements can be used in `do.call("max", ...)` / `do.call("min", ...)`.
.normalize_depth_inputs <- function(x, template, arg) {
  if (inherits(x, "SpatRaster")) {
    x <- list(x)
  } else if (!is.list(x)) {
    x <- as.list(x)
  }

  if (length(x) == 0L) {
    stop("`", arg, "` is empty. Supply at least one numeric value or SpatRaster.", call. = FALSE)
  }

  for (i in seq_along(x)) {
    el <- x[[i]]
    label <- if (length(x) == 1L) {
      paste0("`", arg, "`")
    } else {
      paste0("`", arg, "` element ", i)
    }

    if (inherits(el, "SpatRaster")) {
      if (!terra::same.crs(el, template)) {
        stop(label, " has a different CRS than `template`. ",
             "Project it onto the template grid first.", call. = FALSE)
      }
      aligned <- terra::compareGeom(el, template, crs = FALSE,
                                    stopOnError = FALSE, messages = FALSE)
      if (!isTRUE(aligned)) {
        stop(label, " does not align with `template`; its extent or resolution ",
             "differs. Resample it onto the template grid first.", call. = FALSE)
      }
    } else if (!is.numeric(el)) {
      stop(label, " must be numeric or a SpatRaster, not ", class(el)[1], ".",
           call. = FALSE)
    } else if (anyNA(el)) {
      # A scalar NA carries no information but, under `na.rm = FALSE`, would
      # blank every cell of the output. Reject it here where we can say which
      # argument was at fault instead of reporting an empty envelope later.
      stop(label, " is NA. Drop it, or supply a real depth constraint.",
           call. = FALSE)
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
  is_rast <- vapply(x, function(el) inherits(el, "SpatRaster"), logical(1))

  # skips the reordered do.call when there are no rasters in input list x
  # creates raster with single value across
  if (!any(is_rast)) {
    return(terra::setValues(terra::rast(template[[1]]), do.call(fun, x)))
  }
  # terra has min() and max() functions that can apply to a list that includes
  # both SpatRaster and numeric values. However, it only works when the SpatRaster
  # is first in the list; the base::min is called instead if the first element is 
  # a numeric object. This do.call line takes the reordered input list, with the 
  # SpatRaster objects first. 

  # `na.rm = FALSE` so that cell where constraints / depth raster
  # has NA value keeps it in the output raster. Avoids filling the value from 
  # remaining constraints
  do.call(fun, c(x[is_rast], x[!is_rast], list(na.rm = FALSE)))
}

#' Convert SpatVector or sf to SpatEnvelope
#' 
#' Rasterize polygons onto a template raster grid and assign per-cell depth limits,
#' outputting an SpatEnvelope object. The depths are clamped by depth_min (shallowest) 
#' and depth_max (deepest), and/or by rasters (ex. bathymetry) that are the same 
#' coordinate and resolution as the template. Where a constraint raster is NA at
#' a cell its limit there is unknown, so the cell is dropped rather than falling
#' back to the remaining constraints.
#' 
#' @param polygon sf or SpatVector, for example species ranges or fishery footprints. 
#' @param template SpatRaster that defines the horizontal grid of the SpatEnvelope output. 
#' @param depth_min List. Can contain both numeric and SpatRasters that match CRS, 
#'   resolution, extent of `template`` and contain numeric values. For each cell, the maximum
#'   across the `depth_min` list parameters is used as output SpatEnvelope depth_min layer cell
#'   value. 
#' @param depth_max List. Can contain both numeric and SpatRasters that match CRS, 
#'   resolution, extent of `template`` and contain numeric values. For each cell, the minimum
#'   across the `depth_max` list parameters is used as output SpatEnvelope depth_max layer cell
#'   value. Note that when `depth_min` is exactly `depth_max`, these cell values are dropped and
#'   replaced with NA.
#' 
#' @returns SpatEnvelope, with depth_min and depth_max layers. 
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
    stop("`polygon` needs to be of class SpatVector, sf, or sfc")
  }
  if (inherits(polygon, "sf") || inherits(polygon, "sfc")) {
    polygon <- terra::vect(polygon)
  }
  if (!is(template, "SpatRaster")) {
    stop("`template` needs to be of class SpatRaster")
  }
  # Only the grid is used. Drop any class tag so a voxel can serve as the
  # template without its layers reaching the depth-aware mask() methods.
  template <- .as_plain_raster(template[[1]])
  if(!terra::same.crs(crs(polygon), crs(template))) {
    stop("`polygon` and `template` have different CRS.")
  }

  # check the depth inputs for CRS, resolution, extent matching and reject NA params
  depth_min <- .normalize_depth_inputs(depth_min, template, "depth_min")
  depth_max <- .normalize_depth_inputs(depth_max, template, "depth_max")

  # take first layer in case multi-layer raster
  all_one <- template[[1]]
  # make all cell values 1, since mask leaves NAs alone 
  terra::values(all_one) <- 1
  # intersect between polygon and template
  masked <- terra::mask(all_one, polygon) 
  if(terra::global(masked, "notNA")[1,1] == 0) {
    stop("All output cell values are NA. `polygon` and `template` may not spatially overlap.")
  }

  # Creating depth_min, depth_max rasters that form the SpatEnvelope output
  # - value that narrows the envelope the most is selected
  # - preserves NA values from the constraints, rather than fall back on other constraints
  # deepest of depth_min constraints
  depth_min_rast <- .combine_depths(depth_min, template, "max") %>%
    mask(masked)
  # shallowest of depth_max constraints
  depth_max_rast <- .combine_depths(depth_max, template, "min") %>%
    mask(masked)

  # check for validity, where depth_max is greater than depth_min
  diffs <- depth_max_rast - depth_min_rast
  msk <- ifel(diffs > 0, 1, NA)
  if(terra::global(msk, "notNA")[1,1] == 0) {
    warning("All depth_max values are shallower than depth_min. Check `depth_min` and `depth_max`, you may have swapped these two.")
  }

  # mask cells that are not valid, depth_min is deeper than depth_max
  depth_min_rast <- mask(depth_min_rast, msk)
  depth_max_rast <- mask(depth_max_rast, msk)

  # assign SpatEnvelope class to output
  out <- rast(c(depth_min = depth_min_rast, depth_max = depth_max_rast)) %>% as_envelope()
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
