# The summarise step for voxels (#116): statistics of a SpatVoxel's values
# across its cells and depth levels, weighted by the volume of each cell and
# level, as the per-shape workflow needs after mask(v, shape). It reads no
# files and handles no geometry: the weights come from volume(). The weighted
# statistics are shared with buffer extraction.

#' Summarise a voxel's values across cells and depth levels
#'
#' `summarise_3d()` computes statistics of every non-`NA` value of a
#' [SpatVoxel-class], across all its cells and depth levels at once. It is
#' the summarise step of an analysis: the values inside a shape are kept with
#' [mask()], and `summarise_3d()` turns them into a mean, a standard
#' deviation, a range or a median. Compare [terra::global()], which
#' summarises each layer, here each depth level, on its own.
#'
#' With `weighting = "volume"` (the default), each value is weighted by the
#' volume of its cell *and* depth level, from
#' `volume(x, per_cell = TRUE, per_level = TRUE)`: cells shrink towards the
#' poles, and depth levels range from a few metres thick to 100 m or more.
#' In a water column with thin warm levels near the surface and thick cold
#' ones below, the volume-weighted mean is colder than the equal-weighted one.
#' `weighting = "equal"` counts each value once, so `sd` and `median` are
#' then [stats::sd()] and [stats::median()].
#'
#' With normalised weights \eqn{w_i}, the weighted standard deviation is the
#' square root of \eqn{\sum w_i (v_i - \bar v)^2 / \sum w_i (1 - w_i)}, the
#' sample standard deviation under equal weights and the same formula as in
#' [extract_to_track()]. The weighted median is the value at which the
#' cumulative weight, over the values sorted from smallest to largest,
#' reaches half the total; when it reaches half exactly, the two values
#' either side are averaged, as [stats::median()] does. `min` and `max`
#' ignore the weights.
#'
#' @param x A [SpatVoxel-class], typically `mask(v, shape)` for the values
#'   inside a shape, or a list of voxels, one per time step, as [as_voxel()]
#'   gives for a time window. `x` can also be a one-layer `SpatRaster` of a
#'   variable without depth, as [rast()][rast-SpatProxy] reads from a
#'   [SpatProxy-class], or a list of them, one per time step. Envelopes have
#'   no values, so they aren't accepted; a netCDF file isn't either: describe
#'   it with [as_proxy()] and read one time step with [as_voxel()] or
#'   `rast()`.
#' @param statistic One or more of `"mean"`, `"sd"`, `"min"`, `"max"` and
#'   `"median"`. Each becomes a column of the result, in this order.
#' @param weighting `"volume"` (default) weights each value by the volume of
#'   its cell and depth level, and each value of a one-layer `SpatRaster` by
#'   its cell's area; `"equal"` counts each value once. A row whose
#'   values all have zero volume, such as a voxel with a single level and no
#'   depth edges, or the deepest level under `bounds = "top"`, uses equal
#'   weights.
#' @param by `NULL` (default) gives one row per voxel, so one per time step
#'   for a list; `"depth"` one row per depth level of each voxel, with a
#'   `depth` column. Values from different time steps aren't combined, so a
#'   list shows the change from step to step. `"time"` means the same as
#'   `NULL`, and says in a call on a list that the rows are per time step.
#'   `"all"` gives one row for every value of every voxel in `x` together,
#'   for example all the time steps a time window covers, each value keeping
#'   its own weight; its `time` is `NA`. A one-layer `SpatRaster` has a
#'   single level, so `"depth"` is an error for it.
#' @param bounds As in [volume()]: where the slab edges fall for a voxel that
#'   carries no depth edges; used with `weighting = "volume"` only. Voxels
#'   read with [as_voxel()] from a [SpatProxy-class] carry their edges, so
#'   `bounds` isn't used for them.
#'
#' @returns A data frame with one row per voxel, or per voxel and depth level
#'   with `by = "depth"`, in the order of `x` and then of depth, or one row
#'   with `by = "all"`. Columns:
#'   `time` (the voxel's [terra::time()], UTC, or `NA` when it has none),
#'   `depth` (with `by = "depth"` only, in metres), one column per
#'   `statistic`, and `n`, the number of values the row is built on, so rows
#'   built on too few values can be dropped. A row with no values has `NA`
#'   statistics and `n = 0`.
#'
#' @seealso [mask()] to keep the values inside a shape; [volume()] for the
#'   weights; [as_voxel()] to read one time step of a [SpatProxy-class].
#' @examples
#' # One water column sampled at seven depths: thin warm levels near the
#' # surface and thick cold ones below.
#' r <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0,
#'                  ymax = 1, nlyrs = 7, crs = "EPSG:4326")
#' terra::values(r) <- t(c(20, 19, 18, 17, 16, 8, 4))
#' v <- as_voxel(r, depths = c(0, 10, 20, 30, 40, 100, 200), varname = "temp")
#'
#' # The volume-weighted mean is colder: the deep levels stand for more water.
#' summarise_3d(v, statistic = c("mean", "sd", "min", "max"))
#' summarise_3d(v, statistic = "mean", weighting = "equal")
#'
#' # One row per depth level.
#' summarise_3d(v, by = "depth")
#'
#' # A list of voxels, one per time step, gives one row per step.
#' later <- v + 1
#' terra::time(v) <- rep(as.POSIXct("2020-01-01", tz = "UTC"), 7)
#' terra::time(later) <- rep(as.POSIXct("2020-01-02", tz = "UTC"), 7)
#' summarise_3d(list(v, later), statistic = c("mean", "median"))
#'
#' # Or one row for the values of both steps together.
#' summarise_3d(list(v, later), statistic = "mean", by = "all")
#' @export
summarise_3d <- function(x, statistic = "mean",
                         weighting = c("volume", "equal"), by = NULL,
                         bounds = c("midpoint", "top")) {
  statistic <- match.arg(statistic, c("mean", "sd", "min", "max", "median"),
                         several.ok = TRUE)
  weighting <- match.arg(weighting)
  bounds <- match.arg(bounds)
  if (!is.null(by) && !(rlang::is_string(by) &&
                           by %in% c("depth", "time", "all"))) {
    .abort_summarise_input(paste(
      "`by` must be NULL (one row per voxel), \"depth\", \"time\" or",
      "\"all\" (one row for every voxel together)."
    ))
  }
  # A 2D layer, as rast() reads a variable without depth: a plain one-layer
  # SpatRaster, not one of the 3D classes.
  flat <- function(r) {
    methods::is(r, "SpatRaster") && !methods::is(r, "SpatVoxel") &&
      !methods::is(r, "SpatEnvelope") && terra::nlyr(r) == 1L
  }
  # One voxel is checked as every 3D input is; a list must hold only voxels
  # or 2D layers.
  if (is.list(x)) {
    usable <- vapply(x, function(r) methods::is(r, "SpatVoxel") || flat(r),
                     logical(1))
    if (!length(x) || !all(usable)) {
      .abort_summarise_input(paste(
        "A list `x` must hold one or more SpatVoxels or one-layer",
        "SpatRasters, one per time step, and nothing else. Read each with",
        "as_voxel() or rast(), and keep the values inside a shape with mask()."
      ))
    }
    voxels <- x
  } else {
    if (!flat(x)) .check_3d(x, "x", "SpatVoxel")
    voxels <- list(x)
  }

  # Each voxel's values (cells x depth levels), their weights and its time.
  parts <- lapply(voxels, function(v) {
    two_d <- flat(v)
    if (!two_d) .check_3d(v, "x", "SpatVoxel")
    if (two_d && identical(by, "depth")) {
      .abort_summarise_input(paste(
        "A one-layer SpatRaster has a single level, so `by = \"depth\"`",
        "doesn't apply."
      ))
    }
    # Cells x depth levels. The per-level volume is NA exactly where the
    # voxel is, so the weights line up with the values; a 2D layer is
    # weighted by cell area.
    values <- terra::values(v, mat = TRUE)
    weights <- if (weighting == "equal") {
      array(1, dim(values))
    } else if (two_d) {
      terra::values(.cell_area_km2(v), mat = TRUE)
    } else {
      terra::values(volume(v, per_cell = TRUE, per_level = TRUE,
                           bounds = bounds), mat = TRUE)
    }
    list(values = values, weights = weights, time = .voxel_time(v),
         depths = if (!two_d) depths(v))
  })
  # "all": every value of every voxel in one row, each with its own weight.
  if (identical(by, "all")) {
    parts <- list(list(
      values = unlist(lapply(parts, `[[`, "values")),
      weights = unlist(lapply(parts, `[[`, "weights")),
      time = .POSIXct(NA_real_, tz = "UTC")
    ))
  }
  rows <- lapply(parts, function(p) {
    groups <- if (identical(by, "depth")) {
      lapply(seq_len(ncol(p$values)), function(i) {
        c(list(depth = p$depths[i]),
          .weighted_statistics(p$values[, i], p$weights[, i], statistic))
      })
    } else {
      list(.weighted_statistics(p$values, p$weights, statistic))
    }
    do.call(rbind, lapply(groups, function(g) {
      data.frame(time = p$time, as.data.frame(g), stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows)
  attr(out$time, "tzone") <- "UTC"
  rownames(out) <- NULL
  out
}

# The one producer of `ocean3d_error_summarise_input`: a `by` or a list `x`
# summarise_3d() can't use.
.abort_summarise_input <- function(message) {
  rlang::abort(message, class = "ocean3d_error_summarise_input", call = NULL)
}

# Internal: the statistics of the non-NA `values` with their `weights`, as a
# named list ending in `n`, the number of values. Called per voxel and per
# depth level by summarise_3d().
.weighted_statistics <- function(values, weights, statistic) {
  keep <- !is.na(values)
  values <- as.numeric(values[keep])
  weights <- as.numeric(weights[keep])
  n <- length(values)
  out <- if (n) {
    vapply(statistic, .weighted_statistic, numeric(1), values = values,
           weights = weights)
  } else {
    stats::setNames(rep(NA_real_, length(statistic)), statistic)
  }
  c(as.list(out), list(n = n))
}

# Internal: one weighted statistic of `values` (one or more, none NA) with
# their `weights`: "mean", "sd", "se", "min", "max" or "median". The sd is
# the square root of sum(w (v - mean)^2) / sum(w (1 - w)) with normalised
# weights, the sample sd under equal weights, and the se divides it by the
# square root of the effective number of values, 1 / sum(w^2). The median is
# the value at which the cumulative weight reaches half, averaged with the
# next when it lands exactly on half, as stats::median() does. Weights that
# sum to 0, such as volumes of a zero-thickness slab, count each value once.
# Shared by buffer extraction (.buffer_summarise_samples()) and
# summarise_3d().
.weighted_statistic <- function(values, weights, statistic) {
  if (statistic == "min") return(min(values))
  if (statistic == "max") return(max(values))

  if (sum(weights) <= 0) weights <- rep(1, length(values))
  weights <- weights / sum(weights)
  if (statistic == "median") {
    # A value with no weight can't be the middle, nor half of a middle pair.
    values <- values[weights > 0]
    weights <- weights[weights > 0]
    o <- order(values)
    sorted <- values[o]
    cumulative <- cumsum(weights[o])
    at <- which(cumulative >= 0.5 - 1e-10)[1L]
    if (abs(cumulative[at] - 0.5) < 1e-10 && at < length(values)) {
      return((sorted[at] + sorted[at + 1L]) / 2)
    }
    return(sorted[at])
  }
  average <- sum(values * weights)
  if (statistic == "mean") return(average)
  if (length(values) < 2L) return(NA_real_)

  correction <- sum(weights * (1 - weights))
  if (correction <= 0) return(NA_real_)
  deviation <- sqrt(sum(weights * (values - average)^2) / correction)
  if (statistic == "sd") return(deviation)

  effective_n <- 1 / sum(weights^2)
  deviation / sqrt(effective_n)
}
