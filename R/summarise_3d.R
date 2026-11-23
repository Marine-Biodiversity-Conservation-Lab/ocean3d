# The summarise step for voxels (#116): statistics of a SpatVoxel's values
# across its cells and depth levels, weighted by the volume of each cell and
# level, as the per-shape workflow needs after mask(v, shape). It reads no
# files and handles no geometry: the weights come from volume().

#' Summarise a voxel's values across cells and depth levels
#'
#' `summarise_3d()` computes statistics of every non-`NA` value of a
#' [SpatVoxel-class], across all its cells and depth levels at once. It is
#' the summarise step of an analysis: the values inside a shape are kept with
#' [mask()], and `summarise_3d()` turns them into a mean, a standard
#' deviation, a range or a median. [terra::global()] summarises each layer,
#' which here is each depth level, separately.
#'
#' With `weighting = "volume"` (the default), each value is weighted by the
#' volume of its cell *and* depth level, from
#' `volume(x, per_cell = TRUE, per_level = TRUE)`: cells shrink towards the
#' poles, and depth levels range from a few metres to 100 m or more thick. A
#' mean over a water column with thin warm levels near the surface and thick
#' cold ones below is colder weighted by volume than with equal weights.
#' `weighting = "equal"` counts each value once, so `sd` and `median` are
#' then [stats::sd()] and [stats::median()].
#'
#' The weighted standard deviation is the one [extract_to_buffer()] reports:
#' with normalised weights \eqn{w_i}, the square root of
#' \eqn{\sum w_i (v_i - \bar v)^2 / \sum w_i (1 - w_i)}, which is the sample
#' standard deviation under equal weights. The weighted median is the value
#' at which the cumulative weight, over the values in order, reaches half the
#' total; when it reaches half exactly, the two values either side are
#' averaged, as [stats::median()] does. `min` and `max` ignore the weights.
#'
#' @param x A [SpatVoxel-class], typically `mask(v, shape)` for the values
#'   inside a shape, or a list of voxels, one per time step, as [as_voxel()]
#'   gives for a time window. Envelopes have no values, so they aren't
#'   accepted; a netCDF file isn't either: describe it with [as_proxy()] and
#'   read one time step with [as_voxel()].
#' @param statistic One or more of `"mean"`, `"sd"`, `"min"`, `"max"` and
#'   `"median"`. Each becomes a column of the result, in this order.
#' @param weighting `"volume"` (default) weights each value by the volume of
#'   its cell and depth level; `"equal"` counts each value once.
#' @param by `NULL` (default) gives one row per voxel, so one per time step
#'   for a list; `"depth"` one row per depth level of each voxel, with a
#'   `depth` column. Time steps are never pooled, so a list keeps the change
#'   from step to step. `"time"` means the same as `NULL`.
#' @param bounds As in [volume()]: how the slab edges are placed for a voxel
#'   without depth edges, for `weighting = "volume"`. Voxels read with
#'   [as_voxel()] from a [SpatProxy-class] carry their edges, so it rarely
#'   applies.
#'
#' @returns A data frame with one row per voxel, or per voxel and depth level
#'   with `by = "depth"`, in the order of `x` and then of depth. Columns:
#'   `time` (the voxel's [terra::time()], UTC, or `NA` when it has none),
#'   `depth` (with `by = "depth"` only, in metres), one column per
#'   `statistic`, and `n`, the number of values the row is built on. A row
#'   with no values has `NA` statistics and `n = 0`, so summaries built on
#'   too few values can be dropped.
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
#' @export
summarise_3d <- function(x, statistic = "mean",
                         weighting = c("volume", "equal"), by = NULL,
                         bounds = c("midpoint", "top")) {
  statistic <- match.arg(statistic, c("mean", "sd", "min", "max", "median"),
                         several.ok = TRUE)
  weighting <- match.arg(weighting)
  bounds <- match.arg(bounds)
  if (!is.null(by) && !identical(by, "depth") && !identical(by, "time")) {
    .abort_summarise_input(
      "`by` must be NULL (one row per voxel), \"depth\" or \"time\"."
    )
  }
  voxels <- if (is.list(x)) x else list(x)
  is_voxel <- vapply(voxels, methods::is, logical(1), "SpatVoxel")
  if (!length(voxels) || !all(is_voxel)) {
    got <- if (is.list(x)) {
      "a list with other objects in it"
    } else {
      paste(class(x), collapse = "/")
    }
    .abort_summarise_input(glue::glue(
      "`x` must be a SpatVoxel or a list of SpatVoxels; got {got}. Build one ",
      "with as_voxel(), and keep the values inside a shape with mask()."
    ))
  }

  rows <- lapply(voxels, function(v) {
    .check_3d(v, "x", "SpatVoxel")
    # Cells x depth levels. The per-level volume is NA exactly where the
    # voxel is, so the weights line up with the values.
    values <- terra::values(v, mat = TRUE)
    weights <- if (weighting == "volume") {
      terra::values(volume(v, per_cell = TRUE, per_level = TRUE,
                           bounds = bounds), mat = TRUE)
    } else {
      array(1, dim(values))
    }
    # A date is midnight UTC; a voxel without a time gives NA.
    time <- if (isTRUE(terra::timeInfo(v)$time)) {
      as.POSIXct(terra::time(v)[1L], tz = "UTC")
    } else {
      .POSIXct(NA_real_, tz = "UTC")
    }
    groups <- if (identical(by, "depth")) {
      lapply(seq_len(ncol(values)), function(i) {
        c(list(depth = depths(v)[i]),
          .weighted_statistics(values[, i], weights[, i], statistic))
      })
    } else {
      list(.weighted_statistics(values, weights, statistic))
    }
    do.call(rbind, lapply(groups, function(g) {
      data.frame(time = time, as.data.frame(g), stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows)
  attr(out$time, "tzone") <- "UTC"
  rownames(out) <- NULL
  out
}

# Internal: the statistics of the non-NA `values` with their `weights`, as a
# named list ending in `n`, the number of values. Called per voxel and per
# depth level by summarise_3d(). The weighted sd is the one
# .buffer_statistic() (R/extract_to_buffer.R) reports, so the two agree.
.weighted_statistics <- function(values, weights, statistic) {
  keep <- !is.na(values)
  values <- as.numeric(values[keep])
  weights <- as.numeric(weights[keep])
  n <- length(values)
  out <- stats::setNames(rep(NA_real_, length(statistic)), statistic)
  if (n) {
    w <- weights / sum(weights)
    average <- sum(w * values)
    for (s in statistic) {
      out[[s]] <- switch(
        s,
        mean = average,
        min = min(values),
        max = max(values),
        sd = {
          correction <- sum(w * (1 - w))
          if (n < 2L || correction <= 0) NA_real_ else
            sqrt(sum(w * (values - average)^2) / correction)
        },
        median = {
          o <- order(values)
          sorted <- values[o]
          cumulative <- cumsum(w[o])
          # The first value at which the cumulative weight reaches half;
          # exactly half (within rounding) averages it with the next, as
          # stats::median() does for an even count.
          at <- which(cumulative >= 0.5 - 1e-10)[1L]
          if (abs(cumulative[at] - 0.5) < 1e-10 && at < n) {
            (sorted[at] + sorted[at + 1L]) / 2
          } else {
            sorted[at]
          }
        }
      )
    }
  }
  c(as.list(out), list(n = n))
}

# The one producer of `ocean3d_error_summarise_input`: input summarise_3d()
# can't summarise.
.abort_summarise_input <- function(message) {
  rlang::abort(message, class = "ocean3d_error_summarise_input", call = NULL)
}
