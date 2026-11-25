# Values at many points from a SpatProxy at once (#163): the fast path the
# per-shape loop (#115) lacks. Points are matched to time steps together,
# each step is read once for its points, and each point's cell and level are
# picked with the point queries' rules (.point_levels()), so the values are
# those the loop gives.

#' Values at many points from a SpatProxy
#'
#' `extract()` reads the value of a gridded variable at each point of a
#' [SpatVector3D-class], for many points at once. It gives the values the
#' per-shape loop gives, `mask(as_voxel(crop(proxy, p), time = p), p)` for
#' each point `p`, much faster: the points are matched to time steps
#' together, and each time step is read once, cropped to its points.
#'
#' Each point takes the time step nearest its time, or the step equal to it
#' with `time_match = "exact"`, as [as_voxel()] matches a point. With dates,
#' a point's date matches by UTC day. A source with one time step, such as a
#' World Ocean Atlas climatology, gives that step to every point, whatever
#' its time. On a source with several steps, a point without a time gets
#' `NA`, with a warning.
#'
#' Each point takes the cell that holds it and, for a variable with depth,
#' one depth level, chosen with `method` as [mask()][mask-3d] chooses it for
#' points: `"nearest"`, the nearest level that has a value, as
#' `extract_to_point(method = "nearest")` picks it; `"nearest_level"`, the
#' level nearest the point's depth whatever its value; `"bottom"`, the
#' deepest level with a value; or `"surface"`, the shallowest level. For a
#' variable without depth, each point takes the value in its cell, and
#' `method` doesn't apply.
#'
#' A time step whose points spread over more cells than the read budget
#' allows (the option `ocean3d.max_read_mb`, 1024 MB by default) is read in
#' smaller groups of nearby points, so memory follows the budget; the
#' values don't change.
#'
#' @inheritSection SpatProxy-class Experimental
#'
#' @param x A [SpatProxy-class], from [as_proxy()]: a netCDF file, a series
#'   of files, or a voxel in memory.
#' @param y [SpatVector3D-class] points, from [as_vector3d()].
#' @param method For a variable with depth, the depth level each point
#'   takes: `"nearest"` (default), `"nearest_level"`, `"bottom"` or
#'   `"surface"`, as in [mask()][mask-3d].
#' @param time_match,max_time_diff How a point's time is matched to a time
#'   step: `"nearest"` (default) or `"exact"`, and for nearest matching, the
#'   largest gap in days, as in [as_voxel()]. `NULL` (default) sets no
#'   limit.
#' @param ... Not used.
#'
#' @returns A data frame with one row per point, in the order of `y`:
#'   `feature_id`; `value`; `time`, the date-time of the time step used
#'   (`NA` for a source without dated steps); and `depth`, the depth of the
#'   level used in metres (`NA` for a variable without depth). A point off
#'   the grid, without a matching time step or with no level its `method`
#'   can take has `NA` for `value`.
#'
#' @seealso [mask()][mask-3d] and [as_voxel()] for one shape at a time;
#'   [extract_to_point()].
#' @examples
#' p <- as_proxy(system.file("extdata", "example_3d.nc", package = "ocean3d"))
#' pts <- as_vector3d(data.frame(
#'   lon = c(0.1, 0.9), lat = c(40.2, 40.8), depth = c(5, 60),
#'   date = as.Date(c("2020-01-01", "2020-01-05"))
#' ), time = "date")
#' terra::extract(p, pts)
#'
#' # The deepest level with a value at each point instead.
#' terra::extract(p, pts, method = "bottom")
#' @name extract-SpatProxy
#' @aliases extract,SpatProxy,SpatVector3D-method
#' @export
setMethod(
  "extract", c("SpatProxy", "SpatVector3D"),
  function(x, y, method = c("nearest", "nearest_level", "bottom", "surface"),
           time_match = c("nearest", "exact"), max_time_diff = NULL, ...) {
    .no_extra_args(...)
    .check_geomtype(y, "points")
    has_depth <- !is.null(x@axes$depth)
    rule <- if (missing(method)) "nearest" else .check_level_rule(method)
    # A cell holds one value of a variable without depth, so there is no
    # level to choose.
    if (!has_depth && !missing(method)) .no_extra_args(method = method)
    time_match <- match.arg(time_match)
    .check_max_time_diff(max_time_diff)
    budget <- .read_budget_mb()

    v <- y@vertices
    n <- nrow(v)
    # Each point's time step: 0 for a source without a time axis (read
    # whole), 1 for a source with one step, else the step matched to its
    # time.
    axis <- x@axes$time
    n_steps <- if (is.null(axis)) 1L else x@variable$dim_lengths[[axis$dim]]
    dates <- if (!is.null(axis) && !isTRUE(axis$climatology)) {
      axis$cftime$as_timestamp(asPOSIX = TRUE)
    }
    step <- if (is.null(axis)) {
      rep(0L, n)
    } else if (n_steps == 1L) {
      rep(1L, n)
    } else {
      untimed <- sum(is.na(v$time))
      if (untimed > 0L) {
        have <- if (untimed == 1L) "point has" else "points have"
        rlang::warn(glue::glue(
          "{untimed} of {n} {have} no time, but variable ",
          "`{x@variable$name}` has {n_steps} time steps, so their values ",
          "are NA. Give the points a time in as_vector3d()."
        ), class = "ocean3d_warning_point_time")
      }
      # Dates compare by UTC day, date-times exactly, as for a point given
      # to as_voxel(); several steps on a date's day can't be told apart.
      resolution <- if (is.na(y@time_resolution)) "datetime" else
        y@time_resolution
      tryCatch(
        .match_time_indices(dates, v$time, time_match, max_time_diff,
                            resolution),
        ocean3d_error_repeated_time = function(e) {
          .abort_repeated_time(
            "Several time steps fall on the day of some points' dates", NULL,
            "Give the points date-times (with hours) in as_vector3d()."
          )
        }
      )
    }

    value <- rep(NA_real_, n)
    level_depth <- rep(NA_real_, n)
    # Groups of points read together: those on one time step, split into
    # nearby halves while a group's cells would read more than the budget.
    # A group whose points all share one position can't be split, and its
    # cell always fits.
    groups <- unname(split(seq_len(n)[!is.na(step)], step[!is.na(step)]))
    while (length(groups)) {
      idx <- groups[[1L]]
      groups <- groups[-1L]
      points <- y[idx]
      cropped <- tryCatch(crop(x, points),
                          ocean3d_error_proxy_crop = function(e) NULL)
      if (is.null(cropped)) next # every point of the group is off the grid
      sizes <- cropped@variable$dim_lengths
      dims <- unlist(cropped@variable$axis_dims[c("lon", "lat", "depth")])
      mb <- prod(sizes[dims]) * 8 / 2^20
      lon <- v$lon[idx]
      lat <- v$lat[idx]
      if (mb > budget && (length(unique(lon)) > 1L ||
                            length(unique(lat)) > 1L)) {
        along <- if (diff(range(lon)) >= diff(range(lat))) lon else lat
        first <- along <= stats::median(along)
        # All on the median: split by position in the group instead.
        if (all(first)) first <- seq_along(idx) <= length(idx) %/% 2L
        groups <- c(groups, list(idx[first], idx[!first]))
        next
      }
      s <- step[idx[1L]]
      r <- .proxy_steps(cropped, if (s == 0L) NULL else s, x@variable$name,
                        "nearest", NULL, voxel = has_depth)
      if (has_depth) {
        picked <- .point_levels(points, r, rule)
        on <- which(!is.na(picked$level))
        # remove the 3D class
        layers <- terra::extract(methods::as(r, "SpatRaster"),
                                 picked$cell[on])
        value[idx[on]] <- as.matrix(layers)[cbind(seq_along(on),
                                                  picked$level[on])]
        level_depth[idx[on]] <- depths(r)[picked$level[on]]
      } else {
        cells <- .point_cells(points, r)
        on <- which(!is.na(cells))
        value[idx[on]] <- terra::extract(r, cells[on])[, 1L]
      }
    }
    data.frame(
      feature_id = v$feature_id,
      value = value,
      time = if (is.null(dates)) .POSIXct(rep(NA_real_, n), tz = "UTC") else
        dates[as.integer(step)],
      depth = level_depth,
      stringsAsFactors = FALSE
    )
  }
)
