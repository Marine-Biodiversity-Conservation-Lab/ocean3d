# Values at many points from a SpatProxy at once (#163): the fast path the
# per-shape loop (#115) lacks. Points are matched to time steps together
# (.match_proxy_steps(), as as_voxel() matches them), each step is read once
# for its points, and each point's cell and level are picked with the point
# queries' rules (.point_levels()), so the values are those the loop gives.

#' Values at many points from a SpatProxy
#'
#' `extract()` reads the value of a gridded variable at each point of a
#' [SpatVector3D-class], for many points at once. It gives the same values
#' as calling `mask(as_voxel(crop(proxy, p), time = p), p)` for each point
#' `p`, much faster: the points are matched to time steps together, and
#' each time step is read once, cropped to its points. The article
#' `vignette("workflow-values-in-shapes")` shows it in "Many points at
#' once".
#'
#' Each point takes the time step nearest its time, or the step equal to it
#' with `time_match = "exact"`, as [as_voxel()] matches a point. With dates,
#' a point's date matches by UTC day. A source with one time step, such as a
#' World Ocean Atlas climatology, gives that step to every point, whatever
#' its time. On a source with several steps, a point without a time gets
#' `NA`, with a warning.
#'
#' Each point takes the cell that holds it and, for a variable with depth,
#' one depth level, chosen with `method` as [mask()][mask-3d] does for points:
#' * `"nearest"` (default): the level nearest the point's depth, skipping
#'   levels without a value, as `extract_to_point(method = "nearest")` does.
#' * `"nearest_level"`: the level nearest the point's depth, even if it has
#'   no value.
#' * `"bottom"`: the deepest level with a value.
#' * `"surface"`: the shallowest level.
#'
#' For a variable without depth, each point takes the value in its cell, and
#' `method` doesn't apply.
#'
#' A time step's points are read together while they are close enough. A
#' group of points is read in smaller groups of nearby points when one read
#' would exceed the read budget (the option `ocean3d.max_read_mb`, 1024 MB
#' by default), or would hold more than 10,000 cells and depth levels for
#' each point, so that points far apart aren't read with the water between
#' them. The values are the same.
#'
#' @inheritSection SpatProxy-class Experimental
#'
#' @param x A [SpatProxy-class], from [as_proxy()]: a netCDF file, a series
#'   of files, or a voxel in memory.
#' @param y [SpatVector3D-class] points, from [as_vector3d()].
#' @param method For a variable with depth, the depth level each point
#'   takes: `"nearest"` (default), `"nearest_level"`, `"bottom"` or
#'   `"surface"`, as in [mask()][mask-3d]. `method` sets the depth level;
#'   `time_match` sets the time step.
#' @param time_match,max_time_diff How a point's time is matched to a time
#'   step, as in [as_voxel()]: `"nearest"` (default), the step nearest in
#'   time, or `"exact"`; and for nearest matching, the largest gap in days.
#'   `NULL` (default) sets no limit. These set the time step; `method` sets
#'   the depth level.
#' @param ... Not used.
#'
#' @returns A data frame with one row per point, in the order of `y`:
#'   * `feature_id`, the point's ID;
#'   * `value`, `NA` for a point off the grid, without a matching time step,
#'     or with no level its `method` can take;
#'   * `time`, the date-time of the time step used, `NA` for a source
#'     without dated steps or a point with no matching step. A point off the
#'     grid keeps its step's time;
#'   * `depth`, the depth of the level used in metres, `NA` for a variable
#'     without depth or when `value` is `NA`.
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
    # time, as as_voxel() matches a point.
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
        points_word <- if (n == 1L) "point" else "points"
        verb <- if (untimed == 1L) "has" else "have"
        rlang::warn(glue::glue(
          "{untimed} of {n} {points_word} {verb} no time, but variable ",
          "`{x@variable$name}` has {n_steps} time steps, so their values ",
          "are NA. Give the points a time in as_vector3d()."
        ), class = "ocean3d_warning_point_time")
      }
      resolution <- if (is.na(y@time_resolution)) "datetime" else
        y@time_resolution
      .match_proxy_steps(x, v$time, resolution, time_match, max_time_diff,
                         on_missing = "na")
    }

    # Half a cell each way: a group's crop then holds every cell a point on
    # the group's edge can fall in, as the cells a point's own crop holds,
    # so a point takes the same cell whatever other points share its step.
    half <- vapply(x@axes[c("lon", "lat")], function(a) {
      max(abs(a$upper - a$lower)) / 2
    }, numeric(1))
    # Longitudes in the grid's own frame (0 to 360, -180 to 180, or a
    # regional span): each moved by whole turns to the first value at or
    # east of the grid's western edge. A group's box is built from these,
    # so it never wraps round the globe onto a grid that doesn't (two points
    # far apart on a regional grid), and a cluster across the antimeridian
    # on a global grid is a wide box, which the split below separates.
    west <- min(x@axes$lon$lower, x@axes$lon$upper)
    lon_grid <- v$lon + 360 * ceiling((west - v$lon) / 360)
    # Groups of points read together: those on one time step, split into
    # nearby halves while a group's read would exceed the budget or hold
    # more than `cells_per_point` cells and levels for each point, so points
    # far apart aren't read with the water between them. A group whose
    # points all share one position can't be split, and its cells always
    # fit.
    cells_per_point <- 1e4
    value <- rep(NA_real_, n)
    level_depth <- rep(NA_real_, n)
    groups <- unname(split(seq_len(n)[!is.na(step)], step[!is.na(step)]))
    while (length(groups)) {
      idx <- groups[[1L]]
      groups <- groups[-1L]
      points <- y[idx]
      lon <- lon_grid[idx]
      lat <- v$lat[idx]
      box <- terra::ext(min(lon) - half[["lon"]], max(lon) + half[["lon"]],
                        min(lat) - half[["lat"]], max(lat) + half[["lat"]])
      cropped <- tryCatch(crop(x, box),
                          ocean3d_error_proxy_crop = function(e) NULL)
      if (is.null(cropped)) next # every point of the group is off the grid
      sizes <- cropped@variable$dim_lengths
      dims <- unlist(cropped@variable$axis_dims[c("lon", "lat", "depth")])
      cells <- prod(sizes[dims])
      too_big <- cells * 8 / 2^20 > budget || cells / length(idx) >
        cells_per_point
      if (too_big && nrow(unique(cbind(lon, lat))) > 1L) {
        # Split along the axis the crop spans more cells of.
        by_lon <- sizes[[cropped@variable$axis_dims$lon]] >=
          sizes[[cropped@variable$axis_dims$lat]]
        along <- if (by_lon) lon else lat
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
        at <- .point_cells(points, r)
        on <- which(!is.na(at))
        value[idx[on]] <- terra::extract(r, at[on])[, 1L]
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
