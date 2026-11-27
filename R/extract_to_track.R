# -----------------------------------------------------------------------------
# Track extraction: a convenience wrapper over the three steps (#116)
# -----------------------------------------------------------------------------
#
# The vertex table becomes SpatVector3D lines, sampled into points along each
# track with sample_3d(), and with a radius buffered with buffer_3d()
# (build); each sample's value is read with extract(), or inside its buffer
# with crop(), as_voxel() (rast() for a variable without depth) and mask()
# (intersect); summarise_3d() gives a buffer's mean, and the samples' values
# are summarised per track, weighted by the distance or time each sample
# stands for (summarise). The wrapper calls only exported functions of those
# steps.

#' Summarise environmental values along tracks
#'
#' Samples points along each track, such as a dive or a vessel route, reads
#' the value at each sample, and summarises the samples of each track. It is
#' a convenience wrapper over the three steps of an analysis, and the
#' [workflow article](../articles/workflow-values-in-shapes.html) shows the
#' same steps written out:
#'
#' 1. **Build:** the vertices of each track become a line with
#'    [as_vector3d()], and [sample_3d()] places points along it, with depth
#'    and time interpolated between vertices. With `horizontal_radius_km`,
#'    [buffer_3d()] turns each sample into a buffer.
#' 2. **Intersect:** without a radius, [extract()][extract-SpatProxy] reads
#'    the value at each sample, in the cell that holds it and at the depth
#'    level `method` picks. With a radius, the data for the sample's time
#'    step are read with [crop()] and [as_voxel()] (or [rast()][rast-SpatProxy]
#'    for a variable without depth), and [mask()] keeps the values inside its
#'    buffer: the cells the circle touches, and the levels whose slab of
#'    water overlaps the depth window. [summarise_3d()] gives their mean,
#'    weighted by cell volume.
#' 3. **Summarise:** the statistics of each track's sample values, each
#'    sample weighted by `weighting`.
#'
#' Each sample is matched to the time step nearest its time, or equal to it
#' with `time_match = "exact"`, as [extract()][extract-SpatProxy] matches
#' points. A sample gets `NA` when it is outside the grid, on a cell or level
#' without a value, or has no time step that matches. One warning counts
#' these samples, by reason. The track summaries leave them out, and `n` is
#' the number of samples with a value.
#'
#' With `weighting = "distance"` each sample stands for half the distance to
#' the samples either side (`support_km`), so a stretch sampled more densely
#' doesn't count for more; with `"duration"` for half the time
#' (`support_seconds`). The weighted standard deviation is the one
#' [summarise_3d()] uses. A track whose samples all have zero weight, such
#' as a track that doesn't move with distance weights, gets equal weights.
#' Neighbouring samples are usually correlated, so the standard deviation
#' describes their spread, not the uncertainty of the mean.
#'
#' @param data Vertices of the tracks, in the order they were travelled: a
#'   data frame, matrix, named list or sf POINT object. sf points are
#'   transformed to longitude and latitude. To use sf lines, convert them
#'   first, e.g. with
#'   `as.data.frame(sf::st_coordinates(sf::st_transform(lines, 4326)))` and
#'   `lon_col = "X", lat_col = "Y", track_col = "L1"`.
#' @param nc The data: a netCDF file path, several paths for a series of
#'   files on one grid, an open [ncdf4::nc_open()] connection, or a
#'   [SpatProxy-class] from [as_proxy()].
#' @param var Name of the variable, for a netCDF source. `NULL` detects the
#'   only data variable of the file.
#' @param track_col Column with track IDs. `NULL` (default) reads all of
#'   `data` as one track, with ID `"1"`. Each track needs at least two
#'   vertices.
#' @param sampling,spacing_km,depth_step_m,time_step,time_unit,max_samples
#'   Where the samples go along each track, as in [sample_3d()]: only the
#'   vertices (default without a step), or with a spacing in kilometres or
#'   a depth or time step, points in between as well.
#' @param horizontal_radius_km Radius of a circle around each sample, in
#'   kilometres, below 10000 as in [buffer_3d()]. `0` (default) reads the
#'   value at the sample itself.
#' @param depth_radius_m Distance in metres above and below each sample's
#'   depth, for the buffers. Needed with a radius for a variable with depth,
#'   and used only with a radius.
#' @param method For a variable with depth and no radius, the depth level
#'   each sample takes, as in [extract()][extract-SpatProxy]: `"nearest"`
#'   (default; the level nearest the sample's depth that has a value),
#'   `"nearest_level"` (the nearest level, even if it has no value),
#'   `"bottom"` (the deepest level with a value) or `"surface"`. Leave
#'   `NULL` with a radius; for a variable without depth, leave `NULL` or use
#'   `"2d"`.
#' @param statistic One or more of `"mean"` (default), `"sd"`, `"min"` and
#'   `"max"`, of the samples of each track.
#' @param weighting How each sample is weighted in the track summaries:
#'   `"distance"` (default) by `support_km`, `"duration"` by
#'   `support_seconds` (needs dates), or `"equal"`.
#' @param lon_col,lat_col,depth_col,date_col Vertex columns. The depth is
#'   needed for a variable with depth, unless `method` is `"bottom"` or
#'   `"surface"`, and for `depth_step_m`. The date is needed for data with
#'   dated time steps, for `time_step` and for duration weights. Dates and
#'   date-times are read as in [as_vector3d()].
#' @param lon_dim,lat_dim,depth_dim,time_dim Dimension names for a netCDF
#'   source, as in [as_proxy()].
#' @param time_match,max_time_diff How a sample's time is matched to a time
#'   step, as in [as_voxel()]: the nearest step (default) or an exact one,
#'   and the largest gap in days for a nearest match.
#' @param output_col Name of the value column in `samples`. `NULL` uses the
#'   variable's name.
#'
#' @returns A list of two data frames:
#'   * `samples`: one row per sample, in order along each track, with
#'     `track_id`, `sample_id`, the position (in `lon_col` and `lat_col`),
#'     the depth and date when used (in `depth_col` and `date_col`),
#'     `distance_km` from the track's start, `support_km`,
#'     `support_seconds`, `is_interpolated` (`FALSE` for the track's own
#'     vertices), and the value in `output_col`. With a radius, also
#'     `<output>_n`, the number of values in the sample's buffer.
#'   * `summary`: one row per track, with `track_id`, a column per
#'     statistic, and `n`, the number of samples with a value.
#'
#' @seealso [sample_3d()], [extract_to_point()], [extract_to_buffer()], and
#'   the [workflow article](../articles/workflow-values-in-shapes.html).
#' @examples
#' # Daily temperature on a 10 x 10 grid of 1 degree cells at four depths.
#' path <- tempfile(fileext = ".nc")
#' lon <- ncdf4::ncdim_def("lon", "degrees_east", seq(0.5, 9.5))
#' lat <- ncdf4::ncdim_def("lat", "degrees_north", seq(0.5, 9.5))
#' depth <- ncdf4::ncdim_def("depth", "m", c(5, 15, 30, 60))
#' time <- ncdf4::ncdim_def("time", "days since 2020-01-01", 0:2)
#' temp <- ncdf4::ncvar_def("temp", "degree_Celsius",
#'                          list(lon, lat, depth, time), -9999)
#' nc <- ncdf4::nc_create(path, temp)
#' ncdf4::ncvar_put(nc, temp, seq_len(10 * 10 * 4 * 3) / 10)
#' ncdf4::nc_close(nc)
#'
#' # A dive from the surface to 50 m over two days.
#' dive <- data.frame(lon = c(2.5, 6.5), lat = c(3.5, 5.5), depth = c(0, 50),
#'                    date = as.Date(c("2020-01-01", "2020-01-03")))
#'
#' # The value at a sample at least every 50 km, and the mean along the dive.
#' track <- extract_to_track(dive, path, spacing_km = 50)
#' track$summary
#'
#' # The mean within 30 km and 10 m above and below each sample instead.
#' extract_to_track(dive, path, spacing_km = 50, horizontal_radius_km = 30,
#'                  depth_radius_m = 10)$summary
#' @export
extract_to_track <- function(
  data, nc, var = NULL, track_col = NULL,
  sampling = NULL, spacing_km = NULL, depth_step_m = NULL, time_step = NULL,
  time_unit = c("hours", "minutes", "seconds", "days", "months", "years"),
  max_samples = 100000L,
  horizontal_radius_km = 0, depth_radius_m = 0, method = NULL,
  statistic = "mean", weighting = c("distance", "duration", "equal"),
  lon_col = "lon", lat_col = "lat", depth_col = "depth", date_col = "date",
  lon_dim = NULL, lat_dim = NULL, depth_dim = NULL, time_dim = NULL,
  time_match = c("nearest", "exact"), max_time_diff = NULL,
  output_col = NULL
) {
  time_unit <- match.arg(time_unit)
  statistic <- match.arg(statistic, c("mean", "sd", "min", "max"),
                         several.ok = TRUE)
  weighting <- match.arg(weighting)
  time_match <- match.arg(time_match)
  stop_track <- function(message) {
    rlang::abort(message, class = "ocean3d_error_track_input", call = NULL)
  }
  for (amount in list(horizontal_radius_km, depth_radius_m)) {
    if (!is.numeric(amount) || length(amount) != 1L || !is.finite(amount) ||
          amount < 0) {
      stop_track(paste(
        "`horizontal_radius_km` and `depth_radius_m` must each be one",
        "number, 0 or more."
      ))
    }
  }
  buffered <- horizontal_radius_km > 0
  if (!buffered && depth_radius_m > 0) {
    stop_track("`depth_radius_m` needs a `horizontal_radius_km` above 0.")
  }

  proxy <- if (methods::is(nc, "SpatProxy")) {
    nc
  } else {
    as_proxy(nc, var = var, lon_dim = lon_dim, lat_dim = lat_dim,
             depth_dim = depth_dim, time_dim = time_dim)
  }
  var_name <- proxy@variable$name
  has_depth <- !is.null(proxy@axes$depth)
  if (!has_depth) {
    if (!is.null(method) && !identical(method, "2d")) {
      stop_track(glue::glue(
        "Variable `{var_name}` has no depth, so `method` doesn't apply; ",
        "leave it NULL."
      ))
    }
    method <- NULL
  } else if (buffered) {
    if (!is.null(method)) stop_track("Leave `method` NULL with a radius.")
    if (depth_radius_m == 0) {
      stop_track(glue::glue(
        "Variable `{var_name}` has depth, so buffers need a ",
        "`depth_radius_m` above 0."
      ))
    }
  } else {
    if (is.null(method)) method <- "nearest"
    rules <- c("nearest", "nearest_level", "bottom", "surface")
    if (!rlang::is_string(method) || !method %in% rules) {
      stop_track(glue::glue("`method` must be one of {toString(rules)}."))
    }
  }
  # Dated time steps, even one, are matched to each sample's time; a
  # climatology gives its steps to every sample.
  time_axis <- proxy@axes$time
  timed <- !is.null(time_axis) && !isTRUE(time_axis$climatology)
  needs_depth <- (has_depth && (buffered ||
                                  method %in% c("nearest", "nearest_level"))) ||
    !is.null(depth_step_m)
  needs_time <- timed || weighting == "duration" || !is.null(time_step)

  # One plain table of the vertices, whatever form `data` takes; sf points
  # give their coordinates in WGS84.
  table <- if (inherits(data, "sf")) {
    xy <- sf::st_coordinates(sf::st_transform(sf::st_geometry(data), 4326))
    out_table <- sf::st_drop_geometry(data)
    out_table[[lon_col]] <- xy[, "X"]
    out_table[[lat_col]] <- xy[, "Y"]
    out_table
  } else {
    as.data.frame(data, stringsAsFactors = FALSE)
  }
  wanted <- c(lon_col, lat_col, track_col, if (needs_depth) depth_col,
              if (needs_time) date_col)
  missing_cols <- setdiff(wanted, names(table))
  if (length(missing_cols)) {
    stop_track(glue::glue(
      "`data` needs the column(s) {toString(missing_cols)}."
    ))
  }
  output <- if (is.null(output_col)) var_name else output_col
  if (!rlang::is_string(output) || !nzchar(output)) {
    stop_track("`output_col` must be one non-empty name.")
  }

  # Build: each track becomes a line, sampled into points.
  work <- data.frame(
    track = if (is.null(track_col)) "1" else as.character(table[[track_col]]),
    lon = table[[lon_col]], lat = table[[lat_col]],
    depth = if (needs_depth) table[[depth_col]] else 0,
    stringsAsFactors = FALSE
  )
  if (needs_time) work$time <- table[[date_col]]
  lines <- as_vector3d(work, type = "lines", id = "track",
                       time = if (needs_time) "time")
  samples <- sample_3d(lines, sampling = sampling, spacing_km = spacing_km,
                       depth_step_m = depth_step_m, time_step = time_step,
                       time_unit = time_unit, max_samples = max_samples)
  out <- as.data.frame(samples)
  if (weighting == "duration" && anyNA(out$support_seconds)) {
    stop_track(paste(
      "`weighting = \"duration\"` needs a date at every vertex. Give the",
      "dates, or use \"distance\" or \"equal\"."
    ))
  }

  # Intersect: one value per sample. The time step each sample takes comes
  # from extract(), which matches all the samples at once.
  rule <- if (!is.null(method)) list(method = method)
  found <- do.call(terra::extract, c(list(proxy, samples), rule, list(
    time_match = time_match, max_time_diff = max_time_diff
  )))
  no_step <- if (timed) is.na(found$time_step) else rep(FALSE, nrow(out))
  value <- found$value
  if (buffered) {
    # For a variable without depth the depth window doesn't matter, but
    # mask() keeps no cell of a plain raster for a buffer with no thickness,
    # so it gets one metre.
    circles <- buffer_3d(samples, radius_km = horizontal_radius_km,
                         depth_window = if (has_depth) depth_radius_m else 1)
    value <- rep(NA_real_, nrow(out))
    n <- rep(0L, nrow(out))
    for (i in which(!no_step)) {
      circle <- circles[i]
      # Cropped first: an error raised while S4 evaluates an argument for
      # dispatch loses its class. A buffer off the grid gets NA.
      source <- tryCatch(crop(proxy, circle),
                         ocean3d_error_proxy_crop = function(e) NULL)
      if (is.null(source)) next
      when <- if (timed) found$time_step[i]
      layers <- if (has_depth) {
        as_voxel(source, time = when, time_match = "exact")
      } else {
        terra::rast(source, time = when, time_match = "exact")
      }
      # The buffers keep the sample's time, which the voxel already matches,
      # so mask()'s message that it drops polygon times is muted.
      kept <- suppressMessages(mask(layers, circle))
      summary <- summarise_3d(kept, statistic = "mean")
      value[i] <- summary$mean
      n[i] <- summary$n
    }
  }

  # One warning for every sample that got NA.
  missing_value <- is.na(value) & !no_step
  reasons <- c(
    if (any(no_step)) glue::glue(
      "{sum(no_step)} with no time step that matches"
    ),
    if (any(missing_value)) glue::glue(
      "{sum(missing_value)} outside the grid or without a value"
    )
  )
  if (length(reasons)) {
    rlang::warn(glue::glue(
      "{sum(is.na(value))} of {nrow(out)} samples get NA: ",
      "{glue::glue_collapse(reasons, sep = ' and ')}."
    ), class = "ocean3d_warning_track_samples")
  }

  # Summarise: the statistics of each track's samples, weighted by the
  # distance or time each stands for, as summarise_3d() weights values.
  weights <- switch(weighting, distance = out$support_km,
                    duration = out$support_seconds, equal = rep(1, nrow(out)))
  ids <- unique(out$track_id)
  summary <- do.call(rbind, lapply(ids, function(id) {
    keep <- out$track_id == id & !is.na(value)
    v <- value[keep]
    w <- weights[keep]
    if (sum(w) <= 0) w <- rep(1, length(v))
    w <- w / sum(w)
    average <- sum(w * v)
    correction <- 1 - sum(w^2)
    stats <- list(
      mean = if (length(v)) average else NA_real_,
      sd = if (length(v) > 1L && correction > 0) {
        sqrt(sum(w * (v - average)^2) / correction)
      } else {
        NA_real_
      },
      min = if (length(v)) min(v) else NA_real_,
      max = if (length(v)) max(v) else NA_real_
    )
    data.frame(track_id = id, stats[statistic], n = length(v),
               stringsAsFactors = FALSE)
  }))
  rownames(summary) <- NULL

  # The samples, track and sample number first, with the caller's column
  # names for position, depth and date, and the columns not used dropped.
  first <- c("track_id", "sample_id")
  out <- out[c(first, setdiff(names(out), c(first, "feature_id", "vertex",
                                            if (!needs_depth) "depth",
                                            if (!needs_time) "time")))]
  renamed <- c(lon = lon_col, lat = lat_col, depth = depth_col,
               time = date_col)
  at <- match(names(renamed), names(out))
  names(out)[at[!is.na(at)]] <- renamed[!is.na(at)]
  added <- c(output, if (buffered) paste0(output, "_n"))
  if (any(added %in% names(out))) {
    stop_track(glue::glue(
      "The samples already have a column named ",
      "{toString(intersect(added, names(out)))}. Choose another `output_col`."
    ))
  }
  out[[output]] <- value
  if (buffered) out[[paste0(output, "_n")]] <- n
  list(samples = out, summary = summary)
}
