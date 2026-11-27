# -----------------------------------------------------------------------------
# Buffer extraction: a convenience wrapper over the three steps (#116)
# -----------------------------------------------------------------------------
#
# Each observation becomes a SpatVector3D point and, with a radius, a buffer
# from buffer_3d() (build); the data for its time steps are read with crop()
# and as_voxel(), or rast() for a variable without depth, and the values
# inside the buffer kept with mask() (intersect); summarise_3d() gives one
# statistic per observation (summarise). The wrapper calls only exported
# functions of those steps.

#' Summarise environmental values around observations
#'
#' Calculates one statistic per observation from the values in a buffer
#' around it: a circle of `horizontal_radius_km`, a depth window of
#' `depth_radius_m` above and below the observation, and a time window from
#' `time_before` to `time_after`. It is a convenience wrapper over the three
#' steps of an analysis, and the [workflow
#' article](../articles/workflow-values-in-shapes.html) shows the same steps
#' written out:
#'
#' 1. **Build:** each observation becomes a point with [as_vector3d()], and a
#'    buffer with [buffer_3d()] when `horizontal_radius_km` is positive.
#' 2. **Intersect:** the data for its time steps are read with [crop()] and
#'    [as_voxel()] (or [rast()][rast-SpatProxy] for a variable without
#'    depth), and the values inside the buffer are kept with [mask()].
#' 3. **Summarise:** [summarise_3d()] gives the statistic of those values.
#'
#' A cell or depth level counts in full when the buffer reaches it, as in
#' [mask()]: a cell the circle touches, and a level whose slab of water
#' overlaps the depth window. With `weighting = "weighted"` each value is
#' weighted by the volume of its cell and level (by the cell's area for a
#' variable without depth, or when `method` picks one level per cell), as
#' [summarise_3d()] does.
#'
#' Without a depth window, `method` picks the level in each cell: by default
#' the level nearest the observation's depth that has a value, as
#' [extract_to_point()] does. Without a radius, the buffer is the
#' observation's own cell. Without a time window, the time step is the one
#' `time_match` matches. A time window pools the values of every time step
#' inside it into one statistic, each value keeping its own weight.
#'
#' An observation outside the grid, without a date on data with several time
#' steps, without a depth where one is needed (or with a negative one), or
#' with no time step in its window or within `max_time_diff` gets `NA`, and
#' `<output>_n` is 0. One warning counts these observations, by reason. On
#' data with several time steps on one day, observations need date-times or a
#' time window: a date alone can't pick one of the steps, so the function
#' stops.
#'
#' @param data Observation table, matrix, named list or sf POINT object.
#' @param nc The data: a netCDF file path, several paths for a series of
#'   files on one grid, an open [ncdf4::nc_open()] connection, a
#'   [SpatProxy-class] from [as_proxy()], or a [SpatVoxel-class] in memory.
#' @param var Name of the variable, for a netCDF source. `NULL` detects the
#'   only data variable of the file.
#' @param horizontal_radius_km Radius of the circle around each observation,
#'   in kilometres, measured on the WGS84 ellipsoid, below 10000 as in
#'   [buffer_3d()]. `0` (default) uses the observation's own cell.
#' @param depth_radius_m Distance in metres above and below the observation's
#'   depth. `0` (default) picks one level per cell with `method`.
#' @param time_before,time_after Time before and after the observation, in
#'   `time_unit`. Both `0` (default) read the one time step `time_match`
#'   matches. For observations with dates and no hours, the window covers
#'   whole UTC days, including the last.
#' @param time_unit `"days"` (default), `"hours"`, `"minutes"` or
#'   `"seconds"`.
#' @param statistic One of `"mean"` (default), `"sd"`, `"min"` and `"max"`.
#' @param weighting `"weighted"` (default) weights each value by the volume of
#'   its cell and level with a depth window, and by its cell's area without
#'   one (one level per cell) or for a variable without depth; `"equal"`
#'   counts each value once.
#' @param method How the level in each cell is chosen when there is no depth
#'   window: one of the rules of [mask()] for points, `"nearest"` (default;
#'   the level nearest the observation's depth that has a value),
#'   `"nearest_level"` (the nearest level, even if it has no value),
#'   `"bottom"` or `"surface"`. Leave `NULL` with a depth window. For a
#'   variable without depth, leave `NULL` or use `"2d"`.
#' @param lon_col,lat_col,depth_col,date_col Observation columns. The depth is
#'   needed with a depth window or a `"nearest"` rule, and the date for data
#'   with several time steps. Dates and date-times are read as in
#'   [as_vector3d()].
#' @param lon_dim,lat_dim,depth_dim,time_dim Dimension names for a netCDF
#'   source, as in [as_proxy()].
#' @param time_match,max_time_diff Without a time window, how the time step is
#'   matched, as in [as_voxel()]: the nearest step (default) or an exact one,
#'   and the largest gap in days for a nearest match.
#' @param output_col Name of the output column. `NULL` uses the variable's
#'   name.
#'
#' @returns `data` with two columns added: the statistic, named `output_col`
#'   or after the variable, and `<output>_n`, the number of values it is
#'   built on. sf input keeps its geometry. If `data` already has either
#'   column, the function stops with an error; choose another `output_col`.
#'
#' @seealso [extract_to_point()], [extract_to_area()], and the [workflow
#'   article](../articles/workflow-values-in-shapes.html).
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
#' obs <- data.frame(lon = 5.2, lat = 5.3, depth = 20,
#'                   date = as.Date("2020-01-02"))
#'
#' # The mean within 150 km and 10 m above and below, on the observation's day.
#' extract_to_buffer(obs, path, horizontal_radius_km = 150,
#'                   depth_radius_m = 10)
#'
#' # Pooled over a day either side.
#' extract_to_buffer(obs, path, horizontal_radius_km = 150,
#'                   depth_radius_m = 10, time_before = 1, time_after = 1)
#' @export
extract_to_buffer <- function(
  data, nc, var = NULL,
  horizontal_radius_km = 0, depth_radius_m = 0,
  time_before = 0, time_after = 0,
  time_unit = c("days", "hours", "minutes", "seconds"),
  statistic = c("mean", "sd", "min", "max"),
  weighting = c("weighted", "equal"), method = NULL,
  lon_col = "lon", lat_col = "lat", depth_col = "depth", date_col = "date",
  lon_dim = NULL, lat_dim = NULL, depth_dim = NULL, time_dim = NULL,
  time_match = c("nearest", "exact"), max_time_diff = NULL,
  output_col = NULL
) {
  time_unit <- match.arg(time_unit)
  statistic <- match.arg(statistic)
  weighting <- match.arg(weighting)
  time_match <- match.arg(time_match)
  stop_buffer <- function(message) {
    rlang::abort(message, class = "ocean3d_error_buffer_input", call = NULL)
  }
  for (amount in list(horizontal_radius_km, depth_radius_m, time_before,
                      time_after)) {
    if (!is.numeric(amount) || length(amount) != 1L || !is.finite(amount) ||
          amount < 0) {
      stop_buffer(paste(
        "`horizontal_radius_km`, `depth_radius_m`, `time_before` and",
        "`time_after` must each be one number, 0 or more."
      ))
    }
  }
  windowed <- time_before > 0 || time_after > 0
  if (windowed && !is.null(max_time_diff)) {
    stop_buffer("`max_time_diff` applies without a time window only.")
  }

  proxy <- if (methods::is(nc, "SpatProxy")) {
    nc
  } else {
    as_proxy(nc, var = var, lon_dim = lon_dim, lat_dim = lat_dim,
             depth_dim = depth_dim, time_dim = time_dim)
  }
  has_depth <- !is.null(proxy@axes$depth)
  time_axis <- proxy@axes$time
  steps <- if (is.null(time_axis)) 1L else
    proxy@variable$dim_lengths[[time_axis$dim]]
  timed <- steps > 1L
  if (!has_depth) {
    if (!is.null(method) && !identical(method, "2d")) {
      stop_buffer(glue::glue(
        "Variable `{proxy@variable$name}` has no depth, so `method` doesn't ",
        "apply; leave it NULL."
      ))
    }
    method <- NULL
  } else if (depth_radius_m > 0) {
    if (!is.null(method)) {
      stop_buffer("Leave `method` NULL with a depth window.")
    }
  } else {
    if (is.null(method)) method <- "nearest"
    rules <- c("nearest", "nearest_level", "bottom", "surface")
    if (!rlang::is_string(method) || !method %in% rules) {
      stop_buffer(glue::glue(
        "`method` must be one of {toString(rules)}."
      ))
    }
  }
  needs_depth <- has_depth &&
    (depth_radius_m > 0 || method %in% c("nearest", "nearest_level"))

  # One plain table of the observations, whatever form `data` takes; sf
  # points give their coordinates in WGS84.
  table <- if (inherits(data, "sf")) {
    xy <- sf::st_coordinates(sf::st_transform(data, 4326))
    out_table <- sf::st_drop_geometry(data)
    out_table[[lon_col]] <- xy[, "X"]
    out_table[[lat_col]] <- xy[, "Y"]
    out_table
  } else {
    as.data.frame(data, stringsAsFactors = FALSE)
  }
  wanted <- c(lon_col, lat_col, if (needs_depth) depth_col,
              if (timed) date_col)
  missing_cols <- setdiff(wanted, names(table))
  if (length(missing_cols)) {
    stop_buffer(glue::glue(
      "`data` needs the column(s) {toString(missing_cols)}."
    ))
  }
  output <- if (is.null(output_col)) proxy@variable$name else output_col
  if (!rlang::is_string(output) || !nzchar(output)) {
    stop_buffer("`output_col` must be one non-empty name.")
  }
  added <- c(output, paste0(output, "_n"))
  # A matrix has no names() of its own, so its converted table is checked too.
  taken <- union(names(data), names(table))
  if (any(added %in% taken)) {
    clash <- toString(intersect(added, taken))
    stop_buffer(glue::glue(
      "`data` already has a column named {clash}. Choose another `output_col`."
    ))
  }

  work <- data.frame(lon = table[[lon_col]], lat = table[[lat_col]],
                     depth = if (needs_depth) table[[depth_col]] else 0)
  if (timed) work$time <- table[[date_col]]
  # A depth that can't be one (NA, negative or infinite) leaves its row out,
  # as a missing position or date does.
  usable <- is.finite(work$lon) & is.finite(work$lat) & abs(work$lat) <= 90 &
    is.finite(work$depth) & work$depth >= 0
  if (timed) usable <- usable & !is.na(work$time)
  value <- rep(NA_real_, nrow(work))
  n <- rep(0L, nrow(work))
  points <- if (any(usable)) {
    as_vector3d(work[usable, , drop = FALSE], time = if (timed) "time")
  }
  daily <- timed && !is.null(points) &&
    identical(points@time_resolution, "daily")
  seconds <- .time_unit_seconds(time_unit)
  if (daily && windowed) {
    days <- c(time_before, time_after) * seconds / 86400
    if (any(abs(days - round(days)) > 1e-6)) {
      stop_buffer(paste(
        "`data` has dates without hours, so `time_before` and `time_after`",
        "must be whole days. Give the times with hours to use part of a day."
      ))
    }
  }
  # A date can't pick one of several steps on its day: the caller has to give
  # date-times, as as_voxel() asks.
  if (daily && !windowed) {
    step_days <- floor(as.numeric(
      time_axis$cftime$as_timestamp(asPOSIX = TRUE)
    ) / 86400)
    if (anyDuplicated(step_days)) {
      stop_buffer(paste(
        "Several time steps fall on one day in the data, so a date can't",
        "pick one of them. Give the observations date-times (with hours), or",
        "a time window of whole days to pool them."
      ))
    }
  }
  # Every buffer at once, so buffer_3d() checks the radius before any read.
  circles <- if (horizontal_radius_km > 0 && !is.null(points)) {
    buffer_3d(points, radius_km = horizontal_radius_km,
              depth_window = depth_radius_m)
  }
  # The whole water column, for the cells a circle covers.
  deepest <- if (has_depth) {
    max(proxy@axes$depth$upper, proxy@axes$depth$centres) + 1
  } else {
    1
  }

  off_grid <- 0L
  no_step <- 0L
  for (k in seq_along(which(usable))) {
    row <- which(usable)[k]
    point <- points[k]
    vertex <- as.data.frame(point)
    circle <- if (!is.null(circles)) circles[k]
    at <- if (timed) vertex$time
    when <- if (!timed) {
      NULL
    } else if (windowed && daily) {
      as.Date(at, tz = "UTC") + c(-time_before, time_after)
    } else if (windowed) {
      at + c(-time_before, time_after) * seconds
    } else if (daily) {
      as.Date(at, tz = "UTC")
    } else {
      at
    }
    read <- function(source) {
      if (has_depth) {
        as_voxel(source, time = when, time_match = time_match,
                 max_time_diff = max_time_diff)
      } else {
        terra::rast(source, time = when, time_match = time_match,
                    max_time_diff = max_time_diff)
      }
    }
    # Cropped first: an error raised while S4 evaluates an argument for
    # dispatch loses its class. An observation off the grid, or with no time
    # step in its window or within `max_time_diff`, gets NA and is counted in
    # the warning below; the several-steps-a-day case stopped above.
    source <- tryCatch(crop(proxy, if (is.null(circle)) point else circle),
                       ocean3d_error_proxy_crop = function(e) NULL)
    if (is.null(source)) {
      off_grid <- off_grid + 1L
      next
    }
    layers <- tryCatch(read(source),
                       ocean3d_error_proxy_time = function(e) NULL)
    if (is.null(layers)) {
      no_step <- no_step + 1L
      next
    }
    if (!is.list(layers)) layers <- list(layers)

    # The outline of the buffer over the whole water column: the cells the
    # circle covers, whatever their levels.
    column <- if (!is.null(circle)) {
      outline <- sf::st_sf(geometry = as.data.frame(circle)$geometry)
      as_vector3d(outline, depth_min = 0, depth_max = deepest)
    }
    # The shapes keep the observation's time, which the voxels already
    # match, so mask()'s message that it drops polygon times is muted.
    kept <- suppressMessages(lapply(layers, function(v) {
      if (!has_depth) {
        return(mask(v, if (is.null(column)) point else column))
      }
      if (depth_radius_m > 0) {
        if (!is.null(circle)) return(mask(v, circle))
        # The observation's own cell, as a small square inside it, over the
        # depth window.
        grid <- methods::as(v, "SpatRaster")[[1L]] # remove the 3D class
        filled <- terra::setValues(terra::rast(grid), 1)
        cell <- which(!is.na(terra::values(mask(filled, point))[, 1]))
        if (!length(cell)) return(NULL)
        xy <- terra::xyFromCell(grid, cell)
        half <- terra::res(grid) / 4
        square <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
          xy + c(-1, -1) * half, xy + c(1, -1) * half, xy + c(1, 1) * half,
          xy + c(-1, 1) * half, xy + c(-1, -1) * half
        ))), crs = terra::crs(grid)))
        band <- as_vector3d(square,
                            depth_min = max(0, vertex$depth - depth_radius_m),
                            depth_max = vertex$depth + depth_radius_m)
        return(mask(v, band))
      }
      # One level per cell, picked with `method`: at the observation itself,
      # or at the centre of every cell the circle covers, then flattened to
      # one layer so each cell is weighted by its area.
      at_cells <- if (is.null(column)) {
        point
      } else {
        grid <- methods::as(v, "SpatRaster")[[1L]] # remove the 3D class
        filled <- terra::setValues(terra::rast(grid), 1)
        covered <- which(!is.na(terra::values(mask(filled, column))[, 1]))
        if (!length(covered)) return(NULL)
        xy <- terra::xyFromCell(grid, covered)
        if (!terra::is.lonlat(grid)) {
          xy <- terra::crds(terra::project(
            terra::vect(xy, crs = terra::crs(grid)), "EPSG:4326"
          ))
        }
        as_vector3d(data.frame(lon = xy[, 1L], lat = xy[, 2L],
                               depth = vertex$depth))
      }
      picked <- mask(v, at_cells, method = method)
      terra::app(methods::as(picked, "SpatRaster"), function(z) {
        if (all(is.na(z))) NA_real_ else z[!is.na(z)][1L]
      })
    }))
    kept <- Filter(Negate(is.null), kept)
    if (!length(kept)) next
    summary <- summarise_3d(kept, statistic = statistic, by = "all",
                            weighting = if (weighting == "weighted") "volume"
                            else "equal")
    value[row] <- summary[[statistic]]
    n[row] <- summary$n
  }

  # One warning for every observation that got NA without a value to read.
  reasons <- c(
    if (any(!usable)) glue::glue(
      "{sum(!usable)} without a usable position",
      "{if (needs_depth) ', depth' else ''}{if (timed) ' or date' else ''}"
    ),
    if (off_grid) glue::glue("{off_grid} outside the grid"),
    if (no_step) glue::glue(
      "{no_step} with no time step ",
      "{if (windowed) 'in the window' else 'that matches'}"
    )
  )
  if (length(reasons)) {
    dropped <- sum(!usable) + off_grid + no_step
    rlang::warn(glue::glue(
      "{dropped} of {nrow(work)} observations get NA: ",
      "{glue::glue_collapse(reasons, sep = ', ', last = ' and ')}."
    ), class = "ocean3d_warning_buffer_rows")
  }

  out <- if (inherits(data, "sf") || is.data.frame(data)) data else table
  out[[output]] <- value
  out[[paste0(output, "_n")]] <- n
  out
}
