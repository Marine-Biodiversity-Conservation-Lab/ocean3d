#-------------------------------------------------------------------------------
# Point-based extraction from netCDF files
#-------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------
# These functions are intentionally small and focused. This makes the code
# easier to test and easier to maintain inside a package.
# -----------------------------------------------------------------------------

.check_required_cols <- function(x, cols, x_name = "data") {
  cols <- stats::na.omit(cols)
  missing_cols <- setdiff(cols, names(x))

  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required column(s) in `%s`: %s",
        x_name,
        toString(missing_cols)
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.normalize_points_input <- function(x, lon_col = "lon", lat_col = "lat") {
  # data.frame subclasses such as tibbles already use the internal format.
  if (is.data.frame(x) && !inherits(x, "sf")) {
    return(x)
  }

  # sf objects can provide longitude and latitude directly from POINT geometry.
  if (inherits(x, "sf")) {
    geom_type <- unique(as.character(sf::st_geometry_type(x)))

    if (length(geom_type) != 1L || geom_type != "POINT") {
      stop("`sf` input must contain POINT geometries.", call. = FALSE)
    }

    coords_source <- x

    # Any CRS other than WGS84, geographic or not, is transformed.
    if (sf::st_crs(coords_source) != sf::st_crs(4326)) {
      if (is.na(sf::st_crs(coords_source))) {
        stop(
          "`sf` input must have a geographic CRS or a CRS that can be ",
          "transformed to longitude/latitude.",
          call. = FALSE
        )
      }

      coords_source <- sf::st_transform(coords_source, 4326)
    }

    coords <- sf::st_coordinates(coords_source)
    x[[lon_col]] <- coords[, 1]
    x[[lat_col]] <- coords[, 2]

    return(x)
  }

  # Matrices are converted to a data frame while preserving column names.
  if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop("Matrix input must have column names.", call. = FALSE)
    }

    return(as.data.frame(x))
  }

  # Named lists are converted to a data frame.
  if (is.list(x)) {
    if (is.null(names(x)) || any(names(x) == "")) {
      stop("List input must be named.", call. = FALSE)
    }

    return(as.data.frame(x))
  }

  stop(
    "`data` must be a data frame, tibble, sf object, matrix, or named list.",
    call. = FALSE
  )
}


# Select layers by physical depth, independently of storage order. Returns
# the selected values, the physical depth (m) of each selected layer and its
# index in `profile`. Shared by extract_to_point() and mask() of a voxel by
# points.
.summarise_depth_profile <- function(profile, depth_m, target_depth) {
  index <- c(nearest = NA_integer_, surface = which.min(depth_m),
             seabottom = NA_integer_)
  valid <- which(!is.na(profile))
  if (length(valid)) {
    index["seabottom"] <- valid[which.max(depth_m[valid])]
    if (!is.na(target_depth)) {
      # A tie chooses the shallower valid layer.
      valid <- valid[order(depth_m[valid])]
      index["nearest"] <- valid[which.min(abs(depth_m[valid] - target_depth))]
    }
  }
  list(
    value = stats::setNames(profile[index], names(index)),
    depth_m = stats::setNames(depth_m[index], names(index)),
    index = index
  )
}


# -----------------------------------------------------------------------------
# Main extraction function
# -----------------------------------------------------------------------------


#' Extract values from netCDF files at point observations
#'
#' A convenience wrapper over the build and intersect steps: it turns the
#' observations into [SpatVector3D-class] points with [as_vector3d()],
#' describes the netCDF data with [as_proxy()], and reads the value at every
#' point with [extract()][extract-SpatProxy], which reads each time step once
#' for all the points in it. The same can be done step by step; see the
#' workflow article, `vignette("workflow-values-in-shapes")`.
#'
#' You choose the depth level with `method`, which follows the rules of the
#' same name in `extract()`: `"nearest"` takes the level nearest the
#' observation's depth, skipping levels without a value; `"surface"` takes the
#' shallowest level, even if it has no value; `"bottom"` takes the deepest
#' level with a value. `"2d"` is for a variable without depth. Each observation
#' takes the cell that contains it. Times are matched to the nearest time step
#' by default (see `time_match` and `max_time_diff`): a date by its UTC day, a
#' date-time as an instant.
#'
#' @param data Optional observations: a data frame, tibble, `sf` object with
#'   POINT geometries, numeric matrix with column names, or named list. If
#'   `NULL`, give the coordinates directly with `lon`, `lat`, and `depth` and
#'   `date` as needed. `sf` points need a CRS and are transformed to longitude
#'   and latitude.
#' @param nc The netCDF data: a file path, several paths for a series of files
#'   on one grid (each with its own time steps, such as one file per day), a
#'   connection opened with [ncdf4::nc_open()], or a `SpatProxy` from
#'   [as_proxy()]. Paths and connections are read through [as_proxy()]; a
#'   `SpatProxy` is used as it is, so `var` and the dimension arguments don't
#'   apply to it.
#' @param lon,lat,depth,date Direct input when `data = NULL`: numeric
#'   longitudes, latitudes and depths (metres, positive down), and dates or
#'   date-times. Each has length 1 or the same length as the others.
#'   `depth` is needed for `method = "nearest"`.
#' @param var Name of the variable. `NULL` detects it when the file holds one
#'   data variable.
#' @param method `"2d"` (default) for a variable without depth or with a
#'   single depth level, or `"surface"`, `"bottom"` or `"nearest"` for one
#'   with depth levels (see Details).
#' @param lon_col,lat_col,date_col,depth_col Columns of `data` with the
#'   longitude, latitude, date and depth. `date_col` is needed when the
#'   variable has dated time steps (not a climatology); `depth_col` for
#'   `method = "nearest"`. The coordinate columns must be numeric.
#' @param id_col Optional column of observation identifiers. It is checked to
#'   exist; extraction doesn't use it.
#' @param lon_dim,lat_dim,depth_dim,time_dim Optional dimension names, passed
#'   to [as_proxy()].
#' @param time_match `"nearest"` (default) or `"exact"`: how observation times
#'   are matched to time steps, as in [as_voxel()].
#' @param max_time_diff Optional largest gap, in days, between an observation
#'   and its nearest time step with `time_match = "nearest"`; `NULL` (default)
#'   allows any gap.
#' @param output_col,output_prefix Name of the output column. If both are
#'   given, `output_col` is used. By default, the variable's name, or with
#'   `var` detected, `surface_<var>`, `seabottom_<var>` or `nearest_<var>` for
#'   those methods.
#'
#' @return With `data`, the observations with the value column added. Data
#'   frames, tibbles and `sf` objects keep their class; matrices and named
#'   lists become data frames. With direct input, a numeric vector, one value
#'   per point. A value is `NA` when the observation lies off the grid, has
#'   no coordinates, has no time step within `max_time_diff`, or has no
#'   level its `method` can take.
#'
#' @seealso [extract()][extract-SpatProxy] for the step this wraps;
#'   [extract_to_buffer()] and [extract_to_track()].
#' @examples
#' nc_file <- system.file("extdata", "example_3d.nc", package = "ocean3d")
#'
#' # One observation: the nearest level with a value, on 3 January.
#' extract_to_point(nc = nc_file, lon = 0, lat = 40, depth = 60,
#'                  date = "2020-01-03", var = "temp", method = "nearest")
#'
#' # A table of observations keeps its columns and gains one.
#' obs <- data.frame(id = c("a", "b"), lon = c(0, 1), lat = c(40, 41),
#'                   depth = c(60, 100),
#'                   date = as.Date(c("2020-01-03", "2020-01-03")))
#' extract_to_point(obs, nc_file, var = "temp", method = "bottom")
#' @export
extract_to_point <- function(data = NULL, nc, lon = NULL, lat = NULL,
                             depth = NULL, date = NULL, var = NULL,
                             method = c("2d", "surface", "bottom", "nearest"),
                             lon_col = "lon", lat_col = "lat",
                             date_col = "date", depth_col = "depth",
                             id_col = NULL, lon_dim = NULL, lat_dim = NULL,
                             depth_dim = NULL, time_dim = NULL,
                             time_match = c("nearest", "exact"),
                             max_time_diff = NULL, output_prefix = NULL,
                             output_col = NULL) {
  method <- match.arg(method)
  time_match <- match.arg(time_match)
  stop_input <- function(message) {
    rlang::abort(message, class = "ocean3d_error_point_input", call = NULL)
  }
  if (is.data.frame(nc)) {
    stop_input(paste(
      "`nc` as a table of files is no longer supported. Give one path,",
      "several paths of a series on one grid, a connection or a SpatProxy;",
      "for files on different grids, call once per grid."
    ))
  }
  if (is.list(nc) && !inherits(nc, "ncdf4")) nc <- unlist(nc)
  proxy <- if (methods::is(nc, "SpatProxy")) {
    nc
  } else {
    as_proxy(nc, var = var, lon_dim = lon_dim, lat_dim = lat_dim,
             depth_dim = depth_dim, time_dim = time_dim)
  }
  var_name <- proxy@variable$name
  has_depth <- !is.null(proxy@axes$depth)
  # A single depth level, as in many surface products, reads as 2D.
  if (method == "2d" && length(proxy@axes$depth$centres) > 1L) {
    stop_input(glue::glue(
      "Variable `{var_name}` has several depth levels, so `method = \"2d\"` ",
      "doesn't apply. Use \"surface\", \"bottom\" or \"nearest\"."
    ))
  }
  if (method != "2d" && !has_depth) {
    stop_input(glue::glue(
      "Variable `{var_name}` has no depth axis. Use `method = \"2d\"`."
    ))
  }
  # Dated time steps, even one, are matched to each observation's date; a
  # climatology gives its steps to every observation.
  time_axis <- proxy@axes$time
  timed <- !is.null(time_axis) && !isTRUE(time_axis$climatology)

  # The observations as a plain table of coordinates, depth and time.
  direct <- is.null(data)
  if (direct) {
    inputs <- Filter(Negate(is.null),
                     list(lon = lon, lat = lat, depth = depth, date = date))
    n <- max(lengths(inputs))
    if (!all(lengths(inputs) %in% c(1L, n))) {
      stop_input(paste(
        "`lon`, `lat`, `depth` and `date` must each have length 1 or the",
        "same length."
      ))
    }
    table <- data.frame(lapply(inputs, rep_len, length.out = n))
    lon_col <- "lon"
    lat_col <- "lat"
    depth_col <- "depth"
    date_col <- if (!is.null(date)) "date"
  } else if (inherits(data, "sf")) {
    if (is.na(sf::st_crs(data))) {
      stop_input(paste(
        "`data` has no CRS, so its points can't be placed on the grid. Set",
        "it with `sf::st_set_crs()`, for example to 4326 for longitude and",
        "latitude."
      ))
    }
    xy <- sf::st_coordinates(sf::st_transform(sf::st_geometry(data), 4326))
    table <- sf::st_drop_geometry(data)
    table[[lon_col]] <- xy[, "X"]
    table[[lat_col]] <- xy[, "Y"]
  } else {
    table <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  n <- nrow(table)
  if (timed && is.null(date_col)) {
    stop_input(glue::glue(
      "Variable `{var_name}` has dated time steps, so the observations need ",
      "a date."
    ))
  }
  needed <- c(lon_col, lat_col, if (method == "nearest") depth_col,
              if (timed) date_col, id_col)
  missing_cols <- setdiff(needed, names(table))
  if (length(missing_cols)) {
    stop_input(glue::glue("`data` has no column {toString(missing_cols)}."))
  }
  # Coordinates held as text, as in a matrix with a date column or a column
  # read from a file as text, would otherwise all give NA.
  coordinate_cols <- c(lon_col, lat_col, if (method == "nearest") depth_col)
  text_cols <- coordinate_cols[
    !vapply(table[coordinate_cols], is.numeric, logical(1))
  ]
  if (length(text_cols)) {
    columns_word <- if (length(text_cols) == 1L) "Column" else "Columns"
    stop_input(glue::glue(
      "{columns_word} {toString(paste0('`', text_cols, '`'))} must be ",
      "numeric. A matrix with a date column holds its numbers as text; give ",
      "a data frame instead."
    ))
  }

  output <- if (!is.null(output_col)) {
    output_col
  } else if (!is.null(output_prefix)) {
    output_prefix
  } else if (is.null(var) && method != "2d") {
    paste0(c(surface = "surface_", bottom = "seabottom_",
             nearest = "nearest_")[[method]], var_name)
  } else {
    var_name
  }
  if (!direct && output %in% names(table)) {
    stop_input(glue::glue(
      "`data` already has a column `{output}`; choose another name with ",
      "`output_col`."
    ))
  }

  # Observations without coordinates, or without the depth `nearest` needs,
  # get NA; the rest become points. Depth only matters for `nearest`.
  point_depth <- rep_len(if (method == "nearest") table[[depth_col]] else 0, n)
  usable <- is.finite(table[[lon_col]]) & is.finite(table[[lat_col]]) &
    is.finite(point_depth)
  values <- rep(NA_real_, n)
  if (any(usable)) {
    points <- data.frame(lon = table[[lon_col]][usable],
                         lat = table[[lat_col]][usable],
                         depth = point_depth[usable])
    # The date column keeps its own name, so an unreadable date is reported
    # by the name the caller gave it.
    time_col <- if (timed) {
      if (date_col %in% names(points)) "time" else date_col
    }
    if (timed) points[[time_col]] <- table[[date_col]][usable]
    points <- as_vector3d(points, time = time_col)
    rule <- if (method != "2d") {
      list(method = method)
    } else if (has_depth) {
      list(method = "surface")
    }
    found <- do.call(terra::extract, c(list(proxy, points), rule, list(
      time_match = time_match, max_time_diff = max_time_diff
    )))
    values[usable] <- found$value
  }
  if (direct) return(values)
  data <- if (inherits(data, c("data.frame", "sf"))) data else table
  data[[output]] <- values
  data
}
