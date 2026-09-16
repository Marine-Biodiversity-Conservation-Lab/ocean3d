#-------------------------------------------------------------------------------
# Point-based extraction from netCDF files
#-------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------
# These functions are intentionally small and focused. This makes the code easier
# to test and easier to maintain inside a package.
# -----------------------------------------------------------------------------

.check_required_cols <- function(x, cols, x_name = "data") {
  cols <- stats::na.omit(cols)
  missing_cols <- setdiff(cols, names(x))

  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required column(s) in `%s`: %s",
        x_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.is_ncdf4 <- function(x) {
  inherits(x, "ncdf4")
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

    if (!isTRUE(sf::st_is_longlat(coords_source))) {
      if (is.na(sf::st_crs(coords_source))) {
        stop(
          "`sf` input must have a geographic CRS or a CRS that can be transformed to longitude/latitude.",
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


.recycle_point_inputs <- function(lon, lat, depth = NULL, date = NULL) {
  inputs <- list(lon = lon, lat = lat, depth = depth, date = date)
  lengths <- vapply(inputs[!vapply(inputs, is.null, logical(1))], length, integer(1))
  n <- max(lengths)

  incompatible <- lengths != 1L & lengths != n
  if (any(incompatible)) {
    stop(
      "Direct coordinate inputs must have compatible lengths: each must have length 1 or the same maximum length.",
      call. = FALSE
    )
  }

  recycle <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 1L && n > 1L) rep(x, n) else x
  }

  lapply(inputs, recycle)
}


.flatten_list <- function(x) {
  # Users often pass nested lists when they group files by variable or by year.
  # This helper flattens nested lists while preserving names where possible.
  out <- list()

  add_one <- function(z, nm = NULL) {
    if (is.list(z) && !.is_ncdf4(z) && !is.data.frame(z)) {
      z_names <- names(z)
      for (i in seq_along(z)) {
        child_nm <- if (!is.null(z_names) && nzchar(z_names[i])) z_names[i] else nm
        add_one(z[[i]], child_nm)
      }
    } else {
      out[[length(out) + 1L]] <<- list(object = z, name = nm)
    }
  }

  x_names <- names(x)
  for (i in seq_along(x)) {
    nm <- if (!is.null(x_names) && nzchar(x_names[i])) x_names[i] else NULL
    add_one(x[[i]], nm)
  }

  out
}

.as_netcdf_sources <- function(nc, file_col = "file") {
  # This helper standardises all accepted input formats into the same internal
  # representation. Each source keeps the original object plus a type flag telling
  # the extraction code whether the file needs to be opened and closed.

  if (is.character(nc)) {
    sources <- lapply(nc, function(x) list(source = x, type = "path", name = basename(x)))
    return(sources)
  }

  if (.is_ncdf4(nc)) {
    return(list(list(source = nc, type = "connection", name = NA_character_)))
  }

  if (is.data.frame(nc)) {
    .check_required_cols(nc, file_col, "nc")
    files <- as.character(nc[[file_col]])
    sources <- lapply(files, function(x) list(source = x, type = "path", name = basename(x)))
    return(sources)
  }

  if (is.list(nc)) {
    flat <- .flatten_list(nc)
    sources <- vector("list", length(flat))

    for (i in seq_along(flat)) {
      obj <- flat[[i]]$object
      nm <- flat[[i]]$name

      if (is.character(obj) && length(obj) == 1L) {
        sources[[i]] <- list(
          source = obj,
          type = "path",
          name = if (!is.null(nm)) nm else basename(obj)
        )
      } else if (.is_ncdf4(obj)) {
        sources[[i]] <- list(
          source = obj,
          type = "connection",
          name = if (!is.null(nm)) nm else NA_character_
        )
      } else if (is.data.frame(obj)) {
        .check_required_cols(obj, file_col, "nc")
        files <- as.character(obj[[file_col]])
        sources <- lapply(files, function(x) {
          list(source = x, type = "path", name = if (!is.null(nm)) nm else basename(x))
        })
      } else {
        stop(
          "Unsupported object inside `nc`. Use file paths, opened ncdf4 objects, ",
          "lists of these objects, or a data frame with a file column.",
          call. = FALSE
        )
      }
    }

    return(sources)
  }

  stop(
    "Unsupported `nc` input. Use a file path, vector of file paths, opened ncdf4 object, ",
    "list, or data frame with a file column.",
    call. = FALSE
  )
}

.open_source <- function(source) {
  if (identical(source$type, "path")) {
    if (!file.exists(source$source)) {
      stop("netCDF file does not exist: ", source$source, call. = FALSE)
    }

    return(list(nc = ncdf4::nc_open(source$source), close = TRUE))
  }

  if (identical(source$type, "connection")) {
    return(list(nc = source$source, close = FALSE))
  }

  stop("Unknown internal netCDF source type.", call. = FALSE)
}

.close_source <- function(opened) {
  if (isTRUE(opened$close)) {
    try(ncdf4::nc_close(opened$nc), silent = TRUE)
  }

  invisible(TRUE)
}


.detect_netcdf_var <- function(nc) {
  # When `var = NULL`, the package tries to detect the variable automatically.
  # This is useful for workflows where each netCDF file contains one main
  # environmental variable, for example one file for NPPV and another file for
  # salinity.
  #
  # The detection is deliberately conservative. If more than one variable is
  # present, the function stops and asks the user to provide `var` manually. This
  # avoids extracting the wrong variable by accident.

  vars <- names(nc$var)

  if (length(vars) == 1L) {
    return(vars[1])
  }

  stop(
    "Automatic variable detection requires exactly one variable in each netCDF file. ",
    "Variables found: ",
    paste(vars, collapse = ", "),
    ". Please specify `var` manually.",
    call. = FALSE
  )
}

.resolve_source_var <- function(source, var) {
  # If the user provides `var`, the same variable name is used for all sources.
  # If `var = NULL`, each source is opened briefly and its variable is detected.

  if (!is.null(var)) {
    return(var)
  }

  opened <- .open_source(source)
  on.exit(.close_source(opened), add = TRUE)

  .detect_netcdf_var(opened$nc)
}

.resolve_output_prefix <- function(method, var, output_prefix, auto_var) {
  # If the user provides `output_prefix` or `output_col`, that name is respected.
  # If `var = NULL`, the output name is generated from the detected variable and
  # the extraction method, so multi-variable calls create clear column names:
  #   bottom  -> seabottom_nppv, seabottom_so
  #   surface -> surface_nppv, surface_so
  #   nearest -> nearest_nppv, nearest_so
  #   all     -> surface_*, nearest_* and seabottom_*
  #
  # If the user provided `var`, the previous behaviour is kept.

  if (!is.null(output_prefix)) {
    return(output_prefix)
  }

  if (!isTRUE(auto_var)) {
    return(NULL)
  }

  if (identical(method, "surface")) {
    return(paste0("surface_", var))
  }

  if (identical(method, "bottom")) {
    return(paste0("seabottom_", var))
  }

  if (identical(method, "nearest")) {
    return(paste0("nearest_", var))
  }

  var
}

.match_name <- function(available, candidates, required = TRUE, what = "dimension") {
  hit <- candidates[candidates %in% available]

  if (length(hit) > 0) {
    return(hit[1])
  }

  if (isTRUE(required)) {
    stop(
      sprintf(
        "Could not find %s. Tried: %s. Available names are: %s",
        what,
        paste(candidates, collapse = ", "),
        paste(available, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  NULL
}

.get_dim_values <- function(nc, dim_name) {
  if (is.null(dim_name)) {
    return(NULL)
  }

  if (dim_name %in% names(nc$dim)) {
    return(nc$dim[[dim_name]]$vals)
  }

  # Some files store coordinate values as variables rather than only as dimension
  # values. This fallback covers those files.
  if (dim_name %in% names(nc$var)) {
    return(ncdf4::ncvar_get(nc, dim_name))
  }

  stop("Could not read dimension or coordinate variable: ", dim_name, call. = FALSE)
}

.get_time_origin <- function(nc, time_dim) {
  # netCDF time is commonly stored as numeric values with units such as
  # "days since 1950-01-01". This helper reads that metadata and converts it into
  # a format that R can use.
  att <- try(ncdf4::ncatt_get(nc, time_dim, "units"), silent = TRUE)

  if (inherits(att, "try-error") || is.null(att$value) || is.na(att$value)) {
    return(NULL)
  }

  units <- att$value
  origin <- sub(".*since[[:space:]]+", "", units, ignore.case = TRUE)
  origin <- sub("[[:space:]]+.*$", "", origin)

  multiplier <- if (grepl("seconds since", units, ignore.case = TRUE)) {
    1 / 86400
  } else if (grepl("minutes since", units, ignore.case = TRUE)) {
    1 / 1440
  } else if (grepl("hours since", units, ignore.case = TRUE)) {
    1 / 24
  } else {
    1
  }

  list(origin = origin, multiplier = multiplier, units = units)
}

.convert_time_values <- function(nc, time_dim) {
  if (is.null(time_dim)) {
    return(NULL)
  }

  vals <- .get_dim_values(nc, time_dim)

  if (inherits(vals, "Date")) {
    return(vals)
  }

  if (inherits(vals, "POSIXt")) {
    return(as.Date(vals))
  }

  if (is.character(vals)) {
    return(as.Date(vals))
  }

  info <- .get_time_origin(nc, time_dim)

  if (is.null(info)) {
    # If there is no usable time metadata, the function falls back to treating the
    # values as day offsets from 1970-01-01. This is not ideal, but it keeps the
    # behaviour explicit and avoids silently assuming a product-specific origin.
    return(as.Date(vals, origin = "1970-01-01"))
  }

  as.Date(vals * info$multiplier, origin = info$origin)
}

.nearest_index <- function(values, target) {
  if (length(values) == 0 || all(is.na(values)) || is.na(target)) {
    return(NA_integer_)
  }

  which.min(abs(values - target))
}

.match_time_index <- function(values, target, time_match = c("nearest", "exact"),
                              max_time_diff = NULL) {
  time_match <- match.arg(time_match)

  if (length(values) == 0 || all(is.na(values)) || is.na(target)) {
    return(NA_integer_)
  }

  if (identical(time_match, "exact")) {
    hit <- which(values == target)
    return(if (length(hit) > 0) hit[1] else NA_integer_)
  }

  idx <- which.min(abs(values - target))

  if (!is.null(max_time_diff)) {
    diff_days <- abs(as.numeric(values[idx] - target))
    if (is.na(diff_days) || diff_days > max_time_diff) return(NA_integer_)
  }

  idx
}

.prepare_lon <- function(lon_values, lon_target) {
  # Many global netCDF files store longitude from 0 to 360, whereas observation
  # tables often use -180 to 180. This conversion is applied only when the file
  # clearly uses 0 to 360 and the observation longitude is negative.
  if (is.na(lon_target)) {
    return(lon_target)
  }

  if (min(lon_values, na.rm = TRUE) >= 0 && max(lon_values, na.rm = TRUE) > 180 && lon_target < 0) {
    lon_target <- lon_target + 360
  }

  lon_target
}

.detect_dims <- function(nc,
                         var,
                         lon_dim = NULL,
                         lat_dim = NULL,
                         depth_dim = NULL,
                         time_dim = NULL) {
  if (!var %in% names(nc$var)) {
    stop(
      "Variable `", var, "` not found in netCDF file. Available variables are: ",
      paste(names(nc$var), collapse = ", "),
      call. = FALSE
    )
  }

  available_dims <- names(nc$dim)
  var_dims <- vapply(nc$var[[var]]$dim, function(z) z$name, character(1))

  lon_candidates <- unique(stats::na.omit(c(lon_dim, "lon", "longitude", "x", "nav_lon")))
  lat_candidates <- unique(stats::na.omit(c(lat_dim, "lat", "latitude", "y", "nav_lat")))
  depth_candidates <- unique(stats::na.omit(c(depth_dim, "depth", "deptht", "lev", "level", "z")))
  time_candidates <- unique(stats::na.omit(c(time_dim, "time", "time_counter", "t")))

  lon_name <- .match_name(c(var_dims, available_dims), lon_candidates, TRUE, "longitude dimension")
  lat_name <- .match_name(c(var_dims, available_dims), lat_candidates, TRUE, "latitude dimension")
  depth_name <- .match_name(c(var_dims, available_dims), depth_candidates, FALSE, "depth dimension")
  time_name <- .match_name(c(var_dims, available_dims), time_candidates, FALSE, "time dimension")

  list(
    lon = lon_name,
    lat = lat_name,
    depth = depth_name,
    time = time_name,
    var_dims = var_dims
  )
}

.build_start_count <- function(var_dims, index_list) {
  # ncdf4::ncvar_get() expects indices in the exact order used by the variable in
  # the netCDF file. This helper builds `start` and `count` vectors by matching
  # dimension names rather than assuming a fixed order such as lon-lat-depth-time.
  start <- integer(length(var_dims))
  count <- integer(length(var_dims))

  for (i in seq_along(var_dims)) {
    dim_name <- var_dims[i]
    idx <- index_list[[dim_name]]

    if (is.null(idx)) {
      start[i] <- 1L
      count[i] <- -1L
    } else {
      start[i] <- as.integer(idx$start)
      count[i] <- as.integer(idx$count)
    }
  }

  list(start = start, count = count)
}

.extract_array <- function(nc, var, dims, lon_idx, lat_idx, time_idx = NULL, depth_start = NULL, depth_count = NULL) {
  index_list <- list()
  index_list[[dims$lon]] <- list(start = lon_idx, count = 1L)
  index_list[[dims$lat]] <- list(start = lat_idx, count = 1L)

  if (!is.null(dims$time) && !is.null(time_idx)) {
    index_list[[dims$time]] <- list(start = time_idx, count = 1L)
  }

  if (!is.null(dims$depth) && !is.null(depth_start)) {
    index_list[[dims$depth]] <- list(start = depth_start, count = depth_count)
  }

  sc <- .build_start_count(dims$var_dims, index_list)

  ncdf4::ncvar_get(
    nc = nc,
    varid = var,
    start = sc$start,
    count = sc$count
  )
}

.extract_one_file <- function(data,
                              source,
                              var,
                              method,
                              lon_col,
                              lat_col,
                              date_col,
                              depth_col,
                              id_col,
                              lon_dim,
                              lat_dim,
                              depth_dim,
                              time_dim,
                              time_match,
                              max_time_diff,
                              output_prefix,
                              verbose) {
  opened <- .open_source(source)
  on.exit(.close_source(opened), add = TRUE)
  nc <- opened$nc

  dims <- .detect_dims(
    nc = nc,
    var = var,
    lon_dim = lon_dim,
    lat_dim = lat_dim,
    depth_dim = depth_dim,
    time_dim = time_dim
  )

  if (method != "2d" && is.null(dims$depth)) {
    stop("Method `", method, "` requires a depth dimension, but none was found.", call. = FALSE)
  }

  lon_values <- .get_dim_values(nc, dims$lon)
  lat_values <- .get_dim_values(nc, dims$lat)
  depth_values <- .get_dim_values(nc, dims$depth)
  time_values <- .convert_time_values(nc, dims$time)

  if (!is.null(dims$time) && !is.null(date_col)) {
    date_values <- as.Date(data[[date_col]])
  } else {
    date_values <- rep(NA, nrow(data))
  }

  prefix <- if (is.null(output_prefix)) var else output_prefix

  if (method == "all") {
    out <- data # `all` now returns three summary columns, not one column per depth layer.
    out[[paste0("nearest_", prefix)]] <- NA_real_ # nearest valid depth summary column.
    out[[paste0("surface_", prefix)]] <- NA_real_ # surface layer summary column.
    out[[paste0("seabottom_", prefix)]] <- NA_real_ # deepest valid layer summary column.
  } else {
    out <- data
    out[[prefix]] <- NA_real_
  }

  for (i in seq_len(nrow(data))) {
    if (isTRUE(verbose) && i %% 100 == 0) {
      message("Extracting row ", i, " of ", nrow(data))
    }

    lon_target <- .prepare_lon(lon_values, data[[lon_col]][i])
    lat_target <- data[[lat_col]][i]

    lon_idx <- .nearest_index(lon_values, lon_target)
    lat_idx <- .nearest_index(lat_values, lat_target)

    if (is.na(lon_idx) || is.na(lat_idx)) next

    time_idx <- NULL
    if (!is.null(time_values)) {
      time_idx <- .match_time_index(
        values = time_values,
        target = date_values[i],
        time_match = time_match,
        max_time_diff = max_time_diff
      )
      if (is.na(time_idx)) next
    }

    if (method == "2d") {
      val <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx
      )
      out[[prefix]][i] <- as.numeric(val)[1]
    }

    if (method == "surface") {
      val <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = 1L
      )
      out[[prefix]][i] <- as.numeric(val)[1]
    }

    if (method == "nearest") {
      depth_target <- data[[depth_col]][i]

      # read the full profile before selecting the nearest depth.
      # this allows the function to ignore depth layers that are NA,
      # which usually correspond to layers below the seabed or outside
      # the valid water column for that grid cell.
      profile <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = length(depth_values)
      )
      profile <- as.numeric(profile)
      valid <- which(!is.na(profile))

      if (length(valid) > 0 && !is.na(depth_target)) {
        depth_idx <- valid[which.min(abs(depth_values[valid] - depth_target))]
        out[[prefix]][i] <- profile[depth_idx]
      }
    }

    if (method == "bottom") {
      profile <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = length(depth_values)
      )
      profile <- as.numeric(profile)
      valid <- which(!is.na(profile))
      if (length(valid) > 0) {
        out[[prefix]][i] <- profile[utils::tail(valid, 1)]
      }
    }

    if (method == "all") {
      # `all` now extracts one vertical profile and summarises it into
      # nearest valid depth, surface and bottom columns. This mirrors
      # applying extract3d_nearest(), extract3d_surface() and
      # extract3d_bottom() separately, but avoids reading the same profile
      # three times.
      profile <- .extract_array(
        nc = nc,
        var = var,
        dims = dims,
        lon_idx = lon_idx,
        lat_idx = lat_idx,
        time_idx = time_idx,
        depth_start = 1L,
        depth_count = length(depth_values)
      )
      profile <- as.numeric(profile)
      valid <- which(!is.na(profile))

      out[[paste0("surface_", prefix)]][i] <- profile[1]

      if (length(valid) > 0) {
        out[[paste0("seabottom_", prefix)]][i] <- profile[utils::tail(valid, 1)]

        depth_target <- data[[depth_col]][i]
        if (!is.na(depth_target)) {
          depth_idx <- valid[which.min(abs(depth_values[valid] - depth_target))]
          out[[paste0("nearest_", prefix)]][i] <- profile[depth_idx]
        }
      }
    }
  }

  out
}

.merge_source_outputs <- function(outputs, original_names) {
  # When several files are provided, the function extracts from each file and then
  # combines the results. For ordinary single-variable extraction, the first
  # non-missing value across files is used. This is useful when files represent
  # different time periods or spatial tiles.
  if (length(outputs) == 1L) {
    return(outputs[[1]])
  }

  out <- outputs[[1]][original_names]
  new_cols <- unique(unlist(lapply(outputs, function(x) setdiff(names(x), original_names))))

  for (col in new_cols) {
    values <- rep(NA_real_, nrow(out))

    for (x in outputs) {
      if (!col %in% names(x)) next
      replace <- is.na(values) & !is.na(x[[col]])
      values[replace] <- x[[col]][replace]
    }

    out[[col]] <- values
  }

  out
}

# -----------------------------------------------------------------------------
# Main extraction function
# -----------------------------------------------------------------------------

#' Extract values from netCDFs to one or more point observations
#'
#' This is the general point-extraction interface. It can extract environmental
#' values for one or multiple observations. The user chooses the type of
#' extraction with the `method` argument. The more specific functions
#' `extract2d()`, `extract3d_surface()`, `extract3d_bottom()`,
#' `extract3d_nearest()` and `extract3d_all()` are wrappers around this function.
#'
#' @param data Optional object containing observation points. Can be a data frame,
#'   tibble, `sf` object with POINT geometries, matrix with column names, or named
#'   list. If `NULL`, coordinates can be supplied directly through `lon`, `lat`,
#'   and optionally `depth` and `date`. For `sf` objects, longitude and latitude
#'   are obtained from the POINT geometries; projected geometries are transformed
#'   internally to longitude and latitude when a CRS is available.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of file
#'   paths, opened `ncdf4` object, list of opened objects, or a data frame with a
#'   file column.
#' @param lon Optional numeric longitude value or vector for direct point
#'   extraction when `data = NULL`.
#' @param lat Optional numeric latitude value or vector for direct point
#'   extraction when `data = NULL`.
#' @param depth Optional numeric observation depth value or vector for direct
#'   point extraction. Required when `data = NULL` and `method = "nearest"` or
#'   `method = "all"`.
#' @param date Optional observation date or vector of dates for direct point
#'   extraction. Used when `data = NULL` and the netCDF variable has a time
#'   dimension. Character dates should use `"YYYY-MM-DD"` format. Direct input
#'   vectors may have length 1 or a common compatible length.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param method Extraction method. Options are `"2d"`, `"surface"`, `"bottom"`,
#'   `"nearest"` and `"all"`. With `method = "all"`, the function returns
#'   nearest, surface and bottom outputs together.
#' @param lon_col Name of the longitude column in `data`.
#' @param lat_col Name of the latitude column in `data`.
#' @param date_col Name of the date column in `data`. Required only when the
#'   netCDF variable has a time dimension.
#' @param depth_col Name of the observation-depth column in `data`. Required
#'   only when `method = "nearest"` or `method = "all"`.
#' @param id_col Optional observation identifier column. It is not required for
#'   extraction, but it is checked when provided to help users detect mistakes.
#' @param file_col Column containing file paths when `nc` is a data frame.
#' @param lon_dim Optional name of the longitude dimension in the netCDF file.
#' @param lat_dim Optional name of the latitude dimension in the netCDF file.
#' @param depth_dim Optional name of the depth dimension in the netCDF file.
#' @param time_dim Optional name of the time dimension in the netCDF file.
#' @param time_match Character. How observation dates are matched to netCDF time
#'   values. `"nearest"` (default) selects the closest available date; `"exact"`
#'   requires an exact date match and returns `NA` when no match exists.
#' @param max_time_diff Optional non-negative numeric value giving the maximum
#'   allowed difference, in days, between an observation date and the nearest
#'   netCDF date when `time_match = "nearest"`. If `NULL` (default), no maximum
#'   temporal distance is imposed. Ignored when `time_match = "exact"`.
#' @param output_prefix Optional name used for the extracted output column. If
#'   omitted, `var` is used.
#' @param output_col Optional alias for `output_prefix`. Use this when you want
#'   to define the exact name of the output column directly. If both
#'   `output_col` and `output_prefix` are provided, `output_col` is used. With
#'   `method = "all"`, this name is used as the suffix in `nearest_*`,
#'   `surface_*` and `seabottom_*` columns.
#' @param verbose Logical. If `TRUE`, progress messages are printed every 100
#'   rows.
#'
#' @return If `data` is provided, the input observations with extracted values
#'   appended. Data frames, tibbles and `sf` objects retain their structure;
#'   matrices and named lists are returned as data frames. If `data = NULL`,
#'   direct extraction returns a numeric value for one point and one output, a
#'   numeric vector for multiple points and one output, a named numeric vector
#'   for one point and multiple outputs (for example, `method = "all"`), or a
#'   data frame for multiple points and multiple outputs.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' # Extract the value closest to a single observation
#' extract_to_point(
#'   nc = nc_file,
#'   lon = 0,
#'   lat = 40,
#'   depth = 60,
#'   date = "2020-01-03",
#'   var = "temp",
#'   method = "nearest"
#' )
#'
#' # Extract values for multiple observations
#' extract_to_point(
#'   nc = nc_file,
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   depth = c(60, 100),
#'   date = "2020-01-03",
#'   var = "temp",
#'   method = "nearest"
#' )
#'
#' @export

extract_to_point <- function(data = NULL,
                             nc,
                             lon = NULL,
                             lat = NULL,
                             depth = NULL,
                             date = NULL,
                             var = NULL,
                             method = c("2d", "surface", "bottom", "nearest", "all"),
                             lon_col = "lon",
                             lat_col = "lat",
                             date_col = "date",
                             depth_col = "depth",
                             id_col = NULL,
                             file_col = "file",
                             lon_dim = NULL,
                             lat_dim = NULL,
                             depth_dim = NULL,
                             time_dim = NULL,
                             time_match = c("nearest", "exact"),
                             max_time_diff = NULL,
                             output_prefix = NULL,
                             output_col = NULL,
                             verbose = FALSE) {
  method <- match.arg(method)
  time_match <- match.arg(time_match)

  if (!is.null(max_time_diff) &&
      (!is.numeric(max_time_diff) || length(max_time_diff) != 1L ||
       is.na(max_time_diff) || max_time_diff < 0)) {
    stop("`max_time_diff` must be `NULL` or one non-negative numeric value in days.",
         call. = FALSE)
  }

  direct_input <- is.null(data)

  if (direct_input) {
    if (!is.numeric(lon) || length(lon) == 0L || anyNA(lon)) {
      stop("`lon` must contain numeric values.", call. = FALSE)
    }

    if (!is.numeric(lat) || length(lat) == 0L || anyNA(lat)) {
      stop("`lat` must contain numeric values.", call. = FALSE)
    }

    if (method %in% c("nearest", "all") &&
        (!is.numeric(depth) || length(depth) == 0L || anyNA(depth))) {
      stop(
        "`depth` must contain numeric values for `nearest` and `all` extraction.",
        call. = FALSE
      )
    }

    direct_values <- .recycle_point_inputs(
      lon = lon,
      lat = lat,
      depth = depth,
      date = date
    )

    data <- data.frame(
      lon = direct_values$lon,
      lat = direct_values$lat
    )

    lon_col <- "lon"
    lat_col <- "lat"

    if (!is.null(direct_values$depth)) {
      data$depth <- direct_values$depth
      depth_col <- "depth"
    }

    if (!is.null(direct_values$date)) {
      date_value <- suppressWarnings(
        tryCatch(
          as.Date(direct_values$date),
          error = function(e) rep(as.Date(NA), length(direct_values$date))
        )
      )

      if (anyNA(date_value)) {
        stop(
          "`date` could not be converted to a valid date. Use `YYYY-MM-DD` format, for example `2020-01-03`.",
          call. = FALSE
        )
      }

      data$date <- date_value
      date_col <- "date"
    } else {
      date_col <- NULL
    }
  } else {
    data <- .normalize_points_input(
      data,
      lon_col = lon_col,
      lat_col = lat_col
    )
  }

  if (!is.null(var) && (!is.character(var) || length(var) != 1L || is.na(var))) {
    stop("`var` must be `NULL` or a single character value naming one netCDF variable.", call. = FALSE)
  }

  # allow users to define the output column name with `output_col`.
  # `output_prefix` is kept for backward compatibility with previous code.
  if (!is.null(output_col)) {
    if (!is.character(output_col) || length(output_col) != 1L || is.na(output_col)) {
      stop("`output_col` must be a single character value naming the output column.", call. = FALSE)
    }
    output_prefix <- output_col
  }

  required_cols <- c(lon_col, lat_col)
  if (!is.null(date_col)) required_cols <- c(required_cols, date_col)
  if (method %in% c("nearest", "all")) required_cols <- c(required_cols, depth_col) # `all` needs depth to calculate nearest-depth output.
  if (!is.null(id_col)) required_cols <- c(required_cols, id_col)
  .check_required_cols(data, required_cols, "data")

  sources <- .as_netcdf_sources(nc, file_col = file_col)

  auto_var <- is.null(var)

  outputs <- lapply(sources, function(src) {
    source_var <- .resolve_source_var(src, var)

    source_output_prefix <- .resolve_output_prefix(
      method = method,
      var = source_var,
      output_prefix = output_prefix,
      auto_var = auto_var
    )

    .extract_one_file(
      data = data,
      source = src,
      var = source_var,
      method = method,
      lon_col = lon_col,
      lat_col = lat_col,
      date_col = date_col,
      depth_col = depth_col,
      id_col = id_col,
      lon_dim = lon_dim,
      lat_dim = lat_dim,
      depth_dim = depth_dim,
      time_dim = time_dim,
      time_match = time_match,
      max_time_diff = max_time_diff,
      output_prefix = source_output_prefix,
      verbose = verbose
    )
  })

  out <- .merge_source_outputs(outputs, original_names = names(data))

  if (direct_input) {
    value_cols <- setdiff(names(out), names(data))

    if (length(value_cols) == 1L) {
      return(out[[value_cols]])
    }

    if (nrow(out) == 1L) {
      return(unlist(out[1, value_cols], use.names = TRUE))
    }

    return(out[value_cols])
  }

  out
}

# -----------------------------------------------------------------------------
# Wrapper functions
# -----------------------------------------------------------------------------
# These functions keep explicit names for the most common extraction modes. They
# are easier to discover and easier to document, while the internal workflow stays
# centralised in `extract_to_point()`.
# -----------------------------------------------------------------------------

#' Extract values from a 2D netCDF variable
#'
#' Use this function when the variable has longitude and latitude dimensions,
#' and optionally a time dimension, but no depth dimension.
#'
#' @param data Object containing observation points. Can be a data frame, tibble,
#'   `sf` object with POINT geometries, matrix with column names, or named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of file
#'   paths, opened `ncdf4` object, list of opened objects, or a data frame with a
#'   file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with extracted values appended. Data frames,
#'   tibbles and `sf` objects retain their structure; matrices and named lists
#'   are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_2d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract2d(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export
extract2d <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "2d", ...)
}

#' Extract the surface layer from a 3D netCDF variable
#'
#' The surface layer is defined as the first available depth layer in the
#' netCDF file. This is usually the shallowest layer.
#'
#' @param data Object containing observation points. Can be a data frame, tibble,
#'   `sf` object with POINT geometries, matrix with column names, or named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of file
#'   paths, opened `ncdf4` object, list of opened objects, or a data frame with a
#'   file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with extracted values appended. Data frames,
#'   tibbles and `sf` objects retain their structure; matrices and named lists
#'   are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract3d_surface(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export
extract3d_surface <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "surface", ...)
}

#' Extract the bottom available layer from a 3D netCDF variable
#'
#' For each observation, this function extracts the full vertical profile at
#' the nearest longitude, latitude and time cell, then returns the deepest
#' non-missing value. This is useful when the bathymetry of the netCDF grid
#' means that deeper layers are missing over shallow areas.
#'
#' @param data Object containing observation points. Can be a data frame, tibble,
#'   `sf` object with POINT geometries, matrix with column names, or named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of file
#'   paths, opened `ncdf4` object, list of opened objects, or a data frame with a
#'   file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with extracted values appended. Data frames,
#'   tibbles and `sf` objects retain their structure; matrices and named lists
#'   are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract3d_bottom(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export
extract3d_bottom <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "bottom", ...)
}

#' Extract the nearest valid depth layer from a 3D netCDF variable
#'
#' For each observation, this function reads the full vertical profile and
#' finds the valid, non-missing netCDF depth layer closest to the observation
#' depth. This avoids returning `NA` when the geometrically closest depth layer
#' is invalid at that grid cell, for example below the seabed.
#'
#' @param data Object containing observation points. Can be a data frame, tibble,
#'   `sf` object with POINT geometries, matrix with column names, or named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of file
#'   paths, opened `ncdf4` object, list of opened objects, or a data frame with a
#'   file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with extracted values appended. Data frames,
#'   tibbles and `sf` objects retain their structure; matrices and named lists
#'   are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   depth = c(60, 100),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract3d_nearest(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export

extract3d_nearest <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "nearest", ...)
}

#' Extract nearest, surface and bottom values from a 3D netCDF variable
#'
#' For each observation, this function extracts the full vertical profile at
#' the nearest longitude, latitude and time cell. It then returns three summary
#' columns: nearest valid depth, surface and bottom.
#'
#' @param data Object containing observation points. Can be a data frame, tibble,
#'   `sf` object with POINT geometries, matrix with column names, or named list.
#' @param nc netCDF source. Can be a file path, vector of file paths, list of file
#'   paths, opened `ncdf4` object, list of opened objects, or a data frame with a
#'   file column.
#' @param var Name of the variable to extract from the netCDF file. If `NULL`,
#'   the function tries to detect the variable automatically from each file.
#'   This only works when each netCDF file contains one single variable.
#' @param ... Additional arguments passed to [extract_to_point()].
#'
#' @return The input observations with three extracted columns: `nearest_*`,
#'   `surface_*` and `seabottom_*`. Data frames, tibbles and `sf` objects retain
#'   their structure; matrices and named lists are returned as data frames.
#'
#' @examples
#' nc_file <- system.file(
#'   "extdata", "example_3d.nc",
#'   package = "ocean3d"
#' )
#'
#' observations <- data.frame(
#'   lon = c(0, 1),
#'   lat = c(40, 41),
#'   depth = c(60, 100),
#'   date = as.Date(c("2020-01-01", "2020-01-03"))
#' )
#'
#' extract3d_all(
#'   data = observations,
#'   nc = nc_file,
#'   var = "temp"
#' )
#'
#' @export

extract3d_all <- function(data, nc, var = NULL, ...) {
  extract_to_point(data = data, nc = nc, var = var, method = "all", ...)
}
