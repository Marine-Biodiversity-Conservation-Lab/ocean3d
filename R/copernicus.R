#' Copernicus cache directory
#'
#' Returns the path to the package's persistent cache directory for downloaded
#' Copernicus environmental data. Uses [tools::R_user_dir()] so the location
#' survives across sessions and follows platform conventions.
#'
#' @returns Character. Path to cache directory (created if missing).
#' @export
copernicus_cache_dir <- function() {
  path <- file.path(tools::R_user_dir("sharkabc3d", which = "cache"), "copernicus")
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

#' Clear the Copernicus cache
#'
#' Remove all files downloaded to the default Copernicus cache.
#'
#' @param confirm Logical. Require interactive confirmation. Default `TRUE`.
#'
#' @returns Invisibly, `TRUE` on success.
#' @export
copernicus_cache_clear <- function(confirm = TRUE) {
  path <- copernicus_cache_dir()

  if (confirm && interactive()) {
    ans <- readline(sprintf("Delete all cached Copernicus files in %s? [y/N]: ", path))
    if (!tolower(ans) %in% c("y", "yes")) {
      message("Cancelled.")
      return(invisible(FALSE))
    }
  }

  unlink(path, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

# Internal: one-time interactive consent for writing to persistent cache.
.copernicus_cache_consent <- function(cache_dir) {
  sentinel <- file.path(cache_dir, ".consent")
  if (file.exists(sentinel)) return(invisible(TRUE))

  msg <- paste0(
    "sharkabc3d will cache downloaded Copernicus files in:\n  ",
    cache_dir, "\n",
    "Copernicus environmental datasets can be large.\n",
    "Pass `output_dir` to use a different location, or call ",
    "copernicus_cache_clear() to reclaim space later."
  )

  if (!interactive()) {
    stop(
      msg,
      "\nNon-interactive session: pass `output_dir` explicitly, or run ",
      "copernicus_load() once interactively to record consent.",
      call. = FALSE
    )
  }

  message(msg)
  ans <- readline("Proceed with caching here? [y/N]: ")
  if (!tolower(ans) %in% c("y", "yes")) stop("Cache consent declined.", call. = FALSE)

  file.create(sentinel)
  invisible(TRUE)
}

# Internal: validate a single character argument.
.copernicus_validate_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be a single non-empty character string.", call. = FALSE)
  }
  invisible(TRUE)
}

# Internal: validate geographical bounds.
.copernicus_validate_bbox <- function(xmin, xmax, ymin, ymax) {
  bbox <- list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  supplied <- !vapply(bbox, is.null, logical(1))

  if (!any(supplied)) return(invisible(TRUE))
  if (!all(supplied)) {
    stop("`xmin`, `xmax`, `ymin`, and `ymax` must be supplied together.", call. = FALSE)
  }

  vals <- unlist(bbox, use.names = FALSE)
  if (!is.numeric(vals) || any(!is.finite(vals)) || any(lengths(bbox) != 1)) {
    stop("Geographical bounds must be single finite numeric values.", call. = FALSE)
  }

  if (xmin >= xmax) {
    stop("Longitude bounds must satisfy `xmin < xmax`.", call. = FALSE)
  }

  if (ymin < -90 || ymax > 90 || ymin >= ymax) {
    stop("Latitude bounds must satisfy -90 <= ymin < ymax <= 90.", call. = FALSE)
  }

  invisible(TRUE)
}

# Internal: validate depth range.
.copernicus_validate_depth <- function(depth_min, depth_max) {
  if (is.null(depth_min) && is.null(depth_max)) return(invisible(TRUE))

  depth <- list(depth_min = depth_min, depth_max = depth_max)
  supplied <- !vapply(depth, is.null, logical(1))
  vals <- unlist(depth[supplied], use.names = FALSE)

  if (!is.numeric(vals) || any(!is.finite(vals)) ||
      any(lengths(depth[supplied]) != 1)) {
    stop("Depth bounds must be single finite numeric values.", call. = FALSE)
  }

  if (any(vals < 0)) {
    stop("Depth bounds must be >= 0 metres (positive downward).", call. = FALSE)
  }

  if (all(supplied) && depth_min > depth_max) {
    stop("`depth_min` must be <= `depth_max`.", call. = FALSE)
  }

  invisible(TRUE)
}

# Internal: convert Date/POSIXt/character input to UTC ISO datetime.
.copernicus_datetime <- function(x, name = "datetime") {
  if (is.null(x)) return(NULL)
  
  if (inherits(x, "Date")) {
    x <- as.POSIXct(x, tz = "UTC")
  } else if (inherits(x, "POSIXt")) {
    x <- as.POSIXct(x, tz = "UTC")
  } else if (is.character(x) && length(x) == 1 && !is.na(x)) {
    if (grepl("^\\d{4}-\\d{2}-\\d{2}$", x)) x <- paste0(x, " 00:00:00")
    
    x <- tryCatch(
      suppressWarnings(as.POSIXct(
        gsub("T", " ", sub("Z$", "", x), fixed = TRUE),
        tz = "UTC",
        tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d")
      )),
      error = function(e) as.POSIXct(NA, tz = "UTC")
    )
  } else {
    stop(
      "`", name, "` must be a Date, POSIXt, or single character datetime.",
      call. = FALSE
    )
  }
  
  if (length(x) != 1 || is.na(x)) {
    stop("Could not parse `", name, "`.", call. = FALSE)
  }
  
  format(x, "%Y-%m-%dT%H:%M:%S", tz = "UTC")
}

# Internal: prepare and validate temporal bounds.
.copernicus_time_range <- function(start_datetime, end_datetime) {
  start <- .copernicus_datetime(start_datetime, "start_datetime")
  end <- .copernicus_datetime(end_datetime, "end_datetime")

  if (is.null(start) && !is.null(end)) {
    stop("`start_datetime` is required when `end_datetime` is supplied.", call. = FALSE)
  }

  if (!is.null(start) && is.null(end)) end <- start

  if (!is.null(start) &&
      as.POSIXct(end, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC") <
      as.POSIXct(start, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")) {
    stop("`end_datetime` must be equal to or later than `start_datetime`.",
         call. = FALSE)
  }

  list(start = start, end = end)
}

# Internal: extract paths from Copernicus Marine response objects.
.copernicus_marine_paths <- function(result) {
  get_one <- function(x) {
    path <- tryCatch(as.character(x$file_path), error = function(e) character())
    if (!length(path)) return(character())

    file_names <- tryCatch(
      unlist(x$file_names, use.names = FALSE),
      error = function(e) character()
    )

    if (length(path) == 1 && dir.exists(path) && length(file_names)) {
      candidates <- file.path(path, file_names)
      if (all(file.exists(candidates))) return(candidates)
    }

    path
  }

  if (is.list(result) && is.null(result$file_path)) {
    paths <- unlist(lapply(result, get_one), use.names = FALSE)
  } else {
    paths <- get_one(result)
  }

  paths <- unique(as.character(paths))
  if (!length(paths)) {
    stop("Copernicus Marine completed the request but returned no output path.",
         call. = FALSE)
  }

  paths
}

#' Load environmental data from Copernicus services
#'
#' Acquire environmental data from Copernicus Marine or ECMWF-operated
#' Copernicus Data Stores and save the resulting files locally. Supported
#' sources are Copernicus Marine (`"marine"`), Climate Data Store (`"cds"`),
#' and Atmosphere Data Store (`"ads"`).
#'
#' Copernicus Marine requests use the Python `copernicusmarine` package through
#' `reticulate`. CDS and ADS requests use `ecmwfr`.
#'
#' The common arguments provide a consistent interface for dataset, variable,
#' spatial, temporal, and depth selection. Dataset-specific request parameters
#' can be supplied through `request`. Explicit arguments to `copernicus_load()`
#' take precedence over duplicate entries in `request`.
#'
#' Authentication is managed through the official provider clients. Copernicus
#' Marine credentials are checked automatically and interactive login is offered
#' when needed. CDS and ADS use an ECMWF Personal Access Token (PAT),
#' configured with [ecmwfr::wf_set_key()]. If no PAT is found in an interactive
#' session, the official CDS API setup page is opened and token setup is started.
#'
#' For CERRA requests, geographic subsetting is performed locally when a bounding
#' box is supplied. The provider response is downloaded to a temporary directory,
#' the raster is cropped with `terra`, only the cropped result is retained, and
#' temporary full-domain files are removed automatically.
#'
#' @param source Character. Copernicus service: `"marine"`, `"cds"`, or `"ads"`.
#' @param dataset_id Character. Provider dataset identifier.
#' @param variables Optional character vector of variables to download.
#' @param start_datetime Optional start datetime. A `Date`, `POSIXt`, or
#'   character value. Times are interpreted in UTC.
#' @param end_datetime Optional end datetime. Defaults to `start_datetime`.
#' @param xmin,xmax,ymin,ymax Optional geographic bounding box in decimal
#'   degrees. The ordering used internally for ECMWF `area` is north, west,
#'   south, east.
#' @param depth_min,depth_max Optional depth range in metres, positive downward.
#'   Supported directly for Copernicus Marine. For CDS/ADS, vertical
#'   dimensions such as pressure or model levels should be supplied through
#'   `request`.
#' @param request Named list of provider- or dataset-specific request arguments.
#'   Examples include `product_type`, `pressure_level`, `data_format`,
#'   `download_format`, `service`, or `coordinates_selection_method`.
#' @param output_dir Character. Destination directory. Defaults to a
#'   source-specific subdirectory of [copernicus_cache_dir()].
#' @param filename Optional output filename. For ECMWF services this is passed as
#'   the request `target`.
#' @param split_by Optional character. For Copernicus Marine, split downloaded
#'   subsets by `"variable"`, `"hour"`, `"day"`, `"month"`, or `"year"`.
#' @param concurrent_processes Optional integer. Number of concurrent processes
#'   used by Copernicus Marine when `split_by` is supplied.
#' @param compression Optional integer from 0 to 9 giving the Copernicus Marine
#'   NetCDF compression level.
#' @param dataset_version Optional Copernicus Marine dataset version.
#' @param dataset_part Optional Copernicus Marine dataset part.
#' @param force Logical. Overwrite an identifiable existing output file.
#'   Default `FALSE`.
#' @param quiet Logical. Suppress progress messages where supported.
#'
#' @returns Character vector containing paths to downloaded files.
#'
#' @examples
#' \dontrun{
#' # Copernicus Marine
#' files <- copernicus_load(
#'   source = "marine",
#'   dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
#'   variables = c("thetao", "so"),
#'   start_datetime = "2020-01-01",
#'   end_datetime = "2020-01-31",
#'   xmin = -6, xmax = 10, ymin = 35, ymax = 45,
#'   depth_min = 0, depth_max = 500
#' )
#'
#' # ERA5 through the Climate Data Store
#' files <- copernicus_load(
#'   source = "cds",
#'   dataset_id = "reanalysis-era5-single-levels",
#'   variables = "2m_temperature",
#'   start_datetime = "2021-01-01T12:00:00",
#'   end_datetime = "2021-01-31T12:00:00",
#'   xmin = -2, xmax = 4, ymin = 36, ymax = 43,
#'   request = list(
#'     product_type = "reanalysis",
#'     data_format = "netcdf",
#'     download_format = "unarchived"
#'   ),
#'   filename = "era5_temperature.nc"
#' )
#' }
#'
#' @export
copernicus_load <- function(source = c("marine", "cds", "ads"),
                            dataset_id,
                            variables = NULL,
                            start_datetime = NULL,
                            end_datetime = start_datetime,
                            xmin = NULL, xmax = NULL,
                            ymin = NULL, ymax = NULL,
                            depth_min = NULL, depth_max = NULL,
                            request = list(),
                            output_dir = NULL,
                            filename = NULL,
                            split_by = NULL,
                            concurrent_processes = NULL,
                            compression = NULL,
                            dataset_version = NULL,
                            dataset_part = NULL,
                            force = FALSE,
                            quiet = FALSE) {

  source <- match.arg(source)
  .copernicus_validate_string(dataset_id, "dataset_id")

  if (!is.null(variables) &&
      (!is.character(variables) || !length(variables) || anyNA(variables) ||
       any(!nzchar(variables)))) {
    stop("`variables` must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.list(request) ||
      (length(request) && (is.null(names(request)) || any(!nzchar(names(request)))))) {
    stop("`request` must be a named list.", call. = FALSE)
  }

  if (!is.logical(force) || length(force) != 1 || is.na(force)) {
    stop("`force` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(quiet) || length(quiet) != 1 || is.na(quiet)) {
    stop("`quiet` must be TRUE or FALSE.", call. = FALSE)
  }

  .copernicus_validate_bbox(xmin, xmax, ymin, ymax)
  .copernicus_validate_depth(depth_min, depth_max)
  time <- .copernicus_time_range(start_datetime, end_datetime)

  using_default_cache <- is.null(output_dir)

  if (using_default_cache) {
    cache_root <- copernicus_cache_dir()
    .copernicus_cache_consent(cache_root)
    output_dir <- file.path(cache_root, source)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!is.null(filename)) {
    .copernicus_validate_string(filename, "filename")
    dest <- file.path(output_dir, filename)

    if (file.exists(dest) && !force) {
      if (!quiet) message("Cached: ", filename)
      return(normalizePath(dest, winslash = "/", mustWork = TRUE))
    }

    if (file.exists(dest) && force) unlink(dest, recursive = TRUE, force = TRUE)
  }

  if (source == "marine") {
    .copernicus_load_marine(
      dataset_id = dataset_id, variables = variables,
      start_datetime = time$start, end_datetime = time$end,
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      depth_min = depth_min, depth_max = depth_max,
      request = request, output_dir = output_dir, filename = filename,
      split_by = split_by, concurrent_processes = concurrent_processes,
      compression = compression, dataset_version = dataset_version,
      dataset_part = dataset_part, force = force, quiet = quiet
    )
  } else {
    .copernicus_load_ecmwf(
      source = source, dataset_id = dataset_id, variables = variables,
      start_datetime = time$start, end_datetime = time$end,
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      depth_min = depth_min, depth_max = depth_max,
      request = request, output_dir = output_dir, filename = filename,
      split_by = split_by, concurrent_processes = concurrent_processes,
      compression = compression, dataset_version = dataset_version,
      dataset_part = dataset_part, force = force, quiet = quiet
    )
  }
}

# Internal: prepare the Python dependency used by Copernicus Marine.
.copernicus_marine_module <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Package 'reticulate' is required for `source = \"marine\"`. ",
      "Install it with install.packages('reticulate').",
      call. = FALSE
    )
  }

  # If Python has not yet been initialized, prefer reticulate's managed
  # environment unless the user has explicitly selected another Python through
  # RETICULATE_PYTHON, RETICULATE_PYTHON_ENV, use_python(), or use_virtualenv().
  # This avoids an unrelated pre-existing `r-reticulate` environment taking
  # precedence over the requirements declared below.
  if (!reticulate::py_available(initialize = FALSE)) {
    Sys.setenv(RETICULATE_USE_MANAGED_VENV = "yes")
  }

  reticulate::py_require(
    packages = "copernicusmarine>=2.2,<3",
    python_version = ">=3.10,<3.14"
  )

  tryCatch(
    reticulate::import("copernicusmarine", convert = TRUE),
    error = function(e) {
      py_initialized <- reticulate::py_available(initialize = FALSE)

      stop(
        "Could not import the Python package 'copernicusmarine'. ",
        if (py_initialized) {
          paste0(
            "Python is already initialized in this R session. Restart R and ",
            "call `copernicus_load()` before initializing another Python ",
            "environment, or explicitly select an environment containing ",
            "copernicusmarine with `reticulate::use_virtualenv()` immediately ",
            "after restarting R. "
          )
        } else {
          paste0(
            "reticulate could not provision a compatible managed Python ",
            "environment. "
          )
        },
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

# Internal: validate or interactively configure Copernicus Marine credentials.
.copernicus_marine_auth <- function(cm) {
  valid <- tryCatch(
    isTRUE(cm$login(check_credentials_valid = TRUE)),
    error = function(e) FALSE
  )

  if (valid) return(invisible(TRUE))

  if (!interactive()) {
    stop(
      "Valid Copernicus Marine credentials were not found. ",
      "Run the Copernicus Marine Toolbox login once before using ",
      "`copernicus_load()` non-interactively.",
      call. = FALSE
    )
  }

  message(
    "Copernicus Marine authentication is required.\n",
    "The official Copernicus Marine Toolbox will request your credentials and ",
    "store them in its standard user configuration."
  )

  ok <- tryCatch(
    isTRUE(cm$login()),
    error = function(e) {
      stop(
        "Copernicus Marine authentication failed: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!ok) stop("Copernicus Marine authentication failed.", call. = FALSE)
  invisible(TRUE)
}

# Internal: Copernicus Marine backend.
.copernicus_load_marine <- function(dataset_id, variables,
                                    start_datetime, end_datetime,
                                    xmin, xmax, ymin, ymax,
                                    depth_min, depth_max,
                                    request, output_dir, filename,
                                    split_by, concurrent_processes,
                                    compression, dataset_version,
                                    dataset_part, force, quiet) {

  if (!is.null(split_by)) {
    split_by <- match.arg(split_by, c("variable", "hour", "day", "month", "year"))
  }

  if (!is.null(concurrent_processes)) {
    if (is.null(split_by)) {
      stop("`concurrent_processes` requires `split_by`.", call. = FALSE)
    }
    if (!is.numeric(concurrent_processes) || length(concurrent_processes) != 1 ||
        !is.finite(concurrent_processes) || concurrent_processes < 1 ||
        concurrent_processes != as.integer(concurrent_processes)) {
      stop("`concurrent_processes` must be a positive integer.", call. = FALSE)
    }
  }

  if (!is.null(compression) &&
      (!is.numeric(compression) || length(compression) != 1 ||
       !is.finite(compression) || compression < 0 || compression > 9 ||
       compression != as.integer(compression))) {
    stop("`compression` must be an integer from 0 to 9.", call. = FALSE)
  }

  cm <- .copernicus_marine_module()
  .copernicus_marine_auth(cm)

  common <- list(
    dataset_id = dataset_id,
    output_directory = output_dir,
    overwrite = force,
    skip_existing = !force,
    disable_progress_bar = quiet
  )

  if (!is.null(variables)) common$variables <- as.list(variables)
  if (!is.null(start_datetime)) common$start_datetime <- start_datetime
  if (!is.null(end_datetime)) common$end_datetime <- end_datetime
  if (!is.null(xmin)) common$minimum_longitude <- xmin
  if (!is.null(xmax)) common$maximum_longitude <- xmax
  if (!is.null(ymin)) common$minimum_latitude <- ymin
  if (!is.null(ymax)) common$maximum_latitude <- ymax
  if (!is.null(depth_min)) common$minimum_depth <- depth_min
  if (!is.null(depth_max)) common$maximum_depth <- depth_max
  if (!is.null(filename)) common$output_filename <- filename
  if (!is.null(compression)) common$netcdf_compression_level <- as.integer(compression)
  if (!is.null(dataset_version)) common$dataset_version <- dataset_version
  if (!is.null(dataset_part)) common$dataset_part <- dataset_part

  # Explicit function arguments take precedence over duplicate entries in request.
  args <- utils::modifyList(request, common)

  if (!quiet) {
    message(
      "Requesting Copernicus Marine dataset ", dataset_id,
      if (is.null(split_by)) "..." else paste0(" split by ", split_by, "...")
    )
  }

  result <- tryCatch(
    {
      if (is.null(split_by)) {
        do.call(cm$subset, args)
      } else {
        split_args <- args

        if (split_by == "variable") {
          split_args$on_variables <- TRUE
        } else {
          split_args$on_time <- split_by
        }

        if (!is.null(concurrent_processes)) {
          split_args$concurrent_processes <- as.integer(concurrent_processes)
        }

        do.call(cm$subset_split_on, split_args)
      }
    },
    error = function(e) {
      stop(
        "Copernicus Marine request failed: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  paths <- .copernicus_marine_paths(result)
  missing <- !file.exists(paths) & !dir.exists(paths)

  if (any(missing)) {
    stop(
      "Copernicus Marine returned output path(s) that do not exist:\n  ",
      paste(paths[missing], collapse = "\n  "),
      call. = FALSE
    )
  }

  normalizePath(paths, winslash = "/", mustWork = TRUE)
}

# Internal: validate or interactively configure ECMWF Data Store credentials.
.copernicus_ecmwf_auth <- function() {
  key <- tryCatch(ecmwfr::wf_get_key(), error = function(e) NULL)

  if (is.character(key) && length(key) == 1 && !is.na(key) && nzchar(key)) {
    return(invisible(TRUE))
  }

  if (!interactive()) {
    stop(
      "ECMWF Data Store authentication is required. `ecmwfr >= 2.0.0` uses ",
      "a Personal Access Token (PAT), not the legacy username/password scheme. ",
      "Run `ecmwfr::wf_set_key()` once before using CDS or ADS ",
      "non-interactively.",
      call. = FALSE
    )
  }

  message(
    "ECMWF Data Store authentication is required.\n",
    "CDS and ADS use the same Personal Access Token (PAT).\n",
    "Opening the official CDS API setup page. Log in, copy your PAT, then paste ",
    "it into the authentication dialog.\n",
    "Legacy ERA5/CERRA username/password files are not used."
  )

  try(utils::browseURL("https://cds.climate.copernicus.eu/how-to-api"), silent = TRUE)

  tryCatch(
    ecmwfr::wf_set_key(),
    error = function(e) {
      stop(
        "Could not configure the ECMWF Data Store token: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  key <- tryCatch(ecmwfr::wf_get_key(), error = function(e) NULL)

  if (!is.character(key) || length(key) != 1 || is.na(key) || !nzchar(key)) {
    stop(
      "ECMWF Data Store authentication was not configured successfully. ",
      "Obtain a Personal Access Token from the official CDS API setup page ",
      "and run `ecmwfr::wf_set_key()`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Internal: identify CERRA datasets requiring local geographic cropping.
.copernicus_is_cerra <- function(dataset_id) {
  grepl("^reanalysis-cerra-", dataset_id, ignore.case = TRUE)
}

# Internal: crop CERRA locally while preserving its native Lambert grid.
.copernicus_crop_cerra <- function(input_file, output_file,
                                   xmin, xmax, ymin, ymax,
                                   overwrite = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required to crop CERRA downloads locally.", call. = FALSE)
  }

  r <- tryCatch(
    terra::rast(input_file),
    error = function(e) stop(
      "Could not open the downloaded CERRA file with terra: ",
      conditionMessage(e), call. = FALSE
    )
  )

  bbox <- terra::vect(
    matrix(c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
           ncol = 2, byrow = TRUE),
    type = "polygons",
    crs = "EPSG:4326"
  )

  bbox <- tryCatch(
    terra::project(bbox, terra::crs(r)),
    error = function(e) stop(
      "Could not project the requested bounding box to the CERRA CRS: ",
      conditionMessage(e), call. = FALSE
    )
  )

  cropped <- tryCatch(
    terra::crop(r, bbox),
    error = function(e) stop(
      "Could not crop the CERRA raster to the requested bounding box: ",
      conditionMessage(e), call. = FALSE
    )
  )

  if (terra::ncell(cropped) == 0) {
    stop("The requested bounding box does not overlap the CERRA domain.", call. = FALSE)
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  ext <- tolower(tools::file_ext(output_file))
  filetype <- switch(
    ext,
    grib = "GRIB",
    grb = "GRIB",
    grib2 = "GRIB",
    grb2 = "GRIB",
    nc = "netCDF",
    nc4 = "netCDF",
    tif = "GTiff",
    tiff = "GTiff",
    NULL
  )

  if (is.null(filetype)) {
    stop(
      "Cannot determine an output raster format from `output_file`: ",
      output_file,
      ". Use a supported extension such as .grib, .nc, or .tif.",
      call. = FALSE
    )
  }

  terra::writeRaster(
    cropped,
    output_file,
    overwrite = overwrite,
    filetype = filetype
  )

  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

# Internal: CDS/ADS backend using ecmwfr.
.copernicus_load_ecmwf <- function(source, dataset_id, variables,
                                   start_datetime, end_datetime,
                                   xmin, xmax, ymin, ymax,
                                   depth_min, depth_max,
                                   request, output_dir, filename,
                                   split_by, concurrent_processes,
                                   compression, dataset_version,
                                   dataset_part, force, quiet) {

  if (!requireNamespace("ecmwfr", quietly = TRUE)) {
    stop(
      "Package 'ecmwfr' is required for `source = \"", source, "\"`. ",
      "Install it with install.packages('ecmwfr').",
      call. = FALSE
    )
  }

  .copernicus_ecmwf_auth()

  unsupported <- c(
    split_by = !is.null(split_by),
    concurrent_processes = !is.null(concurrent_processes),
    compression = !is.null(compression),
    dataset_version = !is.null(dataset_version),
    dataset_part = !is.null(dataset_part)
  )

  if (any(unsupported)) {
    stop(
      "The following arguments are specific to Copernicus Marine and cannot ",
      "be used with `source = \"", source, "\"`: ",
      paste(names(unsupported)[unsupported], collapse = ", "),
      ". Supply dataset-specific ECMWF parameters through `request` instead.",
      call. = FALSE
    )
  }

  if (!is.null(depth_min) || !is.null(depth_max)) {
    stop(
      "`depth_min` and `depth_max` are currently supported directly only for ",
      "Copernicus Marine. For ", toupper(source),
      " datasets, specify vertical dimensions such as `pressure_level`, ",
      "`model_level`, or equivalent through `request`.",
      call. = FALSE
    )
  }

  req <- request
  req$dataset_short_name <- dataset_id

  # Explicit common arguments override duplicate request entries.
  if (!is.null(variables)) req$variable <- variables

  if (!is.null(start_datetime)) {
    start_date <- substr(start_datetime, 1, 10)
    end_date <- substr(end_datetime, 1, 10)
    start_time <- substr(start_datetime, 12, 16)
    end_time <- substr(end_datetime, 12, 16)

    req$date <- if (identical(start_date, end_date)) {
      start_date
    } else {
      paste(start_date, end_date, sep = "/")
    }

    if (!identical(start_time, end_time)) {
      stop(
        "`start_datetime` and `end_datetime` contain different clock times. ",
        "For CDS/ADS requests requiring multiple times, supply the time ",
        "selection explicitly through `request = list(time = ...)` and omit ",
        "`start_datetime`/`end_datetime`.",
        call. = FALSE
      )
    }

    req$time <- start_time
  }

  is_cerra <- .copernicus_is_cerra(dataset_id)
  local_crop <- is_cerra && !is.null(xmin)

  if (!is.null(xmin) && !local_crop) req$area <- c(ymax, xmin, ymin, xmax)
  if (!local_crop && !is.null(filename)) req$target <- filename

  if (local_crop && !quiet) {
    message(
      "CERRA dataset detected: downloading the native full domain temporarily ",
      "and cropping locally to the requested bounding box."
    )
  }

  if (!is.null(req$target)) {
    .copernicus_validate_string(req$target, "request$target")
    dest <- file.path(output_dir, req$target)

    if (file.exists(dest) && !force) {
      if (!quiet) message("Cached: ", req$target)
      return(normalizePath(dest, winslash = "/", mustWork = TRUE))
    }

    if (file.exists(dest) && force) unlink(dest, recursive = TRUE, force = TRUE)
  }

  if (!quiet) message("Requesting ", toupper(source), " dataset ", dataset_id, "...")

  if (local_crop) {
    if (is.null(filename)) {
      ext <- if (!is.null(req$data_format) &&
                 tolower(as.character(req$data_format)[1]) %in% c("netcdf", "netcdf4")) ".nc" else ".grib"
      var_tag <- if (!is.null(variables) && length(variables)) {
        gsub("[^A-Za-z0-9._-]+", "_", variables[1])
      } else "cerra"
      filename <- paste0(var_tag, "_cropped", ext)
    }

    .copernicus_validate_string(filename, "filename")
    dest <- file.path(output_dir, filename)

    if (file.exists(dest) && !force) {
      if (!quiet) message("Cached: ", filename)
      return(normalizePath(dest, winslash = "/", mustWork = TRUE))
    }
    if (file.exists(dest) && force) unlink(dest, recursive = TRUE, force = TRUE)

    tmp_dir <- tempfile(pattern = "sharkabc3d_copernicus_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    tmp_ext <- if (!is.null(req$data_format) &&
                  tolower(as.character(req$data_format)[1]) %in% c("netcdf", "netcdf4")) ".nc" else ".grib"
    req$target <- paste0("cerra_full_", format(Sys.time(), "%Y%m%d%H%M%S"), tmp_ext)

    full_path <- tryCatch(
      ecmwfr::wf_request(
        request = req,
        transfer = TRUE,
        path = tmp_dir,
        verbose = !quiet
      ),
      error = function(e) stop(
        toupper(source), " request failed: ", conditionMessage(e), call. = FALSE
      )
    )

    if (!is.character(full_path) || length(full_path) != 1 || !file.exists(full_path)) {
      stop("CERRA request completed but did not return one readable local file.", call. = FALSE)
    }

    cropped <- .copernicus_crop_cerra(
      input_file = full_path,
      output_file = dest,
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      overwrite = force
    )

    if (!quiet) {
      message(sprintf(
        "CERRA local crop complete: %.2f MB full domain -> %.2f MB cropped.",
        file.info(full_path)$size / 1024^2,
        file.info(cropped)$size / 1024^2
      ))
    }

    return(cropped)
  }

  path <- tryCatch(
    ecmwfr::wf_request(
      request = req,
      transfer = TRUE,
      path = output_dir,
      verbose = !quiet
    ),
    error = function(e) {
      stop(
        toupper(source), " request failed: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!is.character(path) || !length(path)) {
    stop(
      toupper(source),
      " request was submitted but no completed local file path was returned.",
      call. = FALSE
    )
  }

  missing <- !file.exists(path)
  if (any(missing)) {
    stop(
      toupper(source), " returned output path(s) that do not exist:\n  ",
      paste(path[missing], collapse = "\n  "),
      call. = FALSE
    )
  }

  normalizePath(path, winslash = "/", mustWork = TRUE)
}
