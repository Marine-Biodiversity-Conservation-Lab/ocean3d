#' WOA cache directory
#'
#' Returns the path to the package's persistent cache directory for downloaded
#' WOA NetCDF files. Uses [tools::R_user_dir()] so the location survives across
#' sessions and follows platform conventions.
#'
#' @returns Character. Path to cache directory (created if missing).
#' @examples
#' # Redirect the cache into a temporary directory so this example does not
#' # write to user filespace. In normal use you would skip this step and let
#' # the cache live in its platform-conventional location.
#' old <- Sys.getenv("R_USER_CACHE_DIR", unset = NA)
#' Sys.setenv(R_USER_CACHE_DIR = tempdir())
#'
#' # The directory is created if it does not already exist
#' cache <- woa_cache_dir()
#' dir.exists(cache)
#'
#' # This is where woa_download() writes when `output_dir` is not supplied
#' list.files(cache, pattern = "\\.nc$")
#'
#' if (is.na(old)) Sys.unsetenv("R_USER_CACHE_DIR") else
#'   Sys.setenv(R_USER_CACHE_DIR = old)
#'
#' @export
woa_cache_dir <- function() {
  path <- file.path(tools::R_user_dir("ocean3d", which = "cache"), "woa")
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

#' Clear the WOA cache
#'
#' Remove all cached WOA NetCDF files.
#'
#' @param confirm Logical. Require interactive confirmation. Default `TRUE`.
#'
#' @returns Invisibly, `TRUE` on success.
#' @examples
#' # Redirect the cache into a temporary directory so this example does not
#' # delete real downloads. In normal use you would skip this step.
#' old <- Sys.getenv("R_USER_CACHE_DIR", unset = NA)
#' Sys.setenv(R_USER_CACHE_DIR = tempdir())
#'
#' cache <- woa_cache_dir()
#' file.create(file.path(cache, "woa23_decav_t00_01.nc"))
#' list.files(cache)
#'
#' # confirm = FALSE deletes without prompting (e.g. from a script)
#' woa_cache_clear(confirm = FALSE)
#' dir.exists(cache)
#'
#' # The cache is recreated empty by the next call that needs it
#' woa_cache_dir()
#'
#' if (is.na(old)) Sys.unsetenv("R_USER_CACHE_DIR") else
#'   Sys.setenv(R_USER_CACHE_DIR = old)
#'
#' \dontrun{
#' # Called with no arguments it asks for confirmation when interactive
#' woa_cache_clear()
#' }
#'
#' @export
woa_cache_clear <- function(confirm = TRUE) {
  path <- woa_cache_dir()
  if (confirm && interactive()) {
    ans <- readline(sprintf("Delete all cached WOA files in %s? [y/N]: ", path))
    if (!tolower(ans) %in% c("y", "yes")) {
      message("Cancelled.")
      return(invisible(FALSE))
    }
  }
  unlink(path, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

# Internal: accepted variable names and the URL/filename components each maps
# to. `code` is the one-letter code used both in WOA filenames and in the
# `{code}_{field}_depth={value}` layer names inside the NetCDFs. Note that
# codes are case-sensitive: "o"/"O" and "i"/"I" are different variables.
.woa_variables <- function() {
  data.frame(
    name = c("temperature", "salinity", "dissolved_oxygen",
             "oxygen_saturation", "AOU", "nitrate", "phosphate",
             "silicate", "density"),
    dir = c("temperature", "salinity", "oxygen",
            "o2sat", "AOU", "nitrate", "phosphate",
            "silicate", "density"),
    code = c("t", "s", "o",
             "O", "A", "n", "p",
             "i", "I"),
    decade = c("decav", "decav", "all",
               "all", "all", "all", "all",
               "all", "decav"),
    stringsAsFactors = FALSE
  )
}

# Internal: accepted statistical fields. `code` is the code used in WOA
# filenames and in the `{variable}_{code}_depth={value}` layer names inside
# the NetCDFs. See the WOA 2023 Product Documentation:
# https://repository.library.noaa.gov/view/noaa/70581
.woa_fields <- function() {
  data.frame(
    code = c("an", "mn", "dd",
             "sd", "se", "oa",
             "gp", "sdo", "sea"),
    description = c("objectively analyzed climatology",
                    "statistical mean",
                    "number of observations",
                    "standard deviation",
                    "standard error",
                    "mean minus climatology",
                    "number of mean values within radius of influence",
                    "objectively analyzed standard deviation",
                    "standard error of the analysis"),
    stringsAsFactors = FALSE
  )
}

# Internal: map variable name to WOA URL/filename components
.woa_variable_spec <- function(variable) {
  specs <- .woa_variables()
  if (!variable %in% specs$name) {
    stop("Unknown variable '", variable,
         "'. Supported: ", toString(specs$name))
  }
  specs[specs$name == variable,]
}

# Internal: map resolution to URL dir + filename code
.woa_resolution_spec <- function(resolution) {
  res <- as.character(resolution)
  specs <- list(
    "0.25" = list(dir = "0.25", code = "04"),
    "1"    = list(dir = "1.00", code = "01"),
    "1.00" = list(dir = "1.00", code = "01"),
    "5"    = list(dir = "5deg", code = "5d")
  )
  if (!res %in% names(specs)) {
    stop("Unknown resolution '", resolution, "'. Supported: 0.25, 1, 5.")
  }
  specs[[res]]
}

# Internal: expand period into one or more 2-digit period codes
.woa_period_codes <- function(period) {
  if (is.character(period) && length(period) == 1) {
    if (period == "annual") return("00")
    if (period == "monthly") return(sprintf("%02d", 1:12))
    if (period == "seasonal") return(sprintf("%02d", 13:16))
  }
  if (is.numeric(period)) {
    if (any(period < 0 | period > 16)) {
      stop("Numeric period must be 0 (annual), 1:12 (monthly), ",
           "or 13:16 (seasonal).")
    }
    return(sprintf("%02d", period))
  }
  stop("period must be 'annual', 'monthly', 'seasonal', or a numeric 0-16.")
}

# Internal: one-time interactive consent for writing to the persistent cache.
# CRAN Repository Policy requires explicit user confirmation before writing to
# user filespace outside tempdir(). A sentinel file records consent so the
# prompt only fires once per cache directory.
.woa_cache_consent <- function(cache_dir) {
  sentinel <- file.path(cache_dir, ".consent")
  if (file.exists(sentinel)) return(invisible(TRUE))

  msg <- paste0(
    "ocean3d will cache downloaded WOA NetCDF files in:\n  ",
    cache_dir, "\n",
    "A full WOA set (temperature + oxygen, annual + monthly, 0.25 deg) ",
    "can exceed 10 GB.\n",
    "Pass `output_dir` to use a different location, or call ",
    "woa_cache_clear() to reclaim space later."
  )

  if (!interactive()) {
    stop(msg,
         "\nNon-interactive session: pass `output_dir` explicitly, or run ",
         "woa_download() once interactively to record consent.")
  }

  message(msg)
  ans <- readline("Proceed with caching here? [y/N]: ")
  if (!tolower(ans) %in% c("y", "yes")) {
    stop("Cache consent declined.")
  }
  file.create(sentinel)
  invisible(TRUE)
}

# Internal: fetch a single URL to `dest` via curl (no artificial timeout,
# resumable). Partial files are removed on failure so they don't masquerade
# as cached on the next call.
.woa_fetch <- function(url, dest, quiet = FALSE) {
  cleanup <- function() if (file.exists(dest)) file.remove(dest)

  tryCatch(
    curl::curl_download(url, destfile = dest, mode = "wb", quiet = quiet),
    error = function(e) {
      cleanup()
      stop("Download failed for ", url, ": ", conditionMessage(e),
           call. = FALSE)
    },
    interrupt = function(e) {
      cleanup()
      stop("Download interrupted.", call. = FALSE)
    }
  )

  if (!file.exists(dest) || file.size(dest) == 0) {
    cleanup()
    stop("Download produced empty file: ", url, call. = FALSE)
  }
  invisible(dest)
}

# Internal: HEAD request for remote file size in MB. Returns NA if the server
# doesn't report Content-Length (unusual for static NCEI files).
.woa_remote_size_mb <- function(url) {
  tryCatch({
    h <- curl::new_handle()
    curl::handle_setopt(h, nobody = TRUE)
    resp <- curl::curl_fetch_memory(url, handle = h)
    hdrs <- curl::parse_headers_list(resp$headers)
    len <- as.numeric(hdrs[["content-length"]])
    if (is.na(len) || len <= 0) NA_real_ else len / 1024^2
  }, error = function(e) NA_real_)
}

#' Download a WOA NetCDF file (with caching)
#'
#' Download World Ocean Atlas 2023 NetCDF files from the NCEI THREDDS server.
#' Files are cached in [woa_cache_dir()] (or a user-supplied `output_dir`) and
#' skipped on subsequent calls unless `force = TRUE`.
#'
#' On first use of the persistent cache in an interactive session, the function
#' prompts for consent to write to user filespace (per CRAN Repository Policy).
#' When the `curl` package is available, the estimated size of each remote
#' file is reported before downloading. A full set of WOA files can be many
#' gigabytes — supply `output_dir` to direct the cache elsewhere (HPC scratch,
#' external drive, etc.).
#'
#' @param variable Character. One of `"temperature"`, `"salinity"`,
#'   `"dissolved_oxygen"`, `"oxygen_saturation"`, `"AOU"`, `"nitrate"`,
#'   `"phosphate"`, `"silicate"`, `"density"`.
#' @param period Character or numeric. `"annual"` (default), `"monthly"`
#'   (all 12 months), `"seasonal"` (4 seasons), or a numeric vector where 0 =
#'   annual, 1:12 = monthly, 13:16 = seasonal.
#' @param resolution Character or numeric. `"0.25"` (default), `"1"`, or `"5"`
#'   degrees.
#' @param decade Character. Decade code (e.g. `"decav"`, `"all"`). Defaults to
#'   the canonical decade for the variable.
#' @param output_dir Character. Destination directory. Defaults to
#'   [woa_cache_dir()].
#' @param force Logical. If `TRUE`, re-download even if the file exists.
#'   Default `FALSE`.
#' @param quiet Logical. Suppress download progress. Default `FALSE`.
#'
#' @returns Character vector of paths to downloaded .nc files.
#' @examples
#' \dontrun{
#' # Annual temperature climatology at 1-degree resolution, cached in a
#' # temporary directory so nothing is written to user filespace
#' tmp <- file.path(tempdir(), "woa")
#' f <- woa_download("temperature", period = "annual", resolution = "1",
#'                   output_dir = tmp)
#' woa_load_nc(f)
#'
#' # All 12 monthly files, downloaded in parallel, into the default cache
#' monthly <- woa_download("temperature", period = "monthly",
#'                         resolution = "1")
#' length(monthly)
#'
#' # Specific periods: 0 = annual, 1:12 = monthly, 13:16 = seasonal
#' winter <- woa_download("dissolved_oxygen", period = 13, resolution = "1")
#'
#' # Other variables use their own canonical decade codes. Note that
#' # 0.25-degree files are large (~1 GB each for salinity).
#' sal <- woa_download("salinity", period = "annual", resolution = "0.25")
#'
#' # Re-download a file that is already cached
#' f <- woa_download("temperature", period = "annual", resolution = "1",
#'                   output_dir = tmp, force = TRUE, quiet = TRUE)
#' }
#'
#' @export
woa_download <- function(variable,
                         period = "annual",
                         resolution = "0.25",
                         decade = NULL,
                         output_dir = NULL,
                         force = FALSE,
                         quiet = FALSE) {
  vspec <- .woa_variable_spec(variable)
  rspec <- .woa_resolution_spec(resolution)
  period_codes <- .woa_period_codes(period)
  if (is.null(decade)) decade <- vspec$decade

  using_default_cache <- is.null(output_dir)
  if (using_default_cache) output_dir <- woa_cache_dir()
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Only gate consent on the persistent default cache — user-supplied paths
  # are an explicit opt-in already.
  if (using_default_cache) .woa_cache_consent(output_dir)
  base <- "https://www.ncei.noaa.gov/thredds-ocean/fileServer/woa23/DATA"
  plan <- lapply(period_codes, function(pc) {
    filename <- sprintf("woa23_%s_%s%s_%s.nc",
                        decade, vspec$code, pc, rspec$code)
    list(
      url = sprintf("%s/%s/netcdf/%s/%s/%s",
                    base, vspec$dir, decade, rspec$dir, filename),
      dest = file.path(output_dir, filename),
      filename = filename
    )
  })

  # Partition into already-cached vs. needs-download. Re-download any cached
  # file whose size disagrees with the remote (i.e., truncated from a prior
  # interrupted run).
  to_fetch <- Filter(Negate(is.null), lapply(plan, function(p) {
    if (!force && file.exists(p$dest) && file.size(p$dest) > 0) {
      remote_mb <- .woa_remote_size_mb(p$url)
      local_mb <- file.size(p$dest) / 1024^2
      if (is.na(remote_mb) || abs(local_mb - remote_mb) < 1) {
        if (!quiet) message("Cached: ", p$filename)
        return(NULL)
      }
      if (!quiet) {
        message(sprintf(
          "Cached file is truncated (%.1f / %.1f MB); re-downloading: %s",
          local_mb, remote_mb, p$filename
        ))
      }
    }
    p
  }))

  if (length(to_fetch) == 1) {
    p <- to_fetch[[1]]
    if (!quiet) {
      size_mb <- .woa_remote_size_mb(p$url)
      size_str <- if (is.na(size_mb)) "unknown size" else
        sprintf("%.1f MB", size_mb)
      message("Downloading (", size_str, "): ", p$filename)
    }
    .woa_fetch(p$url, p$dest, quiet = quiet)
  } else if (length(to_fetch) > 1) {
    .woa_fetch_parallel(to_fetch, quiet = quiet)
  }

  unname(vapply(plan, `[[`, character(1), "dest"))
}

# Internal: download multiple files concurrently via curl::multi_download.
# NCEI THREDDS throttles per-connection rather than per-client, so several
# parallel streams typically give 3-5x aggregate throughput over serial.
.woa_fetch_parallel <- function(items, quiet = FALSE) {
  urls <- vapply(items, `[[`, character(1), "url")
  dests <- vapply(items, `[[`, character(1), "dest")
  filenames <- vapply(items, `[[`, character(1), "filename")

  if (!quiet) {
    message("Downloading ", length(items), " files in parallel:")
    for (fn in filenames) message("  - ", fn)
  }

  cleanup <- function() {
    for (d in dests) if (file.exists(d)) file.remove(d)
  }

  res <- tryCatch(
    curl::multi_download(urls, destfiles = dests,
                         resume = TRUE, progress = !quiet),
    error = function(e) {
      cleanup()
      stop("Parallel download failed: ", conditionMessage(e), call. = FALSE)
    },
    interrupt = function(e) {
      cleanup()
      stop("Download interrupted.", call. = FALSE)
    }
  )

  failed <- is.na(res$success) | !res$success
  if (any(failed)) {
    for (d in dests[failed]) if (file.exists(d)) file.remove(d)
    stop(
      "Download failed for ", sum(failed), " file(s):\n  ",
      paste(filenames[failed], collapse = "\n  "),
      call. = FALSE
    )
  }
  invisible(dests)
}

#' Load a WOA NetCDF file
#'
#' Load a WOA .nc file and select layers for a given statistical field.
#' Returns a SpatRaster with layer names already following the
#' `{variable}_depth={value}` convention used throughout this package
#' (native to WOA NetCDFs).
#'
#' @param file_path Character. Path to a WOA .nc file.
#' @param field Character. Statistical field to select. One of:
#'   `"an"` (objectively analyzed climatology), `"mn"` (statistical mean),
#'   `"dd"` (number of observations), `"sd"` (standard deviation),
#'   `"se"` (standard error), `"oa"` (mean minus climatology),
#'   `"gp"` (number of mean values within radius of influence),
#'   `"sdo"` (objectively analyzed standard deviation),
#'   `"sea"` (standard error of the analysis). See the WOA 2023 Product
#'   Documentation:
#'   <https://repository.library.noaa.gov/view/noaa/70581>
#'
#' @returns SpatRaster with standardized depth layer names.
#' @examples
#' # Build a small WOA-style NetCDF to demonstrate field selection.
#' # Real files come from woa_download().
#' nc <- file.path(tempdir(), "woa_demo.nc")
#' demo <- terra::rast(nrows = 2, ncols = 2, nlyrs = 4,
#'                     xmin = 0, xmax = 2, ymin = 0, ymax = 2,
#'                     crs = "EPSG:4326")
#' names(demo) <- c("t_an_depth=0", "t_an_depth=100",
#'                  "t_sd_depth=0", "t_sd_depth=100")
#' terra::values(demo) <- matrix(1:16, nrow = 4)
#' terra::writeCDF(demo, nc, overwrite = TRUE, split = TRUE)
#'
#' # Objectively analyzed climatology (the default field)
#' temp <- woa_load_nc(nc)
#' names(temp)
#'
#' # Standard deviation layers from the same file
#' names(woa_load_nc(nc, field = "sd"))
#'
#' unlink(nc)
#'
#' \dontrun{
#' # Typical use on a downloaded file
#' f <- woa_download("temperature", period = "annual", resolution = "1")
#' temp <- woa_load_nc(f, field = "an")
#' terra::plot(temp[["t_an_depth=0"]])
#' }
#'
#' @export
woa_load_nc <- function(file_path, field = "an") {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  r <- terra::rast(file_path)
  
  available_fields <- .woa_fields()$code
  if (!field %in% available_fields) {
    stop(
      "field must be one of: ",
      toString(available_fields),
      ". See https://repository.library.noaa.gov/view/noaa/70581"
    )
  }

  # WOA layers take the form of v_field_depth=123
  # v is the abbreviated variable, which has to be one character
  var_abr <- r[[1]] %>%
    names() %>%
    stringr::str_split_i(pattern = stringr::fixed("_"), 1)

  # check that layer names follow convention
  layer_name_split <- r[[1]] %>%
    names() %>%
    stringr::str_split(pattern = stringr::fixed("_"))
  available_vars <- .woa_variables()$code
  if (!layer_name_split[[1]][1] %in% available_vars) {
    stop(
      "check raster layer names. For WOA data, oceanographic variable ",
      "one-letter code must be one of: ",
      toString(available_vars),
      ". See https://repository.library.noaa.gov/view/noaa/70581"
    )
  }
  if (!layer_name_split[[1]][2] %in% available_fields) {
    stop(
      "check raster layer names. For WOA data, field code must be one of: ",
      toString(available_fields),
      ". See https://repository.library.noaa.gov/view/noaa/70581"
    )
  }

  field_pattern <- paste0(var_abr, "_", field, "_depth=")
  selected_names <- r %>%
    names() %>%
    stringr::str_subset(field_pattern)

  # `field` is a valid WOA code, but this particular file may not carry it:
  # NCEI ships some fields (e.g. "sea", "sdo") only for certain variables and
  # periods. Report what the file does contain instead of letting the empty
  # subset fail inside terra with "is.numeric(i) is not TRUE".
  if (length(selected_names) == 0) {
    present <- names(r) %>%
      stringr::str_split_i(pattern = stringr::fixed("_"), 2) %>%
      unique() %>%
      sort()
    stop(
      "No layers for field '", field, "' in ", basename(file_path),
      ". Fields present in this file: ", toString(present),
      "."
    )
  }

  r[[selected_names]]
}
