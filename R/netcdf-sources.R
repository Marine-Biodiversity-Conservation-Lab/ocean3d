# netCDF inputs shared by point, buffer and track extraction. File paths, open
# ncdf4 connections and file tables become one list of sources; each source is
# opened and closed here, and the variable to read is detected when the caller
# doesn't name one.

.is_ncdf4 <- function(x) {
  inherits(x, "ncdf4")
}

# TRUE for a URL such as an OPeNDAP or NCZarr address, which ncdf4 opens
# without a local file.
.is_url <- function(path) grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)

.as_netcdf_sources <- function(nc, file_col = "file") {
  # This helper standardises all accepted input formats into the same
  # internal representation. Each source keeps the original object plus a type
  # flag telling the extraction code whether the file needs to be opened and
  # closed.

  if (is.character(nc)) {
    sources <- lapply(nc, function(x) list(source = x, type = "path"))
    return(sources)
  }

  # ncdf4 objects and data frames are lists too, so both are matched before
  # the list branch below.
  if (.is_ncdf4(nc)) {
    return(list(list(source = nc, type = "connection")))
  }

  if (is.data.frame(nc)) {
    .check_required_cols(nc, file_col, "nc")
    files <- as.character(nc[[file_col]])
    sources <- lapply(files, function(x) list(source = x, type = "path"))
    return(sources)
  }

  # Users often pass nested lists when they group files by variable or by year.
  # Each element is resolved on its own, recursing into nested lists, and the
  # results are joined in order.
  if (is.list(nc)) {
    sources <- lapply(nc, .as_netcdf_sources, file_col = file_col)
    return(Reduce(c, sources, list()))
  }

  stop(
    "Unsupported `nc` input. Use a file path, vector of file paths, ",
    "opened ncdf4 object, list, or data frame with a file column.",
    call. = FALSE
  )
}

.open_source <- function(source) {
  if (identical(source$type, "path")) {
    path <- source$source
    if (!.is_url(path) && !file.exists(path)) {
      rlang::abort(glue::glue("netCDF file does not exist: {path}"),
                   class = "ocean3d_error_netcdf_missing", call = NULL)
    }

    return(list(nc = ncdf4::nc_open(path), close = TRUE))
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
  # salinity. Used by extract_to_point() and extract_to_buffer().
  #
  # Coordinate variables and the variables that CF metadata attaches to them
  # (`bounds`, `climatology`, `coordinates` and `grid_mapping` attributes,
  # e.g. `lat_bnds` or a `crs` variable) are not data, so they are excluded.
  # A grid mapping may use CF's extended form, `"crs: lat lon"`. The
  # detection is deliberately conservative: unless exactly one data variable
  # remains, it stops and asks the user to provide `var` manually, to avoid
  # extracting the wrong variable.

  coordinates <- names(Filter(
    function(x) isTRUE(x$dimvarid$id >= 0), nc$dim
  ))
  metadata <- unlist(lapply(
    unique(c(coordinates, names(nc$var))),
    function(name) {
      atts <- c("bounds", "climatology", "coordinates", "grid_mapping")
      unlist(lapply(atts, function(att) {
        value <- .read_attribute(nc, name, att)
        if (is.character(value)) {
          sub(":$", "", unlist(strsplit(value, "[[:space:]]+")))
        }
      }))
    }
  ))
  candidates <- setdiff(names(nc$var), c(names(nc$dim), metadata))

  if (length(candidates) == 1L) {
    return(candidates)
  }

  found <- if (length(candidates)) toString(candidates) else "none"
  file <- basename(nc$filename)
  rlang::abort(glue::glue(
    "Automatic variable detection requires exactly one data variable in ",
    "each netCDF file (coordinate and bounds variables are ignored).\n",
    "File: {file}\n",
    "Data variables found: {found}.\n",
    "Please specify `var` manually."
  ), class = "ocean3d_error_var_detection", call = NULL)
}
