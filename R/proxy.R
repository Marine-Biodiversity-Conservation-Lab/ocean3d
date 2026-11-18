# The core object (#113, SPEC.md "Direction: a hybrid"): a description of one
# gridded netCDF variable, with the CF metadata applied once when it is built
# and no values read. Axes reuse the shared netCDF helpers in R/netcdf-*.R.

#' @importFrom methods setOldClass show
NULL

setOldClass(c("CFTime", "R6"))

#' A gridded dataset described but not loaded (experimental)
#'
#' A `SpatProxy` describes one variable of a netCDF file without reading its
#' values: where the data are, how they are stored, and each axis in a fixed
#' form. In terra, "proxy" means the same thing: a `SpatVectorProxy`
#' describes vector data on disk. Build one with [as_proxy()].
#'
#' A `SpatProxy` is the planned core data object of ocean3d (SPEC.md, "Data
#' model"). Later versions will read one time step of it into memory as a
#' [SpatVoxel-class].
#'
#' @section Experimental:
#' `SpatProxy` and [as_proxy()] are experimental. Their slots and behaviour
#' may change without notice while the core object is built, and no other
#' function accepts a `SpatProxy` yet.
#'
#' @slot source A list: `type` (`"netcdf"`) and `path`, the full path of the
#'   file. A path rather than an open connection, so the object can be saved
#'   with [base::saveRDS()].
#' @slot variable A list: `name`, `units`, `long_name`, `dims` (the
#'   variable's dimension names in ncdf4 order), `dim_lengths`, `axis_dims`
#'   (the dimension of each axis), `singleton_dims` (other dimensions of
#'   length 1) and `chunksizes` (`NULL` for unchunked storage).
#' @slot axes A list with one element per axis, `NULL` for an axis the
#'   variable doesn't have. `lon`, `lat` and `depth` hold `dim`, `centres`
#'   and cell edges `lower` and `upper`, in storage order, and `edges`: where
#'   the edges come from (`"bounds"`, `"inferred"`, or `"none"` for a single
#'   cell without bounds). `depth` is in metres positive down and also has
#'   `storage` (`"shallow_first"` or `"deep_first"`). `time` holds `dim`,
#'   `cftime` (a CFtime object, with the file's bounds attached when it has
#'   them) and `edges` (`"bounds"` or `"none"`).
#' @slot crs The coordinate reference system, as text.
#'
#' @seealso [as_proxy()]
#' @export
setClass("SpatProxy", slots = c(
  source = "list", variable = "list", axes = "list", crs = "character"
))

setValidity("SpatProxy", function(object) {
  axes <- object@axes
  lengths <- object@variable$dim_lengths
  if (is.null(axes$lon) || is.null(axes$lat)) {
    return("a SpatProxy needs a longitude and a latitude axis")
  }
  gridded <- Filter(Negate(is.null), axes[c("lon", "lat", "depth")])
  for (name in names(gridded)) {
    axis <- gridded[[name]]
    n <- lengths[[axis$dim]]
    if (length(axis$centres) != n || length(axis$lower) != n ||
          length(axis$upper) != n || any(axis$lower > axis$upper)) {
      return(glue::glue(
        "the {name} axis needs one centre, lower and upper edge per cell, ",
        "with lower <= upper"
      ))
    }
  }
  depth <- axes$depth$centres
  if (any(depth < 0) || anyDuplicated(depth)) {
    return("depths must be distinct metres, positive down")
  }
  time <- axes$time$cftime
  if (!is.null(axes$time) &&
        (!inherits(time, "CFTime") || any(diff(time$offsets) <= 0))) {
    return("the time axis must be a strictly increasing CFtime object")
  }
  TRUE
})

#' Describe a gridded netCDF variable without reading its values
#' (experimental)
#'
#' `as_proxy()` opens one netCDF file, reads the metadata of one variable and
#' its coordinates, and returns a [SpatProxy-class]. The CF metadata is
#' applied here, once, so later steps can rely on each axis being in a fixed
#' form. No values of the variable are read.
#'
#' * **Longitude and latitude:** cell edges come from the coordinate's CF
#'   `bounds` variable, or halfway between neighbouring centres.
#' * **Depth:** converted to metres positive down from the coordinate's
#'   `units` and `positive` attributes. Edges come from the CF bounds, or
#'   halfway between centres, cut at the surface. A depth axis stored deep to
#'   shallow is kept in that order and recorded as such.
#' * **Time:** decoded with CFtime. Bounds are attached when the file has
#'   them; otherwise the axis has none (they aren't inferred). The time axis
#'   must be strictly increasing; a backwards or out-of-order axis is an
#'   error that suggests sorting the file, for example with CDO's
#'   `sorttimestamp`. Only the `standard`, `gregorian` and
#'   `proleptic_gregorian` calendars are supported, and climatological time
#'   axes aren't supported yet.
#' * **CRS:** EPSG:4326, unless the variable's `grid_mapping` gives another
#'   geographic CRS as WKT (`crs_wkt` or `spatial_ref`).
#'
#' Regular longitude-latitude grids are supported; curvilinear and projected
#' grids are rejected.
#'
#' @inheritSection SpatProxy-class Experimental
#'
#' @param x A netCDF file path, or a connection opened with
#'   [ncdf4::nc_open()] (netCDF model). A connection is left open.
#' @param var Name of the variable. `NULL` (the default) detects it when the
#'   file holds exactly one data variable.
#' @param lon_dim,lat_dim,depth_dim,time_dim Names of the dimensions for each
#'   axis. `NULL` (the default) detects them from CF metadata (`axis`,
#'   `standard_name`, `units`) or common names.
#'
#' @return A [SpatProxy-class].
#' @seealso [SpatProxy-class]
#' @examples
#' x <- as_proxy(system.file("extdata", "example_3d.nc", package = "ocean3d"))
#' x
#'
#' # Each axis in its fixed form: depth in metres, time as a CFtime object.
#' x@axes$depth$centres
#' x@axes$time$cftime
#' @export
as_proxy <- function(x, var = NULL, lon_dim = NULL, lat_dim = NULL,
                     depth_dim = NULL, time_dim = NULL) {
  if (!.is_ncdf4(x) && !(is.character(x) && length(x) == 1L && !is.na(x))) {
    rlang::abort(
      "`x` must be one netCDF file path or an open ncdf4 connection.",
      class = "ocean3d_error_proxy_input", call = NULL
    )
  }
  type <- if (.is_ncdf4(x)) "connection" else "path"
  opened <- .open_source(list(source = x, type = type))
  on.exit(.close_source(opened), add = TRUE)
  nc <- opened$nc
  if (is.null(var)) var <- .detect_netcdf_var(nc)
  dims <- .detect_dims(nc, var, lon_dim, lat_dim, depth_dim, time_dim)
  info <- nc$var[[var]]

  # Where an axis's edges come from, from the shared readers' `basis`.
  edges_from <- function(basis) {
    switch(basis, explicit_bounds = "bounds",
           singleton_coordinate = "none", "inferred")
  }
  horizontal <- lapply(c(lon = "lon", lat = "lat"), function(axis) {
    values <- as.vector(.get_dim_values(nc, dims[[axis]]))
    out <- .prepare_spatial_axis(nc, dims[[axis]], values, axis == "lon")
    c(list(dim = dims[[axis]]), out[c("centres", "lower", "upper")],
      list(edges = edges_from(out$basis)))
  })

  # The file's full path, or a URL (e.g. OPeNDAP) as it is. A connection's
  # relative path resolves against the working directory, so stop if no such
  # file is there rather than record a wrong path.
  path <- nc$filename
  if (!grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)) {
    if (!file.exists(path)) {
      rlang::abort(glue::glue(
        "Can't find the file `{path}` that this connection was opened ",
        "from. Open it with its full path, or pass the path to as_proxy()."
      ), class = "ocean3d_error_proxy_source", call = NULL)
    }
    path <- normalizePath(path)
  }
  methods::new(
    "SpatProxy",
    source = list(type = "netcdf", path = path),
    variable = list(
      name = var, units = info$units, long_name = info$longname,
      dims = dims$var_dims,
      dim_lengths = stats::setNames(as.integer(info$varsize), dims$var_dims),
      axis_dims = dims[c("lon", "lat", "depth", "time")],
      singleton_dims = dims$singleton_dims,
      # ncdf4 leaves `chunksizes` unset for contiguous storage (storage 1).
      chunksizes = if (identical(as.integer(info$storage), 2L)) {
        stats::setNames(as.integer(info$chunksizes), dims$var_dims)
      }
    ),
    axes = c(horizontal, list(
      depth = if (!is.null(dims$depth)) {
        out <- .read_interval_axis(nc, dims$depth, "depth",
                                   infer_bounds = TRUE, allow_singleton = TRUE)
        n <- length(out$centres)
        c(list(dim = dims$depth), out[c("centres", "lower", "upper")], list(
          edges = edges_from(out$basis),
          storage = if (out$centres[1L] > out$centres[n]) "deep_first" else
            "shallow_first"
        ))
      },
      time = if (!is.null(dims$time)) .proxy_time_axis(nc, dims$time, path)
    )),
    crs = .proxy_crs(nc, var)
  )
}

# The time axis as a CFtime object, with the file's bounds attached when it
# has them. Time must run strictly forwards (SPEC.md, "Time runs forwards
# only"); the core object doesn't reorder it.
.proxy_time_axis <- function(nc, dim, path) {
  file <- basename(path)
  .check_climatology(nc, dim, file)
  offsets <- as.vector(.get_dim_values(nc, dim))
  units <- .read_time_attribute(nc, dim, "units")
  usable <- is.numeric(offsets) && all(is.finite(offsets)) &&
    is.character(units) && length(units) == 1L && nzchar(trimws(units))
  # Without usable units or offsets, the axis builder gives the error.
  if (!usable) .netcdf_time_axis(nc, dim)
  steps <- diff(offsets)
  if (any(steps == 0)) {
    value <- offsets[which(steps == 0)[1L] + 1L]
    .abort_repeated_time(
      glue::glue("Time step {value} ({units}) appears ",
                 "{sum(offsets == value)} times"),
      file,
      paste("ocean3d needs each time step once; remove or combine",
            "the repeated steps first.")
    )
  }
  if (any(steps < 0)) {
    order <- if (all(steps < 0)) "runs backwards" else "is out of order"
    rlang::abort(glue::glue(
      "The time axis `{dim}` in {file} {order}. ocean3d needs time steps ",
      "in increasing order. Write a sorted copy of the file first, for ",
      "example with CDO: `cdo sorttimestamp {file} sorted.nc`."
    ), class = "ocean3d_error_time_order", call = NULL)
  }
  bounds <- .read_bounds(nc, dim, length(offsets))
  if (is.null(bounds)) {
    return(list(dim = dim, cftime = .netcdf_time_axis(nc, dim),
                edges = "none"))
  }
  # Checked as buffer extraction checks them: the bounds calendar, each
  # interval containing its step (in the axis builder), then positive widths
  # and no overlaps. In the file's units, as CFtime can't date every step
  # (e.g. half a month).
  .check_time_bounds_calendar(nc, dim, bounds)
  cftime <- .netcdf_time_axis(nc, dim, bounds$values)
  .axis_intervals(offsets, bounds$values, "time", dim, file = file)
  list(dim = dim, cftime = cftime, edges = "bounds")
}

# EPSG:4326, or the geographic CRS the variable's grid mapping gives as WKT.
.proxy_crs <- function(nc, var) {
  mapping <- .read_time_attribute(nc, var, "grid_mapping")
  if (is.null(mapping) || !mapping %in% names(nc$var)) return("EPSG:4326")
  read <- function(name) {
    att <- ncdf4::ncatt_get(nc, mapping, name)
    if (isTRUE(att$hasatt)) att$value
  }
  wkt <- read("crs_wkt")
  if (is.null(wkt)) wkt <- read("spatial_ref")
  name <- read("grid_mapping_name")
  geographic <- if (is.null(wkt)) {
    identical(name, "latitude_longitude")
  } else {
    isTRUE(tryCatch(terra::is.lonlat(wkt), error = function(e) FALSE))
  }
  if (!geographic) {
    rlang::abort(glue::glue(
      "Variable `{var}` uses the grid mapping `{mapping}`, which isn't a ",
      "geographic longitude-latitude CRS. Only regular longitude-latitude ",
      "grids are supported."
    ), class = "ocean3d_error_grid_mapping", call = NULL)
  }
  if (is.null(wkt)) "EPSG:4326" else wkt
}

#' @rdname SpatProxy-class
#' @param object A `SpatProxy`.
#' @export
setMethod("show", "SpatProxy", function(object) {
  axes <- object@axes
  variable <- object@variable
  extent <- function(x) paste(signif(range(x), 6), collapse = " to ")
  edges <- function(axis) {
    switch(axis$edges, bounds = "edges from bounds",
           inferred = "edges inferred", none = "no edges")
  }
  # Sizes are shown in the fixed axis order, whatever the storage order.
  dims <- unlist(variable$axis_dims)
  sizes <- function(x) paste(x[dims], collapse = " x ")
  chunks <- if (is.null(variable$chunksizes)) "" else
    paste0(", chunks ", sizes(variable$chunksizes))
  time_bounds <- if (identical(axes$time$edges, "bounds")) "with" else "no"
  lines <- c(
    class = "SpatProxy (experimental)",
    source = basename(object@source$path),
    variable = glue::glue(
      "{variable$name} ({variable$units}){chunks}"
    ),
    dimensions = glue::glue(
      "{sizes(variable$dim_lengths)} ({paste(names(dims), collapse = ' x ')})"
    ),
    lon = glue::glue(
      "{extent(axes$lon$centres)}, {length(axes$lon$centres)} cells, ",
      "{edges(axes$lon)}"
    ),
    lat = glue::glue(
      "{extent(axes$lat$centres)}, {length(axes$lat$centres)} cells, ",
      "{edges(axes$lat)}"
    ),
    depth = if (is.null(axes$depth)) "none" else glue::glue(
      "{extent(axes$depth$centres)} m, {length(axes$depth$centres)} levels, ",
      "{sub('_', ' ', axes$depth$storage)}, {edges(axes$depth)}"
    ),
    time = if (is.null(axes$time)) "none" else glue::glue(
      "{paste(axes$time$cftime$range(), collapse = ' to ')}, ",
      "{length(axes$time$cftime)} steps, {axes$time$cftime$cal$name} ",
      "calendar, {time_bounds} bounds"
    ),
    crs = if (grepl("^EPSG:", object@crs)) object@crs else
      terra::crs(object@crs, describe = TRUE)$name
  )
  cat(sprintf("%-10s : %s", names(lines), lines), sep = "\n")
  invisible(object)
})
