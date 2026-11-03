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
#'   length 1), `packing` (`storage` type, `scale_factor`, `add_offset`,
#'   `fill_value` and `missing_value`, `NULL` when absent) and `chunksizes`
#'   (`NULL` for unchunked storage).
#' @slot axes A list with one element per axis, `NULL` for an axis the
#'   variable doesn't have. `lon`, `lat` and `depth` hold `dim`, `centres`
#'   and cell edges `lower` and `upper`, in storage order, and `edges`: where
#'   the edges come from (`"bounds"`, `"inferred"`, or `"none"` for a single
#'   cell without bounds). `lon` also has `convention` (`"0_360"` or
#'   `"-180_180"`). `depth` is in metres positive down and also has
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
#'   `bounds` variable, or halfway between neighbouring centres. The
#'   longitude convention (0 to 360 or -180 to 180) is recorded.
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

  horizontal <- lapply(c(lon = "lon", lat = "lat"), function(axis) {
    values <- as.vector(.get_dim_values(nc, dims[[axis]]))
    out <- .prepare_spatial_axis(nc, dims[[axis]], values, axis == "lon")
    edges <- switch(out$basis, explicit_bounds = "bounds",
                    singleton_coordinate = "none", "inferred")
    c(list(dim = dims[[axis]]), out[c("centres", "lower", "upper")],
      list(edges = edges),
      if (axis == "lon") {
        list(convention = if (any(values > 180)) "0_360" else "-180_180")
      })
  })

  attribute <- function(name) .read_time_attribute(nc, var, name)
  path <- normalizePath(nc$filename)
  methods::new(
    "SpatProxy",
    source = list(type = "netcdf", path = path),
    variable = list(
      name = var, units = info$units, long_name = info$longname,
      dims = dims$var_dims,
      dim_lengths = stats::setNames(as.integer(info$varsize), dims$var_dims),
      axis_dims = dims[c("lon", "lat", "depth", "time")],
      singleton_dims = dims$singleton_dims,
      packing = list(
        storage = info$prec, scale_factor = attribute("scale_factor"),
        add_offset = attribute("add_offset"),
        fill_value = attribute("_FillValue"),
        missing_value = attribute("missing_value")
      ),
      chunksizes = if (!anyNA(info$chunksizes)) {
        stats::setNames(as.integer(info$chunksizes), dims$var_dims)
      }
    ),
    axes = c(horizontal, list(
      depth = if (!is.null(dims$depth)) .proxy_depth_axis(nc, dims$depth),
      time = if (!is.null(dims$time)) .proxy_time_axis(nc, dims$time, path)
    )),
    crs = .proxy_crs(nc, var)
  )
}

# Depth in metres positive down, in storage order, with edges from CF bounds
# or halfway between centres (cut at the surface).
.proxy_depth_axis <- function(nc, dim) {
  values <- as.vector(.get_dim_values(nc, dim))
  centres <- .read_depth_axis(nc, dim, values)$depth_m
  steps <- diff(centres)
  if (length(steps) && !(all(steps > 0) || all(steps < 0))) {
    rlang::abort(glue::glue(
      "The depth axis `{dim}` must run shallow to deep or deep to shallow."
    ), class = "ocean3d_error_depth_order", call = NULL)
  }
  bounds <- .read_bounds(nc, dim, length(values), strict = TRUE)
  if (!is.null(bounds)) {
    edges <- .depth_bounds_m(nc, dim, values, bounds)
    lower <- pmin(edges[1L, ], edges[2L, ])
    upper <- pmax(edges[1L, ], edges[2L, ])
    source <- "bounds"
  } else if (length(centres) > 1L) {
    edges <- .midpoint_edges(centres)
    lower <- pmax(edges$lower, 0)
    upper <- edges$upper
    source <- "inferred"
  } else {
    lower <- upper <- centres
    source <- "none"
  }
  storage <- if (length(steps) && steps[1L] < 0) "deep_first" else
    "shallow_first"
  list(dim = dim, centres = centres, lower = lower, upper = upper,
       edges = source, storage = storage)
}

# The time axis as a CFtime object, with the file's bounds attached when it
# has them. Time must run strictly forwards (SPEC.md, "Time runs forwards
# only"); the core object doesn't reorder it.
.proxy_time_axis <- function(nc, dim, path) {
  file <- basename(path)
  units <- .read_time_attribute(nc, dim, "units")
  if (!is.null(.read_time_attribute(nc, dim, "climatology")) ||
        (is.character(units) && grepl("since\\s+0+-", units))) {
    rlang::abort(glue::glue(
      "The time axis `{dim}` in {file} is a climatology ",
      "(units '{units}'). Climatological time axes aren't supported yet."
    ), class = "ocean3d_error_time_climatology", call = NULL)
  }
  offsets <- as.vector(.get_dim_values(nc, dim))
  steps <- diff(offsets)
  if (any(steps == 0, na.rm = TRUE)) {
    rlang::abort(glue::glue(
      "The time axis `{dim}` in {file} repeats a time step ",
      "(step {which(steps == 0)[1L] + 1L}). ocean3d needs each time step ",
      "once; remove or combine the repeated steps first."
    ), class = "ocean3d_error_time_repeated", call = NULL)
  }
  if (any(steps < 0, na.rm = TRUE)) {
    order <- if (all(steps < 0)) "runs backwards" else "is out of order"
    rlang::abort(glue::glue(
      "The time axis `{dim}` in {file} {order}. ocean3d needs time steps ",
      "in increasing order. Write a sorted copy of the file first, for ",
      "example with CDO: `cdo sorttimestamp {file} sorted.nc`."
    ), class = "ocean3d_error_time_order", call = NULL)
  }
  bounds <- .read_bounds(nc, dim, length(offsets), strict = TRUE)
  list(
    dim = dim, cftime = .netcdf_time_axis(nc, dim, bounds$values),
    edges = if (is.null(bounds)) "none" else "bounds"
  )
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
  packing <- variable$packing
  packed <- if (is.null(packing$scale_factor) &&
                  is.null(packing$add_offset)) "" else ", packed"
  chunks <- if (is.null(variable$chunksizes)) "" else
    paste0(", chunks ", sizes(variable$chunksizes))
  time_bounds <- if (identical(axes$time$edges, "bounds")) "with" else "no"
  lines <- c(
    class = "SpatProxy (experimental)",
    source = basename(object@source$path),
    variable = glue::glue(
      "{variable$name} ({variable$units}), stored as {packing$storage}",
      "{packed}{chunks}"
    ),
    dimensions = glue::glue(
      "{sizes(variable$dim_lengths)} ({paste(names(dims), collapse = ' x ')})"
    ),
    lon = glue::glue(
      "{extent(axes$lon$centres)}, {length(axes$lon$centres)} cells, ",
      "{sub('_', ' to ', axes$lon$convention)}, {edges(axes$lon)}"
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
