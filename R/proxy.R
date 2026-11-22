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
#' model"). [as_voxel()] reads one time step of it into memory as a
#' [SpatVoxel-class], with the proxy's depths, depth edges, CRS and date.
#'
#' @section Experimental:
#' `SpatProxy` and [as_proxy()] are experimental. Their slots and behaviour
#' may change without notice while the core object is built.
#'
#' @slot source A list with `type` and where the data are: for `"netcdf"`,
#'   `path`, the full path of the file (a path rather than an open connection,
#'   so the object can be saved with [base::saveRDS()]); for `"series"`,
#'   `path`, the files in time order, and `steps`, the number of time steps
#'   in each; for `"voxel"`, `voxel`, the [SpatVoxel-class] in memory.
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
#'   them) and `edges` (`"bounds"` or `"none"`). A one-step climatology, such
#'   as a World Ocean Atlas period, has `climatology = TRUE`, its `units` and
#'   no `cftime`.
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
  # The length of an axis's dimension, or NA when it isn't one of the
  # variable's dimensions.
  size <- function(axis) {
    dim <- axis$dim
    if (is.character(dim) && length(dim) == 1L && dim %in% names(lengths)) {
      lengths[[dim]]
    } else {
      NA_integer_
    }
  }
  gridded <- Filter(Negate(is.null), axes[c("lon", "lat", "depth")])
  for (name in names(gridded)) {
    axis <- gridded[[name]]
    n <- size(axis)
    if (is.na(n) || length(axis$centres) != n || length(axis$lower) != n ||
          length(axis$upper) != n ||
          !isTRUE(all(axis$lower <= axis$upper))) {
      return(glue::glue(
        "the {name} axis needs a dimension of the variable and one centre, ",
        "lower and upper edge per cell, with lower <= upper"
      ))
    }
  }
  depth <- axes$depth$centres
  if (!is.null(depth) && (anyNA(depth) || any(depth < 0) ||
                            anyDuplicated(depth))) {
    return("depths must be distinct metres, positive down")
  }
  if (isTRUE(axes$time$climatology)) {
    if (!identical(as.integer(size(axes$time)), 1L)) {
      return("a climatological time axis must have one step")
    }
  } else if (!is.null(axes$time)) {
    time <- axes$time$cftime
    if (!inherits(time, "CFTime") || !isTRUE(all(diff(time$offsets) > 0)) ||
          !identical(length(time$offsets), as.integer(size(axes$time)))) {
      return(paste("the time axis must be a strictly increasing CFtime",
                   "object with one offset per time step"))
    }
  }
  TRUE
})

#' Describe gridded data without reading its values (experimental)
#'
#' `as_proxy()` describes one variable of gridded data and returns a
#' [SpatProxy-class]. The data can be one netCDF file, a series of netCDF
#' files, or a [SpatVoxel-class] already in memory. For a netCDF file it
#' reads the metadata of the variable and its coordinates and applies the CF
#' metadata once, so later steps can rely on each axis being in a fixed form.
#' No values of the variable are read.
#'
#' Besides one file, two sources give the same kind of proxy:
#'
#' * **A series of files** (several paths): the same variable on the same
#'   grid and depths in every file, each with its own dated time axis, such as
#'   one file per day. The files are put in time order and their time axes
#'   joined into one, in the first file's units. A series stops with an error
#'   when files repeat or overlap each other's time steps, when a file has no
#'   dated time axis, or when the files differ in grid (spatial tiles) or in
#'   the variable's units. The time axis keeps its bounds only when every
#'   file has them.
#' * **A [SpatVoxel-class] already in memory**: one time step (the voxel's
#'   [terra::time()] when set), with the voxel's cells, depths and depth
#'   edges; a voxel without edges gets slabs halfway between its depths, as
#'   `bounds = "midpoint"` gives. The time must be a date or a date-time; a
#'   date is read as midnight UTC. A voxel whose time is a year, a month or a
#'   plain number stops with an error: give it a date-time, or remove its
#'   time with `terra::time(v) <- NULL`. [as_voxel()] returns the voxel
#'   itself.
#'
#' For netCDF files, each axis is put in a fixed form:
#'
#' * **Longitude and latitude:** cell edges come from the coordinate's CF
#'   `bounds` variable, or halfway between neighbouring centres.
#' * **Depth:** converted to metres positive down from the coordinate's
#'   `units` and `positive` attributes. Edges come from the CF bounds, or
#'   halfway between centres, cut at the surface and running half a gap
#'   below the deepest centre. A voxel read with [as_voxel()] keeps these
#'   edges. A depth axis stored deep to shallow is kept in that order and
#'   recorded as such.
#' * **Time:** decoded with CFtime. Bounds are attached when the file has
#'   them; otherwise the axis has none (they aren't inferred). The time axis
#'   must be strictly increasing; a backwards or out-of-order axis is an
#'   error that suggests sorting the file, for example with CDO's
#'   `sorttimestamp`. Only the `standard`, `gregorian` and
#'   `proleptic_gregorian` calendars are supported. A climatological time
#'   axis with one step, such as a World Ocean Atlas period, is recorded as a
#'   climatology without a date; climatologies with more steps aren't
#'   supported yet.
#' * **CRS:** EPSG:4326, unless the variable's `grid_mapping` gives another
#'   geographic CRS as WKT (`crs_wkt` or `spatial_ref`).
#'
#' Regular longitude-latitude grids with coordinate values are supported;
#' curvilinear, projected and rotated-pole grids are rejected. Bad cell bounds
#' (cells that don't contain their centre, have no width or overlap) stop with
#' an error, and so does a variable with no values along a dimension (for
#' example an unlimited time dimension with no records).
#'
#' @inheritSection SpatProxy-class Experimental
#'
#' @param x A netCDF file path, several paths for a series of files, a
#'   connection opened with [ncdf4::nc_open()] (left open), or a
#'   [SpatVoxel-class]. `var` and the dimension arguments apply to files.
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
#'
#' # Read one time step into memory as a SpatVoxel.
#' v <- as_voxel(x, time = 1)
#' depths(v)
#' @export
as_proxy <- function(x, var = NULL, lon_dim = NULL, lat_dim = NULL,
                     depth_dim = NULL, time_dim = NULL) {
  if (methods::is(x, "SpatVoxel")) return(.voxel_to_proxy(x))
  paths <- is.character(x) && length(x) >= 1L && !anyNA(x)
  if (!.is_ncdf4(x) && !paths) {
    rlang::abort(paste(
      "`x` must be a netCDF file path, several paths (a series), an open",
      "ncdf4 connection or a SpatVoxel."
    ), class = "ocean3d_error_proxy_input", call = NULL)
  }
  if (paths && length(x) > 1L) {
    return(.series_to_proxy(x, var, lon_dim, lat_dim, depth_dim, time_dim))
  }
  type <- if (.is_ncdf4(x)) "connection" else "path"
  opened <- .open_source(list(source = x, type = type))
  on.exit(.close_source(opened), add = TRUE)
  nc <- opened$nc
  # The reader reopens the file by its path, so write pending changes to it.
  if (isTRUE(nc$writable)) ncdf4::nc_sync(nc)
  if (is.null(var)) var <- .detect_netcdf_var(nc)
  dims <- .detect_dims(nc, var, lon_dim, lat_dim, depth_dim, time_dim)
  info <- nc$var[[var]]
  empty <- dims$var_dims[as.integer(info$varsize) == 0L]
  if (length(empty)) {
    rlang::abort(glue::glue(
      "Variable `{var}` has no values along {toString(empty)}, so there is ",
      "nothing to read."
    ), class = "ocean3d_error_empty_dimension", call = NULL)
  }

  # Where an axis's edges come from, from the shared readers' `basis`.
  edges_from <- function(basis) {
    switch(basis, explicit_bounds = "bounds",
           singleton_coordinate = "none", "inferred")
  }
  horizontal <- lapply(c(lon = "lon", lat = "lat"), function(axis) {
    # Without a coordinate variable, ncdf4 gives the indices 1 to n, which
    # aren't degrees (depth stops the same way in .read_depth_axis()).
    if (!.has_coordinate_variable(nc, dims[[axis]])) {
      rlang::abort(glue::glue(
        "The {axis} dimension `{dims[[axis]]}` has no coordinate variable, ",
        "so its values aren't known. Only regular longitude-latitude grids ",
        "with coordinate values are supported."
      ), class = "ocean3d_error_coordinate_variable", call = NULL)
    }
    out <- .read_interval_axis(nc, dims[[axis]], axis)
    c(list(dim = dims[[axis]]), out[c("centres", "lower", "upper")],
      list(edges = edges_from(out$basis)))
  })

  # The file's full path, or a URL (e.g. OPeNDAP) as it is. A connection's
  # relative path resolves against the working directory, so stop if no such
  # file is there rather than record a wrong path.
  path <- nc$filename
  if (!.is_url(path)) {
    if (!file.exists(path)) {
      rlang::abort(glue::glue(
        "Can't find the file `{path}` that this connection was opened ",
        "from. Open it with its full path, or pass the path to as_proxy()."
      ), class = "ocean3d_error_proxy_source", call = NULL)
    }
    path <- normalizePath(path)
  }
  # Units are NULL when the file gives none, as for a voxel source; ncdf4's
  # own `units` is "" then.
  units <- .read_attribute(nc, var, "units")
  if (!is.null(units) && !nzchar(trimws(units))) units <- NULL
  methods::new(
    "SpatProxy",
    source = list(type = "netcdf", path = path),
    variable = list(
      name = var, units = units, long_name = info$longname,
      dims = dims$var_dims,
      dim_lengths = stats::setNames(as.integer(info$varsize), dims$var_dims),
      axis_dims = dims[c("lon", "lat", "depth", "time")],
      singleton_dims = dims$singleton_dims,
      # ncdf4 leaves `chunksizes` unset for contiguous storage (storage 1),
      # and NA for a classic-format variable with an unlimited dimension.
      chunksizes = if (identical(as.integer(info$storage), 2L) &&
                         !anyNA(info$chunksizes)) {
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
# only"); the core object doesn't reorder it. A one-step climatology is
# recorded without a CFtime object.
.proxy_time_axis <- function(nc, dim, path) {
  file <- basename(path)
  # A one-step climatology, such as a World Ocean Atlas period, has no date
  # CFtime can decode: it is recorded as a climatology, with its units. A
  # climatology with more steps stops here.
  if (.check_climatology(nc, dim, file, one_step_ok = TRUE)) {
    return(list(dim = dim, cftime = NULL, edges = "none", climatology = TRUE,
                units = .read_attribute(nc, dim, "units")))
  }
  offsets <- as.vector(.get_dim_values(nc, dim))
  units <- .read_attribute(nc, dim, "units")
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
  mapping <- .read_attribute(nc, var, "grid_mapping")
  if (is.null(mapping) || !mapping %in% names(nc$var)) return("EPSG:4326")
  read <- function(name) {
    att <- ncdf4::ncatt_get(nc, mapping, name)
    if (isTRUE(att$hasatt)) att$value
  }
  wkt <- read("crs_wkt")
  if (is.null(wkt)) wkt <- read("spatial_ref")
  name <- read("grid_mapping_name")
  geographic <- if (!is.null(name) && !identical(name, "latitude_longitude")) {
    FALSE # e.g. rotated_latitude_longitude, even with a WKT
  } else if (is.null(wkt)) {
    identical(name, "latitude_longitude")
  } else {
    # A rotated pole is a derived geographic CRS, which terra calls lon/lat.
    !grepl("DERIVINGCONVERSION", wkt, fixed = TRUE) &&
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

# A series of netCDF files as one SpatProxy: the same variable on the same
# grid in every file, each with its own time axis, ordered by time and
# joined into one axis in the first file's units (#94, #96). Each file is
# described and checked by as_proxy() first.
.series_to_proxy <- function(paths, var, lon_dim, lat_dim, depth_dim,
                             time_dim) {
  stop_series <- function(problem) {
    rlang::abort(glue::glue("{problem}."),
                 class = "ocean3d_error_proxy_series", call = NULL)
  }
  full <- normalizePath(paths, mustWork = FALSE)
  if (anyDuplicated(full)) {
    stop_series(glue::glue(
      "{basename(full[duplicated(full)][1])} appears more than once in the ",
      "series"
    ))
  }
  parts <- lapply(paths, as_proxy, var = var, lon_dim = lon_dim,
                  lat_dim = lat_dim, depth_dim = depth_dim,
                  time_dim = time_dim)
  file <- function(p) basename(p@source$path)
  first <- parts[[1L]]
  for (p in parts) {
    if (is.null(p@axes$time) || isTRUE(p@axes$time$climatology)) {
      stop_series(glue::glue(
        "{file(p)} has no dated time axis. Every file in a series needs one, ",
        "so its steps can be placed in time order"
      ))
    }
    if (!identical(p@variable$name, first@variable$name)) {
      stop_series(glue::glue(
        "{file(p)} holds `{p@variable$name}`, not ",
        "`{first@variable$name}`. A series is one variable"
      ))
    }
    if (!identical(p@variable$units, first@variable$units)) {
      units <- vapply(list(p, first), function(q) {
        if (is.null(q@variable$units)) "no units" else q@variable$units
      }, character(1))
      stop_series(glue::glue(
        "{file(p)} gives `{p@variable$name}` in {units[1]}, but ",
        "{file(first)} gives it in {units[2]}. A series needs one unit; ",
        "convert the files first"
      ))
    }
    # The same grid: the same dimensions in the same order and CRS, and the
    # same longitude, latitude and depth cells. Time may differ.
    spatial <- setdiff(p@variable$dims, p@variable$axis_dims$time)
    cells <- lapply(list(p, first), function(q) {
      lapply(q@axes[c("lon", "lat", "depth")], function(axis) {
        axis[c("centres", "lower", "upper")]
      })
    })
    same_grid <- identical(p@variable$dims, first@variable$dims) &&
      identical(p@variable$axis_dims, first@variable$axis_dims) &&
      identical(p@crs, first@crs) &&
      identical(p@variable$dim_lengths[spatial],
                first@variable$dim_lengths[spatial]) &&
      isTRUE(all.equal(cells[[1L]], cells[[2L]]))
    if (!same_grid) {
      stop_series(glue::glue(
        "{file(p)} has a different grid from {file(first)}. A series needs ",
        "the same grid and depths in every file; spatial tiles aren't ",
        "supported"
      ))
    }
  }
  first_date <- vapply(parts, function(p) {
    as.numeric(p@axes$time$cftime$as_timestamp(asPOSIX = TRUE)[1L])
  }, numeric(1))
  parts <- parts[order(first_date)]

  # One time axis, in the first file's units.
  axis <- parts[[1L]]@axes$time$cftime
  for (p in parts[-1L]) {
    # CFtime warns about overlapping steps; the errors below say which file.
    joined <- suppressWarnings(axis + .cftime_in(p@axes$time$cftime, axis))
    steps <- diff(joined$offsets)
    if (any(steps == 0)) {
      # A repeated step, as within one file.
      when <- joined$as_timestamp()[which(steps == 0)[1L] + 1L]
      .abort_repeated_time(
        glue::glue("Time step {when} is also in an earlier file of the series"),
        file(p),
        paste("Each time step must come from one file; leave out the file",
              "or the steps it repeats.")
      )
    }
    if (any(steps < 0)) {
      stop_series(glue::glue(
        "{file(p)} overlaps the time steps of an earlier file. Each file's ",
        "steps must come after the previous file's"
      ))
    }
    axis <- joined
  }
  bounded <- all(vapply(parts, function(p) {
    identical(p@axes$time$edges, "bounds")
  }, logical(1)))

  out <- parts[[1L]]
  steps <- vapply(parts, function(p) {
    p@variable$dim_lengths[[p@axes$time$dim]]
  }, integer(1))
  out@source <- list(type = "series",
                     path = vapply(parts, function(p) p@source$path, ""),
                     steps = steps)
  out@variable$dim_lengths[[out@axes$time$dim]] <- sum(steps)
  out@axes$time <- list(dim = out@axes$time$dim, cftime = axis,
                        edges = if (bounded) "bounds" else "none")
  methods::validObject(out)
  out
}

# A CFtime axis in the units of `target`, with its bounds, so the two can be
# joined with `+`. Unchanged when the units already match.
.cftime_in <- function(cftime, target) {
  if (identical(cftime$cal$definition, target$cal$definition) &&
        identical(cftime$cal$name, target$cal$name)) {
    return(cftime)
  }
  out <- CFtime::CFtime(target$cal$definition, target$cal$name,
                        cftime$as_timestamp())
  bounds <- cftime$get_bounds(format = "%Y-%m-%dT%H:%M:%S")
  if (!is.null(bounds)) {
    # The bounds are converted as a flat vector, not in time order.
    offsets <- suppressWarnings(CFtime::CFtime(
      target$cal$definition, target$cal$name, as.vector(bounds)
    )$offsets)
    out$set_bounds(matrix(offsets, nrow = 2L))
  }
  out
}

# A SpatVoxel already in memory as a SpatProxy, so voxels and files go
# through the same steps. One time step: the voxel's date when it has one.
# Cell edges come from the raster; depth edges from the voxel's stored edges,
# or the slabs its default `bounds` ("midpoint") gives.
.voxel_to_proxy <- function(v) {
  v <- as_voxel(v)
  nc <- terra::ncol(v)
  nr <- terra::nrow(v)
  half <- terra::res(v) / 2
  lon <- terra::xFromCol(v, seq_len(nc))
  lat <- terra::yFromRow(v, seq_len(nr))
  d <- depths(v)
  stored <- depth_bounds(v)
  slabs <- .depth_slabs(d, "midpoint", stored)
  dated <- isTRUE(terra::timeInfo(v)$time)
  # terra::time() is a date-time (step "seconds") or a date ("days"); a
  # year, a year and month, a month of the year or a plain number has no
  # date to place on the time axis, so it stops rather than being guessed
  # or dropped.
  step <- if (dated) terra::timeInfo(v)$step
  if (dated && !step %in% c("seconds", "days")) {
    what <- c(years = "a year", yearmonths = "a year and month",
              months = "a month of the year")[step]
    if (is.na(what)) what <- "a plain number"
    rlang::abort(glue::glue(
      "The voxel's time is {what} ({format(terra::time(v)[1L])}), not a ",
      "date or date-time, so it can't be placed on a time axis. Give it a ",
      "date-time with terra::time(), or remove it with ",
      "`terra::time(v) <- NULL`."
    ), class = "ocean3d_error_voxel_time", call = NULL)
  }
  # A date is midnight UTC.
  seconds <- if (dated) {
    as.numeric(as.POSIXct(terra::time(v)[1L], tz = "UTC"))
  }
  dims <- c("lon", "lat", "depth", if (dated) "time")
  lengths <- c(lon = nc, lat = nr, depth = terra::nlyr(v),
               if (dated) c(time = 1L))
  # The variable's name from terra's own field, which as_voxel() sets and
  # terra keeps through renaming, subsetting and arithmetic; layer names are
  # labels and aren't read. "value", as_voxel()'s default, when the voxel
  # has no single variable name.
  name <- unique(terra::varnames(v))
  if (length(name) != 1L || is.na(name) || !nzchar(name)) name <- "value"
  units <- terra::units(v)[1L]
  crs <- terra::crs(v)
  methods::new(
    "SpatProxy",
    source = list(type = "voxel", voxel = v),
    variable = list(
      name = name, units = if (nzchar(units)) units, long_name = NULL,
      dims = dims, dim_lengths = vapply(lengths, as.integer, integer(1)),
      axis_dims = list(lon = "lon", lat = "lat", depth = "depth",
                       time = if (dated) "time"),
      singleton_dims = character(0), chunksizes = NULL
    ),
    axes = list(
      lon = list(dim = "lon", centres = lon, lower = lon - half[1L],
                 upper = lon + half[1L], edges = "bounds"),
      lat = list(dim = "lat", centres = lat, lower = lat - half[2L],
                 upper = lat + half[2L], edges = "bounds"),
      depth = list(dim = "depth", centres = d, lower = slabs$lower,
                   upper = slabs$upper,
                   edges = if (is.null(stored)) "inferred" else "bounds",
                   storage = "shallow_first"),
      time = if (dated) {
        list(dim = "time", edges = "none", cftime = CFtime::CFtime(
          "seconds since 1970-01-01 00:00:00", "standard", seconds
        ))
      }
    ),
    crs = crs
  )
}

# Reads one index window of a SpatProxy's variable. `window` is a named list
# with up to one element per axis (`lon`, `lat`, `depth`, `time`): the first
# and last index to read, in storage order (so for a deep-first depth axis,
# index 1 is the deepest level). An axis left out is read whole; an axis the
# variable doesn't have counts as one of length 1. Returns the values as a 4D
# array in the fixed order lon, lat, depth, time, with length 1 for an axis
# the variable doesn't have. For a netCDF file, one ncdf4 call reads the
# window, and ncdf4 applies the scale and offset and sets cells equal to
# `missing_value` (or `_FillValue` when there is no `missing_value`) to NA. A
# series reads each file's share of the time window and joins them; a voxel
# in memory gives its own values.
.read_proxy_window <- function(x, window = list()) {
  variable <- x@variable
  axes <- c("lon", "lat", "depth", "time")
  if (!is.list(window) || (length(window) && (
    is.null(names(window)) || !all(names(window) %in% axes) ||
      anyDuplicated(names(window))
  ))) {
    rlang::abort(glue::glue(
      "The read window must be a list naming each axis at most once, from ",
      "lon, lat, depth and time."
    ), class = "ocean3d_error_read_window", call = NULL)
  }
  spans <- list()
  for (axis in axes) {
    dim <- variable$axis_dims[[axis]]
    n <- if (is.null(dim)) 1L else variable$dim_lengths[[dim]]
    span <- if (is.null(window[[axis]])) c(1L, n) else window[[axis]]
    if (!is.numeric(span) || length(span) != 2L || anyNA(span) ||
          any(span != round(span)) || span[1L] < 1L || span[2L] > n ||
          span[1L] > span[2L]) {
      rlang::abort(glue::glue(
        "The read window for `{axis}` must be a first and last index, ",
        "whole numbers from 1 to {n}, with first <= last."
      ), class = "ocean3d_error_read_window", call = NULL)
    }
    spans[[axis]] <- as.integer(span)
  }

  source <- x@source
  if (identical(source$type, "voxel")) {
    v <- source$voxel
    # terra's cells run west to east, then north to south: the storage order
    # the proxy recorded.
    values <- array(terra::values(v),
                    c(terra::ncol(v), terra::nrow(v), terra::nlyr(v), 1L))
    keep <- lapply(spans, function(s) seq(s[1L], s[2L]))
    return(values[keep$lon, keep$lat, keep$depth, keep$time, drop = FALSE])
  }
  if (!identical(source$type, "series")) {
    return(.read_file_window(source$path, variable, spans,
                             variable$dim_lengths))
  }
  # A series: each file's share of the time window, joined along time, the
  # slowest-varying axis of the fixed order.
  time_dim <- variable$axis_dims$time
  last <- cumsum(source$steps)
  first <- last - source$steps + 1L
  parts <- list()
  for (i in which(first <= spans$time[2L] & last >= spans$time[1L])) {
    local <- spans
    local$time <- c(max(spans$time[1L], first[i]),
                    min(spans$time[2L], last[i])) - first[i] + 1L
    lengths <- variable$dim_lengths
    lengths[[time_dim]] <- source$steps[i]
    parts[[length(parts) + 1L]] <- .read_file_window(source$path[i], variable,
                                                     local, lengths)
  }
  shape <- vapply(spans, function(s) s[2L] - s[1L] + 1L, integer(1))
  array(unlist(parts, use.names = FALSE), unname(shape))
}

# One file's index window for .read_proxy_window(): `spans` per axis, and
# `dim_lengths` the variable's dimension sizes in that file. The indices are
# only right for the file the proxy described, so a file whose variable is
# gone or has other dimensions stops.
.read_file_window <- function(path, variable, spans, dim_lengths) {
  axes <- c("lon", "lat", "depth", "time")
  index_list <- list()
  shape <- stats::setNames(rep(1L, 4L), axes)
  for (axis in axes) {
    dim <- variable$axis_dims[[axis]]
    if (is.null(dim)) next
    shape[[axis]] <- spans[[axis]][2L] - spans[[axis]][1L] + 1L
    index_list[[dim]] <- list(start = spans[[axis]][1L],
                              count = shape[[axis]])
  }
  sc <- .build_start_count(variable$dims, index_list, variable$singleton_dims)

  file <- basename(path)
  nc <- if (.is_url(path) || file.exists(path)) ncdf4::nc_open(path)
  if (!is.null(nc)) on.exit(ncdf4::nc_close(nc), add = TRUE)
  info <- if (!is.null(nc)) nc$var[[variable$name]]
  changed <- if (is.null(nc)) {
    glue::glue("Can't find {path}")
  } else if (is.null(info)) {
    glue::glue("{file} no longer has the variable `{variable$name}`")
  } else if (!identical(vapply(info$dim, `[[`, "", "name"), variable$dims) ||
               !identical(as.integer(info$varsize), unname(dim_lengths))) {
    glue::glue(
      "Variable `{variable$name}` in {file} no longer has the dimensions ",
      "this SpatProxy describes"
    )
  }
  if (!is.null(changed)) {
    rlang::abort(glue::glue(
      "{changed}. Build the SpatProxy again with as_proxy()."
    ), class = "ocean3d_error_proxy_changed", call = NULL)
  }
  values <- ncdf4::ncvar_get(nc, variable$name, start = sc$start,
                             count = sc$count, collapse_degen = FALSE)
  # Put the axes in fixed order, then the other (length 1) dimensions, and
  # add length 1 for missing axes, which doesn't move any values. Each step
  # copies the window, so skip the ones that change nothing.
  present <- unlist(variable$axis_dims[axes])
  order <- c(match(present, variable$dims), which(!variable$dims %in% present))
  if (!identical(order, seq_along(variable$dims))) {
    if (!identical(dim(values), sc$count)) dim(values) <- sc$count
    values <- aperm(values, order)
  }
  if (!identical(dim(values), unname(shape))) dim(values) <- unname(shape)
  values
}

# TRUE for a URL such as an OPeNDAP or NCZarr address, which ncdf4 opens
# without a local file.
.is_url <- function(path) grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)

# The voxel builder (SPEC.md, "SpatVoxel / SpatEnvelope", Builder): one time
# step of a SpatProxy read into memory as a SpatVoxel, for as_voxel(). The
# raster's extent comes from the cell edges and its CRS from the proxy; depth
# goes in terra::depth() with the proxy's depth edges, and the date of the
# step in terra::time().
.proxy_to_voxel <- function(x, time, varname) {
  # A voxel already in memory is its own only time step.
  if (identical(x@source$type, "voxel")) {
    .proxy_time_step(x, time)
    return(x@source$voxel)
  }
  axes <- x@axes
  if (is.null(axes$depth)) {
    rlang::abort(glue::glue(
      "Variable `{x@variable$name}` has no depth axis, so it isn't a voxel. ",
      "Read 2D data with terra::rast()."
    ), class = "ocean3d_error_proxy_no_depth", call = NULL)
  }
  step <- .proxy_time_step(x, time)
  window <- if (is.null(step)) list() else list(time = c(step, step))
  values <- .read_proxy_window(x, window)

  # terra's cells run west to east, then north to south.
  lon <- .regular_axis(axes$lon, "longitude")
  lat <- .regular_axis(axes$lat, "latitude")
  lon_order <- order(axes$lon$centres)
  lat_order <- order(axes$lat$centres, decreasing = TRUE)
  n_depth <- length(axes$depth$centres)
  r <- terra::rast(nrows = length(lat_order), ncols = length(lon_order),
                   nlyrs = n_depth, xmin = lon$min, xmax = lon$max,
                   ymin = lat$min, ymax = lat$max, crs = x@crs)
  terra::values(r) <- matrix(values[lon_order, lat_order, , 1L, drop = FALSE],
                             ncol = n_depth)

  depth_bounds <- if (!identical(axes$depth$edges, "none")) {
    cbind(axes$depth$lower, axes$depth$upper)
  }
  v <- as_voxel(r, depths = axes$depth$centres, varname = varname,
                depth_bounds = depth_bounds)
  if (!is.null(x@variable$units)) terra::units(v) <- x@variable$units
  if (!is.null(step) && !isTRUE(axes$time$climatology)) {
    date <- axes$time$cftime$as_timestamp(asPOSIX = TRUE)[step]
    terra::time(v) <- rep(date, terra::nlyr(v))
  }
  v
}

# The index of the time step `time` selects: a step number, or a date or
# date-time matched exactly at its own resolution (a date selects that day).
# NULL when the proxy has no time axis. Without `time`, a proxy with one step
# (or a one-step climatology) gives that step, and one with more is an error.
.proxy_time_step <- function(x, time) {
  axis <- x@axes$time
  if (is.null(axis)) {
    if (!is.null(time)) {
      rlang::abort(glue::glue(
        "Variable `{x@variable$name}` has no time axis, so `time` can't ",
        "select a step."
      ), class = "ocean3d_error_proxy_time", call = NULL)
    }
    return(NULL)
  }
  n <- x@variable$dim_lengths[[axis$dim]]
  if (is.null(time)) {
    if (n == 1L) return(1L)
    rlang::abort(glue::glue(
      "Variable `{x@variable$name}` has {n} time steps. Choose one with ",
      "`time`, a step number or a date."
    ), class = "ocean3d_error_proxy_time", call = NULL)
  }
  if (is.numeric(time) && length(time) == 1L && !is.na(time) &&
        time == round(time)) {
    if (time >= 1 && time <= n) return(as.integer(time))
    rlang::abort(glue::glue("`time` must be a step from 1 to {n}."),
                 class = "ocean3d_error_proxy_time", call = NULL)
  }
  if (isTRUE(axis$climatology)) {
    rlang::abort(
      "The time axis is a climatology, so `time` can only be step 1.",
      class = "ocean3d_error_proxy_time", call = NULL
    )
  }
  dates <- axis$cftime$as_timestamp(asPOSIX = TRUE)
  step <- if (length(time) == 1L) {
    .match_time_index(dates, time, time_match = "exact")
  } else {
    NA_integer_
  }
  if (is.na(step)) {
    shown <- paste(format(range(dates), "%Y-%m-%d %H:%M:%S"),
                   collapse = " to ")
    rlang::abort(glue::glue(
      "No time step matches `time` ({toString(format(time))}). The steps run ",
      "from {shown} UTC; give a step number, or a date or date-time that ",
      "matches one."
    ), class = "ocean3d_error_proxy_time", call = NULL)
  }
  step
}

# One horizontal axis as the extent of a regular grid: `min` and `max` from
# the cell edges. Coordinates stored as float32 are rarely exactly regular:
# rounding builds up along the axis (about 0.02% of a cell on Copernicus
# Marine's 1/12 degree grid), so each centre may sit up to 0.1% of a cell from
# where a regular grid puts it. A truly irregular axis can't be a terra raster
# and stops.
.regular_axis <- function(axis, label) {
  centres <- sort(axis$centres)
  lower <- min(axis$lower)
  upper <- max(axis$upper)
  n <- length(centres)
  if (!(upper > lower)) {
    rlang::abort(glue::glue(
      "The {label} axis has a single cell without bounds, so its width isn't ",
      "known. Add cell bounds to the file."
    ), class = "ocean3d_error_irregular_grid", call = NULL)
  }
  step <- (upper - lower) / n
  expected <- lower + step * (seq_len(n) - 0.5)
  tol <- max(4 * .interval_tolerance(axis$lower, axis$upper), 1e-3 * step)
  if (any(abs(centres - expected) > tol)) {
    rlang::abort(glue::glue(
      "The {label} axis isn't regular: its cells differ in width. A ",
      "SpatVoxel needs a regular grid; regrid the data first."
    ), class = "ocean3d_error_irregular_grid", call = NULL)
  }
  list(min = lower, max = upper)
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
    source = switch(
      object@source$type,
      voxel = "SpatVoxel in memory",
      series = glue::glue(
        "{length(object@source$path)} files, ",
        "{basename(object@source$path[1L])} to ",
        "{basename(object@source$path[length(object@source$path)])}"
      ),
      basename(object@source$path)
    ),
    # Units are NULL when the data give none.
    variable = paste0(
      variable$name,
      if (!is.null(variable$units)) glue::glue(" ({variable$units})"),
      chunks
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
    time = if (is.null(axes$time)) {
      "none"
    } else if (isTRUE(axes$time$climatology)) {
      glue::glue("climatology, 1 step ({axes$time$units})")
    } else {
      n_steps <- length(axes$time$cftime)
      steps <- if (n_steps == 1L) "1 step" else glue::glue("{n_steps} steps")
      glue::glue(
        "{paste(axes$time$cftime$range(), collapse = ' to ')}, {steps}, ",
        "{axes$time$cftime$cal$name} calendar, {time_bounds} bounds"
      )
    },
    crs = if (!nzchar(object@crs)) {
      "none"
    } else if (grepl("^EPSG:", object@crs)) {
      object@crs
    } else {
      terra::crs(object@crs, describe = TRUE)$name
    }
  )
  cat(sprintf("%-10s : %s", names(lines), lines), sep = "\n")
  invisible(object)
})
