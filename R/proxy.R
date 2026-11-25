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
  problem <- if (!is.null(depth)) .depth_problem(depth)
  if (!is.null(problem)) return(paste("depths", problem))
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
#' @param x A netCDF file path or URL (e.g. OPeNDAP or NCZarr), several
#'   paths for a series of files, a connection opened with
#'   [ncdf4::nc_open()] (left open), or a [SpatVoxel-class]. `var` and the dimension arguments apply to files.
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
  file_name <- function(p) basename(p@source$path)
  first <- parts[[1L]]
  for (p in parts) {
    if (is.null(p@axes$time) || isTRUE(p@axes$time$climatology)) {
      stop_series(glue::glue(
        "{file_name(p)} has no dated time axis. Every file in a series ",
        "needs one, so its steps can be placed in time order"
      ))
    }
    if (!identical(p@variable$name, first@variable$name)) {
      stop_series(glue::glue(
        "{file_name(p)} holds `{p@variable$name}`, not ",
        "`{first@variable$name}`. A series is one variable"
      ))
    }
    if (!identical(p@variable$units, first@variable$units)) {
      units <- vapply(list(p, first), function(q) {
        u <- q@variable$units
        if (is.null(u)) "without units" else paste("in", u)
      }, character(1))
      stop_series(glue::glue(
        "{file_name(p)} gives `{p@variable$name}` {units[1]}, but ",
        "{file_name(first)} gives it {units[2]}. A series needs one unit; ",
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
        "{file_name(p)} has a different grid from {file_name(first)}. A ",
        "series needs the same grid and depths in every file; spatial tiles ",
        "aren't supported"
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
    joined <- suppressWarnings(
      axis + .cftime_to_units_of(p@axes$time$cftime, axis)
    )
    steps <- diff(joined$offsets)
    if (any(steps == 0)) {
      # A repeated step, as within one file.
      when <- joined$as_timestamp()[which(steps == 0)[1L] + 1L]
      .abort_repeated_time(
        glue::glue("Time step {when} is also in an earlier file of the series"),
        file_name(p),
        paste("Each time step must come from one file; leave out the file",
              "or the steps it repeats.")
      )
    }
    if (any(steps < 0)) {
      stop_series(glue::glue(
        "{file_name(p)} overlaps the time steps of an earlier file. Each ",
        "file's steps must come after the previous file's"
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
.cftime_to_units_of <- function(cftime, target) {
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
  # The voxel's date-time, or NA without one; a year or a month stops there.
  time <- .voxel_time(v)
  dated <- !is.na(time)
  seconds <- if (dated) as.numeric(time)
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
                   # A single level without stored edges has no slab, as
                   # in a file without depth bounds.
                   edges = if (!is.null(stored)) "bounds" else
                     if (length(d) == 1L) "none" else "inferred",
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
  # A proxy cut down with crop() records the file's longitude index of each
  # cell and its latitude offset in the file, and checks the file against
  # the file's own dimension lengths. Across the grid's longitude edge the
  # cells are two runs in the file, each read as its own window.
  file_spans <- spans
  if (!is.null(source$lat_offset)) {
    file_spans$lat <- spans$lat + source$lat_offset
  }
  lon_file <- seq(spans$lon[1L], spans$lon[2L])
  if (!is.null(source$lon_index)) lon_file <- source$lon_index[lon_file]
  runs <- split(lon_file, cumsum(c(1L, diff(lon_file) != 1L)))
  file_lengths <- if (is.null(source$file_lengths)) variable$dim_lengths else
    source$file_lengths
  shape <- vapply(spans, function(s) s[2L] - s[1L] + 1L, integer(1))
  if (length(runs) > 1L) out <- array(NA_real_, unname(shape))
  done <- 0L
  for (run in runs) {
    file_spans$lon <- range(run)
    values <- if (!identical(source$type, "series")) {
      .read_file_window(source$path, variable, file_spans, file_lengths)
    } else {
      # A series: each file's share of the time window, joined along time,
      # the slowest-varying axis of the fixed order.
      time_dim <- variable$axis_dims$time
      last <- cumsum(source$steps)
      first <- last - source$steps + 1L
      parts <- list()
      for (i in which(first <= spans$time[2L] & last >= spans$time[1L])) {
        local <- file_spans
        local$time <- c(max(spans$time[1L], first[i]),
                        min(spans$time[2L], last[i])) - first[i] + 1L
        lengths <- file_lengths
        lengths[[time_dim]] <- source$steps[i]
        parts[[length(parts) + 1L]] <- .read_file_window(
          source$path[i], variable, local, lengths
        )
      }
      array(unlist(parts, use.names = FALSE),
            c(length(run), unname(shape[-1L])))
    }
    if (length(runs) == 1L) return(values)
    out[done + seq_along(run), , , ] <- values
    done <- done + length(run)
  }
  out
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
  opened <- tryCatch(.open_source(list(source = path, type = "path")),
                     ocean3d_error_netcdf_missing = function(e) NULL)
  on.exit(.close_source(opened), add = TRUE)
  nc <- opened$nc
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

# The step reader (SPEC.md, "SpatVoxel / SpatEnvelope", Builder): the time
# steps `time` selects from a SpatProxy, read into memory. Each raster's
# extent comes from the cell edges and its CRS from the proxy, and the date
# of the step goes in terra::time(). Each kind of variable has one reader:
# with `voxel = TRUE` (as_voxel()) a variable with depth becomes SpatVoxels
# with the proxy's depths and depth edges; with `voxel = FALSE` (rast()) a
# variable without depth becomes plain one-layer SpatRasters named after the
# variable. One raster when `time` names one step, a list for a time window
# or a shape. `depth`, a first and last index in storage order, reads only
# those levels, so volume() can measure a large step in groups of layers;
# NULL reads them all. Every source, a voxel in memory too, reads a span
# through the same per-step path, so each group carries the proxy's depth
# edges.
.proxy_steps <- function(x, time, varname, time_match, max_time_diff,
                         depth = NULL, voxel = TRUE) {
  axes <- x@axes
  if (voxel && is.null(axes$depth)) {
    rlang::abort(glue::glue(
      "Variable `{x@variable$name}` has no depth axis, so it isn't a voxel. ",
      "Read it as a plain SpatRaster with terra::rast() on the proxy."
    ), class = "ocean3d_error_proxy_no_depth", call = NULL)
  }
  if (!voxel && !is.null(axes$depth)) {
    rlang::abort(glue::glue(
      "Variable `{x@variable$name}` has a depth axis. Read it with ",
      "as_voxel(), which keeps the depth levels and their edges; ",
      "terra::rast() on a proxy reads variables without depth."
    ), class = "ocean3d_error_proxy_has_depth", call = NULL)
  }
  .check_max_time_diff(max_time_diff)
  selected <- .proxy_time_step(x, time, time_match, max_time_diff)
  if (identical(x@source$type, "voxel") && is.null(depth)) {
    # A voxel already in memory is its own only time step.
    out <- list(x@source$voxel)
  } else {
    if (!is.null(depth)) {
      axes$depth[c("centres", "lower", "upper")] <- lapply(
        axes$depth[c("centres", "lower", "upper")], `[`,
        seq(depth[1L], depth[2L])
      )
    }
    # terra's cells run west to east, then north to south.
    lon <- .regular_extent(axes$lon, "longitude")
    lat <- .regular_extent(axes$lat, "latitude")
    lon_order <- order(axes$lon$centres)
    lat_order <- order(axes$lat$centres, decreasing = TRUE)
    # A variable without depth reads as one layer.
    n_depth <- max(1L, length(axes$depth$centres))
    depth_bounds <- if (!is.null(axes$depth) &&
                          !identical(axes$depth$edges, "none")) {
      cbind(axes$depth$lower, axes$depth$upper)
    }
    # One read per step; a proxy without a time axis is read whole.
    steps <- if (is.null(selected$steps)) list(NULL) else selected$steps
    out <- lapply(seq_along(steps), function(i) {
      step <- steps[[i]]
      window <- if (is.null(step)) list() else list(time = c(step, step))
      window$depth <- depth
      values <- .read_proxy_window(x, window)
      r <- terra::rast(nrows = length(lat_order), ncols = length(lon_order),
                       nlyrs = n_depth, xmin = lon$min, xmax = lon$max,
                       ymin = lat$min, ymax = lat$max, crs = x@crs)
      terra::values(r) <- matrix(
        values[lon_order, lat_order, , 1L, drop = FALSE], ncol = n_depth
      )
      if (voxel) {
        r <- as_voxel(r, depths = axes$depth$centres, varname = varname,
                      depth_bounds = depth_bounds)
      } else {
        names(r) <- varname
        terra::varnames(r) <- varname
      }
      if (!is.null(x@variable$units)) terra::units(r) <- x@variable$units
      if (!is.null(selected$dates)) {
        terra::time(r) <- rep(selected$dates[i], terra::nlyr(r))
      }
      r
    })
  }
  if (selected$list) out else out[[1L]]
}

# The read budget, in MB: the option `ocean3d.max_read_mb`, 1024 by default,
# which caps how much one read of a SpatProxy holds in memory. Read by
# volume() on a proxy and extract() on a proxy, which split their reads to
# fit it.
.read_budget_mb <- function() {
  budget <- getOption("ocean3d.max_read_mb", 1024)
  if (!is.numeric(budget) || length(budget) != 1L || !is.finite(budget) ||
        budget <= 0) {
    rlang::abort(glue::glue(
      "The option `ocean3d.max_read_mb` must be one positive number of ",
      "megabytes; got {format(budget)}."
    ), class = "ocean3d_error_read_budget", call = NULL)
  }
  budget
}

# `max_time_diff` is NULL or one number of days, 0 or more. Checked before
# any data are read, by the step reader and by extract() on a SpatProxy.
.check_max_time_diff <- function(max_time_diff) {
  if (!is.null(max_time_diff) && (!is.numeric(max_time_diff) ||
                                    length(max_time_diff) != 1L ||
                                    is.na(max_time_diff) ||
                                    max_time_diff < 0)) {
    .abort_voxel_input("`max_time_diff` must be one number of days, 0 or more.")
  }
  invisible(max_time_diff)
}

# The time steps `time` selects: `steps` (their indices, or NULL when the
# proxy has no time axis), `dates` (their date-times, or NULL when the axis
# has none: a climatology or no time axis) and `list` (whether as_voxel()
# returns a list). `time` is a step number; a date or date-time matched with
# the nearest or exact rule (`time_match`, `max_time_diff`, as in
# extract_to_point()); a window of two dates or date-times, giving every step
# inside it; or a SpatVector3D with one feature, whose time (a point) or
# interval (a polygon) is used the same way and always gives a list. Without
# `time`, a proxy with one step (or a one-step climatology) gives that step,
# and one with more is an error.
.proxy_time_step <- function(x, time, time_match = "nearest",
                             max_time_diff = NULL) {
  stop_time <- function(message) {
    rlang::abort(message, class = "ocean3d_error_proxy_time", call = NULL)
  }
  axis <- x@axes$time
  n <- if (is.null(axis)) 1L else x@variable$dim_lengths[[axis$dim]]
  dated <- !is.null(axis) && !isTRUE(axis$climatology)
  selected <- function(steps, as_list) {
    list(steps = steps, dates = if (dated && !is.null(steps)) dates[steps],
         list = as_list)
  }
  dates <- if (dated) axis$cftime$as_timestamp(asPOSIX = TRUE)
  as_list <- FALSE
  resolution <- NULL

  # A shape stands for its own time: a point's time, or a polygon's interval.
  # On a source with one step, such as a WOA climatology, that step is read
  # whatever time the shape has.
  if (methods::is(time, "SpatVector3D")) {
    if (length(time) != 1L) {
      stop_time(glue::glue(
        "`time` must be a SpatVector3D with one feature; this one has ",
        "{length(time)}. Pick one with `x[i]`."
      ))
    }
    # Lines have no single time to read; buffer_3d() turns them into
    # polygons, whose interval is a window.
    .check_geomtype(time, c("polygons", "points"))
    as_list <- TRUE
    if (n == 1L) return(selected(if (!is.null(axis)) 1L, as_list))
    resolution <- time@time_resolution
    table <- as.data.frame(time)
    time <- if (time@type == "polygons") {
      c(table$time_start, table$time_end)
    } else {
      table$time
    }
    if (anyNA(time)) {
      stop_time(glue::glue(
        "The shape has no time, but variable `{x@variable$name}` has {n} ",
        "time steps. Give the shape a time in as_vector3d(), or choose the ",
        "steps with `time`."
      ))
    }
  }

  if (is.null(axis)) {
    if (!is.null(time)) {
      stop_time(glue::glue(
        "Variable `{x@variable$name}` has no time axis, so `time` can't ",
        "select a step."
      ))
    }
    return(selected(NULL, as_list))
  }
  if (is.null(time)) {
    if (n == 1L) return(selected(1L, as_list))
    stop_time(glue::glue(
      "Variable `{x@variable$name}` has {n} time steps. Choose one with ",
      "`time`, a step number or a date, or several with a window of two ",
      "dates."
    ))
  }
  if (anyNA(time)) stop_time("`time` must not be NA.")
  if (is.numeric(time)) {
    if (length(time) == 1L && time == round(time) && time >= 1 && time <= n) {
      return(selected(as.integer(time), as_list))
    }
    stop_time(glue::glue(
      "A numeric `time` must be one whole step number from 1 to {n}. For ",
      "several steps, give a window of two dates or date-times."
    ))
  }
  if (isTRUE(axis$climatology)) {
    stop_time("The time axis is a climatology, so `time` can only be step 1.")
  }
  if (!length(time) %in% 1:2) {
    stop_time(paste(
      "`time` must be a step number, one date or date-time, a window of two,",
      "or a SpatVector3D with one feature."
    ))
  }
  # Dates compare by UTC day, date-times exactly; a shape keeps its own
  # resolution. Text is read as extract_to_point() reads observation times,
  # and warns the same way when it has no time zone.
  prepared <- tryCatch(
    .prepare_time_input(time, "auto", input_name = "`time`"),
    error = function(e) {
      reason <- sub("\n.*", "", conditionMessage(e))
      stop_time(glue::glue(
        "`time` couldn't be read as dates or date-times: {reason} Give a ",
        "date (\"2020-01-02\") or a date-time (\"2020-01-02 06:00\"), and ",
        "for a window two of the same kind."
      ))
    }
  )
  .warn_time_zone(prepared, "`time`")
  if (is.null(resolution) || is.na(resolution)) {
    resolution <- prepared$resolution
  }
  shown <- paste(format(range(dates), "%Y-%m-%d %H:%M:%S"),
                 collapse = " to ")
  if (length(time) == 1L) {
    # Several steps on the day of a date tie; a date-time tells them apart.
    step <- tryCatch(
      .match_time_indices(dates, prepared$values, time_match, max_time_diff,
                          resolution),
      ocean3d_error_repeated_time = function(e) {
        stop_time(glue::glue(
          "Several time steps fall on the day of `time` ",
          "({format(prepared$values)}). Give a date-time to read one of ",
          "them, or a window of two dates to read them all."
        ))
      }
    )
    if (is.na(step)) {
      how <- if (time_match == "exact" || is.null(max_time_diff)) {
        "matches"
      } else {
        glue::glue("is within {max_time_diff} days of")
      }
      stop_time(glue::glue(
        "No time step {how} `time` ({format(prepared$values)}). The steps ",
        "run from {shown} UTC; give a step number, a date or date-time, or ",
        "a window of two."
      ))
    }
    return(selected(step, as_list))
  }
  # A window: every step from the start to the end. With dates, the whole
  # of the end day counts, as a polygon's interval does in as_vector3d().
  edges <- as.numeric(.time_comparison_values(prepared$values, resolution))
  if (edges[1L] > edges[2L]) {
    stop_time("A time window's start must not be after its end.")
  }
  at <- as.numeric(dates)
  inside <- if (resolution == "daily") {
    at >= edges[1L] & at < edges[2L] + 86400
  } else {
    at >= edges[1L] & at <= edges[2L]
  }
  if (!any(inside)) {
    window <- paste(format(prepared$values, "%Y-%m-%d %H:%M:%S"),
                    collapse = " to ")
    stop_time(glue::glue(
      "No time step falls in the window {window} UTC. The steps run from ",
      "{shown} UTC; widen the window (for a shape, `time_window` in ",
      "buffer_3d())."
    ))
  }
  selected(which(inside), TRUE)
}

# The extent of a regular horizontal axis: `min` and `max` from its cell
# edges, once the cells are checked to be regular. Coordinates stored as
# float32 are rarely exactly regular: rounding builds up along the axis (about
# 0.02% of a cell on Copernicus Marine's 1/12 degree grid), so each centre may
# sit up to 0.1% of a cell from where a regular grid puts it. A truly
# irregular axis can't be a terra raster and stops.
.regular_extent <- function(axis, label) {
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

#' Read a variable without depth from a SpatProxy
#'
#' `rast()` reads time steps of a [SpatProxy-class] into memory as plain
#' one-layer [terra::SpatRaster]s. Use it for a variable without depth, such
#' as sea-surface temperature or chlorophyll: each step is one layer named
#' after the variable. A variable with depth stops with an error: read it
#' with [as_voxel()], which keeps its depth levels and their edges. So each
#' kind of variable has one reader.
#'
#' The proxy is read as [as_voxel()] reads it. The time axis is decoded with
#' CFtime rather than by terra, which misreads the time of some files, and
#' the date of each step is set in [terra::time()]. A proxy cut down with
#' `crop()` reads only its cells.
#'
#' [mask()] by [SpatVector3D-class] polygons or points works on the result.
#' It uses only their outlines or positions, not their depths.
#' [summarise_3d()] also takes the one-layer raster and weights each cell by
#' its area.
#'
#' @inheritSection SpatProxy-class Experimental
#'
#' @param x A [SpatProxy-class].
#' @param time,time_match,max_time_diff Which time steps to read, as in
#'   [as_voxel()]: a step number, a date or date-time (matched with
#'   `time_match`, by default to the nearest step), a window of two dates or
#'   date-times, or a [SpatVector3D-class] with one feature.
#' @param ... Not used.
#'
#' @returns A one-layer [terra::SpatRaster], or a list of them, one per
#'   step, for a time window or a shape.
#' @seealso [as_voxel()] for a variable with depth.
#' @examples
#' # A 2D temperature field on 1, 3 and 6 January 2020.
#' p <- as_proxy(system.file("extdata", "example_2d.nc", package = "ocean3d"))
#' r <- terra::rast(p, time = "2020-01-03")
#' r
#' terra::time(r)
#'
#' # A window gives one raster per step inside it.
#' length(terra::rast(p, time = c("2020-01-01", "2020-01-03")))
#'
#' # Summarised with each cell weighted by its area.
#' summarise_3d(r, statistic = c("mean", "sd"))
#' @name rast-SpatProxy
#' @aliases rast,SpatProxy-method
#' @export
setMethod("rast", "SpatProxy", function(x, time = NULL,
                                        time_match = c("nearest", "exact"),
                                        max_time_diff = NULL, ...) {
  if (...length()) {
    .abort_voxel_input(paste(
      "rast() on a SpatProxy takes only `time`, `time_match` and",
      "`max_time_diff`."
    ))
  }
  time_match <- match.arg(time_match)
  .proxy_steps(x, time, x@variable$name, time_match, max_time_diff,
               voxel = FALSE)
})

#' @rdname SpatProxy-class
#' @param x A `SpatProxy`.
#' @param y For `crop()`: the area to keep, as a [terra::SpatExtent] or
#'   anything with an extent ([terra::SpatRaster], [terra::SpatVector], an
#'   sf object, or a [SpatVector3D-class] such as a buffer from
#'   [buffer_3d()]). The result keeps the cells whose longitude and latitude
#'   ranges overlap the area's extent, still without reading values, so
#'   [as_voxel()] then reads only those cells from the file. A polygon is cut
#'   to its extent; to keep only the cells inside it, mask the voxel with it
#'   afterwards. An area in another CRS is projected first, and on a
#'   longitude-latitude grid an extent in the other longitude convention (0 to
#'   360 or -180 to 180) finds its cells too. An area that crosses the edge of
#'   the grid's longitudes (the antimeridian on a -180 to 180 grid, or 0
#'   degrees on a 0 to 360 grid) is read from both ends of the grid: the
#'   cropped grid runs on across the edge, for example from 179 to 181
#'   degrees, so [as_voxel()] gives one voxel and [mask()] with the same
#'   shape keeps the right cells. A [SpatVector3D-class] across the
#'   antimeridian has such an extent (see its `ext()`). A `SpatVector3D` with
#'   no features has no extent, so `crop()` stops.
#' @param ... Not used.
#' @export
setMethod("crop", "SpatProxy", function(x, y, ...) {
  if (inherits(y, c("sf", "sfc"))) y <- terra::vect(y)
  # The area's extent in the proxy's CRS, for every source.
  if (!inherits(y, "SpatExtent")) {
    crs <- terra::crs(y)
    if (nzchar(crs) && nzchar(x@crs) && !terra::same.crs(crs, x@crs)) {
      # Points along the box's edges, so the projected extent also reaches
      # where a curved projection (a polar one, say) bends an edge beyond
      # its corners.
      e <- as.vector(terra::ext(y))
      box <- terra::densify(
        terra::as.polygons(terra::ext(y), crs = crs),
        interval = max(e[["xmax"]] - e[["xmin"]], e[["ymax"]] - e[["ymin"]]) /
          100,
        flat = TRUE
      )
      projected <- terra::project(box, x@crs)
      y <- terra::ext(projected)
      # A box across the antimeridian: going round its edge, the projected
      # longitude jumps between 180 and -180 (consecutive points more than
      # 180 degrees apart), so their range spans nearly the globe. Written
      # from 0 to 360 they are continuous, so take that range; crop() then
      # finds the cells, or says the area crosses the grid's edge. A box
      # that covers every longitude has no such jump and keeps its range.
      if (terra::is.lonlat(x@crs)) {
        lon <- terra::crds(projected)[, 1L]
        if (any(abs(diff(c(lon, lon[1L]))) > 180)) {
          east <- ifelse(lon < 0, lon + 360, lon)
          y <- terra::ext(min(east), max(east), terra::ymin(y),
                          terra::ymax(y))
        }
      }
      # A box around a pole: no point on its edges reaches the pole, so the
      # projected extent stops short of it. If a pole lies inside the box,
      # the area reaches that pole at every longitude. A projection that
      # can't place a pole (Mercator, say) can't contain one either.
      if (terra::is.lonlat(x@crs) && !terra::is.lonlat(crs)) {
        for (lat in c(90, -90)) {
          pole <- tryCatch(
            terra::crds(terra::project(terra::vect(cbind(0, lat),
                                                   crs = x@crs), crs)),
            error = function(err) NULL
          )
          inside <- !is.null(pole) && all(is.finite(pole)) &&
            pole[1L] >= e[["xmin"]] && pole[1L] <= e[["xmax"]] &&
            pole[2L] >= e[["ymin"]] && pole[2L] <= e[["ymax"]]
          if (inside) {
            box_lat <- as.vector(y)[c("ymin", "ymax")]
            y <- terra::ext(-180, 180, min(box_lat, lat), max(box_lat, lat))
          }
        }
      }
    }
    if (!inherits(y, "SpatExtent")) y <- terra::ext(y)
  }
  e <- as.vector(y)
  # Longitudes wrap only on a longitude-latitude grid, not in map units.
  lonlat <- nzchar(x@crs) && terra::is.lonlat(x@crs)
  spans <- list(
    lon = .crop_span(x@axes$lon, e[["xmin"]], e[["xmax"]], "longitude",
                     wrap = lonlat),
    lat = .crop_span(x@axes$lat, e[["ymin"]], e[["ymax"]], "latitude")
  )
  # A voxel in memory is cropped as the raster it is, to the edges of the
  # same cells a file would give. Across the grid's longitude edge, each run
  # is cropped on its own, shifted by its 360 degrees and the two joined.
  if (identical(x@source$type, "voxel")) {
    v <- x@source$voxel
    lat <- spans$lat$index
    lat_edges <- range(x@axes$lat$lower[lat], x@axes$lat$upper[lat])
    parts <- lapply(split(spans$lon$index, spans$lon$shift), function(lon) {
      box <- terra::ext(c(range(x@axes$lon$lower[lon], x@axes$lon$upper[lon]),
                          lat_edges))
      shift <- spans$lon$shift[match(lon[1L], spans$lon$index)]
      # remove the 3D class
      terra::shift(terra::crop(methods::as(v, "SpatRaster"), box), dx = shift)
    })
    if (length(parts) == 1L) {
      return(.voxel_to_proxy(terra::crop(v, terra::ext(parts[[1L]]))))
    }
    joined <- terra::merge(parts[[1L]], parts[[2L]])
    terra::time(joined) <- terra::time(v)
    terra::units(joined) <- terra::units(v)
    terra::varnames(joined) <- terra::varnames(v)
    names(joined) <- names(v)
    return(.voxel_to_proxy(as_voxel(joined, depths = depths(v),
                                    varname = terra::varnames(v)[1L],
                                    depth_bounds = depth_bounds(v))))
  }
  # The cropped axes and lengths describe the smaller grid. The file's own
  # lengths, its longitude index of each cell and its latitude offset let the
  # reader find it in the file.
  out <- x
  source <- x@source
  source$file_lengths <- source$file_lengths %||% x@variable$dim_lengths
  lon_index <- source$lon_index %||% seq_along(x@axes$lon$centres)
  source$lon_index <- lon_index[spans$lon$index]
  source$lat_offset <- (source$lat_offset %||% 0L) + spans$lat$index[1L] - 1L
  for (axis in names(spans)) {
    keep <- spans[[axis]]$index
    out@axes[[axis]][c("centres", "lower", "upper")] <-
      lapply(x@axes[[axis]][c("centres", "lower", "upper")],
             function(values) values[keep] + spans[[axis]]$shift)
    out@variable$dim_lengths[[x@variable$axis_dims[[axis]]]] <- length(keep)
  }
  out@source <- source
  methods::validObject(out)
  out
})

# Internal: the cells of a longitude or latitude `axis` whose range overlaps
# `low` to `high`, as `index` (storage indices, in the cropped axis's order)
# and `shift` (degrees to add to each cell's coordinates). With `wrap`, the
# range is also tried 360 degrees either way, so an area in the other
# longitude convention finds its cells. An area whose edge only touches a
# cell, within rounding (as after projecting it), doesn't take that cell. An
# area with no width (a single point) takes the cell that contains it, or
# both cells when it sits on their shared edge. The cells are next to each
# other in the file, or, for an area across the grid's longitude edge (the
# antimeridian of a -180 to 180 grid), one run at each end: the run with the
# higher longitudes comes first, and the other is shifted by 360 degrees so
# the axis stays continuous, in the area's own convention. Its error is
# raised from here only, for both axes of crop().
.crop_span <- function(axis, low, high, label, wrap = FALSE) {
  lo <- pmin(axis$lower, axis$upper)
  hi <- pmax(axis$lower, axis$upper)
  tol <- .interval_tolerance(lo, hi)
  shifts <- if (wrap) c(0, 360, -360) else 0
  hit <- sort(unique(unlist(lapply(shifts, function(s) {
    if (high - low <= tol) {
      which(lo - tol <= low + s & hi + tol >= high + s)
    } else {
      which(hi - tol > low + s & lo + tol < high + s)
    }
  }))))
  if (!length(hit)) {
    rlang::abort(glue::glue("The area doesn't overlap the grid's {label}s."),
                 class = "ocean3d_error_proxy_crop", call = NULL)
  }
  runs <- split(hit, cumsum(c(1L, diff(hit) != 1L)))
  if (length(runs) == 1L) return(list(index = hit, shift = rep(0, length(hit))))
  # Axes run strictly one way, so cells split into two runs only when a
  # wrapped area takes both ends of the grid, across its longitude edge.
  east <- which.max(vapply(runs, function(r) mean(axis$centres[r]), 0))
  east_run <- runs[[east]]
  west_run <- runs[[3L - east]]
  # An area reaching west of the grid keeps its western longitudes.
  west_side <- low < min(lo) - tol
  list(index = c(east_run, west_run),
       shift = c(rep(if (west_side) -360 else 0, length(east_run)),
                 rep(if (west_side) 0 else 360, length(west_run))))
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
