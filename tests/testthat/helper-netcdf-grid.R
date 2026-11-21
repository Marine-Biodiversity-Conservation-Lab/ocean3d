# Synthetic gridded netCDF files for tests. One variable, `temp` by default, on
# lon, lat and optional depth and time axes. Unless `values` is given, its
# values are 1, 2, 3, ... in storage order, so reads can be checked against
# ncdf4::ncvar_get() on the whole variable.
#
# `values` are given in lon, lat, depth, time order whatever `dim_order` is; a
# shorter vector is recycled over the grid. `*_bounds` are 2 x n matrices of
# cell edges in the coordinate's units. Set `depth` or `time` to NULL to leave
# that axis out.
# `dim_order` lists the variable's dimensions in ncdf4 order (fastest-varying
# first); axes left out with NULL are skipped. With `prec = "short"`, the
# stored integers are packed with `scale` and `offset`. `extra_vars` adds
# variables on lon and lat, named as in the list and filled with its values.
# `depth_atts` and `time_atts` are attributes for those axes; `depth_name`
# names the depth dimension.
make_grid_file <- function(lon = c(0.5, 1.5, 2.5), lat = c(10.5, 11.5),
                           depth = c(5, 15), time = c(0, 1),
                           depth_units = "m", positive = NULL,
                           time_units = "days since 2020-01-01",
                           calendar = NULL, lon_bounds = NULL,
                           lat_bounds = NULL, depth_bounds = NULL,
                           time_bounds = NULL, time_atts = list(),
                           dim_order = NULL, prec = "double",
                           chunksizes = NULL, scale = 0.1, offset = 0,
                           values = NULL, var = "temp", var_units = "degC",
                           time_unlim = FALSE, extra_vars = list(),
                           depth_atts = list(), depth_name = "depth",
                           path = tempfile(fileext = ".nc")) {
  axes <- list(
    lon = ncdf4::ncdim_def("lon", "degrees_east", lon),
    lat = ncdf4::ncdim_def("lat", "degrees_north", lat)
  )
  if (!is.null(depth)) {
    axes$depth <- ncdf4::ncdim_def(depth_name, depth_units, depth)
  }
  if (!is.null(time)) {
    axes$time <- ncdf4::ncdim_def("time", time_units, time,
                                  unlim = time_unlim)
  }
  shape <- vapply(axes, function(d) length(d$vals), integer(1))
  if (!is.null(values)) {
    stopifnot(length(values) >= 1L, prod(shape) %% length(values) == 0L)
    values <- array(values, shape)
  }
  if (!is.null(dim_order)) {
    axes <- axes[intersect(dim_order, names(axes))]
    if (!is.null(values)) {
      values <- aperm(values, match(names(axes), names(shape)))
    }
  }
  if (is.null(values)) values <- array(seq_len(prod(shape)), shape[names(axes)])

  bounds <- Filter(Negate(is.null), list(
    lon = lon_bounds, lat = lat_bounds, depth = depth_bounds,
    time = time_bounds
  ))
  nv <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
  bounds_vars <- lapply(names(bounds), function(axis) {
    dim <- axes[[axis]]
    ncdf4::ncvar_def(paste0(dim$name, "_bnds"), dim$units, list(nv, dim),
                     prec = "double")
  })
  extra <- lapply(names(extra_vars), function(name) {
    ncdf4::ncvar_def(name, "1", list(axes$lon, axes$lat), prec = "double")
  })
  missval <- if (prec == "short") -32767L else -9999
  chunked <- !is.null(chunksizes)
  main <- ncdf4::ncvar_def(var, var_units, unname(axes), missval,
                           prec = prec,
                           chunksizes = if (chunked) chunksizes else NA,
                           compression = if (chunked) 4 else NA)
  nc <- ncdf4::nc_create(path, c(list(main), bounds_vars, extra),
                         force_v4 = chunked)

  # Raw values go in before the packing attributes, so ncdf4 stores them as
  # they are.
  ncdf4::ncvar_put(nc, var, values)
  if (prec == "short") {
    ncdf4::ncatt_put(nc, var, "scale_factor", scale)
    ncdf4::ncatt_put(nc, var, "add_offset", offset)
  }
  for (axis in names(bounds)) {
    name <- axes[[axis]]$name
    ncdf4::ncvar_put(nc, paste0(name, "_bnds"), bounds[[axis]])
    ncdf4::ncatt_put(nc, name, "bounds", paste0(name, "_bnds"))
  }
  for (name in names(extra_vars)) {
    ncdf4::ncvar_put(nc, name, array(extra_vars[[name]], shape[1:2]))
  }
  if (!is.null(positive)) depth_atts <- c(list(positive = positive), depth_atts)
  for (att in names(depth_atts)) {
    ncdf4::ncatt_put(nc, depth_name, att, depth_atts[[att]])
  }
  if (!is.null(calendar)) ncdf4::ncatt_put(nc, "time", "calendar", calendar)
  for (att in names(time_atts)) {
    ncdf4::ncatt_put(nc, "time", att, time_atts[[att]])
  }
  ncdf4::nc_close(nc)
  path
}

# make_grid_file() with the track and buffer tests' defaults. Unless `values`
# is given, each value says which cell it is in: 10 * lon + 100 * lat + depth +
# 1000 * time, with time in hours since 2020-01-01. `missing = TRUE` empties
# the first longitude. Other arguments go to make_grid_file(), except
# `positive`, `var_units` and `time_units`, which this wrapper sets.
.make_track_nc <- function(
  path, lon = 0:2, lat = c(0, 1), depth = NULL, time = NULL,
  values = outer(outer(10 * lon, 100 * lat, "+"),
                 outer(depth %||% 0, 1000 * (time %||% 0), "+"), "+"),
  missing = FALSE, ...
) {
  if (missing) values[1L, , , ] <- NA_real_
  make_grid_file(
    lon, lat, depth, time, path = path, values = values,
    positive = if (!is.null(depth)) "down", var_units = "degree_Celsius",
    time_units = "hours since 2020-01-01 00:00:00", ...
  )
}

# Sets attributes, given as `name = value`, on variable `var` of an existing
# file. For attributes make_grid_file() can't set, such as those of a bounds
# variable.
set_atts <- function(f, var, ...) {
  nc <- ncdf4::nc_open(f, write = TRUE)
  on.exit(ncdf4::nc_close(nc))
  atts <- list(...)
  for (name in names(atts)) ncdf4::ncatt_put(nc, var, name, atts[[name]])
}

# A SpatProxy of a make_grid_file() file, and the file, which must outlive the
# proxy: each test removes it at the end.
proxy_and_file <- function(...) {
  f <- make_grid_file(...)
  list(proxy = as_proxy(f, var = "temp"), file = f)
}
