# Synthetic gridded netCDF files for the core object (as_proxy()) and the
# shared axis helpers. One variable, `temp`, on lon, lat and optional depth and
# time axes. Its values are 1, 2, 3, ... in storage order, so reads can be
# checked against ncdf4::ncvar_get() on the whole variable.
#
# `*_bounds` are 2 x n matrices of cell edges in the coordinate's units. Set
# `depth` or `time` to NULL to leave that axis out. `dim_order` lists the
# variable's dimensions in ncdf4 order (fastest-varying first). With `prec =
# "short"`, the stored integers are packed with `scale` and `offset`.
make_grid_file <- function(lon = c(0.5, 1.5, 2.5), lat = c(10.5, 11.5),
                           depth = c(5, 15), time = c(0, 1),
                           depth_units = "m", positive = NULL,
                           time_units = "days since 2020-01-01",
                           calendar = NULL, lon_bounds = NULL,
                           lat_bounds = NULL, depth_bounds = NULL,
                           time_bounds = NULL, time_atts = list(),
                           dim_order = NULL, prec = "double",
                           chunksizes = NULL, scale = 0.1, offset = 0) {
  f <- tempfile(fileext = ".nc")
  axes <- list(
    lon = ncdf4::ncdim_def("lon", "degrees_east", lon),
    lat = ncdf4::ncdim_def("lat", "degrees_north", lat)
  )
  if (!is.null(depth)) {
    axes$depth <- ncdf4::ncdim_def("depth", depth_units, depth)
  }
  if (!is.null(time)) {
    axes$time <- ncdf4::ncdim_def("time", time_units, time)
  }
  if (!is.null(dim_order)) axes <- axes[dim_order]

  bounds <- Filter(Negate(is.null), list(
    lon = lon_bounds, lat = lat_bounds, depth = depth_bounds,
    time = time_bounds
  ))
  nv <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
  bounds_vars <- lapply(names(bounds), function(axis) {
    ncdf4::ncvar_def(paste0(axis, "_bnds"), axes[[axis]]$units,
                     list(nv, axes[[axis]]), prec = "double")
  })
  missval <- if (prec == "short") -32767L else -9999
  chunked <- !is.null(chunksizes)
  var <- ncdf4::ncvar_def("temp", "degC", unname(axes), missval,
                          prec = prec,
                          chunksizes = if (chunked) chunksizes else NA,
                          compression = if (chunked) 4 else NA)
  nc <- ncdf4::nc_create(f, c(list(var), bounds_vars), force_v4 = chunked)

  # Raw values go in before the packing attributes, so ncdf4 stores them as
  # they are.
  shape <- vapply(axes, function(d) length(d$vals), integer(1))
  ncdf4::ncvar_put(nc, "temp", array(seq_len(prod(shape)), shape))
  if (prec == "short") {
    ncdf4::ncatt_put(nc, "temp", "scale_factor", scale)
    ncdf4::ncatt_put(nc, "temp", "add_offset", offset)
  }
  for (axis in names(bounds)) {
    ncdf4::ncvar_put(nc, paste0(axis, "_bnds"), bounds[[axis]])
    ncdf4::ncatt_put(nc, axis, "bounds", paste0(axis, "_bnds"))
  }
  if (!is.null(positive)) ncdf4::ncatt_put(nc, "depth", "positive", positive)
  if (!is.null(calendar)) ncdf4::ncatt_put(nc, "time", "calendar", calendar)
  for (att in names(time_atts)) {
    ncdf4::ncatt_put(nc, "time", att, time_atts[[att]])
  }
  ncdf4::nc_close(nc)
  f
}
