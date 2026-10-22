# Helper: 2x2 voxel over four standard depths, cell 1 has an interior gap
make_multidepth_rast <- function(depths = c(0, 100, 200, 300),
                                 vals = list(c(10, NA,  5, NA),
                                             c(NA, 20,  6, NA),
                                             c(12, 21,  7, NA),
                                             c(NA, NA,  8, NA)),
                                 varname = "temp") {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  lays <- lapply(vals, function(v) terra::setValues(r, v))
  out <- terra::rast(lays)
  names(out) <- paste0(varname, "_depth=", depths)
  out
}
# Helpers for 3D-object tests: a lon/lat template grid, a GEBCO-style
# elevation raster (negative below sea level), and a polygon inside the grid.
make_template <- function(vals = NA) {
  terra::rast(nrows = 5, ncols = 5,
              xmin = 0, xmax = 5, ymin = 0, ymax = 5,
              vals = vals, crs = "EPSG:4326")
}
make_elevation <- function(elev = -500) {
  r <- make_template()
  terra::values(r) <- elev
  r
}
make_polygon <- function() {
  poly <- sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(1, 1), c(4, 1), c(4, 4), c(1, 4), c(1, 1)
    ))),
    crs = "EPSG:4326"
  )
  sf::st_sf(geometry = poly)
}
# Helper: 2x2 footprint, cells 1-3 present, cell 4 absent. Values are
# deliberately varied (and one is 0) because only their non-NA-ness counts.
make_footprint <- function(vals = c(1, 5, 0, NA)) {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::setValues(r, vals)
}
# Helper: an occupancy stack of the shape a vertical profile is handed — one
# logical layer per depth, TRUE where a cell's envelope overlaps the slab that
# depth stands for.
# Built by hand rather than through envelope_to_voxel() so profile tests do not
# depend on how occupancy is derived. Defaults give a 2x2 grid over four
# depths whose cells occupy 1, 2, 4 and 0 levels respectively; the last is the
# cell whose envelope contains none of the depths.
make_occupancy <- function(cells = list(c(TRUE, FALSE, FALSE, FALSE),
                                        c(TRUE, TRUE, FALSE, FALSE),
                                        c(TRUE, TRUE, TRUE, TRUE),
                                        c(FALSE, FALSE, FALSE, FALSE)),
                           depths = c(0, 100, 200, 300)) {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  # A layer holds one value per cell, so the per-cell rows are read down.
  lays <- lapply(seq_along(depths), function(i) {
    terra::setValues(r, vapply(cells, function(c) c[i], logical(1)))
  })
  terra::rast(lays)
}


# Synthetic sources for track extraction
# ---------------------------------------

.make_track_nc <- function(
  path, lon = 0:2, lat = c(0, 1),
  depth = NULL, time = NULL, extra_static = FALSE, missing = FALSE
) {
  x <- ncdf4::ncdim_def("lon", "degrees_east", lon)
  y <- ncdf4::ncdim_def("lat", "degrees_north", lat)
  z <- if (is.null(depth)) {
    NULL
  } else {
    ncdf4::ncdim_def("depth", "m", depth)
  }
  t <- if (is.null(time)) {
    NULL
  } else {
    ncdf4::ncdim_def("time", "hours since 2020-01-01 00:00:00", time)
  }
  dims <- Filter(
    Negate(is.null),
    list(x, y, z, t)
  )
  variables <- list(
    ncdf4::ncvar_def(
      "temp", "degree_Celsius", dims,
      missval = -9999,
      prec = "double"
    )
  )
  if (extra_static) {
    variables <- c(
      variables, list(
        ncdf4::ncvar_def(
          "static", "1", list(x, y),
          prec = "double"
        )
      )
    )
  }
  nc <- ncdf4::nc_create(path, variables)
  on.exit(ncdf4::nc_close(nc))
  zz <- if (is.null(depth)) {
    0
  } else {
    depth
  }
  tt <- if (is.null(time)) {
    0
  } else {
    time
  }
  field <- array(
    NA_real_, c(
      length(lon),
      length(lat),
      length(zz),
      length(tt)
    )
  )
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      for (k in seq_along(zz)) {
        for (l in seq_along(tt)) {
          field[i, j, k, l] <- 10 * lon[i] + 100 * lat[j] + zz[k] + 1000 *
            tt[l]
        }
      }
    }
  }
  if (missing) {
    field[1L, , , ] <- NA_real_
  }
  ncdf4::ncvar_put(
    nc, "temp", array(
      field, vapply(
        dims, function(d) length(d$vals),
        integer(1)
      )
    )
  )
  if (!is.null(z)) {
    ncdf4::ncatt_put(nc, "depth", "positive", "down")
  }
  if (extra_static) {
    ncdf4::ncvar_put(
      nc, "static", matrix(
        7, length(lon),
        length(lat)
      )
    )
  }
  invisible(path)
}


.track_clock <- function(hours = c(0, 1, 2)) {
  as.POSIXct("2020-01-01", tz = "UTC") +
    hours * 3600
}
