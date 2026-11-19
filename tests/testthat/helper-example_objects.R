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
# Helpers for 3D-object tests: a lon/lat template grid and a polygon inside
# the grid.
make_template <- function(vals = NA) {
  terra::rast(nrows = 5, ncols = 5,
              xmin = 0, xmax = 5, ymin = 0, ymax = 5,
              vals = vals, crs = "EPSG:4326")
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
# Helpers for intersect and volume tests. The LAEA projection makes every cell
# 1 km x 1 km, so volumes are hand-computable. With `vals`, make_grid() is a
# single-layer footprint (NA = outside).
make_grid <- function(ncol = 3, nrow = 3, vals = NULL) {
  r <- terra::rast(nrows = nrow, ncols = ncol,
                   xmin = 0, xmax = ncol * 1000,
                   ymin = 0, ymax = nrow * 1000,
                   crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
  if (is.null(vals)) r else terra::setValues(r, vals)
}

# Helper: a SpatEnvelope with per-cell depth limits (NA marks absence).
make_range_rast <- function(depth_min_vals, depth_max_vals,
                            ncol = 3, nrow = 3) {
  dmin <- make_grid(ncol, nrow)
  terra::values(dmin) <- depth_min_vals
  names(dmin) <- "depth_min"

  dmax <- terra::rast(dmin)
  terra::values(dmax) <- depth_max_vals
  names(dmax) <- "depth_max"

  as_envelope(c(dmin, dmax))
}

# Helper: a voxel's stored depth edges for its layers, as a two-column (lower,
# upper) matrix, or NULL when it carries none for every layer.
edges_of <- function(v) {
  b <- depth_bounds(v)
  if (is.null(b)) NULL else unname(as.matrix(b[c("lower", "upper")]))
}

# Helper: a SpatVoxel from a cells x depths matrix of values (NA = unoccupied).
make_voxel <- function(vals, depths, ncol = 3, nrow = 3, varname = "presence") {
  r <- terra::rast(make_grid(ncol, nrow), nlyrs = length(depths))
  terra::values(r) <- vals
  as_voxel(r, depths = depths, varname = varname)
}

# Helper: a temperature field on five levels over make_grid()'s 3x3 grid,
# every cell a distinct value so a misplaced mask would show.
make_field <- function(depths = c(0, 50, 100, 150, 200)) {
  vals <- matrix(seq_len(9 * length(depths)), nrow = 9)
  make_voxel(vals, depths, varname = "temp")
}

# Helper: a polygon covering the centres of the left two columns of
# make_grid()'s 3x3 grid, in the grid's CRS.
make_left_polygon <- function() {
  terra::vect("POLYGON ((0 0, 2000 0, 2000 3000, 0 3000, 0 0))",
              crs = terra::crs(make_grid()))
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


# Helpers for track extraction tests
# ----------------------------------

.track_clock <- function(hours = c(0, 1, 2)) {
  as.POSIXct("2020-01-01", tz = "UTC") +
    hours * 3600
}


.track_test_polygon <- function(xmin, xmax, ymin, ymax) {
  sf::st_sfc(
    sf::st_polygon(
      list(
        rbind(
          c(xmin, ymin), c(xmax, ymin),
          c(xmax, ymax), c(xmin, ymax), c(xmin, ymin)
        )
      )
    ),
    crs = 4326
  )
}
