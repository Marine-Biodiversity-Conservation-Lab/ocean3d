# CF test cases for the data model trial (#112).
#
# Each case writes one small netCDF file and states the expected result:
#   - for #88-#96, the expected output written in the issue;
#   - otherwise, the package's current, tested behaviour.
# Values are chosen so the value at a cell identifies the cell (usually the
# index along the axis being tested).
#
# A case is a list with:
#   id, group, from, what  - labels for the results table
#   build(path)            - writes the file
#   var                    - variable to read
#   query                  - lon, lat and optionally depth and time
#   expected               - value at the query cell (NA_character_ = error
#                            expected); `expected_note` explains it

library(ncdf4)

# 2 x 2 grid with a time axis; temp at time step t equals t (#88-#92 helper).
.mk_time <- function(path, time_vals, units, calendar = NULL, unlim = FALSE,
                     bounds = NULL) {
  lon <- ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdim_def("lat", "degrees_north", c(40, 41))
  tm <- ncdim_def("time", units, time_vals, unlim = unlim)
  vars <- list(ncvar_def("temp", "1", list(lon, lat, tm), -9999,
                         prec = "double"))
  if (!is.null(bounds)) {
    nv <- ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
    vars[[2]] <- ncvar_def("time_bnds", "", list(nv, tm), prec = "double")
  }
  nc <- nc_create(path, vars)
  nt <- length(time_vals)
  ncvar_put(nc, "temp", array(rep(seq_len(nt), each = 4), c(2, 2, nt)),
            count = c(2, 2, nt))
  if (!is.null(calendar)) ncatt_put(nc, "time", "calendar", calendar)
  if (!is.null(bounds)) {
    ncvar_put(nc, "time_bnds", bounds)
    ncatt_put(nc, "time", "bounds", "time_bnds")
  }
  nc_close(nc)
}

# 2 x 2 grid with a vertical axis; value at layer k equals k.
.mk_depth <- function(path, name, values, units, positive = NULL,
                      standard_name = NULL, axis = NULL, bounds = NULL) {
  lon <- ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdim_def("lat", "degrees_north", c(40, 41))
  z <- ncdim_def(name, units, values)
  vars <- list(ncvar_def("temp", "1", list(lon, lat, z), -9999,
                         prec = "double"))
  if (!is.null(bounds)) {
    nv <- ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
    vars[[2]] <- ncvar_def(paste0(name, "_bnds"), units, list(nv, z),
                           prec = "double")
  }
  nc <- nc_create(path, vars)
  nz <- length(values)
  ncvar_put(nc, "temp", array(rep(seq_len(nz), each = 4), c(2, 2, nz)))
  if (!is.null(positive)) ncatt_put(nc, name, "positive", positive)
  if (!is.null(standard_name)) {
    ncatt_put(nc, name, "standard_name", standard_name)
  }
  if (!is.null(axis)) ncatt_put(nc, name, "axis", axis)
  if (!is.null(bounds)) {
    ncvar_put(nc, paste0(name, "_bnds"), bounds)
    ncatt_put(nc, name, "bounds", paste0(name, "_bnds"))
  }
  nc_close(nc)
}

# Longitude axis with one latitude row; value at cell i equals vals[i].
.mk_lon <- function(path, lons, vals = seq_along(lons), bounds = NULL,
                    lats = 40) {
  lon <- ncdim_def("lon", "degrees_east", lons)
  lat <- ncdim_def("lat", "degrees_north", lats)
  vars <- list(ncvar_def("temp", "1", list(lon, lat), -9999, prec = "double"))
  if (!is.null(bounds)) {
    nv <- ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
    vars[[2]] <- ncvar_def("lon_bnds", "degrees_east", list(nv, lon),
                           prec = "double")
  }
  nc <- nc_create(path, vars)
  ncvar_put(nc, "temp", matrix(rep(vals, length(lats)), length(lons)))
  if (!is.null(bounds)) {
    ncvar_put(nc, "lon_bnds", bounds)
    ncatt_put(nc, "lon", "bounds", "lon_bnds")
  }
  nc_close(nc)
}

.utc <- function(x) as.POSIXct(x, tz = "UTC")

cases <- list(
  # Time ----------------------------------------------------------------------
  list(
    id = "T01", group = "time", from = "test fixture",
    what = "days since 2020-01-01, steps 0, 2, 5",
    build = function(p) .mk_time(p, c(0, 2, 5), "days since 2020-01-01"),
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 2,
    expected_note = "2020-01-01, 2020-01-03, 2020-01-06"
  ),
  list(
    id = "T02", group = "time", from = "test (origin hour)",
    what = "hours since 2020-01-01 06:00:00, steps 0, 24, 48",
    build = function(p) {
      .mk_time(p, c(0, 24, 48), "hours since 2020-01-01 06:00:00")
    },
    query = list(lon = 0, lat = 40, time = .utc("2020-01-02 06:00:00")),
    expected = 2, expected_note = "steps at 06:00 UTC"
  ),
  list(
    id = "T03", group = "time", from = "test (calendars)",
    what = "noleap calendar, days since 2000-01-01, steps 58, 59",
    build = function(p) {
      .mk_time(p, c(58, 59), "days since 2000-01-01", calendar = "noleap")
    },
    query = list(lon = 0, lat = 40, time = .utc("2000-03-01")),
    expected = 2,
    expected_note = "2000-02-28, 2000-03-01 (no 29 February)"
  ),
  list(
    id = "T04", group = "time", from = "test (calendars)",
    what = "360_day calendar, days since 2000-01-01, steps 58, 59",
    build = function(p) {
      .mk_time(p, c(58, 59), "days since 2000-01-01", calendar = "360_day")
    },
    query = list(lon = 0, lat = 40, time = "2000-02-30"),
    expected = 2,
    expected_note = "2000-02-29, 2000-02-30 (not valid R dates)"
  ),
  list(
    id = "T05", group = "time", from = "#88",
    what = "months since 1955-01-01 00:00:00, one step (WOA23)",
    build = function(p) .mk_time(p, 0, "months since 1955-01-01 00:00:00"),
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 1, expected_note = "the single step"
  ),
  list(
    id = "T06", group = "time", from = "#88",
    what = "days since 2020-01-01 UTC (GHRSST)",
    build = function(p) .mk_time(p, c(0, 2, 5), "days since 2020-01-01 UTC"),
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 2, expected_note = ""
  ),
  list(
    id = "T07", group = "time", from = "#88",
    what = "hours since 2020-01-01 00:00:0.0 (NCEP, legacy ERA5)",
    build = function(p) {
      .mk_time(p, c(0, 48, 120), "hours since 2020-01-01 00:00:0.0")
    },
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 2, expected_note = ""
  ),
  list(
    id = "T08", group = "time", from = "#88",
    what = "days since 2020-1-1 00:00:00 (PSL ERSST)",
    build = function(p) .mk_time(p, c(0, 2, 5), "days since 2020-1-1 00:00:00"),
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 2, expected_note = ""
  ),
  list(
    id = "T09", group = "time", from = "#88",
    what = "hours since 2020-01-01 00:00, no seconds (NCEI ERSST)",
    build = function(p) {
      .mk_time(p, c(0, 48, 120), "hours since 2020-01-01 00:00")
    },
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 2, expected_note = ""
  ),
  list(
    id = "T10", group = "time", from = "#88",
    what = "days since 0001-01-01, proleptic_gregorian",
    build = function(p) {
      .mk_time(p, c(737424, 737426, 737429), "days since 0001-01-01 00:00:00",
               calendar = "proleptic_gregorian")
    },
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 2, expected_note = "2020-01-01, 2020-01-03, 2020-01-06"
  ),
  list(
    id = "T11", group = "time", from = "#90",
    what = "10-minute steps as fractional days since 1970-01-01",
    build = function(p) {
      base <- as.numeric(as.Date("2020-06-15"))
      .mk_time(p, base + (0:143) / 144, "days since 1970-01-01")
    },
    query = list(lon = 0, lat = 40,
                 time = .utc("2020-06-15") + (0:143) * 600, exact = TRUE),
    expected = 0,
    expected_note = "count of the 144 steps that fail an exact match"
  ),
  list(
    id = "T12", group = "time", from = "#91",
    what = "repeated timestamp: days since 2020-01-01, steps 0, 1, 1, 2",
    build = function(p) .mk_time(p, c(0, 1, 1, 2), "days since 2020-01-01"),
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 4, expected_note = ""
  ),
  list(
    id = "T13", group = "time", from = "#92",
    what = "UNLIMITED time dimension",
    build = function(p) {
      .mk_time(p, c(0, 2, 5), "days since 2020-01-01", unlim = TRUE)
    },
    query = list(lon = 0, lat = 40, time = .utc("2020-01-03")),
    expected = 2,
    expected_note = "#92 itself is a caller-opened ncdf4 connection"
  ),
  list(
    id = "T14", group = "time", from = "test (bounds)",
    what = "monthly means with time_bnds",
    build = function(p) {
      .mk_time(p, c(15.5, 45), "days since 2020-01-01",
               bounds = matrix(c(0, 31, 31, 60), 2))
    },
    query = list(lon = 0, lat = 40, time = .utc("2020-02-01")),
    expected = 2,
    expected_note = "bounds [Jan 1, Feb 1), [Feb 1, Mar 1)"
  ),
  # Depth ---------------------------------------------------------------------
  list(
    id = "D01", group = "depth", from = "test fixture",
    what = "depth (m), positive down, 0, 50, 100; query 60 m",
    build = function(p) .mk_depth(p, "depth", c(0, 50, 100), "m", "down"),
    query = list(lon = 0, lat = 40, depth = 60),
    expected = 2, expected_note = "nearest layer is 50 m"
  ),
  list(
    id = "D02", group = "depth", from = "test (units)",
    what = "depth in cm: 0, 5000, 10000; query 60 m",
    build = function(p) .mk_depth(p, "depth", c(0, 5000, 10000), "cm", "down"),
    query = list(lon = 0, lat = 40, depth = 60),
    expected = 2, expected_note = "5000 cm = 50 m"
  ),
  list(
    id = "D03", group = "depth", from = "test (direction)",
    what = "elevation, positive up: 0, -50, -100; query 60 m",
    build = function(p) {
      .mk_depth(p, "elevation", c(0, -50, -100), "m", "up", "elevation", "Z")
    },
    query = list(lon = 0, lat = 40, depth = 60),
    expected = 2,
    expected_note = "-50 m elevation = 50 m depth (Copernicus Zarr layout)"
  ),
  list(
    id = "D04", group = "depth", from = "test (storage order)",
    what = "depth stored deepest first: 100, 50, 0; query 60 m",
    build = function(p) .mk_depth(p, "depth", c(100, 50, 0), "m", "down"),
    query = list(lon = 0, lat = 40, depth = 60),
    expected = 2, expected_note = "layer 2 is 50 m"
  ),
  list(
    id = "D05", group = "depth", from = "#93",
    what = "altitude axis with axis = Z, standard_name = altitude",
    build = function(p) {
      .mk_depth(p, "altitude", c(0, 10), "m", NULL, "altitude", "Z")
    },
    query = list(lon = 0, lat = 40, depth = 0),
    expected = 1, expected_note = "recognised as the vertical axis"
  ),
  list(
    id = "D06", group = "depth", from = "test (bounds)",
    what = "depth bounds [0,40], [40,100]; centres 5, 55; query 35 m",
    build = function(p) {
      .mk_depth(p, "depth", c(5, 55), "m", "down", "depth", "Z",
                bounds = matrix(c(0, 40, 40, 100), 2))
    },
    query = list(lon = 0, lat = 40, depth = 35),
    expected = 1,
    expected_note = "35 m lies in [0, 40]; nearest centre is layer 2"
  ),
  list(
    id = "D07", group = "depth", from = "Copernicus layout",
    what = "depth with standard_name = depth, axis = Z; query 60 m",
    build = function(p) {
      .mk_depth(p, "depth", c(0, 50, 100), "m", "down", "depth", "Z")
    },
    query = list(lon = 0, lat = 40, depth = 60),
    expected = 2, expected_note = "as D01, with full CF labels"
  ),
  # Horizontal ----------------------------------------------------------------
  list(
    id = "H01", group = "horizontal", from = "test (conventions)",
    what = "longitudes 0..360 (0, 90, 180, 270); query lon -90",
    build = function(p) .mk_lon(p, c(0, 90, 180, 270)),
    query = list(lon = -90, lat = 40),
    expected = 4, expected_note = "-90 = 270"
  ),
  list(
    id = "H02", group = "horizontal", from = "test (conventions)",
    what = "longitudes -180..180 (-180, -90, 0, 90); query lon 270",
    build = function(p) .mk_lon(p, c(-180, -90, 0, 90)),
    query = list(lon = 270, lat = 40),
    expected = 2, expected_note = "270 = -90"
  ),
  list(
    id = "H03", group = "horizontal", from = "test (explicit bounds)",
    what = "lon centres 0, 1 with bounds [-0.2,0.2], [0.2,1.8]; query 0.4",
    build = function(p) {
      .mk_lon(p, c(0, 1), bounds = matrix(c(-0.2, 0.2, 0.2, 1.8), 2))
    },
    query = list(lon = 0.4, lat = 40),
    expected = 2,
    expected_note = "inside cell 2 by bounds; nearest centre is cell 1"
  ),
  list(
    id = "H04", group = "horizontal", from = "#93",
    what = "x labelled latitude/Y, y labelled longitude/X",
    build = function(p) {
      x <- ncdim_def("x", "degrees", c(10, 20))
      y <- ncdim_def("y", "degrees", c(30, 40))
      n <- nc_create(p, ncvar_def("temp", "1", list(x, y), -999))
      ncvar_put(n, "temp", matrix(1:4, 2))
      ncatt_put(n, "x", "standard_name", "latitude")
      ncatt_put(n, "x", "axis", "Y")
      ncatt_put(n, "y", "standard_name", "longitude")
      ncatt_put(n, "y", "axis", "X")
      nc_close(n)
    },
    query = list(lon = 40, lat = 10),
    expected = 3, expected_note = ""
  ),
  list(
    id = "H05", group = "horizontal", from = "#93",
    what = "2D nav_lon/nav_lat with units only (NEMO)",
    build = function(p) {
      x <- ncdim_def("x", "", 1:2)
      y <- ncdim_def("y", "", 1:2)
      n <- nc_create(p, list(
        ncvar_def("temp", "1", list(x, y), -999),
        ncvar_def("nav_lon", "degrees_east", list(x, y), -999),
        ncvar_def("nav_lat", "degrees_north", list(x, y), -999)
      ))
      ncvar_put(n, "temp", matrix(1:4, 2))
      ncvar_put(n, "nav_lon", matrix(c(100, 101, 100.5, 101.5), 2))
      ncvar_put(n, "nav_lat", matrix(c(10, 10.2, 11, 11.2), 2))
      ncatt_put(n, "temp", "coordinates", "nav_lon nav_lat")
      nc_close(n)
    },
    query = list(lon = 100.5, lat = 11),
    expected = 3, error_ok = TRUE,
    expected_note = "3 via nav_lon/nav_lat or a clear error, not x/y as lon/lat"
  ),
  list(
    id = "H06", group = "horizontal", from = "#95",
    what = "coarse longitudes 0, 200, 300; query lon 300",
    build = function(p) .mk_lon(p, c(0, 200, 300), vals = c(10, 20, 30)),
    query = list(lon = 300, lat = 40),
    expected = 30, expected_note = ""
  ),
  list(
    id = "H07", group = "horizontal", from = "test (storage order)",
    what = "latitude stored north to south (41, 40); query lat 40",
    build = function(p) {
      lon <- ncdim_def("lon", "degrees_east", c(0, 1))
      lat <- ncdim_def("lat", "degrees_north", c(41, 40))
      n <- nc_create(p, ncvar_def("temp", "1", list(lon, lat), -9999))
      ncvar_put(n, "temp", matrix(1:4, 2))
      nc_close(n)
    },
    query = list(lon = 0, lat = 40),
    expected = 3, expected_note = "row 2 is lat 40"
  ),
  # Values --------------------------------------------------------------------
  list(
    id = "V01", group = "values", from = "Copernicus layout",
    what = "short with scale_factor 0.001, add_offset 20, _FillValue",
    build = function(p) {
      lon <- ncdim_def("lon", "degrees_east", c(0, 1))
      lat <- ncdim_def("lat", "degrees_north", c(40, 41))
      n <- nc_create(p, ncvar_def("temp", "degrees_C", list(lon, lat),
                                  missval = -32767L, prec = "short"))
      ncvar_put(n, "temp", matrix(c(1234L, -32767L, 0L, 5L), 2))
      ncatt_put(n, "temp", "scale_factor", 0.001)
      ncatt_put(n, "temp", "add_offset", 20)
      nc_close(n)
    },
    query = list(lon = 0, lat = 40),
    expected = 21.234,
    expected_note = "1234 * 0.001 + 20; cell (1, 40) beside it is NA"
  ),
  # Multi-file ----------------------------------------------------------------
  list(
    id = "M01", group = "multi-file", from = "#94",
    what = "time series plus a static file",
    multi = TRUE, expected = 11,
    expected_note = "selection between files is package logic"
  ),
  list(
    id = "M02", group = "multi-file", from = "#96",
    what = "cadence borrowed from an unrelated tile",
    multi = TRUE, expected = NA_real_,
    expected_note = "selection between files is package logic"
  )
)
