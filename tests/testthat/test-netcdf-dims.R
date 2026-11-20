# Axis detection from CF metadata (#93), tested on the shared helpers
# directly: the core object (#113) reuses `.detect_dims()`, while the
# extraction functions that call it may be retired (#116; AGENTS.md, "Tests").
# The files use ordinary dimension coordinates made with ncdim_def(), which
# ncdf4 lists in `nc$dim` rather than `nc$var`.

# A 2 x 2 file whose dimensions are named `x` and `y` and labelled the other
# way round by their `standard_name` and `axis` attributes. `group` puts the
# variable and its dimensions in a netCDF-4 group.
make_swapped_axes_file <- function(group = NULL) {
  name <- function(x) if (is.null(group)) x else paste0(group, "/", x)
  f <- tempfile(fileext = ".nc")
  x <- ncdf4::ncdim_def(name("x"), "degrees", c(10, 20))
  y <- ncdf4::ncdim_def(name("y"), "degrees", c(30, 40))
  nc <- ncdf4::nc_create(
    f, ncdf4::ncvar_def(name("v"), "1", list(x, y), -999),
    force_v4 = !is.null(group)
  )
  ncdf4::ncvar_put(nc, name("v"), matrix(1:4, 2))
  ncdf4::ncatt_put(nc, name("x"), "standard_name", "latitude")
  ncdf4::ncatt_put(nc, name("x"), "axis", "Y")
  ncdf4::ncatt_put(nc, name("y"), "standard_name", "longitude")
  ncdf4::ncatt_put(nc, name("y"), "axis", "X")
  ncdf4::nc_close(nc)
  f
}

# A NEMO-style file: index dimensions `x` and `y`, and 2D `nav_lon` /
# `nav_lat` coordinates named by the variable's `coordinates` attribute.
# `standard_name = FALSE` leaves the coordinates with units only. The
# `coordinates` attribute holds names relative to the variable's group.
make_curvilinear_file <- function(standard_name = FALSE, group = NULL) {
  name <- function(x) if (is.null(group)) x else paste0(group, "/", x)
  f <- tempfile(fileext = ".nc")
  x <- ncdf4::ncdim_def(name("x"), "", 1:2)
  y <- ncdf4::ncdim_def(name("y"), "", 1:2)
  nc <- ncdf4::nc_create(f, list(
    ncdf4::ncvar_def(name("v"), "1", list(x, y), -999),
    ncdf4::ncvar_def(name("nav_lon"), "degrees_east", list(x, y), -999),
    ncdf4::ncvar_def(name("nav_lat"), "degrees_north", list(x, y), -999)
  ), force_v4 = !is.null(group))
  ncdf4::ncvar_put(nc, name("v"), matrix(1:4, 2))
  ncdf4::ncvar_put(nc, name("nav_lon"), matrix(c(10, 20, 10, 20), 2))
  ncdf4::ncvar_put(nc, name("nav_lat"), matrix(c(40, 40, 50, 50), 2))
  ncdf4::ncatt_put(nc, name("v"), "coordinates", "nav_lon nav_lat")
  if (standard_name) {
    ncdf4::ncatt_put(nc, name("nav_lon"), "standard_name", "longitude")
    ncdf4::ncatt_put(nc, name("nav_lat"), "standard_name", "latitude")
  }
  ncdf4::nc_close(nc)
  f
}

detect <- function(f, var = "v", ...) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc))
  .detect_dims(nc, var, ...)
}

# A file with one variable, `temp`, on `dims`: each a list of the dimension's
# name, units and values. `atts` holds attributes for the dimensions, as
# list(<dimension> = list(<attribute> = <value>)).
make_dims_file <- function(dims, atts = list()) {
  f <- tempfile(fileext = ".nc")
  defs <- lapply(dims, function(d) ncdf4::ncdim_def(d[[1]], d[[2]], d[[3]]))
  nc <- ncdf4::nc_create(f, ncdf4::ncvar_def("temp", "1", defs, -9999))
  for (dim in names(atts)) {
    for (name in names(atts[[dim]])) {
      ncdf4::ncatt_put(nc, dim, name, atts[[dim]][[name]])
    }
  }
  ncdf4::nc_close(nc)
  f
}

lon_dim <- list("lon", "degrees_east", c(0, 1))
lat_dim <- list("lat", "degrees_north", c(40, 41))

test_that("CF metadata and dimension names identify the axes", {
  # `args` go to .detect_dims(); `dims` are the axes it should find.
  cases <- list(
    list(label = "axis and standard_name on dimension coordinates",
         file = make_swapped_axes_file(), var = "v",
         dims = list(lon = "y", lat = "x")),
    list(label = "axis attributes on dimensions inside a group",
         file = make_swapped_axes_file(group = "g"), var = "g/v",
         dims = list(lon = "g/y", lat = "g/x")),
    list(label = "an axis attribute alone",
         file = make_dims_file(
           list(list("a", "1", c(0, 1)), list("b", "1", c(40, 41))),
           atts = list(a = list(axis = "X"), b = list(axis = "Y"))
         ),
         dims = list(lon = "a", lat = "b")),
    list(label = "an altitude dimension with axis = 'Z'",
         file = make_dims_file(
           list(lon_dim, lat_dim, list("altitude", "m", c(0, -10))),
           atts = list(altitude = list(standard_name = "altitude", axis = "Z"))
         ),
         dims = list(depth = "altitude")),
    list(label = "valid_time",
         file = make_dims_file(list(
           lon_dim, lat_dim,
           list("valid_time", "days since 2020-01-01", c(0, 1))
         )),
         dims = list(time = "valid_time")),
    list(label = "non-standard names with CF units",
         file = make_dims_file(list(
           list("horizontal_a", "degrees_east", c(0, 1)),
           list("horizontal_b", "degrees_north", c(40, 41)),
           list("forecast_step", "days since 2020-01-01", c(0, 1))
         )),
         dims = list(lon = "horizontal_a", lat = "horizontal_b",
                     time = "forecast_step")),
    list(label = "an explicit dimension",
         file = make_grid_file(depth = NULL, time = NULL),
         args = list(lon_dim = "lon"), dims = list(lon = "lon"))
  )
  for (case in cases) {
    on.exit(unlink(case$file), add = TRUE)
    dims <- do.call(detect, c(list(case$file, case$var %||% "temp"),
                              case$args))
    for (axis in names(case$dims)) {
      expect_identical(dims[[axis]], case$dims[[axis]],
                       label = paste0(case$label, ": ", axis))
    }
  }
})

test_that("grids .detect_dims() can't read are rejected", {
  # `args` go to .detect_dims(). The error contains each of `messages`, as
  # fixed strings if `fixed` is TRUE.
  curvilinear <- "ocean3d_error_curvilinear_grid"
  curvilinear_messages <- c("two-dimensional geographic coordinate `nav_lon`",
                            "regular longitude/latitude grid")
  cases <- list(
    list(label = "2D coordinates identified by units",
         file = make_curvilinear_file(), var = "v", class = curvilinear,
         messages = curvilinear_messages, fixed = TRUE),
    list(label = "2D coordinates identified by standard_name",
         file = make_curvilinear_file(standard_name = TRUE), var = "v",
         class = curvilinear, messages = curvilinear_messages, fixed = TRUE),
    list(label = "2D coordinates inside a group",
         file = make_curvilinear_file(group = "g"), var = "g/v",
         class = curvilinear),
    # An explicit dimension doesn't bypass the curvilinear check.
    list(label = "2D coordinates with explicit dimensions",
         file = make_curvilinear_file(), var = "v",
         args = list(lon_dim = "x", lat_dim = "y"), class = curvilinear),
    list(label = "an unknown non-singleton dimension",
         file = make_dims_file(
           list(lon_dim, lat_dim, list("member", "1", 1:3))
         ),
         messages = "Unrecognised non-singleton dimension"),
    list(label = "projected coordinates",
         file = make_dims_file(list(list("x", "m", c(1000, 2000)),
                                    list("y", "m", c(3000, 4000)))),
         messages = "projected coordinates"),
    list(label = "two longitude dimensions",
         file = make_dims_file(list(
           lon_dim, list("longitude", "degrees_east", c(10, 11)), lat_dim
         )),
         messages = "Ambiguous lon axis"),
    list(label = "an explicit dimension the variable doesn't have",
         file = make_grid_file(depth = NULL, time = NULL),
         args = list(lon_dim = "longitude"),
         messages = "not a dimension of the selected variable")
  )
  for (case in cases) {
    on.exit(unlink(case$file), add = TRUE)
    err <- expect_error(
      do.call(detect, c(list(case$file, case$var %||% "temp"), case$args)),
      class = case$class, label = case$label
    )
    for (message in case$messages) {
      expect_match(conditionMessage(err), message,
                   fixed = isTRUE(case$fixed), label = case$label)
    }
  }
})

# Dimension detection and reading (.detect_dims(), .get_dim_values()) ----------

test_that("dimension helpers read coordinates and report missing names", {
  f <- make_grid_file(lon = c(0, 1))
  on.exit(unlink(f), add = TRUE)
  # Add a `longitude` coordinate variable on the lon dimension.
  nc <- ncdf4::nc_open(f, write = TRUE)
  nc <- ncdf4::ncvar_add(nc, ncdf4::ncvar_def(
    "longitude", "degrees_east", list(nc$dim$lon), prec = "double"
  ))
  ncdf4::ncvar_put(nc, "longitude", c(10, 11))
  ncdf4::nc_close(nc)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_equal(as.vector(.get_dim_values(nc, "longitude")), c(10, 11))
  expect_null(.get_dim_values(nc, NULL))
  expect_error(.get_dim_values(nc, "missing"),
               "Could not read dimension or coordinate variable")
  expect_error(.detect_dims(nc, "missing"),
               "Available variables are: temp, longitude")

})

test_that("dimension detection uses only the selected variable", {
  f <- tempfile(fileext = ".nc")
  on.exit(unlink(f), add = TRUE)

  lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
  lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
  other_lon <- ncdf4::ncdim_def(
    "longitude", "degrees_east", c(10, 11)
  )
  time <- ncdf4::ncdim_def(
    "time", "days since 2020-01-01", c(0, 1)
  )

  static <- ncdf4::ncvar_def(
    "static", "1", list(lon, lat), -9999
  )
  dynamic <- ncdf4::ncvar_def(
    "dynamic", "1", list(other_lon, lat, time), -9999
  )

  nc <- ncdf4::nc_create(f, list(static, dynamic))
  ncdf4::ncvar_put(nc, "static", matrix(1:4, 2, 2))
  ncdf4::ncvar_put(
    nc, "dynamic", array(1:8, c(2, 2, 2))
  )
  ncdf4::nc_close(nc)

  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  dims <- .detect_dims(nc, "static")

  expect_identical(dims$lon, "lon")
  expect_identical(dims$lat, "lat")
  expect_null(dims$time)
  expect_identical(dims$var_dims, c("lon", "lat"))
})

test_that("unknown singleton dimensions are identified", {
  f <- make_dims_file(list(lon_dim, lat_dim, list("member", "1", 1)))
  on.exit(unlink(f), add = TRUE)

  dims <- detect(f, "temp")

  expect_identical(dims$singleton_dims, "member")

  sc <- .build_start_count(
    dims$var_dims,
    list(
      lon = list(start = 1L, count = 1L),
      lat = list(start = 1L, count = 1L)
    ),
    singleton_dims = dims$singleton_dims
  )

  expect_identical(sc$start, c(1L, 1L, 1L))
  expect_identical(sc$count, c(1L, 1L, 1L))
})

test_that("unassigned dimensions cannot be read silently", {
  expect_error(
    .build_start_count(
      c("lon", "lat", "member"),
      list(
        lon = list(start = 1L, count = 1L),
        lat = list(start = 1L, count = 1L)
      )
    ),
    "No extraction index was assigned"
  )
})
