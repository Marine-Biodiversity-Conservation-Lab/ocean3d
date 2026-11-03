# 3D points, lines and polygons (SpatVector3D, #114). The class is defined
# in R/AllClasses.R. Point and line vertices are stored as a table, so a
# line's depth and time at each vertex can't be lost by terra or sf
# operations; polygons keep their outlines in an sf geometry column. sf input
# is read in; there is no conversion back to sf (#114).

#' Build 3D points, lines or polygons
#'
#' Builds a [SpatVector3D-class] from a table with one row per vertex, or from
#' an sf object or `terra::SpatVector` with POINT, LINESTRING, POLYGON or
#' MULTIPOLYGON geometries.
#'
#' @param x A data frame (or tibble, matrix, named list) with one row per
#'   vertex, or an sf object or `terra::SpatVector`. sf points take their
#'   depth and time from columns; sf lines take them from the Z and M
#'   coordinates; polygons take them from `depth_min`, `depth_max` and `time`.
#' @param type `"points"` (default) or `"lines"`, for tables and sf points.
#'   With `"lines"`, the rows that share an `id` form one line, in row order.
#'   sf lines are always lines, and sf polygons always polygons.
#' @param depth Column with the depth of each vertex, in positive metres
#'   downward. Not used for sf lines, whose depth is minus the Z coordinate
#'   (Z is height, as in other GIS tools), or for polygons.
#' @param time Optional column with the time of each vertex: `Date`, `POSIXct`
#'   or text, read as in [extract_to_point()]. For sf lines, the M coordinate
#'   in seconds since 1970-01-01 UTC, if present. For polygons, one column
#'   (each polygon at a single time) or two (the start and end of each
#'   polygon's time interval). With dates, the interval covers whole UTC
#'   days, including the last: 1 to 15 February includes all of 15 February.
#'   The dates are stored as given, at midnight UTC, with `time_resolution`
#'   `"daily"`. With hours, the start and end are exact times.
#' @param id Column with feature IDs. Required for lines from a table; for
#'   points and polygons, each row is a feature and IDs default to the row
#'   numbers.
#' @param lon,lat Longitude and latitude columns of a table.
#' @param depth_min,depth_max Polygons only, and required for them: the top
#'   and bottom of each polygon, in positive metres downward. Each is a
#'   single number for all polygons (e.g. `depth_min = 0, depth_max = 200`)
#'   or the name of a column of `x`.
#'
#' @return A [SpatVector3D-class]. Columns other than the coordinates, depth,
#'   time and ID are kept with each vertex (points and lines) or polygon.
#' @seealso [geomtype()]
#' @examples
#' obs <- data.frame(lon = c(150, 151), lat = c(-20, -21), depth = c(10, 35),
#'                   date = as.Date(c("2020-01-01", "2020-01-02")))
#' as_vector3d(obs, time = "date")
#'
#' dive <- data.frame(tag = "A", lon = c(150, 150.1, 150.2),
#'                    lat = -20, depth = c(0, 120, 15))
#' track <- as_vector3d(dive, type = "lines", id = "tag")
#' geomtype(track)
#'
#' # A protected area from the surface to 200 m.
#' mpa <- sf::st_sf(
#'   name = "reef",
#'   geometry = sf::st_sfc(sf::st_polygon(list(rbind(
#'     c(150, -20), c(151, -20), c(151, -21), c(150, -21), c(150, -20)
#'   ))), crs = 4326)
#' )
#' as_vector3d(mpa, depth_min = 0, depth_max = 200)
#' @export
as_vector3d <- function(x, type = c("points", "lines"), depth = "depth",
                        time = NULL, id = NULL, lon = "lon", lat = "lat",
                        depth_min = NULL, depth_max = NULL) {
  if (inherits(x, "SpatVector")) x <- sf::st_as_sf(x)
  geometry <- if (inherits(x, "sf")) {
    unique(as.character(sf::st_geometry_type(x)))
  }
  if (length(geometry) && all(geometry %in% c("POLYGON", "MULTIPOLYGON"))) {
    return(.vector3d_from_sf_polygons(x, id, time, depth_min, depth_max))
  }
  if (!is.null(depth_min) || !is.null(depth_max)) {
    .abort_vector3d(paste0(
      "`depth_min` and `depth_max` are for polygons. Points and lines ",
      "have a depth at each vertex instead (`depth`)."
    ))
  }
  if (identical(geometry, "LINESTRING")) {
    return(.vector3d_from_sf_lines(x, id))
  }
  type <- match.arg(type)
  table <- .normalize_points_input(x, lon, lat)
  if (inherits(table, "sf")) table <- sf::st_drop_geometry(table)
  .check_required_cols(table, c(lon, lat, depth, time, id), "x")
  .check_numeric_columns(table, c(lon, lat, depth))
  if (type == "lines" && is.null(id)) {
    .abort_vector3d(
      "Lines need `id`: the column that says which line each vertex is on."
    )
  }

  n <- nrow(table)
  feature_id <- if (is.null(id)) {
    as.character(seq_len(n))
  } else {
    as.character(table[[id]])
  }

  if (is.null(time)) {
    times <- .POSIXct(rep(NA_real_, n), tz = "UTC")
    resolution <- NA_character_
  } else {
    prepared <- .prepare_time_input(
      table[[time]], "auto", input_name = glue::glue("`x${time}`")
    )
    times <- prepared$values
    resolution <- prepared$resolution
  }

  vertices <- data.frame(
    feature_id = feature_id,
    vertex = stats::ave(seq_len(n), feature_id, FUN = seq_along),
    lon = as.numeric(table[[lon]]),
    lat = as.numeric(table[[lat]]),
    depth = as.numeric(table[[depth]]),
    time = times,
    stringsAsFactors = FALSE
  )
  extra <- table[setdiff(names(table), c(lon, lat, depth, time, id))]
  vertices <- .add_attributes(vertices, extra)
  # Features in order of first appearance, vertices in row order.
  vertices <- vertices[order(match(feature_id, unique(feature_id)),
                             vertices$vertex), , drop = FALSE]
  .new_vector3d(type, resolution, vertices = vertices)
}

# sf LINESTRING input: depth from minus Z (Z is height), time from M.
.vector3d_from_sf_lines <- function(x, id) {
  x <- .check_sf_input(x)
  coords <- sf::st_coordinates(x)
  if (!"Z" %in% colnames(coords)) {
    .abort_vector3d(paste0(
      "sf lines need a depth at each vertex: give the geometry Z ",
      "coordinates (Z is height, so depth = -Z)."
    ))
  }
  attributes <- sf::st_drop_geometry(x)
  feature_id <- if (is.null(id)) {
    as.character(seq_len(nrow(x)))
  } else {
    .check_required_cols(attributes, id, "x")
    as.character(attributes[[id]])
  }
  if (anyNA(feature_id) || any(!nzchar(feature_id)) ||
        anyDuplicated(feature_id)) {
    .abort_vector3d(glue::glue(
      "Each sf line needs its own non-missing ID in `{id}`; otherwise ",
      "separate lines would be joined into one."
    ))
  }
  row <- coords[, "L1"]
  has_time <- "M" %in% colnames(coords)
  vertices <- data.frame(
    feature_id = feature_id[row],
    vertex = stats::ave(seq_along(row), row, FUN = seq_along),
    lon = coords[, "X"], lat = coords[, "Y"], depth = -coords[, "Z"],
    time = .POSIXct(if (has_time) coords[, "M"] else NA_real_, tz = "UTC"),
    stringsAsFactors = FALSE
  )
  extra <- attributes[setdiff(names(attributes), id)]
  vertices <- .add_attributes(vertices, extra[row, , drop = FALSE])
  .new_vector3d("lines", if (has_time) "datetime" else NA_character_,
                vertices = vertices)
}

# sf POLYGON or MULTIPOLYGON input: one depth range and optional time
# interval per polygon, kept with its outline in an sf geometry column.
.vector3d_from_sf_polygons <- function(x, id, time, depth_min, depth_max) {
  if (is.null(depth_min) || is.null(depth_max)) {
    .abort_vector3d(paste0(
      "Polygons need `depth_min` and `depth_max`: numbers in metres (for ",
      "example depth_min = 0, depth_max = 200) or columns of `x`."
    ))
  }
  if (length(time) > 2L) {
    .abort_vector3d(paste0(
      "`time` for polygons is one column (a single time) or two (the start ",
      "and end of each polygon's interval)."
    ))
  }
  for (value in list(depth_min, depth_max)) {
    if (length(value) != 1L || !(is.numeric(value) || is.character(value))) {
      .abort_vector3d(paste0(
        "`depth_min` and `depth_max` must each be one number or the name ",
        "of a column of `x`."
      ))
    }
  }
  x <- .check_sf_input(x)
  table <- sf::st_drop_geometry(x)
  n <- nrow(table)
  depth_columns <- c(if (is.character(depth_min)) depth_min,
                     if (is.character(depth_max)) depth_max)
  .check_required_cols(table, c(id, time, depth_columns), "x")
  .check_numeric_columns(table, depth_columns)

  if (is.null(time)) {
    start <- end <- .POSIXct(rep(NA_real_, n), tz = "UTC")
    resolution <- NA_character_
  } else {
    # Each column is read on its own: c() on a Date and a POSIXct column
    # would convert the second to the first's class and lose its hours.
    prepared <- lapply(unique(time), function(column) {
      .prepare_time_input(table[[column]], "auto",
                          input_name = glue::glue("`x${column}`"))
    })
    resolution <- unique(vapply(prepared, `[[`, "", "resolution"))
    if (length(resolution) > 1L) {
      .abort_vector3d(glue::glue(
        "`{time[1L]}` and `{time[2L]}` must both be dates or both have ",
        "hours, so the interval has one resolution. Convert one of them ",
        "first, for example with as.Date() or as.POSIXct()."
      ))
    }
    start <- prepared[[1L]]$values
    end <- prepared[[length(prepared)]]$values
  }

  features <- data.frame(
    feature_id = if (is.null(id)) seq_len(n) else table[[id]],
    depth_min = if (is.character(depth_min)) table[[depth_min]] else depth_min,
    depth_max = if (is.character(depth_max)) table[[depth_max]] else depth_max,
    time_start = start,
    time_end = end,
    stringsAsFactors = FALSE
  )
  features$feature_id <- as.character(features$feature_id)
  features$geometry <- sf::st_geometry(x)
  extra <- table[setdiff(names(table), c(id, time, depth_columns))]
  .new_vector3d("polygons", resolution,
                features = .add_attributes(features, extra))
}

# sf lines and polygons: transform any CRS other than WGS84, geographic or
# not, and stop on empty geometries, which would vanish from the vertex
# table or have no outline.
.check_sf_input <- function(x) {
  if (sf::st_crs(x) != sf::st_crs(4326)) {
    if (is.na(sf::st_crs(x))) {
      .abort_vector3d("sf input needs a coordinate reference system.")
    }
    x <- sf::st_transform(x, 4326)
  }
  empty <- which(sf::st_is_empty(x))
  if (length(empty)) {
    rows <- if (length(empty) == 1L) {
      "Row {empty} of `x` has an empty geometry"
    } else {
      "Rows {toString(empty)} of `x` have empty geometries"
    }
    .abort_vector3d(glue::glue(rows, ". Remove the empty rows first."))
  }
  x
}

# as.numeric() on a factor gives its level codes, not its values, so
# coordinate and depth columns must already be numeric.
.check_numeric_columns <- function(table, columns) {
  for (column in columns) {
    if (!is.numeric(table[[column]])) {
      .abort_vector3d(glue::glue(
        "`{column}` must be numeric, but it is ",
        "{class(table[[column]])[1L]}. Convert it first, for example with ",
        "as.numeric(as.character(x)) for a factor."
      ))
    }
  }
}

# Adds the input's other columns, one value per vertex or polygon. A column
# named like one SpatVector3D uses would be lost, so it stops.
.add_attributes <- function(table, extra) {
  clash <- intersect(names(extra), names(table))
  if (length(clash)) {
    .abort_vector3d(glue::glue(
      "Column names {toString(clash)} are used by SpatVector3D. Rename ",
      "them, or pass the column as an argument if that is what it holds ",
      "(for example `time = \"time\"` or `id = \"feature_id\"`)."
    ))
  }
  cbind(table, extra)
}

# The class's validity rules are the only checks; this reports a failure as
# an ocean3d error.
.new_vector3d <- function(type, resolution, vertices = data.frame(),
                          features = data.frame()) {
  rownames(vertices) <- NULL
  rownames(features) <- NULL
  tryCatch(
    methods::new("SpatVector3D", type = type, vertices = vertices,
                 features = features, time_resolution = resolution),
    error = function(e) {
      .abort_vector3d(
        sub("^invalid class .*? object: ", "", conditionMessage(e))
      )
    }
  )
}

# The one producer of `ocean3d_error_vector3d`: input that can't become a
# SpatVector3D. Every check in as_vector3d() and the class's validity
# rules report through it.
.abort_vector3d <- function(message) {
  rlang::abort(message, class = "ocean3d_error_vector3d", call = NULL)
}

#' @rdname SpatVector3D-class
#' @param x A `SpatVector3D`.
#' @param ... Not used.
#' @export
setMethod("geomtype", "SpatVector3D", function(x, ...) x@type)

#' @rdname SpatVector3D-class
#' @export
setMethod("length", "SpatVector3D", function(x) {
  if (x@type == "polygons") return(nrow(x@features))
  length(unique(x@vertices$feature_id))
})

#' @rdname SpatVector3D-class
#' @param object A `SpatVector3D`.
#' @export
setMethod("show", "SpatVector3D", function(object) {
  cat("class      : SpatVector3D\n")
  cat("geometry   :", object@type, "\n")
  cat("features   :", length(object), "\n")
  if (object@type == "polygons") {
    f <- object@features
    depths <- c(min(f$depth_min), max(f$depth_max))
    times <- c(f$time_start, f$time_end)
  } else {
    v <- object@vertices
    cat("vertices   :", nrow(v), "\n")
    depths <- range(v$depth)
    times <- v$time
  }
  cat("depth (m)  :", depths[1L], "to", depths[2L], "\n")
  if (all(is.na(times))) {
    cat("time       : none\n")
  } else {
    shown <- format(range(times, na.rm = TRUE), tz = "UTC",
                    usetz = object@time_resolution == "datetime")
    cat("time       :", shown[1L], "to", shown[2L], "\n")
  }
  invisible(object)
})

#' @rdname SpatVector3D-class
#' @param row.names,optional Not used.
#' @return `as.data.frame()`: the vertex table for points and lines, or the
#'   polygon table, with its sf geometry column, for polygons.
#' @exportS3Method base::as.data.frame
as.data.frame.SpatVector3D <- function(x, row.names = NULL, # nolint
                                       optional = FALSE,
                                       ...) {
  if (x@type == "polygons") x@features else x@vertices
}

#' @rdname volume
#' @export
setMethod("volume", "SpatVector3D", function(x, ...) {
  message <- if (x@type == "polygons") {
    "volume() doesn't support SpatVector3D polygons yet."
  } else {
    glue::glue(
      "These are {x@type}, which have no volume. Volume is defined for ",
      "shapes with a depth range: voxels, envelopes and polygons."
    )
  }
  rlang::abort(message, class = "ocean3d_error_volume_type", call = NULL)
})
