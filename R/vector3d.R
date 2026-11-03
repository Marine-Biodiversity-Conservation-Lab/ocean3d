# 3D points and lines (SpatVector3D, #114). The class is defined in
# R/AllClasses.R. Vertices are stored as a table, so a line's depth and time
# at each vertex can't be lost by terra or sf operations. sf input is read
# in; there is no conversion back to sf (#114).

#' Build 3D points or lines
#'
#' Builds a [SpatVector3D-class] from a table with one row per vertex, or from
#' an sf object with POINT or LINESTRING geometries.
#'
#' @param x A data frame (or tibble, matrix, named list) with one row per
#'   vertex, or an sf object. sf points take their depth and time from
#'   columns; sf lines take them from the Z and M coordinates.
#' @param type `"points"` (default) or `"lines"`, for tables and sf points.
#'   With `"lines"`, the rows that share an `id` form one line, in row order.
#'   sf lines are always lines.
#' @param depth Column with the depth of each vertex, in positive metres
#'   downward. Not used for sf lines, whose depth is minus the Z coordinate
#'   (Z is height, as in other GIS tools).
#' @param time Optional column with the time of each vertex: `Date`, `POSIXct`
#'   or text, read as in [extract_to_point()]. For sf lines, the M coordinate
#'   in seconds since 1970-01-01 UTC, if present.
#' @param id Column with feature IDs. Required for lines from a table; for
#'   points, each row is a feature and IDs default to the row numbers.
#' @param lon,lat Longitude and latitude columns of a table.
#'
#' @return A [SpatVector3D-class]. Columns other than the coordinates, depth,
#'   time and ID are kept with each vertex.
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
#' @export
as_vector3d <- function(x, type = c("points", "lines"), depth = "depth",
                        time = NULL, id = NULL, lon = "lon", lat = "lat") {
  if (inherits(x, "sf") &&
        identical(unique(as.character(sf::st_geometry_type(x))),
                  "LINESTRING")) {
    return(.vector3d_from_sf_lines(x, id))
  }
  type <- match.arg(type)
  table <- .normalize_points_input(x, lon, lat)
  if (inherits(table, "sf")) table <- sf::st_drop_geometry(table)
  .check_required_cols(table, c(lon, lat, depth, time, id), "x")
  # as.numeric() on a factor gives its level codes, not its values.
  for (column in c(lon, lat, depth)) {
    if (!is.numeric(table[[column]])) {
      .abort_vector3d(glue::glue(
        "`{column}` must be numeric, but it is ",
        "{class(table[[column]])[1L]}. Convert it first, for example with ",
        "as.numeric(as.character(x)) for a factor."
      ))
    }
  }
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
  vertices <- .add_vertex_attributes(vertices, extra)
  # Features in order of first appearance, vertices in row order.
  vertices <- vertices[order(match(feature_id, unique(feature_id)),
                             vertices$vertex), , drop = FALSE]
  .new_vector3d(type, vertices, resolution)
}

# sf LINESTRING input: depth from minus Z (Z is height), time from M.
.vector3d_from_sf_lines <- function(x, id) {
  # Any CRS other than WGS84, geographic or not, is transformed.
  if (sf::st_crs(x) != sf::st_crs(4326)) {
    if (is.na(sf::st_crs(x))) {
      .abort_vector3d("sf input needs a coordinate reference system.")
    }
    x <- sf::st_transform(x, 4326)
  }
  # st_coordinates() gives no rows for an empty line, so it would vanish.
  empty <- which(sf::st_is_empty(x))
  if (length(empty)) {
    rows <- if (length(empty) == 1L) {
      "Row {empty} of `x` has an empty geometry"
    } else {
      "Rows {toString(empty)} of `x` have empty geometries"
    }
    .abort_vector3d(glue::glue(
      rows, ". Each line ",
      "needs at least two vertices; remove the empty rows first."
    ))
  }
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
  vertices <- .add_vertex_attributes(vertices, extra[row, , drop = FALSE])
  rownames(vertices) <- NULL
  .new_vector3d("lines", vertices,
                if (has_time) "datetime" else NA_character_)
}

# Adds the input's other columns to the vertex table, one value per vertex.
# A column named like a vertex-table column would be lost, so it stops.
.add_vertex_attributes <- function(vertices, extra) {
  clash <- intersect(names(extra), names(vertices))
  if (length(clash)) {
    .abort_vector3d(glue::glue(
      "Column names {toString(clash)} are used by SpatVector3D. Rename ",
      "them, or for a table, pass the column as an argument if that is what ",
      "it holds (for example `time = \"time\"` or `id = \"feature_id\"`)."
    ))
  }
  cbind(vertices, extra)
}

# The class's validity rules are the only checks; this reports a failure as
# an ocean3d error.
.new_vector3d <- function(type, vertices, resolution) {
  rownames(vertices) <- NULL
  tryCatch(
    methods::new("SpatVector3D", type = type, vertices = vertices,
                 time_resolution = resolution),
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
  length(unique(x@vertices$feature_id))
})

#' @rdname SpatVector3D-class
#' @param object A `SpatVector3D`.
#' @export
setMethod("show", "SpatVector3D", function(object) {
  v <- object@vertices
  cat("class      : SpatVector3D\n")
  cat("geometry   :", object@type, "\n")
  cat("features   :", length(object), "\n")
  cat("vertices   :", nrow(v), "\n")
  cat("depth (m)  :", min(v$depth), "to", max(v$depth), "\n")
  if (all(is.na(v$time))) {
    cat("time       : none\n")
  } else {
    shown <- format(range(v$time, na.rm = TRUE), tz = "UTC",
                    usetz = object@time_resolution == "datetime")
    cat("time       :", shown[1L], "to", shown[2L], "\n")
  }
  invisible(object)
})

#' @rdname SpatVector3D-class
#' @param row.names,optional Not used.
#' @exportS3Method base::as.data.frame
as.data.frame.SpatVector3D <- function(x, row.names = NULL, # nolint
                                       optional = FALSE,
                                       ...) {
  x@vertices
}

#' @rdname volume
#' @export
setMethod("volume", "SpatVector3D", function(x, ...) {
  rlang::abort(glue::glue(
    "These are {x@type}, which have no volume. Volume is defined for shapes ",
    "with a depth range: voxels, envelopes and polygons."
  ), class = "ocean3d_error_volume_type", call = NULL)
})
