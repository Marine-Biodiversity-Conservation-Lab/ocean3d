# 3D points and lines (SpatVector3D, #114). The class is defined in
# R/AllClasses.R. Vertices are stored as a table, so a line's depth and time
# at each vertex can't be lost by terra or sf operations.

#' Build 3D points or lines
#'
#' Builds a [SpatVector3D-class] from a table with one row per vertex, or from
#' an sf object with POINT geometries.
#'
#' @param x A data frame (or tibble, matrix, named list) with one row per
#'   vertex, or an sf object with POINT geometries, which take their depth
#'   and time from columns.
#' @param type `"points"` (default) or `"lines"`, for tables and sf points.
#'   With `"lines"`, the rows that share an `id` form one line, in row order.
#' @param depth Column with the depth of each vertex, in positive metres
#'   downward.
#' @param time Optional column with the time of each vertex: `Date`, `POSIXct`
#'   or text, read as in [extract_to_point()].
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
  type <- match.arg(type)
  table <- .normalize_points_input(x, lon, lat)
  if (inherits(table, "sf")) table <- sf::st_drop_geometry(table)
  .check_required_cols(table, c(lon, lat, depth, time, id), "x")
  if (type == "lines" && is.null(id)) {
    rlang::abort(
      "Lines need `id`: the column that says which line each vertex is on.",
      class = "ocean3d_error_vector3d", call = NULL
    )
  }

  n <- nrow(table)
  feature_id <- if (is.null(id)) {
    as.character(seq_len(n))
  } else {
    as.character(table[[id]])
  }
  if (anyNA(feature_id) || any(!nzchar(feature_id))) {
    rlang::abort("Feature IDs must be non-missing and non-empty.",
                 class = "ocean3d_error_vector3d", call = NULL)
  }
  if (type == "points" && anyDuplicated(feature_id)) {
    rlang::abort("Each point needs its own ID; IDs repeat.",
                 class = "ocean3d_error_vector3d", call = NULL)
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
  clash <- intersect(names(extra), names(vertices))
  if (length(clash)) {
    rlang::abort(glue::glue(
      "Column names {toString(clash)} are used by SpatVector3D. Pass the ",
      "column as an argument if that is what it holds (for example ",
      "`time = \"time\"` or `id = \"feature_id\"`), or rename it."
    ), class = "ocean3d_error_vector3d", call = NULL)
  }
  vertices <- cbind(vertices, extra)
  # Features in order of first appearance, vertices in row order.
  vertices <- vertices[order(match(feature_id, unique(feature_id)),
                             vertices$vertex), , drop = FALSE]
  .new_vector3d(type, vertices, resolution)
}

# The class's validity rules are the only checks; this reports a failure as
# an ocean3d error.
.new_vector3d <- function(type, vertices, resolution) {
  rownames(vertices) <- NULL
  tryCatch(
    methods::new("SpatVector3D", type = type, vertices = vertices,
                 time_resolution = resolution),
    error = function(e) {
      rlang::abort(
        sub("^invalid class .*? object: ", "", conditionMessage(e)),
        class = "ocean3d_error_vector3d", call = NULL
      )
    }
  )
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
