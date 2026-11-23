# show() methods for the raster 3D classes (#61). Without them, printing a
# SpatVoxel or SpatEnvelope fell through to terra's SpatRaster method and
# reported `class : SpatRaster`, so a reader couldn't tell which 3D
# representation they were holding.

# Internal: terra's own description of the grid (size, resolution, extent,
# CRS, sources, names, value ranges), for the show() methods below. Its
# `class` line is replaced by the 3D class, and its `depth` and `time` lines
# by the lines the methods print with the depth edges. Printing must never
# throw, so terra's output is read from a plain copy.
.grid_lines <- function(object) {
  lines <- utils::capture.output(show(methods::as(object, "SpatRaster")))
  lines[!grepl("^(class|depth|time) ", lines)]
}

# Internal: the 3D lines of a show() method, in terra's layout, or a `status`
# line when the object isn't a valid one (a terra operation changed its
# layers), so printing never stops. `lines` is a function so it only reads
# the object once it is known to be valid.
.show_3d <- function(object, class, lines) {
  ok <- methods::validObject(object, test = TRUE)
  shown <- if (isTRUE(ok)) lines() else
    c(status = glue::glue("not a valid {class}: {ok}"))
  shown <- c(class = class, shown)
  cat(sprintf("%-12s: %s", names(shown), shown), .grid_lines(object),
      sep = "\n")
  invisible(object)
}

#' @rdname SpatVoxel-class
#' @param object A `SpatVoxel`. `show()` prints the class, the variable, the
#'   depth levels, whether the voxel carries depth edges (see
#'   [depth_bounds()]), the time step, and then terra's description of the
#'   grid.
#' @export
setMethod("show", "SpatVoxel", function(object) {
  .show_3d(object, "SpatVoxel", function() {
    d <- depths(object)
    n <- length(d)
    levels <- if (n == 1L) "1 level" else glue::glue("{n} levels")
    shown <- if (n <= 6L) toString(d) else glue::glue("{min(d)} to {max(d)}")
    # The variable's name from terra's own field, which as_voxel() sets.
    name <- unique(terra::varnames(object))
    name <- name[!is.na(name) & nzchar(name)]
    c(
      variable = if (length(name)) toString(name) else "none",
      `depth (m)` = glue::glue("{shown} ({levels})"),
      `depth edges` = if (is.null(depth_bounds(object))) {
        "none; `bounds` applies"
      } else {
        "from depth_bounds"
      },
      time = if (isTRUE(terra::timeInfo(object)$time)) {
        format(terra::time(object)[1L], "%Y-%m-%d %H:%M:%S", tz = "UTC",
               usetz = TRUE)
      } else {
        "none"
      }
    )
  })
})

#' @rdname SpatEnvelope-class
#' @param object A `SpatEnvelope`. `show()` prints the class, the depth range
#'   (the shallowest `depth_min` to the deepest `depth_max`), how many cells
#'   the envelope occupies, and then terra's description of the grid.
#' @export
setMethod("show", "SpatEnvelope", function(object) {
  .show_3d(object, "SpatEnvelope", function() {
    # Summarised by terra, so a large grid isn't pulled into R.
    limit <- function(layer, fun) {
      terra::global(.envelope_layer(object, layer), fun, na.rm = TRUE)[1, 1]
    }
    cells <- terra::global(.footprint(object), "sum")[1, 1]
    shallowest <- limit("depth_min", "min")
    c(
      `depth (m)` = if (is.na(shallowest)) "none" else
        glue::glue("{shallowest} to {limit('depth_max', 'max')}"),
      cells = glue::glue("{cells} of {terra::ncell(object)} occupied")
    )
  })
})
