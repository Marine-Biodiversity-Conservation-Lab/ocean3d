# show() methods for the raster 3D classes (#61). Without them, printing a
# SpatVoxel or SpatEnvelope fell through to terra's SpatRaster method and
# reported `class : SpatRaster`, so a reader couldn't tell which 3D
# representation they were holding.

# Internal: the 3D lines of a show() method, in terra's layout, or a `status`
# line when the object isn't a valid one (a terra operation changed its
# layers), so printing never stops. `lines` is a function so it only reads
# the object once it is known to be valid. terra's own description of the
# grid (size, resolution, extent, CRS, sources, names, value ranges) follows,
# read from a plain copy, without the lines named in `drop`, which the 3D
# lines replace.
.show_3d <- function(object, class, lines, drop = "class") {
  ok <- methods::validObject(object, test = TRUE)
  shown <- if (isTRUE(ok)) lines() else
    c(status = glue::glue("not a valid {class}: {ok}"))
  shown <- c(class = class, shown)
  grid <- utils::capture.output(show(methods::as(object, "SpatRaster")))
  pattern <- paste0("^(", paste(drop, collapse = "|"), ") ")
  cat(sprintf("%-12s: %s", names(shown), shown), grid[!grepl(pattern, grid)],
      sep = "\n")
  invisible(object)
}

#' @rdname SpatVoxel-class
#' @param object A `SpatVoxel`. `show()` prints the class, the depth levels,
#'   whether the voxel carries depth edges (see [depth_bounds()]) and the
#'   time, then terra's description of the grid.
#' @export
setMethod("show", "SpatVoxel", function(object) {
  .show_3d(object, "SpatVoxel", drop = c("class", "depth", "time"), function() {
    d <- depths(object)
    n <- length(d)
    levels <- if (n == 1L) "1 level" else glue::glue("{n} levels")
    shown <- if (n <= 6L) toString(d) else glue::glue("{min(d)} to {max(d)}")
    # terra::time() is a date-time, a date, or a number for a year, a year
    # and month or a plain step, which format() shows as it is. Layers with
    # different times show their range.
    time <- range(terra::time(object))
    shown_time <- if (inherits(time, "POSIXt")) {
      format(time, "%Y-%m-%d %H:%M:%S", tz = "UTC", usetz = TRUE)
    } else {
      format(time)
    }
    c(
      `depth (m)` = glue::glue("{shown} ({levels})"),
      `depth edges` = if (is.null(depth_bounds(object))) "none" else
        "from depth_bounds",
      time = if (!isTRUE(terra::timeInfo(object)$time)) {
        "none"
      } else {
        paste(unique(shown_time), collapse = " to ")
      }
    )
  })
})

#' @rdname SpatEnvelope-class
#' @param object A `SpatEnvelope`. `show()` prints the class, the depth range
#'   (the shallowest `depth_min` to the deepest `depth_max`) and how many
#'   cells the envelope occupies, then terra's description of the grid.
#' @export
setMethod("show", "SpatEnvelope", function(object) {
  .show_3d(object, "SpatEnvelope", function() {
    # Both limits in one pass by terra, so a large grid isn't pulled into R.
    # An empty envelope gives non-finite limits. remove the 3D class
    limits <- terra::minmax(methods::as(object, "SpatRaster"), compute = TRUE)
    shallowest <- limits["min", "depth_min"]
    cells <- terra::global(.footprint(object), "sum")[1, 1]
    counts <- format(c(cells, terra::ncell(object)), scientific = FALSE,
                     trim = TRUE)
    c(
      `depth (m)` = if (!is.finite(shallowest)) "none" else
        glue::glue("{shallowest} to {limits['max', 'depth_max']}"),
      cells = glue::glue("{counts[1L]} of {counts[2L]} occupied")
    )
  })
})
