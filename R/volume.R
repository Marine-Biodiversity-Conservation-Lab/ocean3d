# R/volume.R
#
# Volume for the package's 3D representations. The generic is declared in
# R/AllGenerics.R; this file holds the arithmetic.
#
# SpatEnvelope and SpatVoxel both describe a 3D domain over a 2D grid, but they
# store depth in dual roles -- the envelope as a per-cell value, the voxel as
# the layer index -- so the vertical extent is derived differently in each and
# each gets its own method rather than a shared SpatVolume one. See
# [SpatVolume-class].
#
# The internals both this file and R/intersect.R rely on (occupancy, presence,
# cell areas, grid and depth checks) live in R/utils-3d.R. Where two domains
# overlap is intersect_3d()'s job (R/intersect.R); volume() measures the result
# like any other domain.

# Internal: what both methods return from their per-cell volume raster (km³,
# NA outside the domain): the raster itself, or its total. The generic has
# already checked `per_cell`.
.volume_result <- function(vol_rast, per_cell) {
  if (!per_cell) {
    # terra sums an all-NA raster to NA; an empty domain holds no water.
    total <- terra::global(vol_rast, "sum", na.rm = TRUE)[[1]]
    return(if (is.na(total)) 0 else total)
  }
  names(vol_rast) <- "volume"
  methods::as(vol_rast, "SpatRaster")  # remove the 3D class
}

# Internal: the one producer of `ocean3d_error_volume_per_level`. The generic
# checks that `per_level` is a flag and comes with `per_cell`; the envelope,
# polygon and proxy methods name their objects in `what`, with the `reason`
# per-level volume doesn't apply and the `fix`.
.check_per_level <- function(per_level, per_cell, what = NULL,
                             reason = "have no depth levels",
                             fix = paste("Put them on a voxel's levels with",
                                         "envelope_to_voxel() first.")) {
  problem <- if (!rlang::is_bool(per_level)) {
    "`per_level` must be TRUE or FALSE."
  } else if (per_level && !per_cell) {
    paste("`per_level = TRUE` gives the volume of each cell at each depth",
          "level, so it needs `per_cell = TRUE`.")
  } else if (per_level && !is.null(what)) {
    glue::glue("{what} {reason}, so `per_level = TRUE` doesn't apply. {fix}")
  }
  if (!is.null(problem)) {
    rlang::abort(problem, class = "ocean3d_error_volume_per_level",
                 call = NULL)
  }
  invisible(TRUE)
}

#' @rdname volume
#' @export
setMethod("volume", "SpatEnvelope", function(x, per_cell = FALSE,
                                             per_level = FALSE, ...) {
  .no_extra_args(...)
  .check_per_level(per_level, per_cell, "Envelopes")
  .check_3d(x, "x")

  depth_extent <- .envelope_layer(x, "depth_max") -
    .envelope_layer(x, "depth_min")

  # Volume per cell in km³ (depth in m, converted to km). NA depths, outside
  # the domain, give NA.
  .volume_result(.cell_area_km2(x) * (depth_extent / 1000), per_cell)
})

#' @rdname volume
#' @export
setMethod(
  "volume", "SpatVoxel",
  function(x, per_cell = FALSE, per_level = FALSE,
           bounds = c("midpoint", "top")) {
    bounds <- match.arg(bounds)
    .check_3d(x, "x")

    occ <- .voxel_occupancy(x)
    slab <- .depth_slabs(depths(x), bounds, depth_bounds(x))
    # Summing thickness over occupied levels only is what makes an interior
    # gap cost volume, unlike voxel_to_envelope() plus the envelope formula.
    # The raster stays on the left of the operator: terra mishandles
    # `numeric <op> SpatRaster`. The vector is recycled layer-wise.
    occupied_m <- occ * (slab$upper - slab$lower)
    if (per_level) {
      # Each occupied level's volume, NA where a level is unoccupied, as
      # everywhere a shape is absent (SPEC.md, "Absent is NA").
      per_level_km3 <- terra::mask(.cell_area_km2(x) * (occupied_m / 1000),
                                   occ, maskvalues = 0)
      return(as_voxel(per_level_km3, depths = depths(x), varname = "volume",
                      depth_bounds = depth_bounds(x)))
    }
    # The column's thickness times the cell area, in that order, so the result
    # matches an envelope of the same column to the last bit. A cell with no
    # occupied level would read 0; masking makes it NA, as outside an
    # envelope.
    vol_rast <- terra::ifel(.footprint(x),
                            .cell_area_km2(x) * (sum(occupied_m) / 1000), NA)

    .volume_result(vol_rast, per_cell)
  }
)

# SpatVector3D polygons are measured from their outlines, not a grid
# (SPEC.md, "Summarise: volume calculations"): geodesic area times the depth
# range, with terra::expanse(). It agrees with terra::cellSize() for a polygon
# drawn along a cell's edges to within about 0.02% at low latitudes: the
# polygon's edges are geodesics, while cells are bounded by parallels.
# Points and lines have no volume.
#' @rdname volume
#' @export
setMethod("volume", "SpatVector3D", function(x, per_cell = FALSE,
                                             per_level = FALSE, ...) {
  .no_extra_args(...)
  .check_per_level(per_level, per_cell, "Polygons")
  if (x@type != "polygons") {
    rlang::abort(glue::glue(
      "These are {x@type}, which have no volume. Volume is defined for ",
      "shapes with a depth range: voxels, envelopes and polygons. ",
      "buffer_3d() turns {x@type} into polygons."
    ), class = "ocean3d_error_volume_type", call = NULL)
  }
  if (per_cell) {
    rlang::abort(
      "Polygons have no cells, so `per_cell = TRUE` doesn't apply.",
      class = "ocean3d_error_volume_per_cell", call = NULL
    )
  }
  f <- x@features
  if (nrow(f) == 0L) return(0)
  # With a raster surface, thickness times the area inside each polygon,
  # cell by cell on the surface's grid: an estimate at its resolution.
  if (length(x@surfaces)) {
    cells <- .surface_cells(f, lapply(x@surfaces, terra::unwrap),
                            exact = TRUE)
    return(sum(cells$fraction * cells$area_km2 *
                 (cells$bottom - cells$top), na.rm = TRUE) / 1000)
  }
  area_km2 <- terra::expanse(terra::vect(f$geometry), unit = "km")
  sum(area_km2 * (f$depth_max - f$depth_min) / 1000)
})

# Rejected here rather than by S4's default "unable to find an inherited
# method", which says nothing about how to build the right object.
#' @rdname volume
#' @export
setMethod("volume", "SpatRaster",
          function(x, per_cell = FALSE, per_level = FALSE, ...) {
            .check_3d(x, "x")
          })

# A SpatProxy is measured one time step at a time, reading each step through
# .proxy_to_voxel(); a step larger than the read budget is read in groups of
# depth layers. Every group is read with the proxy's depth edges, so each
# layer stands for the same slab whichever layers are read with it, and the
# groups' volumes add up to the whole step's. Only the per-cell result keeps
# a raster per step.
#' @rdname volume
#' @param fun SpatProxy only. Which values count as occupied, as in
#'   [occupied()]: `NULL` (default) counts every non-`NA` value; otherwise a
#'   function of a layer's values returning `TRUE` where the cell is
#'   occupied, for example `function(v) v > 10 & v < 15` for the water
#'   between 10 and 15 °C.
#' @section Memory:
#' `volume()` on a [SpatProxy-class] holds at most one time step in memory.
#' A step larger than the read budget (option `ocean3d.max_read_mb`, 1024 MB
#' by default) is read in groups of depth layers, and their volumes are
#' added. A step takes cells x layers x 8 bytes. Each read holds at least one
#' depth layer of the whole grid, even when that is larger than the budget.
#' Set the option lower on a machine with little memory, or higher for
#' faster reads. With `per_cell = TRUE`, the result itself holds one layer
#' per time step.
#' @section Overlap through time:
#' `intersect_3d()`, `intersects_3d()` and `mask()` don't take a
#' `SpatProxy`. To measure the overlap of a time-varying dataset with a
#' shape, read each step and intersect it:
#' `volume(intersect_3d(as_voxel(x, time = i), shape))` for each step `i`.
#' @export
setMethod("volume", "SpatProxy", function(x, per_cell = FALSE,
                                          per_level = FALSE, fun = NULL,
                                          ...) {
  .no_extra_args(...)
  .check_per_level(per_level, per_cell, "A SpatProxy",
                   reason = "is measured one time step at a time",
                   fix = paste("Read a step with as_voxel() and measure it",
                               "with `per_level = TRUE`."))
  if (!is.null(fun) && !is.function(fun)) {
    rlang::abort("`fun` must be NULL or a function of a layer's values.",
                 class = "ocean3d_error_volume_fun", call = NULL)
  }
  budget <- getOption("ocean3d.max_read_mb", 1024)
  if (!is.numeric(budget) || length(budget) != 1L || !is.finite(budget) ||
        budget <= 0) {
    rlang::abort(glue::glue(
      "The option `ocean3d.max_read_mb` must be one positive number of ",
      "megabytes; got {format(budget)}."
    ), class = "ocean3d_error_read_budget", call = NULL)
  }
  variable <- x@variable
  axis <- x@axes$time
  n_steps <- if (is.null(axis)) 1L else variable$dim_lengths[[axis$dim]]

  # Groups of depth layers (first and last index, storage order) that fit
  # the read budget, at least one layer each; a proxy without a depth axis
  # gets one group of nothing, and .proxy_to_voxel() stops for it.
  depth_dim <- variable$axis_dims$depth
  n_depth <- if (is.null(depth_dim)) 1L else variable$dim_lengths[[depth_dim]]
  n_cells <- prod(variable$dim_lengths[unlist(variable$axis_dims[c("lon",
                                                                   "lat")])])
  per_group <- max(1L, min(n_depth, floor(budget * 1024^2 / (8 * n_cells))))
  groups <- if (is.null(depth_dim)) {
    list(NULL)
  } else {
    lapply(seq(1L, n_depth, by = per_group), function(first) {
      c(first, min(first + per_group - 1L, n_depth))
    })
  }

  # Each step's total (or per-cell raster), and its date from the voxel the
  # builder dated, NA when it has none. Each step is asked for by number.
  times <- rep(NA_real_, n_steps)
  results <- lapply(seq_len(n_steps), function(step) {
    total <- NULL
    for (span in groups) {
      v <- .proxy_to_voxel(x, if (is.null(axis)) NULL else step,
                           variable$name, "nearest", NULL, depth = span)
      times[step] <<- as.numeric(.voxel_time(v))
      if (!is.null(fun)) v <- occupied(v, fun)
      part <- volume(v, per_cell = per_cell)
      # Added up as each group is read, so only the running total is kept.
      # A cell stays NA until a group has water in it, as volume() gives.
      total <- if (is.null(total)) {
        part
      } else if (per_cell) {
        terra::cover(terra::cover(total + part, total), part)
      } else {
        total + part
      }
    }
    total
  })
  times <- .POSIXct(times, tz = "UTC")

  if (!per_cell) {
    return(data.frame(time = times, volume_km3 = unlist(results)))
  }
  out <- terra::rast(results)
  names(out) <- rep("volume", n_steps)
  if (!anyNA(times)) terra::time(out) <- times
  out
})
