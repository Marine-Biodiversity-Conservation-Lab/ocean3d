#' Vertical profiles for voxel building
#'
#' A vertical profile decides how the value a cell carries is spread down the
#' depth levels that cell occupies. [envelope_to_voxel()] takes one through its
#' `profile` argument, as a function — either one of these, or one you write.
#'
#' @section Writing a profile:
#'
#' A profile is called once per [envelope_to_voxel()] call, with three
#' arguments, positionally:
#'
#' \describe{
#'   \item{`ind`}{Logical SpatRaster, one layer per depth, `TRUE` where a
#'     cell's envelope overlaps the slab that depth stands for — from the depth
#'     down to the next one, the deepest standing for itself alone. A level can
#'     therefore be occupied without its own value lying in the envelope. Cells
#'     that are `NA` in the envelope are `NA` here.}
#'   \item{`depths`}{Numeric vector of the depths those layers stand for,
#'     shallow to deep. Always `terra::nlyr(ind)` long.}
#'   \item{`n_depths`}{Single-layer SpatRaster counting the levels each cell
#'     occupies, i.e. `sum(ind)`. Passed in rather than left to be recomputed,
#'     since the caller already needs it.}
#' }
#'
#' It returns the weight that multiplies `ind`, in one of three forms: a single
#' number applying to the whole grid, a single-layer SpatRaster giving one
#' weight per cell, or a `terra::nlyr(ind)`-layer SpatRaster giving one per cell
#' per level. Rasters must be on the grid of `ind`.
#'
#' Because a profile is handed rasters and returns rasters, write it in terra
#' algebra rather than in terms of one cell at a time. Unoccupied cells are
#' masked out by the caller afterwards, so a profile need not zero them — but
#' it should avoid dividing by a zero level count, since `Inf` is harder to read
#' in a partial result than `NA`.
#'
#' A profile that conserves does so *per cell*: the weights it returns sum to 1
#' within each cell, whatever that cell's interval.
#'
#' @param ind Logical SpatRaster of occupancy, one layer per depth.
#' @param depths Numeric vector of depths, shallow to deep, one per layer of
#'   `ind`.
#' @param n_depths Single-layer SpatRaster counting each cell's occupied levels.
#'
#' @returns A weight to multiply `ind` by: a single number, or a SpatRaster on
#'   the grid of `ind` with either one layer or one per depth.
#'
#' @examples
#' fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
#'                   ymax = 1)
#' terra::values(fp) <- c(1, 1)
#' dmax <- terra::setValues(terra::rast(fp), c(100, 300))
#' e <- as_envelope(fp, depth_min = 0, depth_max = dmax)
#' depths <- c(0, 100, 200, 300)
#'
#' # the whole value at every occupied level
#' terra::values(envelope_to_voxel(e, depths, values = 12,
#'                                 profile = profile_flat))
#'
#' # an even share instead, so each cell's levels sum back to its value
#' terra::values(envelope_to_voxel(e, depths, values = 12,
#'                                 profile = profile_equal))
#'
#' # a profile of your own: weight the shallow levels twice as heavily as the
#' # deep ones, still conserving each cell's total
#' surface_weighted <- function(ind, depths, n_depths) {
#'   w <- terra::rast(lapply(depths,
#'                           function(d) ifelse(d <= 100, 2, 1) * ind[[1]]))
#'   w <- w * ind
#'   w / sum(w)
#' }
#' terra::values(envelope_to_voxel(e, depths, values = 12,
#'                                 profile = surface_weighted))
#' @name voxel_profiles
NULL

#' @describeIn voxel_profiles The whole value at every occupied level. This is
#'   what `profile = NULL` means in [envelope_to_voxel()]; it exists as a
#'   function so the no-profile case is a profile like any other rather than a
#'   special case in the caller.
#' @export
profile_flat <- function(ind, depths, n_depths) {
  1
}

#' @describeIn voxel_profiles An even share at each occupied level, so the
#'   levels sum back to the cell's value. The share depends on how many levels
#'   the cell occupies, which varies across the grid wherever the envelope
#'   does: a cell spanning three levels gets a third at each, one spanning ten
#'   gets a tenth.
#' @export
profile_equal <- function(ind, depths, n_depths) {
  terra::ifel(n_depths > 0, 1 / n_depths, NA)
}

# Internal: run a profile and check what it gives back.
#
# `profile` is now arbitrary caller code, so its return is checked before it
# reaches the arithmetic — a wrong-sized raster would otherwise surface as a
# terra recycling error with nothing pointing at the profile as the cause.
.profile_weights <- function(profile, ind, depths, n_depths) {
  w <- withCallingHandlers(
    profile(ind, depths, n_depths),
    error = function(e) {
      stop("`profile` failed when called with (ind, depths, n_depths): ",
           conditionMessage(e), call. = FALSE)
    }
  )

  if (inherits(w, "SpatRaster")) {
    n <- terra::nlyr(w)
    if (n != 1 && n != terra::nlyr(ind)) {
      stop("`profile` must return a SpatRaster with 1 layer or one per depth ",
           "(", terra::nlyr(ind), "); got ", n, ".", call. = FALSE)
    }
    if (!terra::compareGeom(w, ind, stopOnError = FALSE)) {
      stop("`profile` must return a SpatRaster on the same grid (CRS, extent, ",
           "resolution) as `x`.", call. = FALSE)
    }
    return(w)
  }

  if (!is.numeric(w) || length(w) != 1 || is.na(w)) {
    stop("`profile` must return a single non-missing number or a SpatRaster; ",
         "got ", paste(class(w), collapse = "/"), " of length ", length(w),
         ".", call. = FALSE)
  }
  w
}
