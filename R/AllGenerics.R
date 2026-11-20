# R/AllGenerics.R
#
# Generics for the operations that are meaningful on either 3D representation.
# The methods live in R/volume.R and R/intersect.R; only the dispatch contract
# is declared here.

#' @importFrom methods setGeneric setMethod
NULL

#' 3D volume of a rasterized domain
#'
#' Measures the occupied volume of a [SpatEnvelope-class] or a
#' [SpatVoxel-class], as a total over the whole grid or cell by cell. Both
#' representations describe a 3D domain over a 2D grid, but they store depth
#' in dual roles, so each computes the vertical extent differently:
#'
#' \describe{
#'   \item{[SpatEnvelope-class]}{Depth is the cell value, and the interval is
#'     solid by construction. Volume is
#'     `sum(cell_area * (depth_max - depth_min))` over present cells.}
#'   \item{[SpatVoxel-class]}{Depth is the layer index. Each occupied voxel
#'     contributes its cell area times the thickness of the slab its depth
#'     level stands for, so interior gaps are excluded rather than filled.}
#' }
#'
#' A voxel volume is quantized to the levels the voxel was sampled at, and a
#' voxel built by [envelope_to_voxel()] counts a level as occupied on *any*
#' overlap with the envelope. So a range of 0-100 m discretized onto
#' 0, 50, 100, 150, 200 m occupies three levels spanning 0-150 m and reports
#' more water than the envelope it came from. Coarse levels overstate; the
#' error shrinks as the levels tighten around the range. Where an envelope's
#' limits fall exactly on levels under `bounds = "top"`, the two agree.
#'
#' To measure the space two domains share, measure their intersection:
#' `volume(intersect_3d(a, b))`.
#'
#' With one envelope and one voxel, [intersect_3d()] places the envelope on
#' the voxel's depth levels, so the overlap is measured on those levels while
#' `volume()` of the envelope uses its own depth limits. To compare the two on
#' the same levels, measure the envelope on them too:
#' `volume(envelope_to_voxel(a, depths(b), bounds = bounds,`
#' `depth_bounds = depth_bounds(b)), bounds = bounds)`. Pass the same `bounds`
#' to [intersect_3d()] and `volume()`; otherwise the overlap and the input
#' volumes are measured on different slabs and don't add up.
#'
#' @param x A [SpatEnvelope-class] or [SpatVoxel-class]. A bare SpatRaster is
#'   rejected; build one with [as_envelope()], [vect_to_envelope()] or
#'   [as_voxel()] first. A netCDF file isn't accepted: load it as a
#'   `SpatRaster` (e.g. with [woa_load_nc()]) and wrap it with [as_voxel()].
#' @param per_cell `FALSE` (default) for the total volume, `TRUE` for the
#'   volume of each cell.
#' @param ... Arguments for the [SpatVoxel-class] method, which are an error
#'   for an envelope.
#' @param bounds SpatVoxel only. What each depth level stands for vertically:
#'   `"midpoint"` (default) puts the slab edges halfway to each neighbouring
#'   level (the World Ocean Atlas convention);
#'   `"top"` treats the level as the top of its slab, so the slab runs down to
#'   the next level. Under `"top"` the deepest level has no next level and
#'   contributes no volume. Outer edges are clamped to `range(depths)` under
#'   both. Same argument, same meaning, as in [envelope_to_voxel()]. Not used
#'   for a voxel that carries its own depth edges (see [SpatVoxel-class]).
#'
#' @returns With `per_cell = FALSE`, a number: the total volume in km³. With
#'   `per_cell = TRUE`, a single-layer `SpatRaster` named `volume` on the grid
#'   of `x`: the volume of each cell in km³, `NA` where `x` is absent (no
#'   depth interval, or no occupied level). Its sum is the total.
#'
#' @seealso [intersect_3d()] for the space two domains share, which `volume()`
#'   then measures; [intersects_3d()] for whether two domains share any depth
#'   in each cell; [SpatVolume-class] for why the two representations
#'   dispatch separately.
#'
#' @examples
#' fp <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2000,
#'                   ymin = 0, ymax = 2000,
#'                   crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
#' terra::values(fp) <- c(1, 1, 1, NA)
#'
#' # 2.5D: three 1 km² cells, each 0-200 m => 3 * 1 * 0.2 = 0.6 km³.
#' e <- as_envelope(fp, depth_min = 0, depth_max = 200)
#' volume(e)
#'
#' # 3D: the same domain on three levels. The slabs are 50 m, 100 m and 50 m,
#' # so the discretized volume matches exactly here.
#' v <- envelope_to_voxel(e, depths = c(0, 100, 200))
#' volume(v)
#'
#' # An interior gap is volume a voxel can drop and an envelope cannot.
#' gappy <- v
#' gappy[[2]] <- terra::setValues(terra::rast(v[[2]]), NA)
#' gappy <- as_voxel(gappy)
#' volume(gappy)
#' volume(voxel_to_envelope(gappy))
#'
#' # Per cell: the fourth cell is outside the domain, so it is NA.
#' terra::values(volume(e, per_cell = TRUE))
#'
#' # Overlap with a second domain at 100-300 m: intersect, then measure.
#' b <- as_envelope(fp, depth_min = 100, depth_max = 300)
#' ov <- intersect_3d(e, b)
#' volume(ov)
#' terra::values(volume(ov, per_cell = TRUE))
#' @export
setGeneric("volume", function(x, per_cell = FALSE, ...) {
  # Checked once here, before a method computes anything.
  if (!rlang::is_bool(per_cell)) {
    rlang::abort("`per_cell` must be TRUE or FALSE.",
                 class = "ocean3d_error_per_cell", call = NULL)
  }
  standardGeneric("volume")
}, signature = "x")

#' The 3D space two objects share
#'
#' `intersect_3d()` returns the part of space that is inside both `x` and
#' `y`. It answers the question "what do these two share?". To ask only
#' "do they share anything?", use [intersects_3d()]. To keep one object's
#' values where the other is present, use [mask()].
#'
#' At least one of `x` and `y` must be a 3D object: a [SpatEnvelope-class] or
#' a [SpatVoxel-class]. The other can be a 3D object too, or a 2D object. A 2D
#' object restricts the result horizontally and leaves its depths alone.
#'
#' \describe{
#'   \item{Two envelopes}{Cell by cell, the result runs from the deeper
#'     `depth_min` to the shallower `depth_max`. Intervals that only touch
#'     share no water, so that cell is empty (`NA`).}
#'   \item{Two voxels}{Both must be sampled at the same depth levels. A level
#'     is in the result where both voxels occupy it.}
#'   \item{An envelope and a voxel}{The envelope is first placed on the
#'     voxel's depth levels with [envelope_to_voxel()]. The result is a voxel.}
#'   \item{A 3D object and a `SpatRaster`}{The raster is a footprint. Its
#'     non-`NA` cells are inside; its values are not read. It must have one
#'     layer and be on the same grid.}
#'   \item{A 3D object and polygons}{A `SpatVector`, `sf` or `sfc` object. The
#'     polygons are rasterised onto the 3D object's grid. A cell is inside when
#'     its centre is. Polygons in another CRS are projected first.}
#' }
#'
#' The order of `x` and `y` does not matter.
#'
#' The result is a domain, not a field. An envelope result carries the shared
#' depth interval. A voxel result carries presence: `1` where a level is
#' shared, `NA` elsewhere, in layers named `presence_depth=<value>`. The cell
#' values of a voxel input are not carried over. To keep them, use [mask()].
#'
#' @param x,y The two objects. At least one must be a [SpatEnvelope-class] or
#'   a [SpatVoxel-class]. The other can be either of those, a single-layer
#'   `SpatRaster` footprint, or polygons (`SpatVector`, `sf`, `sfc`). Rasters
#'   must share one grid (CRS, extent, resolution). netCDF files aren't
#'   accepted; load them with terra and wrap them with [as_voxel()] first.
#' @param ... Arguments for the voxel methods. They are an error when both
#'   inputs are envelopes.
#' @param bounds Only when one input is an envelope and the other a voxel. How
#'   the envelope is placed on the voxel's depth levels: `"midpoint"`
#'   (default) or `"top"`. See [envelope_to_voxel()]. Not used when the voxel
#'   carries its own depth edges (see [SpatVoxel-class]). To measure the
#'   result, pass the same `bounds` to [volume()].
#'
#' @returns A [SpatEnvelope-class] when both inputs are envelopes, or when one
#'   is an envelope and the other is 2D. Otherwise a [SpatVoxel-class] of
#'   presence. Either is on the grid of the 3D input.
#'
#' @seealso [intersects_3d()] for the yes/no form; [mask()] to keep a voxel's
#'   values inside a domain; [volume()] for the volume of the result, in
#'   total or per cell. When one input is a voxel, the envelope is placed on
#'   its levels, so compare the overlap with `volume(envelope_to_voxel(...))`
#'   rather than with the envelope's own volume.
#'
#' @examples
#' fp <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3000,
#'                   ymin = 0, ymax = 1000,
#'                   crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
#' terra::values(fp) <- c(1, 1, 1)
#'
#' # Two species: one at 0-100 m, one at 50-200 m everywhere.
#' a <- as_envelope(fp, depth_min = 0, depth_max = 100)
#' b <- as_envelope(fp, depth_min = 50, depth_max = 200)
#'
#' # They share 50-100 m in every cell.
#' shared <- intersect_3d(a, b)
#' terra::values(shared)
#' volume(shared)
#'
#' # With a voxel, the result is a presence voxel on the voxel's levels.
#' bv <- envelope_to_voxel(b, depths = c(0, 50, 100, 150, 200))
#' names(intersect_3d(a, bv))
#'
#' # With polygons, the envelope is kept where the polygon is.
#' poly <- terra::vect("POLYGON ((0 0, 2000 0, 2000 1000, 0 1000, 0 0))",
#'                     crs = terra::crs(fp))
#' terra::values(intersect_3d(a, poly))
#' @export
setGeneric("intersect_3d", function(x, y, ...) standardGeneric("intersect_3d"))

#' Do two objects share any 3D space?
#'
#' `intersects_3d()` tests, cell by cell, whether `x` and `y` overlap in 3D.
#' It is the yes/no form of [intersect_3d()] and takes the same inputs. It
#' computes only the presence pattern, with no cell areas and no volumes, so
#' it is the cheap choice for richness and tally maps.
#'
#' Each cell gets one of three answers:
#' \itemize{
#'   \item `TRUE`: both objects are present and their depths overlap.
#'   \item `FALSE`: both are present but their depths do not overlap, or only
#'     one of them is present.
#'   \item `NA`: neither is present. There is nothing to compare.
#' }
#'
#' This is how `terra` answers the 2D question for two rasters, so the result
#' sums and plots like any other boolean layer. Depth intervals that only touch
#' do not overlap. When `y` is a 2D footprint or polygons there is no depth
#' axis to compare, so the test is only whether both are present.
#'
#' The order of `x` and `y` does not matter.
#'
#' @inheritParams intersect_3d
#'
#' @returns A single-layer boolean `SpatRaster` named `intersects`, on the
#'   grid of the 3D input: `TRUE`, `FALSE` or `NA` per cell as above.
#'
#' @seealso [intersect_3d()] for the shared space itself; [mask()] to keep a
#'   voxel's values inside a domain.
#'
#' @examples
#' fp <- terra::rast(nrows = 1, ncols = 4, xmin = 0, xmax = 4000,
#'                   ymin = 0, ymax = 1000,
#'                   crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
#'
#' # Four cells. A is present in cells 1-3 at 0-100 m. B is present in cells
#' # 1, 2 and 4, at 50-200 m in cell 1 and 300-400 m elsewhere.
#' a <- as_envelope(terra::setValues(fp, c(1, 1, 1, NA)),
#'                  depth_min = 0, depth_max = 100)
#' b <- as_envelope(terra::setValues(fp, c(1, 1, NA, 1)),
#'                  depth_min = terra::setValues(fp, c(50, 300, NA, 300)),
#'                  depth_max = terra::setValues(fp, c(200, 400, NA, 400)))
#'
#' # Cell 1: overlap. Cell 2: both present, depths disjoint. Cell 3: A only.
#' # Cell 4: B only.
#' terra::values(intersects_3d(a, b))
#'
#' # Summing a stack of these gives a richness map. na.rm = TRUE counts each
#' # TRUE as 1 and each FALSE as 0.
#' richness <- sum(c(intersects_3d(a, b), intersects_3d(a, a)), na.rm = TRUE)
#' terra::values(richness)
#' @export
setGeneric("intersects_3d",
           function(x, y, ...) standardGeneric("intersects_3d"))
