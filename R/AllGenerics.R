# R/AllGenerics.R
#
# Generics for the operations that are meaningful on either 3D representation.
# The methods live in R/volume.R and R/intersect.R; only the dispatch contract
# is declared here.

#' @importFrom methods setGeneric setMethod
NULL

#' Total 3D volume of a rasterized domain
#'
#' Sums the occupied volume of a [SpatEnvelope-class] or a [SpatVoxel-class]
#' over the whole grid. Both representations describe a 3D domain over a 2D
#' grid, but they store depth in dual roles, so each computes the vertical
#' extent differently:
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
#' @param x A [SpatEnvelope-class] or [SpatVoxel-class]. A bare SpatRaster is
#'   rejected; build one with [as_envelope()], [vect_to_envelope()] or
#'   [as_voxel()] first. A netCDF file isn't accepted: load it as a
#'   `SpatRaster` (e.g. with [woa_load_nc()]) and wrap it with [as_voxel()].
#' @param ... Arguments for the [SpatVoxel-class] method, which are an error
#'   for an envelope.
#' @param bounds SpatVoxel only. What each depth level stands for vertically:
#'   `"top"` (default) treats the level as the top of its slab, so the slab
#'   runs down to the next level; `"midpoint"` puts the level in the middle,
#'   with edges halfway to each neighbour (the World Ocean Atlas convention).
#'   Under `"top"` the deepest level has no next level and contributes no
#'   volume. Outer edges are clamped to `range(depths)` under both. Same
#'   argument, same meaning, as in [envelope_to_voxel()].
#'
#' @returns Numeric of length 1. Total volume in km³.
#'
#' @seealso [calc_volume_overlap()] for the volume two domains share, and
#'   [intersect_3d()] for the shared domain itself; [SpatVolume-class] for why
#'   the two representations dispatch separately.
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
#' # 3D: the same domain on three levels. Under "top" the slabs are 100 m,
#' # 100 m and 0 m, so the discretized volume matches exactly here.
#' v <- envelope_to_voxel(e, depths = c(0, 100, 200))
#' volume(v)
#'
#' # An interior gap is volume a voxel can drop and an envelope cannot.
#' gappy <- v
#' gappy[[2]] <- terra::setValues(terra::rast(v[[2]]), NA)
#' gappy <- as_voxel(gappy)
#' volume(gappy)
#' volume(voxel_to_envelope(gappy))
#' @export
setGeneric("volume", function(x, ...) standardGeneric("volume"))

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
#'   the envelope is placed on the voxel's depth levels: `"top"` (default) or
#'   `"midpoint"`. See [envelope_to_voxel()].
#'
#' @returns A [SpatEnvelope-class] when both inputs are envelopes, or when one
#'   is an envelope and the other is 2D. Otherwise a [SpatVoxel-class] of
#'   presence. Either is on the grid of the 3D input.
#'
#' @seealso [intersects_3d()] for the yes/no form; [mask()] to keep a voxel's
#'   values inside a domain; [volume()] for the volume of the result;
#'   [calc_volume_overlap()], which is built on this function.
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

#' Per-cell 3D volume overlap between two rasterized domains
#'
#' Computes the depth interval and volume of each domain and of their
#' intersection, cell by cell. The intersection is [intersect_3d()]; this
#' function measures it. Accepts any combination of [SpatEnvelope-class] and
#' [SpatVoxel-class]; when the two differ, the envelope is discretized onto
#' the voxel's depth levels with [envelope_to_voxel()] rather than the voxel
#' being collapsed, because [voxel_to_envelope()] fills interior gaps and
#' would overstate the overlap.
#'
#' Two voxels must be sampled at the same depth levels, and all inputs must be
#' on the same grid. To ask only whether two domains overlap, use
#' [intersects_3d()]; it computes no volumes.
#'
#' @param x,y The two domains, each a [SpatEnvelope-class] or
#'   [SpatVoxel-class]. netCDF files aren't accepted; load them with terra
#'   and wrap them with [as_voxel()] first.
#' @param ... Arguments for the voxel methods, which are an error for a pair
#'   of envelopes.
#' @param bounds Voxel methods only. See [volume()]. Also governs how an
#'   envelope is discretized when the two inputs differ.
#'
#' @returns Multi-layer SpatRaster with 9 layers:
#'   \describe{
#'     \item{depth_min_a, depth_max_a}{Depth limits of `x` (m)}
#'     \item{depth_min_b, depth_max_b}{Depth limits of `y` (m)}
#'     \item{depth_min_overlap, depth_max_overlap}{Depth limits of the
#'       intersection (m). NA where the domains do not overlap in depth.}
#'     \item{volume_a, volume_b}{Per-cell volume of each domain (km³)}
#'     \item{volume_overlap}{Per-cell overlap volume (km³). `0` where both are
#'       present but their depths do not intersect, NA where either is
#'       absent.}
#'   }
#'
#'   For voxel inputs the depth layers are *summary bounds* — the shallowest
#'   and deepest occupied level. A voxel with an interior gap is not solid
#'   between them, so read the volume layers, not the depth layers, for
#'   occupied volume. The result is a plain SpatRaster: it is neither an
#'   envelope nor a voxel.
#'
#' @seealso [intersect_3d()] for the shared domain itself, and [volume()] for
#'   its total volume; [intersects_3d()] when only the presence of overlap is
#'   needed.
#'
#' @examples
#' fp <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2000,
#'                   ymin = 0, ymax = 1000,
#'                   crs = "+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m")
#' terra::values(fp) <- c(1, 1)
#'
#' a <- as_envelope(fp, depth_min = 0, depth_max = 100)
#' b <- as_envelope(fp, depth_min = 50, depth_max = 200)
#'
#' ov <- calc_volume_overlap(a, b)
#' names(ov)
#' terra::global(ov[[c("volume_a", "volume_b", "volume_overlap")]], "sum",
#'               na.rm = TRUE)
#'
#' # Mixed input: the envelope is discretized onto the voxel's levels.
#' bv <- envelope_to_voxel(b, depths = c(0, 50, 100, 150, 200))
#' terra::global(calc_volume_overlap(a, bv)[["volume_overlap"]], "sum",
#'               na.rm = TRUE)
#' @export
setGeneric("calc_volume_overlap",
           function(x, y, ...) standardGeneric("calc_volume_overlap"))
