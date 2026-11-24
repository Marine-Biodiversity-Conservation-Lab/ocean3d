# R/AllGenerics.R
#
# Generics for the operations that are meaningful on either 3D representation.
# The methods live in R/volume.R and R/intersect.R; only the dispatch contract
# is declared here.

#' @importFrom methods setGeneric setMethod
NULL

#' 3D volume of a domain
#'
#' Measures the occupied volume of a [SpatEnvelope-class], a
#' [SpatVoxel-class] or [SpatVector3D-class] polygons, as a total, or cell by
#' cell for envelopes and voxels. For a [SpatProxy-class], it gives the volume
#' of each time step (see the "Memory" and "Overlap through time" sections).
#' Envelopes and voxels describe a 3D
#' domain over a 2D grid, but they store depth in dual roles, so each
#' computes the vertical extent differently:
#'
#' \describe{
#'   \item{[SpatEnvelope-class]}{Depth is the cell value, and the interval is
#'     solid by construction. Volume is
#'     `sum(cell_area * (depth_max - depth_min))` over present cells.}
#'   \item{[SpatVoxel-class]}{Depth is the layer index. Each occupied voxel
#'     contributes its cell area times the thickness of the slab its depth
#'     level stands for, so interior gaps are excluded rather than filled.}
#'   \item{[SpatVector3D-class] polygons}{Measured from the outlines, not a
#'     grid: each polygon's area on the WGS84 ellipsoid times
#'     `depth_max - depth_min`. Polygon edges are geodesics, so a polygon
#'     drawn along a grid's cell edges measures slightly differently from
#'     its envelope on that grid, whose cells are bounded by parallels.
#'     The time interval doesn't change the volume. The total adds up the
#'     polygons, so where they overlap the shared space is counted once for
#'     each; the overlap between two sets is `volume(intersect_3d(x, y))`.
#'     `per_cell = TRUE` is an error, since polygons have no cells; put them
#'     on a grid with [vector3d_to_envelope()] first. Points and lines have
#'     no volume, so they are an error; [buffer_3d()] turns them into
#'     polygons.}
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
#' @param x A [SpatEnvelope-class], a [SpatVoxel-class],
#'   [SpatVector3D-class] polygons, or a [SpatProxy-class], which is
#'   measured step by step. A bare SpatRaster is
#'   rejected; build one with [as_envelope()], [vect_to_envelope()] or
#'   [as_voxel()] first. A netCDF file isn't accepted: describe it with
#'   [as_proxy()] and pass the proxy, or read one time step with
#'   [as_voxel()].
#' @param per_cell `FALSE` (default) for the total volume, `TRUE` for the
#'   volume of each cell. Not for [SpatVector3D-class] polygons.
#' @param per_level SpatVoxel only, with `per_cell = TRUE`: `TRUE` for the
#'   volume of each cell *at each depth level*, the weights
#'   [summarise_3d()] uses. Envelopes and polygons have no levels, so it is
#'   an error for them. It is an error for a [SpatProxy-class] too: read a
#'   time step with [as_voxel()] and measure that.
#' @param ... Arguments for the [SpatVoxel-class] method, which are an error
#'   for an envelope or polygons.
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
#'   depth interval, or no occupied level). Its sum is the total. With
#'   `per_level = TRUE` as well, a [SpatVoxel-class] on the depth levels and
#'   edges of `x`, layers `volume_depth=<value>`: the volume of each cell at
#'   each level in km³, `NA` where the level is unoccupied. Summing its
#'   layers gives the per-cell raster.
#'
#'   For a [SpatProxy-class], one result per time step: with
#'   `per_cell = FALSE` a data frame with columns `time` (the step's
#'   date-time in UTC, or `NA` for a one-step climatology or a dataset
#'   without a time axis) and `volume_km3`; with `per_cell = TRUE` a
#'   `SpatRaster` with one
#'   layer per step, named `volume`, with the steps' dates in
#'   [terra::time()].
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
#' # Per cell and level: each occupied level's share of the column.
#' terra::values(volume(v, per_cell = TRUE, per_level = TRUE))
#'
#' # Overlap with a second domain at 100-300 m: intersect, then measure.
#' b <- as_envelope(fp, depth_min = 100, depth_max = 300)
#' ov <- intersect_3d(e, b)
#' volume(ov)
#' terra::values(volume(ov, per_cell = TRUE))
#'
#' # Polygons: a 10 km circle from 10 to 60 m deep, about
#' # pi * 10^2 * 0.05 = 15.7 km³.
#' detection <- as_vector3d(data.frame(lon = 150, lat = -20, depth = 35))
#' volume(buffer_3d(detection, radius_km = 10, depth_window = 25))
#'
#' # Through time: the volume of water warmer than 2 °C in each step.
#' p <- as_proxy(system.file("extdata", "example_3d.nc", package = "ocean3d"))
#' volume(p, fun = function(v) v > 2)
#' @export
setGeneric("volume", function(x, per_cell = FALSE, per_level = FALSE, ...) {
  # Checked once here, before a method computes anything.
  if (!rlang::is_bool(per_cell)) {
    rlang::abort("`per_cell` must be TRUE or FALSE.",
                 class = "ocean3d_error_per_cell", call = NULL)
  }
  .check_per_level(per_level, per_cell)
  standardGeneric("volume")
}, signature = "x")

#' The 3D space two objects share
#'
#' `intersect_3d()` returns the part of space that is inside both `x` and
#' `y`. It answers the question "what do these two share?". To ask only
#' "do they share anything?", use [intersects_3d()]. To keep one object's
#' values where the other is present, use [mask()].
#'
#' At least one of `x` and `y` must be a 3D object: a [SpatEnvelope-class], a
#' [SpatVoxel-class] or [SpatVector3D-class] polygons or points. The other
#' can be a 3D object too, or, with an envelope or voxel, a 2D object. A 2D
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
#'   \item{An envelope or voxel and a `SpatRaster`}{The raster is a
#'     footprint. Its non-`NA` cells are inside; its values are not read. It
#'     must have one layer and be on the same grid.}
#'   \item{An envelope or voxel and 2D polygons}{A `SpatVector`, `sf` or
#'     `sfc` object. The polygons are rasterised onto the 3D object's grid. A
#'     cell is inside when a polygon touches it. Polygons in another CRS are
#'     projected first.}
#'   \item{An envelope and [SpatVector3D-class] polygons}{The polygons are
#'     first put on the envelope's grid with [vector3d_to_envelope()], then
#'     compared as an envelope. Each cell a polygon touches gets its depth
#'     range; where polygons overlap, a cell gets the shallowest top and the
#'     deepest bottom, so a depth gap between two polygons over one cell is
#'     filled.}
#'   \item{A voxel and [SpatVector3D-class] polygons}{The polygons are first
#'     put on the voxel's grid and depth levels with [vector3d_to_voxel()],
#'     one level at a time, then compared as a voxel. A depth gap between two
#'     polygons over one cell stays empty. When the voxel has no depth edges,
#'     you set each level's slab with `bounds`, as in [envelope_to_voxel()].}
#'   \item{An envelope or voxel and [SpatVector3D-class] polygons, either
#'     way}{The polygons' time intervals are dropped, with a message.
#'     Polygons that miss the grid, or a set with no polygons, give an empty
#'     result.}
#'   \item{A voxel and [SpatVector3D-class] points}{A point covers one cell
#'     and one depth level: the cell that holds it and the level whose depth
#'     is nearest the point's depth (a tie goes to the shallower level, as in
#'     [extract_to_point()]). The level is chosen by depth alone, so a point
#'     whose nearest level has no value in the voxel shares nothing there.
#'     When the voxel's depth edges give uneven slabs, the slab that contains
#'     a point can belong to another level; the nearest level is used. The
#'     result is a presence voxel on the voxel's levels, and `bounds` doesn't
#'     apply. Points off the grid share nothing. Points' times aren't
#'     compared, since a voxel shows one time step.}
#'   \item{An envelope and [SpatVector3D-class] points}{In each cell where a
#'     point's depth lies within the envelope's interval (both ends
#'     included), the result holds an interval from the shallowest such point
#'     to the deepest. A lone point gives a zero-thickness interval, so the
#'     result's [volume()] is 0: points have no volume.}
#'   \item{[SpatVector3D-class] points and polygons}{A point is inside a
#'     polygon when [sf::st_intersects()] finds its position in the outline
#'     (on WGS84 coordinates sf uses the sphere by default, where a point
#'     exactly on an edge may not count), its depth is within `depth_min` to
#'     `depth_max` and its time within the polygon's time interval. Depth and
#'     time include both ends (with dates, the interval covers all of its
#'     end day). A point or polygon without a time matches any
#'     time; one side with dates and the other with hours is an error. Each
#'     matching pair gives one point, numbered, with columns `feature_id_x`
#'     and `feature_id_y` naming the pair and the inputs' other columns
#'     dropped. Two sets of points share no space, so they are an error.
#'     Lines aren't supported; [buffer_3d()] turns them into polygons.}
#'   \item{Two sets of [SpatVector3D-class] polygons}{Each pair of polygons,
#'     one from `x` and one from `y`, whose outlines, depth ranges and times
#'     overlap gives one polygon: the shared outline (as
#'     [sf::st_intersection()]), the shared depth range, and the shared time
#'     interval. Depth ranges that only touch share no water. Time intervals
#'     include both ends, so with dates, two that share a day overlap on it,
#'     and with hours, two that meet at one instant overlap at that instant
#'     (so a polygon at a single time can meet another).
#'     A polygon without a time overlaps any time and passes on the other's
#'     interval. One set can't have dates while the other has hours. The
#'     result's feature IDs are numbers; columns `feature_id_x` and
#'     `feature_id_y` say which polygons each came from, and the inputs'
#'     other columns are dropped. Pairs that share nothing give no polygon,
#'     so the result can have none. To measure the overlap, use
#'     `volume(intersect_3d(x, y))`.}
#' }
#'
#' The order of `x` and `y` does not matter.
#'
#' The result is a domain, not a field. An envelope result carries the shared
#' depth interval. A voxel result carries presence: `1` where a level is
#' shared, `NA` elsewhere, in layers named `presence_depth=<value>`. The cell
#' values of a voxel input are not carried over. To keep them, use [mask()].
#'
#' @param x,y The two objects: two of [SpatEnvelope-class],
#'   [SpatVoxel-class] and [SpatVector3D-class] polygons or points, in any
#'   pairing except two sets of points, or an envelope or voxel with a 2D
#'   object (a single-layer `SpatRaster` footprint or 2D polygons:
#'   `SpatVector`, `sf`, `sfc`). Rasters must share one grid (CRS, extent,
#'   resolution). netCDF files aren't accepted; load them with terra and
#'   wrap them with [as_voxel()] first.
#' @param ... Arguments for the voxel methods. They are an error when both
#'   inputs are envelopes, or when one is points and the other an envelope
#'   or a voxel.
#' @param bounds Only when one input is a voxel and the other an envelope or
#'   [SpatVector3D-class] polygons. How the envelope or polygons are placed
#'   on the voxel's depth levels: `"midpoint"` (default) or `"top"`. See
#'   [envelope_to_voxel()]. Not used when the voxel carries its own depth
#'   edges (see [SpatVoxel-class]). To measure the result, pass the same
#'   `bounds` to [volume()].
#'
#' @returns [SpatVector3D-class] polygons when both inputs are polygons, and
#'   [SpatVector3D-class] points when one is points and the other polygons. A
#'   [SpatEnvelope-class] when both inputs are envelopes, or when one
#'   is an envelope and the other is 2D or a [SpatVector3D-class].
#'   Otherwise a [SpatVoxel-class] of
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
#'
#' # With SpatVector3D polygons, the depth ranges are compared too: a polygon
#' # from 80 to 300 m over the first two cells shares 80-100 m with `a`.
#' deep <- as_vector3d(
#'   sf::st_sf(geometry = sf::st_as_sfc(
#'     "POLYGON ((0 0, 2000 0, 2000 1000, 0 1000, 0 0))",
#'     crs = terra::crs(fp)
#'   )),
#'   depth_min = 80, depth_max = 300
#' )
#' terra::values(intersect_3d(a, deep))
#'
#' # A SpatVector3D point covers the cell that holds it and, on a voxel, the
#' # level nearest its depth: here the 50 m level of `bv` in the first cell
#' # (the point is projected onto the grid's CRS).
#' detection <- as_vector3d(data.frame(lon = 0.0045, lat = 0.0045, depth = 60))
#' terra::values(intersect_3d(detection, bv))
#'
#' # Two sets of SpatVector3D polygons give the polygons they share.
#' square <- function(x0) {
#'   sf::st_polygon(list(rbind(c(x0, 0), c(x0 + 2, 0), c(x0 + 2, 2),
#'                             c(x0, 2), c(x0, 0))))
#' }
#' mpa <- as_vector3d(
#'   sf::st_sf(geometry = sf::st_sfc(square(150), crs = 4326)),
#'   depth_min = 0, depth_max = 200
#' )
#' habitat <- as_vector3d(
#'   sf::st_sf(geometry = sf::st_sfc(square(151), crs = 4326)),
#'   depth_min = 100, depth_max = 500
#' )
#' shared <- intersect_3d(mpa, habitat)
#' as.data.frame(shared)[c("feature_id_x", "feature_id_y", "depth_min",
#'                         "depth_max")]
#' volume(shared)
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
#' [SpatVector3D-class] polygons and points have depths, so they are compared
#' as [intersect_3d()] compares them, after being put on the 3D object's
#' grid. A cell is `TRUE` where the voxel has a value at the level nearest a
#' point's depth, or where a point's depth lies within the envelope's
#' interval.
#'
#' Two [SpatVector3D-class] sets (polygons with polygons, or points with
#' polygons) have no grid to answer on, so the answer is a list instead, as
#' [sf::st_intersects()] gives in 2D: for each feature of `x`, the features
#' of `y` it shares 3D space with. A pair is listed exactly when
#' [intersect_3d()] gives it a polygon or point, so outlines that only share
#' an edge, depth ranges that only touch and times that don't overlap aren't
#' listed. Swapping `x` and `y` gives the same pairs listed the other way
#' round.
#'
#' The order of `x` and `y` does not matter.
#'
#' @inheritParams intersect_3d
#'
#' @returns A single-layer boolean `SpatRaster` named `intersects`, on the
#'   grid of the 3D input: `TRUE`, `FALSE` or `NA` per cell as above. For two
#'   [SpatVector3D-class] sets, a list with one element per feature of `x`,
#'   named by its feature ID: the feature IDs of the features of `y` it
#'   shares space with, in the order of `y`, or `character(0)` for none.
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
