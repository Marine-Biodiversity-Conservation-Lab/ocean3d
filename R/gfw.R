# ---------------------------------------------------------------------------
# Global Fishing Watch (GFW) apparent fishing effort utilities.
#
# Ingest is delegated to the `gfwr` package, which queries the GFW 4Wings
# API and returns a long-format tibble of apparent fishing hours per cell,
# already aggregated server-side by the chosen `group_by`. The helpers
# below cover the gap between that tibble and the package's depth-aware
# pipeline:
#
#   - `gfw_effort_to_raster()` rasterises the gfwr tibble onto the study
#     grid as a multi-layer SpatRaster (one layer per group level).
#
# Extending that 2D product into 3D is not a GFW-specific operation: a
# gear's operating window is a depth envelope, and turning an envelope
# into depth layers is `as_envelope()` + `envelope_to_voxel()`. The
# gear → depth-window priors that feed them are assumptions rather than
# data — operating depths vary by region, fleet, and time — so they live
# in the `gfw-fishing-effort-3d` vignette alongside the analysis that owns
# them, with placeholder values users are expected to replace, rather than
# in an exported function here.
# ---------------------------------------------------------------------------

#' Rasterise a GFW effort tibble onto a target grid
#'
#' Turn the long-format apparent-fishing-hours tibble returned by
#' `gfwr::gfw_ais_fishing_hours()` (formerly `gfwr::get_raster()`) into a
#' multi-layer SpatRaster on the package's canonical study grid. Each
#' level of `layer_by` becomes its own layer, named
#' `effort_<level>` (e.g. `effort_drifting_longlines`).
#'
#' The input is expected to carry a cell centroid (`Lat`, `Lon`), a value
#' column (default `"Apparent Fishing Hours"`), and one categorical column
#' matching the API's `group_by` — for example `geartype` or `flag` (note
#' lower-case; this is what `gfwr` actually returns). Records that fall
#' into the same target cell × layer level are aggregated with `fun`.
#'
#' @param effort Data frame. Output of `gfwr::gfw_ais_fishing_hours()` (a
#'   long-format tibble with at minimum `Lat`, `Lon`, a value, and a
#'   grouping column).
#' @param grid SpatRaster. Target grid (extent, resolution, CRS) — typically
#'   the same grid used for species ranges and WOA extraction. If `grid = NULL`, 
#'   assume target grid from input effort data frame. 
#' @param layer_by Character. Column in `effort` whose levels become layers.
#'   `NULL` produces a single-layer total-effort raster. Default
#'   `"geartype"`.
#' @param value Character. Column in `effort` to aggregate. Default
#'   `"Apparent Fishing Hours"`.
#' @param fun Character or function. Aggregation applied to records that
#'   fall into the same cell × layer level. Default `"sum"`.
#'
#' @returns A SpatRaster with one layer per `layer_by` level (or one layer
#'   if `layer_by = NULL`). Layer names follow `effort_<level>`.
#' 
#' @examples
#' \dontrun{
#' jpn_eez_id <- gfwr::gfw_region_id(region = "JPN", region_source = "EEZ")$id
#' effort <- gfwr::gfw_ais_fishing_hours(
#'   spatial_resolution  = "LOW",        # 0.01 deg
#'   temporal_resolution = "YEARLY",      # one bucket per year in [start, end]
#'   start_date          = "2022-01-01",
#'   end_date            = "2022-12-31",
#'   region              = jpn_eez_id,
#'   region_source       = "EEZ",
#'   group_by            = "GEARTYPE"
#' )
#'
#' # Returns data.frame 
#' head(effort)
#'
#' # Create multi-layer SpatRaster from data.frame, one for each geartype
#' effort_by_gear <- gfw_effort_to_raster(
#'   effort   = effort,
#'   layer_by = "geartype",
#'   value    = "Apparent Fishing Hours",
#'   fun      = "sum"
#' )
#' 
#' # Plot drifitng longlines effort
#' terra::plot(effort_by_gear$effort_drifting_longlines)
#' }
#' @export
gfw_effort_to_raster <- function(effort,
                                 grid = NULL,
                                 layer_by = "geartype",
                                 value = "Apparent Fishing Hours",
                                 fun = "sum") {
  effort <- as.data.frame(effort)
  required <- c("Lat", "Lon", value, layer_by)
  missing_cols <- setdiff(required, names(effort))
  if (length(missing_cols) > 0) {
    stop(
      "effort is missing required columns: ",
      toString(missing_cols),
      call. = FALSE
    )
  }

  pts_df <- data.frame(
    x = effort$Lon,
    y = effort$Lat,
    v = effort[[value]]
  )
  if (!is.null(layer_by)) {
    pts_df$grp <- effort[[layer_by]]
  }

  pts <- terra::vect(
    pts_df,
    geom = c("x", "y"),
    crs = "EPSG:4326",
    keepgeom = FALSE
  )

  # rasterize the SpatVector without input grid. 
  # gfwr-derived SpatVector has a regular point 
  # spacing based on the original grid resolution
  if (is.null(grid)) {
    # determine grid resolution from the input effort
    # ascending order unique Lat column values 
    lat_vals <- sort(effort$Lat) %>% unique()
    # determine interval between Lat values
    # round to prevent floating point errors
    lat_intervals <- vapply(
      seq_len(length(lat_vals) - 1L),
      function(i) round((lat_vals[[i + 1L]] - lat_vals[[i]]), 10), 
      numeric(1)
    )
    if(length(unique(lat_intervals)) > 1) {
      stop("Unable to assume grid resolution from effort data frame. Please provide grid.")
    } else {
      grid_res <- unique(lat_intervals)
    }

    grid <- rast(
      # gfwr effort coordinates are for the cell mid-point
      # create grid offset by half of grid_res value
      xmin = (min(effort$Lon) - grid_res/2), 
      ymin = (min(effort$Lat) - grid_res/2), 
      xmax = (max(effort$Lon) + grid_res/2), 
      ymax = (max(effort$Lat) + grid_res/2), 
      res = grid_res, crs = "EPSG:4326"
    )
  }

  if (!terra::same.crs(pts, grid)) {
    pts <- terra::project(pts, terra::crs(grid))
  }

  if (is.null(layer_by)) {
    out <- terra::rasterize(pts, grid, field = "v", fun = fun, background = NA)
    names(out) <- "effort"
  } else {
    out <- terra::rasterize(
      pts, grid,
      field = "v", fun = fun, by = "grp", background = NA
    )
    names(out) <- paste0("effort_", names(out))
  }
  out
}
