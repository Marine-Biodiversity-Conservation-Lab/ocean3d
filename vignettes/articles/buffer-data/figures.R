# Drawing helpers for the track tutorial ---------------------------------------

plot_track_map <- function(samples, vertices) {
  values <- samples$temperature
  palette <- grDevices::hcl.colors(64, "Viridis")
  colour <- if (diff(range(values, na.rm = TRUE)) == 0) {
    rep(palette[32], length(values))
  } else {
    palette[cut(values, 64, labels = FALSE)]
  }
  plot(samples$lon, samples$lat, type = "l", col = "grey60",
       asp = 1 / cos(mean(samples$lat) * pi / 180),
       xlab = "Longitude (degrees east)", ylab = "Latitude (degrees north)")
  points(samples$lon, samples$lat, pch = 16, col = colour, cex = 0.8)
  points(vertices$lon, vertices$lat, pch = 21, bg = "white", cex = 1.4)
  legend("topleft", legend = c("Recorded vertex", "Sampled route"),
         pch = c(21, 16), pt.bg = c("white", NA),
         col = c("black", palette[32]), bty = "n")
}


track_explorer <- function(samples) {
  hover <- paste0(
    "Sample ", samples$sample_id,
    "<br>", ifelse(
      samples$is_interpolated, "Assumed position", "Recorded vertex"
    ),
    "<br>Distance: ", round(samples$distance_km, 2), " km",
    "<br>Requested depth: ", round(samples$depth, 1), " m",
    "<br>UTC: ", format(samples$date, "%Y-%m-%d %H:%M", tz = "UTC"),
    "<br>Selected temperature: ", round(samples$temperature, 2), " degrees C"
  )
  figure <- plotly::plot_ly(
    x = samples$lon, y = samples$lat, z = samples$depth,
    type = "scatter3d", mode = "lines+markers",
    color = samples$temperature, colors = grDevices::hcl.colors(64, "Viridis"),
    text = hover, hoverinfo = "text", marker = list(size = 4), height = 520
  )
  figure <- plotly::layout(figure, scene = list(
    aspectmode = "manual", aspectratio = list(x = 1.4, y = 1, z = 0.9),
    camera = list(eye = list(x = 1.6, y = 1.8, z = 0.8)),
    xaxis = list(title = "Longitude"), yaxis = list(title = "Latitude"),
    zaxis = list(title = "Assumed depth (m)", autorange = "reversed")
  ), margin = list(l = 25, r = 15, b = 20, t = 15))
  figure <- plotly::colorbar(figure, title = "Temperature (C)", thickness = 12)
  figure <- plotly::config(figure, displaylogo = FALSE, responsive = TRUE)
  htmlwidgets::onRender(figure, paste0(
    "function(el,x){el.style.width='100%';",
    "function resize(){window.Plotly.Plots.resize(el);}",
    "window.addEventListener('resize',resize);resize();}"
  ))
}

# Display complete function arguments in related groups for the tutorials.
article_signature <- function(fun, name) {
  groups <- list(
    "Data and sources" = c("data", "nc", "var", "track_col"),
    "Direct observation coordinates" = c("lon", "lat", "depth", "date"),
    "Route sampling" = c("sampling", "spacing_km"),
    "Depth along the route" = c("interpolate_depth", "depth_step_m", "z_as_depth"),
    "Time along the route" = c("interpolate_time", "time_step", "time_unit",
                              "time_step_seconds"),
    "Horizontal buffer" = c("horizontal_radius_km", "lon_cell_radius", "lat_cell_radius"),
    "Depth buffer" = c("depth_radius_m", "depth_layer_radius"),
    "Time buffer" = c("time_before", "time_after", "buffer_time_unit"),
    "Statistics inside each buffer" = c("statistic", "buffer_weighting"),
    "Track summary" = c("statistics", "weighting", "track_min_valid_fraction"),
    "Valid coverage and bounds" = c("min_valid_fraction", "infer_bounds"),
    "Similar-value selection" = c("target_col", "interpolate_target",
                                 "max_target_difference", "target_ties"),
    "Depth selection" = "method",
    "Observation and output columns" = c("lon_col", "lat_col", "depth_col",
                                        "date_col", "output_col", "output_prefix", "id_col"),
    "Source columns and dimensions" = c("file_col", "lon_dim", "lat_dim",
                                       "depth_dim", "time_dim"),
    "Temporal selection" = c("time_match", "time_resolution",
                             "custom_time_interval", "max_time_diff"),
    "Spatial selection" = c("spatial_match", "max_distance_km"),
    "Distances and route limits" = c("distance_model", "max_samples",
                                     "max_segment_km", "max_gap_hours"),
    "Marine mask" = c("land_mask", "mask_type", "marine_values", "mask_domain"),
    "Mask checks" = c("land_policy", "mask_step_km"),
    "Diagnostics and messages" = c("return_match_info", "time_report",
                                   "spatial_report", "depth_report", "verbose", "quiet"),
    "Additional selection options" = "..."
  )
  if (name == "extract_to_buffer") {
    groups[["Time buffer"]] <- c("time_before", "time_after", "time_unit")
    groups[["Time along the route"]] <- character()
    groups[["Statistics inside each buffer"]] <- c("statistic", "weighting")
    groups[["Track summary"]] <- character()
  }
  defaults <- formals(fun)
  groups <- lapply(groups, intersect, y = names(defaults))
  groups <- groups[lengths(groups) > 0L]
  displayed <- unlist(groups, use.names = FALSE)
  stopifnot(!anyDuplicated(displayed), setequal(displayed, names(defaults)))
  lines <- paste0(name, "(")
  for (title in names(groups)) {
    lines <- c(lines, paste0("  # ", title))
    line <- " "
    for (argument in groups[[title]]) {
      value <- paste(deparse(defaults[[argument]], width.cutoff = 500L), collapse = " ")
      entry <- if (value == "" || argument == "...") argument else
        paste(argument, "=", value)
      if (nchar(line) > 1L && nchar(line) + nchar(entry) + 2L > 80L) {
        lines <- c(lines, line)
        line <- " "
      }
      line <- paste0(line, " ", entry, ",")
    }
    lines <- c(lines, line, "")
  }
  lines[length(lines) - 1L] <- sub(",$", "", lines[length(lines) - 1L])
  lines <- c(lines[-length(lines)], ")")
  cat("```r\n", paste(lines, collapse = "\n"), "\n```\n", sep = "")
  invisible(lines)
}


# Explain the selection options forwarded by both track interfaces.
article_selection_notes <- function(arguments) {
  notes <- c(
    file_col = "Source-table path column; default file.",
    lon_dim = "NULL detects the longitude dimension; supply a name to override.",
    lat_dim = "NULL detects the latitude dimension; supply a name to override.",
    depth_dim = "NULL detects the depth dimension; supply a name to override.",
    time_dim = "NULL detects the time dimension; supply a name to override.",
    time_match = "nearest selects the nearest eligible timestamp/label (default); exact requires the requested timestamp/label.",
    max_time_diff = "NULL has no extra temporal-distance limit; a nonnegative number sets a limit in days.",
    custom_time_interval = "Positive difftime warning threshold for custom resolution, which retains full timestamps. NULL by default; does not reject candidates.",
    spatial_match = "within requires a containing cell (default); nearest allows a nearby cell with an explicit max_distance_km.",
    max_distance_km = "Nonnegative distance to cell centre, required for nearest; NULL for within.",
    time_report = "TRUE reports displaced or missing temporal results (default); FALSE hides the report.",
    spatial_report = "NULL reports exceptional spatial selections (default); TRUE adds routine detail; FALSE hides the report.",
    depth_report = "NULL reports exceptional depth selections (default); TRUE adds routine detail; FALSE hides the report.",
    verbose = "FALSE suppresses progress detail (default); TRUE shows it."
  )
  stopifnot(all(arguments %in% names(notes)))
  unname(notes[arguments])
}
