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
    "Depth along the route" = "depth_step_m",
    "Time along the route" = c("time_step", "time_unit"),
    "Horizontal buffer" = "horizontal_radius_km",
    "Depth buffer" = "depth_radius_m",
    "Time buffer" = c("time_before", "time_after"),
    "Statistics inside each buffer" = c("statistic", "weighting"),
    "Track summary" = character(),
    "Depth selection" = "method",
    "Observation and output columns" = c("lon_col", "lat_col", "depth_col",
                                         "date_col", "output_col",
                                         "output_prefix", "id_col"),
    "Source dimensions" = c("lon_dim", "lat_dim", "depth_dim", "time_dim"),
    "Temporal selection" = c("time_match", "max_time_diff"),
    "Route limits" = "max_samples"
  )
  if (name == "extract_to_track") {
    groups[["Statistics inside each buffer"]] <- character()
    groups[["Track summary"]] <- c("statistic", "weighting")
  }
  if (name == "extract_to_buffer") {
    groups[["Time buffer"]] <- c("time_before", "time_after", "time_unit")
    groups[["Time along the route"]] <- character()
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
      value <- paste(deparse(defaults[[argument]], width.cutoff = 500L),
                     collapse = " ")
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
