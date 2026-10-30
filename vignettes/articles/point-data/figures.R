# Drawing helpers for the point tutorial, not exported package functions --------

plot_point_context <- function(example) {
  old <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  cells <- example$cells
  colours <- hcl.colors(64, "Viridis")
  temperatures <- vapply(cells, function(cell) cell$value, numeric(1))
  index <- cut(temperatures, 64, labels = FALSE)
  coordinates <- do.call(rbind, lapply(cells, function(cell) cell$polygon))
  plot(NA, xlim = range(coordinates[, 1]), ylim = range(coordinates[, 2]),
    asp = 1, xlab = "East of observation (km)", ylab = "North of observation (km)",
    main = "Selecting a cell")
  for (i in seq_along(cells)) polygon(cells[[i]]$polygon,
    col = colours[index[i]], border = "grey70")
  selected <- cells[[example$selected_cell]]
  polygon(selected$polygon, border = "black", lwd = 2)
  points(0, 0, pch = 4, col = "red", cex = 1.5, lwd = 2)
  points(selected$centre[1], selected$centre[2], pch = 16)
  legend("topleft", sprintf("%.2f degrees C", range(temperatures)),
    fill = colours[c(1, 64)], bty = "n", cex = 0.65)

  profile <- example$profiles[[2]]$values
  record <- Filter(function(r) r$date == "2026-09-20" && r$target == 60, example$vertical)[[1]]
  plot(profile, example$depths, type = "b", pch = 16, cex = 0.5,
    ylim = rev(range(c(0, example$depths))),
    xlab = "Potential temperature (degrees C)", ylab = "Depth (m)",
    main = "Selecting a layer")
  abline(h = 60, lty = 2, col = "red")
  points(profile[record$nearest], example$depths[record$nearest],
    pch = 21, cex = 1.4, bg = "darkorange")
}


point_explorer <- function(example, view, script) {
  if (view == "horizontal") {
    data <- example[c("cells", "locations", "horizontal")]
    data$locations <- lapply(seq_len(nrow(example$locations)), function(i) {
      as.list(example$locations[i, , drop = FALSE])
    })
  } else if (view == "vertical") {
    data <- example[c("profiles", "depths", "depth_bounds", "vertical")]
    data$cell <- example$cells[[example$selected_cell]]
  } else {
    data <- example[c("temporal", "profiles")]
  }
  data$view <- view
  p <- plotly::plot_ly(x = 0, y = 0, type = "scatter", mode = "markers", height = 500)
  p <- plotly::config(p, displaylogo = FALSE, responsive = TRUE)
  htmlwidgets::onRender(p, paste(readLines(script, warn = FALSE), collapse = "\n"), data = data)
}
