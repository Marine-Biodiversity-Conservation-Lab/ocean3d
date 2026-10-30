# Figure helpers for the buffer article, not exported package functions ---------

plot_buffer_footprint <- function(example) {
  old <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  temperatures <- vapply(example$cells, function(cell) {
    selected <- example$samples[example$samples$cell == cell$id, ]
    weighted.mean(selected$value, selected$weight)
  }, numeric(1))
  colours <- hcl.colors(64, "Viridis")
  index <- cut(temperatures, breaks = 64, labels = FALSE)
  plot(NA, xlim = c(-13, 13), ylim = c(-13, 13), asp = 1,
    xlab = "East of observation (km)", ylab = "North of observation (km)",
    main = "Our horizontal buffer")
  for (i in seq_along(example$cells)) {
    cell <- example$cells[[i]]
    polygon(cell$full, border = "grey75")
    polygon(cell$polygon, col = colours[index[i]], border = "grey45")
  }
  theta <- seq(0, 2 * pi, length.out = 181)
  lines(10 * cos(theta), 10 * sin(theta), lwd = 2)
  points(0, 0, pch = 4, col = "red", cex = 1.5, lwd = 2)
  legend("topleft", sprintf("%.2f degrees C", range(temperatures)),
    fill = colours[c(1, 64)], bty = "n", cex = 0.65)

  depths <- example$depths
  layer_temperature <- vapply(depths$index, function(z) {
    selected <- example$samples[example$samples$depth_index == z, ]
    weighted.mean(selected$value, selected$weight)
  }, numeric(1))
  xr <- range(layer_temperature) + c(-0.15, 0.15)
  plot(NA, xlim = xr, ylim = c(78, 46),
    xlab = "Mean potential temperature (degrees C)", ylab = "Depth (m)",
    main = "Our vertical buffer")
  rect(xr[1], 50, xr[2], 70, col = "grey95", border = NA)
  abline(h = c(50, 70), lty = 2, col = "grey50")
  segments(layer_temperature, depths$lower, layer_temperature, depths$upper,
    lwd = 5, col = "steelblue")
  lines(layer_temperature, depths$centre, type = "b", pch = 16)
  points(mean(xr), 60, pch = 4, col = "red", cex = 1.5, lwd = 2)
}


plot_buffer_weights <- function(example) {
  old <- par(mfrow = c(1, 2), mar = c(4, 5, 3, 1), oma = c(2, 0, 0, 0))
  on.exit(par(old))
  colours <- c("steelblue", "darkorange")
  d <- example$depths
  barplot(rbind(100 * rev(d$fraction), rep(100 / nrow(d), nrow(d))),
    beside = TRUE, horiz = TRUE, col = colours, border = NA,
    names.arg = sprintf("%.2f m", rev(d$centre)), las = 1,
    xlim = c(0, 55), xlab = "Share of vertical support (%)", main = "Depth layers")
  t <- example$times
  barplot(rbind(100 * rev(t$fraction), rep(100 / nrow(t), nrow(t))),
    beside = TRUE, horiz = TRUE, col = colours, border = NA,
    names.arg = substr(rev(t$label), 6, 10), las = 1,
    xlim = c(0, 55), xlab = "Share of temporal support (%)", main = "Timestamps (UTC)")
  legend("bottomright", c("Weighted", "Equal"), fill = colours, bty = "n", cex = 0.85)
}


buffer_explorer <- function(example, script) {
  rows <- function(x) lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
  data <- list(cells = example$cells, depths = rows(example$depths),
    times = rows(example$times), samples = rows(example$samples),
    weighted = example$weighted, equal = example$equal)
  p <- plotly::plot_ly(x = 0, y = 0, z = 60, type = "scatter3d", mode = "markers", height = 560,
    marker = list(size = 6, color = "black"), hoverinfo = "skip")
  p <- plotly::layout(p, margin = list(l = 0, r = 0, t = 10, b = 0))
  p <- plotly::config(p, displaylogo = FALSE, responsive = TRUE)
  htmlwidgets::onRender(p, paste(readLines(script, warn = FALSE), collapse = "\n"), data = data)
}
