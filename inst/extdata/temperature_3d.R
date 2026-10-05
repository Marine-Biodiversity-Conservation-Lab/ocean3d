# Optional 3D visualisation of temperature layers
#
# Educational helper for the ocean3d workshop.
# Uses base R graphics and terra.
#
# Input:
#   x: a SpatVoxel containing environmental data.

plot_temperature_3d <- function(x) {

  # Extract depth levels and geographical coordinates.
  depth_levels <- ocean3d::depths(x)

  lon <- terra::xFromCol(x, seq_len(terra::ncol(x)))
  lat <- terra::yFromRow(x, terra::nrow(x):1)

  # Define a common temperature colour scale.
  temperature_range <- range(
    terra::minmax(x),
    na.rm = TRUE
  )

  palette <- hcl.colors(100, "YlOrRd", rev = TRUE)

  temperature_to_colour <- function(values) {
    indices <- round(
      1 + 99 * (values - temperature_range[1]) /
        diff(temperature_range)
    )

    palette[pmax(1, pmin(100, indices))]
  }

  # Initialise the 3D projection.
  z_bottom <- -max(depth_levels)

  projection <- persp(
    x = range(lon),
    y = range(lat),
    z = matrix(z_bottom, nrow = 2, ncol = 2),
    zlim = c(z_bottom, 0),
    theta = 35,
    phi = 25,
    expand = 0.7,
    col = NA,
    border = NA,
    xlab = "Longitude",
    ylab = "Latitude",
    zlab = "Depth (m)",
    ticktype = "detailed",
    main = "3D temperature layers"
  )

  # Draw the temperature layers from deepest to shallowest.
  for (i in rev(seq_along(depth_levels))) {

    z <- -depth_levels[i]

    values <- terra::as.matrix(x[[i]], wide = TRUE)

    # Arrange rows in increasing latitude.
    values <- values[nrow(values):1, , drop = FALSE]

    for (row in seq_len(nrow(values) - 1)) {
      for (column in seq_len(ncol(values) - 1)) {

        corners <- trans3d(
          x = lon[c(column, column + 1,
                    column + 1, column)],
          y = lat[c(row, row,
                    row + 1, row + 1)],
          z = rep(z, 4),
          pmat = projection
        )

        polygon(
          corners$x,
          corners$y,
          col = temperature_to_colour(
            values[row, column]
          ),
          border = NA
        )
      }
    }
  }

  invisible(NULL)
}
