# Prepare the real-data example for the buffer extraction article ---------------
# Run from the package root, with the original downloaded netCDF:
# source("data-raw/buffer-article.R")
# prepare_buffer_article("path/to/downloaded-temperature.nc")
# No download, credentials or modification of the original file is required.

prepare_buffer_article <- function(source_file, root = ".") {
  stopifnot(file.exists(source_file))
  destination <- file.path(root, "vignettes", "articles", "buffer-data")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  input <- ncdf4::nc_open(source_file)
  on.exit(ncdf4::nc_close(input), add = TRUE)

  coordinates <- lapply(c("longitude", "latitude", "depth", "time"), function(name) {
    input$dim[[name]]$vals
  })
  names(coordinates) <- c("longitude", "latitude", "depth", "time")
  stopifnot(all(vapply(coordinates, length, integer(1)) > 1L))
  limits <- list(c(-0.15, 0.15), c(37.86, 38.14), c(50, 70),
    as.numeric(as.POSIXct(c("2026-09-19", "2026-09-21"), tz = "UTC")))
  # One neighbouring centre on each side preserves the inferred midpoints.
  indices <- Map(function(values, range) {
    inside <- which(values >= range[1] & values <= range[2])
    stopifnot(length(inside) > 0)
    seq.int(max(1L, min(inside) - 1L), min(length(values), max(inside) + 1L))
  }, coordinates, limits)
  dims <- Map(function(name, values, index) {
    ncdf4::ncdim_def(name, input$dim[[name]]$units, values[index])
  }, names(coordinates), coordinates, indices)
  variable <- input$var$thetao
  order <- match(vapply(variable$dim, function(d) d$name, character(1)),
                 names(coordinates))
  values <- ncdf4::ncvar_get(input, "thetao",
    start = vapply(indices[order], min, integer(1)),
    count = lengths(indices[order]), collapse_degen = FALSE)
  values <- aperm(values, match(seq_along(coordinates), order))
  subset_file <- file.path(destination, "temperature-subset.nc")
  output <- ncdf4::nc_create(subset_file, ncdf4::ncvar_def(
    "thetao", "degrees_C", dims, missval = -9999, prec = "double"))
  tryCatch({
    ncdf4::ncvar_put(output, "thetao", values)
    ncdf4::ncatt_put(output, "thetao", "standard_name", "sea_water_potential_temperature")
    ncdf4::ncatt_put(output, "depth", "positive", "down")
    ncdf4::ncatt_put(output, 0, "Conventions", "CF-1.11")
    ncdf4::ncatt_put(output, 0, "title", "Potential temperature subset for the ocean3d buffer tutorial")
    for (name in c("institution", "source", "references", "subset:productId",
                   "subset:datasetId", "subset:date")) {
      attribute <- ncdf4::ncatt_get(input, 0, name)
      if (attribute$hasatt) ncdf4::ncatt_put(output, 0, name, attribute$value)
    }
    ncdf4::ncatt_put(output, 0, "history", paste(
      "Spatial, depth and time subset prepared for an ocean3d tutorial.",
      "Original values and coordinate centres retained; no bounds added."))
    ncdf4::ncatt_put(output, 0, "source_file", basename(source_file))
    ncdf4::ncatt_put(output, 0, "source_md5", unname(tools::md5sum(source_file)))
    ncdf4::ncatt_put(output, 0, "licence",
      "https://marine.copernicus.eu/user-corner/service-commitments-and-licence")
    ncdf4::ncatt_put(output, 0, "credit",
      "E.U. Copernicus Marine Service Information; https://doi.org/10.48670/mds-00375")
  }, finally = ncdf4::nc_close(output))

  midpoint_bounds <- function(x) {
    middle <- (head(x, -1L) + tail(x, -1L)) / 2
    cbind(c(x[1] - (x[2] - x[1]) / 2, middle),
          c(middle, tail(x, 1) + (tail(x, 1) - x[length(x) - 1L]) / 2))
  }
  centres <- Map(function(x, i) x[i], coordinates, indices)
  bounds <- lapply(centres, midpoint_bounds)
  distance_crs <- "+proj=aeqd +lat_0=38 +lon_0=0 +datum=WGS84 +units=m +no_defs"
  area_crs <- "+proj=laea +lat_0=38 +lon_0=0 +datum=WGS84 +units=m +no_defs"
  origin <- sf::st_sfc(sf::st_point(c(0, 0)), crs = distance_crs)
  circle <- sf::st_transform(sf::st_buffer(origin, 10000, nQuadSegs = 180L), area_crs)
  total_area <- as.numeric(sf::st_area(circle))
  rectangle <- function(x, y) {
    # Densify the cell edges before projecting, as for the extraction geometry.
    xx <- seq(x[1], x[2], length.out = max(2L, ceiling(diff(x) / 0.01) + 1L))
    yy <- seq(y[1], y[2], length.out = max(2L, ceiling(diff(y) / 0.01) + 1L))
    ring <- rbind(cbind(xx, y[1]), cbind(x[2], yy), cbind(rev(xx), y[2]),
                  cbind(x[1], rev(yy)))
    pts <- sf::st_transform(sf::st_sfc(lapply(seq_len(nrow(ring)), function(i) {
      sf::st_point(ring[i, ])
    }), crs = 4326), area_crs)
    sf::st_polygon(list(sf::st_coordinates(pts)[, 1:2]))
  }
  cells <- list()
  for (i in seq_along(centres$longitude)) for (j in seq_along(centres$latitude)) {
    full <- sf::st_sfc(rectangle(bounds$longitude[i, ], bounds$latitude[j, ]), crs = area_crs)
    clipped <- suppressWarnings(sf::st_intersection(full, circle))
    if (!length(clipped) || all(sf::st_is_empty(clipped))) next
    area <- as.numeric(sf::st_area(clipped))
    if (area <= 0) next
    cells[[length(cells) + 1L]] <- list(
      id = length(cells) + 1L, x = i, y = j,
      lon = centres$longitude[i], lat = centres$latitude[j],
      area_fraction = area / total_area,
      polygon = unname(sf::st_coordinates(sf::st_simplify(clipped, dTolerance = 5))[, 1:2] / 1000),
      full = unname(sf::st_coordinates(full)[, 1:2] / 1000))
  }
  z <- which(bounds$depth[, 2] > 50 & bounds$depth[, 1] < 70)
  window <- limits[[4]]
  times <- which(bounds$time[, 2] > window[1] & bounds$time[, 1] < window[2])
  depth_table <- data.frame(index = z, centre = centres$depth[z],
    lower = pmax(50, bounds$depth[z, 1]), upper = pmin(70, bounds$depth[z, 2]))
  depth_table$fraction <- (depth_table$upper - depth_table$lower) / 20
  time_table <- data.frame(index = times,
    label = format(as.POSIXct(centres$time[times], origin = "1970-01-01", tz = "UTC"),
                   "%Y-%m-%d %H:%M"),
    lower = pmax(window[1], bounds$time[times, 1]),
    upper = pmin(window[2], bounds$time[times, 2]))
  time_table$fraction <- (time_table$upper - time_table$lower) / diff(window)
  samples <- list()
  for (cell in cells) for (k in seq_len(nrow(depth_table))) {
    for (l in seq_len(nrow(time_table))) {
      samples[[length(samples) + 1L]] <- data.frame(
        cell = cell$id, depth_index = z[k], time_index = times[l],
        value = values[cell$x, cell$y, z[k], times[l]],
        weight = cell$area_fraction * depth_table$fraction[k] * time_table$fraction[l])
    }
  }
  samples <- do.call(rbind, samples)
  observation <- data.frame(lon = 0, lat = 38, depth = 60,
    date = as.POSIXct("2026-09-20", tz = "UTC"))
  arguments <- list(data = observation, var = "thetao", horizontal_radius_km = 10,
    depth_radius_m = 10, time_before = 1, time_after = 1, infer_bounds = TRUE, quiet = TRUE)
  original <- do.call(ocean3d::extract_to_buffer, c(arguments, list(nc = source_file)))
  cropped <- do.call(ocean3d::extract_to_buffer, c(arguments, list(nc = subset_file)))
  # Check the crop and the independently reconstructed figure weights.
  stopifnot(isTRUE(all.equal(original$thetao, cropped$thetao, tolerance = 1e-8)),
    isTRUE(all.equal(sum(samples$value * samples$weight) / sum(samples$weight),
                     original$thetao, tolerance = 1e-7)),
    nrow(samples) == original$thetao_match_info[[1]]$n_valid,
    all(is.finite(samples$value)), abs(sum(samples$weight) - 1) < 1e-7)
  example <- list(cells = cells, depths = depth_table, times = time_table,
    samples = samples, weighted = original$thetao, equal = mean(samples$value),
    observation = observation, source_file = basename(source_file),
    subset_md5 = unname(tools::md5sum(subset_file)))
  saveRDS(example, file.path(destination, "temperature-geometry.rds"), version = 2)
  invisible(example)
}
