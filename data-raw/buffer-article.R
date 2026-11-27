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

  observation <- data.frame(lon = 0, lat = 38, depth = 60,
    date = as.POSIXct("2026-09-20", tz = "UTC"))
  arguments <- list(data = observation, var = "thetao", horizontal_radius_km = 10,
    depth_radius_m = 10, time_before = 1, time_after = 1)
  original <- do.call(ocean3d::extract_to_buffer, c(arguments, list(nc = source_file)))
  cropped <- do.call(ocean3d::extract_to_buffer, c(arguments, list(nc = subset_file)))
  # Check that the crop keeps the article's result.
  stopifnot(isTRUE(all.equal(original$thetao, cropped$thetao, tolerance = 1e-8)),
    original$thetao_n == cropped$thetao_n)
  invisible(cropped)
}
