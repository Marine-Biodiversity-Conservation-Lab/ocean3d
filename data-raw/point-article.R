# Prepare the real-data example for the point extraction article ----------------
# From the package root:
# pkgload::load_all()
# source("data-raw/point-article.R")
# prepare_point_article("path/to/downloaded-temperature.nc")

prepare_point_article <- function(source_file, root = ".") {
  stopifnot(file.exists(source_file))
  destination <- file.path(root, "vignettes", "articles", "point-data")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  input <- ncdf4::nc_open(source_file)
  on.exit(ncdf4::nc_close(input), add = TRUE)
  axis_names <- c("longitude", "latitude", "depth", "time")
  original <- lapply(axis_names, function(name) input$dim[[name]]$vals)
  names(original) <- axis_names
  crop <- function(x, limits) {
    inside <- which(x >= limits[1] & x <= limits[2])
    seq.int(max(1L, min(inside) - 1L), min(length(x), max(inside) + 1L))
  }
  selected <- list(
    longitude = crop(original$longitude, c(-0.15, 0.15)),
    latitude = crop(original$latitude, c(37.86, 38.14)),
    depth = seq_along(original$depth),
    time = which(original$time %in% as.numeric(as.POSIXct(
      c("2026-09-19", "2026-09-20", "2026-09-21"), tz = "UTC")))
  )
  stopifnot(length(selected$time) == 3L)
  centres <- Map(function(x, i) x[i], original, selected)
  read_variable <- function(name, indices) {
    dimensions <- vapply(input$var[[name]]$dim, function(d) d$name, character(1))
    target <- intersect(axis_names, dimensions)
    values <- ncdf4::ncvar_get(input, name,
      start = vapply(indices[dimensions], min, integer(1)),
      count = lengths(indices[dimensions]), collapse_degen = FALSE)
    aperm(values, match(target, dimensions))
  }
  write_subset <- function(file, time_indices, teaching_layers = FALSE) {
    indices <- selected
    indices$time <- time_indices
    dims <- Map(function(name, values, i) {
      ncdf4::ncdim_def(name, input$dim[[name]]$units, values[i])
    }, axis_names, original, indices)
    names(dims) <- axis_names
    variables <- list(ncdf4::ncvar_def("thetao", "degrees_C", dims,
      missval = -9999, prec = "double"))
    if (teaching_layers) variables <- c(variables, list(
      ncdf4::ncvar_def("bottomT", "degrees_C", dims[c("longitude", "latitude", "time")],
        missval = -9999, prec = "double"),
      ncdf4::ncvar_def("cell_id", "1", dims[c("longitude", "latitude")],
        missval = -9999, prec = "integer")
    ))
    output <- ncdf4::nc_create(file, variables)
    tryCatch({
      ncdf4::ncvar_put(output, "thetao", read_variable("thetao", indices))
      ncdf4::ncatt_put(output, "thetao", "standard_name", "sea_water_potential_temperature")
      if (teaching_layers) {
        ncdf4::ncvar_put(output, "bottomT", read_variable("bottomT", indices))
        ncdf4::ncatt_put(output, "bottomT", "long_name", "Original source variable: bottomT")
        identifiers <- matrix(seq_len(length(indices$longitude) * length(indices$latitude)),
          length(indices$longitude), length(indices$latitude))
        ncdf4::ncvar_put(output, "cell_id", identifiers)
        ncdf4::ncatt_put(output, "cell_id", "long_name",
          "Artificial static cell identifiers for this tutorial, not source environmental data")
      }
      ncdf4::ncatt_put(output, "depth", "positive", "down")
      for (name in c("Conventions", "institution", "source", "references",
                     "subset:productId", "subset:datasetId", "subset:date")) {
        attribute <- ncdf4::ncatt_get(input, 0, name)
        if (attribute$hasatt) ncdf4::ncatt_put(output, 0, name, attribute$value)
      }
      ncdf4::ncatt_put(output, 0, "source_file", basename(source_file))
      ncdf4::ncatt_put(output, 0, "source_md5", unname(tools::md5sum(source_file)))
      ncdf4::ncatt_put(output, 0, "history", paste(
        "Subset prepared for the ocean3d point tutorial; original temperature values retained.",
        if (teaching_layers) "Artificial static cell_id added and explicitly labelled." else ""))
      ncdf4::ncatt_put(output, 0, "licence",
        "https://marine.copernicus.eu/user-corner/service-commitments-and-licence")
      ncdf4::ncatt_put(output, 0, "credit",
        "E.U. Copernicus Marine Service Information; https://doi.org/10.48670/mds-00375")
    }, finally = ncdf4::nc_close(output))
  }
  main <- file.path(destination, "temperature-points.nc")
  write_subset(main, selected$time, teaching_layers = TRUE)
  for (i in c(1L, 3L)) write_subset(file.path(destination,
    paste0("temperature-", c("2026-09-19", "2026-09-20", "2026-09-21")[i], ".nc")), selected$time[i])

  # Retain actual messages as well as values for finite interactive examples.
  call_point <- function(data, file = main, method = "nearest", ...) {
    messages <- warnings <- character()
    result <- withCallingHandlers(ocean3d::extract_to_point(data = data, nc = file,
      var = "thetao", method = method, output_col = "temperature", ...),
      message = function(c) {messages <<- c(messages, conditionMessage(c)); invokeRestart("muffleMessage")},
      warning = function(c) {warnings <<- c(warnings, conditionMessage(c)); invokeRestart("muffleWarning")})
    list(result = result, messages = messages, warnings = warnings)
  }
  observation <- data.frame(lon = 0, lat = 38, depth = 60,
    date = as.POSIXct("2026-09-20", tz = "UTC"))
  a <- call_point(observation, source_file, method = "all")$result
  b <- call_point(observation, main, method = "all")$result
  stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-10)))

  midpoint_bounds <- function(x) {
    mid <- (head(x, -1L) + tail(x, -1L)) / 2
    cbind(c(x[1] - diff(x)[1] / 2, mid), c(mid, tail(x, 1) + tail(diff(x), 1) / 2))
  }
  bounds <- lapply(centres, midpoint_bounds)
  projection <- "+proj=aeqd +lat_0=38 +lon_0=0 +datum=WGS84 +units=m +no_defs"
  project <- function(x, y) {
    sf::st_coordinates(sf::st_transform(sf::st_sfc(lapply(seq_along(x), function(i) {
      sf::st_point(c(x[i], y[i]))
    }), crs = 4326), projection))[, 1:2, drop = FALSE] / 1000
  }
  values <- read_variable("thetao", selected)
  base_z <- which.min(abs(centres$depth - 60))
  cells <- list()
  for (i in seq_along(centres$longitude)) for (j in seq_along(centres$latitude)) {
    x <- bounds$longitude[i, ]; y <- bounds$latitude[j, ]
    cells[[length(cells) + 1L]] <- list(id = length(cells) + 1L, i = i, j = j,
      lon = centres$longitude[i], lat = centres$latitude[j],
      centre = unname(as.numeric(project(centres$longitude[i], centres$latitude[j]))),
      polygon = unname(project(c(x[1], x[2], x[2], x[1], x[1]), c(y[1], y[1], y[2], y[2], y[1]))),
      value = values[i, j, base_z, 2])
  }
  locations <- data.frame(key = c("inside", "outside"),
    label = c("Inside the grid", "Just outside the cropped grid"),
    lon = c(0, max(bounds$longitude) + 0.03), lat = c(38, 38))
  horizontal <- list()
  for (i in seq_len(nrow(locations))) for (policy in c("within", "nearest")) {
    for (limit in c(1, 20)) {
      obs <- observation
      obs$lon <- locations$lon[i]; obs$lat <- locations$lat[i]
      candidate_distances <- as.numeric(sf::st_distance(
        sf::st_sfc(sf::st_point(c(obs$lon, obs$lat)), crs = 4326),
        sf::st_sfc(lapply(cells, function(c) sf::st_point(c(c$lon, c$lat))), crs = 4326))) / 1000
      nearest <- which.min(candidate_distances)
      out <- call_point(obs, spatial_match = policy,
        max_distance_km = if (policy == "nearest") limit else NULL)
      horizontal[[length(horizontal) + 1L]] <- list(
        location = locations$key[i], policy = policy, limit = limit,
        point = unname(as.numeric(project(obs$lon, obs$lat))), cell = nearest,
        distance_km = candidate_distances[nearest], value = out$result$temperature,
        messages = out$messages, warnings = out$warnings)
      if (is.finite(out$result$temperature)) stopifnot(isTRUE(all.equal(
        out$result$temperature, cells[[nearest]]$value, tolerance = 1e-10)))
    }
  }
  inside <- horizontal[[1]]$cell
  cell <- cells[[inside]]
  profiles <- lapply(seq_along(centres$time), function(i) list(
    date = format(as.POSIXct(centres$time[i], origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d"),
    values = unname(values[cell$i, cell$j, , i])))
  vertical <- list()
  for (date in seq_along(profiles)) for (target in c(60, 100, 800)) {
    obs <- observation; obs$depth <- target
    obs$date <- as.POSIXct(profiles[[date]]$date, tz = "UTC")
    out <- call_point(obs, method = "all")
    indices <- vapply(c("nearest_temperature", "surface_temperature", "seabottom_temperature"),
      function(name) match(out$result[[name]], profiles[[date]]$values), integer(1))
    stopifnot(!anyNA(indices))
    vertical[[length(vertical) + 1L]] <- list(date = profiles[[date]]$date,
      target = target, nearest = unname(indices[1]), surface = unname(indices[2]), bottom = unname(indices[3]),
      warnings = out$warnings, messages = out$messages)
  }
  requests <- c("2026-09-20 00:00:00", "2026-09-20 12:00:00", "2026-09-25 00:00:00")
  baseline <- vapply(profiles, function(p) p$values[base_z], numeric(1))
  stopifnot(!anyDuplicated(baseline))
  temporal <- list()
  for (request in requests) for (precision in c("datetime", "daily")) {
    for (matching in c("nearest", "exact")) for (limit in c("none", "0.1", "1")) {
      obs <- observation; obs$date <- as.POSIXct(request, tz = "UTC")
      out <- call_point(obs, time_resolution = precision, time_match = matching,
        max_time_diff = if (limit == "none") NULL else as.numeric(limit))
      chosen <- match(out$result$temperature, baseline)
      temporal[[length(temporal) + 1L]] <- list(request = request, precision = precision,
        matching = matching, limit = limit, value = out$result$temperature,
        selected = if (is.na(chosen)) NULL else profiles[[chosen]]$date,
        warnings = out$warnings, messages = out$messages)
    }
  }
  example <- list(cells = cells, locations = locations, horizontal = horizontal,
    profiles = profiles, depths = unname(centres$depth), depth_bounds = unname(bounds$depth),
    vertical = vertical, temporal = temporal, selected_cell = inside,
    source = basename(source_file), main_md5 = unname(tools::md5sum(main)))
  saveRDS(example, file.path(destination, "point-display.rds"), version = 2)
  invisible(example)
}
