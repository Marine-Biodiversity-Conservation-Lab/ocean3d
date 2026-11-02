# Correctness part of the data model trial (#112).
#
# Reads every case in cases.R with four readers and records what each one
# decodes (time labels, vertical values with units and direction, cell
# bounds) and the value it gives at the query cell:
#   stars_ncdf  stars::read_ncdf()   (ncmeta + RNetCDF)
#   stars_mdim  stars::read_mdim()   (GDAL multidimensional API)
#   ncdfCF      ncdfCF::open_ncdf()
#   terra       terra::rast() (GDAL netCDF driver), depth() and time()
#   ocean3d     the package's ncdf4 code (extract_to_point() and its helpers)
#
# For the three external readers the cell is picked from the coordinates the
# reader returns, as returned: no unit conversion, no flip for
# positive = "up", no longitude wrapping beyond +-360. Bounds are used when
# the reader returns them. So a wrong value means the reader's own output
# isn't enough and package code would still be needed on top.
#
# Usage (from the package root): Rscript dev/data-model-trial/correctness.R
#   [lib]   optional library path holding stars, ncdfCF and CFtime
# Writes dev/data-model-trial/results/correctness.csv

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) .libPaths(c(args[[1]], .libPaths()))

suppressPackageStartupMessages({
  library(stars)
  library(ncdfCF)
  devtools::load_all(".", quiet = TRUE)
})
source("dev/data-model-trial/cases.R")

# Coordinates as each reader returns them ------------------------------------

# One entry per dimension, in the order of the value array.
.dim_entry <- function(name, role, centre, start = NULL, end = NULL,
                       units = NULL, positive = NULL) {
  list(name = name, role = role, centre = centre, start = start, end = end,
       units = units, positive = positive)
}

.stars_coords <- function(x) {
  d <- st_dimensions(x)
  xy <- attr(d, "raster")$dimensions
  dims <- lapply(names(d), function(nm) {
    centre <- st_get_dimension_values(x, nm, center = TRUE)
    is_time <- inherits(centre, c("POSIXct", "Date", "PCICt"))
    # Time coordinates are instants: use the stored value, not the centre
    # of the interval stars builds around a regular axis.
    if (is_time) centre <- st_get_dimension_values(x, nm, where = "start")
    # Non-standard calendars: stars keeps a CFtime object in the dimension
    # and returns raw offsets from st_get_dimension_values().
    if (inherits(d[[nm]]$values, "CFTime")) {
      is_time <- TRUE
      centre <- d[[nm]]$values$as_timestamp()
    }
    role <- if (nm == xy[1]) "x" else if (nm == xy[2]) "y" else
      if (is_time) "t" else "z"
    start <- end <- NULL
    if (role != "t") {
      start <- tryCatch(st_get_dimension_values(x, nm, where = "start"),
                        error = function(e) NULL)
      end <- tryCatch(st_get_dimension_values(x, nm, where = "end"),
                      error = function(e) NULL)
      if (identical(start, centre)) start <- end <- NULL
    }
    units <- if (inherits(centre, "units")) units::deparse_unit(centre)
    .dim_entry(nm, role, if (is_time) centre else as.numeric(centre),
               if (!is.null(start)) as.numeric(start),
               if (!is.null(end)) as.numeric(end), units)
  })
  list(dims = dims, values = as.numeric(unclass(x[[1]])),
       shape = dim(x[[1]]))
}

read_stars_ncdf <- function(path, var) {
  .stars_coords(suppressMessages(read_ncdf(path, var = var, proxy = FALSE)))
}

read_stars_mdim <- function(path, var) {
  .stars_coords(suppressWarnings(read_mdim(path, variable = var)))
}

read_ncdfcf <- function(path, var) {
  ds <- open_ncdf(path)
  v <- ds[[var]]
  dims <- lapply(v$axes, function(ax) {
    role <- if (inherits(ax, "CFAxisLongitude")) "x" else
      if (inherits(ax, "CFAxisLatitude")) "y" else
        if (inherits(ax, "CFAxisVertical")) "z" else
          if (inherits(ax, "CFAxisTime")) "t" else
            switch(ax$orientation, X = "x", Y = "y", Z = "z", T = "t", "?")
    centre <- ax$coordinates
    if (role == "t") {
      centre <- ax$time$as_timestamp()
    }
    start <- end <- NULL
    if (!is.null(ax$bounds) && role != "t") {
      b <- ax$bounds$values
      if (is.matrix(b)) {
        start <- b[1, ]
        end <- b[2, ]
      }
    }
    .dim_entry(ax$name, role, centre, start, end,
               ax$attribute("units"), ax$attribute("positive"))
  })
  list(dims = dims, values = as.numeric(v$raw()), shape = v$dim())
}

# terra (GDAL's netCDF driver): one layer per depth and time step, with the
# depth and time of each layer in terra::depth() and terra::time().
read_terra <- function(path, var) {
  r <- suppressWarnings(terra::rast(path, subds = var))
  x <- terra::xFromCol(r, seq_len(terra::ncol(r)))
  y <- terra::yFromRow(r, seq_len(terra::nrow(r)))
  hx <- terra::xres(r) / 2
  hy <- terra::yres(r) / 2
  z_layer <- terra::depth(r)
  t_layer <- terra::time(r)
  has_z <- length(unique(z_layer)) > 1 || !all(z_layer == 0)
  has_t <- !all(is.na(t_layer))
  z <- if (has_z) unique(z_layer) else 0
  t <- if (has_t) unique(t_layer) else NA
  arr <- array(NA_real_, c(length(x), length(y), length(z), length(t)))
  # Cells (row by row from the top) x layers.
  vals <- terra::values(r, mat = TRUE)
  for (l in seq_len(terra::nlyr(r))) {
    zi <- if (has_z) match(z_layer[l], z) else 1L
    ti <- if (has_t) match(t_layer[l], t) else 1L
    arr[, , zi, ti] <- matrix(vals[, l], length(x), length(y))
  }
  dims <- list(
    .dim_entry("x", "x", x, x - hx, x + hx),
    .dim_entry("y", "y", y, y - hy, y + hy)
  )
  if (has_z) dims <- c(dims, list(.dim_entry("z", "z", z)))
  if (has_t) dims <- c(dims, list(.dim_entry("t", "t", t)))
  shape <- c(length(x), length(y), if (has_z) length(z), if (has_t) length(t))
  list(dims = dims, values = as.numeric(arr), shape = shape)
}

# Picking the query cell from returned coordinates ---------------------------

.in_bounds <- function(start, end, q) {
  lo <- pmin(start, end)
  hi <- pmax(start, end)
  which(lo <= q & q < hi)
}

.pick_horizontal <- function(d, q, longitude) {
  if (!is.null(d$start)) {
    shifts <- if (longitude) c(0, -360, 360) else 0
    for (s in shifts) {
      i <- .in_bounds(d$start, d$end, q + s)
      if (length(i)) return(i[1])
    }
    return(NA_integer_)
  }
  dist <- abs(d$centre - q)
  if (longitude) dist <- pmin(dist, abs(d$centre - q - 360),
                              abs(d$centre - q + 360))
  which.min(dist)
}

.as_utc <- function(t) {
  if (inherits(t, "POSIXct")) return(t)
  if (inherits(t, "Date")) return(as.POSIXct(t, tz = "UTC"))
  out <- tryCatch(
    suppressWarnings(as.POSIXct(as.character(t), tz = "UTC",
                                tryFormats = c("%Y-%m-%dT%H:%M:%OS",
                                               "%Y-%m-%d %H:%M:%OS",
                                               "%Y-%m-%d"))),
    error = function(e) NA
  )
  if (anyNA(out)) NULL else out
}

.pick_time <- function(d, q) {
  t <- .as_utc(d$centre)
  if (is.character(q) || is.null(t)) {
    labels <- substr(as.character(d$centre), 1, 10)
    return(match(as.character(q), labels))
  }
  if (isTRUE(attr(q, "exact"))) {
    return(match(as.numeric(q), as.numeric(t)))
  }
  vapply(as.numeric(q), function(v) which.min(abs(as.numeric(t) - v)),
         integer(1))
}

value_at <- function(coords, query) {
  idx <- lapply(coords$dims, function(d) {
    switch(d$role,
      x = .pick_horizontal(d, query$lon, longitude = TRUE),
      y = .pick_horizontal(d, query$lat, longitude = FALSE),
      z = if (is.null(query$depth)) 1L else
        if (!is.null(d$start)) .in_bounds(d$start, d$end, query$depth)[1] else
          which.min(abs(d$centre - query$depth)),
      t = {
        q <- query$time
        if (isTRUE(query$exact)) attr(q, "exact") <- TRUE
        .pick_time(d, q)
      },
      1L
    )
  })
  if (!any(vapply(coords$dims, `[[`, "", "role") == "x")) {
    stop("no longitude axis found", call. = FALSE)
  }
  nt <- max(lengths(idx))
  arr <- array(coords$values, coords$shape)
  out <- vapply(seq_len(nt), function(k) {
    i <- vapply(idx, function(v) v[min(k, length(v))], integer(1))
    if (anyNA(i)) NA_real_ else arr[matrix(i, 1)]
  }, numeric(1))
  if (isTRUE(query$exact)) sum(is.na(out)) else out
}

.describe <- function(coords) {
  parts <- vapply(coords$dims, function(d) {
    if (d$role == "t") {
      lab <- if (inherits(d$centre, "POSIXct")) {
        format(d$centre, "%Y-%m-%d %H:%M:%OS3", tz = "UTC")
      } else {
        as.character(d$centre)
      }
      cls <- class(d$centre)[1]
      return(paste0("t[", cls, "]: ", paste(head(lab, 3), collapse = ", ")))
    }
    out <- paste0(d$role, ": ", paste(head(signif(d$centre, 6), 4),
                                      collapse = ", "))
    if (!is.null(d$units)) out <- paste0(out, " [", d$units, "]")
    if (!is.null(d$positive)) out <- paste0(out, " positive=", d$positive)
    if (!is.null(d$start)) {
      out <- paste0(out, " bounds ", paste(head(signif(d$start, 4), 3),
                                           collapse = ","), "..")
    }
    out
  }, character(1))
  paste(parts, collapse = "; ")
}

# The package -----------------------------------------------------------------

run_ocean3d <- function(path, var, query) {
  method <- if (is.null(query$depth)) "2d" else "nearest"
  date <- query$time
  if (is.character(date)) {
    stop("query time can't be represented as an R date-time", call. = FALSE)
  }
  args <- list(nc = path, var = var, lon = query$lon, lat = query$lat,
               method = method, time_report = FALSE)
  n <- max(1L, length(date))
  if (n > 1) {
    args$lon <- rep(query$lon, n)
    args$lat <- rep(query$lat, n)
  }
  if (!is.null(query$depth)) args$depth <- query$depth
  if (!is.null(date)) args$date <- date
  if (isTRUE(query$exact)) {
    args$time_match <- "exact"
    args$time_resolution <- "datetime"
  }
  value <- suppressMessages(do.call(extract_to_point, args))
  if (is.data.frame(value)) value <- value[[ncol(value)]]
  if (isTRUE(query$exact)) sum(is.na(value)) else value
}

describe_ocean3d <- function(path, var) {
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  dims <- ocean3d:::.detect_dims(nc, var)
  out <- character()
  if (!is.null(dims$time)) {
    t <- ocean3d:::.decode_netcdf_time(nc, dims$time)
    out <- c(out, paste0("t: ", paste(head(format(t, "%Y-%m-%d %H:%M:%OS3",
                                                  tz = "UTC"), 3),
                                      collapse = ", ")))
  }
  if (!is.null(dims$depth)) {
    z <- ocean3d:::.read_depth_axis(nc, dims$depth,
                                    ocean3d:::.get_dim_values(nc, dims$depth))
    out <- c(out, paste0("z (m, down): ",
                         paste(head(signif(z$depth_m, 6), 4),
                               collapse = ", ")))
  }
  paste(c(paste0("x=", dims$lon, " y=", dims$lat), out), collapse = "; ")
}

# Run -------------------------------------------------------------------------

readers <- list(
  stars_ncdf = read_stars_ncdf,
  stars_mdim = read_stars_mdim,
  ncdfCF = read_ncdfcf,
  terra = read_terra
)

.capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(conditionMessage(e), class = "trial_error")
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

.verdict <- function(case, value) {
  if (inherits(value, "trial_error")) {
    return(if (isTRUE(case$error_ok)) "ok (error)" else "error")
  }
  v <- as.numeric(value)
  if (length(v) != 1 || is.na(v)) return("wrong")
  if (isTRUE(all.equal(v, case$expected))) "ok" else "wrong"
}

.clean <- function(x) gsub("[\r\n|]+", " ", substr(x, 1, 160))

rows <- list()
dir <- tempfile("ocean3d-trial-")
dir.create(dir)
for (case in cases) {
  if (isTRUE(case$multi)) next
  path <- file.path(dir, paste0(case$id, ".nc"))
  case$build(path)
  for (r in names(readers)) {
    got <- .capture({
      coords <- readers[[r]](path, "temp")
      structure(value_at(coords, case$query), desc = .describe(coords))
    })
    value <- got$value
    rows[[length(rows) + 1]] <- data.frame(
      id = case$id, group = case$group, from = case$from, what = case$what,
      reader = r, expected = case$expected,
      value = if (inherits(value, "trial_error")) NA else
        paste(signif(as.numeric(value), 6), collapse = ","),
      verdict = .verdict(case, value),
      decoded = if (inherits(value, "trial_error")) "" else
        .clean(attr(value, "desc")),
      message = .clean(paste(c(if (inherits(value, "trial_error")) value,
                               got$warnings), collapse = " / "))
    )
  }
  got <- .capture(run_ocean3d(path, "temp", case$query))
  desc <- .capture(describe_ocean3d(path, "temp"))$value
  value <- got$value
  rows[[length(rows) + 1]] <- data.frame(
    id = case$id, group = case$group, from = case$from, what = case$what,
    reader = "ocean3d", expected = case$expected,
    value = if (inherits(value, "trial_error")) NA else
      paste(signif(as.numeric(value), 6), collapse = ","),
    verdict = .verdict(case, value),
    decoded = if (inherits(desc, "trial_error")) "" else .clean(desc),
    message = .clean(paste(c(if (inherits(value, "trial_error")) value,
                             got$warnings), collapse = " / "))
  )
}

results <- do.call(rbind, rows)

# Verdicts changed after reading the decoded output: the right cell for the
# wrong reason.
overrides <- list(
  list("T11", "stars_mdim", "wrong",
       paste("time not decoded (no calendar attribute): every query fell",
             "back to step 1, so none came back NA")),
  list("D06", "ocean3d", "differs (by design)",
       paste("method = 'nearest' picks the nearest layer centre; depth",
             "bounds are not used")),
  list("D06", "stars_mdim", "wrong",
       paste("right cell by chance: intervals [20, 70), [70, 120) built",
             "from the midpoints of depth_bnds"))
)
for (o in overrides) {
  i <- results$id == o[[1]] & results$reader == o[[2]]
  results$verdict[i] <- o[[3]]
  results$message[i] <- o[[4]]
}
out_dir <- "dev/data-model-trial/results"
dir.create(out_dir, showWarnings = FALSE)
write.csv(results, file.path(out_dir, "correctness.csv"), row.names = FALSE)

wide <- reshape(results[c("id", "reader", "verdict")], idvar = "id",
                timevar = "reader", direction = "wide")
names(wide) <- sub("^verdict\\.", "", names(wide))
print(wide, row.names = FALSE)
