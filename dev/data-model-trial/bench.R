# Benchmark part of the data model trial (#112).
#
# Runs one step with one reader on one dataset, in its own process, so that
# measure.sh can record wall time and peak memory for exactly that step.
#
# Usage (from the package root):
#   dev/data-model-trial/measure.sh <label> \
#     Rscript dev/data-model-trial/bench.R \
#       <lib> <data_dir> <dataset> <reader> <step>
#
# dataset: "synthetic", synthetic_4d.nc from make_large_file.R (4D, 60 GB
#          as doubles), or "oisst", NOAA OISST v2.1 sst.day.mean.2022.nc and
#          .2023.nc (2D + time).
# reader:  ocean3d, stars_ncdf, stars_mdim or ncdfCF.
# step:    load, open, points, series, volume or buffer; also points_native
#          (stars), lazy_subset and lazy_subset_read (ncdfCF).
#
# Prints one line: RESULT <dataset> <reader> <step> <summary of the result>.
# The summaries let the readers' answers be compared with each other.

args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[[1]], .libPaths()))
data_dir <- args[[2]]
dataset <- args[[3]]
reader <- args[[4]]
step <- args[[5]]

suppressPackageStartupMessages({
  library(terra)
  if (reader %in% c("stars_ncdf", "stars_mdim")) library(stars)
  if (reader == "ncdfCF") library(ncdfCF)
  if (reader %in% c("ocean3d", "terra") || step == "volume") {
    devtools::load_all(".", quiet = TRUE)
  }
})

# Datasets --------------------------------------------------------------------

if (dataset == "synthetic") {
  files <- file.path(data_dir, "synthetic_4d.nc")
  var <- "thetao"
} else {
  files <- file.path(data_dir, "oisst",
                     paste0("sst.day.mean.", 2022:2023, ".nc"))
  var <- "sst"
}
has_depth <- dataset == "synthetic"

# 500 observations at random places, depths and days in 2022-2023, avoiding
# the synthetic land boxes.
set.seed(112)
n_obs <- 500
obs <- data.frame(
  lon = runif(n_obs * 2, -179, 179),
  lat = runif(n_obs * 2, -60, 60),
  depth = runif(n_obs * 2, 0, 400),
  date = as.POSIXct("2022-01-01", tz = "UTC") +
    86400 * sample(0:729, n_obs * 2, replace = TRUE)
)
land <- (obs$lon > -100 & obs$lon < -60 & obs$lat > -50 & obs$lat < 50) |
  (obs$lon > 10 & obs$lon < 40 & obs$lat > -30 & obs$lat < 30)
obs <- obs[!land, ][seq_len(n_obs), ]
if (!has_depth) obs$depth <- NULL
# OISST longitudes run 0-360.
if (dataset == "oisst") obs$lon <- obs$lon %% 360

series <- data.frame(
  lon = if (dataset == "oisst") 210.1 else -150.1, lat = 20.1,
  date = as.POSIXct("2022-01-01", tz = "UTC") + 86400 * (0:729)
)
if (has_depth) series$depth <- 0

day <- as.POSIXct("2022-06-15", tz = "UTC")  # volume step
buffer_centre <- list(lon = if (dataset == "oisst") 210 else -150, lat = 20,
                      depth = 50, radius_km = 200, days = 15,
                      date = as.POSIXct(if (dataset == "oisst") "2022-12-31"
                                        else "2022-06-15", tz = "UTC"))

# Shared glue for the external readers -----------------------------------------
# What the package would still write on top of stars or ncdfCF: pick the
# file for a date, the nearest cell on each axis, and the cells in a buffer.

file_for <- function(date) {
  if (length(files) == 1) return(files)
  files[as.integer(format(date, "%Y")) - 2021L]
}

nearest <- function(axis, q, longitude = FALSE) {
  d <- abs(axis - q)
  if (longitude) d <- pmin(d, abs(axis - q - 360), abs(axis - q + 360))
  which.min(d)
}

great_circle_km <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371.0088 * asin(sqrt(pmin(1, a)))
}

# Index windows for the buffer step on one file's axes.
buffer_window <- function(lon_axis, lat_axis, depth_axis, time_axis) {
  b <- buffer_centre
  dlat <- b$radius_km / 111.2
  dlon <- dlat / cos(b$lat * pi / 180)
  ix <- which(abs(((lon_axis - b$lon + 180) %% 360) - 180) <= dlon)
  iy <- which(abs(lat_axis - b$lat) <= dlat)
  iz <- if (!is.null(depth_axis)) {
    which(depth_axis >= b$depth - 50 & depth_axis <= b$depth + 50)
  }
  it <- which(abs(as.numeric(difftime(time_axis, b$date, units = "days"))) <=
                b$days)
  list(x = ix, y = iy, z = iz, t = it)
}

# Equal-weight mean of an [x, y, (z,) t] block over cells inside the radius.
buffer_mean <- function(block, lon_axis, lat_axis) {
  inside <- outer(lon_axis, lat_axis, function(x, y) {
    great_circle_km(buffer_centre$lon, buffer_centre$lat, x, y) <=
      buffer_centre$radius_km
  })
  keep <- array(inside, dim(block))
  mean(block[keep], na.rm = TRUE)
}

habitat <- function(r, depths) {
  names(r) <- paste0(var, "_depth=", depths)
  a <- as_voxel(terra::ifel(r >= 10 & r <= 20, 1, NA), depths = depths)
  b <- as_voxel(terra::ifel(r >= 15 & r <= 25, 1, NA), depths = depths)
  ov <- calc_volume_overlap(a, b)
  c(volume_10_20_km3 = volume(a), volume_15_25_km3 = volume(b),
    overlap_km3 = terra::global(ov[["volume_overlap"]], "sum",
                                na.rm = TRUE)[[1]])
}

fmt <- function(x) {
  if (is.numeric(x) && length(x) > 3) {
    return(sprintf("n=%d mean=%.4f na=%d", length(x), mean(x, na.rm = TRUE),
                   sum(is.na(x))))
  }
  paste(names(x), signif(x, 7), sep = "=", collapse = " ")
}

# Readers ---------------------------------------------------------------------

run <- list()

## ocean3d: the package's ncdf4 code (and terra for volume) -------------------
run$ocean3d <- list(
  load = function() "loaded",
  open = function() {
    out <- vapply(files, function(f) {
      nc <- ncdf4::nc_open(f)
      on.exit(ncdf4::nc_close(nc))
      dims <- ocean3d:::.detect_dims(nc, var)
      length(ocean3d:::.decode_netcdf_time(nc, dims$time))
    }, numeric(1))
    paste("time steps:", sum(out))
  },
  points = function() {
    out <- suppressMessages(extract_to_point(
      obs, nc = files, var = var, method = if (has_depth) "nearest" else "2d",
      time_report = FALSE
    ))
    fmt(out[[ncol(out)]])
  },
  series = function() {
    out <- suppressMessages(extract_to_point(
      series, nc = files, var = var,
      method = if (has_depth) "surface" else "2d", time_report = FALSE
    ))
    fmt(out[[ncol(out)]])
  },
  volume = function() {
    # The current voxel model: terra reads the file, layers named by depth.
    r <- terra::rast(files, subds = var)
    depths <- unique(terra::depth(r))
    t <- which(terra::time(r) == day)
    fmt(habitat(r[[t]], depths))
  },
  buffer = function() {
    b <- buffer_centre
    d <- data.frame(lon = b$lon, lat = b$lat, date = b$date)
    if (has_depth) d$depth <- b$depth
    out <- suppressMessages(extract_to_buffer(
      d, nc = files, var = var, horizontal_radius_km = b$radius_km,
      depth_radius_m = if (has_depth) 50 else 0,
      time_before = b$days, time_after = b$days, statistic = "mean",
      weighting = "equal", infer_bounds = TRUE, quiet = TRUE,
      return_match_info = FALSE
    ))
    fmt(c(mean = out[[ncol(out)]]))
  }
)

## stars::read_ncdf() ---------------------------------------------------------
.stars_axes <- function(p) {
  nms <- names(st_dimensions(p))
  list(
    x = st_get_dimension_values(p, nms[1], center = TRUE),
    y = st_get_dimension_values(p, nms[2], center = TRUE),
    z = if (has_depth) as.numeric(st_get_dimension_values(p, nms[3])),
    t = st_get_dimension_values(p, nms[length(nms)], where = "start")
  )
}
.open_ncdf <- function(f) {
  suppressMessages(read_ncdf(f, var = var, proxy = TRUE))
}
.ncsub_read <- function(f, ix, iy, iz, it) {
  start <- c(min(ix), min(iy), if (has_depth) min(iz), min(it))
  count <- c(length(ix), length(iy), if (has_depth) length(iz), length(it))
  x <- suppressMessages(read_ncdf(f, var = var, ncsub = cbind(start, count),
                                  proxy = FALSE))
  as.numeric(unclass(x[[1]]))
}
.stars_points <- function(d, read_one, open_one) {
  proxies <- lapply(files, open_one)
  axes <- lapply(proxies, .stars_axes)
  vapply(seq_len(nrow(d)), function(i) {
    k <- match(file_for(d$date[i]), files)
    a <- axes[[k]]
    ix <- nearest(a$x, d$lon[i], longitude = TRUE)
    iy <- nearest(a$y, d$lat[i])
    iz <- if (has_depth) nearest(a$z, d$depth[i])
    it <- nearest(as.numeric(a$t), as.numeric(d$date[i]))
    read_one(files[k], ix, iy, iz, it)
  }, numeric(1))
}
.stars_buffer <- function(read_block, open_one) {
  f <- file_for(buffer_centre$date)
  # A window that crosses a file boundary needs both files.
  vals <- c()
  for (f in files) {
    a <- .stars_axes(open_one(f))
    w <- buffer_window(a$x, a$y, a$z, a$t)
    if (!length(w$t)) next
    block <- read_block(f, w)
    dim(block) <- c(length(w$x), length(w$y),
                    if (has_depth) length(w$z), length(w$t))
    inside <- outer(a$x[w$x], a$y[w$y], function(x, y) {
      great_circle_km(buffer_centre$lon, buffer_centre$lat, x, y) <=
        buffer_centre$radius_km
    })
    vals <- c(vals, block[array(inside, dim(block))])
  }
  fmt(c(mean = mean(vals, na.rm = TRUE)))
}
run$stars_ncdf <- list(
  load = function() "loaded",
  open = function() {
    p <- lapply(files, .open_ncdf)
    paste("time steps:", sum(vapply(p, function(x) dim(x)[[length(dim(x))]],
                                    numeric(1))))
  },
  points = function() {
    fmt(.stars_points(obs, function(f, ix, iy, iz, it) {
      .ncsub_read(f, ix, iy, iz, it)
    }, .open_ncdf))
  },
  series = function() {
    # One read per file of the whole column, rather than 730 single cells.
    out <- unlist(lapply(files, function(f) {
      a <- .stars_axes(.open_ncdf(f))
      ix <- nearest(a$x, series$lon[1], longitude = TRUE)
      iy <- nearest(a$y, series$lat[1])
      iz <- if (has_depth) nearest(a$z, 0)
      .ncsub_read(f, ix, iy, iz, seq_along(a$t))
    }))
    fmt(out)
  },
  volume = function() {
    a <- .stars_axes(.open_ncdf(files))
    it <- which(a$t == day)
    # read_ncdf(ncsub = ) with one time step keeps time 1..730 in the
    # returned dimensions, so adrop() fails; slicing the proxy works.
    x <- st_as_stars(.open_ncdf(files)[, , , , it])
    # terra::rast() fails on a stars object whose third dimension is depth
    # ("1:dd[4]: NA/NaN argument"); as(, "SpatRaster") converts it.
    r <- as(adrop(x), "SpatRaster")
    fmt(habitat(r, a$z))
  },
  buffer = function() {
    .stars_buffer(function(f, w) .ncsub_read(f, w$x, w$y, w$z, w$t),
                  .open_ncdf)
  }
)

## stars::read_mdim() ---------------------------------------------------------
.open_mdim <- function(f) {
  suppressWarnings(read_mdim(f, variable = var, proxy = TRUE))
}
.mdim_read <- function(f, ix, iy, iz, it) {
  offset <- c(min(ix), min(iy), if (has_depth) min(iz), min(it)) - 1
  count <- c(length(ix), length(iy), if (has_depth) length(iz), length(it))
  x <- suppressWarnings(read_mdim(f, variable = var, offset = offset,
                                  count = count))
  as.numeric(unclass(x[[1]]))
}
run$stars_mdim <- list(
  load = function() "loaded",
  open = function() {
    p <- lapply(files, .open_mdim)
    paste("time steps:", sum(vapply(p, function(x) dim(x)[[length(dim(x))]],
                                    numeric(1))))
  },
  points = function() fmt(.stars_points(obs, .mdim_read, .open_mdim)),
  series = function() {
    out <- unlist(lapply(files, function(f) {
      a <- .stars_axes(.open_mdim(f))
      ix <- nearest(a$x, series$lon[1], longitude = TRUE)
      iy <- nearest(a$y, series$lat[1])
      iz <- if (has_depth) nearest(a$z, 0)
      .mdim_read(f, ix, iy, iz, seq_along(a$t))
    }))
    fmt(out)
  },
  volume = function() {
    a <- .stars_axes(.open_mdim(files))
    it <- which(as.Date(a$t) == as.Date(day))
    x <- suppressWarnings(read_mdim(
      files, variable = var, offset = c(0, 0, 0, it - 1),
      count = c(length(a$x), length(a$y), length(a$z), 1)
    ))
    # terra::rast() fails on a stars object whose third dimension is depth
    # ("1:dd[4]: NA/NaN argument"); as(, "SpatRaster") converts it.
    r <- as(adrop(x), "SpatRaster")
    fmt(habitat(r, a$z))
  },
  buffer = function() {
    .stars_buffer(function(f, w) .mdim_read(f, w$x, w$y, w$z, w$t),
                  .open_mdim)
  }
)

## ncdfCF ---------------------------------------------------------------------
.cf_open <- function(f) open_ncdf(f)[[var]]
.cf_axes <- function(v) {
  ax <- v$axes
  list(
    x = ax[[1]]$coordinates, y = ax[[2]]$coordinates,
    z = if (has_depth) ax[[3]]$coordinates,
    t = as.POSIXct(ax[[length(ax)]]$time$as_timestamp(), tz = "UTC"),
    names = names(ax)
  )
}
# subset() keeps cells whose coordinate values lie in each range, so a single
# cell is asked for with a range around its own coordinate.
.cf_cell_args <- function(a, ix, iy, iz, it) {
  win <- function(v, i) c(min(v[i]), max(v[i])) + c(-1e-6, 1e-6)
  out <- list(win(a$x, ix), win(a$y, iy))
  if (has_depth) out <- c(out, list(win(a$z, iz)))
  # Time ranges are half-open even with rightmost.closed = TRUE, so the end
  # is moved one second past the last wanted step.
  out <- c(out, list(format(range(a$t[it]) + c(0, 1), "%Y-%m-%dT%H:%M:%S")))
  names(out) <- a$names
  out
}
.cf_read <- function(v, a, ix, iy, iz, it) {
  # Ranges are half-open unless rightmost.closed = TRUE.
  s <- do.call(v$subset, c(.cf_cell_args(a, ix, iy, iz, it),
                           rightmost.closed = TRUE))
  as.numeric(s$raw())
}
run$ncdfCF <- list(
  load = function() "loaded",
  open = function() {
    v <- lapply(files, .cf_open)
    n <- vapply(v, function(x) x$axes[[length(x$axes)]]$length, numeric(1))
    paste("time steps:", sum(n))
  },
  points = function() {
    vs <- lapply(files, .cf_open)
    axes <- lapply(vs, .cf_axes)
    fmt(vapply(seq_len(nrow(obs)), function(i) {
      k <- match(file_for(obs$date[i]), files)
      a <- axes[[k]]
      .cf_read(vs[[k]], a, nearest(a$x, obs$lon[i], TRUE),
               nearest(a$y, obs$lat[i]),
               if (has_depth) nearest(a$z, obs$depth[i]),
               nearest(as.numeric(a$t), as.numeric(obs$date[i])))
    }, numeric(1)))
  },
  series = function() {
    out <- unlist(lapply(files, function(f) {
      v <- .cf_open(f)
      a <- .cf_axes(v)
      .cf_read(v, a, nearest(a$x, series$lon[1], TRUE),
               nearest(a$y, series$lat[1]), if (has_depth) 1L,
               seq_along(a$t))
    }))
    fmt(out)
  },
  volume = function() {
    v <- .cf_open(files)
    a <- .cf_axes(v)
    sel <- list(format(day, "%Y-%m-%d"))
    names(sel) <- a$names[length(a$names)]
    # terra() returns a SpatRasterDataset without a CRS; cell areas (and so
    # volumes) need lon/lat, so the CRS is set here.
    r <- do.call(v$subset, sel)$terra()[[1]]
    terra::crs(r) <- "OGC:CRS84"
    fmt(habitat(r, a$z))
  },
  buffer = function() {
    vals <- c()
    for (f in files) {
      v <- .cf_open(f)
      a <- .cf_axes(v)
      w <- buffer_window(a$x, a$y, a$z, a$t)
      if (!length(w$t)) next
      s <- do.call(v$subset, c(.cf_cell_args(a, w$x, w$y, w$z, w$t),
                               rightmost.closed = TRUE))
      block <- s$raw()
      dim(block) <- c(length(w$x), length(w$y),
                      if (has_depth) length(w$z), length(w$t))
      inside <- outer(a$x[w$x], a$y[w$y], function(x, y) {
        great_circle_km(buffer_centre$lon, buffer_centre$lat, x, y) <=
          buffer_centre$radius_km
      })
      vals <- c(vals, block[array(inside, dim(block))])
    }
    fmt(c(mean = mean(vals, na.rm = TRUE)))
  }
)

# stars' own point extraction on the proxy, for the first two observations
# only: on the synthetic file it reads far more than the needed cells.
.stars_native <- function(open_one) {
  d <- obs[1:2, ]
  k <- match(file_for(d$date[1]), files)
  d <- d[vapply(d$date, file_for, "") == files[k], , drop = FALSE]
  p <- open_one(files[k])
  pts <- st_sf(date = d$date,
               geometry = st_as_sfc(lapply(seq_len(nrow(d)), function(i) {
                 st_point(c(d$lon[i], d$lat[i]))
               }), crs = st_crs(p)))
  e <- st_extract(p, pts, time_column = "date")
  v <- as.numeric(e[[1]])
  sprintf("points=%d values=%d mean=%.4f", nrow(d), length(v),
          mean(v, na.rm = TRUE))
}
# Is ncdfCF's subset() lazy? The whole globe at the top level for every day
# is 757 million cells (6 GB as doubles). lazy_subset only subsets;
# lazy_subset_read then subsets that result again to two days and reads it.
.cf_top_level <- function() {
  v <- .cf_open(files)
  a <- .cf_axes(v)
  sel <- list(c(a$z[1] - 1e-3, a$z[1] + 1e-3))
  names(sel) <- a$names[3]
  list(s = do.call(v$subset, sel), a = a)
}
run$ncdfCF$lazy_subset <- function() {
  s <- .cf_top_level()$s
  paste("subset dims", paste(s$dim(), collapse = "x"))
}
run$ncdfCF$lazy_subset_read <- function() {
  top <- .cf_top_level()
  sel <- list(c("2022-06-15", "2022-06-17"))
  names(sel) <- top$a$names[4]
  s2 <- do.call(top$s$subset, sel)
  v <- as.numeric(s2$raw())
  sprintf("dims %s, read %d values, mean=%.4f",
          paste(s2$dim(), collapse = "x"), length(v), mean(v, na.rm = TRUE))
}

run$stars_ncdf$points_native <- function() .stars_native(.open_ncdf)
run$stars_mdim$points_native <- function() .stars_native(.open_mdim)

result <- tryCatch(run[[reader]][[step]](),
                   error = function(e) paste("ERROR:", conditionMessage(e)))
cat("RESULT", dataset, reader, step, result, "\n", sep = "\t")
