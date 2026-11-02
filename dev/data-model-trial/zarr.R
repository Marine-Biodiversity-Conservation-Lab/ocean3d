# Optional part of the data model trial (#112): the real Copernicus Marine
# global daily temperature product, read as a remote Zarr store.
#
# cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m: 1/12 degree, 50 levels,
# daily from June 2022 (about 1600 days, 2.9 TB as float32). The ARCO stores
# can be read without a login. Two layouts of the same data:
#   timeChunked  chunks (time 1, depth 1, lat 512, lon 2048): maps
#   geoChunked   chunks (time 134, depth 2, lat 32, lon 64): time series
#
# Usage (from the package root):
#   dev/data-model-trial/measure.sh <label> \
#     Rscript dev/data-model-trial/zarr.R <lib> <step>
# step: open, series, map, volume or ncdf4

args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[[1]], .libPaths()))
step <- args[[2]]
suppressPackageStartupMessages({
  library(stars)
  library(terra)
})

base <- paste0(
  "https://s3.waw3-1.cloudferro.com/mdl-arco-%s/arco/",
  "GLOBAL_ANALYSISFORECAST_PHY_001_024/",
  "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m_202406/%s.zarr"
)
# GDAL needs the ZARR:"..." prefix: listing the bucket directory gives 403.
zarr_path <- function(...) sprintf("ZARR:\"/vsicurl/%s\"", sprintf(base, ...))
time_store <- zarr_path("time-012", "timeChunked")
geo_store <- zarr_path("geo-012", "geoChunked")

open_store <- function(store) {
  suppressWarnings(read_mdim(store, variable = "thetao", proxy = TRUE,
                             normalize_path = FALSE))
}

# Dimension order in the store is time, elevation, latitude, longitude;
# stars reverses it to longitude, latitude, elevation, time.
axes <- function(p) {
  d <- st_dimensions(p)
  list(
    lon = st_get_dimension_values(p, "longitude", center = TRUE),
    lat = st_get_dimension_values(p, "latitude", center = TRUE),
    z = as.numeric(st_get_dimension_values(p, "elevation")),
    time = st_get_dimension_values(p, "time"),
    dims = d
  )
}

read_block <- function(store, offset, count) {
  x <- suppressWarnings(read_mdim(store, variable = "thetao",
                                  offset = offset, count = count,
                                  normalize_path = FALSE))
  x
}

out <- switch(step,
  open = {
    p <- open_store(time_store)
    a <- axes(p)
    sprintf("dims %s; time %s .. %s (%s); elevation %s .. %s",
            paste(dim(p), collapse = "x"), format(min(a$time)),
            format(max(a$time)), class(a$time)[1],
            signif(min(a$z), 4), signif(max(a$z), 4))
  },
  series = {
    # Surface temperature at one point for every day, from the geo layout.
    p <- open_store(geo_store)
    a <- axes(p)
    ix <- which.min(abs(a$lon - -150))
    iy <- which.min(abs(a$lat - 20))
    iz <- which.max(a$z)  # elevation closest to 0 = shallowest level
    x <- read_block(geo_store, c(ix, iy, iz, 1) - 1,
                    c(1, 1, 1, length(a$time)))
    v <- as.numeric(x[[1]])
    sprintf("n=%d mean=%.3f na=%d", length(v), mean(v, na.rm = TRUE),
            sum(is.na(v)))
  },
  map = {
    # One day's global surface map, from the time layout.
    p <- open_store(time_store)
    a <- axes(p)
    it <- which(as.Date(a$time) == as.Date("2023-06-15"))
    iz <- which.max(a$z)
    x <- read_block(time_store, c(0, 0, iz - 1, it - 1),
                    c(length(a$lon), length(a$lat), 1, 1))
    v <- as.numeric(x[[1]])
    sprintf("n=%d mean=%.3f na=%d", length(v), mean(v, na.rm = TRUE),
            sum(is.na(v)))
  },
  volume = {
    # Habitat volume and overlap for one day, all 50 levels, 10 x 10 degrees.
    devtools::load_all(".", quiet = TRUE)
    p <- open_store(time_store)
    a <- axes(p)
    it <- which(as.Date(a$time) == as.Date("2023-06-15"))
    ix <- which(a$lon >= -155 & a$lon < -145)
    iy <- which(a$lat >= 15 & a$lat < 25)
    x <- read_block(time_store, c(min(ix), min(iy), 1, it) - 1,
                    c(length(ix), length(iy), length(a$z), 1))
    # The float32 coordinates make stars treat the grid as rectilinear, and
    # as(, "SpatRaster") then refuses it ("only regular rasters can be
    # converted"), so the SpatRaster is built from the array. read_mdim()
    # reports coordinates as cell starts (see cases D02, D04), so the stored
    # values are taken as the cell centres.
    lon <- as.numeric(st_get_dimension_values(x, "longitude", where = "start"))
    lat <- as.numeric(st_get_dimension_values(x, "latitude", where = "start"))
    arr <- array(as.numeric(x[[1]]), dim(x[[1]])[1:3])
    res <- 1 / 12
    r <- terra::rast(
      lapply(seq_len(dim(arr)[3]), function(k) {
        m <- t(arr[, , k])
        if (lat[1] < lat[length(lat)]) m <- m[rev(seq_len(nrow(m))), ]
        terra::rast(m, extent = terra::ext(min(lon) - res / 2,
                                           max(lon) + res / 2,
                                           min(lat) - res / 2,
                                           max(lat) + res / 2),
                    crs = "OGC:CRS84")
      })
    )
    depths <- -a$z  # elevation (positive up) to depth (positive down)
    names(r) <- paste0("thetao_depth=", depths)
    ord <- order(depths)
    r <- r[[ord]]
    depths <- depths[ord]
    v1 <- as_voxel(terra::ifel(r >= 10 & r <= 20, 1, NA), depths = depths)
    v2 <- as_voxel(terra::ifel(r >= 15 & r <= 25, 1, NA), depths = depths)
    ov <- calc_volume_overlap(v1, v2)
    sprintf("volume_10_20_km3=%.1f volume_15_25_km3=%.1f overlap_km3=%.1f",
            volume(v1), volume(v2),
            terra::global(ov[["volume_overlap"]], "sum", na.rm = TRUE)[[1]])
  },
  ncdf4 = {
    # Can netCDF-C (and so ncdf4 / ncdfCF) open the same store via NCZarr?
    url <- paste0(sprintf(base, "time-012", "timeChunked"), "#mode=zarr,s3")
    res <- tryCatch({
      nc <- ncdf4::nc_open(url)
      on.exit(ncdf4::nc_close(nc))
      paste("opened; variables:", toString(names(nc$var)))
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    res
  }
)
cat("RESULT", "zarr", "stars_mdim", step, out, "\n", sep = "\t")
