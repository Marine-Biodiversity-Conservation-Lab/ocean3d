# Write a synthetic 4D daily global file shaped like a Copernicus Marine
# download: lon x lat x depth x time, short integers with scale_factor and
# add_offset, _FillValue on land, deflate, one chunk per depth level and day.
#
# Default size: 1440 x 720 x 10 depths x 730 days. Packed it is 15 GB; read
# into R as doubles it would be 60 GB, more than the 31 GB of RAM on the
# trial machine, so any step that reads it all fails.
#
# Usage: Rscript make_large_file.R <out.nc> [n_days] [n_depth]

args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
n_days <- if (length(args) >= 2) as.integer(args[[2]]) else 730L
n_depth <- if (length(args) >= 3) as.integer(args[[3]]) else 10L

lon <- seq(-179.875, 179.875, by = 0.25)
lat <- seq(-89.875, 89.875, by = 0.25)
depth <- c(0.494, 5.08, 11.4, 21.6, 34.4, 55.8, 92.3, 155.9, 266.0, 453.9)
depth <- depth[seq_len(n_depth)]

d_lon <- ncdf4::ncdim_def("longitude", "degrees_east", lon)
d_lat <- ncdf4::ncdim_def("latitude", "degrees_north", lat)
d_dep <- ncdf4::ncdim_def("depth", "m", depth)
d_tim <- ncdf4::ncdim_def(
  "time", "hours since 1950-01-01", 631152 + 24 * (seq_len(n_days) - 1),
  unlim = FALSE
)
v <- ncdf4::ncvar_def(
  "thetao", "degrees_C", list(d_lon, d_lat, d_dep, d_tim),
  missval = -32767L, prec = "short", compression = 1,
  chunksizes = c(length(lon), length(lat), 1L, 1L)
)
nc <- ncdf4::nc_create(out, v, force_v4 = TRUE)
ncdf4::ncatt_put(nc, "thetao", "scale_factor", 0.001)
ncdf4::ncatt_put(nc, "thetao", "add_offset", 0)
ncdf4::ncatt_put(nc, "thetao", "standard_name",
                 "sea_water_potential_temperature")
ncdf4::ncatt_put(nc, "depth", "positive", "down")
ncdf4::ncatt_put(nc, "depth", "axis", "Z")
ncdf4::ncatt_put(nc, "depth", "standard_name", "depth")
ncdf4::ncatt_put(nc, "time", "calendar", "gregorian")
ncdf4::ncatt_put(nc, "time", "axis", "T")
ncdf4::ncatt_put(nc, "time", "standard_name", "time")
ncdf4::ncatt_put(nc, "longitude", "axis", "X")
ncdf4::ncatt_put(nc, "latitude", "axis", "Y")
ncdf4::ncatt_put(nc, 0, "Conventions", "CF-1.8")

# Smooth field: warm at the equator, cooling with depth, seasonal cycle, plus
# seeded noise of +-0.05 so the file compresses about as well as model output.
# A crude land mask (two rectangles) gives _FillValue cells.
lat_term <- outer(rep(1, length(lon)), 28 * cos(lat * pi / 180))
land <- outer(lon > -100 & lon < -60, lat > -50 & lat < 50) |
  outer(lon > 10 & lon < 40, lat > -30 & lat < 30)
for (t in seq_len(n_days)) {
  season <- 2 * sin(2 * pi * (t - 1) / 365.25)
  for (k in seq_along(depth)) {
    field <- (lat_term + season) * exp(-depth[k] / 500) +
      0.001 * t + 0.01 * k
    set.seed(t * 100 + k)
    field <- field + runif(length(field), -0.05, 0.05)
    field <- round(field / 0.001)
    field[land] <- NA
    ncdf4::ncvar_put(nc, "thetao", field,
                     start = c(1, 1, k, t), count = c(-1, -1, 1, 1))
  }
  if (t %% 50 == 0) message("day ", t, " of ", n_days)
}
ncdf4::nc_close(nc)
message("wrote ", out, " (", round(file.size(out) / 1e9, 2), " GB)")
