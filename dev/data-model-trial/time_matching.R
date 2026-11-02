# #89: matching 500 Date observations to a 10-year daily time axis.
# The package's current matcher (through extract2d()) against CFtime, which
# SPEC.md plans as the time decoder and which stars and ncdfCF both use.
#
# Usage (from the package root):
#   Rscript dev/data-model-trial/time_matching.R <lib>

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) .libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

f <- tempfile(fileext = ".nc")
lon <- ncdf4::ncdim_def("lon", "degrees_east", c(0, 1))
lat <- ncdf4::ncdim_def("lat", "degrees_north", c(40, 41))
tm <- ncdf4::ncdim_def("time", "days since 2010-01-01", 0:3649)
nc <- ncdf4::nc_create(f, ncdf4::ncvar_def("temp", "1", list(lon, lat, tm),
                                           -9999, prec = "double"))
ncdf4::ncvar_put(nc, "temp", array(rep(1:3650, each = 4), c(2, 2, 3650)))
ncdf4::nc_close(nc)

set.seed(1)
d <- data.frame(lon = 0, lat = 40,
                date = as.Date("2010-01-01") + sample(0:3649, 500, TRUE))

t_pkg <- system.time(
  out <- suppressMessages(extract2d(d, f, var = "temp", time_report = FALSE))
)[["elapsed"]]

t_cf <- system.time({
  cf <- CFtime::CFtime("days since 2010-01-01", "standard", 0:3649)
  idx <- CFtime::indexOf(format(d$date), cf)
})[["elapsed"]]

expected <- as.numeric(d$date - as.Date("2010-01-01")) + 1
cat(sprintf("package extract2d(): %.2f s, correct %d of 500\n", t_pkg,
            sum(out$temp == expected)))
cat(sprintf("CFtime::indexOf():   %.3f s, correct %d of 500\n", t_cf,
            sum(idx == expected)))
