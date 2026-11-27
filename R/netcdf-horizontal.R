# Longitude helpers for axes read by .read_interval_axis() and for shapes, in
# either convention (0-360 or -180-180), and the tolerance for a single
# coordinate without bounds.


.normalise_longitude <- function(values) {
  (values + 180) %% 360 - 180
}

# A single horizontal coordinate without bounds (a point time series or a
# transect) has no width. A point matches it when it lies on the coordinate,
# allowing for single-precision rounding of the stored value but no real
# offset. Used by extract() on a SpatProxy.
.singleton_tolerance <- function(centres) {
  pmax(1e-8, abs(centres) * 2^-23)
}


.unwrap_longitude <- function(values) {
  if (length(values) < 2L) return(values)
  original_steps <- diff(values)
  steps <- (original_steps + 180) %% 360 - 180
  ambiguous <- abs(abs(steps) - 180) < 1e-8
  steps[ambiguous] <- original_steps[ambiguous]
  values[1L] + c(0, cumsum(steps))
}
