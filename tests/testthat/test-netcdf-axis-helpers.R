# Shared axis helpers used by buffer extraction and the core object (#113),
# tested directly (AGENTS.md, "Tests").

test_that("depth bounds are converted to metres positive down", {
  # Centres at 5 m and 15 m, stored as heights in centimetres.
  f <- make_grid_file(
    depth = c(-500, -1500), depth_units = "cm", positive = "up",
    depth_bounds = rbind(c(0, -1000), c(-1000, -2000)), time = NULL
  )
  on.exit(unlink(f), add = TRUE)
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  bounds <- .read_bounds(nc, "depth", 2L, strict = TRUE)
  expect_equal(
    .depth_bounds_m(nc, "depth", c(-500, -1500), bounds),
    rbind(c(0, 10), c(10, 20))
  )
})
