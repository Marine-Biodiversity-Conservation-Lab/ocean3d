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

test_that("midpoint edges keep the order of the centres", {
  regular <- .midpoint_edges(c(30, 20, 10))
  expect_equal(regular$lower, c(25, 15, 5))
  expect_equal(regular$upper, c(35, 25, 15))
  expect_true(regular$regular)

  irregular <- .midpoint_edges(c(5, 15, 35))
  expect_equal(irregular$lower, c(0, 10, 25))
  expect_equal(irregular$upper, c(10, 25, 45))
  expect_false(irregular$regular)
})
