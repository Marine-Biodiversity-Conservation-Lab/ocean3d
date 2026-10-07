# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("fetch_species_assessments errors when no input provided", {
  skip_if_not_installed("rredlist")
  expect_error(
    fetch_species_assessments(api_key = Sys.getenv("IUCN_REDLIST_KEY")),
    "Exactly one"
  )
})

test_that("fetch_species_assessments errors when multiple inputs provided", {
  skip_if_not_installed("rredlist")
  expect_error(
    fetch_species_assessments(
      api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
      sis_ids       = 1,
      species_names = "Sphyrna lewini"
    ),
    "Exactly one"
  )
})

test_that("fetch_species_assessments errors when all three inputs provided", {
  skip_if_not_installed("rredlist")
  expect_error(
    fetch_species_assessments(
      api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
      sis_ids       = 1,
      species_names = "Sphyrna lewini",
      group_code    = "sharks_and_rays"
    ),
    "Exactly one"
  )
})

# ---------------------------------------------------------------------------
# group_code path
# ---------------------------------------------------------------------------

# Time consuming test, skip usually
test_that("fetch_species_assessments returns data frame via group_code", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  skip()
  skip_on_ci()
  skip_on_cran()
  vcr::local_cassette("group_code")
  out <- fetch_species_assessments(
    api_key    = Sys.getenv("IUCN_REDLIST_KEY"),
    group_code = "sharks_and_rays"
  )
  expect_s3_class(out, "data.frame")
  expect_gt(nrow(out), 0L)
})

test_that("fetch_species_assessments result has expected columns", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("expected_columns")
  out <- fetch_species_assessments(
    api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
    species_names = "Sphyrna lewini"
  )
  expect_named(out, c(
    "assessment_id", "assessment_date", "sis_id", "scientific_name",
    "kingdom_name", "phylum_name", "class_name", "order_name",
    "family_name", "genus_name", "species_name", "subpopulation_name",
    "red_list_category", "systems_code",
    "upper_depth_limit", "lower_depth_limit", "citation", "url"
  ))
})

# ---------------------------------------------------------------------------
# sis_ids path
# ---------------------------------------------------------------------------

test_that("fetch_species_assessments resolves sis_ids to assessments", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("sis_ids")
  out <- fetch_species_assessments(
    api_key = Sys.getenv("IUCN_REDLIST_KEY"),
    sis_ids = c(44584L, 60191L)
  )
  expect_identical(nrow(out), 2L)
})

test_that("fetch_species_assessments warns and skips bad SIS ID", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("bad_sis_ids")
  expect_warning(
    out <- fetch_species_assessments(
      api_key = Sys.getenv("IUCN_REDLIST_KEY"),
      sis_ids = c(1L, 44584L)
    ),
    "1"
  )
  expect_identical(nrow(out), 1L)
})

test_that("fetch_species_assessments errors when all SIS IDs fail", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("all_sis_ids_fail")
  expect_warning(
    expect_error(
      fetch_species_assessments(
        api_key = Sys.getenv("IUCN_REDLIST_KEY"),
        sis_ids = 99999999L
      ),
      "No Global-scope assessments"
    )
  )
})

# ---------------------------------------------------------------------------
# species_names path
# ---------------------------------------------------------------------------

test_that("fetch_species_assessments resolves species names to assessments", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("species_names")
  out <- fetch_species_assessments(
    api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
    species_names = "Sphyrna lewini"
  )
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 1L)
})

test_that("fetch_species_assessments warns on unparseable species name", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("unparseable_species_name")
  expect_warning(
    out <- fetch_species_assessments(
      api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
      species_names = c("Sphyrna lewini", "Sphyrna")  # second has no epithet
    ),
    "Sphyrna"
  )
  expect_identical(nrow(out), 1L)
})

test_that("fetch_species_assessments warns and skips species name not in API", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("species_not_in_api")
  expect_warning(
    out <- fetch_species_assessments(
      api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
      species_names = c("Sphyrna lewini", "Foo unknownus")
    ),
    "Foo unknownus"
  )
  expect_identical(nrow(out), 1L)
})

# ---------------------------------------------------------------------------
# Output field values
# ---------------------------------------------------------------------------

test_that("fetch_species_assessments populates depth limits from assessment", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("depth_limits")
  out <- fetch_species_assessments(
    api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
    species_names = "Sphyrna lewini"
  )
  expect_false(is.na(out$upper_depth_limit))
  expect_false(is.na(out$lower_depth_limit))
  expect_gte(out$lower_depth_limit, out$upper_depth_limit)
})

test_that("fetch_species_assessments populates systems_code from assessment", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("systems_code")
  out <- fetch_species_assessments(
    api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
    species_names = "Sphyrna lewini"
  )
  expect_type(out$systems_code, "character")
  expect_false(is.na(out$systems_code))
})

test_that("fetch_species_assessments depth fields are numeric or NA", {
  skip_if_not_installed("rredlist")
  skip_if_not_installed("vcr")
  vcr::local_cassette("depth_fields_type")
  out <- fetch_species_assessments(
    api_key       = Sys.getenv("IUCN_REDLIST_KEY"),
    species_names = "Sphyrna lewini"
  )
  expect_true(is.numeric(out$upper_depth_limit) || is.na(out$upper_depth_limit))
  expect_true(is.numeric(out$lower_depth_limit) || is.na(out$lower_depth_limit))
})

# ---------------------------------------------------------------------------
# fill_missing_depths
# ---------------------------------------------------------------------------

test_that("fill_missing_depths swaps reversed upper and lower values", {
  out <- fill_missing_depths(
    upper = c(100, 50),
    lower = c(10, 500),  # first row reversed
    genus = c("Carcharodon", "Sphyrna")
  )
  expect_identical(out$upper_depth, c(10, 50))
  expect_identical(out$lower_depth, c(100, 500))
})

test_that("fill_missing_depths fills NAs with genus means", {
  out <- fill_missing_depths(
    upper = c(0, 10, NA),
    lower = c(100, 200, NA),
    genus = c("Carcharhinus", "Carcharhinus", "Carcharhinus")
  )
  # Mean of non-NA values in the genus: upper = 5, lower = 150
  expect_identical(out$upper_depth[3], 5)
  expect_identical(out$lower_depth[3], 150)
})

test_that("fill_missing_depths leaves NA when entire genus is NA", {
  out <- fill_missing_depths(
    upper = c(NA, NA),
    lower = c(NA, NA),
    genus = c("Raja", "Raja")
  )
  expect_true(all(is.na(out$upper_depth)))
  expect_true(all(is.na(out$lower_depth)))
})

test_that("fill_missing_depths rejects unsupported methods", {
  expect_error(
    fill_missing_depths(1, 10, "Carcharodon", method = "zero"),
    "genus_mean"
  )
})

test_that("fill_missing_depths returns a data frame with expected columns", {
  out <- fill_missing_depths(
    upper = c(0, 10),
    lower = c(100, 200),
    genus = c("A", "B")
  )
  expect_s3_class(out, "data.frame")
  expect_named(out, c("upper_depth", "lower_depth"))
  expect_identical(nrow(out), 2L)
})
