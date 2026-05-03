# tests/testthat/test_utils_paths.R
#
# Direct-assertion tests for the cache/path utilities in R/utils-paths.R:
#   * ardmr_cache_dir()
#   * set_cache_dir()
#   * variant_manifest_path()
#
# All env-var manipulation is wrapped in withr::local_envvar(...) so the
# user's real ARDMR_CACHE_DIR setting is never mutated (each call's
# scope ends with the test_that() block).

skip_if_not_installed("withr")

# ---- ardmr_cache_dir ---------------------------------------------------------

test_that("ardmr_cache_dir() errors when ARDMR_CACHE_DIR is unset", {
  withr::local_envvar(ARDMR_CACHE_DIR = NA)   # NA = unset
  expect_error(ardmr_cache_dir(), "ARDMR_CACHE_DIR is not set")
})

test_that("ardmr_cache_dir() errors when ARDMR_CACHE_DIR is empty string", {
  withr::local_envvar(ARDMR_CACHE_DIR = "")
  expect_error(ardmr_cache_dir(), "ARDMR_CACHE_DIR is not set")
})

test_that("ardmr_cache_dir() returns the value when ARDMR_CACHE_DIR is set", {
  withr::local_envvar(ARDMR_CACHE_DIR = "/tmp/some/cache")
  expect_equal(ardmr_cache_dir(), "/tmp/some/cache")
})

# ---- set_cache_dir -----------------------------------------------------------

test_that("set_cache_dir() errors on missing or empty path", {
  expect_error(set_cache_dir(), "non-empty directory path")
  expect_error(set_cache_dir(""), "non-empty directory path")
})

test_that("set_cache_dir() errors when the directory does not exist", {
  fake <- file.path(tempdir(), "definitely_does_not_exist_xyz_123")
  if (dir.exists(fake)) unlink(fake, recursive = TRUE)
  expect_error(set_cache_dir(fake), "Cache directory does not exist")
})

test_that("set_cache_dir() round-trips a valid directory into ARDMR_CACHE_DIR", {
  withr::local_envvar(ARDMR_CACHE_DIR = NA)
  td <- tempdir()
  ret <- set_cache_dir(td)
  expect_equal(Sys.getenv("ARDMR_CACHE_DIR"), td)
  # Helper returns its input invisibly per docstring.
  expect_equal(ret, td)
})

# ---- variant_manifest_path ---------------------------------------------------

test_that("variant_manifest_path() builds the Neale path", {
  expect_equal(
    variant_manifest_path("neale", "/x"),
    file.path("/x", "neale", "variants.tsv.bgz")
  )
})

test_that("variant_manifest_path() builds the Pan-UKB path", {
  expect_equal(
    variant_manifest_path("panukb", "/x"),
    file.path("/x", "panukb", "full_variant_qc_metrics.txt.bgz")
  )
})

test_that("variant_manifest_path() rejects invalid catalog values", {
  expect_error(variant_manifest_path("invalid", "/x"))
})

test_that("variant_manifest_path() defaults cache_dir to ardmr_cache_dir()", {
  withr::local_envvar(ARDMR_CACHE_DIR = "/from/env")
  expect_equal(
    variant_manifest_path("neale"),
    file.path("/from/env", "neale", "variants.tsv.bgz")
  )
})
