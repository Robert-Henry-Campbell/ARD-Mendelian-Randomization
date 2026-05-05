# tests/testthat/test_coloc_data_fetch.R
# Caching + metadata-derivation invariants for the coloc data fetchers.

skip_if_not_installed("ieugwasr")

cache_dir_test <- function() {
  d <- file.path(tempdir(), paste0("ardmr_coloc_fetch_test_", as.integer(Sys.time())))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

test_that("coloc_panukb_outcome_metadata: case-control phenocode", {
  rec <- tibble::tibble(
    trait_type = "icd10",
    n_cases_EUR = 5000, n_controls_EUR = 400000
  )
  m <- coloc_panukb_outcome_metadata(rec, ancestry = "EUR")
  expect_equal(m$type, "cc")
  expect_equal(m$N, 405000)
  expect_equal(m$ncase, 5000)
  expect_equal(m$ncontrol, 400000)
  expect_equal(m$s, 5000 / 405000)
})

test_that("coloc_panukb_outcome_metadata: continuous trait uses n_cases as N (quant)", {
  rec <- tibble::tibble(
    trait_type = "continuous",
    n_cases_EUR = 173082, n_controls_EUR = NA_real_
  )
  m <- coloc_panukb_outcome_metadata(rec, ancestry = "EUR")
  expect_equal(m$type, "quant")
  expect_equal(m$N, 173082)
  expect_true(is.na(m$ncase))
  expect_true(is.na(m$ncontrol))
  expect_true(is.na(m$s))
})

test_that("coloc_panukb_outcome_metadata: biomarkers behaves like continuous", {
  rec <- tibble::tibble(
    trait_type = "biomarkers",
    n_cases_EUR = 250000, n_controls_EUR = NA_real_
  )
  m <- coloc_panukb_outcome_metadata(rec, ancestry = "EUR")
  expect_equal(m$type, "quant")
  expect_equal(m$N, 250000)
})

test_that("coloc_fetch_metadata does NOT cache when ieugwasr::gwasinfo fails", {
  cd <- cache_dir_test()
  testthat::local_mocked_bindings(
    gwasinfo = function(id) stop("simulated API failure"),
    .package = "ieugwasr"
  )
  m <- coloc_fetch_metadata("ieu-test-fail", cache_dir = cd, verbose = FALSE)
  expect_true(is.na(m$N))
  cache_files <- list.files(file.path(cd, "coloc", "meta"), full.names = TRUE)
  expect_length(cache_files, 0)
})

test_that("coloc_fetch_metadata does cache when N is finite", {
  cd <- cache_dir_test()
  testthat::local_mocked_bindings(
    gwasinfo = function(id) {
      tibble::tibble(id = id, sample_size = 100000, ncase = NA_real_,
                     ncontrol = NA_real_, unit = "SD")
    },
    .package = "ieugwasr"
  )
  m <- coloc_fetch_metadata("ieu-test-good", cache_dir = cd, verbose = FALSE)
  expect_equal(m$N, 100000)
  cache_files <- list.files(file.path(cd, "coloc", "meta"), full.names = TRUE)
  expect_length(cache_files, 1)
})

test_that("coloc_fetch_exposure_region does NOT cache when API call fails", {
  cd <- cache_dir_test()
  testthat::local_mocked_bindings(
    associations = function(variants, id, ...) stop("simulated API failure"),
    .package = "ieugwasr"
  )
  out <- coloc_fetch_exposure_region("ieu-test-x", chr = 1L, start = 1L, end = 1000L,
                                     cache_dir = cd, verbose = FALSE)
  expect_equal(nrow(out), 0)
  cache_files <- list.files(file.path(cd, "coloc", "exposure"), full.names = TRUE)
  expect_length(cache_files, 0)
})

test_that("coloc_fetch_exposure_region caches non-empty success and reuses", {
  cd <- cache_dir_test()
  call_count <- 0L
  testthat::local_mocked_bindings(
    associations = function(variants, id, ...) {
      call_count <<- call_count + 1L
      tibble::tibble(
        rsid = "rs1", chr = 1L, position = 100L,
        ea = "A", nea = "G", beta = 0.05, se = 0.01, eaf = 0.3, p = 1e-8
      )
    },
    .package = "ieugwasr"
  )
  out1 <- coloc_fetch_exposure_region("ieu-test-y", chr = 1L, start = 1L, end = 200L,
                                      cache_dir = cd, verbose = FALSE)
  out2 <- coloc_fetch_exposure_region("ieu-test-y", chr = 1L, start = 1L, end = 200L,
                                      cache_dir = cd, verbose = FALSE)
  expect_equal(nrow(out1), 1)
  expect_equal(nrow(out2), 1)
  expect_equal(call_count, 1L)  # second call hit cache
})
