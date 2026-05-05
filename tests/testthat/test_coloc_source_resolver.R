# tests/testthat/test_coloc_source_resolver.R

skip_if_not_installed("Rsamtools")

mk_snps <- function() {
  tibble::tibble(
    SNP = "rs1", id.exposure = "x",
    effect_allele.exposure = "A", other_allele.exposure = "G",
    beta.exposure = 0.05, se.exposure = 0.01, pval.exposure = 1e-10
  )
}

test_that("nothing supplied -> error", {
  expect_error(
    ardmr:::.resolve_coloc_source(ancestry = "EUR"),
    "Provide one of"
  )
})

test_that("id + sumstats -> error", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  expect_error(
    ardmr:::.resolve_coloc_source(
      exposure_id = "ieu-b-38", exposure_sumstats = vcf, ancestry = "EUR"
    ),
    "both provided"
  )
})

test_that("snps_only without acknowledge_no_coloc -> error", {
  expect_error(
    ardmr:::.resolve_coloc_source(exposure_snps = mk_snps(), ancestry = "EUR"),
    "acknowledge_no_coloc"
  )
})

test_that("snps_only with acknowledge -> mode='snps_only', null fetcher", {
  r <- ardmr:::.resolve_coloc_source(
    exposure_snps = mk_snps(), ancestry = "EUR",
    acknowledge_no_coloc = TRUE, verbose = FALSE
  )
  expect_equal(r$mode, "snps_only")
  expect_null(r$coloc_fetcher)
  expect_null(r$coloc_metadata)
  expect_equal(nrow(r$exposure_snps), 1)
})

test_that("snps + vcf -> mode='snps_vcf' with VCF fetcher", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  r <- ardmr:::.resolve_coloc_source(
    exposure_snps = mk_snps(), exposure_sumstats = vcf,
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(r$mode, "snps_vcf")
  expect_true(is.function(r$coloc_fetcher))
  expect_true(is.list(r$coloc_metadata))
  out <- r$coloc_fetcher(1, 1L, 50000L)
  expect_true(nrow(out) > 0)
})

test_that("snps_vcf with GRCh38 declared -> error", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  expect_error(
    ardmr:::.resolve_coloc_source(
      exposure_snps = mk_snps(), exposure_sumstats = vcf,
      ancestry = "EUR", genome_build = "GRCh38"
    ),
    "GRCh37"
  )
})

test_that("vcf_only -> mode='vcf_only', IVs derived by clumping (mocked)", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  fake_ivs <- tibble::tibble(
    SNP = c("rs1","rs2"), Chr = c(1L, 2L), Pos = c(1000L, 2000L),
    effect_allele.exposure = c("A","C"), other_allele.exposure = c("G","T"),
    beta.exposure = c(0.1, 0.05), se.exposure = c(0.01, 0.01),
    eaf.exposure = c(0.3, 0.4), pval.exposure = c(1e-10, 1e-9),
    samplesize.exposure = c(100000, 100000),
    id.exposure = "mini.bgz", exposure = "mini.bgz"
  )
  testthat::local_mocked_bindings(
    clump_sumstats_to_ivs = function(...) fake_ivs
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_sumstats = vcf, ancestry = "EUR", verbose = FALSE
  )
  expect_equal(r$mode, "vcf_only")
  expect_equal(nrow(r$exposure_snps), 2)
})

test_that("vcf_only -> error when 0 IVs returned", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  testthat::local_mocked_bindings(
    clump_sumstats_to_ivs = function(...) ardmr:::.empty_exposure_tbl()
  )
  expect_error(
    ardmr:::.resolve_coloc_source(
      exposure_sumstats = vcf, ancestry = "EUR", verbose = FALSE
    ),
    "0 IVs"
  )
})

test_that("snps + id -> mode='snps_id' with ieugwasr fetcher", {
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) {
      list(N = 100000, type = "quant", ncase = NA_real_, ncontrol = NA_real_, s = NA_real_)
    }
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_snps = mk_snps(), exposure_id = "ieu-b-38",
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(r$mode, "snps_id")
  expect_true(is.function(r$coloc_fetcher))
  expect_equal(r$coloc_metadata$N, 100000)
})

test_that("ieugwasr (id only) -> extracts instruments + sets up coloc", {
  fake_ivs <- mk_snps()
  fake_ivs$beta.exposure <- 0.1; fake_ivs$se.exposure <- 0.01
  testthat::local_mocked_bindings(
    extract_instruments = function(outcomes, p1, clump, r2, kb) fake_ivs,
    .package = "TwoSampleMR"
  )
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) {
      list(N = 100000, type = "quant", ncase = NA_real_, ncontrol = NA_real_, s = NA_real_)
    }
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_id = "ieu-b-38", ancestry = "EUR", verbose = FALSE
  )
  expect_equal(r$mode, "ieugwasr")
  expect_true(is.function(r$coloc_fetcher))
  expect_equal(nrow(r$exposure_snps), 1)
})

test_that("ieugwasr -> errors when extract_instruments returns nothing", {
  testthat::local_mocked_bindings(
    extract_instruments = function(outcomes, p1, clump, r2, kb) NULL,
    .package = "TwoSampleMR"
  )
  expect_error(
    ardmr:::.resolve_coloc_source(
      exposure_id = "ieu-b-38", ancestry = "EUR", verbose = FALSE
    ),
    "no IVs"
  )
})

test_that("clump_opts overrides default p_threshold/r2/kb/f_threshold", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  observed <- list()
  fake_ivs <- tibble::tibble(
    SNP = "rs1", Chr = 1L, Pos = 1000L,
    effect_allele.exposure = "A", other_allele.exposure = "G",
    beta.exposure = 0.1, se.exposure = 0.01,
    eaf.exposure = 0.3, pval.exposure = 1e-10,
    samplesize.exposure = 100000,
    id.exposure = "mini.bgz", exposure = "mini.bgz"
  )
  testthat::local_mocked_bindings(
    clump_sumstats_to_ivs = function(vcf_path, ancestry, p_threshold, r2, kb,
                                     f_threshold, genome_build, cache_dir, verbose) {
      observed <<- list(p_threshold = p_threshold, r2 = r2, kb = kb, f_threshold = f_threshold)
      fake_ivs
    }
  )
  ardmr:::.resolve_coloc_source(
    exposure_sumstats = vcf, ancestry = "EUR",
    clump_opts = list(p_threshold = 1e-6, r2 = 0.01, kb = 5000L, f_threshold = 5),
    verbose = FALSE
  )
  expect_equal(observed$p_threshold, 1e-6)
  expect_equal(observed$r2, 0.01)
  expect_equal(observed$kb, 5000L)
  expect_equal(observed$f_threshold, 5)
})
