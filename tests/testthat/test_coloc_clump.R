# tests/testthat/test_coloc_clump.R
#
# Unit tests for clump_sumstats_to_ivs(). Skipped when gwasvcf is not
# installed (it is a Suggests dep used only by mode 'vcf_only').

skip_if_not_installed("gwasvcf")
skip_if_not_installed("Rsamtools")

mk_hits_tbl <- function() {
  tibble::tibble(
    rsid = paste0("rs", 1:5),
    seqnames = rep(1L, 5),
    start = c(1000L, 2000L, 3000L, 4000L, 5000L),
    REF = c("A","C","G","T","A"),
    ALT = c("G","T","A","C","G"),
    ES = c(0.1, 0.05, 0.005, 0.08, 0.001),
    SE = c(0.01, 0.01, 0.01, 0.01, 0.01),
    LP = c(20, 8, 0.5, 12, 0.1),
    AF = rep(0.3, 5),
    SS = rep(100000L, 5)
  )
}

mk_fake_ld_ref <- function() {
  list(bfile_prefix = file.path(tempdir(), "fake_bfile"),
       present = TRUE, downloaded = FALSE, aborted = FALSE)
}

test_that("p_threshold filtering applied before clumping", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  observed <- list()
  testthat::local_mocked_bindings(
    query_gwas = function(vcffile, pval, ...) { observed$pval <<- pval; structure(list()) },
    vcf_to_tibble = function(vcf) mk_hits_tbl(),
    .package = "gwasvcf"
  )
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) mk_fake_ld_ref()
  )
  testthat::local_mocked_bindings(
    ld_clump = function(dat, ...) dat,
    .package = "ieugwasr"
  )
  ivs <- clump_sumstats_to_ivs(vcf, ancestry = "EUR", p_threshold = 1e-7,
                               f_threshold = 0, verbose = FALSE)
  expect_equal(observed$pval, 1e-7)
})

test_that("clump params (r2, kb) plumbed to ld_clump", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  observed <- list()
  testthat::local_mocked_bindings(
    query_gwas = function(vcffile, pval, ...) structure(list()),
    vcf_to_tibble = function(vcf) mk_hits_tbl(),
    .package = "gwasvcf"
  )
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) mk_fake_ld_ref()
  )
  testthat::local_mocked_bindings(
    ld_clump = function(dat, clump_r2, clump_kb, ...) {
      observed$r2 <<- clump_r2; observed$kb <<- clump_kb; dat
    },
    .package = "ieugwasr"
  )
  clump_sumstats_to_ivs(vcf, ancestry = "EUR", r2 = 0.005, kb = 1234L,
                       f_threshold = 0, verbose = FALSE)
  expect_equal(observed$r2, 0.005)
  expect_equal(observed$kb, 1234L)
})

test_that("F-statistic threshold filters weak instruments", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  testthat::local_mocked_bindings(
    query_gwas = function(vcffile, pval, ...) structure(list()),
    vcf_to_tibble = function(vcf) mk_hits_tbl(),
    .package = "gwasvcf"
  )
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) mk_fake_ld_ref()
  )
  testthat::local_mocked_bindings(
    ld_clump = function(dat, ...) dat,
    .package = "ieugwasr"
  )
  ivs <- clump_sumstats_to_ivs(vcf, ancestry = "EUR",
                               p_threshold = 1, f_threshold = 50,
                               verbose = FALSE)
  Fstat <- (ivs$beta.exposure / ivs$se.exposure) ^ 2
  expect_true(all(Fstat >= 50))
})

test_that("empty hits returns empty TwoSampleMR-formatted tibble (no error)", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  testthat::local_mocked_bindings(
    query_gwas = function(vcffile, pval, ...) structure(list()),
    vcf_to_tibble = function(vcf) mk_hits_tbl()[0, ],
    .package = "gwasvcf"
  )
  ivs <- clump_sumstats_to_ivs(vcf, ancestry = "EUR", verbose = FALSE)
  expect_equal(nrow(ivs), 0)
  expect_true(all(c("SNP","Chr","Pos","effect_allele.exposure",
                    "other_allele.exposure","beta.exposure","se.exposure",
                    "id.exposure") %in% names(ivs)))
})

test_that("output passes assert_exposure() column requirements", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  testthat::local_mocked_bindings(
    query_gwas = function(vcffile, pval, ...) structure(list()),
    vcf_to_tibble = function(vcf) mk_hits_tbl(),
    .package = "gwasvcf"
  )
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) mk_fake_ld_ref()
  )
  testthat::local_mocked_bindings(
    ld_clump = function(dat, ...) dat,
    .package = "ieugwasr"
  )
  ivs <- clump_sumstats_to_ivs(vcf, ancestry = "EUR", f_threshold = 0,
                               verbose = FALSE)
  expect_silent(assert_exposure(ivs))
})

test_that("missing tabix index errors at validate stage", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  unlink(paste0(vcf, ".tbi"))
  expect_error(
    clump_sumstats_to_ivs(vcf, ancestry = "EUR", verbose = FALSE),
    "Tabix index missing"
  )
})
