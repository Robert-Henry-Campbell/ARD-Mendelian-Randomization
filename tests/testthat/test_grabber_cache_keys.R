# tests/testthat/test_grabber_cache_keys.R

mk_panukb_exposure <- function(snp = "rs1", chr = 1L, pos = 100L,
                               ea = "A", oa = "G") {
  tibble::tibble(
    rsid = snp, panukb_chrom = chr, panukb_pos = pos,
    effect_allele.exposure = ea, other_allele.exposure = oa
  )
}

mk_panukb_mrdf <- function(label = "synth_phenotype") {
  tibble::tibble(
    description    = label,
    aws_link       = "https://invalid.example/synth.tsv.bgz",
    aws_link_tabix = "https://invalid.example/synth.tsv.bgz.tbi"
  )
}

setup_cache <- function() {
  d <- tempfile("ardmr_cache_")
  dir.create(d, recursive = TRUE)
  d
}

test_that("panukb_snp_grabber cache-hit reads new __<iv_hash>.rds filename", {
  skip_if_not_installed("Rsamtools")
  cache_dir <- setup_cache()
  cache_root <- file.path(cache_dir, "panukb_outcome_snps", "EUR")
  dir.create(cache_root, recursive = TRUE)
  exp <- mk_panukb_exposure()
  iv_hash <- ardmr:::.iv_set_hash(exp)
  payload <- tibble::tibble(SNP = "rs1", chr.outcome = 1L, pos.outcome = 100L,
                            effect_allele.outcome = "A", other_allele.outcome = "G",
                            beta.outcome = 0.1, se.outcome = 0.01)
  saveRDS(payload, file.path(cache_root,
                             paste0("synth_phenotype__", iv_hash, ".rds")))
  out <- panukb_snp_grabber(exp, mk_panukb_mrdf(), ancestry = "EUR",
                            cache_dir = cache_dir, verbose = FALSE,
                            force_refresh = FALSE)
  expect_equal(nrow(out$outcome_snps[[1]]), 1L)
  expect_equal(out$outcome_snps[[1]]$SNP, "rs1")
})

test_that("panukb_snp_grabber: two IV sets get separate cache files", {
  skip_if_not_installed("Rsamtools")
  cache_dir <- setup_cache()
  cache_root <- file.path(cache_dir, "panukb_outcome_snps", "EUR")
  dir.create(cache_root, recursive = TRUE)
  expA <- mk_panukb_exposure(snp = "rs1", chr = 1L, pos = 100L)
  expB <- mk_panukb_exposure(snp = "rs2", chr = 2L, pos = 200L)
  hA <- ardmr:::.iv_set_hash(expA)
  hB <- ardmr:::.iv_set_hash(expB)
  expect_false(identical(hA, hB))
  saveRDS(tibble::tibble(SNP = "rs1", chr.outcome=1L, pos.outcome=100L,
                          effect_allele.outcome="A", other_allele.outcome="G",
                          beta.outcome = 0.10, se.outcome = 0.01),
          file.path(cache_root, paste0("synth_phenotype__", hA, ".rds")))
  saveRDS(tibble::tibble(SNP = "rs2", chr.outcome=2L, pos.outcome=200L,
                          effect_allele.outcome="A", other_allele.outcome="G",
                          beta.outcome = 0.05, se.outcome = 0.01),
          file.path(cache_root, paste0("synth_phenotype__", hB, ".rds")))
  outA <- panukb_snp_grabber(expA, mk_panukb_mrdf(), ancestry = "EUR",
                             cache_dir = cache_dir, verbose = FALSE,
                             force_refresh = FALSE)
  outB <- panukb_snp_grabber(expB, mk_panukb_mrdf(), ancestry = "EUR",
                             cache_dir = cache_dir, verbose = FALSE,
                             force_refresh = FALSE)
  expect_equal(outA$outcome_snps[[1]]$SNP, "rs1")
  expect_equal(outB$outcome_snps[[1]]$SNP, "rs2")
  files <- list.files(cache_root)
  expect_equal(sum(grepl("synth_phenotype__", files)), 2L)
})

test_that("panukb_snp_grabber: force_refresh deletes only the matching iv_hash", {
  skip_if_not_installed("Rsamtools")
  cache_dir <- setup_cache()
  cache_root <- file.path(cache_dir, "panukb_outcome_snps", "EUR")
  dir.create(cache_root, recursive = TRUE)
  exp <- mk_panukb_exposure()
  iv_hash <- ardmr:::.iv_set_hash(exp)
  matching <- file.path(cache_root, paste0("synth_phenotype__", iv_hash, ".rds"))
  sibling  <- file.path(cache_root, "synth_phenotype__deadbeef00.rds")
  saveRDS(tibble::tibble(SNP = "rs1"), matching)
  saveRDS(tibble::tibble(SNP = "rs2"), sibling)
  suppressMessages(suppressWarnings(
    panukb_snp_grabber(exp, mk_panukb_mrdf(), ancestry = "EUR",
                       cache_dir = cache_dir, verbose = FALSE,
                       force_refresh = TRUE)
  ))
  expect_true(file.exists(sibling))
  expect_equal(readRDS(sibling)$SNP, "rs2")
  expect_true(!file.exists(matching) || nrow(readRDS(matching)) == 0L)
})

mk_neale_exposure <- function(snp = "rs1", chr = 1L, pos = 100L) {
  tibble::tibble(
    rsid = snp, neale_chrom = chr, neale_pos = pos,
    neale_variant = paste0(chr, ":", pos, ":A:G"),
    effect_allele.exposure = "A", other_allele.exposure = "G"
  )
}

test_that("neale_snp_grabber: two IV sets get separate cache files", {
  skip_if_not_installed("Rsamtools")
  skip_if(!exists("neale_snp_grabber"))
  cache_dir <- setup_cache()
  cache_root <- file.path(cache_dir, "neale_outcome_snps", "unknown")
  dir.create(cache_root, recursive = TRUE)
  expA <- mk_neale_exposure(snp = "rs1", chr = 1L, pos = 100L)
  expB <- mk_neale_exposure(snp = "rs2", chr = 2L, pos = 200L)
  hA <- ardmr:::.iv_set_hash(expA)
  hB <- ardmr:::.iv_set_hash(expB)
  expect_false(identical(hA, hB))
  expect_match(hA, "^[0-9a-f]{10}$")
})
