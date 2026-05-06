# tests/testthat/test_cache_keys.R

mk_iv <- function(snps = c("rs1","rs2","rs3"),
                  chr = c(1L, 2L, 3L),
                  pos = c(100L, 200L, 300L),
                  ea = c("A","C","G"), oa = c("G","T","A"),
                  cols = c("Chr","Pos")) {
  out <- data.frame(SNP = snps, ea = ea, oa = oa, stringsAsFactors = FALSE)
  out[[cols[1]]] <- chr
  out[[cols[2]]] <- pos
  names(out)[names(out) == "ea"] <- "effect_allele.exposure"
  names(out)[names(out) == "oa"] <- "other_allele.exposure"
  out
}

test_that(".iv_set_hash is order-invariant", {
  h1 <- ardmr:::.iv_set_hash(mk_iv())
  h2 <- ardmr:::.iv_set_hash(mk_iv(snps = c("rs3","rs1","rs2"),
                                   chr = c(3L, 1L, 2L),
                                   pos = c(300L, 100L, 200L),
                                   ea  = c("G","A","C"),
                                   oa  = c("A","G","T")))
  expect_equal(h1, h2)
})

test_that(".iv_set_hash differs when chr/pos changes", {
  h1 <- ardmr:::.iv_set_hash(mk_iv())
  h2 <- ardmr:::.iv_set_hash(mk_iv(pos = c(101L, 200L, 300L)))
  expect_false(identical(h1, h2))
})

test_that(".iv_set_hash differs when alleles change", {
  h1 <- ardmr:::.iv_set_hash(mk_iv())
  h2 <- ardmr:::.iv_set_hash(mk_iv(ea = c("T","C","G")))
  expect_false(identical(h1, h2))
})

test_that(".iv_set_hash works with panukb_chrom/panukb_pos columns", {
  h_default <- ardmr:::.iv_set_hash(mk_iv())
  h_panukb  <- ardmr:::.iv_set_hash(mk_iv(cols = c("panukb_chrom","panukb_pos")))
  expect_equal(h_default, h_panukb)
})

test_that(".iv_set_hash works with neale_chrom/neale_pos columns", {
  h_default <- ardmr:::.iv_set_hash(mk_iv())
  h_neale   <- ardmr:::.iv_set_hash(mk_iv(cols = c("neale_chrom","neale_pos")))
  expect_equal(h_default, h_neale)
})

test_that(".iv_set_hash returns short hex of configurable length", {
  h <- ardmr:::.iv_set_hash(mk_iv(), n = 6L)
  expect_equal(nchar(h), 6L)
  expect_match(h, "^[0-9a-f]+$")
})

test_that(".iv_set_hash handles empty input safely", {
  h <- ardmr:::.iv_set_hash(NULL)
  expect_equal(nchar(h), 5L)
  h2 <- ardmr:::.iv_set_hash(data.frame())
  expect_equal(nchar(h2), 5L)
})

test_that(".run_input_hash is deterministic for identical inputs", {
  args <- list(
    exposure_snps = mk_iv(),
    exposure_id = "ieu-a-1",
    sensitivity_enabled = c("a","b","c"),
    sensitivity_pass_min = 6L,
    clump_opts = list(p_threshold = 5e-8, r2 = 0.001),
    coloc_window_kb = 500L,
    coloc_priors = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    coloc_skip_mhc = TRUE,
    multiple_testing_correction = "BH"
  )
  h1 <- do.call(ardmr:::.run_input_hash, args)
  h2 <- do.call(ardmr:::.run_input_hash, args)
  expect_equal(h1, h2)
})

test_that(".run_input_hash is sensitive to clump_opts changes", {
  base <- list(exposure_snps = mk_iv(), exposure_id = "x",
               sensitivity_enabled = c("a","b"))
  h1 <- do.call(ardmr:::.run_input_hash,
                c(base, list(clump_opts = list(p_threshold = 5e-8))))
  h2 <- do.call(ardmr:::.run_input_hash,
                c(base, list(clump_opts = list(p_threshold = 1e-6))))
  expect_false(identical(h1, h2))
})

test_that(".run_input_hash is order-invariant for sensitivity_enabled set", {
  base <- list(exposure_snps = mk_iv(), exposure_id = "x")
  h1 <- do.call(ardmr:::.run_input_hash, c(base, list(sensitivity_enabled = c("a","b","c"))))
  h2 <- do.call(ardmr:::.run_input_hash, c(base, list(sensitivity_enabled = c("c","a","b"))))
  expect_equal(h1, h2)
})

test_that(".run_input_hash differs when exposure_id changes", {
  h1 <- ardmr:::.run_input_hash(exposure_snps = mk_iv(), exposure_id = "ieu-a-1")
  h2 <- ardmr:::.run_input_hash(exposure_snps = mk_iv(), exposure_id = "ieu-a-2")
  expect_false(identical(h1, h2))
})

test_that(".run_input_hash is NULL-safe", {
  h <- ardmr:::.run_input_hash()
  expect_match(h, "^[0-9a-f]{5}__[0-9a-f]{5}$")
})

test_that("ardmr_run_dir composes the expected path", {
  p <- ardmr_run_dir("/cache", "adhd", "both", "EUR", "abc1234567")
  expect_equal(p, "/cache/output/adhd/both/EUR/abc1234567")
})

test_that(".run_input_hash prefix matches .iv_set_hash for the same IVs", {
  iv <- mk_iv()
  rh <- ardmr:::.run_input_hash(exposure_snps = iv, exposure_id = "x")
  expect_equal(substr(rh, 1L, 5L), ardmr:::.iv_set_hash(iv))
})
