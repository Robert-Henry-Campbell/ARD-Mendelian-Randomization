# tests/testthat/test_ld_matrix.R
#
# Pragmatic, offline tests for the LD matrix module in R/ld_matrix.R and
# the population helper in R/utils-ld.R. PLINK + reference-panel-bound
# behaviour is validated indirectly via the single-SNP short-circuit and
# the cache-hit path; full end-to-end PLINK runs are intentionally not
# exercised on CI (binaries + ~1 GB reference panels are out of scope).

skip_if_not_installed("withr")

# ---- ld_pop_from_ancestry ----------------------------------------------------

test_that("ld_pop_from_ancestry() maps canonical labels to 1000G superpops", {
  expect_equal(ld_pop_from_ancestry("EUR"),          "EUR")
  expect_equal(ld_pop_from_ancestry("european"),     "EUR")
  expect_equal(ld_pop_from_ancestry("AFR"),          "AFR")
  expect_equal(ld_pop_from_ancestry("AFRICAN-UKB"),  "AFR")
  expect_equal(ld_pop_from_ancestry("EAS"),          "EAS")
  expect_equal(ld_pop_from_ancestry("east_asian"),   "EAS")
  expect_equal(ld_pop_from_ancestry("SAS"),          "SAS")
  expect_equal(ld_pop_from_ancestry("CSA"),          "SAS")
  expect_equal(ld_pop_from_ancestry("south_asian"),  "SAS")
})

test_that("ld_pop_from_ancestry() warns and defaults to EUR for unknown labels", {
  expect_warning(out <- ld_pop_from_ancestry("MARTIAN"), "Unknown ancestry")
  expect_equal(out, "EUR")
})

test_that("ld_pop_from_ancestry() handles NULL/NA without erroring", {
  expect_warning(p1 <- ld_pop_from_ancestry(NA), "Unknown ancestry")
  expect_equal(p1, "EUR")
  expect_warning(p2 <- ld_pop_from_ancestry(NULL), "Unknown ancestry")
  expect_equal(p2, "EUR")
})

# ---- .ld_cache_key -----------------------------------------------------------

test_that(".ld_cache_key() is invariant to SNP order and duplicates", {
  k1 <- .ld_cache_key(c("rs1","rs2","rs3"), "EUR")
  k2 <- .ld_cache_key(c("rs3","rs1","rs2"), "EUR")
  k3 <- .ld_cache_key(c("rs1","rs2","rs3","rs1"), "EUR")
  expect_equal(k1, k2)
  expect_equal(k1, k3)
})

test_that(".ld_cache_key() distinguishes population", {
  expect_false(identical(
    .ld_cache_key(c("rs1","rs2"), "EUR"),
    .ld_cache_key(c("rs1","rs2"), "AFR")
  ))
})

test_that(".ld_cache_key() distinguishes panel version", {
  expect_false(identical(
    .ld_cache_key(c("rs1","rs2"), "EUR", panel_version = "1kg.v3"),
    .ld_cache_key(c("rs1","rs2"), "EUR", panel_version = "1kg.v4")
  ))
})

test_that(".ld_cache_key() returns a 16-character hex string", {
  k <- .ld_cache_key(c("rs1","rs2"), "EUR")
  expect_true(is.character(k))
  expect_equal(nchar(k), 16L)
  expect_match(k, "^[0-9a-f]{16}$")
})

# ---- .harmonise_ld_signs -----------------------------------------------------

test_that(".harmonise_ld_signs() leaves matrix unchanged when all effect alleles match A1", {
  ld <- matrix(c(1.0, 0.5, 0.2,
                 0.5, 1.0, 0.3,
                 0.2, 0.3, 1.0), 3, 3)
  out <- .harmonise_ld_signs(
    ld_mat = ld,
    bim_alleles      = list(A1 = c("A","C","G"), A2 = c("T","G","A")),
    exposure_alleles = list(effect = c("A","C","G"), other = c("T","G","A")),
    snps = c("rs1","rs2","rs3")
  )
  expect_equal(out$matrix, ld)
  expect_equal(out$kept_idx, 1:3)
  expect_length(out$dropped_snps, 0)
})

test_that(".harmonise_ld_signs() flips sign on a single SNP whose effect allele equals A2", {
  ld <- matrix(c(1.0, 0.5, 0.2,
                 0.5, 1.0, 0.3,
                 0.2, 0.3, 1.0), 3, 3)
  # rs2 has effect/other swapped relative to bim: effect=G (A2), other=C (A1) -> flip rs2 row+col
  out <- .harmonise_ld_signs(
    ld_mat = ld,
    bim_alleles      = list(A1 = c("A","C","G"), A2 = c("T","G","A")),
    exposure_alleles = list(effect = c("A","G","G"), other = c("T","C","A")),
    snps = c("rs1","rs2","rs3")
  )
  expected <- ld
  expected[2, ] <- -expected[2, ]
  expected[, 2] <- -expected[, 2]
  # Diagonal: -1 * 1 * -1 = 1 (preserved); off-diagonals (1,2),(2,3) flipped once
  expect_equal(out$matrix, expected)
  expect_equal(out$matrix[2, 2], 1)
})

test_that(".harmonise_ld_signs() drops SNPs whose alleles match neither A1 nor A2", {
  ld <- matrix(c(1.0, 0.5, 0.2,
                 0.5, 1.0, 0.3,
                 0.2, 0.3, 1.0), 3, 3)
  out <- .harmonise_ld_signs(
    ld_mat = ld,
    bim_alleles      = list(A1 = c("A","C","G"), A2 = c("T","G","A")),
    # rs2 has alleles that don't match the bim at all
    exposure_alleles = list(effect = c("A","X","G"), other = c("T","Y","A")),
    snps = c("rs1","rs2","rs3")
  )
  expect_equal(out$kept_idx, c(1L, 3L))
  expect_equal(out$dropped_snps, "rs2")
  expect_equal(out$dropped_reasons, "allele_mismatch")
  expect_equal(dim(out$matrix), c(2L, 2L))
  expect_equal(out$matrix, ld[c(1, 3), c(1, 3)])
})

test_that(".harmonise_ld_signs() returns an empty matrix when all SNPs mismatch", {
  ld <- matrix(c(1.0, 0.5, 0.5, 1.0), 2, 2)
  out <- .harmonise_ld_signs(
    ld_mat = ld,
    bim_alleles      = list(A1 = c("A","C"), A2 = c("T","G")),
    exposure_alleles = list(effect = c("X","Y"), other = c("Z","W")),
    snps = c("rs1","rs2")
  )
  expect_equal(dim(out$matrix), c(0L, 0L))
  expect_equal(out$dropped_snps, c("rs1","rs2"))
  expect_length(out$kept_idx, 0)
})

# ---- .write_ld_csv -----------------------------------------------------------

test_that(".write_ld_csv() round-trips a small matrix", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  snps <- c("rs1","rs2","rs3","rs4")
  m <- matrix(c(1.0, 0.5, 0.2, 0.1,
                0.5, 1.0, 0.3, 0.0,
                0.2, 0.3, 1.0, 0.4,
                0.1, 0.0, 0.4, 1.0), 4, 4,
              dimnames = list(snps, snps))
  .write_ld_csv(m, snps, tmp)
  back <- utils::read.csv(tmp, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(colnames(back)[1], "SNP")
  expect_equal(back$SNP, snps)
  expect_equal(colnames(back)[-1], snps)
  m_back <- as.matrix(back[, -1])
  storage.mode(m_back) <- "double"
  dimnames(m_back) <- dimnames(m)
  expect_equal(m_back, m, tolerance = 1e-12)
})

test_that(".write_ld_csv() writes NA rows/cols for SNPs missing from the matrix and forces diagonal=1", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  m <- matrix(c(1.0, 0.5, 0.5, 1.0), 2, 2, dimnames = list(c("rs1","rs2"), c("rs1","rs2")))
  full_snps <- c("rs1","rs2","rs3")  # rs3 dropped from matrix
  .write_ld_csv(m, full_snps, tmp)
  back <- utils::read.csv(tmp, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(back$SNP, full_snps)
  expect_equal(back$rs1, c(1.0, 0.5, NA))
  expect_equal(back$rs2, c(0.5, 1.0, NA))
  # Diagonal for rs3 forced to 1 even though rs3 not in matrix
  expect_equal(back$rs3, c(NA, NA, 1.0))
})

test_that(".write_ld_csv() handles 0-SNP and 1-SNP cases", {
  tmp_empty <- withr::local_tempfile(fileext = ".csv")
  .write_ld_csv(matrix(numeric(0), 0, 0), character(0), tmp_empty)
  back0 <- utils::read.csv(tmp_empty, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(back0), 0L)
  expect_equal(colnames(back0), "SNP")

  tmp1 <- withr::local_tempfile(fileext = ".csv")
  m1 <- matrix(1.0, 1, 1, dimnames = list("rs1", "rs1"))
  .write_ld_csv(m1, "rs1", tmp1)
  back1 <- utils::read.csv(tmp1, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(back1$SNP, "rs1")
  expect_equal(back1$rs1, 1.0)
})

# ---- compute_ld_matrix: short-circuit and cache hit --------------------------

minimal_exposure <- function(snps = "rs1", ea = "A", oa = "G") {
  tibble::tibble(
    id.exposure            = "exp1",
    SNP                    = snps,
    beta.exposure          = rep(0.05, length(snps)),
    se.exposure            = rep(0.01, length(snps)),
    effect_allele.exposure = rep(ea, length(snps))[seq_along(snps)],
    other_allele.exposure  = rep(oa, length(snps))[seq_along(snps)],
    pval.exposure          = rep(1e-8, length(snps))
  )
}

test_that("compute_ld_matrix() short-circuits on a single-SNP exposure without invoking PLINK", {
  cache <- withr::local_tempdir()
  withr::local_envvar(ARDMR_CACHE_DIR = cache)
  out_csv <- file.path(cache, "ld.csv")

  res <- compute_ld_matrix(
    exposure_snps = minimal_exposure(),
    ancestry      = "EUR",
    cache_dir     = cache,
    out_csv       = out_csv,
    confirm       = "no",
    verbose       = FALSE
  )
  expect_equal(dim(res$matrix), c(1L, 1L))
  expect_equal(unname(res$matrix[1, 1]), 1.0)
  expect_equal(res$snps_used, "rs1")
  expect_equal(res$ld_pop, "EUR")
  expect_true(file.exists(out_csv))
  expect_true(file.exists(res$cache_path))

  back <- utils::read.csv(out_csv, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(back$SNP, "rs1")
  expect_equal(back$rs1, 1.0)
})

test_that("compute_ld_matrix() returns the cached result on the second call (single-SNP path)", {
  cache <- withr::local_tempdir()
  withr::local_envvar(ARDMR_CACHE_DIR = cache)

  exp <- minimal_exposure()
  res1 <- compute_ld_matrix(
    exposure_snps = exp,
    ancestry      = "EUR",
    cache_dir     = cache,
    confirm       = "no",
    verbose       = FALSE
  )
  expect_false(isTRUE(res1$cached))

  res2 <- compute_ld_matrix(
    exposure_snps = exp,
    ancestry      = "EUR",
    cache_dir     = cache,
    confirm       = "no",
    verbose       = FALSE
  )
  expect_true(isTRUE(res2$cached))
  expect_equal(res1$matrix, res2$matrix)
  expect_equal(res1$cache_path, res2$cache_path)
})

test_that("compute_ld_matrix() force_refresh re-runs even when cache exists", {
  skip_if_not_installed("mockery")

  cache <- withr::local_tempdir()
  withr::local_envvar(ARDMR_CACHE_DIR = cache)
  exp <- minimal_exposure()

  # Seed the cache via a normal call (single-SNP path doesn't touch PLINK).
  res1 <- compute_ld_matrix(
    exposure_snps = exp, ancestry = "EUR", cache_dir = cache,
    confirm = "no", verbose = FALSE
  )
  expect_true(file.exists(res1$cache_path))

  res2 <- compute_ld_matrix(
    exposure_snps = exp, ancestry = "EUR", cache_dir = cache,
    confirm = "no", verbose = FALSE, force_refresh = TRUE
  )
  expect_false(isTRUE(res2$cached))
  expect_equal(res1$matrix, res2$matrix)  # value is unchanged for single-SNP
})

test_that("compute_ld_matrix() rejects empty/all-NA SNP input", {
  cache <- withr::local_tempdir()
  bad <- minimal_exposure(snps = NA_character_)
  expect_error(
    compute_ld_matrix(bad, ancestry = "EUR", cache_dir = cache,
                      confirm = "no", verbose = FALSE),
    "no usable rsids"
  )
})

test_that("compute_ld_matrix() requires a valid exposure_snps argument", {
  cache <- withr::local_tempdir()
  expect_error(
    compute_ld_matrix("not a data frame", ancestry = "EUR", cache_dir = cache),
    "must be a data frame"
  )
  bad <- minimal_exposure()
  bad$pval.exposure <- NULL
  expect_error(
    compute_ld_matrix(bad, ancestry = "EUR", cache_dir = cache,
                      confirm = "no", verbose = FALSE),
    "pval"
  )
})

test_that("compute_ld_matrix() requires an ancestry value", {
  cache <- withr::local_tempdir()
  expect_error(
    compute_ld_matrix(minimal_exposure(), ancestry = "", cache_dir = cache),
    "ancestry"
  )
})
