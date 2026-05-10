# tests/testthat/test_coloc_data_fetch_neale.R
# Behaviour of the Neale coloc fetcher, slice helper, and metadata function.

cache_dir_test <- function() {
  d <- tempfile(pattern = "ardmr_coloc_neale_test_")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

fixture_path <- function() {
  p <- system.file("data", "neale", "sumstats", "TEST_neale_sumstats.tsv", package = "ardmr")
  if (!nzchar(p)) p <- testthat::test_path("..", "..", "inst", "data", "neale", "sumstats", "TEST_neale_sumstats.tsv")
  p
}

# --- metadata ----------------------------------------------------------------

test_that("coloc_neale_outcome_metadata: binary cc trait", {
  rec <- tibble::tibble(
    variable_type = "binary",
    n_cases = 100, n_controls = 1000, n_non_missing = 1100
  )
  m <- coloc_neale_outcome_metadata(rec, sex = "female")
  expect_equal(m$type, "cc")
  expect_equal(m$N, 1100)
  expect_equal(m$ncase, 100)
  expect_equal(m$ncontrol, 1000)
  expect_equal(m$s, 100 / 1100)
})

test_that("coloc_neale_outcome_metadata: categorical with cases is cc", {
  rec <- tibble::tibble(
    variable_type = "categorical",
    n_cases = 149, n_controls = 194025, n_non_missing = 194174
  )
  m <- coloc_neale_outcome_metadata(rec, sex = "female")
  expect_equal(m$type, "cc")
  expect_equal(m$N, 194174)
  expect_equal(m$s, 149 / 194174)
})

test_that("coloc_neale_outcome_metadata: continuous_irnt is quant", {
  rec <- tibble::tibble(
    variable_type = "continuous_irnt",
    n_cases = NA_real_, n_controls = NA_real_, n_non_missing = 5000
  )
  m <- coloc_neale_outcome_metadata(rec, sex = "male")
  expect_equal(m$type, "quant")
  expect_equal(m$N, 5000)
  expect_true(is.na(m$ncase))
  expect_true(is.na(m$ncontrol))
  expect_true(is.na(m$s))
})

test_that("coloc_neale_outcome_metadata: N falls back to ncase + ncontrol when n_non_missing missing", {
  rec <- tibble::tibble(
    variable_type = "binary",
    n_cases = 200, n_controls = 800, n_non_missing = NA_real_
  )
  m <- coloc_neale_outcome_metadata(rec, sex = "female")
  expect_equal(m$N, 1000)
  expect_equal(m$type, "cc")
})

test_that("coloc_neale_outcome_metadata: graceful when columns missing", {
  rec <- tibble::tibble(variable_type = "binary")
  m <- coloc_neale_outcome_metadata(rec, sex = "female")
  expect_true(is.na(m$N))
  expect_true(is.na(m$ncase))
  expect_true(is.na(m$ncontrol))
})

# --- slice helper ------------------------------------------------------------

test_that(".coloc_neale_chr_slice: parses, filters to chr, caches", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()

  s1 <- .coloc_neale_chr_slice(fp, "6", cache_dir = cd, verbose = FALSE)
  expect_true(nrow(s1) >= 8)
  expect_true(all(s1$v_chr == "6"))
  expect_true(is.integer(s1$v_pos))
  expect_true(is.numeric(s1$beta))
  # second call should be a cache hit (no re-read)
  s2 <- .coloc_neale_chr_slice(fp, "6", cache_dir = cd, verbose = FALSE)
  expect_identical(s1, s2)
})

test_that(".coloc_neale_chr_slice: returns chrX as character", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()
  sx <- .coloc_neale_chr_slice(fp, "X", cache_dir = cd, verbose = FALSE)
  expect_equal(nrow(sx), 2L)
  expect_true(all(sx$v_chr == "X"))
})

test_that(".coloc_neale_chr_slice: missing file returns empty tibble (no stop)", {
  cd <- cache_dir_test()
  s <- .coloc_neale_chr_slice(file.path(cd, "does_not_exist.tsv"), "1",
                              cache_dir = cd, verbose = FALSE)
  expect_s3_class(s, "tbl_df")
  expect_equal(nrow(s), 0L)
})

test_that(".coloc_neale_chr_slice: force_refresh re-reads file", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()
  .coloc_neale_chr_slice(fp, "1", cache_dir = cd, verbose = FALSE)
  cache_files <- list.files(file.path(cd, "coloc", "neale_chr_slices"), full.names = TRUE)
  expect_true(length(cache_files) >= 1)
  # corrupt the cache and confirm force_refresh ignores it
  saveRDS(tibble::tibble(), cache_files[[1]])
  s <- .coloc_neale_chr_slice(fp, "1", cache_dir = cd, verbose = FALSE, force_refresh = TRUE)
  expect_true(nrow(s) > 0)
})

# --- regional fetcher --------------------------------------------------------

test_that("coloc_fetch_outcome_neale_region: returns canonical column shape", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()
  rec <- tibble::tibble(File = basename(fp))
  out <- coloc_fetch_outcome_neale_region(
    rec = rec, sex = "female", chr = "6",
    start = 38500000L, end = 39500000L,
    neale_dir = dirname(fp), cache_dir = cd, verbose = FALSE
  )
  expect_named(out, c("SNP", "chr", "pos", "effect_allele", "other_allele",
                      "beta", "se", "eaf", "pval"))
  expect_true(all(out$pos >= 38500000L & out$pos <= 39500000L))
  expect_true(all(out$chr == "6"))
})

test_that("coloc_fetch_outcome_neale_region: drops low_confidence_variant rows", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()
  rec <- tibble::tibble(File = basename(fp))
  out <- coloc_fetch_outcome_neale_region(
    rec = rec, sex = "female", chr = "1",
    start = 4500L, end = 5500L,  # window covers the low_confidence row at 5000
    neale_dir = dirname(fp), cache_dir = cd, verbose = FALSE
  )
  # row "1:5000:G:A" has low_confidence_variant=true and must be filtered out
  expect_false("1:5000_G_A" %in% out$SNP)
})

test_that("coloc_fetch_outcome_neale_region: drops NaN beta/se rows", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()
  rec <- tibble::tibble(File = basename(fp))
  out <- coloc_fetch_outcome_neale_region(
    rec = rec, sex = "female", chr = "1",
    start = 5500L, end = 6500L,  # window covers the NaN row at 6000
    neale_dir = dirname(fp), cache_dir = cd, verbose = FALSE
  )
  # row "1:6000:T:C" has NaN beta/se and must be filtered out
  expect_false("1:6000_T_C" %in% out$SNP)
})

test_that("coloc_fetch_outcome_neale_region: EAF is alt-allele frequency (flips when minor=ref)", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()
  rec <- tibble::tibble(File = basename(fp))
  out <- coloc_fetch_outcome_neale_region(
    rec = rec, sex = "female", chr = "1",
    start = 8500L, end = 9500L,  # window covers row at 9000 where minor_allele=G=ref
    neale_dir = dirname(fp), cache_dir = cd, verbose = FALSE
  )
  # row "1:9000:G:A": ref=G, alt=A, minor_allele=G (=ref), minor_AF=0.30 → eaf = 0.70
  row <- out[out$SNP == "1:9000_G_A", , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_equal(row$effect_allele, "A")
  expect_equal(row$other_allele,  "G")
  expect_equal(row$eaf, 0.70, tolerance = 1e-6)
})

test_that("coloc_fetch_outcome_neale_region: EAF preserved when minor=alt", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()
  rec <- tibble::tibble(File = basename(fp))
  out <- coloc_fetch_outcome_neale_region(
    rec = rec, sex = "female", chr = "1",
    start = 500L, end = 1500L,
    neale_dir = dirname(fp), cache_dir = cd, verbose = FALSE
  )
  # row "1:1000:A:G": ref=A, alt=G, minor_allele=G (=alt), minor_AF=0.10 → eaf = 0.10
  row <- out[out$SNP == "1:1000_A_G", , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_equal(row$effect_allele, "G")
  expect_equal(row$other_allele, "A")
  expect_equal(row$eaf, 0.10, tolerance = 1e-6)
})

test_that("coloc_fetch_outcome_neale_region: missing file returns empty (no stop)", {
  cd <- cache_dir_test()
  rec <- tibble::tibble(File = "no_such_file.tsv.bgz")
  out <- coloc_fetch_outcome_neale_region(
    rec = rec, sex = "female", chr = "6", start = 1L, end = 1e8L,
    neale_dir = cd, cache_dir = cd, verbose = FALSE
  )
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 0L)
})

test_that("coloc_fetch_outcome_neale_region: empty window returns empty tibble", {
  fp <- fixture_path()
  testthat::skip_if(!file.exists(fp), "test fixture missing")
  cd <- cache_dir_test()
  rec <- tibble::tibble(File = basename(fp))
  out <- coloc_fetch_outcome_neale_region(
    rec = rec, sex = "female", chr = "6",
    start = 1e9L, end = 2e9L,  # window with no SNPs
    neale_dir = dirname(fp), cache_dir = cd, verbose = FALSE
  )
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 0L)
})
