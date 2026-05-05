# tests/testthat/test_cache_management.R

build_fake_cache <- function() {
  d <- tempfile("fake_cache_"); dir.create(d, recursive = TRUE)
  for (sub in c("panukb_outcome_snps/EUR", "panukb_outcome_snps/AFR",
                "neale_outcome_snps/both",
                "coloc/exposure", "coloc/outcome", "coloc/meta",
                "ld", "ld_reference", "tabix_index", "hts-cache",
                "output/exp_a/both/EUR/abc1234567",
                "output/exp_b/both/EUR/def9876543")) {
    dir.create(file.path(d, sub), recursive = TRUE)
  }
  saveRDS(tibble::tibble(SNP = "rs1"),
          file.path(d, "panukb_outcome_snps/EUR/synth__aaaaaaaa11.rds"))
  saveRDS(tibble::tibble(SNP = "rs2"),
          file.path(d, "panukb_outcome_snps/EUR/synth__bbbbbbbb22.rds"))
  saveRDS(tibble::tibble(SNP = "rs3"),
          file.path(d, "panukb_outcome_snps/AFR/synth__aaaaaaaa11.rds"))
  saveRDS(tibble::tibble(beta = 0.1),
          file.path(d, "coloc/exposure/key1.rds"))
  saveRDS(list(matrix = matrix(1)),
          file.path(d, "ld/abc.rds"))
  writeLines("manifest", file.path(d, "output/exp_a/both/EUR/abc1234567/run_manifest.json"))
  writeLines("manifest", file.path(d, "output/exp_b/both/EUR/def9876543/run_manifest.json"))
  d
}

test_that("ardmr_cache_summary lists categories with counts", {
  d <- build_fake_cache()
  s <- ardmr_cache_summary(d)
  expect_s3_class(s, "tbl_df")
  expect_true(all(c("panukb_outcome_snps","coloc_exposure","ld_matrix","outputs") %in% s$category))
  expect_equal(s$entries[s$category == "panukb_outcome_snps"], 3L)
  expect_equal(s$entries[s$category == "coloc_exposure"], 1L)
  expect_equal(s$entries[s$category == "ld_matrix"], 1L)
  expect_gt(s$entries[s$category == "outputs"], 0L)
})

test_that("ardmr_cache_summary handles empty categories", {
  d <- build_fake_cache()
  s <- ardmr_cache_summary(d)
  expect_equal(s$entries[s$category == "neale_outcome_snps"], 0L)
})

test_that("ardmr_cache_clean dry_run lists candidates without deleting", {
  d <- build_fake_cache()
  before <- list.files(d, recursive = TRUE)
  cands <- ardmr_cache_clean(d, category = "panukb_outcome_snps", dry_run = TRUE)
  after <- list.files(d, recursive = TRUE)
  expect_equal(before, after)
  expect_equal(length(cands), 3L)
})

test_that("ardmr_cache_clean prunes the requested category only", {
  d <- build_fake_cache()
  before_coloc <- list.files(file.path(d, "coloc/exposure"))
  ardmr_cache_clean(d, category = "panukb_outcome_snps", dry_run = FALSE)
  expect_equal(length(list.files(file.path(d, "panukb_outcome_snps"), recursive = TRUE)), 0L)
  expect_equal(list.files(file.path(d, "coloc/exposure")), before_coloc)
})

test_that("ardmr_cache_clean filters by ancestry on panukb_outcome_snps", {
  d <- build_fake_cache()
  ardmr_cache_clean(d, category = "panukb_outcome_snps",
                    ancestry = "EUR", dry_run = FALSE)
  expect_equal(length(list.files(file.path(d, "panukb_outcome_snps/EUR"))), 0L)
  expect_equal(length(list.files(file.path(d, "panukb_outcome_snps/AFR"))), 1L)
})

test_that("ardmr_cache_clean filters outputs by run_hash", {
  d <- build_fake_cache()
  ardmr_cache_clean(d, category = "outputs",
                    run_hash = "abc1234567", dry_run = FALSE)
  expect_false(dir.exists(file.path(d, "output/exp_a/both/EUR/abc1234567")))
  expect_true(dir.exists(file.path(d, "output/exp_b/both/EUR/def9876543")))
})

test_that("ardmr_cache_clean errors on unknown category", {
  d <- build_fake_cache()
  expect_error(ardmr_cache_clean(d, category = "nonsense"), "Unknown category")
})

test_that("ardmr_cache_clean requires explicit category (no default)", {
  d <- build_fake_cache()
  expect_error(ardmr_cache_clean(d), "category.*required")
})
