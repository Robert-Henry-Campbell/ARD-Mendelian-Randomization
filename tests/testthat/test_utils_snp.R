# tests/testthat/test_utils_snp.R
# Unit tests for the small SNP-shaped helpers in R/utils-snp.R.

test_that("is_rsid recognises rs IDs and rejects everything else", {
  expect_true(ardmr:::is_rsid("rs1234"))
  expect_true(ardmr:::is_rsid("rs1"))
  expect_true(ardmr:::is_rsid("rs999999999"))
  expect_false(ardmr:::is_rsid("chr1:12345"))
  expect_false(ardmr:::is_rsid("RS1234"))     # case-sensitive
  expect_false(ardmr:::is_rsid("rs"))
  expect_false(ardmr:::is_rsid("rs12a"))
  expect_false(ardmr:::is_rsid(NA_character_))
  expect_false(ardmr:::is_rsid(""))
})

test_that("is_rsid is vectorised", {
  v <- c("rs1", "chr1:1", NA, "rs2")
  expect_equal(ardmr:::is_rsid(v), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("palindrome_flag flags A/T and C/G pairs only", {
  expect_true(ardmr:::palindrome_flag("A", "T"))
  expect_true(ardmr:::palindrome_flag("T", "A"))
  expect_true(ardmr:::palindrome_flag("C", "G"))
  expect_true(ardmr:::palindrome_flag("G", "C"))
  expect_false(ardmr:::palindrome_flag("A", "G"))
  expect_false(ardmr:::palindrome_flag("C", "T"))
  expect_false(ardmr:::palindrome_flag("A", "C"))
})

test_that("palindrome_flag is case-insensitive and vectorised", {
  expect_true(ardmr:::palindrome_flag("a", "t"))
  expect_true(ardmr:::palindrome_flag("a", "T"))
  expect_equal(
    ardmr:::palindrome_flag(c("A", "C", "A"), c("T", "T", "G")),
    c(TRUE, FALSE, FALSE)
  )
})

test_that("f_stat returns Wald chi-squared (beta/se)^2", {
  expect_equal(ardmr:::f_stat(0.05, 0.01), 25)
  expect_equal(ardmr:::f_stat(-0.05, 0.01), 25)
  expect_equal(ardmr:::f_stat(0, 0.01), 0)
})

test_that("f_stat is vectorised", {
  expect_equal(
    ardmr:::f_stat(c(0.05, 0.10, 0.02), c(0.01, 0.01, 0.01)),
    c(25, 100, 4)
  )
})

test_that("safe_tophits returns NULL when tophits errors out", {
  testthat::local_mocked_bindings(
    tophits = function(id, p = NULL, pval = NULL, ...) stop("server down"),
    .package = "ieugwasr"
  )
  expect_null(ardmr:::safe_tophits("ieu-a-1", 5e-8))
})

test_that("safe_tophits returns 0-row tibble when tophits returns no rows", {
  testthat::local_mocked_bindings(
    tophits = function(id, p = NULL, pval = NULL, ...) {
      tibble::tibble(rsid = character(), p = numeric())
    },
    .package = "ieugwasr"
  )
  out <- ardmr:::safe_tophits("ieu-a-1", 5e-8)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
})

test_that("safe_tophits synthesises 'p' column when ieugwasr returns 'pval'", {
  testthat::local_mocked_bindings(
    tophits = function(id, p = NULL, pval = NULL, ...) {
      # Simulate older ieugwasr that returns 'pval' instead of 'p'.
      tibble::tibble(rsid = c("rs1", "rs2"), pval = c(1e-9, 1e-10))
    },
    .package = "ieugwasr"
  )
  out <- ardmr:::safe_tophits("ieu-a-1", 5e-8)
  expect_true("p" %in% names(out))
  expect_equal(out$p, c(1e-9, 1e-10))
})
