# tests/testthat/test_utils_data.R
#
# Direct-assertion tests for the small pure utilities in R/utils-data.R:
#   * assert_exposure()
#   * fix_exposure_names()
#   * coerce_lowercase_colnames()
#   * get_pkg_obj()
#
# Fully offline — uses only in-memory tibbles plus packaged .rda data.

# ---- assert_exposure ----------------------------------------------------------

complete_exposure <- function() {
  tibble::tibble(
    id.exposure            = "exp1",
    SNP                    = "rs1",
    beta.exposure          = 0.05,
    se.exposure            = 0.01,
    effect_allele.exposure = "A",
    other_allele.exposure  = "G",
    pval.exposure          = 1e-8
  )
}

test_that("assert_exposure() accepts a valid exposure tibble", {
  expect_invisible(assert_exposure(complete_exposure()))
  expect_true(assert_exposure(complete_exposure()))
})

test_that("assert_exposure() errors and names the missing column(s)", {
  bad <- complete_exposure()
  bad$pval.exposure <- NULL
  expect_error(assert_exposure(bad), "pval\\.exposure")

  bad2 <- complete_exposure()
  bad2$beta.exposure <- NULL
  bad2$se.exposure <- NULL
  err <- tryCatch(assert_exposure(bad2), error = function(e) conditionMessage(e))
  expect_match(err, "beta\\.exposure")
  expect_match(err, "se\\.exposure")
})

# ---- fix_exposure_names ------------------------------------------------------

test_that("fix_exposure_names() lowercases names and adds .exposure suffix", {
  raw <- tibble::tibble(
    SNP           = "rs1",
    Beta          = 0.05,
    SE            = 0.01,
    Effect_Allele = "A",
    Other_Allele  = "G",
    Pval          = 1e-8,
    ID            = "exp1"
  )
  out <- fix_exposure_names(raw)
  expect_setequal(
    names(out),
    c("SNP", "beta.exposure", "se.exposure",
      "effect_allele.exposure", "other_allele.exposure",
      "pval.exposure", "id.exposure")
  )
  expect_equal(out$beta.exposure, 0.05)
})

test_that("fix_exposure_names() preserves uppercase 'SNP' even when input is 'snp' or 'Snp'", {
  for (nm in c("SNP", "snp", "Snp")) {
    df <- tibble::tibble(x = 1)
    names(df)[1] <- nm
    out <- fix_exposure_names(df)
    expect_true("SNP" %in% names(out),
                info = paste("input column name:", nm))
  }
})

test_that("fix_exposure_names() is idempotent on already-suffixed input", {
  already <- tibble::tibble(
    SNP                    = "rs1",
    beta.exposure          = 0.05,
    se.exposure            = 0.01,
    effect_allele.exposure = "A",
    other_allele.exposure  = "G",
    pval.exposure          = 1e-8,
    id.exposure            = "exp1"
  )
  twice <- fix_exposure_names(fix_exposure_names(already))
  expect_equal(names(twice), names(already))
})

test_that("fix_exposure_names() collapses 'exposure.exposure' to 'exposure'", {
  raw <- tibble::tibble(SNP = "rs1", exposure = "Trait A")
  out <- fix_exposure_names(raw)
  expect_true("exposure" %in% names(out))
  expect_false("exposure.exposure" %in% names(out))
})

# ---- coerce_lowercase_colnames -----------------------------------------------

test_that("coerce_lowercase_colnames() lowercases names and preserves values", {
  raw <- data.frame(A = 1:3, B = c("x", "y", "z"), stringsAsFactors = FALSE)
  out <- coerce_lowercase_colnames(raw)
  expect_equal(names(out), c("a", "b"))
  expect_equal(out$a, 1:3)
  expect_equal(out$b, c("x", "y", "z"))
})

# ---- get_pkg_obj -------------------------------------------------------------

test_that("get_pkg_obj() loads a packaged dataset as a tibble", {
  obj <- get_pkg_obj("bothsex_all")
  expect_s3_class(obj, "tbl_df")
  expect_gt(nrow(obj), 0)
  expect_true("ICD10_explo" %in% names(obj))
})

test_that("get_pkg_obj() errors with a message naming a missing object", {
  # utils::data() emits a "data set ... not found" warning before
  # get_pkg_obj()'s own stop() fires. We're testing the latter.
  expect_error(
    suppressWarnings(get_pkg_obj("definitely_not_a_real_dataset_xyz")),
    "definitely_not_a_real_dataset_xyz"
  )
})
