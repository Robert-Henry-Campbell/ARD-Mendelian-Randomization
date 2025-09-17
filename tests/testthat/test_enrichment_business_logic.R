test_that("beta_contrast_global_ARD returns confidence intervals", {
  df <- tibble::tibble(
    results_beta_ivw = c(0.2, 0.15, -0.05, 0.0),
    results_se_ivw = c(0.1, 0.12, 0.2, 0.18),
    ARD_selected = c(TRUE, TRUE, FALSE, FALSE)
  )

  res <- ardmr:::beta_contrast_global_ARD(
    results_df = df,
    use_qc_pass = FALSE,
    min_nsnp = NULL
  )

  expect_true(all(c("ci_low", "ci_high") %in% names(res)))
  expect_type(res$ci_low, "double")
  expect_type(res$ci_high, "double")
  expect_equal(res$ci_low, res$delta_beta - 1.96 * res$se_delta)
  expect_equal(res$ci_high, res$delta_beta + 1.96 * res$se_delta)
})

test_that("beta_contrast_global_ARD NA guard returns CI columns", {
  df_all_ard <- tibble::tibble(
    results_beta_ivw = c(0.2, 0.15),
    results_se_ivw = c(0.1, 0.12),
    ARD_selected = c(TRUE, TRUE)
  )

  res_na <- ardmr:::beta_contrast_global_ARD(
    results_df = df_all_ard,
    use_qc_pass = FALSE,
    min_nsnp = NULL
  )

  expect_true(all(c("ci_low", "ci_high") %in% names(res_na)))
  expect_true(all(is.na(res_na$ci_low)))
  expect_true(all(is.na(res_na$ci_high)))
})
