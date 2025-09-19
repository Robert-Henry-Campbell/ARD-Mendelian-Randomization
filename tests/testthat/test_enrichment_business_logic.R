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

test_that("beta_mean_by_cause computes IVW means and drops non-ARD causes", {
  df <- tibble::tibble(
    results_beta_ivw = c(0.2, 0.1, -0.15, 0.05),
    results_se_ivw = c(0.1, 0.2, 0.1, 0.15),
    results_qc_pass = c(TRUE, TRUE, FALSE, TRUE),
    results_nsnp_after = c(3, 4, 3, 2),
    cause_level_1 = c("Cardio", "Cardio", "Neuro", "Neuro"),
    ARD_selected = c(TRUE, TRUE, FALSE, TRUE)
  )

  res_all <- ardmr:::beta_mean_by_cause(
    results_df = df,
    level = "cause_level_1",
    use_qc_pass = TRUE,
    min_nsnp = 2,
    exposure = "TestExposure",
    ard_only = FALSE
  )

  expect_equal(res_all$exposure, rep("TestExposure", nrow(res_all)))

  cardio_row <- res_all[res_all$cause == "Cardio", , drop = FALSE]
  expect_equal(cardio_row$n, 2L)
  expect_equal(cardio_row$ivw_mean_beta, 0.18, tolerance = 1e-8)
  expect_equal(cardio_row$se_ivw_mean, sqrt(1 / (1/0.1^2 + 1/0.2^2)), tolerance = 1e-8)

  res_ard <- ardmr:::beta_mean_by_cause(
    results_df = df,
    level = "cause_level_1",
    use_qc_pass = TRUE,
    min_nsnp = 2,
    exposure = "TestExposure",
    ard_only = TRUE
  )

  expect_true("Cardio" %in% res_ard$cause)
  expect_false("Neuro" %in% res_ard$cause)
})

test_that("beta_mean_global returns all and ARD summaries", {
  df <- tibble::tibble(
    results_beta_ivw = c(0.2, 0.1, -0.05, 0.0),
    results_se_ivw = c(0.1, 0.2, 0.15, 0.12),
    results_qc_pass = TRUE,
    results_nsnp_after = 3,
    ARD_selected = c(TRUE, TRUE, FALSE, TRUE)
  )

  res_global <- ardmr:::beta_mean_global(
    results_df = df,
    use_qc_pass = TRUE,
    min_nsnp = 2,
    exposure = "TestExposure"
  )

  expect_true(all(c("All Diseases", "Age-Related Diseases") %in% res_global$group))

  all_row <- res_global[res_global$group == "All Diseases", , drop = FALSE]
  ard_row <- res_global[res_global$group == "Age-Related Diseases", , drop = FALSE]

  w_all <- sum(1 / df$results_se_ivw^2)
  mean_all <- sum((1 / df$results_se_ivw^2) * df$results_beta_ivw) / w_all
  expect_equal(all_row$ivw_mean_beta, mean_all)
  expect_equal(all_row$se_ivw_mean, sqrt(1 / w_all))

  df_ard <- df[df$ARD_selected %in% TRUE, ]
  w_ard <- sum(1 / df_ard$results_se_ivw^2)
  mean_ard <- sum((1 / df_ard$results_se_ivw^2) * df_ard$results_beta_ivw) / w_ard
  expect_equal(ard_row$ivw_mean_beta, mean_ard)
  expect_equal(ard_row$se_ivw_mean, sqrt(1 / w_ard))
})
