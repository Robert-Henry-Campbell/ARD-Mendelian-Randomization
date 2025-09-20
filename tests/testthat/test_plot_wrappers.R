test_that("plot_beta_contrast_forest_wrap returns base plot for empty data", {
  empty_tbl <- tibble::tibble(
    cause = character(),
    delta_beta = double(),
    se_delta = double(),
    ci_low = double(),
    ci_high = double(),
    p = double(),
    q = double(),
    sig = logical()
  )

  expect_s3_class(plot_beta_contrast_forest_wrap(empty_tbl), "ggplot")
})

test_that("plot_beta_mean_forest_wrap returns base plot for empty data", {
  empty_tbl <- tibble::tibble(
    cause = character(),
    ivw_mean_beta = double(),
    se_ivw_mean = double(),
    ci_low = double(),
    ci_high = double()
  )

  expect_s3_class(plot_beta_mean_forest_wrap(empty_tbl, exposure_units = "SD"), "ggplot")
})

test_that("run_phenome_mr handles empty ARD-only tables and wraps plots", {
  cache_dir <- tempfile("ardmr_cache_")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  old_cache <- Sys.getenv("ARDMR_CACHE_DIR", unset = NA_character_)
  Sys.setenv(ARDMR_CACHE_DIR = cache_dir)
  on.exit({
    if (is.na(old_cache)) Sys.unsetenv("ARDMR_CACHE_DIR") else Sys.setenv(ARDMR_CACHE_DIR = old_cache)
    unlink(cache_dir, recursive = TRUE)
  }, add = TRUE)

  fake_exposure <- tibble::tibble(
    id.exposure = "fake",
    SNP = "rs123",
    beta.exposure = 0.1,
    se.exposure = 0.01,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-6
  )

  fake_results <- tibble::tibble(
    results_beta_ivw = c(0.2, -0.1, 0.05),
    results_se_ivw = c(0.05, 0.04, 0.06),
    results_p_ivw = c(0.01, 0.02, 0.5),
    results_outcome = c("Outcome A", "Outcome B", "Outcome C"),
    results_qc_pass = TRUE,
    results_nsnp_after = 3L,
    cause_level_1 = c("Cause 1", "Cause 1", "Cause 2"),
    cause_level_2 = c(
      "Very long cause level two label for wrapping",
      "Another extended cause level two description",
      "Third cause level two label"
    ),
    cause_level_3 = c(
      "Extremely long cause level three label used for wrapping",
      "Second cause level three label also long",
      "Third cause level three label for wrapping"
    ),
    ARD_selected = c(FALSE, FALSE, FALSE)
  )

  fake_mr_df <- tibble::tibble(
    cause_level_1 = fake_results$cause_level_1,
    ARD_selected = fake_results$ARD_selected,
    outcome_snps = replicate(nrow(fake_results), tibble::tibble(), simplify = FALSE)
  )

  fake_enrich <- list(
    global_tbl = tibble::tibble(compare_mode = character(), estimate = double()),
    by_cause_tbl = tibble::tibble(level = character(), compare_mode = character())
  )

  res <- testthat::with_mocked_bindings(
    run_phenome_mr(
      exposure = "Test exposure",
      exposure_snps = fake_exposure,
      exposure_units = "SD",
      ancestry = "EUR",
      sex = "both",
      sensitivity_enabled = character(0),
      sensitivity_pass_min = 0,
      Multiple_testing_correction = "BH",
      scatterplot = FALSE,
      snpforestplot = FALSE,
      leaveoneoutplot = FALSE,
      verbose = FALSE,
      confirm = "no"
    ),
    Outcome_setup = function(sex, ancestry) fake_mr_df,
    Variant_manifest_downloader = function(...) NULL,
    exposure_snp_mapper = function(exposure_snps, ...) tibble::as_tibble(exposure_snps),
    panukb_snp_grabber = function(exposure_snps, MR_df, ...) {
      if (!"outcome_snps" %in% names(MR_df)) {
        MR_df$outcome_snps <- replicate(nrow(MR_df), tibble::tibble(), simplify = FALSE)
      }
      MR_df
    },
    mr_business_logic = function(...) list(MR_df = fake_mr_df, results_df = fake_results),
    run_enrichment = function(...) fake_enrich,
    plot_enrichment_global = function(...) ggplot2::ggplot(),
    plot_enrichment_global_signed = function(...) ggplot2::ggplot(),
    plot_enrichment_directional_forest = function(...) ggplot2::ggplot(),
    plot_enrichment_signed_forest = function(...) ggplot2::ggplot(),
    manhattan_plot = function(...) ggplot2::ggplot(),
    volcano_plot = function(...) ggplot2::ggplot()
  )

  expect_s3_class(res$summary_plots$beta$cause_level_1$all_diseases_wrap, "ggplot")
  expect_s3_class(res$summary_plots$beta$cause_level_1$age_related_diseases_wrap, "ggplot")
  expect_s3_class(res$summary_plots$beta_contrast$cause_level_1$cause_vs_rest_all_wrap, "ggplot")
  expect_s3_class(res$summary_plots$beta_contrast$cause_level_1$ARD_vs_nonARD_within_cause_wrap, "ggplot")

  expect_no_error(ggplot2::ggplot_build(res$summary_plots$beta$cause_level_1$all_diseases_wrap))
  expect_no_error(ggplot2::ggplot_build(res$summary_plots$beta$cause_level_1$age_related_diseases_wrap))
})

