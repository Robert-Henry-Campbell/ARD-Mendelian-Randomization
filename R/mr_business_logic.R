#' Harmonise, run MR methods, sensitivity & QC across outcomes
#' @param MR_df tibble with outcome_snps list-column
#' @param exposure_snps exposure instruments
#' @param sensitivity_enabled vector of checks
#' @param sensitivity_pass_min integer threshold
#' @param scatterplot,snpforestplot,leaveoneoutplot logical
#' @param plot_output_dir output dir for optional per-outcome plots
#' @param cache_dir path
#' @param verbose logical
#' @return list(MR_df=updated, results_df=tidy summary)
#' @export
mr_business_logic <- function(
    MR_df, exposure_snps,
    sensitivity_enabled, sensitivity_pass_min,
    scatterplot, snpforestplot, leaveoneoutplot,
    plot_output_dir, cache_dir, verbose = TRUE
) {
  # TODO: for each outcome:
  #  - harmonise exposure/outcome
  #  - run IVW, Egger (slope + intercept, I2_GX), weighted median/mode
  #  - heterogeneity (Q, I2), LOO, Steiger
  #  - compute pass/fail using 8 checks (NA-aware denominator)
  #  - collect per-outcome summary row
  results_df <- tibble::tibble(
    outcome_id = character(),
    icd10 = character(),
    gbd_cause = character(),
    beta_ivw = numeric(),
    se_ivw = numeric(),
    p_ivw = numeric(),
    q_ivw = numeric(),      # adjusted later
    qc_pass = logical(),
    nsnp_after = integer()
  )
  logger::log_info("MR business logic: {nrow(results_df)} outcomes analysed")
  list(MR_df = MR_df, results_df = results_df)
}
