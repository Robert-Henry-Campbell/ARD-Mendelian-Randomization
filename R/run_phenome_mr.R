#' Run phenome-wide MR on ARDs (orchestrator)
#'
#' Builds the outcome plan, maps exposure SNPs, fetches outcome SNP rows,
#' runs MR + sensitivity/QC, and returns results + plots.
#'
#' @param exposure_snps Data frame of exposure instruments (TwoSampleMR-like: rsid, beta, se, effect_allele, other_allele, eaf, etc.).
#' @param ancestry Character (mandatory), e.g. "EUR".
#' @param sex One of "both","male","female". If not "both", you likely use Neale.
#' @param sensitivity_enabled Character vector of checks (default = all 8).
#' @param sensitivity_pass_min Integer threshold to pass QC (default 6).
#' @param Multiple_testing_correction "BH" or "bonferroni" (default "BH").
#' @param scatterplot,snpforestplot,leaveoneoutplot Logical; produce per-outcome diagnostics.
#' @param plot_output_dir Directory to write plots ("" = do not write).
#' @param Neale_GWAS_dir Optional path to Neale sumstats; created if missing.
#' @param cache_dir Cache directory for temporary files (default: user cache).
#' @param verbose Logical; print progress messages.
#'
#' @return A list: MR_df, results_df, manhattan (ggplot), volcano (ggplot)
#' @export
run_phenome_mr <- function(
    exposure_snps,
    ancestry,
    sex = c("both","male","female"),
    sensitivity_enabled = c(
      "egger_intercept","egger_slope_agreement",
      "weighted_median","weighted_mode",
      "steiger_direction","leave_one_out",
      "ivw_Q","ivw_I2"
    ),
    sensitivity_pass_min = 6,
    Multiple_testing_correction = c("BH","bonferroni"),
    scatterplot = FALSE,
    snpforestplot = FALSE,
    leaveoneoutplot = FALSE,
    plot_output_dir = "",
    Neale_GWAS_dir = NULL,
    cache_dir = tools::R_user_dir("ardmr","cache"),
    verbose = TRUE
) {
  # ---- validate args ----
  sex <- match.arg(sex)
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  if (missing(ancestry) || !nzchar(ancestry)) stop("`ancestry` is mandatory.")
  .assert_exposure(exposure_snps)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (nzchar(plot_output_dir) && !dir.exists(plot_output_dir)) {
    dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!is.null(Neale_GWAS_dir) && !dir.exists(Neale_GWAS_dir)) {
    dir.create(Neale_GWAS_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Keep config bundled for helper calls
  cfg <- list(
    ancestry = ancestry,
    sex = sex,
    checks_enabled = sensitivity_enabled,
    checks_pass_min = sensitivity_pass_min,
    mtc = Multiple_testing_correction,
    diag_scatter = scatterplot,
    diag_forest  = snpforestplot,
    diag_loo     = leaveoneoutplot,
    plot_dir = plot_output_dir,
    cache_dir = cache_dir,
    neale_dir = Neale_GWAS_dir,
    verbose = verbose
  )

  if (cfg$verbose) message("1) Outcome setup…")
  MR_df <- outcome_setup(sex = cfg$sex, ancestry = cfg$ancestry, verbose = cfg$verbose)

  if (cfg$verbose) message("2) Map exposure SNPs to provider positions…")
  exposure_snps2 <- exposure_snp_mapper(exposure_snps, sex = cfg$sex, cache_dir = cfg$cache_dir, verbose = cfg$verbose)

  if (cfg$sex == "both") {
    if (cfg$verbose) message("3) Pull Pan-UKB outcome SNP rows…")
    MR_df <- panukb_snp_grabber(exposure_snps2, MR_df, ancestry = cfg$ancestry, cache_dir = cfg$cache_dir, verbose = cfg$verbose)
  } else {
    if (cfg$verbose) message("3a) Ensure Neale GWAS files + tbi present…")
    neale_gwas_checker(MR_df, neale_dir = cfg$neale_dir, verbose = cfg$verbose)
    neale_tbi_maker(neale_dir = cfg$neale_dir, verbose = cfg$verbose)
    if (cfg$verbose) message("3b) Pull Neale outcome SNP rows…")
    MR_df <- neale_snp_grabber(exposure_snps2, MR_df, neale_dir = cfg$neale_dir, cache_dir = cfg$cache_dir, verbose = cfg$verbose)
  }

  if (cfg$verbose) message("4) Run MR + sensitivity/QC…")
  mr_out <- mr_business_logic(
    MR_df = MR_df,
    exposure_snps = exposure_snps2,
    sensitivity_enabled = cfg$checks_enabled,
    sensitivity_pass_min = cfg$checks_pass_min,
    scatterplot = cfg$diag_scatter,
    snpforestplot = cfg$diag_forest,
    leaveoneoutplot = cfg$diag_loo,
    plot_output_dir = cfg$plot_dir,
    cache_dir = cfg$cache_dir,
    verbose = cfg$verbose
  )
  MR_df <- mr_out$MR_df
  results_df <- mr_out$results_df

  if (cfg$verbose) message("5) Build summary plots…")
  manhattan <- manhattan_plot(results_df, Multiple_testing_correction = cfg$mtc)
  volcano   <- volcano_plot(results_df, Multiple_testing_correction = cfg$mtc)

  if (nzchar(cfg$plot_dir)) {
    ggplot2::ggsave(file.path(cfg$plot_dir, "manhattan.png"), manhattan, width = 10, height = 6, dpi = 150)
    ggplot2::ggsave(file.path(cfg$plot_dir, "volcano.png"),   volcano,   width = 10, height = 6, dpi = 150)
  }

  invisible(list(
    MR_df = MR_df,
    results_df = results_df,
    manhattan = manhattan,
    volcano = volcano
  ))
}

# small internal assert (keeps orchestrator readable)
.assert_exposure <- function(x) {
  req <- c("rsid","beta","se","effect_allele","other_allele")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("exposure_snps missing columns: ", paste(miss, collapse=", "))
  invisible(TRUE)
}
