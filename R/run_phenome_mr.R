#R/run_phenome.MR
#' Run phenome-wide MR across phenotypes (orchestrator)
#'
#' Builds the outcome plan, maps exposure SNPs, fetches outcome SNP rows,
#' runs MR + sensitivity/QC, and returns results + plots. Phenotypes flagged by
#' `ARD_selected` denote age-related diseases.
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
#' @param cache_dir Cache directory for temporary files (default:
#'   [ardmr_cache_dir()]).
#' @param logfile Optional path to log file. Defaults to a timestamped file in
#'   `cache_dir` if not supplied.
#' @param verbose Logical; if `FALSE`, only warnings/errors are logged.
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
    scatterplot = TRUE,
    snpforestplot = TRUE,
    leaveoneoutplot = TRUE,
    cache_dir = ardmr_cache_dir(),
    logfile = NULL,
    verbose = TRUE
) {
  # ---- validate args ----
  sex <- match.arg(sex)
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  if (missing(ancestry) || !nzchar(ancestry)) stop("`ancestry` is mandatory.")
  assert_exposure(exposure_snps)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Always derive standard subfolders inside cache_dir
  plot_output_dir <- file.path(cache_dir, "plots")
  Neale_GWAS_dir  <- file.path(cache_dir, "neale_sumstats")
  dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(Neale_GWAS_dir,  recursive = TRUE, showWarnings = FALSE)


  # ---- choose catalog ----
  catalog <- if (sex == "both") {
    "panukb"
  } else {
    "neale"
  }


  logfile <- if (is.null(logfile) || !nzchar(logfile)) {
    file.path(cache_dir, sprintf("ardmr_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
  } else {
    logfile
  }
  setup_logging(logfile, level = if (verbose) "INFO" else "WARN")
  logger::log_info("Logging to {logfile}")

  # Keep config bundled for helper calls
  cfg <- list(
    ancestry = ancestry,
    sex = sex,
    catalog = catalog,
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

  metrics <- list()

  logger::log_info("1) Outcome setup…")
  MR_df <- Outcome_setup(sex = cfg$sex, ancestry = cfg$ancestry)
  metrics$outcomes <- nrow(MR_df)

  metrics$n_ard <- sum(MR_df$ARD_selected)
  metrics$n_non <- metrics$outcomes - metrics$n_ard
  logger::log_info("Outcome setup: {metrics$outcomes} phenotypes loaded ({metrics$n_ard} ARD, {metrics$n_non} non-ARD)")

  logger::log_info("1.1) Variant manifest downloading…")
  Variant_manifest_downloader(catalog = cfg$catalog, cache_dir = cfg$cache_dir, overwrite = FALSE)

  logger::log_info("2) Map exposure SNPs to provider positions…")
  n_before <- nrow(exposure_snps)
  exposure_snps2 <- exposure_snp_mapper(exposure_snps, sex = cfg$sex, cache_dir = cfg$cache_dir, verbose = cfg$verbose)
  metrics$exposure_in <- n_before
  metrics$exposure_mapped <- nrow(exposure_snps2)
  logger::log_info("Exposure mapping: {metrics$exposure_in - metrics$exposure_mapped} filtered; {metrics$exposure_mapped} remain")

  if (cfg$sex == "both") {
    logger::log_info("3) Pull Pan-UKB outcome SNP rows…")
    MR_df <- panukb_snp_grabber(exposure_snps2, MR_df, ancestry = cfg$ancestry, cache_dir = cfg$cache_dir, verbose = cfg$verbose)
  } else {
    logger::log_info("3a) Ensure Neale GWAS files + tbi present…")
    neale_gwas_checker(MR_df, neale_dir = cfg$neale_dir, verbose = cfg$verbose)
    neale_tbi_maker(neale_dir = cfg$neale_dir, verbose = cfg$verbose)
    logger::log_info("3b) Pull Neale outcome SNP rows…")
    MR_df <- neale_snp_grabber(exposure_snps2, MR_df, neale_dir = cfg$neale_dir, cache_dir = cfg$cache_dir, verbose = cfg$verbose)
  }
  metrics$outcome_snps <- sum(lengths(MR_df$outcome_snps))
  logger::log_info("Outcome SNPs: {metrics$outcome_snps} rows fetched")

  logger::log_info("4) Run MR + sensitivity/QC…")
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
  metrics$results <- nrow(results_df)
  logger::log_info("MR results: {metrics$results} outcomes analysed")


  logger::log_info("5) Build summary plots…")
  manhattan <- manhattan_plot(results_df, Multiple_testing_correction = cfg$mtc)
  volcano   <- volcano_plot(results_df, Multiple_testing_correction = cfg$mtc)

  if (nzchar(cfg$plot_dir)) {
    ggplot2::ggsave(file.path(cfg$plot_dir, "manhattan.png"), manhattan, width = 10, height = 6, dpi = 150)
    ggplot2::ggsave(file.path(cfg$plot_dir, "volcano.png"),   volcano,   width = 10, height = 6, dpi = 150)
  }

  summary_tbl <- tibble::tibble(
    stage = c(
      "outcomes","outcomes_ARD","outcomes_nonARD",
      "exposure_in","exposure_mapped","outcome_snps","results"
    ),
    count = c(
      metrics$outcomes, metrics$n_ard, metrics$n_non,
      metrics$exposure_in, metrics$exposure_mapped, metrics$outcome_snps, metrics$results
    )
  )
  logger::log_info(
    "Summary counts:\n{paste(capture.output(print(summary_tbl)), collapse = '\n')}"
  )

  invisible(list(
    MR_df = MR_df,
    results_df = results_df,
    manhattan = manhattan,
    volcano = volcano
  ))
}


