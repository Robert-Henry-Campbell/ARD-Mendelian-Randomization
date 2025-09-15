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
    verbose = TRUE,
    confirm = 'ask',
    force_refresh = FALSE
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

  # ensure logs/ exists inside cache_dir
  log_dir <- file.path(cache_dir, "logs")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)


  logfile <- if (is.null(logfile) || !nzchar(logfile)) {
    file.path(log_dir,sprintf("ardmr_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
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
    verbose = verbose,
    confirm = confirm,
    force_refresh = force_refresh
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
    neale_gwas_checker(MR_df, neale_dir = cfg$neale_dir, verbose = cfg$verbose, confirm = cfg$confirm)
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

  # ---- 5A. MANHATTAN: BH vs Bonf × (all vs ARD-only) ----
  # Helper: ARD-only slice (logical column)
  results_ard_only <- if ("ARD_selected" %in% names(results_df)) {
    results_df[!is.na(results_df$ARD_selected) & results_df$ARD_selected, , drop = FALSE]
  } else {
    results_df[FALSE, , drop = FALSE]
  }

  manhattan_BH_all  <- manhattan_plot(results_df, Multiple_testing_correction = "BH",         exposure = exposure_snps$id.exposure)
  manhattan_BH_ARD  <- manhattan_plot(results_ard_only, Multiple_testing_correction = "BH",   exposure = exposure_snps$id.exposure)
  manhattan_Bonf_all<- manhattan_plot(results_df, Multiple_testing_correction = "bonferroni", exposure = exposure_snps$id.exposure)
  manhattan_Bonf_ARD<- manhattan_plot(results_ard_only, Multiple_testing_correction = "bonferroni", exposure = exposure_snps$id.exposure)

  # ---- 5B. VOLCANO: (placeholder/default for now) ----
  volcano_default <- volcano_plot(results_df, Multiple_testing_correction = cfg$mtc)

  # ---- 6) Enrichment analyses ----
  logger::log_info("6) Enrichment analyses…")
  enrich <- run_enrichment(
    results_df,
    exposure = exposure_snps$id.exposure,
    levels = c("cause_level_1","cause_level_2","cause_level_3"),
    modes  = c("ARD_vs_nonARD_within_cause","cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere"),
    use_qc_pass = TRUE,
    min_nsnp    = 2,
    weight_scheme = "inv_se2",
    exact_max_combn = 1e5,
    mc_B = 100000,
    seed = 1
  )

  # ---- 6A. ENRICHMENT PLOTS ----
  # Uses plotting helpers defined in R/plots.R:
  #   plot_enrichment_global()
  #   plot_enrichment_directional_forest()
  enrichment_global_plot <- plot_enrichment_global(enrich$global_tbl)

  # Per-level × mode forest plots
  cause_levels <- c("cause_level_1","cause_level_2","cause_level_3")
  compare_modes <- c("ARD_vs_nonARD_within_cause","cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere")

  enrichment_cause_plots <- lapply(cause_levels, function(lv) {
    lv_list <- lapply(compare_modes, function(md) {
      plot_enrichment_directional_forest(
        by_cause_tbl = enrich$by_cause_tbl,
        level = lv,
        compare_mode = md
      )
    })
    names(lv_list) <- compare_modes
    lv_list
  })
  names(enrichment_cause_plots) <- cause_levels

  # ---- 7) Assemble hierarchical summary_plots list ----
  summary_plots <- list(
    manhattan = list(
      BH = list(
        all      = manhattan_BH_all,
        ARD_only = manhattan_BH_ARD
      ),
      bonferroni = list(
        all      = manhattan_Bonf_all,
        ARD_only = manhattan_Bonf_ARD
      )
    ),
    volcano = list(
      default = volcano_default
    ),
    enrichment = list(
      global = list(ARD_vs_nonARD = enrichment_global_plot),
      cause_level_1 = enrichment_cause_plots[["cause_level_1"]],
      cause_level_2 = enrichment_cause_plots[["cause_level_2"]],
      cause_level_3 = enrichment_cause_plots[["cause_level_3"]]
    )
  )
  # Assert we made 9 cause-level plots (3 levels × 3 modes)
  n_cause_plots <- sum(vapply(summary_plots$enrichment[c("cause_level_1","cause_level_2","cause_level_3")],
                              function(l) length(l), integer(1)))
  logger::log_info("Enrichment cause-level plots generated: {n_cause_plots} (expected 9)")


  # ---- 8) Save plots mirroring the list structure under cache_dir/plots ----
  save_plot_hierarchy <- function(x, base_dir, path_parts = character(0)) {
    # path_parts: vector of folder names reflecting position in list
    if (inherits(x, "ggplot")) {
      # choose sizes
      subpath <- paste(path_parts, collapse = "/")
      # Defaults
      width <- 6.5; height <- 6.5
      # Manhattan size
      if (grepl("^manhattan", subpath)) {
        width <- 6.5; height <- 6.5
      }
      # Global enrichment (single-row) ~ Nature single-column
      if (grepl("^enrichment/global", subpath)) {
        width <- 3.54; height <- 2.4
      }
      # Cause-level enrichment ~ Nature double-column
      if (grepl("^enrichment/cause_level_", subpath)) {
        width <- 7.2; height <- 5.0
      }
      # Volcano default to Nature double-column for now
      if (grepl("^volcano", subpath)) {
        width <- 7.2; height <- 5.0
      }

      dir_path <- file.path(base_dir, dirname(subpath))
      if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      file_name <- paste0(basename(subpath), ".png")
      file_path <- file.path(dir_path, file_name)
      ggplot2::ggsave(filename = file_path, plot = x, width = width, height = height, dpi = 300)
      return(invisible(NULL))
    }

    if (is.list(x)) {
      # Recurse into lists
      nms <- names(x)
      for (i in seq_along(x)) {
        nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0("item", i)
        save_plot_hierarchy(x[[i]], base_dir, c(path_parts, nm))
      }
      return(invisible(NULL))
    }

    # Not a plot or list: ignore
    invisible(NULL)
  }

  # Save everything
  save_plot_hierarchy(summary_plots, cfg$plot_dir)

  # ---- 9) Keep the originals around too (optional) ----
  manhattan <- manhattan_BH_all
  volcano   <- volcano_default

  # ---- 10) Summary counts (unchanged) ----
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
  logger::log_info("Summary counts:\n{paste(capture.output(print(summary_tbl)), collapse = '\n')}")

  invisible(list(
    MR_df = MR_df,
    results_df = results_df,
    summary_plots = summary_plots,
    enrich = enrich
  ))

}


