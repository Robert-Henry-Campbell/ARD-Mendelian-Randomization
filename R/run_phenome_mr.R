#R/run_phenome.MR
#' Run phenome-wide MR across phenotypes (orchestrator)
#'
#' Builds the outcome plan, maps exposure SNPs, fetches outcome SNP rows,
#' runs MR + sensitivity/QC, and returns results + plots. Phenotypes flagged by
#' `ARD_selected` denote age-related diseases.
#'
#' @param exposure Character (mandatory) label for the exposure; used for
#'   naming output directories and plot titles.
#' @param exposure_snps Data frame of exposure instruments (TwoSampleMR-like:
#'   rsid, beta, se, effect_allele, other_allele, eaf, etc.).
#' @param ancestry Character (mandatory), e.g. "EUR".
#' @param sex One of "both","male","female". If not "both", you likely use Neale.
#' @param sensitivity_enabled Character vector of checks (default = all 8).
#' @param sensitivity_pass_min Integer threshold to pass QC (default 6).
#' @param Multiple_testing_correction "BH" or "bonferroni" (default "BH").
#' @param scatterplot,snpforestplot,leaveoneoutplot Logical; produce per-outcome diagnostics.
#' @param Neale_GWAS_dir Optional path to Neale sumstats; created if missing.
#' @param cache_dir Cache directory for temporary files (default:
#'   [ardmr_cache_dir()]).
#' @param logfile Optional path to log file. Defaults to a timestamped file in
#'   `cache_dir` if not supplied.
#' @param verbose Logical; if `FALSE`, only warnings/errors are logged.
#' @param confirm Character; whether to prompt before downloading required
#'   resources when using Neale data. One of "ask", "yes", or "no".
#' @param force_refresh Logical; if `TRUE`, re-download remote resources even if
#'   they are already cached.
#'
#' @return A list: MR_df, results_df, manhattan (ggplot), volcano (ggplot)
#' @export
run_phenome_mr <- function(
    exposure,
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
  if (missing(exposure)) stop("`exposure` is mandatory.")
  exposure <- as.character(exposure)
  if (!length(exposure)) stop("`exposure` is mandatory.")
  exposure <- exposure[1]
  if (is.na(exposure)) stop("`exposure` is mandatory.")
  exposure <- trimws(exposure)
  if (!nzchar(exposure)) stop("`exposure` is mandatory.")
  if (missing(ancestry) || !nzchar(ancestry)) stop("`ancestry` is mandatory.")
  assert_exposure(exposure_snps)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Always derive standard subfolders inside cache_dir
  .slug <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "NA"
    x <- gsub("[^[:alnum:]]+", "-", x)
    x <- gsub("-+", "-", x)
    x <- gsub("^-|-$", "", x)
    tolower(x)
  }

  plot_output_dir <- file.path(cache_dir, "output", .slug(exposure), sex, ancestry)
  Neale_GWAS_dir  <- file.path(cache_dir, "neale_sumstats")
  dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(Neale_GWAS_dir,  recursive = TRUE, showWarnings = FALSE)


  # ---- choose catalog ----
  catalog <- if (sex == "both") {
    "panukb"
  } else {
    "neale"
  }

  # ensure logs live alongside other run outputs
  log_dir <- plot_output_dir


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
    exposure = exposure,
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
    MR_df <- panukb_snp_grabber(
      exposure_snps2,
      MR_df,
      ancestry = cfg$ancestry,
      cache_dir = cfg$cache_dir,
      verbose = cfg$verbose,
      force_refresh = cfg$force_refresh
    )
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
  results_ard_only <- if ("ARD_selected" %in% names(results_df)) {
    results_df[!is.na(results_df$ARD_selected) & results_df$ARD_selected, , drop = FALSE]
  } else {
    results_df[FALSE, , drop = FALSE]
  }

  manhattan_BH_all   <- manhattan_plot(results_df,        Multiple_testing_correction = "BH",         exposure = exposure)
  manhattan_BH_ARD   <- manhattan_plot(results_ard_only,  Multiple_testing_correction = "BH",         exposure = exposure)
  manhattan_Bonf_all <- manhattan_plot(results_df,        Multiple_testing_correction = "bonferroni", exposure = exposure)
  manhattan_Bonf_ARD <- manhattan_plot(results_ard_only,  Multiple_testing_correction = "bonferroni", exposure = exposure)

  # ---- 5B. VOLCANO ----
  volcano_default <- volcano_plot(results_df, Multiple_testing_correction = cfg$mtc)

  # ---- 6) SES enrichment analyses ----
  logger::log_info("6) Enrichment analyses…")
  enrich <- run_enrichment(
    results_df,
    exposure = exposure,
    levels = c("cause_level_1","cause_level_2","cause_level_3"),
    modes  = c("ARD_vs_nonARD_within_cause","cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere"),
    use_qc_pass = TRUE,
    min_nsnp    = 2,
    weight_scheme = "inv_se2",
    exact_max_combn = 1e5,
    mc_B = 100000,
    seed = 1
  )

  # ---- 6A. SES enrichment plots ----
  enrichment_global_plot_dir    <- plot_enrichment_global(enrich$global_tbl)
  enrichment_global_plot_signed <- plot_enrichment_global_signed(enrich$global_tbl)

  cause_levels  <- c("cause_level_1","cause_level_2","cause_level_3")
  compare_modes <- c("ARD_vs_nonARD_within_cause","cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere")

  enrichment_cause_plots_dir <- lapply(cause_levels, function(lv) {
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
  names(enrichment_cause_plots_dir) <- cause_levels

  enrichment_cause_plots_signed <- lapply(cause_levels, function(lv) {
    lv_list <- lapply(compare_modes, function(md) {
      plot_enrichment_signed_forest(
        by_cause_tbl = enrich$by_cause_tbl,
        level = lv,
        compare_mode = md
      )
    })
    names(lv_list) <- compare_modes
    lv_list
  })
  names(enrichment_cause_plots_signed) <- cause_levels

  # ---- 6C) β-scale contrasts (Δβ) analyses + plots ----
  logger::log_info("6C) Beta-scale contrasts (Δβ)…")

  # Containers for tables & plots
  beta_contrast_tables <- list()
  beta_contrast_plots  <- list()

  # Global ARD vs non-ARD
  tbl_global_bc <- beta_contrast_global_ARD(
    results_df,
    use_qc_pass = TRUE, min_nsnp = 2,
    exposure = exposure,
    Multiple_testing_correction = cfg$mtc, alpha = 0.05
  )
  beta_contrast_tables$global <- list(ARD_vs_nonARD = tbl_global_bc)
  beta_contrast_plots$global  <- list(
    ARD_vs_nonARD = plot_beta_contrast_forest(
      tbl_global_bc,
      title = sprintf("Global Δβ: ARD vs non-ARD (%s)", exposure),
      Multiple_testing_correction = cfg$mtc, alpha = 0.05
    )
  )

  # Cause-levels × modes: cause_vs_rest_all + ARD-based modes
  for (lv in cause_levels) {
    beta_contrast_tables[[lv]] <- list()
    beta_contrast_plots[[lv]]  <- list()

    for (md in compare_modes) {
      tbl <- if (md == "cause_vs_rest_all") {
        # all phenotypes: inside cause vs outside
        beta_contrast_by_cause(
          results_df,
          level = lv,
          use_qc_pass = TRUE, min_nsnp = 2,
          exposure = exposure,
          Multiple_testing_correction = cfg$mtc, alpha = 0.05
        )
      } else {
        # ARD-based scopes mirroring SES contexts
        beta_contrast_by_cause_mode(
          results_df,
          level = lv, compare_mode = md,
          use_qc_pass = TRUE, min_nsnp = 2,
          exposure = exposure,
          Multiple_testing_correction = cfg$mtc, alpha = 0.05
        )
      }

      beta_contrast_tables[[lv]][[md]] <- tbl
      beta_contrast_plots[[lv]][[md]]  <- plot_beta_contrast_forest(
        tbl,
        title = sprintf("Δβ by %s — %s (%s)",
                        gsub("_"," ", lv),
                        .pretty_compare(md),
                        exposure),
        Multiple_testing_correction = cfg$mtc, alpha = 0.05
      )
    }
  }

  # ---- 6D) β-scale IVW means ("beta" analysis) ----
  logger::log_info("6D) Beta analysis (IVW mean MR β)…")

  beta_tables <- list()
  beta_plots  <- list()

  tbl_beta_global <- beta_mean_global(
    results_df,
    use_qc_pass = TRUE, min_nsnp = 2,
    exposure = exposure
  )
  beta_tables$global <- tbl_beta_global
  beta_plots$global  <- list(
    mean_effect = plot_beta_mean_global(
      tbl_beta_global,
      title = sprintf("Global Mean Effect of %s", exposure)
    )
  )

  pretty_level <- function(lv) {
    switch(lv,
           cause_level_1 = "Cause Level 1",
           cause_level_2 = "Cause Level 2",
           cause_level_3 = "Cause Level 3",
           gsub("_", " ", lv)
    )
  }

  for (lv in cause_levels) {
    tbl_all <- beta_mean_by_cause(
      results_df,
      level = lv,
      use_qc_pass = TRUE, min_nsnp = 2,
      exposure = exposure,
      ard_only = FALSE,
      drop_empty = TRUE
    )
    tbl_ard <- beta_mean_by_cause(
      results_df,
      level = lv,
      use_qc_pass = TRUE, min_nsnp = 2,
      exposure = exposure,
      ard_only = TRUE,
      drop_empty = TRUE
    )

    beta_tables[[lv]] <- list(
      all_diseases = tbl_all,
      age_related_diseases = tbl_ard
    )

    beta_plots[[lv]] <- list(
      all_diseases = plot_beta_mean_forest(
        tbl_all,
        title = sprintf("Mean effect of %s on all disease by %s", exposure, pretty_level(lv))
      ),
      age_related_diseases = plot_beta_mean_forest(
        tbl_ard,
        title = sprintf("Mean effect of %s on ARDs by %s", exposure, pretty_level(lv))
      )
    )
  }

  # ---- 7) Assemble hierarchical summary_plots list ----
  summary_plots <- list(
    manhattan = list(
      BH = list(all = manhattan_BH_all, ARD_only = manhattan_BH_ARD),
      bonferroni = list(all = manhattan_Bonf_all, ARD_only = manhattan_Bonf_ARD)
    ),
    volcano = list(default = volcano_default),
    enrichment = list(
      global = list(
        directional = list(ARD_vs_nonARD = enrichment_global_plot_dir),
        signed      = list(ARD_vs_nonARD = enrichment_global_plot_signed)
      ),
      cause_level_1 = list(
        directional = enrichment_cause_plots_dir[["cause_level_1"]],
        signed      = enrichment_cause_plots_signed[["cause_level_1"]]
      ),
      cause_level_2 = list(
        directional = enrichment_cause_plots_dir[["cause_level_2"]],
        signed      = enrichment_cause_plots_signed[["cause_level_2"]]
      ),
      cause_level_3 = list(
        directional = enrichment_cause_plots_dir[["cause_level_3"]],
        signed      = enrichment_cause_plots_signed[["cause_level_3"]]
      )
    ),
    beta = beta_plots,
    beta_contrast = beta_contrast_plots
  )

  # Assert counts
  n_enrich_cause_plots <- sum(vapply(
    summary_plots$enrichment[c("cause_level_1","cause_level_2","cause_level_3")],
    function(l) if (!is.null(l$directional)) length(l$directional) else 0L,
    integer(1)
  ))
  logger::log_info("Enrichment cause-level plots generated: {n_enrich_cause_plots}")

  n_beta_cause_plots <- sum(vapply(
    summary_plots$beta_contrast[cause_levels],
    length,
    integer(1)
  ))
  logger::log_info("β-contrast cause-level plots generated: {n_beta_cause_plots}")

  n_beta_mean_plots <- sum(vapply(
    summary_plots$beta[cause_levels],
    function(l) if (is.list(l)) length(l) else 0L,
    integer(1)
  ))
  logger::log_info("Beta mean cause-level plots generated: {n_beta_mean_plots}")

  # ---- 8) Save plots mirroring the list structure under cfg$plot_dir ----
  save_plot_hierarchy <- function(x, base_dir, path_parts = character(0)) {
    safe_name <- function(s) {
      s <- as.character(s)
      s <- gsub("[/\\?%*:|\"<>]", "_", s)
      s <- gsub("\\s+", "_", s)
      s <- gsub("_+", "_", s)
      s <- sub("^_+", "", s)
      s <- sub("_+$", "", s)
      if (!nzchar(s)) "plot" else s
    }

    if (inherits(x, "ggplot")) {
      subpath <- paste(vapply(path_parts, safe_name, character(1)), collapse = "/")
      width <- 6.5; height <- 6.5

      if (grepl("^manhattan", subpath)) {
        width <- 6.5; height <- 6.5
      }
      if (grepl("^enrichment/global", subpath)) {
        width <- 3.54; height <- 2.4
      }
      if (grepl("^enrichment/cause_level_", subpath)) {
        width <- 7.2; height <- 5.0
      }
      if (grepl("^volcano", subpath)) {
        width <- 7.2; height <- 5.0
      }
      # NEW: sizes for beta_contrast
      if (grepl("^beta/global", subpath)) {
        width <- 3.54; height <- 2.4
      }
      if (grepl("^beta/cause_level_", subpath)) {
        width <- 7.2; height <- 5.0
      }
      if (grepl("^beta_contrast/global", subpath)) {
        width <- 3.54; height <- 2.4
      }
      if (grepl("^beta_contrast/cause_level_", subpath)) {
        width <- 7.2; height <- 5.0
      }

      dir_path  <- file.path(base_dir, dirname(subpath))
      if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      file_stem <- safe_name(basename(subpath))
      file_name <- paste0(file_stem, ".png")
      file_path <- file.path(dir_path, file_name)
      ggplot2::ggsave(filename = file_path, plot = x, width = width, height = height, dpi = 300)
      .ardmr_write_plot_data(x, dir_path = dir_path, base_name = file_stem)
      return(invisible(NULL))
    }

    if (is.list(x)) {
      nms <- names(x)
      for (i in seq_along(x)) {
        nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0("item", i)
        save_plot_hierarchy(x[[i]], base_dir, c(path_parts, nm))
      }
      return(invisible(NULL))
    }

    invisible(NULL)
  }

  # Save all plots
  save_plot_hierarchy(summary_plots, cfg$plot_dir)

  # ---- 8b) Export enrichment tables (CSV) + beta-contrast tables ----
  write_enrichment_tables <- function(enrich, base_dir) {
    out_dir <- file.path(base_dir, "enrichment", "tables")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    global_csv   <- file.path(out_dir, "global.csv")
    by_cause_csv <- file.path(out_dir, "by_cause.csv")
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_csv(enrich$global_tbl,   global_csv)
      readr::write_csv(enrich$by_cause_tbl, by_cause_csv)
    } else {
      utils::write.csv(enrich$global_tbl,   global_csv,   row.names = FALSE)
      utils::write.csv(enrich$by_cause_tbl, by_cause_csv, row.names = FALSE)
    }
    logger::log_info("Enrichment tables written to {out_dir}")
  }

  write_beta_contrast_tables <- function(beta_tables, base_dir) {
    out_dir <- file.path(base_dir, "beta_contrast", "tables")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    safe_name <- function(s) {
      s <- as.character(s)
      s <- gsub("[/\\?%*:|\"<>]", "_", s)
      s <- gsub("\\s+", "_", s)
      s <- gsub("_+", "_", s)
      s <- sub("^_+", "", s)
      s <- sub("_+$", "", s)
      if (!nzchar(s)) "table" else s
    }
    wcsv <- function(df, path) {
      if (!nrow(df)) return(invisible(NULL))
      if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(df, path)
      else utils::write.csv(df, path, row.names = FALSE)
    }
    # Global
    if (!is.null(beta_tables$global$ARD_vs_nonARD)) {
      wcsv(beta_tables$global$ARD_vs_nonARD, file.path(out_dir, "global_ARD_vs_nonARD.csv"))
    }
    # Cause levels × modes
    clv <- c("cause_level_1","cause_level_2","cause_level_3")
    md  <- c("cause_vs_rest_all","ARD_vs_nonARD_within_cause","ARD_in_cause_vs_ARD_elsewhere")
    for (lv in clv) {
      if (!lv %in% names(beta_tables)) next
      for (m in md) {
        if (!m %in% names(beta_tables[[lv]])) next
        fn <- sprintf("%s__%s.csv", safe_name(lv), safe_name(m))
        wcsv(beta_tables[[lv]][[m]], file.path(out_dir, fn))
      }
    }
    logger::log_info("β-contrast tables written to {out_dir}")
  }

  write_beta_tables <- function(beta_tables, base_dir) {
    out_dir <- file.path(base_dir, "beta", "tables")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    safe_name <- function(s) {
      s <- as.character(s)
      s <- gsub("[/\\?%*:|\"<>]", "_", s)
      s <- gsub("\\s+", "_", s)
      s <- gsub("_+", "_", s)
      s <- sub("^_+", "", s)
      s <- sub("_+$", "", s)
      if (!nzchar(s)) "table" else s
    }
    wcsv <- function(df, path) {
      if (!nrow(df)) return(invisible(NULL))
      if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(df, path)
      else utils::write.csv(df, path, row.names = FALSE)
    }

    if (!is.null(beta_tables$global)) {
      wcsv(beta_tables$global, file.path(out_dir, "global.csv"))
    }

    clv <- c("cause_level_1","cause_level_2","cause_level_3")
    scopes <- c("all_diseases", "age_related_diseases")
    for (lv in clv) {
      if (!lv %in% names(beta_tables)) next
      for (sc in scopes) {
        if (!sc %in% names(beta_tables[[lv]])) next
        tbl <- beta_tables[[lv]][[sc]]
        if (!nrow(tbl)) next
        fn <- sprintf("%s__%s.csv", safe_name(lv), safe_name(sc))
        wcsv(tbl, file.path(out_dir, fn))
      }
    }

    logger::log_info("Beta tables written to {out_dir}")
  }

  write_enrichment_tables(enrich, cfg$plot_dir)
  write_beta_contrast_tables(beta_contrast_tables, cfg$plot_dir)
  write_beta_tables(beta_tables, cfg$plot_dir)

  # ---- 9) Keep the originals around too (optional) ----
  manhattan <- manhattan_BH_all
  volcano   <- volcano_default

  # ---- 10) Summary counts (unchanged placeholder) ----
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

  output <- list(
    MR_df = MR_df,
    results_df    = results_df,
    summary_plots = summary_plots,
    enrich        = enrich,
    beta          = beta_tables,
    beta_contrast = beta_contrast_tables
  )

  saveRDS(output, file = file.path(cfg$plot_dir, "results.rds"))
  output
}
