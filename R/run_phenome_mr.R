#R/run_phenome.MR
#' Run phenome-wide MR across phenotypes (orchestrator)
#'
#' Builds the outcome plan, maps exposure SNPs, fetches outcome SNP rows,
#' runs MR + sensitivity/QC, and returns results + plots. Phenotypes flagged by
#' `ARD_selected` denote age-related diseases.
#'
#' @param exposure Character (mandatory) label for the exposure; used for
#'   naming output directories and plot titles.
#' @param exposure_snps Data frame of exposure instruments (TwoSampleMR-like). needs cols:
#'[1] "id.exposure"           "Exposure"                "SNP" (must be rsid format) "Chr"
#'[5] "Pos"                    "effect_allele.exposure" "other_allele.exposure"  "eaf.exposure"
#'[9] "beta.exposure"          "se.exposure"            "palindromic"            "pval.exposure"
#'[13] "mr_keep"                "F"
#' @param exposure_units Character (mandatory) description of the exposure
#'   units for β-scale plots.
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
#'   they are already cached. Defaults to `TRUE`.
#'
#' @return A list: MR_df, results_df, manhattan (ggplot), volcano (ggplot)
#' @export
run_phenome_mr <- function(
    exposure,
    exposure_snps,
    exposure_units,
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
    force_refresh = TRUE
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
  if (missing(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- as.character(exposure_units)
  if (!length(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- exposure_units[1]
  if (is.na(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- trimws(exposure_units)
  if (!nzchar(exposure_units)) stop("`exposure_units` is mandatory.")
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
  purged_plot_dir <- FALSE
  if (isTRUE(force_refresh) && dir.exists(plot_output_dir)) {
    unlink(plot_output_dir, recursive = TRUE, force = TRUE)
    purged_plot_dir <- TRUE
  }
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
  if (isTRUE(purged_plot_dir)) {
    logger::log_info(
      "force_refresh enabled: cleared existing artifacts at {plot_output_dir}"
    )
  }

  # Keep config bundled for helper calls
  cfg <- list(
    ancestry = ancestry,
    exposure = exposure,
    exposure_units = exposure_units,
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

  emit_summary_counts <- function(metrics) {
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
    summary_tbl
  }

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
      path_labels <- vapply(path_parts, safe_name, character(1))
      file_stem_override <- NULL
      if (length(path_labels) >= 1 && identical(path_labels[1], "manhattan_recolor")) {
        path_labels <- c("manhattan", path_labels[-1])
        if (length(path_labels) >= 2) {
          path_labels[length(path_labels)] <- paste0("recolor_", path_labels[length(path_labels)])
        } else {
          path_labels <- c("manhattan", "recolor")
        }
      }
      if (length(path_labels) >= 1 && identical(path_labels[1], "volcano_recolor")) {
        correction <- if (length(path_labels) >= 2) path_labels[2] else "BH"
        slice      <- if (length(path_labels) >= 3) path_labels[3] else "all"
        correction <- safe_name(correction)
        slice      <- safe_name(slice)
        path_labels <- c("volcano", correction, paste0("recolor_", slice))
        file_stem_override <- safe_name(paste0("recolor_", correction, "_", slice))
      }
      subpath <- paste(path_labels, collapse = "/")
      width <- 7.2; height <- 6.5
      plot_data <- attr(x, "ardmr_plot_data", exact = TRUE)
      yfloat_base <- 1.8 #was 1.2
      yfloat_coef <- 0.28

      if (grepl("^manhattan", subpath)) {
        width <- 7.2; height <- 6.5
      }
      if (grepl("^enrichment/global", subpath)) {
        width <- 7.2; height <- 4.0
      }
      if (grepl("^enrichment/cause_level_", subpath)) {
        width <- 7.2; height <- 5.0
      }
      is_enrichment_cause_forest <- length(path_parts) >= 3 &&
        identical(path_parts[1], "enrichment") &&
        grepl("^cause_level_", path_parts[2]) &&
        identical(path_parts[3], "violin_forest")
      if (is_enrichment_cause_forest) {
        observed_df <- NULL
        if (is.list(plot_data) && "observed" %in% names(plot_data)) {
          observed_df <- plot_data$observed
        }
        n_rows <- if (is.data.frame(observed_df)) {
          if ("group_label" %in% names(observed_df)) {
            groups <- unique(observed_df$group_label)
            groups <- groups[!is.na(groups)]
            length(groups)
          } else {
            nrow(observed_df)
          }
        } else {
          0L
        }
        n_rows <- as.integer(n_rows)
        if (!is.finite(n_rows) || n_rows <= 0) n_rows <- 1L
        height <- yfloat_base + (yfloat_coef * 1.3) * n_rows #adjust this coefficient to stretch rows
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
      if (!is.null(file_stem_override)) file_stem <- file_stem_override
      file_name <- paste0(file_stem, ".png")
      file_path <- file.path(dir_path, file_name)
      is_yfloat <- length(path_labels) >= 1 && grepl("_wrap_yfloat$", tail(path_labels, 1))
      if (is_yfloat) {
        main_df <- NULL
        if (is.list(plot_data) && "main" %in% names(plot_data)) {
          main_df <- plot_data$main
        }
        n_rows <- if (is.data.frame(main_df)) nrow(main_df) else 0L
        if (!is.finite(n_rows) || n_rows <= 0) n_rows <- 1L
        height <- yfloat_base + yfloat_coef * n_rows
      }
      if (length(path_labels) >= 1 && identical(path_labels[1], "beta")) {
        return(invisible(NULL))
      }
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
    MR_df <- neale_snp_grabber(
      exposure_snps2,
      MR_df,
      neale_dir = cfg$neale_dir,
      cache_dir = cfg$cache_dir,
      verbose = cfg$verbose,
      force_refresh = cfg$force_refresh
    )
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

  ard_flags <- if ("ARD_selected" %in% names(MR_df)) {
    vapply(MR_df$ARD_selected, isTRUE, logical(1))
  } else {
    logical(0)
  }
  qc_flags <- if ("results_qc_pass" %in% names(MR_df)) {
    vapply(MR_df$results_qc_pass, isTRUE, logical(1))
  } else if (length(ard_flags)) {
    rep(FALSE, length(ard_flags))
  } else {
    logical(0)
  }
  ard_qc_survivors <- ard_flags & qc_flags
  has_ard_qc_survivor <- length(ard_qc_survivors) > 0 && any(ard_qc_survivors)

  if (!has_ard_qc_survivor) {
    n_ard_candidates <- sum(ard_flags)
    warn_msg <- if (n_ard_candidates > 0) {
      sprintf(
        "No age-related diseases passed QC out of %d candidates; skipping downstream analyses.",
        n_ard_candidates
      )
    } else {
      "No age-related diseases were available for analysis; skipping downstream analyses."
    }
    warning(warn_msg, call. = FALSE)
    logger::log_warn("{warn_msg}")

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
      summary_plots = list(),
      enrich        = list(
        global_tbl = tibble::tibble(),
        by_cause_tbl = tibble::tibble()
      ),
      beta          = list()
    )

    saveRDS(output, file = file.path(cfg$plot_dir, "results.rds"))
    logger::log_info("Early exit: wrote results.rds to {cfg$plot_dir} after QC filtering.")
    return(output)
  }

  logger::log_info("5) Build summary plots…")

  # ---- 5A. MANHATTAN: BH vs Bonf × (all vs ARD-only) ----
  results_ard_only <- if ("ARD_selected" %in% names(results_df)) {
    results_df[!is.na(results_df$ARD_selected) & results_df$ARD_selected, , drop = FALSE]
  } else {
    results_df[FALSE, , drop = FALSE]
  }

  # manhattan_BH_all   <- manhattan_plot(results_df,        Multiple_testing_correction = "BH",         exposure = exposure)
  # manhattan_BH_ARD   <- manhattan_plot(results_ard_only,  Multiple_testing_correction = "BH",         exposure = exposure)
  # manhattan_Bonf_all <- manhattan_plot(results_df,        Multiple_testing_correction = "bonferroni", exposure = exposure)
  # manhattan_Bonf_ARD <- manhattan_plot(results_ard_only,  Multiple_testing_correction = "bonferroni", exposure = exposure)

  manhattan_with_names <- manhattan_plot(
    results_df,
    Multiple_testing_correction = cfg$mtc,
    exposure = exposure,
    dot_names = TRUE
  )
  manhattan_without_names <- manhattan_plot(
    results_df,
    Multiple_testing_correction = cfg$mtc,
    exposure = exposure,
    dot_names = FALSE
  )

  manhattan_recolor_BH_all   <- manhattan_plot_recolor(results_df,       Multiple_testing_correction = "BH",         exposure = exposure)
  manhattan_recolor_BH_ARD   <- manhattan_plot_recolor(results_ard_only, Multiple_testing_correction = "BH",         exposure = exposure)
  manhattan_recolor_Bonf_all <- manhattan_plot_recolor(results_df,       Multiple_testing_correction = "bonferroni", exposure = exposure)
  manhattan_recolor_Bonf_ARD <- manhattan_plot_recolor(results_ard_only, Multiple_testing_correction = "bonferroni", exposure = exposure)

  # ---- 5B. VOLCANO ----
  volcano_with_dot_names <- volcano_plot(
    results_df,
    Multiple_testing_correction = cfg$mtc,
    exposure = exposure,
    verbose = cfg$verbose,
    dot_names = TRUE
  )
  volcano_without_dot_names <- volcano_plot(
    results_df,
    Multiple_testing_correction = cfg$mtc,
    exposure = exposure,
    verbose = cfg$verbose,
    dot_names = FALSE
  )
  # volcano_default <- volcano_plot(results_df, Multiple_testing_correction = cfg$mtc)
  volcano_recolor_BH_all   <- volcano_plot_recolor(results_df,       Multiple_testing_correction = "BH",         exposure = exposure)
  volcano_recolor_BH_ARD   <- volcano_plot_recolor(results_ard_only, Multiple_testing_correction = "BH",         exposure = exposure)
  volcano_recolor_Bonf_all <- volcano_plot_recolor(results_df,       Multiple_testing_correction = "bonferroni", exposure = exposure)
  volcano_recolor_Bonf_ARD <- volcano_plot_recolor(results_ard_only, Multiple_testing_correction = "bonferroni", exposure = exposure)

  summary_plots_core <- list(
    manhattan = list(
      with_names = manhattan_with_names,
      without_names = manhattan_without_names
    ),
    manhattan_recolor = list(
      BH = list(all = manhattan_recolor_BH_all, ARD_only = manhattan_recolor_BH_ARD),
      bonferroni = list(all = manhattan_recolor_Bonf_all, ARD_only = manhattan_recolor_Bonf_ARD)
    ),
    volcano = list(
      with_dot_names = volcano_with_dot_names,
      without_dot_names = volcano_without_dot_names
    ),
    volcano_recolor = list(
      BH = list(all = volcano_recolor_BH_all, ARD_only = volcano_recolor_BH_ARD),
      bonferroni = list(all = volcano_recolor_Bonf_all, ARD_only = volcano_recolor_Bonf_ARD)
    )
  )

  # ---- 6) Signed enrichment analyses ----
  logger::log_info("6) Enrichment analyses…")
  partial_output <- NULL
  enrichment_failed <- FALSE
  enrich <- tryCatch(
    run_enrichment(
      results_df,
      exposure = exposure,
      levels = c("cause_level_1","cause_level_2","cause_level_3"),
      modes  = c("cause_vs_rest_all"),
      use_qc_pass = TRUE,
      min_nsnp    = 2,
      weight_scheme = "inv_se2",
      exact_max_combn = 1e5,
      mc_B = 10000,
      seed = 1,
      retain_permutations = TRUE,
      retain_perm_max = 10000
    ),
    error = function(e) {
      warn_msg <- sprintf("Enrichment failed: %s", conditionMessage(e))
      warning(warn_msg, call. = FALSE)
      logger::log_warn(warn_msg)
      summary_plots <- summary_plots_core
      save_plot_hierarchy(summary_plots, cfg$plot_dir)
      emit_summary_counts(metrics)
      output <- list(
        MR_df = MR_df,
        results_df    = results_df,
        summary_plots = summary_plots,
        enrich        = NULL,
        beta          = NULL
      )
      saveRDS(output, file = file.path(cfg$plot_dir, "results.rds"))
      logger::log_info(
        "Early exit: wrote results.rds to {cfg$plot_dir} after enrichment failure."
      )
      partial_output <<- output
      enrichment_failed <<- TRUE
      NULL
    }
  )
  if (isTRUE(enrichment_failed)) {
    return(partial_output)
  }

  logger::log_info("6A) Enrichment violin plots…")
  cause_levels  <- c("cause_level_1","cause_level_2","cause_level_3")
  compare_mode  <- "cause_vs_rest_all"

  enrichment_global_violin_vertical <- plot_enrichment_signed_violin_global(
    enrich$global_tbl,
    orientation = "vertical",
    alpha = 0.05
  )
  enrichment_global_violin_forest <- plot_enrichment_signed_violin_global(
    enrich$global_tbl,
    orientation = "forest",
    alpha = 0.05
  )

  enrichment_cause_violin <- lapply(cause_levels, function(lv) {
    list(
      violin_vertical = plot_enrichment_signed_violin_by_cause(
        enrich$by_cause_tbl,
        level = lv,
        compare_mode = compare_mode,
        orientation = "vertical",
        alpha = 0.05
      ),
      violin_forest = plot_enrichment_signed_violin_by_cause(
        enrich$by_cause_tbl,
        level = lv,
        compare_mode = compare_mode,
        orientation = "forest",
        alpha = 0.05
      )
    )
  })
  names(enrichment_cause_violin) <- cause_levels

  n_enrichment_cause_plots <- sum(vapply(enrichment_cause_violin, length, integer(1)))
  logger::log_info("Enrichment cause-level plots generated: {n_enrichment_cause_plots}")

  # ---- 6C) β-scale contrasts (Δβ) analyses + plots ----
  # logger::log_info("6C) Beta-scale contrasts (Δβ)…")

  # Containers for tables & plots
  # beta_contrast_tables <- list()
  # beta_contrast_plots  <- list()

  # Global ARD vs non-ARD
  # tbl_global_bc <- beta_contrast_global_ARD(
  #   results_df,
  #   use_qc_pass = TRUE, min_nsnp = 2,
  #   exposure = exposure,
  #   Multiple_testing_correction = cfg$mtc, alpha = 0.05
  # )
  # beta_contrast_tables$global <- list(ARD_vs_nonARD = tbl_global_bc)
  # beta_contrast_plots$global  <- list(
  #   ARD_vs_nonARD = plot_beta_contrast_forest(
  #     tbl_global_bc,
  #     title = sprintf("Global Δβ: ARD vs non-ARD (%s)", exposure),
  #     Multiple_testing_correction = cfg$mtc, alpha = 0.05
  #   )
  # )

  # Cause-levels × modes: cause_vs_rest_all + ARD-based modes
  # for (lv in cause_levels) {
  #   beta_contrast_tables[[lv]] <- list()
  #   beta_contrast_plots[[lv]]  <- list()

  #   for (md in compare_modes) {
  #     tbl <- if (md == "cause_vs_rest_all") {
  #       # all phenotypes: inside cause vs outside
  #       beta_contrast_by_cause(
  #         results_df,
  #         level = lv,
  #         use_qc_pass = TRUE, min_nsnp = 2,
  #         exposure = exposure,
  #         Multiple_testing_correction = cfg$mtc, alpha = 0.05
  #       )
  #     } else {
  #       # ARD-based scopes mirroring SES contexts
  #       beta_contrast_by_cause_mode(
  #         results_df,
  #         level = lv, compare_mode = md,
  #         use_qc_pass = TRUE, min_nsnp = 2,
  #         exposure = exposure,
  #         Multiple_testing_correction = cfg$mtc, alpha = 0.05
  #       )
  #     }

  #     beta_contrast_tables[[lv]][[md]] <- tbl
  #     beta_contrast_plots[[lv]][[md]]  <- plot_beta_contrast_forest(
  #       tbl,
  #       title = sprintf("Δβ by %s — %s (%s)",
  #                       gsub("_"," ", lv),
  #                       .pretty_compare(md),
  #                       exposure),
  #       Multiple_testing_correction = cfg$mtc, alpha = 0.05
  #     )
  #     beta_contrast_plots[[lv]][[paste0(md, "_wrap")]] <- plot_beta_contrast_forest_wrap(
  #       tbl,
  #       title = sprintf("Δβ by %s — %s (%s)",
  #                       gsub("_"," ", lv),
  #                       .pretty_compare(md),
  #                       exposure),
  #       Multiple_testing_correction = cfg$mtc, alpha = 0.05
  #     )
  #   }
  # }

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

  beta_effect_scale <- if (identical(cfg$sex, "both")) "odds_ratio" else "absolute_risk"

  beta_plots$global  <- list(
   mean_effect = plot_beta_mean_global(
     tbl_beta_global,
     title = sprintf("Global Mean Effect of %s", exposure),
     exposure_units = exposure_units,
     effect_scale = beta_effect_scale
   )
  )


  pretty_level <- function(lv) {
    pretty <- switch(
      lv,
      cause_level_1 = "cause level 1",
      cause_level_2 = "cause level 2",
      cause_level_3 = "cause level 3",
      gsub("_", " ", lv, fixed = TRUE)
    )
    tolower(pretty)
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
        title = sprintf("Mean effect of %s on all disease by %s", exposure, pretty_level(lv)),
        exposure_units = exposure_units,
        effect_scale = beta_effect_scale
      ),
      age_related_diseases = plot_beta_mean_forest(
        tbl_ard,
        title = sprintf("Mean effect of %s on ARDs by %s", exposure, pretty_level(lv)),
        exposure_units = exposure_units,
        effect_scale = beta_effect_scale
      )
    )

    beta_plots[[lv]][["all_diseases_wrap"]] <- plot_beta_mean_forest_wrap(
      tbl_all,
      title = sprintf("Mean effect of %s on all disease by %s", exposure, pretty_level(lv)),
      exposure_units = exposure_units,
      effect_scale = beta_effect_scale
    )
    beta_plots[[lv]][["all_diseases_wrap_yfloat"]] <- plot_beta_mean_forest_wrap_yfloat(
      tbl_all,
      title = sprintf("Mean effect of %s on all disease by %s", exposure, pretty_level(lv)),
      exposure_units = exposure_units,
      effect_scale = beta_effect_scale
    )
    beta_plots[[lv]][["age_related_diseases_wrap"]] <- plot_beta_mean_forest_wrap(
      tbl_ard,
      title = sprintf("Mean effect of %s on ARDs by %s", exposure, pretty_level(lv)),
      exposure_units = exposure_units,
      effect_scale = beta_effect_scale
    )
    beta_plots[[lv]][["age_related_diseases_wrap_yfloat"]] <- plot_beta_mean_forest_wrap_yfloat(
      tbl_ard,
      title = sprintf("Mean effect of %s on ARDs by %s", exposure, pretty_level(lv)),
      exposure_units = exposure_units,
      effect_scale = beta_effect_scale
    )
  }

  # ---- 7) Assemble hierarchical summary_plots list ----
  summary_plots <- c(
    summary_plots_core,
    list(
      enrichment = list(
        global = list(
          violin_vertical = enrichment_global_violin_vertical,
          violin_forest = enrichment_global_violin_forest
        ),
        cause_level_1 = enrichment_cause_violin[["cause_level_1"]],
        cause_level_2 = enrichment_cause_violin[["cause_level_2"]],
        cause_level_3 = enrichment_cause_violin[["cause_level_3"]]
      ),
      # enrichment = list(
      #   global = list(
      #     directional = list(ARD_vs_nonARD = enrichment_global_plot_dir),
      #     signed      = list(ARD_vs_nonARD = enrichment_global_plot_signed)
      #   ),
      #   cause_level_1 = list(
      #     directional = enrichment_cause_plots_dir[["cause_level_1"]],
      #     signed      = enrichment_cause_plots_signed[["cause_level_1"]]
      #   ),
      #   cause_level_2 = list(
      #     directional = enrichment_cause_plots_dir[["cause_level_2"]],
      #     signed      = enrichment_cause_plots_signed[["cause_level_2"]]
      #   ),
      #   cause_level_3 = list(
      #     directional = enrichment_cause_plots_dir[["cause_level_3"]],
      #     signed      = enrichment_cause_plots_signed[["cause_level_3"]]
      #   )
      # ),
      beta = beta_plots
      # , beta_contrast = beta_contrast_plots
    )
  )

  # Assert counts
  # n_enrich_cause_plots <- sum(vapply(
  #   summary_plots$enrichment[c("cause_level_1","cause_level_2","cause_level_3")],
  #   function(l) if (!is.null(l$directional)) length(l$directional) else 0L,
  #   integer(1)
  # ))
  # logger::log_info("Enrichment cause-level plots generated: {n_enrich_cause_plots}")

  # n_beta_cause_plots <- sum(vapply(
  #   summary_plots$beta_contrast[cause_levels],
  #   length,
  #   integer(1)
  # ))
  # logger::log_info("β-contrast cause-level plots generated: {n_beta_cause_plots}")

  n_beta_mean_plots <- sum(vapply(
    summary_plots$beta[cause_levels],
    function(l) if (is.list(l)) length(l) else 0L,
    integer(1)
  ))
  logger::log_info("Beta mean cause-level plots generated: {n_beta_mean_plots}")

  # ---- 8) Save plots mirroring the list structure under cfg$plot_dir ----

  # Save all plots
  save_plot_hierarchy(summary_plots, cfg$plot_dir)

  # ---- 8b) Export enrichment tables (CSV) + beta-contrast tables ----
  # write_enrichment_tables <- function(enrich, base_dir) {
  #   out_dir <- file.path(base_dir, "enrichment", "tables")
  #   if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  #   global_csv   <- file.path(out_dir, "global.csv")
  #   by_cause_csv <- file.path(out_dir, "by_cause.csv")
  #   if (requireNamespace("readr", quietly = TRUE)) {
  #     readr::write_csv(enrich$global_tbl,   global_csv)
  #     readr::write_csv(enrich$by_cause_tbl, by_cause_csv)
  #   } else {
  #     utils::write.csv(enrich$global_tbl,   global_csv,   row.names = FALSE)
  #     utils::write.csv(enrich$by_cause_tbl, by_cause_csv, row.names = FALSE)
  #   }
  #   logger::log_info("Enrichment tables written to {out_dir}")
  # }

  # write_beta_contrast_tables <- function(beta_tables, base_dir) {
  #   out_dir <- file.path(base_dir, "beta_contrast", "tables")
  #   if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  #   safe_name <- function(s) {
  #     s <- as.character(s)
  #     s <- gsub("[/\\?%*:|\"<>]", "_", s)
  #     s <- gsub("\\s+", "_", s)
  #     s <- gsub("_+", "_", s)
  #     s <- sub("^_+", "", s)
  #     s <- sub("_+$", "", s)
  #     if (!nzchar(s)) "table" else s
  #   }
  #   wcsv <- function(df, path) {
  #     if (!nrow(df)) return(invisible(NULL))
  #     if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(df, path)
  #     else utils::write.csv(df, path, row.names = FALSE)
  #   }
  #   # Global
  #   if (!is.null(beta_tables$global$ARD_vs_nonARD)) {
  #     wcsv(beta_tables$global$ARD_vs_nonARD, file.path(out_dir, "global_ARD_vs_nonARD.csv"))
  #   }
  #   # Cause levels × modes
  #   clv <- c("cause_level_1","cause_level_2","cause_level_3")
  #   md  <- c("cause_vs_rest_all","ARD_vs_nonARD_within_cause","ARD_in_cause_vs_ARD_elsewhere")
  #   for (lv in clv) {
  #     if (!lv %in% names(beta_tables)) next
  #     for (m in md) {
  #       if (!m %in% names(beta_tables[[lv]])) next
  #       fn <- sprintf("%s__%s.csv", safe_name(lv), safe_name(m))
  #       wcsv(beta_tables[[lv]][[m]], file.path(out_dir, fn))
  #     }
  #   }
  #   logger::log_info("β-contrast tables written to {out_dir}")
  # }

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

  # write_enrichment_tables(enrich, cfg$plot_dir)
  # write_beta_contrast_tables(beta_contrast_tables, cfg$plot_dir)
  write_beta_tables(beta_tables, cfg$plot_dir)

  # ---- 9) Keep the originals around too (optional) ----
  # manhattan <- manhattan_BH_all
  # volcano   <- volcano_default

  # ---- 10) Summary counts (unchanged placeholder) ----
  summary_tbl <- emit_summary_counts(metrics)

  output <- list(
    MR_df = MR_df,
    results_df    = results_df,
    summary_plots = summary_plots,

    enrich        = enrich,
    beta          = beta_tables
#    beta_contrast = beta_contrast_tables,
#    exposure_units = exposure_units
  )

  saveRDS(output, file = file.path(cfg$plot_dir, "results.rds"))
  output
}
