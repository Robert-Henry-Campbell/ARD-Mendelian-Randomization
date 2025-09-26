#R/run_phenome.MR
#' Run phenome-wide MR across phenotypes (orchestrator)
#'
#' Builds the outcome plan, maps exposure SNPs, fetches outcome SNP rows,
#' runs MR + sensitivity/QC, and returns results + plots. Phenotypes flagged by
#' `ARD_selected` denote age-related diseases.
#'
#' @param exposure Character (mandatory) label for the exposure; used for
#'   naming output directories and plot titles.
#'   rsid, beta, se, effect_allele, other_allele, eaf, etc.).
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
#'   they are already cached.
#'
#' @return A list: MR_df, results_df, manhattan (ggplot), volcano (ggplot)
#' @export
run_phenome_mr_plotting_only_man_and_volc_only <- function(
    run_output,
    exposure,
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
  if (missing(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- as.character(exposure_units)
  if (!length(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- exposure_units[1]
  if (is.na(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- trimws(exposure_units)
  if (!nzchar(exposure_units)) stop("`exposure_units` is mandatory.")
  if (missing(ancestry) || !nzchar(ancestry)) stop("`ancestry` is mandatory.")

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

#####################################################################################

  #ALL ANALYSIS REMOVED FROM HERE


########################################################################################3


### LOADING UP MY STUFF

  results_df <- run_output$results_df
  MR_df <- run_output$MR_df


#################################################################################
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
      yfloat_base <- 1.8
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
        height <- yfloat_base + (yfloat_coef * 1.3) * n_rows
      }
      if (grepl("^volcano", subpath)) {
        width <- 7.2; height <- 5.0
      }
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

  manhattan_recolor_BH_all <- list(
    with_dot_names = manhattan_plot_recolor(
      results_df,
      Multiple_testing_correction = "BH",
      exposure = exposure,
      dot_names = TRUE
    ),
    without_dot_names = manhattan_plot_recolor(
      results_df,
      Multiple_testing_correction = "BH",
      exposure = exposure,
      dot_names = FALSE
    )
  )
  manhattan_recolor_BH_ARD <- list(
    with_dot_names = manhattan_plot_recolor(
      results_ard_only,
      Multiple_testing_correction = "BH",
      exposure = exposure,
      dot_names = TRUE
    ),
    without_dot_names = manhattan_plot_recolor(
      results_ard_only,
      Multiple_testing_correction = "BH",
      exposure = exposure,
      dot_names = FALSE
    )
  )
  manhattan_recolor_Bonf_all <- list(
    with_dot_names = manhattan_plot_recolor(
      results_df,
      Multiple_testing_correction = "bonferroni",
      exposure = exposure,
      dot_names = TRUE
    ),
    without_dot_names = manhattan_plot_recolor(
      results_df,
      Multiple_testing_correction = "bonferroni",
      exposure = exposure,
      dot_names = FALSE
    )
  )
  manhattan_recolor_Bonf_ARD <- list(
    with_dot_names = manhattan_plot_recolor(
      results_ard_only,
      Multiple_testing_correction = "bonferroni",
      exposure = exposure,
      dot_names = TRUE
    ),
    without_dot_names = manhattan_plot_recolor(
      results_ard_only,
      Multiple_testing_correction = "bonferroni",
      exposure = exposure,
      dot_names = FALSE
    )
  )

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

  volcano_recolor_with_dot_names <- volcano_plot_recolor(
    results_df,
    Multiple_testing_correction = cfg$mtc,
    exposure = exposure,
    verbose = cfg$verbose,
    dot_names = TRUE
  )
  volcano_recolor_without_dot_names <- volcano_plot_recolor(
    results_df,
    Multiple_testing_correction = cfg$mtc,
    exposure = exposure,
    verbose = cfg$verbose,
    dot_names = FALSE
  )

  volcano_recolor_plot_list <- list(
    BH = list(all = volcano_recolor_BH_all, ARD_only = volcano_recolor_BH_ARD),
    bonferroni = list(all = volcano_recolor_Bonf_all, ARD_only = volcano_recolor_Bonf_ARD)
  )
  volcano_recolor_plot_list[[cfg$mtc]][["with_dot_names"]] <- volcano_recolor_with_dot_names
  volcano_recolor_plot_list[[cfg$mtc]][["without_dot_names"]] <- volcano_recolor_without_dot_names

  summary_plots_base <- list(
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
    volcano_recolor = volcano_recolor_plot_list
  )


    summary_plots <- summary_plots_base
    plots_to_save <- list(
      manhattan_recolor = summary_plots$manhattan_recolor,
      volcano_recolor = summary_plots$volcano_recolor
    )
    save_plot_hierarchy(plots_to_save, cfg$plot_dir)

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



  # ---- 7) Assemble hierarchical summary_plots list ----
  summary_plots <- summary_plots_base



  # ---- 8) Save plots mirroring the list structure under cfg$plot_dir ----
  # save_plot_hierarchy(summary_plots, cfg$plot_dir)
  plots_to_save <- summary_plots
  plots_to_save$beta <- NULL
  plots_to_save$manhattan <- NULL
  plots_to_save$volcano <- NULL
  save_plot_hierarchy(plots_to_save, cfg$plot_dir)

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
    summary_plots = summary_plots

    #enrich        = enrich,
    #beta          = beta_tables
    #    beta_contrast = beta_contrast_tables,
    #    exposure_units = exposure_units
  )

  saveRDS(output, file = file.path(cfg$plot_dir, "results.rds"))
  output
}
