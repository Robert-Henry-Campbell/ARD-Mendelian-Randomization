#R/run_phenome.MR
#' Run phenome-wide MR across phenotypes (orchestrator)
#'
#' Builds the outcome plan, maps exposure SNPs, fetches outcome SNP rows,
#' runs MR + sensitivity/QC, and returns results + plots. Phenotypes flagged by
#' `ARD_selected` denote age-related diseases.
#'
#' @param exposure Character (mandatory) label for the exposure; used for
#'   naming output directories and plot titles.
#' @param exposure_snps Optional data frame of pre-filtered exposure instruments
#'   (TwoSampleMR-like). Required columns include `id.exposure`, `SNP`,
#'   `effect_allele.exposure`, `other_allele.exposure`, `beta.exposure`,
#'   `se.exposure`, `pval.exposure`. May be combined with `exposure_sumstats`
#'   (snps used as IVs, sumstats used for coloc).
#' @param exposure_sumstats Optional path to a tabix-indexed GWAS-VCF (MRC IEU
#'   spec) for the exposure. Used as the regional summary-statistics source
#'   for colocalization. Cannot be combined with `exposure_id`.
#' @param exposure_id Optional OpenGWAS study id (e.g. `"ieu-b-38"`). When
#'   supplied alone, IVs are extracted via `TwoSampleMR::extract_instruments`
#'   and the same id is used for coloc.
#' @param exposure_units Character (mandatory) description of the exposure
#'   units for β-scale plots.
#' @param ancestry Character (mandatory), e.g. "EUR".
#' @param sex One of "both","male","female". If not "both", you likely use Neale.
#' @param sensitivity_enabled Character vector of checks (default = all 9 incl. coloc).
#' @param sensitivity_pass_min Integer threshold to pass QC (default 6).
#' @param genome_build Genome build of `exposure_sumstats` (default `"GRCh37"`;
#'   anything else triggers a hard error).
#' @param acknowledge_no_coloc Set to `TRUE` to acknowledge that supplying
#'   `exposure_snps` alone disables coloc. Required in that case.
#' @param clump_opts Named list of preprocessing parameters used by the
#'   unified exposure-SNP preprocessing pipeline (see
#'   [preprocess_exposure_snps()]). Applied to all input modes
#'   (`snps_only`, `snps_vcf`, `vcf_only`, `ieugwasr`, `snps_id`) for
#'   parity. The pipeline runs in this order: p-value -> rsid validate
#'   -> drop indels -> MAF -> INFO -> palindromic -> LD clumping ->
#'   F-stat. Data-quality filters precede clumping so a clump-elected
#'   lead SNP that fails one of them does not cost the entire locus.
#'   Recognised fields:
#'   \describe{
#'     \item{`p_backoff`}{Numeric vector of p-value thresholds tried
#'       strictest-first. Default `c(5e-8)` (single rung); pass a
#'       multi-element vector to enable a backoff ladder
#'       (e.g. `c(5e-8, 5e-7, 5e-6)`). Replaces the deprecated
#'       `p_threshold` scalar (silently translated with a deprecation
#'       warning).}
#'     \item{`r2`, `kb`}{LD clumping params (defaults `0.001`, `10000`).
#'       Clumping is always local against the ancestry-aware 1kg.v3
#'       panel.}
#'     \item{`f_threshold`}{F-statistic minimum (default `10`). Runs
#'       AFTER clumping (per-SNP, no locus loss).}
#'     \item{`maf_min`}{Minimum minor-allele frequency. Default `0.01`
#'       (1%, standard MR practice); applied as `pmin(eaf, 1-eaf) >=
#'       maf_min`. Pass `NULL` to disable. Rows with `NA` EAF are
#'       kept; the step is automatically skipped when the
#'       `eaf.exposure` column is missing.}
#'     \item{`info_min`}{Minimum INFO/imputation-quality score (default
#'       `NULL` = off); auto-detects an `SI`/`Rsq`/`R2`/`INFO`/`info`
#'       column (SI/Rsq/R2 take priority because in real GWAS-VCF data
#'       the `INFO` column is typically the raw `;`-separated VCF INFO
#'       blob, not a numeric score). Rows with `NA` INFO are kept. The
#'       same threshold also applies to coloc-region fetching.}
#'     \item{`drop_indels`}{Drop variants where either allele is non-SNV
#'       (`nchar > 1` or `-`/`D`/`I` codes). Default `TRUE`.}
#'     \item{`drop_palindromic`}{Drop AMBIGUOUS palindromic SNPs
#'       (A/T or C/G pairs AND `eaf in [0.42, 0.58]`).
#'       Strand-resolvable palindromes (extreme EAF) are kept for
#'       downstream `harmonise_data` to align from EAF. The
#'       preprocessor does NOT add a `palindromic` column;
#'       TwoSampleMR's harmonise computes its own. Default `FALSE`.}
#'     \item{`already_clumped`, `already_p_filtered`}{Caller hint that LD
#'       clumping or p-value filtering have already been done upstream;
#'       skips the corresponding local step. Defaults `FALSE`.}
#'     \item{`preprocess`}{Master opt-out: when `FALSE`, no preprocessing
#'       is applied (the input is passed through as-is). Default `TRUE`.}
#'   }
#'   The merged-and-applied list is recorded in `run_manifest.json` as
#'   `clump_opts_resolved`; per-step `n_in`/`n_out` accounting is in
#'   `preprocessing_steps`.
#' @param coloc_window_kb Half-window (kb) around each IV used to build coloc
#'   regions; overlapping windows are merged (default 500).
#' @param coloc_priors Coloc priors `list(p1, p2, p12)` (defaults match coloc).
#' @param coloc_skip_mhc Skip loci overlapping chr6:25-35Mb (default `TRUE`).
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
#'   they are already cached. Defaults to `FALSE` — caches are now keyed by
#'   IV-set hash + run-input hash so cross-run collisions are impossible.
#'   When `TRUE`, only entries matching the current run's keys are deleted
#'   (other exposures' caches and other runs' outputs are preserved).
#' @param n_pheno_limit Optional positive integer. If supplied, truncate the
#'   outcome phenotype list to the first N rows immediately after
#'   [Outcome_setup()]. Intended for fast end-to-end integration tests; the
#'   resulting multiple-testing thresholds, p-values, and enrichment p-values
#'   are NOT scientifically valid. Default `NULL` (no truncation).
#'
#' @return A list: `MR_df`, `results_df`, `manhattan` (ggplot),
#'   `volcano` (ggplot). When preprocessing strips every IV, the
#'   function emits a `logger::log_warn` (not an error), still writes
#'   the run manifest so the user can inspect per-step counts, and
#'   returns
#'   `list(MR_df = tibble(), results_df = tibble(), manhattan = NULL,
#'   volcano = NULL, run_dir = <plot_output_dir>, status = "no_ivs")`.
#' @export
run_phenome_mr <- function(
    exposure,
    exposure_units,
    ancestry,
    exposure_snps = NULL,
    exposure_sumstats = NULL,
    exposure_id = NULL,
    sex = c("both","male","female"),
    sensitivity_enabled = c(
      "egger_intercept","egger_slope_agreement",
      "weighted_median","weighted_mode",
      "steiger_direction","leave_one_out",
      "ivw_Q","ivw_I2","coloc"
    ),
    sensitivity_pass_min = 6,
    Multiple_testing_correction = c("BH","bonferroni"),
    scatterplot = TRUE,
    snpforestplot = TRUE,
    leaveoneoutplot = TRUE,
    genome_build = "GRCh37",
    acknowledge_no_coloc = FALSE,
    clump_opts = list(),
    coloc_window_kb = 500L,
    coloc_priors = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    coloc_skip_mhc = TRUE,
    cache_dir = ardmr_cache_dir(),
    logfile = NULL,
    verbose = TRUE,
    confirm = 'ask',
    force_refresh = FALSE,
    n_pheno_limit = NULL

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
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- resolve coloc source / IV source ----
  resolved <- .resolve_coloc_source(
    exposure_snps      = exposure_snps,
    exposure_sumstats  = exposure_sumstats,
    exposure_id        = exposure_id,
    ancestry           = ancestry,
    genome_build       = genome_build,
    acknowledge_no_coloc = acknowledge_no_coloc,
    clump_opts         = clump_opts,
    cache_dir          = cache_dir,
    confirm            = confirm,
    verbose            = verbose
  )
  exposure_snps       <- resolved$exposure_snps
  clump_opts_resolved <- resolved$clump_opts_resolved
  preprocessing_steps <- resolved$preprocessing_steps
  assert_exposure(exposure_snps)
  # 0-IV is no longer fatal (Job 4): we still build the run dir + write
  # the manifest so the user can inspect what happened, then warn and
  # return an empty result list with status = "no_ivs". The actual
  # warn+return happens after the manifest is written, below.
  no_ivs <- (NROW(exposure_snps) == 0L)
  if (resolved$mode == "snps_only") {
    sensitivity_enabled <- setdiff(sensitivity_enabled, "coloc")
    if (verbose) logger::log_info("coloc-source: mode='snps_only'; coloc disabled for this run")
  } else if (verbose) {
    logger::log_info("coloc-source: mode='{resolved$mode}'")
  }

  # Always derive standard subfolders inside cache_dir
  .slug <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "NA"
    x <- gsub("[^[:alnum:]]+", "-", x)
    x <- gsub("-+", "-", x)
    x <- gsub("^-|-$", "", x)
    tolower(x)
  }

  run_hash <- .run_input_hash(
    exposure_snps               = exposure_snps,
    exposure_id                 = exposure_id,
    exposure_sumstats           = exposure_sumstats,
    sensitivity_enabled         = sensitivity_enabled,
    sensitivity_pass_min        = sensitivity_pass_min,
    clump_opts                  = clump_opts,
    coloc_window_kb             = coloc_window_kb,
    coloc_priors                = coloc_priors,
    coloc_skip_mhc              = coloc_skip_mhc,
    multiple_testing_correction = Multiple_testing_correction
  )
  plot_output_dir <- ardmr_run_dir(cache_dir, .slug(exposure), sex, ancestry, run_hash)
  Neale_GWAS_dir  <- file.path(cache_dir, "neale_sumstats")
  purged_plot_dir <- FALSE
  if (isTRUE(force_refresh) && dir.exists(plot_output_dir)) {
    unlink(plot_output_dir, recursive = TRUE, force = TRUE)
    purged_plot_dir <- TRUE
  }
  dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(Neale_GWAS_dir,  recursive = TRUE, showWarnings = FALSE)

  # ---- write run_manifest.json ----
  iv_hash <- .iv_set_hash(exposure_snps, n = 5L)
  run_manifest <- list(
    run_hash             = run_hash,
    started_at           = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    exposure             = exposure,
    exposure_units       = exposure_units,
    ancestry             = ancestry,
    sex                  = sex,
    n_ivs                = NROW(exposure_snps),
    iv_hash              = iv_hash,
    coloc_source_mode    = resolved$mode,
    exposure_id          = exposure_id,
    exposure_sumstats    = if (is.null(exposure_sumstats)) NULL else basename(as.character(exposure_sumstats)),
    sensitivity_enabled  = sensitivity_enabled,
    sensitivity_pass_min = sensitivity_pass_min,
    clump_opts           = clump_opts,
    clump_opts_resolved  = clump_opts_resolved,
    preprocessing_steps  = preprocessing_steps,
    coloc                = list(window_kb = coloc_window_kb, priors = coloc_priors,
                                skip_mhc = coloc_skip_mhc),
    multiple_testing_correction = Multiple_testing_correction,
    force_refresh        = force_refresh,
    genome_build         = genome_build
  )
  tryCatch(
    jsonlite::write_json(run_manifest,
                         file.path(plot_output_dir, "run_manifest.json"),
                         pretty = TRUE, auto_unbox = TRUE, null = "null"),
    error = function(e) logger::log_warn("Could not write run_manifest.json: {conditionMessage(e)}")
  )

  # ---- 0-IV early return (warn, do not stop) ----
  if (isTRUE(no_ivs)) {
    logger::log_warn(.format_zero_iv_error(preprocessing_steps))
    return(list(
      MR_df      = tibble::tibble(),
      results_df = tibble::tibble(),
      manhattan  = NULL,
      volcano    = NULL,
      run_dir    = plot_output_dir,
      status     = "no_ivs"
    ))
  }

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

  metrics <- list()

  logger::log_info("1) Outcome setup…")
  MR_df <- Outcome_setup(sex = cfg$sex, ancestry = cfg$ancestry)

  if (!is.null(n_pheno_limit)) {
    n_pheno_limit <- suppressWarnings(as.integer(n_pheno_limit))
    if (length(n_pheno_limit) != 1L || is.na(n_pheno_limit) || n_pheno_limit < 1L) {
      stop("`n_pheno_limit` must be a single positive integer or NULL.")
    }
    n_before <- nrow(MR_df)
    MR_df <- head(MR_df, n_pheno_limit)
    logger::log_warn(
      "n_pheno_limit set: truncated outcome list from {n_before} to {nrow(MR_df)} phenotypes. Multiple-testing thresholds and enrichment p-values from this run are NOT scientifically valid."
    )
  }

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

  try(utils::write.csv(exposure_snps2, file.path(cfg$plot_dir, "exposure_snps.csv"), row.names = FALSE), silent = TRUE)

  logger::log_info("2.5) Compute LD matrix for exposure SNPs…")
  ld_csv_path <- file.path(cfg$plot_dir, "exposure_snps_ld_matrix.csv")
  ld_try <- tryCatch(
    compute_ld_matrix(
      exposure_snps = exposure_snps2,
      ancestry      = cfg$ancestry,
      cache_dir     = cfg$cache_dir,
      out_csv       = ld_csv_path,
      confirm       = cfg$confirm,
      force_refresh = cfg$force_refresh,
      verbose       = cfg$verbose
    ),
    error = function(e) e
  )
  if (inherits(ld_try, "error")) {
    logger::log_warn("LD matrix computation failed (non-fatal): {conditionMessage(ld_try)}")
  } else {
    metrics$ld_snps_used <- length(ld_try$snps_used)
  }

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
  coloc_opts <- list()
  if ("coloc" %in% cfg$checks_enabled) {
    if (cfg$catalog == "panukb") {
      fetcher_factory <- function(rec) {
        function(chr, start, end) {
          coloc_fetch_outcome_panukb_region(
            rec = rec, ancestry = cfg$ancestry,
            chr = chr, start = start, end = end,
            cache_dir = cfg$cache_dir, verbose = cfg$verbose
          )
        }
      }
      meta_fn <- function(rec) coloc_panukb_outcome_metadata(rec, ancestry = cfg$ancestry)
    } else if (cfg$catalog == "neale") {
      fetcher_factory <- function(rec) {
        function(chr, start, end) {
          coloc_fetch_outcome_neale_region(
            rec = rec, sex = cfg$sex,
            chr = chr, start = start, end = end,
            neale_dir = cfg$neale_dir, cache_dir = cfg$cache_dir,
            verbose = cfg$verbose,
            force_refresh = isTRUE(cfg$force_refresh)
          )
        }
      }
      meta_fn <- function(rec) coloc_neale_outcome_metadata(rec, sex = cfg$sex)
    } else {
      logger::log_warn("coloc: catalog '{cfg$catalog}' not yet supported; coloc will be skipped per outcome")
      fetcher_factory <- function(rec) NULL
      meta_fn <- function(rec) NULL
    }
    coloc_opts <- list(
      exposure_region_fetcher = resolved$coloc_fetcher,
      exposure_metadata       = resolved$coloc_metadata,
      ancestry                = cfg$ancestry,
      window_kb               = coloc_window_kb,
      priors                  = coloc_priors,
      skip_mhc                = coloc_skip_mhc,
      outcome_fetcher_factory = fetcher_factory,
      outcome_metadata_fn     = meta_fn
    )
    if (is.null(resolved$coloc_fetcher) || is.null(resolved$coloc_metadata)) {
      logger::log_warn("coloc enabled but no exposure source available; coloc will be skipped at runtime")
    }
  }
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
    coloc_opts = coloc_opts,
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

  # ---- 6) Signed enrichment analyses ----
  logger::log_info("6) Enrichment analyses…")
  enrich_try <- tryCatch(
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
    error = function(e) e
  )

  if (inherits(enrich_try, "error")) {
    err_msg <- conditionMessage(enrich_try)
    warn_msg <- sprintf("Enrichment failed: %s", err_msg)
    warning(warn_msg, call. = FALSE)
    logger::log_error("{warn_msg}")

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

    output <- list(
      MR_df = MR_df,
      results_df    = results_df,
      summary_plots = summary_plots,
      enrich        = NULL,
      beta          = NULL
    )

    saveRDS(output, file = file.path(cfg$plot_dir, "results.rds"))
    logger::log_info("Early exit: wrote results.rds to {cfg$plot_dir} after enrichment failure.")
    return(output)
  }

  enrich <- enrich_try

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
  summary_plots <- summary_plots_base
  summary_plots$enrichment <- list(
      global = list(
        violin_vertical = enrichment_global_violin_vertical,
        violin_forest = enrichment_global_violin_forest
      ),
      cause_level_1 = enrichment_cause_violin[["cause_level_1"]],
      cause_level_2 = enrichment_cause_violin[["cause_level_2"]],
      cause_level_3 = enrichment_cause_violin[["cause_level_3"]]
    )
  summary_plots$beta <- beta_plots

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
  # save_plot_hierarchy(summary_plots, cfg$plot_dir)
  plots_to_save <- summary_plots
  plots_to_save$beta <- NULL
  plots_to_save$manhattan <- NULL
  plots_to_save$volcano <- NULL
  save_plot_hierarchy(plots_to_save, cfg$plot_dir)

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
  # write_beta_tables(beta_tables, cfg$plot_dir)

  # ---- 9) Keep the originals around too (optional) ----
  # manhattan <- manhattan_BH_all
  # volcano   <- volcano_default

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
    beta          = beta_tables
#    beta_contrast = beta_contrast_tables,
#    exposure_units = exposure_units
  )

  saveRDS(output, file = file.path(cfg$plot_dir, "results.rds"))
  output
}
