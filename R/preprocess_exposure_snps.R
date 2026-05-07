# R/preprocess_exposure_snps.R
# Unified exposure-SNP preprocessing across all input modes (snps_only,
# snps_vcf, vcf_only, ieugwasr, snps_id). Called by .resolve_coloc_source()
# and (transitively) by clump_sumstats_to_ivs().

#' Apply the unified exposure-SNP preprocessing pipeline.
#'
#' Steps, in order (locked):
#'   1. p-value filter        (default single rung c(5e-8); ladder opt-in;
#'                             skipped if `already_p_filtered = TRUE`)
#'   2. rsid validate + dedupe
#'   3. drop indels           (default ON)
#'   4. MAF filter            (default ON at `maf_min = 0.01`; pass
#'                             `maf_min = NULL` to disable; auto-skipped
#'                             when `eaf.exposure` is absent;
#'                             `pmin(eaf, 1-eaf)`)
#'   5. INFO filter           (default OFF; if `info_min` set AND an INFO
#'                             column is present; auto-detects
#'                             `SI`/`Rsq`/`R2`/`INFO`/`info`)
#'   6. palindromic drop      (default OFF; EAF-aware -- only ambiguous
#'                             palindromes (eaf in [0.42, 0.58]) dropped)
#'   7. LD clumping           (always local against the ancestry-aware
#'                             1kg.v3 panel; skipped if
#'                             `already_clumped = TRUE`)
#'   8. F-statistic           (>= `f_threshold`; runs after clumping
#'                             because it's per-SNP and won't cause
#'                             locus loss)
#'
#' Data-quality filters (3-6) run BEFORE clumping so the clump step
#' picks lead SNPs from a clean universe. F-stat runs AFTER clumping
#' because it's a per-SNP property (no locus loss).
#'
#' Returns a list with the filtered tibble, per-step accounting, and the
#' merged `clump_opts`. A 0-row return is normal when filters strip
#' everything; the function never errors on empty output (caller decides).
#'
#' @param snps_tbl A TwoSampleMR-shaped tibble (must contain `SNP`).
#' @param clump_opts Named list overriding the defaults documented at
#'   `?run_phenome_mr` (`p_backoff`, `r2`, `kb`, `f_threshold`, `maf_min`,
#'   `info_min`, `drop_indels`, `drop_palindromic`, `already_clumped`,
#'   `already_p_filtered`, `preprocess`). For backward compat,
#'   `p_threshold` is accepted and translated to `p_backoff =
#'   c(p_threshold)` with a deprecation warning.
#' @param ancestry Ancestry code, e.g. "EUR".
#' @param cache_dir Cache root; LD reference is read from
#'   `<cache_dir>/ld_reference`.
#' @param confirm One of `"ask"`, `"yes"`, `"no"` -- passed to
#'   [ld_reference_checker()] when the LD panel must be installed for
#'   clumping. `"no"` aborts cleanly if the panel is missing (no silent
#'   download).
#' @param verbose Logical; per-step counts logged when TRUE.
#' @return A list:
#'   * `snps`: filtered tibble (TwoSampleMR shape; may be 0-row).
#'   * `steps`: list of per-step accounting records (each has `step`,
#'     `n_in`, `n_out`, plus step-specific extras like `value`,
#'     `skipped_reason`).
#'   * `resolved_opts`: the merged `defaults <- modifyList(...) -> co`.
#' @keywords internal
#' Canonical defaults for the preprocess pipeline.
#'
#' Single source of truth shared by [preprocess_exposure_snps()] and
#' [.resolve_coloc_source()] so the two cannot drift.
#' @keywords internal
.preprocess_defaults <- function() {
  list(
    p_backoff           = c(5e-8),
    r2                  = 0.001,
    kb                  = 10000L,
    f_threshold         = 10,
    # MAF filter is ON by default at 1% (standard MR practice). Step
    # is skipped automatically when `eaf.exposure` is missing, and
    # rows with NA EAF are kept. Pass `maf_min = NULL` to disable.
    maf_min             = 0.01,
    info_min            = NULL,
    drop_indels         = TRUE,
    drop_palindromic    = FALSE,
    already_clumped     = FALSE,
    already_p_filtered  = FALSE,
    preprocess          = TRUE
  )
}

#' Format a `preprocessing_steps` list into a human-readable error message
#' for the 0-IV case.
#'
#' Lists every recorded step with its `n_in -> n_out` and a marker on the
#' first step that took the count from positive to zero. Falls back to a
#' simpler message when no steps were recorded (e.g. preprocess = FALSE
#' bypass).
#' @keywords internal
.format_zero_iv_error <- function(steps) {
  if (!length(steps)) {
    return(paste0(
      "Preprocessing returned 0 instruments and no preprocessing steps were ",
      "recorded (likely clump_opts$preprocess = FALSE or a bypassed mode). ",
      "Check exposure input or re-enable preprocessing."
    ))
  }
  killed <- NA_character_
  fmt_step <- function(s) {
    desc <- as.character(s$step %||% "?")
    val_str <- ""
    if (!is.null(s$value) && length(s$value) == 1L && !is.na(s$value)) {
      val_str <- sprintf(" (%g)", s$value)
    }
    skip_str <- if (!is.null(s$skipped_reason)) sprintf(" [skipped: %s]", s$skipped_reason) else ""
    n_in  <- if (is.null(s$n_in)) NA_integer_ else as.integer(s$n_in)
    n_out <- if (is.null(s$n_out)) NA_integer_ else as.integer(s$n_out)
    arrow <- ""
    if (is.na(killed) && !is.na(n_in) && !is.na(n_out) && n_in > 0L && n_out == 0L) {
      killed <<- desc
      arrow <- "  <- everything dropped here"
    }
    sprintf("  %s%s: %d -> %d%s%s", desc, val_str, n_in, n_out, skip_str, arrow)
  }
  lines <- c(
    "Preprocessing returned 0 instruments. Per-step counts:",
    vapply(steps, fmt_step, character(1))
  )
  if (!is.na(killed)) {
    lines <- c(lines, sprintf(
      "Likely cause: '%s' step. Adjust clump_opts (e.g. drop_indels = FALSE, looser p_backoff/maf_min/info_min) or check input quality.",
      killed
    ))
  }
  paste(lines, collapse = "\n")
}

#' Translate deprecated `p_threshold` -> `p_backoff`, then merge with defaults.
#'
#' The `p_threshold` deprecation warning fires here; downstream callers can
#' assume `co$p_backoff` is the only pval-source field.
#' @keywords internal
.resolve_clump_opts <- function(clump_opts = list()) {
  if (!is.null(clump_opts$p_threshold)) {
    if (!is.null(clump_opts$p_backoff)) {
      logger::log_warn(
        "clump_opts: both 'p_threshold' (deprecated) and 'p_backoff' supplied; using 'p_backoff' and ignoring 'p_threshold'."
      )
    } else {
      pt <- as.numeric(clump_opts$p_threshold)
      logger::log_warn(
        "clump_opts$p_threshold is deprecated; translating to p_backoff = c({pt}). Use p_backoff directly."
      )
      clump_opts$p_backoff <- pt
    }
    clump_opts$p_threshold <- NULL
  }
  # `keep.null = TRUE` so a user explicitly passing `maf_min = NULL`
  # (or `info_min = NULL`) actually disables the corresponding step
  # rather than silently inheriting the default.
  utils::modifyList(.preprocess_defaults(), clump_opts, keep.null = TRUE)
}

preprocess_exposure_snps <- function(snps_tbl,
                                     clump_opts = list(),
                                     ancestry,
                                     cache_dir = ardmr_cache_dir(),
                                     confirm = "ask",
                                     verbose = TRUE) {
  co <- .resolve_clump_opts(clump_opts)

  if (!is.data.frame(snps_tbl)) {
    stop("preprocess_exposure_snps(): snps_tbl must be a data.frame.", call. = FALSE)
  }
  if (!"SNP" %in% names(snps_tbl)) {
    stop("preprocess_exposure_snps(): snps_tbl must have a 'SNP' column.", call. = FALSE)
  }
  if (missing(ancestry) || is.null(ancestry) || !nzchar(as.character(ancestry)[1])) {
    stop("preprocess_exposure_snps(): `ancestry` is required.", call. = FALSE)
  }

  steps_acc <- list()
  add_step <- function(step, n_in, n_out, ...) {
    rec <- c(list(step = step, n_in = as.integer(n_in), n_out = as.integer(n_out)), list(...))
    steps_acc[[length(steps_acc) + 1L]] <<- rec
  }

  # ---- Step 1: p-value filter (single rung default; ladder opt-in) ----
  n_in <- nrow(snps_tbl)
  if (isTRUE(co$already_p_filtered)) {
    if (verbose) logger::log_info("preprocess: step 1 (p-value) skipped (already_p_filtered)")
    add_step("p_threshold", n_in, n_in, value = NA_real_, skipped_reason = "already_p_filtered")
  } else if (!"pval.exposure" %in% names(snps_tbl)) {
    logger::log_warn("preprocess: step 1 (p-value) skipped: 'pval.exposure' column missing")
    add_step("p_threshold", n_in, n_in, value = NA_real_, skipped_reason = "no pval.exposure column")
  } else {
    pvals <- suppressWarnings(as.numeric(snps_tbl$pval.exposure))
    rungs <- as.numeric(co$p_backoff)
    chosen <- NA_real_
    keep <- rep(FALSE, n_in)
    for (rung in rungs) {
      kk <- !is.na(pvals) & pvals < rung
      if (sum(kk) > 0L) {
        keep <- kk
        chosen <- rung
        break
      }
    }
    if (is.na(chosen)) {
      chosen <- rungs[length(rungs)]
      keep <- rep(FALSE, n_in)
    }
    snps_tbl <- snps_tbl[keep, , drop = FALSE]
    if (verbose) logger::log_info("preprocess: step 1 (p-value) p<{chosen}: {n_in} -> {nrow(snps_tbl)}")
    add_step("p_threshold", n_in, nrow(snps_tbl), value = chosen)
  }

  # ---- Step 2: rsid validate + dedupe (must run before clumping) ----
  n_in <- nrow(snps_tbl)
  if (n_in == 0L) {
    add_step("rsid_validate", 0L, 0L, n_invalid = 0L, n_duplicates = 0L)
  } else {
    valid <- is_rsid(as.character(snps_tbl$SNP))
    n_invalid <- sum(!valid)
    if (n_invalid > 0L) {
      logger::log_warn("preprocess: step 2 dropped {n_invalid} rows with non-rsid SNP IDs")
    }
    snps_tbl <- snps_tbl[valid, , drop = FALSE]

    dups <- duplicated(snps_tbl$SNP)
    n_dup <- sum(dups)
    if (n_dup > 0L) {
      logger::log_warn(
        "preprocess: step 2 collapsed {n_dup} duplicate rsid rows (keeping first occurrence)"
      )
    }
    snps_tbl <- snps_tbl[!dups, , drop = FALSE]
    if (verbose) {
      logger::log_info("preprocess: step 2 (rsid validate+dedupe) {n_in} -> {nrow(snps_tbl)}")
    }
    add_step("rsid_validate", n_in, nrow(snps_tbl),
             n_invalid = as.integer(n_invalid), n_duplicates = as.integer(n_dup))
  }

  # ---- Step 3: drop_indels (data quality, before clumping) ----
  n_in <- nrow(snps_tbl)
  if (!isTRUE(co$drop_indels)) {
    add_step("drop_indels", n_in, n_in, skipped_reason = "drop_indels = FALSE")
  } else if (n_in == 0L) {
    add_step("drop_indels", 0L, 0L)
  } else if (!all(c("effect_allele.exposure", "other_allele.exposure") %in% names(snps_tbl))) {
    logger::log_warn("preprocess: step 3 (indel) skipped: allele columns missing")
    add_step("drop_indels", n_in, n_in, skipped_reason = "allele columns missing")
  } else {
    ea <- toupper(as.character(snps_tbl$effect_allele.exposure))
    oa <- toupper(as.character(snps_tbl$other_allele.exposure))
    indel_codes <- c("-", "D", "I")
    is_indel <- nchar(ea) > 1 | nchar(oa) > 1 | ea %in% indel_codes | oa %in% indel_codes
    is_indel[is.na(is_indel)] <- FALSE
    snps_tbl <- snps_tbl[!is_indel, , drop = FALSE]
    if (verbose) logger::log_info("preprocess: step 3 (drop_indels) {n_in} -> {nrow(snps_tbl)}")
    add_step("drop_indels", n_in, nrow(snps_tbl))
  }

  # ---- Step 4: MAF (data quality, before clumping) ----
  n_in <- nrow(snps_tbl)
  if (is.null(co$maf_min)) {
    add_step("maf", n_in, n_in, value = NA_real_, skipped_reason = "maf_min not set")
  } else if (n_in == 0L) {
    add_step("maf", 0L, 0L, value = co$maf_min)
  } else if (!"eaf.exposure" %in% names(snps_tbl)) {
    logger::log_warn("preprocess: step 4 (MAF) skipped: 'eaf.exposure' column missing")
    add_step("maf", n_in, n_in, value = co$maf_min,
             skipped_reason = "no eaf.exposure column")
  } else {
    eaf <- suppressWarnings(as.numeric(snps_tbl$eaf.exposure))
    maf <- pmin(eaf, 1 - eaf)
    keep <- is.na(maf) | maf >= co$maf_min
    snps_tbl <- snps_tbl[keep, , drop = FALSE]
    if (verbose) {
      logger::log_info("preprocess: step 4 (MAF>={co$maf_min}) {n_in} -> {nrow(snps_tbl)}")
    }
    add_step("maf", n_in, nrow(snps_tbl), value = co$maf_min)
  }

  # ---- Step 5: INFO (data quality, before clumping) ----
  # Priority order: SI/Rsq/R2 (canonical numeric imputation-quality
  # fields) BEFORE INFO/info (which in real GWAS-VCF data is usually
  # the raw `;`-separated VCF INFO blob, not a numeric quality score).
  n_in <- nrow(snps_tbl)
  info_col_candidates <- c("SI", "Rsq", "R2", "INFO", "info")
  if (is.null(co$info_min)) {
    add_step("info", n_in, n_in, value = NA_real_, skipped_reason = "info_min not set")
  } else if (n_in == 0L) {
    add_step("info", 0L, 0L, value = co$info_min)
  } else {
    matched <- info_col_candidates[info_col_candidates %in% names(snps_tbl)]
    if (!length(matched)) {
      logger::log_warn(
        "preprocess: step 5 (INFO) skipped: no INFO column found (looked for {paste(info_col_candidates, collapse=', ')})"
      )
      add_step("info", n_in, n_in, value = co$info_min, skipped_reason = "no INFO column")
    } else {
      info_col <- matched[1]
      info_vals <- suppressWarnings(as.numeric(snps_tbl[[info_col]]))
      keep <- is.na(info_vals) | info_vals >= co$info_min
      snps_tbl <- snps_tbl[keep, , drop = FALSE]
      if (verbose) {
        logger::log_info(
          "preprocess: step 5 (INFO>={co$info_min}, col={info_col}) {n_in} -> {nrow(snps_tbl)}"
        )
      }
      add_step("info", n_in, nrow(snps_tbl), value = co$info_min, column = info_col)
    }
  }

  # ---- Step 6: palindromic drop (EAF-aware, before clumping) ----
  # Only AMBIGUOUS palindromes (A/T or C/G AND eaf in [0.42, 0.58])
  # are dropped. Strand-resolvable palindromes are kept for downstream
  # `harmonise_data` to align from EAF. No `palindromic` column is
  # added to the output -- TwoSampleMR's harmonise computes its own.
  n_in <- nrow(snps_tbl)
  if (!isTRUE(co$drop_palindromic)) {
    add_step("palindromic", n_in, n_in, dropped = FALSE,
             skipped_reason = "drop_palindromic = FALSE")
  } else if (n_in == 0L) {
    add_step("palindromic", 0L, 0L, dropped = TRUE, n_dropped = 0L)
  } else if (!all(c("effect_allele.exposure", "other_allele.exposure") %in% names(snps_tbl))) {
    logger::log_warn("preprocess: step 6 (palindromic) skipped: allele columns missing")
    add_step("palindromic", n_in, n_in, dropped = TRUE,
             skipped_reason = "allele columns missing")
  } else {
    eaf_vec <- if ("eaf.exposure" %in% names(snps_tbl)) snps_tbl$eaf.exposure else rep(NA_real_, n_in)
    ambig <- is_ambiguous_palindrome(snps_tbl$effect_allele.exposure,
                                     snps_tbl$other_allele.exposure,
                                     eaf_vec)
    ambig[is.na(ambig)] <- FALSE
    snps_tbl <- snps_tbl[!ambig, , drop = FALSE]
    if (verbose) {
      logger::log_info(
        "preprocess: step 6 (palindromic, ambiguous-only) {n_in} -> {nrow(snps_tbl)}"
      )
    }
    add_step("palindromic", n_in, nrow(snps_tbl), dropped = TRUE,
             n_dropped = as.integer(sum(ambig)))
  }

  # ---- Step 7: LD clumping (always local, ancestry-aware) ----
  n_in <- nrow(snps_tbl)
  if (isTRUE(co$already_clumped)) {
    if (verbose) logger::log_info("preprocess: step 7 (clump) skipped (already_clumped)")
    add_step("clump", n_in, n_in, r2 = co$r2, kb = co$kb, skipped_reason = "already_clumped")
  } else if (n_in == 0L) {
    add_step("clump", 0L, 0L, r2 = co$r2, kb = co$kb, skipped_reason = "no input rows")
  } else if (!"pval.exposure" %in% names(snps_tbl)) {
    logger::log_warn("preprocess: step 7 (clump) skipped: 'pval.exposure' column missing")
    add_step("clump", n_in, n_in, r2 = co$r2, kb = co$kb, skipped_reason = "no pval.exposure column")
  } else {
    if (!requireNamespace("ieugwasr", quietly = TRUE)) {
      stop("preprocess_exposure_snps(): clumping requires the ieugwasr package.", call. = FALSE)
    }
    ld_pop <- ld_pop_from_ancestry(ancestry)
    ld_ref_dir <- file.path(cache_dir, "ld_reference")
    ref <- ld_reference_checker(pop = ld_pop, ld_ref_dir = ld_ref_dir,
                                verbose = verbose, confirm = confirm)
    plink_bin <- tryCatch(.find_plink_bin(NULL), error = function(e) NULL)
    if (verbose) {
      logger::log_info(
        "preprocess: step 7 (clump) ld_clump on {n_in} SNPs (pop={ld_pop}, r2={co$r2}, kb={co$kb})"
      )
    }
    pvals <- suppressWarnings(as.numeric(snps_tbl$pval.exposure))
    clump_input <- tibble::tibble(rsid = as.character(snps_tbl$SNP), pval = pvals)
    clumped <- tryCatch(
      ieugwasr::ld_clump(clump_input, clump_r2 = co$r2, clump_kb = co$kb,
                         plink_bin = plink_bin, bfile = ref$bfile_prefix, pop = ld_pop),
      error = function(e) {
        logger::log_warn("preprocess: ld_clump failed: {conditionMessage(e)}")
        NULL
      }
    )
    if (is.null(clumped) || !nrow(clumped)) {
      snps_tbl <- snps_tbl[0L, , drop = FALSE]
      add_step("clump", n_in, 0L, r2 = co$r2, kb = co$kb,
               skipped_reason = "ld_clump failed or returned 0 rows")
    } else {
      snps_tbl <- snps_tbl[as.character(snps_tbl$SNP) %in% as.character(clumped$rsid), ,
                           drop = FALSE]
      add_step("clump", n_in, nrow(snps_tbl), r2 = co$r2, kb = co$kb)
    }
    if (verbose) logger::log_info("preprocess: step 7 (clump) {n_in} -> {nrow(snps_tbl)}")
  }

  # ---- Step 8: F-stat (per-SNP; safe to run after clumping) ----
  n_in <- nrow(snps_tbl)
  if (n_in == 0L) {
    add_step("f_stat", 0L, 0L, value = co$f_threshold)
  } else if (!all(c("beta.exposure", "se.exposure") %in% names(snps_tbl))) {
    logger::log_warn("preprocess: step 8 (f_stat) skipped: beta.exposure/se.exposure missing")
    add_step("f_stat", n_in, n_in, value = co$f_threshold,
             skipped_reason = "beta.exposure/se.exposure columns missing")
  } else {
    fs <- f_stat(suppressWarnings(as.numeric(snps_tbl$beta.exposure)),
                 suppressWarnings(as.numeric(snps_tbl$se.exposure)))
    keep <- is.finite(fs) & fs >= co$f_threshold
    snps_tbl <- snps_tbl[keep, , drop = FALSE]
    if (verbose) {
      logger::log_info("preprocess: step 8 (f_stat>={co$f_threshold}) {n_in} -> {nrow(snps_tbl)}")
    }
    add_step("f_stat", n_in, nrow(snps_tbl), value = co$f_threshold)
  }

  if (verbose) {
    logger::log_info("preprocess: complete -- {nrow(snps_tbl)} IVs survived")
  }

  list(
    snps          = tibble::as_tibble(snps_tbl),
    steps         = steps_acc,
    resolved_opts = co
  )
}
