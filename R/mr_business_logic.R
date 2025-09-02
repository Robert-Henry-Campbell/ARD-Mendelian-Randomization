#' Harmonise, run MR methods, sensitivity & QC across outcomes (no auto-format of exposures)
#' @param MR_df tibble with `outcome_snps` list-column (TwoSampleMR-formatted)
#' @param exposure_snps exposure instruments (already asserted/validated; must have *.exposure cols)
#' @param sensitivity_enabled character vector in
#'   c("egger_intercept","egger_slope_agreement","weighted_median","weighted_mode",
#'     "steiger_direction","leave_one_out","ivw_Q","ivw_I2")
#' @param sensitivity_pass_min integer; minimum # of checks that must pass
#' @param scatterplot,snpforestplot,leaveoneoutplot logical; write per-outcome plots if TRUE
#' @param plot_output_dir directory to write plots ("" = off)
#' @param cache_dir path (unused here; kept for parity)
#' @param verbose logical
#' @return list(MR_df=updated, results_df=tidy summary)
#' @export
mr_business_logic <- function(
    MR_df, exposure_snps,
    sensitivity_enabled,
    sensitivity_pass_min,
    scatterplot, snpforestplot, leaveoneoutplot,
    plot_output_dir, cache_dir = ardmr_cache_dir(), verbose = TRUE
) {
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
    stop("mr_business_logic() requires the TwoSampleMR package.", call. = FALSE)
  }
  TS <- asNamespace("TwoSampleMR")

  # ---- guards ----
  need_expo <- c("SNP","effect_allele.exposure","other_allele.exposure",
                 "beta.exposure","se.exposure","id.exposure")
  if (!all(need_expo %in% names(exposure_snps))) {
    stop("Exposure does not contain required *.exposure columns: ",
         paste(setdiff(need_expo, names(exposure_snps)), collapse=", "), call. = FALSE)
  }
  if (!"outcome_snps" %in% names(MR_df)) {
    stop("MR_df has no 'outcome_snps' list-column.", call. = FALSE)
  }

  # ---- helpers ----
  .slug <- function(x) { x <- as.character(x); x[is.na(x)|!nzchar(x)] <- "NA"
  x <- gsub("[^A-Za-z0-9._-]+","_",x); substr(x,1,120) }
  .nz <- function(x) !is.null(x) && length(x) > 0
  .as_num <- function(x) suppressWarnings(as.numeric(x))

  #quick naming of exposure columns
  exposure_snps <- fix_exposure_names(exposure_snps)

  # We won't auto-format exposures any more, but add a label if needed (for plots/titles)
  exposure_fmt <- tibble::as_tibble(exposure_snps)
  if (!"exposure" %in% names(exposure_fmt) && "Exposure" %in% names(exposure_fmt)) {
    exposure_fmt$exposure <- exposure_fmt$Exposure
  }

  # prepare storage
  if (!"mr_summary"      %in% names(MR_df)) MR_df$mr_summary      <- vector("list", nrow(MR_df))
  if (!"heterogeneity"   %in% names(MR_df)) MR_df$heterogeneity   <- vector("list", nrow(MR_df))
  if (!"egger_intercept" %in% names(MR_df)) MR_df$egger_intercept <- vector("list", nrow(MR_df))
  if (!"leaveoneout"     %in% names(MR_df)) MR_df$leaveoneout     <- vector("list", nrow(MR_df))
  if (!"steiger"         %in% names(MR_df)) MR_df$steiger         <- vector("list", nrow(MR_df))
  if (!"harmonised"      %in% names(MR_df)) MR_df$harmonised      <- vector("list", nrow(MR_df))
  if (!"plots"           %in% names(MR_df)) MR_df$plots           <- vector("list", nrow(MR_df))

  results_rows <- vector("list", nrow(MR_df))

  # pre-allocate result columns to silence "Unknown or uninitialised column" warnings
  if (!"qc_pass"      %in% names(MR_df)) MR_df$qc_pass      <- rep(NA, nrow(MR_df))
  if (!"nsnp_before"  %in% names(MR_df)) MR_df$nsnp_before  <- rep(NA_integer_, nrow(MR_df))
  if (!"nsnp_after"   %in% names(MR_df)) MR_df$nsnp_after   <- rep(NA_integer_, nrow(MR_df))
  if (!"ivw_beta"     %in% names(MR_df)) MR_df$ivw_beta     <- rep(NA_real_,    nrow(MR_df))
  if (!"ivw_se"       %in% names(MR_df)) MR_df$ivw_se       <- rep(NA_real_,    nrow(MR_df))
  if (!"ivw_p"        %in% names(MR_df)) MR_df$ivw_p        <- rep(NA_real_,    nrow(MR_df))

  # progress reporting
  total <- nrow(MR_df)
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)

  for (i in seq_len(nrow(MR_df))) {

    #progress bar reporting
    utils::setTxtProgressBar(pb, i)
    if (i %% 10 == 0 || i == total) {
      remaining <- total - i
      logger::log_info("MR business: processed {i}/{total} outcomes; {remaining} to go")
    }


    rec <- MR_df[i, , drop = FALSE]
    outdat <- rec$outcome_snps[[1]]
    if (!.nz(outdat) || nrow(outdat) == 0) {
      if (verbose) logger::log_warn("MR business: outcome row {i} has no outcome_snps; skipping")
      next
    }

    # Outcome label preference: 'outcome' (from format_data), else 'Phenotype', else MR_df$description
    outcome_label <- if ("outcome" %in% names(outdat)) outdat$outcome[1] else
      if ("Phenotype" %in% names(outdat)) outdat$Phenotype[1] else
        if ("description" %in% names(rec)) as.character(rec$description[[1]]) else sprintf("outcome_%03d", i)
    outcome_label <- ifelse(is.na(outcome_label) | !nzchar(outcome_label),
                            sprintf("outcome_%03d", i), outcome_label)

    # 1) Harmonise (action=2 to handle palindromes using EAF where available)
    nsnp_start <- length(unique(exposure_fmt$SNP))
    hdat <- tryCatch(
      TS$harmonise_data(exposure_dat = exposure_fmt, outcome_dat = outdat, action = 2),
      error = function(e) {
        if (verbose) logger::log_warn("Harmonise failed for '{outcome_label}': {conditionMessage(e)}")
        NULL
      }
    )
    if (is.null(hdat) || nrow(hdat) == 0) {
      MR_df$harmonised[[i]] <- tibble::tibble()
      results_rows[[i]] <- tibble::tibble(
        outcome = outcome_label,
        nsnp_before = nsnp_start, nsnp_after = 0L,
        beta_ivw = NA_real_, se_ivw = NA_real_, p_ivw = NA_real_,
        Q_ivw = NA_real_, Q_df_ivw = NA_real_, Q_p_ivw = NA_real_, I2_ivw = NA_real_,
        beta_egger = NA_real_, se_egger = NA_real_, p_egger = NA_real_,
        egger_intercept = NA_real_, egger_intercept_p = NA_real_, I2GX = NA_real_,
        beta_wmed = NA_real_, se_wmed = NA_real_, p_wmed = NA_real_,
        beta_wmode = NA_real_, se_wmode = NA_real_, p_wmode = NA_real_,
        loo_flip = NA, steiger_fail_frac = NA_real_,
        checks_attempted = 0L, checks_passed = 0L, qc_pass = FALSE
      )
      next
    }
    # Respect mr_keep if present
    hdat_use <- if ("mr_keep" %in% names(hdat)) hdat[hdat$mr_keep %in% TRUE, , drop = FALSE] else hdat
    MR_df$harmonised[[i]] <- hdat_use
    nsnp_after <- length(unique(hdat_use$SNP))
    if (nsnp_after < 1) {
      if (verbose) logger::log_warn("No overlapping instruments after harmonisation for '{outcome_label}'")
    }

    # sensitivity bookkeeping
    enabled <- intersect(
      sensitivity_enabled,
      c("egger_intercept","egger_slope_agreement","weighted_median","weighted_mode",
        "steiger_direction","leave_one_out","ivw_Q","ivw_I2")
    )
    checks_attempted <- 0L; checks_passed <- 0L
    check_pass <- function(ok) {
      if (!is.na(ok)) { checks_attempted <<- checks_attempted + 1L; if (isTRUE(ok)) checks_passed <<- checks_passed + 1L }
      ok
    }

    # 1) Main MR
    method_list <- if (nsnp_after == 1) "mr_wald_ratio"
    else c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
    mr_res <- tryCatch(TS$mr(hdat_use, method_list = method_list), error = function(e) NULL)
    MR_df$mr_summary[[i]] <- mr_res

    pick_est <- function(df, pattern) {
      if (is.null(df) || nrow(df) == 0) return(c(b = NA_real_, se = NA_real_, p = NA_real_, nsnp = NA_integer_))
      row <- df[grepl(pattern, df$method, ignore.case = TRUE), , drop = FALSE]
      if (nrow(row) == 0) row <- df[df$method == pattern, , drop = FALSE]
      if (nrow(row) == 0) return(c(b = NA_real_, se = NA_real_, p = NA_real_, nsnp = NA_integer_))
      c(b = .as_num(row$b[1]), se = .as_num(row$se[1]), p = .as_num(row$pval[1]), nsnp = as.integer(row$nsnp[1]))
    }
    ivw   <- pick_est(mr_res, "Inverse variance weighted|mr_ivw")
    egger <- pick_est(mr_res, "Egger|egger_regression")
    wmed  <- pick_est(mr_res, "Weighted median|weighted_median")
    wmode <- pick_est(mr_res, "Weighted mode|weighted_mode")

    # Optional plots
    plot_files <- list()
    if (nzchar(plot_output_dir)) {
      out_dir <- file.path(plot_output_dir, .slug(outcome_label))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      if (scatterplot && !is.null(mr_res) && nrow(mr_res) > 0) {
        p_sc <- tryCatch(TS$mr_scatter_plot(mr_res, hdat_use), error = function(e) NULL)
        if (.nz(p_sc)) {
          f <- file.path(out_dir, "scatter.png")
          try(ggplot2::ggsave(f, p_sc[[1]], width = 6, height = 5, dpi = 150), silent = TRUE)
          plot_files$scatter <- f
        }
      }
      if (snpforestplot && nsnp_after > 1) {
        single <- tryCatch(TS$mr_singlesnp(hdat_use), error = function(e) NULL)
        if (.nz(single)) {
          p_forest <- tryCatch(TS$mr_forest_plot(single), error = function(e) NULL)
          if (.nz(p_forest)) {
            f <- file.path(out_dir, "forest.png")
            try(ggplot2::ggsave(f, p_forest[[1]], width = 6, height = 7, dpi = 150), silent = TRUE)
            plot_files$forest <- f
          }
        }
      }
      if (leaveoneoutplot && nsnp_after > 1) {
        loo_tmp <- tryCatch(TS$mr_leaveoneout(hdat_use), error = function(e) NULL)
        if (.nz(loo_tmp)) {
          p_loo <- tryCatch(TS$mr_leaveoneout_plot(loo_tmp), error = function(e) NULL)
          if (.nz(p_loo)) {
            f <- file.path(out_dir, "leaveoneout.png")
            try(ggplot2::ggsave(f, p_loo[[1]], width = 6, height = 5, dpi = 150), silent = TRUE)
            plot_files$leaveoneout <- f
          }
        }
      }
    }
    MR_df$plots[[i]] <- plot_files

    # 1d) IVW heterogeneity (Q, I2)
    het <- tryCatch(TS$mr_heterogeneity(hdat_use, method_list = c("mr_ivw","mr_egger_regression")), error = function(e) NULL)
    MR_df$heterogeneity[[i]] <- het
    Q_ivw <- Q_df_ivw <- Q_p_ivw <- I2_ivw <- NA_real_
    if (.nz(het)) {
      row_ivw <- het[grepl("Inverse variance weighted|mr_ivw", het$method, ignore.case = TRUE), , drop = FALSE]
      if (nrow(row_ivw)) {
        Q_ivw    <- .as_num(row_ivw$Q[1])
        Q_df_ivw <- .as_num(row_ivw$Q_df[1])
        Q_p_ivw  <- .as_num(row_ivw$Q_pval[1])
        if (!is.na(Q_ivw) && !is.na(Q_df_ivw) && Q_ivw > 0) {
          I2_ivw <- max(0, (Q_ivw - Q_df_ivw) / Q_ivw) * 100
        }
      }
    }
    chk_Q  <- if ("ivw_Q" %in% enabled)  check_pass(!(is.na(Q_p_ivw) || Q_p_ivw < 0.05)) else NA
    chk_I2 <- if ("ivw_I2" %in% enabled) check_pass(!(is.na(I2_ivw) || I2_ivw > 50)) else NA

    # 1e) Egger: slope & intercept; I2GX for precision
    egger_int <- tryCatch(TS$mr_pleiotropy_test(hdat_use), error = function(e) NULL)
    MR_df$egger_intercept[[i]] <- egger_int
    egger_int_est <- egger_int_p <- I2GX <- NA_real_
    if (.nz(egger_int)) {
      egger_int_est <- .as_num(egger_int$estimate[1])
      egger_int_p   <- .as_num(egger_int$pval[1])
    }
    if (nsnp_after > 1 && all(c("beta.exposure","se.exposure") %in% names(hdat_use))) {
      I2GX_raw <- tryCatch(TS$Isq(hdat_use$beta.exposure, hdat_use$se.exposure), error = function(e) NA_real_)
      if (is.list(I2GX_raw) && "I2" %in% names(I2GX_raw)) I2GX_raw <- I2GX_raw$I2
      I2GX <- suppressWarnings(as.numeric(I2GX_raw))
      if (!is.na(I2GX) && I2GX > 1) I2GX <- I2GX / 100
    }
    chk_Eint  <- if ("egger_intercept" %in% enabled) check_pass(!(is.na(egger_int_p) || egger_int_p < 0.05)) else NA
    slope_agree <- (sign(ivw["b"]) == sign(egger["b"]))
    if (is.na(slope_agree)) slope_agree <- NA
    if (!is.na(I2GX) && I2GX < 0.60) {
      chk_Eslope <- if ("egger_slope_agreement" %in% enabled) check_pass(TRUE) else NA
    } else {
      chk_Eslope <- if ("egger_slope_agreement" %in% enabled) check_pass(isTRUE(slope_agree)) else NA
    }

    # 1f/1g) Weighted median/mode checks (sign agreement + |Δβ| <= 1*SE_IVW)
    compare_to_ivw <- function(b_alt) {
      if (is.na(ivw["b"]) || is.na(ivw["se"]) || is.na(b_alt)) return(NA)
      agree  <- sign(b_alt) == sign(ivw["b"])
      within <- abs(b_alt - ivw["b"]) <= as.numeric(ivw["se"])
      isTRUE(agree && within)
    }
    chk_Wmed  <- if ("weighted_median" %in% enabled) check_pass(compare_to_ivw(wmed["b"]))  else NA
    chk_Wmode <- if ("weighted_mode"   %in% enabled) check_pass(compare_to_ivw(wmode["b"])) else NA

    # 2) Leave-one-out: IVW dropping each SNP
    loo <- NULL; chk_LOO <- NA
    if (nsnp_after > 1) {
      loo <- tryCatch(TS$mr_leaveoneout(hdat_use), error = function(e) NULL)
      MR_df$leaveoneout[[i]] <- loo
      if (.nz(loo)) {
        ivw_sign <- sign(ivw["b"])
        if (!is.na(ivw_sign)) {
          flips <- any(sign(.as_num(loo$b)) != ivw_sign, na.rm = TRUE)
          chk_LOO <- if ("leave_one_out" %in% enabled) check_pass(!flips) else NA
        }
      }
    } else {
      MR_df$leaveoneout[[i]] <- tibble::tibble()
    }

    # 3) Steiger directionality
    steiger_dat <- tryCatch(TS$steiger_filtering(hdat_use), error = function(e) NULL)
    MR_df$steiger[[i]] <- if (!is.null(steiger_dat)) steiger_dat else tibble::tibble()
    steiger_fail_frac <- NA_real_
    if (.nz(steiger_dat) && "steiger_dir" %in% names(steiger_dat)) {
      frac_fail <- mean(!as.logical(steiger_dat$steiger_dir), na.rm = TRUE)
      if (is.nan(frac_fail)) frac_fail <- NA_real_
      steiger_fail_frac <- frac_fail
      chk_Steiger <- if ("steiger_direction" %in% enabled) check_pass(!(is.na(frac_fail) || frac_fail > 0.30)) else NA
    } else {
      chk_Steiger <- NA
    }

    # QC decision
    core_viable <- is.finite(as.numeric(ivw["b"])) && is.finite(as.numeric(ivw["se"]))
    qc_ok <- (checks_passed >= sensitivity_pass_min) && core_viable

    # results row
    results_rows[[i]] <- tibble::tibble(
      outcome = outcome_label,
      nsnp_before = nsnp_start,
      nsnp_after  = nsnp_after,
      beta_ivw = as.numeric(ivw["b"]),  se_ivw = as.numeric(ivw["se"]),  p_ivw = as.numeric(ivw["p"]),
      Q_ivw = Q_ivw, Q_df_ivw = Q_df_ivw, Q_p_ivw = Q_p_ivw, I2_ivw = I2_ivw,
      beta_egger = as.numeric(egger["b"]), se_egger = as.numeric(egger["se"]), p_egger = as.numeric(egger["p"]),
      egger_intercept = egger_int_est, egger_intercept_p = egger_int_p, I2GX = I2GX,
      beta_wmed = as.numeric(wmed["b"]),  se_wmed = as.numeric(wmed["se"]),  p_wmed = as.numeric(wmed["p"]),
      beta_wmode = as.numeric(wmode["b"]), se_wmode = as.numeric(wmode["se"]), p_wmode = as.numeric(wmode["p"]),
      loo_flip = if (is.na(chk_LOO)) NA else !chk_LOO,
      steiger_fail_frac = steiger_fail_frac,
      checks_attempted = checks_attempted,
      checks_passed = checks_passed,
      qc_pass = qc_ok
    )

    # compact summary back on MR_df
    MR_df$qc_pass[i]     <- qc_ok
    MR_df$nsnp_before[i] <- nsnp_start
    MR_df$nsnp_after[i]  <- nsnp_after
    MR_df$ivw_beta[i]    <- as.numeric(ivw["b"])
    MR_df$ivw_se[i]      <- as.numeric(ivw["se"])
    MR_df$ivw_p[i]       <- as.numeric(ivw["p"])
  }

  results_df <- dplyr::bind_rows(results_rows)
  if (verbose) logger::log_info("MR business logic: {nrow(results_df)} outcomes analysed")
  list(MR_df = MR_df, results_df = results_df)
}
