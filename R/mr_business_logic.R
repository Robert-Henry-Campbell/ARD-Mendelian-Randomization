#R/mr_business_logic.R
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
    plot_output_dir,
    cache_dir = ardmr_cache_dir(),
    verbose = TRUE,
    test = FALSE
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
  .nz <- function(x) {
    if (is.null(x)) return(FALSE)
    if (is.data.frame(x)) return(nrow(x) > 0)  # require at least one row
    length(x) > 0
  }
  .as_num <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.factor(x)) x <- as.character(x)
    suppressWarnings(as.numeric(x))
  }

  .save_pub_plot <- function(plot, basepath, w_mm = 89, h_mm = 80, dpi = 500) {
    # High-res PNG for quick viewing
    ggplot2::ggsave(paste0(basepath, ".png"), plot,
                    width = w_mm, height = h_mm, units = "mm",
                    dpi = dpi, type = "cairo")         # better antialiasing
    # Vector PDF for submission
    ggplot2::ggsave(paste0(basepath, ".pdf"), plot,
                    width = w_mm, height = h_mm, units = "mm",
                    device = cairo_pdf)
  }


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

  # ---- pre-allocate all results_* columns on MR_df ----
  add_col <- function(nm, prototype) {
    if (!nm %in% names(MR_df)) MR_df[[nm]] <<- rep(prototype, nrow(MR_df))
  }
  add_col("results_outcome",              NA_character_)
  add_col("results_nsnp_before",          NA_integer_)
  add_col("results_nsnp_after",           NA_integer_)
  add_col("results_beta_ivw",             NA_real_)
  add_col("results_se_ivw",               NA_real_)
  add_col("results_log10p_ivw",           NA_real_)
  add_col("results_p_ivw",                NA_real_)
  add_col("results_Q_ivw",                NA_real_)
  add_col("results_Q_df_ivw",             NA_real_)
  add_col("results_Q_p_ivw",              NA_real_)
  add_col("results_I2_ivw",               NA_real_)
  add_col("results_beta_egger",           NA_real_)
  add_col("results_se_egger",             NA_real_)
  add_col("results_p_egger",              NA_real_)
  add_col("results_egger_intercept",      NA_real_)
  add_col("results_egger_intercept_p",    NA_real_)
  add_col("results_I2GX",                 NA_real_)
  add_col("results_beta_wmed",            NA_real_)
  add_col("results_se_wmed",              NA_real_)
  add_col("results_p_wmed",               NA_real_)
  add_col("results_beta_wmode",           NA_real_)
  add_col("results_se_wmode",             NA_real_)
  add_col("results_p_wmode",              NA_real_)
  add_col("results_loo_flip",             NA)            # logical
  add_col("results_steiger_fail_frac",    NA_real_)
  add_col("results_checks_attempted",     NA_integer_)
  add_col("results_checks_passed",        NA_integer_)
  add_col("results_qc_pass",              NA)            # logical

  # just before creating the progress bar
  total <- nrow(MR_df)

  # build the index of rows to run
  idx <- seq_len(total)
  if (isTRUE(test)) idx <- idx[seq_len(min(5L, total))]
  n_to_run <- length(idx)
  # progress bar uses the number we will actually run
  pb <- utils::txtProgressBar(min = 0, max = n_to_run, style = 3)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)


  for (k in seq_along(idx)) {
    i <- idx[k]  # actual MR_df row we’re working on

    utils::setTxtProgressBar(pb, k)
    if (k %% 10 == 0 || k == n_to_run) {
      remaining <- n_to_run - k
      logger::log_info("MR business: processed {k}/{n_to_run} outcomes; {remaining} to go")
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
      # record results_* directly on MR_df
      MR_df$results_outcome[i]           <- outcome_label
      MR_df$results_nsnp_before[i]       <- nsnp_start
      MR_df$results_nsnp_after[i]        <- 0L
      MR_df$results_checks_attempted[i]  <- 0L
      MR_df$results_checks_passed[i]     <- 0L
      MR_df$results_qc_pass[i]           <- FALSE
      next
    }


    # Respect mr_keep if present (drop NAs)
    hdat_use <- if ("mr_keep" %in% names(hdat)) hdat[!is.na(hdat$mr_keep) & hdat$mr_keep, , drop = FALSE] else hdat
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
    method_list <- if (nsnp_after == 1) "mr_wald_ratio" else
      c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
    mr_res <- tryCatch(TS$mr(hdat_use, method_list = method_list), error = function(e) NULL)
    MR_df$mr_summary[[i]] <- mr_res

    pick_est <- function(df, pattern) {
      if (is.null(df) || nrow(df) == 0) return(c(b = NA_real_, se = NA_real_, p = NA_real_, nsnp = NA_integer_))
      row <- df[grepl(pattern, df$method, ignore.case = TRUE), , drop = FALSE]
      if (nrow(row) == 0) row <- df[df$method == pattern, , drop = FALSE]
      if (nrow(row) == 0) return(c(b = NA_real_, se = NA_real_, p = NA_real_, nsnp = NA_integer_))
      c(b = .as_num(row$b[1]), se = .as_num(row$se[1]), p = .as_num(row$pval[1]), nsnp = as.integer(row$nsnp[1]))
    }

    ivw_pattern <- if (nsnp_after == 1) "Wald ratio|mr_wald_ratio" else "Inverse variance weighted|mr_ivw"
    ivw   <- pick_est(mr_res, ivw_pattern)
    egger <- pick_est(mr_res, "Egger|egger_regression")
    wmed  <- pick_est(mr_res, "Weighted median|weighted_median")
    wmode <- pick_est(mr_res, "Weighted mode|weighted_mode")

    # Optional plots
    plot_files <- list()
    if (nzchar(plot_output_dir)) {
      out_dir <- file.path(plot_output_dir, .slug(outcome_label))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      # --- SCATTER ---
      if (scatterplot && !is.null(mr_res) && nrow(mr_res) > 0) {
        p_sc <- tryCatch(TS$mr_scatter_plot(mr_res, hdat_use), error = function(e) NULL)
        if (.nz(p_sc)) {
          # save to a STEM (no extension) to avoid ".png.png" if your helper appends one
          f_sc <- file.path(out_dir, "scatter")
          .save_pub_plot(p_sc[[1]], f_sc, w_mm = 92, h_mm = 100)
          plot_files$scatter <- paste0(f_sc, ".png")
        }
      }

      # --- FOREST ---
      if (snpforestplot && nsnp_after > 1) {
        single <- tryCatch(TS$mr_singlesnp(hdat_use), error = function(e) NULL)
        if (.nz(single)) {
          # shorten pooled label in the data
          single$SNP <- gsub(
            "^All\\s*-\\s*Inverse variance weighted(?:\\s*\\((?:fixed|random) effects\\))?$",
            "All - IVW",
            single$SNP
          )

          p_forest <- tryCatch(TS$mr_forest_plot(single), error = function(e) NULL)
          if (.nz(p_forest)) {
            # also relabel on the axis + wrap long SNP labels so they don't run off
            pf <- p_forest[[1]] +
              ggplot2::scale_y_discrete(
                labels = if (requireNamespace("scales", quietly = TRUE))
                  scales::label_wrap(28) else function(x) x
              ) +
              ggplot2::scale_y_discrete(labels = function(x) gsub(
                "^All\\s*-\\s*Inverse variance weighted(?:\\s*\\((?:fixed|random) effects\\))?$",
                "All - IVW", x
              ))

            f_forest <- file.path(out_dir, "forest")
            # IMPORTANT: save the forest object, not the scatter
            .save_pub_plot(pf, f_forest, h_mm = 80, w_mm = 89)
            plot_files$forest <- paste0(f_forest, ".png")
          }
        }
      }

      # --- LEAVE-ONE-OUT ---
      if (leaveoneoutplot && nsnp_after > 1) {
        loo_tmp <- tryCatch(TS$mr_leaveoneout(hdat_use), error = function(e) NULL)
        if (.nz(loo_tmp)) {
          p_loo <- tryCatch(TS$mr_leaveoneout_plot(loo_tmp), error = function(e) NULL)
          if (.nz(p_loo)) {
            f_loo <- file.path(out_dir, "leaveoneout")  # stem (no extension)
            .save_pub_plot(p_loo[[1]], f_loo)           # keep your defaults
            plot_files$leaveoneout <- paste0(f_loo, ".png")
          }
        }
      }
    }
    MR_df$plots[[i]] <- plot_files

    # 1d) IVW heterogeneity (Q, I2)
    #testing a potential issue with null het shrinking row size
    het <- if (nsnp_after > 1) {
      tryCatch(TS$mr_heterogeneity(hdat_use, method_list = c("mr_ivw","mr_egger_regression")),
               error = function(e) NULL)
    } else {
      NULL
    }

    # ensure the stored value is never a raw NULL (which would delete the element)
    if (is.null(het)) het <- tibble::tibble()

    # always assign list-cols with [i] <- list(...)
    MR_df$heterogeneity[i] <- list(het)


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
    chk_Q  <- if (nsnp_after > 1 && "ivw_Q" %in% enabled)  check_pass(!(is.na(Q_p_ivw) || Q_p_ivw < 0.05)) else NA
    chk_I2 <- if (nsnp_after > 1 && "ivw_I2" %in% enabled) check_pass(!(is.na(I2_ivw) || I2_ivw > 50)) else NA

    # 1e) Egger: slope & intercept; I2GX for precision
    egger_int <- tryCatch(TS$mr_pleiotropy_test(hdat_use), error = function(e) NULL)
    if (is.null(egger_int)) egger_int <- tibble::tibble()   # ensure tibble, not NULL
    MR_df$egger_intercept[[i]] <- egger_int

    egger_int_est <- NA_real_
    egger_int_p   <- NA_real_
    I2GX          <- NA_real_

    # use the columns you actually have: egger_intercept, pval
    if (nrow(egger_int) > 0) {
      egger_int_est <- .as_num(egger_int$egger_intercept[1])
      egger_int_p   <- .as_num(egger_int$pval[1])
    }

    if (nsnp_after > 1 && all(c("beta.exposure","se.exposure") %in% names(hdat_use))) {
      I2GX_raw <- tryCatch(TS$Isq(hdat_use$beta.exposure, hdat_use$se.exposure), error = function(e) NA_real_)
      if (is.list(I2GX_raw) && "I2" %in% names(I2GX_raw)) I2GX_raw <- I2GX_raw$I2
      I2GX <- suppressWarnings(as.numeric(I2GX_raw))
      if (!is.na(I2GX) && I2GX >= 1) I2GX <- I2GX / 100
    }
    chk_Eint  <- if ("egger_intercept" %in% enabled) check_pass(!(is.na(egger_int_p) || egger_int_p < 0.05)) else NA
    slope_agree <- (sign(ivw["b"]) == sign(egger["b"]))
    if (!is.na(I2GX) && I2GX < 0.60) {
      chk_Eslope <- if ("egger_slope_agreement" %in% enabled) check_pass(TRUE) else NA
    } else if (!is.na(slope_agree)) {
      chk_Eslope <- if ("egger_slope_agreement" %in% enabled) check_pass(slope_agree) else NA
    } else {
      chk_Eslope <- NA
    }


    # 1f/1g) Weighted median/mode checks (sign agreement + |Δβ| <= 1*SE_IVW)
    compare_to_ivw <- function(b_alt) {
      if (is.na(ivw["b"]) || is.na(ivw["se"]) || is.na(b_alt)) return(NA)
      agree  <- sign(b_alt) == sign(ivw["b"])
      within <- abs(b_alt - ivw["b"]) <= as.numeric(ivw["se"])
      if (is.na(agree) || is.na(within)) return(NA)
      agree && within
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
    qc_ok <- if (nsnp_after == 1) {
      core_viable
    } else {
      (checks_passed >= sensitivity_pass_min) && core_viable
    }

    # write all results_* fields directly on MR_df
    MR_df$results_outcome[i]            <- outcome_label
    MR_df$results_nsnp_before[i]        <- nsnp_start
    MR_df$results_nsnp_after[i]         <- nsnp_after
    b  <- as.numeric(ivw["b"])
    se <- as.numeric(ivw["se"])
    z  <- b / se
    # two-sided log10(p); still NEGATIVE for small p
    log10p <- (stats::pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10)
    MR_df$results_beta_ivw[i]       <- b
    MR_df$results_se_ivw[i]         <- se
    MR_df$results_log10p_ivw[i]     <- log10p
    MR_df$results_p_ivw[i]          <- 10^log10p
    MR_df$results_Q_ivw[i]              <- Q_ivw
    MR_df$results_Q_df_ivw[i]           <- Q_df_ivw
    MR_df$results_Q_p_ivw[i]            <- Q_p_ivw
    MR_df$results_I2_ivw[i]             <- I2_ivw
    MR_df$results_beta_egger[i]         <- as.numeric(egger["b"])
    MR_df$results_se_egger[i]           <- as.numeric(egger["se"])
    MR_df$results_p_egger[i]            <- as.numeric(egger["p"])
    MR_df$results_egger_intercept[i]    <- egger_int_est
    MR_df$results_egger_intercept_p[i]  <- egger_int_p
    MR_df$results_I2GX[i]               <- I2GX
    MR_df$results_beta_wmed[i]          <- as.numeric(wmed["b"])
    MR_df$results_se_wmed[i]            <- as.numeric(wmed["se"])
    MR_df$results_p_wmed[i]             <- as.numeric(wmed["p"])
    MR_df$results_beta_wmode[i]         <- as.numeric(wmode["b"])
    MR_df$results_se_wmode[i]           <- as.numeric(wmode["se"])
    MR_df$results_p_wmode[i]            <- as.numeric(wmode["p"])
    MR_df$results_loo_flip[i]           <- if (is.na(chk_LOO)) NA else !chk_LOO
    MR_df$results_steiger_fail_frac[i]  <- steiger_fail_frac
    MR_df$results_checks_attempted[i]   <- checks_attempted
    MR_df$results_checks_passed[i]      <- checks_passed
    MR_df$results_qc_pass[i]            <- qc_ok

  }

  # derive results_df from MR_df's results_* columns (+ carry group/meta)
  results_cols <- grep("^results_", names(MR_df), value = TRUE)
  group_meta <- intersect(
    c("Cause Name", "Cause Hierarchy Level", "ICD10",
      "cause_level_1", "cause_level_2", "cause_level_3", "cause_level_4",
      "icd10_explo", "pheno_sex", "description", "plots","ARD_selected"),
    names(MR_df)
  )
  results_df <- tibble::as_tibble(
    dplyr::bind_cols(
      MR_df[, group_meta, drop = FALSE],
      MR_df[, results_cols, drop = FALSE]
    )
  )

  if (verbose) logger::log_info("MR business logic: {nrow(results_df)} outcomes analysed")
  list(MR_df = MR_df, results_df = results_df)

}
