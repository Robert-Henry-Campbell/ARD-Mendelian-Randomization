# R/plots.R (UNTESTED)

#' Manhattan plot over outcomes (IVW p-values on QC-pass outcomes)
#' @param results_df tidy results (must contain p_ivw; ideally qc_pass)
#' @param Multiple_testing_correction "BH" or "bonferroni"
#' @param alpha significance level for corrections (default 0.05)
#' @param verbose if TRUE, emit step-by-step logs via {logger}
#' @export
manhattan_plot <- function(results_df,
                           Multiple_testing_correction = c("BH","bonferroni"),
                           alpha = 0.05,
                           verbose = TRUE) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)

  if (is.null(results_df) || !nrow(results_df)) {
    if (verbose) logger::log_warn("Manhattan: input results_df is NULL or has 0 rows; returning placeholder plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no results)"))
  }

  if (verbose) logger::log_info("Manhattan: starting with {nrow(results_df)} rows; MTC = {Multiple_testing_correction}, alpha = {alpha}")

  df <- tibble::as_tibble(results_df)

  # keep only QC-pass (if column missing, assume all pass)
  if (!"qc_pass" %in% names(df)) {
    if (verbose) logger::log_info("Manhattan: 'qc_pass' column not present; assuming all pass.")
    df$qc_pass <- TRUE
  }
  n_before <- nrow(df)
  df <- dplyr::filter(df, .data$qc_pass %in% TRUE)
  n_after  <- nrow(df)
  if (verbose) logger::log_info("Manhattan: filtered to QC-pass => {n_after}/{n_before} rows remain.")

  if (!nrow(df) || !"p_ivw" %in% names(df)) {
    if (!"p_ivw" %in% names(df) && verbose) {
      logger::log_warn("Manhattan: required column 'p_ivw' missing after filtering; returning placeholder plot.")
    } else if (verbose) {
      logger::log_warn("Manhattan: no QC-pass rows; returning placeholder plot.")
    }
    return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no QC-pass results)"))
  }

  # group label for coloring (use whatever is available)
  group_source <- NULL
  group_col <- dplyr::coalesce(
    df$cause_level_3 %||% NULL,
    df$cause_level_2 %||% NULL,
    df$gbd_cause     %||% NULL
  )
  if (!is.null(df$cause_level_3 %||% NULL)) group_source <- "cause_level_3" else
    if (!is.null(df$cause_level_2 %||% NULL)) group_source <- "cause_level_2" else
      if (!is.null(df$gbd_cause %||% NULL))   group_source <- "gbd_cause" else
        group_source <- "constant('All')"

  if (is.null(group_col)) group_col <- factor("All")
  df$group <- as.factor(group_col)

  if (verbose) {
    logger::log_info("Manhattan: grouping by {group_source}; {length(levels(df$group))} groups.")
  }

  # positions (deterministic within group, then outcome label if present)
  df <- dplyr::arrange(df, .data$group, dplyr::coalesce(.data$outcome, ""))
  df$idx <- seq_len(nrow(df))

  # -log10 p
  df$logp <- -log10(pmax(df$p_ivw, .Machine$double.xmin))
  if (verbose) {
    rng <- range(df$logp, finite = TRUE)
    logger::log_info("Manhattan: computed -log10(p); range = [{round(rng[1], 3)}, {round(rng[2], 3)}].")
  }

  # Multiple testing handling
  if (Multiple_testing_correction == "BH") {
    df$q_bh <- stats::p.adjust(df$p_ivw, method = "BH")
    df$sig  <- df$q_bh < alpha
    thr_y   <- NA_real_
    if (verbose) {
      n_sig <- sum(df$sig %in% TRUE, na.rm = TRUE)
      logger::log_info("Manhattan: BH adjustment done; {n_sig} hits at q < {alpha}. No single threshold line drawn.")
    }
  } else {
    m       <- nrow(df)
    alpha_b <- alpha / max(1L, m)
    df$sig  <- df$p_ivw < alpha_b
    thr_y   <- -log10(alpha_b)
    if (verbose) {
      n_sig <- sum(df$sig %in% TRUE, na.rm = TRUE)
      logger::log_info("Manhattan: Bonferroni with m={m}; alpha_b={signif(alpha_b, 3)}; threshold y={round(thr_y, 3)}; {n_sig} hits.")
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$logp, colour = .data$group)) +
    ggplot2::geom_point(alpha = 0.9) +
    { if (!is.na(thr_y)) ggplot2::geom_hline(yintercept = thr_y, linetype = "dashed") } +
    ggplot2::labs(
      x = "Outcome (grouped)", y = expression(-log[10](p[IVW])),
      colour = "Cause"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )

  if (verbose) logger::log_info("Manhattan: plot object constructed; returning ggplot.")
  p
}

#' Volcano plot over outcomes
#' @param results_df tidy results (needs beta_ivw, se_ivw, p_ivw; ideally qc_pass, nsnp_after)
#' @param Multiple_testing_correction "BH" or "bonferroni"
#' @param alpha significance level for corrections (default 0.05)
#' @param label_top_n integer; label top N significant hits if ggrepel available
#' @export
volcano_plot <- function(results_df,
                         Multiple_testing_correction = c("BH","bonferroni"),
                         alpha = 0.05,
                         label_top_n = 20) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)

  if (is.null(results_df) || !nrow(results_df)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "Volcano (no results)"))
  }
  df <- tibble::as_tibble(results_df)
  need <- c("beta_ivw","se_ivw","p_ivw")
  if (!all(need %in% names(df))) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "Volcano (missing beta/se/p)"))
  }

  # computed axes
  df$z_ivw   <- df$beta_ivw / df$se_ivw
  df$logp    <- -log10(pmax(df$p_ivw, .Machine$double.xmin))
  if (!"qc_pass" %in% names(df)) df$qc_pass <- TRUE
  if (!"nsnp_after" %in% names(df)) df$nsnp_after <- NA_integer_

  # group (colour)
  group_col <- dplyr::coalesce(
    df$cause_level_3 %||% NULL,
    df$cause_level_2 %||% NULL,
    df$gbd_cause     %||% NULL
  )
  if (is.null(group_col)) group_col <- factor("All")
  df$group <- as.factor(group_col)

  # flags for outline ring
  het_fail   <- (!is.na(df$Q_p_ivw) & df$Q_p_ivw < 0.05) | (!is.na(df$I2_ivw) & df$I2_ivw > 50)
  egger_fail <- !is.na(df$egger_intercept_p) & df$egger_intercept_p < 0.05
  df$flag_any <- het_fail | egger_fail

  # multiple-testing significance
  if (Multiple_testing_correction == "BH") {
    df$q_bh <- stats::p.adjust(df$p_ivw, method = "BH")
    df$sig  <- df$q_bh < alpha
    thr_y   <- NA_real_
  } else {
    m       <- sum(!is.na(df$p_ivw))
    alpha_b <- alpha / max(1L, m)
    df$sig  <- df$p_ivw < alpha_b
    thr_y   <- -log10(alpha_b)
  }

  # base plot: shape 21 so we can have outline ring for flagged points
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$z_ivw, y = .data$logp,
                 fill = .data$group,
                 alpha = .data$qc_pass,
                 size = .data$nsnp_after)
  ) +
    ggplot2::geom_point(shape = 21, colour = NA) +
    # add outline for flagged
    ggplot2::geom_point(
      data = dplyr::filter(df, .data$flag_any %in% TRUE),
      ggplot2::aes(x = .data$z_ivw, y = .data$logp),
      inherit.aes = FALSE, shape = 21, fill = NA, colour = "black", stroke = 0.8
    ) +
    { if (!is.na(thr_y)) ggplot2::geom_hline(yintercept = thr_y, linetype = "dashed") } +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.5, 4), guide = ggplot2::guide_legend(title = "NSNP")) +
    ggplot2::labs(
      x = expression(Z == beta[IVW] / SE[IVW]),
      y = expression(-log[10](p[IVW])),
      fill = "Cause",
      title = "Volcano of MR IVW results"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  # Label top N significant hits if ggrepel present
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    lab_df <- dplyr::mutate(df, label = dplyr::coalesce(.data$outcome, NA_character_))
    lab_df <- dplyr::filter(lab_df, .data$sig %in% TRUE, !is.na(.data$label))
    if (nrow(lab_df)) {
      lab_df <- lab_df %>%
        dplyr::arrange(dplyr::desc(abs(.data$z_ivw))) %>%
        dplyr::slice_head(n = min(label_top_n, n()))
      p <- p + ggrepel::geom_text_repel(
        data = lab_df,
        ggplot2::aes(x = .data$z_ivw, y = .data$logp, label = .data$label),
        max.overlaps = 100, box.padding = 0.3, point.padding = 0.2, size = 3
      )
    }
  }

  p
}

# helper for coalescing with NULLs (base R compatible)
`%||%` <- function(a, b) if (!is.null(a)) a else b
