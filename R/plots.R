#' Manhattan plot over outcomes (IVW p-values on QC-pass outcomes)
#' L2-only brackets; right-tip anchored labels; x-axis numbers hidden
#' Colors: grey = protective (β<0), black = risk (β≥0)
#' Legends:
#'   - Always: Effect direction (grey/black)
#'   - BH:     Significance (red ring = BH significant) and colour legend title becomes "FDR < 0.05 significant results:"
#'   - Bonf.:  Significance (dashed line = Bonferroni threshold)
#'
#' @param results_df tidy results (must contain results_p_ivw; ideally results_qc_pass)
#' @param Multiple_testing_correction "BH" or "bonferroni"
#' @param qc_pass_only logical; if TRUE, plot only QC-pass rows
#' @param alpha significance level for corrections (default 0.05)
#' @param verbose if TRUE, emit step-by-step logs via {logger}
#' @param left_nuge numeric; horizontal leftward nudge of L2 labels in x data units (default 1)
#' @param tree_down numeric; vertical nudge (down if negative) for L2 bracket + labels (default -0.08)
#' @param label_down numeric; extra vertical nudge (down if negative) for labels only (default -0.02)
#' @export
manhattan_plot <- function(results_df,
                           Multiple_testing_correction = c("BH","bonferroni"),
                           qc_pass_only = TRUE,
                           alpha = 0.05,
                           verbose = TRUE,
                           left_nuge = 1,
                           tree_down = -0.08,
                           label_down = -0.02,
                           exposure = NULL){
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)

  if (is.null(results_df) || !nrow(results_df)) {
    if (verbose) logger::log_warn("Manhattan: input results_df is NULL or has 0 rows; returning placeholder plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no results)"))
  }

  if (verbose) logger::log_info("Manhattan: starting with {nrow(results_df)} rows; MTC = {Multiple_testing_correction}, alpha = {alpha}")

  df <- tibble::as_tibble(results_df)

  # Alias AFTER df exists
  if (!"gbd_cause" %in% names(df) && "Cause Name" %in% names(df)) {
    df$gbd_cause <- df[["Cause Name"]]
  }

  # ---------------- ARD-only filter ----------------
  #if (!"ARD_selected" %in% names(df)) {
  #  if (verbose) logger::log_warn("Manhattan: 'ARD_selected' column missing; cannot restrict to ARDs.")
  #  return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (ARD_selected column missing)"))
  #}
  #if (!is.logical(df$ARD_selected)) {
  #  df$ARD_selected <- df$ARD_selected %in% c(TRUE, "TRUE", "True", "true", 1, "1", "T")
  #}
  #n_before_ard <- nrow(df)
  #df <- dplyr::filter(df, .data$ARD_selected %in% TRUE)
  #if (verbose) logger::log_info("Manhattan: restricted to ARD-selected => {nrow(df)}/{n_before_ard} rows.")
  #if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no ARD-selected rows)"))
  # -------------------------------------------------

  #----------helpers---------
  # human-friendly labels for the compare modes
  .pretty_compare <- function(x) {
    rec <- c(
      "cause_vs_rest_all"              = "Within-cause effect size vs all-cause effect size",
      "ARD_vs_nonARD_within_cause"     = "ARD vs non-ARD (within cause)",
      "ARD_in_cause_vs_ARD_elsewhere"  = "ARD in cause vs ARD elsewhere"
    )
    out <- rec[as.character(x)]
    ifelse(is.na(out), as.character(x), out)
  }

  # cause_level_*  -> integer level number
  .level_number <- function(level) {
    as.integer(gsub("\\D", "", as.character(level)))
  }

  # Bonferroni two-sided z threshold for m tested causes at alpha
  .z_thr_bonf <- function(m, alpha = 0.05) {
    if (!is.finite(m) || m < 1) return(Inf)
    stats::qnorm(1 - alpha/(2*m))
  }
  #--------


  # ---- QC-pass (optional filter) ----
  if (!"results_qc_pass" %in% names(df)) {
    if (verbose) logger::log_info("Manhattan: 'results_qc_pass' not present; assuming all pass.")
    df$results_qc_pass <- TRUE
  } else if (!is.logical(df$results_qc_pass)) {
    df$results_qc_pass <- df$results_qc_pass %in% c(TRUE,"TRUE","True","true",1,"1","T")
  }

  n_before <- nrow(df)
  if (isTRUE(qc_pass_only)) {
    df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
    if (verbose) logger::log_info("Manhattan: filtered to QC-pass => {nrow(df)}/{n_before} rows remain.")
  } else {
    if (verbose) logger::log_info("Manhattan: qc_pass_only=FALSE; keeping all {n_before} rows (non-QC will be shown).")
  }
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no rows after QC filter)"))

  # Required columns
  if (!"results_p_ivw" %in% names(df)) {
    if (verbose) logger::log_warn("Manhattan: required column 'results_p_ivw' missing; returning placeholder plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (missing results_p_ivw)"))
  }
  if (!"results_beta_ivw" %in% names(df)) {
    stop("manhattan_plot(): requires 'results_beta_ivw' to colour points by effect direction.", call. = FALSE)
  }

  # Compute -log10(p)
  df <- dplyr::filter(df, !is.na(.data$results_p_ivw))
  df$logp <- dplyr::if_else(
    !is.na(df$results_log10p_ivw),
    -df$results_log10p_ivw,
    -log10(pmax(df$results_p_ivw, .Machine$double.xmin))
  )
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no finite p-values)"))

  # =========== GROUPING & X POSITIONS (L1→L2) ===========
  df <- dplyr::arrange(
    df,
    dplyr::coalesce(.data$cause_level_1, "~"),
    dplyr::coalesce(.data$cause_level_2, "~"),
    dplyr::coalesce(.data$results_outcome, "")
  )
  df$idx <- seq_len(nrow(df))

  l2_map <- df |>
    dplyr::group_by(.data$cause_level_1, .data$cause_level_2) |>
    dplyr::summarise(start = min(.data$idx), end = max(.data$idx), .groups = "drop") |>
    dplyr::arrange(.data$cause_level_1, .data$start) |>
    dplyr::mutate(center = (start + end)/2, alt_id = dplyr::row_number())

  # ========== MULTIPLE TESTING ==========
  if (Multiple_testing_correction == "BH") {
    df$q_bh <- stats::p.adjust(df$results_p_ivw, method = "BH")
    df$sig  <- df$q_bh < alpha
    thr_y   <- NA_real_   # BH has no single flat threshold line
  } else {
    m       <- nrow(df)
    alpha_b <- alpha / max(1L, m)
    df$sig  <- df$results_p_ivw < alpha_b
    thr_y   <- -log10(alpha_b)  # Bonferroni flat line
  }

  # ========== EFFECT DIRECTION COLOURING ==========
  df$effect_dir <- ifelse(is.finite(df$results_beta_ivw) & df$results_beta_ivw < 0,
                          "Protective (β<0)", "Risk (β≥0)")
  df$effect_dir <- factor(df$effect_dir, levels = c("Protective (β<0)", "Risk (β≥0)"))

  # ========== L2 BRACKETS ONLY ==========
  y_max   <- max(df$logp, na.rm = TRUE)
  base_y  <- 0.02 + tree_down
  cap     <- 0.15
  lab_gap <- 0.08

  seg_l2 <- l2_map |> dplyr::transmute(x = start, xend = end, y = base_y, yend = base_y)
  teeth_l2 <- dplyr::bind_rows(
    l2_map |> dplyr::transmute(x = start, xend = start, y = base_y, yend = base_y - cap/2),
    l2_map |> dplyr::transmute(x = end,   xend = end,   y = base_y, yend = base_y - cap/2)
  )
  lab_l2 <- l2_map |> dplyr::transmute(
    x = center - left_nuge,
    y = base_y - lab_gap/2 + label_down,
    label = dplyr::coalesce(cause_level_2, "NA")
  )

  # -------- legend title: BH gets custom text ----------
  legend_title <- if (Multiple_testing_correction == "BH") {
    "FDR < 0.05 significant results:"
  } else {
    "Effect direction:"
  }

  # --- title uses the explicit `exposure` argument ---
  exposure_lab <- tryCatch(as.character(exposure)[1], error = function(e) NA_character_)
  if (is.null(exposure_lab) || is.na(exposure_lab) || !nzchar(exposure_lab)) exposure_lab <- "exposure"
  plot_title <- sprintf("ARD-wide MR IVW of %s by GBD Cause", exposure_lab)

  # --- END: dynamic title ---

  # ========== BASE PLOT (effect direction legend) ==========
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$logp)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$effect_dir), alpha = 0.9, size = 1.6) +
    ggplot2::scale_colour_manual(
      name   = legend_title,
      values = c("Protective (β<0)" = "grey60", "Risk (β≥0)" = "black")
    ) +
    ggplot2::labs(x = NULL, y = expression(-log[10](p[IVW])), title = plot_title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(t = 10, r = 10, b = 80, l = 10),
      legend.position    = "top",
      legend.title       = ggplot2::element_text(size = 9),
      legend.text        = ggplot2::element_text(size = 9)
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1))

  # --- Significance legend (BH: red ring) ---
  if (Multiple_testing_correction == "BH") {
    p <- p + ggplot2::geom_point(
      data = dplyr::filter(df, .data$sig %in% TRUE),
      mapping = ggplot2::aes(x = .data$idx, y = .data$logp, shape = "BH significant"),
      inherit.aes = FALSE,
      shape = 21, fill = NA, colour = "firebrick",
      stroke = 1.4, size = 2.8, show.legend = TRUE
    ) +
      ggplot2::scale_shape_manual(
        name   = "Significance:",
        values = c("BH significant" = 21),
        guide  = ggplot2::guide_legend(
          order = 2,
          override.aes = list(colour = "firebrick", fill = NA, size = 3.2, stroke = 1.6)
        )
      )
  } else {
    # --- Significance legend (Bonferroni: dashed line) ---
    if (!is.na(thr_y)) {
      p <- p +
        ggplot2::geom_hline(
          data = data.frame(yint = thr_y),
          ggplot2::aes(yintercept = yint, linetype = "Bonferroni threshold"),
          linewidth = 0.4
        ) +
        ggplot2::scale_linetype_manual(
          name   = "Significance:",
          values = c("Bonferroni threshold" = "dashed"),
          guide  = ggplot2::guide_legend(order = 2)
        )
    }
  }

  # --- label ALL significant points with phenotype name ---
  if (!"results_outcome" %in% names(df)) {
    stop("manhattan_plot(): 'results_outcome' column is required to label significant points.", call. = FALSE)
  }
  sig_df <- dplyr::filter(
    df, .data$sig %in% TRUE, !is.na(.data$results_outcome),
    is.finite(.data$idx), is.finite(.data$logp)
  )
  if (nrow(sig_df)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = sig_df,
        ggplot2::aes(x = .data$idx, y = .data$logp, label = .data$results_outcome),
        max.overlaps = Inf, box.padding = 0.25, point.padding = 0.15,
        min.segment.length = 0, size = 2.7, seed = 123,
        segment.size = 0.2
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = sig_df,
        ggplot2::aes(x = .data$idx, y = .data$logp, label = .data$results_outcome),
        size = 2.7, vjust = -0.25, check_overlap = TRUE
      )
    }
  }

  # draw the L2 bracket and its labels last (near the bottom)
  p <- p +
    ggplot2::geom_segment(data = seg_l2,   ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    ggplot2::geom_segment(data = teeth_l2, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    ggplot2::geom_text(
      data = lab_l2,
      ggplot2::aes(x = x, y = y, label = label),
      angle = 45, hjust = 1, vjust = 1.15, size = 3
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.1, y_max * 1.05), clip = "off")

  if (verbose) logger::log_info("Manhattan: plot object constructed; returning ggplot.")
  p
}

# helper for coalescing with NULLs (base R compatible)
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Directional enrichment forest (protective vs risk)
#'
#' Draws a diverging "forest/dumbbell" plot from the output of
#' `enrichment_by_cause_directional()` or a filtered slice of
#' `run_enrichment(...)$by_cause_tbl`.
#'
#' @param by_cause_tbl Tibble with columns: level, cause, compare_mode,
#'        SES_prot, SES_risk, q_prot, q_risk, n_pos, n_neg, n_total, exposure.
#' @param level One of "cause_level_1","cause_level_2","cause_level_3".
#' @param compare_mode One of "ARD_vs_nonARD_within_cause",
#'        "cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere".
#' @param title Optional plot title; if NULL a default is constructed.
#' @param subtitle Optional subtitle; if NULL a default is constructed.
#' @param y_lab Y axis label (default "Cause").
#' @param ring_q FDR threshold for red ring highlight (default 0.05).
#' @return A ggplot object.
#' @export
plot_enrichment_directional_forest <- function(
    by_cause_tbl,
    level,
    compare_mode,
    title = NULL,
    subtitle = NULL,
    y_lab = "Cause",
    ring_q = 0.05
) {
  stopifnot(all(c("level","cause","compare_mode","SES_prot","SES_risk","q_prot","q_risk") %in% names(by_cause_tbl)))

  df <- by_cause_tbl[by_cause_tbl$level == level & by_cause_tbl$compare_mode == compare_mode, , drop = FALSE]
  if (nrow(df) == 0) {
    warning("No rows to plot for level=", level, " & compare_mode=", compare_mode)
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  # order by max |SES|
  df$._ord <- pmax(abs(df$SES_prot), abs(df$SES_risk), na.rm = TRUE)
  df <- df[order(df$._ord, decreasing = FALSE), , drop = FALSE]
  df$cause_f <- factor(df$cause, levels = df$cause)

  seg <- data.frame(cause_f = df$cause_f, xmin = df$SES_prot, xmax = df$SES_risk)

  # ===== titles / subtitles =====
  exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NULL
  if (is.null(title)) {
    title <- sprintf("Cause Level %d Enrichment Analysis for %s",
                     .level_number(level),
                     ifelse(is.null(exp_label) || !nzchar(exp_label), "exposure", exp_label))
  }
  if (is.null(subtitle)) {
    subtitle <- .pretty_compare(unique(as.character(df$compare_mode)))
  }

  # ===== bubble sizes (+2 pt), no black ring =====
  if ("n_pos" %in% names(df) && any(is.finite(df$n_pos))) {
    rng <- range(df$n_pos, na.rm = TRUE)
    pts <- if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
      scales::rescale(df$n_pos, to = c(3,6), from = rng)
    } else rep(4.5, nrow(df))
  } else {
    pts <- rep(4.5, nrow(df))
  }
  pts <- pmax(2.5, pts + 2)  # +2 points everywhere

  # ===== Bonferroni ±z* lines (two-sided) computed per-plot =====
  m      <- nrow(df)
  z_thr  <- .z_thr_bonf(m, alpha = 0.05)

  ggplot2::ggplot() +
    # significance guide: ±z* (Bonferroni)
    ggplot2::geom_vline(xintercept = c(-z_thr, z_thr),
                        ggplot2::aes(linetype = "Bonferroni (α=0.05, 2-sided)"),
                        linewidth = 0.5, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    ggplot2::scale_linetype_manual(
      name   = "Significance:",
      values = c("Bonferroni (α=0.05, 2-sided)" = "dashed"),
      guide  = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xmin, xend = xmax, y = cause_f, yend = cause_f),
                          linewidth = 0.6, alpha = 0.5) +
    # protective dots (no ring)
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = SES_prot, y = cause_f),
      shape = 16, size = pts, colour = "grey60"
    ) +
    # risk dots (no ring)
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = SES_risk, y = cause_f),
      shape = 16, size = pts, colour = "black"
    ) +
    ggplot2::scale_x_continuous("SES z-score (− protective  ←  0  →  risk +)",
                                expand = ggplot2::expansion(mult = c(0.05,0.1))) +
    ggplot2::ylab(y_lab) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 11),
      legend.text  = ggplot2::element_text(size = 11)
    )
}

#' Cause-level signed enrichment forest (single-dot)
#'
#' Draws a one-axis forest using the signed SES from enrichment.
#' Negative values = protective; positive = risk.
#'
#' Expects by_cause_tbl to contain at least:
#'   level, cause, compare_mode, SES_signed, q_signed
#' If you used different names (e.g., SES_comb/q_comb), rename below.
#'
#' @param by_cause_tbl Tibble with signed enrichment columns.
#' @param level One of "cause_level_1","cause_level_2","cause_level_3".
#' @param compare_mode One of "ARD_vs_nonARD_within_cause",
#'        "cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere".
#' @param title,subtitle Optional strings.
#' @param y_lab Y axis label.
#' @param ring_q FDR cutoff; red ring when q_signed < ring_q.
#' @export
plot_enrichment_signed_forest <- function(
    by_cause_tbl,
    level,
    compare_mode,
    title = NULL,
    subtitle = NULL,
    y_lab = "Cause",
    ring_q = 0.05
) {
  need <- c("level","cause","compare_mode","SES_signed","q_signed")
  stopifnot(all(need %in% names(by_cause_tbl)))

  df <- by_cause_tbl[by_cause_tbl$level == level & by_cause_tbl$compare_mode == compare_mode, , drop = FALSE]
  if (nrow(df) == 0) {
    warning("No rows to plot for level=", level, " & compare_mode=", compare_mode)
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  # order by |SES_signed|
  df$._ord <- abs(df$SES_signed)
  df <- df[order(df$._ord, decreasing = FALSE), , drop = FALSE]
  df$cause_f <- factor(df$cause, levels = df$cause)

  # titles
  exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NULL
  if (is.null(title)) {
    title <- sprintf("Cause Level %d Enrichment Analysis for %s",
                     .level_number(level),
                     ifelse(is.null(exp_label) || !nzchar(exp_label), "exposure", exp_label))
  }
  if (is.null(subtitle)) {
    subtitle <- .pretty_compare(unique(as.character(df$compare_mode)))
  }

  # sizes (+2), no ring
  if ("n_pos" %in% names(df) && any(is.finite(df$n_pos))) {
    rng <- range(df$n_pos, na.rm = TRUE)
    pts <- if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
      scales::rescale(df$n_pos, to = c(3,6), from = rng)
    } else rep(4.5, nrow(df))
  } else {
    pts <- rep(4.5, nrow(df))
  }
  pts <- pmax(2.5, pts + 2)

  # direction fill (kept for legend), but no stroke ring
  df$dir <- ifelse(df$SES_signed < 0, "Protective", "Risk")

  # per-plot Bonferroni ±z*
  m     <- nrow(df)
  z_thr <- .z_thr_bonf(m, alpha = 0.05)

  seg <- data.frame(cause_f = df$cause_f, xmin = 0, xmax = df$SES_signed)

  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = c(-z_thr, z_thr),
                        ggplot2::aes(linetype = "Bonferroni (α=0.05, 2-sided)"),
                        linewidth = 0.5, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    ggplot2::scale_linetype_manual(
      name   = "Significance:",
      values = c("Bonferroni (α=0.05, 2-sided)" = "dashed"),
      guide  = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xmin, xend = xmax, y = cause_f, yend = cause_f),
                          linewidth = 0.6, alpha = 0.5) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = SES_signed, y = cause_f, colour = dir),
      shape = 16, size = pts
    ) +
    ggplot2::scale_colour_manual(
      name   = "Direction:",
      values = c("Protective" = "grey60", "Risk" = "black")
    ) +
    ggplot2::scale_x_continuous("SES z-score (− protective  ←  0  →  risk +)",
                                expand = ggplot2::expansion(mult = c(0.05,0.1))) +
    ggplot2::ylab(y_lab) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 11),
      legend.text  = ggplot2::element_text(size = 11)
    )
}


#' Global ARD vs non-ARD directional summary (protective vs risk)
#'
#' One-row dumbbell-like plot for the global enrichment.
#'
#' @param global_tbl Output of `enrichment_global_directional()`.
#' @param title Optional title.
#' @param ring_q FDR threshold for red ring.
#' @return ggplot.
#' @export
plot_enrichment_global <- function(global_tbl, title = NULL, ring_q = 0.05) {
  stopifnot(all(c("SES_prot","SES_risk","q_prot","q_risk") %in% names(global_tbl)))
  df <- global_tbl[1, , drop = FALSE]
  df$y <- factor("ARD vs non-ARD", levels = "ARD vs non-ARD")
  seg <- data.frame(y = df$y, xmin = df$SES_prot, xmax = df$SES_risk)

  exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NULL
  if (is.null(title)) title <- sprintf("Global Enrichment Analysis for %s",
                                       ifelse(is.null(exp_label) || !nzchar(exp_label), "exposure", exp_label))

  z_thr <- .z_thr_bonf(m = 1, alpha = 0.05)

  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = c(-z_thr, z_thr),
                        ggplot2::aes(linetype = "Bonferroni (α=0.05, 2-sided)"),
                        linewidth = 0.5, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    ggplot2::scale_linetype_manual(
      name   = "Significance:",
      values = c("Bonferroni (α=0.05, 2-sided)" = "dashed"),
      guide  = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y),
                          linewidth = 0.6, alpha = 0.5) +
    ggplot2::geom_point(data = df, ggplot2::aes(x = SES_prot, y = y),
                        shape = 16, size = 7, colour = "grey60") +
    ggplot2::geom_point(data = df, ggplot2::aes(x = SES_risk, y = y),
                        shape = 16, size = 7, colour = "black") +
    ggplot2::scale_x_continuous("SES z-score (− protective  ←  0  →  risk +)",
                                expand = ggplot2::expansion(mult = c(0.05,0.1))) +
    ggplot2::labs(title = title, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 11),
      legend.text  = ggplot2::element_text(size = 11)
    )
}

#' Global ARD vs non-ARD signed summary (single-dot)
#'
#' One-row forest for the global signed SES.
#'
#' Expects global_tbl to contain at least: SES_signed, q_signed.
#' If you used different names (e.g., SES_comb/q_comb), rename below.
#'
#' @param global_tbl Output of enrichment_global_* with signed columns.
#' @param title Optional title.
#' @param ring_q FDR threshold for red ring.
#' @export
plot_enrichment_global_signed <- function(global_tbl, title = NULL, ring_q = 0.05) {
  need <- c("SES_signed","q_signed")
  stopifnot(all(need %in% names(global_tbl)))

  df <- global_tbl[1, , drop = FALSE]
  df$y <- factor("ARD vs non-ARD", levels = "ARD vs non-ARD")
  df$dir <- ifelse(df$SES_signed < 0, "Protective", "Risk")
  seg <- data.frame(y = df$y, xmin = 0, xmax = df$SES_signed)

  exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NULL
  if (is.null(title)) title <- sprintf("Global Enrichment Analysis for %s",
                                       ifelse(is.null(exp_label) || !nzchar(exp_label), "exposure", exp_label))

  z_thr <- .z_thr_bonf(m = 1, alpha = 0.05)

  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = c(-z_thr, z_thr),
                        ggplot2::aes(linetype = "Bonferroni (α=0.05, 2-sided)"),
                        linewidth = 0.5, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    ggplot2::scale_linetype_manual(
      name   = "Significance:",
      values = c("Bonferroni (α=0.05, 2-sided)" = "dashed"),
      guide  = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y),
                          linewidth = 0.6, alpha = 0.5) +
    ggplot2::geom_point(
      data = df, ggplot2::aes(x = SES_signed, y = y, colour = dir),
      shape = 16, size = 7
    ) +
    ggplot2::scale_colour_manual(
      name   = "Direction:",
      values = c("Protective" = "grey60", "Risk" = "black")
    ) +
    ggplot2::scale_x_continuous("SES z-score (− protective  ←  0  →  risk +)",
                                expand = ggplot2::expansion(mult = c(0.05,0.1))) +
    ggplot2::labs(title = title, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 11),
      legend.text  = ggplot2::element_text(size = 11)
    )
}

#' Volcano plot (IVW): effect size vs significance
#'
#' Colours: grey = protective (β<0), black = risk (β≥0).
#' Significance legend:
#'   - BH: red ring on significant points (q < alpha).
#'   - Bonferroni: dashed horizontal line at alpha/m.
#'
#' Uses results_log10p_ivw if available (expected to be log10(p), i.e. negative),
#' otherwise falls back to -log10(results_p_ivw) with underflow protection.
#'
#' @param results_df Tidy MR results (needs results_beta_ivw and results_p_ivw or results_log10p_ivw).
#' @param Multiple_testing_correction "BH" or "bonferroni".
#' @param qc_pass_only If TRUE (default), drop non-QC rows if results_qc_pass exists.
#' @param alpha Significance level (default 0.05).
#' @param exposure Optional exposure label for the title.
#' @param label_top_n How many most-significant points to label (default 25). Use Inf to label all significant.
#' @param verbose If TRUE, log a few steps via {logger} if available.
#' @export
volcano_plot <- function(results_df,
                         Multiple_testing_correction = c("BH","bonferroni"),
                         qc_pass_only = TRUE,
                         alpha = 0.05,
                         exposure = NULL,
                         label_top_n = 25,
                         verbose = TRUE) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)

  if (is.null(results_df) || !nrow(results_df)) {
    if (verbose && requireNamespace("logger", quietly = TRUE)) logger::log_warn("Volcano: no rows.")
    return(ggplot2::ggplot() + ggplot2::labs(title = "Volcano (no results)"))
  }

  df <- tibble::as_tibble(results_df)

  # Optional QC filter
  if ("results_qc_pass" %in% names(df) && isTRUE(qc_pass_only)) {
    if (!is.logical(df$results_qc_pass)) {
      df$results_qc_pass <- df$results_qc_pass %in% c(TRUE,"TRUE","True","true",1,"1","T")
    }
    df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  }

  # Required columns
  if (!all(c("results_beta_ivw","results_p_ivw") %in% names(df)) &&
      !"results_log10p_ivw" %in% names(df)) {
    stop("volcano_plot(): need 'results_beta_ivw' and either 'results_p_ivw' or 'results_log10p_ivw'.", call. = FALSE)
  }

  # y = -log10(p) using your negative log10(p) column if present
  df <- dplyr::mutate(
    df,
    logp = dplyr::if_else(
      !is.na(.data$results_log10p_ivw),
      -.data$results_log10p_ivw,                                # results_log10p_ivw is log10(p) (negative) -> negate
      -log10(pmax(.data$results_p_ivw, .Machine$double.xmin))   # fallback
    )
  )

  # Keep finite x/y
  df <- dplyr::filter(df, is.finite(.data$results_beta_ivw), is.finite(.data$logp))
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "Volcano (no finite rows)"))
  }

  # Multiple testing
  if (Multiple_testing_correction == "BH") {
    df$q_bh <- stats::p.adjust(df$results_p_ivw, method = "BH")
    df$sig  <- df$q_bh < alpha
    thr_y   <- NA_real_
  } else {
    m       <- nrow(df)
    alpha_b <- alpha / max(1L, m)
    df$sig  <- df$results_p_ivw < alpha_b
    thr_y   <- -log10(alpha_b)
  }

  # Effect direction
  df$effect_dir <- ifelse(df$results_beta_ivw < 0, "Protective (β<0)", "Risk (β≥0)")
  df$effect_dir <- factor(df$effect_dir, levels = c("Protective (β<0)","Risk (β≥0)"))

  # Title
  exposure_lab <- tryCatch(as.character(exposure)[1], error = function(e) NA_character_)
  if (is.null(exposure_lab) || is.na(exposure_lab) || !nzchar(exposure_lab)) exposure_lab <- "exposure"
  plot_title <- sprintf("Phenome-wide MR IVW of %s — volcano", exposure_lab)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$results_beta_ivw, y = .data$logp)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$effect_dir), alpha = 0.9, size = 1.7) +
    ggplot2::scale_colour_manual(
      name   = "Effect direction:",
      values = c("Protective (β<0)" = "grey60", "Risk (β≥0)" = "black")
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    ggplot2::labs(
      x = expression(beta[IVW]),
      y = expression(-log[10](p[IVW])),
      title = plot_title
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.title    = ggplot2::element_text(size = 9),
      legend.text     = ggplot2::element_text(size = 9)
    )

  # Significance layer / legend
  if (Multiple_testing_correction == "BH") {
    p <- p +
      ggplot2::geom_point(
        data = dplyr::filter(df, .data$sig %in% TRUE),
        mapping = ggplot2::aes(x = .data$results_beta_ivw, y = .data$logp, shape = "BH significant"),
        inherit.aes = FALSE,
        shape = 21, fill = NA, colour = "firebrick",
        stroke = 1.2, size = 2.8, show.legend = TRUE
      ) +
      ggplot2::scale_shape_manual(
        name   = "Significance:",
        values = c("BH significant" = 21),
        guide  = ggplot2::guide_legend(
          order = 2,
          override.aes = list(colour = "firebrick", fill = NA, size = 3.2, stroke = 1.6)
        )
      )
  } else {
    if (!is.na(thr_y)) {
      p <- p +
        ggplot2::geom_hline(
          yintercept = thr_y,
          ggplot2::aes(linetype = "Bonferroni threshold"),
          linewidth = 0.4
        ) +
        ggplot2::scale_linetype_manual(
          name   = "Significance:",
          values = c("Bonferroni threshold" = "dashed"),
          guide  = ggplot2::guide_legend(order = 2)
        )
    }
  }

  # Labels for most significant points (if results_outcome exists)
  if ("results_outcome" %in% names(df)) {
    lab_df <- dplyr::filter(df, .data$sig %in% TRUE & !is.na(.data$results_outcome))
    if (nrow(lab_df)) {
      # rank by significance (higher -log10 p first)
      lab_df <- lab_df[order(lab_df$logp, decreasing = TRUE), , drop = FALSE]
      if (is.finite(label_top_n)) lab_df <- utils::head(lab_df, label_top_n)

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = lab_df,
          ggplot2::aes(x = .data$results_beta_ivw, y = .data$logp, label = .data$results_outcome),
          max.overlaps = Inf, box.padding = 0.25, point.padding = 0.15,
          min.segment.length = 0, size = 2.7, seed = 123, segment.size = 0.2
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = lab_df,
          ggplot2::aes(x = .data$results_beta_ivw, y = .data$logp, label = .data$results_outcome),
          size = 2.7, vjust = -0.25, check_overlap = TRUE
        )
      }
    }
  }

  p
}

