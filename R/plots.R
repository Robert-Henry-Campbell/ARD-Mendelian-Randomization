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
#' @param dot_names logical; if `TRUE`, annotate significant points with their
#'   outcome names (default `TRUE`).
#' @export
manhattan_plot <- function(results_df,
                           Multiple_testing_correction = c("BH","bonferroni"),
                           qc_pass_only = TRUE,
                           alpha = 0.05,
                           verbose = TRUE,
                           left_nuge = 1,
                           tree_down = -0.08,
                           label_down = -0.02,
                           dot_names = TRUE,
                           exposure = NULL){
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  dot_names <- isTRUE(dot_names)

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
  if (dot_names && nrow(sig_df)) {
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
  p <- .ardmr_attach_plot_data(
    p,
    main = df,
    significant = sig_df,
    l2_map = l2_map,
    segments = seg_l2,
    teeth = teeth_l2,
    labels = lab_l2
  )
  p
}

#' Manhattan plot with recoloured significance highlighting
#'
#' Generates a Manhattan-style plot over outcomes, ordering phenotypes within
#' each cause-level-2 block from smallest to largest -log10(p) and colouring
#' points by significance and effect direction. Non-significant points are shown
#' in dark grey, significant protective associations (β < 0) in blue, and
#' significant risk associations (β ≥ 0) in red.
#'
#' @inheritParams manhattan_plot
#' @export
manhattan_plot_recolor <- function(results_df,
                                   Multiple_testing_correction = c("BH","bonferroni"),
                                   qc_pass_only = TRUE,
                                   alpha = 0.05,
                                   verbose = TRUE,
                                   left_nuge = 1,
                                   tree_down = -0.08,
                                   label_down = -0.02,
                                   exposure = NULL) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)

  if (is.null(results_df) || !nrow(results_df)) {
    if (verbose) logger::log_warn("Manhattan (recolor): input results_df is NULL or has 0 rows; returning placeholder plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (recolor, no results)"))
  }

  if (verbose) {
    logger::log_info(
      "Manhattan (recolor): starting with {nrow(results_df)} rows; MTC = {Multiple_testing_correction}, alpha = {alpha}"
    )
  }

  df <- tibble::as_tibble(results_df)

  if (!"gbd_cause" %in% names(df) && "Cause Name" %in% names(df)) {
    df$gbd_cause <- df[["Cause Name"]]
  }

  if (!"results_qc_pass" %in% names(df)) {
    if (verbose) logger::log_info("Manhattan (recolor): 'results_qc_pass' not present; assuming all pass.")
    df$results_qc_pass <- TRUE
  } else if (!is.logical(df$results_qc_pass)) {
    df$results_qc_pass <- df$results_qc_pass %in% c(TRUE,"TRUE","True","true",1,"1","T")
  }

  n_before <- nrow(df)
  if (isTRUE(qc_pass_only)) {
    df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
    if (verbose) logger::log_info("Manhattan (recolor): filtered to QC-pass => {nrow(df)}/{n_before} rows remain.")
  } else {
    if (verbose) logger::log_info("Manhattan (recolor): qc_pass_only=FALSE; keeping all {n_before} rows (non-QC will be shown).")
  }
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (recolor, no rows after QC filter)"))

  if (!"results_p_ivw" %in% names(df)) {
    if (verbose) logger::log_warn("Manhattan (recolor): required column 'results_p_ivw' missing; returning placeholder plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (recolor, missing results_p_ivw)"))
  }
  if (!"results_beta_ivw" %in% names(df)) {
    stop("manhattan_plot_recolor(): requires 'results_beta_ivw' to colour points by effect direction.", call. = FALSE)
  }

  df <- dplyr::filter(df, !is.na(.data$results_p_ivw))
  df$logp <- dplyr::if_else(
    !is.na(df$results_log10p_ivw),
    -df$results_log10p_ivw,
    -log10(pmax(df$results_p_ivw, .Machine$double.xmin))
  )
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (recolor, no finite p-values)"))

  df <- dplyr::arrange(
    df,
    dplyr::coalesce(.data$cause_level_1, "~"),
    dplyr::coalesce(.data$cause_level_2, "~"),
    .data$logp,
    dplyr::coalesce(.data$results_outcome, "")
  )
  df$idx <- seq_len(nrow(df))

  l2_map <- df |>
    dplyr::group_by(.data$cause_level_1, .data$cause_level_2) |>
    dplyr::summarise(start = min(.data$idx), end = max(.data$idx), .groups = "drop") |>
    dplyr::arrange(.data$cause_level_1, .data$start) |>
    dplyr::mutate(center = (start + end)/2, alt_id = dplyr::row_number())

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

  df <- df |>
    dplyr::mutate(
      colour_group = dplyr::case_when(
        .data$sig %in% TRUE & is.finite(.data$results_beta_ivw) & .data$results_beta_ivw < 0 ~
          "Significant protective (β<0)",
        .data$sig %in% TRUE & is.finite(.data$results_beta_ivw) & .data$results_beta_ivw >= 0 ~
          "Significant risk (β≥0)",
        TRUE ~ "Not significant"
      )
    )

  colour_levels <- c("Not significant", "Significant protective (β<0)", "Significant risk (β≥0)")
  present_levels <- colour_levels[colour_levels %in% unique(df$colour_group)]
  df$colour_group <- factor(df$colour_group, levels = present_levels)

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

  legend_title <- if (Multiple_testing_correction == "BH") {
    "Significance & direction (FDR < 0.05)"
  } else {
    NULL
  }

  exposure_lab <- tryCatch(as.character(exposure)[1], error = function(e) NA_character_)
  if (is.null(exposure_lab) || is.na(exposure_lab) || !nzchar(exposure_lab)) exposure_lab <- "exposure"
  plot_title <- sprintf("ARD-wide MR IVW of %s by GBD Cause", exposure_lab)

  colour_values <- c(
    "Not significant" = "#4B4B4B",
    "Significant protective (β<0)" = "#1f78b4",
    "Significant risk (β≥0)" = "#e31a1c"
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$logp)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$colour_group), alpha = 0.9, size = 1.8) +
    ggplot2::scale_colour_manual(name = legend_title, values = colour_values, drop = FALSE) +
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

  if (Multiple_testing_correction == "bonferroni" && !is.na(thr_y)) {
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

  if (!"results_outcome" %in% names(df)) {
    stop("manhattan_plot_recolor(): 'results_outcome' column is required to label significant points.", call. = FALSE)
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

  p <- p +
    ggplot2::geom_segment(data = seg_l2,   ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    ggplot2::geom_segment(data = teeth_l2, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    ggplot2::geom_text(
      data = lab_l2,
      ggplot2::aes(x = x, y = y, label = label),
      angle = 45, hjust = 1, vjust = 1.15, size = 3
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.1, y_max * 1.05), clip = "off")

  if (verbose) logger::log_info("Manhattan (recolor): plot object constructed; returning ggplot.")
  p <- .ardmr_attach_plot_data(
    p,
    main = df,
    significant = sig_df,
    l2_map = l2_map,
    segments = seg_l2,
    teeth = teeth_l2,
    labels = lab_l2
  )
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
# ========== Directional (two-dot) forest ==========
plot_enrichment_directional_forest <- function(
    by_cause_tbl, level, compare_mode,
    title = NULL, subtitle = NULL, y_lab = "Cause",
    Multiple_testing_correction = c("bonferroni","BH"),
    alpha = 0.05
) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  stopifnot(all(c("level","cause","compare_mode","SES_prot","SES_risk","q_prot","q_risk") %in% names(by_cause_tbl)))

  df <- by_cause_tbl[by_cause_tbl$level == level & by_cause_tbl$compare_mode == compare_mode, , drop = FALSE]
  if (nrow(df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())

  # order & factors
  df$._ord <- pmax(abs(df$SES_prot), abs(df$SES_risk), na.rm = TRUE)
  df       <- df[order(df$._ord, decreasing = FALSE), , drop = FALSE]
  df$cause_f <- factor(df$cause, levels = df$cause)

  # titles
  exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NULL
  if (is.null(title)) {
    title <- sprintf("Cause Level %d Enrichment Analysis for %s",
                     .level_number(level),
                     ifelse(is.null(exp_label) || !nzchar(exp_label), "exposure", exp_label))
  }
  if (is.null(subtitle)) subtitle <- .pretty_compare(unique(as.character(df$compare_mode)))

  seg <- data.frame(cause_f = df$cause_f, xmin = df$SES_prot, xmax = df$SES_risk)
  m   <- nrow(df)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xmin, xend = xmax, y = cause_f, yend = cause_f),
                          linewidth = 0.6, alpha = 0.5) +
    # base points (fixed size, no stroke)
    ggplot2::geom_point(data = df, ggplot2::aes(x = SES_prot, y = cause_f),
                        shape = 16, size = 4, colour = "grey60") +
    ggplot2::geom_point(data = df, ggplot2::aes(x = SES_risk,  y = cause_f),
                        shape = 16, size = 4, colour = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7)

  if (Multiple_testing_correction == "bonferroni") {
    zthr <- .z_thr_bonf(m, alpha = alpha)
    p <- p +
      ggplot2::geom_vline(xintercept = c(-zthr, zthr),
                          ggplot2::aes(linetype = sprintf("Bonferroni (α=%.2f, 2-sided)", alpha)),
                          linewidth = 0.5, alpha = 0.85) +
      ggplot2::scale_linetype_manual(
        name = "Significance:", values = setNames("dashed", sprintf("Bonferroni (α=%.2f, 2-sided)", alpha))
      )
  } else {
    # BH: draw firebrick rings around significant points
    ring_df_prot <- df[is.finite(df$q_prot) & df$q_prot < alpha, , drop = FALSE]
    ring_df_risk <- df[is.finite(df$q_risk) & df$q_risk < alpha, , drop = FALSE]
    p <- p +
      ggplot2::geom_point(data = ring_df_prot,
                          ggplot2::aes(x = SES_prot, y = cause_f, shape = "BH significant (q<α)"),
                          size = 5, stroke = 1.2, colour = "firebrick", fill = NA) +
      ggplot2::geom_point(data = ring_df_risk,
                          ggplot2::aes(x = SES_risk,  y = cause_f, shape = "BH significant (q<α)"),
                          size = 5, stroke = 1.2, colour = "firebrick", fill = NA) +
      ggplot2::scale_shape_manual(
        name = "Significance:", values = c("BH significant (q<α)" = 21)
      )
  }

  p <- p +
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
  p <- .ardmr_attach_plot_data(p, main = df, segments = seg)
  p
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
# ========== Signed (single-dot) forest ==========
# ========== Signed (single-dot) forest ==========
plot_enrichment_signed_forest <- function(
    by_cause_tbl, level, compare_mode,
    title = NULL, subtitle = NULL, y_lab = "Cause",
    Multiple_testing_correction = c("bonferroni","BH"),
    alpha = 0.05
) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  need <- c("level","cause","compare_mode","SES_signed","q_signed")
  stopifnot(all(need %in% names(by_cause_tbl)))

  df <- by_cause_tbl[by_cause_tbl$level == level & by_cause_tbl$compare_mode == compare_mode, , drop = FALSE]
  if (nrow(df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())

  df$._ord  <- abs(df$SES_signed)
  df        <- df[order(df$._ord, decreasing = FALSE), , drop = FALSE]
  df$cause_f <- factor(df$cause, levels = df$cause)
  df$dir     <- ifelse(df$SES_signed < 0, "Protective", "Risk")

  exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NULL
  if (is.null(title)) {
    title <- sprintf("Cause Level %d Enrichment Analysis for %s",
                     .level_number(level),
                     ifelse(is.null(exp_label) || !nzchar(exp_label), "exposure", exp_label))
  }
  if (is.null(subtitle)) subtitle <- .pretty_compare(unique(as.character(df$compare_mode)))

  seg <- data.frame(cause_f = df$cause_f, xmin = 0, xmax = df$SES_signed)
  m   <- nrow(df)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xmin, xend = xmax, y = cause_f, yend = cause_f),
                          linewidth = 0.6, alpha = 0.5) +
    ggplot2::geom_point(
      data = df, ggplot2::aes(x = SES_signed, y = cause_f, colour = dir),
      shape = 16, size = 4
    ) +
    ggplot2::scale_colour_manual(
      name = "Direction:", values = c("Protective" = "grey60", "Risk" = "black")
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7)

  if (Multiple_testing_correction == "bonferroni") {
    zthr <- .z_thr_bonf(m, alpha = alpha)
    p <- p +
      ggplot2::geom_vline(xintercept = c(-zthr, zthr),
                          ggplot2::aes(linetype = sprintf("Bonferroni (α=%.2f, 2-sided)", alpha)),
                          linewidth = 0.5, alpha = 0.85) +
      ggplot2::scale_linetype_manual(
        name = "Significance:", values = setNames("dashed", sprintf("Bonferroni (α=%.2f, 2-sided)", alpha))
      )
  } else {
    ring_df <- df[is.finite(df$q_signed) & df$q_signed < alpha, , drop = FALSE]
    p <- p +
      ggplot2::geom_point(
        data = ring_df,
        ggplot2::aes(x = SES_signed, y = cause_f, shape = "BH significant (q<α)"),
        size = 5, stroke = 1.2, colour = "firebrick", fill = NA
      ) +
      ggplot2::scale_shape_manual(
        name = "Significance:", values = c("BH significant (q<α)" = 21)
      )
  }

  p <- p +
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
  p <- .ardmr_attach_plot_data(p, main = df, segments = seg)
  p
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
# ========== Global (signed) one-dot summary ==========
# ========== Global (signed) one-dot summary ==========
plot_enrichment_global_signed <- function(
    global_tbl, title = NULL,
    Multiple_testing_correction = c("bonferroni","BH"),
    alpha = 0.05
) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  need <- c("SES_signed","q_signed")
  stopifnot(all(need %in% names(global_tbl)))

  df <- global_tbl[1, , drop = FALSE]
  df$y   <- factor("ARD vs non-ARD", levels = "ARD vs non-ARD")
  df$dir <- ifelse(df$SES_signed < 0, "Protective", "Risk")
  seg    <- data.frame(y = df$y, xmin = 0, xmax = df$SES_signed)

  exp_label <- if ("exposure" %in% names(df)) df$exposure[1] else NULL
  if (is.null(title)) title <- sprintf("Global Enrichment Analysis for %s",
                                       ifelse(is.null(exp_label) || !nzchar(exp_label), "exposure", exp_label))

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y),
                          linewidth = 0.6, alpha = 0.5) +
    ggplot2::geom_point(
      data = df, ggplot2::aes(x = SES_signed, y = y, colour = dir),
      shape = 16, size = 7
    ) +
    ggplot2::scale_colour_manual(name = "Direction:", values = c("Protective" = "grey60", "Risk" = "black")) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7)

  if (Multiple_testing_correction == "bonferroni") {
    zthr <- .z_thr_bonf(1, alpha = alpha)
    p <- p +
      ggplot2::geom_vline(xintercept = c(-zthr, zthr),
                          ggplot2::aes(linetype = sprintf("Bonferroni (α=%.2f, 2-sided)", alpha)),
                          linewidth = 0.5, alpha = 0.85) +
      ggplot2::scale_linetype_manual(
        name = "Significance:", values = setNames("dashed", sprintf("Bonferroni (α=%.2f, 2-sided)", alpha))
      )
  } else {
    ring_df <- df[is.finite(df$q_signed) & df$q_signed < alpha, , drop = FALSE]
    p <- p +
      ggplot2::geom_point(
        data = ring_df,
        ggplot2::aes(x = SES_signed, y = y, shape = "BH significant (q<α)"),
        size = 8, stroke = 1.3, colour = "firebrick", fill = NA
      ) +
      ggplot2::scale_shape_manual(
        name = "Significance:", values = c("BH significant (q<α)" = 21)
      )
  }

  p <- p +
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
  p <- .ardmr_attach_plot_data(p, main = df, segment = seg)
  p
}


# ---------- Signed enrichment violins ----------

.enrichment_violin_prepare <- function(df, group_col, extra_cols = character(), group_levels = NULL) {
  df <- tibble::as_tibble(df)
  if (!nrow(df)) {
    return(list(draws = tibble::tibble(), observed = tibble::tibble(), levels = character()))
  }
  if (!group_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in data frame.", group_col), call. = FALSE)
  }
  if (!"perm_stat" %in% names(df)) {
    stop("Input data frame must contain a 'perm_stat' list-column of permutation draws.", call. = FALSE)
  }
  missing_extra <- setdiff(extra_cols, names(df))
  if (length(missing_extra)) {
    stop(sprintf("Missing columns required for plotting: %s", paste(missing_extra, collapse = ", ")), call. = FALSE)
  }

  df$group_label <- as.character(df[[group_col]])
  keep <- !is.na(df$group_label)
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) {
    return(list(draws = tibble::tibble(), observed = tibble::tibble(), levels = character()))
  }

  if (!is.null(group_levels)) {
    levels_use <- unique(as.character(group_levels))
  } else {
    ses_vals <- df$SES_signed
    ord <- order(ses_vals, df$group_label, na.last = NA)
    if (!length(ord)) ord <- order(df$group_label)
    levels_use <- unique(df$group_label[ord])
  }
  levels_use <- levels_use[!is.na(levels_use)]
  df$group_f <- factor(df$group_label, levels = levels_use)
  df$row_id <- seq_len(nrow(df))

  draws_list <- lapply(seq_len(nrow(df)), function(i) {
    draws <- df$perm_stat[[i]]
    if (is.null(draws) || !length(draws)) return(NULL)
    mu <- df$mean_signed[i]
    sd <- df$sd_signed[i]
    ses <- if (is.finite(sd) && sd > 0) (draws - mu) / sd else rep(0, length(draws))
    tibble::tibble(
      row_id = df$row_id[i],
      group_label = df$group_label[i],
      group_f = df$group_f[i],
      Tg = draws,
      Tg_ses = ses
    )
  })

  draws_df <- dplyr::bind_rows(draws_list)
  if (nrow(draws_df)) {
    meta_cols <- unique(c("row_id", extra_cols))
    meta_cols <- meta_cols[meta_cols %in% names(df)]
    if (length(meta_cols) > 1L) {
      meta_df <- df[, meta_cols, drop = FALSE]
      draws_df <- dplyr::left_join(draws_df, meta_df, by = "row_id")
    }
    draws_df$row_id <- NULL
  } else {
    draws_df <- tibble::tibble()
  }

  obs_cols <- unique(c(group_col, extra_cols,
                       "stat_signed","mean_signed","sd_signed","SES_signed",
                       "p_signed","q_signed","n_perm_signed","exact_signed","perm_scope"))
  obs_cols <- obs_cols[obs_cols %in% names(df)]
  obs_df <- df[, obs_cols, drop = FALSE]
  obs_df$group_label <- df$group_label
  obs_df$group_f <- df$group_f

  list(draws = draws_df, observed = obs_df, levels = levels_use)
}

.build_enrichment_violin_plot <- function(
    draw_df,
    obs_df,
    orientation = c("vertical","forest"),
    alpha = 0.05,
    violin_width = 1.25
) {
  orientation <- match.arg(orientation)
  if (!nrow(obs_df)) {
    placeholder <- ggplot2::ggplot() +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::geom_text(
        ggplot2::aes(x = 0, y = 0, label = "No data available"),
        colour = "grey40"
      ) +
      ggplot2::xlab(NULL) + ggplot2::ylab(NULL)
    return(placeholder)
  }

  obs_df <- tibble::as_tibble(obs_df)
  draw_df <- tibble::as_tibble(draw_df)

  alpha_formatted <- sprintf("%.2f", alpha)
  alpha_formatted <- sub("^0\\.", ".", alpha_formatted)
  sig_labels <- c(
    sprintf("q < %s (BH)", alpha_formatted),
    sprintf("q ≥ %s (BH)", alpha_formatted)
  )
  obs_df$signif_label <- ifelse(
    is.finite(obs_df$q_signed) & obs_df$q_signed < alpha,
    sig_labels[1], sig_labels[2]
  )
  obs_df$signif_label <- factor(obs_df$signif_label, levels = sig_labels)
  obs_df$signif_label <- base::droplevels(obs_df$signif_label)

  colour_map <- stats::setNames(c("firebrick", "grey30"), sig_labels)
  active_levels <- levels(obs_df$signif_label)
  if (length(active_levels)) {
    colour_map <- colour_map[active_levels]
  }

  if (orientation == "vertical") {
    base <- ggplot2::ggplot(draw_df, ggplot2::aes(x = group_f, y = Tg_ses)) +
      ggplot2::geom_violin(
        fill = "grey85",
        colour = "grey60",
        alpha = 0.9,
        trim = FALSE,
        width = violin_width
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
      ggplot2::geom_segment(
        data = obs_df,
        ggplot2::aes(x = group_f, xend = group_f, y = 0, yend = SES_signed, colour = signif_label),
        linewidth = 0.6, alpha = 0.7
      ) +
      ggplot2::geom_point(
        data = obs_df,
        ggplot2::aes(x = group_f, y = SES_signed, colour = signif_label),
        size = 2.5
      ) +
      ggplot2::labs(x = NULL, y = expression(Standardized~T[G]~"(SES)")) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        legend.position = "top",
        legend.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 10)
      )
  } else {
    base <- ggplot2::ggplot(draw_df, ggplot2::aes(y = group_f, x = Tg_ses)) +
      ggplot2::geom_violin(
        fill = "grey85",
        colour = "grey60",
        alpha = 0.9,
        trim = FALSE,
        width = violin_width
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
      ggplot2::geom_segment(
        data = obs_df,
        ggplot2::aes(y = group_f, yend = group_f, x = 0, xend = SES_signed, colour = signif_label),
        linewidth = 0.6, alpha = 0.7
      ) +
      ggplot2::geom_point(
        data = obs_df,
        ggplot2::aes(y = group_f, x = SES_signed, colour = signif_label),
        size = 2.5
      ) +
      ggplot2::labs(y = NULL, x = expression(Standardized~T[G]~"(SES)")) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position = "top",
        legend.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 10)
      )
  }

  base +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = colour_map,
      drop = TRUE
    )
}

#' Global signed enrichment violin plot
#'
#' @param global_tbl Output of [enrichment_global_directional()] (or similar)
#'   with retained permutation draws.
#' @param orientation Either "vertical" (group on x-axis) or "forest"
#'   (group on y-axis, violin oriented horizontally).
#' @param alpha BH significance threshold for colouring the observed statistic.
#' @param title Optional plot title.
#' @param subtitle Optional subtitle.
#' @export
plot_enrichment_signed_violin_global <- function(
    global_tbl,
    orientation = c("vertical","forest"),
    alpha = 0.05,
    title = NULL,
    subtitle = NULL
) {
  orientation <- match.arg(orientation)
  df <- tibble::as_tibble(global_tbl)
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = title %||% "Global enrichment (no data)"))
  }

  df$plot_group <- "ARD vs non-ARD"
  prep <- .enrichment_violin_prepare(
    df,
    group_col = "plot_group",
    extra_cols = c("scope","exposure","pheno_sex","n_total","n_ard","n_non","perm_scope")
  )

  draws <- prep$draws
  obs <- prep$observed

  if (!nrow(draws)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = title %||% "Global enrichment", subtitle = "No permutation draws available"))
  }

  exp_label <- tryCatch(unique(na.omit(obs$exposure))[1], error = function(e) NULL)
  if (is.null(exp_label) || !nzchar(exp_label)) exp_label <- "exposure"
  if (is.null(title)) {
    title <- sprintf("Global enrichment of %s", exp_label)
  }
  if (is.null(subtitle)) {
    subtitle <- "ARD vs non-ARD, two-sided permutation test"
  }

  p <- .build_enrichment_violin_plot(draws, obs, orientation = orientation, alpha = alpha) +
    ggplot2::labs(title = title, subtitle = subtitle)

  if (orientation == "vertical") {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))
  } else {
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
  }

  .ardmr_attach_plot_data(p, draws = draws, observed = obs)
}

#' Global signed enrichment violin plot across ARD compare groups
#'
#' @param combined_global_tbl Data frame containing one row per group with
#'   retained permutation draws and associated metadata. Expected to include
#'   columns `perm_stat`, `group_axis_label`, and `group_order` in addition to
#'   the enrichment outputs.
#' @param alpha BH significance threshold for colouring the observed statistic.
#' @param title Optional plot title.
#' @param subtitle Optional subtitle.
plot_enrichment_signed_violin_global_compare <- function(
    combined_global_tbl,
    alpha = 0.05,
    title = NULL,
    subtitle = NULL
) {
  df <- tibble::as_tibble(combined_global_tbl)
  if (!nrow(df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::labs(title = title %||% "Global enrichment (no data)")
    )
  }

  required_cols <- c("perm_stat", "group_axis_label", "group_order")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols)) {
    stop(
      sprintf(
        "Combined global enrichment table is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  df <- df[order(df$group_order, seq_len(nrow(df))), , drop = FALSE]
  df$plot_group <- df$group_axis_label
  axis_levels <- unique(df$group_axis_label)
  axis_levels <- axis_levels[!is.na(axis_levels)]
  if (!length(axis_levels)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::labs(title = title %||% "Global enrichment (no groups)")
    )
  }

  prep <- .enrichment_violin_prepare(
    df,
    group_col = "plot_group",
    extra_cols = c(
      "group_id", "group_order", "group_axis_label", "sex", "ancestry",
      "scope", "exposure", "pheno_sex", "n_total", "n_ard", "n_non",
      "perm_scope"
    ),
    group_levels = axis_levels
  )

  draws <- prep$draws
  obs <- prep$observed

  if (!nrow(draws)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::labs(
          title = title %||% "Global enrichment (compare)",
          subtitle = "No permutation draws available"
        )
    )
  }

  exp_label <- tryCatch(unique(na.omit(obs$exposure))[1], error = function(e) NULL)
  if (is.null(exp_label) || !nzchar(exp_label)) exp_label <- "exposure"
  if (is.null(title)) {
    title <- sprintf("Global enrichment of %s across groups", exp_label)
  }
  if (is.null(subtitle)) {
    subtitle <- "ARD vs non-ARD, two-sided permutation test"
  }

  p <- .build_enrichment_violin_plot(
    draws,
    obs,
    orientation = "forest",
    alpha = alpha,
    violin_width = 0.75
  ) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 27)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))

  .ardmr_attach_plot_data(p, draws = draws, observed = obs)
}

#' Cause-level signed enrichment violin plot
#'
#' @param by_cause_tbl Table returned by [enrichment_by_cause_directional()].
#' @param level Cause level identifier (e.g., "cause_level_1").
#' @param compare_mode Comparison mode (default "cause_vs_rest_all").
#' @inheritParams plot_enrichment_signed_violin_global
#' @export
plot_enrichment_signed_violin_by_cause <- function(
    by_cause_tbl,
    level,
    compare_mode = "cause_vs_rest_all",
    orientation = c("vertical","forest"),
    alpha = 0.05,
    title = NULL,
    subtitle = NULL
) {
  orientation <- match.arg(orientation)
  df <- tibble::as_tibble(by_cause_tbl)
  df <- df[df$level == level & df$compare_mode == compare_mode, , drop = FALSE]
  if (!nrow(df)) {
    pretty_level <- gsub("_", " ", level)
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = title %||% sprintf("%s enrichment (no data)", pretty_level)))
  }

  df$plot_group <- df$cause
  prep <- .enrichment_violin_prepare(
    df,
    group_col = "plot_group",
    extra_cols = c("level","cause","compare_mode","exposure","pheno_sex","n_total","n_pos","n_neg","perm_scope")
  )

  draws <- prep$draws
  obs <- prep$observed

  if (!nrow(draws)) {
    pretty_level <- sprintf("Cause level %d", .level_number(level))
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = title %||% sprintf("%s enrichment", pretty_level), subtitle = "No permutation draws available"))
  }

  exp_label <- tryCatch(unique(na.omit(obs$exposure))[1], error = function(e) NULL)
  if (is.null(exp_label) || !nzchar(exp_label)) exp_label <- "exposure"
  if (is.null(title)) {
    title <- sprintf("Cause level %d signed enrichment of %s", .level_number(level), exp_label)
  }
  if (is.null(subtitle)) {
    subtitle <- "Within-cause effect size vs all-cause effect size, two-sided permutation test"
  }

  p <- .build_enrichment_violin_plot(draws, obs, orientation = orientation, alpha = alpha) +
    ggplot2::labs(title = title, subtitle = subtitle)

  if (orientation == "vertical") {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9))
  } else {
    p <- p +
      ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 27)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
  }

  .ardmr_attach_plot_data(p, draws = draws, observed = obs)
}


#' Heatmap of signed enrichment p-values across compare groups
#'
#' @param heatmap_df Data frame containing one row per cause × group cell.
#'   Expected columns include: `cause`, `group_label`, `p_signed`,
#'   `SES_signed`, `significant`, `direction`, `label`, `cause_order`,
#'   `group_order`.
#' @param level Cause level identifier (e.g., "cause_level_1").
#' @param scope Scope string ("all_diseases" or "age_related_diseases").
#' @param exposure Optional exposure label for the title.
#' @param title,subtitle Optional strings overriding the defaults.
#' @keywords internal
plot_enrichment_group_heatmap <- function(
    heatmap_df,
    level,
    scope,
    exposure = NULL,
    title = NULL,
    subtitle = NULL
) {
  df <- tibble::as_tibble(heatmap_df)
  required <- c(
    "cause","group_label","p_signed","SES_signed","q_value",
    "significant","direction","label","cause_order","group_order"
  )
  if (!all(required %in% names(df))) {
    stop("plot_enrichment_group_heatmap(): heatmap_df missing required columns.", call. = FALSE)
  }
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  df$cause <- as.character(df$cause)
  df$group_label <- as.character(df$group_label)
  cause_levels <- unique(df$cause[order(df$cause_order, df$cause)])
  group_levels <- unique(df$group_label[order(df$group_order, df$group_label)])
  df$cause_f <- factor(df$cause, levels = rev(cause_levels))
  df$group_f <- factor(df$group_label, levels = group_levels)
  df$significant <- df$significant %in% TRUE

  missing_cells <- !is.finite(df$p_signed) | !is.finite(df$q_value)

  dir_vals <- ifelse(df$direction %in% c("protective","risk"), df$direction, "neutral")
  dir_vals[!df$significant] <- "neutral"
  dir_vals[missing_cells] <- NA_character_
  df$fill_state <- factor(dir_vals, levels = c("protective","neutral","risk"))

  palette <- c(
    protective = "#c6dbef",
    neutral = "#e0e0e0",
    risk = "#fcbba1"
  )

  exp_label <- tryCatch(as.character(exposure)[1], error = function(e) NA_character_)
  if (is.null(exp_label) || is.na(exp_label) || !nzchar(exp_label)) {
    exp_label <- "exposure"
  }

  if (is.null(title)) {
    lvl_num <- suppressWarnings(.level_number(level))
    if (is.na(lvl_num) || !is.finite(lvl_num)) {
      title <- sprintf("Enrichment heatmap for %s", exp_label)
    } else {
      title <- sprintf("Cause level %d enrichment heatmap for %s", lvl_num, exp_label)
    }
  }

  if (is.null(subtitle)) {
    scope_label <- switch(
      scope,
      age_related_diseases = "Age-related diseases",
      all_diseases = "All diseases",
      scope
    )
    subtitle <- sprintf("%s scope · cause vs rest (BH-adjusted q-values)", scope_label)
  }

  df$label <- as.character(df$label)
  df$label[is.na(df$label) | !nzchar(df$label)] <- NA_character_

  p <- ggplot2::ggplot(df, ggplot2::aes(x = group_f, y = cause_f, fill = fill_state)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2, colour = "#1a1a1a", na.rm = TRUE) +
    ggplot2::scale_fill_manual(
      values = palette,
      limits = c("protective","neutral","risk"),
      breaks = c("protective","neutral","risk"),
      labels = c(
        "Protective (BH q-val < .05)",
        "Not significant",
        "Risk (BH q-val < .05)"
      ),
      name = NULL,
      drop = FALSE,
      na.value = "#ffffff",
      na.translate = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(hjust = 1),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 9)
    )

  .ardmr_attach_plot_data(p, main = df)
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
#' @param dot_names Logical; if TRUE, label significant points with their outcome names.
#' @export
volcano_plot <- function(results_df,
                         Multiple_testing_correction = c("BH","bonferroni"),
                         qc_pass_only = TRUE,
                         alpha = 0.05,
                         exposure = NULL,
                         label_top_n = 25,
                         verbose = TRUE,
                         dot_names = TRUE) {
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
  plot_title <- sprintf("Phenome-wide MR of %s", exposure_lab)

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
  lab_df <- tibble::tibble()
  if (isTRUE(dot_names) && "results_outcome" %in% names(df)) {
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

  p <- .ardmr_attach_plot_data(p, main = df, labelled = lab_df)
  p
}


#' Volcano plot with recoloured significance highlighting
#'
#' Generates a volcano plot of MR IVW estimates, colouring each point by
#' significance (FDR or Bonferroni) and effect direction. Non-significant points
#' are shown in dark grey, significant protective associations (β < 0) in blue,
#' and significant risk associations (β ≥ 0) in red.
#'
#' @inheritParams volcano_plot
#' @param dot_names Logical; if `TRUE`, annotate significant outcomes on the
#'   plot (default `TRUE`).
#' @export
volcano_plot_recolor <- function(results_df,
                                 Multiple_testing_correction = c("BH","bonferroni"),
                                 qc_pass_only = TRUE,
                                 alpha = 0.05,
                                 exposure = NULL,
                                 dot_names = TRUE,
                                 label_top_n = 25,
                                 verbose = TRUE) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  dot_names <- isTRUE(dot_names)

  if (is.null(results_df) || !nrow(results_df)) {
    if (verbose && requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn("Volcano (recolor): no rows.")
    }
    return(ggplot2::ggplot() + ggplot2::labs(title = "Volcano (recolor, no results)"))
  }

  df <- tibble::as_tibble(results_df)

  if ("results_qc_pass" %in% names(df) && isTRUE(qc_pass_only)) {
    if (!is.logical(df$results_qc_pass)) {
      df$results_qc_pass <- df$results_qc_pass %in% c(TRUE,"TRUE","True","true",1,"1","T")
    }
    df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  }

  if (!all(c("results_beta_ivw","results_p_ivw") %in% names(df)) &&
      !"results_log10p_ivw" %in% names(df)) {
    stop("volcano_plot_recolor(): need 'results_beta_ivw' and either 'results_p_ivw' or 'results_log10p_ivw'.",
         call. = FALSE)
  }

  df <- dplyr::mutate(
    df,
    logp = dplyr::if_else(
      !is.na(.data$results_log10p_ivw),
      -.data$results_log10p_ivw,
      -log10(pmax(.data$results_p_ivw, .Machine$double.xmin))
    )
  )

  df <- dplyr::filter(df, is.finite(.data$results_beta_ivw), is.finite(.data$logp))
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "Volcano (recolor, no finite rows)"))
  }

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

  df <- df |>
    dplyr::mutate(
      colour_group = dplyr::case_when(
        .data$sig %in% TRUE & is.finite(.data$results_beta_ivw) & .data$results_beta_ivw < 0 ~
          "Significant protective (β<0)",
        .data$sig %in% TRUE & is.finite(.data$results_beta_ivw) & .data$results_beta_ivw >= 0 ~
          "Significant risk (β≥0)",
        TRUE ~ "Not significant"
      )
    )

  colour_levels <- c("Not significant", "Significant protective (β<0)", "Significant risk (β≥0)")
  present_levels <- colour_levels[colour_levels %in% unique(df$colour_group)]
  df$colour_group <- factor(df$colour_group, levels = present_levels)

  exposure_lab <- tryCatch(as.character(exposure)[1], error = function(e) NA_character_)
  if (is.null(exposure_lab) || is.na(exposure_lab) || !nzchar(exposure_lab)) exposure_lab <- "exposure"
  plot_title <- sprintf("Phenome-wide MR of %s", exposure_lab)

  legend_title <- if (Multiple_testing_correction == "BH") {
    "Significance & direction (FDR < 0.05)"
  } else {
    "Significance & direction:"
  }

  colour_values <- c(
    "Not significant" = "#4B4B4B",
    "Significant protective (β<0)" = "#1f78b4",
    "Significant risk (β≥0)" = "#e31a1c"
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$results_beta_ivw, y = .data$logp)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$colour_group), alpha = 0.9, size = 1.7) +
    ggplot2::scale_colour_manual(
      name = legend_title,
      values = colour_values,
      drop = FALSE
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

  if (Multiple_testing_correction == "bonferroni" && !is.na(thr_y)) {
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

  lab_df <- tibble::tibble()
  if (isTRUE(dot_names) && "results_outcome" %in% names(df)) {
    lab_df <- dplyr::filter(df, .data$sig %in% TRUE & !is.na(.data$results_outcome))
    if (nrow(lab_df)) {
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

  p <- .ardmr_attach_plot_data(p, main = df, labelled = lab_df)
  p
}


#' Global ARD vs non-ARD summary (two-dot)
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
# ========== Global (directional) one-row summary ==========
# ========== Global (directional) one-row summary ==========
plot_enrichment_global <- function(
    global_tbl, title = NULL,
    Multiple_testing_correction = c("bonferroni","BH"),
    alpha = 0.05
) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  stopifnot(all(c("SES_prot","SES_risk","q_prot","q_risk") %in% names(global_tbl)))
  df <- global_tbl[1, , drop = FALSE]
  df$y <- factor("ARD vs non-ARD", levels = "ARD vs non-ARD")
  seg  <- data.frame(y = df$y, xmin = df$SES_prot, xmax = df$SES_risk)

  exp_label <- if ("exposure" %in% names(df)) df$exposure[1] else NULL
  if (is.null(title)) title <- sprintf("Global Enrichment Analysis for %s",
                                       ifelse(is.null(exp_label) || !nzchar(exp_label), "exposure", exp_label))

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y),
                          linewidth = 0.6, alpha = 0.5) +
    ggplot2::geom_point(data = df, ggplot2::aes(x = SES_prot, y = y),
                        shape = 16, size = 7, colour = "grey60") +
    ggplot2::geom_point(data = df, ggplot2::aes(x = SES_risk, y = y),
                        shape = 16, size = 7, colour = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7)

  if (Multiple_testing_correction == "bonferroni") {
    zthr <- .z_thr_bonf(1, alpha = alpha)
    p <- p +
      ggplot2::geom_vline(xintercept = c(-zthr, zthr),
                          ggplot2::aes(linetype = sprintf("Bonferroni (α=%.2f, 2-sided)", alpha)),
                          linewidth = 0.5, alpha = 0.85) +
      ggplot2::scale_linetype_manual(
        name = "Significance:", values = setNames("dashed", sprintf("Bonferroni (α=%.2f, 2-sided)", alpha))
      )
  } else {
    ring_df_prot <- df[is.finite(df$q_prot) & df$q_prot < alpha, , drop = FALSE]
    ring_df_risk <- df[is.finite(df$q_risk) & df$q_risk < alpha, , drop = FALSE]
    p <- p +
      ggplot2::geom_point(data = ring_df_prot,
                          ggplot2::aes(x = SES_prot, y = y, shape = "BH significant (q<α)"),
                          size = 8, stroke = 1.3, colour = "firebrick", fill = NA) +
      ggplot2::geom_point(data = ring_df_risk,
                          ggplot2::aes(x = SES_risk,  y = y, shape = "BH significant (q<α)"),
                          size = 8, stroke = 1.3, colour = "firebrick", fill = NA) +
      ggplot2::scale_shape_manual(
        name = "Significance:", values = c("BH significant (q<α)" = 21)
      )
  }

  p <- p +
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
  p <- .ardmr_attach_plot_data(p, main = df, segments = seg)
  p
}

# ---------- forest for Δβ (inside vs outside) ----------
#' Forest of Δβ (inside - outside) by cause with multiplicity marks
#'
#' X-axis = Δβ on the original MR beta scale (e.g., log-OR per 1 SD exposure).
#' Points = Δβ; whiskers = 95% CI; vertical line at 0.
#' BH: red rings around q<alpha; Bonferroni: dashed legend line & bold labels optional.
#'
#' @param beta_tbl Output of beta_contrast_by_cause().
#' @param title,subtitle Optional strings.
#' @param Multiple_testing_correction "BH" or "bonferroni" (must match the table).
#' @param alpha Significance level used for the table.
#' @return ggplot
#' @export
plot_beta_contrast_forest <- function(
    beta_tbl,
    title = NULL,
    subtitle = "Within-cause mean β vs outside all other causes",
    Multiple_testing_correction = c("BH","bonferroni"),
    alpha = 0.05
) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  df <- tibble::as_tibble(beta_tbl)
  stopifnot(all(c("cause","delta_beta","se_delta","ci_low","ci_high","p","q","sig") %in% names(df)))
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void())

  # order by |Δβ|
  df$._ord <- abs(df$delta_beta)
  df <- df[order(df$._ord, decreasing = FALSE), , drop = FALSE]
  df$cause_f <- factor(df$cause, levels = df$cause)

  # default title from exposure if present
  if (is.null(title)) {
    exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NA_character_
    title <- sprintf("Cause-level contrast of mean MR β for %s",
                     ifelse(is.na(exp_label) || !nzchar(exp_label), "exposure", exp_label))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(y = cause_f)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ci_low, xmax = ci_high),
      height = 0.2, linewidth = 0.6, alpha = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = delta_beta),
      shape = 16, size = 2.8, colour = "black"
    ) +
    ggplot2::labs(
      x = expression(Delta*beta ~ "(inside − outside, log-OR per 1 SD exposure)"),
      y = "Cause",
      title = title, subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 11),
      legend.text  = ggplot2::element_text(size = 11)
    )

  if (Multiple_testing_correction == "BH") {
    p <- p +
      ggplot2::geom_point(
        data = df[df$sig %in% TRUE, , drop = FALSE],
        ggplot2::aes(x = delta_beta, y = cause_f, shape = "BH significant (q<α)"),
        size = 4.2, stroke = 1.1, colour = "firebrick", fill = NA
      ) +
      ggplot2::scale_shape_manual(
        name = "Significance:",
        values = c("BH significant (q<α)" = 21)
      )
  } else {
    # Bonferroni: show legend item + (optional) bold significant labels
    lab_col <- ifelse(df$sig, "black", "grey30")
    p <- p +
      ggplot2::scale_linetype_manual(
        name = "Significance:", values = c("Bonferroni threshold" = "dashed")
      ) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(colour = lab_col))
  }

  p <- .ardmr_attach_plot_data(p, main = df)
  p
}

#' Wrapped-axis companion to `plot_beta_contrast_forest()`
#'
#' Applies a 27-character wrap to the displayed cause labels while keeping
#' the attached data unchanged.
#'
#' @inheritParams plot_beta_contrast_forest
#' @return ggplot object with wrapped y-axis labels.
#' @export
plot_beta_contrast_forest_wrap <- function(
    beta_tbl,
    title = NULL,
    subtitle = "Within-cause mean β vs outside all other causes",
    Multiple_testing_correction = c("BH","bonferroni"),
    alpha = 0.05
) {
  df <- tibble::as_tibble(beta_tbl)
  p <- plot_beta_contrast_forest(
    beta_tbl = df,
    title = title,
    subtitle = subtitle,
    Multiple_testing_correction = Multiple_testing_correction,
    alpha = alpha
  )

  plot_data <- attr(p, "ardmr_plot_data", exact = TRUE)
  if (!nrow(df)) {
    if (!is.null(plot_data)) {
      attr(p, "ardmr_plot_data") <- plot_data
    }
    return(p)
  }
  p_wrapped <- p + ggplot2::scale_y_discrete(
    labels = function(x) stringr::str_wrap(x, width = 27)
  )
  if (!is.null(plot_data)) {
    attr(p_wrapped, "ardmr_plot_data") <- plot_data
  }
  p_wrapped
}


# helper to format the β-scale axis label using exposure units
.format_beta_axis_label <- function(exposure_units, effect_scale = c("odds_ratio", "absolute_risk")) {
  effect_scale <- match.arg(effect_scale)
  if (missing(exposure_units)) stop("`exposure_units` is required.", call. = FALSE)
  units_label <- as.character(exposure_units)
  if (!length(units_label)) stop("`exposure_units` is required.", call. = FALSE)
  units_label <- units_label[1]
  if (is.na(units_label)) stop("`exposure_units` is required.", call. = FALSE)
  units_label <- trimws(units_label)
  if (!nzchar(units_label)) stop("`exposure_units` must be a non-empty string.", call. = FALSE)
  label_template <- switch(
    effect_scale,
    odds_ratio = "Mean Δ odds ratio per %s",
    absolute_risk = "Mean Δ absolute risk per %s"
  )
  sprintf(label_template, units_label)
}


# ---------- forest for IVW mean β ----------
#' Forest plot of IVW mean MR beta by cause
#'
#' @param beta_tbl Tibble with columns cause, ivw_mean_beta, se_ivw_mean,
#'   ci_low, ci_high.
#' @param title Optional plot title.
#' @param subtitle Optional subtitle string.
#' @param exposure_units Character string describing the exposure units used
#'   for the β-scale axis label.
#' @param effect_scale Character string specifying whether the β-axis should be
#'   expressed as odds ratios ("odds_ratio") or absolute risk
#'   ("absolute_risk").
plot_beta_mean_forest <- function(
    beta_tbl,
    title = NULL,
    subtitle = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_tbl)
  required_cols <- c("cause", "ivw_mean_beta", "se_ivw_mean", "ci_low", "ci_high")
  stopifnot(all(required_cols %in% names(df)))
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void())

  effect_scale <- match.arg(effect_scale)
  axis_label <- .format_beta_axis_label(exposure_units, effect_scale = effect_scale)

  df$._ord <- abs(df$ivw_mean_beta)
  df <- df[order(df$._ord, decreasing = FALSE), , drop = FALSE]
  df$cause_f <- factor(df$cause, levels = df$cause)

  if (is.null(title)) {
    exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NA_character_
    title <- sprintf("IVW mean MR β by cause (%s)",
                     ifelse(is.na(exp_label) || !nzchar(exp_label), "exposure", exp_label))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(y = cause_f)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ci_low, xmax = ci_high),
      height = 0.2, linewidth = 0.6, alpha = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = ivw_mean_beta),
      shape = 16, size = 2.8, colour = "black"
    ) +
    ggplot2::labs(
      x = axis_label,
      y = "Cause",
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6))
    )

  p <- .ardmr_attach_plot_data(p, main = df)
  p
}

#' Wrapped-axis companion to `plot_beta_mean_forest()`
#'
#' Applies a 27-character wrap to the displayed cause labels while keeping
#' the attached data unchanged.
#'
#' @inheritParams plot_beta_mean_forest
#' @return ggplot object with wrapped y-axis labels.
#' @export
plot_beta_mean_forest_wrap <- function(
    beta_tbl,
    title = NULL,
    subtitle = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_tbl)
  effect_scale <- match.arg(effect_scale)
  p <- plot_beta_mean_forest(
    beta_tbl = df,
    title = title,
    subtitle = subtitle,
    exposure_units = exposure_units,
    effect_scale = effect_scale
  )

  plot_data <- attr(p, "ardmr_plot_data", exact = TRUE)
  if (!nrow(df)) {
    if (!is.null(plot_data)) {
      attr(p, "ardmr_plot_data") <- plot_data
    }
    return(p)
  }
  p_wrapped <- p + ggplot2::scale_y_discrete(
    labels = function(x) stringr::str_wrap(x, width = 27)
  )
  if (!is.null(plot_data)) {
    attr(p_wrapped, "ardmr_plot_data") <- plot_data
  }
  p_wrapped
}

#' Wrapped-axis companion to `plot_beta_mean_forest()` with floating height
#'
#' Identical to `plot_beta_mean_forest_wrap()` but intended for saving without
#' enforcing a fixed output height.
#'
#' @inheritParams plot_beta_mean_forest
#' @return ggplot object with wrapped y-axis labels.
#' @export
plot_beta_mean_forest_wrap_yfloat <- function(
    beta_tbl,
    title = NULL,
    subtitle = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_tbl)
  effect_scale <- match.arg(effect_scale)
  p <- plot_beta_mean_forest(
    beta_tbl = df,
    title = title,
    subtitle = subtitle,
    exposure_units = exposure_units,
    effect_scale = effect_scale
  )

  plot_data <- attr(p, "ardmr_plot_data", exact = TRUE)
  if (!nrow(df)) {
    if (!is.null(plot_data)) {
      attr(p, "ardmr_plot_data") <- plot_data
    }
    return(p)
  }
  p_wrapped <- p + ggplot2::scale_y_discrete(
    labels = function(x) stringr::str_wrap(x, width = 27)
  )
  if (!is.null(plot_data)) {
    attr(p_wrapped, "ardmr_plot_data") <- plot_data
  }
  p_wrapped
}

#' Global forest plot of IVW mean MR beta (all diseases vs ARDs)
#'
#' @param beta_global_tbl Tibble from beta_mean_global().
#' @param title Optional title.
#' @param exposure_units Character string describing the exposure units used
#'   for the β-scale axis label.
#' @param effect_scale Character string specifying whether the β-axis should be
#'   expressed as odds ratios ("odds_ratio") or absolute risk
#'   ("absolute_risk").
plot_beta_mean_global <- function(
    beta_global_tbl,
    title = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_global_tbl)
  required_cols <- c("group", "ivw_mean_beta", "se_ivw_mean", "ci_low", "ci_high")
  stopifnot(all(required_cols %in% names(df)))
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void())

  effect_scale <- match.arg(effect_scale)
  axis_label <- .format_beta_axis_label(exposure_units, effect_scale = effect_scale)

  df$._ord <- ifelse(is.finite(df$ivw_mean_beta), abs(df$ivw_mean_beta), Inf)
  df <- df[order(df$._ord, decreasing = FALSE), , drop = FALSE]
  df$group_f <- factor(df$group, levels = df$group)

  if (is.null(title)) {
    exp_label <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure))[1] else NA_character_
    title <- sprintf("Global IVW mean MR β (%s)",
                     ifelse(is.na(exp_label) || !nzchar(exp_label), "exposure", exp_label))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(y = group_f)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    ggplot2::geom_errorbarh(
      data = df[is.finite(df$ci_low) & is.finite(df$ci_high), , drop = FALSE],
      ggplot2::aes(xmin = ci_low, xmax = ci_high),
      height = 0.25, linewidth = 0.6, alpha = 0.7
    ) +
    ggplot2::geom_point(
      data = df[is.finite(df$ivw_mean_beta), , drop = FALSE],
      ggplot2::aes(x = ivw_mean_beta),
      shape = 16, size = 3.2, colour = "black"
    ) +
    ggplot2::labs(
      x = axis_label,
      y = NULL,
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6))
    )

  p <- .ardmr_attach_plot_data(p, main = df)
  p
}

#' @keywords internal
#' @param exposure_units Character string describing the exposure units used
#'   for the β-scale axis label.
plot_beta_mean_global_compare <- function(
    beta_global_tbl,
    title = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_global_tbl)
  required_cols <- c("display_label", "ivw_mean_beta", "ci_low", "ci_high")
  stopifnot(all(required_cols %in% names(df)))
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void())

  effect_scale <- match.arg(effect_scale)
  axis_label <- .format_beta_axis_label(exposure_units, effect_scale = effect_scale)

  df$axis_label <- df$display_label
  levels_vec <- rev(df$display_label)
  levels_vec <- unique(levels_vec)
  df$axis_f <- factor(df$display_label, levels = levels_vec)
  axis_labels <- stats::setNames(levels_vec, levels_vec)

  if (is.null(title)) {
    exposure_vals <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure)) else character()
    exposure_lab <- tryCatch(as.character(exposure_vals)[1], error = function(e) NA_character_)
    if (is.na(exposure_lab) || !nzchar(exposure_lab)) {
      exposure_lab <- "the exposure"
    }
    title <- sprintf("Mean effect of %s by group", exposure_lab)
  }

  err_data <- df[is.finite(df$ci_low) & is.finite(df$ci_high), , drop = FALSE]
  point_data <- df[is.finite(df$ivw_mean_beta), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(y = .data$axis_f)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    {if (nrow(err_data)) ggplot2::geom_errorbarh(
      data = err_data,
      ggplot2::aes(xmin = .data$ci_low, xmax = .data$ci_high),
      height = 0.25,
      linewidth = 0.6,
      alpha = 0.7
    ) else NULL} +
    {if (nrow(point_data)) ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = .data$ivw_mean_beta),
      shape = 16,
      size = 3.2,
      colour = "black"
    ) else NULL} +
    ggplot2::labs(
      x = axis_label,
      y = NULL,
      title = title
    ) +
    ggplot2::scale_y_discrete(labels = axis_labels) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6))
    )

  p <- .ardmr_attach_plot_data(p, main = df)
  p
}

#' Wrapped-axis companion to `plot_beta_mean_global_compare()`
#'
#' Applies a 27-character wrap to the displayed group labels while keeping
#' the attached data unchanged.
#'
#' @inheritParams plot_beta_mean_global_compare
#' @keywords internal
plot_beta_mean_global_compare_wrap <- function(
    beta_global_tbl,
    title = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_global_tbl)
  effect_scale <- match.arg(effect_scale)
  p <- if (missing(title)) {
    plot_beta_mean_global_compare(
      beta_global_tbl = df,
      exposure_units = exposure_units,
      effect_scale = effect_scale
    )
  } else {
    plot_beta_mean_global_compare(
      beta_global_tbl = df,
      title = title,
      exposure_units = exposure_units,
      effect_scale = effect_scale
    )
  }

  plot_data <- attr(p, "ardmr_plot_data", exact = TRUE)
  if (!nrow(df)) {
    if (!is.null(plot_data)) {
      attr(p, "ardmr_plot_data") <- plot_data
    }
    return(p)
  }

  levels_vec <- rev(df$display_label)
  levels_vec <- unique(levels_vec)
  axis_labels <- stats::setNames(levels_vec, levels_vec)
  wrapped_labels <- stringr::str_wrap(unname(axis_labels), width = 27)
  names(wrapped_labels) <- names(axis_labels)

  p_wrapped <- p + ggplot2::scale_y_discrete(
    limits = levels_vec,
    labels = wrapped_labels
  )
  if (!is.null(plot_data)) {
    attr(p_wrapped, "ardmr_plot_data") <- plot_data
  }
  p_wrapped
}

#' @keywords internal
plot_beta_mean_global_compare_wrap_yfloat <- function(
    beta_global_tbl,
    title = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_global_tbl)
  effect_scale <- match.arg(effect_scale)
  p <- if (missing(title)) {
    plot_beta_mean_global_compare(
      beta_global_tbl = df,
      exposure_units = exposure_units,
      effect_scale = effect_scale
    )
  } else {
    plot_beta_mean_global_compare(
      beta_global_tbl = df,
      title = title,
      exposure_units = exposure_units,
      effect_scale = effect_scale
    )
  }

  plot_data <- attr(p, "ardmr_plot_data", exact = TRUE)
  if (!nrow(df)) {
    if (!is.null(plot_data)) {
      attr(p, "ardmr_plot_data") <- plot_data
    }
    return(p)
  }

  levels_vec <- rev(df$display_label)
  levels_vec <- unique(levels_vec)
  axis_labels <- stats::setNames(levels_vec, levels_vec)
  wrapped_labels <- stringr::str_wrap(unname(axis_labels), width = 27)
  names(wrapped_labels) <- names(axis_labels)

  p_wrapped <- p + ggplot2::scale_y_discrete(
    limits = levels_vec,
    labels = wrapped_labels
  )
  if (!is.null(plot_data)) {
    attr(p_wrapped, "ardmr_plot_data") <- plot_data
  }
  p_wrapped
}

#' @keywords internal
#' @param exposure_units Character string describing the exposure units used
#'   for the β-scale axis label.
#' @param effect_scale Character string specifying whether the β-axis should be
#'   expressed as odds ratios ("odds_ratio") or absolute risk
#'   ("absolute_risk").
plot_beta_mean_cause_compare <- function(
    beta_cause_tbl,
    title = NULL,
    subtitle = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_cause_tbl)
  required_cols <- c("axis_label", "axis_id", "is_separator", "ivw_mean_beta", "ci_low", "ci_high")
  stopifnot(all(required_cols %in% names(df)))
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void())

  effect_scale <- match.arg(effect_scale)
  axis_label <- .format_beta_axis_label(exposure_units, effect_scale = effect_scale)

  df$is_separator <- df$is_separator %in% TRUE
  axis_levels <- rev(df$axis_id)
  axis_levels <- unique(axis_levels)
  df$axis_f <- factor(df$axis_id, levels = axis_levels)
  axis_labels <- stats::setNames(df$axis_label, df$axis_id)
  axis_labels <- axis_labels[axis_levels]

  if (is.null(title)) {
    exposure_vals <- if ("exposure" %in% names(df)) unique(na.omit(df$exposure)) else character()
    exposure_lab <- tryCatch(as.character(exposure_vals)[1], error = function(e) NA_character_)
    if (is.na(exposure_lab) || !nzchar(exposure_lab)) {
      exposure_lab <- "the exposure"
    }

    level_vals <- if ("level" %in% names(df)) as.character(df$level) else character()
    level_idx <- which(!is.na(level_vals) & nzchar(level_vals))
    level_value <- if (length(level_idx)) level_vals[level_idx[1]] else NA_character_
    if (!is.na(level_value) && nzchar(level_value)) {
      cause_level <- trimws(gsub("_", " ", tolower(level_value)))
    } else {
      cause_level <- "cause level"
    }

    scope_vals <- if ("scope" %in% names(df)) as.character(df$scope) else character()
    scope_idx <- which(!is.na(scope_vals) & nzchar(scope_vals))
    scope_value <- if (length(scope_idx)) scope_vals[scope_idx[1]] else NA_character_

    if (identical(scope_value, "all_diseases")) {
      title <- sprintf("Mean effect of %s on all diseases by %s", exposure_lab, cause_level)
    } else if (identical(scope_value, "age_related_diseases")) {
      title <- sprintf("Mean effect of %s on ARDs by %s", exposure_lab, cause_level)
    } else {
      scope_label <- if (!is.na(scope_value) && nzchar(scope_value)) {
        if (scope_value == "age_related_diseases") {
          "Age-Related Diseases"
        } else if (scope_value == "all_diseases") {
          "All Diseases"
        } else {
          scope_value
        }
      } else {
        "All Diseases"
      }
      title <- sprintf("%s IVW mean MR β (%s)", tools::toTitleCase(cause_level), scope_label)
    }
  }

  err_data <- df[!df$is_separator & is.finite(df$ci_low) & is.finite(df$ci_high), , drop = FALSE]
  point_data <- df[!df$is_separator & is.finite(df$ivw_mean_beta), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(y = .data$axis_f)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    {if (nrow(err_data)) ggplot2::geom_errorbarh(
      data = err_data,
      ggplot2::aes(xmin = .data$ci_low, xmax = .data$ci_high),
      height = 0.2,
      linewidth = 0.6,
      alpha = 0.7
    ) else NULL} +
    {if (nrow(point_data)) ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = .data$ivw_mean_beta),
      shape = 16,
      size = 2.8,
      colour = "black"
    ) else NULL} +
    ggplot2::labs(
      x = axis_label,
      y = NULL,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::scale_y_discrete(labels = axis_labels) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6))
    )

  p <- .ardmr_attach_plot_data(p, main = df)
  p
}

#' Wrapped-axis companion to `plot_beta_mean_cause_compare()`
#'
#' Applies a 27-character wrap to the displayed cause labels while keeping
#' the attached data unchanged.
#'
#' @inheritParams plot_beta_mean_cause_compare
#' @keywords internal
plot_beta_mean_cause_compare_wrap <- function(
    beta_cause_tbl,
    title = NULL,
    subtitle = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_cause_tbl)
  effect_scale <- match.arg(effect_scale)
  p <- plot_beta_mean_cause_compare(
    beta_cause_tbl = df,
    title = title,
    subtitle = subtitle,
    exposure_units = exposure_units,
    effect_scale = effect_scale
  )

  plot_data <- attr(p, "ardmr_plot_data", exact = TRUE)
  if (!nrow(df)) {
    if (!is.null(plot_data)) {
      attr(p, "ardmr_plot_data") <- plot_data
    }
    return(p)
  }

  df$is_separator <- df$is_separator %in% TRUE
  axis_levels <- rev(df$axis_id)
  axis_levels <- unique(axis_levels)
  axis_labels <- stats::setNames(df$axis_label, df$axis_id)
  axis_labels <- axis_labels[axis_levels]
  wrapped_labels <- stringr::str_wrap(unname(axis_labels), width = 27)
  names(wrapped_labels) <- names(axis_labels)

  p_wrapped <- p + ggplot2::scale_y_discrete(
    limits = axis_levels,
    labels = wrapped_labels
  )
  if (!is.null(plot_data)) {
    attr(p_wrapped, "ardmr_plot_data") <- plot_data
  }
  p_wrapped
}

#' @keywords internal
plot_beta_mean_cause_compare_wrap_yfloat <- function(
    beta_cause_tbl,
    title = NULL,
    subtitle = NULL,
    exposure_units,
    effect_scale = c("odds_ratio", "absolute_risk")
) {
  df <- tibble::as_tibble(beta_cause_tbl)
  effect_scale <- match.arg(effect_scale)
  p <- plot_beta_mean_cause_compare(
    beta_cause_tbl = df,
    title = title,
    subtitle = subtitle,
    exposure_units = exposure_units,
    effect_scale = effect_scale
  )

  plot_data <- attr(p, "ardmr_plot_data", exact = TRUE)
  if (!nrow(df)) {
    if (!is.null(plot_data)) {
      attr(p, "ardmr_plot_data") <- plot_data
    }
    return(p)
  }

  df$is_separator <- df$is_separator %in% TRUE
  axis_levels <- rev(df$axis_id)
  axis_levels <- unique(axis_levels)
  axis_labels <- stats::setNames(df$axis_label, df$axis_id)
  axis_labels <- axis_labels[axis_levels]
  wrapped_labels <- stringr::str_wrap(unname(axis_labels), width = 27)
  names(wrapped_labels) <- names(axis_labels)

  p_wrapped <- p + ggplot2::scale_y_discrete(
    limits = axis_levels,
    labels = wrapped_labels
  )
  if (!is.null(plot_data)) {
    attr(p_wrapped, "ardmr_plot_data") <- plot_data
  }
  p_wrapped
}

