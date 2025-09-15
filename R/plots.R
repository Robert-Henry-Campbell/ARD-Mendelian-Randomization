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
                           exposure = exposure_snps2$id.exposure){
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
  if (!"ARD_selected" %in% names(df)) {
    if (verbose) logger::log_warn("Manhattan: 'ARD_selected' column missing; cannot restrict to ARDs.")
    return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (ARD_selected column missing)"))
  }
  if (!is.logical(df$ARD_selected)) {
    df$ARD_selected <- df$ARD_selected %in% c(TRUE, "TRUE", "True", "true", 1, "1", "T")
  }
  n_before_ard <- nrow(df)
  df <- dplyr::filter(df, .data$ARD_selected %in% TRUE)
  if (verbose) logger::log_info("Manhattan: restricted to ARD-selected => {nrow(df)}/{n_before_ard} rows.")
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no ARD-selected rows)"))
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

