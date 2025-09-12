# R/plots.R (UNTESTED)

#' Manhattan plot over outcomes (IVW p-values on QC-pass outcomes)
#' @param results_df tidy results (must contain results_p_ivw; ideally results_qc_pass)
#' @param Multiple_testing_correction "BH" or "bonferroni"
#' @param alpha significance level for corrections (default 0.05)
#' @param verbose if TRUE, emit step-by-step logs via {logger}
#' @export
# helpers (keep these once at top-level)
col_or_null <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NULL
`%||%`      <- function(a, b) if (!is.null(a)) a else b

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

  # Alias AFTER df exists
  if (!"gbd_cause" %in% names(df) && "Cause Name" %in% names(df)) {
    df$gbd_cause <- df[["Cause Name"]]
  }

  # QC-pass
  if (!"results_qc_pass" %in% names(df)) {
    if (verbose) logger::log_info("Manhattan: 'results_qc_pass' column not present; assuming all pass.")
    df$results_qc_pass <- TRUE
  }
  n_before <- nrow(df)
  df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  n_after  <- nrow(df)
  if (verbose) logger::log_info("Manhattan: filtered to QC-pass => {n_after}/{n_before} rows remain.")

  # Required p
  if (!nrow(df) || !"results_p_ivw" %in% names(df)) {
    if (!"results_p_ivw" %in% names(df) && verbose) {
      logger::log_warn("Manhattan: required column 'results_p_ivw' missing after filtering; returning placeholder plot.")
    } else if (verbose) {
      logger::log_warn("Manhattan: no QC-pass rows; returning placeholder plot.")
    }
    return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no QC-pass results)"))
  }

  # Compute -log10(p)
  df <- dplyr::filter(df, !is.na(.data$results_p_ivw))
  df$logp <- -log10(pmax(df$results_p_ivw, .Machine$double.xmin))
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::labs(title = "Manhattan (no finite p-values)"))

  # ==============================
  # GROUPING & X POSITIONS (L1→L2)
  # ==============================
  # Use L3 only for colouring if you ever want; for layout we stop at L2 as requested.
  # Arrange deterministically: L1 → L2 → outcome label
  df <- dplyr::arrange(
    df,
    dplyr::coalesce(.data$cause_level_1, "~"),
    dplyr::coalesce(.data$cause_level_2, "~"),
    dplyr::coalesce(.data$results_outcome, "")
  )
  df$idx <- seq_len(nrow(df))

  # Build block maps for L2 (within L1) and for L1
  l2_map <- df |>
    dplyr::group_by(.data$cause_level_1, .data$cause_level_2) |>
    dplyr::summarise(start = min(.data$idx), end = max(.data$idx), .groups = "drop") |>
    dplyr::arrange(.data$cause_level_1, .data$start) |>
    dplyr::mutate(center = (start + end)/2,
                  alt_id = dplyr::row_number()) # for alternating colours

  l1_map <- df |>
    dplyr::group_by(.data$cause_level_1) |>
    dplyr::summarise(start = min(.data$idx), end = max(.data$idx), .groups = "drop") |>
    dplyr::arrange(.data$start) |>
    dplyr::mutate(center = (start + end)/2)

  # Alternate black/grey by L2 blocks
  df <- df |>
    dplyr::left_join(l2_map |> dplyr::select(.data$cause_level_1, .data$cause_level_2, .data$alt_id),
                     by = c("cause_level_1","cause_level_2")) |>
    dplyr::mutate(alt_col = ifelse(.data$alt_id %% 2L == 1L, "grey20", "grey65"))

  # =======================
  # MULTIPLE TESTING LINES
  # =======================
  if (Multiple_testing_correction == "BH") {
    df$q_bh <- stats::p.adjust(df$results_p_ivw, method = "BH")
    df$sig  <- df$q_bh < alpha
    thr_y   <- NA_real_
    if (verbose) logger::log_info("Manhattan: BH done; {sum(df$sig, na.rm=TRUE)} hits at q < {alpha}.")
  } else {
    m       <- nrow(df)
    alpha_b <- alpha / max(1L, m)
    df$sig  <- df$results_p_ivw < alpha_b
    thr_y   <- -log10(alpha_b)
    if (verbose) logger::log_info("Manhattan: Bonferroni m={m}; ythr={round(thr_y,3)}; hits={sum(df$sig, na.rm=TRUE)}.")
  }

  # ============================
  # BRACKET (“TREE”) ANNOTATIONS
  # ============================
  # We'll allocate a little negative y-space to draw brackets and labels.
  y_max <- max(df$logp, na.rm = TRUE)
  y_l2  <- -0.8   # L2 bracket baseline
  y_l1  <- -1.6   # L1 bracket baseline
  cap   <- 0.35   # bracket end "teeth" height
  lab_gap <- 0.25 # label offset below the bracket

  # Build segments for horizontal bars and vertical "teeth"
  seg_l2 <- l2_map |>
    dplyr::transmute(x = start, xend = end, y = y_l2, yend = y_l2)
  teeth_l2 <- dplyr::bind_rows(
    l2_map |> dplyr::transmute(x = start, xend = start, y = y_l2, yend = y_l2 - cap),
    l2_map |> dplyr::transmute(x = end,   xend = end,   y = y_l2, yend = y_l2 - cap)
  )
  lab_l2 <- l2_map |>
    dplyr::transmute(x = center, y = y_l2 - lab_gap,
                     label = dplyr::coalesce(cause_level_2, "NA"))

  seg_l1 <- l1_map |>
    dplyr::transmute(x = start, xend = end, y = y_l1, yend = y_l1)
  teeth_l1 <- dplyr::bind_rows(
    l1_map |> dplyr::transmute(x = start, xend = start, y = y_l1, yend = y_l1 - cap),
    l1_map |> dplyr::transmute(x = end,   xend = end,   y = y_l1, yend = y_l1 - cap)
  )
  lab_l1 <- l1_map |>
    dplyr::transmute(x = center, y = y_l1 - lab_gap,
                     label = dplyr::coalesce(cause_level_1, "NA"))

  # =================
  # BUILD THE PLOT
  # =================
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$logp)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$alt_col), alpha = 0.9) +
    ggplot2::scale_colour_identity(guide = "none") +
    { if (!is.na(thr_y)) ggplot2::geom_hline(yintercept = thr_y, linetype = "dashed") } +
    ggplot2::labs(
      x = NULL,
      y = expression(-log[10](p[IVW])),
      title = "Manhattan of MR IVW results (grouped L1→L2)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 60, l = 10) # extra bottom room for brackets
    ) +
    # Brackets & labels for L2 (near axis)
    ggplot2::geom_segment(data = seg_l2, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    ggplot2::geom_segment(data = teeth_l2, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    ggplot2::geom_text(data = lab_l2, ggplot2::aes(x = x, y = y, label = label), vjust = 1, size = 3) +
    # Brackets & labels for L1 (below)
    ggplot2::geom_segment(data = seg_l1, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.5) +
    ggplot2::geom_segment(data = teeth_l1, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.5) +
    ggplot2::geom_text(data = lab_l1, ggplot2::aes(x = x, y = y, label = label), vjust = 1, fontface = "bold", size = 3.3) +
    # Make room below 0 for the “brackets”
    ggplot2::coord_cartesian(ylim = c(min(y_l1 - 0.8, 0), y_max * 1.05), clip = "off")

  if (verbose) logger::log_info("Manhattan: plot object constructed; returning ggplot.")
  p
}


#' Volcano plot over outcomes
#' @param results_df tidy results (needs results_beta_ivw, results_se_ivw, results_p_ivw; ideally results_qc_pass, results_nsnp_after)
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
  need <- c("results_beta_ivw","results_se_ivw","results_p_ivw")
  if (!all(need %in% names(df))) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "Volcano (missing beta/se/p)"))
  }

  # computed axes
  df$z_ivw   <- df$results_beta_ivw / df$results_se_ivw
  df$logp    <- -log10(pmax(df$results_p_ivw, .Machine$double.xmin))
  if (!"results_qc_pass" %in% names(df)) df$results_qc_pass <- TRUE
  if (!"results_nsnp_after" %in% names(df)) df$results_nsnp_after <- NA_integer_

  # group (colour)
  group_col <- dplyr::coalesce(
    df$cause_level_3 %||% NULL,
    df$cause_level_2 %||% NULL,
    df$`Cause Name`     %||% NULL
  )
  if (is.null(group_col)) group_col <- factor("All")
  df$group <- as.factor(group_col)

  # flags for outline ring
  het_fail   <- (!is.na(df$results_Q_p_ivw) & df$results_Q_p_ivw < 0.05) | (!is.na(df$results_I2_ivw) & df$results_I2_ivw > 50)
  egger_fail <- !is.na(df$results_egger_intercept_p) & df$results_egger_intercept_p < 0.05
  df$flag_any <- het_fail | egger_fail

  # multiple-testing significance
  if (Multiple_testing_correction == "BH") {
    df$q_bh <- stats::p.adjust(df$results_p_ivw, method = "BH")
    df$sig  <- df$q_bh < alpha
    thr_y   <- NA_real_
  } else {
    m       <- sum(!is.na(df$results_p_ivw))
    alpha_b <- alpha / max(1L, m)
    df$sig  <- df$results_p_ivw < alpha_b
    thr_y   <- -log10(alpha_b)
  }

  # base plot: shape 21 so we can have outline ring for flagged points
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$z_ivw, y = .data$logp,
                 fill = .data$group,
                 alpha = .data$results_qc_pass,
                 size = .data$results_nsnp_after)
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
    lab_df <- dplyr::mutate(df, label = dplyr::coalesce(.data$results_outcome, NA_character_))
    lab_df <- dplyr::filter(lab_df, .data$sig %in% TRUE, !is.na(.data$label))

    if (nrow(lab_df)) {
      lab_df <- dplyr::arrange(lab_df, dplyr::desc(abs(.data$z_ivw)))
      n_lab  <- min(label_top_n, nrow(lab_df))
      lab_df <- dplyr::slice_head(lab_df, n = n_lab)  # n must be a constant integer

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
