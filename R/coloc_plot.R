# R/coloc_plot.R
# Coloc-specific plot helpers: stacked locus zoom + wrappers around coloc's
# native plot_dataset() and sensitivity().

#' Stacked locus-zoom plot (exposure on top, outcome on bottom)
#'
#' @param exp_tbl tibble with `pos`, `pval`, `SNP`.
#' @param out_tbl tibble with `pos`, `pval`, `SNP`.
#' @param ivs character vector of IV rsIDs to highlight.
#' @param title plot title.
#' @return A patchwork composition (ggplot/patchwork).
#' @keywords internal
#' @export
coloc_stacked_plot <- function(exp_tbl, out_tbl, ivs = character(), title = "") {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("coloc_stacked_plot(): requires patchwork.", call. = FALSE)
  }
  mk_panel <- function(df, ylab, hl_ivs) {
    df <- df[is.finite(df$pos) & is.finite(df$pval) & df$pval > 0, , drop = FALSE]
    df$mlog10p <- -log10(df$pval)
    df$is_iv <- df$SNP %in% hl_ivs
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pos, y = .data$mlog10p)) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$is_iv), size = 0.6, alpha = 0.7) +
      ggplot2::scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#D55E00"),
                                   guide = "none") +
      ggplot2::labs(x = NULL, y = ylab) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    iv_df <- df[df$is_iv, , drop = FALSE]
    if (nrow(iv_df)) {
      p <- p + ggplot2::geom_point(data = iv_df,
                                   ggplot2::aes(x = .data$pos, y = .data$mlog10p),
                                   shape = 18, colour = "#D55E00", size = 2.4)
    }
    p
  }
  top <- mk_panel(exp_tbl, expression(-log[10] * "(p)" ~ "exposure"), ivs)
  bot <- mk_panel(out_tbl, expression(-log[10] * "(p)" ~ "outcome"), ivs) +
    ggplot2::labs(x = "position (bp)")
  if (nzchar(title)) top <- top + ggplot2::ggtitle(title)
  patchwork::wrap_plots(top, bot, ncol = 1)
}

.coloc_save_pub_plot <- function(plot, basepath, w_mm = 130, h_mm = 110, dpi = 400) {
  d <- dirname(basepath)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(basepath, ".png"), plot,
                  width = w_mm, height = h_mm, units = "mm", dpi = dpi, type = "cairo")
  ggplot2::ggsave(paste0(basepath, ".pdf"), plot,
                  width = w_mm, height = h_mm, units = "mm", device = cairo_pdf)
}

.coloc_save_base_plot <- function(fun, basepath, w_in = 5.1, h_in = 4.0) {
  d <- dirname(basepath)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(paste0(basepath, ".png"), width = w_in, height = h_in,
                 units = "in", res = 400, type = "cairo")
  ok_png <- tryCatch({ fun(); TRUE }, error = function(e) FALSE)
  grDevices::dev.off()
  if (!ok_png) try(unlink(paste0(basepath, ".png")), silent = TRUE)
  grDevices::cairo_pdf(paste0(basepath, ".pdf"), width = w_in, height = h_in)
  ok_pdf <- tryCatch({ fun(); TRUE }, error = function(e) FALSE)
  grDevices::dev.off()
  if (!ok_pdf) try(unlink(paste0(basepath, ".pdf")), silent = TRUE)
  invisible(ok_png || ok_pdf)
}
