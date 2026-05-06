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

#' Render coloc::sensitivity() with proper trait names and wider margins
#'
#' Wraps `coloc::sensitivity()` so that:
#'   * the `"trait 1"` / `"trait 2"` panel titles are replaced with the actual
#'     exposure / outcome labels (text-wrapped at `wrap` characters), and
#'   * the y-axis labels in the prior/posterior probability panels stop
#'     overlapping the legend (achieved with a wider PNG/PDF and a slightly
#'     larger left margin).
#'
#' Implementation note: the installed `coloc` (v5.x) hardcodes the panel
#' titles and the `par(mar = c(3, 3, 2, 1))` left-margin via base graphics.
#' To override without forking the package we briefly `assignInNamespace()`
#' a wrapper around `graphics::title` (translating `"trait 1"` /
#' `"trait 2"` to the user's labels) and `graphics::par` (substituting a
#' larger left margin only when the exact `mar = c(3, 3, 2, 1)` call comes
#' through). Both are restored on exit.
#'
#' @param abf_res Output of `coloc::coloc.abf()`.
#' @param basepath File-path stem (no extension); `.png` and `.pdf` are
#'   written.
#' @param exp_label,out_label Character. Names of the exposure and outcome
#'   traits. Long values are wrapped at `wrap` characters across multiple
#'   lines. NULL leaves the default `"trait 1"` / `"trait 2"`.
#' @param wrap Integer. Wrap width in characters for trait labels.
#' @param rule Sensitivity rule passed straight to `coloc::sensitivity()`.
#' @param w_in,h_in Plot dimensions in inches.
#' @param left_mar Substituted `par(mar)[2]` value (in lines) to give the
#'   prior/posterior y-axis labels breathing room.
#' @return Invisibly returns the `basepath` stem.
#' @keywords internal
#' @export
coloc_plot_sensitivity <- function(abf_res, basepath,
                                   exp_label = NULL, out_label = NULL,
                                   wrap = 30, rule = "H4 > 0.8",
                                   w_in = 8.5, h_in = 6.0,
                                   left_mar = 4.4) {
  if (!requireNamespace("coloc", quietly = TRUE)) {
    stop("coloc_plot_sensitivity(): requires the coloc package.", call. = FALSE)
  }
  wrap_one <- function(s) {
    if (is.null(s)) return(NULL)
    s <- as.character(s)
    if (length(s) == 0L || is.na(s[1]) || !nzchar(s[1])) return(NULL)
    paste(strwrap(s[1], width = wrap), collapse = "\n")
  }
  exp_t <- wrap_one(exp_label)
  out_t <- wrap_one(out_label)

  # Build a "shadow" environment that intercepts title() and par() calls
  # made inside coloc::sensitivity()'s body, then attach it as that
  # function's environment. Lexical scoping means coloc::sensitivity will
  # resolve title/par against this shadow first (giving us our overrides),
  # while everything else still resolves up the chain to the original
  # coloc namespace. assignInNamespace() can't touch base packages, so this
  # is the cleanest non-fork way to override the hardcoded "trait 1" /
  # "trait 2" titles and the cramped par(mar = c(3,3,2,1)) left-margin.
  shadow <- new.env(parent = environment(coloc::sensitivity))
  shadow$title <- function(main = NULL, ...) {
    if (!is.null(main) && is.character(main) && length(main) >= 1L) {
      m <- as.character(main)[1]
      if (startsWith(m, "trait 1") && !is.null(exp_t)) main <- exp_t
      else if (startsWith(m, "trait 2") && !is.null(out_t)) main <- out_t
    }
    graphics::title(main = main, ...)
  }
  shadow$par <- function(...) {
    args <- list(...)
    nms  <- names(args)
    if (length(nms) && "mar" %in% nms &&
        length(args$mar) == 4L &&
        isTRUE(all(suppressWarnings(as.numeric(args$mar)) == c(3, 3, 2, 1)))) {
      args$mar <- c(3, left_mar, 2, 1)
    }
    do.call(graphics::par, args)
  }
  shadowed_sensitivity <- coloc::sensitivity
  environment(shadowed_sensitivity) <- shadow

  call_fn <- function() shadowed_sensitivity(abf_res, rule = rule)
  .coloc_save_base_plot(call_fn, basepath, w_in = w_in, h_in = h_in)
  invisible(basepath)
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
