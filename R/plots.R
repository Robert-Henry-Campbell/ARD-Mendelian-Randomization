#' Manhattan plot over outcomes
#' @param results_df tidy results
#' @param Multiple_testing_correction "BH" or "bonferroni"
#' @export
manhattan_plot <- function(results_df, Multiple_testing_correction = c("BH","bonferroni")) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  # TODO: compute p_adj (if BH) on QC-pass; draw threshold/colouring
  ggplot2::ggplot(results_df, ggplot2::aes(x = seq_len(nrow(results_df)), y = -log10(p_ivw))) +
    ggplot2::geom_point()
}

#' Volcano plot over outcomes
#' @param results_df tidy results
#' @param Multiple_testing_correction "BH" or "bonferroni"
#' @export
volcano_plot <- function(results_df, Multiple_testing_correction = c("BH","bonferroni")) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  # TODO: color by GBD cause, alpha by qc_pass, etc.
  ggplot2::ggplot(results_df, ggplot2::aes(x = beta_ivw, y = -log10(p_ivw))) +
    ggplot2::geom_point()
}
