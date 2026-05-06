# R/utils-snp.R
# Small SNP-shaped helpers shared across the package. Lifted out of
# run_ieugwasr_ard_compare() closures so preprocess_exposure_snps() and
# .resolve_coloc_source() can reuse them without duplication.

#' @keywords internal
is_rsid <- function(x) !is.na(x) & grepl("^rs\\d+$", x)

#' @keywords internal
palindrome_flag <- function(a1, a2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

#' @keywords internal
# Wald-chi-squared statistic; named f_stat for community convention.
f_stat <- function(beta, se) (beta / se) ^ 2

#' @keywords internal
# Robust ieugwasr::tophits across versions (p vs pval arg). Returns NULL on
# error; otherwise a data.frame with a `p` column (synthesised from `pval`
# when missing) so callers don't have to branch.
safe_tophits <- function(id, thr) {
  th <- try(ieugwasr::tophits(id = id, p = thr), silent = TRUE)
  if (inherits(th, "try-error") || !is.data.frame(th)) {
    th <- try(ieugwasr::tophits(id = id, pval = thr), silent = TRUE)
  }
  if (inherits(th, "try-error") || !is.data.frame(th)) return(NULL)
  if (!nrow(th)) return(th[0, , drop = FALSE])
  if (!"p" %in% names(th) && "pval" %in% names(th)) th$p <- th$pval
  th
}
