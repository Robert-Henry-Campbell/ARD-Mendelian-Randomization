#' Ensure Neale sumstats exist locally (download if needed)
#' @param MR_df outcomes plan (neale rows)
#' @param neale_dir directory with .bgz sumstats
#' @param verbose logical
#' @export
neale_gwas_checker <- function(MR_df, neale_dir, verbose = TRUE) {
  # TODO: check presence per outcome; optionally download; log actions
  invisible(TRUE)
}

#' Generate .tbi for all bgz files (no-op if present)
#' @param neale_dir directory
#' @param verbose logical
#' @export
neale_tbi_maker <- function(neale_dir, verbose = TRUE) {
  # TODO: run tabix indexing if missing
  invisible(TRUE)
}

#' Fetch outcome SNP rows from Neale sumstats for each outcome
#' @param exposure_snps mapped snps (must include neale_pos when implemented)
#' @param MR_df tibble (neale rows)
#' @param neale_dir directory with sumstats
#' @param cache_dir path
#' @param verbose logical
#' @return MR_df with list-column `outcome_snps`
#' @export
neale_snp_grabber <- function(exposure_snps, MR_df, neale_dir, cache_dir, verbose = TRUE) {
  # TODO: iterate MR_df[provider=="neale"], extract SNP rows; attach rsid; format for TwoSampleMR
  MR_df$outcome_snps <- vector("list", nrow(MR_df))
  MR_df
}
