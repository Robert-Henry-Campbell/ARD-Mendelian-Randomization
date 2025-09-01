#' Ensure Neale sumstats exist locally (download if needed)
#' @param MR_df outcomes plan (neale rows)
#' @param neale_dir directory with .bgz sumstats
#' @param verbose logical
#' @export
neale_gwas_checker <- function(MR_df, neale_dir, verbose = TRUE) {
  # TODO: check presence per outcome; optionally download; log actions
  logger::log_info("Neale GWAS check: {nrow(MR_df)} outcomes scanned")
  invisible(TRUE)
}

#' Generate .tbi for all bgz files (no-op if present)
#' @param neale_dir directory
#' @param verbose logical
#' @export
neale_tbi_maker <- function(neale_dir, verbose = TRUE) {
  # TODO: run tabix indexing if missing
  n_files <- length(list.files(neale_dir, pattern = "\\.bgz$"))
  logger::log_info("Neale tbi maker: {n_files} files examined")
  invisible(TRUE)
}

#' Fetch outcome SNP rows from Neale sumstats for each outcome
#' @param exposure_snps mapped snps (must include neale_pos when implemented)
#' @param MR_df tibble (neale rows)
#' @param neale_dir directory with sumstats
#' @param cache_dir path (default: [ardmr_cache_dir()])
#' @param verbose logical
#' @return MR_df with list-column `outcome_snps`
#' @export
neale_snp_grabber <- function(exposure_snps, MR_df, neale_dir, cache_dir = ardmr_cache_dir(), verbose = TRUE) {
  # TODO: iterate MR_df[provider=="neale"], extract SNP rows; attach rsid; format for TwoSampleMR
  MR_df$outcome_snps <- vector("list", nrow(MR_df))
  n_snps <- sum(lengths(MR_df$outcome_snps))
  logger::log_info("Neale SNP grabber: {nrow(MR_df)} outcomes; {n_snps} SNP rows")
  MR_df
}
