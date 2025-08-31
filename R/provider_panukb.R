#' Fetch outcome SNP rows from Pan-UKB sumstats for each outcome
#' @param exposure_snps mapped snps (must include panukb_pos when implemented)
#' @param MR_df tibble with provider info (panukb rows)
#' @param ancestry character
#' @param cache_dir path
#' @param verbose logical
#' @return MR_df with a new list-column `outcome_snps` (TwoSampleMR-formatted)
#' @export
panukb_snp_grabber <- function(exposure_snps, MR_df, ancestry, cache_dir, verbose = TRUE) {
  # TODO: iterate MR_df[provider=="panukb"], pull only needed rows per phenotype
  MR_df$outcome_snps <- vector("list", nrow(MR_df))
  n_snps <- sum(lengths(MR_df$outcome_snps))
  logger::log_info("Pan-UKB SNP grabber: {nrow(MR_df)} outcomes; {n_snps} SNP rows")
  MR_df
}
