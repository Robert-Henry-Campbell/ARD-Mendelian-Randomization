#' Map RSIDs to provider-specific positions and append columns
#' @param exposure_snps data frame with rsid
#' @param sex "both" -> panukb mapping; otherwise neale mapping
#' @param cache_dir path for crosswalk caches
#' @param verbose logical
#' @return exposure_snps with new col "panukb_pos" or "neale_pos"
#' @export
exposure_snp_mapper <- function(exposure_snps, sex, cache_dir, verbose = TRUE) {
  n_before <- nrow(exposure_snps)
  # TODO: load variant manifest (provider) and map rsid -> chr:pos (GRCh37)
  # Return unchanged for now (stub)
  out <- exposure_snps
  n_after <- nrow(out)
  logger::log_info("Exposure mapping: {n_before - n_after} filtered; {n_after} remain")
  out
}
