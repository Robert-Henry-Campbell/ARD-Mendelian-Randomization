#' Map RSIDs to provider-specific positions and append columns
#' @param exposure_snps data frame with rsid
#' @param sex "both" -> panukb mapping; otherwise neale mapping
#' @param cache_dir path for crosswalk caches
#' @param verbose logical
#' @return exposure_snps with new col "panukb_pos" or "neale_pos"
#' @export
exposure_snp_mapper <- function(exposure_snps, sex, cache_dir, verbose = TRUE) {
  # TODO: load variant manifest (provider) and map rsid -> chr:pos (GRCh37)
  # Return unchanged for now (stub)
  exposure_snps
}
