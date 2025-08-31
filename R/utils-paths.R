# R/utils-paths.R
#' @keywords internal
ardmr_cache_dir <- function() {
  tools::R_user_dir("ardmr", which = "cache")
}

#' @keywords internal
variant_manifest_path <- function(catalog, cache_dir = ardmr_cache_dir()) {
  stopifnot(catalog %in% c("neale","panukb"))
  file.path(
    cache_dir, catalog,
    if (catalog == "neale") "variants.tsv.bgz" else "full_variant_qc_metrics.txt.bgz"
  )
}
