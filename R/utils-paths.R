# R/utils-paths.R
#' Resolve the ARDMR cache directory from the environment
#'
#' Looks for the `ARDMR_CACHE_DIR` environment variable and returns its
#' value. This must point to a writable directory used for caching large
#' files (variant manifests, etc.).
#'
#' @return String path to the cache directory.
#' @keywords internal
ardmr_cache_dir <- function() {
  cache <- Sys.getenv("ARDMR_CACHE_DIR", unset = "")
  if (!nzchar(cache)) {
    stop(
      paste(
        "Environment variable ARDMR_CACHE_DIR is not set.",
        "Set it via set_cache_dir('/path/to/cache') or",
        "Sys.setenv(ARDMR_CACHE_DIR = '/path/to/cache')."
      ),
      call. = FALSE
    )
  }
  cache
}

#' Set the ARDMR cache directory
#'
#' Convenience helper to set `ARDMR_CACHE_DIR` for the current R session.
#'
#' @param path Existing directory to use for caching.
#' @return Invisibly returns `path`.
#' @export
set_cache_dir <- function(path) {
  if (missing(path) || !nzchar(path)) {
    stop("`path` must be a non-empty directory path.", call. = FALSE)
  }
  if (!dir.exists(path)) {
    stop(sprintf("Cache directory does not exist: %s", path), call. = FALSE)
  }
  Sys.setenv(ARDMR_CACHE_DIR = path)
  invisible(path)
}

#' @keywords internal
variant_manifest_path <- function(catalog, cache_dir = ardmr_cache_dir()) {
  stopifnot(catalog %in% c("neale","panukb"))
  file.path(
    cache_dir, catalog,
    if (catalog == "neale") "variants.tsv.bgz" else "full_variant_qc_metrics.txt.bgz"
  )
}
