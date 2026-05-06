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

#' @keywords internal
.iv_chr_pos <- function(snps) {
  if (is.null(snps) || !NROW(snps)) {
    return(data.frame(chr = character(), pos = integer(),
                      ea = character(), oa = character()))
  }
  pcols <- list(c("Chr", "Pos"),
                c("panukb_chrom", "panukb_pos"),
                c("neale_chrom", "neale_pos"),
                c("chr.exposure", "pos.exposure"),
                c("chr", "pos"))
  for (pc in pcols) {
    if (all(pc %in% names(snps))) {
      ea <- if ("effect_allele.exposure" %in% names(snps)) toupper(as.character(snps$effect_allele.exposure))
            else if ("effect_allele" %in% names(snps)) toupper(as.character(snps$effect_allele))
            else rep("?", NROW(snps))
      oa <- if ("other_allele.exposure" %in% names(snps)) toupper(as.character(snps$other_allele.exposure))
            else if ("other_allele" %in% names(snps)) toupper(as.character(snps$other_allele))
            else rep("?", NROW(snps))
      return(data.frame(
        chr = as.character(snps[[pc[1]]]),
        pos = suppressWarnings(as.integer(snps[[pc[2]]])),
        ea = ea, oa = oa,
        stringsAsFactors = FALSE
      ))
    }
  }
  data.frame(chr = character(), pos = integer(),
             ea = character(), oa = character())
}

#' Hash an IV set (chr+pos+alleles) into a short stable identifier.
#'
#' Order-invariant. Returns the first `n` characters of an xxhash64 over the
#' sorted `chr:pos:ea:oa` tuples. Used to disambiguate cached fetches that
#' depend on which IVs were queried.
#'
#' @keywords internal
.iv_set_hash <- function(snps, n = 10L) {
  df <- .iv_chr_pos(snps)
  if (!nrow(df)) return(strrep("0", n))
  tuples <- sprintf("%s:%d:%s:%s", df$chr, df$pos, df$ea, df$oa)
  tuples <- sort(unique(tuples))
  substr(digest::digest(tuples, algo = "xxhash64", serialize = FALSE), 1L, n)
}

#' Hash the inputs that determine a run's *outputs* into a short identifier.
#'
#' Used to give each unique run its own subfolder under
#' `<cache>/output/<exposure>/<sex>/<ancestry>/<run_hash>/` so different
#' parameter sets don't overwrite each other's plots and CSVs.
#'
#' @keywords internal
.run_input_hash <- function(exposure_snps = NULL,
                            exposure_id = NULL,
                            exposure_sumstats = NULL,
                            sensitivity_enabled = NULL,
                            sensitivity_pass_min = NULL,
                            clump_opts = list(),
                            coloc_window_kb = NULL,
                            coloc_priors = NULL,
                            coloc_skip_mhc = NULL,
                            multiple_testing_correction = NULL,
                            n = 10L) {
  norm_set <- function(x) if (is.null(x)) NULL else sort(unique(as.character(x)))
  norm_list <- function(x) if (is.null(x) || !length(x)) NULL else x[order(names(x))]
  payload <- list(
    iv_hash = if (!is.null(exposure_snps)) .iv_set_hash(exposure_snps, n = 16L) else NULL,
    exposure_id = if (is.null(exposure_id)) NULL else as.character(exposure_id),
    exposure_sumstats = if (is.null(exposure_sumstats)) NULL else basename(as.character(exposure_sumstats)),
    sensitivity_enabled = norm_set(sensitivity_enabled),
    sensitivity_pass_min = if (is.null(sensitivity_pass_min)) NULL else as.integer(sensitivity_pass_min),
    clump_opts = norm_list(clump_opts),
    coloc_window_kb = if (is.null(coloc_window_kb)) NULL else as.integer(coloc_window_kb),
    coloc_priors = norm_list(coloc_priors),
    coloc_skip_mhc = if (is.null(coloc_skip_mhc)) NULL else isTRUE(coloc_skip_mhc),
    multiple_testing_correction = if (is.null(multiple_testing_correction)) NULL else as.character(multiple_testing_correction)
  )
  substr(digest::digest(payload, algo = "xxhash64"), 1L, n)
}

#' @keywords internal
ardmr_run_dir <- function(cache_dir, exposure_slug, sex, ancestry, run_hash) {
  file.path(cache_dir, "output", exposure_slug, sex, ancestry, run_hash)
}
