# R/Variant_manifest_downloader.R

#' Download (if needed) the provider variant manifest (GRCh37) for RSID -> chr:pos mapping
#'
#' @description
#' Downloads and caches the large **variant manifest** required to map RSIDs to
#' genomic positions for a given catalog (`"neale"` or `"panukb"`). Files are saved under:
#'
#' - `<cache_dir>/neale/variants.tsv.bgz`
#' - `<cache_dir>/panukb/full_variant_qc_metrics.txt.bgz`
#'
#' If the file already exists, nothing is downloaded.
#'
#' Uses **curl** to show a console progress bar. Make sure you have internet access and
#' enough disk space (~hundreds of MB).
#'
#' @param catalog Character. One of `"neale"`, `"panukb"`.
#' @param data_root Character. Root directory to store data
#'   (default: ardmr_cache_dir()).
#' @param overwrite Logical. If `TRUE`, re-download even if the file exists.
#' @return (Invisibly) the local file path to the downloaded (or existing) manifest.
#' @export
Variant_manifest_downloader <- function(
    catalog,
    data_root = ardmr_cache_dir(),
    overwrite = FALSE
) {
  allowed <- c("neale", "panukb")
  if (missing(catalog) || !nzchar(catalog)) {
    stop("[catalog] must be provided ('neale' or 'panukb').", call. = FALSE)
  }
  if (!catalog %in% allowed) {
    stop(sprintf("[catalog] must be one of: %s", paste(allowed, collapse = ", ")), call. = FALSE)
  }

  # URL map only (no paths here)
  url_map <- list(
    neale  = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz",
    panukb = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz"
  )
  src_url <- url_map[[catalog]]

  # Single source of truth for destination path
  dest <- variant_manifest_path(catalog, cache_dir = data_root)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(dest) && !overwrite) {
    message(sprintf("[%s] Variant manifest already present: %s", catalog, dest))
    return(invisible(dest))
  }

  message(sprintf(
    "[%s] Downloading variant manifest to %s\n  URL: %s\n  This is a large file; please expect a 5–25 min wait.",
    catalog, dest, src_url
  ))

  # Optional size estimate
  size_str <- tryCatch(.http_content_length_human(src_url), error = function(e) NA_character_)
  if (is.character(size_str) && nzchar(size_str)) {
    message(sprintf("  Estimated size: ~%s", size_str))
  }

  # Download with progress bar (curl preferred)
  ok <- tryCatch({
    if (!requireNamespace("curl", quietly = TRUE)) {
      utils::download.file(src_url, destfile = dest, mode = "wb", method = "auto", quiet = FALSE)
    } else {
      curl::curl_download(src_url, destfile = dest, mode = "wb")
    }
    TRUE
  }, error = function(e) {
    message(sprintf("Download failed: %s", conditionMessage(e)))
    FALSE
  })

  if (!ok || !file.exists(dest)) {
    stop(sprintf("[%s] Failed to download variant manifest to %s", catalog, dest), call. = FALSE)
  }

  message(sprintf("[%s] Download complete: %s", catalog, dest))
  invisible(dest)
}

# --- helpers --------------------------------------------------------------

# Best-effort HEAD request to report Content-Length (human readable)
.http_content_length_human <- function(url) {
  if (!requireNamespace("curl", quietly = TRUE)) return(NA_character_)
  h <- curl::new_handle(nobody = TRUE, customrequest = "HEAD")
  resp <- curl::curl_fetch_memory(url, handle = h)
  if (resp$status_code >= 400) return(NA_character_)
  headers <- rawToChar(resp$headers)
  m <- regexec("Content-Length:\\s*([0-9]+)", headers, ignore.case = TRUE)
  mt <- regmatches(headers, m)[[1]]
  if (length(mt) < 2) return(NA_character_)
  bytes <- as.numeric(mt[2])
  .format_bytes(bytes)
}

.format_bytes <- function(x) {
  if (is.na(x)) return(NA_character_)
  units <- c("B","KB","MB","GB","TB")
  if (x <= 0) return(paste0(x, " B"))
  e <- floor(log(x, 1024))
  e <- max(min(e, length(units) - 1), 0)
  out <- x / (1024 ^ e)
  sprintf("%.2f %s", out, units[e + 1])
}
