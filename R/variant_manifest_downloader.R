# R/Variant_manifest_downloader.R UNTESTED

#' Download (if needed) the provider variant manifest (GRCh37) for RSID -> chr:pos mapping
#'
#' @description
#' Downloads and caches the large **variant manifest** required to map RSIDs to
#' genomic positions for a given catalog (`"neale"` or `"panukb"`). Files are saved under:
#'
#' - `data/neale/variants.tsv.bgz`
#' - `data/panukb/full_variant_qc_metrics.txt.bgz`
#'
#' If the file already exists, nothing is downloaded.
#'
#' Uses **curl** to show a console progress bar. Make sure you have internet access and
#' enough disk space (~hundreds of MB).
#'
#' @param catalog Character. One of `"neale"`, `"panukb"`.
#' @param data_root Character. Root directory to store data (default `"data"`).
#' @param overwrite Logical. If `TRUE`, re-download even if the file exists.
#' @return (Invisibly) the local file path to the downloaded (or existing) manifest.
#' @export
Variant_manifest_downloader <- function(catalog, data_root = "data", overwrite = FALSE) {
  allowed <- c("neale", "panukb")
  if (missing(catalog) || !nzchar(catalog)) {
    stop("[catalog] must be provided ('neale' or 'panukb').", call. = FALSE)
  }
  if (!catalog %in% allowed) {
    stop(sprintf("[catalog] must be one of: %s", paste(allowed, collapse = ", ")), call. = FALSE)
  }

  # URLs & destinations
  cfg <- list(
    neale = list(
      url  = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz",
      dir  = file.path(data_root, "neale"),
      file = "variants.tsv.bgz"
    ),
    panukb = list(
      url  = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz",
      dir  = file.path(data_root, "panukb"),
      file = "full_variant_qc_metrics.txt.bgz"
    )
  )[[catalog]]

  dir.create(cfg$dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(cfg$dir, cfg$file)

  if (file.exists(dest) && !overwrite) {
    message(sprintf("[%s] Variant manifest already present: %s", catalog, dest))
    return(invisible(dest))
  }

  # Inform the user this can be large/slow
  message(sprintf(
    "[%s] Downloading variant manifest to %s\n  URL: %s\n  This is a large file; please expect a long wait. Showing progress bar...",
    catalog, dest, cfg$url
  ))

  # Try to fetch content length for a size estimate (optional best-effort)
  size_str <- tryCatch(.http_content_length_human(cfg$url), error = function(e) NA_character_)
  if (is.character(size_str) && nzchar(size_str)) {
    message(sprintf("  Estimated size: ~%s", size_str))
  }

  # Download with progress bar (curl)
  ok <- tryCatch({
    if (!requireNamespace("curl", quietly = TRUE)) {
      utils::download.file(cfg$url, destfile = dest, mode = "wb", method = "auto", quiet = FALSE)
    } else {
      curl::curl_download(cfg$url, destfile = dest, mode = "wb")  # shows progress bar
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
  # Parse Content-Length
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
