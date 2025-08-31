#' Map RSIDs to provider-specific positions and append columns
#' @param exposure_snps data.frame/tibble with RSIDs (`rsid` or `SNP`)
#' @param sex "both" → Pan-UKB; "male"/"female" → Neale
#' @param cache_dir cache root holding the downloaded variant manifest
#' @param verbose logical; log counts
#' @return exposure_snps with `panukb_chrom/pos` or `neale_chrom/pos`
#' @export
exposure_snp_mapper <- function(exposure_snps, sex, cache_dir, verbose = TRUE) {
  stopifnot(is.data.frame(exposure_snps))
  sex <- match.arg(tolower(sex), c("both","male","female"))
  catalog <- if (sex == "both") "panukb" else "neale"
  chr_col <- paste0(catalog, "_chrom")
  pos_col <- paste0(catalog, "_pos")

  exp <- tibble::as_tibble(exposure_snps)
  if (!("rsid" %in% names(exp))) {
    if ("SNP" %in% names(exp)) exp$rsid <- exp$SNP
    else stop("exposure_snps must contain 'rsid' or 'SNP'.")
  }
  exp$rsid <- .normalize_rsid(exp$rsid)
  exp <- dplyr::filter(exp, !is.na(rsid) & rsid != "")
  exp <- dplyr::distinct(exp, rsid, .keep_all = TRUE)
  n_before <- nrow(exp)

  mani_path <- variant_manifest_path(catalog = catalog, cache_dir = cache_dir)
  if (!file.exists(mani_path)) {
    stop(sprintf("Variant manifest not found at %s. Run Variant_manifest_downloader('%s') first.",
                 mani_path, catalog), call. = FALSE)
  }



  # Estimate memory: rule of thumb ~3x file size in RAM for fread
  mani_size <- file.info(mani_path)$size
  need_ram  <- max(2.5 * 1024^3, mani_size * 3)  # at least 2.5 GB or 3× file size

  # Try to query available system RAM (optional; skip if ps not present)
  free_ram <- NA_real_
  if (requireNamespace("ps", quietly = TRUE)) {
    sm <- try(ps::ps_system_memory(), silent = TRUE)
    if (!inherits(sm, "try-error") && is.list(sm) && "available" %in% names(sm)) {
      free_ram <- sm[["available"]]
    }
  }

  if (!is.na(free_ram) && free_ram < need_ram) {
    warning(sprintf(
      "Low free memory (%.1f GB available). Reading %s may fail; need ~%.1f GB.",
      free_ram/1024^3, basename(mani_path), need_ram/1024^3
    ))
  } else {
    message(sprintf(
      "Loading variant manifest (%s). File size ~%.2f GB. Expect a wait (minutes)…",
      catalog, mani_size/1024^3
    ))
  }



  hdr <- data.table::fread(mani_path, nrows = 0, showProgress = FALSE)
  cols <- names(hdr)
  pick_col <- function(cands, cols) { hit <- intersect(cands, cols); if (length(hit)) hit[[1]] else NA_character_ }
  rsid_name <- pick_col(c("rsid","RSID","variant_id","snp","SNP"), cols)
  chr_name  <- pick_col(c("chrom","chr","chromosome","CHR"), cols)
  pos_name  <- pick_col(c("pos","position","bp","POS"), cols)
  if (anyNA(c(rsid_name, chr_name, pos_name))) {
    stop(sprintf("[%s] manifest missing rsid/chrom/pos columns. Available: %s",
                 catalog, paste(cols, collapse = ", ")), call. = FALSE)
  }

  mani <- data.table::fread(
    mani_path,
    select = c(rsid_name, chr_name, pos_name),
    showProgress = FALSE
  ) |> tibble::as_tibble()
  names(mani) <- c("rsid","chrom","pos")
  mani$rsid <- .normalize_rsid(mani$rsid)
  mani <- dplyr::mutate(mani,
                        chrom = suppressWarnings(as.integer(chrom)),
                        pos   = suppressWarnings(as.integer(pos))
  )
  mani <- dplyr::distinct(mani, rsid, .keep_all = TRUE)
  mani <- dplyr::rename(mani, !!chr_col := chrom, !!pos_col := pos)

  out <- dplyr::inner_join(exp, mani, by = "rsid")
  n_after <- nrow(out)
  n_missing <- n_before - n_after

  if (verbose) {
    if (!requireNamespace("logger", quietly = TRUE)) {
      message(sprintf("Exposure mapping (%s): %d mapped, %d not found", catalog, n_after, n_missing))
    } else {
      logger::log_info("Exposure mapping ({catalog}): {n_after} mapped; {n_missing} not found")
    }
  }
  out
}

.normalize_rsid <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  add_rs <- !is.na(x) & grepl("^\\d+$", x)
  x[add_rs] <- paste0("rs", x[add_rs])
  x <- sub("^RS", "rs", x)
  x
}
