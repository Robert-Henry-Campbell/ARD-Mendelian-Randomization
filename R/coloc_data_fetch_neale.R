# R/coloc_data_fetch_neale.R
# Regional summary-statistic fetcher for Neale Lab UK Biobank sumstats.
# Mirrors the panukb pair in R/coloc_data_fetch.R at the interface level
# (same factory contract, same return shape) but reads local .tsv.bgz files
# (no tabix; chr/pos must be parsed from the single `variant` column).

#' Read and cache a per-chromosome slice of a Neale .tsv.bgz file.
#'
#' Neale sumstats have ~13.8M rows and no chr/pos columns; the chr/pos info
#' is encoded in the single `variant = chr:pos:ref:alt` column. We must read
#' the whole file once to extract any chromosome's rows. This helper writes
#' just the requested chromosome's parsed slice to disk and returns it; later
#' calls for the same (file, chr) pair hit cache. A request for a different
#' chromosome of the same file pays the read cost again.
#'
#' @param file_path Absolute path to a Neale `.tsv.bgz` file.
#' @param chr Chromosome to extract (character; e.g. `"6"` or `"X"`).
#' @param cache_dir Cache root.
#' @param verbose Logical.
#' @param force_refresh Logical; if TRUE delete the chr's cache before re-reading.
#' @return Tibble with `v_chr, v_pos, ref, alt, minor_allele, minor_AF, low_confidence_variant, beta, se, pval`.
#' @keywords internal
.coloc_neale_chr_slice <- function(file_path, chr,
                                   cache_dir = ardmr_cache_dir(),
                                   verbose = TRUE,
                                   force_refresh = FALSE) {
  chr <- as.character(chr)
  key <- .coloc_region_key("neale_chr_slice", basename(file_path), chr)
  cache_path <- .coloc_cache_path(cache_dir, "neale_chr_slices", key)

  if (isTRUE(force_refresh) && file.exists(cache_path)) try(unlink(cache_path), silent = TRUE)
  if (file.exists(cache_path)) {
    if (verbose) logger::log_info("coloc[neale]: cached slice {basename(file_path)} chr{chr}")
    return(readRDS(cache_path))
  }

  if (!file.exists(file_path)) {
    logger::log_warn("coloc[neale]: file missing: {file_path}; returning empty slice")
    return(tibble::tibble())
  }

  if (verbose) {
    logger::log_info(
      "coloc[neale]: parsing {basename(file_path)} for chr{chr} (~30-60s for ~13.8M rows)"
    )
  }
  t0 <- Sys.time()
  dt <- tryCatch(
    data.table::fread(
      file_path,
      select = c("variant", "minor_allele", "minor_AF",
                 "low_confidence_variant", "beta", "se", "pval"),
      colClasses = c(low_confidence_variant = "character"),
      showProgress = FALSE,
      data.table = TRUE
    ),
    error = function(e) { logger::log_warn("coloc[neale]: fread failed for {basename(file_path)}: {conditionMessage(e)}"); NULL }
  )
  if (is.null(dt) || !nrow(dt)) return(tibble::tibble())

  parts <- data.table::tstrsplit(dt$variant, ":", fixed = TRUE, type.convert = FALSE)
  if (length(parts) != 4L) {
    logger::log_warn("coloc[neale]: unexpected variant format in {basename(file_path)}; returning empty slice")
    return(tibble::tibble())
  }
  dt$v_chr <- parts[[1]]
  dt$v_pos <- parts[[2]]
  dt$ref   <- parts[[3]]
  dt$alt   <- parts[[4]]
  dt$variant <- NULL

  keep <- dt$v_chr == chr
  slice <- dt[keep, , drop = FALSE]
  rm(dt, parts, keep); invisible(gc(verbose = FALSE))
  slice$v_pos    <- suppressWarnings(as.integer(slice$v_pos))
  slice$minor_AF <- suppressWarnings(as.numeric(slice$minor_AF))
  slice$beta     <- suppressWarnings(as.numeric(slice$beta))
  slice$se       <- suppressWarnings(as.numeric(slice$se))
  slice$pval     <- suppressWarnings(as.numeric(slice$pval))
  slice <- tibble::as_tibble(slice)

  tmp <- paste0(cache_path, ".tmp")
  saveRDS(slice, tmp)
  file.rename(tmp, cache_path)

  if (verbose) {
    logger::log_info("coloc[neale]: slice {basename(file_path)} chr{chr} -> {nrow(slice)} rows in {round(as.numeric(Sys.time() - t0, units='secs'), 1)}s")
  }
  slice
}

#' Fetch outcome regional summary statistics from a Neale Lab .tsv.bgz file
#'
#' @param rec One-row tibble/data.frame from MR_df, must hold a `File` column.
#' @param sex Outcome sex stratum (`"male"` or `"female"`); informational here
#'   (file selection happens upstream in `Outcome_setup()`).
#' @param chr,start,end region coords (chr coerced to character).
#' @param neale_dir Directory containing Neale `.tsv.bgz` files.
#' @param cache_dir Cache root.
#' @param verbose Logical.
#' @param force_refresh Logical; if TRUE invalidate the per-chr slice cache.
#' @return Tibble with `SNP, chr, pos, effect_allele, other_allele, beta, se, eaf, pval`.
#' @keywords internal
#' @export
coloc_fetch_outcome_neale_region <- function(rec, sex, chr, start, end,
                                             neale_dir,
                                             cache_dir = ardmr_cache_dir(),
                                             verbose = TRUE,
                                             force_refresh = FALSE) {
  empty <- tibble::tibble(SNP = character(), chr = character(), pos = integer(),
                          effect_allele = character(), other_allele = character(),
                          beta = double(), se = double(), eaf = double(), pval = double())
  if (is.null(rec) || !("File" %in% names(rec))) {
    logger::log_warn("coloc[neale]: rec missing 'File' column; skipping")
    return(empty)
  }
  fname <- as.character(rec$File[[1]])
  if (!nzchar(fname) || is.na(fname)) {
    logger::log_warn("coloc[neale]: empty File field; skipping")
    return(empty)
  }
  file_path <- file.path(neale_dir, basename(fname))
  if (!file.exists(file_path)) {
    logger::log_warn("coloc[neale]: file not found on disk: {file_path}; skipping")
    return(empty)
  }

  chr_chr <- as.character(chr)
  start_i <- as.integer(start)
  end_i   <- as.integer(end)

  slice <- .coloc_neale_chr_slice(file_path, chr_chr, cache_dir = cache_dir,
                                  verbose = verbose, force_refresh = force_refresh)
  if (!nrow(slice)) return(empty)

  win <- slice[!is.na(slice$v_pos) & slice$v_pos >= start_i & slice$v_pos <= end_i, , drop = FALSE]
  if (!nrow(win)) return(empty)

  win <- win[is.na(win$low_confidence_variant) | win$low_confidence_variant != "true", , drop = FALSE]
  if (!nrow(win)) return(empty)

  ref <- toupper(as.character(win$ref))
  alt <- toupper(as.character(win$alt))
  ma  <- toupper(as.character(win$minor_allele))
  eaf <- ifelse(ma == alt, win$minor_AF, 1 - win$minor_AF)

  out <- tibble::tibble(
    chr           = win$v_chr,
    pos           = win$v_pos,
    effect_allele = alt,
    other_allele  = ref,
    beta          = win$beta,
    se            = win$se,
    eaf           = eaf,
    pval          = win$pval
  )
  out$SNP <- paste0(out$chr, ":", out$pos, "_", out$other_allele, "_", out$effect_allele)
  out <- dplyr::relocate(out, "SNP", .before = 1)
  out <- out[!is.na(out$beta) & !is.na(out$se) & out$se > 0 & !is.na(out$pval), , drop = FALSE]
  out
}

#' Derive Neale outcome metadata (N, trait type, case/control) from MR_df row.
#'
#' Mirrors `coloc_panukb_outcome_metadata()`. Pulls from columns merged into
#' `MR_df` by `Outcome_setup()` from the Neale per-sex phenotype manifest
#' (`neale_male_manifest` / `neale_female_manifest`).
#'
#' @param rec One-row tibble/data.frame from MR_df.
#' @param sex Outcome sex stratum (informational).
#' @return list(N, type, ncase, ncontrol, s).
#' @keywords internal
#' @export
coloc_neale_outcome_metadata <- function(rec, sex) {
  pull <- function(nms) {
    for (n in nms) if (n %in% names(rec)) {
      v <- suppressWarnings(as.numeric(rec[[n]][[1]]))
      if (!is.na(v)) return(v)
    }
    NA_real_
  }
  ncase    <- pull(c("n_cases"))
  ncontrol <- pull(c("n_controls"))
  N        <- pull(c("n_non_missing"))
  vt <- if ("variable_type" %in% names(rec)) tolower(as.character(rec$variable_type[[1]])) else NA_character_

  is_quant <- !is.na(vt) && vt %in% c("continuous_irnt", "continuous_raw", "ordinal")
  if (is_quant) {
    if (is.na(N) && !is.na(ncase) && !is.na(ncontrol)) N <- ncase + ncontrol
    type <- "quant"
    ncase <- NA_real_; ncontrol <- NA_real_; s <- NA_real_
  } else {
    if (is.na(N) && !is.na(ncase) && !is.na(ncontrol)) N <- ncase + ncontrol
    type <- if (!is.na(ncase) && ncase > 0) "cc" else "quant"
    s <- if (!is.na(ncase) && !is.na(ncontrol) && (ncase + ncontrol) > 0) ncase / (ncase + ncontrol) else NA_real_
  }
  list(N = N, type = type, ncase = ncase, ncontrol = ncontrol, s = s)
}
