# R/panukb_snp_grabber.R

#' Fetch outcome SNP rows from Pan-UKB sumstats for each outcome (R-based Tabix)
#'
#' Uses Rsamtools/Rhtslib to remotely query tabix-indexed BGZF files over HTTP(S),
#' removing the need for a system 'tabix' in PATH.
#'
#' @param exposure_snps mapped SNPs (must include `panukb_chrom` and `panukb_pos`)
#' @param MR_df tibble with provider info (Pan-UKB rows and a column containing a sumstats URL or wget cmd)
#' @param ancestry character (e.g. "EUR")
#' @param cache_dir (unused here but kept for parity; default:
#'   [ardmr_cache_dir()])
#' @param verbose logical
#' @return MR_df with a new list-column `outcome_snps` (TwoSampleMR-formatted)
#' @export
#' @importFrom Rsamtools TabixFile headerTabix scanTabix open.TabixFile close.TabixFile
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
panukb_snp_grabber <- function(exposure_snps, MR_df, ancestry, cache_dir = ardmr_cache_dir(), verbose = TRUE) {
  stopifnot(is.data.frame(exposure_snps), is.data.frame(MR_df))
  ancestry <- match.arg(toupper(ancestry), c("AFR","AMR","CSA","EAS","EUR","MID"))

  # ---- package checks (Bioconductor) ----
  if (!requireNamespace("Rsamtools", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    stop(
      "This function requires Bioconductor packages Rsamtools, GenomicRanges, and IRanges.\n",
      "Install via: if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); ",
      "BiocManager::install(c('Rsamtools','GenomicRanges','IRanges'))",
      call. = FALSE
    )
  }

  # ---- require mapping columns ----
  req_cols <- c("panukb_chrom","panukb_pos","rsid")
  missing <- setdiff(req_cols, names(exposure_snps))
  if (length(missing)) {
    stop("exposure_snps missing columns: ", paste(missing, collapse = ", "),
         "\nDid you run exposure_snp_mapper() with sex='both'?", call. = FALSE)
  }

  # Minimal, deduped exposure lookup for joining rsid later
  exp_lu <- tibble::as_tibble(exposure_snps) |>
    dplyr::select(rsid, panukb_chrom, panukb_pos) |>
    dplyr::distinct()

  # ---- discover URL column in MR_df ----
  url_cols <- c("tabix", "tabix_url", "url", "wget", "download_cmd", "download")
  url_cols <- intersect(url_cols, names(MR_df))
  if (!length(url_cols)) {
    char_cols <- names(MR_df)[vapply(MR_df, is.character, logical(1))]
    url_cols <- char_cols
  }

  .extract_url <- function(x) {
    m <- regexpr("(https?://\\S+|s3://\\S+)", x, perl = TRUE)
    if (m[1] > 0) {
      url <- regmatches(x, m)
      url <- sub("[\"')\\s].*$", "", url)
      return(url)
    }
    if (is.character(x) && grepl("^(https?|s3)://", x) && grepl("\\.(bgz|gz)(\\?|$)", x)) return(x)
    NA_character_
  }

  # ---- helpers ----

  # prefer ancestry-suffixed columns when standardizing
  .pick_ancestry_col <- function(nms, base) {
    cand <- c(paste0(base, "_", ancestry), paste0(base, ".", ancestry), base, toupper(base), tolower(base))
    hit <- intersect(cand, nms)
    if (length(hit)) hit[[1]] else NA_character_
  }

  # programmatic standardization to TwoSampleMR-like columns
  .standardize_two_sample <- function(df, ancestry) {
    nms <- names(df)
    pick <- function(nms, base, anc = ancestry) {
      cand <- c(paste0(base, "_", anc), paste0(base, ".", anc), base, toupper(base), tolower(base))
      hit <- intersect(cand, nms)
      if (length(hit)) hit[[1]] else NA_character_
    }

    chr_col  <- pick(nms, "chrom"); if (is.na(chr_col)) chr_col <- pick(nms, "chr"); if (is.na(chr_col)) chr_col <- pick(nms, "chromosome"); if (is.na(chr_col)) chr_col <- pick(nms, "CHR")
    pos_col  <- pick(nms, "pos");   if (is.na(pos_col)) pos_col <- pick(nms, "position"); if (is.na(pos_col)) pos_col <- pick(nms, "bp"); if (is.na(pos_col)) pos_col <- pick(nms, "POS")
    ea_col   <- pick(nms, "effect_allele"); if (is.na(ea_col)) ea_col <- pick(nms, "alt"); if (is.na(ea_col)) ea_col <- pick(nms, "A1"); if (is.na(ea_col)) ea_col <- pick(nms, "EA")
    oa_col   <- pick(nms, "other_allele");  if (is.na(oa_col)) oa_col <- pick(nms, "ref"); if (is.na(oa_col)) oa_col <- pick(nms, "A2"); if (is.na(oa_col)) oa_col <- pick(nms, "OA")
    beta_col <- pick(nms, "beta")
    se_col   <- pick(nms, "se")
    p_col    <- pick(nms, "p"); if (is.na(p_col)) p_col <- pick(nms, "pval"); if (is.na(p_col)) p_col <- pick(nms, "P")

    out <- tibble::as_tibble(df)
    ren <- c(chrom = chr_col, pos = pos_col,
             effect_allele = ea_col, other_allele = oa_col,
             beta = beta_col, se = se_col, pval = p_col)
    ren <- ren[!is.na(ren)]
    if (length(ren)) out <- dplyr::rename(out, !!!ren)

    # normalize chrom as integer (strip any "chr")
    if ("chrom" %in% names(out)) {
      out$chrom <- sub("^chr", "", as.character(out$chrom))
      out$chrom <- suppressWarnings(as.integer(out$chrom))
    }
    if ("pos" %in% names(out)) {
      out$pos <- suppressWarnings(as.integer(out$pos))
    }

    keep <- intersect(c("chrom","pos","effect_allele","other_allele","beta","se","pval"), names(out))
    if (length(keep)) out <- dplyr::select(out, dplyr::all_of(keep))
    out
  }

  # Build GRanges once (numeric and chr-prefixed variants)
  gr_numeric <- GenomicRanges::GRanges(
    seqnames = as.character(exp_lu$panukb_chrom),
    ranges   = IRanges::IRanges(start = exp_lu$panukb_pos, end = exp_lu$panukb_pos)
  )
  gr_chr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", as.character(exp_lu$panukb_chrom)),
    ranges   = IRanges::IRanges(start = exp_lu$panukb_pos, end = exp_lu$panukb_pos)
  )

  MR_df$outcome_snps <- vector("list", nrow(MR_df))

  pb <- utils::txtProgressBar(min = 0, max = nrow(MR_df), style = 3)
  on.exit(close(pb), add = TRUE)

  for (i in seq_len(nrow(MR_df))) {
    utils::setTxtProgressBar(pb, i)

    # find a usable URL field on this row
    rec <- MR_df[i, , drop = FALSE]
    url <- NA_character_
    for (col in url_cols) {
      val <- rec[[col]][[1]]
      if (is.character(val) && nzchar(val)) {
        cand <- .extract_url(val)
        if (!is.na(cand)) { url <- cand; break }
      }
    }
    if (is.na(url)) {
      if (verbose) logger::log_warn("Pan-UKB row {i}: no sumstats URL/wget found; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      next
    }

    # Open remote Tabix (expects <url>.tbi to exist)
    tf <- Rsamtools::TabixFile(url, index = paste0(url, ".tbi"))
    ok_open <- TRUE
    tryCatch(Rsamtools::open.TabixFile(tf), error = function(e) { ok_open <<- FALSE })
    if (!ok_open) {
      if (verbose) logger::log_warn("Pan-UKB row {i}: failed to open Tabix at {url}; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      next
    }
    on.exit(try(Rsamtools::close.TabixFile(tf), silent = TRUE), add = TRUE)

    # Grab header (to recover column names)
    hdr <- try(Rsamtools::headerTabix(tf), silent = TRUE)
    header_line <- NA_character_
    if (!inherits(hdr, "try-error") && length(hdr)) {
      header_line <- sub("^#", "", utils::tail(hdr, 1))
    }

    # Query numeric seqnames first; if empty, retry with chr-prefixed
    res <- Rsamtools::scanTabix(tf, param = gr_numeric)
    if (all(lengths(res) == 0L)) {
      res <- Rsamtools::scanTabix(tf, param = gr_chr)
    }
    lines <- unlist(res, use.names = FALSE)

    if (!length(lines)) {
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      next
    }

    # Parse to data.frame (prefer attaching header if available)
    if (is.character(header_line) && nzchar(header_line)) {
      txt <- paste(c(header_line, lines), collapse = "\n")
      panukb_tmp <- data.table::fread(text = txt, header = TRUE, showProgress = FALSE) |> tibble::as_tibble()
    } else {
      # Fallback: try to infer headerless table
      panukb_tmp <- data.table::fread(text = paste(lines, collapse = "\n"),
                                      header = FALSE, showProgress = FALSE, fill = TRUE) |> tibble::as_tibble()
    }

    # Standardize columns (ancestry-aware), then join RSIDs by chrom/pos
    panukb_std <- .standardize_two_sample(panukb_tmp, ancestry)

    if (all(c("chrom","pos") %in% names(panukb_std))) {
      panukb_std <- dplyr::inner_join(
        panukb_std,
        exp_lu,
        by = c("chrom" = "panukb_chrom", "pos" = "panukb_pos")
      )
      if ("rsid" %in% names(panukb_std)) {
        panukb_std <- dplyr::relocate(panukb_std, rsid, .before = 1)
        names(panukb_std)[names(panukb_std) == "rsid"] <- "SNP"
      }
    } else if (verbose) {
      logger::log_warn("Pan-UKB row {i}: could not detect chrom/pos to join RSIDs")
    }

    MR_df$outcome_snps[[i]] <- panukb_std

    # Close handle for this file
    try(Rsamtools::close.TabixFile(tf), silent = TRUE)
  }

  total_rows <- sum(vapply(MR_df$outcome_snps, nrow, integer(1)))
  if (verbose) {
    logger::log_info("Pan-UKB SNP grabber (Rsamtools): {nrow(MR_df)} outcomes; {total_rows} SNP rows fetched")
  }

  MR_df
}
