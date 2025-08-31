#' Fetch outcome SNP rows from Pan-UKB sumstats for each outcome
#'
#' @param exposure_snps mapped SNPs (must include `panukb_chrom` and `panukb_pos`)
#' @param MR_df tibble with provider info (Pan-UKB rows and a column containing a sumstats URL or wget cmd)
#' @param ancestry character (e.g. "EUR")
#' @param cache_dir (unused here but kept for parity)
#' @param verbose logical
#' @return MR_df with a new list-column `outcome_snps` (TwoSampleMR-formatted)
#' @export
panukb_snp_grabber <- function(exposure_snps, MR_df, ancestry, cache_dir, verbose = TRUE) {
  stopifnot(is.data.frame(exposure_snps), is.data.frame(MR_df))
  ancestry <- match.arg(toupper(ancestry), c("AFR","AMR","CSA","EAS","EUR","MID"))

  # Require tabix
  if (nzchar(Sys.which("tabix")) == FALSE) {
    stop("tabix not found in PATH. Please install htslib (tabix/bgzip) and try again.", call. = FALSE)
  }

  # Require mapping columns
  req_cols <- c("panukb_chrom","panukb_pos","rsid")
  missing <- setdiff(req_cols, names(exposure_snps))
  if (length(missing)) {
    stop("exposure_snps missing columns: ", paste(missing, collapse = ", "),
         "\nDid you run exposure_snp_mapper() with sex='both'?", call. = FALSE)
  }

  # Minimal, deduped exposure lookup for joining rsid later
  exp_lu <- exposure_snps |>
    dplyr::select(rsid, panukb_chrom, panukb_pos) |>
    dplyr::distinct()

  # Discover a column in MR_df that holds a usable URL or a wget command
  url_cols <- c("tabix", "tabix_url", "url", "wget", "download_cmd", "download")
  url_cols <- intersect(url_cols, names(MR_df))
  if (!length(url_cols)) {
    # Fallback: try ANY character column that looks like an http/s3 and ends with .bgz or .gz
    char_cols <- names(MR_df)[vapply(MR_df, is.character, logical(1))]
    url_cols <- char_cols
  }

  .extract_url <- function(x) {
    # if it's a wget line, grab the first http/https/s3 URL
    m <- regexpr("(https?://\\S+|s3://\\S+)", x, perl = TRUE)
    if (m[1] > 0) {
      url <- regmatches(x, m)
      # strip trailing characters (quotes, ) , etc.)
      url <- sub("[\"')\\s].*$", "", url)
      return(url)
    }
    # already a bare URL?
    if (is.character(x) && grepl("^(https?|s3)://", x) &&
        grepl("\\.(bgz|gz)(\\?|$)", x)) {
      return(x)
    }
    NA_character_
  }

  # Helper: choose ancestry-specific columns if available
  .pick_ancestry_col <- function(nms, base) {
    # e.g. prefer beta_EUR over beta
    cand <- c(paste0(base, "_", ancestry), paste0(base, ".", ancestry), base, toupper(base), tolower(base))
    hit <- intersect(cand, nms)
    if (length(hit)) hit[[1]] else NA_character_
  }

  # Helper: standardize fetched table to TwoSampleMR-ish
  .standardize_two_sample <- function(df, ancestry) {
    nms <- names(df)

    chr_col <- dplyr::coalesce(
      .pick_ancestry_col(nms, "chrom"), .pick_ancestry_col(nms, "chr"), .pick_ancestry_col(nms, "chromosome"), .pick_ancestry_col(nms, "CHR")
    )
    pos_col <- dplyr::coalesce(
      .pick_ancestry_col(nms, "pos"), .pick_ancestry_col(nms, "position"), .pick_ancestry_col(nms, "bp"), .pick_ancestry_col(nms, "POS")
    )
    ea_col  <- dplyr::coalesce(.pick_ancestry_col(nms, "effect_allele"), .pick_ancestry_col(nms, "alt"), .pick_ancestry_col(nms, "A1"), .pick_ancestry_col(nms, "EA"))
    oa_col  <- dplyr::coalesce(.pick_ancestry_col(nms, "other_allele"),  .pick_ancestry_col(nms, "ref"), .pick_ancestry_col(nms, "A2"), .pick_ancestry_col(nms, "OA"))
    beta_col<- .pick_ancestry_col(nms, "beta")
    se_col  <- .pick_ancestry_col(nms, "se")
    p_col   <- .pick_ancestry_col(nms, "p") %||% .pick_ancestry_col(nms, "pval") %||% .pick_ancestry_col(nms, "P")

    out <- df
    # best-effort rename if present
    ren <- c()
    if (!is.na(chr_col))  ren[chr_col]  <- "chrom"
    if (!is.na(pos_col))  ren[pos_col]  <- "pos"
    if (!is.na(ea_col))   ren[ea_col]   <- "effect_allele"
    if (!is.na(oa_col))   ren[oa_col]   <- "other_allele"
    if (!is.na(beta_col)) ren[beta_col] <- "beta"
    if (!is.na(se_col))   ren[se_col]   <- "se"
    if (!is.na(p_col))    ren[p_col]    <- "pval"
    if (length(ren)) out <- dplyr::rename(out, dplyr::all_of(ren))

    # Keep only the essentials if present
    keep <- intersect(c("chrom","pos","effect_allele","other_allele","beta","se","pval"), names(out))
    if (length(keep)) {
      out <- dplyr::select(out, dplyr::all_of(keep))
    }
    tibble::as_tibble(out)
  }

  # Build regions vector once
  regions <- paste0(exp_lu$panukb_chrom, ":", exp_lu$panukb_pos, "-", exp_lu$panukb_pos)

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

    # tabix each region; rbind results
    # (If Pan-UKB needs chromosome prefix 'chr', adjust here)
    fetched <- lapply(seq_along(regions), function(j) {
      reg <- regions[j]
      cmd <- sprintf("tabix -h %s %s", shQuote(url), shQuote(reg))
      dt <- try(data.table::fread(cmd = cmd, showProgress = FALSE), silent = TRUE)
      if (inherits(dt, "try-error")) return(NULL)
      as.data.frame(dt)
    })
    fetched <- fetched[!vapply(fetched, is.null, logical(1))]
    if (!length(fetched)) {
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      next
    }
    panukb_tmp <- tibble::as_tibble(data.table::rbindlist(fetched, fill = TRUE))

    # Try to standardize columns (ancestry-aware picks)
    panukb_std <- .standardize_two_sample(panukb_tmp, ancestry)

    # Append rsid by joining back on chrom/pos vs panukb_chrom/panukb_pos
    join_by <- NULL
    # Detect names to join on in standardized output
    if (all(c("chrom","pos") %in% names(panukb_std))) {
      join_by <- c("chrom" = "panukb_chrom", "pos" = "panukb_pos")
      panukb_std <- dplyr::inner_join(panukb_std, exp_lu, by = join_by)
      # move rsid first if present
      if ("rsid" %in% names(panukb_std)) {
        panukb_std <- dplyr::relocate(panukb_std, rsid, .before = 1)
        names(panukb_std)[names(panukb_std) == "rsid"] <- "SNP"
      }
    } else {
      # fallback: leave as-is
      if (verbose) logger::log_warn("Pan-UKB row {i}: could not detect chrom/pos to join RSIDs")
    }

    MR_df$outcome_snps[[i]] <- panukb_std
  }

  total_rows <- sum(vapply(MR_df$outcome_snps, nrow, integer(1)))
  if (verbose) {
    logger::log_info("Pan-UKB SNP grabber: {nrow(MR_df)} outcomes; {total_rows} SNP rows fetched")
  }

  MR_df
}
