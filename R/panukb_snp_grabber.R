# R/panukb_snp_grabber.R

#' Fetch outcome SNP rows from Pan-UKB sumstats for each outcome (R-based Tabix)
#'
#' Uses Rsamtools/Rhtslib to remotely query tabix-indexed BGZF files over HTTP(S),
#' removing the need for a system 'tabix' in PATH. Standardizes output with
#' TwoSampleMR::format_data(type = "outcome").
#'
#' @param exposure_snps mapped SNPs (must include `panukb_chrom` and `panukb_pos`)
#' @param MR_df tibble with provider info (must contain 'aws_link' and 'aws_link_tabix')
#' @param ancestry character (AFR, AMR, CSA, EAS, EUR, MID)
#' @param cache_dir (unused here; kept for parity)
#' @param verbose logical
#' @return MR_df with a new list-column `outcome_snps` (TwoSampleMR-formatted)
#' @export
#' @importFrom Rsamtools TabixFile headerTabix scanTabix open.TabixFile close.TabixFile
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr select distinct inner_join relocate rename all_of mutate transmute
#' @importFrom rlang .data
panukb_snp_grabber <- function(exposure_snps, MR_df, ancestry, cache_dir = ardmr_cache_dir(), verbose = TRUE) {
  stopifnot(is.data.frame(exposure_snps), is.data.frame(MR_df))
  ancestry <- match.arg(toupper(ancestry), c("AFR","AMR","CSA","EAS","EUR","MID"))

  # ---- Bioconductor deps ----
  if (!requireNamespace("Rsamtools", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    stop(
      "Requires Bioconductor packages Rsamtools, GenomicRanges, IRanges.\n",
      "Install: if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); ",
      "BiocManager::install(c('Rsamtools','GenomicRanges','IRanges'))",
      call. = FALSE
    )
  }

  # ---- exposure mapping cols ----
  req_cols <- c("panukb_chrom","panukb_pos","rsid")
  missing <- setdiff(req_cols, names(exposure_snps))
  if (length(missing)) {
    stop("exposure_snps missing columns: ", paste(missing, collapse = ", "),
         "\nDid you run exposure_snp_mapper() with sex='both'?", call. = FALSE)
  }

  exp_lu <- tibble::as_tibble(exposure_snps) |>
    dplyr::select(rsid, panukb_chrom, panukb_pos, effect_allele) |>
    dplyr::distinct()

  # ---- enforce link columns in MR_df ----
  if (!all(c("aws_link", "aws_link_tabix") %in% names(MR_df))) {
    stop("MR_df must contain columns 'aws_link' and 'aws_link_tabix'.", call. = FALSE)
  }

  .clean_url <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- sub("^(['\"])(.*)\\1$", "\\2", x, perl = TRUE)  # strip quotes
    # optional: normalise s3://bucket/key -> https://bucket.s3.amazonaws.com/key
    if (grepl("^s3://", x)) x <- sub("^s3://([^/]+)/(.*)$", "https://\\1.s3.amazonaws.com/\\2", x)
    x
  }

  # also require an outcome label column
  if (!"description" %in% names(MR_df)) {
    stop("MR_df must contain a 'description' column to label outcomes.", call. = FALSE)
  }

  # ---- expected Pan-UKB column order (headerless tabix output) ----
  base_need <- c("chr","pos","ref","alt")
  suf <- paste0("_", ancestry)
  anc_need <- paste0(c("af","beta","se","neglog10_pval","low_confidence"), suf)
  need <- c(base_need, anc_need)

  # ---- direct htslib's remote cache into our cache_dir ----
  cache_dir <- normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
  hts_cache <- file.path(cache_dir, "hts-cache")
  dir.create(hts_cache, recursive = TRUE, showWarnings = FALSE)
  old_hts_cache <- Sys.getenv("HTS_CACHE_DIR", unset = NA)
  Sys.setenv(HTS_CACHE_DIR = hts_cache)
  on.exit({
    if (is.na(old_hts_cache)) Sys.unsetenv("HTS_CACHE_DIR") else Sys.setenv(HTS_CACHE_DIR = old_hts_cache)
  }, add = TRUE)

  # cache directory for downloaded .tbi files (so they don't land in getwd())
  idx_cache_dir <- file.path(cache_dir, "tabix_index")
  dir.create(idx_cache_dir, recursive = TRUE, showWarnings = FALSE)


  # ---- setup RDS caching per-phenotype SNPs ----
  cache_root <- file.path(cache_dir, "panukb_outcome_snps", ancestry)
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  .slug <- function(x) {
    x <- as.character(x); x[is.na(x) | !nzchar(x)] <- "NA"
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    substr(x, 1, 120)
  }
  .close_tf <- function(tf) try(Rsamtools::close.TabixFile(tf), silent = TRUE)

  # ---- build GRanges once ----
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

    # explicit data and index URLs
    rec <- MR_df[i, , drop = FALSE]
    outcome_label <- as.character(rec$description[[1]])
    if (is.na(outcome_label) || !nzchar(outcome_label)) {
      outcome_label <- sprintf("panUKB_%03d", i)
    }
    cache_file <- file.path(cache_root, paste0(.slug(outcome_label), ".rds"))

    # resume from cache if available
    if (file.exists(cache_file)) {
      if (verbose) logger::log_info("Pan-UKB row {i}: using cached result {basename(cache_file)}")
      MR_df$outcome_snps[[i]] <- readRDS(cache_file)
      next
    }

    url <- .clean_url(rec$aws_link[[1]])
    idx <- .clean_url(rec$aws_link_tabix[[1]])

    if (!is.character(url) || !nzchar(url) || !grepl("^(https?|s3)://", url)) {
      if (verbose) logger::log_warn("Pan-UKB row {i}: invalid or missing 'aws_link'; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
      next
    }
    if (!is.character(idx) || !nzchar(idx) || !grepl("^(https?|s3)://", idx)) {
      if (verbose) logger::log_warn("Pan-UKB row {i}: invalid or missing 'aws_link_tabix'; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
      next
    }

    # Open remote tabix
    # ---- ensure the .tbi is cached locally (avoid Rsamtools writing to getwd) ----
    idx_hash  <- tryCatch(digest::digest(idx, algo = "xxhash64"),
                          error = function(e) digest::digest(idx))
    idx_local <- file.path(idx_cache_dir, paste0(.slug(outcome_label), "_", idx_hash, ".tbi"))

    if (!file.exists(idx_local)) {
      if (verbose) logger::log_info("Pan-UKB row {i}: caching tabix index => {basename(idx_local)}")
      ok_dl <- tryCatch({
        utils::download.file(idx, destfile = idx_local, mode = "wb", quiet = !verbose)
        TRUE
      }, error = function(e) FALSE)
      if (!ok_dl || !file.exists(idx_local)) {
        if (verbose) logger::log_warn("Pan-UKB row {i}: failed to download index {idx}; skipping")
        MR_df$outcome_snps[[i]] <- tibble::tibble()
        next
      }
    }

    # Open tabix using remote data URL + local index path
    tf <- Rsamtools::TabixFile(file = url, index = idx_local)
    ok_open <- TRUE
    tryCatch(Rsamtools::open.TabixFile(tf), error = function(e) { ok_open <<- FALSE })

    if (!ok_open) {
      if (verbose) logger::log_warn("Pan-UKB row {i}: failed to open Tabix (data={url}, index={idx}); skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
      next
    }

    # query all regions at once; if empty, retry with 'chr' prefix
    res <- Rsamtools::scanTabix(tf, param = gr_numeric)
    if (all(lengths(res) == 0L)) res <- Rsamtools::scanTabix(tf, param = gr_chr)
    lines <- unlist(res, use.names = FALSE)

    if (!length(lines)) {
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
      .close_tf(tf)
      next
    }

    # ---- parse headerless and assign the expected names ----
    panukb_tmp <- data.table::fread(
      text = paste(lines, collapse = "\n"),
      header = FALSE, fill = TRUE, showProgress = FALSE
    )

    # Must have at least the columns we expect; if more, drop extras; if fewer, skip
    if (ncol(panukb_tmp) < length(need)) {
      if (verbose) logger::log_warn("Pan-UKB row {i}: expected >= {length(need)} columns, found {ncol(panukb_tmp)}; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
      .close_tf(tf)
      next
    }
    if (ncol(panukb_tmp) > length(need)) {
      panukb_tmp <- panukb_tmp[, seq_len(length(need))]
    }
    data.table::setnames(panukb_tmp, need)
    panukb_tmp <- tibble::as_tibble(panukb_tmp)

    # ---- standardise for MR (EA=ALT; EAF=af_<ANC>; pval is -log10) ----
    pval_col <- paste0("neglog10_pval", suf)
    beta_col <- paste0("beta", suf)
    se_col   <- paste0("se", suf)
    eaf_col  <- paste0("af", suf)
    lowq_col <- paste0("low_confidence", suf)

    panukb_std <- panukb_tmp |>
      dplyr::transmute(
        chr = as.integer(sub("^chr", "", as.character(.data$chr))),
        pos = as.integer(.data$pos),
        effect_allele = toupper(.data$alt),
        other_allele  = toupper(.data$ref),
        beta          = suppressWarnings(as.numeric(.data[[beta_col]])),
        se            = suppressWarnings(as.numeric(.data[[se_col]])),
        eaf           = suppressWarnings(as.numeric(.data[[eaf_col]])),
        p_log10       = suppressWarnings(as.numeric(.data[[pval_col]])),
        low_confidence = tolower(as.character(.data[[lowq_col]])) %in% c("true","t","1")
      )

    # ---- join RSIDs by chr:pos using your exposure mapping ----
    panukb_std <- dplyr::inner_join(
      panukb_std,
      exp_lu,
      by = c("chr" = "panukb_chrom", "pos" = "panukb_pos", "effect_allele")
    )

    if (!"rsid" %in% names(panukb_std)) {
      if (verbose) logger::log_warn("Pan-UKB row {i}: RSID join failed; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
      .close_tf(tf)
      next
    }

    panukb_std <- dplyr::relocate(panukb_std, rsid, .before = 1)
    names(panukb_std)[names(panukb_std) == "rsid"] <- "SNP"
    panukb_std <- dplyr::distinct(panukb_std, SNP, effect_allele, .keep_all = TRUE)
    panukb_std$chrpos <- paste0(panukb_std$chr, ":", panukb_std$pos)
    panukb_std$Phenotype <- outcome_label

    # ---- TwoSampleMR formatting (log10 p-values) ----
    if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
      panukb_fmt <- TwoSampleMR::format_data(
        dat               = panukb_std,
        type              = "outcome",
        snp_col           = "SNP",
        beta_col          = "beta",
        se_col            = "se",
        effect_allele_col = "effect_allele",
        other_allele_col  = "other_allele",
        eaf_col           = "eaf",
        pval_col          = "p_log10",
        chr_col           = "chr",
        pos_col           = "pos",
        log_pval          = TRUE,
        phenotype_col     = "Phenotype"
      )
      panukb_fmt$chrpos <- panukb_std$chrpos
      panukb_fmt$low_confidence <- panukb_std$low_confidence
      MR_df$outcome_snps[[i]] <- tibble::as_tibble(panukb_fmt)
    } else {
      panukb_std$pval <- 10^(-panukb_std$p_log10)
      MR_df$outcome_snps[[i]] <- tibble::as_tibble(
        dplyr::select(panukb_std, SNP, chr, pos, effect_allele, other_allele, beta, se, eaf, pval, chrpos, low_confidence)
      )
    }

    # cache formatted result for this phenotype and close tabix
    try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
    .close_tf(tf)
  }

  total_rows <- sum(vapply(MR_df$outcome_snps, nrow, integer(1)))
  if (verbose) {
    logger::log_info("Pan-UKB SNP grabber (Rsamtools + format_data): {nrow(MR_df)} outcomes; {total_rows} SNP rows fetched")
  }
  MR_df
}
