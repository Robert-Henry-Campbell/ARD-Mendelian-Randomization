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
panukb_snp_grabber <- function(exposure_snps, MR_df, ancestry, cache_dir = ardmr_cache_dir(), verbose = TRUE, force_refresh = FALSE) {
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
  req_cols <- c("panukb_chrom","panukb_pos","rsid","effect_allele.exposure")
  if (!"effect_allele.exposure" %in% names(exposure_snps) &&
      "effect_allele" %in% names(exposure_snps)) {
    exposure_snps <- dplyr::rename(exposure_snps, effect_allele.exposure = effect_allele)
  }
  missing <- setdiff(req_cols, names(exposure_snps))
  if (length(missing)) {
    stop("exposure_snps missing columns: ", paste(missing, collapse = ", "),
         "\nDid you run exposure_snp_mapper() with sex='both'?", call. = FALSE)
  }

  exposure_snps <- exposure_snps |>
    dplyr::mutate(effect_allele.exposure = toupper(.data$effect_allele.exposure))
  bad_alleles <- unique(exposure_snps$effect_allele.exposure[
    !exposure_snps$effect_allele.exposure %in% c("A", "C", "G", "T")
  ])
  if (length(bad_alleles)) {
    stop(
      "effect_allele.exposure contains non-ACGT values: ",
      paste(bad_alleles, collapse = ", "),
      call. = FALSE
    )
  }

  exp_lu <- tibble::as_tibble(exposure_snps) |>
    dplyr::select(rsid, panukb_chrom, panukb_pos) |>
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


  # ---- header utilities (EXACT matching) --------------------------------------
  .get_sumstats_header <- function(tf, data_url, verbose = TRUE) {
    hdr_lines <- tryCatch(Rsamtools::headerTabix(tf), error = function(e) character())
    hdr_lines <- unlist(hdr_lines, use.names = FALSE)

    hdr_line <- NA_character_
    if (length(hdr_lines)) {
      cand <- hdr_lines[grepl("\t", hdr_lines, fixed = TRUE)]
      if (length(cand)) hdr_line <- utils::tail(cand, 1)
    }

    warn_msg <- NULL
    if (is.na(hdr_line) || !nzchar(hdr_line)) {
      con <- NULL
      hdr_line <- tryCatch({
        con <- gzcon(url(data_url, open = "rb"))
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        withCallingHandlers(
          readLines(con, n = 1),
          warning = function(w) {
            warn_msg <<- conditionMessage(w)
            invokeRestart("muffleWarning")
          }
        )
      }, error = function(e) {
        warn_msg <<- conditionMessage(e)
        NA_character_
      })
    }

    if (length(hdr_line) > 1) hdr_line <- hdr_line[[1]]

    if (is.null(hdr_line) || !length(hdr_line) || is.na(hdr_line) || !nzchar(hdr_line)) {
      warn_text <- if (!is.null(warn_msg) && nzchar(warn_msg)) warn_msg else "no warning provided"
      reason <- sprintf(
        "failed to retrieve header from %s via readLines; warning: %s",
        data_url,
        warn_text
      )
      return(list(fields = NULL, error = reason))
    }

    hdr_line <- sub("^#+", "", hdr_line)
    fields <- strsplit(hdr_line, "\t", fixed = TRUE)[[1]]
    list(fields = trimws(fields), error = NULL)
  }

  # Build an EXACT, case-sensitive mapping for the columns we require.
  # We REQUIRE: chr, pos, ref, alt, and the ancestry-suffixed
  # (Optionally low_confidence_<ANC> if present; if not, set NA.)

  # 1) Replace your current .build_header_mapping_exact() helper with this:
  .build_header_mapping_eaf_controls_fallback <- function(header_fields, ancestry) {
    anc <- toupper(ancestry)  # Pan-UKB uses upper-case suffixes

    # required (exact)
    req <- c("chr", "pos", "ref", "alt",
             paste0("beta_", anc),
             paste0("se_", anc),
             paste0("neglog10_pval_", anc))
    idx_req <- match(req, header_fields)
    if (any(is.na(idx_req))) {
      stop(
        "Header missing required columns: ",
        paste(req[is.na(idx_req)], collapse = ", "),
        "\nHeader was: ", paste(header_fields, collapse = " | "),
        call. = FALSE
      )
    }

    # EAF sources (prefer af_<ANC>, else af_controls_<ANC>)
    idx_af_primary  <- match(paste0("af_", anc), header_fields)
    idx_af_controls <- match(paste0("af_controls_", anc), header_fields)
    af_idx <- if (!is.na(idx_af_primary)) idx_af_primary else idx_af_controls
    af_src <- if (!is.na(idx_af_primary)) "primary" else if (!is.na(idx_af_controls)) "controls" else "absent"

    list(
      chr      = match("chr", header_fields),
      pos      = match("pos", header_fields),
      ref      = match("ref", header_fields),
      alt      = match("alt", header_fields),
      af       = af_idx,                                # may be NA -> we'll handle below
      beta     = match(paste0("beta_", anc), header_fields),
      se       = match(paste0("se_", anc), header_fields),
      p_log10  = match(paste0("neglog10_pval_", anc), header_fields),
      low_conf = match(paste0("low_confidence_", anc), header_fields),  # may be NA
      af_src   = af_src,
      p_is_log10 = TRUE
    )
  }

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

    # -- optional hard reset of cache for this phenotype
    if (isTRUE(force_refresh) && file.exists(cache_file)) {
      if (verbose) logger::log_info("Pan-UKB row {i}: force_refresh -> deleting cache {basename(cache_file)}")
      try(unlink(cache_file), silent = TRUE)
    }

    # resume from cache if available
    if (!isTRUE(force_refresh) && file.exists(cache_file)) {
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

    # ---- read lines and <- columns by the actual header for this GWAS ----
    header_info <- .get_sumstats_header(tf, url, verbose = verbose)
    if (!is.null(header_info$error)) {
      reason <- header_info$error
      if (verbose) {
        logger::log_warn(
          "Pan-UKB row {i}: unable to download header ({reason}); skipping phenotype",
          reason = reason
        )
      }
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
      .close_tf(tf)
      next
    }
    header_fields <- header_info$fields
    if (verbose) logger::log_info("Pan-UKB row {i}: detected header fields: {paste(header_fields, collapse=' | ')}")

    panukb_tmp <- data.table::fread(
      text = paste(lines, collapse = "\n"),
      header = FALSE, fill = TRUE, showProgress = FALSE
    )

    # Build mapping from header (EXACT names); apply to the fetched rows
    map <- .build_header_mapping_eaf_controls_fallback(header_fields, ancestry)
    if (verbose && identical(map$af_src, "controls")) {
      logger::log_warn("Pan-UKB row {i}: using af_controls_{ancestry} for EAF.")
    }
    if (verbose && identical(map$af_src, "absent")) {
      logger::log_warn("Pan-UKB row {i}: no af_{ancestry} or af_controls_{ancestry}; EAF will be NA (more palindromes may drop).")
    }

    # --- PART 3: tolerant subset/rename that skips NA af/low_conf and backfills if absent ---
    keep_idx <- c(map$chr, map$pos, map$ref, map$alt,
                  if (!is.na(map$af)) map$af else NULL,
                  map$beta, map$se, map$p_log10,
                  if (!is.na(map$low_conf)) map$low_conf else NULL)
    panukb_tmp <- panukb_tmp[, keep_idx, with = FALSE]

    newn <- c("chr","pos","ref","alt",
              if (!is.na(map$af)) "af" else NULL,
              "beta","se","p_log10",
              if (!is.na(map$low_conf)) "low_conf" else NULL)
    data.table::setnames(panukb_tmp, newn)
    panukb_tmp <- tibble::as_tibble(panukb_tmp)

    # Backfill absent columns to keep schema stable for downstream code
    if (!"af" %in% names(panukb_tmp))       panukb_tmp$af       <- NA_real_
    if (!"low_conf" %in% names(panukb_tmp)) panukb_tmp$low_conf <- NA


    # ---- standardise for MR (EA=ALT; EAF=af; pval as -log10) ----
    panukb_std <- panukb_tmp |>
      dplyr::transmute(
        chr = suppressWarnings(as.integer(sub("^chr", "", as.character(.data$chr)))),
        pos = suppressWarnings(as.integer(.data$pos)),
        effect_allele = toupper(as.character(.data$alt)),
        other_allele  = toupper(as.character(.data$ref)),
        beta          = suppressWarnings(as.numeric(.data$beta)),
        se            = suppressWarnings(as.numeric(.data$se)),
        eaf           = suppressWarnings(as.numeric(.data$af)),
        p_log10       = suppressWarnings(as.numeric(.data$p_log10)),  # already -log10(p)
        low_confidence = dplyr::case_when(
          is.null(.data$low_conf) ~ NA,                              # column absent
          TRUE ~ tolower(as.character(.data$low_conf)) %in% c("true","t","1")
        )
      )




    # ---- join RSIDs by chr:pos using your exposure mapping ----
    panukb_std <- dplyr::inner_join(
      panukb_std,
      exp_lu,
      by = c("chr" = "panukb_chrom", "pos" = "panukb_pos")
    )

    # If multiple outcome rows map to the same rsID (multi-allelic etc.), keep the strongest line
    panukb_std <- panukb_std |>
      dplyr::arrange(.data$rsid, dplyr::desc(.data$p_log10)) |>
      dplyr::distinct(.data$rsid, .keep_all = TRUE)

    panukb_std <- dplyr::relocate(panukb_std, rsid, .before = 1)
    names(panukb_std)[names(panukb_std) == "rsid"] <- "SNP"
    panukb_std$chrpos    <- paste0(panukb_std$chr, ":", panukb_std$pos)
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

    # Sanity: p from -log10 vs beta/se agree?
    p2s   <- 10^(-panukb_std$p_log10)          # two-sided p from the raw -log10(p)
    z_p   <- stats::qnorm(p2s/2, lower.tail = FALSE)
    z_bse <- abs(panukb_fmt$beta.outcome / panukb_fmt$se.outcome)
    bad   <- which(is.finite(z_p) & is.finite(z_bse) & abs(z_p - z_bse) > 3)

    if (length(bad)) {
      stop(sprintf("Outcome '%s': %d SNPs show beta/se vs p inconsistency (check units/columns).",
                   outcome_label, length(bad)))
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
