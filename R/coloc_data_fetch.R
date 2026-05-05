# R/coloc_data_fetch.R
# Regional summary-statistic fetchers for colocalization analysis.

.coloc_cache_path <- function(cache_dir, kind, key) {
  d <- file.path(cache_dir, "coloc", kind)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, paste0(key, ".rds"))
}

.coloc_region_key <- function(...) digest::digest(list(...), algo = "xxhash64")

.coloc_clean_url <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(['\"])(.*)\\1$", "\\2", x, perl = TRUE)
  if (grepl("^s3://", x)) x <- sub("^s3://([^/]+)/(.*)$", "https://\\1.s3.amazonaws.com/\\2", x)
  x
}

#' Fetch exposure regional summary statistics from OpenGWAS
#'
#' @param exposure_id OpenGWAS study id (e.g. `"ieu-b-38"`).
#' @param chr Chromosome (numeric or character).
#' @param start,end 1-based inclusive region coordinates.
#' @param cache_dir Cache directory.
#' @param verbose Logical.
#' @return Tibble with columns `SNP, chr, pos, effect_allele, other_allele, beta, se, eaf, pval`.
#' @keywords internal
#' @export
coloc_fetch_exposure_region <- function(exposure_id, chr, start, end,
                                        cache_dir = ardmr_cache_dir(),
                                        verbose = TRUE) {
  if (!requireNamespace("ieugwasr", quietly = TRUE)) {
    stop("coloc_fetch_exposure_region(): requires ieugwasr.", call. = FALSE)
  }
  key <- .coloc_region_key("exposure", exposure_id, as.character(chr), as.integer(start), as.integer(end))
  cache_path <- .coloc_cache_path(cache_dir, "exposure", key)
  if (file.exists(cache_path)) {
    if (verbose) logger::log_info("coloc: cached exposure region {chr}:{start}-{end} for {exposure_id}")
    return(readRDS(cache_path))
  }
  region <- sprintf("%s:%d-%d", as.character(chr), as.integer(start), as.integer(end))
  if (verbose) logger::log_info("coloc: ieugwasr::associations id={exposure_id} region={region}")
  res <- tryCatch(
    ieugwasr::associations(variants = region, id = exposure_id),
    error = function(e) {
      logger::log_warn("ieugwasr::associations failed for {exposure_id} {region}: {conditionMessage(e)}")
      NULL
    }
  )
  fetch_failed <- is.null(res)
  out <- .coloc_standardise_exposure(res)
  if (!fetch_failed && nrow(out) > 0) saveRDS(out, cache_path)
  out
}

.coloc_or <- function(a, b) if (is.null(a)) b else a

.coloc_standardise_exposure <- function(res) {
  empty <- tibble::tibble(SNP = character(), chr = integer(), pos = integer(),
                          effect_allele = character(), other_allele = character(),
                          beta = double(), se = double(), eaf = double(), pval = double())
  if (is.null(res) || NROW(res) == 0) return(empty)
  res <- tibble::as_tibble(res)
  pick <- function(nms) {
    for (n in nms) if (n %in% names(res)) return(res[[n]])
    NULL
  }
  snp <- pick(c("rsid", "SNP", "name"))
  ch  <- pick(c("chr", "chromosome"))
  ps  <- pick(c("position", "pos", "bp"))
  ea  <- pick(c("ea", "effect_allele"))
  nea <- pick(c("nea", "other_allele", "non_effect_allele"))
  b   <- pick(c("beta", "b"))
  s   <- pick(c("se", "standard_error"))
  ef  <- pick(c("eaf", "af", "maf"))
  p   <- pick(c("p", "pval", "p_value"))
  if (is.null(snp) || is.null(b) || is.null(s)) return(empty)
  out <- tibble::tibble(
    SNP           = as.character(snp),
    chr           = suppressWarnings(as.integer(sub("^chr", "", as.character(.coloc_or(ch, NA))))),
    pos           = suppressWarnings(as.integer(.coloc_or(ps, NA))),
    effect_allele = toupper(as.character(.coloc_or(ea, NA))),
    other_allele  = toupper(as.character(.coloc_or(nea, NA))),
    beta          = suppressWarnings(as.numeric(b)),
    se            = suppressWarnings(as.numeric(s)),
    eaf           = suppressWarnings(as.numeric(.coloc_or(ef, NA))),
    pval          = suppressWarnings(as.numeric(.coloc_or(p, NA)))
  )
  out[!is.na(out$SNP) & !is.na(out$beta) & !is.na(out$se) & out$se > 0, , drop = FALSE]
}

#' Fetch exposure metadata (N, type, ncase/ncontrol) from OpenGWAS
#' @keywords internal
#' @export
coloc_fetch_metadata <- function(exposure_id, cache_dir = ardmr_cache_dir(), verbose = TRUE) {
  if (!requireNamespace("ieugwasr", quietly = TRUE)) {
    stop("coloc_fetch_metadata(): requires ieugwasr.", call. = FALSE)
  }
  key <- .coloc_region_key("meta", exposure_id)
  cache_path <- .coloc_cache_path(cache_dir, "meta", key)
  if (file.exists(cache_path)) return(readRDS(cache_path))
  info <- tryCatch(ieugwasr::gwasinfo(exposure_id), error = function(e) NULL)
  out <- list(N = NA_real_, type = "quant", ncase = NA_real_, ncontrol = NA_real_, s = NA_real_)
  if (!is.null(info) && NROW(info) > 0) {
    info <- as.list(tibble::as_tibble(info)[1, ])
    n_total <- suppressWarnings(as.numeric(.coloc_or(.coloc_or(info$sample_size, info$n), NA)))
    ncase <- suppressWarnings(as.numeric(.coloc_or(info$ncase, NA)))
    ncontrol <- suppressWarnings(as.numeric(.coloc_or(info$ncontrol, NA)))
    is_cc <- (!is.na(ncase) && ncase > 0) ||
             (!is.null(info$unit) && grepl("logOR|log\\(OR\\)|log\\s*odds", info$unit, ignore.case = TRUE))
    out$type <- if (is_cc) "cc" else "quant"
    out$ncase <- ncase
    out$ncontrol <- ncontrol
    out$N <- if (is.na(n_total) && !is.na(ncase) && !is.na(ncontrol)) ncase + ncontrol else n_total
    out$s <- if (!is.na(ncase) && !is.na(ncontrol) && (ncase + ncontrol) > 0) ncase / (ncase + ncontrol) else NA_real_
  }
  if (is.finite(out$N)) saveRDS(out, cache_path)
  out
}

#' Fetch outcome regional summary statistics from a Pan-UKB tabix file
#'
#' @param rec One-row tibble/data.frame from the MR_df, holding `aws_link` and `aws_link_tabix`.
#' @param ancestry Pan-UKB ancestry code (uppercase: AFR, AMR, CSA, EAS, EUR, MID).
#' @param chr,start,end region coords.
#' @param cache_dir cache directory.
#' @param verbose logical.
#' @return Tibble with `SNP, chr, pos, effect_allele, other_allele, beta, se, eaf, pval`.
#' @keywords internal
#' @export
coloc_fetch_outcome_panukb_region <- function(rec, ancestry, chr, start, end,
                                              cache_dir = ardmr_cache_dir(),
                                              verbose = TRUE) {
  if (!requireNamespace("Rsamtools", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    stop("coloc_fetch_outcome_panukb_region(): requires Bioc Rsamtools/GenomicRanges/IRanges.", call. = FALSE)
  }
  ancestry <- toupper(trimws(as.character(ancestry)))
  url <- .coloc_clean_url(rec$aws_link[[1]])
  idx <- .coloc_clean_url(rec$aws_link_tabix[[1]])
  if (!nzchar(url) || !nzchar(idx)) stop("coloc: missing aws_link / aws_link_tabix on rec", call. = FALSE)

  url_hash <- digest::digest(url, algo = "xxhash64")
  key <- .coloc_region_key("panukb_outcome", url_hash, ancestry, as.character(chr), as.integer(start), as.integer(end))
  cache_path <- .coloc_cache_path(cache_dir, "outcome", key)
  if (file.exists(cache_path)) {
    if (verbose) logger::log_info("coloc: cached panukb outcome region {chr}:{start}-{end}")
    return(readRDS(cache_path))
  }

  hts_cache <- file.path(cache_dir, "hts-cache")
  if (!dir.exists(hts_cache)) dir.create(hts_cache, recursive = TRUE, showWarnings = FALSE)
  old_hts <- Sys.getenv("HTS_CACHE_DIR", unset = NA)
  Sys.setenv(HTS_CACHE_DIR = hts_cache)
  on.exit({ if (is.na(old_hts)) Sys.unsetenv("HTS_CACHE_DIR") else Sys.setenv(HTS_CACHE_DIR = old_hts) }, add = TRUE)

  idx_dir <- file.path(cache_dir, "tabix_index")
  if (!dir.exists(idx_dir)) dir.create(idx_dir, recursive = TRUE, showWarnings = FALSE)
  idx_local <- file.path(idx_dir, paste0(url_hash, ".tbi"))
  if (!file.exists(idx_local)) {
    ok <- tryCatch({ utils::download.file(idx, idx_local, mode = "wb", quiet = !verbose); TRUE },
                   error = function(e) { logger::log_warn("coloc: failed to download tabix index: {conditionMessage(e)}"); FALSE })
    if (!ok || !file.exists(idx_local)) stop("coloc: could not fetch Pan-UKB tabix index", call. = FALSE)
  }

  tf <- Rsamtools::TabixFile(file = url, index = idx_local)
  on.exit(try(Rsamtools::close.TabixFile(tf), silent = TRUE), add = TRUE)
  Rsamtools::open.TabixFile(tf)

  gr_num <- GenomicRanges::GRanges(as.character(chr), IRanges::IRanges(start = as.integer(start), end = as.integer(end)))
  gr_chr <- GenomicRanges::GRanges(paste0("chr", chr), IRanges::IRanges(start = as.integer(start), end = as.integer(end)))
  scan_ok <- FALSE
  res <- tryCatch({ x <- Rsamtools::scanTabix(tf, param = gr_num); scan_ok <- TRUE; x },
                  error = function(e) NULL)
  if (is.null(res) || !length(res) || all(lengths(res) == 0)) {
    res2 <- tryCatch({ x <- Rsamtools::scanTabix(tf, param = gr_chr); scan_ok <- TRUE; x },
                     error = function(e) NULL)
    if (!is.null(res2)) res <- res2
  }
  lines <- if (!is.null(res)) unlist(res, use.names = FALSE) else character()

  empty <- tibble::tibble(SNP = character(), chr = integer(), pos = integer(),
                          effect_allele = character(), other_allele = character(),
                          beta = double(), se = double(), eaf = double(), pval = double())
  if (!length(lines)) {
    if (scan_ok) saveRDS(empty, cache_path)
    return(empty)
  }

  hdr <- tryCatch(Rsamtools::headerTabix(tf), error = function(e) NULL)
  hdr_fields <- character()
  if (!is.null(hdr)) {
    hl <- unlist(hdr, use.names = FALSE)
    cand <- hl[grepl("\t", hl, fixed = TRUE)]
    if (length(cand)) hdr_fields <- strsplit(sub("^#+", "", utils::tail(cand, 1)), "\t", fixed = TRUE)[[1]]
  }
  if (!length(hdr_fields)) stop("coloc: could not read Pan-UKB header for column mapping", call. = FALSE)

  af_idx <- match(paste0("af_", ancestry), hdr_fields)
  if (is.na(af_idx)) af_idx <- match(paste0("af_controls_", ancestry), hdr_fields)
  map <- list(
    chr     = match("chr", hdr_fields),
    pos     = match("pos", hdr_fields),
    ref     = match("ref", hdr_fields),
    alt     = match("alt", hdr_fields),
    af      = af_idx,
    beta    = match(paste0("beta_", ancestry), hdr_fields),
    se      = match(paste0("se_", ancestry), hdr_fields),
    p_log10 = match(paste0("neglog10_pval_", ancestry), hdr_fields)
  )
  if (any(is.na(c(map$chr, map$pos, map$ref, map$alt, map$beta, map$se, map$p_log10)))) {
    stop(sprintf("coloc: Pan-UKB header missing required columns for ancestry %s", ancestry), call. = FALSE)
  }

  panukb_tmp <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE,
                                  fill = TRUE, showProgress = FALSE)
  out <- tibble::tibble(
    chr           = suppressWarnings(as.integer(sub("^chr", "", as.character(panukb_tmp[[map$chr]])))),
    pos           = suppressWarnings(as.integer(panukb_tmp[[map$pos]])),
    effect_allele = toupper(as.character(panukb_tmp[[map$alt]])),
    other_allele  = toupper(as.character(panukb_tmp[[map$ref]])),
    beta          = suppressWarnings(as.numeric(panukb_tmp[[map$beta]])),
    se            = suppressWarnings(as.numeric(panukb_tmp[[map$se]])),
    eaf           = if (!is.na(map$af)) suppressWarnings(as.numeric(panukb_tmp[[map$af]]))
                    else rep(NA_real_, nrow(panukb_tmp)),
    pval          = 10^(-suppressWarnings(as.numeric(panukb_tmp[[map$p_log10]])))
  )
  out$SNP <- paste0(out$chr, ":", out$pos, "_", out$other_allele, "_", out$effect_allele)
  out <- dplyr::relocate(out, "SNP", .before = 1)
  out <- out[!is.na(out$beta) & !is.na(out$se) & out$se > 0, , drop = FALSE]
  saveRDS(out, cache_path)
  out
}

#' Derive Pan-UKB outcome metadata (N, trait type, case/control) from MR_df row.
#' @keywords internal
#' @export
coloc_panukb_outcome_metadata <- function(rec, ancestry) {
  ancestry <- toupper(trimws(as.character(ancestry)))
  pull <- function(nms) {
    for (n in nms) if (n %in% names(rec)) {
      v <- suppressWarnings(as.numeric(rec[[n]][[1]]))
      if (!is.na(v)) return(v)
    }
    NA_real_
  }
  ncase <- pull(c(paste0("n_cases_", ancestry), "n_cases", "ncase"))
  ncontrol <- pull(c(paste0("n_controls_", ancestry), "n_controls", "ncontrol"))
  N <- pull(c(paste0("n_samples_", ancestry), "n_samples", "sample_size"))
  trait_type <- if ("trait_type" %in% names(rec)) tolower(as.character(rec$trait_type[[1]])) else NA_character_
  is_quant <- !is.na(trait_type) && trait_type %in% c("continuous", "biomarkers")
  if (is_quant) {
    if (is.na(N)) N <- ncase
    type <- "quant"
    ncase <- NA_real_; ncontrol <- NA_real_; s <- NA_real_
  } else {
    if (is.na(N) && !is.na(ncase) && !is.na(ncontrol)) N <- ncase + ncontrol
    type <- if (!is.na(ncase) && ncase > 0) "cc" else "quant"
    s <- if (!is.na(ncase) && !is.na(ncontrol) && (ncase + ncontrol) > 0) ncase / (ncase + ncontrol) else NA_real_
  }
  list(N = N, type = type, ncase = ncase, ncontrol = ncontrol, s = s)
}
