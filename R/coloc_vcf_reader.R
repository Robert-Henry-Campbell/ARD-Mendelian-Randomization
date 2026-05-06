# R/coloc_vcf_reader.R

#' @keywords internal
.assert_grch37 <- function(declared) {
  if (!identical(toupper(as.character(declared)), "GRCH37")) {
    stop("ardmr only supports GRCh37 sumstats; declared genome_build='", declared,
         "'. Lift coordinates to GRCh37 before re-running.", call. = FALSE)
  }
}

#' Validate a GWAS-VCF for use as the coloc exposure source.
#' @export
validate_gwasvcf <- function(vcf_path, genome_build = "GRCh37") {
  .assert_grch37(genome_build)
  if (!file.exists(vcf_path)) stop("VCF not found: ", vcf_path, call. = FALSE)
  if (!file.exists(paste0(vcf_path, ".tbi"))) {
    stop("Tabix index missing for ", vcf_path,
         "; create with Rsamtools::indexTabix() or `tabix -p vcf <file>`.", call. = FALSE)
  }
  hdr <- tryCatch(Rsamtools::headerTabix(vcf_path), error = function(e) NULL)
  if (is.null(hdr) || !length(hdr$header)) {
    stop("VCF header unreadable: ", vcf_path, call. = FALSE)
  }
  fmt <- grep("^##FORMAT=<ID=", hdr$header, value = TRUE)
  keys <- sub("^##FORMAT=<ID=([^,]+),.*$", "\\1", fmt)
  need <- c("ES", "SE")
  miss <- setdiff(need, keys)
  if (length(miss)) {
    stop("VCF missing required FORMAT fields: ", paste(miss, collapse = ", "),
         "; not GWAS-VCF.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.read_gwasvcf_lines <- function(vcf_path, chr, start, end) {
  gr_num <- GenomicRanges::GRanges(as.character(chr),
                                   IRanges::IRanges(as.integer(start), as.integer(end)))
  res <- tryCatch(Rsamtools::scanTabix(vcf_path, param = gr_num), error = function(e) NULL)
  if (is.null(res) || !length(res) || all(lengths(res) == 0)) {
    gr_chr <- GenomicRanges::GRanges(paste0("chr", chr),
                                     IRanges::IRanges(as.integer(start), as.integer(end)))
    res <- tryCatch(Rsamtools::scanTabix(vcf_path, param = gr_chr), error = function(e) NULL)
  }
  if (is.null(res)) character() else unlist(res, use.names = FALSE)
}

#' @keywords internal
.empty_gwasvcf_tbl <- function() {
  tibble::tibble(SNP = character(), chr = integer(), pos = integer(),
                 effect_allele = character(), other_allele = character(),
                 beta = double(), se = double(), eaf = double(), pval = double(),
                 n = double(), ncase = double())
}

#' @keywords internal
.parse_gwasvcf_lines <- function(lines) {
  if (!length(lines)) return(.empty_gwasvcf_tbl())
  parts <- strsplit(lines, "\t", fixed = TRUE)
  ncol <- vapply(parts, length, integer(1))
  parts <- parts[ncol >= 10]
  if (!length(parts)) return(.empty_gwasvcf_tbl())
  m <- do.call(rbind, lapply(parts, function(x) x[1:10]))
  fmt_keys <- strsplit(m[, 9], ":", fixed = TRUE)
  fmt_vals <- strsplit(m[, 10], ":", fixed = TRUE)
  field <- function(key) {
    vapply(seq_along(fmt_keys), function(i) {
      idx <- match(key, fmt_keys[[i]])
      if (is.na(idx) || idx > length(fmt_vals[[i]])) NA_character_
      else fmt_vals[[i]][idx]
    }, character(1))
  }
  ES <- field("ES"); SE <- field("SE"); LP <- field("LP")
  AF <- field("AF"); SS <- field("SS"); NC <- field("NC"); SI <- field("SI")
  tibble::tibble(
    SNP           = m[, 3],
    chr           = m[, 1],
    pos           = suppressWarnings(as.integer(m[, 2])),
    REF           = m[, 4],
    ALT           = m[, 5],
    INFO          = m[, 8],
    SI            = suppressWarnings(as.numeric(SI)),
    ES = ES, SE = SE, LP = LP, AF = AF, SS = SS, NC = NC
  )
}

#' @keywords internal
.split_multiallelic <- function(df) {
  if (!nrow(df)) return(df)
  alt_lists <- strsplit(df$ALT, ",", fixed = TRUE)
  n_alts <- lengths(alt_lists)
  if (all(n_alts == 1)) return(df)
  idx <- rep(seq_len(nrow(df)), times = n_alts)
  out <- df[idx, , drop = FALSE]
  out$ALT <- unlist(alt_lists, use.names = FALSE)
  alt_pos <- unlist(lapply(n_alts, seq_len), use.names = FALSE)
  for (col in c("ES", "SE", "LP", "AF")) {
    if (col %in% names(df)) {
      v_lists <- strsplit(as.character(df[[col]]), ",", fixed = TRUE)
      v_lists <- v_lists[idx]
      out[[col]] <- mapply(function(v, k) if (length(v) >= k) v[k] else NA_character_,
                           v_lists, alt_pos, USE.NAMES = FALSE)
    }
  }
  out
}

#' Read a GWAS-VCF region and return a canonical exposure-region tibble.
#' @export
read_gwasvcf_region <- function(vcf_path, chr, start, end,
                                genome_build = "GRCh37",
                                drop_indels = TRUE,
                                info_min = 0.8,
                                verbose = TRUE) {
  .assert_grch37(genome_build)
  raw <- .parse_gwasvcf_lines(.read_gwasvcf_lines(vcf_path, chr, start, end))
  if (!nrow(raw)) return(.empty_gwasvcf_tbl())
  raw <- .split_multiallelic(raw)
  if (drop_indels) {
    is_indel <- nchar(raw$REF) > 1 | nchar(raw$ALT) > 1 | raw$REF == "." | raw$ALT == "."
    n_drop <- sum(is_indel, na.rm = TRUE)
    if (n_drop && verbose) logger::log_info("coloc-vcf: dropped {n_drop} indel/non-SNV rows")
    raw <- raw[!is_indel, , drop = FALSE]
  }
  out <- tibble::tibble(
    SNP           = raw$SNP,
    chr           = suppressWarnings(as.integer(sub("^chr", "", raw$chr))),
    pos           = raw$pos,
    effect_allele = toupper(raw$ALT),
    other_allele  = toupper(raw$REF),
    beta          = suppressWarnings(as.numeric(raw$ES)),
    se            = suppressWarnings(as.numeric(raw$SE)),
    eaf           = suppressWarnings(as.numeric(raw$AF)),
    pval          = 10 ^ (-suppressWarnings(as.numeric(raw$LP))),
    n             = suppressWarnings(as.numeric(raw$SS)),
    ncase         = suppressWarnings(as.numeric(raw$NC))
  )
  out <- out[is.finite(out$beta) & is.finite(out$se) & out$se > 0, , drop = FALSE]
  if ("SI" %in% names(raw) && any(is.finite(raw$SI))) {
    si <- suppressWarnings(as.numeric(raw$SI))[is.finite(out$beta)]
    keep <- !is.finite(si) | si >= info_min
    n_low <- sum(!keep, na.rm = TRUE)
    if (n_low && verbose) logger::log_info("coloc-vcf: dropped {n_low} low-INFO (<{info_min}) variants")
    out <- out[keep, , drop = FALSE]
  }
  big <- abs(out$beta) > 10
  if (any(big, na.rm = TRUE)) {
    logger::log_warn("coloc-vcf: {sum(big, na.rm=TRUE)} SNPs have |beta|>10 (raw OR not log(OR)?)")
  }
  out
}

#' Read aggregate metadata (N, ncase, type) from a GWAS-VCF.
#' @export
read_gwasvcf_metadata <- function(vcf_path, genome_build = "GRCh37") {
  validate_gwasvcf(vcf_path, genome_build)
  hdr <- Rsamtools::headerTabix(vcf_path)$header
  seqs <- grep("^##contig=<ID=", hdr, value = TRUE)
  chrs <- sub("^##contig=<ID=([^,>]+).*$", "\\1", seqs)
  chrs <- chrs[grepl("^(chr)?(1[0-9]?|2[0-2]?|[3-9])$", chrs)]
  res <- character()
  for (c0 in head(chrs, 5)) {
    res <- .read_gwasvcf_lines(vcf_path, sub("^chr", "", c0), 1L, 250000000L)
    if (length(res)) break
  }
  raw <- .parse_gwasvcf_lines(head(res, 200))
  if (!nrow(raw)) {
    return(list(N = NA_real_, type = "quant", ncase = NA_real_, ncontrol = NA_real_, s = NA_real_))
  }
  N  <- stats::median(suppressWarnings(as.numeric(raw$SS)), na.rm = TRUE)
  NC <- stats::median(suppressWarnings(as.numeric(raw$NC)), na.rm = TRUE)
  is_cc <- is.finite(NC) && NC > 0
  list(
    N        = if (is.finite(N)) N else NA_real_,
    type     = if (is_cc) "cc" else "quant",
    ncase    = if (is_cc) NC else NA_real_,
    ncontrol = if (is_cc && is.finite(N)) N - NC else NA_real_,
    s        = if (is_cc && is.finite(N) && N > 0) NC / N else NA_real_
  )
}
