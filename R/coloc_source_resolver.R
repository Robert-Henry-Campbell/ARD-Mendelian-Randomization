# R/coloc_source_resolver.R

#' Resolve the (exposure_snps, exposure_sumstats, exposure_id) trio into the
#' inputs the rest of the pipeline expects: a tibble of IVs plus a coloc
#' fetcher closure and metadata list.
#'
#' Modes:
#'   - "snps_only"   : SNPs supplied, no coloc; requires acknowledge_no_coloc.
#'   - "snps_vcf"    : SNPs as IVs; VCF for coloc.
#'   - "vcf_only"    : Derive IVs by clumping VCF; VCF for coloc.
#'   - "ieugwasr"    : Extract IVs from OpenGWAS id; same id for coloc.
#'   - "snps_id"     : SNPs as IVs (already filtered) + id for coloc only.
#'                     Allowed but warns; primarily for ard_compare's
#'                     PheWAS-filtered flow.
#' @keywords internal
.resolve_coloc_source <- function(exposure_snps = NULL,
                                  exposure_sumstats = NULL,
                                  exposure_id = NULL,
                                  ancestry,
                                  genome_build = "GRCh37",
                                  acknowledge_no_coloc = FALSE,
                                  clump_opts = list(),
                                  cache_dir = ardmr_cache_dir(),
                                  verbose = TRUE) {
  has_snps <- !is.null(exposure_snps) && NROW(exposure_snps) > 0
  has_vcf  <- !is.null(exposure_sumstats) && nzchar(as.character(exposure_sumstats))
  has_id   <- !is.null(exposure_id) && nzchar(as.character(exposure_id))

  if (has_id && has_vcf) {
    stop("exposure_id and exposure_sumstats are both provided; pick one.", call. = FALSE)
  }
  if (!has_snps && !has_vcf && !has_id) {
    stop("Provide one of: exposure_snps, exposure_sumstats (VCF), or exposure_id.",
         call. = FALSE)
  }

  defaults <- list(p_threshold = 5e-8, r2 = 0.001, kb = 10000L, f_threshold = 10)
  co <- utils::modifyList(defaults, clump_opts)

  # Mode: snps_only
  if (has_snps && !has_vcf && !has_id) {
    if (!isTRUE(acknowledge_no_coloc)) {
      stop("exposure_snps supplied without exposure_sumstats / exposure_id: ",
           "coloc cannot run. Pass acknowledge_no_coloc = TRUE to proceed without coloc.",
           call. = FALSE)
    }
    return(list(mode = "snps_only", exposure_snps = exposure_snps,
                coloc_fetcher = NULL, coloc_metadata = NULL))
  }

  # Mode: snps_vcf
  if (has_snps && has_vcf) {
    validate_gwasvcf(exposure_sumstats, genome_build)
    return(list(
      mode = "snps_vcf", exposure_snps = exposure_snps,
      coloc_fetcher = .vcf_fetcher(exposure_sumstats, genome_build, verbose),
      coloc_metadata = read_gwasvcf_metadata(exposure_sumstats, genome_build)
    ))
  }

  # Mode: vcf_only
  if (has_vcf) {
    validate_gwasvcf(exposure_sumstats, genome_build)
    snps <- clump_sumstats_to_ivs(
      vcf_path = exposure_sumstats, ancestry = ancestry,
      p_threshold = co$p_threshold, r2 = co$r2, kb = co$kb,
      f_threshold = co$f_threshold, genome_build = genome_build,
      cache_dir = cache_dir, verbose = verbose
    )
    if (!nrow(snps)) {
      stop("clump_sumstats_to_ivs() returned 0 IVs at p<", co$p_threshold,
           ", r2<", co$r2, ", F>=", co$f_threshold, call. = FALSE)
    }
    return(list(
      mode = "vcf_only", exposure_snps = snps,
      coloc_fetcher = .vcf_fetcher(exposure_sumstats, genome_build, verbose),
      coloc_metadata = read_gwasvcf_metadata(exposure_sumstats, genome_build)
    ))
  }

  # Mode: snps_id (advanced — PheWAS-filtered IVs + id for coloc)
  if (has_snps && has_id) {
    if (verbose) {
      logger::log_warn(
        "Both exposure_snps and exposure_id provided: using snps as IVs, id for coloc only."
      )
    }
    return(list(
      mode = "snps_id", exposure_snps = exposure_snps,
      coloc_fetcher = .ieugwasr_fetcher(exposure_id, cache_dir, verbose),
      coloc_metadata = coloc_fetch_metadata(exposure_id, cache_dir, verbose)
    ))
  }

  # Mode: ieugwasr (id only)
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
    stop("Mode 'ieugwasr' requires the TwoSampleMR package.", call. = FALSE)
  }
  if (verbose) logger::log_info("coloc-source: extract_instruments({exposure_id})")
  ivs <- tryCatch(
    TwoSampleMR::extract_instruments(
      outcomes = exposure_id, p1 = co$p_threshold,
      clump = TRUE, r2 = co$r2, kb = co$kb
    ),
    error = function(e) { logger::log_warn("extract_instruments failed: {conditionMessage(e)}"); NULL }
  )
  if (is.null(ivs) || !nrow(ivs)) {
    stop("extract_instruments(", exposure_id, ") returned no IVs at p<",
         co$p_threshold, call. = FALSE)
  }
  if (all(c("beta.exposure", "se.exposure") %in% names(ivs))) {
    Fstat <- (ivs$beta.exposure / ivs$se.exposure) ^ 2
    ivs <- ivs[is.finite(Fstat) & Fstat >= co$f_threshold, , drop = FALSE]
  }
  list(
    mode = "ieugwasr", exposure_snps = tibble::as_tibble(ivs),
    coloc_fetcher = .ieugwasr_fetcher(exposure_id, cache_dir, verbose),
    coloc_metadata = coloc_fetch_metadata(exposure_id, cache_dir, verbose)
  )
}

.vcf_fetcher <- function(vcf_path, genome_build, verbose) {
  force(vcf_path); force(genome_build); force(verbose)
  function(chr, start, end) {
    read_gwasvcf_region(vcf_path, chr, start, end,
                        genome_build = genome_build, verbose = verbose)
  }
}

.ieugwasr_fetcher <- function(exposure_id, cache_dir, verbose) {
  force(exposure_id); force(cache_dir); force(verbose)
  function(chr, start, end) {
    coloc_fetch_exposure_region(exposure_id, chr, start, end,
                                cache_dir = cache_dir, verbose = verbose)
  }
}
