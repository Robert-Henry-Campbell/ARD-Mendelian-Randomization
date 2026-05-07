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
#'
#' All modes route their IV tibble through [preprocess_exposure_snps()] for
#' a single, unified exposure-SNP preprocessing pipeline. The `ieugwasr`
#' branch handles the p-value backoff ladder at the source-fetch layer
#' (fetching unclumped via `extract_instruments(clump = FALSE)` or
#' `tophits` per rung) and passes `already_p_filtered = TRUE,
#' already_clumped = FALSE` to the preprocessor. Local LD clumping then
#' runs against the ancestry-aware 1kg.v3 panel inside
#' `preprocess_exposure_snps`. The `snps_*` modes leave both flags
#' FALSE so the preprocessor walks the ladder itself by re-filtering
#' the input tibble.
#'
#' @return list(mode, exposure_snps, coloc_fetcher, coloc_metadata,
#'   clump_opts_resolved, preprocessing_steps).
#'   `preprocessing_steps` is a list (may be empty when
#'   `clump_opts$preprocess = FALSE` or for `vcf_only` whose
#'   `clump_sumstats_to_ivs` runs the preprocessor internally and
#'   returns just the SNPs).
#' @keywords internal
.resolve_coloc_source <- function(exposure_snps = NULL,
                                  exposure_sumstats = NULL,
                                  exposure_id = NULL,
                                  ancestry,
                                  genome_build = "GRCh37",
                                  acknowledge_no_coloc = FALSE,
                                  clump_opts = list(),
                                  cache_dir = ardmr_cache_dir(),
                                  confirm = "ask",
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

  # Single source of truth for defaults + p_threshold deprecation.
  co <- .resolve_clump_opts(clump_opts)

  # Clumping always runs locally against the ancestry-aware 1kg.v3
  # panel (Job 4). The previous prefer_server_clump branch was removed
  # because OpenGWAS's server-side clumping uses its default (EUR) LD
  # panel regardless of the requested ancestry.

  # Helper: invoke preprocess (or bypass when co$preprocess = FALSE) and
  # build the per-mode return list.
  #
  # `already_p_filtered` / `already_clumped` are RESOLVER-SIDE hints
  # (facts about what the resolver itself just did). They are merged
  # with the user's clump_opts via OR: skip the corresponding local
  # step if EITHER the resolver did it (per-mode hint TRUE) OR the
  # caller signalled it was done upstream (clump_opts$already_X = TRUE).
  # Previously these args UNCONDITIONALLY overwrote the user's values,
  # so ard_compare's `already_clumped = TRUE` was silently discarded
  # for snps_id mode and the unified preprocessor re-clumped already-
  # clumped SNPs.
  finalize <- function(mode, snps, coloc_fetcher, coloc_metadata,
                       already_p_filtered = FALSE, already_clumped = FALSE) {
    if (!isTRUE(co$preprocess)) {
      return(list(
        mode                = mode,
        exposure_snps       = tibble::as_tibble(snps),
        coloc_fetcher       = coloc_fetcher,
        coloc_metadata      = coloc_metadata,
        clump_opts_resolved = co,
        preprocessing_steps = list()
      ))
    }
    mode_opts <- co
    mode_opts$already_p_filtered <-
      isTRUE(co$already_p_filtered) || isTRUE(already_p_filtered)
    mode_opts$already_clumped    <-
      isTRUE(co$already_clumped)    || isTRUE(already_clumped)
    pp <- preprocess_exposure_snps(
      snps_tbl   = snps,
      clump_opts = mode_opts,
      ancestry   = ancestry,
      cache_dir  = cache_dir,
      confirm    = confirm,
      verbose    = verbose
    )
    list(
      mode                = mode,
      exposure_snps       = pp$snps,
      coloc_fetcher       = coloc_fetcher,
      coloc_metadata      = coloc_metadata,
      clump_opts_resolved = pp$resolved_opts,
      preprocessing_steps = pp$steps
    )
  }

  # Mode: snps_only
  if (has_snps && !has_vcf && !has_id) {
    if (!isTRUE(acknowledge_no_coloc)) {
      stop("exposure_snps supplied without exposure_sumstats / exposure_id: ",
           "coloc cannot run. Pass acknowledge_no_coloc = TRUE to proceed without coloc.",
           call. = FALSE)
    }
    return(finalize("snps_only", exposure_snps, NULL, NULL))
  }

  # Mode: snps_vcf
  if (has_snps && has_vcf) {
    validate_gwasvcf(exposure_sumstats, genome_build)
    return(finalize(
      "snps_vcf", exposure_snps,
      coloc_fetcher  = .vcf_fetcher(exposure_sumstats, genome_build,
                                    co$info_min, verbose),
      coloc_metadata = read_gwasvcf_metadata(exposure_sumstats, genome_build)
    ))
  }

  # Mode: vcf_only
  # clump_sumstats_to_ivs (Job 2b refactor) now delegates to
  # preprocess_exposure_snps internally. Pass the resolved clump_opts so
  # the unified pipeline (with the same defaults this resolver uses) runs
  # against the VCF-derived IVs. The vcf_only branch does NOT layer an
  # additional preprocess pass on top -- preprocessing_steps is empty
  # here because it lives inside clump_sumstats_to_ivs.
  if (has_vcf) {
    validate_gwasvcf(exposure_sumstats, genome_build)
    snps <- clump_sumstats_to_ivs(
      vcf_path     = exposure_sumstats,
      ancestry     = ancestry,
      clump_opts   = co,
      genome_build = genome_build,
      cache_dir    = cache_dir,
      confirm      = confirm,
      verbose      = verbose
    )
    return(list(
      mode                = "vcf_only",
      exposure_snps       = tibble::as_tibble(snps),
      coloc_fetcher       = .vcf_fetcher(exposure_sumstats, genome_build,
                                         co$info_min, verbose),
      coloc_metadata      = read_gwasvcf_metadata(exposure_sumstats, genome_build),
      clump_opts_resolved = co,
      preprocessing_steps = list()
    ))
  }

  # Mode: snps_id (advanced -- PheWAS-filtered IVs + id for coloc)
  if (has_snps && has_id) {
    if (verbose) {
      logger::log_warn(
        "Both exposure_snps and exposure_id provided: using snps as IVs, id for coloc only."
      )
    }
    return(finalize(
      "snps_id", exposure_snps,
      coloc_fetcher  = .ieugwasr_fetcher(exposure_id, cache_dir, verbose),
      coloc_metadata = coloc_fetch_metadata(exposure_id, cache_dir, verbose)
    ))
  }

  # Mode: ieugwasr (id only)
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
    stop("Mode 'ieugwasr' requires the TwoSampleMR package.", call. = FALSE)
  }
  ivs <- .ieugwasr_extract_with_backoff(
    exposure_id = exposure_id, co = co, verbose = verbose
  )
  if (is.null(ivs) || !nrow(ivs)) ivs <- .empty_exposure_tbl()
  return(finalize(
    "ieugwasr", ivs,
    coloc_fetcher      = .ieugwasr_fetcher(exposure_id, cache_dir, verbose),
    coloc_metadata     = coloc_fetch_metadata(exposure_id, cache_dir, verbose),
    already_p_filtered = TRUE,
    already_clumped    = FALSE  # local clump always runs in preprocess
  ))
}

# ---- ieugwasr extraction with p-value backoff at the source layer --------

#' Fetch unclumped instruments from OpenGWAS, walking the p_backoff
#' ladder strictest-first.
#'
#' Always fetches WITHOUT server-side clumping so the unified
#' preprocessor can clump locally against the ancestry-aware 1kg.v3
#' panel. Falls back from `extract_instruments(clump = FALSE)` to
#' `safe_tophits()` when the former fails.
#' @keywords internal
.ieugwasr_extract_with_backoff <- function(exposure_id, co, verbose = TRUE) {
  rungs <- as.numeric(co$p_backoff)
  for (rung in rungs) {
    if (verbose) {
      logger::log_info(
        "coloc-source: extract_instruments({exposure_id}) at p<{rung} (clump=FALSE; local clump downstream)"
      )
    }
    ivs <- tryCatch(
      TwoSampleMR::extract_instruments(
        outcomes = exposure_id, p1 = rung, clump = FALSE
      ),
      error = function(e) {
        logger::log_warn(
          "extract_instruments failed at p<{rung}: {conditionMessage(e)}; falling back to tophits()"
        )
        NULL
      }
    )
    if (is.null(ivs) || !NROW(ivs)) {
      th <- safe_tophits(exposure_id, rung)
      if (!is.null(th) && nrow(th)) ivs <- .tophits_to_tsmr(th, exposure_id)
    }
    if (!is.null(ivs) && NROW(ivs) > 0L) {
      if (verbose) {
        logger::log_info(
          "coloc-source: {nrow(ivs)} IVs from {exposure_id} at p<{rung} (pre-clump)"
        )
      }
      return(tibble::as_tibble(ivs))
    }
  }
  if (verbose) {
    logger::log_warn(
      "coloc-source: no IVs from {exposure_id} across p_backoff = ({paste(rungs, collapse=', ')})"
    )
  }
  NULL
}

#' Convert an `ieugwasr::tophits` tibble to TwoSampleMR-shaped exposure SNPs.
#' @keywords internal
.tophits_to_tsmr <- function(th, exposure_id) {
  pick <- function(nms, default = NA) {
    hit <- intersect(nms, names(th))
    if (!length(hit)) return(rep(default, nrow(th)))
    th[[hit[1]]]
  }
  tibble::tibble(
    SNP                    = as.character(pick(c("rsid", "SNP"))),
    id.exposure            = as.character(exposure_id),
    exposure               = as.character(exposure_id),
    effect_allele.exposure = toupper(as.character(pick(c("ea", "effect_allele")))),
    other_allele.exposure  = toupper(as.character(pick(c("nea", "other_allele")))),
    beta.exposure          = suppressWarnings(as.numeric(pick(c("beta", "es")))),
    se.exposure            = suppressWarnings(as.numeric(pick(c("se")))),
    eaf.exposure           = suppressWarnings(as.numeric(pick(c("eaf", "af")))),
    pval.exposure          = suppressWarnings(as.numeric(pick(c("p", "pval")))),
    samplesize.exposure    = suppressWarnings(as.numeric(pick(c("n", "samplesize"))))
  )
}

.vcf_fetcher <- function(vcf_path, genome_build, info_min, verbose) {
  force(vcf_path); force(genome_build); force(info_min); force(verbose)
  effective_info_min <- if (is.null(info_min)) 0.8 else info_min
  function(chr, start, end) {
    read_gwasvcf_region(vcf_path, chr, start, end,
                        genome_build = genome_build,
                        info_min = effective_info_min,
                        verbose = verbose)
  }
}

.ieugwasr_fetcher <- function(exposure_id, cache_dir, verbose) {
  force(exposure_id); force(cache_dir); force(verbose)
  function(chr, start, end) {
    coloc_fetch_exposure_region(exposure_id, chr, start, end,
                                cache_dir = cache_dir, verbose = verbose)
  }
}
