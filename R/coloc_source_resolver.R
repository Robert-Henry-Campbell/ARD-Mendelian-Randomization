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
#' (since `extract_instruments`/`tophits` bound the data the preprocessor
#' sees) and passes `already_p_filtered = TRUE` (and, when
#' `prefer_server_clump = TRUE`, `already_clumped = TRUE`) downstream.
#' The `snps_*` modes leave both flags FALSE so the preprocessor walks
#' the ladder itself by re-filtering the input tibble.
#'
#' @return list(mode, exposure_snps, coloc_fetcher, coloc_metadata,
#'   clump_opts_resolved, preprocessing_steps).
#'   `preprocessing_steps` is a list (may be empty when
#'   `clump_opts$preprocess = FALSE` or for vcf_only's interim state until
#'   Job 2b refactors `clump_sumstats_to_ivs`).
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

  # ld_pop alignment note: when prefer_server_clump = TRUE for a non-EUR
  # ancestry, OpenGWAS does the clumping with its own (typically EUR) LD
  # panel. Local steps still use ld_pop_from_ancestry(ancestry). Surface the
  # discrepancy so users see it in their run logs.
  if (verbose && has_id && !has_vcf && isTRUE(co$prefer_server_clump)) {
    anc_norm <- toupper(trimws(as.character(ancestry)[1]))
    if (!identical(anc_norm, "EUR") && !identical(anc_norm, "EUROPEAN")) {
      logger::log_info(
        "coloc-source: ancestry='{ancestry}' but prefer_server_clump=TRUE -- server clumping uses its default (typically EUR) LD panel."
      )
    }
  }

  # Helper: invoke preprocess (or bypass when co$preprocess = FALSE) and
  # build the per-mode return list.
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
    mode_opts$already_p_filtered <- already_p_filtered
    mode_opts$already_clumped    <- already_clumped
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
      coloc_fetcher  = .vcf_fetcher(exposure_sumstats, genome_build, verbose),
      coloc_metadata = read_gwasvcf_metadata(exposure_sumstats, genome_build)
    ))
  }

  # Mode: vcf_only
  # Interim state: clump_sumstats_to_ivs still does its own indel/clump/
  # F-stat work. Job 2b will refactor it to delegate to
  # preprocess_exposure_snps. For now, plumb scalars from co (using the
  # strictest p_backoff rung) and return without an extra preprocess pass.
  if (has_vcf) {
    validate_gwasvcf(exposure_sumstats, genome_build)
    snps <- clump_sumstats_to_ivs(
      vcf_path = exposure_sumstats, ancestry = ancestry,
      p_threshold = co$p_backoff[1], r2 = co$r2, kb = co$kb,
      f_threshold = co$f_threshold, genome_build = genome_build,
      cache_dir = cache_dir, verbose = verbose
    )
    return(list(
      mode                = "vcf_only",
      exposure_snps       = tibble::as_tibble(snps),
      coloc_fetcher       = .vcf_fetcher(exposure_sumstats, genome_build, verbose),
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
    already_clumped    = isTRUE(co$prefer_server_clump)
  ))
}

# ---- ieugwasr extraction with p-value backoff at the source layer --------

#' @keywords internal
.ieugwasr_extract_with_backoff <- function(exposure_id, co, verbose = TRUE) {
  rungs <- as.numeric(co$p_backoff)
  for (rung in rungs) {
    ivs <- if (isTRUE(co$prefer_server_clump)) {
      if (verbose) {
        logger::log_info(
          "coloc-source: extract_instruments({exposure_id}) at p<{rung} (server clump)"
        )
      }
      tryCatch(
        TwoSampleMR::extract_instruments(
          outcomes = exposure_id, p1 = rung, clump = TRUE, r2 = co$r2, kb = co$kb
        ),
        error = function(e) {
          logger::log_warn("extract_instruments failed at p<{rung}: {conditionMessage(e)}")
          NULL
        }
      )
    } else {
      if (verbose) {
        logger::log_info(
          "coloc-source: tophits({exposure_id}) at p<{rung} (local clump downstream)"
        )
      }
      th <- safe_tophits(exposure_id, rung)
      if (is.null(th) || !nrow(th)) NULL else .tophits_to_tsmr(th, exposure_id)
    }
    if (!is.null(ivs) && NROW(ivs) > 0L) {
      if (verbose) {
        logger::log_info(
          "coloc-source: {nrow(ivs)} IVs from {exposure_id} at p<{rung}"
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
