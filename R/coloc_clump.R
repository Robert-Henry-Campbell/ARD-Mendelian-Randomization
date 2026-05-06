# R/coloc_clump.R

#' @keywords internal
.empty_exposure_tbl <- function() {
  tibble::tibble(
    SNP = character(), Chr = integer(), Pos = integer(),
    effect_allele.exposure = character(), other_allele.exposure = character(),
    beta.exposure = double(), se.exposure = double(),
    eaf.exposure = double(), pval.exposure = double(),
    samplesize.exposure = double(),
    id.exposure = character(), exposure = character()
  )
}

#' Derive TwoSampleMR-formatted exposure IVs from a GWAS-VCF.
#'
#' Queries the VCF for hits at the loosest p-value rung in `p_backoff`,
#' converts to TwoSampleMR shape, and routes through
#' [preprocess_exposure_snps()] for the unified pipeline (p-value backoff,
#' LD clumping, rsid-validation/dedupe, indel/palindromic, F-stat, MAF,
#' INFO). The clumping always runs locally against the 1kg.v3 panel;
#' `already_clumped` and `already_p_filtered` are forced to FALSE because
#' VCF-derived IVs come straight from the file with no upstream filtering.
#'
#' @param vcf_path Path to a tabix-indexed GWAS-VCF (MRC IEU spec).
#' @param ancestry Ancestry code (used to derive the LD population).
#' @param p_threshold,r2,kb,f_threshold Legacy scalar args kept for
#'   backward compat. When `clump_opts` is supplied, its values override
#'   these scalars.
#' @param genome_build Genome build of the VCF (only "GRCh37" supported).
#' @param cache_dir Cache root.
#' @param clump_opts Named list overriding the scalar args; same fields as
#'   documented at `?run_phenome_mr` (`p_backoff`, `r2`, `kb`,
#'   `f_threshold`, `maf_min`, `info_min`, `drop_indels`,
#'   `drop_palindromic`).
#' @param confirm One of `"ask"`, `"yes"`, `"no"`; passed to
#'   [ld_reference_checker()] when the LD panel must be installed for
#'   clumping.
#' @param verbose Logical.
#' @export
clump_sumstats_to_ivs <- function(vcf_path, ancestry,
                                  p_threshold = 5e-8, r2 = 0.001, kb = 10000L,
                                  f_threshold = 10,
                                  genome_build = "GRCh37",
                                  cache_dir = ardmr_cache_dir(),
                                  clump_opts = list(),
                                  confirm = "ask",
                                  verbose = TRUE) {
  validate_gwasvcf(vcf_path, genome_build)
  if (!requireNamespace("gwasvcf", quietly = TRUE)) {
    stop("clump_sumstats_to_ivs(): requires gwasvcf (mrcieu/gwasvcf).", call. = FALSE)
  }

  # Build effective clump_opts: scalar args provide defaults, clump_opts
  # entries override. Then run through .resolve_clump_opts() to merge in
  # the canonical preprocess defaults (maf_min, info_min, drop_indels, ...).
  effective_opts <- utils::modifyList(
    list(p_backoff = c(p_threshold), r2 = r2, kb = kb, f_threshold = f_threshold),
    clump_opts
  )
  co <- .resolve_clump_opts(effective_opts)

  # VCF-derived IVs come from a fresh query, so the unified pipeline must
  # NOT skip clumping or p-filter regardless of what the caller passed.
  co$already_clumped    <- FALSE
  co$already_p_filtered <- FALSE

  # Query at the loosest rung so the preprocessor has headroom to walk the
  # backoff ladder locally.
  p_loosest <- max(as.numeric(co$p_backoff))
  if (verbose) logger::log_info("coloc-clump: scanning VCF for hits at p < {p_loosest}")
  vcf <- gwasvcf::query_gwas(vcf_path, pval = p_loosest)
  hits <- gwasvcf::vcf_to_tibble(vcf)
  if (!NROW(hits)) {
    logger::log_warn("coloc-clump: no SNPs at p < {p_loosest}; returning empty.")
    return(.empty_exposure_tbl())
  }

  # Convert to TwoSampleMR shape FIRST so preprocess_exposure_snps sees
  # the canonical column names. The unified preprocessor then handles
  # indel drop, clumping, F-stat, palindromic, MAF, INFO.
  pval_num <- 10 ^ (-suppressWarnings(as.numeric(hits$LP)))
  tag <- sub("\\.vcf(\\.b?gz)?$", "", basename(vcf_path))
  snps_tbl <- tibble::tibble(
    SNP                    = as.character(hits$rsid),
    Chr                    = suppressWarnings(as.integer(sub("^chr", "", as.character(hits$seqnames)))),
    Pos                    = suppressWarnings(as.integer(hits$start)),
    effect_allele.exposure = toupper(as.character(hits$ALT)),
    other_allele.exposure  = toupper(as.character(hits$REF)),
    beta.exposure          = suppressWarnings(as.numeric(hits$ES)),
    se.exposure            = suppressWarnings(as.numeric(hits$SE)),
    eaf.exposure           = suppressWarnings(as.numeric(hits$AF)),
    pval.exposure          = pval_num,
    samplesize.exposure    = suppressWarnings(as.numeric(hits$SS)),
    id.exposure            = tag,
    exposure               = tag
  )
  # Propagate any INFO-shaped column so preprocess_exposure_snps' INFO
  # step can see it. The preprocessor accepts any of INFO/info/Rsq/R2/SI;
  # gwasvcf::vcf_to_tibble typically yields `SI` (the FORMAT field).
  info_candidates <- c("INFO", "info", "Rsq", "R2", "SI")
  matched_info <- info_candidates[info_candidates %in% names(hits)]
  if (length(matched_info)) {
    nm <- matched_info[1]
    snps_tbl[[nm]] <- suppressWarnings(as.numeric(hits[[nm]]))
  }
  # Drop rows without finite beta/se -- preprocess_exposure_snps' F-stat
  # step would also drop these, but doing it here keeps the clump step's
  # input universe sane (and gives the F-stat step a clean dataset).
  ok <- is.finite(snps_tbl$beta.exposure) & is.finite(snps_tbl$se.exposure) &
        snps_tbl$se.exposure > 0
  snps_tbl <- snps_tbl[ok, , drop = FALSE]
  if (!nrow(snps_tbl)) return(.empty_exposure_tbl())

  pp <- preprocess_exposure_snps(
    snps_tbl   = snps_tbl,
    clump_opts = co,
    ancestry   = ancestry,
    cache_dir  = cache_dir,
    confirm    = confirm,
    verbose    = verbose
  )
  pp$snps
}
