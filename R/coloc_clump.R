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
#' Applies p-value thresholding, LD clumping (against the local 1kg.v3 panel),
#' and an F-statistic filter.
#'
#' @export
clump_sumstats_to_ivs <- function(vcf_path, ancestry,
                                  p_threshold = 5e-8, r2 = 0.001, kb = 10000L,
                                  f_threshold = 10,
                                  genome_build = "GRCh37",
                                  cache_dir = ardmr_cache_dir(),
                                  verbose = TRUE) {
  validate_gwasvcf(vcf_path, genome_build)
  if (!requireNamespace("gwasvcf", quietly = TRUE)) {
    stop("clump_sumstats_to_ivs(): requires gwasvcf (mrcieu/gwasvcf).", call. = FALSE)
  }
  if (verbose) logger::log_info("coloc-clump: scanning VCF for hits at p < {p_threshold}")
  vcf <- gwasvcf::query_gwas(vcf_path, pval = p_threshold)
  hits <- gwasvcf::vcf_to_tibble(vcf)
  if (!NROW(hits)) {
    logger::log_warn("coloc-clump: no SNPs at p < {p_threshold}; returning empty.")
    return(.empty_exposure_tbl())
  }
  pval_num <- 10 ^ (-suppressWarnings(as.numeric(hits$LP)))
  std <- tibble::tibble(
    rsid = as.character(hits$rsid),
    pval = pval_num,
    beta = suppressWarnings(as.numeric(hits$ES)),
    se   = suppressWarnings(as.numeric(hits$SE)),
    eaf  = suppressWarnings(as.numeric(hits$AF)),
    chr  = suppressWarnings(as.integer(sub("^chr", "", as.character(hits$seqnames)))),
    pos  = suppressWarnings(as.integer(hits$start)),
    REF  = as.character(hits$REF),
    ALT  = as.character(hits$ALT),
    SS   = suppressWarnings(as.numeric(hits$SS))
  )
  std <- std[is.finite(std$beta) & is.finite(std$se) & std$se > 0 &
             nchar(std$REF) == 1 & nchar(std$ALT) == 1, , drop = FALSE]
  if (!nrow(std)) return(.empty_exposure_tbl())
  ld_pop <- ld_pop_from_ancestry(ancestry)
  ld_ref_dir <- file.path(cache_dir, "ld_reference")
  ref <- ld_reference_checker(pop = ld_pop, ld_ref_dir = ld_ref_dir,
                              verbose = verbose, confirm = "no")
  plink_bin <- tryCatch(.find_plink_bin(NULL), error = function(e) NULL)
  if (verbose) logger::log_info("coloc-clump: ld_clump on {nrow(std)} SNPs (pop={ld_pop}, r2={r2}, kb={kb})")
  clumped <- tryCatch(
    ieugwasr::ld_clump(
      tibble::tibble(rsid = std$rsid, pval = std$pval),
      clump_r2 = r2, clump_kb = kb, plink_bin = plink_bin,
      bfile = ref$bfile_prefix, pop = ld_pop
    ),
    error = function(e) { logger::log_warn("ld_clump failed: {conditionMessage(e)}"); NULL }
  )
  if (is.null(clumped) || !NROW(clumped)) return(.empty_exposure_tbl())
  std <- std[std$rsid %in% clumped$rsid, , drop = FALSE]
  Fstat <- (std$beta / std$se) ^ 2
  std <- std[is.finite(Fstat) & Fstat >= f_threshold, , drop = FALSE]
  if (verbose) logger::log_info("coloc-clump: {nrow(std)} IVs after F>={f_threshold}")
  tag <- sub("\\.vcf(\\.b?gz)?$", "", basename(vcf_path))
  tibble::tibble(
    SNP                    = std$rsid,
    Chr                    = std$chr,
    Pos                    = std$pos,
    effect_allele.exposure = toupper(std$ALT),
    other_allele.exposure  = toupper(std$REF),
    beta.exposure          = std$beta,
    se.exposure            = std$se,
    eaf.exposure           = std$eaf,
    pval.exposure          = std$pval,
    samplesize.exposure    = std$SS,
    id.exposure            = tag,
    exposure               = tag
  )
}
