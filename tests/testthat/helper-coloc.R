# tests/testthat/helper-coloc.R
#
# Synthetic regional summary statistics for coloc tests. Builds two
# datasets at the same chr:pos grid, with a configurable causal SNP.
# All SNPs are non-palindromic (A/G or C/T) so the harmoniser keeps them.

make_coloc_synthetic_region <- function(chr = 1L, start = 1L, end = 100000L,
                                        n_snps = 50L,
                                        causal_pos = NULL,
                                        causal_z = 5,
                                        noise_z = 0.5,
                                        N = 100000L,
                                        seed = 42L) {
  set.seed(seed)
  pos <- as.integer(round(seq(start, end, length.out = n_snps)))
  ea <- rep(c("A", "C"), length.out = n_snps)
  oa <- rep(c("G", "T"), length.out = n_snps)
  se <- rep(0.01, n_snps)
  z  <- rnorm(n_snps, sd = noise_z)
  if (!is.null(causal_pos)) {
    idx <- which.min(abs(pos - causal_pos))
    z[idx] <- causal_z
  }
  beta <- z * se
  pval <- 2 * stats::pnorm(-abs(z))
  tibble::tibble(
    SNP           = paste0("rs", chr, "_", pos),
    chr           = as.integer(chr),
    pos           = pos,
    effect_allele = ea,
    other_allele  = oa,
    beta          = beta,
    se            = se,
    eaf           = rep(0.30, n_snps),
    pval          = pval
  )
}

make_coloc_synthetic_iv_hdat <- function(snp_id, chr, pos,
                                         beta_exp, se_exp,
                                         beta_out, se_out) {
  tibble::tibble(
    SNP                   = snp_id,
    chr.outcome           = as.integer(chr),
    pos.outcome           = as.integer(pos),
    effect_allele.exposure = "A",
    other_allele.exposure  = "G",
    effect_allele.outcome  = "A",
    other_allele.outcome   = "G",
    beta.exposure         = beta_exp,
    se.exposure           = se_exp,
    pval.exposure         = 2 * stats::pnorm(-abs(beta_exp / se_exp)),
    beta.outcome          = beta_out,
    se.outcome            = se_out,
    pval.outcome          = 2 * stats::pnorm(-abs(beta_out / se_out)),
    eaf.exposure          = 0.30,
    eaf.outcome           = 0.30,
    samplesize.exposure   = 100000L,
    samplesize.outcome    = 50000L,
    mr_keep               = TRUE,
    id.exposure           = "synth",
    id.outcome            = "synth_out",
    exposure              = "Synth",
    outcome               = "SynthOut"
  )
}

make_coloc_synthetic_exposure_snps <- function(snp_ids, chrs, positions) {
  tibble::tibble(
    SNP                    = snp_ids,
    Chr                    = as.integer(chrs),
    Pos                    = as.integer(positions),
    effect_allele.exposure = rep("A", length(snp_ids)),
    other_allele.exposure  = rep("G", length(snp_ids)),
    beta.exposure          = rep(0.05, length(snp_ids)),
    se.exposure            = rep(0.01, length(snp_ids)),
    pval.exposure          = rep(5e-7, length(snp_ids)),
    eaf.exposure           = rep(0.30, length(snp_ids)),
    samplesize.exposure    = rep(100000L, length(snp_ids)),
    id.exposure            = "synth"
  )
}
