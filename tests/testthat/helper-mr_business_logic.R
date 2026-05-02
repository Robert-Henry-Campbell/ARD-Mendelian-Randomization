# tests/testthat/helper-mr_business_logic.R
#
# Builders for synthetic inputs to mr_business_logic().
# Used by both the golden-fixture regenerator (data-raw/) and the
# regression test (test_mr_business_logic.R).
#
# All values are hand-picked to:
#   * survive TwoSampleMR::harmonise_data(action = 2) without palindrome drops
#     (no A/T or G/C effect/other allele pairs at any SNP)
#   * yield finite, non-degenerate IVW / Egger / weighted-median / weighted-mode
#     estimates
#   * produce a mix of pass/fail sensitivity checks

make_synthetic_exposure <- function() {
  tibble::tibble(
    SNP                    = c("rs1", "rs2", "rs3"),
    effect_allele.exposure = c("A", "C", "G"),
    other_allele.exposure  = c("G", "T", "A"),
    beta.exposure          = c(0.05, -0.03, 0.12),
    se.exposure            = c(0.010, 0.015, 0.008),
    pval.exposure          = c(5e-7, 4.5e-2, 1e-30),
    eaf.exposure           = c(0.35, 0.42, 0.28),
    samplesize.exposure    = c(100000L, 100000L, 100000L),
    id.exposure            = "synth_exposure",
    exposure               = "Synthetic exposure"
  )
}

# Build a TwoSampleMR-formatted outcome tibble for a given subset of the
# synthetic exposure SNPs, with provided per-SNP outcome betas/SEs.
make_synthetic_outcome <- function(snps, beta_out, se_out,
                                   pval_out = NULL, label = "Synth_outcome") {
  stopifnot(length(snps) == length(beta_out),
            length(snps) == length(se_out))
  all_alleles <- list(
    rs1 = list(ea = "A", oa = "G", chr = 1L, pos = 100000L, eaf = 0.35),
    rs2 = list(ea = "C", oa = "T", chr = 1L, pos = 200000L, eaf = 0.42),
    rs3 = list(ea = "G", oa = "A", chr = 2L, pos = 150000L, eaf = 0.28)
  )
  if (!all(snps %in% names(all_alleles))) {
    stop("Unknown SNP requested: ",
         paste(setdiff(snps, names(all_alleles)), collapse = ", "))
  }
  if (is.null(pval_out)) {
    z <- beta_out / se_out
    pval_out <- 2 * stats::pnorm(-abs(z))
  }
  tibble::tibble(
    SNP                  = snps,
    chr.outcome          = vapply(snps, function(s) all_alleles[[s]]$chr, integer(1)),
    pos.outcome          = vapply(snps, function(s) all_alleles[[s]]$pos, integer(1)),
    effect_allele.outcome = vapply(snps, function(s) all_alleles[[s]]$ea, character(1)),
    other_allele.outcome  = vapply(snps, function(s) all_alleles[[s]]$oa, character(1)),
    beta.outcome         = beta_out,
    se.outcome           = se_out,
    eaf.outcome          = vapply(snps, function(s) all_alleles[[s]]$eaf, numeric(1)),
    pval.outcome         = pval_out,
    samplesize.outcome   = rep(50000L, length(snps)),
    outcome              = label,
    Phenotype            = label,
    id.outcome           = "synth_outcome",
    mr_keep.outcome      = TRUE
  )
}

# Top-level builder: returns a list with $exposure_snps and $MR_df.
# MR_df has 3 rows exercising the main code paths in mr_business_logic():
#   row 1: 3-SNP outcome  -> full IVW/Egger/WMed/WMode + heterogeneity + LOO + Steiger
#   row 2: 1-SNP outcome  -> Wald-ratio branch (no heterogeneity, no LOO)
#   row 3: empty outcome  -> early-skip branch
make_synthetic_mr_inputs <- function() {
  exposure_snps <- make_synthetic_exposure()

  out_3snp <- make_synthetic_outcome(
    snps     = c("rs1", "rs2", "rs3"),
    beta_out = c( 0.080, -0.060,  0.150),
    se_out   = c( 0.012,  0.018,  0.010),
    label    = "Synth_3snp"
  )

  out_1snp <- make_synthetic_outcome(
    snps     = "rs1",
    beta_out = 0.090,
    se_out   = 0.014,
    label    = "Synth_1snp"
  )

  out_empty <- tibble::tibble()

  MR_df <- tibble::tibble(
    description  = c("Synth_3snp", "Synth_1snp", "Synth_empty"),
    outcome_snps = list(out_3snp, out_1snp, out_empty)
  )

  list(exposure_snps = exposure_snps, MR_df = MR_df)
}

# Canonical mr_business_logic() invocation used by both regenerator and test.
SYNTH_SENSITIVITY_ENABLED <- c(
  "egger_intercept", "egger_slope_agreement",
  "weighted_median", "weighted_mode",
  "steiger_direction", "leave_one_out",
  "ivw_Q", "ivw_I2"
)
SYNTH_SENSITIVITY_PASS_MIN <- 4L
SYNTH_SEED <- 20260501L

run_synthetic_mr_business_logic <- function(inputs = make_synthetic_mr_inputs()) {
  set.seed(SYNTH_SEED)
  mr_business_logic(
    MR_df                = inputs$MR_df,
    exposure_snps        = inputs$exposure_snps,
    sensitivity_enabled  = SYNTH_SENSITIVITY_ENABLED,
    sensitivity_pass_min = SYNTH_SENSITIVITY_PASS_MIN,
    scatterplot          = FALSE,
    snpforestplot        = FALSE,
    leaveoneoutplot      = FALSE,
    plot_output_dir      = "",
    cache_dir            = tempdir(),
    verbose              = FALSE,
    test                 = FALSE
  )
}
