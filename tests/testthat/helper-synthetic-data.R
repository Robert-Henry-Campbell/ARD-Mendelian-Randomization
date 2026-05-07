# tests/testthat/helper-synthetic-data.R
#
# Shared synthetic-input helpers used across multiple test files.
# testthat auto-sources every helper-*.R at the top of every test run,
# so constants and builders defined here are visible everywhere without
# explicit source() calls.
#
# Anything tightly coupled to a single function-under-test belongs in
# its own helper-<function>.R file; this file is for cross-test fixtures.

# ---- Shared RNG seed ----
# Used by every test that calls a function exercising stochastic code
# (e.g. mr_weighted_mode bootstrap, enrichment Monte Carlo permutation).
# Single source of truth so seeds cannot silently drift across tests.
SYNTH_SEED <- 20260501L

# ---- Synthetic exposure-SNP builder ----
# A flexible TwoSampleMR-shaped tibble used by test_preprocess_exposure_snps.R
# and (later) other tests that need a known IV input. Defaults produce one
# clean SNV; override individual columns to construct edge cases (indels,
# palindromes, NA EAF, low INFO, weak F-stat, non-rsid IDs, duplicates).
#
# `n` rows are generated with default values; supply any of the column
# vectors to override. Recycled to length `n`.
make_synthetic_exposure_snps <- function(
    n = 1L,
    SNP                    = sprintf("rs%d", seq_len(n)),
    id.exposure            = "synth",
    exposure               = "synth",
    Chr                    = 1L,
    Pos                    = seq.int(1000L, length.out = n, by = 1000L),
    effect_allele.exposure = rep("A", n),
    other_allele.exposure  = rep("G", n),
    beta.exposure          = rep(0.1, n),
    se.exposure            = rep(0.01, n),
    eaf.exposure           = rep(0.30, n),
    pval.exposure          = rep(1e-10, n),
    samplesize.exposure    = rep(100000L, n),
    INFO                   = NULL
) {
  out <- tibble::tibble(
    SNP                    = SNP,
    id.exposure            = id.exposure,
    exposure               = exposure,
    Chr                    = Chr,
    Pos                    = Pos,
    effect_allele.exposure = effect_allele.exposure,
    other_allele.exposure  = other_allele.exposure,
    beta.exposure          = beta.exposure,
    se.exposure            = se.exposure,
    eaf.exposure           = eaf.exposure,
    pval.exposure          = pval.exposure,
    samplesize.exposure    = samplesize.exposure
  )
  if (!is.null(INFO)) out$INFO <- INFO
  out
}

# ---- Synthetic enrichment input ----
# Chains the mr_business_logic golden fixture: loads its results_df,
# augments with the columns enrichment_business_logic requires
# (ARD_selected, cause_level_1/2/3, pheno_sex). The mr fixture has 3 rows;
# enrichment will auto-filter the Synth_empty row (results_qc_pass = NA),
# leaving 2 rows -- 1 ARD vs 1 non-ARD, both in cause "A". Small input
# but lands on the exact-enumeration permutation path (choose(2,1) = 2),
# so output is deterministic without RNG.

make_synthetic_enrichment_input <- function() {
  fpath <- testthat::test_path("fixtures", "mr_business_logic_golden.rds")
  if (!file.exists(fpath)) {
    stop("Required upstream fixture not found: ", fpath,
         ". Source tests/testthat/fixtures/regenerate_mr_business_logic_golden.R first.",
         call. = FALSE)
  }
  rdf <- readRDS(fpath)$results_df

  # Augment with enrichment-required columns.
  # Row 1 (Synth_3snp):  ARD,     cause_level_1 = "A"
  # Row 2 (Synth_1snp):  non-ARD, cause_level_1 = "A"
  # Row 3 (Synth_empty): NA across the board (filtered out by use_qc_pass)
  rdf$ARD_selected  <- c(TRUE, FALSE, NA)
  rdf$cause_level_1 <- c("A", "A", NA_character_)
  rdf$cause_level_2 <- c("A1", "A1", NA_character_)
  rdf$cause_level_3 <- c("A1a", "A1b", NA_character_)
  rdf$pheno_sex     <- c("both_sex", "both_sex", NA_character_)

  tibble::as_tibble(rdf)
}

# Canonical enrichment_business_logic invocation. Both the regenerator and
# the test call this exact function so the snapshot and live runs stay in
# lockstep. min_nsnp = 1 is required because the upstream mr fixture
# includes a 1-SNP outcome (Synth_1snp) which the enrichment default
# (min_nsnp = 2) would filter out.
run_synthetic_enrichment <- function(input = make_synthetic_enrichment_input()) {
  set.seed(SYNTH_SEED)
  list(
    run_enrichment = run_enrichment(
      input,
      min_nsnp = 1,
      seed = SYNTH_SEED
    ),
    enrichment_global_directional = enrichment_global_directional(
      input,
      min_nsnp = 1,
      seed = SYNTH_SEED
    ),
    enrichment_by_cause_directional = list(
      cause_vs_rest_all = enrichment_by_cause_directional(
        input, level = "cause_level_1",
        compare_mode = "cause_vs_rest_all",
        min_nsnp = 1, seed = SYNTH_SEED
      ),
      ARD_vs_nonARD_within_cause = enrichment_by_cause_directional(
        input, level = "cause_level_1",
        compare_mode = "ARD_vs_nonARD_within_cause",
        min_nsnp = 1, seed = SYNTH_SEED
      ),
      ARD_in_cause_vs_ARD_elsewhere = enrichment_by_cause_directional(
        input, level = "cause_level_1",
        compare_mode = "ARD_in_cause_vs_ARD_elsewhere",
        min_nsnp = 1, seed = SYNTH_SEED
      )
    ),
    beta_contrast_by_cause = beta_contrast_by_cause(
      input, level = "cause_level_1", min_nsnp = 1
    ),
    beta_contrast_global_ARD = beta_contrast_global_ARD(
      input, min_nsnp = 1
    ),
    beta_mean_by_cause = beta_mean_by_cause(
      input, level = "cause_level_1", min_nsnp = 1
    ),
    beta_mean_global = beta_mean_global(
      input, min_nsnp = 1
    )
  )
}
