# =============================================================================
# ardmr template: genome-wide MR  WITHOUT  coloc (no GWAS-VCF)
# =============================================================================
#
# Use this when:
#   * You're running phenome-wide MR for a polygenic exposure (e.g. BMI, LDL-C)
#     using a curated set of genome-wide-significant lead SNPs.
#   * You don't have an exposure GWAS-VCF (or don't want to run coloc).
#
# No package patch needed for this variant -- ard_compare() automatically
# sets acknowledge_no_coloc=TRUE when no exposure_sumstats / exposure_id
# is supplied per group.
#
# Catalog constraint (R/ard_compare.R:247-251):
#   sex = "both" forces the panukb catalog; sex = "male"/"female" forces the
#   neale catalog. They cannot be mixed in one ard_compare() call. So a
#   "both, female, male" comparison needs TWO ard_compare() calls.
# -----------------------------------------------------------------------------

setwd("PATH/TO/ardmr")          # path to package source
devtools::load_all()
# devtools::test()              # uncomment to verify nothing's broken first

Sys.setenv(ARDMR_CACHE_DIR = "PATH/TO/CACHE")     # tens-of-GB scratch dir
Sys.setenv(OPENGWAS_JWT    = "YOUR_JWT_HERE")     # only needed if any group uses ieu_id

# ---- 1. Inputs --------------------------------------------------------------
# CSV must contain TwoSampleMR-style columns required by assert_exposure():
#   SNP, effect_allele.exposure, other_allele.exposure,
#   beta.exposure, se.exposure, pval.exposure, id.exposure
snps_both   <- read.csv("PATH/TO/exposure_snps_both.csv")
snps_female <- read.csv("PATH/TO/exposure_snps_female.csv")
snps_male   <- read.csv("PATH/TO/exposure_snps_male.csv")

# Genome-wide MR clump opts: trust the upstream LD-clump + p-filter that
# produced your CSV (matches ard_compare's default). Apply on-top
# data-quality filters here.
clump_opts_gw <- list(
  already_clumped    = TRUE,
  already_p_filtered = TRUE,
  f_threshold        = 10,
  maf_min            = 0.01,
  drop_indels        = TRUE,
  drop_palindromic   = FALSE          # TwoSampleMR's harmonise handles these
)
# If your CSV is NOT pre-clumped, set already_clumped=FALSE,
# already_p_filtered=FALSE, p_backoff=c(5e-8), r2=0.001, kb=10000.

# Sensitivity panel WITHOUT coloc (8 of 9 checks).
sens_no_coloc <- c("egger_intercept","egger_slope_agreement",
                   "weighted_median","weighted_mode",
                   "steiger_direction","leave_one_out",
                   "ivw_Q","ivw_I2")

# ---- 2. Run 1: sex = "both" (panukb catalog) --------------------------------
ard_compare(
  exposure       = "EXPOSURE LABEL",       # e.g. "BMI"
  exposure_units = "EXPOSURE UNITS",       # e.g. "SD BMI"
  groups = list(
    eur_both = list(
      sex           = "both",
      ancestry      = "EUR",
      exposure_snps = snps_both
      # no exposure_sumstats -> acknowledge_no_coloc = TRUE automatically
    )
  ),
  sensitivity_enabled         = sens_no_coloc,
  sensitivity_pass_min        = 5,         # of 8; standard genome-wide MR
  Multiple_testing_correction = "BH",
  clump_opts                  = clump_opts_gw,
  force_refresh               = FALSE,
  verbose                     = TRUE,
  confirm                     = "no"
)

# ---- 3. Run 2: sex-stratified male + female (neale catalog) -----------------
ard_compare(
  exposure       = "EXPOSURE LABEL",
  exposure_units = "EXPOSURE UNITS",
  groups = list(
    eur_female = list(
      sex           = "female",
      ancestry      = "EUR",
      exposure_snps = snps_female
    ),
    eur_male = list(
      sex           = "male",
      ancestry      = "EUR",
      exposure_snps = snps_male
    )
  ),
  sensitivity_enabled         = sens_no_coloc,
  sensitivity_pass_min        = 5,
  Multiple_testing_correction = "BH",
  clump_opts                  = clump_opts_gw,
  force_refresh               = FALSE,
  verbose                     = TRUE,
  confirm                     = "no"
)

# Notes:
# - Without coloc, the run is faster and needs no GWAS-VCF, but you lose the
#   single-causal-variant evidence layer. Steiger directionality + Egger +
#   weighted median/mode carry the QC weight.
# - To enable coloc later: drop in a tabix-indexed GWAS-VCF per group via
#   `exposure_sumstats = "..."`, switch the sensitivity panel to include
#   "coloc", and apply the package patch (see the with-coloc template).
# - Non-EUR ancestries are only allowed with sex = "both".
