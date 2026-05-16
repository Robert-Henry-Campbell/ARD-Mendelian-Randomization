# =============================================================================
# ardmr template: cis-gene / drug-target MR  WITH  coloc (per-group GWAS-VCF)
# =============================================================================
#
# Use this when:
#   * You're running cis-MR for a single drug target / gene region.
#   * You have a tabix-indexed GWAS-VCF (MRC-IEU spec, GRCh37) for the
#     exposure and want coloc as a sensitivity check.
#
# Input pattern: sex-stratified discovery SNPs + sex-stratified GWAS-VCFs.
# If you only have sex-combined inputs, you can reuse them across groups
# but it's a study-design oddity (note in your methods section).
#
# Catalog constraint (R/ard_compare.R:247-251):
#   sex = "both" forces the panukb catalog; sex = "male"/"female" forces the
#   neale catalog. They cannot be mixed in one ard_compare() call. So a
#   "both, female, male" comparison needs TWO ard_compare() calls.
#
# *** Required local patch to ard_compare.R (coloc-enabled runs only) ***
#   In R/ard_compare.R, inside the run_phenome_mr() call (~line 326):
#     - Add `exposure_sumstats = info$exposure_sumstats,`
#     - Update `acknowledge_no_coloc` to:
#         is.null(info$exposure_id) && is.null(info$ieu_id) && is.null(info$exposure_sumstats)
#   Without this patch the per-group VCF is silently dropped and coloc skipped.
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

# Tabix-indexed GWAS-VCF paths (a .tbi must sit next to each .vcf.gz).
vcf_both   <- "PATH/TO/exposure_sumstats_both.vcf.gz"
vcf_female <- "PATH/TO/exposure_sumstats_female.vcf.gz"
vcf_male   <- "PATH/TO/exposure_sumstats_male.vcf.gz"

# Cis-MR (drug-target) clump opts: more permissive than genome-wide.
# Schmidt AF 2020 Nat Comms / Burgess cis-MR guidance: r2 ~0.1-0.3, kb 250-1000.
clump_opts_cis <- list(
  already_clumped    = FALSE,
  already_p_filtered = FALSE,
  p_backoff          = c(5e-8, 1e-6),   # strict first; fall back if too few IVs
  r2                 = 0.1,             # 0.1-0.3 typical for cis-MR
  kb                 = 1000,            # cis-window (~±500 kb of gene)
  f_threshold        = 10,
  maf_min            = 0.01
)

# Full sensitivity panel (ard_compare's default omits coloc -- add it back).
sens_full <- c("egger_intercept","egger_slope_agreement",
               "weighted_median","weighted_mode",
               "steiger_direction","leave_one_out",
               "ivw_Q","ivw_I2","coloc")

# ---- 2. Run 1: sex = "both" (panukb catalog) --------------------------------
ard_compare(
  exposure       = "EXPOSURE LABEL",       # e.g. "GLP1R (BMI Path)"
  exposure_units = "EXPOSURE UNITS",       # e.g. "SD BMI"
  groups = list(
    eur_both = list(
      sex               = "both",
      ancestry          = "EUR",
      exposure_snps     = snps_both,
      exposure_sumstats = vcf_both
    )
  ),
  sensitivity_enabled         = sens_full,
  sensitivity_pass_min        = 6,         # of 9; cis-MR with few IVs -> 4-5
  Multiple_testing_correction = "BH",
  clump_opts                  = clump_opts_cis,
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
      sex               = "female",
      ancestry          = "EUR",
      exposure_snps     = snps_female,
      exposure_sumstats = vcf_female
    ),
    eur_male = list(
      sex               = "male",
      ancestry          = "EUR",
      exposure_snps     = snps_male,
      exposure_sumstats = vcf_male
    )
  ),
  sensitivity_enabled         = sens_full,
  sensitivity_pass_min        = 6,
  Multiple_testing_correction = "BH",
  clump_opts                  = clump_opts_cis,
  force_refresh               = FALSE,
  verbose                     = TRUE,
  confirm                     = "no"
)

# Notes:
# - Per-group `exposure_sumstats` (VCF) and `ieu_id` (OpenGWAS) are mutually
#   exclusive (R/coloc_source_resolver.R:48-50). Pick one per group.
# - Non-EUR ancestries are only allowed with sex = "both"
#   (R/ard_compare.R:261-268).
# - For drug-target MR, coloc PP.H4 > 0.8 is the most informative sensitivity
#   anchor when IV counts are small; MR-Egger / median / mode get noisy.
