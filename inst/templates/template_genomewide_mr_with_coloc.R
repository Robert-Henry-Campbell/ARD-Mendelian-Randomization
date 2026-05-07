# =============================================================================
# ardmr template: genome-wide MR  WITH  coloc (per-group GWAS-VCF)
# =============================================================================
#
# Use this when:
#   * You're running phenome-wide MR for a polygenic exposure (e.g. BMI, LDL-C)
#     using a curated set of genome-wide-significant lead SNPs.
#   * You have a tabix-indexed GWAS-VCF (MRC-IEU spec, GRCh37) for the
#     exposure and want coloc as a sensitivity check at each IV's locus.
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

# Full sensitivity panel (ard_compare's default omits coloc -- add it back).
sens_full <- c("egger_intercept","egger_slope_agreement",
               "weighted_median","weighted_mode",
               "steiger_direction","leave_one_out",
               "ivw_Q","ivw_I2","coloc")

# ---- 2. Run 1: sex = "both" (panukb catalog) --------------------------------
ard_compare(
  exposure       = "EXPOSURE LABEL",       # e.g. "BMI"
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
  sensitivity_pass_min        = 6,         # of 9; standard genome-wide MR
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
  clump_opts                  = clump_opts_gw,
  force_refresh               = FALSE,
  verbose                     = TRUE,
  confirm                     = "no"
)

# Notes:
# - Per-group `exposure_sumstats` (VCF) and `ieu_id` (OpenGWAS) are mutually
#   exclusive (R/coloc_source_resolver.R:48-50). Pick one per group.
# - Non-EUR ancestries are only allowed with sex = "both"
#   (R/ard_compare.R:261-268).
# - Coloc with hundreds of IVs scales -- per-IV ±500kb windows are merged,
#   but expect a few minutes per outcome. For very large IV sets consider
#   restricting to nominally-significant outcomes (already gated by default,
#   see README "Per-outcome gating").
