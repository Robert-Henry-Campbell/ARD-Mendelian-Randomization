# =============================================================================
# ardmr template: genome-wide phenome-wide MR via OpenGWAS / IEU GWAS
# (run_ieugwasr_ard_compare)
# =============================================================================
#
# Use this when:
#   * Your exposures live in OpenGWAS (you have ieu_ids, not local sumstats).
#   * You want phenome-wide MR with built-in PheWAS instrument screening
#     and coloc-as-sensitivity (coloc is *always on* in this path; the
#     OpenGWAS sumstats are used as the coloc source via ieu_id).
#
# How this differs from the CSV+VCF templates:
#   - Instruments are extracted by `get_instruments()` from OpenGWAS, then
#     LD-clumped + p-filtered upstream. The r2/kb/p_threshold/p_backoff
#     args below govern that upstream pass.
#   - Coloc uses each group's `ieu_id` for the regional sumstats; you do
#     NOT supply a GWAS-VCF here.
#   - No package patch needed (this function bypasses the
#     exposure_sumstats path).
#
# *** AFK / batch warning re: prompt_for_units ***
#   If `prompt_for_units = TRUE` and any ieu_id has missing units in the CSV
#   AND `gwasinfo()` returns no unit, the run BLOCKS forever at readline().
#   This template defaults to FALSE so unattended runs fail fast and clean.
#   FILL IN `exposure_units` FOR EVERY ROW IN YOUR CSV before running AFK.
#
# Required CSV columns (one row per exposure x sex x ancestry combination):
#   ieu_id                      OpenGWAS study id, e.g. "ieu-b-40"
#   exposure                    Display label, e.g. "BMI"
#   exposure_group              Groups rows that should be compared
#                                 side-by-side in one ard_compare() call.
#                                 Within a group, sex is allowed to vary
#                                 (both vs male+female) provided you
#                                 follow the panukb/neale catalog rule.
#   sex                         "both" | "male" | "female"
#   ancestry                    "EUR" (sex-stratified must be EUR)
#   multiple_testing_correction "BH" | "bonferroni"
#   confirm                     "yes" | "no" | "ask"  (Neale-download prompt
#                                 behaviour for that row)
# Optional CSV columns:
#   exposure_units              e.g. "SD BMI". STRONGLY recommended for
#                                 unattended runs (see warning above).
#   phenoscanner_exclusions     Comma-separated regex patterns;
#                                 instruments associated (PheWAS or
#                                 PhenoScanner) with phenotypes matching
#                                 any pattern at p <= phewas_pval are
#                                 dropped. Empty = skip the filter.
#   samplesize                  Per-row override when gwasinfo() returns
#                                 missing/zero N (common for some pQTLs).
#                                 Without this, coloc skips with
#                                 "exposure N unknown".
# -----------------------------------------------------------------------------

setwd("PATH/TO/ardmr")          # path to package source
devtools::load_all()
# devtools::test()              # uncomment for a green-bar check first

Sys.setenv(ARDMR_CACHE_DIR = "PATH/TO/CACHE")     # tens-of-GB scratch dir
Sys.setenv(OPENGWAS_JWT    = "YOUR_JWT_HERE")     # mandatory for this path

# ---- 1. Inputs --------------------------------------------------------------
csv_path <- "PATH/TO/exposures.csv"

# ---- 2. Run -----------------------------------------------------------------
result <- run_ieugwasr_ard_compare(
  csv_path             = csv_path,

  # auth + IO
  # cache_dir          = ardmr_cache_dir(),  # default: ARDMR_CACHE_DIR env var
  # jwt                = Sys.getenv("OPENGWAS_JWT"),  # default

  # interactivity
  prompt_for_units     = FALSE,        # FALSE = fail fast on missing units
                                       # (safe for AFK). Flip to TRUE only
                                       # for hand-held interactive runs.

  # instrument extraction (genome-wide standard)
  p_threshold          = 5e-8,         # genome-wide significance
  p_backoff            = c(5e-8, 5e-7, 5e-6),
                                       # batch-friendly relaxation ladder;
                                       # accepts progressively weaker IVs
                                       # for exposures with no 5e-8 hits.
                                       # Set to NULL for strict single rung.
  r2                   = 0.001,        # strict LD pruning (genome-wide std)
  kb                   = 10000,        # 10 Mb clump window
  f_threshold          = 10,           # per-SNP F-stat min

  # PheWAS / Phenoscanner instrument-exclusion screening
  phewas_pval          = 5e-8,         # acts only if CSV has
                                       # phenoscanner_exclusions patterns
  phenoscanner         = FALSE,        # Phenoscanner v2 API deprecated
                                       # (per README); leave FALSE
  # phenoscanner_pval  = 5e-8,         # only if phenoscanner = TRUE

  # MR sensitivity / multiple testing
  sensitivity_pass_min = 6,            # of 9 (panel always incl. coloc)

  # cache / debug
  force_refresh        = FALSE,        # TRUE only when re-running from scratch
  n_pheno_limit        = NULL          # set to e.g. 5 for a smoke test;
                                       # results not scientifically valid
)

# `result` is a list with: processed exposure_group names, any failures,
# and (if phenoscanner=TRUE) phenoscanner summary tibble + written file paths.
