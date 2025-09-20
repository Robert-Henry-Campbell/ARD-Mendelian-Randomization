# ARDMR

**One-exposure → many-phenotype Mendelian Randomization.** Automates outcome resolution (ICD-10 ↔ GBD), SNP row extraction from Pan-UKB / Neale, harmonisation, a battery of **8 sensitivity/QC checks**, multiple testing correction, and scalable visualisation (volcano, Manhattan). Phenotype tables span age-related diseases (ARDs) and other traits, with the `ARD_selected` column flagging the ARD subset. Designed to **ship fast**; runtime efficiency is secondary.

---

## Install

Before developing or running tests, install package dependencies:

```bash
Rscript -e "source('R/install_deps.R'); install_deps()"
```

```r
# requires: R >= 4.1, tabix in PATH, bgzip (.tbi) support
# and CRAN/BioC deps: data.table, dplyr, ggplot2, TwoSampleMR, etc.

# devtools::install_github("YOUR_ORG/ardmr")
```

**System tools required**

* `tabix` and `bgzip` (htslib) in your PATH
* `wget` or `curl` for downloads (Neale / remote Pan-UKB slices)

## Cache directory

ARDMR caches large reference files (e.g., variant manifests) under a user
specified directory. Set the location via the `ARDMR_CACHE_DIR` environment
variable before running any functions:

```bash
export ARDMR_CACHE_DIR=/path/to/cache
```

Or within R:

```r
ardmr::set_cache_dir("/path/to/cache")
```

---

## Quick start

```r
library(ardmr)

# 1) Exposure SNPs (TwoSampleMR-style; RSIDs required)
exposure <- tibble::tibble(
  SNP = c("rs123","rs456","rs789"),
  beta = c(0.05, -0.03, 0.02),
  se   = c(0.01,  0.01, 0.01),
  effect_allele = c("A","C","T"),
  other_allele  = c("G","T","C")
)

# 2) Run phenome-wide MR
res <- run_phenome_mr(
  exposure = "LDL-C",
  exposure_snps = exposure,
  exposure_units = "1-SD LDL-C",
  sex = "both",                # "both" => Pan-UKB, "male"/"female" => Neale
  ancestry = "EUR",
  sensitivity_pass_min = 6,
  multiple_testing_correction = "BH",
  neale_gwas_dir = NULL        # required if sex != "both"
)

res$results_df      # tidy results across outcomes
res$volcano         # ggplot
res$manhattan       # ggplot
```

All run-specific plots, diagnostics, and logs are written to
`<cache_dir>/output/<exposure>/<sex>/<ancestry>`.

`run_phenome_mr()` builds `MR_df` from phenotypes defined by the Global Burden of Disease (GBD) cause ontology. Age-related diseases and other outcomes share this structure and are distinguished by an `ARD_selected` flag, so the instructions above apply uniformly across phenotypes.`

---

## Package workflow (modules)

Below, each script is listed with **Inputs → Actions → Outputs**, using the project’s wording wherever possible.
We also mark **Status** for the current v0:

* ✅ **Complete (v0)**
* 🟡 **Partial (MVP in place; refine/extend)**
* ⏳ **Planned / later**

---

### 0) `load_data`  — *manifest loader*

**Purpose:** Loads the packaged / cached manifest objects.

* **Inputs**

  * (implicit) packaged files or cached copies
* **Actions**

  * Load the **Pan-UKB pheno manifest df**
  * Load the **Neale pheno manifest df**, including:

    * Neale **male** phenotype supplement
    * Neale **female** phenotype supplement
* **Outputs**

  * In-memory tibbles used by `Outcome_setup()`
* **Status:** 🟡 Partial — manifests are loaded via helper functions (e.g., `load_panukb_manifest()`, `load_neale_manifest()`). Wire to your packaged files in `inst/extdata/` and cache to RDS.

---

### 1) `Run_phenome_mr.R` — *orchestrator; runs the whole show*

* **Inputs**

  * `exposure_snps`, `sex`, `ancestry`
  * sensitivity settings, plotting flags, `neale_gwas_dir`, `cache_dir`
* **Actions**

  * Calls `Outcome_setup()` to build **MR\_df** (phenotypes + provider metadata)
  * Maps exposure RSIDs to provider-specific positions via `Exposure_SNP_mapper()`
  * Fetches outcome SNP rows:

    * `Panukb_SNP_grabber()` if `sex == "both"`
    * `Neale_GWAS_checker()` → `Neale_tbi_maker()` → `Neale_SNP_grabber()` if sex-specific
  * Runs **MR + sensitivity/QC** via `MR_business_logic()`
  * Builds **Volcano** and **Manhattan** plots
* **Outputs**

  * `list(MR_df, results_df, volcano, manhattan)`
* **Status:** ✅ Complete (v0)

---

### 2) `Outcome_setup.R`

**i. Defines `MR_df`**

* **Inputs (Params)**

  * `sex` (both/male/female)
  * `ancestry` (EUR, AFR, …)
* **Actions**

  1. **Sets `MR_df` equal to relevant sex using the GBD-derived phenotype tibble.**
     *In practice:* filter the tibble to the requested sex and attach a **provider** (Pan-UKB for `both`; Neale for `male`/`female`).  
     All phenotypes originate from the Global Burden of Disease cause ontology and include an `ARD_selected` flag distinguishing age-related diseases from other outcomes.  
     *(Your earlier phrase “define `MR_df` as one of `bothsex_ARD`, `male_ARD`, `female_ARD`” is captured here via filtering / provider selection.)*
  2. **Log phenotypes** at **ICD-10** and **GBD cause** level.
  3. **Map `ICD10_explo`** in `MR_df` to the **provider pheno manifest** (Pan-UKB if sex=both, Neale otherwise).

     * Log **mapped vs. unmapped** rows
     * **Drop unmapped** (no available GWAS)
     * `MR_df` now contains phenotype data + **GWAS manifest** data per row, with ARD and non-ARD phenotypes distinguished by `ARD_selected`
* **Checks**

  1. If `sex` is *female* or *male*, require `ancestry == "EUR"`.
  2. Validate `sex` and `ancestry` are approved options.
  3. (Hook for additional checks.)
* **Outputs**

  * `MR_df`
* **Status:** ✅ Complete (v0) — provider choice + join + logs + basic checks

---

### 3) `Exposure_SNP_mapper.R`

* **Inputs**

  * `exposure_snps` (must have: `SNP`, `beta`, `se`, `effect_allele`, `other_allele`)
  * `sex`
* **Outputs**

  * `exposure_snps` **augmented** with provider-specific positions
* **Actions**

  1. **Check columns** on `exposure_snps`.
  2. If `sex = "both"` (Pan-UKB):

     * Ensure/Download **Pan-UKB variants manifest** *(big)*
     * Map RSID → **chr\:pos** using manifest
     * Append position as `panukb_chr`, `panukb_pos`
  3. If `sex = "male" | "female"` (Neale):

     * Ensure/Download **Neale variants manifest** *(big)*
     * Map RSID → **chr\:pos**
     * Append as `neale_chr`, `neale_pos`
* **Status:** ✅ Complete (v0) — stubs load manifests, mapping + columns are enforced

---

### 4) `Panukb_SNP_grabber.R`

* **Inputs**

  * `Exposure_snps`, `MR_df`, `Ancestry`
* **Outputs**

  * `MR_df` with list-column `outcome_snps` (TwoSampleMR outcome schema)
* **Actions**

  1. Validate `exposure_snps` columns (`panukb_chr`, `panukb_pos` present).
  2. **For each phenotype** (row) in `MR_df`, use manifest metadata (download/`tabix` URLs):

     1. For **each chrom\:pos** in `exposure_snps`, **tabix-slice** that row from remote sumstats.
        *Notify counts/size/progress (planned enhancement).*
     2. **Drop other ancestries** (keep requested `ancestry`).
     3. **Append RSIDs** by joining back on chr\:pos.
     4. **Transform to TwoSampleMR outcome format**, preserving true RSID and chr\:pos.
     5. **Store** the prepared tibble in `MR_df$outcome_snps`.
* **Status:** 🟡 Partial — MVP implemented (remote tabix slice + formatting); progress/meters + robust ancestry handling can be extended.

---

### 5) `Neale_GWAS_checker.R`

* **Inputs**

  * `MR_df`, `Neale_GWAS_dir`
* **Outputs**

  * None (side-effects: ensures files present)
* **Actions**

  1. If `neale_gwas_dir` exists, **check** whether required **.tsv.bgz** exist for all phenotypes:

     * If missing → **download** via **manifest commands**.
     * If present → skip.
  2. Ensure **.tbi** indices exist (or leave to `Neale_tbi_maker`).
  3. If dir missing → **create**, **estimate total download** (assume \~500 MB per phenotype), prompt user, and download with progress.
     *Note: non-sex-specific MR can be done with minimal download (Pan-UKB remote slicing).*
  4. **Print & log** actions taken.
* **Status:** 🟡 Partial — checks & downloads wired; add size estimation, prompts, and structured logs as next polish.

---

### 6) `Neale_tbi_maker.R`

* **Inputs**

  * `Neale_GWAS_dir`
* **Outputs**

  * None (side-effects: `.tbi` created)
* **Actions**

  1. **Find all `.bgz`** in dir.
  2. **Generate `.tbi`** for each (`tabix -s <seq> -b <start> -e <end>`).
  3. **Log** what was indexed.
* **Status:** ✅ Complete (v0) — column settings may need provider-specific tweaks.

---

### 7) `Neale_SNP_grabber.R`

* **Inputs**

  * `Exposure_snps`, `MR_df`, `Neale_GWAS_dir`
* **Outputs**

  * `MR_df` with list-column `outcome_snps` (TwoSampleMR outcome schema)
* **Actions**

  1. For every row in `MR_df`:

     1. **Locate** the `.bgz` in `Neale_GWAS_dir` (via manifest mapping; typically by `local_filename`, or by `ICD10_explo` prefix if needed).
     2. For every SNP in `exposure_snps`, **tabix-slice** the matching **chr\:pos** row.
     3. **Assemble** rows into `neale_temp_df`.
     4. **Add RSID** column by joining back on chr\:pos.
     5. **Transform** to TwoSampleMR outcome format.
     6. **Store** the tibble in `MR_df$outcome_snps`.
* **Status:** ✅ Complete (v0)

---

### 8) `MR_business_logic.R`

* **Inputs**

  * `MR_df`, `Exposure_snps`
  * `Sensitivity_pass_min`, `Sensitivity analyses`
  * Plot flags: `scatterplot`, `snpforestplot`, `leaveoneoutplot`
* **Outputs**

  * Updated `MR_df`
  * Tidy `results_df` (per-outcome summary)
* **Actions**

  * **For each phenotype** in `MR_df`:

    1. **Harmonize** exposure/outcome; **log** `results_nsnp_before` → `results_nsnp_after`.
       Use MAF/EAF where available to deal with palindromes.
    2. **Main MR** (IVW, MR-Egger, Weighted median, Weighted mode).
       Optional per-outcome plots:

       * `mr_scatter_plot()`
       * `mr_forest_plot()` (per-SNP IVW)
       * `mr_leaveoneout_plot()`
    3. **Heterogeneity:** collect **Cochran’s Q** and **I²**.
       *Flag if Q p<.05 **OR** I²>50% (not sufficient alone to drop).*
    4. **MR-Egger:** slope±SE±p and **intercept p**.

       * Intercept p<.05 → **possible directional pleiotropy**
       * **Slope sign** disagree with IVW → **flag**
       * If **I²\_GX < .6**, **downgrade** Egger slope disagreements (low precision)
    5. **Weighted median/mode:** record estimates; **flag** if sign discordance or |Δβ| > 1×SE\_IVW.
    6. **LOO IVW:** refit, **flag** if any single-SNP **flips** the IVW direction.
    7. **Steiger:** compute **fraction** with r²(outcome) > r²(exposure); **flag** if >30%.
  * **Sensitivity rule (8 checks; NA-aware):**

    * **Core viability:** IVW returns **finite** β and SE
    * Require **≥ `Sensitivity_pass_min`** of the enabled checks to **pass**
    * Analyses that **fail to run** are excluded from the denominator but **logged**
  * **Emit** tidy `results_df` with key columns (β\_IVW, SE\_IVW, p\_IVW, Q/I², Egger intercept p, I²\_GX, LOO flip flag, Steiger fraction, nsnp before/after, counts of checks applied/passed, QC gate)
* **Status:** ✅ Complete (v0) — implements your full pass/fail logic with NA-aware denominator

---

### 9) `Manhattan_plot.R`

* **Inputs**

  * `Results_df` (or `MR_df` post-run)
  * `Multiple_testing_correction` (BH/Bonferroni)
* **Outputs**

  * Manhattan-style plot across outcomes
* **Actions**

  * Plot **−log10(p\_IVW)** for **QC-pass** outcomes
  * Draw **Bonferroni** threshold line if requested; BH visualised via discovery colouring
* **Status:** ✅ Complete (v0)

---

### 10) `Enrichment_plot` *(do later)*

* **Purpose**

  * Hypergeometric over-representation of **GBD causes** among discoveries (`q<0.05`)
* **Status:** ⏳ Planned — code hooks exist; you already have similar logic elsewhere.

---

### 11) `Volcano_plot.R`

* **Inputs**

  * `Results_df` (or `MR_df` post-run)
  * `Multiple_testing_correction` (BH/Bonferroni)
* **Output**

  * Volcano plot
* **Attributes**

  * **x** = effect direction/size (uses **Z = β\_IVW / SE\_IVW** for comparability)
  * **y** = −log10(p\_IVW) (or −log10(q) visually)
  * **Color** = GBD cause (level 2/3)
  * **Alpha** = QC pass (opaque) vs fail (faint)
  * **Size** = `results_nsnp_after`
  * **Shape/ring** = residual issues (e.g., heterogeneity fail or Egger intercept p<0.05) *(optional)*
  * **Labels** = top-N by |Z| or smallest q (`ggrepel`)
  * **Thresholds**

    * Bonferroni line at −log10(α\_Bonf)
    * BH/FDR: colour points by `q<0.05` (no single line)
  * **Variants**

    * Mirror volcano by GBD cause (small multiples)
    * QC overlay / influence-aware LOO ring
    * Volcano + enrichment strip
* **Status:** ✅ Complete (v0) — core volcano in place; add optional overlays as desired

---

## Multiple testing

* Default on **QC-pass** set: `BH` or `bonferroni`
* You can expose `mtc_population = c("results_qc_pass","all_tested")` later if reviewers want both.

---

## Data expectations

* **Genome build:** GRCh37. RSID→pos crosswalks tied to provider **variant manifests**.
* **Exposure:** RSIDs required; alleles, β, SE. EAF / N optional but improves harmonisation (palindromes) and Steiger.
* **Outcomes:** Provider sumstats accessible via:

  * **Pan-UKB:** remote tabix URL + manifested column map
  * **Neale:** local `.tsv.bgz` + `.tbi` in `neale_gwas_dir`

### Packaged manifests

The repository bundles several `.rda` files (in `data/`) consumed by
`Outcome_setup()`:

 - `bothsex_ARD.rda`, `female_ARD.rda`, `male_ARD.rda` – GBD-derived phenotype tables
   for both sexes, females, and males; each row includes `ARD_selected` to identify
   age-related diseases.
- `neale_female_manifest.rda`, `neale_male_manifest.rda` – Neale sex-specific
  phenotype manifests with case and control counts.
- `neale_file_manifest.rda` – Neale download information for each GWAS.
- `panukb_pheno_manifest.rda` – Pan-UKB phenotype manifest with metadata and
  download links.

---

## QC checklist (8)

1. **IVW core viability** (finite β & SE)
2. **Heterogeneity** (Q p≥.05 **and** I² ≤ 50%)
3. **Egger intercept** p≥.05
4. **Egger slope sign** agrees with IVW (unless **I²\_GX < .6**)
5. **I²\_GX ≥ .6** (precision for Egger)
6. **Weighted median** \~ sign & |Δβ| ≤ 1×SE\_IVW
7. **Weighted mode** \~ sign & |Δβ| ≤ 1×SE\_IVW
8. **LOO flip** absent & **Steiger** fraction bad ≤ 30%  *(two checks recorded; denominator shrinks if any step not runnable)*

Pass gate: **≥ `sensitivity_pass_min`** of enabled checks.

---

## Logging & outputs

* `results_df` CSV and summary tables saved under
  `<cache_dir>/output/<exposure>/<sex>/<ancestry>`
* Optional per-outcome harmonised RDS (hook available)
* Console and file logs stored alongside run outputs

---

## Roadmap

* **Provider robustness:** richer column maps; auto-detect sequence/position columns per file
* **Progress bars & size estimates** for downloads and remote tabix pulls
* **Enrichment\_plot:** finalise & document universe choice (`"all_tested"` vs `"all_ards"`)
* **Vignette:** end-to-end demo with tiny fixtures
* **Caching:** persistent crosswalks and harmonised objects for resumable runs

---

## Citation / dependencies

Built on **TwoSampleMR** and standard tidyverse/htslib tools. Please cite those projects when publishing results produced with ARDMR.

---

### Contact

Maintainer: **Robert Campbell** · Oxford University
Issues/PRs welcome.
