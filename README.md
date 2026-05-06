# ARDMR

Phenome-wide Mendelian randomization for age-related disease research.

## Overview

`ardmr` orchestrates the end-to-end workflow required to run one-exposure →
phenome Mendelian randomization (MR) analyses across age-related diseases (ARDs)
and related traits. The package combines curated phenotype manifests, automated
variant lookups, Mendelian randomization with multiple sensitivity checks, and
publication-ready visualisations. Outcomes are grouped using the Global Burden
of Disease (GBD) cause hierarchy and flagged when they belong to the core ARD
set so that the same codepath can support ARD-only or pan-phenome analyses.

The primary entry point is [`run_phenome_mr()`](./R/run_phenome_mr.R), which:

1. Loads the packaged phenotype catalog for the requested sex and ancestry.
2. Downloads or reuses variant manifests and maps exposure SNPs to provider
   coordinates.
3. Computes a pairwise LD matrix between the supplied exposure instruments
   against the 1000 Genomes Phase 3 reference panel and writes
   `exposure_snps_ld_matrix.csv` to the run output directory.
4. Retrieves outcome SNP rows either from Pan-UK Biobank (sex = "both") or local
   Neale Lab summary statistics (sex = "male"/"female").
5. Runs inverse-variance weighted MR along with up to nine supporting
   sensitivity and quality-control checks (Cochran's Q, I², MR-Egger,
   weighted median/mode, Steiger directionality, leave-one-out, and
   Bayesian colocalisation via `coloc.abf` / `coloc.susie`).
6. Produces tidy result tables and ggplot2-based volcano and Manhattan plots and
   writes run artefacts to disk under a per-run hash subfolder so different
   parameterisations don't overwrite each other.

## Key capabilities

- **Automated phenotype manifests.** Packaged datasets (`bothsex_all`,
  `male_all`, `female_all`) link ICD-10 codes to GBD causes, Pan-UKB phenocodes,
  and Neale Lab identifiers while tracking which outcomes are considered ARDs.
- **Provider-aware SNP mapping.** `exposure_snp_mapper()` joins user-supplied
  instruments to Pan-UKB or Neale variant manifests and handles multi-allelic
  loci.
- **Automatic exposure-SNP LD matrix.** `compute_ld_matrix()` runs inside
  `run_phenome_mr()` (step 3 above) and writes a per-run
  `exposure_snps_ld_matrix.csv` of pairwise correlations against 1000G v3.
  Useful for diagnosing IV correlation structure post hoc.
- **Remote-friendly Pan-UKB queries.** `panukb_snp_grabber()` streams tabix
  slices over HTTPS via Bioconductor's `Rsamtools`, so a system `tabix`
  executable is not required.
- **Neale GWAS management.** `neale_gwas_checker()` verifies large Neale Lab
  files locally, optionally downloading missing datasets with resumable logging.
- **QC-rich MR engine.** `mr_business_logic()` harmonises, runs multiple MR
  estimators, and records up to nine sensitivity checks (including Bayesian
  colocalisation) with configurable pass thresholds.
- **Bayesian colocalisation.** `coloc_business_logic()` runs `coloc::coloc.abf`
  on a ±500 kb window around each IV and `coloc::coloc.susie` whenever a
  per-locus LD matrix is computable. Per-locus artefacts (CSV summaries,
  stacked locus-zoom plots, sensitivity curves) are written under
  `<run_dir>/<outcome_slug>/coloc/<chrX_start_end>/`. See [Colocalisation](#colocalisation).
- **Plotting + enrichment.** Helpers such as `plots.R`,
  `run_phenome_mr_plotting_only()`, and `enrichment_business_logic.R` generate
  reusable visuals and enrichment summaries for downstream reporting.

## Installation

### Prerequisites

- R ≥ 4.1
- Write access to a cache directory with tens of gigabytes of free space. Neale
  Lab summary statistics can exceed 400 GB if all phenotypes are downloaded.
- Optional (recommended) command line tools: `wget` or `curl` for faster Neale
  downloads. `tabix` is not required because remote queries use `Rsamtools`.

### Install package dependencies

The repository ships with an installer helper. Run it once after cloning to
install CRAN and Bioconductor dependencies listed in `DESCRIPTION`, including
`MendelianRandomization`, `jsonlite`, and other core tooling:

```bash
Rscript -e "source('R/install_deps.R'); install_deps()"
```

### Install the package locally

To install the development version in your R library:

```r
# install.packages("pak") # if needed
pak::pkg_install("aging-research-discovery/ARD-Mendelian-Randomization")
# or: devtools::install_github("aging-research-discovery/ARD-Mendelian-Randomization")
```

During installation Bioconductor packages (`Rsamtools`, `GenomicRanges`,
`IRanges`) will also be installed via the `Remotes` field.

## Configure caching

Large reference files (variant manifests, downloaded summary statistics, run
outputs) are stored under a user-controlled cache root. Set the location before
calling any package functions. Two common approaches are:

```bash
export ARDMR_CACHE_DIR=/absolute/path/to/cache
```

or inside R:

```r
ardmr::set_cache_dir("/absolute/path/to/cache")
```

`ardmr_cache_dir()` will throw an error if the environment variable is missing,
ensuring that expensive downloads never default to a temporary directory.

## Preparing exposure instruments

`run_phenome_mr()` expects harmonised exposure instruments in a
TwoSampleMR-compatible format. At minimum the following columns are required:

| column name             | description                                           |
|------------------------|-------------------------------------------------------|
| `id.exposure`          | identifier for the exposure GWAS                      |
| `SNP`                  | RSID-formatted variant identifier                     |
| `beta.exposure`        | effect size for the exposure                          |
| `se.exposure`          | standard error for `beta.exposure`                    |
| `effect_allele.exposure` | effect allele (A/C/G/T)                            |
| `other_allele.exposure`  | non-effect allele                                   |
| `pval.exposure`        | p-value associated with the exposure effect estimate  |

Additional metadata columns (sample size, eaf, etc.) are retained and
propagated when available. Use `exposure_snp_mapper()` prior to MR runs if you
need access to provider-specific chromosome/position columns.

## Quick start

```r
library(ardmr)

# Tell ardmr where to store manifests, downloads, and outputs
set_cache_dir("/data/ardmr-cache")

# Exposure instruments (example with LDL-C)
exposure <- tibble::tibble(
  id.exposure = "ldl.gwas",
  Exposure = "LDL cholesterol",
  SNP = c("rs12345", "rs67890", "rs24680"),
  beta.exposure = c(0.05, -0.03, 0.02),
  se.exposure   = c(0.01, 0.01, 0.01),
  effect_allele.exposure = c("A", "C", "T"),
  other_allele.exposure  = c("G", "T", "C"),
  pval.exposure = c(1e-8, 5e-9, 2e-7)
)

results <- run_phenome_mr(
  exposure = "LDL cholesterol",
  exposure_snps = exposure,
  exposure_id = "ieu-b-110",  # OpenGWAS id; required for colocalisation
  exposure_units = "1-SD increase",
  ancestry = "EUR",
  sex = "both",                  # "both" uses Pan-UKB, "male"/"female" use Neale
  sensitivity_pass_min = 6,      # require six of the nine checks
  Multiple_testing_correction = "BH"
)

# To run without coloc (no exposure_id / exposure_sumstats), pass
# `acknowledge_no_coloc = TRUE`. Coloc will be removed from the sensitivity
# panel for the run and a single INFO line will note this in the run log.

# Inspect outputs
results$results_df   # tidy phenome-wide MR summary; includes results_coloc_* cols
results$volcano      # ggplot object; can be customised further
results$manhattan    # ggplot object for -log10(p) across outcomes
```

Each run writes a timestamped log file, diagnostic plots, and CSV exports under
a per-run-hash subfolder:

```
<cache_dir>/output/<exposure>/<sex>/<ancestry>/<run_hash>/
```

The run hash is derived from the inputs (IV set, sensitivity panel, clump opts,
coloc params, multiple-testing rule), so different parameter choices land in
distinct folders and cannot overwrite each other. It has the form
`<iv_hash>__<params_hash>`; the prefix matches the suffix on cached

outcome-SNP files (`<phenotype>__<iv_hash>.rds`) so you can trace lineage at a
glance. Variant lookups continue to be cached under provider-specific 
subdirectories so subsequent runs are faster.

### Quick test runs (`n_pheno_limit`)

For end-to-end smoke tests and CI integration, pass `n_pheno_limit = N` to
truncate the outcome list to the first `N` phenotypes immediately after
`Outcome_setup()`:

```r
run_phenome_mr(..., n_pheno_limit = 10)   # only the first 10 outcomes
```

Multiple-testing thresholds, p-values, and enrichment p-values from a
truncated run are **not scientifically valid** — use this knob for
debugging and pipeline regression-tests only. A `WARN` line is emitted in
the run log so a truncated run is impossible to mistake for a full one.
The `dev/integration_test_ieugwasr.R` script demonstrates the recommended
usage.

## Running multiple analyses concurrently

`run_phenome_mr()` does **not** lock its cache directory and writes most cache
and output files directly to their destination paths (no atomic rename). Two
runs that share a cache can therefore corrupt each other's intermediates and
outputs unless one of the safe configurations below applies.

### Safe

- **Each concurrent run uses a distinct `ARDMR_CACHE_DIR`.** Fully isolated;
  no shared state.
- **Shared cache, fully warm.** If the variant manifests, 1000G LD reference
  panel, tabix indices, and per-IV outcome-SNP RDS files already exist on
  disk, the cache writers all early-return on `file.exists()` and the runs
  only read.
- **Shared cache, warm, *and* each run has a different input hash** (different
  IV set, sex, ancestry, sensitivity panel, clumping options, coloc params,
  or multiple-testing rule). Output folders segregate by `<run_hash>` so
  `results.rds`, `run_manifest.json`, plots, and per-run logs do not collide.
- **The `testthat` suite** (`devtools::test()`). Every test scopes
  `ARDMR_CACHE_DIR` to its own `tempfile()` cache via `withr::local_envvar`,
  so parallel test execution is safe.

### Not safe

- **Cold cache + concurrency.** Two simultaneous runs that both need to
  download the variant manifest, the 1000G LD reference panel, or a tabix
  index will race on the same destination file. There is no atomic rename
  and no lock, so a partial download can leave a silently corrupt file that
  every later run treats as a cache hit.
- **Two runs with identical inputs sharing a cache.** Both produce the same
  `<run_hash>` and write to the same `results.rds`, `run_manifest.json`,
  per-IV outcome RDS, LD-matrix RDS, and per-run log paths. Last writer wins;
  an interrupted writer leaves a corrupt file.
- **Parallelism inside a single R session** (e.g. `future`, `furrr`,
  `parallel::mclapply`). The Pan-UKB grabber and the coloc data fetcher both
  set `HTS_CACHE_DIR` via `Sys.setenv()`, which is process-global and races
  between workers. The structured logger is also process-global and gets
  reset by every call to `setup_logging()`. Use separate `Rscript`
  invocations instead.

### Recommended workflow for parallel runs

1. Pre-warm the cache once, serially, by running a single small job (or
   issuing the manifest / LD-reference / tabix downloads directly).
2. Launch each parallel job as a separate `Rscript` process.
3. Either give each job its own `ARDMR_CACHE_DIR`, or guarantee that each
   job's input hash differs (different exposure, IV set, sex/ancestry, or
   parameter sweep value).

## Quality-control checks

Up to nine sensitivity diagnostics are calculated for every outcome when data
permit:

1. Cochran's Q (heterogeneity) p-value ≥ 0.05
2. Cochran's Q I² ≤ 50%
3. MR-Egger intercept p-value ≥ 0.05
4. MR-Egger slope sign agrees with IVW (unless I²<sub>GX</sub> < 0.6)
5. I²<sub>GX</sub> ≥ 0.6
6. Weighted median estimate within one IVW standard error
7. Weighted mode estimate within one IVW standard error
8. Leave-one-out analysis shows no sign-flip and Steiger filtering failures ≤ 30%
9. Bayesian colocalisation PP.H4 > 0.80 in at least one IV-anchored locus
   (`coloc.abf`, with `coloc.susie` as a multi-causal-variant supplement when
   LD is computable). Only computed for outcomes that pass the gating rules
   described in [Colocalisation](#colocalisation).

`sensitivity_pass_min` defines how many of the enabled checks an outcome must
pass to be considered QC-approved. This flag feeds into the default multiple
testing correction and plot aesthetics. Coloc is included in the default
panel; pass `sensitivity_enabled` without `"coloc"` to disable it for a run.

## Colocalisation

Bayesian colocalisation is run as one of the sensitivity checks via
[`coloc_business_logic()`](./R/coloc_business_logic.R). It builds a ±500 kb
window around each IV (overlapping windows are merged), pulls regional
summary statistics for both the exposure and the outcome, and applies
`coloc::coloc.abf`. Where in-sample-style LD is computable against the
1000G v3 panel, `coloc::coloc.susie` is also fitted to relax the
single-causal-variant assumption. The MHC region (chr6:25–35 Mb) is
skipped by default; toggle with `coloc_skip_mhc = FALSE`.

### Required exposure data

Coloc needs a **regional summary-statistics source** for the exposure, in
addition to the instrument table. Either:

- `exposure_id = "ieu-b-110"` — an OpenGWAS study id; regions are pulled via
  `ieugwasr::associations()` (requires `OPENGWAS_JWT`), or
- `exposure_sumstats = "/path/to.vcf.gz"` — a tabix-indexed
  [GWAS-VCF](https://github.com/MRCIEU/gwasvcf) (MRC-IEU spec) on disk
  (GRCh37 only; FORMAT fields `ES` and `SE` required).

If you supply only `exposure_snps` and neither `exposure_id` nor
`exposure_sumstats`, `run_phenome_mr()` errors and asks you to set
`acknowledge_no_coloc = TRUE`; doing so removes coloc from the sensitivity
panel for that run.

### Per-outcome gating

To keep coloc time bounded on phenome-wide runs, coloc is **only run on
outcomes that are nominally significant in IVW (raw p < 0.05)**. The first
analysed outcome with at least one IV after harmonisation is exempted —
it always runs as a mandatory smoke test, so a clean run log will always
contain at least one `coloc: ENTER coloc_business_logic ...` line.
Outcomes skipped by the gate get a clear log entry:

```
coloc: skipping '<NAME>' (IVW p=<value> >= 0.050, not first analysed outcome)
```

and **no `coloc/` subfolder** under their outcome directory. This is by
design — empty folders are not "missing artefacts".

### Per-locus artefacts

For every outcome that passes the gate, a `coloc/<chrX_start_end>/`
subfolder is written for each merged locus. The locus containing the
IV with the lowest exposure p-value (computed after MHC removal) is
prefixed `1-LEAD-` (e.g. `coloc/1-LEAD-chr5_1234567_1244567/`) so it
sorts to the top and is easy to spot when reviewing a run. Each folder
contains:

- `coloc_abf_summary.csv` — PP.H0 through PP.H4 plus the top SNP
- `coloc_susie_summary.csv` — SuSiE credible-set table (when `coloc.susie` ran)
- `stacked_manhattan.{png,pdf}` — exposure (top) and outcome (bottom) regional
  −log₁₀(p) plots with IV positions highlighted
- `coloc_dataset_exposure.{png,pdf}` and `coloc_dataset_outcome.{png,pdf}` —
  `coloc::plot_dataset()` Manhattan for each trait
- `coloc_sensitivity.{png,pdf}` — `coloc::sensitivity()` curve for the p12 prior
- `coloc_aligned_data.csv` — the harmonised regional table fed to coloc

Per-outcome columns are added to the returned `MR_df` and `results_df`:

| column | meaning |
|---|---|
| `results_coloc_n_loci_total`  | merged loci tested for that outcome |
| `results_coloc_n_loci_pass`   | loci with PP.H4 > 0.80 |
| `results_coloc_max_PPH4`      | max PP.H4 across all loci (abf and susie) |
| `results_coloc_method`        | "susie" / "abf" / "none" — which method's PP.H4 cleared 0.80 |

### Run-level diagnostics

After every run, a `COLOC_DIAGNOSTICS.txt` file is written at the run root
summarising how many outcomes were eligible, attempted, returned loci,
and produced artefact folders, plus a per-skip-reason breakdown:

```
COLOC TELEMETRY (mr_business_logic)
  outcomes_total              : 396
  coloc_eligible_by_gate      : 23   (first-outcome smoke test or IVW p < 0.050)
  coloc_attempted             : 23
  coloc_with_loci             : 21
  coloc_artefacts_dir         : 21
  skip-reason breakdown:
    not_nominally_significant              373
    coloc_business_logic_returned_empty      2
```

If a run produces no coloc artefacts, this file (plus `coloc:` log lines
from `coloc_business_logic`) is your first stop.

### Tuning knobs

| arg | default | meaning |
|---|---|---|
| `coloc_window_kb` | `500L` | half-window around each IV |
| `coloc_priors`    | `list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)` | coloc.abf priors |
| `coloc_skip_mhc`  | `TRUE` | skip chr6:25–35 Mb (long-range LD breaks coloc.abf) |

## Working with outcome providers

### Pan-UK Biobank (sex = "both")

- Outcomes are fetched on the fly using the AWS HTTPS endpoints recorded in the
  packaged `panukb_pheno_manifest` dataset.
- `panukb_snp_grabber()` caches decompressed tabix indices and streams the
  required SNP rows only, keeping disk usage manageable.
- Outcomes that do not list the requested ancestry in their population column
  are automatically excluded.

### Neale Lab (sex = "male" / "female")

- Sex-specific manifests (`neale_male_manifest`, `neale_female_manifest`) are
  joined with the Neale file manifest to resolve download URLs.
- `neale_gwas_checker()` ensures the corresponding `.tsv.bgz` files exist under
  `<cache_dir>/neale_sumstats`. Set `confirm = "yes"` to allow large batch
  downloads in non-interactive sessions.
- `neale_snp_grabber()` performs exact variant matches within the downloaded
  files and caches per-outcome SNP tables to avoid repeated parsing.

## Packaged data resources

The `data/` directory ships several tibble-formatted datasets that can be
inspected with `utils::data()` or `get_pkg_obj()`:

- `bothsex_all`, `male_all`, `female_all`: outcome catalogs with ARD flags,
  ICD-10 codes, and GBD metadata.
- `panukb_pheno_manifest`: Pan-UKB phenotype manifest, including AWS URLs and
  ancestry coverage.
- `neale_male_manifest`, `neale_female_manifest`: Neale phenotype tables with
  case/control counts.
- `neale_file_manifest`: file-level metadata (download links, checksums) for
  Neale sumstats.

## Additional utilities

- `run_phenome_mr_plotting_only()` regenerates plots from a previous run without
  re-downloading or re-running MR analyses.
- `ard_compare()` and `run_ieugwasr_ard_compare()` contrast ARDMR findings
  against IEU OpenGWAS or other external MR results.
- `compute_ld_matrix()` computes pairwise correlations between exposure
  instruments using PLINK 1.9 against the bundled 1000G v3 reference. Called
  automatically inside `run_phenome_mr()` (writes `exposure_snps_ld_matrix.csv`
  to the run directory); also callable on its own.
- `coloc_fetch_exposure_region()`, `coloc_fetch_outcome_panukb_region()`,
  `coloc_panukb_outcome_metadata()` — building blocks for the coloc data
  fetchers; useful if you want to drive coloc directly without the
  `run_phenome_mr()` wrapper.
- `validate_gwasvcf()` and `read_gwasvcf_region()` — light wrappers for the
  GWAS-VCF format used as the on-disk exposure-sumstats source.
- `clump_sumstats_to_ivs()` — when `exposure_sumstats` is supplied without
  `exposure_snps`, derives IVs by p-thresholding, LD-clumping (against 1000G
  v3) and F-filtering the VCF (`gwasvcf` package required).
- `enrichment_business_logic.R` performs over-representation testing of GBD
  causes among significant discoveries.
- `logging.R` sets up structured logging via the `logger` package so runs can be
  audited later.
- `dev/integration_test_ieugwasr.R` is a minimal end-to-end smoke test using
  `n_pheno_limit` and a tiny exposure CSV; safe to run for ~5 outcomes to
  confirm the pipeline is intact after a code change.

### Phenoscanner deprecation and IEU OpenGWAS PheWAS support

The optional Phenoscanner lookups used for post-hoc instrument screening are now
**deprecated** because the public Phenoscanner v2 API is no longer maintained
reliably. Calls remain available for backward compatibility via
`phenoscanner = TRUE`, but new workflows should avoid the dependency and rely on
locally documented exclusion patterns instead.

`run_ieugwasr_ard_compare()` now bundles a PheWAS helper that queries
`ieugwasr::phewas()` to screen instruments across the IEU OpenGWAS catalogue.
This integrates with the existing grouped comparison pipeline so that exposure
CSV manifests can automatically record PheWAS hits, filter variants, and write
audit trails alongside the MR outputs.

## Development and testing

Run the unit test suite after making changes:

```bash
Rscript -e "devtools::test()"
```

The repository uses `testthat` (edition 3). Tests rely on small fixtures so they
can execute without downloading large external data files.

When contributing code:

- Follow the tidyverse style conventions used throughout `R/`.
- Ensure new features respect the cache directory contract and avoid writing to
  the working directory.
- Update documentation and examples as needed.

## License

The package is released under the MIT License. See [`LICENSE.md`](./LICENSE.md)
for details.

## Getting help

Open an issue or pull request if you encounter bugs, have questions about the
pipeline, or wish to propose enhancements. Contributions that improve
reproducibility, portability, or test coverage are especially welcome.
