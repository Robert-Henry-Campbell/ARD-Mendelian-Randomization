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
3. Retrieves outcome SNP rows either from Pan-UK Biobank (sex = "both") or local
   Neale Lab summary statistics (sex = "male"/"female").
4. Runs inverse-variance weighted MR along with eight supporting sensitivity and
   quality-control checks.
5. Produces tidy result tables and ggplot2-based volcano and Manhattan plots and
   writes run artefacts to disk for reproducibility.

## Key capabilities

- **Automated phenotype manifests.** Packaged datasets (`bothsex_all`,
  `male_all`, `female_all`) link ICD-10 codes to GBD causes, Pan-UKB phenocodes,
  and Neale Lab identifiers while tracking which outcomes are considered ARDs.
- **Provider-aware SNP mapping.** `exposure_snp_mapper()` joins user-supplied
  instruments to Pan-UKB or Neale variant manifests and handles multi-allelic
  loci.
- **Remote-friendly Pan-UKB queries.** `panukb_snp_grabber()` streams tabix
  slices over HTTPS via Bioconductor's `Rsamtools`, so a system `tabix`
  executable is not required.
- **Neale GWAS management.** `neale_gwas_checker()` verifies large Neale Lab
  files locally, optionally downloading missing datasets with resumable logging.
- **QC-rich MR engine.** `mr_business_logic()` harmonises, runs multiple MR
  estimators, and records eight sensitivity checks with configurable pass
  thresholds.
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
  exposure_units = "1-SD increase",
  ancestry = "EUR",
  sex = "both",             # "both" uses Pan-UKB, "male"/"female" use Neale
  sensitivity_pass_min = 6,  # require six of the eight checks
  Multiple_testing_correction = "BH"
)

# Inspect outputs
results$results_df   # tidy phenome-wide MR summary
results$volcano      # ggplot object; can be customised further
results$manhattan    # ggplot object for -log10(p) across outcomes
```

Each run writes a timestamped log file, diagnostic plots, and CSV exports to:

```
<cache_dir>/output/<exposure>/<sex>/<ancestry>/
```

and caches variant lookups under provider-specific subdirectories so subsequent
runs are faster.

## Quality-control checks

Eight sensitivity diagnostics are calculated for every outcome when data permit:

1. Cochran's Q (heterogeneity) p-value ≥ 0.05
2. Cochran's Q I² ≤ 50%
3. MR-Egger intercept p-value ≥ 0.05
4. MR-Egger slope sign agrees with IVW (unless I²<sub>GX</sub> < 0.6)
5. I²<sub>GX</sub> ≥ 0.6
6. Weighted median estimate within one IVW standard error
7. Weighted mode estimate within one IVW standard error
8. Leave-one-out analysis shows no sign-flip and Steiger filtering failures ≤ 30%

`sensitivity_pass_min` defines how many of the enabled checks an outcome must
pass to be considered QC-approved. This flag feeds into the default multiple
testing correction and plot aesthetics.

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
- `ard_compare()` and `ard_compare_grouped_ieugwasr()` contrast ARDMR findings
  against IEU OpenGWAS or other external MR results.
- `enrichment_business_logic.R` performs over-representation testing of GBD
  causes among significant discoveries.
- `logging.R` sets up structured logging via the `logger` package so runs can be
  audited later.

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
