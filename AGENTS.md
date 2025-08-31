\# AGENTS.md – Guidance for AI Coding Agents



This file guides AI agents (like OpenAI Codex) working on the \*\*ARDMR\*\* package. Think of it as a README for AI, not humans — follow it precisely to produce consistent, high-quality code.



---



\## Project Overview

\- \*\*Purpose:\*\* Automate Mendelian Randomization (one-exposure to many-outcomes) across age-related diseases (ARDs), with automatic SNP selection, sensitivity/QC checks, multiple testing, and visualization.

\- \*\*Language:\*\* R with tidyverse + TwoSampleMR.



\### Key modules:

1\. `Run\_phenome\_mr.R` – orchestrator; runs the full pipeline, outputs results + plots.

2\. `Outcome\_setup.R` – builds MR\_df by mapping ARD ICD10s to GWAS manifests (Pan-UKB or Neale).

3\. `Exposure\_SNP\_mapper.R` – maps exposure RSIDs to chr:pos via variant manifests.

4\. `Panukb\_SNP\_grabber.R` – fetches outcome SNP data via remote tabix slices.

5\. `Neale\_GWAS\_checker.R` – ensures Neale sumstats files are present/indexed locally.

6\. `Neale\_tbi\_maker.R` – indexes Neale bgz files with .tbi.

7\. `Neale\_SNP\_grabber.R` – fetches outcome SNP data from local bgz files.

8\. `MR\_business\_logic.R` – harmonizes, runs MR + sensitivity checks (8-check QC), outputs results.

9\. `Volcano\_plot.R` \& `Manhattan\_plot.R` – visualizes results across outcomes.



---



\## Build \& Test



\*\*To run the full pipeline end-to-end (on real or fixture data):\*\*

```bash

Rscript -e "devtools::load\_all(); run\_phenome\_mr(exposure\_snps = …, sex = 'both', ancestry = 'EUR', plot\_output\_dir = 'output/')"




## Packaged data files

The package ships with several `.rda` files (located in `data/`) used by `Outcome_setup()`:

- `bothsex_ARD.rda`, `female_ARD.rda`, `male_ARD.rda` – ARD phenotype tables for both sexes, females, and males.
- `neale_female_manifest.rda`, `neale_male_manifest.rda` – Neale sex-specific phenotype manifests containing case and control counts.
- `neale_file_manifest.rda` – Neale download information for each GWAS.
- `panukb_pheno_manifest.rda` – Pan-UKB phenotype manifest with metadata and download links.
