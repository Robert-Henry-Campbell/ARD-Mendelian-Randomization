# ardmr (development version)

## Unified exposure-SNP preprocessing

Every input mode (`snps_only`, `snps_vcf`, `vcf_only`, `ieugwasr`,
`snps_id`) now routes its IV tibble through a single
[`preprocess_exposure_snps()`](R/preprocess_exposure_snps.R) pipeline
called from inside [`.resolve_coloc_source()`](R/coloc_source_resolver.R).
Previously, each mode applied a different (and partial) subset of
preprocessing steps; cross-mode IV sets were not directly comparable.

The pipeline runs (in this fixed order): p-value filter with backoff
ladder → LD clumping → rsid validation + dedupe → indel drop →
palindromic flag (and optional drop) → F-statistic filter → MAF filter
→ INFO filter. Each step records `n_in` / `n_out` counts in the new
`run_manifest.json$preprocessing_steps` field.

The full set of `clump_opts` knobs (with their defaults) is documented
at `?run_phenome_mr`. Highlights:

- **`p_backoff = c(5e-8, 5e-7, 5e-6)`** -- a ladder tried strictest-first.
  When fewer than one SNP survives the strictest rung, the next-loosest
  rung is tried. Replaces the deprecated `p_threshold` scalar (silently
  translated to `p_backoff = c(p_threshold)` with a one-time deprecation
  warning).
- **`maf_min = NULL`** -- minor-allele frequency floor; computed as
  `pmin(eaf, 1 - eaf)`. Off by default. Rows with NA EAF are kept.
- **`info_min = NULL`** -- INFO/imputation-quality floor; auto-detects
  any of `INFO`/`info`/`Rsq`/`R2`/`SI`. Off by default. Rows with NA
  INFO are kept. Same threshold applies to coloc-region fetching.
- **`drop_indels = TRUE`** -- detects multi-character or `-`/`D`/`I`
  alleles.
- **`drop_palindromic = FALSE`** -- a `palindromic` boolean column is
  added regardless; only dropped when this flag is TRUE.
- **`prefer_server_clump = TRUE`** -- in `ieugwasr` mode, OpenGWAS does
  the clumping (via `extract_instruments`) and the local clump step is
  skipped. Set to FALSE to use `tophits` + local clump instead.
- **`already_clumped`, `already_p_filtered`** -- caller hints to skip
  the corresponding local step. Used internally by ard_compare for its
  pre-curated input.
- **`preprocess = TRUE`** -- master opt-out. When FALSE, the input is
  passed through as-is.

### Removed / changed

- `clump_opts$p_threshold` is **deprecated**; pass `p_backoff` instead.
  Existing configs continue to work via auto-translation, with a
  deprecation warning at first use per run.
- The `ld_panel` field that previously appeared in some examples has
  been removed (it was a phantom -- not implemented anywhere). Local
  clumping always uses the 1kg.v3 panel.
- The CSV-driven flow `run_ieugwasr_ard_compare()` retains its
  per-row server-side instrument extraction (so PheWAS / PhenoScanner
  filtering can run on the SNPs before they reach `ard_compare()`),
  but now uses the canonical preprocess defaults for the p-value
  ladder and tells the unified preprocessor downstream to skip
  redundant clump and p-filter via `already_clumped = TRUE` and
  `already_p_filtered = TRUE`.

## Bug fixes

- **`.iv_set_hash()`**: previously called `digest::digest(serialize =
  FALSE)` on a character vector, which only hashed the first element.
  Two IV sets that shared their lexicographically-first
  `chr:pos:ea:oa` tuple aliased to the same `iv_hash` -- the run
  output directory and per-outcome cache files would silently collide
  across distinct IV sets that differed only in their tail.
  Fixed by collapsing tuples to one string before hashing. Existing
  tests continued to pass because they varied the first tuple.

## Cache-wipe advisory

Two changes in this release alter cache keys for *some* runs:

1. `iv_hash` will change for any IV set that was previously aliasing
   to a different one due to the `.iv_set_hash` bug above. This is
   rare in practice (requires shared first tuple) but means stale
   per-outcome cache files in `panukb_outcome_snps/` /
   `neale_outcome_snps/` that were pinned by the old (buggy) hash
   should be wiped.
2. `params_hash` will change for any run whose `clump_opts` previously
   omitted fields that defaults now make explicit (e.g. ard_compare's
   per-row `clump_opts` previously omitted `f_threshold`).

To clear scoped caches without touching raw downloads, run any of:

```r
source("dev/delete_cache_outcome_snps.R")  # per-outcome IV-keyed RDS files
source("dev/delete_cache_coloc.R")         # per-locus coloc artefacts
source("dev/delete_cache_ld.R")            # locally computed LD matrices
```
