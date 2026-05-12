# =============================================================================
# ardmr template: warm the cache (one-shot setup downloader)
# =============================================================================
#
# Use this once after pointing ARDMR_CACHE_DIR at a fresh scratch directory,
# to pre-fetch every large remote resource the orchestrator
# (run_phenome_mr / ard_compare / run_ieugwasr_ard_compare) will otherwise
# download lazily during a real run.
#
# What this script downloads
#   1. Variant manifests (RSID -> chr:pos lookup), one per catalog:
#        <cache>/neale/variants.tsv.bgz                       (~250 MB)
#        <cache>/panukb/full_variant_qc_metrics.txt.bgz       (~250 MB)
#      Source of truth: Variant_manifest_downloader().
#
#   2. 1kg.v3 LD reference bundle (one .tgz contains all 5 superpopulations):
#        <cache>/ld_reference/{EUR,AFR,EAS,SAS,AMR}.{bed,bim,fam}
#        ~1.5 GB compressed download, ~5 GB extracted.
#      Source of truth: ld_reference_checker().
#
#   3. (Optional, gated) Neale sex-stratified sumstats:
#        <cache>/neale_sumstats/*.tsv.bgz                     (~500 GB FULL SET)
#      Only needed if you intend to run sex = "male" or "female".
#      The set actually pulled is the union of files referenced by the
#      male + female outcome plans for ancestry = "EUR" (the only ancestry
#      Outcome_setup() permits for sex-stratified runs).
#
#   4. (Optional, gated) Pan-UKB tabix indices (.tbi) for every mappable
#      Pan-UKB outcome at the requested ancestry:
#        <cache>/tabix_index/*.tbi                            (a few GB total)
#      Cuts the per-phenotype "first hit" latency in sex = "both" runs.
#
# What this script does NOT download
#   * OpenGWAS / ieugwasr data: there is no static cache to warm; these are
#     fetched per query at runtime and require OPENGWAS_JWT.
#   * GWAS-VCFs for coloc (exposure_sumstats): user-supplied per run.
#   * Pan-UKB outcome SNP rows (.rds in panukb_outcome_snps/): IV-set keyed,
#     populated lazily at run time.
#
# Cost: with both optional sections OFF this is ~2 GB of network traffic and
# a few minutes. With Neale ON expect many hours and ~500 GB on disk.
# -----------------------------------------------------------------------------

#setwd("/mnt/sdg/robert/ardmr/ARD-Mendelian-Randomization")          # path to package source
devtools::load_all()

Sys.setenv(ARDMR_CACHE_DIR = "/mnt/sdg/robert/ardmr/ardmr_cache")     # tens-of-GB scratch dir

# ---- toggles ----------------------------------------------------------------
download_neale_sumstats     <- TRUE      # ~500 GB, many hours
download_panukb_tabix_index <- TRUE      # ~few GB, ~30 min
panukb_ancestry             <- "EUR"      # ancestry to pre-fetch tabix .tbi for

cache_dir <- ardmr_cache_dir()
message(sprintf("Cache root: %s", normalizePath(cache_dir, mustWork = FALSE)))

# ---- 1. Variant manifests (both catalogs) -----------------------------------
Variant_manifest_downloader("neale",  cache_dir = cache_dir, overwrite = FALSE)
Variant_manifest_downloader("panukb", cache_dir = cache_dir, overwrite = FALSE)

# ---- 2. 1kg.v3 LD reference bundle (all superpopulations in one go) ---------
# The bundle ships EUR/AFR/EAS/SAS/AMR together; one call extracts all of them
# and a second call for any other pop is a no-op once the .bed/.bim/.fam exist.
ld_ref_dir <- file.path(cache_dir, "ld_reference")
ld_reference_checker(
  pop         = "EUR",
  ld_ref_dir  = ld_ref_dir,
  verbose     = TRUE,
  confirm     = "yes"
)
# Sanity-confirm the other pops were extracted from the same bundle. These
# calls do NOT redownload; they just verify the .bed/.bim/.fam are in place.
for (pop in c("AFR", "EAS", "SAS")) {
  ld_reference_checker(pop = pop, ld_ref_dir = ld_ref_dir,
                       verbose = TRUE, confirm = "no")
}

# ---- 3. (Optional) Neale sex-stratified sumstats ----------------------------
if (isTRUE(download_neale_sumstats)) {
  neale_dir <- file.path(cache_dir, "neale_sumstats")
  dir.create(neale_dir, recursive = TRUE, showWarnings = FALSE)

  # Build the union of files referenced by the male + female outcome plans
  # at ancestry = "EUR" (Outcome_setup forbids any other ancestry here).
  mr_male   <- Outcome_setup(sex = "male",   ancestry = "EUR")
  mr_female <- Outcome_setup(sex = "female", ancestry = "EUR")
  mr_neale  <- dplyr::bind_rows(mr_male, mr_female)
  mr_neale  <- dplyr::distinct(mr_neale, File, .keep_all = TRUE)

  neale_gwas_checker(
    MR_df       = mr_neale,
    neale_dir   = neale_dir,
    verbose     = TRUE,
    confirm     = "yes",
    estimate_gb = 500
  )
}

# ---- 4. (Optional) Pan-UKB tabix indices for mappable phenotypes ------------
if (isTRUE(download_panukb_tabix_index)) {
  idx_cache_dir <- file.path(cache_dir, "tabix_index")
  dir.create(idx_cache_dir, recursive = TRUE, showWarnings = FALSE)

  mr_panukb <- Outcome_setup(sex = "both", ancestry = panukb_ancestry)
  if (!all(c("aws_link_tabix", "description") %in% names(mr_panukb))) {
    stop("Outcome_setup() did not return aws_link_tabix/description; cannot pre-cache .tbi.",
         call. = FALSE)
  }

  .slug <- function(x) {
    x <- as.character(x); x[is.na(x) | !nzchar(x)] <- "NA"
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    substr(x, 1, 120)
  }
  .clean_url <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^(['\"])(.*)\\1$", "\\2", x, perl = TRUE)
    if (grepl("^s3://", x)) {
      x <- sub("^s3://([^/]+)/(.*)$", "https://\\1.s3.amazonaws.com/\\2", x)
    }
    x
  }

  old_to <- getOption("timeout"); options(timeout = max(3600, old_to))
  on.exit(options(timeout = old_to), add = TRUE)

  urls   <- vapply(mr_panukb$aws_link_tabix, .clean_url, character(1))
  labels <- as.character(mr_panukb$description)
  keep   <- nzchar(urls) & grepl("^(https?|s3)://", urls) &
            nzchar(labels) & !is.na(labels)
  urls   <- urls[keep]; labels <- labels[keep]

  total <- length(urls)
  message(sprintf("Pan-UKB tabix index: pre-caching .tbi for %d phenotypes (%s)",
                  total, panukb_ancestry))

  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)

  n_present <- 0L; n_downloaded <- 0L; n_failed <- 0L
  for (k in seq_len(total)) {
    idx_hash  <- tryCatch(digest::digest(unname(urls[k]), algo = "xxhash64"),
                          error = function(e) digest::digest(urls[k]))
    dest <- file.path(idx_cache_dir,
                      paste0(.slug(labels[k]), "_", idx_hash, ".tbi"))
    if (file.exists(dest)) { n_present <- n_present + 1L; utils::setTxtProgressBar(pb, k); next }
    ok <- tryCatch({
      utils::download.file(urls[k], destfile = dest, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (ok && file.exists(dest)) n_downloaded <- n_downloaded + 1L
    else { n_failed <- n_failed + 1L; if (file.exists(dest)) try(unlink(dest, force = TRUE), silent = TRUE) }
    utils::setTxtProgressBar(pb, k)
  }
  message(sprintf("\nPan-UKB tabix index: present=%d, downloaded=%d, failed=%d",
                  n_present, n_downloaded, n_failed))
}

# ---- summary ----------------------------------------------------------------
print(ardmr_cache_summary(cache_dir))
