# Delete the panukb_outcome_snps and neale_outcome_snps caches under
# $ARDMR_CACHE_DIR. These hold IV-hashed query results (.rds) that the
# pipeline regenerates on next run. Does NOT touch raw server downloads
# (manifests, sumstats, ld_reference, tabix_index, hts-cache).
#
# Usage: source("dev/delete_cache_outcome_snps.R") from the package root.

stopifnot(file.exists("DESCRIPTION"))
cache_dir <- Sys.getenv("ARDMR_CACHE_DIR", unset = "")
stopifnot(nzchar(cache_dir), dir.exists(cache_dir))

for (t in file.path(cache_dir, c("panukb_outcome_snps", "neale_outcome_snps"))) {
  if (dir.exists(t)) {
    n <- length(list.files(t, recursive = TRUE, all.files = TRUE, no.. = TRUE))
    unlink(t, recursive = TRUE, force = TRUE)
    message(sprintf("Deleted %s (%d files)", t, n))
  } else {
    message(sprintf("Skipped (does not exist): %s", t))
  }
}
