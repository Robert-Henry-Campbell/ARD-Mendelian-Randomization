# Delete the ld/ cache under $ARDMR_CACHE_DIR (LD matrices computed locally
# by PLINK). Does NOT touch ld_reference/ -- those are 1000G PLINK panels
# downloaded from MRC-IEU and never modified, expensive to re-fetch.
#
# Usage: source("dev/delete_cache_ld.R") from the package root.

stopifnot(file.exists("DESCRIPTION"))
cache_dir <- Sys.getenv("ARDMR_CACHE_DIR", unset = "")
stopifnot(nzchar(cache_dir), dir.exists(cache_dir))

for (t in file.path(cache_dir, "ld")) {
  if (dir.exists(t)) {
    n <- length(list.files(t, recursive = TRUE, all.files = TRUE, no.. = TRUE))
    unlink(t, recursive = TRUE, force = TRUE)
    message(sprintf("Deleted %s (%d files)", t, n))
  } else {
    message(sprintf("Skipped (does not exist): %s", t))
  }
}
