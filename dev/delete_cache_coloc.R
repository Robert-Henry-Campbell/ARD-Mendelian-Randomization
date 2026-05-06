# Delete the coloc/ cache under $ARDMR_CACHE_DIR (covers exposure/, outcome/,
# meta/). Region-keyed .rds files are regenerated from OpenGWAS / Pan-UKB / VCF
# on next coloc run. Does NOT touch raw server downloads.
#
# Usage: source("dev/delete_cache_coloc.R") from the package root.

stopifnot(file.exists("DESCRIPTION"))
cache_dir <- Sys.getenv("ARDMR_CACHE_DIR", unset = "")
stopifnot(nzchar(cache_dir), dir.exists(cache_dir))

for (t in file.path(cache_dir, "coloc")) {
  if (dir.exists(t)) {
    n <- length(list.files(t, recursive = TRUE, all.files = TRUE, no.. = TRUE))
    unlink(t, recursive = TRUE, force = TRUE)
    message(sprintf("Deleted %s (%d files)", t, n))
  } else {
    message(sprintf("Skipped (does not exist): %s", t))
  }
}
