# Manual integration test for run_ieugwasr_ard_compare()
# Truncates outcomes to 10 phenotypes via n_pheno_limit so the run finishes
# in minutes, not hours. NOT a science run -- p-values, BH/Bonferroni flags,
# enrichment tests, and plot thresholds are computed against N=10 and are
# meaningless. Use this only to catch crashes / regressions in the pipeline.
#
# Usage: source("dev/integration_test_ieugwasr.R") from the package root.
rm(list = ls())
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJoLnJvYmVydC5jYW1wYmVsbEBnbWFpbC5jb20iLCJpYXQiOjE3Nzc5Nzc3OTIsImV4cCI6MTc3OTE4NzM5Mn0.n--pL9YObsrKDfMPEOMoXQAlyPcQ-CuQcQiQBeMCZ_gBfrXuunNcEuyyqjAGeDP0NG_knFZQtzYPPz2zqCSpe_zRuPi2MKy64nJy2ZJVy0nYhokUWk7henI8autyo_g2k-Z1_J56cU1CaoKDTxTerNTbySwmY4KY71FkK3QYXE222jZErUivC79M_4DUEz4Vq0wlzXmWSwEAcfv2B_j4CLjEVHbRbRMSTm03RTmZK-4KSPlXAPLDAwSWbvH8s7XAzNKC4Os6-qzmg9txHEy1GX4AyJwnrSX8xNWP9ze_HCRDETdsFodXUb2LRBHOamby8f96YMbhwwcWsxkMDhxp5Q")
#set cache dir
Sys.setenv(ARDMR_CACHE_DIR = 'C:\\Users\\Robert\\Downloads\\ardmr')


stopifnot(file.exists("DESCRIPTION"))  # must run from package root

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required. install.packages('devtools').")
}

jwt <- Sys.getenv("IEU_JWT", unset = "")
if (!nzchar(jwt)) jwt <- Sys.getenv("OPENGWAS_JWT", unset = "")
if (!nzchar(jwt)) {
  stop("Set IEU_JWT or OPENGWAS_JWT env var before running this script.")
}

devtools::load_all(".", quiet = TRUE)

csv_path <- "G:\\My Drive\\Documents\\0Oxford_main\\ARD paper\\3_ARD_MR\\exposures to run on\\test_phenos\\tnfA_goodtests.csv"
Sys.setenv(ARDMR_CACHE_DIR = 'C:\\Users\\Robert\\Downloads\\ardmr')


message("=== integration test ===")
message("csv_path : ", csv_path)
message("cache_dir: ", ardmr_cache_dir())
message("n_pheno_limit: 10")
message("--- starting run ---")

t0 <- Sys.time()
out <- run_ieugwasr_ard_compare(
  csv_path         = csv_path,
  cache_dir        = ardmr_cache_dir(),
  prompt_for_units = FALSE,
  phenoscanner     = FALSE,
  force_refresh    = FALSE,
  n_pheno_limit    = 3
)
t1 <- Sys.time()
elapsed <- difftime(t1, t0, units = "mins")
message(sprintf("--- finished in %.2f min ---", as.numeric(elapsed)))

# ---- post-run sanity checks ----
slug <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "NA"
  x <- gsub("[^[:alnum:]]+", "-", x)
  x <- gsub("-+", "-", x)
  x <- gsub("^-|-$", "", x)
  tolower(x)
}

exposure_slug <- slug("LDL_cholesterol")
group_dir <- file.path(cache_dir, "output", exposure_slug, "both", "EUR")
results_rds <- file.path(group_dir, "results.rds")

message("\n=== sanity checks ===")
message("group_dir exists           : ", dir.exists(group_dir))
message("results.rds exists         : ", file.exists(results_rds))
if (file.exists(results_rds)) {
  res <- readRDS(results_rds)
  res_df <- if (is.data.frame(res)) res else res$results_df
  if (!is.null(res_df)) {
    message(sprintf("results_df rows (expect 10): %d", nrow(res_df)))
  }
}
message("\nDone. Inspect logs in: ", cache_dir)
