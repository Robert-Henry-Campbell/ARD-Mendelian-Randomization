pkg_path  <- Sys.getenv("ARDMR_PKG_PATH",  unset = NA_character_)
cache_dir <- Sys.getenv("ARDMR_CACHE_DIR", unset = NA_character_)
jwt       <- Sys.getenv("OPENGWAS_JWT",    unset = NA_character_)

if (is.na(pkg_path)  || !nzchar(pkg_path))  stop("ARDMR_PKG_PATH not set")
if (is.na(cache_dir) || !nzchar(cache_dir)) stop("ARDMR_CACHE_DIR not set")
if (is.na(jwt)       || !nzchar(jwt))       stop("OPENGWAS_JWT not set")

setwd(pkg_path)
devtools::load_all()

csv_path <- file.path("inst", "extdata", "19_batches",
                      "batch_17_methylation_telomere_imaging_GWAS.csv")
stopifnot(file.exists(csv_path))

result <- run_ieugwasr_ard_compare(
  csv_path             = csv_path,
  prompt_for_units     = FALSE,
  p_threshold          = 5e-8,
  p_backoff            = NULL,
  r2                   = 0.001,
  kb                   = 10000,
  f_threshold          = 10,
  phewas_pval          = 5e-8,
  phenoscanner         = FALSE,
  sensitivity_pass_min = 5,
  force_refresh        = FALSE,
  n_pheno_limit        = NULL
)

dir.create(file.path(cache_dir, "logs_19batch"), showWarnings = FALSE, recursive = TRUE)
saveRDS(result, file.path(cache_dir, "logs_19batch",
                          sprintf("batch_%02d_result.rds", 17L)))
