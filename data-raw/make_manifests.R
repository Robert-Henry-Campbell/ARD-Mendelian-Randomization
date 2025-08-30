# data-raw/make_manifests.R

library(readr)
library(usethis)

# 1) Raw files live OUTSIDE the package tarball.
#    Put them under data-raw/ and keep data-raw/ in .Rbuildignore (default).
data_dir <- "G:/My Drive/Documents/0Oxford_main/ARD paper/ardmr/data-raw"

# Helper: assert file exists
req <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  path
}

# ---- Pan-UKB (CSV) ----
panukb_csv <- file.path(data_dir, "PanUKB.phenotype.manifest.trimmed.csv")
req(panukb_csv)
panukb_pheno_manifest <- read_csv(panukb_csv, show_col_types = FALSE)

# ---- Neale manifests ----
neale_file_csv   <- file.path(data_dir, "Neale.File.Manifest.Release20180731.trimmed.csv")
neale_male_bgz   <- file.path(data_dir, "neale.phenotypes.male.tsv.bgz")
neale_female_bgz <- file.path(data_dir, "neale.phenotypes.female.tsv.bgz")  # fixed extension

req(neale_file_csv)
req(neale_male_bgz)
req(neale_female_bgz)

# File manifest (CSV)
neale_file_manifest <- read_csv(neale_file_csv, show_col_types = FALSE)

# Phenotype manifests (TSV, bgzip)
# If these are huge, consider selecting columns to reduce package size, e.g.:
# cols_keep <- c("phenocode","description","n_cases","n_controls", "sex", "category")
# neale_male_manifest   <- read_tsv(neale_male_bgz,   col_select = all_of(cols_keep), show_col_types = FALSE)
# neale_female_manifest <- read_tsv(neale_female_bgz, col_select = all_of(cols_keep), show_col_types = FALSE)

neale_male_manifest   <- read_tsv(neale_male_bgz,   show_col_types = FALSE)
neale_female_manifest <- read_tsv(neale_female_bgz, show_col_types = FALSE)

# ---- Save as compressed R objects under data/ ----
use_data(
  panukb_pheno_manifest,
  neale_file_manifest,
  neale_male_manifest,
  neale_female_manifest,
  overwrite = TRUE,
  compress  = "xz"
)

# ---- Report sizes on disk ----
rda_files <- list.files("data", pattern = "\\.rda$", full.names = TRUE)
info <- file.info(rda_files)
sizes <- data.frame(
  file   = basename(rownames(info)),
  sizeMB = round(info$size / 1024^2, 2)
)
print(sizes)
