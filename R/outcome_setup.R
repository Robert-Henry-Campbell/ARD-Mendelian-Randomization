# R/Outcome_setup.R

#' Build MR_df (outcome plan) by mapping ARD phenotypes to provider manifests
#'
#' @description
#' Creates the outcome plan **MR_df** for a phenome-wide MR run by:
#' 1) selecting the ARD phenotype tibble for the requested `sex`
#' 2) logging counts at ICD-10 and GBD cause level
#' 3) mapping `ICD10_explo` to the appropriate provider manifest
#'    - Pan-UKB if `sex == "both"`
#'    - Neale if `sex %in% c("male","female")`
#' 4) dropping unmapped phenotypes (no available GWAS)
#' 5) returning a row per outcome with ARD + provider metadata
#'
#' @param sex Character. One of `"both"`, `"male"`, `"female"`.
#' @param ancestry Character. One of allowed ancestries (e.g., `"EUR"`, `"AFR"`, `"AMR"`, `"EAS"`, `"SAS"`).
#'                 **Constraint:** if `sex` is `"male"` or `"female"`, `ancestry` must be `"EUR"` (Neale).
#' @param verbose Logical. Print summary logs.
#'
#' @return A tibble `MR_df` with one row per outcome. Contains ARD fields
#' (e.g., `ICD10_explo`, `gbd_cause`) plus provider fields (e.g., `provider`,
#' `provider_key`, `download_url`/`local_filename`, column maps, etc.).
#' @export
Outcome_setup <- function(sex, ancestry, verbose = TRUE) {
  # ---- validations ----
  check_choice(sex, SEX_OPTIONS, "sex")
  check_choice(ancestry, ANCESTRY_OPTIONS, "ancestry")
  if (sex %in% c("male", "female") && ancestry != "EUR") {
    stop("Sex-specific (male/female) runs currently require ancestry == 'EUR' (Neale).", call. = FALSE)
  }

  # ---- select ARD tibble for the requested sex ----
  # Expect packaged data objects: bothsex_ARD, male_ARD, female_ARD
  ard_df <- .select_ard_by_sex(sex)

  # Canonicalize expected ARD columns
  # Required: ICD10_explo, gbd_cause
  # Allow some common variants and rename if present
  rename_map <- c(
    "ICD10" = "ICD10_explo",
    "icd10" = "ICD10_explo",
    "icd10_explo" = "ICD10_explo",
    "GBD_cause" = "gbd_cause",
    "gbd" = "gbd_cause"
  )
  for (from in names(rename_map)) {
    to <- rename_map[[from]]
    if (from %in% names(ard_df) && !to %in% names(ard_df)) {
      ard_df <- dplyr::rename(ard_df, !!to := .data[[from]])
    }
  }

  assert_cols(ard_df, c("ICD10_explo", "gbd_cause"), "ARD phenotype tibble")

  # ---- basic logs ----
  if (verbose) {
    n_rows <- nrow(ard_df)
    n_icd  <- dplyr::n_distinct(ard_df$ICD10_explo)
    n_gbd  <- dplyr::n_distinct(ard_df$gbd_cause)
    log_info("Outcome_setup: %d phenotypes | %d unique ICD10_explo | %d unique GBD causes", n_rows, n_icd, n_gbd)
  }

  # ---- choose provider & load manifest ----
  provider <- if (sex == "both") "panukb" else "neale"

  if (provider == "panukb") {
    man <- load_panukb_manifest()
    # Expected minimal columns in manifest:
    #   icd10_explo, provider_key, download_url OR tabix_url,
    #   chr_col, pos_col, ea_col, oa_col, beta_col, se_col, p_col, (optional ancestry column(s))
    assert_cols(man, c("icd10_explo", "provider_key"), "Pan-UKB manifest")
    by_cols <- c("ICD10_explo" = "icd10_explo")
  } else {
    man <- load_neale_manifest()
    # Expected minimal columns in manifest:
    #   icd10_explo, provider_key, local_filename OR download_cmd,
    #   chr_col, pos_col, ea_col, oa_col, beta_col, se_col, p_col
    assert_cols(man, c("icd10_explo", "provider_key"), "Neale manifest")
    by_cols <- c("ICD10_explo" = "icd10_explo")
  }

  # ---- join ARD to manifest on ICD10_explo ----
  # Perform a case-insensitive join by normalizing the join key
  ard_key <- ard_df |>
    dplyr::mutate(.key_ard = tolower(.data$ICD10_explo))

  man_key <- man |>
    dplyr::mutate(.key_man = tolower(.data[[unname(by_cols)] ]))

  joined <- ard_key |>
    dplyr::left_join(
      dplyr::select(man_key, dplyr::all_of(c(".key_man", names(man)))),
      by = c(".key_ard" = ".key_man")
    ) |>
    dplyr::select(-.key_ard)

  # ---- mapping diagnostics & drop unmapped ----
  mapped_flag <- !is.na(joined$provider_key)
  n_total  <- nrow(joined)
  n_mapped <- sum(mapped_flag, na.rm = TRUE)
  n_unmap  <- n_total - n_mapped

  if (verbose) {
    log_info("Provider = %s | mapped = %d / %d | unmapped = %d", provider, n_mapped, n_total, n_unmap)
    if (n_unmap > 0) {
      unmapped_preview <- joined$ICD10_explo[which(!mapped_flag)][seq_len(min(5, n_unmap))]
      log_info("Unmapped (first up to 5): %s", paste(unmapped_preview, collapse = ", "))
    }
  }

  MR_df <- dplyr::filter(joined, !!mapped_flag)

  # ---- finalize shape & provider column ----
  MR_df <- MR_df |>
    dplyr::mutate(provider = provider) |>
    # create a stable outcome_id if not present
    dplyr::mutate(outcome_id = dplyr::coalesce(.data$provider_key, .data$ICD10_explo)) |>
    dplyr::relocate(outcome_id, provider, provider_key, ICD10_explo, gbd_cause)

  # ---- extra checks (extensible) ----
  # Ensure manifest core columns exist post-join (soft assert: warn rather than stop)
  expected_cols <- c("chr_col", "pos_col", "ea_col", "oa_col", "beta_col", "se_col", "p_col")
  missing_manifest_cols <- setdiff(expected_cols, names(MR_df))
  if (length(missing_manifest_cols) && verbose) {
    log_info("Warning: Missing expected manifest columns in MR_df: %s", paste(missing_manifest_cols, collapse = ", "))
  }

  MR_df
}

# ---- internal helpers ----

# Select ARD tibble by sex. Expects packaged data objects named:
# bothsex_ARD, male_ARD, female_ARD (with columns including ICD10_explo, gbd_cause).
.select_ard_by_sex <- function(sex) {
  get_pkg_obj <- function(name) {
    obj <- get0(name, envir = asNamespace(utils::packageName()), inherits = FALSE)
    if (is.null(obj)) stop(sprintf("Packaged data object '%s' not found. Add it to data/ via data-raw/.", name), call. = FALSE)
    dplyr::as_tibble(obj)
  }
  switch(
    sex,
    "both"   = get_pkg_obj("bothsex_ARD"),
    "male"   = get_pkg_obj("male_ARD"),
    "female" = get_pkg_obj("female_ARD"),
    stop("Invalid sex option.", call. = FALSE)
  )
}

# Note: load_panukb_manifest() and load_neale_manifest() are expected to be
# defined elsewhere (e.g., in load_data helpers). They must return tibbles
# with at least:
#   icd10_explo, provider_key, and the column map fields used downstream:
#   chr_col, pos_col, ea_col, oa_col, beta_col, se_col, p_col,
#   plus either (download_url/tabix_url) for Pan-UKB or (local_filename/download_cmd) for Neale.
