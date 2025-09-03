#R/outcome_setup.R
#' Build the outcome plan (`MR_df`) by mapping ARD phenotypes to GWAS manifests
#'
#' @description
#' Selects the appropriate age-related disease (ARD) phenotype table based on
#' `sex`, logs counts at the ICD-10 and cause levels, then maps the
#' `ICD10_explo` column to the relevant GWAS manifest. Pan-UKB is used when
#' `sex == "both"`; otherwise Neale manifests are consulted. Phenotypes with no
#' available GWAS are dropped.
#'
#' @param sex Character. One of `"both"`, `"male"`, or `"female"`.
#' @param ancestry Character ancestry code. If `sex` is `"male"` or
#'   `"female"`, this must be `"EUR"`.
#'
#' @return A tibble `MR_df` combining the ARD data with GWAS manifest
#'   metadata.
#' @export
Outcome_setup <- function(sex, ancestry) {
  # ---- validations ---------------------------------------------------------
  sex <- match.arg(tolower(sex), c("both", "male", "female"))
  ancestry <- match.arg(toupper(ancestry), c("AFR","AMR","CSA","EAS","EUR","MID"))

  if (sex %in% c("male", "female") && ancestry != "EUR") {
    stop("Sex-specific (male/female) runs currently require ancestry == 'EUR'", call. = FALSE)
  }

  # ---- select ARD tibble ---------------------------------------------------
  ard_df <- switch(
    sex,
    "both"   = get_pkg_obj("bothsex_ARD"),
    "male"   = get_pkg_obj("male_ARD"),
    "female" = get_pkg_obj("female_ARD")
  )

  logger::log_info(
    "ARD phenotypes| {nrow(ard_df)} rows; {dplyr::n_distinct(ard_df$ICD10_explo)} unique ICD10; {dplyr::n_distinct(ard_df$cause_level_3)} unique causes"
  )

  # ---- load manifest -------------------------------------------------------
  if (sex == "both") {
    manifest <- get_pkg_obj("panukb_pheno_manifest")
    join_col <- "phenocode"
  } else {
    file_man <- get_pkg_obj("neale_file_manifest")
    sex_man <- if (sex == "male") get_pkg_obj("neale_male_manifest") else get_pkg_obj("neale_female_manifest")
    manifest <- dplyr::left_join(sex_man, file_man, by = c("phenotype" = "Phenotype Code"))
    join_col <- "phenotype"
  }

  # ---- map ICD10 to manifest ----------------------------------------------
  MR_df <- dplyr::left_join(ard_df, manifest, by = c("ICD10_explo" = join_col))

  mapped <- !is.na(MR_df[[ncol(MR_df)]]) #hand editted
  n_total <- nrow(MR_df)
  n_mapped <- sum(mapped)
  logger::log_info("Mapped {n_mapped} of {n_total} phenotypes; {n_total - n_mapped} without available GWAS")

  MR_df <- MR_df[!is.na(MR_df[[ncol(MR_df)]]), ] #for some reason using mapped isn't ok here

  MR_df
}

