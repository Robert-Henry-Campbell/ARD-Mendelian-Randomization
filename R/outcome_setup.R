#' Build MR_df (outcome plan) from packaged manifests
#' @param sex "both","male","female"
#' @param ancestry character
#' @param verbose logical
#' @return tibble MR_df (one row per outcome), with columns: provider, provider_key, icd10, gbd_cause, file_path (if known), etc.
#' @export
outcome_setup <- function(sex, ancestry, verbose = TRUE) {
  # TODO: validate sex/ancestry; join your ARD→ICD10→GBD map with panukb/neale manifests
  # For now, return a tiny placeholder tibble
  tibble::tibble(
    outcome_id = character(),
    provider = character(),      # "panukb" or "neale"
    provider_key = character(),  # e.g., file stem or phenocode
    icd10 = character(),
    gbd_cause = character(),
    file_path = character()
  )
}
