# R/utils-ld.R
#' Map a free-form ancestry label to a 1000G LD superpopulation code
#'
#' Normalises a user-supplied ancestry string (e.g. `"EUR"`, `"european"`,
#' `"south_asian"`) to one of the four 1000G superpopulation codes used by
#' MRC-IEU LD reference panels: `"EUR"`, `"AFR"`, `"EAS"`, or `"SAS"`. Unknown
#' values fall back to `"EUR"` with a warning so downstream LD operations can
#' still proceed.
#'
#' @param ancestry Character scalar; the ancestry label to map.
#' @return One of `"EUR"`, `"AFR"`, `"EAS"`, `"SAS"`.
#' @export
#' @examples
#' ld_pop_from_ancestry("EUR")          # "EUR"
#' ld_pop_from_ancestry("south_asian")  # "SAS"
ld_pop_from_ancestry <- function(ancestry) {
  if (is.null(ancestry) || length(ancestry) == 0L || is.na(ancestry[[1]])) {
    label <- if (is.null(ancestry)) "NULL" else as.character(ancestry[[1]])
    warning(sprintf("Unknown ancestry '%s'; defaulting LD pop to EUR for clumping.", label))
    return("EUR")
  }
  anc0 <- toupper(trimws(as.character(ancestry[[1]])))
  if (anc0 %in% c("EUR","EUROPEAN")) return("EUR")
  if (anc0 %in% c("AFR","AFRICAN","AFRICAN-UKB")) return("AFR")
  if (anc0 %in% c("EAS","EAST_ASIAN","EAST-ASIAN","EAST ASIAN")) return("EAS")
  if (anc0 %in% c("CSA","SAS","SOUTH_ASIAN","CENTRAL_SOUTH_ASIAN","SOUTH ASIAN")) return("SAS")
  warning(sprintf("Unknown ancestry '%s'; defaulting LD pop to EUR for clumping.", ancestry))
  "EUR"
}
