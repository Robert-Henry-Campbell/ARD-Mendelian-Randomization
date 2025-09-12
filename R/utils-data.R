# R/utils-data.R

#' Internal: fetch a dataset shipped in the package
#'
#' @param name String, name of the dataset (without quotes).
#' @param pkg Package name (default = "ardmr").
#'
#' @return A tibble with the dataset contents.
#' @keywords internal
get_pkg_obj <- function(name, pkg = "ardmr") {
  e <- new.env(parent = emptyenv())
  utils::data(list = name, package = pkg, envir = e)
  if (!exists(name, envir = e, inherits = FALSE)) {
    stop(sprintf("Data object '%s' not found in package '%s'", name, pkg), call. = FALSE)
  }
  get(name, envir = e, inherits = FALSE) |> dplyr::as_tibble()
}

# small internal assert (keeps orchestrator readable)
assert_exposure <- function(x) {
  req <- c("id.exposure", "SNP","beta.exposure","se.exposure","effect_allele.exposure","other_allele.exposure","pval.exposure")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("exposure_snps missing columns: ", paste(miss, collapse=", "))
  invisible(TRUE)
}


fix_exposure_names <- function(df) {
  stopifnot(is.data.frame(df))
  nms  <- names(df)
  nmsL <- tolower(nms)

  # keep these unsuffixed (everything else gets '.exposure' if missing)
  keep <- c("snp","exposure")

  nms_fixed <- ifelse(
    nmsL %in% keep,
    nmsL,
    ifelse(grepl("\\.exposure$", nmsL), nmsL, paste0(nmsL, ".exposure"))
  )

  # guard against accidental "exposure.exposure"
  nms_fixed <- sub("^exposure\\.exposure$", "exposure", nms_fixed)

  # de-dup just in case
  nms_fixed <- make.unique(nms_fixed, sep = "_")

  # enforce uppercase for the unsuffixed SNP column
  nms_fixed[nms_fixed == "snp"] <- "SNP"

  names(df) <- nms_fixed
  tibble::as_tibble(df)
}


coerce_lowercase_colnames <- function(df) {
  names(df) <- tolower(names(df))
  return(df)
}

