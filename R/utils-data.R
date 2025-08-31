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
