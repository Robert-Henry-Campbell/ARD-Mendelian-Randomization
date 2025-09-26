#' Install package dependencies
#'
#' Ensures the required packages for this project are installed.
#' This includes core MR libraries such as `MendelianRandomization`
#' and support tooling like `jsonlite` as declared in `DESCRIPTION`.
#' If the `remotes` package is missing, it will be installed first.
#'
#' @keywords internal
install_deps <- function() {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_deps()
  invisible(NULL)
}
