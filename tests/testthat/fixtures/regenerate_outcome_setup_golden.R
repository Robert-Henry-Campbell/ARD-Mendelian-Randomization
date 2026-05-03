# tests/testthat/fixtures/regenerate_outcome_setup_golden.R
#
# Regenerate the golden snapshot fixtures used by
# tests/testthat/test_outcome_setup.R.
#
# Produces three fixtures (one per sex) at ancestry = "EUR":
#   * outcome_setup_golden_both.rds
#   * outcome_setup_golden_male.rds
#   * outcome_setup_golden_female.rds
#
# RUN MANUALLY. testthat does NOT auto-source files inside
# tests/testthat/fixtures/, so this script will not run during
# devtools::test() / R CMD check. Re-run only when you have INTENTIONALLY
# changed Outcome_setup() or the packaged manifest .rda files and want to
# pin the new behaviour.
#
# Usage from an R session at the project root:
#   source("tests/testthat/fixtures/regenerate_outcome_setup_golden.R")

stopifnot(
  basename(getwd()) == "ardmr" || file.exists("DESCRIPTION")
)

devtools::load_all(".", quiet = TRUE)

fixture_dir <- file.path("tests", "testthat", "fixtures")
if (!dir.exists(fixture_dir)) dir.create(fixture_dir, recursive = TRUE)

write_one <- function(sex_arg) {
  obj <- Outcome_setup(sex = sex_arg, ancestry = "EUR")
  attr(obj, "regenerated_on") <- Sys.time()
  attr(obj, "R_version")      <- R.version.string
  attr(obj, "sex_arg")        <- sex_arg
  attr(obj, "ancestry_arg")   <- "EUR"
  fname <- sprintf("outcome_setup_golden_%s.rds", sex_arg)
  fpath <- file.path(fixture_dir, fname)
  saveRDS(obj, fpath)
  message(
    "Wrote ", fpath,
    " (", format(structure(file.size(fpath), class = "object_size"),
                 units = "auto"),
    "; ", nrow(obj), " rows, ", ncol(obj), " cols)"
  )
  invisible(obj)
}

write_one("both")
write_one("male")
write_one("female")

message("R ", R.version.string)
