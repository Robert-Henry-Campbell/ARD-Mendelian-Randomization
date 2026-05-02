# tests/testthat/fixtures/regenerate_mr_business_logic_golden.R
#
# Regenerate the golden snapshot fixture used by
# tests/testthat/test_mr_business_logic.R.
#
# RUN MANUALLY. testthat does NOT auto-source files inside
# tests/testthat/fixtures/, so this script will not run during
# devtools::test() / R CMD check. Re-run only when you have INTENTIONALLY
# changed mr_business_logic()'s output (e.g. behavioural fix, TwoSampleMR
# upgrade) and want to pin the new behaviour.
#
# Usage from an R session at the project root:
#   source("tests/testthat/fixtures/regenerate_mr_business_logic_golden.R")

stopifnot(
  basename(getwd()) == "ardmr" || file.exists("DESCRIPTION")
)

devtools::load_all(".", quiet = TRUE)
source(file.path("tests", "testthat", "helper-mr_business_logic.R"))

inputs <- make_synthetic_mr_inputs()
golden <- run_synthetic_mr_business_logic(inputs)

# Record provenance so future drift is diagnosable.
attr(golden, "regenerated_on")      <- Sys.time()
attr(golden, "TwoSampleMR_version") <- as.character(utils::packageVersion("TwoSampleMR"))
attr(golden, "R_version")           <- R.version.string
attr(golden, "synth_seed")          <- SYNTH_SEED

fixture_dir  <- file.path("tests", "testthat", "fixtures")
fixture_file <- file.path(fixture_dir, "mr_business_logic_golden.rds")
if (!dir.exists(fixture_dir)) dir.create(fixture_dir, recursive = TRUE)
saveRDS(golden, fixture_file)

message("Wrote ", fixture_file,
        " (", format(structure(file.size(fixture_file), class = "object_size"),
                     units = "auto"), ")")
message("TwoSampleMR ", attr(golden, "TwoSampleMR_version"),
        " | ", attr(golden, "R_version"))
