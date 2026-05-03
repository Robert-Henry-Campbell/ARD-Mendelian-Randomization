# tests/testthat/fixtures/regenerate_enrichment_business_logic_golden.R
#
# Regenerate the golden snapshot fixture used by
# tests/testthat/test_enrichment_business_logic.R.
#
# Snapshots the output of all 7 exported enrichment_business_logic
# entry points, run on a synthetic input chained from the
# mr_business_logic golden fixture (see
# tests/testthat/helper-synthetic-data.R::make_synthetic_enrichment_input()).
#
# RUN MANUALLY. testthat does NOT auto-source files inside
# tests/testthat/fixtures/, so this script will not run during
# devtools::test() / R CMD check. Re-run only when you have INTENTIONALLY
# changed enrichment_business_logic, the upstream mr_business_logic
# golden fixture, or the synthetic-input augmentation.
#
# Usage from an R session at the project root:
#   source("tests/testthat/fixtures/regenerate_enrichment_business_logic_golden.R")

stopifnot(
  basename(getwd()) == "ardmr" || file.exists("DESCRIPTION")
)

devtools::load_all(".", quiet = TRUE)
source(file.path("tests", "testthat", "helper-synthetic-data.R"))

# Upstream dependency: the mr_business_logic golden fixture must exist,
# because make_synthetic_enrichment_input() reads its results_df.
mr_fixture <- file.path("tests", "testthat", "fixtures",
                        "mr_business_logic_golden.rds")
if (!file.exists(mr_fixture)) {
  stop("Upstream fixture missing: ", mr_fixture,
       ". Source tests/testthat/fixtures/regenerate_mr_business_logic_golden.R first.",
       call. = FALSE)
}

input  <- make_synthetic_enrichment_input()
golden <- run_synthetic_enrichment(input)

# Provenance for future drift diagnosis.
attr(golden, "regenerated_on")     <- Sys.time()
attr(golden, "R_version")          <- R.version.string
attr(golden, "synth_seed")         <- SYNTH_SEED
attr(golden, "upstream_mr_fixture")<- mr_fixture
attr(golden, "upstream_mr_R_version") <-
  attr(readRDS(mr_fixture), "R_version")

fixture_dir  <- file.path("tests", "testthat", "fixtures")
fixture_file <- file.path(fixture_dir, "enrichment_business_logic_golden.rds")
if (!dir.exists(fixture_dir)) dir.create(fixture_dir, recursive = TRUE)
saveRDS(golden, fixture_file)

message(
  "Wrote ", fixture_file,
  " (", format(structure(file.size(fixture_file), class = "object_size"),
               units = "auto"), ")"
)
message("R ", R.version.string, " | seed ", SYNTH_SEED)
