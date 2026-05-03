# tests/testthat/test_enrichment_business_logic.R
#
# Golden-snapshot regression test for the seven exported entry points of
# enrichment_business_logic. Synthetic input is chained from the
# mr_business_logic golden fixture (its results_df, augmented with
# ARD_selected + cause_level_* + pheno_sex columns).
#
# With 2 valid rows after qc filtering (1 ARD + 1 non-ARD in cause "A"),
# the permutation tests land on the exact-enumeration path
# (choose(2, 1) = 2 permutations), so output is bit-exact reproducible
# without RNG concerns.
#
# If a snapshot fails after an INTENTIONAL change to enrichment_business_logic
# or to the upstream mr_business_logic golden fixture, regenerate by
# sourcing tests/testthat/fixtures/regenerate_enrichment_business_logic_golden.R
# and commit the new .rds.

skip_if_not_installed("dplyr")

fixture_path <- testthat::test_path("fixtures",
                                    "enrichment_business_logic_golden.rds")

test_that("golden fixture exists (run tests/testthat/fixtures/regenerate_enrichment_business_logic_golden.R if missing)", {
  expect_true(
    file.exists(fixture_path),
    info = paste0("Missing fixture: ", fixture_path)
  )
})

if (!file.exists(fixture_path)) return(invisible())

golden <- readRDS(fixture_path)
out    <- run_synthetic_enrichment()

# ---- Top-level structure -----------------------------------------------------

test_that("run_synthetic_enrichment() returns the expected list of seven entries", {
  expect_named(out, c(
    "run_enrichment", "enrichment_global_directional",
    "enrichment_by_cause_directional",
    "beta_contrast_by_cause", "beta_contrast_global_ARD",
    "beta_mean_by_cause", "beta_mean_global"
  ))
  expect_named(out, names(golden))
})

# ---- run_enrichment orchestrator --------------------------------------------

test_that("run_enrichment matches the golden fixture", {
  expect_named(out$run_enrichment, c("global_tbl", "by_cause_tbl"))
  expect_equal(out$run_enrichment$global_tbl,
               golden$run_enrichment$global_tbl,
               tolerance = 1e-6)
  expect_equal(out$run_enrichment$by_cause_tbl,
               golden$run_enrichment$by_cause_tbl,
               tolerance = 1e-6)
})

# ---- Direct calls per entry point -------------------------------------------

test_that("enrichment_global_directional matches the golden fixture", {
  expect_equal(out$enrichment_global_directional,
               golden$enrichment_global_directional,
               tolerance = 1e-6)
})

test_that("enrichment_by_cause_directional matches the golden fixture across all 3 modes", {
  expect_named(out$enrichment_by_cause_directional, names(golden$enrichment_by_cause_directional))
  for (md in names(golden$enrichment_by_cause_directional)) {
    expect_equal(out$enrichment_by_cause_directional[[md]],
                 golden$enrichment_by_cause_directional[[md]],
                 tolerance = 1e-6,
                 info = paste("compare_mode =", md))
  }
})

test_that("beta_contrast_by_cause matches the golden fixture (empty-result shape)", {
  # With only one cause in the synthetic input, this function legitimately
  # returns an empty tibble with the documented column structure. The
  # snapshot pins both the empty-ness AND the column shape.
  expect_equal(out$beta_contrast_by_cause,
               golden$beta_contrast_by_cause,
               tolerance = 1e-6)
  expect_equal(nrow(out$beta_contrast_by_cause), 0L)
})

test_that("beta_contrast_global_ARD matches the golden fixture", {
  expect_equal(out$beta_contrast_global_ARD,
               golden$beta_contrast_global_ARD,
               tolerance = 1e-6)
})

test_that("beta_mean_by_cause matches the golden fixture", {
  expect_equal(out$beta_mean_by_cause,
               golden$beta_mean_by_cause,
               tolerance = 1e-6)
})

test_that("beta_mean_global matches the golden fixture", {
  expect_equal(out$beta_mean_global,
               golden$beta_mean_global,
               tolerance = 1e-6)
})

# ---- Determinism --------------------------------------------------------------

test_that("rerunning produces bit-identical output (exact-enumeration path)", {
  out2 <- run_synthetic_enrichment()
  expect_identical(out$run_enrichment$global_tbl,
                   out2$run_enrichment$global_tbl)
  expect_identical(out$enrichment_global_directional,
                   out2$enrichment_global_directional)
})
