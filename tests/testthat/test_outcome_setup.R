# tests/testthat/test_outcome_setup.R
#
# Golden-snapshot regression test for Outcome_setup() — maps ARD ICD10
# codes to GWAS manifests using packaged .rda data only. Fully offline.
#
# Output is fully deterministic (pure joins + filters, no RNG), so
# tolerance can be tight.
#
# If a snapshot fails after an INTENTIONAL change to Outcome_setup() or
# to the packaged manifests, regenerate by sourcing
#   tests/testthat/fixtures/regenerate_outcome_setup_golden.R
# and commit the new .rds files.

skip_if_not_installed("logger")

fixture_path <- function(sex_arg) {
  testthat::test_path("fixtures",
                      sprintf("outcome_setup_golden_%s.rds", sex_arg))
}

# ---- Fixture existence ------------------------------------------------------

test_that("golden fixtures exist (run tests/testthat/fixtures/regenerate_outcome_setup_golden.R if missing)", {
  for (sx in c("both", "male", "female")) {
    expect_true(file.exists(fixture_path(sx)),
                info = paste("Missing fixture for sex =", sx))
  }
})

if (!all(vapply(c("both", "male", "female"),
                function(sx) file.exists(fixture_path(sx)),
                logical(1)))) {
  return(invisible())
}

# ---- Snapshot equality per sex branch ---------------------------------------

run_outcome_setup <- function(sex_arg) {
  Outcome_setup(sex = sex_arg, ancestry = "EUR")
}

for (sx in c("both", "male", "female")) {
  local({
    sex_arg <- sx
    test_that(paste0("Outcome_setup('", sex_arg, "', 'EUR') matches the golden fixture"),
              {
                golden <- readRDS(fixture_path(sex_arg))
                out    <- run_outcome_setup(sex_arg)

                # Drop the regenerator's provenance attributes before comparing,
                # so a fresh run (which has none) is comparable to the .rds.
                attr(golden, "regenerated_on") <- NULL
                attr(golden, "R_version")      <- NULL
                attr(golden, "sex_arg")        <- NULL
                attr(golden, "ancestry_arg")   <- NULL

                expect_equal(nrow(out), nrow(golden))
                expect_equal(ncol(out), ncol(golden))
                expect_named(out, names(golden))
                expect_equal(out, golden, tolerance = 1e-8)
              })
  })
}

# ---- Validation paths -------------------------------------------------------

test_that("Outcome_setup() rejects sex-specific runs at non-EUR ancestry", {
  expect_error(Outcome_setup("male", "AFR"), "EUR")
  expect_error(Outcome_setup("female", "AFR"), "EUR")
})

test_that("Outcome_setup() rejects unknown sex / ancestry args", {
  expect_error(Outcome_setup("xyz", "EUR"))      # match.arg on sex
  expect_error(Outcome_setup("both", "ZZZ"))     # match.arg on ancestry
})

# ---- Structural invariants --------------------------------------------------

test_that("Outcome_setup() output has logical ARD_selected and unique ICD10_explo", {
  for (sx in c("both", "male", "female")) {
    out <- run_outcome_setup(sx)
    expect_type(out$ARD_selected, "logical")
    expect_equal(length(unique(out$ICD10_explo)), nrow(out),
                 info = paste("duplicate ICD10_explo for sex =", sx))
  }
})

test_that("Outcome_setup() routes to Pan-UKB for sex='both' and Neale for sex-specific", {
  # Distinguishing columns introduced by each manifest's join (the join key
  # column itself is consumed by left_join's `by =` mapping):
  #   Pan-UKB manifest contributes e.g. `aws_link`, `aws_link_tabix`.
  #   Neale manifest contributes e.g. `Phenotype Description`, `AWS File`.
  out_both   <- run_outcome_setup("both")
  out_male   <- run_outcome_setup("male")
  out_female <- run_outcome_setup("female")

  expect_true("aws_link" %in% names(out_both),
              info = "sex='both' should route through Pan-UKB manifest")
  expect_false("Phenotype Description" %in% names(out_both),
               info = "sex='both' should NOT carry Neale manifest columns")

  expect_true("Phenotype Description" %in% names(out_male),
              info = "sex='male' should route through Neale manifest")
  expect_false("aws_link" %in% names(out_male),
               info = "sex='male' should NOT carry Pan-UKB manifest columns")

  expect_true("Phenotype Description" %in% names(out_female),
              info = "sex='female' should route through Neale manifest")
  expect_false("aws_link" %in% names(out_female),
               info = "sex='female' should NOT carry Pan-UKB manifest columns")
})
