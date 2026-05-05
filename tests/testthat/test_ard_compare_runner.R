# tests/testthat/test_ard_compare_runner.R
#
# Structural test: the ieugwasr runner must enable "coloc" when it
# constructs the sensitivity_enabled list it passes to ard_compare().
# This guards against regression of the silent-skip bug where the runner
# hardcoded an 8-element list that omitted "coloc".

test_that("run_ieugwasr_ard_compare passes 'coloc' in sensitivity_enabled", {
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  src <- paste(readLines(src_path), collapse = "\n")
  block <- regmatches(
    src,
    regexpr("sensitivity_enabled\\s*=\\s*c\\([^)]+\\)", src)
  )
  expect_true(length(block) == 1L && nzchar(block),
              info = "sensitivity_enabled = c(...) literal not found")
  expect_match(block, "\"coloc\"",
               info = paste("Found:", block))
})

test_that("ieugwasr runner's sensitivity list is a superset of run_phenome_mr default", {
  default <- eval(formals(run_phenome_mr)$sensitivity_enabled)
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  src <- paste(readLines(src_path), collapse = "\n")
  block <- regmatches(
    src,
    regexpr("sensitivity_enabled\\s*=\\s*c\\([^)]+\\)", src)
  )
  vals <- regmatches(block, gregexpr('"[^"]+"', block))[[1]]
  vals <- gsub('"', '', vals)
  expect_true(all(default %in% vals),
              info = paste("Missing:", paste(setdiff(default, vals), collapse = ", ")))
})
