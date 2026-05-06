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

# Regression for the bug where exposures with NA sample_size in OpenGWAS
# (notably some pQTL records, e.g. prot-c-3722_49_2 = TNF-alpha) silently
# disabled coloc with "exposure N unknown; skipping". The fix added an
# optional `samplesize` CSV column that build_exposure_snps prefers over
# fetch_samplesize().
test_that("run_ieugwasr_ard_compare exposes samplesize override codepath", {
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  src <- paste(readLines(src_path), collapse = "\n")

  # build_exposure_snps must accept samplesize_override as a formal
  expect_match(
    src,
    "build_exposure_snps\\s*<-\\s*function\\s*\\([^)]*samplesize_override",
    info = "build_exposure_snps must accept a samplesize_override argument"
  )

  # build_exposure_snps must use samplesize_override to populate n_total
  expect_match(
    src,
    "n_override\\s*<-\\s*suppressWarnings\\(as\\.numeric\\(samplesize_override\\)\\)",
    info = "build_exposure_snps must coerce and inspect samplesize_override"
  )

  # the per-row map step must thread r$samplesize into the call
  expect_match(
    src,
    "build_exposure_snps\\([^)]*samplesize_override\\s*=\\s*r\\$samplesize",
    info = "the per-row purrr::map() must forward r$samplesize to build_exposure_snps"
  )

  # the CSV-reading step must default `samplesize` to NA_real_ if absent
  expect_match(
    src,
    'if\\s*\\(!"samplesize"\\s*%in%\\s*names\\(expos\\)\\)\\s*expos\\$samplesize\\s*<-\\s*NA_real_',
    info = "CSV must accept an optional `samplesize` column with NA fallback"
  )
})

# Regression for the bug where ard_compare() silently dropped ieu_id /
# exposure_id / exposure_sumstats from groups_info[[i]], causing every
# downstream run_phenome_mr() call to receive exposure_id = NULL and
# acknowledge_no_coloc = TRUE -> coloc disabled for the run.
test_that("ard_compare forwards ieu_id from groups to run_phenome_mr", {
  captured <- list()
  testthat::local_mocked_bindings(
    run_phenome_mr = function(exposure_snps, exposure_id = NULL,
                              acknowledge_no_coloc = FALSE, ...) {
      captured$exposure_id          <<- exposure_id
      captured$acknowledge_no_coloc <<- acknowledge_no_coloc
      list(MR_df = tibble::tibble(), results_df = tibble::tibble())
    }
  )
  exposure_snps <- tibble::tibble(
    SNP                    = c("rs1", "rs2", "rs3"),
    effect_allele.exposure = c("A", "C", "G"),
    other_allele.exposure  = c("G", "T", "A"),
    beta.exposure          = c(0.05, -0.03, 0.12),
    se.exposure            = c(0.010, 0.015, 0.008),
    pval.exposure          = c(5e-7, 4.5e-2, 1e-30),
    eaf.exposure           = c(0.35, 0.42, 0.28),
    samplesize.exposure    = c(100000L, 100000L, 100000L),
    id.exposure            = "synth_exposure",
    exposure               = "Synthetic exposure"
  )
  invisible(tryCatch(
    ard_compare(
      exposure       = "TEST",
      exposure_units = "1-SD",
      groups         = list(eur_both = list(
        sex           = "both",
        ancestry      = "EUR",
        exposure_snps = exposure_snps,
        ieu_id        = "ieu-b-110"
      )),
      cache_dir      = tempdir(),
      verbose        = FALSE,
      force_refresh  = TRUE
    ),
    error = function(e) NULL  # post-call assembly may error; we only care about the captured args
  ))
  expect_equal(captured$exposure_id, "ieu-b-110")
  expect_false(isTRUE(captured$acknowledge_no_coloc))
})

# Structural test: pre-filter and review CSV snapshots must include the
# run_hash so different parameter sets don't silently overwrite each other.
test_that("phewas/phenoscanner snapshot CSVs are suffixed with run_hash", {
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  src <- paste(readLines(src_path), collapse = "\n")

  # Old hardcoded literals must be gone.
  for (lit in c("phewas_pre_filter\\.csv",
                "phenoscanner_pre_filter\\.csv",
                "phewas_review\\.csv",
                "phenoscanner_review\\.csv")) {
    expect_false(grepl(paste0('"', lit, '"'), src),
                 info = paste0("Bare '", gsub("\\\\", "", lit),
                               "' literal still present."))
  }

  # New suffixed forms must be present.
  for (stem in c("phewas_pre_filter",
                 "phenoscanner_pre_filter",
                 "phewas_review",
                 "phenoscanner_review")) {
    pattern <- paste0('sprintf\\("', stem, '_%s\\.csv", run_hash\\)')
    expect_match(src, pattern, fixed = FALSE,
                 info = paste0("Missing sprintf form for ", stem))
  }

  # Both filter functions must declare a run_hash parameter.
  expect_match(src,
               "apply_phewas_filters\\s*<-\\s*function\\([^)]*run_hash[^)]*\\)",
               fixed = FALSE,
               info = "apply_phewas_filters() missing run_hash parameter")
  expect_match(src,
               "apply_phenoscanner_filters\\s*<-\\s*function\\([^)]*run_hash[^)]*\\)",
               fixed = FALSE,
               info = "apply_phenoscanner_filters() missing run_hash parameter")

  # The runner must compute run_hash via .run_input_hash() before the calls.
  expect_match(src, "run_hash\\s*<-\\s*\\.run_input_hash\\(", fixed = FALSE)
})
