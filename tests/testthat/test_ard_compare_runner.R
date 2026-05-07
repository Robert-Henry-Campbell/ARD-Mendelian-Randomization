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

# ---- Job 3: clump_opts plumbing ------------------------------------------

mk_curated_snps <- function() {
  tibble::tibble(
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
}

run_ardc_capture <- function(clump_opts_arg = NULL) {
  captured <- new.env()
  testthat::local_mocked_bindings(
    run_phenome_mr = function(exposure_snps, exposure_id = NULL,
                              clump_opts = list(),
                              acknowledge_no_coloc = FALSE, ...) {
      captured$clump_opts <- clump_opts
      list(MR_df = tibble::tibble(), results_df = tibble::tibble())
    },
    .env = parent.frame()
  )
  args <- list(
    exposure       = "TEST",
    exposure_units = "1-SD",
    groups         = list(eur_both = list(
      sex           = "both",
      ancestry      = "EUR",
      exposure_snps = mk_curated_snps(),
      ieu_id        = "ieu-b-110"
    )),
    cache_dir      = tempdir(),
    verbose        = FALSE,
    force_refresh  = TRUE
  )
  if (!is.null(clump_opts_arg)) args$clump_opts <- clump_opts_arg
  invisible(tryCatch(do.call(ard_compare, args), error = function(e) NULL))
  captured
}

test_that("ard_compare default clump_opts marks SNPs as already_clumped + already_p_filtered", {
  cap <- run_ardc_capture()
  expect_true(isTRUE(cap$clump_opts$already_clumped))
  expect_true(isTRUE(cap$clump_opts$already_p_filtered))
})

test_that("ard_compare's clump_opts arg is forwarded verbatim to run_phenome_mr", {
  user_opts <- list(maf_min = 0.05,
                    already_clumped = TRUE,
                    already_p_filtered = TRUE,
                    drop_palindromic = TRUE)
  cap <- run_ardc_capture(clump_opts_arg = user_opts)
  expect_equal(cap$clump_opts$maf_min, 0.05)
  expect_true(isTRUE(cap$clump_opts$already_clumped))
  expect_true(isTRUE(cap$clump_opts$already_p_filtered))
  expect_true(isTRUE(cap$clump_opts$drop_palindromic))
})

test_that("ard_compare exposes clump_opts as a formal parameter", {
  expect_true("clump_opts" %in% names(formals(ard_compare)))
  default <- eval(formals(ard_compare)$clump_opts)
  expect_true(isTRUE(default$already_clumped))
  expect_true(isTRUE(default$already_p_filtered))
})

# ---- Job 3: structural checks on ard_compare_grouped_ieugwasr.R ---------

test_that("inline literal p-value backoff ladder no longer present in ieugwasr runner", {
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  src <- paste(readLines(src_path), collapse = "\n")
  # No hardcoded c(p_threshold, 5e-7, 5e-6) literal anywhere (Job 3
  # removed the original; Job 4 does NOT silently re-introduce it).
  expect_false(grepl("c\\(p_threshold,\\s*5e-7,\\s*5e-6\\)", src),
               info = "Inline c(p_threshold, 5e-7, 5e-6) literal must be gone.")
  # Job 4: default behaviour must be a single rung c(p_threshold).
  expect_match(src, "if\\s*\\(is\\.null\\(p_backoff\\)\\)\\s*c\\(p_threshold\\)",
               info = "Runner must default to a single-rung ladder when p_backoff is NULL.")
})

# ---- Job 4: p_backoff parameter and local-clumping refactor --------------

test_that("run_ieugwasr_ard_compare has a p_backoff formal parameter (default NULL)", {
  expect_true("p_backoff" %in% names(formals(run_ieugwasr_ard_compare)))
  expect_null(eval(formals(run_ieugwasr_ard_compare)$p_backoff))
})

test_that("get_instruments uses local clumping via preprocess_exposure_snps (Job 4)", {
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  lines <- readLines(src_path)
  # Find the get_instruments closure body and assert the new local-clump call.
  start <- grep("get_instruments\\s*<-\\s*function", lines)
  expect_true(length(start) >= 1L, info = "get_instruments closure not found")
  # Walk forward until matching close brace.
  depth <- 0L; end <- NA_integer_
  for (i in seq.int(start[1], length(lines))) {
    chars <- strsplit(lines[i], "")[[1]]
    for (ch in chars) {
      if (ch == "{") depth <- depth + 1L
      if (ch == "}") {
        depth <- depth - 1L
        if (depth == 0L) { end <- i; break }
      }
    }
    if (!is.na(end)) break
  }
  expect_false(is.na(end), info = "Could not find end of get_instruments body")
  block <- paste(lines[start[1]:end], collapse = "\n")
  # Must NOT use the old server-clump default (clump = TRUE) anywhere.
  expect_false(grepl("clump\\s*=\\s*TRUE", block),
               info = "get_instruments must not server-clump (Job 4 always-local).")
  # Must call extract_instruments with clump = FALSE.
  expect_match(block, "clump\\s*=\\s*FALSE",
               info = "get_instruments must call extract_instruments(clump = FALSE).")
  # Must delegate the post-fetch work to preprocess_exposure_snps.
  expect_match(block, "preprocess_exposure_snps\\s*\\(",
               info = "get_instruments must call preprocess_exposure_snps for local clump + filters.")
  # Must mark the p-filter as already-done (server fetched the rung).
  expect_match(block, "already_p_filtered\\s*=\\s*TRUE",
               info = "get_instruments must pass already_p_filtered = TRUE.")
  # Must request local clumping (already_clumped = FALSE).
  expect_match(block, "already_clumped\\s*=\\s*FALSE",
               info = "get_instruments must pass already_clumped = FALSE for local clump.")
})

test_that("per-row clump_opts in run_ieugwasr_ard_compare includes p_backoff (Job 4)", {
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  lines <- readLines(src_path)
  # Find the .run_input_hash( call site inside the runner and walk
  # forward until the matching close paren, then assert p_backoff
  # appears inside its argument block.
  open_idx <- grep("\\.run_input_hash\\(", lines)
  expect_true(length(open_idx) >= 1L,
              info = ".run_input_hash( call site not found")
  start <- open_idx[1]
  depth <- 0L; end <- NA_integer_
  for (i in seq.int(start, length(lines))) {
    chars <- strsplit(lines[i], "")[[1]]
    for (ch in chars) {
      if (ch == "(") depth <- depth + 1L
      if (ch == ")") {
        depth <- depth - 1L
        if (depth == 0L) { end <- i; break }
      }
    }
    if (!is.na(end)) break
  }
  expect_false(is.na(end), info = "Could not find matching ) for .run_input_hash( call")
  block <- paste(lines[start:end], collapse = "\n")
  expect_match(block, "p_backoff\\s*=",
               info = "p_backoff must appear in the per-row clump_opts list passed to .run_input_hash().")
})

test_that("per-row clump_opts in run_ieugwasr_ard_compare includes f_threshold", {
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  src <- paste(readLines(src_path), collapse = "\n")
  # The .run_input_hash call's clump_opts list must include f_threshold
  # so the cache key reflects every parameter actually used downstream.
  pattern <- "clump_opts\\s*=\\s*list\\([^)]*f_threshold\\s*="
  expect_match(src, pattern, info = "f_threshold must appear in the per-row clump_opts list.")
})

test_that("run_ieugwasr_ard_compare passes already_clumped/already_p_filtered when calling ard_compare", {
  src_path <- testthat::test_path("..", "..", "R", "ard_compare_grouped_ieugwasr.R")
  skip_if_not(file.exists(src_path), "Source file not found")
  lines <- readLines(src_path)
  # Find the line that opens the ard_compare() call inside the runner.
  open_idx <- grep("compare_info\\s*<-\\s*ard_compare\\(", lines)
  expect_true(length(open_idx) >= 1L,
              info = "compare_info <- ard_compare( call site not found")
  # Walk forward until the matching close paren; everything in between is
  # the call's argument block.
  start <- open_idx[1]
  depth <- 0L
  end <- NA_integer_
  for (i in seq.int(start, length(lines))) {
    chars <- strsplit(lines[i], "")[[1]]
    for (ch in chars) {
      if (ch == "(") depth <- depth + 1L
      if (ch == ")") {
        depth <- depth - 1L
        if (depth == 0L) { end <- i; break }
      }
    }
    if (!is.na(end)) break
  }
  expect_false(is.na(end), info = "Could not find matching ) for ard_compare( call")
  block <- paste(lines[start:end], collapse = "\n")
  expect_match(block, "clump_opts\\s*=",
               info = "ard_compare() call must thread a clump_opts argument")
  expect_match(block, "already_clumped\\s*=\\s*TRUE",
               info = "ard_compare() call's clump_opts must set already_clumped = TRUE")
  expect_match(block, "already_p_filtered\\s*=\\s*TRUE",
               info = "ard_compare() call's clump_opts must set already_p_filtered = TRUE")
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
