# tests/testthat/test_preprocess_exposure_snps.R
# Per-step unit tests for the unified exposure-SNP preprocessing pipeline.
# Each test isolates one step by configuring `clump_opts` so the other
# steps are no-ops (already_p_filtered/already_clumped, default-OFF
# thresholds, etc.).

# ---- shared helpers -------------------------------------------------------

# clump_opts that turns every step into a no-op except the one under test.
# Defaults: skip p-backoff, skip clump, F-stat unreachable, MAF/INFO off,
# indel/palindromic at their pipeline defaults.
mk_skip_opts <- function(...) {
  utils::modifyList(
    list(
      already_p_filtered = TRUE,
      already_clumped    = TRUE,
      f_threshold        = -Inf,
      drop_indels        = FALSE,
      drop_palindromic   = FALSE,
      maf_min            = NULL,
      info_min           = NULL
    ),
    list(...)
  )
}

# Capture logger output for a single call. Restores the appender on exit.
.simple_console_appender <- function(lines) cat(lines, sep = "\n", file = stderr())
with_captured_logs <- function(expr) {
  captured <- character()
  logger::log_appender(function(lines) captured <<- c(captured, lines))
  on.exit(logger::log_appender(.simple_console_appender), add = TRUE)
  val <- force(expr)
  list(value = val, logs = captured)
}

# Find a recorded step by name (returns NULL if not found).
get_step <- function(res, step_name) {
  hit <- Filter(function(x) identical(x$step, step_name), res$steps)
  if (!length(hit)) NULL else hit[[1]]
}

# ---- p-value backoff ------------------------------------------------------

test_that("p-backoff: all rows pass at strictest rung -> value = 5e-8", {
  snps <- make_synthetic_exposure_snps(
    n = 3L, pval.exposure = c(1e-10, 1e-9, 1e-12)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(already_p_filtered = FALSE),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "p_threshold")
  expect_equal(step$value, 5e-8)
  expect_equal(step$n_in, 3L)
  expect_equal(step$n_out, 3L)
  expect_equal(nrow(res$snps), 3L)
})

test_that("p-backoff: strictest empty -> walks to next rung", {
  snps <- make_synthetic_exposure_snps(
    n = 3L, pval.exposure = c(1e-7, 5e-7, 9e-7)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(already_p_filtered = FALSE),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "p_threshold")
  expect_equal(step$value, 5e-7)
  expect_equal(step$n_in, 3L)
  expect_equal(step$n_out, 1L)  # only 1e-7 passes p < 5e-7
})

test_that("p-backoff: all rungs empty -> 0-row tibble, no error", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, pval.exposure = c(1e-3, 1e-2)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(already_p_filtered = FALSE),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "p_threshold")
  expect_equal(step$value, 5e-6)  # last rung tried
  expect_equal(step$n_out, 0L)
  expect_equal(nrow(res$snps), 0L)
})

test_that("p-backoff: already_p_filtered = TRUE skips the step entirely", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, pval.exposure = c(1e-3, 1e-2)  # would all be dropped
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(already_p_filtered = TRUE),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "p_threshold")
  expect_equal(step$skipped_reason, "already_p_filtered")
  expect_equal(step$n_out, 2L)  # nothing dropped
  expect_equal(nrow(res$snps), 2L)
})

test_that("p-backoff: single-element ladder does not cascade", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, pval.exposure = c(1e-7, 1e-6)
  )
  res <- preprocess_exposure_snps(
    snps,
    clump_opts = mk_skip_opts(already_p_filtered = FALSE,
                              p_backoff = c(5e-8)),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "p_threshold")
  expect_equal(step$value, 5e-8)
  expect_equal(step$n_out, 0L)
})

# ---- clump ----------------------------------------------------------------

test_that("clump: already_clumped = TRUE -> ld_clump never called", {
  snps <- make_synthetic_exposure_snps(n = 2L)
  call_count <- 0L
  testthat::local_mocked_bindings(
    ld_clump = function(...) { call_count <<- call_count + 1L; tibble::tibble() },
    .package = "ieugwasr"
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(already_clumped = TRUE),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(call_count, 0L)
  step <- get_step(res, "clump")
  expect_equal(step$skipped_reason, "already_clumped")
  expect_equal(step$n_out, 2L)
})

test_that("clump: already_clumped = FALSE + LD ref present -> ld_clump called once", {
  snps <- make_synthetic_exposure_snps(n = 3L)
  call_count <- 0L
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) {
      list(bfile_prefix = file.path(tempdir(), "fake_bfile"),
           present = TRUE, downloaded = FALSE, aborted = FALSE)
    }
  )
  testthat::local_mocked_bindings(
    .find_plink_bin = function(plink_bin = NULL) NULL
  )
  testthat::local_mocked_bindings(
    ld_clump = function(dat, ...) {
      call_count <<- call_count + 1L
      # Keep first 2 SNPs to mimic clumping
      dat[seq_len(min(2L, nrow(dat))), , drop = FALSE]
    },
    .package = "ieugwasr"
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(already_clumped = FALSE),
    ancestry = "EUR", cache_dir = tempdir(), verbose = FALSE
  )
  expect_equal(call_count, 1L)
  step <- get_step(res, "clump")
  expect_equal(step$n_out, 2L)
  expect_equal(nrow(res$snps), 2L)
})

test_that("clump: confirm = 'no' + LD reference absent -> aborts cleanly", {
  snps <- make_synthetic_exposure_snps(n = 1L)
  testthat::local_mocked_bindings(
    ld_reference_checker = function(pop, ld_ref_dir, verbose, confirm) {
      if (identical(confirm, "no")) {
        stop("ld_reference_checker(): panel missing and confirm='no' -- aborting.",
             call. = FALSE)
      }
      list(bfile_prefix = "fake", present = TRUE)
    }
  )
  expect_error(
    preprocess_exposure_snps(
      snps, clump_opts = mk_skip_opts(already_clumped = FALSE),
      ancestry = "EUR", cache_dir = tempdir(),
      confirm = "no", verbose = FALSE
    ),
    "panel missing"
  )
})

# ---- rsid validate + dedupe -----------------------------------------------

test_that("rsid_validate: non-rsid IDs dropped + warning", {
  snps <- make_synthetic_exposure_snps(
    n = 3L, SNP = c("rs1", "chr1:12345", "rs3")
  )
  expect_warning(
    res <- preprocess_exposure_snps(
      snps, clump_opts = mk_skip_opts(),
      ancestry = "EUR", verbose = FALSE
    ),
    NA  # logger::log_warn does not raise an R warning, just emits a log line
  )
  step <- get_step(res, "rsid_validate")
  expect_equal(step$n_invalid, 1L)
  expect_equal(step$n_out, 2L)
  expect_setequal(res$snps$SNP, c("rs1", "rs3"))
})

test_that("rsid_validate: duplicate rsids collapsed to first occurrence", {
  snps <- make_synthetic_exposure_snps(
    n = 4L, SNP = c("rs1", "rs2", "rs1", "rs3"),
    beta.exposure = c(0.10, 0.20, 0.99, 0.30)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "rsid_validate")
  expect_equal(step$n_duplicates, 1L)
  expect_equal(step$n_out, 3L)
  # First occurrence kept: rs1 -> beta 0.10 (not 0.99)
  rs1_row <- res$snps[res$snps$SNP == "rs1", ]
  expect_equal(rs1_row$beta.exposure, 0.10)
})

# ---- drop_indels (B3) -----------------------------------------------------

test_that("drop_indels: AT/-/D/I dropped, SNV kept", {
  snps <- make_synthetic_exposure_snps(
    n = 5L,
    SNP = paste0("rs", 1:5),
    effect_allele.exposure = c("AT", "-", "D", "I", "A"),
    other_allele.exposure  = c("A",  "A", "A", "A", "G")
  )
  res <- preprocess_exposure_snps(
    snps,
    clump_opts = mk_skip_opts(drop_indels = TRUE),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "drop_indels")
  expect_equal(step$n_in, 5L)
  expect_equal(step$n_out, 1L)
  expect_equal(res$snps$SNP, "rs5")
})

test_that("drop_indels: indel in other_allele.exposure also caught", {
  snps <- make_synthetic_exposure_snps(
    n = 2L,
    SNP = c("rs1", "rs2"),
    effect_allele.exposure = c("A", "A"),
    other_allele.exposure  = c("AT", "G")
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(drop_indels = TRUE),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(get_step(res, "drop_indels")$n_out, 1L)
  expect_equal(res$snps$SNP, "rs2")
})

test_that("drop_indels = FALSE leaves indels in place", {
  snps <- make_synthetic_exposure_snps(
    n = 2L,
    SNP = c("rs1", "rs2"),
    effect_allele.exposure = c("AT", "A"),
    other_allele.exposure  = c("A",  "G")
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(drop_indels = FALSE),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "drop_indels")
  expect_equal(step$skipped_reason, "drop_indels = FALSE")
  expect_equal(step$n_out, 2L)
})

# ---- palindromic ----------------------------------------------------------

test_that("palindromic: A/T pair flagged TRUE; non-pal flagged FALSE", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, SNP = c("rs1", "rs2"),
    effect_allele.exposure = c("A", "A"),
    other_allele.exposure  = c("T", "G")
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(drop_palindromic = FALSE),
    ancestry = "EUR", verbose = FALSE
  )
  expect_true("palindromic" %in% names(res$snps))
  expect_equal(res$snps$palindromic, c(TRUE, FALSE))
  expect_equal(nrow(res$snps), 2L)  # not dropped
  step <- get_step(res, "palindromic")
  expect_equal(step$n_flagged, 1L)
  expect_false(step$dropped)
})

test_that("palindromic: drop_palindromic = TRUE removes flagged rows", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, SNP = c("rs1", "rs2"),
    effect_allele.exposure = c("C", "A"),
    other_allele.exposure  = c("G", "G")
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(drop_palindromic = TRUE),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "palindromic")
  expect_true(step$dropped)
  expect_equal(step$n_in, 2L)
  expect_equal(step$n_out, 1L)
  expect_equal(res$snps$SNP, "rs2")
})

# ---- F-stat ---------------------------------------------------------------

test_that("f_stat: above threshold kept, below dropped", {
  snps <- make_synthetic_exposure_snps(
    n = 3L,
    SNP = c("rs1", "rs2", "rs3"),
    beta.exposure = c(0.05, 0.005, 0.10),  # F = 25, 0.25, 100
    se.exposure   = c(0.01, 0.01,  0.01)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(f_threshold = 10),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "f_stat")
  expect_equal(step$value, 10)
  expect_equal(step$n_out, 2L)
  expect_setequal(res$snps$SNP, c("rs1", "rs3"))
})

test_that("f_stat: NA/Inf dropped", {
  snps <- make_synthetic_exposure_snps(
    n = 3L,
    SNP = c("rs1", "rs2", "rs3"),
    beta.exposure = c(0.10, NA_real_, 0.10),
    se.exposure   = c(0.01, 0.01,     0)  # division by zero -> Inf
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(f_threshold = 10),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(get_step(res, "f_stat")$n_out, 1L)
  expect_equal(res$snps$SNP, "rs1")
})

# ---- MAF (B1) -------------------------------------------------------------

test_that("MAF: eaf = 0.95 with maf_min = 0.05 is KEPT (pmin semantics)", {
  snps <- make_synthetic_exposure_snps(n = 1L, eaf.exposure = 0.95)
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(maf_min = 0.05),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(get_step(res, "maf")$n_out, 1L)
  expect_equal(nrow(res$snps), 1L)
})

test_that("MAF: eaf = 0.01 with maf_min = 0.05 is DROPPED", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, SNP = c("rs1","rs2"),
    eaf.exposure = c(0.01, 0.30)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(maf_min = 0.05),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(get_step(res, "maf")$n_out, 1L)
  expect_equal(res$snps$SNP, "rs2")
})

test_that("MAF: NA EAF rows are KEPT", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, SNP = c("rs1","rs2"),
    eaf.exposure = c(NA_real_, 0.30)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(maf_min = 0.05),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(get_step(res, "maf")$n_out, 2L)
})

test_that("MAF: maf_min = NULL skips the step", {
  snps <- make_synthetic_exposure_snps(n = 1L, eaf.exposure = 0.001)
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(maf_min = NULL),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "maf")
  expect_equal(step$skipped_reason, "maf_min not set")
  expect_equal(step$n_out, 1L)
})

# ---- INFO -----------------------------------------------------------------

test_that("INFO: column 'INFO' present, value below threshold -> dropped", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, SNP = c("rs1","rs2"), INFO = c(0.5, 0.95)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(info_min = 0.8),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "info")
  expect_equal(step$column, "INFO")
  expect_equal(step$n_out, 1L)
  expect_equal(res$snps$SNP, "rs2")
})

test_that("INFO: alternate column name 'Rsq' is detected", {
  snps <- make_synthetic_exposure_snps(n = 2L, SNP = c("rs1","rs2"))
  snps$Rsq <- c(0.5, 0.95)
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(info_min = 0.8),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "info")
  expect_equal(step$column, "Rsq")
  expect_equal(step$n_out, 1L)
  expect_equal(res$snps$SNP, "rs2")
})

test_that("INFO: 'info_score' column NOT in detect list -> step skipped", {
  snps <- make_synthetic_exposure_snps(n = 2L, SNP = c("rs1","rs2"))
  snps$info_score <- c(0.1, 0.9)  # would be dropped if detected
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(info_min = 0.8),
    ancestry = "EUR", verbose = FALSE
  )
  step <- get_step(res, "info")
  expect_equal(step$skipped_reason, "no INFO column")
  expect_equal(step$n_out, 2L)
})

test_that("INFO: column missing + info_min set -> warning logged, step skipped", {
  snps <- make_synthetic_exposure_snps(n = 2L)
  cap <- with_captured_logs(
    preprocess_exposure_snps(
      snps, clump_opts = mk_skip_opts(info_min = 0.8),
      ancestry = "EUR", verbose = FALSE
    )
  )
  expect_true(any(grepl("INFO.*skipped", cap$logs)))
  step <- get_step(cap$value, "info")
  expect_equal(step$skipped_reason, "no INFO column")
  expect_equal(step$n_out, 2L)
})

test_that("INFO: NA value + info_min set -> row KEPT", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, SNP = c("rs1","rs2"), INFO = c(NA_real_, 0.95)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(info_min = 0.8),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(get_step(res, "info")$n_out, 2L)
})

# ---- order regression ------------------------------------------------------

test_that("order regression: indel runs AFTER clump (locked order)", {
  # An indel SNP would be "preserved" if it survived clumping, then dropped
  # by the indel filter. With the locked order (clump -> indel), the indel
  # SNP is present until step 4, then dropped. If the order were reversed
  # (indel -> clump), the indel would be dropped earlier and the recorded
  # `n_in` for clump would be smaller. We assert step records here.
  snps <- make_synthetic_exposure_snps(
    n = 3L,
    SNP = c("rs1", "rs2", "rs3"),
    effect_allele.exposure = c("AT", "A", "A"),  # rs1 is an indel
    other_allele.exposure  = c("A",  "G", "C")
  )
  res <- preprocess_exposure_snps(
    snps,
    clump_opts = mk_skip_opts(already_clumped = TRUE, drop_indels = TRUE),
    ancestry = "EUR", verbose = FALSE
  )
  rsid_step  <- get_step(res, "rsid_validate")
  indel_step <- get_step(res, "drop_indels")
  expect_equal(rsid_step$n_in,  3L)
  expect_equal(rsid_step$n_out, 3L)
  expect_equal(indel_step$n_in,  3L)  # all 3 reach the indel step
  expect_equal(indel_step$n_out, 2L)  # rs1 dropped
})

# ---- p_threshold deprecation ----------------------------------------------

test_that("deprecation: p_threshold translated to p_backoff = c(p_threshold)", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, pval.exposure = c(1e-8, 1e-6)
  )
  cap <- with_captured_logs(
    preprocess_exposure_snps(
      snps,
      clump_opts = list(p_threshold = 1e-7,
                        already_clumped = TRUE,
                        f_threshold = -Inf,
                        drop_indels = FALSE),
      ancestry = "EUR", verbose = FALSE
    )
  )
  res <- cap$value
  expect_true(any(grepl("p_threshold is deprecated", cap$logs)))
  expect_equal(res$resolved_opts$p_backoff, 1e-7)
  expect_null(res$resolved_opts$p_threshold)
  step <- get_step(res, "p_threshold")
  expect_equal(step$value, 1e-7)
  expect_equal(step$n_out, 1L)  # only 1e-8 < 1e-7
})

test_that("deprecation: both p_threshold AND p_backoff -> p_backoff wins, warn", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, pval.exposure = c(1e-9, 1e-6)
  )
  cap <- with_captured_logs(
    preprocess_exposure_snps(
      snps,
      clump_opts = list(p_threshold = 1e-7, p_backoff = c(5e-8),
                        already_clumped = TRUE,
                        f_threshold = -Inf,
                        drop_indels = FALSE),
      ancestry = "EUR", verbose = FALSE
    )
  )
  res <- cap$value
  expect_true(any(grepl("both 'p_threshold'", cap$logs)))
  expect_equal(res$resolved_opts$p_backoff, 5e-8)
  step <- get_step(res, "p_threshold")
  expect_equal(step$value, 5e-8)
})

# ---- manifest schema -------------------------------------------------------

test_that("manifest schema: 8 step records, each has step+n_in+n_out", {
  snps <- make_synthetic_exposure_snps(n = 1L)
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(),
    ancestry = "EUR", verbose = FALSE
  )
  expect_length(res$steps, 8L)
  expected_steps <- c("p_threshold", "clump", "rsid_validate", "drop_indels",
                      "palindromic", "f_stat", "maf", "info")
  observed_steps <- vapply(res$steps, function(x) x$step, character(1))
  expect_equal(observed_steps, expected_steps)
  for (s in res$steps) {
    expect_true(all(c("step", "n_in", "n_out") %in% names(s)))
    expect_type(s$n_in,  "integer")
    expect_type(s$n_out, "integer")
  }
})

test_that("manifest schema: resolved_opts is a superset of user clump_opts", {
  snps <- make_synthetic_exposure_snps(n = 1L)
  user_opts <- list(maf_min = 0.05, already_clumped = TRUE,
                    already_p_filtered = TRUE, drop_indels = FALSE,
                    f_threshold = -Inf)
  res <- preprocess_exposure_snps(
    snps, clump_opts = user_opts,
    ancestry = "EUR", verbose = FALSE
  )
  for (k in names(user_opts)) {
    expect_equal(res$resolved_opts[[k]], user_opts[[k]],
                 info = sprintf("user-supplied key %s not preserved", k))
  }
  # And defaults present
  expect_equal(res$resolved_opts$r2, 0.001)
  expect_equal(res$resolved_opts$kb, 10000L)
})

# ---- 0-row return ----------------------------------------------------------

test_that("0-row input returns 0-row tibble (does not error)", {
  snps <- make_synthetic_exposure_snps(n = 1L)[0, ]
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(nrow(res$snps), 0L)
  expect_length(res$steps, 8L)
})

test_that("filters that strip everything return 0-row tibble (does not error)", {
  snps <- make_synthetic_exposure_snps(
    n = 2L, SNP = c("rs1","rs2"),
    eaf.exposure = c(0.001, 0.002)
  )
  res <- preprocess_exposure_snps(
    snps, clump_opts = mk_skip_opts(maf_min = 0.05),
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(nrow(res$snps), 0L)
  expect_equal(get_step(res, "maf")$n_out, 0L)
})

# ---- input validation ------------------------------------------------------

test_that("snps_tbl missing 'SNP' column -> hard error", {
  bad <- tibble::tibble(rsid = "rs1", pval.exposure = 1e-9)
  expect_error(
    preprocess_exposure_snps(bad, ancestry = "EUR", verbose = FALSE),
    "must have a 'SNP' column"
  )
})

test_that("ancestry required", {
  snps <- make_synthetic_exposure_snps(n = 1L)
  expect_error(
    preprocess_exposure_snps(snps, ancestry = NULL, verbose = FALSE),
    "ancestry"
  )
})
