# tests/testthat/test_run_hash.R
# Verifies run_phenome_mr writes outputs into a per-run hash subfolder,
# emits run_manifest.json with the new clump_opts_resolved /
# preprocessing_steps fields (Job 2c), and errors helpfully when
# preprocessing strips every IV. Heavy stages are mocked inline so the
# test stays offline.

skip_if_not_installed("jsonlite")

mk_synth_exposure_snps <- function(beta = 0.05, pos = 100L) {
  tibble::tibble(
    SNP = "rs1", id.exposure = "synth", exposure = "synth",
    Chr = 1L, Pos = pos,
    effect_allele.exposure = "A", other_allele.exposure = "G",
    beta.exposure = beta, se.exposure = 0.01,
    pval.exposure = 5e-7, eaf.exposure = 0.3,
    samplesize.exposure = 100000L
  )
}

mk_args <- function() list(
  exposure       = "synth_exposure",
  exposure_units = "log(OR)",
  ancestry       = "EUR",
  sex            = "both",
  exposure_id    = "ieu-a-test",
  sensitivity_enabled = c("egger_intercept","weighted_median","ivw_Q","ivw_I2"),
  sensitivity_pass_min = 4L,
  scatterplot = FALSE, snpforestplot = FALSE, leaveoneoutplot = FALSE,
  confirm = "no", verbose = FALSE
)

# Default preprocess stub: pass-through (preserves the input as the IV set).
.default_preprocess_stub <- function(snps_tbl, clump_opts = list(),
                                     ancestry, cache_dir, confirm = "ask",
                                     verbose = TRUE) {
  # Job 4 locked order: p -> rsid -> indel -> MAF -> INFO -> palindromic
  # -> clump -> F-stat. The default stub passes everything through.
  list(
    snps          = tibble::as_tibble(snps_tbl),
    steps         = list(
      list(step = "p_threshold",  n_in = NROW(snps_tbl), n_out = NROW(snps_tbl), value = 5e-8),
      list(step = "rsid_validate",n_in = NROW(snps_tbl), n_out = NROW(snps_tbl)),
      list(step = "drop_indels",  n_in = NROW(snps_tbl), n_out = NROW(snps_tbl)),
      list(step = "maf",          n_in = NROW(snps_tbl), n_out = NROW(snps_tbl), value = NA_real_, skipped_reason = "maf_min not set"),
      list(step = "info",         n_in = NROW(snps_tbl), n_out = NROW(snps_tbl), value = NA_real_, skipped_reason = "info_min not set"),
      list(step = "palindromic",  n_in = NROW(snps_tbl), n_out = NROW(snps_tbl), dropped = FALSE, skipped_reason = "drop_palindromic = FALSE"),
      list(step = "clump",        n_in = NROW(snps_tbl), n_out = NROW(snps_tbl), r2 = 0.001, kb = 10000L),
      list(step = "f_stat",       n_in = NROW(snps_tbl), n_out = NROW(snps_tbl), value = 10)
    ),
    resolved_opts = ardmr:::.resolve_clump_opts(clump_opts)
  )
}

run_with_mocks <- function(args, snps = NULL, preprocess_stub = NULL) {
  cache_dir <- tempfile("ardmr_run_")
  dir.create(cache_dir, recursive = TRUE)
  args$cache_dir <- cache_dir
  args$exposure_snps <- if (is.null(snps)) mk_synth_exposure_snps() else snps
  pp <- if (is.null(preprocess_stub)) .default_preprocess_stub else preprocess_stub
  testthat::local_mocked_bindings(
    Outcome_setup = function(sex, ancestry) {
      tibble::tibble(description = "synth_outcome", ARD_selected = FALSE,
                     aws_link = "https://invalid.example/x.tsv.bgz",
                     aws_link_tabix = "https://invalid.example/x.tsv.bgz.tbi")
    },
    Variant_manifest_downloader = function(...) invisible(NULL),
    exposure_snp_mapper = function(exposure_snps, ...) {
      exposure_snps$panukb_chrom <- 1L; exposure_snps$panukb_pos <- 100L
      exposure_snps$rsid <- exposure_snps$SNP; exposure_snps
    },
    compute_ld_matrix = function(...) list(matrix = matrix(1.0, 1L, 1L), snps_used = "rs1"),
    panukb_snp_grabber = function(exposure_snps, MR_df, ...) {
      MR_df$outcome_snps <- list(tibble::tibble()); MR_df
    },
    mr_business_logic = function(MR_df, exposure_snps, ...) {
      list(MR_df = MR_df, results_df = tibble::tibble())
    },
    preprocess_exposure_snps = pp,
    .package = "ardmr",
    .env = parent.frame()
  )
  invisible(do.call(run_phenome_mr, args))
  cache_dir
}

read_runhash <- function(cache_dir) {
  out_root <- file.path(cache_dir, "output", "synth-exposure", "both", "EUR")
  list.dirs(out_root, recursive = FALSE, full.names = FALSE)
}

read_manifest <- function(cache_dir) {
  rh <- read_runhash(cache_dir)
  manifest_path <- file.path(cache_dir, "output", "synth-exposure", "both", "EUR",
                             rh, "run_manifest.json")
  jsonlite::read_json(manifest_path, simplifyVector = FALSE)
}

# ---- existing run_hash invariants (regression) ----------------------------

test_that("run_phenome_mr writes outputs into a <run_hash>/ subfolder", {
  cache_dir <- run_with_mocks(mk_args())
  rh <- read_runhash(cache_dir)
  expect_equal(length(rh), 1L)
  expect_match(rh, "^[0-9a-f]{5}__[0-9a-f]{5}$")
})

test_that("run_manifest.json is written with run_hash and key inputs", {
  cache_dir <- run_with_mocks(mk_args())
  rh <- read_runhash(cache_dir)
  manifest_path <- file.path(cache_dir, "output", "synth-exposure", "both", "EUR",
                             rh, "run_manifest.json")
  expect_true(file.exists(manifest_path))
  m <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  expect_equal(m$run_hash, rh)
  expect_equal(m$exposure, "synth_exposure")
  expect_equal(m$ancestry, "EUR")
  expect_equal(m$sex, "both")
  expect_equal(m$exposure_id, "ieu-a-test")
  expect_match(m$iv_hash, "^[0-9a-f]{5}$")
  expect_equal(substr(m$run_hash, 1L, 5L), m$iv_hash)
})

test_that("identical inputs produce identical run_hash", {
  c1 <- run_with_mocks(mk_args())
  c2 <- run_with_mocks(mk_args())
  expect_equal(read_runhash(c1), read_runhash(c2))
})

test_that("changing clump_opts -> different run_hash", {
  c1 <- run_with_mocks(mk_args())
  args2 <- mk_args(); args2$clump_opts <- list(p_threshold = 1e-6)
  c2 <- run_with_mocks(args2)
  expect_false(identical(read_runhash(c1), read_runhash(c2)))
})

test_that("changing exposure IVs -> different run_hash", {
  c1 <- run_with_mocks(mk_args(), mk_synth_exposure_snps(pos = 100L))
  c2 <- run_with_mocks(mk_args(), mk_synth_exposure_snps(pos = 200L))
  expect_false(identical(read_runhash(c1), read_runhash(c2)))
})

# ---- new manifest fields (Job 2c) ----------------------------------------

test_that("manifest includes clump_opts_resolved as a superset of user clump_opts", {
  args <- mk_args(); args$clump_opts <- list(maf_min = 0.05)
  cache_dir <- run_with_mocks(args)
  m <- read_manifest(cache_dir)
  expect_true("clump_opts_resolved" %in% names(m))
  # User-supplied key preserved
  expect_equal(m$clump_opts_resolved$maf_min, 0.05)
  # Defaults filled in
  expect_equal(m$clump_opts_resolved$r2, 0.001)
  expect_equal(m$clump_opts_resolved$kb, 10000L)
  expect_equal(m$clump_opts_resolved$f_threshold, 10)
  expect_equal(unlist(m$clump_opts_resolved$p_backoff), c(5e-8))   # Job 4 default
  expect_true(isTRUE(m$clump_opts_resolved$drop_indels))
  expect_false(isTRUE(m$clump_opts_resolved$drop_palindromic))
  expect_true(isTRUE(m$clump_opts_resolved$preprocess))
})

test_that("manifest preserves the raw clump_opts field separately from clump_opts_resolved", {
  args <- mk_args(); args$clump_opts <- list(maf_min = 0.05)
  cache_dir <- run_with_mocks(args)
  m <- read_manifest(cache_dir)
  # `clump_opts` is the raw user input (single field set, no defaults)
  expect_equal(names(m$clump_opts), "maf_min")
  expect_equal(m$clump_opts$maf_min, 0.05)
})

test_that("manifest includes preprocessing_steps with 8 entries each having step/n_in/n_out", {
  cache_dir <- run_with_mocks(mk_args())
  m <- read_manifest(cache_dir)
  expect_true("preprocessing_steps" %in% names(m))
  expect_length(m$preprocessing_steps, 8L)
  expected_steps <- c("p_threshold", "rsid_validate", "drop_indels",
                      "maf", "info", "palindromic", "clump", "f_stat")
  observed <- vapply(m$preprocessing_steps, function(x) x$step, character(1))
  expect_equal(observed, expected_steps)
  for (s in m$preprocessing_steps) {
    expect_true(all(c("step", "n_in", "n_out") %in% names(s)),
                info = sprintf("step '%s' missing required keys", s$step %||% "?"))
  }
})

test_that("preprocessing_steps roundtrips through JSON without losing structure", {
  cache_dir <- run_with_mocks(mk_args())
  m <- read_manifest(cache_dir)
  # Step `palindromic` carries a logical `dropped` field
  pal <- Filter(function(x) identical(x$step, "palindromic"), m$preprocessing_steps)[[1]]
  expect_false(isTRUE(pal$dropped))
  # Step `clump` carries r2/kb numerics
  clump <- Filter(function(x) identical(x$step, "clump"), m$preprocessing_steps)[[1]]
  expect_equal(clump$r2, 0.001)
  expect_equal(clump$kb, 10000L)
  # Skipped steps record their reason
  maf <- Filter(function(x) identical(x$step, "maf"), m$preprocessing_steps)[[1]]
  expect_equal(maf$skipped_reason, "maf_min not set")
})

# ---- deprecation: p_threshold in raw, p_backoff in resolved ---------------

test_that("deprecated p_threshold: raw retains p_threshold, resolved has p_backoff", {
  args <- mk_args(); args$clump_opts <- list(p_threshold = 1e-7)
  cache_dir <- run_with_mocks(args)
  m <- read_manifest(cache_dir)
  expect_equal(m$clump_opts$p_threshold, 1e-7)  # raw
  expect_equal(unlist(m$clump_opts_resolved$p_backoff), 1e-7)  # translated
  expect_null(m$clump_opts_resolved$p_threshold)
})

# ---- preprocess = FALSE bypass --------------------------------------------

test_that("preprocess = FALSE: empty preprocessing_steps in manifest, preprocessor never called", {
  call_count <- 0L
  bypass_stub <- function(snps_tbl, ...) {
    call_count <<- call_count + 1L
    list(snps = tibble::as_tibble(snps_tbl), steps = list(),
         resolved_opts = ardmr:::.resolve_clump_opts(list()))
  }
  args <- mk_args(); args$clump_opts <- list(preprocess = FALSE)
  # The resolver's finalize() bypasses preprocess when preprocess = FALSE,
  # so the stub should not be called even though we mock it.
  cache_dir <- run_with_mocks(args, preprocess_stub = bypass_stub)
  m <- read_manifest(cache_dir)
  expect_equal(call_count, 0L)
  expect_length(m$preprocessing_steps, 0L)
  expect_false(isTRUE(m$clump_opts_resolved$preprocess))
})

# ---- iv_hash uses POST-preprocess SNPs ------------------------------------

test_that("preprocessing that drops a SNP changes iv_hash (post-preprocess)", {
  drop_first_stub <- function(snps_tbl, ...) {
    keep <- snps_tbl[seq_len(max(0L, nrow(snps_tbl) - 1L)), , drop = FALSE]
    list(snps = tibble::as_tibble(keep), steps = list(),
         resolved_opts = ardmr:::.resolve_clump_opts(list()))
  }
  snps_two <- rbind(
    mk_synth_exposure_snps(pos = 100L),
    mk_synth_exposure_snps(pos = 200L)
  )
  c_full <- run_with_mocks(mk_args(), snps = snps_two)
  c_drop <- run_with_mocks(mk_args(), snps = snps_two, preprocess_stub = drop_first_stub)
  rh_full <- read_runhash(c_full)
  rh_drop <- read_runhash(c_drop)
  iv_full <- substr(rh_full, 1L, 5L)
  iv_drop <- substr(rh_drop, 1L, 5L)
  expect_false(identical(iv_full, iv_drop))
})

# ---- 0-IV preprocessing failure ------------------------------------------

# Job 4 changed 0-IV from a fatal stop() to a logger::log_warn() +
# early return. The capture helper below collects logger output for
# warning-text assertions; the function still returns its list.
.simple_console_appender <- function(lines) cat(lines, sep = "\n", file = stderr())
with_captured_logs <- function(expr) {
  captured <- character()
  prev_thr <- logger::log_threshold()
  logger::log_threshold(logger::TRACE)
  logger::log_appender(function(lines) captured <<- c(captured, lines))
  on.exit({
    logger::log_appender(.simple_console_appender)
    logger::log_threshold(prev_thr)
  }, add = TRUE)
  val <- force(expr)
  list(value = val, logs = captured)
}

# Run with mocks AND capture logs. Returns list(cache_dir, logs).
run_with_mocks_capture <- function(args, snps = NULL, preprocess_stub = NULL) {
  cap <- with_captured_logs(run_with_mocks(args, snps, preprocess_stub))
  list(cache_dir = cap$value, logs = cap$logs)
}

test_that("0-IV preprocessing -> warning lists per-step counts (no error)", {
  zero_stub <- function(snps_tbl, ...) {
    n <- NROW(snps_tbl)
    list(
      snps  = tibble::as_tibble(snps_tbl)[0, , drop = FALSE],
      steps = list(
        list(step = "p_threshold",   n_in = n, n_out = n, value = 5e-8),
        list(step = "rsid_validate", n_in = n, n_out = n),
        list(step = "drop_indels",   n_in = n, n_out = 0L),
        list(step = "maf",           n_in = 0L, n_out = 0L, value = NA_real_, skipped_reason = "maf_min not set"),
        list(step = "info",          n_in = 0L, n_out = 0L, value = NA_real_, skipped_reason = "info_min not set"),
        list(step = "palindromic",   n_in = 0L, n_out = 0L, dropped = FALSE, skipped_reason = "drop_palindromic = FALSE"),
        list(step = "clump",         n_in = 0L, n_out = 0L, r2 = 0.001, kb = 10000L),
        list(step = "f_stat",        n_in = 0L, n_out = 0L, value = 10)
      ),
      resolved_opts = ardmr:::.resolve_clump_opts(list())
    )
  }
  cap <- run_with_mocks_capture(mk_args(), preprocess_stub = zero_stub)
  expect_true(any(grepl("Preprocessing returned 0 instruments", cap$logs)))
})

test_that("0-IV warning identifies the killing step", {
  zero_stub <- function(snps_tbl, ...) {
    n <- NROW(snps_tbl)
    list(
      snps  = tibble::as_tibble(snps_tbl)[0, , drop = FALSE],
      steps = list(
        list(step = "p_threshold",   n_in = n, n_out = n, value = 5e-8),
        list(step = "rsid_validate", n_in = n, n_out = n),
        list(step = "drop_indels",   n_in = n, n_out = n),
        list(step = "maf",           n_in = n, n_out = n, value = NA_real_, skipped_reason = "maf_min not set"),
        list(step = "info",          n_in = n, n_out = n, value = NA_real_, skipped_reason = "info_min not set"),
        list(step = "palindromic",   n_in = n, n_out = n, dropped = FALSE, skipped_reason = "drop_palindromic = FALSE"),
        list(step = "clump",         n_in = n, n_out = n, r2 = 0.001, kb = 10000L),
        list(step = "f_stat",        n_in = n, n_out = 0L, value = 10)
      ),
      resolved_opts = ardmr:::.resolve_clump_opts(list())
    )
  }
  cap <- run_with_mocks_capture(mk_args(), preprocess_stub = zero_stub)
  expect_true(any(grepl("f_stat", cap$logs)))
})

test_that("0-IV warning includes per-step counts in n_in -> n_out form", {
  zero_stub <- function(snps_tbl, ...) {
    n <- NROW(snps_tbl)
    list(
      snps  = tibble::as_tibble(snps_tbl)[0, , drop = FALSE],
      steps = list(
        list(step = "p_threshold",   n_in = n, n_out = 0L, value = 5e-8),
        list(step = "rsid_validate", n_in = 0L, n_out = 0L),
        list(step = "drop_indels",   n_in = 0L, n_out = 0L),
        list(step = "maf",           n_in = 0L, n_out = 0L, value = NA_real_, skipped_reason = "maf_min not set"),
        list(step = "info",          n_in = 0L, n_out = 0L, value = NA_real_, skipped_reason = "info_min not set"),
        list(step = "palindromic",   n_in = 0L, n_out = 0L, dropped = FALSE, skipped_reason = "drop_palindromic = FALSE"),
        list(step = "clump",         n_in = 0L, n_out = 0L, r2 = 0.001, kb = 10000L),
        list(step = "f_stat",        n_in = 0L, n_out = 0L, value = 10)
      ),
      resolved_opts = ardmr:::.resolve_clump_opts(list())
    )
  }
  cap <- run_with_mocks_capture(mk_args(), preprocess_stub = zero_stub)
  expect_true(any(grepl("1 -> 0", cap$logs)))  # the killing step's n_in -> n_out
})

test_that("0-IV with empty preprocessing_steps -> bypass-mode warning", {
  bypass_zero_stub <- function(snps_tbl, ...) {
    list(snps = tibble::as_tibble(snps_tbl)[0, , drop = FALSE],
         steps = list(),
         resolved_opts = ardmr:::.resolve_clump_opts(list()))
  }
  args <- mk_args(); args$clump_opts <- list(preprocess = FALSE)
  empty_snps <- mk_synth_exposure_snps()[0, ]
  cap <- run_with_mocks_capture(args, snps = empty_snps,
                                preprocess_stub = bypass_zero_stub)
  expect_true(any(grepl("no preprocessing steps were recorded", cap$logs)))
})

# ---- 0-IV: return list shape + manifest still written --------------------

run_with_mocks_return <- function(args, snps = NULL, preprocess_stub = NULL) {
  # Identical to run_with_mocks but returns the function result instead
  # of just the cache_dir, so the caller can inspect the result list.
  cache_dir <- tempfile("ardmr_run_")
  dir.create(cache_dir, recursive = TRUE)
  args$cache_dir <- cache_dir
  args$exposure_snps <- if (is.null(snps)) mk_synth_exposure_snps() else snps
  pp <- if (is.null(preprocess_stub)) .default_preprocess_stub else preprocess_stub
  testthat::local_mocked_bindings(
    Outcome_setup = function(sex, ancestry) {
      tibble::tibble(description = "synth_outcome", ARD_selected = FALSE,
                     aws_link = "https://invalid.example/x.tsv.bgz",
                     aws_link_tabix = "https://invalid.example/x.tsv.bgz.tbi")
    },
    Variant_manifest_downloader = function(...) invisible(NULL),
    exposure_snp_mapper = function(exposure_snps, ...) {
      exposure_snps$panukb_chrom <- 1L; exposure_snps$panukb_pos <- 100L
      exposure_snps$rsid <- exposure_snps$SNP; exposure_snps
    },
    compute_ld_matrix = function(...) list(matrix = matrix(1.0, 1L, 1L), snps_used = "rs1"),
    panukb_snp_grabber = function(exposure_snps, MR_df, ...) {
      MR_df$outcome_snps <- list(tibble::tibble()); MR_df
    },
    mr_business_logic = function(MR_df, exposure_snps, ...) {
      list(MR_df = MR_df, results_df = tibble::tibble())
    },
    preprocess_exposure_snps = pp,
    .package = "ardmr",
    .env = parent.frame()
  )
  res <- do.call(run_phenome_mr, args)
  list(result = res, cache_dir = cache_dir)
}

test_that("0-IV: returned list has status = 'no_ivs' + empty MR/results tables + NULL plots", {
  zero_stub <- function(snps_tbl, ...) {
    n <- NROW(snps_tbl)
    list(
      snps  = tibble::as_tibble(snps_tbl)[0, , drop = FALSE],
      steps = list(
        list(step = "p_threshold",   n_in = n, n_out = 0L, value = 5e-8),
        list(step = "rsid_validate", n_in = 0L, n_out = 0L),
        list(step = "drop_indels",   n_in = 0L, n_out = 0L),
        list(step = "maf",           n_in = 0L, n_out = 0L, value = NA_real_, skipped_reason = "maf_min not set"),
        list(step = "info",          n_in = 0L, n_out = 0L, value = NA_real_, skipped_reason = "info_min not set"),
        list(step = "palindromic",   n_in = 0L, n_out = 0L, dropped = FALSE, skipped_reason = "drop_palindromic = FALSE"),
        list(step = "clump",         n_in = 0L, n_out = 0L, r2 = 0.001, kb = 10000L),
        list(step = "f_stat",        n_in = 0L, n_out = 0L, value = 10)
      ),
      resolved_opts = ardmr:::.resolve_clump_opts(list())
    )
  }
  out <- suppressWarnings(run_with_mocks_return(mk_args(), preprocess_stub = zero_stub))
  expect_equal(out$result$status, "no_ivs")
  expect_equal(nrow(out$result$MR_df), 0L)
  expect_equal(nrow(out$result$results_df), 0L)
  expect_null(out$result$manhattan)
  expect_null(out$result$volcano)
  expect_true(nzchar(out$result$run_dir))
})

test_that("0-IV: manifest is still written so the user can inspect the steps", {
  zero_stub <- function(snps_tbl, ...) {
    n <- NROW(snps_tbl)
    list(
      snps  = tibble::as_tibble(snps_tbl)[0, , drop = FALSE],
      steps = list(
        list(step = "p_threshold",   n_in = n, n_out = 0L, value = 5e-8),
        list(step = "rsid_validate", n_in = 0L, n_out = 0L),
        list(step = "drop_indels",   n_in = 0L, n_out = 0L),
        list(step = "maf",           n_in = 0L, n_out = 0L, value = NA_real_, skipped_reason = "maf_min not set"),
        list(step = "info",          n_in = 0L, n_out = 0L, value = NA_real_, skipped_reason = "info_min not set"),
        list(step = "palindromic",   n_in = 0L, n_out = 0L, dropped = FALSE, skipped_reason = "drop_palindromic = FALSE"),
        list(step = "clump",         n_in = 0L, n_out = 0L, r2 = 0.001, kb = 10000L),
        list(step = "f_stat",        n_in = 0L, n_out = 0L, value = 10)
      ),
      resolved_opts = ardmr:::.resolve_clump_opts(list())
    )
  }
  out <- suppressWarnings(run_with_mocks_return(mk_args(), preprocess_stub = zero_stub))
  manifest_path <- file.path(out$result$run_dir, "run_manifest.json")
  expect_true(file.exists(manifest_path))
  m <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  expect_equal(m$n_ivs, 0L)
  expect_length(m$preprocessing_steps, 8L)
})

# ---- .format_zero_iv_error unit tests (the helper itself) ----------------

test_that(".format_zero_iv_error: empty steps -> bypass message", {
  msg <- ardmr:::.format_zero_iv_error(list())
  expect_match(msg, "no preprocessing steps were recorded")
  expect_match(msg, "preprocess = FALSE")
})

test_that(".format_zero_iv_error: lists every step with n_in -> n_out", {
  steps <- list(
    list(step = "p_threshold", n_in = 100L, n_out = 50L, value = 5e-8),
    list(step = "clump",       n_in = 50L,  n_out = 0L,  r2 = 0.001, kb = 10000L)
  )
  msg <- ardmr:::.format_zero_iv_error(steps)
  expect_match(msg, "Preprocessing returned 0 instruments")
  expect_match(msg, "p_threshold")
  expect_match(msg, "100 -> 50")
  expect_match(msg, "50 -> 0")
})

test_that(".format_zero_iv_error: marks the first step that took count to 0", {
  steps <- list(
    list(step = "p_threshold", n_in = 100L, n_out = 100L, value = 5e-8),
    list(step = "drop_indels", n_in = 100L, n_out = 0L),
    list(step = "f_stat",      n_in = 0L,   n_out = 0L,   value = 10)
  )
  msg <- ardmr:::.format_zero_iv_error(steps)
  expect_match(msg, "drop_indels.*<- everything dropped here")
  expect_match(msg, "Likely cause: 'drop_indels'")
  # `f_stat` (which takes 0 -> 0) should NOT be marked as the killer
  killer_count <- length(regmatches(msg, gregexpr("everything dropped here", msg))[[1]])
  expect_equal(killer_count, 1L)
})

test_that(".format_zero_iv_error: handles steps with NULL value gracefully", {
  steps <- list(
    list(step = "rsid_validate", n_in = 10L, n_out = 0L),
    list(step = "maf", n_in = 0L, n_out = 0L, value = NULL,
         skipped_reason = "maf_min not set")
  )
  msg <- ardmr:::.format_zero_iv_error(steps)
  expect_match(msg, "rsid_validate")
  expect_match(msg, "skipped: maf_min not set")
})

# ---- confirm propagation ---------------------------------------------------

test_that("run_phenome_mr's `confirm` arg is plumbed through to the resolver", {
  observed <- new.env()
  testthat::local_mocked_bindings(
    .resolve_coloc_source = function(exposure_snps, exposure_sumstats, exposure_id,
                                     ancestry, genome_build, acknowledge_no_coloc,
                                     clump_opts, cache_dir, confirm = "ask",
                                     verbose = TRUE) {
      observed$confirm <- confirm
      list(
        mode                = "snps_id",
        exposure_snps       = mk_synth_exposure_snps(),
        coloc_fetcher       = function(chr, start, end) tibble::tibble(),
        coloc_metadata      = list(N = 100000, type = "quant", ncase = NA_real_,
                                   ncontrol = NA_real_, s = NA_real_),
        clump_opts_resolved = ardmr:::.resolve_clump_opts(clump_opts),
        preprocessing_steps = list()
      )
    }
  )
  # Use the existing mock chain so the rest of the pipeline runs.
  args <- mk_args(); args$confirm <- "yes"
  run_with_mocks(args)
  expect_equal(observed$confirm, "yes")
})
