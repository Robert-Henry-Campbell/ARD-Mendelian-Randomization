# tests/testthat/test_coloc_source_resolver.R
# Tests for .resolve_coloc_source(). Covers (a) input validation and mode
# dispatch, (b) per-mode wiring of preprocess_exposure_snps with the right
# already_p_filtered/already_clumped flags, (c) source-layer p-value
# backoff for the ieugwasr branch, (d) preprocess = FALSE bypass, (e)
# clump_opts merge / p_threshold deprecation translation.

skip_if_not_installed("Rsamtools")

mk_snps <- function() {
  tibble::tibble(
    SNP = "rs1", id.exposure = "x", exposure = "x",
    effect_allele.exposure = "A", other_allele.exposure = "G",
    beta.exposure = 0.05, se.exposure = 0.01,
    eaf.exposure = 0.30, pval.exposure = 1e-10,
    samplesize.exposure = 100000
  )
}

mk_fake_meta <- function() {
  list(N = 100000, type = "quant",
       ncase = NA_real_, ncontrol = NA_real_, s = NA_real_)
}

# Recording mock for preprocess_exposure_snps -- replays the input as the
# output (so per-mode behaviour can be observed without needing the LD
# reference).
mk_preprocess_recorder <- function(captured_env, return_snps = NULL) {
  function(snps_tbl, clump_opts = list(), ancestry,
           cache_dir = ardmr_cache_dir(), confirm = "ask", verbose = TRUE) {
    captured_env$call_count    <- (captured_env$call_count %||% 0L) + 1L
    captured_env$last_co       <- clump_opts
    captured_env$last_snps_in  <- snps_tbl
    captured_env$last_ancestry <- ancestry
    out <- if (is.null(return_snps)) snps_tbl else return_snps
    list(
      snps          = tibble::as_tibble(out),
      steps         = list(list(step = "p_threshold", n_in = NROW(snps_tbl),
                                n_out = NROW(out), value = NA_real_)),
      resolved_opts = ardmr:::.resolve_clump_opts(clump_opts)
    )
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- input validation -----------------------------------------------------

test_that("nothing supplied -> error", {
  expect_error(
    ardmr:::.resolve_coloc_source(ancestry = "EUR"),
    "Provide one of"
  )
})

test_that("id + sumstats -> error", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  expect_error(
    ardmr:::.resolve_coloc_source(
      exposure_id = "ieu-b-38", exposure_sumstats = vcf, ancestry = "EUR"
    ),
    "both provided"
  )
})

test_that("snps_only without acknowledge_no_coloc -> error", {
  expect_error(
    ardmr:::.resolve_coloc_source(exposure_snps = mk_snps(), ancestry = "EUR"),
    "acknowledge_no_coloc"
  )
})

# ---- snps_only ------------------------------------------------------------

test_that("snps_only: preprocess called with already_clumped=FALSE, already_p_filtered=FALSE", {
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_snps = mk_snps(), ancestry = "EUR",
    acknowledge_no_coloc = TRUE, verbose = FALSE
  )
  expect_equal(r$mode, "snps_only")
  expect_null(r$coloc_fetcher)
  expect_null(r$coloc_metadata)
  expect_equal(cap$call_count, 1L)
  expect_false(isTRUE(cap$last_co$already_clumped))
  expect_false(isTRUE(cap$last_co$already_p_filtered))
  expect_true("clump_opts_resolved" %in% names(r))
  expect_true("preprocessing_steps" %in% names(r))
  expect_true(is.list(r$preprocessing_steps))
})

# ---- snps_vcf -------------------------------------------------------------

test_that("snps_vcf: preprocess called once + VCF fetcher attached", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_snps = mk_snps(), exposure_sumstats = vcf,
    ancestry = "EUR", verbose = FALSE
  )
  expect_equal(r$mode, "snps_vcf")
  expect_true(is.function(r$coloc_fetcher))
  expect_true(is.list(r$coloc_metadata))
  expect_equal(cap$call_count, 1L)
  expect_false(isTRUE(cap$last_co$already_clumped))
  expect_false(isTRUE(cap$last_co$already_p_filtered))
  out <- r$coloc_fetcher(1, 1L, 50000L)
  expect_true(nrow(out) > 0)
})

test_that("snps_vcf with GRCh38 declared -> error", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  expect_error(
    ardmr:::.resolve_coloc_source(
      exposure_snps = mk_snps(), exposure_sumstats = vcf,
      ancestry = "EUR", genome_build = "GRCh38"
    ),
    "GRCh37"
  )
})

# ---- snps_id --------------------------------------------------------------

test_that("snps_id: preprocess called once + ieugwasr fetcher attached", {
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) mk_fake_meta()
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_snps = mk_snps(), exposure_id = "ieu-b-38",
    ancestry = "EUR", cache_dir = tempdir(), verbose = FALSE
  )
  expect_equal(r$mode, "snps_id")
  expect_true(is.function(r$coloc_fetcher))
  expect_equal(r$coloc_metadata$N, 100000)
  expect_equal(cap$call_count, 1L)
  expect_false(isTRUE(cap$last_co$already_clumped))
  expect_false(isTRUE(cap$last_co$already_p_filtered))
})

# ---- vcf_only -------------------------------------------------------------

test_that("vcf_only: clump_sumstats_to_ivs called once with merged opts", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  fake_ivs <- tibble::tibble(
    SNP = c("rs1","rs2"), Chr = c(1L, 2L), Pos = c(1000L, 2000L),
    effect_allele.exposure = c("A","C"), other_allele.exposure = c("G","T"),
    beta.exposure = c(0.1, 0.05), se.exposure = c(0.01, 0.01),
    eaf.exposure = c(0.3, 0.4), pval.exposure = c(1e-10, 1e-9),
    samplesize.exposure = c(100000, 100000),
    id.exposure = "mini.bgz", exposure = "mini.bgz"
  )
  observed <- new.env()
  testthat::local_mocked_bindings(
    clump_sumstats_to_ivs = function(vcf_path, ancestry, p_threshold, r2, kb,
                                     f_threshold, genome_build, cache_dir, verbose) {
      observed$count       <- (observed$count %||% 0L) + 1L
      observed$p_threshold <- p_threshold
      observed$r2          <- r2
      observed$kb          <- kb
      observed$f_threshold <- f_threshold
      fake_ivs
    }
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_sumstats = vcf, ancestry = "EUR", verbose = FALSE
  )
  expect_equal(r$mode, "vcf_only")
  expect_equal(observed$count, 1L)
  expect_equal(observed$p_threshold, 5e-8)  # strictest p_backoff rung
  expect_equal(observed$r2, 0.001)
  expect_equal(observed$kb, 10000L)
  expect_equal(observed$f_threshold, 10)
  expect_equal(nrow(r$exposure_snps), 2L)
  expect_true("clump_opts_resolved" %in% names(r))
  expect_true("preprocessing_steps" %in% names(r))
})

test_that("vcf_only: 0-row return is NOT an error (downstream decides)", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  testthat::local_mocked_bindings(
    clump_sumstats_to_ivs = function(...) ardmr:::.empty_exposure_tbl()
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_sumstats = vcf, ancestry = "EUR", verbose = FALSE
  )
  expect_equal(r$mode, "vcf_only")
  expect_equal(nrow(r$exposure_snps), 0L)
})

# ---- ieugwasr (server clump) ----------------------------------------------

test_that("ieugwasr + prefer_server_clump=TRUE: extract_instruments called, preprocess flagged", {
  fake_ivs <- mk_snps()
  testthat::local_mocked_bindings(
    extract_instruments = function(outcomes, p1, clump, r2, kb, ...) fake_ivs,
    .package = "TwoSampleMR"
  )
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) mk_fake_meta()
  )
  ieugwasr_called <- 0L
  testthat::local_mocked_bindings(
    ld_clump = function(...) { ieugwasr_called <<- ieugwasr_called + 1L; tibble::tibble() },
    .package = "ieugwasr"
  )
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_id = "ieu-b-38", ancestry = "EUR",
    cache_dir = tempdir(), verbose = FALSE
  )
  expect_equal(r$mode, "ieugwasr")
  expect_equal(cap$call_count, 1L)
  expect_true(isTRUE(cap$last_co$already_p_filtered))
  expect_true(isTRUE(cap$last_co$already_clumped))
  expect_equal(ieugwasr_called, 0L)  # no double-clump
})

test_that("ieugwasr + prefer_server_clump=FALSE: safe_tophits called, only p flagged", {
  fake_top <- tibble::tibble(
    rsid = "rs1", ea = "A", nea = "G",
    beta = 0.05, se = 0.01, eaf = 0.3, p = 1e-10, n = 100000
  )
  testthat::local_mocked_bindings(
    tophits = function(id, p = NULL, pval = NULL, ...) fake_top,
    .package = "ieugwasr"
  )
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) mk_fake_meta()
  )
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_id = "ieu-b-38", ancestry = "EUR",
    clump_opts = list(prefer_server_clump = FALSE),
    cache_dir = tempdir(), verbose = FALSE
  )
  expect_equal(r$mode, "ieugwasr")
  expect_equal(cap$call_count, 1L)
  expect_true(isTRUE(cap$last_co$already_p_filtered))
  expect_false(isTRUE(cap$last_co$already_clumped))
  # The IVs sent to preprocess should be the .tophits_to_tsmr-converted
  # tibble (TwoSampleMR shape).
  expect_true("effect_allele.exposure" %in% names(cap$last_snps_in))
  expect_equal(cap$last_snps_in$SNP, "rs1")
  expect_equal(cap$last_snps_in$effect_allele.exposure, "A")
})

# ---- ieugwasr backoff at source -------------------------------------------

test_that("ieugwasr backoff: 0 at strictest, hits at next rung -> called twice", {
  call_log <- new.env(); call_log$pthr <- numeric(); call_log$count <- 0L
  testthat::local_mocked_bindings(
    extract_instruments = function(outcomes, p1, clump, r2, kb, ...) {
      call_log$count <- call_log$count + 1L
      call_log$pthr  <- c(call_log$pthr, p1)
      if (p1 < 5e-7) NULL else mk_snps()
    },
    .package = "TwoSampleMR"
  )
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) mk_fake_meta()
  )
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_id = "ieu-b-38", ancestry = "EUR",
    cache_dir = tempdir(), verbose = FALSE
  )
  expect_equal(r$mode, "ieugwasr")
  expect_equal(call_log$count, 2L)
  expect_equal(call_log$pthr, c(5e-8, 5e-7))
})

test_that("ieugwasr backoff: 0 at all rungs -> 0 IVs, no error", {
  call_log <- new.env(); call_log$count <- 0L
  testthat::local_mocked_bindings(
    extract_instruments = function(outcomes, p1, clump, r2, kb, ...) {
      call_log$count <- call_log$count + 1L
      NULL
    },
    .package = "TwoSampleMR"
  )
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) mk_fake_meta()
  )
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_id = "ieu-b-38", ancestry = "EUR",
    cache_dir = tempdir(), verbose = FALSE
  )
  expect_equal(r$mode, "ieugwasr")
  expect_equal(call_log$count, 3L)  # all three rungs tried
  expect_equal(nrow(r$exposure_snps), 0L)
})

# ---- preprocess = FALSE bypass --------------------------------------------

test_that("preprocess = FALSE: preprocess_exposure_snps never called, input passed through", {
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_snps = mk_snps(), ancestry = "EUR",
    acknowledge_no_coloc = TRUE,
    clump_opts = list(preprocess = FALSE),
    verbose = FALSE
  )
  expect_equal(r$mode, "snps_only")
  expect_null(cap$call_count)
  expect_equal(nrow(r$exposure_snps), 1L)
  expect_equal(r$preprocessing_steps, list())
  expect_false(isTRUE(r$clump_opts_resolved$preprocess))
})

# ---- clump_opts merge -----------------------------------------------------

test_that("defaults merge: user maf_min preserved, defaults filled in", {
  cap <- new.env()
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(cap)
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_snps = mk_snps(), ancestry = "EUR",
    acknowledge_no_coloc = TRUE,
    clump_opts = list(maf_min = 0.05),
    verbose = FALSE
  )
  expect_equal(r$clump_opts_resolved$maf_min, 0.05)
  expect_equal(r$clump_opts_resolved$r2, 0.001)
  expect_equal(r$clump_opts_resolved$kb, 10000L)
  expect_equal(r$clump_opts_resolved$f_threshold, 10)
  expect_equal(r$clump_opts_resolved$p_backoff, c(5e-8, 5e-7, 5e-6))
})

# ---- p_threshold deprecation through the resolver ------------------------

.simple_console_appender <- function(lines) cat(lines, sep = "\n", file = stderr())
with_captured_logs <- function(expr) {
  captured <- character()
  prev_thr <- logger::log_threshold()
  logger::log_threshold(logger::TRACE)  # capture INFO and below; reset on exit
  logger::log_appender(function(lines) captured <<- c(captured, lines))
  on.exit({
    logger::log_appender(.simple_console_appender)
    logger::log_threshold(prev_thr)
  }, add = TRUE)
  val <- force(expr)
  list(value = val, logs = captured)
}

test_that("deprecation: clump_opts$p_threshold translated to p_backoff via resolver", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  observed <- new.env()
  testthat::local_mocked_bindings(
    clump_sumstats_to_ivs = function(vcf_path, ancestry, p_threshold, r2, kb,
                                     f_threshold, genome_build, cache_dir, verbose) {
      observed$p_threshold <- p_threshold
      observed$r2          <- r2
      observed$kb          <- kb
      observed$f_threshold <- f_threshold
      tibble::tibble(
        SNP = "rs1", Chr = 1L, Pos = 1000L,
        effect_allele.exposure = "A", other_allele.exposure = "G",
        beta.exposure = 0.1, se.exposure = 0.01,
        eaf.exposure = 0.3, pval.exposure = 1e-10,
        samplesize.exposure = 100000,
        id.exposure = "mini.bgz", exposure = "mini.bgz"
      )
    }
  )
  cap <- with_captured_logs(
    ardmr:::.resolve_coloc_source(
      exposure_sumstats = vcf, ancestry = "EUR",
      clump_opts = list(p_threshold = 1e-6, r2 = 0.01, kb = 5000L, f_threshold = 5),
      verbose = FALSE
    )
  )
  expect_true(any(grepl("p_threshold is deprecated", cap$logs)))
  # The strictest (and only) p_backoff rung is 1e-6 after translation.
  expect_equal(observed$p_threshold, 1e-6)
  expect_equal(observed$r2, 0.01)
  expect_equal(observed$kb, 5000L)
  expect_equal(observed$f_threshold, 5)
  expect_equal(cap$value$clump_opts_resolved$p_backoff, 1e-6)
  expect_null(cap$value$clump_opts_resolved$p_threshold)
})

# ---- ld_pop alignment notice ----------------------------------------------

test_that("ld_pop alignment: non-EUR ancestry + prefer_server_clump=TRUE logs note", {
  testthat::local_mocked_bindings(
    extract_instruments = function(outcomes, p1, ...) mk_snps(),
    .package = "TwoSampleMR"
  )
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) mk_fake_meta()
  )
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(new.env())
  )
  cap <- with_captured_logs(
    ardmr:::.resolve_coloc_source(
      exposure_id = "ieu-b-38", ancestry = "AFR",
      cache_dir = tempdir(), verbose = TRUE
    )
  )
  expect_true(any(grepl("server clumping uses its default", cap$logs)))
})

test_that("ld_pop alignment: EUR ancestry does NOT log the note", {
  testthat::local_mocked_bindings(
    extract_instruments = function(outcomes, p1, ...) mk_snps(),
    .package = "TwoSampleMR"
  )
  testthat::local_mocked_bindings(
    coloc_fetch_metadata = function(exposure_id, ...) mk_fake_meta()
  )
  testthat::local_mocked_bindings(
    preprocess_exposure_snps = mk_preprocess_recorder(new.env())
  )
  cap <- with_captured_logs(
    ardmr:::.resolve_coloc_source(
      exposure_id = "ieu-b-38", ancestry = "EUR",
      cache_dir = tempdir(), verbose = TRUE
    )
  )
  expect_false(any(grepl("server clumping uses its default", cap$logs)))
})
