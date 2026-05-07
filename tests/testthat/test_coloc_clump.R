# tests/testthat/test_coloc_clump.R
#
# Unit tests for clump_sumstats_to_ivs() AFTER its Job 2b refactor: the
# function now does VCF query -> reshape to TwoSampleMR shape -> delegate
# to preprocess_exposure_snps() for the unified pipeline (p-backoff,
# clump, rsid, indel, palindromic, F-stat, MAF, INFO).
#
# We use mockery::stub to patch gwasvcf::* and validate_gwasvcf so the
# tests run without gwasvcf installed. Tests that exercise the full
# preprocess pipeline additionally mock the LD machinery via
# local_mocked_bindings.

skip_if_not_installed("mockery")
skip_if_not_installed("ieugwasr")

# ---- shared helpers -------------------------------------------------------

# Builder for a synthetic VCF "hits" tibble in the shape returned by
# gwasvcf::vcf_to_tibble. Defaults: 5 SNVs at distinct positions.
mk_hits_tbl <- function(n = 5,
                        rsid     = paste0("rs", seq_len(n)),
                        seqnames = rep(1L, n),
                        start    = seq.int(1000L, length.out = n, by = 1000L),
                        REF      = c("A","C","G","T","A")[seq_len(n)],
                        ALT      = c("G","T","A","C","G")[seq_len(n)],
                        ES       = c(0.1, 0.05, 0.005, 0.08, 0.001)[seq_len(n)],
                        SE       = rep(0.01, n),
                        LP       = c(20, 8, 0.5, 12, 0.1)[seq_len(n)],
                        AF       = rep(0.3, n),
                        SS       = rep(100000L, n)) {
  tibble::tibble(
    rsid = rsid, seqnames = seqnames, start = start,
    REF = REF, ALT = ALT,
    ES = ES, SE = SE, LP = LP, AF = AF, SS = SS
  )
}

mk_fake_ld_ref <- function() {
  list(bfile_prefix = file.path(tempdir(), "fake_bfile"),
       present = TRUE, downloaded = FALSE, aborted = FALSE)
}

# Build a stubbed clump_sumstats_to_ivs that bypasses gwasvcf and
# validate_gwasvcf. Returns the modified function -- callers run it
# directly. `hits` controls what gwasvcf::vcf_to_tibble returns;
# `obs_env` (optional) captures the pval passed to query_gwas.
make_clump_fn <- function(hits, obs_env = NULL) {
  fn <- clump_sumstats_to_ivs
  mockery::stub(fn, "validate_gwasvcf", function(...) invisible(TRUE))
  mockery::stub(fn, "gwasvcf::query_gwas", function(vcffile, pval, ...) {
    if (!is.null(obs_env)) obs_env$pval <- pval
    structure(list())
  })
  mockery::stub(fn, "gwasvcf::vcf_to_tibble", function(vcf) hits)
  # Also stub the requireNamespace check so it returns TRUE without gwasvcf.
  mockery::stub(fn, "requireNamespace", function(pkg, ...) TRUE)
  fn
}

# Mock the LD machinery so preprocess_exposure_snps can clump offline.
mock_ld_stack <- function(ld_clump_fn = function(dat, ...) dat) {
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) mk_fake_ld_ref(),
    .find_plink_bin = function(plink_bin = NULL) NULL,
    .env = parent.frame()
  )
  testthat::local_mocked_bindings(
    ld_clump = ld_clump_fn,
    .package = "ieugwasr",
    .env = parent.frame()
  )
}

# ---- public-signature regression -----------------------------------------

test_that("legacy scalar interface (p_threshold/r2/kb/f_threshold) still works", {
  fn <- make_clump_fn(mk_hits_tbl())
  mock_ld_stack()
  ivs <- fn(
    "fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
    p_threshold = 1e-7, r2 = 0.001, kb = 10000L,
    f_threshold = 0, verbose = FALSE
  )
  expect_s3_class(ivs, "tbl_df")
  expect_silent(assert_exposure(ivs))
})

test_that("query_gwas is called with the LOOSEST p_backoff rung", {
  obs <- new.env()
  fn <- make_clump_fn(mk_hits_tbl(), obs_env = obs)
  mock_ld_stack()
  fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
     clump_opts = list(p_backoff = c(5e-8, 5e-7, 5e-6)),
     verbose = FALSE)
  expect_equal(obs$pval, 5e-6)
})

test_that("scalar p_threshold becomes a single-rung ladder", {
  obs <- new.env()
  fn <- make_clump_fn(mk_hits_tbl(), obs_env = obs)
  mock_ld_stack()
  fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
     p_threshold = 1e-9, verbose = FALSE)
  expect_equal(obs$pval, 1e-9)
})

test_that("clump_opts overrides scalar args (r2, kb)", {
  obs_clump <- new.env()
  fn <- make_clump_fn(mk_hits_tbl())
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) mk_fake_ld_ref(),
    .find_plink_bin = function(plink_bin = NULL) NULL
  )
  testthat::local_mocked_bindings(
    ld_clump = function(dat, clump_r2, clump_kb, ...) {
      obs_clump$r2 <- clump_r2
      obs_clump$kb <- clump_kb
      dat
    },
    .package = "ieugwasr"
  )
  fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
     r2 = 0.001, kb = 10000L,                       # scalar args
     clump_opts = list(r2 = 0.005, kb = 1234L),     # overrides
     verbose = FALSE)
  expect_equal(obs_clump$r2, 0.005)
  expect_equal(obs_clump$kb, 1234L)
})

# ---- VCF -> tibble -> preprocess flow -------------------------------------

test_that("preprocess_exposure_snps invoked with already_clumped/p_filtered = FALSE", {
  cap <- new.env()
  fn <- make_clump_fn(mk_hits_tbl())
  mockery::stub(fn, "preprocess_exposure_snps",
    function(snps_tbl, clump_opts = list(), ancestry, cache_dir,
             confirm = "ask", verbose = TRUE) {
      cap$call_count    <- (if (is.null(cap$call_count)) 0L else cap$call_count) + 1L
      cap$last_co       <- clump_opts
      cap$last_snps     <- snps_tbl
      cap$last_ancestry <- ancestry
      list(snps = snps_tbl, steps = list(), resolved_opts = clump_opts)
    })
  fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
     clump_opts = list(p_backoff = c(5e-8), already_clumped = TRUE,
                        already_p_filtered = TRUE),
     verbose = FALSE)
  expect_equal(cap$call_count, 1L)
  # Forced inside clump_sumstats_to_ivs regardless of caller flags
  expect_false(isTRUE(cap$last_co$already_clumped))
  expect_false(isTRUE(cap$last_co$already_p_filtered))
  # The TwoSampleMR-shape conversion happened in clump_sumstats_to_ivs
  expect_true(all(c("SNP", "effect_allele.exposure", "other_allele.exposure",
                    "beta.exposure", "se.exposure", "pval.exposure",
                    "id.exposure") %in% names(cap$last_snps)))
  expect_equal(cap$last_ancestry, "EUR")
})

# ---- empty / degenerate inputs --------------------------------------------

test_that("empty VCF query -> empty TwoSampleMR-formatted tibble (no error)", {
  fn <- make_clump_fn(mk_hits_tbl(n = 5)[0, ])
  mock_ld_stack()
  ivs <- fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(), verbose = FALSE)
  expect_equal(nrow(ivs), 0)
  expect_true(all(c("SNP","Chr","Pos","effect_allele.exposure",
                    "other_allele.exposure","beta.exposure","se.exposure",
                    "id.exposure") %in% names(ivs)))
})

test_that("rows with non-finite beta/se are dropped before reaching preprocess", {
  hits <- mk_hits_tbl(n = 3)
  hits$ES[2] <- NA_real_       # this row should be dropped pre-preprocess
  hits$SE[3] <- 0              # se = 0 -> dropped pre-preprocess
  cap <- new.env()
  fn <- make_clump_fn(hits)
  mockery::stub(fn, "preprocess_exposure_snps",
    function(snps_tbl, clump_opts = list(), ancestry, cache_dir,
             confirm = "ask", verbose = TRUE) {
      cap$nrow_in <- nrow(snps_tbl)
      list(snps = snps_tbl, steps = list(), resolved_opts = clump_opts)
    })
  fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(), verbose = FALSE)
  expect_equal(cap$nrow_in, 1L)
})

# ---- F-stat / clumping happen INSIDE preprocess (end-to-end) -------------

test_that("F-statistic threshold filters weak instruments (via preprocess)", {
  fn <- make_clump_fn(mk_hits_tbl())
  mock_ld_stack()
  ivs <- fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
            p_threshold = 1, f_threshold = 50, verbose = FALSE)
  Fstat <- (ivs$beta.exposure / ivs$se.exposure) ^ 2
  expect_true(all(Fstat >= 50))
})

test_that("output passes assert_exposure() column requirements", {
  fn <- make_clump_fn(mk_hits_tbl())
  mock_ld_stack()
  ivs <- fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(), f_threshold = 0,
            verbose = FALSE)
  expect_silent(assert_exposure(ivs))
})

# ---- new ordering: indel/MAF dropped BEFORE clump (Job 4) ----------------

test_that("indels are dropped BEFORE clumping (regression for new order)", {
  # rs1 is an indel; rs2, rs3 are SNVs. Job 4 order drops indels first
  # so the indel must NEVER reach the clump library.
  hits <- tibble::tibble(
    rsid = c("rs1","rs2","rs3"), seqnames = rep(1L, 3),
    start = c(1000L, 2000L, 3000L),
    REF = c("A","C","G"), ALT = c("AT","T","A"),     # rs1 is an indel
    ES  = c(0.20, 0.05, 0.04), SE = rep(0.01, 3),
    LP  = c(50, 12, 11),
    AF  = rep(0.3, 3), SS = rep(100000L, 3)
  )
  fn <- make_clump_fn(hits)
  observed_clump_input <- new.env()
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) mk_fake_ld_ref(),
    .find_plink_bin = function(plink_bin = NULL) NULL
  )
  testthat::local_mocked_bindings(
    ld_clump = function(dat, ...) {
      observed_clump_input$rsids <- as.character(dat$rsid)
      dat
    },
    .package = "ieugwasr"
  )
  fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
     p_threshold = 1, f_threshold = 0, verbose = FALSE)
  # rs1 (the indel) must not appear in the clump library; only rs2, rs3.
  expect_false("rs1" %in% observed_clump_input$rsids)
  expect_setequal(observed_clump_input$rsids, c("rs2", "rs3"))
})

test_that("low-MAF SNPs are dropped BEFORE clumping when maf_min set", {
  hits <- tibble::tibble(
    rsid = c("rs1","rs2","rs3"), seqnames = rep(1L, 3),
    start = c(1000L, 2000L, 3000L),
    REF = c("A","C","G"), ALT = c("G","T","A"),
    ES  = c(0.20, 0.05, 0.04), SE = rep(0.01, 3),
    LP  = c(50, 12, 11),
    AF  = c(0.005, 0.30, 0.40),  # rs1 is rare (MAF = 0.005)
    SS  = rep(100000L, 3)
  )
  fn <- make_clump_fn(hits)
  observed_clump_input <- new.env()
  testthat::local_mocked_bindings(
    ld_reference_checker = function(...) mk_fake_ld_ref(),
    .find_plink_bin = function(plink_bin = NULL) NULL
  )
  testthat::local_mocked_bindings(
    ld_clump = function(dat, ...) {
      observed_clump_input$rsids <- as.character(dat$rsid)
      dat
    },
    .package = "ieugwasr"
  )
  fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
     p_threshold = 1, f_threshold = 0,
     clump_opts = list(maf_min = 0.05),
     verbose = FALSE)
  # rs1 (low MAF) must not reach the clump library.
  expect_false("rs1" %in% observed_clump_input$rsids)
  expect_setequal(observed_clump_input$rsids, c("rs2", "rs3"))
})

test_that("SNV clump-elected lead survives end-to-end", {
  # Single SNV survives all filters and clumping.
  hits <- tibble::tibble(
    rsid = "rs2", seqnames = 1L, start = 2000L,
    REF = "C", ALT = "T",
    ES  = 0.05, SE = 0.01, LP = 12, AF = 0.30, SS = 100000L
  )
  fn <- make_clump_fn(hits)
  mock_ld_stack(ld_clump_fn = function(dat, ...) dat)
  ivs <- fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
            p_threshold = 1, f_threshold = 0, verbose = FALSE)
  expect_equal(nrow(ivs), 1L)
  expect_equal(ivs$SNP, "rs2")
})

# ---- INFO filter at IV stage ----------------------------------------------

test_that("INFO filter applies when clump_opts$info_min is set and column present", {
  hits <- mk_hits_tbl(n = 3)
  hits$INFO <- c(0.30, 0.95, 0.50)  # rs1 below 0.8, rs2 above, rs3 below 0.8
  fn <- make_clump_fn(hits)
  mock_ld_stack()
  ivs <- fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
            p_threshold = 1, f_threshold = 0,
            clump_opts = list(info_min = 0.8),
            verbose = FALSE)
  expect_equal(nrow(ivs), 1L)
  expect_equal(ivs$SNP, "rs2")
})

test_that("INFO column missing -> info_min has no effect, no error", {
  fn <- make_clump_fn(mk_hits_tbl(n = 3))  # no INFO column
  mock_ld_stack()
  ivs <- fn("fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
            p_threshold = 1, f_threshold = 0,
            clump_opts = list(info_min = 0.8),
            verbose = FALSE)
  expect_equal(nrow(ivs), 3L)
})

# ---- INFO propagation to coloc-region fetcher (resolver-level) -----------

test_that("info_min in clump_opts propagates to .vcf_fetcher closure", {
  observed <- new.env()
  testthat::local_mocked_bindings(
    read_gwasvcf_region = function(vcf_path, chr, start, end,
                                   genome_build = "GRCh37",
                                   drop_indels = TRUE,
                                   info_min = 0.8,
                                   verbose = TRUE) {
      observed$info_min <- info_min
      tibble::tibble(SNP = "rs1", chr = 1L, pos = 100L,
                     effect_allele = "A", other_allele = "G",
                     beta = 0.1, se = 0.01, eaf = 0.3, pval = 1e-9,
                     n = 100000, ncase = NA_real_)
    },
    clump_sumstats_to_ivs = function(...) {
      tibble::tibble(SNP = "rs1", Chr = 1L, Pos = 100L,
                     effect_allele.exposure = "A", other_allele.exposure = "G",
                     beta.exposure = 0.1, se.exposure = 0.01,
                     eaf.exposure = 0.3, pval.exposure = 1e-9,
                     samplesize.exposure = 100000,
                     id.exposure = "x", exposure = "x")
    },
    validate_gwasvcf = function(...) invisible(TRUE),
    read_gwasvcf_metadata = function(...) {
      list(N = 100000, type = "quant", ncase = NA_real_,
           ncontrol = NA_real_, s = NA_real_)
    }
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_sumstats = "fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(),
    clump_opts = list(info_min = 0.5),
    verbose = FALSE
  )
  out <- r$coloc_fetcher(1, 1L, 50000L)
  expect_equal(observed$info_min, 0.5)
})

test_that(".vcf_fetcher uses default info_min = 0.8 when clump_opts$info_min unset", {
  observed <- new.env()
  testthat::local_mocked_bindings(
    read_gwasvcf_region = function(vcf_path, chr, start, end,
                                   genome_build = "GRCh37",
                                   drop_indels = TRUE,
                                   info_min = 0.8,
                                   verbose = TRUE) {
      observed$info_min <- info_min
      tibble::tibble(SNP = "rs1", chr = 1L, pos = 100L,
                     effect_allele = "A", other_allele = "G",
                     beta = 0.1, se = 0.01, eaf = 0.3, pval = 1e-9,
                     n = 100000, ncase = NA_real_)
    },
    clump_sumstats_to_ivs = function(...) {
      tibble::tibble(SNP = "rs1", Chr = 1L, Pos = 100L,
                     effect_allele.exposure = "A", other_allele.exposure = "G",
                     beta.exposure = 0.1, se.exposure = 0.01,
                     eaf.exposure = 0.3, pval.exposure = 1e-9,
                     samplesize.exposure = 100000,
                     id.exposure = "x", exposure = "x")
    },
    validate_gwasvcf = function(...) invisible(TRUE),
    read_gwasvcf_metadata = function(...) {
      list(N = 100000, type = "quant", ncase = NA_real_,
           ncontrol = NA_real_, s = NA_real_)
    }
  )
  r <- ardmr:::.resolve_coloc_source(
    exposure_sumstats = "fake.vcf.bgz", ancestry = "EUR", cache_dir = tempdir(), verbose = FALSE
  )
  out <- r$coloc_fetcher(1, 1L, 50000L)
  expect_equal(observed$info_min, 0.8)
})

# ---- preserved errors -----------------------------------------------------

test_that("missing tabix index errors at validate stage", {
  # Use a real fixture for this test (Rsamtools is installed).
  skip_if_not_installed("Rsamtools")
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  unlink(paste0(vcf, ".tbi"))
  expect_error(
    clump_sumstats_to_ivs(vcf, ancestry = "EUR", verbose = FALSE),
    "Tabix index missing"
  )
})
