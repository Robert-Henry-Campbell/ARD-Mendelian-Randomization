# tests/testthat/test_run_hash.R
# Verifies run_phenome_mr writes outputs into a per-run hash subfolder and
# emits run_manifest.json. Heavy stages are mocked inline so the test stays
# offline.

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

run_with_mocks <- function(args, snps = NULL) {
  cache_dir <- tempfile("ardmr_run_")
  dir.create(cache_dir, recursive = TRUE)
  args$cache_dir <- cache_dir
  args$exposure_snps <- if (is.null(snps)) mk_synth_exposure_snps() else snps
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
    # Snps_id mode now routes through preprocess_exposure_snps (Job 2a).
    # Stub it to pass-through so the test stays offline (no LD reference).
    preprocess_exposure_snps = function(snps_tbl, clump_opts = list(),
                                        ancestry, cache_dir, confirm = "ask",
                                        verbose = TRUE) {
      list(snps = tibble::as_tibble(snps_tbl), steps = list(),
           resolved_opts = ardmr:::.resolve_clump_opts(clump_opts))
    },
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
