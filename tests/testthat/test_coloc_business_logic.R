# tests/testthat/test_coloc_business_logic.R
#
# Offline tests for coloc_business_logic(). External calls are bypassed by
# passing exposure_region_fetcher / outcome_region_fetcher closures and
# the exposure_metadata / outcome_metadata lists directly.

skip_if_not_installed("coloc")

cache_dir_test <- function() {
  d <- file.path(tempdir(), "ardmr_coloc_test")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

QUANT_META <- function(N) list(N = N, type = "quant",
                               ncase = NA_real_, ncontrol = NA_real_, s = NA_real_)

test_that("locus merging: two IVs within window collapse to one locus", {
  iv_pos <- tibble::tibble(SNP = c("rs1","rs2"), chr = c("1","1"), pos = c(50000L, 60000L))
  loci <- ardmr:::.coloc_build_loci(iv_pos, window_kb = 100L)
  expect_equal(nrow(loci), 1)
  expect_equal(loci$chr[1], "1")
  expect_lt(loci$start[1], 50000L)
  expect_gt(loci$end[1], 60000L)
})

test_that("locus merging: two distant IVs yield two loci", {
  iv_pos <- tibble::tibble(SNP = c("rs1","rs2"), chr = c("1","1"), pos = c(50000L, 5000000L))
  loci <- ardmr:::.coloc_build_loci(iv_pos, window_kb = 100L)
  expect_equal(nrow(loci), 2)
})

test_that("locus merging: separate chromosomes yield separate loci", {
  iv_pos <- tibble::tibble(SNP = c("rs1","rs2"), chr = c("1","2"), pos = c(50000L, 50000L))
  loci <- ardmr:::.coloc_build_loci(iv_pos, window_kb = 100L)
  expect_equal(nrow(loci), 2)
})

test_that(".coloc_iv_positions returns pval column from pval.exposure", {
  hdat <- make_coloc_synthetic_iv_hdat("rs1", 1L, 100L, 0.08, 0.01, 0.06, 0.01)
  iv_pos <- ardmr:::.coloc_iv_positions(hdat, tibble::tibble(SNP = character()))
  expect_true("pval" %in% names(iv_pos))
  expect_equal(iv_pos$pval[1], hdat$pval.exposure[1])
})

test_that(".coloc_build_loci attaches min_pval_exp per locus", {
  iv_pos <- tibble::tibble(
    SNP  = c("rs1","rs2","rs3"),
    chr  = c("1","1","2"),
    pos  = c(50000L, 60000L, 50000L),
    pval = c(1e-12, 1e-5, 1e-8)
  )
  loci <- ardmr:::.coloc_build_loci(iv_pos, window_kb = 100L)
  expect_true("min_pval_exp" %in% names(loci))
  loc1 <- loci[loci$chr == "1", , drop = FALSE]
  loc2 <- loci[loci$chr == "2", , drop = FALSE]
  expect_equal(loc1$min_pval_exp, 1e-12)
  expect_equal(loc2$min_pval_exp, 1e-8)
})

test_that("harmonisation flips beta when alleles swap", {
  exp_region <- tibble::tibble(
    SNP = c("rs1","rs2"), chr = c(1L,1L), pos = c(100L,200L),
    effect_allele = c("A","A"), other_allele = c("G","G"),
    beta = c(0.1, 0.05), se = c(0.01, 0.01),
    eaf = c(0.3, 0.3), pval = c(1e-10, 1e-3)
  )
  out_region <- tibble::tibble(
    SNP = c("a","b"), chr = c(1L,1L), pos = c(100L,200L),
    effect_allele = c("G","A"), other_allele = c("A","G"),
    beta = c(0.08, 0.04), se = c(0.01, 0.01),
    eaf = c(0.7, 0.3), pval = c(1e-8, 1e-3)
  )
  harm <- ardmr:::.coloc_harmonise(exp_region, out_region)
  expect_equal(nrow(harm), 2)
  expect_equal(harm$beta_out[harm$pos == 100L], -0.08, tolerance = 1e-9)
  expect_equal(harm$beta_out[harm$pos == 200L],  0.04, tolerance = 1e-9)
})

test_that("harmonisation drops mismatched-allele SNPs", {
  exp_region <- tibble::tibble(
    SNP = "rs1", chr = 1L, pos = 100L,
    effect_allele = "A", other_allele = "G",
    beta = 0.1, se = 0.01, eaf = 0.3, pval = 1e-10
  )
  out_region <- tibble::tibble(
    SNP = "x", chr = 1L, pos = 100L,
    effect_allele = "C", other_allele = "T",
    beta = 0.05, se = 0.01, eaf = 0.3, pval = 1e-3
  )
  harm <- ardmr:::.coloc_harmonise(exp_region, out_region)
  expect_equal(nrow(harm), 0)
})

test_that("coloc_business_logic: shared causal variant yields PP.H4 > 0.80", {
  exp_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, N = 100000L, seed = 1L)
  out_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 6, N = 50000L, seed = 1L)
  iv_snp <- exp_region$SNP[which.min(abs(exp_region$pos - 100000L))]
  iv_pos <- exp_region$pos[which.min(abs(exp_region$pos - 100000L))]
  hdat_use <- make_coloc_synthetic_iv_hdat(iv_snp, 1L, iv_pos, 0.08, 0.01, 0.06, 0.01)
  exposure_snps <- make_coloc_synthetic_exposure_snps(iv_snp, 1L, iv_pos)
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = function(chr, start, end) exp_region,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = function(chr, start, end) out_region,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = "", cache_dir = cache_dir_test(), verbose = FALSE
  )
  expect_equal(nrow(res), 1)
  expect_gt(res$pp_h4_abf[1], 0.80)
  expect_equal(res$method_passed[1], "abf")
})

test_that("coloc_business_logic: no shared signal yields low PP.H4", {
  exp_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, N = 100000L, seed = 2L)
  out_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L,  30000L, 8, N =  50000L, seed = 3L)
  iv_snp <- exp_region$SNP[which.min(abs(exp_region$pos - 100000L))]
  iv_pos <- exp_region$pos[which.min(abs(exp_region$pos - 100000L))]
  hdat_use <- make_coloc_synthetic_iv_hdat(iv_snp, 1L, iv_pos, 0.08, 0.01, 0.005, 0.01)
  exposure_snps <- make_coloc_synthetic_exposure_snps(iv_snp, 1L, iv_pos)
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = function(chr, start, end) exp_region,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = function(chr, start, end) out_region,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = "", cache_dir = cache_dir_test(), verbose = FALSE
  )
  expect_equal(nrow(res), 1)
  expect_lt(res$pp_h4_abf[1], 0.80)
  expect_equal(res$method_passed[1], "none")
})

test_that("coloc_business_logic: aggregation rule passes if ANY locus PP.H4 > 0.80", {
  shared_exp <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, seed = 4L)
  shared_out <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 6, seed = 4L)
  null_exp   <- make_coloc_synthetic_region(2L, 1L, 200000L, 60L, 100000L, 8, seed = 5L)
  null_out   <- make_coloc_synthetic_region(2L, 1L, 200000L, 60L,  30000L, 8, seed = 6L)
  iv1 <- shared_exp$SNP[which.min(abs(shared_exp$pos - 100000L))]
  iv2 <-   null_exp$SNP[which.min(abs(null_exp$pos   - 100000L))]
  hdat_use <- dplyr::bind_rows(
    make_coloc_synthetic_iv_hdat(iv1, 1L, 100000L, 0.08, 0.01, 0.06, 0.01),
    make_coloc_synthetic_iv_hdat(iv2, 2L, 100000L, 0.08, 0.01, 0.005, 0.01)
  )
  exposure_snps <- make_coloc_synthetic_exposure_snps(c(iv1, iv2), c(1L, 2L), c(100000L, 100000L))
  exp_fetcher <- function(chr, start, end) if (as.character(chr) == "1") shared_exp else null_exp
  out_fetcher <- function(chr, start, end) if (as.character(chr) == "1") shared_out else null_out
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = exp_fetcher,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = out_fetcher,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = "", cache_dir = cache_dir_test(), verbose = FALSE
  )
  expect_equal(nrow(res), 2)
  expect_true(any(res$pp_h4_abf > 0.80))
  expect_true(any(res$pp_h4_abf < 0.80))
})

test_that("coloc_business_logic: empty IV positions returns empty tibble", {
  res <- coloc_business_logic(
    hdat_use = tibble::tibble(),
    exposure_snps = tibble::tibble(SNP = character()),
    exposure_region_fetcher = function(chr, start, end) tibble::tibble(),
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = function(chr, start, end) tibble::tibble(),
    outcome_metadata = QUANT_META(50000L),
    window_kb = 100L, plot_dir = "", cache_dir = cache_dir_test(), verbose = FALSE
  )
  expect_equal(nrow(res), 0)
})

test_that("coloc artefacts are written to plot_dir", {
  exp_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, seed = 7L)
  out_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 6, seed = 7L)
  iv_snp <- exp_region$SNP[which.min(abs(exp_region$pos - 100000L))]
  hdat_use <- make_coloc_synthetic_iv_hdat(iv_snp, 1L, 100000L, 0.08, 0.01, 0.06, 0.01)
  exposure_snps <- make_coloc_synthetic_exposure_snps(iv_snp, 1L, 100000L)
  plot_root <- file.path(tempdir(), paste0("ardmr_coloc_plots_", as.integer(Sys.time())))
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = function(chr, start, end) exp_region,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = function(chr, start, end) out_region,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = plot_root, cache_dir = cache_dir_test(), verbose = FALSE
  )
  expect_equal(nrow(res), 1)
  locus_dir <- res$coloc_dir[1]
  expect_true(file.exists(file.path(locus_dir, "coloc_abf_summary.csv")))
  expect_true(file.exists(file.path(locus_dir, "coloc_aligned_data.csv")))
  expect_true(file.exists(file.path(locus_dir, "stacked_manhattan.png")) ||
              file.exists(file.path(locus_dir, "stacked_manhattan.pdf")))
})

test_that("single locus folder is prefixed with 1-LEAD-", {
  exp_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, N = 100000L, seed = 1L)
  out_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 6, N = 50000L, seed = 1L)
  iv_snp <- exp_region$SNP[which.min(abs(exp_region$pos - 100000L))]
  iv_p   <- exp_region$pos[which.min(abs(exp_region$pos - 100000L))]
  hdat_use <- make_coloc_synthetic_iv_hdat(iv_snp, 1L, iv_p, 0.08, 0.01, 0.06, 0.01)
  exposure_snps <- make_coloc_synthetic_exposure_snps(iv_snp, 1L, iv_p)
  plot_root <- file.path(tempdir(), paste0("ardmr_coloc_lead1_", as.integer(Sys.time())))
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = function(chr, start, end) exp_region,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = function(chr, start, end) out_region,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = plot_root,
    cache_dir = cache_dir_test(), verbose = FALSE
  )
  expect_equal(nrow(res), 1)
  dirs <- list.dirs(plot_root, recursive = FALSE, full.names = FALSE)
  expect_equal(length(dirs), 1L)
  expect_match(dirs[1], "^1-LEAD-chr1_")
  expect_false(grepl("^1-LEAD-", res$locus_id[1]))
  expect_match(basename(res$coloc_dir[1]), "^1-LEAD-")
})

test_that("multi-locus run prefixes only the lead-IV locus folder", {
  shared_exp <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, seed = 4L)
  shared_out <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 6, seed = 4L)
  null_exp   <- make_coloc_synthetic_region(2L, 1L, 200000L, 60L, 100000L, 8, seed = 5L)
  null_out   <- make_coloc_synthetic_region(2L, 1L, 200000L, 60L,  30000L, 8, seed = 6L)
  iv1 <- shared_exp$SNP[which.min(abs(shared_exp$pos - 100000L))]
  iv2 <-   null_exp$SNP[which.min(abs(null_exp$pos   - 100000L))]
  hdat_use <- dplyr::bind_rows(
    make_coloc_synthetic_iv_hdat(iv1, 1L, 100000L, 0.20, 0.01, 0.06, 0.01),
    make_coloc_synthetic_iv_hdat(iv2, 2L, 100000L, 0.05, 0.02, 0.005, 0.01)
  )
  exposure_snps <- make_coloc_synthetic_exposure_snps(c(iv1, iv2), c(1L, 2L),
                                                     c(100000L, 100000L))
  exp_fetcher <- function(chr, start, end) if (as.character(chr) == "1") shared_exp else null_exp
  out_fetcher <- function(chr, start, end) if (as.character(chr) == "1") shared_out else null_out
  plot_root <- file.path(tempdir(), paste0("ardmr_coloc_lead2_", as.integer(Sys.time())))
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = exp_fetcher,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = out_fetcher,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = plot_root,
    cache_dir = cache_dir_test(), verbose = FALSE
  )
  expect_equal(nrow(res), 2)
  dirs <- list.dirs(plot_root, recursive = FALSE, full.names = FALSE)
  lead_dirs <- grep("^1-LEAD-", dirs, value = TRUE)
  expect_equal(length(lead_dirs), 1L)
  expect_match(lead_dirs[1], "^1-LEAD-chr1_")
  expect_false(any(grepl("^1-LEAD-", res$locus_id)))
  lead_row  <- res[res$chr == "1", , drop = FALSE]
  other_row <- res[res$chr == "2", , drop = FALSE]
  expect_match(basename(lead_row$coloc_dir[1]), "^1-LEAD-")
  expect_false(grepl("^1-LEAD-", basename(other_row$coloc_dir[1])))
})

test_that("plot_dir = '' leaves locus_id unprefixed and coloc_dir NA", {
  exp_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, seed = 1L)
  out_region <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 6, seed = 1L)
  iv_snp <- exp_region$SNP[which.min(abs(exp_region$pos - 100000L))]
  iv_p   <- exp_region$pos[which.min(abs(exp_region$pos - 100000L))]
  hdat_use <- make_coloc_synthetic_iv_hdat(iv_snp, 1L, iv_p, 0.08, 0.01, 0.06, 0.01)
  exposure_snps <- make_coloc_synthetic_exposure_snps(iv_snp, 1L, iv_p)
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = function(chr, start, end) exp_region,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = function(chr, start, end) out_region,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = "",
    cache_dir = cache_dir_test(), verbose = FALSE
  )
  expect_equal(nrow(res), 1)
  expect_false(grepl("^1-LEAD-", res$locus_id[1]))
  expect_true(is.na(res$coloc_dir[1]))
})

test_that("MHC removal reassigns 1-LEAD- to next-lowest-p-value surviving locus", {
  iv_mhc  <- "rs_mhc"
  iv_chr1 <- "rs_chr1"
  iv_chr3 <- "rs_chr3"
  hdat_use <- dplyr::bind_rows(
    make_coloc_synthetic_iv_hdat(iv_mhc,  6L, 30000000L, 0.30, 0.01, 0.06, 0.01),
    make_coloc_synthetic_iv_hdat(iv_chr1, 1L,   100000L, 0.20, 0.01, 0.06, 0.01),
    make_coloc_synthetic_iv_hdat(iv_chr3, 3L,   100000L, 0.05, 0.02, 0.06, 0.01)
  )
  exposure_snps <- make_coloc_synthetic_exposure_snps(
    c(iv_mhc, iv_chr1, iv_chr3), c(6L, 1L, 3L),
    c(30000000L, 100000L, 100000L)
  )
  region_chr1 <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, seed = 11L)
  region_chr3 <- make_coloc_synthetic_region(3L, 1L, 200000L, 60L, 100000L, 8, seed = 12L)
  out_chr1    <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 6, seed = 11L)
  out_chr3    <- make_coloc_synthetic_region(3L, 1L, 200000L, 60L,  30000L, 8, seed = 13L)
  exp_fetcher <- function(chr, start, end) if (as.character(chr) == "1") region_chr1 else region_chr3
  out_fetcher <- function(chr, start, end) if (as.character(chr) == "1") out_chr1    else out_chr3
  plot_root <- file.path(tempdir(), paste0("ardmr_coloc_mhc_lead_", as.integer(Sys.time())))
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = exp_fetcher,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = out_fetcher,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = plot_root,
    cache_dir = cache_dir_test(), verbose = FALSE,
    skip_mhc = TRUE
  )
  expect_equal(nrow(res), 2)
  dirs <- list.dirs(plot_root, recursive = FALSE, full.names = FALSE)
  lead_dirs <- grep("^1-LEAD-", dirs, value = TRUE)
  expect_equal(length(lead_dirs), 1L)
  expect_match(lead_dirs[1], "^1-LEAD-chr1_")
})

.simple_console_appender <- function(lines) {
  cat(lines, sep = "\n", file = stderr())
}

.run_mr_with_coloc <- function(coloc_opts) {
  inputs <- make_synthetic_mr_inputs()
  captured <- character()
  logger::log_appender(function(lines) captured <<- c(captured, lines))
  on.exit(logger::log_appender(.simple_console_appender), add = TRUE)
  set.seed(SYNTH_SEED)
  out <- mr_business_logic(
    MR_df = inputs$MR_df, exposure_snps = inputs$exposure_snps,
    sensitivity_enabled = c(SYNTH_SENSITIVITY_ENABLED, "coloc"),
    sensitivity_pass_min = SYNTH_SENSITIVITY_PASS_MIN,
    scatterplot = FALSE, snpforestplot = FALSE, leaveoneoutplot = FALSE,
    plot_output_dir = "", cache_dir = tempdir(),
    coloc_opts = coloc_opts, verbose = TRUE
  )
  list(out = out, log = paste(captured, collapse = "\n"))
}

test_that("mr_business_logic warns when outcome N is NA (no silent skip)", {
  r <- .run_mr_with_coloc(list(
    exposure_region_fetcher = function(chr, start, end) tibble::tibble(),
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_fetcher_factory = function(rec) function(chr, start, end) tibble::tibble(),
    outcome_metadata_fn     = function(rec) {
      list(N = NA_real_, type = "quant", ncase = NA_real_, ncontrol = NA_real_, s = NA_real_)
    }
  ))
  expect_match(r$log, "outcome N unknown")
  expect_true(all(is.na(r$out$MR_df$results_coloc_max_PPH4)))
})

test_that("mr_business_logic warns when outcome fetcher returns NULL", {
  r <- .run_mr_with_coloc(list(
    exposure_region_fetcher = function(chr, start, end) tibble::tibble(),
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_fetcher_factory = function(rec) NULL,
    outcome_metadata_fn     = function(rec) QUANT_META(50000L)
  ))
  expect_match(r$log, "outcome fetcher unavailable")
})

test_that("coloc gating: non-significant outcomes skipped, first analysed outcome runs as smoke test", {
  exposure_snps <- make_synthetic_exposure()
  # Three outcomes, all designed to give VERY non-significant IVW (large SEs).
  out1 <- make_synthetic_outcome(
    snps     = c("rs1", "rs2", "rs3"),
    beta_out = c(0.001, -0.001, 0.001),
    se_out   = c(0.50,  0.50,   0.50),
    label    = "OutNullA"
  )
  out2 <- make_synthetic_outcome(
    snps     = c("rs1", "rs2", "rs3"),
    beta_out = c(0.0005, -0.0005, 0.0005),
    se_out   = c(0.50,   0.50,    0.50),
    label    = "OutNullB"
  )
  out3 <- make_synthetic_outcome(
    snps     = c("rs1", "rs2", "rs3"),
    beta_out = c(0.0008, -0.0009, 0.0007),
    se_out   = c(0.50,   0.50,    0.50),
    label    = "OutNullC"
  )
  MR_df <- tibble::tibble(
    description  = c("OutNullA", "OutNullB", "OutNullC"),
    outcome_snps = list(out1, out2, out3)
  )
  call_log <- character()
  testthat::local_mocked_bindings(
    coloc_business_logic = function(hdat_use, ...) {
      out_label <- if ("outcome" %in% names(hdat_use)) hdat_use$outcome[1] else "?"
      call_log <<- c(call_log, as.character(out_label))
      tibble::tibble(
        locus_id = "chr1_x", chr = "1", start = 1L, end = 1L,
        n_ivs_in_locus = 1L, n_snps_used = 5L,
        pp_h4_abf = 0.5, pp_h4_susie_max = NA_real_,
        method_passed = "none", coloc_dir = NA_character_
      )
    }
  )
  set.seed(SYNTH_SEED)
  out <- mr_business_logic(
    MR_df = MR_df, exposure_snps = exposure_snps,
    sensitivity_enabled = c(SYNTH_SENSITIVITY_ENABLED, "coloc"),
    sensitivity_pass_min = SYNTH_SENSITIVITY_PASS_MIN,
    scatterplot = FALSE, snpforestplot = FALSE, leaveoneoutplot = FALSE,
    plot_output_dir = "", cache_dir = tempdir(),
    coloc_opts = list(
      exposure_region_fetcher = function(chr, start, end) tibble::tibble(),
      exposure_metadata = QUANT_META(100000L),
      ancestry = "EUR",
      outcome_fetcher_factory = function(rec) function(chr, start, end) tibble::tibble(),
      outcome_metadata_fn     = function(rec) QUANT_META(50000L)
    ),
    verbose = FALSE
  )
  expect_equal(length(call_log), 1L)
  expect_equal(call_log[1], "OutNullA")
  # Smoke-test outcome has coloc-tbl populated; later outcomes do not.
  expect_false(is.na(out$MR_df$results_coloc_max_PPH4[1]))
  expect_true(is.na(out$MR_df$results_coloc_max_PPH4[2]))
  expect_true(is.na(out$MR_df$results_coloc_max_PPH4[3]))
})

test_that("coloc gating: nominally significant outcomes also trigger coloc", {
  exposure_snps <- make_synthetic_exposure()
  out_null <- make_synthetic_outcome(
    snps     = c("rs1", "rs2", "rs3"),
    beta_out = c(0.001, -0.001, 0.001),
    se_out   = c(0.50,  0.50,   0.50),
    label    = "OutNull"
  )
  out_sig <- make_synthetic_outcome(
    snps     = c("rs1", "rs2", "rs3"),
    beta_out = c(0.080, -0.060, 0.150),
    se_out   = c(0.012, 0.018,  0.010),
    label    = "OutSig"
  )
  MR_df <- tibble::tibble(
    description  = c("OutNull", "OutSig"),
    outcome_snps = list(out_null, out_sig)
  )
  call_log <- character()
  testthat::local_mocked_bindings(
    coloc_business_logic = function(hdat_use, ...) {
      out_label <- if ("outcome" %in% names(hdat_use)) hdat_use$outcome[1] else "?"
      call_log <<- c(call_log, as.character(out_label))
      tibble::tibble(
        locus_id = "chr1_x", chr = "1", start = 1L, end = 1L,
        n_ivs_in_locus = 1L, n_snps_used = 5L,
        pp_h4_abf = 0.5, pp_h4_susie_max = NA_real_,
        method_passed = "none", coloc_dir = NA_character_
      )
    }
  )
  set.seed(SYNTH_SEED)
  out <- mr_business_logic(
    MR_df = MR_df, exposure_snps = exposure_snps,
    sensitivity_enabled = c(SYNTH_SENSITIVITY_ENABLED, "coloc"),
    sensitivity_pass_min = SYNTH_SENSITIVITY_PASS_MIN,
    scatterplot = FALSE, snpforestplot = FALSE, leaveoneoutplot = FALSE,
    plot_output_dir = "", cache_dir = tempdir(),
    coloc_opts = list(
      exposure_region_fetcher = function(chr, start, end) tibble::tibble(),
      exposure_metadata = QUANT_META(100000L),
      ancestry = "EUR",
      outcome_fetcher_factory = function(rec) function(chr, start, end) tibble::tibble(),
      outcome_metadata_fn     = function(rec) QUANT_META(50000L)
    ),
    verbose = FALSE
  )
  # Both should run: first (smoke) + significant
  expect_equal(length(call_log), 2L)
  expect_setequal(call_log, c("OutNull", "OutSig"))
})

test_that("skip_mhc drops chr6:25-35Mb loci by default", {
  iv_snp_mhc <- "rs_mhc"; iv_snp_chr1 <- "rs_chr1"
  hdat_use <- dplyr::bind_rows(
    make_coloc_synthetic_iv_hdat(iv_snp_mhc, 6L, 30000000L, 0.08, 0.01, 0.06, 0.01),
    make_coloc_synthetic_iv_hdat(iv_snp_chr1, 1L, 100000L, 0.08, 0.01, 0.06, 0.01)
  )
  exposure_snps <- make_coloc_synthetic_exposure_snps(
    c(iv_snp_mhc, iv_snp_chr1), c(6L, 1L), c(30000000L, 100000L)
  )
  exp_chr1 <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 8, seed = 8L)
  out_chr1 <- make_coloc_synthetic_region(1L, 1L, 200000L, 60L, 100000L, 6, seed = 8L)
  exp_fetcher <- function(chr, start, end) if (as.character(chr) == "1") exp_chr1 else exp_chr1
  out_fetcher <- function(chr, start, end) if (as.character(chr) == "1") out_chr1 else out_chr1
  res <- coloc_business_logic(
    hdat_use = hdat_use, exposure_snps = exposure_snps,
    exposure_region_fetcher = exp_fetcher,
    exposure_metadata = QUANT_META(100000L),
    ancestry = "EUR",
    outcome_region_fetcher = out_fetcher,
    outcome_metadata = QUANT_META(50000L),
    window_kb = 200L, plot_dir = "", cache_dir = cache_dir_test(), verbose = FALSE,
    skip_mhc = TRUE
  )
  expect_equal(nrow(res), 1)
  expect_equal(res$chr[1], "1")
})
