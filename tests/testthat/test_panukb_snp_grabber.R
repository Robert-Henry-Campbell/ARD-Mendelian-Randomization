test_that("panukb_snp_grabber skips phenotype when header download fails", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_if_not_installed("glue")

  exposure_snps <- tibble::tibble(
    rsid = "rs123",
    panukb_chrom = 1L,
    panukb_pos = 12345L,
    effect_allele.exposure = "A"
  )

  MR_df <- tibble::tibble(
    aws_link = "https://example.com/data.bgz",
    aws_link_tabix = "https://example.com/data.bgz.tbi",
    description = "Mock phenotype"
  )

  cache_dir <- tempfile("ardmr-cache-")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  warn_calls <- character()

  res <- with_mocked_bindings(
    panukb_snp_grabber(
      exposure_snps = exposure_snps,
      MR_df = MR_df,
      ancestry = "EUR",
      cache_dir = cache_dir,
      verbose = TRUE
    ),
    Rsamtools::TabixFile = function(file, index) list(file = file, index = index),
    Rsamtools::open.TabixFile = function(x) invisible(x),
    Rsamtools::headerTabix = function(tf) character(),
    Rsamtools::scanTabix = function(tf, param) list("chr1\t12345\tA\tG\t0.1\t0.2\t1.3"),
    Rsamtools::close.TabixFile = function(x) invisible(NULL),
    gzcon = function(con) "mock_con",
    url = function(...) "mock_url",
    readLines = function(con, n) {
      warning("mock download failure", call. = FALSE)
      NA_character_
    },
    close = function(con) invisible(NULL),
    digest::digest = function(x, algo = "xxhash64") "digesthash",
    utils::download.file = function(url, destfile, mode, quiet) {
      file.create(destfile)
      0
    },
    logger::log_info = function(msg, ..., .envir = parent.frame()) invisible(NULL),
    logger::log_warn = function(msg, ..., .envir = parent.frame()) {
      warn_calls <<- c(warn_calls, as.character(glue::glue(msg, ..., .envir = .envir)))
      invisible(NULL)
    }
  )

  expect_true(tibble::is_tibble(res$outcome_snps[[1]]))
  expect_identical(nrow(res$outcome_snps[[1]]), 0L)
  expect_true(any(grepl("unable to download header", warn_calls, fixed = TRUE)))
  expect_true(any(grepl("mock download failure", warn_calls, fixed = TRUE)))
})
