# tests/testthat/test_coloc_vcf_reader.R

skip_if_not_installed("Rsamtools")

test_that("validate_gwasvcf accepts a well-formed GWAS-VCF", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  expect_true(validate_gwasvcf(vcf, genome_build = "GRCh37"))
})

test_that("validate_gwasvcf rejects GRCh38 declared build", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  expect_error(validate_gwasvcf(vcf, genome_build = "GRCh38"),
               "only supports GRCh37")
})

test_that("validate_gwasvcf errors when tabix index is missing", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  unlink(paste0(vcf, ".tbi"))
  expect_error(validate_gwasvcf(vcf, genome_build = "GRCh37"),
               "Tabix index missing")
})

test_that("validate_gwasvcf errors on missing required FORMAT fields", {
  dir <- tempfile("badvcf_"); dir.create(dir, recursive = TRUE)
  bad <- file.path(dir, "bad.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"af\">",
    "##contig=<ID=1>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tstudy1",
    "1\t1000\trs1\tA\tG\t.\t.\t.\tAF\t0.3"
  ), bad)
  bgz <- Rsamtools::bgzip(bad, dest = paste0(bad, ".bgz"), overwrite = TRUE)
  Rsamtools::indexTabix(bgz, format = "vcf")
  expect_error(validate_gwasvcf(bgz, genome_build = "GRCh37"),
               "missing required FORMAT")
})

test_that("read_gwasvcf_region returns canonical schema", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  r <- read_gwasvcf_region(vcf, chr = 1, start = 1L, end = 50000L,
                           genome_build = "GRCh37", verbose = FALSE)
  expect_named(r,
               c("SNP","chr","pos","effect_allele","other_allele",
                 "beta","se","eaf","pval","n","ncase"),
               ignore.order = TRUE)
  expect_true(all(r$chr == 1L))
  expect_true(all(is.finite(r$beta) & is.finite(r$se) & r$se > 0))
})

test_that("read_gwasvcf_region drops indels by default", {
  vcf <- make_gwasvcf_fixture(include_indel = TRUE, include_multiallelic = FALSE)
  r <- read_gwasvcf_region(vcf, chr = 1, start = 1L, end = 1000000L,
                           genome_build = "GRCh37", verbose = FALSE)
  expect_false("rs_ind" %in% r$SNP)
  expect_true(all(nchar(r$effect_allele) == 1L))
  expect_true(all(nchar(r$other_allele) == 1L))
})

test_that("read_gwasvcf_region splits multiallelic rows", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = TRUE)
  r <- read_gwasvcf_region(vcf, chr = 1, start = 1L, end = 50000L,
                           genome_build = "GRCh37", verbose = FALSE)
  multi <- r[r$SNP == "rs_multi", , drop = FALSE]
  expect_equal(nrow(multi), 2)
  expect_setequal(multi$effect_allele, c("A", "C"))
  expect_setequal(round(multi$beta, 4), c(0.05, 0.08))
})

test_that("read_gwasvcf_region tabix-restricts to the requested window", {
  vcf <- make_gwasvcf_fixture(hits_per_chr = 60L, include_indel = FALSE, include_multiallelic = FALSE)
  r <- read_gwasvcf_region(vcf, chr = 1, start = 5000L, end = 12000L,
                           genome_build = "GRCh37", verbose = FALSE)
  expect_true(all(r$pos >= 5000L & r$pos <= 12000L))
  expect_true(all(r$chr == 1L))
})

test_that("read_gwasvcf_metadata returns quant when no NC field", {
  vcf <- make_gwasvcf_fixture(include_NC = FALSE,
                              include_indel = FALSE, include_multiallelic = FALSE)
  m <- read_gwasvcf_metadata(vcf, genome_build = "GRCh37")
  expect_equal(m$type, "quant")
  expect_equal(m$N, 100000)
  expect_true(is.na(m$ncase))
})

test_that("read_gwasvcf_metadata flags cc when NC present", {
  vcf <- make_gwasvcf_fixture(include_NC = TRUE,
                              include_indel = FALSE, include_multiallelic = FALSE)
  m <- read_gwasvcf_metadata(vcf, genome_build = "GRCh37")
  expect_equal(m$type, "cc")
  expect_equal(m$ncase, 5000)
  expect_equal(m$ncontrol, 95000)
  expect_equal(m$s, 0.05, tolerance = 1e-9)
})

test_that("read_gwasvcf_region fills per-SNP N from SS", {
  vcf <- make_gwasvcf_fixture(include_indel = FALSE, include_multiallelic = FALSE)
  r <- read_gwasvcf_region(vcf, chr = 1, start = 1L, end = 50000L,
                           genome_build = "GRCh37", verbose = FALSE)
  expect_true(all(r$n == 100000))
})
