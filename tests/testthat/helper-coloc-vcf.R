# tests/testthat/helper-coloc-vcf.R
# Build a tiny GWAS-VCF on disk (bgzip+tabix) for the VCF reader tests.

make_gwasvcf_fixture <- function(dir = tempfile("vcf_"),
                                 include_indel = TRUE,
                                 include_multiallelic = TRUE,
                                 include_NC = FALSE,
                                 chrs = c("1","2"),
                                 hits_per_chr = 30L) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  vcf_path <- file.path(dir, "mini.vcf")
  fmt_keys <- c("ES","SE","LP","AF","SS")
  if (include_NC) fmt_keys <- c(fmt_keys, "NC")
  fmt_str <- paste(fmt_keys, collapse = ":")
  hdr <- c(
    "##fileformat=VCFv4.2",
    paste0("##FORMAT=<ID=", fmt_keys,
           ',Number=A,Type=Float,Description="', fmt_keys, '">'),
    paste0("##contig=<ID=", chrs, ">"),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tstudy1"
  )
  rows <- character()
  set.seed(1L)
  for (ch in chrs) {
    for (i in seq_len(hits_per_chr)) {
      pos <- as.integer(i * 1000L)
      ea  <- c("A","C","G","T")[(i %% 4L) + 1L]
      oa  <- c("G","T","A","C")[(i %% 4L) + 1L]
      es  <- round(stats::rnorm(1L, sd = 0.05), 4L)
      se  <- 0.01
      lp  <- abs(round(es / se, 3L))
      af  <- 0.3
      ss  <- 100000L
      vals <- c(es, se, lp, af, ss)
      if (include_NC) vals <- c(vals, 5000L)
      sample_str <- paste(vals, collapse = ":")
      rows <- c(rows, paste(ch, pos, paste0("rs_", ch, "_", pos), ea, oa,
                            ".", ".", ".", fmt_str, sample_str, sep = "\t"))
    }
  }
  if (include_indel) {
    rows <- c(rows, paste("1", "999999", "rs_ind", "A", "AT",
                          ".", ".", ".", fmt_str,
                          if (include_NC) "0.2:0.05:5:0.1:100000:5000" else "0.2:0.05:5:0.1:100000",
                          sep = "\t"))
  }
  if (include_multiallelic) {
    vals_ma <- if (include_NC) "0.05,0.08:0.01,0.01:5,7:0.2,0.1:100000,100000:5000,5000"
               else            "0.05,0.08:0.01,0.01:5,7:0.2,0.1:100000,100000"
    rows <- c(rows, paste("1", "5000", "rs_multi", "G", "A,C",
                          ".", ".", ".", fmt_str, vals_ma, sep = "\t"))
  }
  rows <- rows[order(
    suppressWarnings(as.integer(sub("\t.*$", "", rows))),
    suppressWarnings(as.integer(sub("^[^\t]+\t([0-9]+)\t.*$", "\\1", rows)))
  )]
  writeLines(c(hdr, rows), vcf_path)
  bgz <- Rsamtools::bgzip(vcf_path, dest = paste0(vcf_path, ".bgz"), overwrite = TRUE)
  Rsamtools::indexTabix(bgz, format = "vcf")
  bgz
}
