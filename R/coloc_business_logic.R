# R/coloc_business_logic.R
# Colocalization (coloc.abf + coloc.susie) as a robustness check.

#' Run coloc.abf (+ coloc.susie when LD computable) per IV-region.
#'
#' @param hdat_use Harmonised exposure-outcome tibble (post-`harmonise_data`,
#'   `mr_keep` filtered).
#' @param exposure_snps Original exposure SNPs tibble; used as a fallback for
#'   IV positions if `hdat_use` lacks `chr.outcome`/`pos.outcome`.
#' @param exposure_region_fetcher `function(chr, start, end)` returning a
#'   tibble of exposure SNPs with `chr, pos, effect_allele, other_allele,
#'   beta, se, eaf, pval`.
#' @param exposure_metadata `list(N, type, ncase, ncontrol, s)`.
#' @param ancestry Ancestry code (1000G + Pan-UKB style, e.g. `"EUR"`).
#' @param outcome_region_fetcher Same signature as `exposure_region_fetcher`.
#' @param outcome_metadata Same shape as `exposure_metadata`.
#' @param window_kb Half-window around each IV (default 500).
#' @param plot_dir Directory for per-locus artefacts (`""` = off).
#' @param cache_dir Cache directory.
#' @param priors Coloc priors as `list(p1, p2, p12)`.
#' @param skip_mhc Skip loci overlapping chr6:25Mb-35Mb (default `TRUE`).
#' @param susie_always Try coloc.susie whenever LD is computable, regardless
#'   of IV count in the locus (default `TRUE`).
#' @param verbose Logical.
#' @return Tibble (one row per locus) with `locus_id, chr, start, end,
#'   n_ivs_in_locus, n_snps_used, pp_h4_abf, pp_h4_susie_max, method_passed,
#'   coloc_dir`.
#' @export
coloc_business_logic <- function(hdat_use,
                                 exposure_snps,
                                 exposure_region_fetcher,
                                 exposure_metadata,
                                 ancestry,
                                 outcome_region_fetcher,
                                 outcome_metadata,
                                 window_kb = 500L,
                                 plot_dir = "",
                                 cache_dir = ardmr_cache_dir(),
                                 priors = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
                                 skip_mhc = TRUE,
                                 susie_always = TRUE,
                                 verbose = TRUE) {
  if (!requireNamespace("coloc", quietly = TRUE)) {
    stop("coloc_business_logic(): requires the coloc package.", call. = FALSE)
  }
  if (is.null(exposure_region_fetcher) || is.null(exposure_metadata)) {
    if (verbose) logger::log_warn("coloc: exposure fetcher/metadata missing; skipping")
    return(.coloc_empty_result_tbl())
  }

  iv_pos <- .coloc_iv_positions(hdat_use, exposure_snps)
  if (nrow(iv_pos) == 0) {
    if (verbose) logger::log_warn("coloc: no IV positions available; skipping")
    return(.coloc_empty_result_tbl())
  }

  loci <- .coloc_build_loci(iv_pos, window_kb)
  if (skip_mhc) {
    mhc <- loci$chr == "6" & loci$start <= 35000000L & loci$end >= 25000000L
    n_mhc <- sum(mhc, na.rm = TRUE)
    if (n_mhc) {
      if (verbose) logger::log_info("coloc: skipping {n_mhc} MHC locus/loci (chr6:25-35Mb)")
      loci <- loci[!mhc, , drop = FALSE]
    }
  }
  if (verbose) logger::log_info("coloc: {nrow(loci)} merged loci across {nrow(iv_pos)} IVs (window={window_kb} kb)")
  if (!nrow(loci)) return(.coloc_empty_result_tbl())

  exp_meta <- exposure_metadata
  if (is.null(exp_meta$N) || is.na(exp_meta$N)) {
    if ("samplesize.exposure" %in% names(exposure_snps)) {
      n_exp <- suppressWarnings(as.numeric(exposure_snps$samplesize.exposure[1]))
      if (!is.na(n_exp) && n_exp > 0) exp_meta$N <- n_exp
    }
  }
  if (is.null(exp_meta$N) || is.na(exp_meta$N)) {
    if (verbose) logger::log_warn("coloc: exposure N unknown; skipping")
    return(.coloc_empty_result_tbl())
  }

  results <- vector("list", nrow(loci))
  for (i in seq_len(nrow(loci))) {
    res_i <- tryCatch(
      .coloc_run_one_locus(
        locus = loci[i, , drop = FALSE],
        iv_pos = iv_pos,
        exposure_region_fetcher = exposure_region_fetcher,
        exp_meta = exp_meta,
        ancestry = ancestry,
        outcome_region_fetcher = outcome_region_fetcher,
        outcome_metadata = outcome_metadata,
        plot_dir = plot_dir,
        cache_dir = cache_dir,
        priors = priors,
        susie_always = susie_always,
        verbose = verbose
      ),
      error = function(e) {
        logger::log_warn("coloc: locus {i} failed: {conditionMessage(e)}")
        NULL
      }
    )
    if (!is.null(res_i)) results[[i]] <- res_i
  }
  results <- results[!vapply(results, is.null, logical(1))]
  if (!length(results)) return(.coloc_empty_result_tbl())
  dplyr::bind_rows(results)
}

.coloc_empty_result_tbl <- function() {
  tibble::tibble(
    locus_id        = character(),
    chr             = character(),
    start           = integer(),
    end             = integer(),
    n_ivs_in_locus  = integer(),
    n_snps_used     = integer(),
    pp_h4_abf       = double(),
    pp_h4_susie_max = double(),
    method_passed   = character(),
    coloc_dir       = character()
  )
}

.coloc_iv_positions <- function(hdat_use, exposure_snps) {
  empty <- tibble::tibble(SNP = character(), chr = character(), pos = integer())
  if (is.null(hdat_use) || NROW(hdat_use) == 0) return(empty)
  out <- empty
  if (all(c("chr.outcome", "pos.outcome") %in% names(hdat_use))) {
    out <- tibble::tibble(
      SNP = as.character(hdat_use$SNP),
      chr = as.character(hdat_use$chr.outcome),
      pos = suppressWarnings(as.integer(hdat_use$pos.outcome))
    )
  } else {
    pcols <- list(c("Chr", "Pos"),
                  c("chr.exposure", "pos.exposure"),
                  c("panukb_chrom", "panukb_pos"))
    for (pc in pcols) {
      if (all(pc %in% names(exposure_snps))) {
        idx <- match(hdat_use$SNP, exposure_snps$SNP)
        ok <- !is.na(idx)
        out <- tibble::tibble(
          SNP = as.character(hdat_use$SNP[ok]),
          chr = as.character(exposure_snps[[pc[1]]][idx[ok]]),
          pos = suppressWarnings(as.integer(exposure_snps[[pc[2]]][idx[ok]]))
        )
        break
      }
    }
  }
  out <- out[!is.na(out$chr) & nzchar(out$chr) & !is.na(out$pos), , drop = FALSE]
  out$chr <- sub("^chr", "", out$chr)
  out
}

.coloc_build_loci <- function(iv_pos, window_kb) {
  win <- as.integer(window_kb) * 1000L
  iv_pos <- iv_pos[order(iv_pos$chr, iv_pos$pos), , drop = FALSE]
  iv_pos$lo <- pmax(1L, iv_pos$pos - win)
  iv_pos$hi <- iv_pos$pos + win

  rows <- list()
  for (ch in unique(iv_pos$chr)) {
    sub <- iv_pos[iv_pos$chr == ch, , drop = FALSE]
    sub <- sub[order(sub$lo), , drop = FALSE]
    cur_lo <- sub$lo[1]; cur_hi <- sub$hi[1]; cur_ivs <- sub$SNP[1]
    for (k in seq_len(nrow(sub))[-1]) {
      if (sub$lo[k] <= cur_hi) {
        cur_hi <- max(cur_hi, sub$hi[k])
        cur_ivs <- c(cur_ivs, sub$SNP[k])
      } else {
        rows[[length(rows) + 1]] <- tibble::tibble(chr = ch, start = cur_lo, end = cur_hi,
                                                   ivs = list(cur_ivs))
        cur_lo <- sub$lo[k]; cur_hi <- sub$hi[k]; cur_ivs <- sub$SNP[k]
      }
    }
    rows[[length(rows) + 1]] <- tibble::tibble(chr = ch, start = cur_lo, end = cur_hi,
                                               ivs = list(cur_ivs))
  }
  if (!length(rows)) {
    return(tibble::tibble(chr = character(), start = integer(), end = integer(), ivs = list()))
  }
  dplyr::bind_rows(rows)
}

.coloc_harmonise <- function(exp_region, out_region) {
  exp_region$chr <- suppressWarnings(as.integer(sub("^chr", "", as.character(exp_region$chr))))
  out_region$chr <- suppressWarnings(as.integer(sub("^chr", "", as.character(out_region$chr))))
  exp_use <- tibble::tibble(
    chr = exp_region$chr, pos = exp_region$pos, snp_rs = as.character(exp_region$SNP),
    ea_exp = toupper(as.character(exp_region$effect_allele)),
    oa_exp = toupper(as.character(exp_region$other_allele)),
    beta_exp = exp_region$beta, se_exp = exp_region$se, eaf_exp = exp_region$eaf,
    pval_exp = exp_region$pval
  )
  out_use <- tibble::tibble(
    chr = out_region$chr, pos = out_region$pos,
    ea_out = toupper(as.character(out_region$effect_allele)),
    oa_out = toupper(as.character(out_region$other_allele)),
    beta_out = out_region$beta, se_out = out_region$se, eaf_out = out_region$eaf,
    pval_out = out_region$pval
  )
  j <- dplyr::inner_join(exp_use, out_use, by = c("chr", "pos"))
  if (nrow(j) == 0) return(j)

  same <- (j$ea_exp == j$ea_out) & (j$oa_exp == j$oa_out)
  flip <- (j$ea_exp == j$oa_out) & (j$oa_exp == j$ea_out)
  keep <- same | flip
  j <- j[keep, , drop = FALSE]
  flip_kept <- flip[keep]
  if (any(flip_kept)) {
    j$beta_out[flip_kept] <- -j$beta_out[flip_kept]
    if (any(!is.na(j$eaf_out[flip_kept]))) {
      j$eaf_out[flip_kept] <- 1 - j$eaf_out[flip_kept]
    }
  }

  pal <- (j$ea_exp == "A" & j$oa_exp == "T") | (j$ea_exp == "T" & j$oa_exp == "A") |
         (j$ea_exp == "C" & j$oa_exp == "G") | (j$ea_exp == "G" & j$oa_exp == "C")
  ambig <- pal & !is.na(j$eaf_exp) & j$eaf_exp >= 0.42 & j$eaf_exp <= 0.58
  j <- j[!ambig, , drop = FALSE]

  j <- j[is.finite(j$beta_exp) & is.finite(j$beta_out) &
         j$se_exp > 0 & j$se_out > 0, , drop = FALSE]

  j <- j[!duplicated(paste(j$chr, j$pos)), , drop = FALSE]
  j
}

.coloc_run_one_locus <- function(locus, iv_pos,
                                 exposure_region_fetcher, exp_meta, ancestry,
                                 outcome_region_fetcher, outcome_metadata,
                                 plot_dir, cache_dir, priors, susie_always, verbose) {
  locus_label <- sprintf("chr%s_%d_%d", locus$chr, locus$start, locus$end)
  locus_dir <- if (nzchar(plot_dir)) file.path(plot_dir, locus_label) else ""

  ivs_here <- iv_pos$SNP[
    iv_pos$chr == locus$chr & iv_pos$pos >= locus$start & iv_pos$pos <= locus$end
  ]
  n_ivs_here <- length(ivs_here)

  exp_region <- exposure_region_fetcher(locus$chr, locus$start, locus$end)
  if (is.null(exp_region) || nrow(exp_region) == 0) {
    if (verbose) logger::log_warn("coloc: empty exposure region {locus_label}; skipping")
    return(NULL)
  }
  out_region <- outcome_region_fetcher(locus$chr, locus$start, locus$end)
  if (is.null(out_region) || nrow(out_region) == 0) {
    if (verbose) logger::log_warn("coloc: empty outcome region {locus_label}; skipping")
    return(NULL)
  }

  harm <- .coloc_harmonise(exp_region, out_region)
  if (nrow(harm) < 5) {
    if (verbose) logger::log_warn("coloc: locus {locus_label} has <5 SNPs after harmonisation; skipping")
    return(NULL)
  }

  maf_exp <- pmin(harm$eaf_exp, 1 - harm$eaf_exp)
  maf_out <- pmin(harm$eaf_out, 1 - harm$eaf_out)
  maf_exp[is.na(maf_exp) | maf_exp <= 0] <- 0.05
  maf_out[is.na(maf_out) | maf_out <= 0] <- 0.05
  snp_ids <- paste0("chr", harm$chr, "_", harm$pos)

  D1 <- list(beta = harm$beta_exp, varbeta = harm$se_exp^2, MAF = maf_exp,
             N = exp_meta$N, type = exp_meta$type, snp = snp_ids, position = harm$pos)
  if (identical(exp_meta$type, "cc") && !is.na(exp_meta$s)) D1$s <- exp_meta$s
  D2 <- list(beta = harm$beta_out, varbeta = harm$se_out^2, MAF = maf_out,
             N = outcome_metadata$N, type = outcome_metadata$type,
             snp = snp_ids, position = harm$pos)
  if (identical(outcome_metadata$type, "cc") && !is.na(outcome_metadata$s)) {
    D2$s <- outcome_metadata$s
  }

  set.seed(20260501L)
  abf_res <- tryCatch(
    coloc::coloc.abf(D1, D2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12),
    error = function(e) { logger::log_warn("coloc.abf failed at {locus_label}: {conditionMessage(e)}"); NULL }
  )
  pp_h4_abf <- if (!is.null(abf_res) && !is.null(abf_res$summary))
    suppressWarnings(as.numeric(abf_res$summary["PP.H4.abf"])) else NA_real_

  susie_res <- NULL
  pp_h4_susie_max <- NA_real_
  try_susie <- requireNamespace("susieR", quietly = TRUE) &&
               (susie_always || n_ivs_here >= 2)
  if (try_susie) {
    ld_obj <- tryCatch(
      .coloc_locus_ld(harm, ancestry, cache_dir, verbose),
      error = function(e) { logger::log_warn("coloc LD failed at {locus_label}: {conditionMessage(e)}"); NULL }
    )
    if (!is.null(ld_obj) && !is.null(ld_obj$matrix) && nrow(ld_obj$matrix) >= 5) {
      ld_canon <- ld_obj$canon  # named character vector: rsid -> chr_pos_id
      common <- intersect(snp_ids, ld_canon)
      if (length(common) >= 5) {
        rs_keep <- names(ld_canon)[ld_canon %in% common]
        ld_sub <- ld_obj$matrix[rs_keep, rs_keep, drop = FALSE]
        canon_keep <- ld_canon[rs_keep]
        order_idx <- match(canon_keep, snp_ids)
        keep_mask <- !is.na(order_idx)
        rs_keep <- rs_keep[keep_mask]
        canon_keep <- canon_keep[keep_mask]
        ld_sub <- ld_sub[rs_keep, rs_keep, drop = FALSE]
        dimnames(ld_sub) <- list(canon_keep, canon_keep)
        sel <- match(canon_keep, snp_ids)
        sub_d <- function(D) {
          dd <- list(beta = D$beta[sel], varbeta = D$varbeta[sel], MAF = D$MAF[sel],
                     N = D$N, type = D$type, snp = D$snp[sel],
                     position = D$position[sel], LD = ld_sub)
          if (!is.null(D$s)) dd$s <- D$s
          dd
        }
        D1s <- sub_d(D1); D2s <- sub_d(D2)
        set.seed(20260501L)
        susie_res <- tryCatch(coloc::coloc.susie(D1s, D2s),
                              error = function(e) { logger::log_warn("coloc.susie failed at {locus_label}: {conditionMessage(e)}"); NULL })
        if (!is.null(susie_res) && !is.null(susie_res$summary) && NROW(susie_res$summary) > 0) {
          v <- suppressWarnings(as.numeric(susie_res$summary$PP.H4.abf))
          v <- v[is.finite(v)]
          if (length(v)) pp_h4_susie_max <- max(v)
        }
      } else if (verbose) {
        logger::log_warn("coloc.susie at {locus_label}: <5 SNPs overlap with LD panel; skipped")
      }
    }
  }

  if (nzchar(locus_dir)) {
    .coloc_save_locus_artefacts(locus_dir, locus_label, harm, ivs_here,
                                D1, D2, abf_res, susie_res)
  }

  method_passed <- if (!is.na(pp_h4_susie_max) && pp_h4_susie_max > 0.80) "susie"
                   else if (!is.na(pp_h4_abf) && pp_h4_abf > 0.80) "abf"
                   else "none"

  tibble::tibble(
    locus_id        = locus_label,
    chr             = as.character(locus$chr),
    start           = as.integer(locus$start),
    end             = as.integer(locus$end),
    n_ivs_in_locus  = as.integer(n_ivs_here),
    n_snps_used     = as.integer(nrow(harm)),
    pp_h4_abf       = pp_h4_abf,
    pp_h4_susie_max = pp_h4_susie_max,
    method_passed   = method_passed,
    coloc_dir       = if (nzchar(locus_dir)) locus_dir else NA_character_
  )
}

.coloc_locus_ld <- function(harm, ancestry, cache_dir, verbose) {
  rs <- harm$snp_rs
  has_rs <- grepl("^rs[0-9]+$", rs)
  if (sum(has_rs) < 5) return(NULL)
  exp_for_ld <- tibble::tibble(
    SNP = rs[has_rs],
    id.exposure = "coloc_region",
    beta.exposure = harm$beta_exp[has_rs],
    se.exposure   = harm$se_exp[has_rs],
    pval.exposure = harm$pval_exp[has_rs],
    effect_allele.exposure = harm$ea_exp[has_rs],
    other_allele.exposure  = harm$oa_exp[has_rs]
  )
  ld <- tryCatch(
    compute_ld_matrix(exposure_snps = exp_for_ld, ancestry = ancestry,
                      cache_dir = cache_dir, confirm = "no",
                      force_refresh = FALSE, verbose = verbose,
                      include_dropped = FALSE),
    error = function(e) NULL
  )
  if (is.null(ld) || is.null(ld$matrix) || nrow(ld$matrix) == 0) return(NULL)
  used <- rownames(ld$matrix)
  canon_lookup <- paste0("chr", harm$chr[has_rs], "_", harm$pos[has_rs])
  names(canon_lookup) <- exp_for_ld$SNP
  canon <- canon_lookup[used]
  list(matrix = ld$matrix, canon = canon)
}

.coloc_save_locus_artefacts <- function(locus_dir, locus_label, harm, ivs_here,
                                        D1, D2, abf_res, susie_res) {
  if (!dir.exists(locus_dir)) dir.create(locus_dir, recursive = TRUE, showWarnings = FALSE)

  exp_tbl <- tibble::tibble(SNP = harm$snp_rs, pos = harm$pos, pval = harm$pval_exp)
  out_tbl <- tibble::tibble(SNP = harm$snp_rs, pos = harm$pos, pval = harm$pval_out)

  stack_p <- tryCatch(
    coloc_stacked_plot(exp_tbl, out_tbl, ivs = ivs_here, title = locus_label),
    error = function(e) NULL
  )
  if (!is.null(stack_p)) {
    .coloc_save_pub_plot(stack_p, file.path(locus_dir, "stacked_manhattan"))
  }

  .coloc_save_base_plot(function() coloc::plot_dataset(D1),
                        file.path(locus_dir, "coloc_dataset_exposure"))
  .coloc_save_base_plot(function() coloc::plot_dataset(D2),
                        file.path(locus_dir, "coloc_dataset_outcome"))

  if (!is.null(abf_res)) {
    s <- as.numeric(abf_res$summary)
    abf_summary <- tibble::tibble(
      locus_id  = locus_label,
      nsnps     = abf_res$summary["nsnps"],
      PP.H0     = unname(abf_res$summary["PP.H0.abf"]),
      PP.H1     = unname(abf_res$summary["PP.H1.abf"]),
      PP.H2     = unname(abf_res$summary["PP.H2.abf"]),
      PP.H3     = unname(abf_res$summary["PP.H3.abf"]),
      PP.H4     = unname(abf_res$summary["PP.H4.abf"])
    )
    if (!is.null(abf_res$results) && "snp" %in% names(abf_res$results)) {
      top <- abf_res$results[order(-abf_res$results$SNP.PP.H4)[1], , drop = FALSE]
      abf_summary$top_snp <- as.character(top$snp)
    }
    utils::write.csv(abf_summary, file.path(locus_dir, "coloc_abf_summary.csv"), row.names = FALSE)
    .coloc_save_base_plot(function() coloc::sensitivity(abf_res, rule = "H4 > 0.8"),
                          file.path(locus_dir, "coloc_sensitivity"))
  }

  if (!is.null(susie_res) && !is.null(susie_res$summary) && NROW(susie_res$summary) > 0) {
    susie_tbl <- tibble::as_tibble(susie_res$summary)
    susie_tbl$locus_id <- locus_label
    utils::write.csv(susie_tbl, file.path(locus_dir, "coloc_susie_summary.csv"),
                     row.names = FALSE)
  }

  utils::write.csv(harm, file.path(locus_dir, "coloc_aligned_data.csv"), row.names = FALSE)
  invisible(NULL)
}

.coloc_panukb_outcome_fetcher <- function(rec, ancestry, cache_dir, verbose) {
  function(chr, start, end) {
    coloc_fetch_outcome_panukb_region(
      rec = rec, ancestry = ancestry, chr = chr, start = start, end = end,
      cache_dir = cache_dir, verbose = verbose
    )
  }
}
