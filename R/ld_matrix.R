# R/ld_matrix.R
# Pairwise LD (signed r) matrix between exposure SNPs, computed locally with
# PLINK against MRC-IEU 1000G reference panels. Mirrors the
# download/cache/confirm pattern of neale_gwas_checker() in R/provider_neale.R.

#' Ensure a 1000G LD reference panel exists locally (download if needed)
#'
#' Checks for `<ld_ref_dir>/<POP>.bed/.bim/.fam` for the requested
#' superpopulation. If they are missing, downloads the MRC-IEU 1000G
#' reference *bundle* (`http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz`,
#' ~1.5 GB compressed; ~5 GB extracted) which contains all five
#' superpopulations (EUR, AFR, EAS, SAS, AMR) and extracts them in place.
#' Subsequent calls for any other population in the bundle reuse the
#' already-extracted files without re-downloading.
#'
#' Mirrors the download pattern used by [neale_gwas_checker()]: bumps the
#' download timeout, prompts the user before fetching anything large,
#' retries with exponential backoff, and removes partial files on failure.
#'
#' @param pop One of `"EUR"`, `"AFR"`, `"EAS"`, `"SAS"`.
#' @param ld_ref_dir Directory to hold the panel files. Created if missing.
#' @param verbose Logical.
#' @param confirm One of `"ask"`, `"yes"`, `"no"`. `"ask"` (default) prompts
#'   interactively; `"yes"` proceeds silently; `"no"` aborts if the panel is
#'   missing.
#' @param estimate_mb Approximate compressed download size for the user
#'   prompt (informational; defaults to 1500, the full 1kg.v3 bundle size).
#' @return Invisibly, a list with `bfile_prefix`, `present`, `downloaded`,
#'   `aborted`.
#' @export
ld_reference_checker <- function(pop, ld_ref_dir,
                                 verbose = TRUE,
                                 confirm = c("ask", "yes", "no"),
                                 estimate_mb = 1500) {
  confirm <- match.arg(confirm)
  pop <- toupper(trimws(as.character(pop)))
  if (!pop %in% c("EUR", "AFR", "EAS", "SAS")) {
    stop("ld_reference_checker(): pop must be one of EUR, AFR, EAS, SAS.", call. = FALSE)
  }
  if (missing(ld_ref_dir) || !nzchar(ld_ref_dir)) {
    stop("ld_reference_checker(): ld_ref_dir is required.", call. = FALSE)
  }
  if (!dir.exists(ld_ref_dir)) dir.create(ld_ref_dir, recursive = TRUE, showWarnings = FALSE)

  old_to <- getOption("timeout")
  options(timeout = max(3600, old_to))
  on.exit(options(timeout = old_to), add = TRUE)

  bfile_prefix <- file.path(ld_ref_dir, pop)
  parts <- paste0(bfile_prefix, c(".bed", ".bim", ".fam"))
  min_bed <- 50 * 1024 * 1024  # 50 MB sanity for .bed
  if (all(file.exists(parts))) {
    bed_sz <- suppressWarnings(file.info(parts[1])$size)
    if (!is.na(bed_sz) && bed_sz >= min_bed) {
      if (verbose) logger::log_info("LD reference present for {pop} at '{normalizePath(bfile_prefix, mustWork = FALSE)}'")
      return(invisible(list(bfile_prefix = bfile_prefix, present = TRUE,
                            downloaded = FALSE, aborted = FALSE)))
    }
    if (verbose) logger::log_warn("LD reference for {pop} looks truncated (.bed = {format(bed_sz, big.mark = ',')} bytes); will redownload bundle.")
    try(unlink(parts, force = TRUE), silent = TRUE)
  }

  url  <- "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
  dest <- file.path(ld_ref_dir, "1kg.v3.tgz")

  logger::log_warn(
    "LD reference panel for {pop} is missing; downloading ~{estimate_mb} MB MRC-IEU 1kg.v3 bundle (covers all 5 superpopulations)."
  )

  proceed <- FALSE
  if (confirm == "yes") {
    proceed <- TRUE
  } else if (confirm == "no") {
    proceed <- FALSE
  } else {
    if (interactive()) {
      ans <- readline(prompt = sprintf(
        "About to download ~%d MB (%s, contains EUR/AFR/EAS/SAS/AMR) into '%s'. Proceed? [y/N]: ",
        as.integer(estimate_mb), basename(url), normalizePath(ld_ref_dir, mustWork = FALSE)
      ))
      proceed <- grepl("^[Yy]", ans)
    } else {
      stop("LD reference download required but session is non-interactive and confirm='ask'. ",
           "Re-run with confirm='yes' to proceed, or 'no' to skip.", call. = FALSE)
    }
  }
  if (!proceed) {
    logger::log_warn("User declined LD reference download for {pop}.")
    return(invisible(list(bfile_prefix = bfile_prefix, present = FALSE,
                          downloaded = FALSE, aborted = TRUE)))
  }

  retries <- 3L
  min_tgz <- 500 * 1024 * 1024  # 500 MB compressed sanity (bundle is ~1.5 GB)
  ok <- FALSE
  for (r in seq_len(retries)) {
    ok <- tryCatch({
      utils::download.file(url, destfile = dest, mode = "wb", quiet = !verbose)
      TRUE
    }, error = function(e) {
      logger::log_warn("LD bundle download attempt {r}/{retries} failed: {conditionMessage(e)}")
      FALSE
    })
    if (ok && file.exists(dest)) {
      sz <- suppressWarnings(file.info(dest)$size)
      if (!is.na(sz) && sz >= min_tgz) break
      logger::log_warn("Downloaded {basename(dest)} but size looks small ({format(sz, big.mark = ',')} bytes); retrying ({r}/{retries})")
      ok <- FALSE
    }
    if (!ok && r < retries) Sys.sleep(2^r)
  }
  if (!ok || !file.exists(dest)) {
    if (file.exists(dest)) try(unlink(dest, force = TRUE), silent = TRUE)
    stop(sprintf("Failed to download LD reference bundle from %s", url), call. = FALSE)
  }

  ext_ok <- tryCatch({
    utils::untar(dest, exdir = ld_ref_dir)
    TRUE
  }, error = function(e) {
    logger::log_warn("untar({dest}) failed: {conditionMessage(e)}")
    FALSE
  })
  try(unlink(dest, force = TRUE), silent = TRUE)  # remove .tgz after extract

  # The 1kg.v3 bundle may extract files either flat (<ld_ref_dir>/EUR.bed)
  # or nested under a top-level directory (<ld_ref_dir>/1kg.v3/EUR.bed).
  # Normalise to flat layout so bfile_prefix is stable.
  if (!all(file.exists(parts))) {
    found <- list.files(
      ld_ref_dir,
      pattern = "^(EUR|AFR|EAS|SAS|AMR)\\.(bed|bim|fam)$",
      recursive = TRUE, full.names = TRUE
    )
    for (f in found) {
      target <- file.path(ld_ref_dir, basename(f))
      if (!identical(normalizePath(f, mustWork = FALSE),
                     normalizePath(target, mustWork = FALSE))) {
        try(file.rename(f, target), silent = TRUE)
      }
    }
    # Clean up any now-empty extraction subdirectory.
    sub_dirs <- setdiff(list.dirs(ld_ref_dir, recursive = FALSE), ld_ref_dir)
    for (d in sub_dirs) {
      if (length(list.files(d, recursive = TRUE)) == 0L) {
        try(unlink(d, recursive = TRUE, force = TRUE), silent = TRUE)
      }
    }
  }
  if (!all(file.exists(parts))) {
    stop(sprintf("LD reference extraction for %s did not produce expected files at %s.",
                 pop, normalizePath(ld_ref_dir, mustWork = FALSE)), call. = FALSE)
  }

  bed_sz <- suppressWarnings(file.info(parts[1])$size)
  if (is.na(bed_sz) || bed_sz < min_bed) {
    stop(sprintf("LD reference .bed for %s looks truncated (%s bytes).",
                 pop, format(bed_sz, big.mark = ",")), call. = FALSE)
  }
  logger::log_info("LD reference {pop} ready at '{normalizePath(bfile_prefix, mustWork = FALSE)}' (.bed = {format(bed_sz, big.mark = ',')} bytes)")
  invisible(list(bfile_prefix = bfile_prefix, present = TRUE,
                 downloaded = TRUE, aborted = FALSE))
}


#' Compute a pairwise LD (signed r) matrix between exposure SNPs
#'
#' Wraps PLINK 1.9/2.0 to compute the signed r LD matrix between the SNPs in
#' `exposure_snps`, using the MRC-IEU 1000G reference panel for the
#' superpopulation derived from `ancestry`. Results are cached on disk by
#' (sorted SNPs + population + panel version), and optionally written as a
#' CSV at `out_csv`. Sign is harmonised so that a positive entry corresponds
#' to the exposure effect alleles being correlated in the reference panel.
#'
#' Note: `force_refresh = TRUE` invalidates the LD-matrix RDS cache but does
#' **not** redownload the (large, stable) reference panel.
#'
#' @param exposure_snps Data frame of exposure instruments. Must contain
#'   `SNP`, `effect_allele.exposure`, `other_allele.exposure` (validated via
#'   [assert_exposure()]).
#' @param ancestry Character ancestry label (mapped to a 1000G superpop via
#'   [ld_pop_from_ancestry()]).
#' @param cache_dir Cache root (default [ardmr_cache_dir()]). RDS caches go to
#'   `<cache_dir>/ld/`; the reference panel goes to `<ld_ref_dir>` (default
#'   `<cache_dir>/ld_reference`).
#' @param out_csv Optional path for the CSV output. If `NULL`, no CSV is
#'   written and only the cache + return value are produced.
#' @param confirm One of `"ask"`, `"yes"`, `"no"` (passed to
#'   [ld_reference_checker()] when the panel is missing).
#' @param force_refresh Logical; if `TRUE`, ignore any RDS cache hit and
#'   recompute (the reference panel is left in place).
#' @param verbose Logical.
#' @param plink_bin Optional explicit path to the PLINK binary. If `NULL`,
#'   resolution order is `ARDMR_PLINK_BIN` env var,
#'   `genetics.binaRies::get_plink_binary()`, then `Sys.which("plink"|"plink2")`.
#' @param ld_ref_dir Optional directory holding the reference panel. Defaults
#'   to `file.path(cache_dir, "ld_reference")`.
#' @param include_dropped Logical; if `TRUE` (default), SNPs not represented
#'   in the panel (or with allele mismatches) appear as NA rows/columns in
#'   the CSV so the dimensionality matches the input. If `FALSE`, they are
#'   omitted.
#' @return Invisibly, a list with `matrix` (numeric matrix with SNP
#'   row/colnames), `snps_used`, `snps_dropped` (data frame of `SNP, reason`),
#'   `ld_pop`, `csv_path`, `cache_path`.
#' @export
compute_ld_matrix <- function(exposure_snps,
                              ancestry,
                              cache_dir       = ardmr_cache_dir(),
                              out_csv         = NULL,
                              confirm         = c("ask", "yes", "no"),
                              force_refresh   = FALSE,
                              verbose         = TRUE,
                              plink_bin       = NULL,
                              ld_ref_dir      = NULL,
                              include_dropped = TRUE) {
  confirm <- match.arg(confirm)
  if (missing(exposure_snps) || !is.data.frame(exposure_snps)) {
    stop("compute_ld_matrix(): `exposure_snps` must be a data frame.", call. = FALSE)
  }
  if (missing(ancestry) || !nzchar(as.character(ancestry))) {
    stop("compute_ld_matrix(): `ancestry` is required.", call. = FALSE)
  }
  assert_exposure(exposure_snps)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  exp <- tibble::as_tibble(exposure_snps)
  exp$SNP <- .normalize_rsid(exp$SNP)
  exp <- exp[!is.na(exp$SNP) & nzchar(exp$SNP), , drop = FALSE]
  exp$effect_allele.exposure <- toupper(trimws(as.character(exp$effect_allele.exposure)))
  exp$other_allele.exposure  <- toupper(trimws(as.character(exp$other_allele.exposure)))
  exp <- exp[!duplicated(exp$SNP), , drop = FALSE]

  if (nrow(exp) == 0) {
    stop("compute_ld_matrix(): no usable rsids in exposure_snps after normalization.", call. = FALSE)
  }

  ld_pop <- ld_pop_from_ancestry(ancestry)
  panel_version <- "1kg.v3"

  cache_subdir <- file.path(cache_dir, "ld")
  if (!dir.exists(cache_subdir)) dir.create(cache_subdir, recursive = TRUE, showWarnings = FALSE)
  cache_key <- .ld_cache_key(exp$SNP, ld_pop, panel_version)
  cache_path <- file.path(cache_subdir, paste0(cache_key, ".rds"))

  used_cache <- FALSE
  result <- NULL
  if (!isTRUE(force_refresh) && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) e)
    if (inherits(cached, "list") && !is.null(cached$matrix)) {
      if (verbose) logger::log_info("LD: using cached LD matrix at '{basename(cache_path)}'")
      result <- cached
      used_cache <- TRUE
    } else {
      if (verbose) logger::log_warn("LD: cache at '{basename(cache_path)}' was unusable; recomputing.")
      try(unlink(cache_path, force = TRUE), silent = TRUE)
    }
  }

  if (is.null(result)) {
    if (nrow(exp) == 1L) {
      if (verbose) logger::log_info("LD: single-SNP exposure '{exp$SNP[1]}' — writing 1x1 identity without invoking PLINK.")
      mat <- matrix(1.0, 1L, 1L, dimnames = list(exp$SNP, exp$SNP))
      result <- list(
        matrix = mat,
        snps_input    = exp$SNP,
        snps_used     = exp$SNP,
        snps_dropped  = data.frame(SNP = character(0), reason = character(0),
                                   stringsAsFactors = FALSE),
        ld_pop        = ld_pop,
        panel_version = panel_version,
        computed_at   = Sys.time()
      )
    } else {
      if (is.null(ld_ref_dir) || !nzchar(ld_ref_dir)) {
        ld_ref_dir <- file.path(cache_dir, "ld_reference")
      }
      ref <- ld_reference_checker(
        pop = ld_pop, ld_ref_dir = ld_ref_dir,
        verbose = verbose, confirm = confirm
      )
      if (isTRUE(ref$aborted)) {
        stop("compute_ld_matrix(): user declined LD reference download; cannot compute matrix.",
             call. = FALSE)
      }
      bfile_prefix <- ref$bfile_prefix

      plink_bin_resolved <- .find_plink_bin(plink_bin)
      if (verbose) logger::log_info("LD: PLINK at '{plink_bin_resolved}'")

      bim_path <- paste0(bfile_prefix, ".bim")
      bim <- data.table::fread(
        bim_path, header = FALSE, sep = "\t",
        col.names = c("chr", "snp", "cm", "bp", "A1", "A2"),
        showProgress = FALSE
      )
      bim$snp <- as.character(bim$snp)
      bim$A1 <- toupper(as.character(bim$A1))
      bim$A2 <- toupper(as.character(bim$A2))

      in_panel <- exp$SNP %in% bim$snp
      not_in_panel <- exp$SNP[!in_panel]
      exp_panel <- exp[in_panel, , drop = FALSE]

      if (nrow(exp_panel) == 0L) {
        logger::log_warn("LD: 0 of {length(exp$SNP)} exposure SNPs found in {ld_pop} panel; returning empty matrix.")
        mat <- matrix(numeric(0), 0L, 0L)
        result <- list(
          matrix = mat,
          snps_input    = exp$SNP,
          snps_used     = character(0),
          snps_dropped  = data.frame(
            SNP = exp$SNP, reason = "not_in_panel",
            stringsAsFactors = FALSE
          ),
          ld_pop        = ld_pop,
          panel_version = panel_version,
          computed_at   = Sys.time()
        )
      } else if (nrow(exp_panel) == 1L) {
        if (verbose) logger::log_info("LD: only 1 SNP overlaps with panel; writing 1x1 identity.")
        mat <- matrix(1.0, 1L, 1L, dimnames = list(exp_panel$SNP, exp_panel$SNP))
        dropped <- data.frame(
          SNP = not_in_panel, reason = rep("not_in_panel", length(not_in_panel)),
          stringsAsFactors = FALSE
        )
        result <- list(
          matrix = mat,
          snps_input    = exp$SNP,
          snps_used     = exp_panel$SNP,
          snps_dropped  = dropped,
          ld_pop        = ld_pop,
          panel_version = panel_version,
          computed_at   = Sys.time()
        )
      } else {
        snps_file  <- tempfile("ld_snps_", tmpdir = cache_subdir, fileext = ".txt")
        out_prefix <- tempfile("ld_plink_", tmpdir = cache_subdir)
        on.exit({
          for (suf in c("", ".ld", ".log", ".snplist", ".nosex")) {
            f <- paste0(out_prefix, suf)
            if (file.exists(f)) try(unlink(f, force = TRUE), silent = TRUE)
          }
          if (file.exists(snps_file)) try(unlink(snps_file, force = TRUE), silent = TRUE)
        }, add = TRUE)

        writeLines(exp_panel$SNP, snps_file)

        plink_out <- .run_plink_r_square(
          plink_bin    = plink_bin_resolved,
          bfile_prefix = bfile_prefix,
          snps_file    = snps_file,
          out_prefix   = out_prefix,
          verbose      = verbose
        )
        raw_mat       <- plink_out$matrix
        snps_in_panel <- plink_out$snps   # order matches matrix rows/cols

        bim_lu <- bim[match(snps_in_panel, bim$snp), c("snp", "A1", "A2")]
        ord <- match(snps_in_panel, exp$SNP)
        exp_lu <- exp[ord, c("SNP", "effect_allele.exposure", "other_allele.exposure"), drop = FALSE]

        harm <- .harmonise_ld_signs(
          ld_mat           = raw_mat,
          bim_alleles      = list(A1 = bim_lu$A1, A2 = bim_lu$A2),
          exposure_alleles = list(
            effect = exp_lu$effect_allele.exposure,
            other  = exp_lu$other_allele.exposure
          ),
          snps = snps_in_panel
        )

        kept_snps <- snps_in_panel[harm$kept_idx]
        mat <- harm$matrix
        if (length(kept_snps) > 0L) {
          dimnames(mat) <- list(kept_snps, kept_snps)
        }

        dropped <- rbind(
          data.frame(SNP = not_in_panel,
                     reason = rep("not_in_panel", length(not_in_panel)),
                     stringsAsFactors = FALSE),
          data.frame(SNP = harm$dropped_snps, reason = harm$dropped_reasons,
                     stringsAsFactors = FALSE)
        )

        if (nrow(dropped) > 0L) {
          logger::log_warn("LD: dropped {nrow(dropped)} SNP(s); see sidecar for reasons.")
        }
        if ("mr_keep" %in% names(exp)) {
          mk <- exp$mr_keep
          n_mrk_false <- if (is.logical(mk)) sum(!is.na(mk) & !mk) else 0L
          if (n_mrk_false > 0L) {
            logger::log_info("LD: {n_mrk_false} input SNP(s) had mr_keep=FALSE; included anyway in LD matrix.")
          }
        }

        result <- list(
          matrix = mat,
          snps_input    = exp$SNP,
          snps_used     = kept_snps,
          snps_dropped  = dropped,
          ld_pop        = ld_pop,
          panel_version = panel_version,
          computed_at   = Sys.time()
        )
      }
    }

    try(saveRDS(result, cache_path), silent = TRUE)
  }

  csv_path <- NULL
  if (!is.null(out_csv) && nzchar(out_csv)) {
    out_dir <- dirname(out_csv)
    if (nzchar(out_dir) && !dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    snps_for_csv <- if (isTRUE(include_dropped)) {
      result$snps_input
    } else {
      result$snps_used
    }
    .write_ld_csv(result$matrix, snps_for_csv, out_csv)
    csv_path <- out_csv
    if (nrow(result$snps_dropped) > 0L) {
      sidecar <- file.path(
        dirname(out_csv),
        sub("\\.csv$", "_dropped.csv", basename(out_csv))
      )
      utils::write.csv(result$snps_dropped, sidecar, row.names = FALSE)
    }
    if (verbose) logger::log_info("LD: CSV written to '{out_csv}' ({length(snps_for_csv)} SNPs)")
  }

  invisible(list(
    matrix       = result$matrix,
    snps_used    = result$snps_used,
    snps_dropped = result$snps_dropped,
    ld_pop       = result$ld_pop,
    csv_path     = csv_path,
    cache_path   = cache_path,
    cached       = used_cache
  ))
}


# ---- Internal helpers ------------------------------------------------------

.find_plink_bin <- function(plink_bin = NULL) {
  if (!is.null(plink_bin) && nzchar(plink_bin)) {
    if (!file.exists(plink_bin)) {
      stop(sprintf("PLINK binary not found at '%s'.", plink_bin), call. = FALSE)
    }
    return(plink_bin)
  }
  env_bin <- Sys.getenv("ARDMR_PLINK_BIN", unset = "")
  if (nzchar(env_bin) && file.exists(env_bin)) return(env_bin)

  if (requireNamespace("genetics.binaRies", quietly = TRUE)) {
    pb <- tryCatch(genetics.binaRies::get_plink_binary(), error = function(e) "")
    if (nzchar(pb) && file.exists(pb)) return(pb)
  }
  for (cand in c("plink", "plink2")) {
    pb <- Sys.which(cand)[[1]]
    if (nzchar(pb) && file.exists(pb)) return(unname(pb))
  }
  stop(
    paste(
      "Could not locate a PLINK binary. Options:",
      "  (1) install.packages('genetics.binaRies', repos = 'https://mrcieu.r-universe.dev')",
      "  (2) install PLINK 1.9 or 2.0 system-wide and ensure it is on PATH",
      "  (3) set ARDMR_PLINK_BIN to the binary path",
      "  (4) pass plink_bin = '/path/to/plink' to compute_ld_matrix()",
      sep = "\n"
    ),
    call. = FALSE
  )
}

.ld_cache_key <- function(snps, pop, panel_version = "1kg.v3") {
  payload <- list(
    snps = sort(unique(as.character(snps))),
    pop  = as.character(pop),
    ver  = as.character(panel_version)
  )
  hash <- digest::digest(payload, algo = "sha256")
  substr(hash, 1, 16)
}

.run_plink_r_square <- function(plink_bin, bfile_prefix, snps_file, out_prefix, verbose = TRUE) {
  args <- c(
    "--bfile",         bfile_prefix,
    "--extract",       snps_file,
    "--r",             "square",
    "--write-snplist",
    "--out",           out_prefix
  )
  status <- tryCatch(
    system2(plink_bin, args = args,
            stdout = if (verbose) "" else FALSE,
            stderr = if (verbose) "" else FALSE),
    error = function(e) {
      stop(sprintf("PLINK invocation failed: %s", conditionMessage(e)), call. = FALSE)
    }
  )

  ld_file       <- paste0(out_prefix, ".ld")
  snplist_file  <- paste0(out_prefix, ".snplist")
  log_file      <- paste0(out_prefix, ".log")

  if (!file.exists(ld_file) || !file.exists(snplist_file)) {
    log_tail <- if (file.exists(log_file)) {
      paste(utils::tail(readLines(log_file, warn = FALSE), 30), collapse = "\n")
    } else "(no PLINK log file found)"
    stop(sprintf(
      "PLINK did not produce expected outputs (status=%s).\n--- log tail ---\n%s",
      status, log_tail
    ), call. = FALSE)
  }

  mat <- as.matrix(data.table::fread(ld_file, header = FALSE, showProgress = FALSE))
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- NA_real_

  snps <- readLines(snplist_file, warn = FALSE)
  snps <- trimws(snps)
  snps <- snps[nzchar(snps)]

  if (nrow(mat) != length(snps) || ncol(mat) != length(snps)) {
    stop(sprintf("PLINK output dimensions mismatch: matrix %dx%d vs %d SNPs.",
                 nrow(mat), ncol(mat), length(snps)), call. = FALSE)
  }

  list(matrix = mat, snps = snps)
}

.harmonise_ld_signs <- function(ld_mat, bim_alleles, exposure_alleles, snps) {
  A1 <- toupper(bim_alleles$A1)
  A2 <- toupper(bim_alleles$A2)
  ea <- toupper(exposure_alleles$effect)
  oa <- toupper(exposure_alleles$other)

  signs <- rep(NA_real_, length(snps))
  reasons <- rep(NA_character_, length(snps))

  match_a1 <- ea == A1 & oa == A2
  match_a2 <- ea == A2 & oa == A1

  signs[match_a1] <- 1
  signs[match_a2] <- -1
  reasons[!match_a1 & !match_a2] <- "allele_mismatch"

  kept_idx <- which(!is.na(signs))
  dropped_idx <- which(is.na(signs))

  if (length(kept_idx) == 0L) {
    return(list(
      matrix         = matrix(numeric(0), 0L, 0L),
      kept_idx       = integer(0),
      dropped_idx    = dropped_idx,
      dropped_snps   = snps[dropped_idx],
      dropped_reasons = reasons[dropped_idx]
    ))
  }

  D <- diag(signs[kept_idx], nrow = length(kept_idx))
  sub_mat <- ld_mat[kept_idx, kept_idx, drop = FALSE]
  out <- D %*% sub_mat %*% D

  list(
    matrix          = out,
    kept_idx        = kept_idx,
    dropped_idx     = dropped_idx,
    dropped_snps    = snps[dropped_idx],
    dropped_reasons = reasons[dropped_idx]
  )
}

.write_ld_csv <- function(ld_mat, snps, path) {
  snps <- as.character(snps)
  n <- length(snps)
  out <- matrix(NA_real_, n, n, dimnames = list(snps, snps))
  if (n > 0L && !is.null(rownames(ld_mat)) && !is.null(colnames(ld_mat))) {
    common <- intersect(rownames(ld_mat), snps)
    if (length(common) > 0L) {
      out[common, common] <- ld_mat[common, common]
    }
  } else if (n > 0L && nrow(ld_mat) == n && ncol(ld_mat) == n) {
    out[] <- ld_mat
  }
  if (n > 0L) {
    diag_vals <- diag(out)
    diag_vals[is.na(diag_vals)] <- 1.0
    diag(out) <- diag_vals
  }

  df <- data.frame(SNP = snps, stringsAsFactors = FALSE)
  if (n > 0L) {
    body <- as.data.frame(out, stringsAsFactors = FALSE)
    colnames(body) <- snps
    df <- cbind(df, body)
  }
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}

