#' Ensure Neale sumstats exist locally (download if needed)
#' @param MR_df outcomes plan (must include columns "File" and "AWS File")
#' @param neale_dir directory where .tsv.bgz files should live
#' @param verbose logical
#' @param confirm one of "ask", "yes", "no".
#'   - "ask" (default): prompt the user if downloads are needed (interactive only).
#'   - "yes": proceed without prompting.
#'   - "no": do not download; return immediately.
#' @param estimate_gb numeric ~ how much space the full set could take (informational only).
#' @return (invisible) list with counts and `aborted` flag
#' @export
neale_gwas_checker <- function(MR_df, neale_dir, verbose = TRUE,
                               confirm = c("ask","yes","no"),
                               estimate_gb = 500) {
  confirm <- match.arg(confirm)

  # --- guards ---
  if (!all(c("File", "AWS File") %in% names(MR_df))) {
    stop("neale_gwas_checker(): MR_df must contain columns 'File' and 'AWS File'.", call. = FALSE)
  }
  if (!dir.exists(neale_dir)) dir.create(neale_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- change timeout to longer than 60 seconds
  old_to <- getOption("timeout")
  options(timeout = max(3600, old_to))   # e.g. 60 mins
  on.exit(options(timeout = old_to), add = TRUE)

  .clean <- function(x) {
    x <- as.character(x); x <- trimws(x)
    sub("^(['\"])(.*)\\1$", "\\2", x, perl = TRUE)
  }

  files <- .clean(MR_df$File)
  urls  <- .clean(MR_df[["AWS File"]])

  keep <- nzchar(files) & !is.na(files)
  files <- files[keep]
  urls  <- urls[keep]

  split_urls <- split(urls, files)
  uniq_files <- names(split_urls)
  uniq_urls  <- vapply(split_urls, function(v) {
    w <- v[nzchar(v) & !is.na(v)]
    if (length(w)) w[[1]] else NA_character_
  }, character(1))

  total <- length(uniq_files)
  if (total == 0) {
    logger::log_info("Neale GWAS check: no files listed in MR_df$File; nothing to do.")
    return(invisible(list(n_total = 0L, n_present = 0L, n_downloaded = 0L, n_failed = 0L, aborted = FALSE)))
  }

  dest_paths <- file.path(neale_dir, basename(uniq_files))
  exists     <- file.exists(dest_paths)
  n_present  <- sum(exists)
  n_missing  <- total - n_present

  logger::log_warn(
    "Neale GWAS check: If files are not already present, downloading the full set may require ~{estimate_gb} GB and many hours."
  )
  logger::log_info("Neale GWAS check: {total} unique files to verify in '{normalizePath(neale_dir, mustWork = FALSE)}' (present={n_present}, missing={n_missing})")

  # --- preflight confirmation ---
  if (n_missing > 0) {
    proceed <- FALSE
    if (confirm == "yes") {
      proceed <- TRUE
    } else if (confirm == "no") {
      proceed <- FALSE
    } else {  # "ask"
      if (interactive()) {
        ans <- readline(prompt = sprintf(
          "About to download up to ~%d GB across %d missing files into '%s'. Proceed? [y/N]: ",
          as.integer(estimate_gb), n_missing, normalizePath(neale_dir, mustWork = FALSE)
        ))
        proceed <- grepl("^[Yy]", ans)
      } else {
        stop("Downloads required but session is non-interactive and confirm='ask'. ",
             "Re-run with confirm='yes' to proceed, or 'no' to skip.", call. = FALSE)
      }
    }
    if (!proceed) {
      logger::log_warn("User declined downloads; skipping. (present={n_present}, missing={n_missing})")
      return(invisible(list(n_total = total, n_present = n_present, n_downloaded = 0L,
                            n_failed = n_missing, aborted = TRUE)))
    }
  } else {
    logger::log_info("All required files already present; no downloads needed.")
    return(invisible(list(n_total = total, n_present = n_present, n_downloaded = 0L,
                          n_failed = 0L, aborted = FALSE)))
  }

  # --- downloads (only for missing) ---
  n_downloaded <- 0L; n_failed <- 0L
  min_size  <- 1024 * 1024     # 1 MB sanity
  retries   <- 3L

  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)

  for (k in seq_len(total)) {
    f   <- basename(uniq_files[k])
    url <- uniq_urls[k]
    dest <- file.path(neale_dir, f)

    .log_progress <- function() {
      logger::log_info("Progress: processed {k}/{total}; present={n_present}, downloaded={n_downloaded}, failed={n_failed}; remaining={total - k}")
    }

    if (file.exists(dest)) {
      utils::setTxtProgressBar(pb, k); .log_progress(); next
    }

    if (!nzchar(url) || is.na(url)) {
      n_failed <- n_failed + 1L
      logger::log_warn("[{k}/{total}] missing AWS URL for {f}; skipping")
      utils::setTxtProgressBar(pb, k); .log_progress(); next
    }

    logger::log_info("[{k}/{total}] downloading: {f}")
    ok <- FALSE
    for (r in seq_len(retries)) {
      ok <- tryCatch({
        utils::download.file(url, destfile = dest, mode = "wb", quiet = !verbose)
        TRUE
      }, error = function(e) {
        logger::log_warn("Attempt {r}/{retries} failed for {f}: {conditionMessage(e)}")
        FALSE
      })
      if (ok && file.exists(dest)) {
        sz <- suppressWarnings(file.info(dest)$size)
        if (!is.na(sz) && sz >= min_size) break
        logger::log_warn("Downloaded {f} but size looks small ({format(sz, big.mark = ',')} bytes); retrying ({r}/{retries})")
        ok <- FALSE
      }
      if (!ok && r < retries) Sys.sleep(2^r)
    }

    if (ok && file.exists(dest)) {
      n_downloaded <- n_downloaded + 1L
      logger::log_info("Downloaded: {f} ({format(file.info(dest)$size, big.mark = ',')} bytes)")
    } else {
      if (file.exists(dest)) try(unlink(dest), silent = TRUE)  # clean partial
      n_failed <- n_failed + 1L
      logger::log_warn("Failed to download: {f}")
    }


    utils::setTxtProgressBar(pb, k); .log_progress()
  }

  logger::log_info("Neale GWAS check complete: total={total}, present={n_present}, downloaded={n_downloaded}, failed/skipped={n_failed}")
  invisible(list(n_total = total, n_present = n_present, n_downloaded = n_downloaded,
                 n_failed = n_failed, aborted = FALSE))
}
#' @export
neale_tbi_maker <- function(...) { ... } #unimplemented feature

#' Fetch outcome SNP rows from Neale sumstats for each outcome (no tabix; exact variant matches)
#'
#' @param exposure_snps MAPPED exposure SNPs. Must include:
#'   - rsid
#'   - neale_variant  (string "chr:pos:ref:alt")
#'   - effect_allele.exposure (preferred) or effect_allele (for logs / future harmonization)
#' @param MR_df Outcome plan (Neale rows). Must include:
#'   - File         (bgz filename in neale_dir)
#'   - description  (preferred phenotype label; else ICD10_explo)
#' @param neale_dir Directory containing Neale .tsv.bgz files
#' @param cache_dir Cache root (default: [ardmr_cache_dir()]); results cached per outcome
#' @param verbose Logical
#' @param force_refresh Logical; if TRUE, delete per-outcome cached RDS before reading
#'   (default inherits cfg$force_refresh if available, else FALSE).
#' @return MR_df with list-column `outcome_snps` (TwoSampleMR-formatted if available)
#' @export
neale_snp_grabber <- function(
    exposure_snps, MR_df, neale_dir,
    cache_dir = ardmr_cache_dir(),
    verbose = TRUE,
    force_refresh = get0("cfg", inherits = TRUE, ifnotfound = list(force_refresh = FALSE))$force_refresh
) {
  stopifnot(is.data.frame(exposure_snps), is.data.frame(MR_df))

  # ---- strict input checks ---------------------------------------------------
  need_exp <- c("rsid","neale_variant")
  miss_exp <- setdiff(need_exp, names(exposure_snps))
  if (length(miss_exp)) {
    stop("neale_snp_grabber(): exposure_snps missing required columns: ", paste(miss_exp, collapse=", "), call. = FALSE)
  }
  need_mr <- c("File")
  miss_mr <- setdiff(need_mr, names(MR_df))
  if (length(miss_mr)) {
    stop("neale_snp_grabber(): MR_df missing required columns: ", paste(miss_mr, collapse=", "), call. = FALSE)
  }

  # ---- helpers ---------------------------------------------------------------
  .slug <- function(x) {
    x <- as.character(x); x[is.na(x) | !nzchar(x)] <- "NA"
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    substr(x, 1, 120)
  }
  .parse_variant <- function(v) {
    # Expect "chr:pos:ref:alt"; return tibble(chr, pos, ref, alt)
    parts <- strsplit(v, ":", fixed = TRUE)
    ok <- vapply(parts, length, integer(1)) == 4L
    if (!all(ok)) {
      bad <- paste(utils::head(v[!ok], 3), collapse = ", ")
      stop("Variant parse error: expected 'chr:pos:ref:alt'. Examples of bad values: ", bad, call. = FALSE)
    }
    tibble::tibble(
      chr = suppressWarnings(as.integer(vapply(parts, `[[`, "", 1))),
      pos = suppressWarnings(as.integer(vapply(parts, `[[`, "", 2))),
      ref = toupper(vapply(parts, `[[`, "", 3)),
      alt = toupper(vapply(parts, `[[`, "", 4))
    )
  }
  .to_bool <- function(x) {
    tx <- tolower(as.character(x))
    ifelse(tx %in% c("true","t","1"), TRUE,
           ifelse(tx %in% c("false","f","0"), FALSE, NA))
  }

  # ---- cache per outcome -----------------------------------------------------
  sex_sub <- if ("Sex" %in% names(MR_df)) {
    s <- unique(tolower(na.omit(MR_df$Sex)))
    if (length(s) == 1 && nzchar(s)) s else "unknown"
  } else "unknown"

  cache_root <- file.path(cache_dir, "neale_outcome_snps", sex_sub)
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)

  # (optional) log once if forcing refresh
  if (isTRUE(force_refresh)) {
    logger::log_info("Neale SNP grabber: force_refresh = TRUE (rebuilding cached per-outcome SNP tables).")
  }

  # ---- build variant lookup once (rsid <-> neale_variant) --------------------
  exp <- tibble::as_tibble(exposure_snps) |>
    dplyr::transmute(
      rsid = .normalize_rsid(.data$rsid),
      neale_variant = trimws(as.character(.data$neale_variant))
    ) |>
    dplyr::filter(!is.na(.data$rsid) & .data$rsid != "" &
                    !is.na(.data$neale_variant) & .data$neale_variant != "") |>
    dplyr::distinct(neale_variant, .keep_all = TRUE)

  if (nrow(exp) == 0) {
    logger::log_warn("Neale SNP grabber: no usable exposure SNPs with 'rsid' + 'neale_variant'. Returning empty results.")
    MR_df$outcome_snps <- replicate(nrow(MR_df), tibble::tibble(), simplify = FALSE)
    return(MR_df)
  }

  # normalized, de-duplicated targets (and keep mapping to rsid)
  target_set <- unique(trimws(as.character(exp$neale_variant)))
  names(target_set) <- exp$rsid

  # ---- iterate files/outcomes -----------------------------------------------
  MR_df$outcome_snps <- vector("list", nrow(MR_df))

  pb <- utils::txtProgressBar(min = 0, max = nrow(MR_df), style = 3)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)

  for (i in seq_len(nrow(MR_df))) {
    utils::setTxtProgressBar(pb, i)

    rec <- MR_df[i, , drop = FALSE]
    fname <- as.character(rec$File[[1]])
    if (!nzchar(fname) || is.na(fname)) {
      logger::log_warn("Neale row {i}: empty 'File' field; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      next
    }

    fpath <- file.path(neale_dir, basename(fname))
    label <- if ("description" %in% names(MR_df) && isTRUE(nzchar(rec$description[[1]]))) {
      as.character(rec$description[[1]])
    } else if ("ICD10_explo" %in% names(MR_df) && isTRUE(nzchar(rec$ICD10_explo[[1]]))) {
      as.character(rec$ICD10_explo[[1]])
    } else {
      tools::file_path_sans_ext(basename(fname))
    }

    cache_file <- file.path(cache_root, paste0(.slug(label), ".rds"))

    # optional hard reset for this phenotype
    if (isTRUE(force_refresh) && file.exists(cache_file)) {
      if (verbose) logger::log_info("Neale row {i}: force_refresh -> deleting cache {basename(cache_file)}")
      try(unlink(cache_file), silent = TRUE)
    }

    # reuse from cache if present (only when NOT forcing refresh)
    if (!isTRUE(force_refresh) && file.exists(cache_file)) {
      if (verbose) logger::log_info("Neale row {i}: using cached result {basename(cache_file)}")
      MR_df$outcome_snps[[i]] <- readRDS(cache_file)
      next
    }

    if (!file.exists(fpath)) {
      logger::log_warn("Neale row {i}: file missing => {basename(fpath)}; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      next
    }

    # ---- read header, enforce required columns ------------------------------
    hdr <- tryCatch(data.table::fread(fpath, nrows = 0, showProgress = FALSE),
                    error = function(e) e)
    if (inherits(hdr, "error")) {
      logger::log_warn("Neale row {i}: failed to read header for {basename(fpath)}: {conditionMessage(hdr)}; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      next
    }
    cols <- names(hdr)
    req <- c("variant","minor_allele","minor_AF","beta","se","pval","low_confidence_variant")
    missing <- setdiff(req, cols)
    if (length(missing)) {
      stop(sprintf("Neale row %d: required columns missing in %s: %s. Available: %s",
                   i, basename(fpath), paste(missing, collapse=", "), paste(cols, collapse=", ")), call. = FALSE)
    }

    # ---- read only required columns; filter to target variants --------------
    dt <- tryCatch(
      data.table::fread(
        fpath,
        select     = req,
        colClasses = c(variant = "character", minor_allele = "character"),
        showProgress = FALSE
      ),
      error = function(e) e
    )
    if (inherits(dt, "error")) {
      logger::log_warn("Neale row {i}: fread failed for {basename(fpath)}: {conditionMessage(dt)}; skipping")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      next
    }

    # --- make sure the key columns exist and are clean ---
    if (!"variant" %in% names(dt)) {
      stop("Internal: 'variant' column missing after fread for ", basename(fpath), call. = FALSE)
    }
    dt[["variant"]]      <- trimws(gsub("\r$", "", as.character(dt[["variant"]])))
    dt[["minor_allele"]] <- toupper(trimws(as.character(dt[["minor_allele"]])))

    # base-safe subsetting (works for data.table or data.frame)
    dt <- dt[dt[["variant"]] %in% target_set, , drop = FALSE]

    if (nrow(dt) == 0) {
      if (verbose) logger::log_info("Neale row {i}: 0 matches in {basename(fpath)} for {length(target_set)} requested variants")
      MR_df$outcome_snps[[i]] <- tibble::tibble()
      try(saveRDS(MR_df$outcome_snps[[i]], cache_file), silent = TRUE)
      next
    }

    # ---- parse variant; compute EAF for ALT (effect allele) -----------------
    var_parsed <- .parse_variant(dt$variant)

    # join rsid back via sanitized neale_variant
    rs_lu <- tibble::tibble(variant = target_set, rsid = names(target_set))

    neale_std <- tibble::as_tibble(dt) |>
      dplyr::bind_cols(var_parsed) |>
      dplyr::left_join(rs_lu, by = "variant") |>
      dplyr::transmute(
        SNP            = .data$rsid,
        chr            = .data$chr,
        pos            = .data$pos,
        effect_allele  = toupper(.data$alt),               # beta is w.r.t ALT
        other_allele   = toupper(.data$ref),
        beta           = suppressWarnings(as.numeric(.data$beta)),
        se             = suppressWarnings(as.numeric(.data$se)),
        # effect allele frequency = ALT AF = minor_AF if ALT==minor_allele else 1 - minor_AF
        eaf            = dplyr::case_when(
          is.na(.data$minor_allele) | is.na(.data$minor_AF) ~ NA_real_,
          toupper(.data$minor_allele) == toupper(.data$alt) ~ suppressWarnings(as.numeric(.data$minor_AF)),
          TRUE ~ 1 - suppressWarnings(as.numeric(.data$minor_AF))
        ),
        pval           = suppressWarnings(as.numeric(.data$pval)),
        low_confidence = .to_bool(.data$low_confidence_variant)
      ) |>
      dplyr::arrange(.data$SNP, .data$pval) |>
      dplyr::distinct(.data$SNP, .keep_all = TRUE)

    # add annotation fields
    neale_std$chrpos    <- paste0(neale_std$chr, ":", neale_std$pos)
    neale_std$Phenotype <- label

    # ---- TwoSampleMR formatting (plain p-values) ----------------------------
    if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
      neale_fmt <- TwoSampleMR::format_data(
        dat               = neale_std,
        type              = "outcome",
        snp_col           = "SNP",
        beta_col          = "beta",
        se_col            = "se",
        effect_allele_col = "effect_allele",
        other_allele_col  = "other_allele",
        eaf_col           = "eaf",
        pval_col          = "pval",
        chr_col           = "chr",
        pos_col           = "pos",
        log_pval          = FALSE,
        phenotype_col     = "Phenotype"
      )
      neale_fmt$chrpos         <- neale_std$chrpos
      neale_fmt$low_confidence <- neale_std$low_confidence
      out_tbl <- tibble::as_tibble(neale_fmt)
    } else {
      out_tbl <- tibble::as_tibble(
        dplyr::select(neale_std, SNP, chr, pos, effect_allele, other_allele, beta, se, eaf, pval, chrpos, low_confidence, Phenotype)
      )
    }

    # ---- sanity: p vs beta/se consistency (warn, don't stop) ----------------
    p2s   <- suppressWarnings(as.numeric(out_tbl$pval))
    z_p   <- stats::qnorm(p2s/2, lower.tail = FALSE)
    z_bse <- abs(suppressWarnings(as.numeric(out_tbl$beta)) /
                   suppressWarnings(as.numeric(out_tbl$se)))
    bad <- which(is.finite(z_p) & is.finite(z_bse) & abs(z_p - z_bse) > 3)
    if (length(bad)) {
      logger::log_warn("Neale row {i} ('{label}'): {length(bad)} SNPs show beta/se vs p inconsistency; keeping but flagging.")
    }

    MR_df$outcome_snps[[i]] <- out_tbl
    try(saveRDS(out_tbl, cache_file), silent = TRUE)

    if (verbose) {
      logger::log_info("Neale row {i}: {nrow(out_tbl)} SNP rows for '{label}' from {basename(fpath)}")
    }
  }

  total_rows <- sum(vapply(MR_df$outcome_snps, nrow, integer(1)))
  logger::log_info("Neale SNP grabber: {nrow(MR_df)} outcomes processed; {total_rows} SNP rows harvested")
  MR_df
}
