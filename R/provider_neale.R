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
