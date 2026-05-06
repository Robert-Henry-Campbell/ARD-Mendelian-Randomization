# R/cache_management.R
# Helpers to inspect and selectively prune the ardmr cache directory.

.cache_categories <- function() {
  list(
    panukb_outcome_snps = list(path = "panukb_outcome_snps", pattern = "\\.rds$"),
    neale_outcome_snps  = list(path = "neale_outcome_snps",  pattern = "\\.rds$"),
    coloc_exposure      = list(path = "coloc/exposure",      pattern = "\\.rds$"),
    coloc_outcome       = list(path = "coloc/outcome",       pattern = "\\.rds$"),
    coloc_meta          = list(path = "coloc/meta",          pattern = "\\.rds$"),
    ld_matrix           = list(path = "ld",                  pattern = "\\.rds$"),
    ld_reference        = list(path = "ld_reference",        pattern = "\\.(bed|bim|fam|tgz)$"),
    tabix_index         = list(path = "tabix_index",         pattern = "\\.tbi$"),
    hts_cache           = list(path = "hts-cache",           pattern = ""),
    outputs             = list(path = "output",              pattern = "")
  )
}

.cache_files <- function(cache_dir, category) {
  cats <- .cache_categories()
  if (!category %in% names(cats)) {
    stop("Unknown cache category: ", category, call. = FALSE)
  }
  base <- file.path(cache_dir, cats[[category]]$path)
  if (!dir.exists(base)) return(character())
  list.files(base, pattern = cats[[category]]$pattern,
             recursive = TRUE, full.names = TRUE)
}

#' Summarise the ardmr cache.
#'
#' @param cache_dir Cache root.
#' @return Tibble with one row per category: entries, size_mb, oldest, newest.
#' @export
ardmr_cache_summary <- function(cache_dir = ardmr_cache_dir()) {
  cats <- names(.cache_categories())
  rows <- lapply(cats, function(cat) {
    files <- .cache_files(cache_dir, cat)
    if (!length(files)) {
      return(tibble::tibble(category = cat, entries = 0L, size_mb = 0,
                            oldest = as.POSIXct(NA), newest = as.POSIXct(NA)))
    }
    info <- file.info(files)
    tibble::tibble(
      category = cat,
      entries  = length(files),
      size_mb  = round(sum(info$size, na.rm = TRUE) / 1024^2, 2),
      oldest   = min(info$mtime, na.rm = TRUE),
      newest   = max(info$mtime, na.rm = TRUE)
    )
  })
  dplyr::bind_rows(rows)
}

#' Selectively prune cached files.
#'
#' By default this is a dry-run (returns the paths it *would* delete).
#' Set `dry_run = FALSE` to actually delete.
#'
#' @param cache_dir Cache root.
#' @param category One of the names in [ardmr_cache_summary()] or `"all"`.
#'   Required.
#' @param exposure_id If supplied, restrict to coloc files for this id.
#' @param ancestry If supplied, restrict outcome-SNP files to this ancestry
#'   (`panukb_outcome_snps` only).
#' @param run_hash If supplied, restrict outputs to this run hash.
#' @param dry_run Logical; default `TRUE`. If `TRUE`, only list candidates.
#' @return Character vector of paths deleted (or, if `dry_run`, would-delete).
#' @export
ardmr_cache_clean <- function(cache_dir = ardmr_cache_dir(),
                              category = NULL,
                              exposure_id = NULL,
                              ancestry = NULL,
                              run_hash = NULL,
                              dry_run = TRUE) {
  if (is.null(category)) {
    stop("`category` is required (one of: ",
         paste(c(names(.cache_categories()), "all"), collapse = ", "),
         ", or 'all').", call. = FALSE)
  }
  if (identical(category, "all")) {
    if (isTRUE(dry_run)) {
      return(unlist(lapply(names(.cache_categories()),
                           function(c) .cache_files(cache_dir, c)), use.names = FALSE))
    }
    cats <- names(.cache_categories())
  } else {
    if (!category %in% names(.cache_categories())) {
      stop("Unknown category: ", category, call. = FALSE)
    }
    cats <- category
  }

  candidates <- character()
  for (cat in cats) {
    files <- .cache_files(cache_dir, cat)
    if (!length(files)) next
    if (cat == "panukb_outcome_snps" && !is.null(ancestry)) {
      files <- files[grepl(sprintf("/%s/", ancestry), files, fixed = TRUE)]
    }
    if (cat == "outputs" && !is.null(run_hash)) {
      files <- files[grepl(sprintf("/%s(/|$)", run_hash), files)]
    }
    if (cat %in% c("coloc_exposure", "coloc_meta") && !is.null(exposure_id)) {
      keep <- vapply(files, function(p) {
        v <- tryCatch(readRDS(p), error = function(e) NULL)
        if (is.null(v)) FALSE
        else identical(attr(v, "exposure_id", exact = TRUE), exposure_id)
      }, logical(1))
      files <- files[keep]
    }
    candidates <- c(candidates, files)
  }
  if (isTRUE(dry_run)) return(candidates)
  for (f in candidates) try(unlink(f, force = TRUE), silent = TRUE)
  if (identical(category, "outputs") || identical(category, "all")) {
    out_dir <- file.path(cache_dir, "output")
    if (dir.exists(out_dir)) {
      empty_dirs <- list.dirs(out_dir, recursive = TRUE, full.names = TRUE)
      empty_dirs <- empty_dirs[order(-nchar(empty_dirs))]
      for (d in empty_dirs) {
        if (!length(list.files(d, all.files = TRUE, no.. = TRUE))) {
          try(unlink(d, recursive = TRUE), silent = TRUE)
        }
      }
    }
  }
  candidates
}
