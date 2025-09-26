#' Build instruments from IEU GWAS, validate units, assemble groups, and run ard_compare()
#'
#' @param csv_path Path to a CSV with columns: ieu_id, exposure, exposure_group,
#'   sex, ancestry, multiple_testing_correction, confirm, and optionally
#'   exposure_units and phenoscanner_exclusions
#' @param cache_dir Output/cache directory used by ard_compare/run_phenome_mr
#' @param jwt IEU JWT string; defaults to IEU_JWT/OPENGWAS_JWT env vars
#' @param prompt_for_units If TRUE, interactively prompt for missing units; otherwise stop on missing units
#' @param p_threshold Genome-wide significance threshold for instrument extraction (default 5e-8)
#' @param phenoscanner_pval P-value threshold used when querying Phenoscanner for
#'   exclusion filtering (default 5e-8)
#' @param r2 LD clumping r^2 threshold (default 0.001)
#' @param kb Clumping window in kb (default 10000)
#' @param f_threshold Per-variant F-stat minimum (default 10)
#' @param force_refresh If TRUE, force ard_compare() to refresh cached results
#'
#' @return (invisible) list with the processed exposure group names, any
#'   failures, the combined PhenoScanner summary tibble, and the written
#'   PhenoScanner summary file paths.
#' @export
run_ieugwasr_ard_compare <- function(
    csv_path,
    cache_dir      = ardmr_cache_dir(),

    jwt            = {
      jwt_env <- Sys.getenv("IEU_JWT", unset = "")
      if (!nzchar(jwt_env)) jwt_env <- Sys.getenv("OPENGWAS_JWT", unset = "")
      jwt_env
    },

    prompt_for_units = interactive(),
    p_threshold    = 5e-8,
    phenoscanner_pval = 5e-8,
    r2             = 0.001,
    kb             = 10000,
    f_threshold    = 10,
    force_refresh  = TRUE,
    sensitivity_pass_min = 6
) {
  # ---- packages ----
  pkgs <- c(
    "TwoSampleMR","ieugwasr","dplyr","readr","purrr",
    "stringr","tibble","MendelianRandomization","jsonlite"
  )
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Please install required packages: ", paste(miss, collapse = ", "))

  # ---- auth ----
  if (!nzchar(jwt)) stop("IEU_JWT or OPENGWAS_JWT env var not set and no 'jwt' provided.")
  Sys.setenv(OPENGWAS_JWT = jwt)

  # ---- helpers ----
  is_rsid <- function(x) !is.na(x) & grepl("^rs\\d+$", x)

  palindrome_flag <- function(a1, a2) {
    a1 <- toupper(a1); a2 <- toupper(a2)
    (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C")
  }

  f_stat <- function(beta, se) (beta/se)^2

  safe_chr <- function(x) {
    x <- trimws(as.character(x))
    x[!nzchar(x) | is.na(x)] <- NA_character_
    x
  }

  slugify <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "NA"
    x <- gsub("[^[:alnum:]]+", "-", x)
    x <- gsub("-+", "-", x)
    x <- gsub("^-|-$", "", x)
    tolower(x)
  }

  parse_phenoscanner_patterns <- function(val) {
    val <- safe_chr(val)
    val <- val[!is.na(val) & nzchar(val)]
    if (!length(val)) {
      return(list(skip = TRUE, patterns = list(), raw = character(0)))
    }
    val <- val[1]
    parts <- unlist(strsplit(val, ","), use.names = FALSE)
    parts <- stringr::str_trim(parts)
    parts <- unique(parts[nzchar(parts)])
    if (!length(parts)) {
      return(list(skip = TRUE, patterns = list(), raw = character(0)))
    }
    list(
      skip = FALSE,
      patterns = lapply(parts, function(p) stringr::regex(p, ignore_case = TRUE)),
      raw = parts
    )
  }

  match_ci_col <- function(df, choices) {
    if (!is.data.frame(df) || !ncol(df)) return(NA_character_)
    nm <- names(df)
    nm_lower <- tolower(nm)
    for (choice in choices) {
      idx <- which(nm_lower == tolower(choice))[1]
      if (length(idx) && !is.na(idx)) {
        return(nm[idx])
      }
    }
    NA_character_
  }

  phenoscanner_query_batches <- function(rsids, pval, exposure_label, ieu_id, max_attempts = 3, base_sleep = 1) {
    rsids <- safe_chr(rsids)
    rsids <- unique(rsids[!is.na(rsids)])
    if (!length(rsids)) {
      return(list(per_snp = list(), errors = list()))
    }

    chunks <- split(rsids, ceiling(seq_along(rsids) / 100))
    per_snp <- stats::setNames(vector("list", length(rsids)), rsids)
    errors <- list()

    for (i in seq_along(chunks)) {
      chunk <- chunks[[i]]
      attempt <- 1L
      success <- FALSE
      last_error <- NULL

      while (attempt <= max_attempts && !success) {
        res <- try(
          MendelianRandomization::phenoscanner(
            snpquery = chunk,
            catalogue = "GWAS",
            pvalue = pval,
            proxies = "None"
          ),
          silent = TRUE
        )

        if (!inherits(res, "try-error")) {
          success <- TRUE
          assoc_tbl <- tryCatch(
            tibble::as_tibble(res$results),
            error = function(e) tibble::tibble()
          )
          snp_col <- match_ci_col(assoc_tbl, c("snp", "rsid", "rs"))
          if (is.na(snp_col) && ncol(assoc_tbl) >= 1) {
            snp_col <- names(assoc_tbl)[1]
          }

          if (nrow(assoc_tbl) > 0) {
            message(sprintf(
              "[PhenoScanner][%s|%s] Chunk %d/%d returned %d association rows.",
              exposure_label, ieu_id, i, length(chunks), nrow(assoc_tbl)
            ))
          } else {
            message(sprintf(
              "[PhenoScanner][%s|%s] Chunk %d/%d returned no hits.",
              exposure_label, ieu_id, i, length(chunks)
            ))
          }

          split_tbl <- if (!is.na(snp_col) && nrow(assoc_tbl)) {
            split(assoc_tbl, assoc_tbl[[snp_col]])
          } else {
            list()
          }

          for (snp in chunk) {
            if (length(split_tbl)) {
              tbl <- split_tbl[[snp]]
              if (!is.null(tbl)) {
                per_snp[[snp]] <- list(status = "ok", table = tibble::as_tibble(tbl))
              } else {
                per_snp[[snp]] <- list(status = "no_hits", table = tibble::tibble())
              }
            } else {
              per_snp[[snp]] <- list(status = "no_hits", table = tibble::tibble())
            }
          }
        } else {
          err_cond <- attr(res, "condition")
          if (inherits(err_cond, "condition")) {
            last_error <- conditionMessage(err_cond)
          } else {
            last_error <- as.character(res)
          }
          message(sprintf(
            "[PhenoScanner][%s|%s] Chunk %d/%d attempt %d failed: %s",
            exposure_label, ieu_id, i, length(chunks), attempt, last_error
          ))
          if (attempt < max_attempts) {
            Sys.sleep(base_sleep * 2^(attempt - 1))
          }
          attempt <- attempt + 1L
        }
      }

      if (!success) {
        message(sprintf(
          "[PhenoScanner][%s|%s] Chunk %d/%d failed after %d attempts.",
          exposure_label, ieu_id, i, length(chunks), max_attempts
        ))
        for (snp in chunk) {
          per_snp[[snp]] <- list(status = "failed", table = NULL, error = last_error)
          errors[[snp]] <- last_error
        }
      }
    }

    list(per_snp = per_snp, errors = errors)
  }

  write_phenoscanner_report <- function(path, tbl) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(tbl, path)
    path
  }

  collapse_unique <- function(x, sep = " | ") {
    vals <- safe_chr(x)
    vals <- unique(vals[!is.na(vals) & nzchar(vals)])
    if (!length(vals)) {
      return(NA_character_)
    }
    paste(vals, collapse = sep)
  }

  apply_phenoscanner_filters <- function(snps_df, exposure_label, ieu_id, sex, ancestry,
                                         patterns_info, phenoscanner_pval, cache_dir) {
    total_snps <- nrow(snps_df)
    summary <- list(
      ieu_id = ieu_id,
      exposure = exposure_label,
      sex = as.character(sex),
      ancestry = as.character(ancestry),
      total_snps = total_snps,
      matches_removed = 0L,
      failures = 0L,
      survivors = total_snps,
      patterns = if (length(patterns_info$raw)) paste(patterns_info$raw, collapse = ", ") else NA_character_
    )

    summary$phenoscanner_skipped <- FALSE
    summary$skip_reason <- NA_character_

    output_dir <- file.path(cache_dir, "output", slugify(exposure_label), sex, ancestry)
    pre_snapshot_path <- write_phenoscanner_report(
      file.path(output_dir, "phenoscanner_pre_filter.csv"),
      snps_df
    )
    summary$pre_filter_snapshot <- pre_snapshot_path

    if (!total_snps || isTRUE(patterns_info$skip)) {
      reasons <- character(0)
      if (!total_snps) reasons <- c(reasons, "no_snps")
      if (isTRUE(patterns_info$skip)) reasons <- c(reasons, "no_patterns")
      if (length(reasons)) summary$skip_reason <- paste(reasons, collapse = "; ")
      summary$phenoscanner_skipped <- TRUE
      if (isTRUE(patterns_info$skip)) {
        message(sprintf(
          "[PhenoScanner][%s|%s] No exclusion patterns supplied; skipping PhenoScanner filter.",
          exposure_label, ieu_id
        ))
      }
      summary$report_path <- NA_character_
      return(list(snps = snps_df, summary = summary, report_path = NA_character_))
    }

    query_res <- phenoscanner_query_batches(snps_df$SNP, phenoscanner_pval, exposure_label, ieu_id)
    per_snp <- query_res$per_snp
    if (!length(per_snp)) {
      per_snp <- stats::setNames(vector("list", length(snps_df$SNP)), snps_df$SNP)
      for (snp in snps_df$SNP) {
        per_snp[[snp]] <- list(status = "no_hits", table = tibble::tibble())
      }
    }

    snp_order <- snps_df$SNP
    status_vec <- character(length(snp_order))
    remove_vec <- logical(length(snp_order))
    matched_pattern_vec <- rep(NA_character_, length(snp_order))
    assoc_json_vec <- rep(NA_character_, length(snp_order))

    failures <- character(0)

    for (idx in seq_along(snp_order)) {
      snp <- snp_order[idx]
      info <- per_snp[[snp]]
      if (is.null(info)) {
        info <- list(status = "no_hits", table = tibble::tibble())
      }
      status <- info$status
      if (length(status) > 1) status <- status[1]
      if (is.null(status) || !nzchar(status)) status <- "no_hits"
      status_vec[idx] <- status

      if (identical(status, "failed")) {
        failures <- c(failures, snp)
        message(sprintf(
          "[PhenoScanner][%s|%s] Query failed for %s after retries.",
          exposure_label, ieu_id, snp
        ))
        assoc_json_vec[idx] <- NA_character_
        next
      }

      tbl <- info$table
      if (!inherits(tbl, "tbl_df")) {
        tbl <- tryCatch(tibble::as_tibble(tbl), error = function(e) tibble::tibble())
      }

      p_col <- match_ci_col(tbl, c("p", "pvalue", "p_gwas"))
      if (!is.na(p_col)) {
        pvals <- suppressWarnings(as.numeric(tbl[[p_col]]))
        keep_idx <- is.na(pvals) | pvals <= phenoscanner_pval
        tbl <- tbl[keep_idx, , drop = FALSE]
      }

      trait_col <- match_ci_col(tbl, c("trait", "phenotype", "outcome"))
      traits <- if (!is.na(trait_col)) safe_chr(tbl[[trait_col]]) else character(0)
      traits <- traits[!is.na(traits)]

      matched_flags <- if (length(traits) && length(patterns_info$patterns)) {
        vapply(seq_along(patterns_info$patterns), function(j) {
          any(stringr::str_detect(traits, patterns_info$patterns[[j]]), na.rm = TRUE)
        }, logical(1))
      } else {
        rep(FALSE, length(patterns_info$patterns))
      }

      if (length(traits) == 0 && !identical(status, "no_hits")) {
        message(sprintf(
          "[PhenoScanner][%s|%s] %s returned associations but no trait column detected.",
          exposure_label, ieu_id, snp
        ))
      }

      if (any(matched_flags)) {
        first_match <- patterns_info$raw[which(matched_flags)[1]]
        matched_pattern_vec[idx] <- first_match
        remove_vec[idx] <- TRUE
        message(sprintf(
          "[PhenoScanner][%s|%s] Forbidden trait match (%s) detected for %s — marked for removal.",
          exposure_label, ieu_id, first_match, snp
        ))
      } else if (identical(status, "no_hits") || !nrow(tbl)) {
        message(sprintf(
          "[PhenoScanner][%s|%s] %s has no forbidden traits (no hits).",
          exposure_label, ieu_id, snp
        ))
      } else {
        message(sprintf(
          "[PhenoScanner][%s|%s] %s has no forbidden traits; keeping.",
          exposure_label, ieu_id, snp
        ))
      }

      assoc_json_vec[idx] <- jsonlite::toJSON(tbl, dataframe = "rows", auto_unbox = TRUE, na = "null")
    }

    summary$failures <- length(unique(failures))

    if (length(failures)) {
      for (snp in failures) {
        decision <- "keep"
        if (interactive()) {
          repeat {
            ans <- readline(sprintf(
              "PhenoScanner query failed for %s. Keep (k), drop (d), or abort (a)? ",
              snp
            ))
            ans <- tolower(trimws(ans))
            if (!nzchar(ans) || ans %in% c("k", "keep")) {
              decision <- "keep"
              break
            }
            if (ans %in% c("d", "drop")) {
              decision <- "drop"
              break
            }
            if (ans %in% c("a", "abort")) {
              stop("User aborted due to PhenoScanner query failure.")
            }
            message("Please enter 'k', 'd', or 'a'.")
          }
        } else {
          message(sprintf(
            "[PhenoScanner][%s|%s] Non-interactive mode; keeping %s despite failed query.",
            exposure_label, ieu_id, snp
          ))
        }

        idx <- match(snp, snp_order)
        status_vec[idx] <- "failed"
        assoc_json_vec[idx] <- NA_character_

        if (identical(decision, "drop")) {
          remove_vec[idx] <- TRUE
          message(sprintf(
            "[PhenoScanner][%s|%s] %s dropped after user confirmation despite query failure.",
            exposure_label, ieu_id, snp
          ))
        } else {
          remove_vec[idx] <- FALSE
          message(sprintf(
            "[PhenoScanner][%s|%s] %s retained after query failure.",
            exposure_label, ieu_id, snp
          ))
        }
      }
    }

    survivors <- snps_df[!snps_df$SNP %in% snp_order[remove_vec], , drop = FALSE]

    summary$matches_removed <- sum(remove_vec)
    summary$survivors <- nrow(survivors)

    output_dir <- file.path(cache_dir, "output", slugify(exposure_label), sex, ancestry)
    report_tbl <- tibble::tibble(
      SNP = snp_order,
      query_status = status_vec,
      matched_pattern = matched_pattern_vec,
      associations_json = assoc_json_vec,
      remove = remove_vec
    )
    report_path <- write_phenoscanner_report(
      file.path(output_dir, "phenoscanner_review.csv"),
      report_tbl
    )
    message(sprintf(
      "[PhenoScanner][%s|%s] Review saved to %s",
      exposure_label, ieu_id, report_path
    ))
    summary$report_path <- report_path

    list(snps = survivors, summary = summary, report_path = report_path)
  }

  # pick first existing column name from a list
  pick_colname <- function(df, candidates, required = TRUE, context = "column") {
    nm <- candidates[candidates %in% names(df)][1]
    if (isTRUE(required) && is.na(nm)) {
      stop("Could not find any of: ", paste(candidates, collapse=", "), " for ", context)
    }
    nm
  }

  normalize_mtc <- function(x) {
    x <- safe_chr(x)
    if (!length(x)) return(NA_character_)
    out <- vapply(x, function(val) {
      if (is.na(val)) stop("multiple_testing_correction cannot be missing in the CSV.")
      lv <- tolower(val)
      if (lv %in% c("bh","fdr")) {
        "BH"
      } else if (lv %in% c("bonferroni","bonf")) {
        "bonferroni"
      } else if (val %in% c("BH","bonferroni")) {
        val
      } else {
        stop(sprintf("Unsupported multiple_testing_correction value '%s'. Use 'BH' or 'bonferroni'.", val))
      }
    }, character(1))
    unname(out)
  }

  normalize_confirm <- function(x) {
    x <- safe_chr(x)
    if (!length(x)) return(NA_character_)
    out <- vapply(x, function(val) {
      if (is.na(val)) stop("confirm cannot be missing in the CSV.")
      lv <- tolower(val)
      if (lv %in% c("ask","yes","no")) lv else stop(sprintf("Unsupported confirm value '%s'. Use 'ask', 'yes', or 'no'.", val))
    }, character(1))
    unname(out)
  }

  # LD-pop mapping for clumping
  ld_pop_from_ancestry <- function(anc) {
    anc0 <- toupper(trimws(ifelse(is.null(anc) || is.na(anc), "", as.character(anc))))
    if (anc0 %in% c("EUR","EUROPEAN")) return("EUR")
    if (anc0 %in% c("AFR","AFRICAN","AFRICAN-UKB")) return("AFR")
    if (anc0 %in% c("EAS","EAST_ASIAN","EAST-ASIAN","EAST ASIAN")) return("EAS")
    if (anc0 %in% c("CSA","SAS","SOUTH_ASIAN","CENTRAL_SOUTH_ASIAN","SOUTH ASIAN")) return("SAS")
    warning(sprintf("Unknown ancestry '%s'; defaulting LD pop to EUR for clumping.", anc))
    "EUR"
  }

  # light retry/backoff for flaky API calls
  with_backoff <- function(expr, tries = 3, base_sleep = 0.5) {
    for (i in seq_len(tries)) {
      out <- try(force(expr), silent = TRUE)
      if (!inherits(out, "try-error")) return(out)
      if (i < tries) Sys.sleep(base_sleep * 2^(i-1))
    }
    out
  }

  # units via ieugwasr::gwasinfo()
  fetch_units <- function(ieu_id) {
    info <- try(ieugwasr::gwasinfo(ieu_id), silent = TRUE)
    if (!inherits(info,"try-error") && "unit" %in% names(info)) {
      u <- safe_chr(info$unit[1]); if (!is.na(u)) return(u)
    }
    NA_character_
  }

  # chr/pos annotator: reuse if present, else variants_rsid(), else associations()
  annotate_chrpos <- function(dat) {
    # reuse if already present
    if (all(c("chr.exposure","pos.exposure") %in% names(dat))) {
      return(dat |>
        dplyr::mutate(
          Chr = as.character(.data$chr.exposure),
          Pos = as.numeric(.data$pos.exposure)
        ))
    }
    if (all(c("chr","pos") %in% names(dat))) {
      return(dat |>
        dplyr::mutate(
          Chr = as.character(.data$chr),
          Pos = as.numeric(.data$pos)
        ))
    }

    if (!"SNP" %in% names(dat)) stop("annotate_chrpos(): 'SNP' column is missing.")
    snps <- unique(dat$SNP)

    ann <- try(ieugwasr::variants_rsid(snps), silent = TRUE)
    if (!inherits(ann, "try-error") && all(c("name","chr","pos") %in% names(ann))) {
      ann <- ann |>
        dplyr::transmute(
          SNP = .data$name,
          Chr = as.character(.data$chr),
          Pos = as.numeric(.data$pos)
        )
      return(dplyr::left_join(dat, ann, by = "SNP"))
    }

    look_one_study <- function(study_id, rs) {
      out <- ieugwasr::associations(variants = rs, id = study_id)
      if (!nrow(out)) return(out[0, ])
      out |>
        dplyr::select(rsid, chr, pos) |>
        dplyr::distinct()
    }

    study_id <- if ("id.exposure" %in% names(dat)) dat$id.exposure[1] else NA_character_
    if (is.na(study_id)) {
      warning("No id.exposure in dat; cannot fetch chr/pos. Returning NA Chr/Pos.")
      return(dat |>
        dplyr::mutate(Chr = NA_character_, Pos = as.numeric(NA)))
    }

    chunks <- split(snps, ceiling(seq_along(snps)/500))
    ann_list <- lapply(chunks, function(rs) try(look_one_study(study_id, rs), silent = TRUE))
    ann_list <- ann_list[!vapply(ann_list, inherits, logical(1), "try-error")]
    if (length(ann_list)) {
      ann2 <- dplyr::bind_rows(ann_list) |>
        dplyr::transmute(
          SNP = .data$rsid,
          Chr = as.character(.data$chr),
          Pos = as.numeric(.data$pos)
        )
      return(dplyr::left_join(dat, ann2, by = "SNP"))
    }

    dat |>
      dplyr::mutate(Chr = NA_character_, Pos = as.numeric(NA))
  }

  build_exposure_snps <- function(raw, exposure_label, ieu_id) {
    # resolve column names that differ by version
    beta_nm <- pick_colname(raw, c("beta","beta.exposure"), TRUE, "beta")
    se_nm   <- pick_colname(raw, c("se","se.exposure"), TRUE, "se")
    ea_nm   <- pick_colname(raw, c("effect_allele","effect_allele.exposure"), TRUE, "effect_allele")
    oa_nm   <- pick_colname(raw, c("other_allele","other_allele.exposure"), TRUE, "other_allele")
    eaf_nm  <- pick_colname(raw, c("eaf","eaf.exposure"), TRUE, "eaf")
    pval_nm <- pick_colname(raw, c("pval","pval.exposure"), TRUE, "pval")

    if (!"SNP" %in% names(raw)) stop("Expected 'SNP' in instrument table.")
    if (any(!is_rsid(raw$SNP))) {
      bad <- unique(raw$SNP[!is_rsid(raw$SNP)])
      stop(sprintf("Non-rsID variant IDs for %s: e.g. %s", ieu_id, paste(utils::head(bad,5), collapse=", ")))
    }

    out <- raw |>
      dplyr::mutate(
        id.exposure            = ieu_id,
        exposure               = exposure_label,
        Exposure               = exposure_label,  # optional alias
        effect_allele.exposure = .data[[ea_nm]],
        other_allele.exposure  = .data[[oa_nm]],
        eaf.exposure           = .data[[eaf_nm]],
        beta.exposure          = .data[[beta_nm]],
        se.exposure            = .data[[se_nm]],
        pval.exposure          = .data[[pval_nm]],
        palindromic            = palindrome_flag(.data[[ea_nm]], .data[[oa_nm]]),
        mr_keep                = TRUE,
        F                      = f_stat(.data[[beta_nm]], .data[[se_nm]])
      ) |>
      dplyr::select(
        id.exposure, exposure, Exposure, SNP, Chr, Pos,
        effect_allele.exposure, other_allele.exposure, eaf.exposure,
        beta.exposure, se.exposure, palindromic, pval.exposure,
        mr_keep, F
      ) |>
      dplyr::arrange(.data$SNP)

    need <- c("id.exposure","exposure","SNP","Chr","Pos",
              "effect_allele.exposure","other_allele.exposure","eaf.exposure",
              "beta.exposure","se.exposure","palindromic","pval.exposure",
              "mr_keep","F")
    stopifnot(all(need %in% names(out)))
    out
  }

  # ---- p-value backoff ladder ----
  p_backoff <- c(p_threshold, 5e-7, 5e-6) |> unique()

  # instrument extraction + F filtering with ancestry-aware fallback + p backoff
  # ALWAYS clumps before returning
  get_instruments <- function(ieu_id, ancestry) {
    ld_pop <- ld_pop_from_ancestry(ancestry)

    # helper: robust tophits across ieugwasr versions (p vs pval)
    safe_tophits <- function(id, thr) {
      th <- try(ieugwasr::tophits(id = id, p = thr), silent = TRUE)
      if (inherits(th, "try-error") || !is.data.frame(th)) {
        th <- try(ieugwasr::tophits(id = id, pval = thr), silent = TRUE)
      }
      if (inherits(th, "try-error") || !is.data.frame(th)) return(NULL)
      if (!nrow(th)) return(th[0, ])
      if (!"p" %in% names(th) && "pval" %in% names(th)) th$p <- th$pval
      th
    }

    # helper: self-clump a table that has rsids + p-values, then subset original
    # helper: self-clump a table that has rsids + p-values, then subset original
    self_clump <- function(tbl, rs_col = "SNP", p_col = NULL) {
      # pick a p-value column if not provided
      if (is.null(p_col)) {
        p_col <- pick_colname(tbl, c("p","pval","pval.exposure","pval.outcome"), TRUE, "p-value")
      }

      # build library for clumping
      lib <- tbl |>
        dplyr::transmute(rsid = .data[[rs_col]], p = .data[[p_col]]) |>
        dplyr::filter(is.finite(.data$p), !is.na(.data$rsid)) |>
        dplyr::distinct()

      # if 0 or 1 SNPs, just return the input (no clumping possible/needed)
      if (!nrow(lib) || nrow(lib) == 1L) {
        return(tbl)
      }

      cl <- with_backoff(ieugwasr::ld_clump(lib, pop = ld_pop, r2 = r2, kb = kb))
      if (inherits(cl, "try-error") || !is.data.frame(cl) || !nrow(cl)) {
        warning(sprintf("ld_clump() failed/empty for %s (pop=%s). Proceeding UNCLUMPED.", ieu_id, ld_pop))
        return(tbl)
      }

      # correct key mapping: keep rows in tbl whose rs_col matches cl$rsid
      dplyr::semi_join(tbl, cl, by = setNames("rsid", rs_col))
    }


    for (pthr in p_backoff) {
      message(sprintf("[ID=%s] Trying instrument extraction at p<=%.0e …", ieu_id, pthr))

      # ---- 1) Standard path: already clumped by API ----
      ins <- with_backoff(
        TwoSampleMR::extract_instruments(outcomes = ieu_id, p1 = pthr, clump = TRUE, r2 = r2, kb = kb)
      )

      # If that fails or returns 0 rows, try unclumped then self-clump
      if (inherits(ins, "try-error") || !nrow(ins)) {
        if (inherits(ins, "try-error")) {
          message(sprintf("[ID=%s] extract_instruments(clump=TRUE) errored at p<=%.0e; retrying clump=FALSE + self-clump …", ieu_id, pthr))
        } else {
          message(sprintf("[ID=%s] No SNPs via extract_instruments(clump=TRUE) at p<=%.0e; retrying clump=FALSE + self-clump …", ieu_id, pthr))
        }
        ins <- with_backoff(
          TwoSampleMR::extract_instruments(outcomes = ieu_id, p1 = pthr, clump = FALSE)
        )
        if (!inherits(ins, "try-error") && nrow(ins)) {
          pcol_try <- if ("pval" %in% names(ins)) "pval" else if ("pval.exposure" %in% names(ins)) "pval.exposure" else NULL
          ins <- self_clump(ins, rs_col = "SNP", p_col = pcol_try)
        }
      }

      # ---- 2) Fallback to tophits + self-clump ----
      if (inherits(ins, "try-error") || !nrow(ins)) {
        message(sprintf("[ID=%s] Falling back to ieugwasr::tophits() at p<=%.0e …", ieu_id, pthr))
        th <- safe_tophits(ieu_id, pthr)
        if (is.null(th)) {
          message(sprintf("[ID=%s] tophits error at p<=%.0e; continuing backoff …", ieu_id, pthr))
          next
        }
        if (!nrow(th)) {
          message(sprintf("[ID=%s] No SNPs at p<=%.0e; continuing backoff …", ieu_id, pthr))
          next
        }
        thc <- self_clump(th, rs_col = "rsid", p_col = "p")
        ins <- dplyr::transmute(
          thc,
          SNP           = .data$rsid,
          effect_allele = .data$ea,
          other_allele  = .data$oa,
          eaf           = .data$eaf,
          beta          = .data$beta,
          se            = .data$se,
          pval          = .data$p
        )
      }

      # ---- 3) F filter (robust to column name variants) ----
      beta_nm <- if ("beta" %in% names(ins)) "beta" else if ("beta.exposure" %in% names(ins)) "beta.exposure" else NA_character_
      se_nm   <- if ("se"   %in% names(ins)) "se"   else if ("se.exposure"   %in% names(ins)) "se.exposure"   else NA_character_
      if (is.na(beta_nm) || is.na(se_nm)) {
        message(sprintf("[ID=%s] Missing beta/se at p<=%.0e; continuing backoff …", ieu_id, pthr))
        next
      }
      ins$F_tmp <- (ins[[beta_nm]] / ins[[se_nm]])^2
      ins <- dplyr::filter(ins, is.finite(.data$F_tmp), .data$F_tmp >= f_threshold) |>
        dplyr::select(-.data$F_tmp)
      if (!nrow(ins)) {
        message(sprintf("[ID=%s] All SNPs dropped by F-stat (>= %g) at p<=%.0e; continuing backoff …", ieu_id, f_threshold, pthr))
        next
      }

      # ---- 4) Annotate chr/pos and return ----
      ins <- annotate_chrpos(ins)
      message(sprintf("[ID=%s] SUCCESS at p<=%.0e: %d clumped SNPs retained (F>=%.1f).", ieu_id, pthr, nrow(ins), f_threshold))
      attr(ins, "p_used") <- pthr
      return(ins)
    }

    # If we reach here, nothing found across the ladder
    message(sprintf("[ID=%s] FAILED after backoff (5e-08 → 5e-07 → 5e-06). Skipping.", ieu_id))
    return(tibble::tibble())
  }

  # ---- read CSV ----
  expos <- readr::read_csv(csv_path, show_col_types = FALSE)
  req_cols <- c(
    "ieu_id","exposure","exposure_group","sex","ancestry",
    "multiple_testing_correction","confirm"
  )
  miss_cols <- setdiff(req_cols, names(expos))
  if (length(miss_cols)) stop("CSV missing required columns: ", paste(miss_cols, collapse = ", "))
  if (!"exposure_units" %in% names(expos)) expos$exposure_units <- NA_character_
  if (!"phenoscanner_exclusions" %in% names(expos)) expos$phenoscanner_exclusions <- NA_character_

  expos <- expos |>
    dplyr::mutate(
      exposure = safe_chr(exposure),
      exposure_group = safe_chr(exposure_group),
      sex = safe_chr(sex),
      ancestry = safe_chr(ancestry),
      exposure_units = safe_chr(exposure_units),
      phenoscanner_exclusions = {
        val <- safe_chr(phenoscanner_exclusions)
        if (length(val)) {
          idx <- !is.na(val)
          val[idx] <- stringr::str_squish(val[idx])
          val[!is.na(val) & !nzchar(val)] <- NA_character_
        }
        val
      },
      multiple_testing_correction = normalize_mtc(multiple_testing_correction),
      confirm = normalize_confirm(confirm)
    )

  if (any(is.na(expos$exposure))) stop("'exposure' column contains missing values.")

  # ---- units (preflight) ----
  unique_ids <- unique(expos$ieu_id)
  units_map <- setNames(rep(NA_character_, length(unique_ids)), unique_ids)
  for (id in unique_ids) {
    csv_units <- unique(expos$exposure_units[expos$ieu_id == id])
    csv_units <- csv_units[!is.na(csv_units) & csv_units != ""]
    if (length(csv_units) > 1L) {
      stop(sprintf("Exposure units differ for ieu_id '%s' in the CSV: %s", id, paste(csv_units, collapse = " | ")))
    }
    if (length(csv_units) == 1L) {
      units_map[[id]] <- csv_units
    } else {
      units_map[[id]] <- fetch_units(id)
    }
  }

  need_units <- names(units_map)[is.na(units_map) | units_map == ""]
  if (length(need_units)) {
    if (!prompt_for_units) {
      stop("Missing units for: ", paste(need_units, collapse = ", "),
           ". Re-run with prompt_for_units=TRUE or provide a units mapping.")
    }
    message("Please enter units for the following IEU GWAS (press Enter to confirm each):")
    for (id in need_units) {
      u <- readline(prompt = sprintf("Units for %s: ", id))
      u <- stringr::str_trim(u)
      if (!nzchar(u)) stop(sprintf("Units for %s are required.", id))
      units_map[[id]] <- u
    }
  }

  # ---- build per-row instrument tables + collect meta ----
  message("Fetching & preparing instruments…")
  extraction_meta <- list()  # one entry per attempted row
  phenoscanner_summaries <- list()

  expos_prepared <- expos |>
    dplyr::mutate(.row_id = dplyr::row_number()) |>
    dplyr::group_split(.row_id, .keep = TRUE) |>
    purrr::map(function(rr) {
      r <- dplyr::slice(rr, 1)
      ieu_id   <- r$ieu_id
      exp_name <- r$exposure

      ins <- get_instruments(ieu_id, ancestry = r$ancestry)
      p_used <- attr(ins, "p_used")
      status <- if (nrow(ins)) "success" else "fail"

      meta_entry <- list(
        ieu_id   = ieu_id,
        exposure = exp_name,
        p_used   = if (length(p_used)) p_used else NA_real_,
        status   = status
      )

      if (!nrow(ins)) {
        extraction_meta[[length(extraction_meta)+1]] <<- meta_entry
        warning(sprintf("No instruments retained for %s; skipping.", ieu_id))
        return(NULL)
      }

      snps_df <- build_exposure_snps(ins, exposure_label = exp_name, ieu_id = ieu_id)

      patterns_info <- parse_phenoscanner_patterns(r$phenoscanner_exclusions)
      filter_res <- apply_phenoscanner_filters(
        snps_df,
        exposure_label = exp_name,
        ieu_id = ieu_id,
        sex = r$sex,
        ancestry = r$ancestry,
        patterns_info = patterns_info,
        phenoscanner_pval = phenoscanner_pval,
        cache_dir = cache_dir
      )

      snps_df <- filter_res$snps
      filter_summary <- filter_res$summary
      filter_summary$exposure_group <- r$exposure_group
      phenoscanner_summaries[[length(phenoscanner_summaries)+1]] <<- filter_summary

      meta_entry$phenoscanner_total <- filter_summary$total_snps
      meta_entry$phenoscanner_removed <- filter_summary$matches_removed
      meta_entry$phenoscanner_failures <- filter_summary$failures
      meta_entry$phenoscanner_survivors <- filter_summary$survivors

      mtc_val     <- normalize_mtc(r$multiple_testing_correction)[1]
      confirm_val <- normalize_confirm(r$confirm)[1]

      if (!nrow(snps_df)) {
        extraction_meta[[length(extraction_meta)+1]] <<- meta_entry
        warning(sprintf(
          "All instruments removed for %s after PhenoScanner filtering; skipping exposure.",
          ieu_id
        ))
        return(NULL)
      }

      extraction_meta[[length(extraction_meta)+1]] <<- meta_entry

      list(
        exposure_group = r$exposure_group,
        sex            = match.arg(as.character(r$sex), c("both","male","female")),
        ancestry       = as.character(r$ancestry),
        exposure       = exp_name,
        exposure_units = units_map[[ieu_id]],
        exposure_snps  = snps_df,
        phenoscanner_exclusions = r$phenoscanner_exclusions,
        phenoscanner_pval = phenoscanner_pval,
        multiple_testing_correction = mtc_val,
        confirm        = confirm_val
      )
    }) |>
    purrr::compact()

  phenoscanner_summary_df <-
    if (length(phenoscanner_summaries)) {
      dplyr::bind_rows(phenoscanner_summaries) |>
        dplyr::mutate(
          total_snps = suppressWarnings(as.integer(.data$total_snps)),
          matches_removed = suppressWarnings(as.integer(.data$matches_removed)),
          failures = suppressWarnings(as.integer(.data$failures)),
          survivors = suppressWarnings(as.integer(.data$survivors)),
          phenoscanner_skipped = as.logical(.data$phenoscanner_skipped),
          skip_reason = safe_chr(.data$skip_reason),
          pre_filter_snapshot = safe_chr(.data$pre_filter_snapshot),
          report_path = safe_chr(.data$report_path)
        )
    } else {
      tibble::tibble()
    }

  # ---- summary + user interrupt BEFORE ard_compare() ----
  meta_df <- dplyr::bind_rows(extraction_meta)
  if (is.null(meta_df) || !nrow(meta_df)) stop("No exposure rows were attempted; nothing to do.")

  meta_df <- meta_df |>
    dplyr::mutate(
      status = as.character(.data$status),
      p_used = suppressWarnings(as.numeric(.data$p_used)),
      phenoscanner_total = suppressWarnings(as.integer(.data$phenoscanner_total)),
      phenoscanner_removed = suppressWarnings(as.integer(.data$phenoscanner_removed)),
      phenoscanner_failures = suppressWarnings(as.integer(.data$phenoscanner_failures)),
      phenoscanner_survivors = suppressWarnings(as.integer(.data$phenoscanner_survivors))
    )

  tried_n <- nrow(meta_df)
  succ_5e8 <- sum(meta_df$status == "success" & !is.na(meta_df$p_used) & meta_df$p_used <= 5e-8)
  succ_5e7 <- sum(meta_df$status == "success" & !is.na(meta_df$p_used) & meta_df$p_used > 5e-8 & meta_df$p_used <= 5e-7)
  succ_5e6 <- sum(meta_df$status == "success" & !is.na(meta_df$p_used) & meta_df$p_used > 5e-7 & meta_df$p_used <= 5e-6)
  failed_n <- sum(meta_df$status == "fail")

  failed_wh <- meta_df |>
    dplyr::filter(.data$status == "fail") |>
    dplyr::distinct(.data$ieu_id, .data$exposure)

  summary_tbl <- tibble::tibble(
    total_attempted   = tried_n,
    extracted_at_5e_08 = succ_5e8,
    extracted_at_5e_07 = succ_5e7,
    extracted_at_5e_06 = succ_5e6,
    failed_total       = failed_n
  )

  message("\n=== Instrument Extraction Summary ===")
  print(summary_tbl)
  if (failed_n) {
    message("\nFailed studies (no instruments after backoff):")
    print(failed_wh, n = Inf)
  }

  if (nrow(phenoscanner_summary_df)) {
    message("\n=== PhenoScanner Filter Summary ===")
    print(phenoscanner_summary_df, n = Inf)
  }

  # user interrupt (interactive sessions only)
  if (interactive()) {
    ans <- readline("Proceed to ard_compare() with the successfully extracted exposures? [y/N]: ")
    if (!nzchar(ans) || !tolower(substr(ans,1,1)) %in% c("y")) {
      stop("User aborted before ard_compare().")
    }
  } else {
    message("Non-interactive session detected; proceeding without prompt.")
  }



  if (!length(expos_prepared)) stop("No valid rows after instrument preparation.")

  # ---- assemble groups and run ard_compare() per exposure_group ----
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  by_group <- split(expos_prepared, vapply(expos_prepared, `[[`, "", "exposure_group"))

  processed <- character(0)
  failed_records <- list()
  phenoscanner_summary_paths <- list()

  for (grp_name in names(by_group)) {
    entries <- by_group[[grp_name]]

    # exposure label (use first; warn if multiple)
    exp_label <- unique(vapply(entries, `[[`, "", "exposure"))
    if (length(exp_label) > 1L) {
      warning(sprintf("Multiple 'exposure' labels in group '%s'. Using the first: %s", grp_name, exp_label[1]))
    }
    exposure_label <- exp_label[1]

    # check units agree across entries
    u_vals <- unique(vapply(entries, `[[`, "", "exposure_units"))
    if (length(u_vals) > 1L) {
      stop(sprintf("Units differ within exposure_group '%s': %s", grp_name, paste(u_vals, collapse = " | ")))
    }
    exposure_units <- u_vals[1]

    mtc_vals <- unique(vapply(entries, `[[`, "", "multiple_testing_correction"))
    if (length(mtc_vals) > 1L) {
      stop(sprintf("multiple_testing_correction values differ within exposure_group '%s': %s",
                   grp_name, paste(mtc_vals, collapse = " | ")))
    }
    multiple_testing_correction <- mtc_vals[1]

    confirm_vals <- unique(vapply(entries, `[[`, "", "confirm"))
    if (length(confirm_vals) > 1L) {
      stop(sprintf("confirm values differ within exposure_group '%s': %s",
                   grp_name, paste(confirm_vals, collapse = " | ")))
    }
    confirm_val <- confirm_vals[1]

    # build groups list
    groups_list <- lapply(entries, function(e) {
      list(
        sex = e$sex,
        ancestry = e$ancestry,
        exposure_snps = e$exposure_snps,
        phenoscanner_exclusions = e$phenoscanner_exclusions,
        phenoscanner_pval = e$phenoscanner_pval
      )
    })
    names(groups_list) <- paste0(
      vapply(entries, `[[`, "", "sex"), "_", vapply(entries, `[[`, "", "ancestry")
    )

    message(sprintf("\n=== ard_compare(): group '%s' ===", grp_name))
    message(sprintf("Exposure: %s | Units: %s", exposure_label, exposure_units))

    result <- tryCatch(
      {
        compare_info <- ard_compare(
          exposure                    = exposure_label,
          exposure_units              = exposure_units,
          groups                      = groups_list,
          sensitivity_enabled         = c(
            "egger_intercept","egger_slope_agreement",
            "weighted_median","weighted_mode",
            "steiger_direction","leave_one_out",
            "ivw_Q","ivw_I2"
          ),
          sensitivity_pass_min        = sensitivity_pass_min,
          Multiple_testing_correction = multiple_testing_correction,
          scatterplot                 = TRUE,
          snpforestplot               = TRUE,
          leaveoneoutplot             = TRUE,
          cache_dir                   = cache_dir,
          logfile                     = NULL,
          verbose                     = TRUE,
          confirm                     = confirm_val,
          force_refresh               = force_refresh
        )

        list(ok = TRUE, compare = compare_info)
      },
      error = function(err) {
        err_msg <- conditionMessage(err)
        warning(
          sprintf(
            "ard_compare() failed for exposure_group '%s' (exposure: %s). Error: %s",
            grp_name,
            exposure_label,
            err_msg
          ),
          call. = FALSE
        )

        list(
          ok = FALSE,
          error_message = err_msg
        )
      }
    )

    if (isTRUE(result$ok)) {
      processed <- c(processed, grp_name)

      compare_root <- NA_character_
      compare_logfile <- NA_character_
      if (is.list(result$compare)) {
        cr <- safe_chr(result$compare$compare_root)
        if (length(cr)) compare_root <- cr[1]
        cl <- safe_chr(result$compare$compare_logfile)
        if (length(cl)) compare_logfile <- cl[1]
      }

      group_summary <-
        if (nrow(phenoscanner_summary_df)) {
          phenoscanner_summary_df |>
            dplyr::filter(.data$exposure_group == grp_name)
        } else {
          tibble::tibble()
        }

      if (nrow(group_summary)) {
        stratum_rows <- group_summary |>
          dplyr::mutate(
            level = "stratum",
            sex_ancestry = paste(.data$sex, .data$ancestry, sep = ":"),
            skipped_queries = dplyr::if_else(
              dplyr::coalesce(.data$phenoscanner_skipped, FALSE),
              1L,
              0L
            ),
            compare_root = compare_root,
            compare_logfile = compare_logfile
          )

        total_row <- tibble::tibble(
          level = "group_total",
          exposure_group = grp_name,
          exposure = collapse_unique(group_summary$exposure),
          sex = NA_character_,
          ancestry = NA_character_,
          sex_ancestry = collapse_unique(stratum_rows$sex_ancestry),
          total_snps = sum(group_summary$total_snps, na.rm = TRUE),
          matches_removed = sum(group_summary$matches_removed, na.rm = TRUE),
          failures = sum(group_summary$failures, na.rm = TRUE),
          survivors = sum(group_summary$survivors, na.rm = TRUE),
          phenoscanner_skipped = any(dplyr::coalesce(group_summary$phenoscanner_skipped, FALSE)),
          skipped_queries = sum(stratum_rows$skipped_queries, na.rm = TRUE),
          skip_reason = collapse_unique(group_summary$skip_reason[dplyr::coalesce(group_summary$phenoscanner_skipped, FALSE)]),
          patterns = collapse_unique(group_summary$patterns),
          pre_filter_snapshot = collapse_unique(group_summary$pre_filter_snapshot),
          report_path = collapse_unique(group_summary$report_path),
          compare_root = compare_root,
          compare_logfile = compare_logfile
        )

        summary_cols <- c(
          "level","exposure_group","exposure","sex","ancestry","sex_ancestry",
          "total_snps","matches_removed","failures","survivors",
          "phenoscanner_skipped","skipped_queries","skip_reason","patterns",
          "pre_filter_snapshot","report_path","compare_root","compare_logfile"
        )

        stratum_rows <- stratum_rows |>
          dplyr::select(dplyr::all_of(summary_cols))

        total_row <- total_row |>
          dplyr::select(dplyr::all_of(summary_cols))

        summary_to_write <- dplyr::bind_rows(total_row, stratum_rows)

        summary_dir <- compare_root
        if (is.na(summary_dir) || !nzchar(summary_dir)) {
          summary_dir <- file.path(cache_dir, "output", slugify(exposure_label), "compare", "summary_fallback")
          dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
          warning(sprintf(
            "compare_root missing for exposure_group '%s'; saving PhenoScanner summary to %s.",
            grp_name, summary_dir
          ), call. = FALSE)
        } else {
          dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
        }

        summary_path <- file.path(
          summary_dir,
          sprintf("phenoscanner_summary_%s.csv", slugify(grp_name))
        )
        readr::write_csv(summary_to_write, summary_path)
        message(sprintf(
          "[PhenoScanner][%s] Summary saved to %s",
          grp_name,
          summary_path
        ))
        phenoscanner_summary_paths[[grp_name]] <- summary_path
      } else {
        message(sprintf(
          "[PhenoScanner][%s] No PhenoScanner summary rows available; skipping summary file.",
          grp_name
        ))
      }
    } else {
      failed_records <- append(
        failed_records,
        list(
          list(
            exposure_group = grp_name,
            exposure = exposure_label,
            error = result$error_message
          )
        )
      )
    }
  }

  failure_df <-
    if (length(failed_records)) {
      tibble::tibble(
        exposure_group = vapply(failed_records, `[[`, character(1), "exposure_group"),
        exposure = vapply(failed_records, `[[`, character(1), "exposure"),
        error = vapply(failed_records, `[[`, character(1), "error")
      )
    } else {
      tibble::tibble()
    }

  if (nrow(failure_df)) {
    message("\n=== ard_compare() failures ===")
    print(failure_df, n = Inf)
  }

  invisible(list(
    processed_groups = processed,
    failed_groups = failure_df,
    phenoscanner_summary = phenoscanner_summary_df,
    phenoscanner_summary_paths = if (length(phenoscanner_summary_paths)) {
      unlist(phenoscanner_summary_paths, use.names = TRUE)
    } else {
      character(0)
    }
  ))
}
