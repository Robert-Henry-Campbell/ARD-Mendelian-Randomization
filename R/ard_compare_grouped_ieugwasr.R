#' Build instruments from IEU GWAS, validate units, assemble groups, and run ard_compare()
#'
#' @param csv_path Path to a CSV with columns: ieu_id, exposure, exposure_group,
#'   sex, ancestry, multiple_testing_correction, confirm, and optionally
#'   exposure_units
#' @param cache_dir Output/cache directory used by ard_compare/run_phenome_mr
#' @param jwt IEU JWT string; defaults to IEU_JWT/OPENGWAS_JWT env vars
#' @param prompt_for_units If TRUE, interactively prompt for missing units; otherwise stop on missing units
#' @param p_threshold Genome-wide significance threshold for instrument extraction (default 5e-8)
#' @param r2 LD clumping r^2 threshold (default 0.001)
#' @param kb Clumping window in kb (default 10000)
#' @param f_threshold Per-variant F-stat minimum (default 10)
#'
#' @return (invisible) list with names of processed exposure_groups
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
    r2             = 0.001,
    kb             = 10000,
    f_threshold    = 10
) {
  # ---- packages ----
  pkgs <- c("TwoSampleMR","ieugwasr","dplyr","readr","purrr","stringr","tibble")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Please install required packages: ", paste(miss, collapse = ", "))

  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)

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
      return(dat |> mutate(Chr = as.character(.data$chr.exposure),
                           Pos = as.numeric(.data$pos.exposure)))
    }
    if (all(c("chr","pos") %in% names(dat))) {
      return(dat |> mutate(Chr = as.character(.data$chr),
                           Pos = as.numeric(.data$pos)))
    }

    if (!"SNP" %in% names(dat)) stop("annotate_chrpos(): 'SNP' column is missing.")
    snps <- unique(dat$SNP)

    ann <- try(ieugwasr::variants_rsid(snps), silent = TRUE)
    if (!inherits(ann, "try-error") && all(c("name","chr","pos") %in% names(ann))) {
      ann <- ann |>
        transmute(SNP = .data$name,
                  Chr = as.character(.data$chr),
                  Pos = as.numeric(.data$pos))
      return(left_join(dat, ann, by = "SNP"))
    }

    look_one_study <- function(study_id, rs) {
      out <- ieugwasr::associations(variants = rs, id = study_id)
      if (!nrow(out)) return(out[0, ])
      out |> select(rsid, chr, pos) |> distinct()
    }

    study_id <- if ("id.exposure" %in% names(dat)) dat$id.exposure[1] else NA_character_
    if (is.na(study_id)) {
      warning("No id.exposure in dat; cannot fetch chr/pos. Returning NA Chr/Pos.")
      return(dat |> mutate(Chr = NA_character_, Pos = as.numeric(NA)))
    }

    chunks <- split(snps, ceiling(seq_along(snps)/500))
    ann_list <- lapply(chunks, function(rs) try(look_one_study(study_id, rs), silent = TRUE))
    ann_list <- ann_list[!vapply(ann_list, inherits, logical(1), "try-error")]
    if (length(ann_list)) {
      ann2 <- bind_rows(ann_list) |>
        transmute(SNP = .data$rsid,
                  Chr = as.character(.data$chr),
                  Pos = as.numeric(.data$pos))
      return(left_join(dat, ann2, by = "SNP"))
    }

    dat |> mutate(Chr = NA_character_, Pos = as.numeric(NA))
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

    out <- raw %>%
      mutate(
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
      ) %>%
      select(
        id.exposure, exposure, Exposure, SNP, Chr, Pos,
        effect_allele.exposure, other_allele.exposure, eaf.exposure,
        beta.exposure, se.exposure, palindromic, pval.exposure,
        mr_keep, F
      ) %>%
      arrange(SNP)

    need <- c("id.exposure","exposure","SNP","Chr","Pos",
              "effect_allele.exposure","other_allele.exposure","eaf.exposure",
              "beta.exposure","se.exposure","palindromic","pval.exposure",
              "mr_keep","F")
    stopifnot(all(need %in% names(out)))
    out
  }

  # instrument extraction + F filtering (robust to column names)
  get_instruments <- function(ieu_id) {
    # helper to compute F with whatever column names we get back
    pick_colname <- function(df, candidates) {
      nm <- candidates[candidates %in% names(df)][1]
      if (is.na(nm)) NA_character_ else nm
    }

    # 1) Try the normal call (your settings)
    ins <- try(
      TwoSampleMR::extract_instruments(
        outcomes = ieu_id, p1 = p_threshold, clump = TRUE, r2 = r2, kb = kb
      ),
      silent = TRUE
    )

    if (inherits(ins, "try-error")) {
      message(sprintf("[extract_instruments] error for %s. Retrying with clump=FALSE…", ieu_id))
      # 2) Retry without clumping (sometimes the clumping step is where it breaks)
      ins <- try(
        TwoSampleMR::extract_instruments(
          outcomes = ieu_id, p1 = p_threshold, clump = FALSE
        ),
        silent = TRUE
      )
      if (inherits(ins, "try-error")) {
        message(sprintf("[extract_instruments] still error for %s. Probing dataset availability…", ieu_id))
        # 3) Probe that the study ID is valid & accessible
        info_ok <- !inherits(try(ieugwasr::gwasinfo(ieu_id), silent = TRUE), "try-error")
        # Try to pull a single variant to see if the dataset is reachable at all
        probe <- try(
          TwoSampleMR::extract_outcome_data(variants = "rs56116432", outcomes = ieu_id),
          silent = TRUE
        )
        if (!info_ok) {
          warning(sprintf("Study '%s' not found via gwasinfo() -> likely bad ID. Skipping.", ieu_id))
        } else if (inherits(probe, "try-error")) {
          warning(sprintf("Study '%s' exists but OpenGWAS access failed (JWT/permission/network). Skipping.", ieu_id))
        } else {
          warning(sprintf(
            "Study '%s' is accessible, but instrument endpoint errored (API quirk). Treating as failure & skipping.",
            ieu_id
          ))
        }
        return(tibble::tibble())  # graceful skip
      }
    }

    # If we get here, we have a data.frame (maybe empty = 'no hits')
    if (!nrow(ins)) {
      message(sprintf("No instruments for %s at p<=%g (pre-clump). Returning 0 rows.", ieu_id, p_threshold))
      return(ins[0, ])
    }

    # F filter (handle column name variants)
    beta_nm <- pick_colname(ins, c("beta","beta.exposure"))
    se_nm   <- pick_colname(ins, c("se","se.exposure"))
    if (!nzchar(beta_nm) || !nzchar(se_nm)) {
      warning(sprintf("Missing beta/se columns in instruments for %s; skipping.", ieu_id))
      return(tibble::tibble())
    }
    ins$F_tmp <- (ins[[beta_nm]] / ins[[se_nm]])^2
    ins <- dplyr::filter(ins, is.finite(.data$F_tmp), .data$F_tmp >= f_threshold) |>
      dplyr::select(-.data$F_tmp)
    if (!nrow(ins)) {
      message(sprintf("All instruments for %s dropped by F-stat >= %g. Returning 0 rows.", ieu_id, f_threshold))
      return(ins[0, ])
    }

    # Chr/Pos annotation (your robust helper)
    ins <- annotate_chrpos(ins)
    ins
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

  expos <- expos %>%
    mutate(
      exposure = safe_chr(exposure),
      exposure_group = safe_chr(exposure_group),
      sex = safe_chr(sex),
      ancestry = safe_chr(ancestry),
      exposure_units = safe_chr(exposure_units),
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

  # ---- build per-row instrument tables ----
  message("Fetching & preparing instruments…")
  expos_prepared <- expos %>%
    mutate(.row_id = dplyr::row_number()) %>%
    group_split(.row_id, .keep = TRUE) %>%   # modern dplyr arg
    purrr::map(function(rr) {
      r <- dplyr::slice(rr, 1)
      ieu_id   <- r$ieu_id
      exp_name <- r$exposure

      ins <- get_instruments(ieu_id)
      if (!nrow(ins)) {
        warning(sprintf("No instruments retained for %s; skipping.", ieu_id))
        return(NULL)
      }

      snps_df <- build_exposure_snps(ins, exposure_label = exp_name, ieu_id = ieu_id)

      mtc_val    <- normalize_mtc(r$multiple_testing_correction)[1]
      confirm_val <- normalize_confirm(r$confirm)[1]

      list(
        exposure_group = r$exposure_group,
        sex            = match.arg(as.character(r$sex), c("both","male","female")),
        ancestry       = as.character(r$ancestry),
        exposure       = exp_name,
        exposure_units = units_map[[ieu_id]],
        exposure_snps  = snps_df,
        multiple_testing_correction = mtc_val,
        confirm        = confirm_val
      )
    }) %>%
    purrr::compact()

  if (!length(expos_prepared)) stop("No valid rows after instrument preparation.")

  # ---- assemble groups and run ard_compare() per exposure_group ----
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  by_group <- split(expos_prepared, vapply(expos_prepared, `[[`, "", "exposure_group"))

  processed <- character(0)

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
        exposure_snps = e$exposure_snps
      )
    })
    names(groups_list) <- paste0(
      vapply(entries, `[[`, "", "sex"), "_", vapply(entries, `[[`, "", "ancestry")
    )

    message(sprintf("\n=== ard_compare(): group '%s' ===", grp_name))
    message(sprintf("Exposure: %s | Units: %s", exposure_label, exposure_units))

    ard_compare(
      exposure                    = exposure_label,
      exposure_units              = exposure_units,
      groups                      = groups_list,
      sensitivity_enabled         = c(
        "egger_intercept","egger_slope_agreement",
        "weighted_median","weighted_mode",
        "steiger_direction","leave_one_out",
        "ivw_Q","ivw_I2"
      ),
      sensitivity_pass_min        = 6,
      Multiple_testing_correction = multiple_testing_correction,
      scatterplot                 = TRUE,
      snpforestplot               = TRUE,
      leaveoneoutplot             = TRUE,
      cache_dir                   = cache_dir,
      logfile                     = NULL,
      verbose                     = TRUE,
      confirm                     = confirm_val,
      force_refresh               = FALSE
    )

    processed <- c(processed, grp_name)
  }

  invisible(list(processed_groups = processed))
}
