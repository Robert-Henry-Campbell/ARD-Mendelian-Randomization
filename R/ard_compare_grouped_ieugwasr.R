#' Build instruments from IEU GWAS, validate units, assemble groups, and run ard_compare()
#'
#' @param csv_path Path to a CSV with columns: ieu_id, Exposure, exposure_group, sex, ancestry
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
  # See ieugwasr authentication guide: https://mrcieu.github.io/ieugwasr/articles/guide_authentication.html
  Sys.setenv(OPENGWAS_JWT = jwt)

  # ---- helpers ----
  is_rsid <- function(x) !is.na(x) & grepl("^rs\\d+$", x)
  palindrome_flag <- function(a1, a2) {
    a1 <- toupper(a1); a2 <- toupper(a2)
    (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C")
  }
  f_stat <- function(beta, se) (beta/se)^2
  safe_chr <- function(x) { x <- as.character(x); x[!nzchar(x)|is.na(x)] <- NA_character_; x }

  fetch_units <- function(ieu_id) {
    info <- try(ieugwasr::gwasinfo(ieu_id), silent = TRUE)
    if (!inherits(info,"try-error") && "unit" %in% names(info)) {
      u <- safe_chr(info$unit[1]); if (!is.na(u)) return(u)
    }
    tr <- try(TwoSampleMR::get_trait(ieu_id), silent = TRUE)
    if (!inherits(tr,"try-error") && "unit" %in% names(tr)) {
      u <- safe_chr(tr$unit[1]); if (!is.na(u)) return(u)
    }
    NA_character_
  }

  annotate_chrpos <- function(dat) {
    snps <- unique(dat$SNP)
    anns <- ieugwasr::variants(snps)
    anns <- dplyr::transmute(anns, SNP = rsid, Chr = as.character(chr), Pos = as.numeric(pos))
    dplyr::left_join(dat, anns, by = "SNP")
  }

  build_exposure_snps <- function(raw, exposure_label, ieu_id) {
    if (any(!is_rsid(raw$SNP))) {
      bad <- unique(raw$SNP[!is_rsid(raw$SNP)])
      stop(sprintf("Non-rsID variant IDs for %s: e.g. %s", ieu_id, paste(utils::head(bad,5), collapse=", ")))
    }
    out <- raw %>%
      mutate(
        id.exposure            = ieu_id,
        Exposure               = exposure_label,
        effect_allele.exposure = effect_allele,
        other_allele.exposure  = other_allele,
        eaf.exposure           = eaf,
        beta.exposure          = beta,
        se.exposure            = se,
        pval.exposure          = pval,
        palindromic            = palindrome_flag(effect_allele, other_allele),
        mr_keep                = TRUE,
        F                      = f_stat(beta.exposure, se.exposure)
      ) %>%
      select(
        id.exposure, Exposure, SNP, Chr, Pos,
        effect_allele.exposure, other_allele.exposure, eaf.exposure,
        beta.exposure, se.exposure, palindromic, pval.exposure,
        mr_keep, F
      ) %>%
      arrange(SNP)

    need <- c("id.exposure","Exposure","SNP","Chr","Pos",
              "effect_allele.exposure","other_allele.exposure","eaf.exposure",
              "beta.exposure","se.exposure","palindromic","pval.exposure",
              "mr_keep","F")
    stopifnot(all(need %in% names(out)))
    out
  }

  get_instruments <- function(ieu_id, p1 = p_threshold, r2 = r2, kb = kb) {
    ins <- TwoSampleMR::extract_instruments(outcomes = ieu_id, p1 = p1, clump = TRUE, r2 = r2, kb = kb)
    if (!nrow(ins)) return(ins[0, ])
    ins <- annotate_chrpos(ins)
    ins <- ins %>%
      mutate(F_tmp = f_stat(beta, se)) %>%
      filter(is.finite(F_tmp), F_tmp >= f_threshold) %>%
      select(-F_tmp)
    ins
  }

  # ---- read CSV ----
  expos <- readr::read_csv(csv_path, show_col_types = FALSE)
  req_cols <- c("ieu_id","Exposure","exposure_group","sex","ancestry")
  miss_cols <- setdiff(req_cols, names(expos))
  if (length(miss_cols)) stop("CSV missing required columns: ", paste(miss_cols, collapse = ", "))

  # ---- units (preflight) ----
  unique_ids <- unique(expos$ieu_id)
  units_map <- setNames(rep(NA_character_, length(unique_ids)), unique_ids)
  for (id in unique_ids) units_map[[id]] <- fetch_units(id)

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
    group_split(.row_id, keep = TRUE) %>%
    purrr::map(function(rr) {
      r <- rr[[1]]
      ieu_id   <- r$ieu_id
      exp_name <- r$Exposure

      ins <- get_instruments(ieu_id)
      if (!nrow(ins)) {
        warning(sprintf("No instruments retained for %s; skipping.", ieu_id))
        return(NULL)
      }

      snps_df <- build_exposure_snps(ins, exposure_label = exp_name, ieu_id = ieu_id)

      list(
        exposure_group = r$exposure_group,
        sex            = match.arg(as.character(r$sex), c("both","male","female")),
        ancestry       = as.character(r$ancestry),
        exposure       = exp_name,
        exposure_units = units_map[[ieu_id]],
        exposure_snps  = snps_df
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
      warning(sprintf("Multiple 'Exposure' labels in group '%s'. Using the first: %s", grp_name, exp_label[1]))
    }
    exposure_label <- exp_label[1]

    # check units agree across entries
    u_vals <- unique(vapply(entries, `[[`, "", "exposure_units"))
    if (length(u_vals) > 1L) {
      stop(sprintf("Units differ within exposure_group '%s': %s", grp_name, paste(u_vals, collapse = " | ")))
    }
    exposure_units <- u_vals[1]

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
      Multiple_testing_correction = "BH",
      scatterplot                 = TRUE,
      snpforestplot               = TRUE,
      leaveoneoutplot             = TRUE,
      cache_dir                   = cache_dir,
      logfile                     = NULL,
      verbose                     = TRUE,
      confirm                     = "yes",
      force_refresh               = FALSE
    )

    processed <- c(processed, grp_name)
  }

  invisible(list(processed_groups = processed))
}
