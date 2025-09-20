#' Run ARD compare across multiple sex/ancestry groups
#'
#' `ard_compare()` orchestrates repeated calls to [run_phenome_mr()] for a
#' collection of sex/ancestry groups and then assembles beta-scale summary
#' tables and forests that juxtapose the groups side-by-side. Existing beta
#' outputs are reused when available.
#'
#' @inheritParams run_phenome_mr
#' @param groups A list of named lists. Each inner list must contain
#'   `sex`, `ancestry`, and `exposure_snps` entries describing one group. The
#'   `sex` value must be one of `"both"`, `"male"`, or `"female"`. The
#'   `ancestry` code must match the inputs accepted by
#'   [run_phenome_mr()]. `exposure_snps` should contain the TwoSampleMR-style
#'   instrument data for that group.
#'
#' @return Invisibly returns `NULL` after writing the combined compare outputs
#'   to disk.
#' @export
ard_compare <- function(
    exposure,
    exposure_units,
    groups,
    sensitivity_enabled = c(
      "egger_intercept","egger_slope_agreement",
      "weighted_median","weighted_mode",
      "steiger_direction","leave_one_out",
      "ivw_Q","ivw_I2"
    ),
    sensitivity_pass_min = 6,
    Multiple_testing_correction = c("BH","bonferroni"),
    scatterplot = TRUE,
    snpforestplot = TRUE,
    leaveoneoutplot = TRUE,
    cache_dir = ardmr_cache_dir(),
    logfile = NULL,
    verbose = TRUE,
    confirm = "ask",
    force_refresh = FALSE
) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  if (missing(exposure)) stop("`exposure` is mandatory.")
  exposure <- as.character(exposure)
  if (!length(exposure)) stop("`exposure` is mandatory.")
  exposure <- exposure[1]
  if (is.na(exposure) || !nzchar(exposure)) stop("`exposure` is mandatory.")
  exposure <- trimws(exposure)
  if (!nzchar(exposure)) stop("`exposure` is mandatory.")
  if (missing(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- as.character(exposure_units)
  if (!length(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- exposure_units[1]
  if (is.na(exposure_units)) stop("`exposure_units` is mandatory.")
  exposure_units <- trimws(exposure_units)
  if (!nzchar(exposure_units)) stop("`exposure_units` is mandatory.")

  if (!is.list(groups) || !length(groups)) {
    stop("`groups` must be a non-empty list of group specifications.")
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  slug <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "NA"
    x <- gsub("[^[:alnum:]]+", "-", x)
    x <- gsub("-+", "-", x)
    x <- gsub("^-|-$", "", x)
    tolower(x)
  }

  shorttimestamp <- format(Sys.time(), "%m%d_%H%M")
  exposure_slug <- slug(exposure)
  compare_root <- file.path(cache_dir, "output", exposure_slug, "compare", shorttimestamp)
  beta_compare_dir <- file.path(compare_root, "beta_compare")
  if (!dir.exists(beta_compare_dir)) {
    dir.create(beta_compare_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ancestry_label <- function(sex, ancestry) {
    sx <- tolower(trimws(sex))
    anc <- toupper(trimws(ancestry))
    if (!sx %in% c("both","male","female")) {
      stop(sprintf("Unsupported sex value '%s'.", sex), call. = FALSE)
    }
    if (identical(sx, "male") || identical(sx, "female")) {
      if (!identical(anc, "EUR")) {
        stop(sprintf("Sex-specific group '%s' must use EUR ancestry (got '%s').", sex, ancestry), call. = FALSE)
      }
      return(if (identical(sx, "male")) "Male" else "Female")
    }
    if (!identical(anc, "EUR") && !identical(sx, "both")) {
      stop("Sex must be 'both' when ancestry is not EUR.", call. = FALSE)
    }
    rec <- c(
      "EUR" = "European",
      "AFR" = "African",
      "AMR" = "Admixed American",
      "CSA" = "Central/South Asian",
      "EAS" = "East Asian",
      "SAS" = "South Asian",
      "MID" = "Middle Eastern"
    )
    lbl <- rec[[anc]]
    if (is.null(lbl)) {
      lbl <- tools::toTitleCase(tolower(anc))
    }
    lbl
  }

  read_table <- function(path, required = TRUE) {
    if (!file.exists(path)) {
      if (!required) {
        return(tibble::tibble())
      }
      stop(sprintf("Expected table '%s' was not found.", path), call. = FALSE)
    }
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
    } else {
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }

  ensure_tibble <- function(df) {
    tibble::as_tibble(df)
  }

  groups_info <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    g <- groups[[i]]
    if (!is.list(g)) {
      stop(sprintf("Group %d must be a list with named elements.", i), call. = FALSE)
    }
    nms <- names(g)
    if (is.null(nms) || !length(nms)) {
      stop(sprintf("Group %d must provide named elements (sex, ancestry, exposure_snps).", i), call. = FALSE)
    }
    nms_lower <- tolower(trimws(nms))
    g_named <- stats::setNames(g, nms_lower)
    req <- c("sex","ancestry","exposure_snps")
    missing_req <- setdiff(req, names(g_named))
    if (length(missing_req)) {
      stop(sprintf("Group %d missing required fields: %s", i, paste(missing_req, collapse = ", ")), call. = FALSE)
    }
    sex <- g_named[["sex"]]
    ancestry <- g_named[["ancestry"]]
    exposure_snps <- g_named[["exposure_snps"]]
    if (!length(sex)) stop(sprintf("Group %d: sex is required.", i), call. = FALSE)
    if (!length(ancestry)) stop(sprintf("Group %d: ancestry is required.", i), call. = FALSE)
    sex <- tolower(as.character(sex)[1])
    ancestry <- toupper(as.character(ancestry)[1])
    if (!nzchar(sex) || !nzchar(ancestry)) {
      stop(sprintf("Group %d: sex and ancestry must be non-empty strings.", i), call. = FALSE)
    }
    if (!is.data.frame(exposure_snps)) {
      stop(sprintf("Group %d: exposure_snps must be a data frame.", i), call. = FALSE)
    }
    assert_exposure(exposure_snps)
    display <- ancestry_label(sex, ancestry)
    group_dir <- file.path(cache_dir, "output", exposure_slug, sex, ancestry)
    global_csv <- file.path(group_dir, "beta", "tables", "global.csv")
    groups_info[[i]] <- list(
      index = i,
      sex = sex,
      ancestry = ancestry,
      exposure_snps = exposure_snps,
      display = display,
      group_dir = group_dir,
      global_csv = global_csv,
      label = display,
      id = sprintf("group_%02d", i)
    )
  }

  sex_values <- vapply(groups_info, function(info) info$sex, character(1))
  if (any(sex_values == "both") && !all(sex_values == "both")) {
    stop(
      "comparisons of panukb and neale not allowed: check sex values to ensure one catalog only",
      call. = FALSE
    )
  }

  # ensure sex/both consistency for non-EUR ancestries
  for (info in groups_info) {
    if (!identical(info$ancestry, "EUR") && !identical(info$sex, "both")) {
      stop(sprintf(
        "Group %d uses ancestry '%s' but sex is '%s'. Only sex = 'both' allowed for non-EUR ancestries.",
        info$index, info$ancestry, info$sex
      ), call. = FALSE)
    }
  }

  if (isTRUE(verbose)) {
    logger::log_info("Assembling ARD compare outputs for {length(groups_info)} groups into {compare_root}")
  }

  # run or reuse each group
  for (info in groups_info) {
    if (file.exists(info$global_csv)) {
      if (isTRUE(verbose)) {
        logger::log_info("Reusing existing beta tables for sex={info$sex}, ancestry={info$ancestry} ({info$global_csv}).")
      }
      next
    }
    if (isTRUE(verbose)) {
      logger::log_info("Launching run_phenome_mr for sex={info$sex}, ancestry={info$ancestry}.")
    }
    run_phenome_mr(
      exposure = exposure,
      exposure_snps = info$exposure_snps,
      exposure_units = exposure_units,
      ancestry = info$ancestry,
      sex = info$sex,
      sensitivity_enabled = sensitivity_enabled,
      sensitivity_pass_min = sensitivity_pass_min,
      Multiple_testing_correction = Multiple_testing_correction,
      scatterplot = scatterplot,
      snpforestplot = snpforestplot,
      leaveoneoutplot = leaveoneoutplot,
      cache_dir = cache_dir,
      logfile = logfile,
      verbose = verbose,
      confirm = confirm,
      force_refresh = force_refresh
    )
  }

  load_group_tables <- function(info) {
    base <- file.path(info$group_dir, "beta", "tables")
    list(
      global = ensure_tibble(read_table(file.path(base, "global.csv"), required = TRUE)),
      cause_level_1 = list(
        all_diseases = ensure_tibble(read_table(file.path(base, "cause_level_1__all_diseases.csv"), required = FALSE)),
        age_related_diseases = ensure_tibble(read_table(file.path(base, "cause_level_1__age_related_diseases.csv"), required = FALSE))
      ),
      cause_level_2 = list(
        all_diseases = ensure_tibble(read_table(file.path(base, "cause_level_2__all_diseases.csv"), required = FALSE)),
        age_related_diseases = ensure_tibble(read_table(file.path(base, "cause_level_2__age_related_diseases.csv"), required = FALSE))
      ),
      cause_level_3 = list(
        all_diseases = ensure_tibble(read_table(file.path(base, "cause_level_3__all_diseases.csv"), required = FALSE)),
        age_related_diseases = ensure_tibble(read_table(file.path(base, "cause_level_3__age_related_diseases.csv"), required = FALSE))
      )
    )
  }

  group_tables <- lapply(groups_info, load_group_tables)

  # combine global tables
  combine_global <- function(group_tables, groups_info) {
    pieces <- list()
    for (i in seq_along(group_tables)) {
      tbl <- group_tables[[i]]$global
      info <- groups_info[[i]]
      if (!nrow(tbl)) next
      tbl$scope <- dplyr::case_when(
        tolower(tbl$group) == "age-related diseases" ~ "age_related_diseases",
        tolower(tbl$group) == "all diseases" ~ "all_diseases",
        TRUE ~ "other"
      )
      tbl$sex <- info$sex
      tbl$ancestry <- info$ancestry
      tbl$group_display <- info$display
      tbl$group_order <- info$index
      tbl$group_id <- info$id
      tbl$display_label <- paste0(tbl$group, ": ", info$display)
      tbl$scope_rank <- dplyr::case_when(
        tolower(tbl$group) == "age-related diseases" ~ 1L,
        tolower(tbl$group) == "all diseases" ~ 2L,
        TRUE ~ 3L
      )
      pieces[[length(pieces) + 1L]] <- tbl
    }
    if (!length(pieces)) {
      return(tibble::tibble())
    }
    out <- dplyr::bind_rows(pieces)
    out <- dplyr::arrange(out, .data$group_order, .data$scope_rank)
    out$scope_rank <- NULL
    out
  }

  global_combined <- combine_global(group_tables, groups_info)

  make_cause_dataset <- function(level, scope, group_tables, groups_info) {
    rows <- list()
    for (i in seq_along(group_tables)) {
      tbl <- group_tables[[i]][[level]][[scope]]
      info <- groups_info[[i]]
      if (!nrow(tbl)) next
      tbl$sex <- info$sex
      tbl$ancestry <- info$ancestry
      tbl$group_display <- info$display
      tbl$group_order <- info$index
      tbl$group_id <- info$id
      tbl$scope <- scope
      rows[[length(rows) + 1L]] <- tbl
    }
    if (!length(rows)) {
      return(list(table = tibble::tibble(), plot = tibble::tibble()))
    }
    combined <- dplyr::bind_rows(rows)
    combined <- dplyr::arrange(combined, .data$cause, .data$group_order)
    combined$display_label <- paste0(combined$cause, ": ", combined$group_display)

    causes <- unique(combined$cause)
    plot_rows <- list()
    for (idx in seq_along(causes)) {
      cz <- causes[[idx]]
      cz_rows <- combined[combined$cause == cz, , drop = FALSE]
      if (!nrow(cz_rows)) next
      cz_rows$axis_label <- cz_rows$display_label
      cz_rows$axis_id <- paste0("cause", idx, "_", seq_len(nrow(cz_rows)))
      cz_rows$is_separator <- FALSE
      plot_rows[[length(plot_rows) + 1L]] <- cz_rows
      if (idx < length(causes)) {
        sep <- tibble::tibble(
          level = cz_rows$level[1],
          cause = cz,
          exposure = if ("exposure" %in% names(cz_rows)) cz_rows$exposure[1] else NA_character_,
          n = NA_integer_,
          ivw_mean_beta = NA_real_,
          se_ivw_mean = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          sex = NA_character_,
          ancestry = NA_character_,
          group_display = NA_character_,
          group_order = NA_integer_,
          group_id = NA_character_,
          scope = scope,
          display_label = NA_character_,
          axis_label = "",
          axis_id = paste0("sep_", idx),
          is_separator = TRUE
        )
        missing_cols <- setdiff(names(cz_rows), names(sep))
        if (length(missing_cols)) {
          for (col in missing_cols) {
            sep[[col]] <- NA
          }
        }
        sep <- sep[names(cz_rows)]
        plot_rows[[length(plot_rows) + 1L]] <- sep
      }
    }
    plot_df <- dplyr::bind_rows(plot_rows)
    list(table = combined, plot = plot_df)
  }

  cause_levels <- c("cause_level_1","cause_level_2","cause_level_3")
  scopes <- c("all_diseases","age_related_diseases")
  cause_datasets <- list()
  for (lvl in cause_levels) {
    cause_datasets[[lvl]] <- list()
    for (sc in scopes) {
      cause_datasets[[lvl]][[sc]] <- make_cause_dataset(lvl, sc, group_tables, groups_info)
    }
  }

  # write combined tables
  tables_dir <- file.path(beta_compare_dir, "tables")
  if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(df, path) {
    if (!nrow(df)) return(invisible(NULL))
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_csv(df, path)
    } else {
      utils::write.csv(df, path, row.names = FALSE)
    }
  }

  if (nrow(global_combined)) {
    write_csv(global_combined, file.path(tables_dir, "global.csv"))
  }
  for (lvl in cause_levels) {
    for (sc in scopes) {
      tbl <- cause_datasets[[lvl]][[sc]]$table
      if (!nrow(tbl)) next
      fn <- sprintf("%s__%s.csv", lvl, sc)
      write_csv(tbl, file.path(tables_dir, fn))
    }
  }

  save_plot <- function(plot, dir_path, file_stem, n_rows, width = 6.0) {
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    height <- max(2.4, 1.2 + 0.28 * n_rows)
    file_path <- file.path(dir_path, paste0(file_stem, ".png"))
    ggplot2::ggsave(filename = file_path, plot = plot, width = width, height = height, dpi = 300)
    .ardmr_write_plot_data(plot, dir_path = dir_path, base_name = file_stem)
    invisible(NULL)
  }

  # build and save plots
  if (nrow(global_combined)) {
    out_dir <- file.path(beta_compare_dir, "global")
    plot_global <- plot_beta_mean_global_compare(global_combined, exposure_units = exposure_units)
    save_plot(plot_global, out_dir, "mean_effect", nrow(global_combined))
    plot_global_wrap <- plot_beta_mean_global_compare_wrap(global_combined, exposure_units = exposure_units)
    save_plot(plot_global_wrap, out_dir, "mean_effect_wrap", nrow(global_combined))
  }

  for (lvl in cause_levels) {
    for (sc in scopes) {
      plot_df <- cause_datasets[[lvl]][[sc]]$plot
      if (!nrow(plot_df)) next
      out_dir <- file.path(beta_compare_dir, lvl, sc)
      n_rows <- nrow(plot_df)
      plot_obj <- plot_beta_mean_cause_compare(plot_df, exposure_units = exposure_units)
      save_plot(plot_obj, out_dir, "mean_effect", n_rows)
      plot_obj_wrap <- plot_beta_mean_cause_compare_wrap(plot_df, exposure_units = exposure_units)
      save_plot(plot_obj_wrap, out_dir, "mean_effect_wrap", n_rows)
    }
  }

  invisible(NULL)
}
