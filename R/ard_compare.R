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
    confirm = "yes",
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

  catalog <- if (groups[[1]]$sex == "both") "panukb" else "neale"
  # or by name: catalog <- if (groups[["eur_both"]]$sex == "both") "panukb" else "neale"

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
  if (!dir.exists(compare_root)) {
    dir.create(compare_root, recursive = TRUE, showWarnings = FALSE)
  }
  compare_logfile <- file.path(
    compare_root,
    sprintf("ard_compare_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  compare_log_level <- if (isTRUE(verbose)) "INFO" else "WARN"
  setup_logging(compare_logfile, level = compare_log_level)
  restore_compare_logging <- function() {
    if (nzchar(compare_logfile)) {
      logger::log_appender(logger::appender_tee(compare_logfile))
    } else {
      logger::log_appender(logger::appender_console())
    }
    logger::log_layout(logger::layout_glue_colors)
    logger::log_threshold(compare_log_level)
  }
  beta_compare_dir <- file.path(compare_root, "beta_compare")
  if (!dir.exists(beta_compare_dir)) {
    dir.create(beta_compare_dir, recursive = TRUE, showWarnings = FALSE)
  }

  enrich_compare_dir <- file.path(compare_root, "enrichment_compare")

  if (!dir.exists(enrich_compare_dir)) {
    dir.create(enrich_compare_dir, recursive = TRUE, showWarnings = FALSE)
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

  ancestry_display <- function(ancestry) {
    ancestry_label("both", ancestry)
  }

  sex_display <- function(sex) {
    sx <- tolower(trimws(sex))
    switch(
      sx,
      both = "both sexes",
      male = "male",
      female = "female",
      stop(sprintf("Unsupported sex value '%s'.", sex), call. = FALSE)
    )
  }

  group_axis_label <- function(sex, ancestry) {
    paste0(ancestry_display(ancestry), ", ", sex_display(sex))
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
      results_rds = file.path(group_dir, "results.rds"),
      label = display,
      id = sprintf("group_%02d", i),
      axis_label = group_axis_label(sex, ancestry)
    )
  }

  sex_values <- vapply(groups_info, function(info) info$sex, character(1))
  if (any(sex_values == "both") && !all(sex_values == "both")) {
    stop(
      "comparisons of panukb and neale not allowed: check sex values to ensure one catalog only",
      call. = FALSE
    )
  }

  compare_effect_scale <- if (length(sex_values) && all(sex_values == "both")) {
    "odds_ratio"
  } else {
    "absolute_risk"
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
    restore_compare_logging()
    logger::log_info("Assembling ARD compare outputs for {length(groups_info)} groups into {compare_root}")
  }

  # run or reuse each group
  for (info in groups_info) {
    if (file.exists(info$global_csv)) {
      if (isTRUE(verbose)) {
        restore_compare_logging()
        logger::log_info("Reusing existing beta tables for sex={info$sex}, ancestry={info$ancestry} ({info$global_csv}).")
      }
      next
    }
    if (isTRUE(verbose)) {
      restore_compare_logging()
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

  load_group_enrichment <- function(info) {
    rds_path <- file.path(info$group_dir, "results.rds")
    if (!file.exists(rds_path)) {
      return(tibble::tibble())
    }

    res <- tryCatch(readRDS(rds_path), error = function(e) NULL)
    if (!is.list(res) || is.null(res$enrich)) {
      return(tibble::tibble())
    }
    enrich <- res$enrich
    if (!is.list(enrich) || !is.data.frame(enrich$by_cause_tbl)) {
      return(tibble::tibble())
    }

    tibble::as_tibble(enrich$by_cause_tbl)
  }

  group_tables <- lapply(groups_info, load_group_tables)
  group_enrichment <- lapply(groups_info, load_group_enrichment)

  load_group_enrichment <- function(info) {
    results_path <- info$results_rds
    if (!file.exists(results_path)) {
      if (isTRUE(verbose)) {
        restore_compare_logging()
        logger::log_warn(
          "Global enrichment results not found for sex={info$sex}, ancestry={info$ancestry} ({results_path})."
        )
      }
      return(tibble::tibble())
    }
    enr <- tryCatch(readRDS(results_path), error = function(e) e)
    if (inherits(enr, "error")) {
      if (isTRUE(verbose)) {
        restore_compare_logging()
        logger::log_warn(
          "Failed to read enrichment results for sex={info$sex}, ancestry={info$ancestry}: {enr$message}"
        )
      }
      return(tibble::tibble())
    }
    if (!is.list(enr) || !"enrich" %in% names(enr)) {
      if (isTRUE(verbose)) {
        restore_compare_logging()
        logger::log_warn(
          "Enrichment component missing in results for sex={info$sex}, ancestry={info$ancestry}."
        )
      }
      return(tibble::tibble())
    }
    global_tbl <- tryCatch(enr$enrich$global_tbl, error = function(e) NULL)
    if (!is.data.frame(global_tbl) || !nrow(global_tbl)) {
      return(tibble::tibble())
    }
    df <- tibble::as_tibble(global_tbl)
    df$group_id <- info$id
    df$group_order <- info$index
    df$group_axis_label <- info$axis_label
    df$group_display <- info$display
    df$sex <- info$sex
    df$ancestry <- info$ancestry
    df
  }

  enrichment_global_combined <- dplyr::bind_rows(lapply(groups_info, load_group_enrichment))
  if (nrow(enrichment_global_combined)) {
    enrichment_global_combined <- dplyr::arrange(enrichment_global_combined, .data$group_order)
  }

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

  format_q_label <- function(q) {
    out <- rep("", length(q))
    ok <- is.finite(q)
    if (any(ok)) {
      hi <- ok & q >= 0.001
      lo <- ok & !hi
      out[hi] <- sprintf("%.3f", q[hi])
      out[lo] <- sprintf("%.2e", q[lo])
    }
    out
  }

  build_heatmap_dataset <- function(
      level,
      scope,
      cause_datasets,
      group_enrichment,
      groups_info,
      exposure_label,
      compare_mode = "cause_vs_rest_all",
      p_threshold = 0.05
  ) {
    cause_tbl <- cause_datasets[[level]][[scope]]$table
    if (!is.data.frame(cause_tbl) || !nrow(cause_tbl)) {
      return(tibble::tibble())
    }

    causes <- unique(cause_tbl$cause)
    causes <- causes[!is.na(causes)]
    causes <- trimws(as.character(causes))
    causes <- causes[nzchar(causes)]
    if (!length(causes)) {
      return(tibble::tibble())
    }

    cause_order <- stats::setNames(seq_along(causes), causes)

    rows <- vector("list", length(groups_info))
    for (i in seq_along(groups_info)) {
      info <- groups_info[[i]]
      enrich_tbl <- group_enrichment[[i]]
      if (is.data.frame(enrich_tbl) && nrow(enrich_tbl) &&
          all(c("level","compare_mode","cause") %in% names(enrich_tbl))) {
        enrich_tbl <- enrich_tbl[
          enrich_tbl$level == level &
            enrich_tbl$compare_mode == compare_mode,
          ,
          drop = FALSE
        ]
      } else {
        enrich_tbl <- tibble::tibble()
      }

      if (nrow(enrich_tbl)) {
        enrich_tbl$cause <- trimws(as.character(enrich_tbl$cause))
      }

      match_idx <- if (nrow(enrich_tbl)) match(causes, enrich_tbl$cause) else rep(NA_integer_, length(causes))
      pull_vals <- function(col, default = NA_real_) {
        if (!nrow(enrich_tbl) || is.null(enrich_tbl[[col]])) {
          rep(default, length(causes))
        } else {
          vals <- rep(default, length(causes))
          ok <- !is.na(match_idx)
          if (any(ok)) {
            vals[ok] <- enrich_tbl[[col]][match_idx[ok]]
          }
          vals
        }
      }

      p_vals <- pull_vals("p_signed")
      ses_vals <- pull_vals("SES_signed")
      q_vals <- rep(NA_real_, length(p_vals))
      finite_idx <- which(is.finite(p_vals))
      if (length(finite_idx)) {
        q_vals[finite_idx] <- stats::p.adjust(p_vals[finite_idx], method = "BH")
      }

      rows[[i]] <- tibble::tibble(
        level = level,
        scope = scope,
        compare_mode = compare_mode,
        cause = causes,
        cause_order = unname(cause_order[causes]),
        group_id = info$id,
        group_label = info$display,
        group_order = info$index,
        p_signed = p_vals,
        SES_signed = ses_vals,
        q_value = q_vals,
        significant = is.finite(q_vals) & (q_vals < p_threshold),
        direction = dplyr::case_when(
          is.finite(ses_vals) & ses_vals < 0 ~ "protective",
          is.finite(ses_vals) & ses_vals > 0 ~ "risk",
          TRUE ~ "neutral"
        ),
        label = format_q_label(q_vals),
        exposure = exposure_label
      )
    }

    out <- dplyr::bind_rows(rows)
    if (!nrow(out)) {
      return(tibble::tibble())
    }
    out$direction[out$significant %in% FALSE] <- "neutral"
    out
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

  heatmap_datasets <- list()
  for (lvl in cause_levels) {
    heatmap_datasets[[lvl]] <- list()
    for (sc in scopes) {
      heatmap_datasets[[lvl]][[sc]] <- build_heatmap_dataset(
        level = lvl,
        scope = sc,
        cause_datasets = cause_datasets,
        group_enrichment = group_enrichment,
        groups_info = groups_info,
        exposure_label = exposure
      )
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

  save_plot <- function(plot, dir_path, file_stem, n_rows, width = 6.0, height = NULL) {
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    file_path <- file.path(dir_path, paste0(file_stem, ".png"))
    is_yfloat <- is.character(file_stem) && grepl("_wrap_yfloat$", file_stem)
    if (is.null(height)) {
      height <- if (is_yfloat) {
        1.2 + 0.28 * max(1, n_rows)
      } else {
        max(2.4, 1.2 + 0.28 * n_rows)
      }
    }
    ggsave_args <- list(
      filename = file_path,
      plot = plot,
      width = width,
      height = height,
      dpi = 300
    )
    do.call(ggplot2::ggsave, ggsave_args)
    .ardmr_write_plot_data(plot, dir_path = dir_path, base_name = file_stem)
    invisible(NULL)
  }

  # build and save plots
  if (nrow(global_combined)) {
    out_dir <- file.path(beta_compare_dir, "global")
    plot_global <- plot_beta_mean_global_compare(
      global_combined,
      exposure_units = exposure_units,
      effect_scale = compare_effect_scale
    )
    # save_plot(plot_global, out_dir, "mean_effect", nrow(global_combined))
    plot_global_wrap <- plot_beta_mean_global_compare_wrap(
      global_combined,
      exposure_units = exposure_units,
      effect_scale = compare_effect_scale
    )
    # save_plot(plot_global_wrap, out_dir, "mean_effect_wrap", nrow(global_combined))
    plot_global_wrap_yfloat <- plot_beta_mean_global_compare_wrap_yfloat(
      global_combined,
      exposure_units = exposure_units,
      effect_scale = compare_effect_scale
    )
    # save_plot(plot_global_wrap_yfloat, out_dir, "mean_effect_wrap_yfloat", nrow(global_combined))
  }

  for (lvl in cause_levels) {
    for (sc in scopes) {
      plot_df <- cause_datasets[[lvl]][[sc]]$plot
      if (!nrow(plot_df)) next
      out_dir <- file.path(beta_compare_dir, lvl, sc)
      n_rows <- nrow(plot_df)
      plot_obj <- plot_beta_mean_cause_compare(
        plot_df,
        exposure_units = exposure_units,
        effect_scale = compare_effect_scale
      )
      # save_plot(plot_obj, out_dir, "mean_effect", n_rows)
      plot_obj_wrap <- plot_beta_mean_cause_compare_wrap(
        plot_df,
        exposure_units = exposure_units,
        effect_scale = compare_effect_scale
      )
      # save_plot(plot_obj_wrap, out_dir, "mean_effect_wrap", n_rows)
      plot_obj_wrap_yfloat <- plot_beta_mean_cause_compare_wrap_yfloat(
        plot_df,
        exposure_units = exposure_units,
        effect_scale = compare_effect_scale
      )
      # save_plot(plot_obj_wrap_yfloat, out_dir, "mean_effect_wrap_yfloat", n_rows)

      heat_df <- heatmap_datasets[[lvl]][[sc]]
      if (is.data.frame(heat_df) && nrow(heat_df)) {
        heat_plot <- plot_enrichment_group_heatmap(
          heat_df,
          level = lvl,
          scope = sc,
          exposure = exposure
        )
        n_causes <- length(unique(heat_df$cause))
        heat_out_dir <- file.path(enrich_compare_dir, lvl, sc)
        save_plot(heat_plot, heat_out_dir, "group_heatmap", max(1, n_causes), width = 7.2)
      }
    }
  }

  enrichment_compare_dir <- file.path(compare_root, "enrichment_compare")
  if (nrow(enrichment_global_combined)) {
    n_groups <- enrichment_global_combined$group_axis_label
    n_groups <- n_groups[!is.na(n_groups)]
    n_groups <- length(unique(n_groups))
    n_groups <- if (is.finite(n_groups) && n_groups > 0) n_groups else 1L
    plot_enrichment_compare <- plot_enrichment_signed_violin_global_compare(
      enrichment_global_combined,
      alpha = 0.05
    )
    save_plot(
      plot_enrichment_compare,
      file.path(enrichment_compare_dir, "global"),
      "violin_forest",
      n_groups,
      width = 7.2,
      height = 1.8 + 0.4 * max(1, n_groups)
    )
  } else if (isTRUE(verbose)) {
    restore_compare_logging()
    logger::log_warn("No enrichment data available for compare violin forest plot.")
  }

  invisible(NULL)
}
