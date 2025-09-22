# R/enrichment_business_logic.R

#' Threshold-free directional enrichment via label permutations (business logic)
#'
#' Computes threshold-free, directional enrichment using **label permutations**
#' on a signed evidence score \code{s = sign(beta_ivw) * (-log10 p_ivw)}.
#'
#' Provides:
#' \itemize{
#'   \item \strong{Global}: ARD vs non-ARD across the full phenome
#'         (returns a single two-sided signed test).
#'   \item \strong{Cause-level (L1/L2/L3)}: comparison modes per cause,
#'         notably \code{"cause_vs_rest_all"} (cause vs rest with global permutations).
#' }
#'
#' For each test we report the signed statistic \code{stat_signed}, its permutation
#' mean/SD, the standardized enrichment score \code{SES_signed}, the two-sided
#' permutation p-value \code{p_signed}, the BH-adjusted q-value \code{q_signed}, and
#' (optionally) a retained vector of permutation draws for visualization.
#'
#' @section Expected `results_df` columns:
#' Required: \code{results_outcome} (chr), \code{results_beta_ivw} (dbl),
#' \code{results_p_ivw} (dbl in (0,1]), \code{ARD_selected} (logical),
#' \code{cause_level_1/2/3} (chr; NA allowed).
#'
#' Strongly recommended: \code{results_qc_pass} (lgl), \code{results_nsnp_after} (int),
#' \code{results_se_ivw} (dbl; enables weights), \code{pheno_sex} (chr in
#' {"male","female","both_sex"}).
#'
#' @section Inference:
#' Non-parametric; no normality assumed. Null = exchangeability of labels under the
#' permutation scope (global / within-cause / among-ARD-only). BH is applied per
#' tested family (per tail and for the signed test).
#'
#' @name enrichment_business_logic
#' @keywords enrichment MR permutation ARD
#' @importFrom dplyr filter mutate transmute group_by ungroup arrange bind_rows summarise across n distinct
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom stats sd p.adjust
NULL

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(rlang)
})

# ---------- internal helpers ----------

# @keywords internal
.as_int01 <- function(x) if (is.logical(x)) as.integer(x) else as.integer(x != 0)

# @keywords internal
.sanitize_p <- function(p, floor = 1e-300) { p2 <- p; p2[!is.finite(p2) | p2 <= 0] <- floor; p2 }

#' Signed IVW evidence: s = sign(beta) * (-log10 p)
#' @keywords internal
signed_ivw_score <- function(beta, p) {
  s <- -log10(.sanitize_p(p))
  s[!is.finite(s)] <- NA_real_
  sign(beta) * s
}

#' Choose weights ("none" or "inv_se2")
#' @keywords internal
.choose_weights <- function(se, scheme = c("inv_se2","none")) {
  scheme <- match.arg(scheme)
  if (scheme == "inv_se2") {
    w <- 1 / (se^2); w[!is.finite(w)] <- NA_real_
  } else w <- rep(1, length(se))
  w
}

#' Tail contributions (protective keeps negative; risk keeps positive)
#' @keywords internal
.tail_contrib <- function(signed_score, weights, tail = c("protective","risk")) {
  tail <- match.arg(tail)
  contrib <- if (tail == "protective") weights * pmax(0, -signed_score) else weights * pmax(0, signed_score)
  contrib[is.na(contrib)] <- 0
  contrib
}

#' Core permutation engine for non-negative tail sums (one-sided p)
#' @keywords internal
.perm_test_sum <- function(contrib, y, exact_max_combn = 1e5, mc_B = 10000, seed = NULL) {
  y <- .as_int01(y); n <- length(y); a <- sum(y)
  if (n == 0L || a == 0L || a == n) {
    return(list(stat_obs=0, p=1, mean_perm=0, sd_perm=0, ses=0, n_perm=0L, exact=FALSE))
  }
  stat_obs <- sum(contrib[y == 1L])
  if (sum(contrib) <= 0) {
    return(list(stat_obs=stat_obs, p=1, mean_perm=0, sd_perm=0, ses=0, n_perm=0L, exact=FALSE))
  }
  total_combn <- suppressWarnings(choose(n, a))
  do_exact <- is.finite(total_combn) && (total_combn <= exact_max_combn) && (n <= 30)

  if (do_exact) {
    idx <- utils::combn(n, a)
    m <- ncol(idx); aint <- a
    chunk <- 5000L; acc <- numeric(m); j <- 1L
    while (j <= m) {
      k <- min(m, j + chunk - 1L)
      sel <- idx[, j:k, drop = FALSE]
      acc[j:k] <- colSums(matrix(contrib[sel], nrow = aint))
      j <- k + 1L
    }
    stat_perm <- acc
  } else {
    if (!is.null(seed)) set.seed(seed)
    aint <- as.integer(a)
    idx <- replicate(mc_B, sample.int(n, size = aint, replace = FALSE))
    stat_perm <- colSums(matrix(contrib[idx], nrow = aint))
  }

  n_perm <- length(stat_perm)
  p <- (1 + sum(stat_perm >= stat_obs)) / (1 + n_perm)
  mu <- mean(stat_perm); sdv <- stats::sd(stat_perm)
  ses <- if (sdv > 0) (stat_obs - mu) / sdv else 0
  list(stat_obs=stat_obs, p=p, mean_perm=mu, sd_perm=sdv, ses=ses, n_perm=n_perm, exact=do_exact)
}

#' Core permutation engine for signed sums (two-sided p)
#' @keywords internal
.perm_test_signed <- function(
    contrib_signed,
    y,
    exact_max_combn = 1e5,
    mc_B = 10000,
    seed = NULL,
    retain_perm = FALSE,
    retain_perm_max = 10000
) {
  # contrib_signed can be negative or positive (w * s)
  y <- .as_int01(y); n <- length(y); a <- sum(y)
  if (n == 0L || a == 0L || a == n) {
    return(list(stat_obs=0, p=1, mean_perm=0, sd_perm=0, ses=0, n_perm=0L, exact=FALSE))
  }
  stat_obs <- sum(contrib_signed[y == 1L])

  total_combn <- suppressWarnings(choose(n, a))
  do_exact <- is.finite(total_combn) && (total_combn <= exact_max_combn) && (n <= 30)

  if (do_exact) {
    idx <- utils::combn(n, a)
    m <- ncol(idx); aint <- a
    chunk <- 5000L; acc <- numeric(m); j <- 1L
    while (j <= m) {
      k <- min(m, j + chunk - 1L)
      sel <- idx[, j:k, drop = FALSE]
      acc[j:k] <- colSums(matrix(contrib_signed[sel], nrow = aint))
      j <- k + 1L
    }
    stat_perm <- acc
  } else {
    if (!is.null(seed)) set.seed(seed)
    aint <- as.integer(a)
    idx <- replicate(mc_B, sample.int(n, size = aint, replace = FALSE))
    stat_perm <- colSums(matrix(contrib_signed[idx], nrow = aint))
  }

  n_perm <- length(stat_perm)
  mu <- mean(stat_perm); sdv <- stats::sd(stat_perm)
  # two-sided permutation p using absolute deviation from mean
  p <- (1 + sum(abs(stat_perm - mu) >= abs(stat_obs - mu))) / (1 + n_perm)
  ses <- if (sdv > 0) (stat_obs - mu) / sdv else 0

  perm_keep <- NULL
  if (isTRUE(retain_perm)) {
    if (!is.null(stat_perm)) {
      keep_n <- min(n_perm, if (is.null(retain_perm_max)) n_perm else as.integer(retain_perm_max))
      if (keep_n < 0L || !is.finite(keep_n)) keep_n <- n_perm
      if (keep_n > 0L) {
        perm_keep <- stat_perm[seq_len(keep_n)]
      } else {
        perm_keep <- numeric(0)
      }
    } else {
      perm_keep <- numeric(0)
    }
  }

  list(
    stat_obs = stat_obs,
    p = p,
    mean_perm = mu,
    sd_perm = sdv,
    ses = ses,
    n_perm = n_perm,
    exact = do_exact,
    stat_perm = perm_keep
  )
}

# Try to capture a default exposure label (from exposure_snps$id.exposure if present in caller)
#' @keywords internal
.default_exposure <- function() {
  val <- tryCatch({
    es <- get("exposure_snps", envir = parent.frame(), inherits = TRUE)
    es$id.exposure
  }, error = function(e) NA_character_)
  if (length(val) == 1) as.character(val) else if (length(val) > 1) paste0(unique(as.character(val)), collapse = "; ") else NA_character_
}

# ---------- public API (call these from run_phenome_mr) ----------

#' Run all enrichment analyses for one exposure
#'
#' Orchestrates:
#' \itemize{
#'   \item Global ARD vs non-ARD (tails + single signed).
#'   \item Cause-level L1/L2/L3 with the three comparison modes (tails + signed).
#' }
#'
#' @param results_df Tidy MR results for the exposure.
#' @param exposure Character exposure label. Defaults to \code{exposure_snps$id.exposure} if discoverable.
#' @param levels Character vector among \code{c("cause_level_1","cause_level_2","cause_level_3")}.
#' @param modes Character vector among the three compare modes (see details).
#' @param use_qc_pass Logical; if TRUE and present, filter \code{results_qc_pass==TRUE}.
#' @param min_nsnp Integer; if present, filter \code{results_nsnp_after >= min_nsnp}.
#' @param weight_scheme Weighting scheme; defaults to \code{"inv_se2"}
#'   (1/SE^2, requires \code{results_se_ivw}) with \code{"none"}
#'   available as a fallback.
#' @param exact_max_combn Maximum \eqn{choose(n,a)} for exact enumeration; otherwise Monte Carlo.
#' @param mc_B Monte Carlo permutations if not exact.
#' @param seed RNG seed.
#' @param retain_permutations Logical; if TRUE, keep up to \code{retain_perm_max}
#'   permutation draws per test (for plotting).
#' @param retain_perm_max Maximum number of permutation draws to retain when
#'   \code{retain_permutations = TRUE}.
#' @return A list with: \code{$global_tbl} and \code{$by_cause_tbl} (all levels × modes concatenated).
#' @export
run_enrichment <- function(
    results_df,
    exposure = .default_exposure(),
    levels = c("cause_level_1","cause_level_2","cause_level_3"),
    modes  = c("cause_vs_rest_all"),
    use_qc_pass = TRUE,
    min_nsnp = 2,
    weight_scheme = c("inv_se2","none"),
    exact_max_combn = 1e5,
    mc_B = 10000,
    seed = NULL,
    retain_permutations = FALSE,
    retain_perm_max = 10000
) {
  weight_scheme <- match.arg(weight_scheme)
  levels <- match.arg(levels, c("cause_level_1","cause_level_2","cause_level_3"), several.ok = TRUE)
  modes <- match.arg(modes,
                     c("ARD_vs_nonARD_within_cause","cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere"),
                     several.ok = TRUE)

  # --- global ---
  global_tbl <- enrichment_global_directional(
    results_df = results_df,
    exposure   = exposure,
    use_qc_pass = use_qc_pass,
    min_nsnp    = min_nsnp,
    weight_scheme = weight_scheme,
    exact_max_combn = exact_max_combn,
    mc_B = mc_B,
    seed = seed,
    retain_permutations = retain_permutations,
    retain_perm_max = retain_perm_max
  )

  # --- cause-level across requested levels × modes ---
  by_cause <- list()
  for (lv in levels) {
    for (md in modes) {
      tmp <- enrichment_by_cause_directional(
        results_df = results_df,
        level = lv,
        compare_mode = md,
        exposure = exposure,
        use_qc_pass = use_qc_pass,
        min_nsnp = min_nsnp,
        weight_scheme = weight_scheme,
        exact_max_combn = exact_max_combn,
        mc_B = mc_B,
        seed = seed,
        retain_permutations = retain_permutations,
        retain_perm_max = retain_perm_max
      )
      by_cause[[paste(lv, md, sep = "::")]] <- tmp
    }
  }
  by_cause_tbl <- dplyr::bind_rows(by_cause)

  list(global_tbl = global_tbl, by_cause_tbl = by_cause_tbl)
}

#' Global ARD vs non-ARD signed enrichment (two-sided)
#'
#' @inheritParams run_enrichment
#' @return One-row tibble with counts, signed statistic, permutation summary,
#'   and (optionally) retained permutation draws.
#' @export
enrichment_global_directional <- function(
    results_df,
    exposure = .default_exposure(),
    use_qc_pass = TRUE,
    min_nsnp = 2,
    weight_scheme = c("inv_se2","none"),
    exact_max_combn = 1e5,
    mc_B = 10000,
    seed = NULL,
    retain_permutations = FALSE,
    retain_perm_max = 10000
) {
  weight_scheme <- match.arg(weight_scheme)

  df <- results_df
  if (use_qc_pass && "results_qc_pass" %in% names(df)) df <- dplyr::filter(df, .data$results_qc_pass)
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df)) df <- dplyr::filter(df, .data$results_nsnp_after >= !!min_nsnp)

  stopifnot(all(c("results_beta_ivw","results_p_ivw","ARD_selected") %in% names(df)))
  s <- signed_ivw_score(df$results_beta_ivw, df$results_p_ivw)
  w <- if ("results_se_ivw" %in% names(df)) .choose_weights(df$results_se_ivw, weight_scheme) else .choose_weights(rep(NA_real_, nrow(df)), "none")
  y <- .as_int01(df$ARD_selected)

  # signed (two-sided)
  c_signed <- w * s
  ts <- .perm_test_signed(
    c_signed,
    y,
    exact_max_combn = exact_max_combn,
    mc_B = mc_B,
    seed = if (is.null(seed)) NULL else seed,
    retain_perm = retain_permutations,
    retain_perm_max = retain_perm_max
  )

  sex_summary <- if ("pheno_sex" %in% names(df)) paste0(sort(unique(na.omit(as.character(df$pheno_sex)))), collapse = ";") else NA_character_

  perm_vals <- if (isTRUE(retain_permutations)) ts$stat_perm else NULL
  if (is.null(perm_vals)) perm_vals <- numeric(0)

  out <- tibble::tibble(
    scope   = "GLOBAL",
    exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
    pheno_sex = sex_summary,
    n_total = length(y),
    n_ard   = sum(y),
    n_non   = length(y) - sum(y),
    # signed
    stat_signed = ts$stat_obs,
    mean_signed = ts$mean_perm,
    sd_signed = ts$sd_perm,
    p_signed = ts$p,
    SES_signed = ts$ses,
    n_perm_signed = ts$n_perm,
    exact_signed = ts$exact,
    perm_scope = "global",
    perm_stat = list(perm_vals)
  )

  out |>
    dplyr::mutate(q_signed = stats::p.adjust(p_signed, method = "BH"))
}

#' Cause-level signed enrichment (two-sided) for one level × one mode
#' 
#' @param level One of \code{"cause_level_1"}, \code{"cause_level_2"}, \code{"cause_level_3"}.
#' @param compare_mode One of:
#'   \code{"ARD_vs_nonARD_within_cause"}, \code{"cause_vs_rest_all"},
#'   \code{"ARD_in_cause_vs_ARD_elsewhere"}.
#' @inheritParams run_enrichment
#' @return Tibble: one row per cause with signed statistic, permutation summary,
#'   and (optionally) retained permutation draws.
#' @export
enrichment_by_cause_directional <- function(
    results_df,
    level = c("cause_level_2","cause_level_1","cause_level_3"),
    compare_mode = c("ARD_vs_nonARD_within_cause","cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere"),
    exposure = .default_exposure(),
    use_qc_pass = TRUE,
    min_nsnp = 2,
    weight_scheme = c("inv_se2","none"),
    exact_max_combn = 1e5,
    mc_B = 10000,
    seed = NULL,
    retain_permutations = FALSE,
    retain_perm_max = 10000
) {
  level <- match.arg(level); compare_mode <- match.arg(compare_mode); weight_scheme <- match.arg(weight_scheme)

  df0 <- results_df
  if (use_qc_pass && "results_qc_pass" %in% names(df0)) df0 <- dplyr::filter(df0, .data$results_qc_pass)
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df0)) df0 <- dplyr::filter(df0, .data$results_nsnp_after >= !!min_nsnp)
  df0 <- dplyr::filter(df0, !is.na(.data[[level]]))

  stopifnot(all(c("results_beta_ivw","results_p_ivw") %in% names(df0)))
  df0 <- dplyr::mutate(
    df0,
    signed_score = signed_ivw_score(results_beta_ivw, results_p_ivw),
    w = if ("results_se_ivw" %in% names(df0)) .choose_weights(results_se_ivw, weight_scheme) else 1
  )

  run_one <- function(df_all, cause_label, mode, idx) {
    if (mode == "ARD_vs_nonARD_within_cause") {
      g <- dplyr::filter(df_all, .data[[level]] == cause_label)
      if (!"ARD_selected" %in% names(g)) return(NULL)
      y <- .as_int01(g$ARD_selected)
      scope_counts <- c(n_total = nrow(g), n_pos = sum(y), n_neg = nrow(g) - sum(y))
      perm_scope <- "within_cause"
      s <- g$signed_score; w <- g$w
    } else if (mode == "cause_vs_rest_all") {
      g <- df_all
      y <- .as_int01(g[[level]] == cause_label)
      scope_counts <- c(n_total = nrow(g), n_pos = sum(y), n_neg = nrow(g) - sum(y))
      perm_scope <- "global"
      s <- g$signed_score; w <- g$w
    } else { # ARD_in_cause_vs_ARD_elsewhere
      if (!"ARD_selected" %in% names(df_all)) return(NULL)
      g <- dplyr::filter(df_all, .data$ARD_selected == TRUE)
      y <- .as_int01(g[[level]] == cause_label)
      scope_counts <- c(n_total = nrow(g), n_pos = sum(y), n_neg = nrow(g) - sum(y))
      perm_scope <- "among_ARD_only"
      s <- g$signed_score; w <- g$w
    }

    if (length(y) == 0L || sum(y) == 0L || sum(y) == length(y)) {
      return(tibble::tibble(
        level = level, cause = cause_label, compare_mode = mode,
        exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
        pheno_sex = if ("pheno_sex" %in% names(g)) paste0(sort(unique(na.omit(as.character(g$pheno_sex)))), collapse = ";") else NA_character_,
        n_total = scope_counts["n_total"], n_pos = scope_counts["n_pos"], n_neg = scope_counts["n_neg"],
        stat_signed = NA_real_,
        mean_signed = NA_real_,
        sd_signed = NA_real_,
        p_signed = NA_real_,
        SES_signed = NA_real_,
        n_perm_signed = NA_integer_,
        exact_signed = NA,
        perm_scope = perm_scope,
        perm_stat = list(numeric(0))
      ))
    }

    # signed (two-sided)
    c_signed <- w * s
    seed_use <- if (is.null(seed)) NULL else seed + as.integer(idx)
    ts <- .perm_test_signed(
      c_signed,
      y,
      exact_max_combn = exact_max_combn,
      mc_B = mc_B,
      seed = seed_use,
      retain_perm = retain_permutations,
      retain_perm_max = retain_perm_max
    )

    perm_vals <- if (isTRUE(retain_permutations)) ts$stat_perm else NULL
    if (is.null(perm_vals)) perm_vals <- numeric(0)

    tibble::tibble(
      level = level, cause = cause_label, compare_mode = mode,
      exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
      pheno_sex = if ("pheno_sex" %in% names(g)) paste0(sort(unique(na.omit(as.character(g$pheno_sex)))), collapse = ";") else NA_character_,
      n_total = scope_counts["n_total"], n_pos = scope_counts["n_pos"], n_neg = scope_counts["n_neg"],
      # signed
      stat_signed = ts$stat_obs,
      mean_signed = ts$mean_perm,
      sd_signed = ts$sd_perm,
      p_signed = ts$p,
      SES_signed = ts$ses,
      n_perm_signed = ts$n_perm,
      exact_signed = ts$exact,
      perm_scope = perm_scope,
      perm_stat = list(perm_vals)
    )
  }

  causes <- sort(unique(df0[[level]]))
  out <- dplyr::bind_rows(lapply(seq_along(causes), function(i) run_one(df0, causes[[i]], compare_mode, i)))

  # BH adjustment across the causes within this level/mode
  out |>
    dplyr::mutate(q_signed = stats::p.adjust(p_signed, method = "BH"))
}

# ---------- beta-contrast (inside vs outside) ----------
#' Inside-vs-outside contrast of mean MR betas by cause
#'
#' For each cause (at a chosen level), compute inverse-variance weighted means
#' of beta_IVW inside the cause and outside the cause, then the difference
#' Δβ = β̄_in - β̄_out, its SE, z, p, and BH/Bonferroni q across causes.
#'
#' Assumes all betas are on a common scale (e.g., log-OR per 1 SD exposure).
#'
#' @param results_df Tidy MR results with columns: results_beta_ivw, results_se_ivw,
#'        results_qc_pass (optional), cause_level_* (chosen via `level`).
#' @param level "cause_level_1","cause_level_2","cause_level_3".
#' @param use_qc_pass If TRUE, keep results_qc_pass==TRUE when present.
#' @param min_nsnp Optional integer filter on results_nsnp_after.
#' @param exposure Optional label for reporting.
#' @param Multiple_testing_correction "BH" or "bonferroni" for q-values.
#' @param alpha Significance level (default 0.05).
#' @return Tibble with one row per cause: n_in, n_out, beta_in, beta_out,
#'         delta_beta, se_delta, z, p, q, ci_low, ci_high.
#' @export
beta_contrast_by_cause <- function(
    results_df,
    level = c("cause_level_1","cause_level_2","cause_level_3"),
    use_qc_pass = TRUE,
    min_nsnp = 2,
    exposure = NULL,
    Multiple_testing_correction = c("BH","bonferroni"),
    alpha = 0.05
) {
  level <- match.arg(level)
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)

  df <- tibble::as_tibble(results_df)
  stopifnot(all(c("results_beta_ivw","results_se_ivw", level) %in% names(df)))

  empty_template <- tibble::tibble(
    level = character(), cause = character(), exposure = character(),
    n_in = integer(), n_out = integer(),
    beta_in = double(), beta_out = double(),
    delta_beta = double(), se_delta = double(),
    z = double(), p = double(),
    ci_low = double(), ci_high = double(),
    q = double(),
    sig = logical()
  )

  if (use_qc_pass && "results_qc_pass" %in% names(df)) {
    df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  }
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df)) {
    df <- dplyr::filter(df, .data$results_nsnp_after >= !!min_nsnp)
  }
  df <- dplyr::filter(df, !is.na(.data[[level]]),
                      is.finite(.data$results_beta_ivw),
                      is.finite(.data$results_se_ivw),
                      .data$results_se_ivw > 0)

  if (!nrow(df)) {
    return(empty_template)
  }

  causes <- sort(unique(df[[level]]))

  one <- function(cz) {
    g <- df
    in_i  <- which(g[[level]] == cz)
    out_i <- which(g[[level]] != cz)
    if (!length(in_i) || !length(out_i)) {
      return(NULL)
    }
    b_in  <- g$results_beta_ivw[in_i]
    s_in  <- g$results_se_ivw[in_i]
    w_in  <- 1/(s_in^2)

    b_out <- g$results_beta_ivw[out_i]
    s_out <- g$results_se_ivw[out_i]
    w_out <- 1/(s_out^2)

    beta_in  <- sum(w_in * b_in)  / sum(w_in)
    beta_out <- sum(w_out * b_out)/ sum(w_out)

    var_in   <- 1 / sum(w_in)
    var_out  <- 1 / sum(w_out)
    se_delta <- sqrt(var_in + var_out)

    delta <- beta_in - beta_out
    z     <- if (se_delta > 0) delta / se_delta else NA_real_
    p     <- if (is.finite(z)) 2*stats::pnorm(-abs(z)) else NA_real_

    ci_low  <- delta - 1.96*se_delta
    ci_high <- delta + 1.96*se_delta

    tibble::tibble(
      level = level, cause = cz,
      exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
      n_in = length(in_i), n_out = length(out_i),
      beta_in = beta_in, beta_out = beta_out,
      delta_beta = delta, se_delta = se_delta,
      z = z, p = p, ci_low = ci_low, ci_high = ci_high
    )
  }

  out <- dplyr::bind_rows(lapply(causes, one))
  if (!nrow(out)) return(empty_template)

  out$q <- if (Multiple_testing_correction == "BH") {
    stats::p.adjust(out$p, method = "BH")
  } else {
    # Bonferroni family-wise across causes
    pmin(1, out$p * nrow(out))
  }

  dplyr::mutate(
    out,
    sig = if (Multiple_testing_correction == "BH") q < alpha else p < alpha / nrow(out)
  )
}

# One-row: global ARD vs non-ARD contrast of mean β (IVW)
beta_contrast_global_ARD <- function(
    results_df,
    use_qc_pass = TRUE, min_nsnp = 2,
    exposure = NULL,
    Multiple_testing_correction = c("BH","bonferroni"),
    alpha = 0.05
) {
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  df <- tibble::as_tibble(results_df)
  stopifnot(all(c("results_beta_ivw","results_se_ivw","ARD_selected") %in% names(df)))

  if (use_qc_pass && "results_qc_pass" %in% names(df)) df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df)) df <- dplyr::filter(df, .data$results_nsnp_after >= !!min_nsnp)
  df <- dplyr::filter(df, is.finite(.data$results_beta_ivw), is.finite(.data$results_se_ivw), .data$results_se_ivw > 0)

  w <- 1/(df$results_se_ivw^2)
  ard <- df$ARD_selected %in% TRUE
  if (!any(ard) || all(ard)) {
    return(tibble::tibble(level="GLOBAL", cause="ARD vs non-ARD", exposure=exposure,
                          n_in=sum(ard), n_out=sum(!ard),
                          beta_in=NA_real_, beta_out=NA_real_,
                          delta_beta=NA_real_, se_delta=NA_real_,
                          ci_low=NA_real_, ci_high=NA_real_,
                          z=NA_real_, p=NA_real_, q=NA_real_, sig=NA))
  }

  beta_in  <- sum(w[ard]*df$results_beta_ivw[ard])  / sum(w[ard])
  beta_out <- sum(w[!ard]*df$results_beta_ivw[!ard])/ sum(w[!ard])
  var_in   <- 1/sum(w[ard]); var_out <- 1/sum(w[!ard])
  se_delta <- sqrt(var_in + var_out)
  delta    <- beta_in - beta_out
  ci_low   <- delta - 1.96 * se_delta
  ci_high  <- delta + 1.96 * se_delta
  z        <- delta/se_delta
  p        <- 2*stats::pnorm(-abs(z))
  q        <- if (Multiple_testing_correction=="BH") p else p # single test => same
  tibble::tibble(level="GLOBAL", cause="ARD vs non-ARD", exposure=exposure,
                 n_in=sum(ard), n_out=sum(!ard),
                 beta_in=beta_in, beta_out=beta_out,
                 delta_beta=delta, se_delta=se_delta,
                 ci_low=ci_low, ci_high=ci_high,
                 z=z, p=p, q=q, sig=(p < alpha))
}

# Generalized β-contrast by cause & compare_mode (mirrors enrichment modes)
beta_contrast_by_cause_mode <- function(
    results_df,
    level = c("cause_level_1","cause_level_2","cause_level_3"),
    compare_mode = c("cause_vs_rest_all","ARD_vs_nonARD_within_cause","ARD_in_cause_vs_ARD_elsewhere"),
    use_qc_pass = TRUE, min_nsnp = 2,
    exposure = NULL,
    Multiple_testing_correction = c("BH","bonferroni"),
    alpha = 0.05
) {
  level <- match.arg(level); compare_mode <- match.arg(compare_mode)
  Multiple_testing_correction <- match.arg(Multiple_testing_correction)
  df <- tibble::as_tibble(results_df)

  stopifnot(all(c("results_beta_ivw","results_se_ivw", level) %in% names(df)))

  empty_template <- tibble::tibble(
    level = character(), cause = character(), compare_mode = character(),
    exposure = character(),
    n_in = integer(), n_out = integer(),
    beta_in = double(), beta_out = double(),
    delta_beta = double(), se_delta = double(),
    z = double(), p = double(),
    q = double(), sig = logical(),
    ci_low = double(), ci_high = double()
  )
  if (use_qc_pass && "results_qc_pass" %in% names(df)) df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df)) df <- dplyr::filter(df, .data$results_nsnp_after >= !!min_nsnp)
  df <- dplyr::filter(df, !is.na(.data[[level]]), is.finite(.data$results_beta_ivw),
                      is.finite(.data$results_se_ivw), .data$results_se_ivw > 0)

  if (!nrow(df)) return(empty_template)

  causes <- sort(unique(df[[level]]))
  ivw_delta <- function(b, se, g) {
    w <- 1/(se^2)
    beta_in  <- sum(w[g]*b[g]) / sum(w[g])
    beta_out <- sum(w[!g]*b[!g]) / sum(w[!g])
    var_in <- 1/sum(w[g]); var_out <- 1/sum(w[!g])
    se_delta <- sqrt(var_in + var_out)
    delta <- beta_in - beta_out
    z <- delta/se_delta; p <- 2*stats::pnorm(-abs(z))
    list(beta_in=beta_in, beta_out=beta_out, delta=delta, se_delta=se_delta, z=z, p=p)
  }

  rows <- lapply(causes, function(cz) {
    gdf <- df
    if (compare_mode == "ARD_vs_nonARD_within_cause") {
      if (!"ARD_selected" %in% names(gdf)) return(NULL)
      gdf <- dplyr::filter(gdf, .data[[level]] == cz)
      grp <- gdf$ARD_selected %in% TRUE
    } else if (compare_mode == "ARD_in_cause_vs_ARD_elsewhere") {
      if (!"ARD_selected" %in% names(gdf)) return(NULL)
      gdf <- dplyr::filter(gdf, .data$ARD_selected %in% TRUE)
      grp <- gdf[[level]] == cz
    } else { # cause_vs_rest_all
      grp <- gdf[[level]] == cz
    }
    if (!length(grp) || all(grp) || !any(grp)) return(NULL)
    res <- ivw_delta(gdf$results_beta_ivw, gdf$results_se_ivw, grp)
    tibble::tibble(
      level = level, cause = cz, compare_mode = compare_mode,
      exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
      n_in = sum(grp), n_out = sum(!grp),
      beta_in = res$beta_in, beta_out = res$beta_out,
      delta_beta = res$delta, se_delta = res$se_delta,
      z = res$z, p = res$p
    )
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(empty_template)

  out$q <- if (Multiple_testing_correction == "BH") {
    stats::p.adjust(out$p, method = "BH")
  } else {
    pmin(1, out$p * nrow(out))
  }
  out$sig <- if (Multiple_testing_correction == "BH") out$q < alpha else out$p < alpha / nrow(out)
  out$ci_low  <- out$delta_beta - 1.96*out$se_delta
  out$ci_high <- out$delta_beta + 1.96*out$se_delta
  out
}


# ---------- beta means (IVW) ----------
#' Inverse-variance weighted mean MR beta by cause
#'
#' Computes the IVW mean of MR beta estimates within each cause at a
#' specified cause level. The standard error is derived from the inverse of
#' the summed weights (1/SE^2) and 95% confidence intervals use ±1.96×SE.
#'
#' @param results_df Tidy MR results with columns results_beta_ivw,
#'   results_se_ivw, the chosen cause level column, and optionally
#'   ARD_selected.
#' @param level One of "cause_level_1", "cause_level_2", "cause_level_3".
#' @param use_qc_pass If TRUE, keep rows with results_qc_pass == TRUE when the
#'   column exists.
#' @param min_nsnp Optional integer filter on results_nsnp_after.
#' @param exposure Optional character exposure label recorded in the output.
#' @param ard_only If TRUE, restrict to rows with ARD_selected == TRUE.
#' @param drop_empty If TRUE, omit causes that have zero rows after filtering
#'   (useful for ARD-only summaries).
#' @return Tibble with columns level, cause, exposure, n, ivw_mean_beta,
#'   se_ivw_mean, ci_low, ci_high.
beta_mean_by_cause <- function(
    results_df,
    level = c("cause_level_1","cause_level_2","cause_level_3"),
    use_qc_pass = TRUE,
    min_nsnp = 2,
    exposure = NULL,
    ard_only = FALSE,
    drop_empty = TRUE
) {
  level <- match.arg(level)
  df <- tibble::as_tibble(results_df)

  required_cols <- c("results_beta_ivw", "results_se_ivw", level)
  stopifnot(all(required_cols %in% names(df)))

  if (use_qc_pass && "results_qc_pass" %in% names(df)) {
    df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  }
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df)) {
    df <- dplyr::filter(df, .data$results_nsnp_after >= !!min_nsnp)
  }

  df <- dplyr::filter(
    df,
    !is.na(.data[[level]]),
    is.finite(.data$results_beta_ivw),
    is.finite(.data$results_se_ivw),
    .data$results_se_ivw > 0
  )

  if (ard_only) {
    if (!"ARD_selected" %in% names(df)) {
      return(tibble::tibble(
        level = character(), cause = character(), exposure = character(),
        n = integer(), ivw_mean_beta = double(), se_ivw_mean = double(),
        ci_low = double(), ci_high = double()
      ))
    }
    df <- dplyr::filter(df, .data$ARD_selected %in% TRUE)
  }

  if (!nrow(df)) {
    return(tibble::tibble(
      level = character(), cause = character(), exposure = character(),
      n = integer(), ivw_mean_beta = double(), se_ivw_mean = double(),
      ci_low = double(), ci_high = double()
    ))
  }

  causes <- sort(unique(df[[level]]))

  calc_ivw <- function(beta, se) {
    w <- 1 / (se^2)
    mean_beta <- sum(w * beta) / sum(w)
    se_mean <- sqrt(1 / sum(w))
    list(
      mean = mean_beta,
      se = se_mean,
      ci_low = mean_beta - 1.96 * se_mean,
      ci_high = mean_beta + 1.96 * se_mean
    )
  }

  rows <- lapply(causes, function(cz) {
    g <- df[df[[level]] == cz, , drop = FALSE]
    if (!nrow(g)) {
      if (drop_empty) return(NULL)
      return(tibble::tibble(
        level = level, cause = cz,
        exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
        n = 0L,
        ivw_mean_beta = NA_real_,
        se_ivw_mean = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      ))
    }
    stats <- calc_ivw(g$results_beta_ivw, g$results_se_ivw)
    tibble::tibble(
      level = level,
      cause = cz,
      exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
      n = nrow(g),
      ivw_mean_beta = stats$mean,
      se_ivw_mean = stats$se,
      ci_low = stats$ci_low,
      ci_high = stats$ci_high
    )
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(tibble::tibble(
      level = character(),
      cause = character(),
      exposure = character(),
      n = integer(),
      ivw_mean_beta = double(),
      se_ivw_mean = double(),
      ci_low = double(),
      ci_high = double()
    ))
  }

  out <- dplyr::bind_rows(rows)
  out
}

#' Global IVW mean MR beta for all diseases and ARDs
#'
#' @inheritParams beta_mean_by_cause
#' @return Tibble with rows describing "All Diseases" and "Age-Related Diseases"
#'   when available.
beta_mean_global <- function(
    results_df,
    use_qc_pass = TRUE,
    min_nsnp = 2,
    exposure = NULL
) {
  df <- tibble::as_tibble(results_df)
  stopifnot(all(c("results_beta_ivw", "results_se_ivw") %in% names(df)))

  if (use_qc_pass && "results_qc_pass" %in% names(df)) {
    df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  }
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df)) {
    df <- dplyr::filter(df, .data$results_nsnp_after >= !!min_nsnp)
  }
  df <- dplyr::filter(
    df,
    is.finite(.data$results_beta_ivw),
    is.finite(.data$results_se_ivw),
    .data$results_se_ivw > 0
  )

  calc_ivw <- function(beta, se) {
    w <- 1 / (se^2)
    mean_beta <- sum(w * beta) / sum(w)
    se_mean <- sqrt(1 / sum(w))
    list(
      mean = mean_beta,
      se = se_mean,
      ci_low = mean_beta - 1.96 * se_mean,
      ci_high = mean_beta + 1.96 * se_mean
    )
  }

  rows <- list()

  if (nrow(df)) {
    stats_all <- calc_ivw(df$results_beta_ivw, df$results_se_ivw)
    rows[["all"]] <- tibble::tibble(
      group = "All Diseases",
      exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
      n = nrow(df),
      ivw_mean_beta = stats_all$mean,
      se_ivw_mean = stats_all$se,
      ci_low = stats_all$ci_low,
      ci_high = stats_all$ci_high
    )
  }

  if ("ARD_selected" %in% names(df)) {
    ard_df <- df[df$ARD_selected %in% TRUE, , drop = FALSE]
    if (nrow(ard_df)) {
      stats_ard <- calc_ivw(ard_df$results_beta_ivw, ard_df$results_se_ivw)
      rows[["ard"]] <- tibble::tibble(
        group = "Age-Related Diseases",
        exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
        n = nrow(ard_df),
        ivw_mean_beta = stats_ard$mean,
        se_ivw_mean = stats_ard$se,
        ci_low = stats_ard$ci_low,
        ci_high = stats_ard$ci_high
      )
    } else {
      rows[["ard"]] <- tibble::tibble(
        group = "Age-Related Diseases",
        exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
        n = 0L,
        ivw_mean_beta = NA_real_,
        se_ivw_mean = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    }
  }

  if (!length(rows)) {
    return(tibble::tibble(
      group = character(),
      exposure = character(),
      n = integer(),
      ivw_mean_beta = double(),
      se_ivw_mean = double(),
      ci_low = double(),
      ci_high = double()
    ))
  }

  out <- dplyr::bind_rows(rows)
  out
}



