# R/enrichment_business_logic.R

#' Threshold-free directional enrichment via label permutations (business logic)
#'
#' Computes threshold-free, directional enrichment using **label permutations**
#' on a signed evidence score \code{s = sign(beta_ivw) * (-log10 p_ivw)}.
#'
#' Provides:
#' \itemize{
#'   \item \strong{Global}: ARD vs non-ARD across the full phenome
#'         (returns protective/risk tails \emph{and} a single signed test).
#'   \item \strong{Cause-level (L1/L2/L3)}: three comparison modes, each per cause:
#'         \enumerate{
#'           \item \code{"ARD_vs_nonARD_within_cause"}: permute within the cause.
#'           \item \code{"cause_vs_rest_all"}: permute globally (cause vs rest).
#'           \item \code{"ARD_in_cause_vs_ARD_elsewhere"}: permute among ARDs only.
#'         }
#'         Each returns protective/risk tails \emph{and} a single signed test.
#' }
#'
#' For each test we report:
#' \itemize{
#'   \item Tail-specific one-sided p-values \code{p_prot}, \code{p_risk}, BH \code{q_*}, and
#'         standardized enrichment scores \code{SES_prot}, \code{SES_risk}.
#'   \item A \strong{single signed} statistic with two-sided p-value:
#'         \code{stat_signed}, \code{p_signed}, \code{q_signed}, \code{SES_signed}
#'         (negative = protective enrichment; positive = risk enrichment).
#' }
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
.choose_weights <- function(se, scheme = c("none","inv_se2")) {
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
.perm_test_sum <- function(contrib, y, exact_max_combn = 1e5, mc_B = 100000, seed = NULL) {
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
.perm_test_signed <- function(contrib_signed, y, exact_max_combn = 1e5, mc_B = 100000, seed = NULL) {
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
  list(stat_obs=stat_obs, p=p, mean_perm=mu, sd_perm=sdv, ses=ses, n_perm=n_perm, exact=do_exact)
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
#' @param weight_scheme \code{"none"} or \code{"inv_se2"} (needs \code{results_se_ivw}).
#' @param exact_max_combn Maximum \eqn{choose(n,a)} for exact enumeration; otherwise Monte Carlo.
#' @param mc_B Monte Carlo permutations if not exact.
#' @param seed RNG seed.
#' @return A list with: \code{$global_tbl} and \code{$by_cause_tbl} (all levels × modes concatenated).
#' @export
run_enrichment <- function(
    results_df,
    exposure = .default_exposure(),
    levels = c("cause_level_1","cause_level_2","cause_level_3"),
    modes  = c("ARD_vs_nonARD_within_cause","cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere"),
    use_qc_pass = TRUE,
    min_nsnp = 2,
    weight_scheme = c("none","inv_se2"),
    exact_max_combn = 1e5,
    mc_B = 100000,
    seed = NULL
) {
  weight_scheme <- match.arg(weight_scheme)

  # --- global ---
  global_tbl <- enrichment_global_directional(
    results_df = results_df,
    exposure   = exposure,
    use_qc_pass = use_qc_pass,
    min_nsnp    = min_nsnp,
    weight_scheme = weight_scheme,
    exact_max_combn = exact_max_combn,
    mc_B = mc_B,
    seed = seed
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
        seed = seed
      )
      by_cause[[paste(lv, md, sep = "::")]] <- tmp
    }
  }
  by_cause_tbl <- dplyr::bind_rows(by_cause)

  list(global_tbl = global_tbl, by_cause_tbl = by_cause_tbl)
}

#' Global ARD vs non-ARD directional enrichment (tails + signed)
#'
#' @inheritParams run_enrichment
#' @return One-row tibble with counts; tail stats (prot/risk) and a single signed stat.
#' @export
enrichment_global_directional <- function(
    results_df,
    exposure = .default_exposure(),
    use_qc_pass = TRUE,
    min_nsnp = 2,
    weight_scheme = c("none","inv_se2"),
    exact_max_combn = 1e5,
    mc_B = 100000,
    seed = NULL
) {
  weight_scheme <- match.arg(weight_scheme)

  df <- results_df
  if (use_qc_pass && "results_qc_pass" %in% names(df)) df <- dplyr::filter(df, .data$results_qc_pass)
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df)) df <- dplyr::filter(df, .data$results_nsnp_after >= !!min_nsnp)

  stopifnot(all(c("results_beta_ivw","results_p_ivw","ARD_selected") %in% names(df)))
  s <- signed_ivw_score(df$results_beta_ivw, df$results_p_ivw)
  w <- if ("results_se_ivw" %in% names(df)) .choose_weights(df$results_se_ivw, weight_scheme) else .choose_weights(rep(NA_real_, nrow(df)), "none")
  y <- .as_int01(df$ARD_selected)

  # tails (one-sided)
  c_prot <- .tail_contrib(s, w, "protective"); tp <- .perm_test_sum(c_prot, y, exact_max_combn, mc_B, seed)
  c_risk <- .tail_contrib(s, w, "risk");       tr <- .perm_test_sum(c_risk, y, exact_max_combn, mc_B, if (is.null(seed)) NULL else seed + 1L)

  # signed (two-sided)
  c_signed <- w * s
  ts <- .perm_test_signed(c_signed, y, exact_max_combn, mc_B, if (is.null(seed)) NULL else seed + 2L)

  sex_summary <- if ("pheno_sex" %in% names(df)) paste0(sort(unique(na.omit(as.character(df$pheno_sex)))), collapse = ";") else NA_character_

  out <- tibble::tibble(
    scope   = "GLOBAL",
    exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
    pheno_sex = sex_summary,
    n_total = length(y),
    n_ard   = sum(y),
    n_non   = length(y) - sum(y),
    # tails
    stat_prot = tp$stat_obs, p_prot = tp$p, SES_prot = tp$ses,
    stat_risk = tr$stat_obs, p_risk = tr$p, SES_risk = tr$ses,
    n_perm_prot = tp$n_perm, exact_prot = tp$exact,
    n_perm_risk = tr$n_perm, exact_risk = tr$exact,
    # signed
    stat_signed = ts$stat_obs, p_signed = ts$p, SES_signed = ts$ses,
    n_perm_signed = ts$n_perm, exact_signed = ts$exact
  )

  out |>
    dplyr::mutate(
      q_prot   = p.adjust(p_prot,   method = "BH"),
      q_risk   = p.adjust(p_risk,   method = "BH"),
      q_signed = p.adjust(p_signed, method = "BH")
    )
}

#' Cause-level directional enrichment (tails + signed) for one level × one mode
#'
#' @param level One of \code{"cause_level_1"}, \code{"cause_level_2"}, \code{"cause_level_3"}.
#' @param compare_mode One of:
#'   \code{"ARD_vs_nonARD_within_cause"}, \code{"cause_vs_rest_all"},
#'   \code{"ARD_in_cause_vs_ARD_elsewhere"}.
#' @inheritParams run_enrichment
#' @return Tibble: one row per cause with tail stats and a single signed stat.
#' @export
enrichment_by_cause_directional <- function(
    results_df,
    level = c("cause_level_2","cause_level_1","cause_level_3"),
    compare_mode = c("ARD_vs_nonARD_within_cause","cause_vs_rest_all","ARD_in_cause_vs_ARD_elsewhere"),
    exposure = .default_exposure(),
    use_qc_pass = TRUE,
    min_nsnp = 2,
    weight_scheme = c("none","inv_se2"),
    exact_max_combn = 1e5,
    mc_B = 100000,
    seed = NULL
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

  run_one <- function(df_all, cause_label, mode) {
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
        # tails
        stat_prot = NA_real_, p_prot = NA_real_, SES_prot = NA_real_,
        stat_risk = NA_real_, p_risk = NA_real_, SES_risk = NA_real_,
        n_perm_prot = NA_integer_, exact_prot = NA, n_perm_risk = NA_integer_, exact_risk = NA,
        # signed
        stat_signed = NA_real_, p_signed = NA_real_, SES_signed = NA_real_,
        n_perm_signed = NA_integer_, exact_signed = NA,
        perm_scope = perm_scope
      ))
    }

    # tails (one-sided)
    c_prot <- .tail_contrib(s, w, "protective")
    c_risk <- .tail_contrib(s, w, "risk")
    tp <- .perm_test_sum(c_prot, y, exact_max_combn, mc_B, seed)
    tr <- .perm_test_sum(c_risk, y, exact_max_combn, mc_B, if (is.null(seed)) NULL else seed + 1L)

    # signed (two-sided)
    c_signed <- w * s
    ts <- .perm_test_signed(c_signed, y, exact_max_combn, mc_B, if (is.null(seed)) NULL else seed + 2L)

    tibble::tibble(
      level = level, cause = cause_label, compare_mode = mode,
      exposure = if (is.null(exposure)) NA_character_ else as.character(exposure),
      pheno_sex = if ("pheno_sex" %in% names(g)) paste0(sort(unique(na.omit(as.character(g$pheno_sex)))), collapse = ";") else NA_character_,
      n_total = scope_counts["n_total"], n_pos = scope_counts["n_pos"], n_neg = scope_counts["n_neg"],
      # tails
      stat_prot = tp$stat_obs, p_prot = tp$p, SES_prot = tp$ses,
      stat_risk = tr$stat_obs, p_risk = tr$p, SES_risk = tr$ses,
      n_perm_prot = tp$n_perm, exact_prot = tp$exact,
      n_perm_risk = tr$n_perm, exact_risk = tr$exact,
      # signed
      stat_signed = ts$stat_obs, p_signed = ts$p, SES_signed = ts$ses,
      n_perm_signed = ts$n_perm, exact_signed = ts$exact,
      perm_scope = perm_scope
    )
  }

  causes <- sort(unique(df0[[level]]))
  out <- dplyr::bind_rows(lapply(causes, function(cz) run_one(df0, cz, compare_mode)))

  # BH within this family for all three tests
  out |>
    dplyr::mutate(
      q_prot   = p.adjust(p_prot,   method = "BH"),
      q_risk   = p.adjust(p_risk,   method = "BH"),
      q_signed = p.adjust(p_signed, method = "BH")
    )
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
    return(tibble::tibble(
      level = character(), cause = character(), exposure = character(),
      n_in = integer(), n_out = integer(),
      beta_in = double(), beta_out = double(),
      delta_beta = double(), se_delta = double(),
      z = double(), p = double(), q = double(),
      ci_low = double(), ci_high = double()
    ))
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
  if (!nrow(out)) return(out)

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
  if (use_qc_pass && "results_qc_pass" %in% names(df)) df <- dplyr::filter(df, .data$results_qc_pass %in% TRUE)
  if (!is.null(min_nsnp) && "results_nsnp_after" %in% names(df)) df <- dplyr::filter(df, .data$results_nsnp_after >= !!min_nsnp)
  df <- dplyr::filter(df, !is.na(.data[[level]]), is.finite(.data$results_beta_ivw),
                      is.finite(.data$results_se_ivw), .data$results_se_ivw > 0)

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
  if (!nrow(out)) return(out)

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



