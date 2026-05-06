# dev/replot_coloc_sensitivity.R
#
# Re-renders a coloc sensitivity plot (PNG + PDF) from the per-locus CSVs
# the pipeline writes to <run_dir>/<outcome_slug>/coloc/<locus>/. Useful
# for tweaking plot dimensions / trait labels / margins without re-running
# the whole pipeline.
#
# The script reads two CSVs from this dev/ directory:
#
#   dev/coloc_aligned_data.csv  -- per-SNP aligned betas + SEs for both traits
#   dev/coloc_abf_summary.csv   -- locus-level PP.H0..PP.H4 + top SNP
#
# It reconstructs the two coloc datasets D1 / D2, re-runs coloc::coloc.abf,
# and feeds the result into the same `coloc_plot_sensitivity()` helper that
# the pipeline calls in production. The PNG + PDF land back in dev/.
#
# Edit the trait labels / N / type below before running. The defaults match
# the TNF-alpha (Sun et al. 2018 INTERVAL pQTL) -> A04 Other bacterial
# intestinal infections smoke-test that produced the example CSVs.
#
# Usage: source("dev/replot_coloc_sensitivity.R") from the package root.

stopifnot(file.exists("DESCRIPTION"))
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required. install.packages('devtools').")
}
devtools::load_all(".", quiet = TRUE)

# ---- inputs (edit if you have other CSVs in dev/) ---------------------
aligned_csv <- "dev/coloc_aligned_data.csv"
summary_csv <- "dev/coloc_abf_summary.csv"
out_stem    <- "dev/coloc_sensitivity_replot"

# Trait labels (replaces "trait 1" / "trait 2" in the panel titles).
exp_label <- "TNF-alpha (prot-c-3722_49_2)"
out_label <- "A04 Other bacterial intestinal infections"

# coloc inputs that the aligned-data CSV doesn't preserve. Update these
# when re-plotting from a different exposure / outcome.
exp_N    <- 3301L          # Sun et al. 2018 INTERVAL TNF-alpha pQTL
exp_type <- "quant"
out_N    <- 420531L        # from the run log: Pan-UKB EUR n_total for A04
out_type <- "cc"
out_s    <- 196 / 420531   # n_cases_EUR / n_total for A04 -- adjust if known

# ---- load ------------------------------------------------------------
aligned <- utils::read.csv(aligned_csv, stringsAsFactors = FALSE)
sum_tbl <- utils::read.csv(summary_csv, stringsAsFactors = FALSE)
stopifnot(nrow(aligned) >= 5L)
message(sprintf("Loaded %d SNPs from %s", nrow(aligned), aligned_csv))
message(sprintf("Locus: %s; cached PP.H4 = %.4f",
                sum_tbl$locus_id[1], sum_tbl$PP.H4[1]))

# ---- reconstruct coloc datasets --------------------------------------
# Mirror the MAF-floor and snp-id construction from .coloc_run_one_locus().
maf_exp <- pmin(aligned$eaf_exp, 1 - aligned$eaf_exp)
maf_out <- pmin(aligned$eaf_out, 1 - aligned$eaf_out)
maf_exp[is.na(maf_exp) | maf_exp <= 0] <- 0.05
maf_out[is.na(maf_out) | maf_out <= 0] <- 0.05
snp_ids <- paste0("chr", aligned$chr, "_", aligned$pos)

D1 <- list(
  beta = aligned$beta_exp, varbeta = aligned$se_exp ^ 2, MAF = maf_exp,
  N = exp_N, type = exp_type, snp = snp_ids, position = aligned$pos
)
D2 <- list(
  beta = aligned$beta_out, varbeta = aligned$se_out ^ 2, MAF = maf_out,
  N = out_N, type = out_type, snp = snp_ids, position = aligned$pos
)
if (identical(out_type, "cc") && is.finite(out_s)) D2$s <- out_s

# ---- re-run coloc.abf ------------------------------------------------
set.seed(20260501L)
abf_res <- coloc::coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
new_pp_h4 <- as.numeric(abf_res$summary["PP.H4.abf"])
message(sprintf("Re-computed PP.H4 = %.4f (cached was %.4f)",
                new_pp_h4, sum_tbl$PP.H4[1]))

# ---- render via the exact pipeline plotter ---------------------------
coloc_plot_sensitivity(
  abf_res,
  basepath  = out_stem,
  exp_label = exp_label,
  out_label = out_label
)
message("Wrote: ", paste0(out_stem, c(".png", ".pdf"), collapse = " + "))
