# tests/testthat/helper-synthetic-data.R
#
# Shared synthetic-input helpers used across multiple test files.
# testthat auto-sources every helper-*.R at the top of every test run,
# so constants and builders defined here are visible everywhere without
# explicit source() calls.
#
# Anything tightly coupled to a single function-under-test belongs in
# its own helper-<function>.R file; this file is for cross-test fixtures.

# ---- Shared RNG seed ----
# Used by every test that calls a function exercising stochastic code
# (e.g. mr_weighted_mode bootstrap, enrichment Monte Carlo permutation).
# Single source of truth so seeds cannot silently drift across tests.
SYNTH_SEED <- 20260501L

# ---- Synthetic results_df builder (PR2 placeholder) ----
# make_synthetic_results_df() will be added in PR2 to support tests for
# enrichment_business_logic. It will produce a tibble matching the shape
# of mr_business_logic()'s results_df output, augmented with the
# ARD_selected and cause_level_* columns enrichment expects.
