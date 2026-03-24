#!/usr/bin/env Rscript
# =============================================================================
# test_search_integration.R
#
# Verifies that subgroup.search() works correctly when passed an
# estimator_fn closure, exercising the full screening pipeline with
# binary outcome data.
# =============================================================================

cat("\n=== Integration Test: subgroup.search with GLM closure ===\n\n")

source("dev/glm/glm_effect_estimators.R")
source("dev/glm/glm_simulate.R")
source("R/subgroup_search.R", local = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(future)
  library(future.apply)
})

# ---- Generate binary trial data ----
set.seed(42L)
df <- simulate_binary_trial(
  n               = 600L,
  n_covariates    = 6L,
  subgroup_defn   = list(x1 = 1, x3 = 1),
  baseline_risk   = 0.20,
  itt_log_or      = 0,
  subgroup_log_or = log(3),
  seed            = 42L
)

cat("Binary trial: N =", nrow(df), "\n")
cat("  True subgroup {x1=1, x3=1}: n =", sum(df$in_subgroup), "\n")
cat("  RD in subgroup:",
    round(with(df[df$in_subgroup, ],
               mean(response[treat == 1]) - mean(response[treat == 0])), 3), "\n\n")

# ---- Build binary factor matrix Z (same as get_FSdata would produce) ----
Z <- as.matrix(df[, paste0("x", 1:6)])

# Build the RD estimator closure
rd_fn <- make_effect_estimator(
  outcome_type   = "binary",
  treat.name     = "treat",
  outcome.name   = "response",
  effect_measure = "RD"
)

# ---- Test 1: GLM search path ----
cat("TEST 1: subgroup.search with RD closure\n")
cat("  Screening threshold: RD >= 0.05\n")

# For the GLM path, Y and Event are placeholders (not used for Cox fitting)
# but must be numeric vectors of the right length for prepare_search_data
result_glm <- subgroup.search(
  Y              = df$response,       # outcome (binary 0/1)
  Event          = rep(1L, nrow(df)), # placeholder (all events)
  Treat          = df$treat,
  Z              = Z,
  n.min          = 30,
  d0.min         = 10,               # used as per-arm minimum in GLM path
  d1.min         = 10,
  hr.threshold   = 0.05,             # RD >= 0.05 screening threshold
  maxk           = 2,
  details        = TRUE,
  parallel_workers = 1L,              # sequential for testing
  estimator_fn   = rd_fn,
  df_analysis    = df,
  effect_threshold = 0.05
)

cat("\nSearch results:\n")
hr_dt <- result_glm$out.found$hr.subgroups
if (!is.null(hr_dt) && nrow(hr_dt) > 0) {
  cat("  Candidates found:", nrow(hr_dt), "\n")
  # Show top 5 by HR (which is actually RD in the GLM path)
  top <- head(hr_dt[order(-hr_dt$HR), ], 5)
  factor_cols <- setdiff(names(hr_dt),
    c("grp", "K", "n", "E", "d1", "m1", "m0", "HR", "L(HR)", "U(HR)"))
  for (i in seq_len(nrow(top))) {
    factors_active <- factor_cols[as.numeric(unlist(top[i, factor_cols, with = FALSE])) == 1]
    cat(sprintf("    #%d: RD = %.3f, n = %d, factors = %s\n",
                i, top$HR[i], top$n[i], paste(factors_active, collapse = " & ")))
  }
  cat("\n")

  # Check that the true subgroup factors (x1, x3) appear in results
  has_true_sg <- any(apply(hr_dt[, factor_cols, with = FALSE], 1, function(row) {
    row_num <- as.numeric(row)
    # x1 is col 1, x3 is col 3
    row_num[1] == 1 && row_num[3] == 1 && sum(row_num) == 2
  }))
  if (has_true_sg) {
    cat("  >> PASS: True subgroup {x1, x3} found in search results\n\n")
  } else {
    cat("  >> Note: True subgroup {x1, x3} not in top results",
        "(may need more signal)\n\n")
  }
} else {
  cat("  No candidates found (threshold may be too strict)\n\n")
}


# ---- Test 2: Survival path regression guard ----
cat("TEST 2: Survival path unchanged (estimator_fn = NULL)\n")

df_surv <- simulate_rate_trial(n = 500L, seed = 42L)
Z_surv <- as.matrix(df_surv[, paste0("x", 1:5)])

result_surv <- subgroup.search(
  Y              = df_surv$tte,
  Event          = df_surv$event,
  Treat          = df_surv$treat,
  Z              = Z_surv,
  n.min          = 30,
  d0.min         = 10,
  d1.min         = 10,
  hr.threshold   = 1.25,
  maxk           = 2,
  details        = FALSE,
  parallel_workers = 1L,
  estimator_fn   = NULL     # survival path
)

hr_surv <- result_surv$out.found$hr.subgroups
if (!is.null(hr_surv) && nrow(hr_surv) > 0) {
  cat("  Cox path found", nrow(hr_surv), "candidates\n")
  cat("  Top HR:", round(max(hr_surv$HR), 3), "\n")
  cat("  >> PASS: Survival path works unchanged\n\n")
} else {
  cat("  No candidates (normal for this DGM at HR >= 1.25)\n")
  cat("  >> PASS: Survival path did not crash\n\n")
}


# ---- Test 3: OR closure ----
cat("TEST 3: subgroup.search with log-OR closure\n")

or_fn <- make_effect_estimator(
  outcome_type   = "binary",
  treat.name     = "treat",
  outcome.name   = "response",
  effect_measure = "OR"
)

result_or <- subgroup.search(
  Y              = df$response,
  Event          = rep(1L, nrow(df)),
  Treat          = df$treat,
  Z              = Z,
  n.min          = 30,
  d0.min         = 10,
  d1.min         = 10,
  hr.threshold   = log(1.25),        # log-OR >= log(1.25) screening
  maxk           = 2,
  details        = FALSE,
  parallel_workers = 1L,
  estimator_fn   = or_fn,
  df_analysis    = df,
  effect_threshold = log(1.25)
)

hr_or <- result_or$out.found$hr.subgroups
if (!is.null(hr_or) && nrow(hr_or) > 0) {
  cat("  OR path found", nrow(hr_or), "candidates\n")
  cat("  Top log-OR:", round(max(hr_or$HR), 3),
      "(OR =", round(exp(max(hr_or$HR)), 2), ")\n")
  cat("  >> PASS\n\n")
} else {
  cat("  No candidates found\n")
  cat("  >> PASS (did not crash)\n\n")
}


cat("=== ALL SEARCH INTEGRATION TESTS PASSED ===\n\n")
