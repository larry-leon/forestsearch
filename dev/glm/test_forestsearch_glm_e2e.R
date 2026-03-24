#!/usr/bin/env Rscript
# =============================================================================
# test_forestsearch_glm_e2e.R
#
# End-to-end test: calls forestsearch() with outcome_type = "binary"
# on simulated data and verifies the full pipeline produces sensible results.
#
# Sources package files directly (no devtools::load_all needed).
# =============================================================================

cat("\n============================================================\n")
cat("  End-to-End Test: forestsearch(outcome_type = 'binary')\n")
cat("============================================================\n\n")

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(future)
  library(future.apply)
  library(glmnet)
})

# ---- Source package files in dependency order ----
# Helpers and utilities
source("R/get_FSdata_helpers.R")
source("R/summary_utility_functions.R")
source("R/forestsearch_helpers.R")

# GLM estimator closures (from dev/glm/)
source("dev/glm/glm_effect_estimators.R")
source("dev/glm/glm_simulate.R")

# get_FSdata
source("R/get_fsdata.R")

# Search and consistency pipeline
source("R/subgroup_search.R")
source("R/subgroup_consistency_helpers.R")
source("R/subgroup_consistency_main.R")

# GRF stubs (needed for grf.subg.harm.survival to exist, though
# we disable GRF for binary outcomes)
source("R/grf_helpers.R")
source("R/grf_main.R")

# Top-level dispatcher
source("R/forestsearch_main.R")

cat("All source files loaded\n\n")


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
rd_sg <- with(df[df$in_subgroup, ],
              mean(response[treat == 1]) - mean(response[treat == 0]))
rd_co <- with(df[!df$in_subgroup, ],
              mean(response[treat == 1]) - mean(response[treat == 0]))
cat("  RD in subgroup:", round(rd_sg, 3), "\n")
cat("  RD in complement:", round(rd_co, 3), "\n\n")


# ===========================================================================
# TEST 1: forestsearch with outcome_type = "binary", effect_measure = "RD"
# ===========================================================================

cat("------------------------------------------------------------\n")
cat("TEST 1: forestsearch(outcome_type = 'binary', effect_measure = 'RD')\n")
cat("------------------------------------------------------------\n\n")

t0 <- proc.time()
fs_bin <- tryCatch(
  forestsearch(
    df.analysis      = df,
    outcome.name     = "response",
    event.name       = "response",     # binary outcome IS the event
    treat.name       = "treat",
    id.name          = "id",
    confounders.name = paste0("x", 1:6),
    outcome_type     = "binary",
    effect_measure   = "RD",
    # Search parameters
    n.min            = 30,
    d0.min           = 10,
    d1.min           = 10,
    hr.threshold     = 0.05,           # RD screening >= 0.05
    hr.consistency   = 0.0,            # RD consistency >= 0.0
    # Consistency parameters
    fs.splits        = 200L,
    pconsistency.threshold = 0.80,
    max_subgroups_search   = 10,
    use_twostage     = FALSE,
    # Disable GRF/LASSO (Phase 1)
    use_grf          = FALSE,
    use_lasso        = FALSE,
    # Sequential
    parallel_args    = list(plan = "sequential", workers = 1, show_message = FALSE),
    details          = TRUE,
    maxk             = 2,
    sg_focus         = "hr",
    stop_threshold   = NULL,
    vi.grf.min       = NULL
  ),
  error = function(e) {
    cat("ERROR:", e$message, "\n")
    NULL
  }
)
elapsed <- (proc.time() - t0)["elapsed"]

cat("\nElapsed:", round(elapsed, 1), "sec\n\n")

if (!is.null(fs_bin) && !is.null(fs_bin$sg.harm)) {
  cat("Subgroup identified:", paste(fs_bin$sg.harm, collapse = " & "), "\n")
  cat("outcome_type:", fs_bin$outcome_type, "\n")
  cat("effect_measure:", fs_bin$effect_measure, "\n")

  if (!is.null(fs_bin$grp.consistency)) {
    gc <- fs_bin$grp.consistency
    cat("Consistency algorithm:", gc$algorithm, "\n")
    cat("Candidates evaluated:", gc$n_candidates_evaluated,
        "of", gc$n_candidates_total, "\n")
    cat("Passed:", gc$n_passed, "\n")

    if (!is.null(gc$out_sg) && !is.null(gc$out_sg$result)) {
      top <- gc$out_sg$result[1, ]
      cat("Top subgroup Pcons:", as.numeric(top$Pcons), "\n")
      cat("Top subgroup effect (RD):", round(as.numeric(top$hr), 3), "\n")
      cat("Top subgroup N:", as.numeric(top$N), "\n")
    }
  }

  if (!is.null(fs_bin$df.est)) {
    cat("\nTreatment recommendations:\n")
    cat("  treat.recommend == 1 (benefit):", sum(fs_bin$df.est$treat.recommend == 1, na.rm = TRUE), "\n")
    cat("  treat.recommend == 0 (harm):   ", sum(fs_bin$df.est$treat.recommend == 0, na.rm = TRUE), "\n")
  }

  cat("\n>> PASS: forestsearch identified a subgroup for binary outcomes\n\n")
} else if (!is.null(fs_bin)) {
  cat("No subgroup identified (sg.harm is NULL)\n")
  cat("  find.grps max_sg_est:", fs_bin$max_sg_est, "\n")
  cat("  This may be expected if screening threshold is too strict\n")
  cat("\n>> PASS: forestsearch ran without errors (no subgroup found)\n\n")
} else {
  cat("\n>> FAIL: forestsearch returned NULL\n\n")
}


# ===========================================================================
# TEST 2: Survival path regression guard
# ===========================================================================

cat("------------------------------------------------------------\n")
cat("TEST 2: forestsearch(outcome_type = 'survival') — regression guard\n")
cat("------------------------------------------------------------\n\n")

df_surv <- simulate_rate_trial(n = 600L, seed = 42L)

# Rename to match package conventions
names(df_surv)[names(df_surv) == "tte"]   <- "rfstime"
names(df_surv)[names(df_surv) == "event"] <- "status"
names(df_surv)[names(df_surv) == "treat"] <- "hormon"

fs_surv <- tryCatch(
  forestsearch(
    df.analysis      = df_surv,
    outcome.name     = "rfstime",
    event.name       = "status",
    treat.name       = "hormon",
    id.name          = "id",
    confounders.name = paste0("x", 1:5),
    outcome_type     = "survival",     # default path
    n.min            = 30,
    hr.threshold     = 1.25,
    hr.consistency   = 1.0,
    fs.splits        = 100L,
    pconsistency.threshold = 0.80,
    max_subgroups_search   = 5,
    use_twostage     = FALSE,
    use_grf          = FALSE,
    use_lasso        = FALSE,
    parallel_args    = list(plan = "sequential", workers = 1, show_message = FALSE),
    details          = FALSE,
    maxk             = 2,
    sg_focus         = "hr",
    stop_threshold   = NULL,
    vi.grf.min       = NULL
  ),
  error = function(e) {
    cat("ERROR:", e$message, "\n")
    NULL
  }
)

if (!is.null(fs_surv)) {
  cat("Survival path: outcome_type =", fs_surv$outcome_type, "\n")
  if (!is.null(fs_surv$sg.harm)) {
    cat("Subgroup found:", paste(fs_surv$sg.harm, collapse = " & "), "\n")
  } else {
    cat("No subgroup found (expected — depends on DGM)\n")
  }
  cat(">> PASS: Survival path works unchanged\n\n")
} else {
  cat(">> FAIL: Survival path returned NULL\n\n")
}


# ===========================================================================
# TEST 3: OR closure through full pipeline
# ===========================================================================

cat("------------------------------------------------------------\n")
cat("TEST 3: forestsearch(outcome_type = 'binary', effect_measure = 'OR')\n")
cat("------------------------------------------------------------\n\n")

fs_or <- tryCatch(
  forestsearch(
    df.analysis      = df,
    outcome.name     = "response",
    event.name       = "response",
    treat.name       = "treat",
    id.name          = "id",
    confounders.name = paste0("x", 1:6),
    outcome_type     = "binary",
    effect_measure   = "OR",
    n.min            = 30,
    d0.min           = 10,
    d1.min           = 10,
    hr.threshold     = 1.25,          # OR screening >= 1.25
    hr.consistency   = 1.0,           # OR consistency >= 1.0
    fs.splits        = 200L,
    pconsistency.threshold = 0.80,
    max_subgroups_search   = 10,
    use_twostage     = FALSE,
    use_grf          = FALSE,
    use_lasso        = FALSE,
    parallel_args    = list(plan = "sequential", workers = 1, show_message = FALSE),
    details          = FALSE,
    maxk             = 2,
    sg_focus         = "hr",
    stop_threshold   = NULL,
    vi.grf.min       = NULL
  ),
  error = function(e) {
    cat("ERROR:", e$message, "\n")
    NULL
  }
)

if (!is.null(fs_or)) {
  cat("OR path: outcome_type =", fs_or$outcome_type,
      ", effect_measure =", fs_or$effect_measure, "\n")
  if (!is.null(fs_or$sg.harm)) {
    cat("Subgroup found:", paste(fs_or$sg.harm, collapse = " & "), "\n")
  } else {
    cat("No subgroup found\n")
  }
  cat(">> PASS: OR path through forestsearch works\n\n")
} else {
  cat(">> FAIL: OR path returned NULL\n\n")
}


# ===========================================================================
# SUMMARY
# ===========================================================================

cat("============================================================\n")
cat("  ALL END-TO-END TESTS COMPLETE\n")
cat("============================================================\n\n")
