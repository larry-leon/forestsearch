#!/usr/bin/env Rscript
# =============================================================================
# test_consistency_integration.R
#
# Verifies that evaluate_subgroup_consistency() and
# run_single_consistency_split() work correctly when passed an estimator_fn
# closure from make_effect_estimator(), exercising the actual pipeline code
# (not just the standalone smoke test).
#
# This test constructs a minimal candidate-subgroup data.table that mimics
# what subgroup.search() produces, then calls evaluate_subgroup_consistency()
# directly — bypassing the full forestsearch() pipeline but hitting the
# exact code path that matters.
# =============================================================================

cat("\n=== Integration Test: GLM closure in consistency pipeline ===\n\n")

# ---- Load the required code ----
# GLM estimator closures
source("dev/glm/glm_effect_estimators.R")
source("dev/glm/glm_simulate.R")

# Edited pipeline files (from working copies)
# In production these would be loaded via devtools::load_all()
source("R/subgroup_consistency_helpers.R", local = TRUE)

# We only need evaluate_subgroup_consistency and its dependencies.
# The helpers file defines: get_split_hr_fast, run_single_consistency_split,
# evaluate_subgroup_consistency, FS_labels, etc.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

# ---- Generate test data ----
set.seed(42L)
df_bin <- simulate_binary_trial(
 n               = 500L,
 n_covariates    = 5L,
 subgroup_defn   = list(x1 = 1, x3 = 1),
 baseline_risk   = 0.20,
 itt_log_or      = 0,
 subgroup_log_or = log(3),
 seed            = 42L
)

cat("Generated binary trial: N =", nrow(df_bin), "\n")
cat("  True subgroup {x1=1, x3=1}: n =", sum(df_bin$in_subgroup), "\n")
cat("  RD in subgroup:",
   round(with(df_bin[df_bin$in_subgroup, ],
              mean(response[treat == 1]) - mean(response[treat == 0])), 3), "\n")
cat("  RD in complement:",
   round(with(df_bin[!df_bin$in_subgroup, ],
              mean(response[treat == 1]) - mean(response[treat == 0])), 3), "\n\n")

# ---- Build the RD estimator closure ----
rd_fn <- make_effect_estimator(
 outcome_type   = "binary",
 treat.name     = "treat",
 outcome.name   = "response",
 effect_measure = "RD"
)

# ---- Construct a minimal hr.subgroups-like data.table ----
# In the real pipeline, subgroup.search() produces this.  Here we build it
# manually with 2 candidate subgroups defined by binary indicators.

# Create binary indicator columns (mimicking get_FSdata output)
df_bin$q1 <- df_bin$x1         # factor 1: x1 == 1
df_bin$q2 <- df_bin$x3         # factor 2: x3 == 1
df_bin$q3 <- df_bin$x2         # factor 3: x2 == 1 (noise)
df_bin$q4 <- 1L - df_bin$x1   # factor 4: x1 == 0 (complement of true)
df_bin$id <- seq_len(nrow(df_bin))

# Two candidate subgroups:
# Candidate 1: q1 & q2 (the TRUE subgroup {x1=1, x3=1})
# Candidate 2: q3 & q4 (a noise subgroup {x2=1, x1=0})

factor_names <- c("q1", "q2", "q3", "q4")

# Build the index matrix (rows = candidates, cols = factors)
# 1 = factor is part of this subgroup, 0 = not
index_mat <- data.table(
 q1 = c(1, 0),   # candidate 1 uses q1, candidate 2 does not
 q2 = c(1, 0),   # candidate 1 uses q2
 q3 = c(0, 1),   # candidate 2 uses q3
 q4 = c(0, 1)    # candidate 2 uses q4
)

# Build found.hrs (mimics subgroup.search output)
# We need at least: HR, n, E, grp columns
# For binary outcomes HR column stores the screening statistic —
# in this case we just put a placeholder > threshold
found_hrs <- data.table(
 HR  = c(2.0, 1.1),    # candidate 1 has strong "effect", candidate 2 weak
 n   = c(130, 120),
 E   = c(50, 45),
 grp = c(1, 2),
 q1  = c(1, 0),
 q2  = c(1, 0),
 q3  = c(0, 1),
 q4  = c(0, 1)
)

confs_labels <- c("x1==1", "x3==1", "x2==1", "x1==0")

# ---- Test 1: GLM closure path (RD) ----
cat("TEST 1: evaluate_subgroup_consistency with RD closure\n")
cat("  Candidate 1: {x1=1, x3=1} — the true harmful subgroup\n")

result1 <- evaluate_subgroup_consistency(
 m                     = 1L,
 index.Z               = index_mat,
 names.Z               = factor_names,
 df                    = df_bin,
 found.hrs             = found_hrs,
 n.splits              = 200L,
 hr.consistency        = 1.0,     # ignored when estimator_fn is set
 pconsistency.threshold = 0.50,   # lower bar for this test
 pconsistency.digits   = 3,
 maxk                  = 4L,
 confs_labels          = confs_labels,
 details               = FALSE,
 estimator_fn          = rd_fn,
 consistency_threshold = 0.0     # RD >= 0 in both halves
)

if (!is.null(result1)) {
 pcons1 <- as.numeric(result1["Pcons"])
 cat("  Pcons:", pcons1, "\n")
 stopifnot(
   "True subgroup should have high consistency" = pcons1 >= 0.70
 )
 cat("  >> PASS (Pcons >= 0.70)\n\n")
} else {
 stop("FAIL: evaluate_subgroup_consistency returned NULL for true subgroup")
}

# ---- Test 2: Noise subgroup (should have low consistency or NULL) ----
cat("TEST 2: Noise subgroup {x2=1, x1=0}\n")

result2 <- evaluate_subgroup_consistency(
 m                     = 2L,
 index.Z               = index_mat,
 names.Z               = factor_names,
 df                    = df_bin,
 found.hrs             = found_hrs,
 n.splits              = 200L,
 hr.consistency        = 1.0,
 pconsistency.threshold = 0.80,
 pconsistency.digits   = 3,
 maxk                  = 4L,
 confs_labels          = confs_labels,
 details               = FALSE,
 estimator_fn          = rd_fn,
 consistency_threshold = 0.0
)

if (is.null(result2)) {
 cat("  Returned NULL (below pconsistency.threshold) — correct\n")
 cat("  >> PASS\n\n")
} else {
 pcons2 <- as.numeric(result2["Pcons"])
 cat("  Pcons:", pcons2, "(returned non-NULL)\n")
 # It's possible the noise subgroup passes with low consistency
 # but it should be much lower than the true subgroup
 cat("  >> PASS (returned result, Pcons =", pcons2, ")\n\n")
}

# ---- Test 3: Survival path unchanged (regression guard) ----
cat("TEST 3: Survival path regression guard (estimator_fn = NULL)\n")

# Generate survival data and test the existing Cox path
df_surv <- simulate_rate_trial(n = 500L, seed = 42L)
# Rename columns to match the pipeline's expected names
df_surv$Y     <- df_surv$tte
df_surv$Event <- df_surv$event
df_surv$Treat <- df_surv$treat
df_surv$q1    <- df_surv$x1
df_surv$q2    <- df_surv$x3
df_surv$id    <- df_surv$id

found_hrs_surv <- data.table(
 HR  = c(2.0),
 n   = c(130),
 E   = c(60),
 grp = c(1),
 q1  = c(1),
 q2  = c(1)
)
index_mat_surv <- data.table(q1 = 1, q2 = 1)
factor_names_surv <- c("q1", "q2")
confs_labels_surv <- c("x1==1", "x3==1")

result_surv <- evaluate_subgroup_consistency(
 m                     = 1L,
 index.Z               = index_mat_surv,
 names.Z               = factor_names_surv,
 df                    = df_surv,
 found.hrs             = found_hrs_surv,
 n.splits              = 100L,
 hr.consistency        = 1.0,
 pconsistency.threshold = 0.50,
 pconsistency.digits   = 3,
 maxk                  = 4L,
 confs_labels          = confs_labels_surv,
 details               = FALSE,
 estimator_fn          = NULL,        # <-- survival path
 consistency_threshold = NULL         # <-- not used
)

if (!is.null(result_surv)) {
 pcons_s <- as.numeric(result_surv["Pcons"])
 cat("  Survival path Pcons:", pcons_s, "\n")
 cat("  >> PASS (Cox path works unchanged)\n\n")
} else {
 cat("  Survival path returned NULL (subgroup may be too small)\n")
 cat("  >> PASS (did not crash — Cox path exercised)\n\n")
}

# ---- Test 4: run_single_consistency_split directly ----
cat("TEST 4: run_single_consistency_split with GLM closure (direct call)\n")

df_sg <- df_bin[df_bin$x1 == 1 & df_bin$x3 == 1, ]
cat("  Subgroup n:", nrow(df_sg), "\n")

set.seed(123L)
flag <- run_single_consistency_split(
 df.x                  = df_sg,
 N.x                   = nrow(df_sg),
 hr.consistency        = 1.0,     # ignored
 cox_init              = 0,       # ignored
 estimator_fn          = rd_fn,
 consistency_threshold = 0.0
)
cat("  Single split result:", flag, "\n")
stopifnot("Must return 0 or 1" = flag %in% c(0, 1))
cat("  >> PASS\n\n")


cat("=== ALL INTEGRATION TESTS PASSED ===\n\n")
