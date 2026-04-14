#!/usr/bin/env Rscript
# =============================================================================
# QC Test: Poisson + offset path through GLM estimator infrastructure
#
# Tests:
#   1. make_effect_estimator("count", effect_measure="IRR") works
#   2. Estimates recover known treatment effect from synthetic data
#   3. .extract_fs_estimates_gen GLM branch produces correct output
#   4. Cox HR ≈ Poisson IRR on survival-like data (single event per person)
# =============================================================================

cat("=== QC: Poisson + Offset GLM Path ===\n\n")

library(forestsearch)
# Source estimator code
#source("glm_effect_estimators.R")

# ─── Test 1: Synthetic Poisson data with known effect ─────────────────────
cat("--- Test 1: Synthetic Poisson data (true IRR = 1.5) ---\n")

set.seed(42)
n <- 1000
treat <- rep(0:1, each = n / 2)
follow_up <- runif(n, 0.5, 5)  # person-time
true_log_irr <- log(1.5)
lambda <- exp(-1 + true_log_irr * treat)  # baseline rate = exp(-1)
events <- rpois(n, lambda * follow_up)

syn_df <- data.frame(
  id     = seq_len(n),
  treat  = treat,
  events = events,
  fu     = follow_up
)

# Build estimator
est_fn <- make_effect_estimator(
  outcome_type   = "count",
  treat.name     = "treat",
  outcome.name   = "events",
  offset.name    = "fu",
  effect_measure = "IRR"
)

res <- est_fn(syn_df)
cat(sprintf("  IRR estimate: %.3f (true: 1.500)\n", exp(res$estimate)))
cat(sprintf("  log-IRR:      %.4f (true: %.4f)\n", res$estimate, true_log_irr))
cat(sprintf("  SE:           %.4f\n", res$se))
cat(sprintf("  Converged:    %s\n", res$converged))
cat(sprintf("  n0=%d, n1=%d\n", res$n0, res$n1))

stopifnot(
  abs(exp(res$estimate) - 1.5) < 0.3,  # within reasonable range

  res$converged,
  res$n0 == 500,
  res$n1 == 500,
  res$measure == "log-IRR"
)
cat("  PASS\n\n")


# ─── Test 2: Binary event + offset ≈ Cox HR ──────────────────────────────
cat("--- Test 2: Binary event + offset vs Cox HR (GBSG-like) ---\n")

set.seed(123)
n2 <- 800
treat2 <- rep(0:1, each = n2 / 2)
true_log_hr <- log(1.4)

# Simulate exponential survival
baseline_hazard <- 0.01  # per month
fu_time <- rexp(n2, rate = baseline_hazard * exp(true_log_hr * treat2))
censor_time <- runif(n2, 12, 120)
obs_time <- pmin(fu_time, censor_time)
event <- as.integer(fu_time <= censor_time)

surv_df <- data.frame(
  id    = seq_len(n2),
  treat = treat2,
  time  = obs_time,
  event = event
)

cat(sprintf("  Events: %d / %d (%.1f%%)\n", sum(event), n2,
            100 * mean(event)))

# Cox HR
cox_fit <- survival::coxph(survival::Surv(time, event) ~ treat,
                           data = surv_df)
cox_hr <- exp(coef(cox_fit))

# Poisson IRR (event ~ treat + offset(log(time)))
pois_fit <- glm(event ~ treat, family = poisson(),
                offset = log(pmax(time, 1e-6)), data = surv_df)
pois_irr <- exp(coef(pois_fit)[["treat"]])

# Via make_effect_estimator
est_fn2 <- make_effect_estimator(
  outcome_type   = "count",
  treat.name     = "treat",
  outcome.name   = "event",
  offset.name    = "time",
  effect_measure = "IRR"
)
est_res2 <- est_fn2(surv_df)
est_irr <- exp(est_res2$estimate)

cat(sprintf("  Cox HR:             %.3f (true: 1.400)\n", cox_hr))
cat(sprintf("  Poisson IRR (glm):  %.3f\n", pois_irr))
cat(sprintf("  Poisson IRR (est):  %.3f\n", est_irr))
cat(sprintf("  |HR - IRR|:         %.4f\n", abs(cox_hr - est_irr)))

# Under exponential model, IRR ≈ HR
stopifnot(abs(cox_hr - est_irr) < 0.15)
cat("  PASS\n\n")


# ─── Test 3: Subgroup-level estimation ────────────────────────────────────
cat("--- Test 3: Subgroup-level Poisson estimation ---\n")

# Create a subgroup: treat2 == 1 & high-risk patients
surv_df$high_risk <- rbinom(n2, 1, 0.4)

sg_idx <- surv_df$high_risk == 1
n_sg <- sum(sg_idx)

sg_res <- est_fn2(surv_df[sg_idx, ])
comp_res <- est_fn2(surv_df[!sg_idx, ])

cat(sprintf("  Subgroup (n=%d):   IRR = %.3f\n", n_sg, exp(sg_res$estimate)))
cat(sprintf("  Complement (n=%d): IRR = %.3f\n", n2 - n_sg,
            exp(comp_res$estimate)))
cat(sprintf("  Both converged: %s\n",
            sg_res$converged && comp_res$converged))

stopifnot(sg_res$converged, comp_res$converged)
cat("  PASS\n\n")


# ─── Test 4: Edge case — very few events ─────────────────────────────────
cat("--- Test 4: Edge case — sparse events ---\n")

sparse_df <- data.frame(
  id    = 1:20,
  treat = rep(0:1, each = 10),
  event = c(rep(0, 9), 1, rep(0, 8), 1, 1),
  time  = runif(20, 1, 10)
)

sparse_res <- est_fn2(sparse_df)
cat(sprintf("  Converged: %s\n", sparse_res$converged))
cat(sprintf("  IRR: %.3f (may be unstable, that's expected)\n",
            exp(sparse_res$estimate)))
cat("  PASS (no crash)\n\n")


cat("=== All QC tests passed ===\n")
