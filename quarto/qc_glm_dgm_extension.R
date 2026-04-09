#!/usr/bin/env Rscript
# =============================================================================
# QC: generate_glm_dgm / simulate_from_glm_dgm / calibrate_glm_interaction
# for binary, continuous, and count (Poisson + offset) outcomes
# =============================================================================

cat("=== QC: GLM DGM Extension (binary / continuous / count) ===\n\n")


library(forestsearch)

# Source the three files under test
# source("generate_glm_dgm.R")
# source("simulate_from_glm_dgm.R")
# source("calibrate_glm_interaction.R")
# # Source dependencies from the latest codebase
# source("latest_R/R/generate_aft_dgm_helpers.R")

# ─── Synthetic data common to all tests ──────────────────────────────────
set.seed(42)
n <- 800
syn <- data.frame(
  id    = seq_len(n),
  treat = rep(0:1, each = n / 2),
  z1    = as.factor(sample(0:1, n, replace = TRUE)),
  z2    = as.factor(sample(0:1, n, replace = TRUE)),
  z3    = as.factor(sample(0:1, n, replace = TRUE)),
  x_cont = rnorm(n, 50, 10)
)

# ═══════════════════════════════════════════════════════════════════════════
# TEST 1: Binary outcome (existing functionality — regression test)
# ═══════════════════════════════════════════════════════════════════════════
cat("--- Test 1: Binary outcome (OR) ---\n")

# Generate binary outcome with treatment effect
p_base <- plogis(-1 + 0.5 * as.integer(as.character(syn$z1)) +
                   0.3 * as.integer(as.character(syn$treat)))
syn$y_bin <- rbinom(n, 1, p_base)

dgm_bin <- generate_glm_dgm(
  data          = syn,
  factor_vars   = c("z1", "z2", "z3"),
  outcome_var   = "y_bin",
  treatment_var = "treat",
  outcome_type  = "binary",
  effect_measure = "OR",
  subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model         = "alt",
  k_inter       = 1.0,
  n_super       = 3000L,
  verbose       = TRUE
)

stopifnot(
  inherits(dgm_bin, "glm_dgm"),
  dgm_bin$outcome_type == "binary",
  dgm_bin$effect_measure == "OR",
  "p0" %in% names(dgm_bin$df_super),
  "p1" %in% names(dgm_bin$df_super),
  !is.na(dgm_bin$hazard_ratios$harm_subgroup),
  !is.na(dgm_bin$hazard_ratios$overall)
)

# Simulate from the DGM
df_sim_bin <- simulate_from_glm_dgm(dgm_bin, n = 500, seed = 1)
stopifnot(
  all(df_sim_bin$y_sim %in% c(0L, 1L)),
  "treat_sim" %in% names(df_sim_bin),
  nrow(df_sim_bin) == 500
)
cat("  Binary DGM: OR(Q) =", round(dgm_bin$hazard_ratios$harm_subgroup, 3), "\n")
cat("  Simulated: n =", nrow(df_sim_bin),
    ", event rate =", round(mean(df_sim_bin$y_sim), 3), "\n")
cat("  PASS\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# TEST 2: Continuous outcome (MD)
# ═══════════════════════════════════════════════════════════════════════════
cat("--- Test 2: Continuous outcome (MD) ---\n")

syn$y_cont <- 50 + 3 * as.integer(as.character(syn$treat)) +
  2 * as.integer(as.character(syn$z1)) + rnorm(n, sd = 5)

dgm_cont <- generate_glm_dgm(
  data          = syn,
  factor_vars   = c("z1", "z2", "z3"),
  outcome_var   = "y_cont",
  treatment_var = "treat",
  outcome_type  = "continuous",
  effect_measure = "MD",
  subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model         = "alt",
  k_inter       = 5.0,   # additional 5-unit shift for Q under treatment

  n_super       = 3000L,
  verbose       = TRUE
)

stopifnot(
  inherits(dgm_cont, "glm_dgm"),
  dgm_cont$outcome_type == "continuous",
  "mu0" %in% names(dgm_cont$df_super),
  "mu1" %in% names(dgm_cont$df_super),
  !is.null(dgm_cont$model_params$sigma),
  dgm_cont$model_params$sigma > 0
)

# Check that Q subgroup has larger MD than Qc
cat("  MD(Q)  =", round(dgm_cont$hazard_ratios$harm_subgroup, 3), "\n")
cat("  MD(Qc) =", round(dgm_cont$hazard_ratios$no_harm_subgroup, 3), "\n")
cat("  MD(ITT)=", round(dgm_cont$hazard_ratios$overall, 3), "\n")
cat("  sigma  =", round(dgm_cont$model_params$sigma, 3), "\n")

stopifnot(
  dgm_cont$hazard_ratios$harm_subgroup >
    dgm_cont$hazard_ratios$no_harm_subgroup
)

# Simulate
df_sim_cont <- simulate_from_glm_dgm(dgm_cont, n = 500, seed = 2)
stopifnot(
  is.numeric(df_sim_cont$y_sim),
  nrow(df_sim_cont) == 500
)
cat("  Simulated: mean(y_sim) =", round(mean(df_sim_cont$y_sim), 2),
    ", sd =", round(sd(df_sim_cont$y_sim), 2), "\n")
cat("  PASS\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# TEST 3: Count outcome WITHOUT offset (simple Poisson)
# ═══════════════════════════════════════════════════════════════════════════
cat("--- Test 3: Count outcome without offset (IRR) ---\n")

lambda_base <- exp(1.0 + 0.3 * as.integer(as.character(syn$treat)) +
                    0.2 * as.integer(as.character(syn$z1)))
syn$y_count <- rpois(n, lambda_base)

dgm_count <- generate_glm_dgm(
  data          = syn,
  factor_vars   = c("z1", "z2", "z3"),
  outcome_var   = "y_count",
  treatment_var = "treat",
  outcome_type  = "count",
  effect_measure = "IRR",
  subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model         = "alt",
  k_inter       = 0.5,   # multiplicative: rate * exp(0.5) ≈ 1.65x
  n_super       = 3000L,
  verbose       = TRUE
)

stopifnot(
  inherits(dgm_count, "glm_dgm"),
  dgm_count$outcome_type == "count",
  "mu0" %in% names(dgm_count$df_super),
  "mu1" %in% names(dgm_count$df_super),
  all(dgm_count$df_super$mu0 > 0),
  all(dgm_count$df_super$mu1 > 0)
)

cat("  IRR(Q)  =", round(dgm_count$hazard_ratios$harm_subgroup, 3), "\n")
cat("  IRR(Qc) =", round(dgm_count$hazard_ratios$no_harm_subgroup, 3), "\n")

# Simulate
df_sim_count <- simulate_from_glm_dgm(dgm_count, n = 500, seed = 3)
stopifnot(
  all(df_sim_count$y_sim >= 0),
  all(df_sim_count$y_sim == floor(df_sim_count$y_sim))
)
cat("  Simulated: mean(y_sim) =", round(mean(df_sim_count$y_sim), 2), "\n")
cat("  PASS\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# TEST 4: Count outcome WITH offset (Poisson rate model)
# ═══════════════════════════════════════════════════════════════════════════
cat("--- Test 4: Count outcome with offset (Poisson rate model) ---\n")

syn$follow_up <- runif(n, 0.5, 5.0)
rate_base <- exp(-0.5 + 0.4 * as.integer(as.character(syn$treat)))
syn$y_events <- rpois(n, rate_base * syn$follow_up)

dgm_rate <- generate_glm_dgm(
  data          = syn,
  factor_vars   = c("z1", "z2", "z3"),
  outcome_var   = "y_events",
  treatment_var = "treat",
  outcome_type  = "count",
  effect_measure = "IRR",
  offset_var    = "follow_up",
  subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model         = "alt",
  k_inter       = 0.7,
  n_super       = 3000L,
  verbose       = TRUE
)

stopifnot(
  inherits(dgm_rate, "glm_dgm"),
  dgm_rate$model_params$offset_var == "follow_up",
  "mu0" %in% names(dgm_rate$df_super),
  "mu1" %in% names(dgm_rate$df_super)
)

cat("  IRR(Q)  =", round(dgm_rate$hazard_ratios$harm_subgroup, 3), "\n")
cat("  IRR(Qc) =", round(dgm_rate$hazard_ratios$no_harm_subgroup, 3), "\n")
cat("  IRR(ITT)=", round(dgm_rate$hazard_ratios$overall, 3), "\n")

# Verify Q has larger IRR than Qc (k_inter > 0 increases rate for Q)
stopifnot(
  dgm_rate$hazard_ratios$harm_subgroup >
    dgm_rate$hazard_ratios$no_harm_subgroup
)

# Simulate and verify counts are non-negative integers
df_sim_rate <- simulate_from_glm_dgm(dgm_rate, n = 500, seed = 4)
stopifnot(
  all(df_sim_rate$y_sim >= 0),
  all(df_sim_rate$y_sim == floor(df_sim_rate$y_sim)),
  "follow_up" %in% names(df_sim_rate)
)
cat("  Simulated: mean events =", round(mean(df_sim_rate$y_sim), 2),
    ", mean FU =", round(mean(df_sim_rate$follow_up), 2), "\n")
cat("  PASS\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# TEST 5: Null model (no interaction) — all outcome types
# ═══════════════════════════════════════════════════════════════════════════
cat("--- Test 5: Null model (all types) ---\n")

for (otype in c("binary", "continuous", "count")) {
  ovar <- switch(otype, binary = "y_bin", continuous = "y_cont", count = "y_count")
  em   <- switch(otype, binary = "OR", continuous = "MD", count = "IRR")

  dgm_null <- generate_glm_dgm(
    data          = syn,
    factor_vars   = c("z1", "z2", "z3"),
    outcome_var   = ovar,
    treatment_var = "treat",
    outcome_type  = otype,
    effect_measure = em,
    subgroup_vars = c("z1", "z2"),
    subgroup_cuts = list(z1 = 1L, z2 = 1L),
    model         = "null",
    k_inter       = 999,  # should be forced to 0 under null
    n_super       = 2000L,
    verbose       = FALSE
  )

  eff_Q  <- dgm_null$hazard_ratios$harm_subgroup
  eff_Qc <- dgm_null$hazard_ratios$no_harm_subgroup

  # Under null, Q and Qc effects should be approximately equal
  if (em %in% c("OR", "IRR")) {
    ratio <- eff_Q / eff_Qc
    ok <- abs(ratio - 1) < 0.3
  } else {
    ok <- abs(eff_Q - eff_Qc) < 3
  }

  cat(sprintf("  %s null: Effect(Q)=%.3f, Effect(Qc)=%.3f, similar=%s\n",
      otype, eff_Q, eff_Qc, ok))
  stopifnot(ok)
}
cat("  PASS\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# TEST 6: Calibration (continuous and count)
# ═══════════════════════════════════════════════════════════════════════════
cat("--- Test 6: calibrate_glm_interaction (continuous) ---\n")

dgm_cal_cont <- calibrate_glm_interaction(
  data          = syn,
  factor_vars   = c("z1", "z2", "z3"),
  outcome_var   = "y_cont",
  treatment_var = "treat",
  target_effect = 8.0,   # target MD = 8 in Q
  outcome_type  = "continuous",
  effect_measure = "MD",
  subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  k_inter_range = c(0, 15),
  grid_step     = 0.5,
  n_super       = 3000L,
  verbose       = TRUE
)

cat("  Calibrated MD(Q) =", round(dgm_cal_cont$hazard_ratios$harm_subgroup, 3),
    "(target: 8.0)\n")
stopifnot(abs(dgm_cal_cont$hazard_ratios$harm_subgroup - 8.0) < 1.5)
cat("  PASS\n\n")


cat("--- Test 6b: calibrate_glm_interaction (count with offset) ---\n")

dgm_cal_rate <- calibrate_glm_interaction(
  data          = syn,
  factor_vars   = c("z1", "z2", "z3"),
  outcome_var   = "y_events",
  treatment_var = "treat",
  target_effect = 2.0,   # target IRR = 2.0 in Q
  outcome_type  = "count",
  effect_measure = "IRR",
  offset_var    = "follow_up",
  subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  k_inter_range = c(0, 3),
  grid_step     = 0.1,
  n_super       = 3000L,
  verbose       = TRUE
)

cat("  Calibrated IRR(Q) =", round(dgm_cal_rate$hazard_ratios$harm_subgroup, 3),
    "(target: 2.0)\n")
stopifnot(abs(dgm_cal_rate$hazard_ratios$harm_subgroup - 2.0) < 0.5)
cat("  PASS\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# TEST 7: Edge case — print method
# ═══════════════════════════════════════════════════════════════════════════
cat("--- Test 7: print.glm_dgm ---\n")
print(dgm_rate)
cat("  PASS\n\n")


cat("=== All QC tests passed ===\n")
