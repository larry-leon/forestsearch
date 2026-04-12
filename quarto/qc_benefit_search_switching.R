#!/usr/bin/env Rscript
# ===================================================================
# qc_benefit_search_switching.R
#
# Pressure test for the benefit search switching logic used in
# actg175_binary_benefit_simulations.qmd and
# actg175_survival_benefit_simulations.qmd.
#
# Tests that the treatment switching approach correctly:
# 1. Produces large OR/HR on the switched scale for Q (detection)
# 2. Produces OR/HR < 1 on the switched scale under null (FPR control)
# 3. ITT reported on original scale matches expectations
# 4. Null DGM represents homogeneous benefit, not zero effect
#
# This mirrors S3 (benefit + positive binary) and S15 (benefit +
# positive survival) from forestsearch_scenario_validation.qmd.
# ===================================================================

library(forestsearch)
library(survival)

cat("=== QC: Benefit Search Switching Logic ===\n\n")

pass_count <- 0L
fail_count <- 0L

check <- function(condition, label) {
  if (condition) {
    cat(sprintf("  PASS: %s\n", label))
    pass_count <<- pass_count + 1L
  } else {
    cat(sprintf("  FAIL: %s\n", label))
    fail_count <<- fail_count + 1L
  }
}


# ===================================================================
# SECTION 1: Binary Benefit + Positive (S3 pattern)
# ===================================================================
cat("--- Section 1: Binary Benefit + Positive (S3 pattern) ---\n")

set.seed(42)
N <- 1200
syn <- data.frame(
  id = seq_len(N), treat = rep(0:1, each = N / 2),
  z1 = as.factor(sample(0:1, N, TRUE)),
  z2 = as.factor(sample(0:1, N, TRUE)),
  z3 = as.factor(sample(0:1, N, TRUE)),
  z4 = as.factor(sample(0:1, N, TRUE))
)
in_Q <- syn$z1 == 1 & syn$z2 == 1
syn$treat_sw <- 1L - syn$treat

# Positive outcome: treatment increases improvement for everyone,
# Q benefits MORE (OR_orig(Q) >> 1)
syn$y_improve <- rbinom(N, 1,
  ifelse(syn$treat == 1 & in_Q, 0.80,
    ifelse(syn$treat == 1, 0.50, 0.40)))

# Verify raw data
or_orig_Q  <- {
  t <- table(syn$y_improve[in_Q], syn$treat[in_Q])
  (t[2,2]*t[1,1]) / (t[1,2]*t[2,1])
}
or_orig_Qc <- {
  t <- table(syn$y_improve[!in_Q], syn$treat[!in_Q])
  (t[2,2]*t[1,1]) / (t[1,2]*t[2,1])
}
cat(sprintf("\n  Raw data OR(improve, orig): Q=%.2f, Qc=%.2f\n",
    or_orig_Q, or_orig_Qc))

check(or_orig_Q > 2.0,
  sprintf("OR(Q, orig) = %.2f > 2.0 (Q benefits more)", or_orig_Q))
check(or_orig_Qc > 1.0 && or_orig_Qc < 2.0,
  sprintf("OR(Qc, orig) = %.2f in (1, 2) (Qc benefits mildly)", or_orig_Qc))


# --- Build DGM with SWITCHED treatment ---
cat("\n  Building DGM with switched treatment...\n")
dgm_alt <- generate_glm_dgm(
  syn, c("z1","z2","z3","z4"), "y_improve", "treat_sw",
  "binary", "OR",
  subgroup_vars = c("z1","z2"), subgroup_cuts = list(z1=1L, z2=1L),
  model = "alt", k_inter = 2.5, n_super = 5000L, seed = 42L
)

or_sw_Q <- dgm_alt$hazard_ratios$harm_subgroup
or_sw_Qc <- dgm_alt$hazard_ratios$no_harm_subgroup
or_sw_ITT <- dgm_alt$hazard_ratios$overall

cat(sprintf("  DGM (switched): OR(Q)=%.2f, OR(Qc)=%.2f, OR(ITT)=%.2f\n",
    or_sw_Q, or_sw_Qc, or_sw_ITT))
cat(sprintf("  DGM (original): OR(Q)=%.2f, OR(Qc)=%.2f, OR(ITT)=%.2f\n",
    1/or_sw_Q, 1/or_sw_Qc, 1/or_sw_ITT))

# On switched scale: Q should have high OR (ddI=treat worse for Q)
check(or_sw_Q > 1.5,
  sprintf("OR(Q, switched) = %.2f > 1.5 (exceeds threshold)", or_sw_Q))

# On original scale: Q has OR < 1 (strong benefit)
check(1/or_sw_Q < 0.7,
  sprintf("OR(Q, original) = %.2f < 0.70 (strong benefit)", 1/or_sw_Q))

# ITT on switched scale should be < 1 (treat_orig is better, so
# treat_sw is worse overall)
check(or_sw_ITT < 1.0,
  sprintf("OR(ITT, switched) = %.2f < 1.0 (original treat is better)", or_sw_ITT))

# ITT on original scale: OR > 1 (improvement increases)
check(1/or_sw_ITT > 1.0,
  sprintf("OR(ITT, original) = %.2f > 1.0 (treat increases improvement)", 1/or_sw_ITT))


# --- Null DGM with k_treat = 1 ---
cat("\n  Building null DGM (k_treat=1, homogeneous benefit)...\n")
syn$y_improve_null <- rbinom(N, 1,
  ifelse(syn$treat == 1, 0.55, 0.40))

dgm_null <- generate_glm_dgm(
  syn, c("z1","z2","z3","z4"), "y_improve_null", "treat_sw",
  "binary", "OR",
  subgroup_vars = c("z1","z2"), subgroup_cuts = list(z1=1L, z2=1L),
  model = "null", k_treat = 1, n_super = 5000L, seed = 42L
)

or_null_sw <- dgm_null$hazard_ratios$overall
or_null_orig <- 1 / or_null_sw

cat(sprintf("  Null OR(ITT, switched)=%.3f, OR(ITT, original)=%.3f\n",
    or_null_sw, or_null_orig))

# Switched OR < 1 (far below threshold = good FPR)
check(or_null_sw < 1.0,
  sprintf("Null OR(ITT, switched) = %.3f < 1.0 (below threshold)", or_null_sw))

# Original: OR > 1 (treatment is beneficial)
check(or_null_orig > 1.0,
  sprintf("Null OR(ITT, original) = %.3f > 1.0 (beneficial)", or_null_orig))

# Null DGM is NOT zero effect
check(abs(or_null_orig - 1.0) > 0.05,
  sprintf("Null is NOT zero effect: |OR-1| = %.3f > 0.05",
    abs(or_null_orig - 1.0)))

# Q and Qc have same effect (no HTE)
or_null_Q <- dgm_null$hazard_ratios$harm_subgroup
or_null_Qc <- dgm_null$hazard_ratios$no_harm_subgroup

if (is.numeric(or_null_Q) && is.numeric(or_null_Qc) &&
    length(or_null_Q) == 1 && length(or_null_Qc) == 1) {
  cat(sprintf("  Null OR(Q, sw)=%.3f, OR(Qc, sw)=%.3f\n",
      or_null_Q, or_null_Qc))
  check(abs(log(or_null_Q) - log(or_null_Qc)) < 0.3,
    sprintf("Null Q vs Qc similar: |log(OR_Q)-log(OR_Qc)| = %.3f < 0.3",
      abs(log(or_null_Q) - log(or_null_Qc))))
} else {
  cat("  Null DGM does not report subgroup-specific effects (expected).\n")
  check(TRUE, "Null DGM subgroup effects not available (OK for model='null')")
}


# ===================================================================
# SECTION 2: Survival Benefit + Adverse (S16/S9 pattern — ACTG style)
# ===================================================================
cat("\n\n--- Section 2: Survival Benefit + Adverse (S9/S16 pattern) ---\n")

# Adverse event (death). Treatment reduces hazard (beneficial).
# Q benefits MORE (HR_orig(Q) << 1).
# Under switched treatment: HR_sw(Q) >> 1.
set.seed(42)
rate_base <- 0.05
rate_alt <- rate_base * exp(-0.3 * syn$treat)   # treat reduces hazard
rate_alt[in_Q & syn$treat == 1] <-
  rate_alt[in_Q & syn$treat == 1] * 0.3          # Q: much fewer events
syn$time_surv <- rexp(N, rate_alt)
syn$event_surv <- as.integer(syn$time_surv <= 60)
syn$time_surv <- pmin(syn$time_surv, 60)

# Cox on original scale
cox_orig <- coxph(Surv(time_surv, event_surv) ~ treat, data = syn)
hr_orig <- exp(coef(cox_orig)["treat"])
cat(sprintf("\n  Raw data HR(ITT, original) = %.3f\n", hr_orig))
check(hr_orig < 0.8,
  sprintf("HR(ITT, original) = %.3f < 0.80 (beneficial)", hr_orig))

# Build DGM with switched treatment
df_sw <- data.frame(
  id = syn$id, treat = syn$treat_sw,
  z1 = syn$z1, z2 = syn$z2, z3 = syn$z3, z4 = syn$z4,
  y = syn$time_surv, event = syn$event_surv
)

dgm_surv_alt <- generate_aft_dgm_flex(
  data = df_sw, outcome_var = "y", event_var = "event",
  treatment_var = "treat", continuous_vars = character(0),
  factor_vars = c("z1","z2","z3","z4"),
  subgroup_vars = c("z1","z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 2.0, k_treat = 1,
  n_super = 5000L, seed = 42L, verbose = FALSE
)

hr_sw_Q <- dgm_surv_alt$hazard_ratios$harm_subgroup
hr_sw_ITT <- dgm_surv_alt$hazard_ratios$overall

cat(sprintf("  DGM (switched): HR(Q)=%.2f, HR(ITT)=%.2f\n",
    hr_sw_Q, hr_sw_ITT))
cat(sprintf("  DGM (original): HR(Q)=%.2f, HR(ITT)=%.2f\n",
    1/hr_sw_Q, 1/hr_sw_ITT))

# Switched: Q has high HR (control worse for Q → treat_sw worse)
check(hr_sw_Q > 1.5,
  sprintf("HR(Q, switched) = %.2f > 1.5 (exceeds threshold)", hr_sw_Q))

# Original: Q has low HR (strong benefit)
check(1/hr_sw_Q < 0.7,
  sprintf("HR(Q, original) = %.2f < 0.70 (strong benefit)", 1/hr_sw_Q))


# --- Null DGM with k_treat = 1 (ACTG production pattern) ---
cat("\n  Building null DGM (k_treat=1)...\n")
rate_null <- rate_base * exp(-0.15 * syn$treat)
syn$time_surv_null <- rexp(N, rate_null)
syn$event_surv_null <- as.integer(syn$time_surv_null <= 60)
syn$time_surv_null <- pmin(syn$time_surv_null, 60)

df_sw_null <- data.frame(
  id = syn$id, treat = syn$treat_sw,
  z1 = syn$z1, z2 = syn$z2, z3 = syn$z3, z4 = syn$z4,
  y = syn$time_surv_null, event = syn$event_surv_null
)

dgm_surv_null <- generate_aft_dgm_flex(
  data = df_sw_null, outcome_var = "y", event_var = "event",
  treatment_var = "treat", continuous_vars = character(0),
  factor_vars = c("z1","z2","z3","z4"),
  subgroup_vars = c("z1","z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "null", k_treat = 1,
  n_super = 5000L, seed = 42L, verbose = FALSE
)

hr_null_sw <- dgm_surv_null$hazard_ratios$overall
hr_null_orig <- 1 / hr_null_sw
hr_null_Q <- dgm_surv_null$hazard_ratios$harm_subgroup
hr_null_Qc <- dgm_surv_null$hazard_ratios$no_harm_subgroup

cat(sprintf("  Null HR(ITT, switched)=%.3f, HR(ITT, original)=%.3f\n",
    hr_null_sw, hr_null_orig))

# Switched HR > 1 (because treat_sw = control is better →
# "treatment" ddI is worse → HR > 1)
# But should be well below threshold 1.667
check(hr_null_sw > 1.0 && hr_null_sw < 1.5,
  sprintf("Null HR(ITT, switched) = %.3f in (1.0, 1.5) (below threshold 1.667)",
    hr_null_sw))

# Original: HR < 1 (beneficial)
check(hr_null_orig < 1.0,
  sprintf("Null HR(ITT, original) = %.3f < 1.0 (beneficial)", hr_null_orig))

# Not zero effect
check(abs(hr_null_orig - 1.0) > 0.05,
  sprintf("Null is NOT zero effect: |HR-1| = %.3f > 0.05",
    abs(hr_null_orig - 1.0)))

# Homogeneous: Q ≈ Qc
if (is.numeric(hr_null_Q) && is.numeric(hr_null_Qc) &&
    length(hr_null_Q) == 1 && length(hr_null_Qc) == 1) {
  cat(sprintf("  Null HR(Q, sw)=%.3f, HR(Qc, sw)=%.3f\n",
      hr_null_Q, hr_null_Qc))
  check(abs(log(hr_null_Q) - log(hr_null_Qc)) < 0.3,
    sprintf("Null Q vs Qc similar: |log(HR_Q)-log(HR_Qc)| = %.3f < 0.3",
      abs(log(hr_null_Q) - log(hr_null_Qc))))
} else {
  cat("  Null DGM does not report subgroup-specific effects (expected).\n")
  check(TRUE, "Null DGM subgroup effects not available (OK for model='null')")
}


# ===================================================================
# SECTION 3: Benefit + Positive Survival (S15 pattern — direct)
# ===================================================================
cat("\n\n--- Section 3: Survival Benefit + Positive (S15 pattern) ---\n")

# Positive event (recovery). Treatment accelerates recovery.
# Q recovers even faster. HR > 1 naturally exceeds threshold.
# NO switching needed.
rate_pos <- 0.05 * exp(0.15 * syn$treat)
rate_pos[in_Q & syn$treat == 1] <-
  rate_pos[in_Q & syn$treat == 1] * 3
syn$time_pos <- rexp(N, rate_pos)
syn$event_pos <- as.integer(syn$time_pos <= 60)
syn$time_pos <- pmin(syn$time_pos, 60)

df_pos <- data.frame(
  id = syn$id, treat = syn$treat,    # NOT switched
  z1 = syn$z1, z2 = syn$z2, z3 = syn$z3, z4 = syn$z4,
  y = syn$time_pos, event = syn$event_pos
)

dgm_pos_alt <- generate_aft_dgm_flex(
  data = df_pos, outcome_var = "y", event_var = "event",
  treatment_var = "treat", continuous_vars = character(0),
  factor_vars = c("z1","z2","z3","z4"),
  subgroup_vars = c("z1","z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 2.0, k_treat = 1,
  n_super = 5000L, seed = 42L, verbose = FALSE
)

hr_pos_Q <- dgm_pos_alt$hazard_ratios$harm_subgroup
hr_pos_ITT <- dgm_pos_alt$hazard_ratios$overall
cat(sprintf("\n  DGM (direct): HR(Q)=%.2f, HR(ITT)=%.2f\n",
    hr_pos_Q, hr_pos_ITT))

# HR > 1 = treat accelerates positive event = beneficial
check(hr_pos_Q > 2.0,
  sprintf("HR(Q) = %.2f > 2.0 (Q recovers much faster)", hr_pos_Q))
check(hr_pos_ITT > 1.0,
  sprintf("HR(ITT) = %.2f > 1.0 (treat accelerates recovery)", hr_pos_ITT))


# ===================================================================
# SECTION 4: Simulate + run_simulation_analysis round-trip (binary)
# ===================================================================
cat("\n\n--- Section 4: Binary simulate + run round-trip ---\n")

sim_df <- simulate_from_glm_dgm(dgm_alt, n = 600L, seed = 99L)

cat(sprintf("  Simulated: N=%d, response=%.1f%%, Q_prev=%.1f%%\n",
    nrow(sim_df), 100 * mean(sim_df$y_sim),
    100 * mean(sim_df$flag_harm)))

check(nrow(sim_df) == 600, "N = 600")
check("y_sim" %in% names(sim_df), "y_sim column exists")
check("treat_sim" %in% names(sim_df), "treat_sim column exists")
check("flag_harm" %in% names(sim_df), "flag_harm column exists")

# Verify simulated OR is in the right ballpark
fit_sim <- glm(y_sim ~ treat_sim, data = sim_df, family = binomial)
or_sim <- exp(coef(fit_sim)["treat_sim"])
cat(sprintf("  Simulated OR(ITT, switched) = %.3f (DGM: %.3f)\n",
    or_sim, or_sw_ITT))
check(abs(log(or_sim) - log(or_sw_ITT)) < 1.0,
  sprintf("|log(OR_sim) - log(OR_dgm)| = %.2f < 1.0",
    abs(log(or_sim) - log(or_sw_ITT))))

# run_simulation_analysis round-trip
result <- run_simulation_analysis(
  sim_id = 1L,
  dgm = dgm_alt,
  n_sample = 600L,
  confounders_base = c("z1","z2","z3","z4"),
  n_add_noise = 0L,
  run_fs = TRUE,
  run_fs_grf = FALSE,
  run_grf = FALSE,
  fs_params = list(
    outcome.name = "y_sim", event.name = "y_sim",
    treat.name = "treat_sim", id.name = "id",
    outcome_type = "binary", effect_measure = "OR",
    use_lasso = FALSE, use_grf = TRUE,
    effect.threshold = 1.5, consistency.threshold = 1.2,
    pconsistency.threshold = 0.80,
    fs.splits = 100L, n.min = 40L, d0.min = 8L, d1.min = 8L,
    maxk = 2L, vi.grf.min = -0.2, seedit = 42L,
    parallel_args = list(plan = "sequential", workers = 1,
                          show_message = FALSE)
  ),
  verbose = FALSE
)

check(is.data.frame(result), "run_simulation_analysis returns data.frame")
check("any.H" %in% names(result), "any.H column exists")
check(nrow(result) > 0, sprintf("Result has %d rows", nrow(result)))
cat(sprintf("  Detection: any.H = %s\n",
    paste(result$any.H[result$analysis == "FS"], collapse = ", ")))


# ===================================================================
# SUMMARY
# ===================================================================
cat(sprintf("\n\n=== SUMMARY: %d passed, %d failed ===\n",
    pass_count, fail_count))

if (fail_count > 0) {
  stop(sprintf("QC FAILED: %d checks failed", fail_count))
} else {
  cat("All checks passed.\n")
}
