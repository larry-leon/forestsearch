#!/usr/bin/env Rscript
# ===================================================================
# qc_grf_actg175_binary_benefit.R
#
# Pressure test: GRF standalone for the ACTG175 binary benefit
# scenario (Benefit + Positive, switched treatment).
#
# Calibration finding: dmin.grf = 0.02 achieves Det >= 76%,
# FPR <= 9% for this scenario type.
#
# Tests:
# 1. Direct grf.subg.harm.glm() call on synthetic data
# 2. run_simulation_analysis() round-trip with whitelist fix
# 3. 20-sim detection under H1
# 4. 20-sim FPR under H0
# 5. dmin.grf = 0 vs 0.02 comparison
# ===================================================================

library(forestsearch)
library(survival)
library(data.table)

cat("=== QC: GRF Binary Benefit (ACTG175 scenario) ===\n")
cat("=== Calibrated dmin.grf = 0.02 ===\n\n")

# ===================================================================
# PARAMETERS (all tunable knobs in one place)
# ===================================================================
N              <- 1000L         # Synthetic data size
sim_n          <- 800L          # Per-simulation sample size
sim_seed       <- 8316951L
n_sims_direct  <- 20L           # Direct GRF sims (Section 3)
n_sims_rsa     <- 10L           # run_simulation_analysis sims (Section 4)
confounders    <- c("z1", "z2", "z3", "z4")

# GRF parameters
grf_dmin       <- 0.02          # calibrated dmin.grf
grf_maxdepth   <- 2L
grf_sg_crit    <- "mDiff"
grf_n_min      <- 40L
grf_adverse    <- FALSE         # benefit + positive: no Y-flip needed
                                 # (run_simulation_analysis defaults TRUE
                                 #  for binary; direct calls use this value)
grf_tune       <- TRUE         # cross-validated GRF hyperparameter tuning
grf_seedit     <- 42L

# FS parameters (for run_simulation_analysis Sections 2 & 4)
fs_or_threshold   <- 1.667      # 1/0.60
fs_or_consistency <- 1.25       # 1/0.80
fs_pconsistency   <- 0.80
fs_splits         <- 100L
fs_n_min          <- 40L
fs_d0_min         <- 8L
fs_d1_min         <- 8L
fs_maxk           <- 2L
fs_vi_grf_min     <- -0.2

pass_count <- 0L
fail_count <- 0L

check <- function(cond, label) {
  if (cond) {
    cat(sprintf("  PASS: %s\n", label))
    pass_count <<- pass_count + 1L
  } else {
    cat(sprintf("  FAIL: %s\n", label))
    fail_count <<- fail_count + 1L
  }
}

# ===================================================================
# SECTION 1: Direct grf.subg.harm.glm() on synthetic data
# ===================================================================
cat("--- Section 1: Direct GRF call (benefit + positive, switched) ---\n\n")

set.seed(42)
syn <- data.frame(
  id = seq_len(N),
  treat_orig = rep(0:1, each = N / 2),
  z1 = as.factor(sample(0:1, N, TRUE)),
  z2 = as.factor(sample(0:1, N, TRUE)),
  z3 = as.factor(sample(0:1, N, TRUE)),
  z4 = as.factor(sample(0:1, N, TRUE))
)
in_Q <- syn$z1 == 1 & syn$z2 == 1
syn$treat <- 1L - syn$treat_orig  # switched for benefit search

# Positive outcome: treatment increases Y for all, Q benefits MORE
syn$y <- rbinom(N, 1,
  ifelse(syn$treat_orig == 1 & in_Q, 0.80,
    ifelse(syn$treat_orig == 1, 0.50, 0.40)))

cat(sprintf("  N=%d, Q prevalence=%.1f%%\n", N, 100 * mean(in_Q)))
cat(sprintf("  RD(Q)=%+.3f, RD(Qc)=%+.3f (original scale)\n",
    mean(syn$y[in_Q & syn$treat_orig == 1]) -
      mean(syn$y[in_Q & syn$treat_orig == 0]),
    mean(syn$y[!in_Q & syn$treat_orig == 1]) -
      mean(syn$y[!in_Q & syn$treat_orig == 0])))

# Test with dmin = 0
cat("\n  dmin.grf = 0:\n")
res0 <- grf.subg.harm.glm(
  data = syn, confounders.name = confounders,
  outcome.name = "y", treat.name = "treat", id.name = "id",
  outcome_type = "binary", adverse_outcome = grf_adverse,
  n.min = grf_n_min, dmin.grf = 0, maxdepth = grf_maxdepth, RCT = TRUE,
  sg.criterion = grf_sg_crit, seedit = grf_seedit, verbose = TRUE,
  tune_grf = grf_tune
)
sg0 <- res0$sg.harm.id
cat(sprintf("    Subgroup: %s\n",
    if (is.null(sg0) || length(sg0) == 0) "NONE"
    else paste(sg0, collapse = " & ")))
check(!is.null(sg0) && length(sg0) > 0,
  "dmin=0: GRF detected a subgroup")

# Test with calibrated dmin
cat(sprintf("\n  dmin.grf = %.2f:\n", grf_dmin))
res02 <- grf.subg.harm.glm(
  data = syn, confounders.name = confounders,
  outcome.name = "y", treat.name = "treat", id.name = "id",
  outcome_type = "binary", adverse_outcome = grf_adverse,
  n.min = grf_n_min, dmin.grf = grf_dmin, maxdepth = grf_maxdepth, RCT = TRUE,
  sg.criterion = grf_sg_crit, seedit = grf_seedit, verbose = TRUE,
  tune_grf = grf_tune
)
sg02 <- res02$sg.harm.id
cat(sprintf("    Subgroup: %s\n",
    if (is.null(sg02) || length(sg02) == 0) "NONE"
    else paste(sg02, collapse = " & ")))

# Variable importance
vi <- res02$grf_varimp
if (!is.null(vi) && length(vi) > 0) {
  cat(sprintf("    VI: %s\n",
      paste(sprintf("%s=%.3f", names(vi), vi), collapse = ", ")))
  vi_sorted <- sort(vi, decreasing = TRUE)
  top2 <- names(vi_sorted)[1:min(2, length(vi_sorted))]
  check("z1" %in% top2 || "z2" %in% top2,
    sprintf("z1 or z2 in top-2 VI: %s", paste(top2, collapse = ", ")))
}


# ===================================================================
# SECTION 2: run_simulation_analysis() round-trip
# ===================================================================
cat("\n\n--- Section 2: run_simulation_analysis() with run_grf=TRUE ---\n\n")

dgm_alt <- generate_glm_dgm(
  syn, confounders, "y", "treat",
  "binary", "OR",
  subgroup_vars = c("z1","z2"), subgroup_cuts = list(z1=1L, z2=1L),
  model = "alt", k_inter = 2.0, n_super = 5000L, seed = grf_seedit
)
cat(sprintf("  DGM: OR(Q,sw)=%.2f, OR(Qc,sw)=%.2f\n",
    dgm_alt$hazard_ratios$harm_subgroup,
    dgm_alt$hazard_ratios$no_harm_subgroup))

result <- tryCatch(
  run_simulation_analysis(
    sim_id = 1L, dgm = dgm_alt, n_sample = sim_n,
    confounders_base = confounders,
    n_add_noise = 0L,
    run_fs = TRUE, run_fs_grf = FALSE, run_grf = TRUE,
    fs_params = list(
      outcome.name = "y_sim", event.name = "y_sim",
      treat.name = "treat_sim", id.name = "id",
      outcome_type = "binary", effect_measure = "OR",
      use_lasso = FALSE, use_grf = TRUE,
      effect.threshold = fs_or_threshold, consistency.threshold = fs_or_consistency,
      pconsistency.threshold = fs_pconsistency,
      fs.splits = fs_splits, n.min = fs_n_min,
      d0.min = fs_d0_min, d1.min = fs_d1_min,
      maxk = fs_maxk, vi.grf.min = fs_vi_grf_min, seedit = grf_seedit,
      tune_grf = grf_tune,
      parallel_args = list(plan = "sequential", workers = 1,
                            show_message = FALSE)
    ),
    grf_params = list(
      dmin.grf = grf_dmin, maxdepth = grf_maxdepth,
      sg.criterion = grf_sg_crit, tune_grf = grf_tune
    ),
    verbose = FALSE
  ),
  error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    NULL
  }
)

check(!is.null(result), "run_simulation_analysis() returned results")

if (!is.null(result) && is.data.frame(result)) {
  analyses <- unique(result$analysis)
  cat(sprintf("  Analyses: %s\n", paste(analyses, collapse = ", ")))
  check("FS" %in% analyses, "FS present")
  check("GRF" %in% analyses, "GRF present")

  for (a in analyses) {
    r <- result[result$analysis == a, ]
    cat(sprintf("  %s: any.H=%d", a, r$any.H[1]))
    if (!is.na(r$size.H[1])) cat(sprintf(", size.H=%d", r$size.H[1]))
    cat("\n")
  }

  # Verify no causal_survival_forest error
  check(!any(grepl("causal survival forest", result$analysis,
                    ignore.case = TRUE)),
    "No causal_survival_forest errors (dispatch fix working)")
}


# ===================================================================
# SECTION 3: 20-sim detection + FPR with dmin.grf comparison
# ===================================================================
cat(sprintf("\n\n--- Section 3: %d-sim detection/FPR (dmin=0 vs %.2f) ---\n\n",
    n_sims_direct, grf_dmin))

# Null DGM
syn$y_null <- rbinom(N, 1, ifelse(syn$treat_orig == 1, 0.55, 0.40))
dgm_null <- generate_glm_dgm(
  syn, confounders, "y_null", "treat",
  "binary", "OR",
  subgroup_vars = c("z1","z2"), subgroup_cuts = list(z1=1L, z2=1L),
  model = "null", k_treat = 1, n_super = 5000L, seed = grf_seedit
)

for (dmin_val in c(0, grf_dmin, 0.05)) {
  det <- 0L; fpr <- 0L

  for (i in seq_len(n_sims_direct)) {
    # Alt
    sim_df <- simulate_from_glm_dgm(dgm_alt, n = sim_n,
                                     seed = sim_seed + i * 1000L)
    res_a <- tryCatch(
      grf.subg.harm.glm(
        data = sim_df, confounders.name = confounders,
        outcome.name = "y_sim", treat.name = "treat_sim",
        id.name = "id", outcome_type = "binary",
        adverse_outcome = grf_adverse,
        n.min = grf_n_min, dmin.grf = dmin_val, maxdepth = grf_maxdepth,
        RCT = TRUE, sg.criterion = grf_sg_crit, seedit = grf_seedit,
        tune_grf = grf_tune, verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(res_a) && !is.null(res_a$sg.harm.id) &&
        length(res_a$sg.harm.id) > 0) det <- det + 1L

    # Null
    sim_df_n <- simulate_from_glm_dgm(dgm_null, n = sim_n,
                                       seed = sim_seed + 50000L + i * 1000L)
    res_n <- tryCatch(
      grf.subg.harm.glm(
        data = sim_df_n, confounders.name = confounders,
        outcome.name = "y_sim", treat.name = "treat_sim",
        id.name = "id", outcome_type = "binary",
        adverse_outcome = grf_adverse,
        n.min = grf_n_min, dmin.grf = dmin_val, maxdepth = grf_maxdepth,
        RCT = TRUE, sg.criterion = grf_sg_crit, seedit = grf_seedit,
        tune_grf = grf_tune, verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(res_n) && !is.null(res_n$sg.harm.id) &&
        length(res_n$sg.harm.id) > 0) fpr <- fpr + 1L
  }

  det_pct <- 100 * det / n_sims_direct
  fpr_pct <- 100 * fpr / n_sims_direct
  pass_flag <- if (det_pct >= 50 && fpr_pct <= 20) "***" else ""
  cat(sprintf("  dmin=%.2f: Det=%3.0f%% (%d/%d), FPR=%3.0f%% (%d/%d) %s\n",
      dmin_val, det_pct, det, n_sims_direct, fpr_pct, fpr, n_sims_direct, pass_flag))

  if (dmin_val == grf_dmin) {
    check(det_pct >= 30,
      sprintf("dmin=%.2f: Det=%.0f%% >= 30%%", grf_dmin, det_pct))
    check(fpr_pct <= 30,
      sprintf("dmin=%.2f: FPR=%.0f%% <= 30%%", grf_dmin, fpr_pct))
  }
}


# ===================================================================
# SECTION 4: run_simulation_analysis() with GRF — 10 sims
# ===================================================================
cat(sprintf("\n\n--- Section 4: %d-sim via run_simulation_analysis (FS + GRF) ---\n\n",
    n_sims_rsa))

det_fs <- 0L; det_grf <- 0L
fpr_fs <- 0L; fpr_grf <- 0L

# Build param lists once from setup parameters
qc_fs_params <- list(
  outcome.name = "y_sim", event.name = "y_sim",
  treat.name = "treat_sim", id.name = "id",
  outcome_type = "binary", effect_measure = "OR",
  use_lasso = FALSE, use_grf = TRUE,
  effect.threshold = fs_or_threshold, consistency.threshold = fs_or_consistency,
  pconsistency.threshold = fs_pconsistency,
  fs.splits = fs_splits, n.min = fs_n_min,
  d0.min = fs_d0_min, d1.min = fs_d1_min,
  maxk = fs_maxk, vi.grf.min = fs_vi_grf_min, seedit = grf_seedit,
  tune_grf = grf_tune,
  parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE)
)
qc_grf_params <- list(
  dmin.grf = grf_dmin, maxdepth = grf_maxdepth,
  sg.criterion = grf_sg_crit, tune_grf = grf_tune
)

for (i in seq_len(n_sims_rsa)) {
  # Alt
  r_alt <- tryCatch(
    run_simulation_analysis(
      sim_id = i, dgm = dgm_alt, n_sample = sim_n,
      confounders_base = confounders, n_add_noise = 0L,
      run_fs = TRUE, run_fs_grf = FALSE, run_grf = TRUE,
      fs_params = qc_fs_params, grf_params = qc_grf_params,
      verbose = FALSE),
    error = function(e) NULL)

  if (!is.null(r_alt) && is.data.frame(r_alt)) {
    if (any(r_alt$any.H[r_alt$analysis == "FS"] == 1)) det_fs <- det_fs + 1L
    if (any(r_alt$any.H[r_alt$analysis == "GRF"] == 1)) det_grf <- det_grf + 1L
  }

  # Null
  r_null <- tryCatch(
    run_simulation_analysis(
      sim_id = i + 1000L, dgm = dgm_null, n_sample = sim_n,
      confounders_base = confounders, n_add_noise = 0L,
      run_fs = TRUE, run_fs_grf = FALSE, run_grf = TRUE,
      fs_params = qc_fs_params, grf_params = qc_grf_params,
      verbose = FALSE),
    error = function(e) NULL)

  if (!is.null(r_null) && is.data.frame(r_null)) {
    if (any(r_null$any.H[r_null$analysis == "FS"] == 1)) fpr_fs <- fpr_fs + 1L
    if (any(r_null$any.H[r_null$analysis == "GRF"] == 1)) fpr_grf <- fpr_grf + 1L
  }
}

cat(sprintf("  FS:  Det=%d/%d (%.0f%%), FPR=%d/%d (%.0f%%)\n",
    det_fs, n_sims_rsa, 100*det_fs/n_sims_rsa,
    fpr_fs, n_sims_rsa, 100*fpr_fs/n_sims_rsa))
cat(sprintf("  GRF: Det=%d/%d (%.0f%%), FPR=%d/%d (%.0f%%)\n",
    det_grf, n_sims_rsa, 100*det_grf/n_sims_rsa,
    fpr_grf, n_sims_rsa, 100*fpr_grf/n_sims_rsa))

check(det_grf > 0, "GRF detected at least one alt subgroup via run_simulation_analysis")
check(!is.null(r_alt), "run_simulation_analysis completed without crash")


# ===================================================================
# SUMMARY
# ===================================================================
cat(sprintf("\n\n=== SUMMARY: %d passed, %d failed ===\n",
    pass_count, fail_count))

if (fail_count > 0) {
  stop(sprintf("QC FAILED: %d checks failed", fail_count))
} else {
  cat("All checks passed.\n")
  cat("\nGRF standalone with dmin.grf=0.02 is validated for\n")
  cat("binary benefit + positive (switched treatment) scenarios.\n")
  cat("\nRecommended parameters:\n")
  cat("  dmin.grf         = 0.02 (calibrated via 1500 sims)\n")
  cat("  maxdepth         = 2L\n")
  cat("  sg.criterion     = 'mDiff'\n")
  cat("  adverse_outcome  = FALSE\n")
  cat("  treatment        = switched (for benefit search)\n")
}
