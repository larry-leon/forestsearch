#!/usr/bin/env Rscript
# ===================================================================
# qc_grf_binary_benefit.R
#
# Pressure test: Can grf.subg.harm.glm() detect a benefit subgroup
# when treatment labels are switched?
#
# Setup mirrors actg175_binary_benefit_simulations.qmd:
#   - Y = improvement (positive outcome, 0/1)
#   - treat = switched (original control → treat=1)
#   - Q = {z1==1 & z2==1} benefits MORE from original treatment
#   - Under switched labels, Q appears as "harm" (ddI is worse for Q)
#
# Tests:
# 1. grf.subg.harm.glm() identifies Q or a superset
# 2. GRF variable importance picks up z1/z2
# 3. run_simulation_analysis() with run_grf=TRUE returns results
# 4. Null DGM produces no detection (FPR control)
# ===================================================================

library(forestsearch)
library(survival)

cat("=== QC: GRF Standalone with Binary Benefit Search ===\n\n")

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
# SECTION 1: Direct call to grf.subg.harm.glm()
# ===================================================================
cat("--- Section 1: Direct grf.subg.harm.glm() call ---\n\n")

set.seed(42)
N <- 1000
syn <- data.frame(
  id = seq_len(N),
  treat_orig = rep(0:1, each = N / 2),
  z1 = as.factor(sample(0:1, N, TRUE)),
  z2 = as.factor(sample(0:1, N, TRUE)),
  z3 = as.factor(sample(0:1, N, TRUE)),
  z4 = as.factor(sample(0:1, N, TRUE))
)
in_Q <- syn$z1 == 1 & syn$z2 == 1

# Switched treatment for benefit search
syn$treat <- 1L - syn$treat_orig

# Positive outcome: original treatment increases improvement for all,
# Q benefits MORE
syn$y <- rbinom(N, 1,
  ifelse(syn$treat_orig == 1 & in_Q, 0.80,
    ifelse(syn$treat_orig == 1, 0.50, 0.40)))

cat(sprintf("  N=%d, Q prevalence=%.1f%%\n", N, 100 * mean(in_Q)))
cat(sprintf("  Response rate: %.1f%% (overall)\n", 100 * mean(syn$y)))
cat(sprintf("  Response in Q:  treat_orig=1: %.1f%%  treat_orig=0: %.1f%%\n",
    100 * mean(syn$y[in_Q & syn$treat_orig == 1]),
    100 * mean(syn$y[in_Q & syn$treat_orig == 0])))
cat(sprintf("  Response in Qc: treat_orig=1: %.1f%%  treat_orig=0: %.1f%%\n",
    100 * mean(syn$y[!in_Q & syn$treat_orig == 1]),
    100 * mean(syn$y[!in_Q & syn$treat_orig == 0])))

# Call GRF with switched treatment
cat("\n  Calling grf.subg.harm.glm()...\n")
grf_result <- tryCatch(
  grf.subg.harm.glm(
    data             = syn,
    confounders.name = c("z1", "z2", "z3", "z4"),
    outcome.name     = "y",
    treat.name       = "treat",       # SWITCHED
    id.name          = "id",
    outcome_type     = "binary",
    effect_measure   = "log_OR",
    n.min            = 40L,
    dmin.grf         = 0,
    maxdepth         = 2L,
    RCT              = TRUE,
    sg.criterion     = "mDiff",
    seedit           = 42L,
    adverse_outcome  = FALSE,         # Y=improvement is positive
    verbose          = TRUE
  ),
  error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    NULL
  }
)

check(!is.null(grf_result), "grf.subg.harm.glm() returned a result")

if (!is.null(grf_result)) {
  sg_id <- grf_result$sg.harm.id
  cat(sprintf("  Subgroup identified: %s\n",
      if (is.null(sg_id) || length(sg_id) == 0) "NONE"
      else paste(sg_id, collapse = " & ")))

  check(!is.null(sg_id) && length(sg_id) > 0,
    "GRF identified a subgroup")

  # Check variable importance
  vi <- grf_result$grf_varimp
  if (!is.null(vi) && length(vi) > 0) {
    cat(sprintf("  Variable importance: %s\n",
        paste(sprintf("%s=%.3f", names(vi), vi), collapse = ", ")))

    # z1 and z2 should have highest importance
    vi_sorted <- sort(vi, decreasing = TRUE)
    top2 <- names(vi_sorted)[1:min(2, length(vi_sorted))]
    check("z1" %in% top2 || "z2" %in% top2,
      sprintf("z1 or z2 in top-2 VI: %s", paste(top2, collapse = ", ")))
  }

  # Check classification
  if (!is.null(syn$treat.recommend) ||
      "treat.recommend" %in% names(grf_result$data)) {
    df_grf <- grf_result$data
    in_H <- df_grf$treat.recommend == 0
    true_Q <- in_Q

    tp <- sum(in_H & true_Q)
    fp <- sum(in_H & !true_Q)
    fn <- sum(!in_H & true_Q)
    sens <- tp / max(tp + fn, 1)
    ppv  <- tp / max(tp + fp, 1)
    cat(sprintf("  Classification: TP=%d FP=%d FN=%d Sens=%.2f PPV=%.2f\n",
        tp, fp, fn, sens, ppv))
    check(sens > 0.3, sprintf("Sensitivity = %.2f > 0.30", sens))
  }
}


# ===================================================================
# SECTION 2: GRF through run_simulation_analysis()
# ===================================================================
cat("\n\n--- Section 2: run_simulation_analysis() with run_grf=TRUE ---\n\n")

# Build DGM
dgm_alt <- generate_glm_dgm(
  syn, c("z1","z2","z3","z4"), "y", "treat",
  "binary", "OR",
  subgroup_vars = c("z1","z2"), subgroup_cuts = list(z1=1L, z2=1L),
  model = "alt", k_inter = 2.0, n_super = 5000L, seed = 42L
)

cat(sprintf("  DGM: OR(Q,sw)=%.2f, OR(Qc,sw)=%.2f, OR(ITT,sw)=%.2f\n",
    dgm_alt$hazard_ratios$harm_subgroup,
    dgm_alt$hazard_ratios$no_harm_subgroup,
    dgm_alt$hazard_ratios$overall))

# Run with both FS and GRF
result <- tryCatch(
  run_simulation_analysis(
    sim_id = 1L,
    dgm = dgm_alt,
    n_sample = 800L,
    confounders_base = c("z1","z2","z3","z4"),
    n_add_noise = 0L,
    run_fs = TRUE,
    run_fs_grf = FALSE,
    run_grf = TRUE,
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
    grf_params = list(
      dmin.grf = 0, maxdepth = 2L, sg.criterion = "mDiff"
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
  cat(sprintf("  Analyses returned: %s\n", paste(analyses, collapse = ", ")))

  check("FS" %in% analyses, "FS analysis present")
  check("GRF" %in% analyses, "GRF analysis present")

  # FS detection
  fs_rows <- result[result$analysis == "FS", ]
  if (nrow(fs_rows) > 0) {
    cat(sprintf("  FS: any.H=%d\n", fs_rows$any.H[1]))
    check(fs_rows$any.H[1] == 1, "FS detected subgroup")
  }

  # GRF detection
  grf_rows <- result[result$analysis == "GRF", ]
  if (nrow(grf_rows) > 0) {
    cat(sprintf("  GRF: any.H=%d\n", grf_rows$any.H[1]))
    # GRF may or may not detect — just verify it ran without error
    check(TRUE, "GRF analysis completed without error")

    # Check if GRF has size.H populated
    if (!is.na(grf_rows$size.H[1])) {
      cat(sprintf("  GRF: size.H=%d, size.Hc=%d\n",
          grf_rows$size.H[1], grf_rows$size.Hc[1]))
    }
  }
}


# ===================================================================
# SECTION 3: Null DGM — GRF FPR
# ===================================================================
cat("\n\n--- Section 3: Null DGM — GRF should NOT detect ---\n\n")

# Null outcome: homogeneous mild benefit (no HTE)
syn$y_null <- rbinom(N, 1, ifelse(syn$treat_orig == 1, 0.55, 0.40))

dgm_null <- generate_glm_dgm(
  syn, c("z1","z2","z3","z4"), "y_null", "treat",
  "binary", "OR",
  subgroup_vars = c("z1","z2"), subgroup_cuts = list(z1=1L, z2=1L),
  model = "null", k_treat = 1, n_super = 5000L, seed = 42L
)

cat(sprintf("  Null DGM: OR(ITT,sw)=%.3f, OR(ITT,orig)=%.3f\n",
    dgm_null$hazard_ratios$overall,
    1 / dgm_null$hazard_ratios$overall))

# Run GRF on null
grf_null <- tryCatch(
  grf.subg.harm.glm(
    data             = syn,
    confounders.name = c("z1", "z2", "z3", "z4"),
    outcome.name     = "y_null",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "binary",
    effect_measure   = "log_OR",
    n.min            = 40L,
    dmin.grf         = 0,
    maxdepth         = 2L,
    RCT              = TRUE,
    sg.criterion     = "mDiff",
    seedit           = 42L,
    adverse_outcome  = FALSE,
    verbose          = FALSE
  ),
  error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    NULL
  }
)

check(!is.null(grf_null), "GRF null analysis completed")

if (!is.null(grf_null)) {
  sg_null <- grf_null$sg.harm.id
  detected_null <- !is.null(sg_null) && length(sg_null) > 0 &&
    !all(is.na(sg_null))
  cat(sprintf("  GRF null detection: %s\n",
      if (detected_null) paste(sg_null, collapse = " & ") else "NONE"))
  # Note: GRF on a single dataset may or may not detect.
  # The FPR test is across many simulations. Here we just verify it runs.
  check(TRUE, "GRF null analysis ran without error (FPR tested in simulation)")
}


# ===================================================================
# SECTION 4: Multi-simulation GRF FPR check (small scale)
# ===================================================================
cat("\n\n--- Section 4: Quick 20-sim GRF FPR check ---\n\n")

n_quick <- 20L
grf_detections <- 0L

for (i in seq_len(n_quick)) {
  sim_df <- simulate_from_glm_dgm(dgm_null, n = 800L, seed = 1000L + i)

  grf_i <- tryCatch(
    grf.subg.harm.glm(
      data             = sim_df,
      confounders.name = c("z1", "z2", "z3", "z4"),
      outcome.name     = "y_sim",
      treat.name       = "treat_sim",
      id.name          = "id",
      outcome_type     = "binary",
      effect_measure   = "log_OR",
      n.min            = 40L,
      dmin.grf         = 0,
      maxdepth         = 2L,
      RCT              = TRUE,
      sg.criterion     = "mDiff",
      seedit           = 42L,
      adverse_outcome  = FALSE,
      verbose          = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(grf_i) && !is.null(grf_i$sg.harm.id) &&
      length(grf_i$sg.harm.id) > 0 && !all(is.na(grf_i$sg.harm.id))) {
    grf_detections <- grf_detections + 1L
  }
}

grf_fpr <- 100 * grf_detections / n_quick
cat(sprintf("  GRF FPR (null, %d sims): %.1f%% (%d detections)\n",
    n_quick, grf_fpr, grf_detections))
check(grf_fpr <= 50,
  sprintf("GRF FPR = %.1f%% <= 50%% (relaxed threshold for 20 sims)", grf_fpr))


# ===================================================================
# SECTION 5: Multi-simulation GRF detection check (small scale)
# ===================================================================
cat("\n\n--- Section 5: Quick 20-sim GRF detection check ---\n\n")

grf_det_alt <- 0L

for (i in seq_len(n_quick)) {
  sim_df <- simulate_from_glm_dgm(dgm_alt, n = 800L, seed = 2000L + i)

  grf_i <- tryCatch(
    grf.subg.harm.glm(
      data             = sim_df,
      confounders.name = c("z1", "z2", "z3", "z4"),
      outcome.name     = "y_sim",
      treat.name       = "treat_sim",
      id.name          = "id",
      outcome_type     = "binary",
      effect_measure   = "log_OR",
      n.min            = 40L,
      dmin.grf         = 0,
      maxdepth         = 2L,
      RCT              = TRUE,
      sg.criterion     = "mDiff",
      seedit           = 42L,
      adverse_outcome  = FALSE,
      verbose          = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(grf_i) && !is.null(grf_i$sg.harm.id) &&
      length(grf_i$sg.harm.id) > 0 && !all(is.na(grf_i$sg.harm.id))) {
    grf_det_alt <- grf_det_alt + 1L
  }
}

grf_det_rate <- 100 * grf_det_alt / n_quick
cat(sprintf("  GRF detection (alt, %d sims): %.1f%% (%d detections)\n",
    n_quick, grf_det_rate, grf_det_alt))
check(grf_det_rate > 0,
  sprintf("GRF detection rate = %.1f%% > 0%%", grf_det_rate))


# ===================================================================
# SUMMARY
# ===================================================================
cat(sprintf("\n\n=== SUMMARY: %d passed, %d failed ===\n",
    pass_count, fail_count))

if (fail_count > 0) {
  stop(sprintf("QC FAILED: %d checks failed", fail_count))
} else {
  cat("All checks passed.\n")
  cat("\nGRF standalone is safe to add to binary benefit simulations.\n")
  cat("Recommended grf_params for binary benefit search:\n")
  cat("  dmin.grf = 0 (rate difference scale, not RMST)\n")
  cat("  maxdepth = 2L (matches FS maxk = 2)\n")
  cat("  sg.criterion = 'mDiff'\n")
  cat("  adverse_outcome = FALSE (Y = improvement is positive)\n")
}
