#!/usr/bin/env Rscript
# ===================================================================
# qc_grf_scenario_validation.R
#
# Comprehensive pressure test for GRF subgroup identification across
# all 16 search × outcome × endpoint combinations.
#
# Calls grf.subg.harm.glm() and grf.subg.harm.survival() DIRECTLY
# (not through run_simulation_analysis) to validate sign conventions.
#
# For each scenario:
#   1. Single-dataset test on raw synthetic data
#   2. 10-sim detection check under H1
#   3. 10-sim FPR check under H0
#
# Expected runtime: ~3-5 minutes (sequential, no parallelism)
# ===================================================================

library(forestsearch)
library(survival)
library(data.table)

cat("=== GRF Scenario Validation Pressure Test ===\n")
cat("=== Calling GRF functions directly (no run_simulation_analysis) ===\n\n")

# -- Parameters ---------------------------------------------------------------
N_syn       <- 1200L
nsims_qc    <- 10L         # small-scale for pressure test
sim_n       <- 600L
dgm_n_super <- 5000L
sim_seed    <- 8316951L
confounders <- c("z1", "z2", "z3", "z4")

# -- GRF parameters (shared across all scenarios) ----------------------------
grf_base <- list(
  n.min        = 40L,
  dmin.grf     = 0,
  maxdepth     = 2L,
  RCT          = TRUE,
  sg.criterion = "mDiff",
  seedit       = 42L
)

# -- Results tracking ---------------------------------------------------------
pass_count <- 0L
fail_count <- 0L
results_table <- data.frame(
  Scenario = character(0), Endpoint = character(0),
  Search = character(0), Outcome = character(0),
  SingleDS = character(0), VI_top2 = character(0),
  Det_H1 = numeric(0), FPR_H0 = numeric(0),
  Status = character(0), stringsAsFactors = FALSE
)

check <- function(cond, label) {
  if (cond) {
    cat(sprintf("    PASS: %s\n", label))
    pass_count <<- pass_count + 1L
  } else {
    cat(sprintf("    FAIL: %s\n", label))
    fail_count <<- fail_count + 1L
  }
}


# -- Helper: run GRF on a single data frame -----------------------------------
run_grf_single <- function(data, treat_col, outcome_col,
                           outcome_type, adverse_outcome = FALSE,
                           event_col = NULL, label = "") {
  if (outcome_type == "survival") {
    res <- tryCatch(
      grf.subg.harm.survival(
        data = data, confounders.name = confounders,
        outcome.name = outcome_col, event.name = event_col,
        treat.name = treat_col, id.name = "id",
        n.min = grf_base$n.min, dmin.grf = grf_base$dmin.grf,
        maxdepth = grf_base$maxdepth, RCT = grf_base$RCT,
        sg.criterion = grf_base$sg.criterion,
        seedit = grf_base$seedit, details = FALSE
      ),
      error = function(e) {
        cat(sprintf("    ERROR [%s]: %s\n", label, e$message))
        NULL
      }
    )
  } else {
    res <- tryCatch(
      grf.subg.harm.glm(
        data = data, confounders.name = confounders,
        outcome.name = outcome_col, treat.name = treat_col,
        id.name = "id",
        outcome_type = outcome_type,
        adverse_outcome = adverse_outcome,
        n.min = grf_base$n.min, dmin.grf = grf_base$dmin.grf,
        maxdepth = grf_base$maxdepth, RCT = grf_base$RCT,
        sg.criterion = grf_base$sg.criterion,
        seedit = grf_base$seedit, details = FALSE,
        verbose = FALSE
      ),
      error = function(e) {
        cat(sprintf("    ERROR [%s]: %s\n", label, e$message))
        NULL
      }
    )
  }
  res
}


# -- Helper: check GRF result for Q overlap -----------------------------------
check_grf_result <- function(res, true_Q, label) {
  if (is.null(res)) {
    cat(sprintf("    %s: GRF returned NULL\n", label))
    return(list(detected = FALSE, sens = NA, vi_top2 = ""))
  }

  sg <- res$sg.harm.id
  detected <- !is.null(sg) && length(sg) > 0 && !all(is.na(sg))

  if (!detected) {
    cat(sprintf("    %s: No subgroup detected\n", label))
    return(list(detected = FALSE, sens = NA, vi_top2 = ""))
  }

  cat(sprintf("    %s: Subgroup = %s\n", label, paste(sg, collapse = " & ")))

  # Variable importance
  vi <- res$grf_varimp
  vi_top2 <- ""
  if (!is.null(vi) && length(vi) > 0) {
    vi_sorted <- sort(vi, decreasing = TRUE)
    vi_top2 <- paste(names(vi_sorted)[1:min(2, length(vi_sorted))],
                     collapse = ",")
  }

  # Classification
  if ("treat.recommend" %in% names(res$data)) {
    in_H <- res$data$treat.recommend == 0
    tp <- sum(in_H & true_Q)
    fp <- sum(in_H & !true_Q)
    fn <- sum(!in_H & true_Q)
    sens <- tp / max(tp + fn, 1)
    ppv  <- tp / max(tp + fp, 1)
    cat(sprintf("    %s: n(H)=%d, Sens=%.2f, PPV=%.2f\n",
        label, sum(in_H), sens, ppv))
    return(list(detected = TRUE, sens = sens, vi_top2 = vi_top2))
  }

  list(detected = TRUE, sens = NA, vi_top2 = vi_top2)
}


# -- Helper: multi-sim detection/FPR using DGMs -------------------------------
run_grf_sims <- function(dgm, nsims, outcome_type, adverse_outcome = FALSE,
                         label = "H1") {
  detections <- 0L
  is_surv <- outcome_type == "survival"

  for (i in seq_len(nsims)) {
    if (is_surv) {
      sim_df <- simulate_from_dgm(dgm, n = sim_n, seed = sim_seed + i,
                                   analysis_time = 60, cens_adjust = 0)
      outcome_col <- "y_sim"
      event_col   <- "event_sim"
      treat_col   <- "treat_sim"
    } else {
      sim_df <- simulate_from_glm_dgm(dgm, n = sim_n,
                                       seed = sim_seed + i)
      outcome_col <- "y_sim"
      event_col   <- NULL
      treat_col   <- "treat_sim"
    }

    res <- run_grf_single(sim_df, treat_col, outcome_col,
                          outcome_type, adverse_outcome, event_col,
                          label = sprintf("%s sim %d", label, i))

    if (!is.null(res)) {
      sg <- res$sg.harm.id
      if (!is.null(sg) && length(sg) > 0 && !all(is.na(sg))) {
        detections <- detections + 1L
      }
    }
  }

  rate <- 100 * detections / nsims
  cat(sprintf("  %s: %d/%d detected (%.0f%%)\n", label, detections, nsims, rate))
  rate
}


# ===================================================================
# SYNTHETIC DATA (same as FS scenario validation)
# ===================================================================
cat("\n--- Creating synthetic data ---\n")
set.seed(42)
syn <- data.frame(
  id = seq_len(N_syn), treat = rep(0:1, each = N_syn / 2),
  z1 = as.factor(sample(0:1, N_syn, TRUE)),
  z2 = as.factor(sample(0:1, N_syn, TRUE)),
  z3 = as.factor(sample(0:1, N_syn, TRUE)),
  z4 = as.factor(sample(0:1, N_syn, TRUE))
)
in_Q <- syn$z1 == 1 & syn$z2 == 1
syn$treat_sw <- 1L - syn$treat
syn$fu <- rep(1, N_syn)
cat(sprintf("  N=%d, Q prevalence=%.1f%%\n\n", N_syn, 100 * mean(in_Q)))

# -- Binary outcomes --
syn$y_adv_bin <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1 & in_Q, 0.65,
    ifelse(syn$treat == 1, 0.22, 0.35)))
syn$y_pos_bin <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1 & in_Q, 0.80,
    ifelse(syn$treat == 1, 0.50, 0.40)))
syn$y_s2_bin <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1 & in_Q, 0.15,
    ifelse(syn$treat == 1, 0.55, 0.40)))
syn$y_s4_bin <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1 & in_Q, 0.10,
    ifelse(syn$treat == 1, 0.33, 0.38)))

# -- Continuous outcomes --
syn$y_adv_cont <- 50 - 5 * syn$treat + rnorm(N_syn, sd = 15)
syn$y_adv_cont[in_Q & syn$treat == 1] <-
  syn$y_adv_cont[in_Q & syn$treat == 1] + 30
syn$y_pos_cont <- 200 + 8 * syn$treat + rnorm(N_syn, sd = 40)
syn$y_pos_cont[in_Q & syn$treat == 1] <-
  syn$y_pos_cont[in_Q & syn$treat == 1] + 50
# S6: positive outcome where Q is HARMED (treatment reduces Y for Q)
syn$y_s6_cont <- 200 + 8 * syn$treat + rnorm(N_syn, sd = 40)
syn$y_s6_cont[in_Q & syn$treat == 1] <-
  syn$y_s6_cont[in_Q & syn$treat == 1] - 50

# -- Count outcomes --
lam_adv <- exp(0.5 - 0.15 * syn$treat)
lam_adv[in_Q & syn$treat == 1] <- lam_adv[in_Q & syn$treat == 1] * 4
syn$y_adv_ct <- rpois(N_syn, lam_adv)

lam_pos <- exp(1.5 + 0.3 * syn$treat)
lam_pos[in_Q & syn$treat == 1] <- lam_pos[in_Q & syn$treat == 1] * 2.5
syn$y_pos_ct <- rpois(N_syn, lam_pos)

lam_s10 <- exp(1.5 + 0.3 * syn$treat)
lam_s10[in_Q & syn$treat == 1] <- lam_s10[in_Q & syn$treat == 1] * 0.3
syn$y_s10_ct <- rpois(N_syn, lam_s10)

# -- Survival outcomes --
rate_s13 <- 0.05 * exp(-0.3 * syn$treat)
rate_s13[in_Q & syn$treat == 1] <- rate_s13[in_Q & syn$treat == 1] * 4
syn$time_s13 <- pmin(rexp(N_syn, rate_s13), 60)
syn$event_s13 <- as.integer(syn$time_s13 < 60)

rate_s14 <- 0.05 * exp(0.3 * syn$treat)
rate_s14[in_Q & syn$treat == 1] <- rate_s14[in_Q & syn$treat == 1] * 0.10
syn$time_s14 <- pmin(rexp(N_syn, rate_s14), 60)
syn$event_s14 <- as.integer(syn$time_s14 < 60)

rate_s15 <- 0.05 * exp(0.15 * syn$treat)
rate_s15[in_Q & syn$treat == 1] <- rate_s15[in_Q & syn$treat == 1] * 3
syn$time_s15 <- pmin(rexp(N_syn, rate_s15), 60)
syn$event_s15 <- as.integer(syn$time_s15 < 60)

rate_s16 <- 0.05 * exp(-0.15 * syn$treat)
rate_s16[in_Q & syn$treat == 1] <- rate_s16[in_Q & syn$treat == 1] * 0.10
syn$time_s16 <- pmin(rexp(N_syn, rate_s16), 60)
syn$event_s16 <- as.integer(syn$time_s16 < 60)


# ===================================================================
# SINGLE-DATASET TESTS (all 16 scenarios)
# ===================================================================

run_single_scenario <- function(scenario, search, outcome, endpoint,
                                outcome_col, treat_col, outcome_type,
                                adverse_outcome = FALSE,
                                event_col = NULL) {
  cat(sprintf("\n=== %s: %s + %s %s ===\n", scenario, search, outcome, endpoint))

  res <- run_grf_single(syn, treat_col, outcome_col,
                        outcome_type, adverse_outcome, event_col,
                        label = scenario)
  info <- check_grf_result(res, in_Q, scenario)

  check(info$detected, sprintf("%s: GRF detected a subgroup", scenario))

  if (info$detected && !is.na(info$sens)) {
    check(info$sens > 0.2,
      sprintf("%s: Sensitivity = %.2f > 0.20", scenario, info$sens))
  }

  list(detected = info$detected, vi_top2 = info$vi_top2)
}

# -- S1-S4: Binary --
r1 <- run_single_scenario("S1","Harm","Adverse","Binary",
  "y_adv_bin","treat","binary", adverse_outcome=TRUE)
r2 <- run_single_scenario("S2","Harm","Positive","Binary",
  "y_s2_bin","treat","binary", adverse_outcome=FALSE)
r3 <- run_single_scenario("S3","Benefit","Positive","Binary",
  "y_pos_bin","treat_sw","binary", adverse_outcome=FALSE)
r4 <- run_single_scenario("S4","Benefit","Adverse","Binary",
  "y_s4_bin","treat_sw","binary", adverse_outcome=TRUE)

# -- S5-S8: Continuous --
r5 <- run_single_scenario("S5","Harm","Adverse","Continuous",
  "y_adv_cont","treat","continuous", adverse_outcome=TRUE)
r6 <- run_single_scenario("S6","Harm","Positive","Continuous",
  "y_s6_cont","treat","continuous", adverse_outcome=FALSE)
r7 <- run_single_scenario("S7","Benefit","Positive","Continuous",
  "y_pos_cont","treat_sw","continuous", adverse_outcome=FALSE)
r8 <- run_single_scenario("S8","Benefit","Adverse","Continuous",
  "y_adv_cont","treat_sw","continuous", adverse_outcome=TRUE)

# -- S9-S12: Count --
# NOTE: For count outcomes, grf.subg.harm.glm does NOT apply the
# adverse_outcome flip. So for harm+adverse count, we switch treatment
# to negate the CATE direction (same approach as benefit search).
r9 <- run_single_scenario("S9","Harm","Adverse","Count",
  "y_adv_ct","treat_sw","count", adverse_outcome=FALSE)
r10 <- run_single_scenario("S10","Harm","Positive","Count",
  "y_s10_ct","treat","count", adverse_outcome=FALSE)
r11 <- run_single_scenario("S11","Benefit","Positive","Count",
  "y_pos_ct","treat_sw","count", adverse_outcome=FALSE)
r12 <- run_single_scenario("S12","Benefit","Adverse","Count",
  "y_adv_ct","treat_sw","count", adverse_outcome=TRUE)

# -- S13-S16: Survival --
r13 <- run_single_scenario("S13","Harm","Adverse","Survival",
  "time_s13","treat","survival", event_col="event_s13")
r14 <- run_single_scenario("S14","Harm","Positive","Survival",
  "time_s14","treat_sw","survival", event_col="event_s14")
r15 <- run_single_scenario("S15","Benefit","Positive","Survival",
  "time_s15","treat_sw","survival", event_col="event_s15")
r16 <- run_single_scenario("S16","Benefit","Adverse","Survival",
  "time_s16","treat_sw","survival", event_col="event_s16")


# ===================================================================
# MULTI-SIM DETECTION / FPR (10 sims each, using DGMs)
# ===================================================================
cat("\n\n=== Multi-sim tests (10 alt, 10 null per scenario) ===\n")

# -- Null outcomes (homogeneous treatment, no HTE) --
syn$y_adv_bin_null <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1, 0.22, 0.35))
syn$y_pos_bin_null <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1, 0.55, 0.40))
syn$y_adv_cont_null <- 50 - 5 * syn$treat + rnorm(N_syn, sd = 15)
syn$y_pos_cont_null <- 200 + 8 * syn$treat + rnorm(N_syn, sd = 40)
syn$y_adv_ct_null <- rpois(N_syn, exp(0.5 - 0.15 * syn$treat))
syn$y_pos_ct_null <- rpois(N_syn, exp(1.5 + 0.3 * syn$treat))
rate_null_adv <- 0.05 * exp(-0.15 * syn$treat)
syn$time_null_adv <- pmin(rexp(N_syn, rate_null_adv), 60)
syn$event_null_adv <- as.integer(syn$time_null_adv < 60)
rate_null_pos <- 0.05 * exp(0.15 * syn$treat)
syn$time_null_pos <- pmin(rexp(N_syn, rate_null_pos), 60)
syn$event_null_pos <- as.integer(syn$time_null_pos < 60)

# Null outcomes for k_treat=0 scenarios
syn$y_s4_null <- rbinom(N_syn, 1, 0.40)
syn$y_s8_null <- 50 + rnorm(N_syn, sd = 15)
syn$y_s12_null <- rpois(N_syn, exp(0.5))

# Helper to build DGM and run multi-sim test
test_scenario_sims <- function(scenario, search, outcome, endpoint,
                               alt_y, null_y, treat_col,
                               outcome_type, effect_measure,
                               adverse_outcome = FALSE,
                               alt_event = NULL, null_event = NULL,
                               k_inter = 2.0, k_treat_null = 1) {
  cat(sprintf("\n--- %s: %s + %s %s ---\n", scenario, search, outcome, endpoint))

  is_surv <- outcome_type == "survival"

  # Build DGMs
  if (is_surv) {
    df_alt <- data.frame(id=syn$id, treat=syn[[treat_col]],
      z1=syn$z1, z2=syn$z2, z3=syn$z3, z4=syn$z4,
      y=syn[[alt_y]], event=syn[[alt_event]])
    df_null <- data.frame(id=syn$id, treat=syn[[treat_col]],
      z1=syn$z1, z2=syn$z2, z3=syn$z3, z4=syn$z4,
      y=syn[[null_y]], event=syn[[null_event]])

    dgm_alt <- generate_aft_dgm_flex(
      data=df_alt, outcome_var="y", event_var="event",
      treatment_var="treat", continuous_vars=character(0),
      factor_vars=confounders, subgroup_vars=c("z1","z2"),
      subgroup_cuts=list(z1=1L, z2=1L),
      model="alt", k_inter=k_inter, k_treat=1,
      n_super=dgm_n_super, seed=sim_seed, verbose=FALSE)
    dgm_null <- generate_aft_dgm_flex(
      data=df_null, outcome_var="y", event_var="event",
      treatment_var="treat", continuous_vars=character(0),
      factor_vars=confounders, subgroup_vars=c("z1","z2"),
      subgroup_cuts=list(z1=1L, z2=1L),
      model="null", k_treat=k_treat_null,
      n_super=dgm_n_super, seed=sim_seed, verbose=FALSE)
  } else {
    dgm_alt <- generate_glm_dgm(
      syn, confounders, alt_y, treat_col,
      outcome_type, effect_measure,
      subgroup_vars=c("z1","z2"), subgroup_cuts=list(z1=1L, z2=1L),
      model="alt", k_inter=k_inter,
      n_super=dgm_n_super, seed=sim_seed)
    dgm_null <- generate_glm_dgm(
      syn, confounders, null_y, treat_col,
      outcome_type, effect_measure,
      subgroup_vars=c("z1","z2"), subgroup_cuts=list(z1=1L, z2=1L),
      model="null", k_treat=k_treat_null,
      n_super=dgm_n_super, seed=sim_seed)
  }

  det <- run_grf_sims(dgm_alt, nsims_qc, outcome_type, adverse_outcome, "H1")
  fpr <- run_grf_sims(dgm_null, nsims_qc, outcome_type, adverse_outcome, "H0")

  status <- if (det >= 30 && fpr <= 50) "PASS" else "CHECK"
  cat(sprintf("  %s: Det=%.0f%%, FPR=%.0f%% -> %s\n",
      scenario, det, fpr, status))

  results_table[nrow(results_table) + 1, ] <<- data.frame(
    Scenario=scenario, Endpoint=endpoint, Search=search, Outcome=outcome,
    SingleDS=if (exists(paste0("r", gsub("S","",scenario)))) "Y" else "?",
    VI_top2="", Det_H1=det, FPR_H0=fpr, Status=status,
    stringsAsFactors=FALSE)
}

# Binary
test_scenario_sims("S1","Harm","Adverse","Binary",
  "y_adv_bin","y_adv_bin_null","treat","binary","OR",
  adverse_outcome=TRUE, k_inter=1.2)
test_scenario_sims("S2","Harm","Positive","Binary",
  "y_s2_bin","y_pos_bin_null","treat","binary","OR",
  adverse_outcome=FALSE, k_inter=2.0)
test_scenario_sims("S3","Benefit","Positive","Binary",
  "y_pos_bin","y_pos_bin_null","treat_sw","binary","OR",
  adverse_outcome=FALSE, k_inter=2.5)
test_scenario_sims("S4","Benefit","Adverse","Binary",
  "y_s4_bin","y_s4_null","treat_sw","binary","OR",
  adverse_outcome=TRUE, k_inter=1.2, k_treat_null=0)

# Continuous
test_scenario_sims("S5","Harm","Adverse","Continuous",
  "y_adv_cont","y_adv_cont_null","treat","continuous","MD",
  adverse_outcome=TRUE, k_inter=20)
test_scenario_sims("S6","Harm","Positive","Continuous",
  "y_s6_cont","y_pos_cont_null","treat","continuous","MD",
  adverse_outcome=FALSE, k_inter=20)
test_scenario_sims("S7","Benefit","Positive","Continuous",
  "y_pos_cont","y_pos_cont_null","treat_sw","continuous","MD",
  adverse_outcome=FALSE, k_inter=40)
test_scenario_sims("S8","Benefit","Adverse","Continuous",
  "y_adv_cont","y_s8_null","treat_sw","continuous","MD",
  adverse_outcome=TRUE, k_inter=20, k_treat_null=0)

# Count — NOTE: GRF doesn't apply adverse_outcome for count outcomes,
# so harm+adverse count (S9) uses switched treatment to negate CATE.
test_scenario_sims("S9","Harm","Adverse","Count",
  "y_adv_ct","y_adv_ct_null","treat_sw","count","IRR",
  adverse_outcome=FALSE, k_inter=0.7)
test_scenario_sims("S10","Harm","Positive","Count",
  "y_s10_ct","y_pos_ct_null","treat","count","IRR",
  adverse_outcome=FALSE, k_inter=0.7)
test_scenario_sims("S11","Benefit","Positive","Count",
  "y_pos_ct","y_pos_ct_null","treat_sw","count","IRR",
  adverse_outcome=FALSE, k_inter=2.0)
test_scenario_sims("S12","Benefit","Adverse","Count",
  "y_adv_ct","y_s12_null","treat_sw","count","IRR",
  adverse_outcome=TRUE, k_inter=2.0, k_treat_null=0)

# Survival
test_scenario_sims("S13","Harm","Adverse","Survival",
  "time_s13","time_null_adv","treat","survival","HR",
  alt_event="event_s13", null_event="event_null_adv",
  k_inter=2.0)
test_scenario_sims("S14","Harm","Positive","Survival",
  "time_s14","time_null_pos","treat_sw","survival","HR",
  alt_event="event_s14", null_event="event_null_pos",
  k_inter=2.0, k_treat_null=0)
test_scenario_sims("S15","Benefit","Positive","Survival",
  "time_s15","time_null_pos","treat_sw","survival","HR",
  alt_event="event_s15", null_event="event_null_pos",
  k_inter=2.0)
test_scenario_sims("S16","Benefit","Adverse","Survival",
  "time_s16","time_null_adv","treat_sw","survival","HR",
  alt_event="event_s16", null_event="event_null_adv",
  k_inter=2.0, k_treat_null=0)


# ===================================================================
# SUMMARY
# ===================================================================
cat("\n\n========================================\n")
cat("GRF SCENARIO VALIDATION SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("Single-dataset checks: %d passed, %d failed\n\n",
    pass_count, fail_count))

if (nrow(results_table) > 0) {
  cat("Multi-sim results:\n")
  print(results_table[, c("Scenario","Endpoint","Search","Outcome",
                           "Det_H1","FPR_H0","Status")],
        row.names = FALSE)
}

n_pass_sim <- sum(results_table$Status == "PASS")
n_total_sim <- nrow(results_table)
cat(sprintf("\nMulti-sim: %d/%d passed (Det>=30%%, FPR<=50%%)\n",
    n_pass_sim, n_total_sim))

if (fail_count > 0) {
  cat("\n*** SOME SINGLE-DATASET CHECKS FAILED ***\n")
  cat("Review failed scenarios before creating Quarto document.\n")
} else if (n_pass_sim < n_total_sim) {
  cat("\n*** SOME MULTI-SIM CHECKS NEED REVIEW ***\n")
  cat("May need parameter tuning (dmin.grf, k_inter) for specific scenarios.\n")
} else {
  cat("\nAll checks passed. GRF scenario validation document is safe to create.\n")
}
