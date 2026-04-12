#!/usr/bin/env Rscript
# ===================================================================
# qc_grf_dmin_calibration_binary.R
#
# Systematic evaluation of dmin.grf for binary endpoints.
#
# PHASE 1: Identify correct GRF configuration for each binary scenario
#   Test all combinations of (treat_switching, adverse_outcome)
#   on single datasets to find which detect Q correctly.
#
# PHASE 2: Sweep dmin.grf for working configurations
#   For each working config, evaluate detection and FPR across
#   a grid of dmin.grf values using DGM-based simulations.
#
# PHASE 3: Vary effect size
#   Repeat PHASE 2 at weak/moderate/strong effect sizes to
#   map the detection-FPR tradeoff.
#
# Expected runtime: ~15-30 minutes (sequential)
# ===================================================================

library(forestsearch)
library(survival)
library(data.table)

cat("==========================================================\n")
cat("  GRF Binary Endpoint: dmin.grf Calibration Evaluation\n")
cat("==========================================================\n\n")

# -- Parameters ---------------------------------------------------------------
N_syn       <- 1200L
nsims       <- 20L
sim_n       <- 800L
dgm_n_super <- 5000L
sim_seed    <- 8316951L
confounders <- c("z1", "z2", "z3", "z4")

# dmin.grf grid (on the rate difference scale for binary outcomes)
# Clinically: 0.05 = 5pp difference, 0.10 = 10pp, 0.20 = 20pp
dmin_grid <- c(0, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30)

# Effect size multipliers for k_inter
# These control OR(Q) in the DGM: larger k_inter = stronger effect
k_inter_levels <- c(weak = 0.8, moderate = 1.5, strong = 2.5)


# ===================================================================
# PHASE 1: Configuration Discovery
# ===================================================================
cat("===========================================================\n")
cat("  PHASE 1: Identify correct GRF configuration per scenario\n")
cat("===========================================================\n\n")

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
cat(sprintf("  N=%d, Q prevalence=%.1f%%\n\n", N_syn, 100 * mean(in_Q)))

# Create four binary outcome scenarios
# (A) Harm + Adverse: treatment INCREASES bad outcome for Q
syn$y_harm_adv <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1 & in_Q, 0.65,
    ifelse(syn$treat == 1, 0.22, 0.35)))

# (B) Harm + Positive: treatment REDUCES good outcome for Q
syn$y_harm_pos <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1 & in_Q, 0.15,
    ifelse(syn$treat == 1, 0.55, 0.40)))

# (C) Benefit + Positive: treatment INCREASES good outcome MORE for Q
syn$y_ben_pos <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1 & in_Q, 0.80,
    ifelse(syn$treat == 1, 0.50, 0.40)))

# (D) Benefit + Adverse: treatment REDUCES bad outcome MORE for Q
syn$y_ben_adv <- rbinom(N_syn, 1,
  ifelse(syn$treat == 1 & in_Q, 0.10,
    ifelse(syn$treat == 1, 0.33, 0.38)))

# Test all 4 combinations of (switching, adverse_outcome) for each scenario
test_config <- function(scenario, y_col, treat_col, adverse_outcome) {
  res <- tryCatch(
    grf.subg.harm.glm(
      data = syn, confounders.name = confounders,
      outcome.name = y_col, treat.name = treat_col, id.name = "id",
      outcome_type = "binary", adverse_outcome = adverse_outcome,
      n.min = 40L, dmin.grf = 0, maxdepth = 2L, RCT = TRUE,
      sg.criterion = "mDiff", seedit = 42L, verbose = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(res) || is.null(res$sg.harm.id) ||
      length(res$sg.harm.id) == 0 || all(is.na(res$sg.harm.id))) {
    return(data.frame(Scenario = scenario, treat = treat_col,
      adverse_outcome = adverse_outcome, detected = FALSE,
      sens = NA, ppv = NA, n_H = NA, subgroup = "",
      stringsAsFactors = FALSE))
  }

  in_H <- res$data$treat.recommend == 0
  tp <- sum(in_H & in_Q); fp <- sum(in_H & !in_Q)
  fn <- sum(!in_H & in_Q)
  sens <- tp / max(tp + fn, 1)
  ppv  <- tp / max(tp + fp, 1)

  data.frame(Scenario = scenario, treat = treat_col,
    adverse_outcome = adverse_outcome, detected = TRUE,
    sens = round(sens, 2), ppv = round(ppv, 2), n_H = sum(in_H),
    subgroup = paste(res$sg.harm.id, collapse = " & "),
    stringsAsFactors = FALSE)
}

# Test all configurations
configs <- list()
for (scen in list(
  list(name = "Harm+Adverse", y = "y_harm_adv"),
  list(name = "Harm+Positive", y = "y_harm_pos"),
  list(name = "Benefit+Positive", y = "y_ben_pos"),
  list(name = "Benefit+Adverse", y = "y_ben_adv")
)) {
  for (tr in c("treat", "treat_sw")) {
    for (ao in c(FALSE, TRUE)) {
      configs[[length(configs) + 1]] <- test_config(
        scen$name, scen$y, tr, ao)
    }
  }
}

config_df <- do.call(rbind, configs)

# Display results
cat("Configuration discovery results:\n")
cat("(Sens > 0.50 = correct subgroup found)\n\n")
print(config_df[, c("Scenario", "treat", "adverse_outcome",
                     "detected", "sens", "ppv", "subgroup")],
      row.names = FALSE)

# Identify best config per scenario
cat("\n\nBest configuration per scenario (highest sensitivity):\n")
best_configs <- list()
for (scen in unique(config_df$Scenario)) {
  sub <- config_df[config_df$Scenario == scen & config_df$detected, ]
  if (nrow(sub) == 0) {
    cat(sprintf("  %s: NO WORKING CONFIG FOUND\n", scen))
    next
  }
  best <- sub[which.max(sub$sens), ]
  best_configs[[scen]] <- best
  cat(sprintf("  %s: treat=%s, adverse_outcome=%s -> Sens=%.2f, PPV=%.2f, SG=%s\n",
      scen, best$treat, best$adverse_outcome,
      best$sens, best$ppv, best$subgroup))
}


# ===================================================================
# PHASE 2: dmin.grf Sweep for Working Configurations
# ===================================================================
cat("\n\n===========================================================\n")
cat(sprintf("  PHASE 2: dmin.grf sweep (%d sims per point, N=%d)\n",
    nsims, sim_n))
cat("===========================================================\n\n")

# Helper: run GRF sims at a specific dmin.grf value
run_dmin_eval <- function(dgm, nsims, dmin_val, outcome_type,
                          adverse_outcome = FALSE) {
  detections <- 0L
  for (i in seq_len(nsims)) {
    sim_df <- simulate_from_glm_dgm(dgm, n = sim_n,
                                     seed = sim_seed + i * 1000L)
    res <- tryCatch(
      grf.subg.harm.glm(
        data = sim_df, confounders.name = confounders,
        outcome.name = "y_sim", treat.name = "treat_sim",
        id.name = "id", outcome_type = "binary",
        adverse_outcome = adverse_outcome,
        n.min = 40L, dmin.grf = dmin_val, maxdepth = 2L,
        RCT = TRUE, sg.criterion = "mDiff", seedit = 42L,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(res) && !is.null(res$sg.harm.id) &&
        length(res$sg.harm.id) > 0 && !all(is.na(res$sg.harm.id))) {
      detections <- detections + 1L
    }
  }
  100 * detections / nsims
}

# Build DGMs for each working scenario at each effect size
sweep_results <- data.frame(
  Scenario = character(0), Effect_Size = character(0),
  k_inter = numeric(0), OR_Q = numeric(0),
  dmin = numeric(0), Det_H1 = numeric(0), FPR_H0 = numeric(0),
  stringsAsFactors = FALSE
)

# Define scenario DGM specs using best configs
scenario_specs <- list(
  list(
    name = "Harm+Adverse",
    y_alt = "y_harm_adv",
    y_null = NULL,  # will create
    treat_col = NULL,  # from best_config
    adverse_outcome = NULL,
    null_rates = c(treat = 0.22, ctrl = 0.35)  # homogeneous beneficial
  ),
  list(
    name = "Harm+Positive",
    y_alt = "y_harm_pos",
    y_null = NULL,
    treat_col = NULL,
    adverse_outcome = NULL,
    null_rates = c(treat = 0.55, ctrl = 0.40)
  ),
  list(
    name = "Benefit+Positive",
    y_alt = "y_ben_pos",
    y_null = NULL,
    treat_col = NULL,
    adverse_outcome = NULL,
    null_rates = c(treat = 0.55, ctrl = 0.40)
  ),
  list(
    name = "Benefit+Adverse",
    y_alt = "y_ben_adv",
    y_null = NULL,
    treat_col = NULL,
    adverse_outcome = NULL,
    null_rates = c(treat = 0.22, ctrl = 0.35)
  )
)

# Fill in best configs
for (spec in scenario_specs) {
  if (spec$name %in% names(best_configs)) {
    bc <- best_configs[[spec$name]]
    spec$treat_col <- bc$treat
    spec$adverse_outcome <- bc$adverse_outcome
  }
}

cat(sprintf("dmin.grf grid: %s\n", paste(dmin_grid, collapse = ", ")))
cat(sprintf("Effect sizes: %s\n",
    paste(sprintf("%s (k=%.1f)", names(k_inter_levels), k_inter_levels),
          collapse = ", ")))
cat(sprintf("Sims per point: %d alt + %d null\n\n", nsims, nsims))

for (spec in scenario_specs) {
  scen_name <- spec$name
  if (!scen_name %in% names(best_configs)) {
    cat(sprintf("--- Skipping %s (no working config) ---\n\n", scen_name))
    next
  }

  bc <- best_configs[[scen_name]]
  treat_col <- bc$treat
  adverse_out <- bc$adverse_outcome

  cat(sprintf("--- %s (treat=%s, adverse_outcome=%s) ---\n",
      scen_name, treat_col, adverse_out))

  # Create null outcome
  null_rates <- spec$null_rates
  y_null_col <- paste0("y_null_", gsub("\\+", "", scen_name))
  syn[[y_null_col]] <- rbinom(N_syn, 1,
    ifelse(syn[[treat_col]] == 1, null_rates["treat"], null_rates["ctrl"]))

  for (es_name in names(k_inter_levels)) {
    k_val <- k_inter_levels[es_name]

    # Build alt DGM
    dgm_alt <- tryCatch(
      generate_glm_dgm(
        syn, confounders, spec$y_alt, treat_col,
        "binary", "OR",
        subgroup_vars = c("z1", "z2"),
        subgroup_cuts = list(z1 = 1L, z2 = 1L),
        model = "alt", k_inter = k_val,
        n_super = dgm_n_super, seed = sim_seed
      ),
      error = function(e) { cat("  DGM error:", e$message, "\n"); NULL }
    )

    # Build null DGM
    dgm_null <- tryCatch(
      generate_glm_dgm(
        syn, confounders, y_null_col, treat_col,
        "binary", "OR",
        subgroup_vars = c("z1", "z2"),
        subgroup_cuts = list(z1 = 1L, z2 = 1L),
        model = "null", k_treat = 1,
        n_super = dgm_n_super, seed = sim_seed
      ),
      error = function(e) { cat("  DGM error:", e$message, "\n"); NULL }
    )

    if (is.null(dgm_alt) || is.null(dgm_null)) next

    or_Q <- dgm_alt$hazard_ratios$harm_subgroup
    cat(sprintf("  %s (k=%.1f): OR(Q)=%.2f\n", es_name, k_val, or_Q))

    for (dmin_val in dmin_grid) {
      det <- run_dmin_eval(dgm_alt, nsims, dmin_val, "binary", adverse_out)
      fpr <- run_dmin_eval(dgm_null, nsims, dmin_val, "binary", adverse_out)

      sweep_results[nrow(sweep_results) + 1, ] <- data.frame(
        Scenario = scen_name, Effect_Size = es_name,
        k_inter = k_val, OR_Q = round(or_Q, 2),
        dmin = dmin_val, Det_H1 = det, FPR_H0 = fpr,
        stringsAsFactors = FALSE
      )

      cat(sprintf("    dmin=%.2f: Det=%.0f%%, FPR=%.0f%%%s\n",
          dmin_val, det, fpr,
          if (det >= 50 && fpr <= 15) " ***" else ""))
    }
    cat("\n")
  }
  cat("\n")
}


# ===================================================================
# PHASE 3: Summary and Recommendations
# ===================================================================
cat("\n===========================================================\n")
cat("  RESULTS SUMMARY\n")
cat("===========================================================\n\n")

if (nrow(sweep_results) > 0) {
  # Print full results table
  cat("Full sweep results:\n")
  print(sweep_results, row.names = FALSE)

  # Find optimal dmin per scenario × effect size
  cat("\n\nOptimal dmin.grf per scenario (Det >= 50% AND FPR <= 15%):\n")
  for (scen in unique(sweep_results$Scenario)) {
    for (es in unique(sweep_results$Effect_Size)) {
      sub <- sweep_results[sweep_results$Scenario == scen &
                            sweep_results$Effect_Size == es, ]
      good <- sub[sub$Det_H1 >= 50 & sub$FPR_H0 <= 15, ]
      if (nrow(good) > 0) {
        # Pick smallest dmin that achieves the target
        best <- good[which.min(good$dmin), ]
        cat(sprintf("  %s (%s, OR=%.2f): dmin=%.2f (Det=%.0f%%, FPR=%.0f%%)\n",
            scen, es, best$OR_Q, best$dmin, best$Det_H1, best$FPR_H0))
      } else {
        cat(sprintf("  %s (%s): No dmin achieves Det>=50%% AND FPR<=15%%\n",
            scen, es))
        # Show best tradeoff
        if (nrow(sub) > 0) {
          # Score: maximize detection, penalize FPR
          sub$score <- sub$Det_H1 - 2 * sub$FPR_H0
          best_trade <- sub[which.max(sub$score), ]
          cat(sprintf("    Best tradeoff: dmin=%.2f (Det=%.0f%%, FPR=%.0f%%)\n",
              best_trade$dmin, best_trade$Det_H1, best_trade$FPR_H0))
        }
      }
    }
  }

  # Summary statistics
  cat("\n\nFPR by dmin.grf (averaged across scenarios and effect sizes):\n")
  for (d in dmin_grid) {
    sub <- sweep_results[sweep_results$dmin == d, ]
    if (nrow(sub) > 0) {
      cat(sprintf("  dmin=%.2f: mean FPR=%.0f%%, mean Det=%.0f%%\n",
          d, mean(sub$FPR_H0), mean(sub$Det_H1)))
    }
  }
}

cat("\n\nPhase 1 best configs (for reference):\n")
for (scen in names(best_configs)) {
  bc <- best_configs[[scen]]
  cat(sprintf("  %s: treat=%s, adverse_outcome=%s, Sens=%.2f\n",
      scen, bc$treat, bc$adverse_outcome, bc$sens))
}

cat("\n=== Done ===\n")
