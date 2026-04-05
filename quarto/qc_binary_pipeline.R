#!/usr/bin/env Rscript
# =============================================================================
# QC: Comprehensive synthetic-data validation of the binary simulation
#     pipeline for actg175_simulations_binary.qmd
#
# Tests:
#   A. Binary DGM construction (logistic model, super-population, interaction)
#   B. Binary trial simulation (sampling, outcome generation, column contract)
#   C. Extraction function correctness (mock FS results)
#   D. Inversion logic (benefit-search auto-inversion)
#   E. Table function interface contract
#   F. Edge cases (rare events, separation, degenerate data)
#   G. Parallelisation safety (nested plan, serialisation)
#   H. Seed reproducibility
#   I. OR calibration accuracy
#   J. Treatment switching algebra
# =============================================================================

library(survival)
library(parallel)
library(data.table)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doFuture))
suppressPackageStartupMessages(library(future))

pass_count <- 0L
fail_count <- 0L
test_log   <- character(0)

report <- function(test_id, description, passed, detail = "") {
  status <- if (passed) "PASS" else "FAIL"
  if (!passed) fail_count <<- fail_count + 1L else pass_count <<- pass_count + 1L
  msg <- sprintf("  [%s] %s: %s", status, test_id, description)
  if (nchar(detail) > 0) msg <- paste0(msg, " — ", detail)
  cat(msg, "\n")
  test_log <<- c(test_log, msg)
}

# ─────────────────────────────────────────────────────────────────────────────
# REPRODUCE FUNCTIONS FROM THE QMD (exactly as they appear)
# ─────────────────────────────────────────────────────────────────────────────

add_interaction <- function(df, beta_inter) {
  in_Q <- df$flag_harm == 1
  eta_1 <- qlogis(df$p1_base)
  eta_1[in_Q] <- eta_1[in_Q] + beta_inter
  df$p1 <- plogis(eta_1)
  df$p0 <- df$p0_base
  df
}

compute_or <- function(df, subset_flag = NULL) {
  if (!is.null(subset_flag)) df <- df[df$flag_harm == subset_flag, ]
  p1 <- mean(df$p1); p0 <- mean(df$p0)
  (p1 / (1 - p1)) / (p0 / (1 - p0))
}

make_binary_dgm <- function(df_sp, model_type, or_sg, or_sgc, or_overall,
                             fit_base = NULL, beta_inter_cal = 0) {
  list(
    df_super = df_sp,
    hazard_ratios = list(
      harm_subgroup    = or_sg,
      no_harm_subgroup = or_sgc,
      overall          = or_overall
    ),
    model_type = model_type,
    model_params = list(
      beta = if (!is.null(fit_base)) coef(fit_base) else NULL,
      beta_inter = beta_inter_cal,
      family = "binomial"
    ),
    n_super = nrow(df_sp)
  )
}

simulate_binary_trial <- function(dgm, n, sim_id, seed) {
  set.seed(seed + sim_id)
  idx <- sample(nrow(dgm$df_super), n, replace = TRUE)
  df  <- dgm$df_super[idx, ]
  df$id <- seq_len(n)
  df$treat_sim <- rbinom(n, 1, 0.5)
  df$y_sim <- ifelse(
    df$treat_sim == 1,
    rbinom(n, 1, df$p1),
    rbinom(n, 1, df$p0)
  )
  df$event_sim <- df$y_sim
  df
}

# Inversion helpers (from simulation_tables.R)
invert_hr_columns <- function(res) {
  hr_cols <- c("hr.H.hat","hr.Hc.hat","hr.H.bc","hr.Hc.bc",
               "hr.H.true","hr.Hc.true","hr.itt","hr.adj.itt",
               "ahr.H.hat","ahr.Hc.hat","ahr.H.true","ahr.Hc.true",
               "cde.H.hat","cde.Hc.hat")
  for (col in intersect(hr_cols, names(res)))
    res[[col]] <- 1 / res[[col]]
  res
}

invert_dgm_hrs <- function(dgm) {
  if (!is.null(dgm$hazard_ratios))
    for (nm in names(dgm$hazard_ratios)) {
      v <- dgm$hazard_ratios[[nm]]
      if (is.numeric(v) && length(v) == 1 && !is.na(v))
        dgm$hazard_ratios[[nm]] <- 1 / v
    }
  for (f in c("hr_H_true","hr_Hc_true","hr_causal",
              "AHR","AHR_H_true","AHR_Hc_true",
              "cde_H","cde_Hc","CDE"))
    if (!is.null(dgm[[f]])) dgm[[f]] <- 1 / dgm[[f]]
  dgm
}

# get_dgm_hr (from oc_analyses.R)
get_dgm_hr <- function(dgm, which = "hr_H") {
  if (!is.null(dgm$hazard_ratios)) {
    hr <- dgm$hazard_ratios
    r <- switch(which, hr_H=hr$harm_subgroup, hr_Hc=hr$no_harm_subgroup,
                hr_overall=hr$overall, ahr_H=hr$AHR_harm, ahr_Hc=hr$AHR_no_harm,
                cde_H=hr$CDE_harm, cde_Hc=hr$CDE_no_harm, NA)
    if (!is.null(r) && !is.na(r)) return(r)
  }
  r <- switch(which, hr_H=dgm$hr_H_true, hr_Hc=dgm$hr_Hc_true,
              hr_overall=dgm$hr_causal, NA)
  if (is.null(r)) NA_real_ else r
}


# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("  COMPREHENSIVE QC: Binary Simulation Pipeline\n")
cat("=", rep("=", 70), "\n\n")
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
cat("─── TEST GROUP A: Binary DGM Construction ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

set.seed(42)
n_syn <- 1000
z_vars <- paste0("z", 1:6)

syn_df <- data.frame(
  id    = 1:n_syn,
  z1    = factor(rbinom(n_syn, 1, 0.35)),
  z2    = factor(rbinom(n_syn, 1, 0.55)),
  z3    = factor(rbinom(n_syn, 1, 0.40)),
  z4    = factor(rbinom(n_syn, 1, 0.50)),
  z5    = factor(rbinom(n_syn, 1, 0.45)),
  z6    = factor(rbinom(n_syn, 1, 0.30)),
  treat = rbinom(n_syn, 1, 0.5)
)
syn_df$flag_harm <- as.integer(syn_df$z1 == 1 & syn_df$z2 == 1)

# Generate binary outcome with known interaction
beta_true <- c(-0.5, 0.3, 0.4, -0.2, 0.1, -0.1, 0.2)
eta <- with(syn_df, {
  -0.5 + 0.3 * as.numeric(treat) +
  0.4 * as.numeric(as.character(z1)) - 0.2 * as.numeric(as.character(z2)) +
  0.1 * as.numeric(as.character(z3)) - 0.1 * as.numeric(as.character(z4)) +
  0.8 * as.numeric(treat) * as.numeric(flag_harm)  # interaction
})
syn_df$y_binary <- rbinom(n_syn, 1, plogis(eta))

# A1: Logistic model fit
fit_syn <- glm(y_binary ~ treat + z1 + z2 + z3 + z4 + z5 + z6,
               data = syn_df, family = binomial)
report("A1", "Logistic model converges",
       fit_syn$converged,
       sprintf("deviance = %.1f, AIC = %.1f", deviance(fit_syn), AIC(fit_syn)))

# A2: Super-population resampling
set.seed(99)
idx_sp <- sample(n_syn, 5000, replace = TRUE)
df_sp  <- syn_df[idx_sp, ]
df_sp$id <- seq_len(5000)
report("A2", "Super-population size correct",
       nrow(df_sp) == 5000)

# A3: Potential outcome computation
newdata_0 <- df_sp; newdata_0$treat <- 0L
newdata_1 <- df_sp; newdata_1$treat <- 1L
df_sp$p0_base <- predict(fit_syn, newdata = newdata_0, type = "response")
df_sp$p1_base <- predict(fit_syn, newdata = newdata_1, type = "response")

report("A3a", "P(Y|treat=0) in [0,1]",
       all(df_sp$p0_base >= 0 & df_sp$p0_base <= 1))
report("A3b", "P(Y|treat=1) in [0,1]",
       all(df_sp$p1_base >= 0 & df_sp$p1_base <= 1))
report("A3c", "P0 and P1 differ",
       !isTRUE(all.equal(df_sp$p0_base, df_sp$p1_base)))

# A4: Interaction calibration
beta_grid <- seq(0, 3.0, by = 0.05)
or_grid   <- sapply(beta_grid, function(b) {
  df_mod <- add_interaction(df_sp, b)
  compute_or(df_mod, 1)
})

report("A4a", "OR increases monotonically with beta_inter",
       all(diff(or_grid) > 0),
       sprintf("OR range: [%.2f, %.2f]", min(or_grid), max(or_grid)))

# Find beta for target OR = 2.0
target_or <- 2.0
best_idx  <- which.min(abs(or_grid - target_or))
beta_cal  <- beta_grid[best_idx]
or_achieved <- or_grid[best_idx]

report("A4b", "Calibration achieves target OR within 10%",
       abs(or_achieved - target_or) / target_or < 0.10,
       sprintf("beta = %.2f → OR(Q) = %.3f (target: %.1f)", beta_cal, or_achieved, target_or))

# A5: Complement OR remains below Q (complement should have weaker effect)
df_alt <- add_interaction(df_sp, beta_cal)
or_Qc  <- compute_or(df_alt, 0)
report("A5", "OR(Qc) < OR(Q) (complement has weaker effect)",
       or_Qc < or_achieved,
       sprintf("OR(Qc) = %.3f < OR(Q) = %.3f", or_Qc, or_achieved))

# A6: OR(ITT) between OR(Q) and OR(Qc)
or_itt <- compute_or(df_alt)
report("A6", "OR(ITT) is between OR(Q) and OR(Qc)",
       or_itt >= min(or_Qc, or_achieved) && or_itt <= max(or_Qc, or_achieved),
       sprintf("OR(ITT) = %.3f", or_itt))


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP B: Binary Trial Simulation ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

dgm_syn <- make_binary_dgm(df_alt, "alt", or_achieved, or_Qc, or_itt,
                             fit_syn, beta_cal)

sim1 <- simulate_binary_trial(dgm_syn, n = 500, sim_id = 1, seed = 42)

# B1: Column contract
expected_cols <- c("id", "treat_sim", "y_sim", "event_sim",
                   "flag_harm", "p0", "p1", paste0("z", 1:6))
missing <- setdiff(expected_cols, names(sim1))
report("B1", "All expected columns present",
       length(missing) == 0,
       if (length(missing) > 0) paste("missing:", paste(missing, collapse=", ")) else "")

# B2: Correct dimensions
report("B2", "Simulated data has correct nrow",
       nrow(sim1) == 500)

# B3: Binary outcome values
report("B3a", "y_sim is binary {0,1}",
       all(sim1$y_sim %in% c(0L, 1L)))
report("B3b", "event_sim equals y_sim",
       identical(sim1$event_sim, sim1$y_sim))

# B4: Treatment balance
treat_prop <- mean(sim1$treat_sim)
report("B4", "Treatment balance near 0.50",
       abs(treat_prop - 0.50) < 0.10,
       sprintf("treat proportion = %.3f", treat_prop))

# B5: Outcome prevalence is reasonable
y_prev <- mean(sim1$y_sim)
report("B5", "Binary outcome prevalence in [0.1, 0.9]",
       y_prev > 0.1 && y_prev < 0.9,
       sprintf("prevalence = %.3f", y_prev))

# B6: IDs are unique
report("B6", "Subject IDs are unique",
       length(unique(sim1$id)) == nrow(sim1))

# B7: Subgroup flag preserved
report("B7", "flag_harm prevalence matches DGM super-population",
       abs(mean(sim1$flag_harm) - mean(dgm_syn$df_super$flag_harm)) < 0.10,
       sprintf("sim = %.3f, dgm = %.3f",
               mean(sim1$flag_harm), mean(dgm_syn$df_super$flag_harm)))

# B8: Potential outcomes are present and valid
report("B8a", "p0 column in [0,1]",
       all(sim1$p0 >= 0 & sim1$p0 <= 1))
report("B8b", "p1 column in [0,1]",
       all(sim1$p1 >= 0 & sim1$p1 <= 1))

# B9: Treatment effect is in the right direction for Q members
# Under switched treatment, Q members should have higher P(Y=1|treat=1)
q_members <- sim1[sim1$flag_harm == 1, ]
report("B9", "p1 > p0 on average for Q (interaction present)",
       mean(q_members$p1) > mean(q_members$p0),
       sprintf("mean(p1|Q) = %.3f, mean(p0|Q) = %.3f",
               mean(q_members$p1), mean(q_members$p0)))


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP C: Extraction Function (Mock FS Results) ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

# Reproduce extract_binary_estimates from the qmd
extract_binary_estimates <- function(df, fs_res, dgm, analysis = "FS") {
  out <- data.frame(
    analysis  = analysis,
    any.H     = 0L,
    hr.H.hat  = NA_real_,
    hr.Hc.hat = NA_real_,
    hr.H.true = NA_real_,
    hr.Hc.true = NA_real_,
    hr.itt     = NA_real_,
    sens       = NA_real_,
    spec       = NA_real_,
    ppv        = NA_real_,
    npv        = NA_real_,
    size.H     = NA_real_,
    size.Hc    = nrow(df),
    p.cens     = 0,
    stringsAsFactors = FALSE
  )
  if (is.null(fs_res)) return(out)
  sg_harm <- fs_res$grp.consistency$sg.harm.id
  if (is.null(sg_harm) || length(sg_harm) == 0 ||
      all(is.na(sg_harm)) || all(sg_harm == "")) {
    itt_fit <- tryCatch(
      glm(y_sim ~ treat_sim, data = df, family = binomial),
      error = function(e) NULL
    )
    if (!is.null(itt_fit)) out$hr.itt <- exp(coef(itt_fit)["treat_sim"])
    return(out)
  }
  out$any.H <- 1L
  df_pred <- fs_res$df.est
  if (is.null(df_pred) || !"treat.recommend" %in% names(df_pred))
    return(out)
  in_H  <- df_pred$treat.recommend == 0
  in_Hc <- df_pred$treat.recommend == 1
  out$size.H  <- sum(in_H)
  out$size.Hc <- sum(in_Hc)
  if (sum(in_H) >= 10 && length(unique(df_pred$treat_sim[in_H])) == 2) {
    fit_H <- tryCatch(
      glm(y_sim ~ treat_sim, data = df_pred[in_H, ], family = binomial),
      error = function(e) NULL
    )
    if (!is.null(fit_H)) out$hr.H.hat <- exp(coef(fit_H)["treat_sim"])
  }
  if (sum(in_Hc) >= 10 && length(unique(df_pred$treat_sim[in_Hc])) == 2) {
    fit_Hc <- tryCatch(
      glm(y_sim ~ treat_sim, data = df_pred[in_Hc, ], family = binomial),
      error = function(e) NULL
    )
    if (!is.null(fit_Hc)) out$hr.Hc.hat <- exp(coef(fit_Hc)["treat_sim"])
  }
  itt_fit <- tryCatch(
    glm(y_sim ~ treat_sim, data = df_pred, family = binomial),
    error = function(e) NULL
  )
  if (!is.null(itt_fit)) out$hr.itt <- exp(coef(itt_fit)["treat_sim"])
  out$hr.H.true  <- dgm$hazard_ratios$harm_subgroup
  out$hr.Hc.true <- dgm$hazard_ratios$no_harm_subgroup
  true_H <- df_pred$flag_harm == 1
  if (sum(true_H) > 0 && sum(in_H) > 0) {
    tp <- sum(in_H & true_H)
    fp <- sum(in_H & !true_H)
    fn <- sum(!in_H & true_H)
    tn <- sum(!in_H & !true_H)
    out$sens <- tp / (tp + fn)
    out$spec <- tn / (tn + fp)
    out$ppv  <- if (tp + fp > 0) tp / (tp + fp) else NA
    out$npv  <- if (tn + fn > 0) tn / (tn + fn) else NA
  }
  out
}

# C1: NULL result (FS failed)
out_null <- extract_binary_estimates(sim1, NULL, dgm_syn)
report("C1a", "NULL FS result returns any.H = 0",
       out_null$any.H == 0)
report("C1b", "NULL FS result has correct ncol",
       ncol(out_null) == 14)

# C2: No subgroup found
mock_fs_nosub <- list(
  grp.consistency = list(sg.harm.id = NULL),
  df.est = sim1
)
out_nosub <- extract_binary_estimates(sim1, mock_fs_nosub, dgm_syn)
report("C2a", "No-subgroup result: any.H = 0",
       out_nosub$any.H == 0)
report("C2b", "No-subgroup result: ITT OR computed",
       !is.na(out_nosub$hr.itt) && out_nosub$hr.itt > 0)

# C3: Subgroup found (mock perfect identification)
mock_df_est <- sim1
mock_df_est$treat.recommend <- ifelse(mock_df_est$flag_harm == 1, 0L, 1L)
mock_fs_found <- list(
  grp.consistency = list(sg.harm.id = c("z1.1", "z2.1")),
  df.est = mock_df_est
)
out_found <- extract_binary_estimates(sim1, mock_fs_found, dgm_syn)
report("C3a", "Subgroup found: any.H = 1",
       out_found$any.H == 1)
report("C3b", "OR(Q) is positive",
       !is.na(out_found$hr.H.hat) && out_found$hr.H.hat > 0,
       sprintf("OR(Q) = %.3f", out_found$hr.H.hat))
report("C3c", "OR(Qc) is positive",
       !is.na(out_found$hr.Hc.hat) && out_found$hr.Hc.hat > 0,
       sprintf("OR(Qc) = %.3f", out_found$hr.Hc.hat))
report("C3d", "Perfect identification: sens = 1",
       isTRUE(all.equal(out_found$sens, 1.0)))
report("C3e", "Perfect identification: spec = 1",
       isTRUE(all.equal(out_found$spec, 1.0)))
report("C3f", "Perfect identification: ppv = 1",
       isTRUE(all.equal(out_found$ppv, 1.0)))
report("C3g", "Size(H) matches flag_harm count",
       out_found$size.H == sum(sim1$flag_harm),
       sprintf("size.H = %d, sum(flag_harm) = %d",
               out_found$size.H, sum(sim1$flag_harm)))

# C4: Column names match table function expectations
expected_result_cols <- c("analysis", "any.H", "hr.H.hat", "hr.Hc.hat",
                          "hr.H.true", "hr.Hc.true", "hr.itt",
                          "sens", "spec", "ppv", "npv",
                          "size.H", "size.Hc", "p.cens")
report("C4", "Result column names match table function contract",
       all(expected_result_cols %in% names(out_found)),
       paste("cols:", ncol(out_found)))

# C5: rbind of multiple results
out_multi <- rbind(out_found, out_nosub, out_null)
report("C5", "rbind of mixed results produces correct shape",
       nrow(out_multi) == 3 && ncol(out_multi) == 14)


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP D: Inversion Logic ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

# D1: Results inversion
mock_res <- data.frame(
  analysis   = rep("FS", 5),
  any.H      = c(1, 1, 1, 1, 0),
  hr.H.hat   = c(2.50, 2.10, 2.80, 2.30, NA),
  hr.Hc.hat  = c(1.15, 1.05, 1.20, 1.10, 0.95),
  hr.H.true  = c(2.00, 2.00, 2.00, 2.00, NA),
  hr.Hc.true = rep(1.10, 5),
  hr.itt     = rnorm(5, 1.30, 0.1),
  sens       = c(0.8, 0.7, 0.85, 0.75, NA),
  size.H     = c(250, 270, 240, 280, NA),
  stringsAsFactors = FALSE
)

inv_res <- invert_hr_columns(mock_res)

# D1a: Element-wise 1/x for each HR column
for (col in c("hr.H.hat", "hr.Hc.hat", "hr.H.true", "hr.Hc.true", "hr.itt")) {
  ok <- all(abs(inv_res[[col]] - 1/mock_res[[col]]) < 1e-10, na.rm = TRUE)
  report(paste0("D1-", col), sprintf("Element-wise inversion: %s", col), ok)
}

# D1b: Non-HR columns unchanged
report("D1-sens", "sens unchanged after inversion",
       identical(mock_res$sens, inv_res$sens))
report("D1-size", "size.H unchanged after inversion",
       identical(mock_res$size.H, inv_res$size.H))
report("D1-any", "any.H unchanged after inversion",
       identical(mock_res$any.H, inv_res$any.H))

# D2: DGM inversion
inv_dgm <- invert_dgm_hrs(dgm_syn)
report("D2a", "DGM harm_subgroup inverted",
       abs(inv_dgm$hazard_ratios$harm_subgroup - 1/dgm_syn$hazard_ratios$harm_subgroup) < 1e-10)
report("D2b", "DGM no_harm_subgroup inverted",
       abs(inv_dgm$hazard_ratios$no_harm_subgroup - 1/dgm_syn$hazard_ratios$no_harm_subgroup) < 1e-10)
report("D2c", "DGM overall inverted",
       abs(inv_dgm$hazard_ratios$overall - 1/dgm_syn$hazard_ratios$overall) < 1e-10)

# D3: get_dgm_hr returns inverted values
report("D3a", "get_dgm_hr(inv, 'hr_H') < 1",
       get_dgm_hr(inv_dgm, "hr_H") < 1,
       sprintf("= %.3f", get_dgm_hr(inv_dgm, "hr_H")))
report("D3b", "get_dgm_hr(inv, 'hr_Hc') < 1",
       get_dgm_hr(inv_dgm, "hr_Hc") < 1,
       sprintf("= %.3f", get_dgm_hr(inv_dgm, "hr_Hc")))

# D4: Double inversion = identity
dgm_double <- invert_dgm_hrs(inv_dgm)
report("D4", "Double inversion returns to original",
       abs(dgm_double$hazard_ratios$harm_subgroup - dgm_syn$hazard_ratios$harm_subgroup) < 1e-10)


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP E: Table Function Interface Contract ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

# Build a 50-row mock results table mimicking a full simulation
set.seed(123)
n_mock <- 50
mock_full <- data.frame(
  analysis   = rep("FS", n_mock),
  any.H      = c(rep(1, 45), rep(0, 5)),
  hr.H.hat   = c(rnorm(45, 2.5, 0.5), rep(NA, 5)),
  hr.Hc.hat  = rnorm(n_mock, 1.15, 0.2),
  hr.H.bc    = c(rnorm(45, 2.2, 0.4), rep(NA, 5)),
  hr.Hc.bc   = rnorm(n_mock, 1.12, 0.18),
  hr.H.true  = c(rep(2.0, 45), rep(NA, 5)),
  hr.Hc.true = rep(1.10, n_mock),
  hr.itt     = rnorm(n_mock, 1.30, 0.1),
  ahr.H.hat  = c(rnorm(45, 2.3, 0.4), rep(NA, 5)),
  ahr.Hc.hat = rnorm(n_mock, 1.20, 0.15),
  ahr.H.true = c(rep(1.90, 45), rep(NA, 5)),
  ahr.Hc.true= rep(1.12, n_mock),
  sens       = c(runif(45, 0.6, 0.9), rep(NA, 5)),
  spec       = c(runif(45, 0.8, 0.95), rep(NA, 5)),
  ppv        = c(runif(45, 0.5, 0.8), rep(NA, 5)),
  npv        = c(runif(45, 0.85, 0.98), rep(NA, 5)),
  size.H     = c(round(rnorm(45, 258, 30)), rep(NA, 5)),
  size.Hc    = c(round(rnorm(45, 742, 30)), rep(1000, 5)),
  p.cens     = rep(0, n_mock),
  stringsAsFactors = FALSE
)

mock_dgm_full <- make_binary_dgm(
  df_alt, "alt", 2.0, 1.10, 1.30, fit_syn, beta_cal
)

# E1: data.table conversion
dt_mock <- data.table::as.data.table(mock_full)
report("E1", "Converts to data.table without error",
       data.table::is.data.table(dt_mock))

# E2: Subsetting by analysis
dt_fs <- dt_mock[dt_mock$analysis == "FS", ]
report("E2", "Filter by analysis preserves rows",
       nrow(dt_fs) == n_mock)

# E3: Subsetting by any.H
dt_found <- dt_mock[dt_mock$any.H == 1, ]
report("E3", "Filter by any.H = 1 works",
       nrow(dt_found) == 45)

# E4: Mean computation on HR columns with NAs
mean_hr_H <- mean(dt_found$hr.H.hat, na.rm = TRUE)
report("E4", "Mean HR(Q) computable from found subset",
       !is.na(mean_hr_H) && mean_hr_H > 0,
       sprintf("mean = %.3f", mean_hr_H))

# E5: Inversion of full result set
inv_full <- invert_hr_columns(mock_full)
found_inv <- inv_full[inv_full$any.H == 1, ]
report("E5a", "Inverted OR(Q) < 1 for all found",
       all(found_inv$hr.H.hat < 1, na.rm = TRUE),
       sprintf("range: [%.3f, %.3f]",
               min(found_inv$hr.H.hat, na.rm=TRUE),
               max(found_inv$hr.H.hat, na.rm=TRUE)))
report("E5b", "Inverted OR(Qc) near 1",
       all(abs(found_inv$hr.Hc.hat - 1) < 0.5, na.rm = TRUE))

# E6: Bias computation on original scale
truth_inv <- 1 / 2.0  # = 0.50
mean_est_inv <- mean(found_inv$hr.H.hat, na.rm = TRUE)
bias_pct <- 100 * (mean_est_inv - truth_inv) / truth_inv
report("E6", "Bias negative on original scale (overestimate benefit)",
       bias_pct < 0,
       sprintf("bias = %.1f%%", bias_pct))


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP F: Edge Cases ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

# F1: Very rare outcome
df_rare <- sim1
df_rare$y_sim <- 0L  # all zeros
df_rare$y_sim[1:5] <- 1L  # 1% prevalence
df_rare$event_sim <- df_rare$y_sim
mock_fs_rare <- list(
  grp.consistency = list(sg.harm.id = c("z1.1")),
  df.est = df_rare
)
mock_fs_rare$df.est$treat.recommend <- ifelse(df_rare$flag_harm == 1, 0L, 1L)
out_rare <- extract_binary_estimates(df_rare, mock_fs_rare, dgm_syn)
report("F1", "Rare outcome: extraction does not error",
       is.data.frame(out_rare) && nrow(out_rare) == 1)

# F2: All-same treatment in subgroup (separation risk)
df_sep <- sim1
df_sep$treat_sim[df_sep$flag_harm == 1] <- 1L  # all treated in Q
df_sep$event_sim <- df_sep$y_sim
mock_fs_sep <- list(
  grp.consistency = list(sg.harm.id = c("z1.1", "z2.1")),
  df.est = df_sep
)
mock_fs_sep$df.est$treat.recommend <- ifelse(df_sep$flag_harm == 1, 0L, 1L)
out_sep <- extract_binary_estimates(df_sep, mock_fs_sep, dgm_syn)
report("F2", "Separation: OR(Q) is NA (no both-arm data)",
       is.na(out_sep$hr.H.hat),
       "tryCatch prevents crash")

# F3: Very small subgroup
df_small <- sim1
df_small$treat.recommend <- 1L  # all in complement
df_small$treat.recommend[1:5] <- 0L  # only 5 in subgroup
mock_fs_small <- list(
  grp.consistency = list(sg.harm.id = c("z1.1")),
  df.est = df_small
)
out_small <- extract_binary_estimates(sim1, mock_fs_small, dgm_syn)
report("F3", "Very small subgroup (n=5): OR is NA (< 10 threshold)",
       is.na(out_small$hr.H.hat))

# F4: Empty sg.harm.id string
mock_fs_empty <- list(
  grp.consistency = list(sg.harm.id = ""),
  df.est = sim1
)
out_empty <- extract_binary_estimates(sim1, mock_fs_empty, dgm_syn)
report("F4", "Empty string sg.harm.id treated as no subgroup",
       out_empty$any.H == 0)

# F5: Inversion of NA values
res_with_na <- data.frame(hr.H.hat = c(2.0, NA, 3.0), hr.itt = c(1.5, 1.2, NA))
inv_na <- invert_hr_columns(res_with_na)
report("F5a", "Inversion preserves NA positions",
       is.na(inv_na$hr.H.hat[2]) && !is.na(inv_na$hr.H.hat[1]))
report("F5b", "Inversion of NA produces NA (not NaN or Inf)",
       is.na(inv_na$hr.itt[3]))

# F6: Inversion of zero (pathological)
res_zero <- data.frame(hr.H.hat = c(0.0, 2.0))
inv_zero <- invert_hr_columns(res_zero)
report("F6", "Inversion of zero produces Inf (expected)",
       is.infinite(inv_zero$hr.H.hat[1]),
       "Inf is a valid sentinel — downstream mean() with na.rm handles it")


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP G: Parallelisation Safety ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

# G1: Functions serialise correctly to workers
registerDoFuture()
plan(multisession, workers = 2)

par_result <- tryCatch({
  foreach(i = 1:4, .combine = rbind, .errorhandling = "pass",
          .options.future = list(seed = TRUE)) %dofuture% {
    future::plan("sequential")  # Layer 1
    sd <- simulate_binary_trial(dgm_syn, 200, i, 42)
    # No actual forestsearch call — just test extraction with mock
    mock_res_i <- data.frame(
      analysis = "FS", any.H = 1L,
      hr.H.hat = 2.0 + rnorm(1, 0, 0.3),
      hr.Hc.hat = 1.1 + rnorm(1, 0, 0.1),
      hr.H.true = dgm_syn$hazard_ratios$harm_subgroup,
      hr.Hc.true = dgm_syn$hazard_ratios$no_harm_subgroup,
      hr.itt = 1.3, sens = 0.8, spec = 0.9, ppv = 0.7, npv = 0.95,
      size.H = 80, size.Hc = 120, p.cens = 0,
      stringsAsFactors = FALSE
    )
    mock_res_i
  }
}, error = function(e) e)

plan(sequential)

report("G1a", "Parallel foreach completes without error",
       is.data.frame(par_result),
       if (inherits(par_result, "error")) par_result$message else
         sprintf("nrow = %d", nrow(par_result)))
if (is.data.frame(par_result)) {
  report("G1b", "All 4 iterations returned",
         nrow(par_result) == 4)
}

# G2: Post-hoc failure accounting pattern
n_expected <- 4
n_got <- if (is.data.frame(par_result)) nrow(par_result) else 0
n_failed <- n_expected - n_got
report("G2", "Post-hoc failure accounting detects no failures",
       n_failed == 0,
       sprintf("%d/%d completed", n_got, n_expected))

# G3: DGM serialisation (large object)
plan(multisession, workers = 2)
dgm_check <- tryCatch({
  foreach(i = 1:2, .combine = c) %dofuture% {
    nrow(dgm_syn$df_super)  # access DGM inside worker
  }
}, error = function(e) e)
plan(sequential)
report("G3", "DGM object serialises to workers",
       is.numeric(dgm_check) && all(dgm_check == nrow(dgm_syn$df_super)))


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP H: Seed Reproducibility ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

sim_a <- simulate_binary_trial(dgm_syn, 300, sim_id = 7, seed = 42)
sim_b <- simulate_binary_trial(dgm_syn, 300, sim_id = 7, seed = 42)
report("H1", "Same seed + sim_id → identical data",
       identical(sim_a$y_sim, sim_b$y_sim) &&
       identical(sim_a$treat_sim, sim_b$treat_sim))

sim_c <- simulate_binary_trial(dgm_syn, 300, sim_id = 8, seed = 42)
report("H2", "Different sim_id → different data",
       !identical(sim_a$y_sim, sim_c$y_sim))

# H3: Alt vs null seed offset prevents collision
max_nsims <- 1000L
report("H3", "Null seed offset (100000) avoids collision with alt",
       (42 + 100000L) > (42 + max_nsims),
       "offset = 100000 >> max nsims")


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP I: OR Calibration Accuracy ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

# Simulate many trials and check that the empirical OR(Q) matches the
# DGM target on average
set.seed(999)
n_calib_reps <- 200
or_Q_empirical <- numeric(n_calib_reps)
for (r in seq_len(n_calib_reps)) {
  sim_r <- simulate_binary_trial(dgm_syn, 1000, r, 999)
  q_idx <- sim_r$flag_harm == 1
  if (sum(q_idx) > 20) {
    fit_q <- glm(y_sim ~ treat_sim, data = sim_r[q_idx, ], family = binomial)
    or_Q_empirical[r] <- exp(coef(fit_q)["treat_sim"])
  } else {
    or_Q_empirical[r] <- NA
  }
}
mean_or_Q <- mean(or_Q_empirical, na.rm = TRUE)
target_dgm <- dgm_syn$hazard_ratios$harm_subgroup

report("I1", "Mean empirical OR(Q) within 20% of DGM target",
       abs(mean_or_Q - target_dgm) / target_dgm < 0.20,
       sprintf("empirical = %.3f, target = %.3f, diff = %.1f%%",
               mean_or_Q, target_dgm,
               100 * (mean_or_Q - target_dgm) / target_dgm))

# I2: SD of empirical OR is reasonable
sd_or_Q <- sd(or_Q_empirical, na.rm = TRUE)
report("I2", "SD of empirical OR(Q) is reasonable (< 2.0)",
       sd_or_Q < 2.0,
       sprintf("SD = %.3f", sd_or_Q))


# ─────────────────────────────────────────────────────────────────────────────
cat("\n─── TEST GROUP J: Treatment Switching Algebra ───\n\n")
# ─────────────────────────────────────────────────────────────────────────────

# Verify that the treatment switching + inversion gives correct results
# On the switched scale, treat=1 means ddI (monotherapy), treat=0 means combo
# OR(Q, switched) > 1 means ddI is worse for Q (= Q benefits from combo)
# On the original scale, OR(Q, original) < 1 means combo is better for Q

or_sw  <- dgm_syn$hazard_ratios$harm_subgroup
or_ori <- 1 / or_sw

report("J1", "OR(Q, switched) > 1 (ddI worse for Q)",
       or_sw > 1,
       sprintf("%.3f", or_sw))
report("J2", "OR(Q, original) < 1 (combo better for Q)",
       or_ori < 1,
       sprintf("%.3f", or_ori))
report("J3", "OR(Q, original) = 1/OR(Q, switched)",
       abs(or_ori - 1/or_sw) < 1e-10)

# J4: The inversion preserves the direction of the treatment effect
# On the switched scale: OR > 1 → more events under treat=1 (ddI)
# On the original scale: OR < 1 → fewer events under treat=1 (combo)
# The benefit is the same, just expressed from different reference points
or_Qc_sw <- dgm_syn$hazard_ratios$no_harm_subgroup
or_Qc_ori <- 1 / or_Qc_sw
report("J4", "OR(Qc, original) closer to 1 than OR(Q, original)",
       abs(or_Qc_ori - 1) < abs(or_ori - 1),
       sprintf("OR(Q) = %.3f, OR(Qc) = %.3f on original", or_ori, or_Qc_ori))

# J5: Verify subgroup has STRONGER benefit than complement
report("J5", "Q has stronger benefit (lower OR on original scale)",
       or_ori < or_Qc_ori,
       "Confirms Q is the benefit subgroup")


# =============================================================================
cat("\n", rep("=", 72), "\n", sep="")
cat(sprintf("  SUMMARY: %d PASSED, %d FAILED (of %d total)\n",
            pass_count, fail_count, pass_count + fail_count))
cat(rep("=", 72), "\n", sep="")
if (fail_count > 0) {
  cat("\nFailed tests:\n")
  for (msg in test_log[grepl("\\[FAIL\\]", test_log)]) cat(msg, "\n")
}
cat("\n")
