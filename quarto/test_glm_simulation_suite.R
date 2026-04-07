# =============================================================================
# Smoke Tests: GLM Simulation Suite
# =============================================================================
#
# Usage:
#   devtools::load_all()          # or devtools::install()
#   source("tests/test_glm_simulation_suite.R")
#
# Requires: speff2trial (for ACTG175 data)
#
# Test groups:
#   1. generate_glm_dgm -- construction, structure, effects
#   2. simulate_from_glm_dgm -- sampling, outcomes, reproducibility
#   3. calibrate_glm_interaction -- grid search, convergence
#   4. Integration -- DGM -> simulate -> forestsearch pipeline
#   5. Edge cases -- degenerate data, boundary conditions
# =============================================================================

cat("\n")
cat("=============================================================\n")
cat("  GLM Simulation Suite -- Smoke Tests\n")
cat("=============================================================\n\n")

library(speff2trial)

n_pass <- 0L
n_fail <- 0L

pass <- function(name) {
  n_pass <<- n_pass + 1L
  cat(sprintf("  PASS: %s\n", name))
}

fail <- function(name, msg = "") {
  n_fail <<- n_fail + 1L
  cat(sprintf("  FAIL: %s %s\n", name, msg))
}

check <- function(name, expr) {
  result <- tryCatch({
    force(expr)
    TRUE
  }, error = function(e) {
    fail(name, paste("--", e$message))
    FALSE
  })
  if (result) pass(name)
  invisible(result)
}


# =============================================================================
# Prepare ACTG175 data (shared across all tests)
# =============================================================================

cat("--- Preparing ACTG175 data ---\n")

actg_df <- subset(ACTG175, arms %in% c(1, 3))
actg_df$id <- seq_len(nrow(actg_df))
actg_df$treat_orig <- ifelse(actg_df$arms == 1, 1L, 0L)
actg_df$treat <- 1L - actg_df$treat_orig
actg_df$y_binary <- as.integer(actg_df$cd420 > actg_df$cd40)
actg_df <- actg_df[!is.na(actg_df$cd420), ]

actg_df$z1  <- as.factor(ifelse(actg_df$age > 34, 1L, 0L))
actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
actg_df$z4  <- as.factor(ifelse(
  actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
actg_df$z5  <- as.factor(ifelse(
  actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
actg_df$z6  <- as.factor(ifelse(
  actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))

z_vars <- paste0("z", 1:6)
sg_vars <- c("z1", "z2")
sg_cuts <- list(z1 = 1L, z2 = 1L)

cat(sprintf("  N = %d, Y-bar = %.3f, Q-prevalence = %.3f\n\n",
    nrow(actg_df), mean(actg_df$y_binary),
    mean(actg_df$z1 == 1 & actg_df$z2 == 1)))


# =============================================================================
# GROUP 1: generate_glm_dgm
# =============================================================================

cat("--- Group 1: generate_glm_dgm ---\n")

# 1.1 Alt model returns correct class
check("1.1 alt model returns glm_dgm class", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    outcome_type = "binary", subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, model = "alt", k_inter = 1,
    n_super = 1000L, seed = 42L, verbose = FALSE
  )
  stopifnot(inherits(dgm, "glm_dgm"))
})

# 1.2 Required fields present
check("1.2 required fields present", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  stopifnot(all(c("df_super", "hazard_ratios", "outcome_type",
                   "effect_measure", "model_params", "subgroup_info",
                   "model_type", "n_super", "factor_vars") %in% names(dgm)))
})

# 1.3 df_super has potential outcome columns
check("1.3 df_super has p0, p1, flag_harm", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  stopifnot(all(c("p0", "p1", "flag_harm") %in% names(dgm$df_super)))
})

# 1.4 Potential outcomes are valid probabilities
check("1.4 p0, p1 in [0, 1]", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  stopifnot(
    all(dgm$df_super$p0 >= 0 & dgm$df_super$p0 <= 1),
    all(dgm$df_super$p1 >= 0 & dgm$df_super$p1 <= 1)
  )
})

# 1.5 hazard_ratios has required fields
check("1.5 hazard_ratios structure", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  hr <- dgm$hazard_ratios
  stopifnot(
    all(c("harm_subgroup", "no_harm_subgroup", "overall") %in% names(hr)),
    is.numeric(hr$harm_subgroup),
    is.numeric(hr$no_harm_subgroup),
    is.numeric(hr$overall),
    hr$harm_subgroup > 0,
    hr$no_harm_subgroup > 0,
    hr$overall > 0
  )
})

# 1.6 n_super respected
check("1.6 n_super = 2000 produces 2000 rows", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 2000L, seed = 42L
  )
  stopifnot(nrow(dgm$df_super) == 2000L)
})

# 1.7 Null model: Q and Qc effects approximately equal
check("1.7 null model: OR(Q) ~= OR(Qc)", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    model = "null", n_super = 5000L, seed = 42L
  )
  diff <- abs(dgm$hazard_ratios$harm_subgroup -
              dgm$hazard_ratios$no_harm_subgroup)
  stopifnot(diff < 0.05)
})

# 1.8 Null model: flag_harm is still populated (for classification)
check("1.8 null model preserves flag_harm", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    model = "null", n_super = 2000L, seed = 42L
  )
  stopifnot(
    "flag_harm" %in% names(dgm$df_super),
    sum(dgm$df_super$flag_harm) > 0
  )
})

# 1.9 Alt model with k_inter > 1 creates separation
check("1.9 alt model: OR(Q) != OR(Qc) when k_inter = 3", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    model = "alt", k_inter = 3, n_super = 2000L, seed = 42L
  )
  diff <- abs(dgm$hazard_ratios$harm_subgroup -
              dgm$hazard_ratios$no_harm_subgroup)
  stopifnot(diff > 0.1)
})

# 1.10 Single-level factors dropped silently
check("1.10 single-level factors dropped", {
  # z4 (karnof) is single-level in this ACTG subset
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = paste0("z", 1:6),
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  stopifnot(!"z4" %in% dgm$factor_vars)
  stopifnot(length(dgm$factor_vars) < 6)
})

# 1.11 Seed reproducibility
check("1.11 same seed = identical DGM", {
  dgm_a <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 99L
  )
  dgm_b <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 99L
  )
  stopifnot(
    identical(dgm_a$hazard_ratios, dgm_b$hazard_ratios),
    identical(dgm_a$df_super$p0, dgm_b$df_super$p0),
    identical(dgm_a$df_super$p1, dgm_b$df_super$p1)
  )
})

# 1.12 get_dgm_hr() compatibility
check("1.12 get_dgm_hr compatibility", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  hr_H  <- get_dgm_hr(dgm, "hr_H")
  hr_Hc <- get_dgm_hr(dgm, "hr_Hc")
  hr_o  <- get_dgm_hr(dgm, "hr_overall")
  stopifnot(!is.na(hr_H), !is.na(hr_Hc), !is.na(hr_o))
})

# 1.13 print method works
check("1.13 print.glm_dgm runs without error", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  invisible(capture.output(print(dgm)))
})


# =============================================================================
# GROUP 2: simulate_from_glm_dgm
# =============================================================================

cat("\n--- Group 2: simulate_from_glm_dgm ---\n")

dgm_test <- generate_glm_dgm(
  data = actg_df, factor_vars = z_vars,
  outcome_var = "y_binary", treatment_var = "treat",
  subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
  model = "alt", k_inter = 2, n_super = 2000L, seed = 42L
)

# 2.1 Returns data.frame with correct N
check("2.1 returns N = 500", {
  df <- simulate_from_glm_dgm(dgm_test, n = 500, seed = 1)
  stopifnot(is.data.frame(df), nrow(df) == 500)
})

# 2.2 Required columns present
check("2.2 required columns: id, y_sim, treat_sim, flag_harm", {
  df <- simulate_from_glm_dgm(dgm_test, n = 500, seed = 1)
  stopifnot(all(c("id", "y_sim", "treat_sim", "flag_harm") %in% names(df)))
})

# 2.3 y_sim is binary
check("2.3 y_sim is 0/1", {
  df <- simulate_from_glm_dgm(dgm_test, n = 500, seed = 1)
  stopifnot(all(df$y_sim %in% c(0L, 1L)))
})

# 2.4 treat_sim is binary
check("2.4 treat_sim is 0/1", {
  df <- simulate_from_glm_dgm(dgm_test, n = 500, seed = 1)
  stopifnot(all(df$treat_sim %in% c(0L, 1L)))
})

# 2.5 Treatment allocation near 50%
check("2.5 treatment near 50% (rand_ratio = 1)", {
  df <- simulate_from_glm_dgm(dgm_test, n = 2000, seed = 1)
  p_treat <- mean(df$treat_sim)
  stopifnot(p_treat > 0.45, p_treat < 0.55)
})

# 2.6 Unequal randomisation
check("2.6 rand_ratio = 2 gives ~67% treatment", {
  df <- simulate_from_glm_dgm(dgm_test, n = 3000, seed = 1,
                                rand_ratio = 2)
  p_treat <- mean(df$treat_sim)
  stopifnot(p_treat > 0.60, p_treat < 0.73)
})

# 2.7 n = NULL returns full super-population
check("2.7 n = NULL returns n_super rows", {
  df <- simulate_from_glm_dgm(dgm_test, n = NULL, seed = 1)
  stopifnot(nrow(df) == dgm_test$n_super)
})

# 2.8 Seed reproducibility
check("2.8 same seed = identical simulation", {
  a <- simulate_from_glm_dgm(dgm_test, n = 200, seed = 77)
  b <- simulate_from_glm_dgm(dgm_test, n = 200, seed = 77)
  stopifnot(
    identical(a$y_sim, b$y_sim),
    identical(a$treat_sim, b$treat_sim),
    identical(a$id, b$id)
  )
})

# 2.9 Different seeds produce different data
check("2.9 different seeds = different outcomes", {
  a <- simulate_from_glm_dgm(dgm_test, n = 200, seed = 1)
  b <- simulate_from_glm_dgm(dgm_test, n = 200, seed = 2)
  stopifnot(!identical(a$y_sim, b$y_sim))
})

# 2.10 IDs are sequential
check("2.10 IDs are 1:n", {
  df <- simulate_from_glm_dgm(dgm_test, n = 300, seed = 1)
  stopifnot(identical(df$id, 1:300))
})

# 2.11 Factor vars preserved in simulated data
check("2.11 factor variables present in output", {
  df <- simulate_from_glm_dgm(dgm_test, n = 100, seed = 1)
  fv <- dgm_test$factor_vars
  stopifnot(all(fv %in% names(df)))
})

# 2.12 Simulated OR roughly matches DGM ITT
check("2.12 simulated ITT OR consistent with DGM (large n)", {
  df <- simulate_from_glm_dgm(dgm_test, n = 5000, seed = 1)
  fit <- glm(y_sim ~ treat_sim, data = df, family = binomial)
  or_sim <- exp(coef(fit)["treat_sim"])
  or_dgm <- dgm_test$hazard_ratios$overall
  rel_err <- abs(or_sim - or_dgm) / or_dgm
  stopifnot(rel_err < 0.15)
})


# =============================================================================
# GROUP 3: calibrate_glm_interaction
# =============================================================================

cat("\n--- Group 3: calibrate_glm_interaction ---\n")

# 3.1 Calibration returns glm_dgm
check("3.1 returns glm_dgm class", {
  dgm <- calibrate_glm_interaction(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    target_effect = 2.0, subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, grid_step = 0.2,
    n_super = 1000L, seed = 42L
  )
  stopifnot(inherits(dgm, "glm_dgm"))
})

# 3.2 Calibrated OR(Q) near target
check("3.2 calibrated OR(Q) within 0.15 of target = 2.0", {
  dgm <- calibrate_glm_interaction(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    target_effect = 2.0, subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, grid_step = 0.1,
    n_super = 2000L, seed = 42L
  )
  or_Q <- dgm$hazard_ratios$harm_subgroup
  stopifnot(abs(or_Q - 2.0) < 0.15)
})

# 3.3 Calibration for OR < 1
check("3.3 calibration for target OR = 0.5", {
  dgm <- calibrate_glm_interaction(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    target_effect = 0.5, subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, grid_step = 0.1,
    n_super = 2000L, seed = 42L
  )
  or_Q <- dgm$hazard_ratios$harm_subgroup
  stopifnot(abs(or_Q - 0.5) < 0.15)
})

# 3.4 Calibrated DGM is simulable
check("3.4 calibrated DGM produces valid simulations", {
  dgm <- calibrate_glm_interaction(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    target_effect = 2.0, subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, grid_step = 0.2,
    n_super = 1000L, seed = 42L
  )
  df <- simulate_from_glm_dgm(dgm, n = 500, seed = 1)
  stopifnot(
    nrow(df) == 500,
    all(df$y_sim %in% c(0L, 1L)),
    all(c("flag_harm", "treat_sim") %in% names(df))
  )
})

# 3.5 Finer grid gives closer result
check("3.5 finer grid improves precision", {
  dgm_coarse <- calibrate_glm_interaction(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    target_effect = 1.5, subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, grid_step = 0.5,
    n_super = 2000L, seed = 42L
  )
  dgm_fine <- calibrate_glm_interaction(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    target_effect = 1.5, subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, grid_step = 0.05,
    n_super = 2000L, seed = 42L
  )
  err_coarse <- abs(dgm_coarse$hazard_ratios$harm_subgroup - 1.5)
  err_fine   <- abs(dgm_fine$hazard_ratios$harm_subgroup - 1.5)
  stopifnot(err_fine <= err_coarse)
})


# =============================================================================
# GROUP 4: Integration with run_simulation_analysis
# =============================================================================

cat("\n--- Group 4: run_simulation_analysis integration ---\n")

# 4.1 Class dispatch: glm_dgm triggers simulate_from_glm_dgm
check("4.1 inherits(dgm, 'glm_dgm') is TRUE", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  stopifnot(inherits(dgm, "glm_dgm"))
})

# 4.2 run_simulation_analysis dispatches correctly
check("4.2 run_simulation_analysis works with glm_dgm", {
  dgm <- calibrate_glm_interaction(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    target_effect = 2.0, subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, grid_step = 0.2,
    n_super = 2000L, seed = 42L
  )
  result <- run_simulation_analysis(
    sim_id = 1, dgm = dgm, n_sample = 500,
    confounders_base = dgm$factor_vars,
    run_fs = TRUE, run_fs_grf = FALSE, run_grf = FALSE,
    fs_params = list(
      outcome_type = "binary",
      effect_measure = "OR",
      hr.threshold = 1.667,
      hr.consistency = 1.25,
      pconsistency.threshold = 0.90,
      fs.splits = 100L,
      n.min = 40,
      d0.min = 10,
      d1.min = 10,
      maxk = 2,
      use_lasso = FALSE,
      use_grf = TRUE,
      vi.grf.min = -0.2,
      parallel_args = list(
        plan = "sequential", workers = 1, show_message = FALSE
      )
    ),
    verbose = FALSE
  )
  stopifnot(
    is.data.frame(result) || data.table::is.data.table(result),
    "analysis" %in% names(result),
    "any.H" %in% names(result),
    "hr.H.hat" %in% names(result),
    nrow(result) >= 1
  )
})

# 4.3 GLM params auto-injected from DGM
check("4.3 outcome_type auto-injected when not in fs_params", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    n_super = 1000L, seed = 42L
  )
  # Minimal fs_params -- outcome_type should come from DGM
  result <- run_simulation_analysis(
    sim_id = 1, dgm = dgm, n_sample = 300,
    confounders_base = dgm$factor_vars,
    run_fs = TRUE, run_fs_grf = FALSE, run_grf = FALSE,
    fs_params = list(
      hr.threshold = 1.5,
      hr.consistency = 1.25,
      fs.splits = 50L,
      n.min = 30,
      maxk = 1,
      use_lasso = FALSE,
      use_grf = TRUE,
      parallel_args = list(
        plan = "sequential", workers = 1, show_message = FALSE
      )
    ),
    verbose = FALSE
  )
  stopifnot(nrow(result) >= 1)
})

# 4.4 Result has valid classification columns
check("4.4 classification metrics present", {
  dgm <- calibrate_glm_interaction(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    target_effect = 2.0, subgroup_vars = sg_vars,
    subgroup_cuts = sg_cuts, grid_step = 0.5,
    n_super = 2000L, seed = 42L
  )
  result <- run_simulation_analysis(
    sim_id = 1, dgm = dgm, n_sample = 800,
    confounders_base = dgm$factor_vars,
    run_fs = TRUE, run_fs_grf = FALSE, run_grf = FALSE,
    fs_params = list(
      hr.threshold = 1.667,
      hr.consistency = 1.25,
      fs.splits = 100L,
      n.min = 40,
      maxk = 2,
      use_lasso = FALSE,
      use_grf = TRUE,
      parallel_args = list(
        plan = "sequential", workers = 1, show_message = FALSE
      )
    ),
    verbose = FALSE
  )
  stopifnot(all(c("sens", "spec", "ppv", "npv") %in% names(result)))
})


# =============================================================================
# GROUP 5: Edge cases
# =============================================================================

cat("\n--- Group 5: Edge cases ---\n")

# 5.1 No subgroup vars (ITT-only DGM)
check("5.1 generate_glm_dgm with no subgroup", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = NULL, subgroup_cuts = NULL,
    model = "null", n_super = 500L, seed = 42L
  )
  stopifnot(
    inherits(dgm, "glm_dgm"),
    all(dgm$df_super$flag_harm == 0)
  )
})

# 5.2 Minimal factor set (one variable)
check("5.2 single factor variable", {
  dgm <- generate_glm_dgm(
    data = actg_df, factor_vars = "z1",
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
    n_super = 500L, seed = 42L
  )
  stopifnot(inherits(dgm, "glm_dgm"))
  df <- simulate_from_glm_dgm(dgm, n = 100, seed = 1)
  stopifnot(nrow(df) == 100)
})

# 5.3 Invalid outcome_type
check("5.3 invalid outcome_type raises error", {
  err <- tryCatch(
    generate_glm_dgm(
      data = actg_df, factor_vars = z_vars,
      outcome_var = "y_binary", treatment_var = "treat",
      outcome_type = "survival"
    ),
    error = function(e) TRUE
  )
  stopifnot(isTRUE(err))
})

# 5.4 Non-binary outcome raises error
check("5.4 non-binary outcome raises error", {
  actg_df$y_cont <- actg_df$cd40
  err <- tryCatch(
    generate_glm_dgm(
      data = actg_df, factor_vars = z_vars,
      outcome_var = "y_cont", treatment_var = "treat",
      outcome_type = "binary"
    ),
    error = function(e) TRUE
  )
  stopifnot(isTRUE(err))
  actg_df$y_cont <- NULL
})

# 5.5 Missing column raises error
check("5.5 missing outcome column raises error", {
  err <- tryCatch(
    generate_glm_dgm(
      data = actg_df, factor_vars = z_vars,
      outcome_var = "nonexistent", treatment_var = "treat"
    ),
    error = function(e) TRUE
  )
  stopifnot(isTRUE(err))
})

# 5.6 simulate_from_glm_dgm rejects non-glm_dgm
check("5.6 simulate rejects wrong class", {
  err <- tryCatch(
    simulate_from_glm_dgm(list(a = 1), n = 100),
    error = function(e) TRUE
  )
  stopifnot(isTRUE(err))
})

# 5.7 Null + alt DGMs produce rbind-compatible results
check("5.7 null + alt results are rbind-compatible", {
  dgm_a <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    model = "alt", k_inter = 2, n_super = 1000L, seed = 42L
  )
  dgm_n <- generate_glm_dgm(
    data = actg_df, factor_vars = z_vars,
    outcome_var = "y_binary", treatment_var = "treat",
    subgroup_vars = sg_vars, subgroup_cuts = sg_cuts,
    model = "null", n_super = 1000L, seed = 42L
  )
  r_a <- run_simulation_analysis(
    sim_id = 1, dgm = dgm_a, n_sample = 300,
    confounders_base = dgm_a$factor_vars,
    run_fs = TRUE, run_fs_grf = FALSE, run_grf = FALSE,
    fs_params = list(
      hr.threshold = 1.5, fs.splits = 50L,
      n.min = 30, maxk = 1,
      use_lasso = FALSE, use_grf = TRUE,
      parallel_args = list(plan = "sequential", workers = 1,
                           show_message = FALSE)
    ), verbose = FALSE
  )
  r_n <- run_simulation_analysis(
    sim_id = 1, dgm = dgm_n, n_sample = 300,
    confounders_base = dgm_n$factor_vars,
    run_fs = TRUE, run_fs_grf = FALSE, run_grf = FALSE,
    fs_params = list(
      hr.threshold = 1.5, fs.splits = 50L,
      n.min = 30, maxk = 1,
      use_lasso = FALSE, use_grf = TRUE,
      parallel_args = list(plan = "sequential", workers = 1,
                           show_message = FALSE)
    ), verbose = FALSE
  )
  combined <- rbind(r_a, r_n)
  stopifnot(nrow(combined) == 2)
})


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=============================================================\n")
cat(sprintf("  %d PASSED, %d FAILED, %d total\n",
    n_pass, n_fail, n_pass + n_fail))
if (n_fail == 0) {
  cat("  All tests passed.\n")
} else {
  cat("  *** FAILURES DETECTED ***\n")
}
cat("=============================================================\n\n")

if (n_fail > 0) stop(sprintf("%d test(s) failed", n_fail))
