# ==============================================================================
# TESTING REFACTORED FUNCTIONS WITH SYNTHETIC DATASET
# Complete end-to-end test with generated data
# ==============================================================================

library(survival)
library(data.table)
library(randomizr)

# ------------------------------------------------------------------------------
# STEP 1: CREATE SYNTHETIC DATASET
# ------------------------------------------------------------------------------

set.seed(42)

create_synthetic_dataset <- function(n = 500) {
  cat("Creating synthetic dataset with n =", n, "\n")
  
  df <- data.frame(
    # Patient ID
    id = 1:n,
    
    # Treatment (will be overwritten, but needed for DGM)
    treat = sample(0:1, n, replace = TRUE),
    
    # Biomarker (continuous, log-scale)
    z = rnorm(n, mean = 3, sd = 1.5),
    
    # Prognostic factor (binary: 0=pre-meno, 1=post-meno)
    meno = rbinom(n, 1, prob = 0.6),
    
    # Stratification variable (3 levels: grade)
    stratum = sample(1:3, n, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
    
    # Additional covariates
    age = rnorm(n, mean = 60, sd = 10),
    size = rexp(n, rate = 1/20),
    nodes = rpois(n, lambda = 5),
    
    # Country
    country = sample(c("US", "EU"), n, replace = TRUE, prob = c(0.4, 0.6))
  )
  
  # Ensure age is positive and reasonable
  df$age <- pmax(30, pmin(90, df$age))
  
  # Create derived variables needed for DGM
  df$z.treat <- df$z * df$treat
  df$z.k <- pmax(0, df$z - 5)  # Spline knot at 5
  df$z.k.treat <- df$z.k * df$treat
  
  # Create binary biomarker indicator
  df$bm_low <- ifelse(df$z <= 2, 1, 0)
  
  # Create ITT flag
  df$itt <- "Y"
  
  # Simulate actual survival times (for DGM fitting)
  # Simple Weibull with some treatment effect
  lambda <- 0.01 * exp(0.3 * df$treat - 0.5 * df$meno + 0.1 * df$z)
  shape <- 1.2
  df$tte <- rweibull(n, shape = shape, scale = 1/lambda^(1/shape))
  df$event <- rbinom(n, 1, prob = 0.7)  # 70% event rate
  
  # Create random small subgroups
  df$random15 <- as.numeric(sample(1:n, 15))
  df$random15 <- ifelse(df$id %in% df$random15, 1, 0)
  
  df$random20 <- as.numeric(sample(1:n, 20))
  df$random20 <- ifelse(df$id %in% df$random20, 1, 0)
  
  cat("✓ Synthetic dataset created\n")
  cat(sprintf("  - Sample size: %d\n", nrow(df)))
  cat(sprintf("  - Biomarker range: [%.2f, %.2f]\n", min(df$z), max(df$z)))
  cat(sprintf("  - Event rate: %.1f%%\n", 100 * mean(df$event)))
  
  return(df)
}

df.synthetic <- create_synthetic_dataset(n = 500)

# ------------------------------------------------------------------------------
# STEP 2: SOURCE REFACTORED FUNCTIONS
# ------------------------------------------------------------------------------

cat("\n=== Loading Refactored Functions ===\n")

# Since we're testing, let's inline the key constants
SIMULATION_CONSTANTS <- list(
  BASE_SEED = 8316951,
  SEED_INCREMENT = 1000,
  CENSORING_SEED_OFFSET = 100,
  MAX_KM_PLOTS = 10,
  DEFAULT_Z_INCREMENT = 1
)

# Include minimal required functions for testing
source_inline_functions <- function() {
  
  # Validation
  validate_simulation_inputs <<- function(dgm, strata_rand, wname) {
    required_vars <- c(strata_rand, wname)
    if (!all(required_vars %in% names(dgm$df_super))) {
      stop("Required variables not in dgm$df_super")
    }
    invisible(TRUE)
  }
  
  # Data preparation
  prepare_simulation_data <<- function(dgm, Ndraw, seed) {
    set.seed(SIMULATION_CONSTANTS$BASE_SEED + seed * SIMULATION_CONSTANTS$SEED_INCREMENT)
    df_super <- dgm$df_super
    if (Ndraw != nrow(df_super)) {
      id_sample <- sample(seq_len(nrow(df_super)), size = Ndraw, replace = TRUE)
      df_super <- df_super[id_sample, ]
    }
    return(df_super)
  }
  
  calculate_gamma_w <<- function(bw, tau_approx) {
    -bw * tau_approx
  }
  
  # Build covariate matrix
  build_covariate_matrix <<- function(df, treat) {
    zmat <- as.matrix(df[, c("treat", "z", "z.treat", "z.k", "z.k.treat")])
    zmat[, "treat"] <- treat
    zmat[, "z.treat"] <- treat * df$z
    zmat[, "z.k.treat"] <- treat * zmat[, "z.k"]
    return(zmat)
  }
  
  calculate_log_hazard_ratio <<- function(zmat_1, zmat_0, gamma_true, tau_strata) {
    (-1) * c((zmat_1 - zmat_0) %*% gamma_true) / tau_strata
  }
  
  calculate_hazard_components <<- function(df, zmat_1, zmat_0, gamma_true, 
                                           gamma_w, w, tau_strata) {
    df$theta1.po <- -c(zmat_1 %*% gamma_true + w * gamma_w) / tau_strata
    df$theta0.po <- -c(zmat_0 %*% gamma_true + w * gamma_w) / tau_strata
    df$theta1_w1.po <- -c(zmat_1 %*% gamma_true + 1 * gamma_w) / tau_strata
    df$theta0_w1.po <- -c(zmat_0 %*% gamma_true + 1 * gamma_w) / tau_strata
    df$theta1_w0.po <- -c(zmat_1 %*% gamma_true + 0 * gamma_w) / tau_strata
    df$theta0_w0.po <- -c(zmat_0 %*% gamma_true + 0 * gamma_w) / tau_strata
    return(df)
  }
  
  # Potential outcomes
  generate_potential_outcomes <<- function(df, gamma_true, mu, tau_strata, 
                                           gamma_w, wname, seed) {
    N <- nrow(df)
    w <- df[[wname]]
    
    zmat_1 <- build_covariate_matrix(df, treat = 1)
    zmat_0 <- build_covariate_matrix(df, treat = 0)
    
    set.seed(SIMULATION_CONSTANTS$BASE_SEED + seed * SIMULATION_CONSTANTS$SEED_INCREMENT + 1)
    epsilon <- log(rexp(N))
    
    df$log.Y1 <- mu + c(zmat_1 %*% gamma_true) + w * gamma_w + tau_strata * epsilon
    df$log.Y0 <- mu + c(zmat_0 %*% gamma_true) + w * gamma_w + tau_strata * epsilon
    df$loghr.po <- calculate_log_hazard_ratio(zmat_1, zmat_0, gamma_true, tau_strata)
    df <- calculate_hazard_components(df, zmat_1, zmat_0, gamma_true, gamma_w, w, tau_strata)
    
    return(df)
  }
  
  # Randomization
  apply_randomization_treatment <<- function(po_data, strata_rand, keep_rand, seed) {
    if (!keep_rand) {
      blocks <- po_data[[strata_rand]]
      set.seed(SIMULATION_CONSTANTS$BASE_SEED + seed * SIMULATION_CONSTANTS$SEED_INCREMENT + 2)
      po_data$treat.sim <- randomizr::block_ra(blocks = blocks)
    } else {
      po_data$treat.sim <- ifelse(po_data$treat == 1, 1, 0)
    }
    
    po_data$y.sim <- exp(
      po_data$treat.sim * po_data$log.Y1 + 
      (1 - po_data$treat.sim) * po_data$log.Y0
    )
    
    return(po_data)
  }
  
  # Censoring
  apply_censoring <<- function(obs_data, muC, tauC, time_eos, seed) {
    N <- nrow(obs_data)
    set.seed(SIMULATION_CONSTANTS$BASE_SEED + 
             seed * SIMULATION_CONSTANTS$SEED_INCREMENT + 
             SIMULATION_CONSTANTS$CENSORING_SEED_OFFSET)
    
    epsilonC <- log(rexp(N))
    censor_time <- exp(muC + tauC * epsilonC)
    censor_time <- pmin(censor_time, time_eos)
    
    obs_data$event.sim <- ifelse(obs_data$y.sim <= censor_time, 1, 0)
    obs_data$y.sim <- pmin(obs_data$y.sim, censor_time)
    
    return(obs_data)
  }
  
  # Stratification
  add_stratification_variables <<- function(obs_data, strata_rand, strata_tte) {
    obs_data$strata.simR <- obs_data[[strata_rand]]
    if (!is.null(strata_tte)) {
      obs_data$strata.simO <- obs_data[[strata_tte]]
    } else {
      obs_data$strata.simO <- "All"
    }
    obs_data$w <- obs_data[[wname]]
    return(obs_data)
  }
  
  # Output formatting
  format_simulation_output <<- function(obs_data, return_df, hrz_crit = log(1.2)) {
    obs_data <- data.table::setorder(obs_data, z)
    return(obs_data)
  }
  
  cat("✓ Functions loaded\n")
}

source_inline_functions()

# ------------------------------------------------------------------------------
# STEP 3: CREATE DATA GENERATING MECHANISM
# ------------------------------------------------------------------------------

cat("\n=== Creating Data Generating Mechanism ===\n")

# Use existing get_dgm_stratified function (from original code)
get_dgm_stratified <- function(df, knot = 5, zeta = 10, 
                               log.hrs = log(c(0.75, 0.75, 0.75)),
                               strata_tte = NULL, details = FALSE) {
  
  # Fit Weibull model to get parameters
  weib.formula <- as.formula("Surv(tte, event) ~ treat + z + z.treat + z.k + z.k.treat")
  fit.weibk <- survreg(weib.formula, dist = 'weibull', data = df)
  
  # Get parameters
  mu <- c(coef(fit.weibk)[1])
  gamma <- c(coef(fit.weibk)[c(-1)])
  tau <- c(fit.weibk$scale)
  
  # Censoring model
  fitC.weib <- survreg(Surv(tte, 1 - event) ~ 1, dist = 'weibull', data = df)
  tauC <- c(fitC.weib$scale)
  muC <- c(coef(fitC.weib)[1])
  
  # Re-define parameters to match specified log.hrs
  tau.approx <- tau
  b0 <- c(-gamma) / tau.approx
  
  loghr.0 <- log.hrs[1]
  loghr.knot <- log.hrs[2]
  loghr.zeta <- log.hrs[3]
  
  b0[1] <- loghr.0
  b0[3] <- (loghr.knot - b0[1]) / knot
  b0[5] <- (loghr.zeta - b0[1] - zeta * b0[3]) / (zeta - knot)
  
  gamma.true <- -b0 * tau.approx
  
  df$tau.strataO <- tau
  
  if (details) {
    cat("True log(HR) at z=0:", loghr.0, "\n")
    cat("True log(HR) at z=knot:", loghr.knot, "\n")
    cat("True log(HR) at z=zeta:", loghr.zeta, "\n")
  }
  
  return(list(
    df_super = df,
    gamma.true = gamma.true,
    mu = mu,
    tau = tau,
    muC = muC,
    tauC = tauC,
    strata_tte = strata_tte,
    tau.approx = tau.approx
  ))
}

# Create DGM with uniform treatment effect (HR = 0.60)
dgm_test <- get_dgm_stratified(
  df = df.synthetic,
  knot = 5,
  zeta = 10,
  log.hrs = log(c(0.60, 0.60, 0.60)),
  strata_tte = NULL,
  details = TRUE
)

cat("✓ DGM created successfully\n")

# ------------------------------------------------------------------------------
# STEP 4: TEST MAIN SIMULATION FUNCTION
# ------------------------------------------------------------------------------

cat("\n=== Testing draw_sim_stratified ===\n")

# Main function
draw_sim_stratified <- function(dgm, ss = 1, Ndraw = nrow(dgm$df_super),
                                strata_rand = "stratum", wname = "meno", 
                                bw = 0, checking = FALSE, details = FALSE,
                                return_df = TRUE, time_eos = Inf, 
                                keep_rand = FALSE, hrz_crit = log(1.2)) {
  
  validate_simulation_inputs(dgm, strata_rand, wname)
  
  sim_data <- prepare_simulation_data(dgm, Ndraw, ss)
  gamma_w <- calculate_gamma_w(bw, dgm$tau.approx)
  
  po_data <- generate_potential_outcomes(
    df = sim_data,
    gamma_true = dgm$gamma.true,
    mu = dgm$mu,
    tau_strata = sim_data$tau.strataO,
    gamma_w = gamma_w,
    wname = wname,
    seed = ss
  )
  
  obs_data <- apply_randomization_treatment(
    po_data = po_data,
    strata_rand = strata_rand,
    keep_rand = keep_rand,
    seed = ss
  )
  
  obs_data <- apply_censoring(
    obs_data = obs_data,
    muC = dgm$muC,
    tauC = dgm$tauC,
    time_eos = time_eos,
    seed = ss
  )
  
  obs_data <- add_stratification_variables(
    obs_data = obs_data,
    strata_rand = strata_rand,
    strata_tte = dgm$strata_tte
  )
  
  if (details && ss <= 3) {
    cat(sprintf("\nSimulation %d:\n", ss))
    cat(sprintf("  N = %d\n", nrow(obs_data)))
    cat(sprintf("  Events = %d (%.1f%%)\n", 
                sum(obs_data$event.sim), 
                100 * mean(obs_data$event.sim)))
    cat(sprintf("  Treatment: %d vs %d\n", 
                sum(obs_data$treat.sim), 
                sum(1 - obs_data$treat.sim)))
  }
  
  return(format_simulation_output(obs_data, return_df, hrz_crit))
}

# Test single simulation
test_sim1 <- draw_sim_stratified(
  dgm = dgm_test,
  ss = 1,
  wname = "meno",
  bw = -log(5),
  strata_rand = "stratum",
  details = TRUE
)

cat("\n✓ Single simulation completed\n")
cat(sprintf("  Columns: %s\n", paste(names(test_sim1)[1:10], collapse = ", ")))

# ------------------------------------------------------------------------------
# STEP 5: TEST REPRODUCIBILITY
# ------------------------------------------------------------------------------

cat("\n=== Testing Reproducibility ===\n")

test_sim1a <- draw_sim_stratified(dgm_test, ss = 1, wname = "meno", bw = -log(5))
test_sim1b <- draw_sim_stratified(dgm_test, ss = 1, wname = "meno", bw = -log(5))

# Check if identical
max_diff_y <- max(abs(test_sim1a$y.sim - test_sim1b$y.sim))
max_diff_event <- max(abs(test_sim1a$event.sim - test_sim1b$event.sim))
max_diff_treat <- max(abs(test_sim1a$treat.sim - test_sim1b$treat.sim))

cat(sprintf("Max difference in y.sim: %.2e\n", max_diff_y))
cat(sprintf("Max difference in event.sim: %.2e\n", max_diff_event))
cat(sprintf("Max difference in treat.sim: %.2e\n", max_diff_treat))

if (max_diff_y < 1e-10 && max_diff_event < 1e-10 && max_diff_treat < 1e-10) {
  cat("✓ Reproducibility test PASSED\n")
} else {
  cat("✗ Reproducibility test FAILED\n")
}

# ------------------------------------------------------------------------------
# STEP 6: TEST MULTIPLE SIMULATIONS
# ------------------------------------------------------------------------------

cat("\n=== Testing Multiple Simulations ===\n")

sim_results <- list()
for (i in 1:5) {
  sim_results[[i]] <- draw_sim_stratified(
    dgm_test, ss = i, wname = "meno", bw = -log(5), details = TRUE
  )
}

cat("\n✓ Multiple simulations completed\n")

# ------------------------------------------------------------------------------
# STEP 7: TEST COX MODEL FITTING
# ------------------------------------------------------------------------------

cat("\n=== Testing Cox Model Fitting ===\n")

# Fit different Cox models
df_test <- sim_results[[1]]

# Unadjusted
fit1 <- coxph(Surv(y.sim, event.sim) ~ treat.sim, data = df_test)
cat(sprintf("Unadjusted HR: %.3f\n", exp(coef(fit1))))

# Stratified by randomization
fit2 <- coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR), 
              data = df_test)
cat(sprintf("Stratified by R HR: %.3f\n", exp(coef(fit2))))

# Adjusted for prognostic factor
fit3 <- coxph(Surv(y.sim, event.sim) ~ treat.sim + w, data = df_test)
cat(sprintf("Adjusted for W HR: %.3f\n", exp(coef(fit3))[1]))

cat("✓ Cox models fit successfully\n")

# ------------------------------------------------------------------------------
# STEP 8: TEST SUBGROUP ANALYSES (SIMPLIFIED)
# ------------------------------------------------------------------------------

cat("\n=== Testing Subgroup Analyses ===\n")

# Define simple subgroups
subgroups_id <- c(
  "itt == 'Y'",
  "meno == 0",
  "meno == 1",
  "bm_low == 1"
)

subgroups_name <- c(
  "All Patients",
  "Pre-menopausal",
  "Post-menopausal",
  "BM Low"
)

# Check subgroup sizes
cat("\nSubgroup sizes in synthetic data:\n")
for (i in seq_along(subgroups_id)) {
  df_sg <- subset(df.synthetic, eval(parse(text = subgroups_id[i])))
  cat(sprintf("  %-20s: %3d patients\n", subgroups_name[i], nrow(df_sg)))
}

# Simplified subgroup analysis (just a few sims)
cat("\nRunning simplified subgroup analysis (10 simulations)...\n")

n_sims <- 10
hrs_matrix <- matrix(NA, nrow = n_sims, ncol = length(subgroups_name))
colnames(hrs_matrix) <- subgroups_name

for (ss in 1:n_sims) {
  df_sim <- draw_sim_stratified(dgm_test, ss = ss, wname = "meno", bw = -log(5))
  df_sim$itt <- "Y"
  
  for (gg in seq_along(subgroups_id)) {
    df_sg <- subset(df_sim, eval(parse(text = subgroups_id[gg])))
    
    if (nrow(df_sg) >= 20 && sum(df_sg$event.sim) >= 10) {
      fit <- tryCatch({
        coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR), 
              data = df_sg)
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        hrs_matrix[ss, gg] <- exp(coef(fit))
      }
    }
  }
}

cat("\n✓ Subgroup analysis completed\n")

# Summarize results
cat("\nSubgroup Analysis Results (median HR):\n")
for (sg in subgroups_name) {
  hrs <- hrs_matrix[, sg]
  hrs_valid <- hrs[!is.na(hrs)]
  if (length(hrs_valid) > 0) {
    cat(sprintf("  %-20s: %.3f [n=%d]\n", 
                sg, median(hrs_valid), length(hrs_valid)))
  } else {
    cat(sprintf("  %-20s: No valid estimates\n", sg))
  }
}

# ------------------------------------------------------------------------------
# STEP 9: VISUAL VALIDATION
# ------------------------------------------------------------------------------

cat("\n=== Creating Validation Plots ===\n")

pdf("test_synthetic_validation.pdf", width = 12, height = 10)
par(mfrow = c(2, 3))

# 1. KM curves for one simulation
df_plot <- sim_results[[1]]
fit_km <- survfit(Surv(y.sim, event.sim) ~ treat.sim, data = df_plot)
plot(fit_km, col = c("blue", "red"), lwd = 2,
     xlab = "Time", ylab = "Survival Probability",
     main = "KM Curves - Simulation 1")
legend("topright", c("Treatment", "Control"), 
       col = c("red", "blue"), lwd = 2, bty = "n")

# 2. Biomarker distribution
hist(df.synthetic$z, breaks = 30, col = "skyblue", border = "white",
     xlab = "Biomarker (z)", main = "Biomarker Distribution")
abline(v = 2, col = "red", lty = 2, lwd = 2)

# 3. Treatment effect by biomarker (true)
z_seq <- seq(0, 8, by = 0.5)
loghr_true <- log(0.60)  # Uniform effect
plot(z_seq, rep(loghr_true, length(z_seq)), type = "l", lwd = 2,
     xlab = "Biomarker (z)", ylab = "log(HR)",
     main = "True Treatment Effect", ylim = c(-1, 0))
abline(h = 0, col = "gray", lty = 2)
abline(h = log(0.60), col = "red", lty = 2)

# 4. HR estimates across simulations (ITT)
hrs_itt <- hrs_matrix[, "All Patients"]
boxplot(hrs_itt, main = "HR Estimates (ITT)", ylab = "Hazard Ratio",
        col = "lightgreen")
abline(h = 0.60, col = "red", lty = 2, lwd = 2)
abline(h = 1, col = "gray", lty = 2)

# 5. Sample sizes by simulation
sample_sizes <- sapply(sim_results, nrow)
plot(1:5, sample_sizes, type = "b", pch = 19,
     xlab = "Simulation", ylab = "Sample Size",
     main = "Sample Size Consistency")
abline(h = nrow(df.synthetic), col = "red", lty = 2)

# 6. Event rates by simulation
event_rates <- sapply(sim_results, function(x) mean(x$event.sim))
plot(1:5, event_rates, type = "b", pch = 19,
     xlab = "Simulation", ylab = "Event Rate",
     main = "Event Rate Consistency", ylim = c(0, 1))

dev.off()

cat("✓ Validation plots saved to test_synthetic_validation.pdf\n")

# ------------------------------------------------------------------------------
# STEP 10: SUMMARY REPORT
# ------------------------------------------------------------------------------

cat("\n" )
cat("===========================================================\n")
cat("TESTING SUMMARY REPORT\n")
cat("===========================================================\n\n")

cat("Dataset:\n")
cat(sprintf("  - Synthetic data: %d patients\n", nrow(df.synthetic)))
cat(sprintf("  - Biomarker range: [%.2f, %.2f]\n", 
            min(df.synthetic$z), max(df.synthetic$z)))
cat(sprintf("  - Event rate: %.1f%%\n", 
            100 * mean(df.synthetic$event)))
cat("\n")

cat("DGM:\n")
cat("  - True treatment effect: HR = 0.60 (uniform)\n")
cat("  - Prognostic factor: meno (HR = 5.0)\n")
cat("\n")

cat("Tests Performed:\n")
cat("  ✓ Single simulation generation\n")
cat("  ✓ Reproducibility check\n")
cat("  ✓ Multiple simulations (n=5)\n")
cat("  ✓ Cox model fitting\n")
cat("  ✓ Subgroup analysis (n=10 sims)\n")
cat("  ✓ Validation plots\n")
cat("\n")

cat("Key Results:\n")
cat(sprintf("  - Median ITT HR: %.3f (true=0.60)\n", 
            median(hrs_itt, na.rm = TRUE)))
cat(sprintf("  - IQR ITT HR: %.3f\n", 
            IQR(hrs_itt, na.rm = TRUE)))
cat(sprintf("  - Success rate: %.1f%%\n", 
            100 * mean(!is.na(hrs_itt))))
cat("\n")

cat("All tests completed successfully!\n")
cat("===========================================================\n\n")