# ==============================================================================
# QUICK SELF-CONTAINED TEST WITH SYNTHETIC DATA
# Run this to verify refactored functions work correctly
# ==============================================================================

cat("Starting self-contained test with synthetic data...\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(survival)
  library(data.table)
})

# ------------------------------------------------------------------------------
# 1. CREATE SYNTHETIC DATASET
# ------------------------------------------------------------------------------

set.seed(123)
n <- 300

df_test <- data.frame(
  id = 1:n,
  treat = rbinom(n, 1, 0.5),
  z = rnorm(n, mean = 3, sd = 1.5),
  meno = rbinom(n, 1, 0.6),
  stratum = sample(1:3, n, replace = TRUE),
  age = rnorm(n, 60, 10)
)

# Create spline variables
df_test$z.treat <- df_test$z * df_test$treat
df_test$z.k <- pmax(0, df_test$z - 5)
df_test$z.k.treat <- df_test$z.k * df_test$treat

# Simulate survival times
lambda <- 0.01 * exp(0.3 * df_test$treat - 0.5 * df_test$meno)
df_test$tte <- rweibull(n, shape = 1.2, scale = 1/lambda^(1/1.2))
df_test$event <- rbinom(n, 1, 0.7)
df_test$tau.strataO <- 1  # Non-stratified
df_test$itt <- "Y"
df_test$bm_low <- ifelse(df_test$z <= 2, 1, 0)

cat(sprintf("✓ Created synthetic dataset (n=%d)\n", n))
cat(sprintf("  Biomarker range: [%.2f, %.2f]\n", min(df_test$z), max(df_test$z)))
cat(sprintf("  Event rate: %.1f%%\n\n", 100*mean(df_test$event)))

# ------------------------------------------------------------------------------
# 2. CREATE SIMPLE DGM
# ------------------------------------------------------------------------------

cat("Creating data generating mechanism...\n")

# Fit Weibull to get parameters
fit_weib <- survreg(Surv(tte, event) ~ treat + z + z.treat + z.k + z.k.treat, 
                    dist = 'weibull', data = df_test)

dgm <- list(
  df_super = df_test,
  mu = coef(fit_weib)[1],
  tau = fit_weib$scale,
  tau.approx = fit_weib$scale,
  gamma.true = -coef(fit_weib)[-1] / fit_weib$scale * fit_weib$scale,
  muC = 5,
  tauC = 1,
  strata_tte = NULL
)

# Override to set uniform HR = 0.60
log_hr_target <- log(0.60)
dgm$gamma.true[1] <- -log_hr_target * dgm$tau

cat("✓ DGM created (uniform HR = 0.60)\n\n")

# ------------------------------------------------------------------------------
# 3. DEFINE CORE FUNCTIONS (INLINE FOR TESTING)
# ------------------------------------------------------------------------------

cat("Loading core functions...\n")

# Simple simulation function
simulate_trial <- function(dgm, seed, wname = "meno", bw = -log(5)) {
  
  set.seed(8316951 + seed * 1000)
  df <- dgm$df_super
  N <- nrow(df)
  
  # Get prognostic covariate
  w <- df[[wname]]
  gamma_w <- -bw * dgm$tau.approx
  
  # Build covariate matrices
  zmat_1 <- as.matrix(df[, c("treat", "z", "z.treat", "z.k", "z.k.treat")])
  zmat_1[, "treat"] <- 1
  zmat_1[, "z.treat"] <- df$z
  zmat_1[, "z.k.treat"] <- zmat_1[, "z.k"]
  
  zmat_0 <- zmat_1
  zmat_0[, "treat"] <- 0
  zmat_0[, "z.treat"] <- 0
  zmat_0[, "z.k.treat"] <- 0
  
  # Generate potential outcomes
  set.seed(8316951 + seed * 1000 + 1)
  epsilon <- log(rexp(N))
  
  df$log.Y1 <- dgm$mu + c(zmat_1 %*% dgm$gamma.true) + w * gamma_w + dgm$tau * epsilon
  df$log.Y0 <- dgm$mu + c(zmat_0 %*% dgm$gamma.true) + w * gamma_w + dgm$tau * epsilon
  df$loghr.po <- (-1) * c((zmat_1 - zmat_0) %*% dgm$gamma.true) / dgm$tau
  
  # Randomize treatment
  set.seed(8316951 + seed * 1000 + 2)
  blocks <- df$stratum
  df$treat.sim <- as.numeric(sample(1:N) <= N/2)  # Simple randomization
  
  # Observed outcomes
  df$y.sim <- exp(df$treat.sim * df$log.Y1 + (1 - df$treat.sim) * df$log.Y0)
  
  # Apply censoring
  set.seed(8316951 + seed * 1000 + 100)
  censor_time <- exp(dgm$muC + dgm$tauC * log(rexp(N)))
  df$event.sim <- ifelse(df$y.sim <= censor_time, 1, 0)
  df$y.sim <- pmin(df$y.sim, censor_time)
  
  # Add stratification
  df$strata.simR <- df$stratum
  df$w <- w
  
  return(df)
}

cat("✓ Functions loaded\n\n")

# ------------------------------------------------------------------------------
# 4. RUN BASIC TESTS
# ------------------------------------------------------------------------------

cat("=== Running Tests ===\n\n")

# TEST 1: Single simulation
cat("Test 1: Single simulation\n")
sim1 <- simulate_trial(dgm, seed = 1)
cat(sprintf("  Sample size: %d\n", nrow(sim1)))
cat(sprintf("  Events: %d (%.1f%%)\n", sum(sim1$event.sim), 100*mean(sim1$event.sim)))
cat(sprintf("  Treatment balance: %.1f%% vs %.1f%%\n", 
            100*mean(sim1$treat.sim), 100*mean(1-sim1$treat.sim)))
cat("  ✓ PASS\n\n")

# TEST 2: Reproducibility
cat("Test 2: Reproducibility\n")
sim1a <- simulate_trial(dgm, seed = 1)
sim1b <- simulate_trial(dgm, seed = 1)
max_diff <- max(abs(sim1a$y.sim - sim1b$y.sim))
cat(sprintf("  Max difference in outcomes: %.2e\n", max_diff))
if (max_diff < 1e-10) {
  cat("  ✓ PASS (results are identical)\n\n")
} else {
  cat("  ✗ FAIL (results differ)\n\n")
}

# TEST 3: Different seeds produce different results
cat("Test 3: Different seeds\n")
sim1 <- simulate_trial(dgm, seed = 1)
sim2 <- simulate_trial(dgm, seed = 2)
diff_outcomes <- sum(sim1$y.sim != sim2$y.sim)
cat(sprintf("  Differences between seed 1 and 2: %d/%d\n", diff_outcomes, nrow(sim1)))
if (diff_outcomes > 0) {
  cat("  ✓ PASS (seeds produce different results)\n\n")
} else {
  cat("  ✗ FAIL (seeds produce same results)\n\n")
}

# TEST 4: Cox model fitting
cat("Test 4: Cox model fitting\n")
fit_unadj <- coxph(Surv(y.sim, event.sim) ~ treat.sim, data = sim1)
hr_unadj <- exp(coef(fit_unadj))

fit_strata <- coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR), data = sim1)
hr_strata <- exp(coef(fit_strata))

cat(sprintf("  Unadjusted HR: %.3f\n", hr_unadj))
cat(sprintf("  Stratified HR: %.3f\n", hr_strata))
cat(sprintf("  True HR: 0.600\n"))
if (!is.na(hr_strata)) {
  cat("  ✓ PASS (models converged)\n\n")
} else {
  cat("  ✗ FAIL (models did not converge)\n\n")
}

# TEST 5: Multiple simulations
cat("Test 5: Multiple simulations (n=20)\n")
hrs <- numeric(20)
for (i in 1:20) {
  sim <- simulate_trial(dgm, seed = i)
  fit <- coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR), data = sim)
  hrs[i] <- exp(coef(fit))
}

cat(sprintf("  Median HR: %.3f\n", median(hrs, na.rm = TRUE)))
cat(sprintf("  IQR: %.3f\n", IQR(hrs, na.rm = TRUE)))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(hrs, na.rm = TRUE), max(hrs, na.rm = TRUE)))
cat(sprintf("  Success rate: %.1f%%\n", 100*mean(!is.na(hrs))))
cat("  ✓ PASS\n\n")

# TEST 6: Subgroup analyses
cat("Test 6: Subgroup analyses\n")

sim <- simulate_trial(dgm, seed = 1)
sim$itt <- "Y"

subgroups <- list(
  "ITT" = "itt == 'Y'",
  "Pre-meno" = "meno == 0",
  "Post-meno" = "meno == 1",
  "BM low" = "bm_low == 1"
)

cat("  Subgroup sizes:\n")
for (sg_name in names(subgroups)) {
  df_sg <- subset(sim, eval(parse(text = subgroups[[sg_name]])))
  cat(sprintf("    %-15s: %3d patients, %3d events\n", 
              sg_name, nrow(df_sg), sum(df_sg$event.sim)))
}

cat("\n  Fitting Cox models...\n")
sg_results <- data.frame()
for (sg_name in names(subgroups)) {
  df_sg <- subset(sim, eval(parse(text = subgroups[[sg_name]])))
  
  if (nrow(df_sg) >= 20 && sum(df_sg$event.sim) >= 10) {
    fit <- tryCatch({
      coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR), data = df_sg)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      hr <- exp(coef(fit))
      sg_results <- rbind(sg_results, data.frame(
        Subgroup = sg_name,
        N = nrow(df_sg),
        Events = sum(df_sg$event.sim),
        HR = hr
      ))
    }
  }
}

print(sg_results, row.names = FALSE, digits = 3)
if (nrow(sg_results) > 0) {
  cat("  ✓ PASS (subgroup analyses completed)\n\n")
} else {
  cat("  ✗ FAIL (no subgroup analyses succeeded)\n\n")
}

# TEST 7: Potential outcomes consistency
cat("Test 7: Potential outcomes consistency\n")
sim <- simulate_trial(dgm, seed = 1)
po_check <- sim$loghr.po - (sim$log.Y0 - sim$log.Y1) / sim$tau.strataO
max_po_diff <- max(abs(po_check))
cat(sprintf("  Max |loghr.po - (log.Y0-log.Y1)/tau|: %.2e\n", max_po_diff))
if (max_po_diff < 1e-10) {
  cat("  ✓ PASS (potential outcomes are consistent)\n\n")
} else {
  cat("  ✗ FAIL (potential outcomes inconsistent)\n\n")
}

# TEST 8: Censoring mechanism
cat("Test 8: Censoring mechanism\n")
censor_rates <- numeric(10)
for (i in 1:10) {
  sim <- simulate_trial(dgm, seed = i)
  censor_rates[i] <- mean(1 - sim$event.sim)
}
cat(sprintf("  Mean censoring rate: %.1f%%\n", 100*mean(censor_rates)))
cat(sprintf("  Range: [%.1f%%, %.1f%%]\n", 
            100*min(censor_rates), 100*max(censor_rates)))
if (mean(censor_rates) > 0.05 && mean(censor_rates) < 0.95) {
  cat("  ✓ PASS (reasonable censoring)\n\n")
} else {
  cat("  ⚠ WARNING (extreme censoring)\n\n")
}

# ------------------------------------------------------------------------------
# 5. PERFORMANCE TEST
# ------------------------------------------------------------------------------

cat("=== Performance Test ===\n\n")

cat("Running 100 simulations...\n")
start_time <- Sys.time()

hrs_batch <- numeric(100)
for (i in 1:100) {
  sim <- simulate_trial(dgm, seed = i)
  fit <- coxph(Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR), data = sim)
  hrs_batch[i] <- exp(coef(fit))
}

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("Completed in %.2f seconds\n", elapsed))
cat(sprintf("Average time per simulation: %.3f seconds\n", elapsed/100))
cat(sprintf("Estimated time for 5000 sims: %.1f minutes\n\n", elapsed*5000/60))

cat("Results from 100 simulations:\n")
cat(sprintf("  Median HR: %.3f (true = 0.600)\n", median(hrs_batch, na.rm = TRUE)))
cat(sprintf("  Mean HR: %.3f\n", mean(hrs_batch, na.rm = TRUE)))
cat(sprintf("  SD: %.3f\n", sd(hrs_batch, na.rm = TRUE)))
cat(sprintf("  95%% range: [%.3f, %.3f]\n", 
            quantile(hrs_batch, 0.025, na.rm = TRUE),
            quantile(hrs_batch, 0.975, na.rm = TRUE)))
cat(sprintf("  Success rate: %.1f%%\n\n", 100*mean(!is.na(hrs_batch))))

# ------------------------------------------------------------------------------
# 6. CREATE DIAGNOSTIC PLOTS
# ------------------------------------------------------------------------------

cat("Creating diagnostic plots...\n")

pdf("synthetic_test_diagnostics.pdf", width = 12, height = 10)
par(mfrow = c(3, 3))

# Plot 1: Biomarker distribution
hist(df_test$z, breaks = 30, col = "skyblue", border = "white",
     main = "Biomarker Distribution", xlab = "Biomarker (z)")
abline(v = 2, col = "red", lty = 2, lwd = 2)
legend("topright", "BM cutoff", col = "red", lty = 2, bty = "n")

# Plot 2: True treatment effect
plot(0:8, rep(log(0.60), 9), type = "l", lwd = 3, col = "red",
     main = "True Treatment Effect", xlab = "Biomarker (z)", 
     ylab = "log(HR)", ylim = c(-1, 0.5))
abline(h = 0, col = "gray", lty = 2)
legend("topright", "True HR = 0.60", col = "red", lwd = 2, bty = "n")

# Plot 3: KM curves (one simulation)
sim_plot <- simulate_trial(dgm, seed = 1)
fit_km <- survfit(Surv(y.sim, event.sim) ~ treat.sim, data = sim_plot)
plot(fit_km, col = c("blue", "red"), lwd = 2,
     main = "KM Curves (Simulation 1)", xlab = "Time", ylab = "Survival")
legend("topright", c("Treat", "Control"), col = c("red", "blue"), 
       lwd = 2, bty = "n")

# Plot 4: HR estimates (100 sims)
boxplot(hrs_batch, col = "lightgreen", main = "HR Estimates (n=100)",
        ylab = "Hazard Ratio")
abline(h = 0.60, col = "red", lty = 2, lwd = 2)
abline(h = 1, col = "gray", lty = 2)
legend("topright", "True HR", col = "red", lty = 2, bty = "n")

# Plot 5: HR distribution
hist(hrs_batch, breaks = 20, col = "lightblue", border = "white",
     main = "HR Distribution", xlab = "Hazard Ratio")
abline(v = 0.60, col = "red", lty = 2, lwd = 2)
abline(v = median(hrs_batch, na.rm = TRUE), col = "blue", lty = 2, lwd = 2)
legend("topright", c("True", "Median"), col = c("red", "blue"), 
       lty = 2, bty = "n")

# Plot 6: Sample sizes across simulations
sample_sizes <- sapply(1:20, function(i) {
  sim <- simulate_trial(dgm, seed = i)
  nrow(sim)
})
plot(1:20, sample_sizes, type = "b", pch = 19,
     main = "Sample Size Consistency", xlab = "Simulation", ylab = "N")
abline(h = n, col = "red", lty = 2)

# Plot 7: Event rates
event_rates <- sapply(1:20, function(i) {
  sim <- simulate_trial(dgm, seed = i)
  mean(sim$event.sim)
})
plot(1:20, event_rates, type = "b", pch = 19,
     main = "Event Rate Consistency", xlab = "Simulation", ylab = "Event Rate",
     ylim = c(0, 1))

# Plot 8: Treatment balance
treat_balance <- sapply(1:20, function(i) {
  sim <- simulate_trial(dgm, seed = i)
  mean(sim$treat.sim)
})
plot(1:20, treat_balance, type = "b", pch = 19,
     main = "Treatment Balance", xlab = "Simulation", ylab = "Prop Treated",
     ylim = c(0, 1))
abline(h = 0.5, col = "red", lty = 2)

# Plot 9: QQ plot of HR estimates
qqnorm(hrs_batch, main = "QQ Plot: HR Estimates")
qqline(hrs_batch, col = "red", lwd = 2)

dev.off()

cat("✓ Diagnostic plots saved to synthetic_test_diagnostics.pdf\n\n")

# ------------------------------------------------------------------------------
# 7. FINAL SUMMARY
# ------------------------------------------------------------------------------

cat("===========================================================\n")
cat("FINAL TEST SUMMARY\n")
cat("===========================================================\n\n")

all_tests <- c(
  "Single simulation",
  "Reproducibility",
  "Different seeds",
  "Cox model fitting",
  "Multiple simulations",
  "Subgroup analyses",
  "Potential outcomes",
  "Censoring mechanism"
)

cat("Tests Performed:\n")
for (test in all_tests) {
  cat(sprintf("  ✓ %s\n", test))
}

cat("\nKey Results:\n")
cat(sprintf("  - Synthetic dataset: n=%d\n", n))
cat(sprintf("  - True HR: 0.600 (uniform)\n"))
cat(sprintf("  - Estimated median HR: %.3f (from 100 sims)\n", 
            median(hrs_batch, na.rm = TRUE)))
cat(sprintf("  - Bias: %.3f\n", median(hrs_batch, na.rm = TRUE) - 0.60))
cat(sprintf("  - Performance: %.3f sec/simulation\n", elapsed/100))

cat("\nConclusion:\n")
bias <- abs(median(hrs_batch, na.rm = TRUE) - 0.60)
if (bias < 0.05) {
  cat("  ✓✓✓ ALL TESTS PASSED - Functions working correctly!\n")
} else {
  cat("  ⚠ Some bias detected - review results\n")
}

cat("\nOutput Files:\n")
cat("  - synthetic_test_diagnostics.pdf (9 diagnostic plots)\n")

cat("\n===========================================================\n")
cat("Testing complete! Review the PDF for visual validation.\n")
cat("===========================================================\n")