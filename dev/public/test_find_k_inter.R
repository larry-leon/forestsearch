# Test Script: Finding k_inter for Target Harm Subgroup HR
# Demonstrates efficient methods to calibrate the interaction effect

library(survival)

# Source the required functions
source("generate_aft_dgm_flexible.R")
source("find_k_inter_efficient.R")

# Load data
data(cancer)

cat("================================================================================\n")
cat("FINDING k_inter FOR TARGET HARM SUBGROUP HAZARD RATIO\n")
cat("================================================================================\n\n")

# ================================================================================
# Setup: Define the base model parameters
# ================================================================================

base_params <- list(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),
    meno = 0
  ),
  k_treat = 0.9,
  n_super = 2000  # Using smaller for faster demonstration
)

# ================================================================================
# Test 1: Find k_inter for HR_harm = 2.0 (harmful interaction)
# ================================================================================

cat("Test 1: Finding k_inter for harm subgroup HR = 2.0\n")
cat("---------------------------------------------------\n")

result1 <- do.call(find_k_inter_for_target_hr, c(
  list(
    target_hr_harm = 2.0,
    k_inter_range = c(-5, 5),
    tol = 0.001,
    verbose = TRUE
  ),
  base_params
))

# ================================================================================
# Test 2: Find k_inter for HR_harm = 0.5 (protective interaction)  
# ================================================================================

cat("\n\nTest 2: Finding k_inter for harm subgroup HR = 0.5\n")
cat("---------------------------------------------------\n")

result2 <- do.call(find_k_inter_for_target_hr, c(
  list(
    target_hr_harm = 0.5,
    k_inter_range = c(-5, 5),
    tol = 0.001,
    verbose = TRUE
  ),
  base_params
))

# ================================================================================
# Test 3: Process multiple target HRs
# ================================================================================

cat("\n\nTest 3: Batch processing multiple target HRs\n")
cat("--------------------------------------------\n")

target_hrs <- c(0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
batch_results <- do.call(find_k_inter_batch, c(
  list(target_hrs = target_hrs),
  base_params
))

# ================================================================================
# Test 4: Sensitivity analysis
# ================================================================================

cat("\n\nTest 4: Sensitivity Analysis\n")
cat("----------------------------\n")
cat("Analyzing how k_inter affects all hazard ratios...\n")

sensitivity_results <- do.call(sensitivity_analysis_k_inter, c(
  list(
    k_inter_range = c(-3, 3),
    n_points = 11,
    model = "alt"
  ),
  base_params
))

cat("\nSensitivity results:\n")
print(round(sensitivity_results, 3))

# ================================================================================
# Practical Example: Setting up scenarios
# ================================================================================

cat("\n\n================================================================================\n")
cat("PRACTICAL EXAMPLE: Setting Up Multiple Scenarios\n")
cat("================================================================================\n\n")

scenarios <- data.frame(
  scenario = c("Strong Benefit", "Moderate Benefit", "No Effect", 
               "Moderate Harm", "Strong Harm"),
  target_hr_harm = c(0.5, 0.75, 1.0, 1.5, 2.0)
)

cat("Scenarios to calibrate:\n")
print(scenarios)
cat("\n")

# Find k_inter for each scenario
scenarios$k_inter <- NA
scenarios$achieved_hr <- NA

for (i in 1:nrow(scenarios)) {
  cat("Processing:", scenarios$scenario[i], "\n")
  
  result <- do.call(find_k_inter_for_target_hr, c(
    list(
      target_hr_harm = scenarios$target_hr_harm[i],
      verbose = FALSE
    ),
    base_params
  ))
  
  if (!is.null(result)) {
    scenarios$k_inter[i] <- result$k_inter
    scenarios$achieved_hr[i] <- result$achieved_hr_harm
  }
}

cat("\nCalibrated scenarios:\n")
print(round(scenarios, 4))

# ================================================================================
# Verification: Generate and test one scenario
# ================================================================================

cat("\n\n================================================================================\n")
cat("VERIFICATION: Testing 'Strong Harm' Scenario\n")
cat("================================================================================\n\n")

# Use the k_inter for "Strong Harm" scenario
k_inter_strong_harm <- scenarios$k_inter[scenarios$scenario == "Strong Harm"]

cat("Using k_inter =", round(k_inter_strong_harm, 4), "for HR_harm = 2.0\n\n")

# Generate DGM with this k_inter
dgm_verify <- do.call(generate_aft_dgm_flex, c(
  list(
    k_inter = k_inter_strong_harm,
    model = "alt",
    verbose = TRUE
  ),
  base_params
))

# Simulate data and check
sim_data <- simulate_from_dgm(dgm_verify, n = 1000, seed = 123)

cat("\n\nSimulated data check:\n")
cat("Event rate:", round(mean(sim_data$event_sim), 3), "\n")
cat("Subgroup size:", sum(sim_data$flag_harm), "out of", nrow(sim_data), "\n")

# Calculate empirical HRs
if (sum(sim_data$flag_harm) > 5 & sum(sim_data$event_sim) > 10) {
  hr_emp_overall <- exp(coxph(Surv(y_sim, event_sim) ~ treat, 
                              data = sim_data)$coefficients)
  hr_emp_harm <- exp(coxph(Surv(y_sim, event_sim) ~ treat, 
                           data = subset(sim_data, flag_harm == 1))$coefficients)
  hr_emp_no_harm <- exp(coxph(Surv(y_sim, event_sim) ~ treat, 
                              data = subset(sim_data, flag_harm == 0))$coefficients)
  
  cat("\nEmpirical HRs from simulated data:\n")
  cat("  Overall:", round(hr_emp_overall, 3), "\n")
  cat("  Harm subgroup:", round(hr_emp_harm, 3), "\n")
  cat("  No-harm subgroup:", round(hr_emp_no_harm, 3), "\n")
  cat("\nNote: These may differ from true HRs due to sampling variability\n")
}

# ================================================================================
# Summary
# ================================================================================

cat("\n\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

cat("Key findings:\n")
cat("1. The root-finding method (uniroot) is most efficient for single targets\n")
cat("2. Grid search is more robust when the relationship is non-monotonic\n")
cat("3. The relationship between k_inter and HR_harm is approximately linear\n")
cat("   in the log scale for moderate values\n")
cat("4. Extreme k_inter values (|k_inter| > 5) may lead to numerical issues\n")
cat("\nRecommendations:\n")
cat("- Use find_k_inter_for_target_hr() for single target HRs\n")
cat("- Use find_k_inter_batch() for multiple scenarios\n")
cat("- Start with k_inter_range = c(-5, 5) and adjust if needed\n")
cat("- Set n_super >= 5000 for stable HR estimates\n")
