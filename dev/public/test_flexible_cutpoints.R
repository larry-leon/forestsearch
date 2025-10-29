# Demonstration of Flexible Subgroup Cutpoint Specifications
# Shows various ways to define subgroups using the enhanced generate_aft_dgm_flex function

library(survival)
source("generate_aft_dgm_flexible.R")

# Load GBSG data
data(cancer)

cat("================================================================================\n")
cat("DEMONSTRATION: Flexible Subgroup Cutpoint Specifications\n")
cat("================================================================================\n\n")

# ================================================================================
# Example 1: Quantile-based cutpoints
# ================================================================================

cat("Example 1: Quantile-based cutpoints\n")
cat("------------------------------------\n")

dgm1 <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "pgr"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),   # er <= 25th percentile
    pgr = list(type = "quantile", value = 0.50)   # pgr <= median (50th percentile)
  ),
  model = "alt",
  n_super = 1000,
  verbose = TRUE
)

cat("\n")

# ================================================================================
# Example 2: Function-based cutpoints
# ================================================================================

cat("\nExample 2: Function-based cutpoints\n")
cat("------------------------------------\n")

dgm2 <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("age", "size"),
  subgroup_cuts = list(
    age = list(type = "function", fun = median),   # age <= median
    size = list(type = "function", fun = mean)     # size <= mean
  ),
  model = "alt",
  n_super = 500,
  verbose = TRUE
)

cat("\n")

# ================================================================================
# Example 3: Mixed specifications
# ================================================================================

cat("\nExample 3: Mixed specifications (fixed, quantile, range)\n")
cat("---------------------------------------------------------\n")

dgm3 <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "age", "meno"),
  subgroup_cuts = list(
    er = 20,                                      # Fixed value: er <= 20
    age = list(type = "range", min = 40, max = 60),  # Range: 40 <= age <= 60
    meno = 0                                      # Factor: meno == 0 (pre-menopausal)
  ),
  model = "alt",
  n_super = 500,
  verbose = TRUE
)

cat("\n")

# ================================================================================
# Example 4: Custom function for complex criteria
# ================================================================================

cat("\nExample 4: Custom function for complex criteria\n")
cat("------------------------------------------------\n")

# Define a custom function that creates a complex subgroup
custom_er_function <- function(x) {
  # Complex criteria: in the bottom 30% OR above 90th percentile
  q30 <- quantile(x, 0.30, na.rm = TRUE)
  q90 <- quantile(x, 0.90, na.rm = TRUE)
  return(x <= q30 | x >= q90)
}

dgm4 <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("er", "pgr"),
  factor_vars = c("meno"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er"),
  subgroup_cuts = list(
    er = list(type = "custom", fun = custom_er_function)
  ),
  model = "alt",
  n_super = 500,
  verbose = TRUE
)

cat("\n")

# ================================================================================
# Example 5: Finding quantile for target proportion
# ================================================================================

cat("\nExample 5: Finding quantile for target subgroup proportion\n")
cat("-----------------------------------------------------------\n")

# Find the ER quantile that gives us exactly 12.5% of patients in subgroup
target_result <- find_quantile_for_proportion(
  data = gbsg,
  var_name = "er",
  target_prop = 0.125,
  direction = "less"
)

cat("To achieve 12.5% in subgroup (er <= cutpoint):\n")
cat("  Required quantile:", round(target_result$quantile, 3), "\n")
cat("  Cutpoint value:", round(target_result$cutpoint, 2), "\n")

# Use this in the model
dgm5 <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("er", "pgr"),
  factor_vars = c("meno"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = target_result$quantile),
    meno = 0  # pre-menopausal
  ),
  model = "alt",
  n_super = 500,
  verbose = TRUE
)

cat("\n")

# ================================================================================
# Example 6: Greater than cutpoints
# ================================================================================

cat("\nExample 6: Greater than cutpoints\n")
cat("----------------------------------\n")

dgm6 <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("nodes", "size"),
  factor_vars = c("grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("nodes", "grade"),
  subgroup_cuts = list(
    nodes = list(type = "greater", quantile = 0.75),  # nodes > 75th percentile
    grade = 3                                          # grade == 3 (poor)
  ),
  model = "alt",
  n_super = 500,
  verbose = TRUE
)

# ================================================================================
# Simulate from one of the DGMs
# ================================================================================

cat("\n================================================================================\n")
cat("Simulating data from Example 1 DGM\n")
cat("================================================================================\n")

sim_data <- simulate_from_dgm(
  dgm = dgm1,
  n = 700,
  rand_ratio = 1,
  max_follow = 84,
  seed = 123
)

cat("\nSimulated data summary:\n")
cat("  Sample size:", nrow(sim_data), "\n")
cat("  Event rate:", round(mean(sim_data$event_sim), 3), "\n")
cat("  Subgroup size:", sum(sim_data$flag_harm), "\n")
cat("  Subgroup proportion:", round(mean(sim_data$flag_harm), 3), "\n")

# Calculate observed hazard ratios
hr_overall <- exp(coxph(Surv(y_sim, event_sim) ~ treat, data = sim_data)$coefficients)
cat("\nObserved overall HR:", round(hr_overall, 3), "\n")

if (sum(sim_data$flag_harm) > 0) {
  hr_subgroup <- exp(coxph(Surv(y_sim, event_sim) ~ treat, 
                          data = subset(sim_data, flag_harm == 1))$coefficients)
  cat("Observed subgroup HR:", round(hr_subgroup, 3), "\n")
}

cat("\n================================================================================\n")
cat("Summary of Flexible Cutpoint Options\n")
cat("================================================================================\n\n")

cat("1. FIXED VALUE:\n")
cat("   subgroup_cuts = list(er = 20)\n\n")

cat("2. QUANTILE:\n")
cat("   subgroup_cuts = list(er = list(type = 'quantile', value = 0.25))\n\n")

cat("3. FUNCTION (median, mean, etc.):\n")
cat("   subgroup_cuts = list(er = list(type = 'function', fun = median))\n\n")

cat("4. RANGE:\n")
cat("   subgroup_cuts = list(age = list(type = 'range', min = 40, max = 60))\n\n")

cat("5. GREATER THAN:\n")
cat("   subgroup_cuts = list(nodes = list(type = 'greater', quantile = 0.75))\n\n")

cat("6. CUSTOM FUNCTION:\n")
cat("   subgroup_cuts = list(er = list(type = 'custom', fun = your_function))\n\n")

cat("7. MULTIPLE VALUES (for categorical):\n")
cat("   subgroup_cuts = list(grade = list(type = 'multiple', values = c(2, 3)))\n\n")

cat("This flexibility allows you to:\n")
cat("- Define subgroups based on quantiles (as in your original request)\n")
cat("- Achieve target subgroup proportions\n")
cat("- Create complex subgroup definitions\n")
cat("- Mix different cutpoint types for different variables\n")
