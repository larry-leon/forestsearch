# Comprehensive GBSG Synthetic Data Generation
# This script provides multiple methods to generate mock data mimicking the GBSG dataset

library(survival)
library(MASS)

# Load the original GBSG dataset
data(cancer)

# ================================================================================
# METHOD 1: PARAMETRIC SIMULATION (Using distributions fitted to original data)
# ================================================================================

generate_gbsg_parametric <- function(n = 686, seed = 123) {
  set.seed(seed)
  
  # Extract statistics from original data
  original_stats <- list(
    # Continuous variables statistics
    age_mean = mean(gbsg$age), age_sd = sd(gbsg$age),
    size_mean = mean(gbsg$size), size_sd = sd(gbsg$size),
    nodes_mean = mean(gbsg$nodes), nodes_sd = sd(gbsg$nodes),
    pgr_mean = mean(gbsg$pgr), pgr_sd = sd(gbsg$pgr),
    er_mean = mean(gbsg$er), er_sd = sd(gbsg$er),
    rfstime_mean = mean(gbsg$rfstime), rfstime_sd = sd(gbsg$rfstime),
    
    # Categorical variables proportions
    meno_prop = mean(gbsg$meno),
    grade_probs = prop.table(table(gbsg$grade)),
    hormon_prop = mean(gbsg$hormon),
    status_prop = mean(gbsg$status)
  )
  
  # Generate synthetic data
  synthetic_data <- data.frame(
    # ID variable
    pid = 1:n,
    
    # Age: normal distribution, bounded between 21 and 80
    age = pmin(80, pmax(21, round(rnorm(n, original_stats$age_mean, original_stats$age_sd)))),
    
    # Menopausal status: binary
    meno = rbinom(n, 1, original_stats$meno_prop),
    
    # Tumor size: log-normal distribution to handle skewness
    size = pmin(120, pmax(3, round(exp(rnorm(n, 
                                             log(original_stats$size_mean), 
                                             log(1 + original_stats$size_sd/original_stats$size_mean)))))),
    
    # Grade: categorical (1, 2, 3)
    grade = sample(1:3, n, replace = TRUE, prob = original_stats$grade_probs),
    
    # Number of nodes: negative binomial for count data
    nodes = pmin(51, pmax(1, rnbinom(n, size = 2, mu = original_stats$nodes_mean))),
    
    # Progesterone receptor: highly skewed, use gamma distribution
    pgr = round(pmax(0, rgamma(n, shape = 0.5, scale = original_stats$pgr_mean * 2))),
    
    # Estrogen receptor: similar to pgr
    er = round(pmax(0, rgamma(n, shape = 0.6, scale = original_stats$er_mean * 1.67))),
    
    # Hormone therapy: binary
    hormon = rbinom(n, 1, original_stats$hormon_prop),
    
    # Recurrence-free survival time: Weibull distribution
    rfstime = pmin(2659, pmax(8, round(rweibull(n, shape = 1.5, 
                                                scale = original_stats$rfstime_mean * 1.13)))),
    
    # Status: binary (event indicator)
    status = rbinom(n, 1, original_stats$status_prop)
  )
  
  # Adjust correlations between variables
  # Age tends to correlate with menopausal status
  synthetic_data$meno[synthetic_data$age < 45] <- ifelse(runif(sum(synthetic_data$age < 45)) < 0.2, 1, 0)
  synthetic_data$meno[synthetic_data$age > 55] <- ifelse(runif(sum(synthetic_data$age > 55)) < 0.9, 1, 0)
  
  # Larger tumors tend to have more nodes
  high_size <- synthetic_data$size > median(synthetic_data$size)
  synthetic_data$nodes[high_size] <- pmax(synthetic_data$nodes[high_size], 
                                          rpois(sum(high_size), lambda = 5))
  
  # Event more likely with more nodes and larger size
  risk_score <- scale(synthetic_data$nodes) + scale(synthetic_data$size) - 
                scale(synthetic_data$pgr) - scale(synthetic_data$er)
  event_prob <- plogis(risk_score / 2)
  synthetic_data$status <- rbinom(n, 1, event_prob)
  
  # Adjust survival times based on status
  synthetic_data$rfstime[synthetic_data$status == 1] <- 
    pmin(synthetic_data$rfstime[synthetic_data$status == 1],
         rexp(sum(synthetic_data$status == 1), rate = 1/800))
  
  return(synthetic_data)
}

# ================================================================================
# METHOD 2: COPULA-BASED SIMULATION (Preserves correlation structure)
# ================================================================================

generate_gbsg_copula <- function(n = 686, seed = 123) {
  set.seed(seed)
  
  # Get correlation matrix of continuous variables
  cont_vars <- c("age", "size", "nodes", "pgr", "er", "rfstime")
  cor_mat <- cor(gbsg[, cont_vars], method = "spearman")
  
  # Generate correlated uniform variables using Gaussian copula
  library(MASS)
  normal_vars <- mvrnorm(n, mu = rep(0, length(cont_vars)), Sigma = cor_mat)
  uniform_vars <- pnorm(normal_vars)
  
  # Transform to marginal distributions matching original data
  synthetic_data <- data.frame(
    pid = 1:n,
    age = round(quantile(gbsg$age, uniform_vars[,1])),
    size = round(quantile(gbsg$size, uniform_vars[,2])),
    nodes = pmax(1, round(quantile(gbsg$nodes, uniform_vars[,3]))),
    pgr = pmax(0, round(quantile(gbsg$pgr, uniform_vars[,4]))),
    er = pmax(0, round(quantile(gbsg$er, uniform_vars[,5]))),
    rfstime = pmax(8, round(quantile(gbsg$rfstime, uniform_vars[,6])))
  )
  
  # Add categorical variables with dependencies
  synthetic_data$meno <- ifelse(synthetic_data$age > 50, 
                                rbinom(n, 1, 0.8), 
                                rbinom(n, 1, 0.2))
  
  synthetic_data$grade <- cut(uniform_vars[,2] + runif(n, -0.3, 0.3), 
                              breaks = c(0, 0.12, 0.76, 1), 
                              labels = FALSE)
  
  synthetic_data$hormon <- rbinom(n, 1, mean(gbsg$hormon))
  
  # Status depends on prognostic factors
  linear_predictor <- -0.5 + 0.02 * synthetic_data$size + 
                      0.05 * synthetic_data$nodes - 
                      0.001 * synthetic_data$pgr - 
                      0.001 * synthetic_data$er
  synthetic_data$status <- rbinom(n, 1, plogis(linear_predictor))
  
  # Reorder columns to match original
  synthetic_data <- synthetic_data[, names(gbsg)]
  
  return(synthetic_data)
}

# ================================================================================
# METHOD 3: BOOTSTRAP RESAMPLING WITH PERTURBATION
# ================================================================================

generate_gbsg_bootstrap <- function(n = 686, seed = 123, noise_level = 0.1) {
  set.seed(seed)
  
  # Bootstrap sample from original data
  boot_indices <- sample(nrow(gbsg), n, replace = TRUE)
  synthetic_data <- gbsg[boot_indices, ]
  
  # Add noise to continuous variables
  continuous_vars <- c("age", "size", "nodes", "pgr", "er", "rfstime")
  
  for (var in continuous_vars) {
    # Calculate noise based on variable's standard deviation
    noise_sd <- sd(gbsg[[var]]) * noise_level
    noise <- rnorm(n, 0, noise_sd)
    
    # Add noise and ensure values stay within original bounds
    synthetic_data[[var]] <- synthetic_data[[var]] + noise
    synthetic_data[[var]] <- pmax(min(gbsg[[var]]), 
                                  pmin(max(gbsg[[var]]), 
                                       round(synthetic_data[[var]])))
  }
  
  # Slightly perturb categorical variables
  cat_vars <- c("meno", "grade", "hormon", "status")
  for (var in cat_vars) {
    # Randomly flip a small percentage of values
    flip_prob <- noise_level / 2
    flip_indices <- which(runif(n) < flip_prob)
    
    if (var == "grade") {
      # For grade (1-3), randomly change to adjacent value
      for (i in flip_indices) {
        current <- synthetic_data[[var]][i]
        if (current == 1) synthetic_data[[var]][i] <- 2
        else if (current == 3) synthetic_data[[var]][i] <- 2
        else synthetic_data[[var]][i] <- sample(c(1, 3), 1)
      }
    } else {
      # For binary variables, flip the value
      synthetic_data[[var]][flip_indices] <- 1 - synthetic_data[[var]][flip_indices]
    }
  }
  
  # Reset row names and pid
  rownames(synthetic_data) <- NULL
  synthetic_data$pid <- 1:n
  
  return(synthetic_data)
}

# ================================================================================
# COMPARISON AND VALIDATION FUNCTIONS
# ================================================================================

compare_datasets <- function(original, synthetic, dataset_name = "Synthetic") {
  cat("\n================================================================================\n")
  cat("COMPARISON:", dataset_name, "vs Original GBSG\n")
  cat("================================================================================\n")
  
  # Compare dimensions
  cat("\nDimensions:\n")
  cat("Original:", nrow(original), "x", ncol(original), "\n")
  cat(dataset_name, ":", nrow(synthetic), "x", ncol(synthetic), "\n")
  
  # Compare summary statistics for continuous variables
  cat("\nContinuous Variables Comparison:\n")
  cat("--------------------------------\n")
  continuous_vars <- c("age", "size", "nodes", "pgr", "er", "rfstime")
  
  comparison <- data.frame(
    Variable = continuous_vars,
    Original_Mean = round(sapply(continuous_vars, function(v) mean(original[[v]])), 2),
    Synthetic_Mean = round(sapply(continuous_vars, function(v) mean(synthetic[[v]])), 2),
    Original_SD = round(sapply(continuous_vars, function(v) sd(original[[v]])), 2),
    Synthetic_SD = round(sapply(continuous_vars, function(v) sd(synthetic[[v]])), 2),
    Mean_Diff_Pct = round(100 * abs(sapply(continuous_vars, function(v) 
      (mean(synthetic[[v]]) - mean(original[[v]])) / mean(original[[v]]))), 1)
  )
  print(comparison)
  
  # Compare categorical variables
  cat("\nCategorical Variables Comparison:\n")
  cat("----------------------------------\n")
  cat_vars <- c("meno", "grade", "hormon", "status")
  
  for (var in cat_vars) {
    cat("\n", var, ":\n", sep = "")
    orig_tab <- table(original[[var]])
    synth_tab <- table(synthetic[[var]])
    
    # Ensure same levels
    all_levels <- union(names(orig_tab), names(synth_tab))
    orig_prop <- prop.table(orig_tab)[all_levels]
    synth_prop <- prop.table(synth_tab)[all_levels]
    
    comp_df <- data.frame(
      Level = all_levels,
      Original = round(orig_prop, 3),
      Synthetic = round(synth_prop, 3)
    )
    print(comp_df, row.names = FALSE)
  }
  
  # Correlation comparison
  cat("\nCorrelation Matrix Comparison (key variables):\n")
  cat("-----------------------------------------------\n")
  key_vars <- c("age", "size", "nodes", "rfstime")
  orig_cor <- cor(original[, key_vars])
  synth_cor <- cor(synthetic[, key_vars])
  
  cat("Original Correlations:\n")
  print(round(orig_cor, 3))
  cat("\nSynthetic Correlations:\n")
  print(round(synth_cor, 3))
  cat("\nMaximum Absolute Correlation Difference:", 
      round(max(abs(orig_cor - synth_cor)), 3), "\n")
}

# ================================================================================
# GENERATE AND COMPARE ALL METHODS
# ================================================================================

cat("Generating synthetic GBSG datasets using different methods...\n")

# Generate using each method
synth_parametric <- generate_gbsg_parametric(n = 686, seed = 123)
synth_copula <- generate_gbsg_copula(n = 686, seed = 123)
synth_bootstrap <- generate_gbsg_bootstrap(n = 686, seed = 123, noise_level = 0.1)

# Compare each synthetic dataset with original
compare_datasets(gbsg, synth_parametric, "Parametric Simulation")
compare_datasets(gbsg, synth_copula, "Copula-based Simulation")
compare_datasets(gbsg, synth_bootstrap, "Bootstrap with Perturbation")

# Save the synthetic datasets
write.csv(synth_parametric, "gbsg_synthetic_parametric.csv", row.names = FALSE)
write.csv(synth_copula, "gbsg_synthetic_copula.csv", row.names = FALSE)
write.csv(synth_bootstrap, "gbsg_synthetic_bootstrap.csv", row.names = FALSE)

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")
cat("Three synthetic datasets have been generated and saved:\n")
cat("1. gbsg_synthetic_parametric.csv - Using parametric distributions\n")
cat("2. gbsg_synthetic_copula.csv - Using copula to preserve correlations\n")
cat("3. gbsg_synthetic_bootstrap.csv - Using bootstrap with perturbation\n\n")

cat("Recommendations:\n")
cat("- Parametric method: Good for general simulation, fast, interpretable\n")
cat("- Copula method: Best for preserving correlation structure\n")
cat("- Bootstrap method: Closest to original data, good for small perturbations\n\n")

cat("Note: The 'modgo' package mentioned does not appear to be a standard R package.\n")
cat("These methods provide comprehensive alternatives for synthetic data generation.\n")
