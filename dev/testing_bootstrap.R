# ============================================================================
# FORESTSEARCH BOOTSTRAP TEST WITH SYNTHETIC DATA
# ============================================================================

library(survival)
library(data.table)
library(foreach)
library(doFuture)

# Set seed for reproducibility
set.seed(8316951)

# ============================================================================
# STEP 1: Generate Synthetic Survival Data
# ============================================================================

generate_synthetic_survival <- function(n = 300, p_confounders = 6) {
  # Treatment assignment (0 = control, 1 = treatment)
  treat <- rbinom(n, 1, 0.5)

  # Generate confounders
  z1 <- rbinom(n, 1, 0.4)  # Binary confounder 1
  z2 <- rbinom(n, 1, 0.3)  # Binary confounder 2
  z3 <- rnorm(n, 50, 10)   # Continuous confounder (age)
  z4 <- rnorm(n, 25, 5)    # Continuous confounder (BMI)
  z5 <- rbinom(n, 1, 0.5)  # Binary confounder 3
  z6 <- rnorm(n, 100, 15)  # Continuous confounder

  # Define subgroup based on z1 and z3
  # Subgroup H: z1=1 AND z3 > 50 (harm group - questionable benefit)
  subgroup_H <- (z1 == 1) & (z3 > 50)

  # Generate log hazard ratios
  # In H: treatment has minimal/no benefit (HR ~ 1.0)
  # In H^c: treatment has benefit (HR ~ 0.6)
  log_hr <- ifelse(subgroup_H,
                   rnorm(n, log(0.95), 0.1),  # Near null effect in H
                   rnorm(n, log(0.60), 0.15)) # Benefit in H^c

  # Generate survival times using Weibull distribution
  # Baseline hazard influenced by confounders
  baseline_hazard <- exp(0.01 * z3 + 0.02 * z4 - 0.3 * z5)

  # Scale parameter for Weibull
  lambda <- baseline_hazard * exp(-log_hr * treat)
  shape <- 1.5

  # Generate event times
  tte <- rweibull(n, shape = shape, scale = 1/lambda)

  # Generate censoring times (administrative censoring)
  censor_time <- runif(n, 12, 36)

  # Observed time and event indicator
  time_obs <- pmin(tte, censor_time)
  event <- as.numeric(tte <= censor_time)

  # Create data frame
  df <- data.frame(
    id = 1:n,
    tte = time_obs,
    event = event,
    treat = treat,
    z1 = z1,
    z2 = z2,
    z3 = z3,
    z4 = z4,
    z5 = z5,
    z6 = z6,
    true_subgroup = as.numeric(subgroup_H),
    true_loghr = log_hr
  )

  return(df)
}

cat("Generating synthetic survival data...\n")
df_synthetic <- generate_synthetic_survival(n = 300, p_confounders = 6)

cat("\nData Summary:\n")
cat("Sample size:", nrow(df_synthetic), "\n")
cat("Events:", sum(df_synthetic$event),
    sprintf("(%.1f%%)\n", 100*mean(df_synthetic$event)))
cat("Treatment arm:", sum(df_synthetic$treat),
    sprintf("(%.1f%%)\n", 100*mean(df_synthetic$treat)))
cat("True subgroup H size:", sum(df_synthetic$true_subgroup),
    sprintf("(%.1f%%)\n", 100*mean(df_synthetic$true_subgroup)))

# ============================================================================
# STEP 2: Run ForestSearch to Identify Subgroups
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("RUNNING FORESTSEARCH\n")
cat(strrep("=", 70) %+% "\n\n")

# Source required ForestSearch functions
# (In practice, this would be: library(forestsearch))

fs_result <- forestsearch(
  df.analysis = df_synthetic,
  outcome.name = "tte",
  event.name = "event",
  treat.name = "treat",
  id.name = "id",
  confounders.name = c("z1", "z2", "z3", "z4", "z5", "z6"),
  use_lasso = FALSE,
  use_grf = TRUE,
  grf_depth = 2,
  dmin.grf = 20,
  n.min = 40,
  hr.threshold = 1.10,
  hr.consistency = 1.0,
  maxk = 2,
  fs.splits = 500,
  pconsistency.threshold = 0.70,
  max.minutes = 2,
  details = TRUE,
  plot.grf = FALSE,
  plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1)
)

# Check if subgroup was found
if (is.null(fs_result$sg.harm)) {
  cat("\n*** NO SUBGROUP FOUND - Cannot proceed with bootstrap ***\n")
  stop("ForestSearch did not identify a subgroup")
}

cat("\n*** SUBGROUP FOUND ***\n")
cat("Subgroup definition:", fs_result$sg.harm, "\n")
cat("Sample size in H:", sum(fs_result$df.est$treat.recommend == 0), "\n")
cat("Sample size in H^c:", sum(fs_result$df.est$treat.recommend == 1), "\n")

# ============================================================================
# STEP 3: Run Bootstrap Analysis
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("RUNNING BOOTSTRAP ANALYSIS\n")
cat(strrep("=", 70) %+% "\n\n")

# Bootstrap configuration
nb_boots <- 100  # Use 100 for testing (500-1000 for production)

cat("Bootstrap configuration:\n")
cat("  Number of bootstrap samples:", nb_boots, "\n")
cat("  Parallel plan: multisession\n")
cat("  Workers: 2\n\n")

# Run bootstrap
bootstrap_result <- forestsearch_bootstrap_dofuture(
  fs.est = fs_result,
  nb_boots = nb_boots,
  details = TRUE,
  show_three = TRUE,  # Show details for first 3 iterations
  parallel_args = list(
    plan = "multisession",
    workers = 2,
    show_message = TRUE
  )
)

# ============================================================================
# STEP 4: Examine Bootstrap Results
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("BOOTSTRAP RESULTS\n")
cat(strrep("=", 70) %+% "\n\n")

# Bootstrap success rate
boot_results <- bootstrap_result$results
success_rate <- sum(!is.na(boot_results$H_biasadj_2)) / nb_boots

cat("Bootstrap Success Rate:",
    sprintf("%.1f%% (%d/%d)\n",
            100*success_rate,
            sum(!is.na(boot_results$H_biasadj_2)),
            nb_boots))

# Summary of search times
cat("\nSearch Time Summary (minutes):\n")
cat("  Mean:", sprintf("%.2f\n", mean(boot_results$tmins_search, na.rm=TRUE)))
cat("  Median:", sprintf("%.2f\n", median(boot_results$tmins_search, na.rm=TRUE)))
cat("  Range:", sprintf("[%.2f, %.2f]\n",
                        min(boot_results$tmins_search, na.rm=TRUE),
                        max(boot_results$tmins_search, na.rm=TRUE)))

# Subgroup estimates
cat("\n" %+% strrep("-", 70) %+% "\n")
cat("SUBGROUP H (Questionable) ESTIMATES\n")
cat(strrep("-", 70) %+% "\n")
cat("Unadjusted:         ", bootstrap_result$SG_CIs$H_raw, "\n")
cat("Bias-corrected:     ", bootstrap_result$SG_CIs$H_bc, "\n")

cat("\n" %+% strrep("-", 70) %+% "\n")
cat("SUBGROUP H^c (Recommend) ESTIMATES\n")
cat(strrep("-", 70) %+% "\n")
cat("Unadjusted:         ", bootstrap_result$SG_CIs$Hc_raw, "\n")
cat("Bias-corrected:     ", bootstrap_result$SG_CIs$Hc_bc, "\n")

# Summary table
cat("\n" %+% strrep("-", 70) %+% "\n")
cat("SUMMARY TABLE\n")
cat(strrep("-", 70) %+% "\n")
print(bootstrap_result$FSsg_tab)

# ============================================================================
# STEP 5: Diagnostic Plots and Validation
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("DIAGNOSTIC CHECKS\n")
cat(strrep("=", 70) %+% "\n\n")

# Check Ystar matrix dimensions
cat("Ystar matrix check:\n")
cat("  Dimensions:", dim(bootstrap_result$Ystar_mat), "\n")
cat("  Expected: [", nb_boots, "x", nrow(fs_result$df.est), "]\n")
cat("  Match:",
    ifelse(nrow(bootstrap_result$Ystar_mat) == nb_boots &&
             ncol(bootstrap_result$Ystar_mat) == nrow(fs_result$df.est),
           "✓ PASS", "✗ FAIL"), "\n\n")

# Check for NA patterns
cat("Missing data patterns:\n")
cat("  H_biasadj_1 NAs:", sum(is.na(boot_results$H_biasadj_1)), "\n")
cat("  H_biasadj_2 NAs:", sum(is.na(boot_results$H_biasadj_2)), "\n")
cat("  Hc_biasadj_1 NAs:", sum(is.na(boot_results$Hc_biasadj_1)), "\n")
cat("  Hc_biasadj_2 NAs:", sum(is.na(boot_results$Hc_biasadj_2)), "\n\n")

# Distribution of estimates
cat("Distribution of bias-corrected estimates (H):\n")
cat("  Mean:", sprintf("%.3f\n", mean(boot_results$H_biasadj_2, na.rm=TRUE)))
cat("  SD:", sprintf("%.3f\n", sd(boot_results$H_biasadj_2, na.rm=TRUE)))
cat("  Range:", sprintf("[%.3f, %.3f]\n",
                        min(boot_results$H_biasadj_2, na.rm=TRUE),
                        max(boot_results$H_biasadj_2, na.rm=TRUE)))

# Visualizations
if (require(ggplot2, quietly = TRUE)) {

  # 1. Distribution of bias-corrected estimates
  p1 <- ggplot(boot_results, aes(x = H_biasadj_2)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = bootstrap_result$H_estimates$H2[1],
               color = "red", linetype = "dashed", size = 1) +
    labs(title = "Bootstrap Distribution: Subgroup H",
         subtitle = paste("Based on", nb_boots, "bootstrap samples"),
         x = "Bias-Corrected Log(HR)",
         y = "Count") +
    theme_minimal()

  print(p1)

  # 2. Search time vs success
  p2 <- ggplot(boot_results, aes(x = tmins_search,
                                 y = as.numeric(!is.na(H_biasadj_2)))) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = "Bootstrap Success vs Search Time",
         x = "Search Time (minutes)",
         y = "Success (1=found, 0=not found)") +
    theme_minimal()

  print(p2)
}

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("BOOTSTRAP TEST COMPLETE\n")
cat(strrep("=", 70) %+% "\n")

# Return summary
list(
  data_summary = list(
    n = nrow(df_synthetic),
    events = sum(df_synthetic$event),
    true_H_size = sum(df_synthetic$true_subgroup)
  ),
  fs_result = list(
    sg_found = !is.null(fs_result$sg.harm),
    sg_definition = fs_result$sg.harm,
    H_size = sum(fs_result$df.est$treat.recommend == 0),
    Hc_size = sum(fs_result$df.est$treat.recommend == 1)
  ),
  bootstrap_summary = list(
    nb_boots = nb_boots,
    success_rate = success_rate,
    mean_search_time = mean(boot_results$tmins_search, na.rm=TRUE),
    H_estimate = bootstrap_result$SG_CIs$H_bc,
    Hc_estimate = bootstrap_result$SG_CIs$Hc_bc
  )
)
