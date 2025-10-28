# ==============================================================================
# COMPLETE EXAMPLE: Using Refactored Simulation Functions
# Demonstrates full workflow from setup to analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------------------

library(survival)
library(data.table)
library(randomizr)
library(ggplot2)

# Source refactored functions
source("sim_weibspline_functions_v2.R")
source("sim_subgroup_analyses_v2.R")
source("sim_subgroup_summarytables_v1.R")  # Keep existing table functions

# Load case study data
df.case <- survival::gbsg

# ------------------------------------------------------------------------------
# STEP 1: PREPARE CASE STUDY DATA
# ------------------------------------------------------------------------------

cat("=== Step 1: Preparing Case Study Data ===\n")

# Transform variables
df.case <- within(df.case, {
  treat <- hormon
  event <- status
  tte <- rfstime / 30.4375  # Convert to months
  z <- log(er + 1)          # Log-transform biomarker
  stratum <- grade          # Stratification for randomization
  bm_low <- ifelse(z <= log(2), 1, 0)  # Binary biomarker indicator
})

# Create additional subgroup variables for analysis
set.seed(8316951)
df.case$random15 <- as.numeric(cumsum(rbinom(nrow(df.case), 1, 0.2)) <= 15)
df.case$random20 <- as.numeric(cumsum(rbinom(nrow(df.case), 1, 0.2)) <= 20)
df.case$country <- sample(c("US", "EU"), nrow(df.case), replace = TRUE, prob = c(0.3, 0.7))

cat(sprintf("Dataset prepared: %d patients\n", nrow(df.case)))
cat(sprintf("Biomarker range: [%.2f, %.2f]\n", min(df.case$z), max(df.case$z)))

# ------------------------------------------------------------------------------
# STEP 2: DEFINE DATA GENERATING MECHANISM (DGM)
# ------------------------------------------------------------------------------

cat("\n=== Step 2: Defining Data Generating Mechanism ===\n")

# Scenario 1: Uniform treatment effect (null hypothesis for biomarker)
cat("\nScenario 1: Uniform HR = 0.60\n")
dgm_uniform <- get_dgm_stratified(
  df = df.case,
  knot = 5,
  zeta = 10,
  log.hrs = log(c(0.60, 0.60, 0.60)),  # Constant effect
  strata_tte = NULL,
  details = FALSE
)

cat("  ✓ DGM created (uniform treatment effect)\n")

# Scenario 2: Biomarker-driven effects
cat("\nScenario 2: Biomarker-driven effects\n")
dgm_biomarker <- get_dgm_stratified(
  df = df.case,
  knot = 5,
  zeta = 10,
  log.hrs = log(c(0.90, 0.70, 0.50)),  # Increasing benefit with biomarker
  strata_tte = NULL,
  details = TRUE  # Show treatment effect profile
)

cat("  ✓ DGM created (biomarker-driven effects)\n")

# ------------------------------------------------------------------------------
# STEP 3: EXPLORE POPULATION-LEVEL EFFECTS
# ------------------------------------------------------------------------------

cat("\n=== Step 3: Exploring Population-Level Effects ===\n")

# Generate large sample to approximate population parameters
pop_summary_uniform <- draw_sim_stratified(
  dgm = dgm_uniform,
  ss = 1,
  Ndraw = 10000,
  wname = "meno",
  bw = -log(5),  # Strong prognostic effect
  strata_rand = "stratum",
  checking = FALSE,
  details = FALSE,
  return_df = FALSE
)

pop_summary_biomarker <- draw_sim_stratified(
  dgm = dgm_biomarker,
  ss = 1,
  Ndraw = 10000,
  wname = "meno",
  bw = -log(5),
  strata_rand = "stratum",
  return_df = FALSE
)

cat("\nPopulation summaries (Uniform HR=0.60):\n")
cat(sprintf("  Overall AHR: %.3f\n", pop_summary_uniform$AHR))
cat(sprintf("  AHR (Post-meno): %.3f\n", pop_summary_uniform$AHR_W1))
cat(sprintf("  AHR (Pre-meno): %.3f\n", pop_summary_uniform$AHR_W0))

cat("\nPopulation summaries (Biomarker-driven):\n")
cat(sprintf("  Overall AHR: %.3f\n", pop_summary_biomarker$AHR))
cat(sprintf("  AHR (Post-meno): %.3f\n", pop_summary_biomarker$AHR_W1))
cat(sprintf("  AHR (Pre-meno): %.3f\n", pop_summary_biomarker$AHR_W0))
cat(sprintf("  Optimal cutpoint: z = %.2f\n", pop_summary_biomarker$cut_optimal))
cat(sprintf("  AHR (Optimal subgroup): %.3f\n", pop_summary_biomarker$AHR_optimal))

# Plot AHR profiles
pdf("output/ahr_profiles.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot_AHRs(pop_summary_uniform, df.case)
title("Uniform Treatment Effect", outer = TRUE, line = -1)

plot_AHRs(pop_summary_biomarker, df.case)
title("Biomarker-Driven Effect", outer = TRUE, line = -1)
dev.off()

cat("  ✓ AHR profiles saved to output/ahr_profiles.pdf\n")

# ------------------------------------------------------------------------------
# STEP 4: GENERATE EXAMPLE SIMULATIONS
# ------------------------------------------------------------------------------

cat("\n=== Step 4: Generating Example Simulations ===\n")

# Generate 3 example simulations
set.seed(123)
example_sims <- list()

for (i in 1:3) {
  example_sims[[i]] <- draw_sim_stratified(
    dgm = dgm_uniform,
    ss = i,
    wname = "meno",
    bw = -log(5),
    strata_rand = "stratum",
    checking = (i == 1),  # Check first simulation
    details = TRUE,
    return_df = TRUE
  )
  
  cat(sprintf("  Simulation %d: N=%d, Events=%d (%.1f%%)\n",
              i, nrow(example_sims[[i]]),
              sum(example_sims[[i]]$event.sim),
              100 * mean(example_sims[[i]]$event.sim)))
}

# Plot KM curves for examples
pdf("output/example_km_curves.pdf", width = 12, height = 12)
par(mfrow = c(2, 2))

for (i in 1:3) {
  df_plot <- example_sims[[i]]
  
  # Simplified KM plot (requires your KM plotting function)
  fit <- survfit(Surv(y.sim, event.sim) ~ treat.sim, data = df_plot)
  plot(fit, col = c("blue", "red"), lwd = 2,
       xlab = "Months", ylab = "Survival Probability",
       main = sprintf("Simulation %d", i))
  legend("topright", c("Treatment", "Control"), 
         col = c("red", "blue"), lwd = 2, bty = "n")
}
dev.off()

cat("  ✓ Example KM curves saved to output/example_km_curves.pdf\n")

# ------------------------------------------------------------------------------
# STEP 5: DEFINE SUBGROUPS FOR ANALYSIS
# ------------------------------------------------------------------------------

cat("\n=== Step 5: Defining Subgroups ===\n")

subgroups_id <- c(
  "itt == 'Y'",
  "bm_low == 0",
  "bm_low == 1",
  "meno == 0",
  "meno == 1",
  "country == 'US'",
  "country == 'EU'",
  "age <= 65",
  "age > 65",
  "size <= median(size)",
  "size > median(size)",
  "grade == 1",
  "grade == 2",
  "grade == 3",
  "random15 == 1",
  "random20 == 1"
)

subgroups_name <- c(
  "All Patients",
  "BM non-low",
  "BM low",
  "Pre-menopausal",
  "Post-menopausal",
  "US",
  "EU",
  "Age ≤65",
  "Age >65",
  "Size ≤median",
  "Size >median",
  "Grade 1",
  "Grade 2",
  "Grade 3",
  "Random n=15",
  "Random n=20"
)

cat(sprintf("  Defined %d subgroups for analysis\n", length(subgroups_name)))

# Examine subgroup sizes in case study data
cat("\nSubgroup sizes in case study data:\n")
for (i in seq_along(subgroups_id)) {
  df_sg <- subset(df.case, eval(parse(text = subgroups_id[i])))
  cat(sprintf("  %-20s: %3d patients\n", subgroups_name[i], nrow(df_sg)))
}

# ------------------------------------------------------------------------------
# STEP 6: RUN SUBGROUP ANALYSES (Small Test)
# ------------------------------------------------------------------------------

cat("\n=== Step 6: Running Subgroup Analyses (Test Run) ===\n")

# Test with small number of simulations first
cat("Running 100 test simulations...\n")

test_results <- get_SGanalyses(
  dgm = dgm_uniform,
  subgroups_name = subgroups_name,
  subgroups_id = subgroups_id,
  sims = 100,
  wname = "meno",
  bw = -log(5),
  Ndraw = nrow(dgm_uniform$df_super),
  bmcut = log(2),
  outfile = "output/test_results.Rdata",
  verbose = TRUE
)

cat("  ✓ Test run complete\n")

# Quick summary of test results
summary_test <- summarize_subgroup_results(test_results)
print(head(summary_test, 10))

# ------------------------------------------------------------------------------
# STEP 7: RUN FULL SIMULATION STUDY
# ------------------------------------------------------------------------------

cat("\n=== Step 7: Running Full Simulation Study ===\n")

# Run with full simulation count
cat("Running 5000 simulations (this may take several minutes)...\n")

full_results <- get_SGanalyses(
  dgm = dgm_uniform,
  subgroups_name = subgroups_name,
  subgroups_id = subgroups_id,
  sims = 5000,
  wname = "meno",
  bw = -log(5),
  bmcut = log(2),
  outfile = "output/full_results_uniform.Rdata",
  verbose = TRUE
)

cat("  ✓ Full simulation complete\n")

# ------------------------------------------------------------------------------
# STEP 8: ANALYZE AND VISUALIZE RESULTS
# ------------------------------------------------------------------------------

cat("\n=== Step 8: Analyzing Results ===\n")

# Generate comprehensive summary
summary_full <- summarize_subgroup_results(full_results)

# Save summary table
write.csv(summary_full, "output/subgroup_summary.csv", row.names = FALSE)
cat("  ✓ Summary table saved to output/subgroup_summary.csv\n")

# Create forest plots for each analysis type
library(forestploter)

# Analysis 1: Stratified by randomization factors (sR)
cat("\nGenerating forest plots...\n")

SG_table_sR <- getSG_dfhrONE(
  z = apply(full_results$ns, 2, mean),
  x1 = full_results$hrs1,
  x2 = full_results$hrs2,
  x3 = full_results$hrs3,
  ubx1 = full_results$ubs1,
  ubx2 = full_results$ubs2,
  ubx3 = full_results$ubs3,
  ubx4 = full_results$ubs4,
  ubx5 = full_results$ubs5,
  ubx6 = full_results$ubs6,
  analysisx1 = "sR",
  which_sgs = subgroups_name,
  alpha = 0.01
)

# Create forest plot
dt <- SG_table_sR
dt$Subgroup <- ifelse(is.na(dt$est), dt$Subgroup, paste0("   ", dt$Subgroup))
dt# ==============================================================================
# COMPLETE EXAMPLE: Using Refactored Simulation Functions
# Demonstrates full workflow from setup to analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------------------

library(survival)
library(data.table)
library(randomizr)
library(ggplot2)

# Source refactored functions
source("sim_weibspline_functions_v2.R")
source("sim_subgroup_analyses_v2.R")
source("sim_subgroup_summarytables_v1.R")  # Keep existing table functions

# Load case study data
df.case <- survival::gbsg

# ------------------------------------------------------------------------------
# STEP 1: PREPARE CASE STUDY DATA
# ------------------------------------------------------------------------------

cat("=== Step 1: Preparing Case Study Data ===\n")

# Transform variables
df.case <- within(df.case, {
  treat <- hormon
  event <- status
  tte <- rfstime / 30.4375  # Convert to months
  z <- log(er + 1)          # Log-transform biomarker
  stratum <- grade          # Stratification for randomization
  bm_low <- ifelse(z <= log(2), 1, 0)  # Binary biomarker indicator
})

# Create additional subgroup variables for analysis
set.seed(8316951)
df.case$random15 <- as.numeric(cumsum(rbinom(nrow(df.case), 1, 0.2)) <= 15)
df.case$random20 <- as.numeric(cumsum(rbinom(nrow(df.case), 1, 0.2)) <= 20)
df.case$country <- sample(c("US", "EU"), nrow(df.case), replace = TRUE, prob = c(0.3, 0.7))

cat(sprintf("Dataset prepared: %d patients\n", nrow(df.case)))
cat(sprintf("Biomarker range: [%.2f, %.2f]\n", min(df.case$z), max(df.case$z)))

# ------------------------------------------------------------------------------
# STEP 2: DEFINE DATA GENERATING MECHANISM (DGM)
# ------------------------------------------------------------------------------

cat("\n=== Step 2: Defining Data Generating Mechanism ===\n")

# Scenario 1: Uniform treatment effect (null hypothesis for biomarker)
cat("\nScenario 1: Uniform HR = 0.60\n")
dgm_uniform <- get_dgm_stratified(
  df = df.case,
  knot = 5,
  zeta = 10,
  log.hrs = log(c(0.60, 0.60, 0.60)),  # Constant effect
  strata_tte = NULL,
  details = FALSE
)

cat("  ✓ DGM created (uniform treatment effect)\n")

# Scenario 2: Biomarker-driven effects
cat("\nScenario 2: Biomarker-driven effects\n")
dgm_biomarker <- get_dgm_stratified(
  df = df.case,
  knot = 5,
  zeta = 10,
  log.hrs = log(c(0.90, 0.70, 0.50)),  # Increasing benefit with biomarker
  strata_tte = NULL,
  details = TRUE  # Show treatment effect profile
)

cat("  ✓ DGM created (biomarker-driven effects)\n")

# ------------------------------------------------------------------------------
# STEP 3: EXPLORE POPULATION-LEVEL EFFECTS
# ------------------------------------------------------------------------------

cat("\n=== Step 3: Exploring Population-Level Effects ===\n")

# Generate large sample to approximate population parameters
pop_summary_uniform <- draw_sim_stratified(
  dgm = dgm_uniform,
  ss = 1,
  Ndraw = 10000,
  wname = "meno",
  bw = -log(5),  # Strong prognostic effect
  strata_rand = "stratum",
  checking = FALSE,
  details = FALSE,
  return_df = FALSE
)

pop_summary_biomarker <- draw_sim_stratified(
  dgm = dgm_biomarker,
  ss = 1,
  Ndraw = 10000,
  wname = "meno",
  bw = -log(5),
  strata_rand = "stratum",
  return_df = FALSE
)

cat("\nPopulation summaries (Uniform HR=0.60):\n")
cat(sprintf("  Overall AHR: %.3f\n", pop_summary_uniform$AHR))
cat(sprintf("  AHR (Post-meno): %.3f\n", pop_summary_uniform$AHR_W1))
cat(sprintf("  AHR (Pre-meno): %.3f\n", pop_summary_uniform$AHR_W0))

cat("\nPopulation summaries (Biomarker-driven):\n")
cat(sprintf("  Overall AHR: %.3f\n", pop_summary_biomarker$AHR))
cat(sprintf("  AHR (Post-meno): %.3f\n", pop_summary_biomarker$AHR_W1))
cat(sprintf("  AHR (Pre-meno): %.3f\n", pop_summary_biomarker$AHR_W0))
cat(sprintf("  Optimal cutpoint: z = %.2f\n", pop_summary_biomarker$cut_optimal))
cat(sprintf("  AHR (Optimal subgroup): %.3f\n", pop_summary_biomarker$AHR_optimal))

# Plot AHR profiles
pdf("output/ahr_profiles.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot_AHRs(pop_summary_uniform, df.case)
title("Uniform Treatment Effect", outer = TRUE, line = -1)

plot_AHRs(pop_summary_biomarker, df.case)
title("Biomarker-Driven Effect", outer = TRUE, line = -1)
dev.off()

cat("  ✓ AHR profiles saved to output/ahr_profiles.pdf\n")

# ------------------------------------------------------------------------------
# STEP 4: GENERATE EXAMPLE SIMULATIONS
# ------------------------------------------------------------------------------

cat("\n=== Step 4: Generating Example Simulations ===\n")

# Generate 3 example simulations
set.seed(123)
example_sims <- list()

for (i in 1:3) {
  example_sims[[i]] <- draw_sim_stratified(
    dgm = dgm_uniform,
    ss = i,
    wname = "meno",
    bw = -log(5),
    strata_rand = "stratum",
    checking = (i == 1),  # Check first simulation
    details = TRUE,
    return_df = TRUE
  )
  
  cat(sprintf("  Simulation %d: N=%d, Events=%d (%.1f%%)\n",
              i, nrow(example_sims[[i]]),
              sum(example_sims[[i]]$event.sim),
              100 * mean(example_sims[[i]]$event.sim)))
}

# Plot KM curves for examples
pdf("output/example_km_curves.pdf", width = 12, height = 12)
par(mfrow = c(2, 2))

for (i in 1:3) {
  df_plot <- example_sims[[i]]
  
  # Simplified KM plot (requires your KM plotting function)
  fit <- survfit(Surv(y.sim, event.sim) ~ treat.sim, data = df_plot)
  plot(fit, col = c("blue", "red"), lwd = 2,
       xlab = "Months", ylab = "Survival Probability",
       main = sprintf("Simulation %d", i))
  legend("topright", c("Treatment", "Control"), 
         col = c("red", "blue"), lwd = 2, bty = "n")
}
dev.off()

cat("  ✓ Example KM curves saved to output/example_km_curves.pdf\n")

# ------------------------------------------------------------------------------
# STEP 5: DEFINE SUBGROUPS FOR ANALYSIS
# ------------------------------------------------------------------------------

cat("\n=== Step 5: Defining Subgroups ===\n")

subgroups_id <- c(
  "itt == 'Y'",
  "bm_low == 0",
  "bm_low == 1",
  "meno == 0",
  "meno == 1",
  "country == 'US'",
  "country == 'EU'",
  "age <= 65",
  "age > 65",
  "size <= median(size)",
  "size > median(size)",
  "grade == 1",
  "grade == 2",
  "grade == 3",
  "random15 == 1",
  "random20 == 1"
)

subgroups_name <- c(
  "All Patients",
  "BM non-low",
  "BM low",
  "Pre-menopausal",
  "Post-menopausal",
  "US",
  "EU",
  "Age ≤65",
  "Age >65",
  "Size ≤median",
  "Size >median",
  "Grade 1",
  "Grade 2",
  "Grade 3",
  "Random n=15",
  "Random n=20"
)

cat(sprintf("  Defined %d subgroups for analysis\n", length(subgroups_name)))

# Examine subgroup sizes in case study data
cat("\nSubgroup sizes in case study data:\n")
for (i in seq_along(subgroups_id)) {
  df_sg <- subset(df.case, eval(parse(text = subgroups_id[i])))
  cat(sprintf("  %-20s: %3d patients\n", subgroups_name[i], nrow(df_sg)))
}

# ------------------------------------------------------------------------------
# STEP 6: RUN SUBGROUP ANALYSES (Small Test)
# ------------------------------------------------------------------------------

cat("\n=== Step 6: Running Subgroup Analyses (Test Run) ===\n")

# Test with small number of simulations first
cat("Running 100 test simulations...\n")

test_results <- get_SGanalyses(
  dgm = dgm_uniform,
  subgroups_name = subgroups_name,
  subgroups_id = subgroups_id,
  sims = 100,
  wname = "meno",
  bw = -log(5),
  Ndraw = nrow(dgm_uniform$df_super),
  bmcut = log(2),
  outfile = "output/test_results.Rdata",
  verbose = TRUE
)

cat("  ✓ Test run complete\n")

# Quick summary of test results
summary_test <- summarize_subgroup_results(test_results)
print(head(summary_test, 10))

# ------------------------------------------------------------------------------
# STEP 7: RUN FULL SIMULATION STUDY
# ------------------------------------------------------------------------------

cat("\n=== Step 7: Running Full Simulation Study ===\n")

# Run with full simulation count
cat("Running 5000 simulations (this may take several minutes)...\n")

full_results <- get_SGanalyses(
  dgm = dgm_uniform,
  subgroups_name = subgroups_name,
  subgroups_id = subgroups_id,
  sims = 5000,
  wname = "meno",
  bw = -log(5),
  bmcut = log(2),
  outfile = "output/full_results_uniform.Rdata",
  verbose = TRUE
)

cat("  ✓ Full simulation complete\n")

# ------------------------------------------------------------------------------
# STEP 8: ANALYZE AND VISUALIZE RESULTS
# ------------------------------------------------------------------------------

 ` <- paste(rep(" ", 20), collapse = " ")
dt# ==============================================================================
# COMPLETE EXAMPLE: Using Refactored Simulation Functions
# Demonstrates full workflow from setup to analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------------------

library(survival)
library(data.table)
library(randomizr)
library(ggplot2)

# Source refactored functions
source("sim_weibspline_functions_v2.R")
source("sim_subgroup_analyses_v2.R")
source("sim_subgroup_summarytables_v1.R")  # Keep existing table functions

# Load case study data
df.case <- survival::gbsg

# ------------------------------------------------------------------------------
# STEP 1: PREPARE CASE STUDY DATA
# ------------------------------------------------------------------------------

cat("=== Step 1: Preparing Case Study Data ===\n")

# Transform variables
df.case <- within(df.case, {
  treat <- hormon
  event <- status
  tte <- rfstime / 30.4375  # Convert to months
  z <- log(er + 1)          # Log-transform biomarker
  stratum <- grade          # Stratification for randomization
  bm_low <- ifelse(z <= log(2), 1, 0)  # Binary biomarker indicator
})

# Create additional subgroup variables for analysis
set.seed(8316951)
df.case$random15 <- as.numeric(cumsum(rbinom(nrow(df.case), 1, 0.2)) <= 15)
df.case$random20 <- as.numeric(cumsum(rbinom(nrow(df.case), 1, 0.2)) <= 20)
df.case$country <- sample(c("US", "EU"), nrow(df.case), replace = TRUE, prob = c(0.3, 0.7))

cat(sprintf("Dataset prepared: %d patients\n", nrow(df.case)))
cat(sprintf("Biomarker range: [%.2f, %.2f]\n", min(df.case$z), max(df.case$z)))

# ------------------------------------------------------------------------------
# STEP 2: DEFINE DATA GENERATING MECHANISM (DGM)
# ------------------------------------------------------------------------------

cat("\n=== Step 2: Defining Data Generating Mechanism ===\n")

# Scenario 1: Uniform treatment effect (null hypothesis for biomarker)
cat("\nScenario 1: Uniform HR = 0.60\n")
dgm_uniform <- get_dgm_stratified(
  df = df.case,
  knot = 5,
  zeta = 10,
  log.hrs = log(c(0.60, 0.60, 0.60)),  # Constant effect
  strata_tte = NULL,
  details = FALSE
)

cat("  ✓ DGM created (uniform treatment effect)\n")

# Scenario 2: Biomarker-driven effects
cat("\nScenario 2: Biomarker-driven effects\n")
dgm_biomarker <- get_dgm_stratified(
  df = df.case,
  knot = 5,
  zeta = 10,
  log.hrs = log(c(0.90, 0.70, 0.50)),  # Increasing benefit with biomarker
  strata_tte = NULL,
  details = TRUE  # Show treatment effect profile
)

cat("  ✓ DGM created (biomarker-driven effects)\n")

# ------------------------------------------------------------------------------
# STEP 3: EXPLORE POPULATION-LEVEL EFFECTS
# ------------------------------------------------------------------------------

cat("\n=== Step 3: Exploring Population-Level Effects ===\n")

# Generate large sample to approximate population parameters
pop_summary_uniform <- draw_sim_stratified(
  dgm = dgm_uniform,
  ss = 1,
  Ndraw = 10000,
  wname = "meno",
  bw = -log(5),  # Strong prognostic effect
  strata_rand = "stratum",
  checking = FALSE,
  details = FALSE,
  return_df = FALSE
)

pop_summary_biomarker <- draw_sim_stratified(
  dgm = dgm_biomarker,
  ss = 1,
  Ndraw = 10000,
  wname = "meno",
  bw = -log(5),
  strata_rand = "stratum",
  return_df = FALSE
)

cat("\nPopulation summaries (Uniform HR=0.60):\n")
cat(sprintf("  Overall AHR: %.3f\n", pop_summary_uniform$AHR))
cat(sprintf("  AHR (Post-meno): %.3f\n", pop_summary_uniform$AHR_W1))
cat(sprintf("  AHR (Pre-meno): %.3f\n", pop_summary_uniform$AHR_W0))

cat("\nPopulation summaries (Biomarker-driven):\n")
cat(sprintf("  Overall AHR: %.3f\n", pop_summary_biomarker$AHR))
cat(sprintf("  AHR (Post-meno): %.3f\n", pop_summary_biomarker$AHR_W1))
cat(sprintf("  AHR (Pre-meno): %.3f\n", pop_summary_biomarker$AHR_W0))
cat(sprintf("  Optimal cutpoint: z = %.2f\n", pop_summary_biomarker$cut_optimal))
cat(sprintf("  AHR (Optimal subgroup): %.3f\n", pop_summary_biomarker$AHR_optimal))

# Plot AHR profiles
pdf("output/ahr_profiles.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot_AHRs(pop_summary_uniform, df.case)
title("Uniform Treatment Effect", outer = TRUE, line = -1)

plot_AHRs(pop_summary_biomarker, df.case)
title("Biomarker-Driven Effect", outer = TRUE, line = -1)
dev.off()

cat("  ✓ AHR profiles saved to output/ahr_profiles.pdf\n")

# ------------------------------------------------------------------------------
# STEP 4: GENERATE EXAMPLE SIMULATIONS
# ------------------------------------------------------------------------------

cat("\n=== Step 4: Generating Example Simulations ===\n")

# Generate 3 example simulations
set.seed(123)
example_sims <- list()

for (i in 1:3) {
  example_sims[[i]] <- draw_sim_stratified(
    dgm = dgm_uniform,
    ss = i,
    wname = "meno",
    bw = -log(5),
    strata_rand = "stratum",
    checking = (i == 1),  # Check first simulation
    details = TRUE,
    return_df = TRUE
  )
  
  cat(sprintf("  Simulation %d: N=%d, Events=%d (%.1f%%)\n",
              i, nrow(example_sims[[i]]),
              sum(example_sims[[i]]$event.sim),
              100 * mean(example_sims[[i]]$event.sim)))
}

# Plot KM curves for examples
pdf("output/example_km_curves.pdf", width = 12, height = 12)
par(mfrow = c(2, 2))

for (i in 1:3) {
  df_plot <- example_sims[[i]]
  
  # Simplified KM plot (requires your KM plotting function)
  fit <- survfit(Surv(y.sim, event.sim) ~ treat.sim, data = df_plot)
  plot(fit, col = c("blue", "red"), lwd = 2,
       xlab = "Months", ylab = "Survival Probability",
       main = sprintf("Simulation %d", i))
  legend("topright", c("Treatment", "Control"), 
         col = c("red", "blue"), lwd = 2, bty = "n")
}
dev.off()

cat("  ✓ Example KM curves saved to output/example_km_curves.pdf\n")

# ------------------------------------------------------------------------------
# STEP 5: DEFINE SUBGROUPS FOR ANALYSIS
# ------------------------------------------------------------------------------

cat("\n=== Step 5: Defining Subgroups ===\n")

subgroups_id <- c(
  "itt == 'Y'",
  "bm_low == 0",
  "bm_low == 1",
  "meno == 0",
  "meno == 1",
  "country == 'US'",
  "country == 'EU'",
  "age <= 65",
  "age > 65",
  "size <= median(size)",
  "size > median(size)",
  "grade == 1",
  "grade == 2",
  "grade == 3",
  "random15 == 1",
  "random20 == 1"
)

subgroups_name <- c(
  "All Patients",
  "BM non-low",
  "BM low",
  "Pre-menopausal",
  "Post-menopausal",
  "US",
  "EU",
  "Age ≤65",
  "Age >65",
  "Size ≤median",
  "Size >median",
  "Grade 1",
  "Grade 2",
  "Grade 3",
  "Random n=15",
  "Random n=20"
)

cat(sprintf("  Defined %d subgroups for analysis\n", length(subgroups_name)))

# Examine subgroup sizes in case study data
cat("\nSubgroup sizes in case study data:\n")
for (i in seq_along(subgroups_id)) {
  df_sg <- subset(df.case, eval(parse(text = subgroups_id[i])))
  cat(sprintf("  %-20s: %3d patients\n", subgroups_name[i], nrow(df_sg)))
}

# ------------------------------------------------------------------------------
# STEP 6: RUN SUBGROUP ANALYSES (Small Test)
# ------------------------------------------------------------------------------

cat("\n=== Step 6: Running Subgroup Analyses (Test Run) ===\n")

# Test with small number of simulations first
cat("Running 100 test simulations...\n")

test_results <- get_SGanalyses(
  dgm = dgm_uniform,
  subgroups_name = subgroups_name,
  subgroups_id = subgroups_id,
  sims = 100,
  wname = "meno",
  bw = -log(5),
  Ndraw = nrow(dgm_uniform$df_super),
  bmcut = log(2),
  outfile = "output/test_results.Rdata",
  verbose = TRUE
)

cat("  ✓ Test run complete\n")

# Quick summary of test results
summary_test <- summarize_subgroup_results(test_results)
print(head(summary_test, 10))

# ------------------------------------------------------------------------------
# STEP 7: RUN FULL SIMULATION STUDY
# ------------------------------------------------------------------------------

cat("\n=== Step 7: Running Full Simulation Study ===\n")

# Run with full simulation count
cat("Running 5000 simulations (this may take several minutes)...\n")

full_results <- get_SGanalyses(
  dgm = dgm_uniform,
  subgroups_name = subgroups_name,
  subgroups_id = subgroups_id,
  sims = 5000,
  wname = "meno",
  bw = -log(5),
  bmcut = log(2),
  outfile = "output/full_results_uniform.Rdata",
  verbose = TRUE
)

cat("  ✓ Full simulation complete\n")

# ------------------------------------------------------------------------------
# STEP 8: ANALYZE AND VISUALIZE RESULTS
# ------------------------------------------------------------------------------

HR (99% ECI)` <- ifelse(
  is.na(dt$se), "",
  sprintf("%.2f (%.2f to %.2f)", dt$est, dt$low, dt$hi)
)

tm <- forest_theme(
  base_size = 10,
  refline_gp = gpar(col = "red"),
  footnote_gp = gpar(col = "#636363", fontface = "italic")
)

pdf("output/forest_plot_sR.pdf", width = 10, height = 12)
p <- forest(
  dt[, c(1:4, 9:10)],
  est = dt$est,
  lower = dt$low,
  upper = dt$hi,
  sizes = dt$se / 2,
  ci_column = 5,
  ref_line = 0.60,
  arrow_lab = c("Favors Treatment", "Favors Control"),
  xlim = c(0.3, 1.2),
  ticks_at = c(0.4, 0.6, 0.8, 1.0),
  title = "Hazard Ratio Estimates (Stratified by Randomization)",
  footnote = "Analysis: Cox model stratified by randomization factors",
  theme = tm
)
plot(p)
dev.off()

cat("  ✓ Forest plot saved to output/forest_plot_sR.pdf\n")

# ------------------------------------------------------------------------------
# STEP 9: EXAMINE EXTREME SUBGROUP RESULTS
# ------------------------------------------------------------------------------

cat("\n=== Step 9: Examining Extreme Results ===\n")

# Extract HR estimates for all subgroups (Analysis 1: stratified by R)
hrs_matrix <- full_results$hrs1

# Find extreme results
extreme_results <- data.frame()

for (i in seq_along(subgroups_name)) {
  hrs <- hrs_matrix[, i]
  hrs_valid <- hrs[!is.na(hrs)]
  
  if (length(hrs_valid) > 0) {
    extreme_results <- rbind(extreme_results, data.frame(
      Subgroup = subgroups_name[i],
      Mean_N = mean(full_results$ns[, i], na.rm = TRUE),
      Median_HR = median(hrs_valid),
      Q025_HR = quantile(hrs_valid, 0.025),
      Q975_HR = quantile(hrs_valid, 0.975),
      Min_HR = min(hrs_valid),
      Max_HR = max(hrs_valid),
      Prop_Extreme_Low = mean(hrs_valid < 0.4),
      Prop_Extreme_High = mean(hrs_valid > 1.0),
      N_Success = length(hrs_valid)
    ))
  }
}

# Print extreme results
cat("\nSubgroups with potential for extreme HR estimates:\n")
extreme_subset <- extreme_results[
  extreme_results$Prop_Extreme_Low > 0.10 | 
  extreme_results$Prop_Extreme_High > 0.10,
]

if (nrow(extreme_subset) > 0) {
  print(extreme_subset, row.names = FALSE)
} else {
  cat("  No subgroups showed >10% extreme results\n")
}

# Save extreme results
write.csv(extreme_results, "output/extreme_results_summary.csv", row.names = FALSE)
cat("\n  ✓ Extreme results saved to output/extreme_results_summary.csv\n")

# ------------------------------------------------------------------------------
# STEP 10: COMPARE ANALYSES WITH DIFFERENT ADJUSTMENTS
# ------------------------------------------------------------------------------

cat("\n=== Step 10: Comparing Analysis Approaches ===\n")

# Focus on ITT population
itt_comparison <- data.frame(
  Analysis = c("Unadjusted", "Strata R", "Strata W", "Strata BM", "Strata W+R", "Adj W+age+Strata R"),
  Median_HR = c(
    median(full_results$hrs2[, "All Patients"], na.rm = TRUE),
    median(full_results$hrs1[, "All Patients"], na.rm = TRUE),
    median(full_results$hrs3[, "All Patients"], na.rm = TRUE),
    median(full_results$hrs4[, "All Patients"], na.rm = TRUE),
    median(full_results$hrs5[, "All Patients"], na.rm = TRUE),
    median(full_results$hrs6[, "All Patients"], na.rm = TRUE)
  ),
  IQR_HR = c(
    IQR(full_results$hrs2[, "All Patients"], na.rm = TRUE),
    IQR(full_results$hrs1[, "All Patients"], na.rm = TRUE),
    IQR(full_results$hrs3[, "All Patients"], na.rm = TRUE),
    IQR(full_results$hrs4[, "All Patients"], na.rm = TRUE),
    IQR(full_results$hrs5[, "All Patients"], na.rm = TRUE),
    IQR(full_results$hrs6[, "All Patients"], na.rm = TRUE)
  )
)

cat("\nITT Analysis Comparison:\n")
print(itt_comparison, row.names = FALSE)

# Visualize comparison
pdf("output/analysis_comparison.pdf", width = 10, height = 6)
par(mfrow = c(1, 1))

boxplot(
  list(
    Unadj = full_results$hrs2[, "All Patients"],
    StrataR = full_results$hrs1[, "All Patients"],
    StrataW = full_results$hrs3[, "All Patients"],
    StrataBM = full_results$hrs4[, "All Patients"],
    StrataWR = full_results$hrs5[, "All Patients"],
    AdjWage = full_results$hrs6[, "All Patients"]
  ),
  main = "ITT HR Estimates: Comparison of Analysis Approaches",
  ylab = "Hazard Ratio",
  xlab = "Analysis Method",
  col = rainbow(6, alpha = 0.5),
  las = 2
)
abline(h = 0.60, col = "red", lty = 2, lwd = 2)
legend("topright", "True HR = 0.60", col = "red", lty = 2, lwd = 2, bty = "n")

dev.off()

cat("  ✓ Analysis comparison saved to output/analysis_comparison.pdf\n")

# ------------------------------------------------------------------------------
# STEP 11: SENSITIVITY ANALYSIS - PARALLEL EXECUTION
# ------------------------------------------------------------------------------

cat("\n=== Step 11: Running Sensitivity Analysis (Parallel) ===\n")

# Test parallel execution with different prognostic effect strengths
if (requireNamespace("parallel", quietly = TRUE)) {
  
  cat("Running parallel simulations with varying prognostic effects...\n")
  
  # Weak prognostic effect (bw = -log(2))
  cat("  - Weak prognostic effect (HR=2.0)...\n")
  results_weak <- get_SGanalyses_parallel(
    dgm = dgm_uniform,
    subgroups_name = c("All Patients", "Post-menopausal"),
    subgroups_id = c("itt == 'Y'", "meno == 1"),
    sims = 1000,
    wname = "meno",
    bw = -log(2),
    n_cores = 4,
    outfile = "output/sensitivity_weak.Rdata"
  )
  
  # Strong prognostic effect (bw = -log(5))
  cat("  - Strong prognostic effect (HR=5.0)...\n")
  results_strong <- get_SGanalyses_parallel(
    dgm = dgm_uniform,
    subgroups_name = c("All Patients", "Post-menopausal"),
    subgroups_id = c("itt == 'Y'", "meno == 1"),
    sims = 1000,
    wname = "meno",
    bw = -log(5),
    n_cores = 4,
    outfile = "output/sensitivity_strong.Rdata"
  )
  
  # Compare
  sensitivity_comparison <- data.frame(
    Scenario = c("Weak (HR=2.0)", "Strong (HR=5.0)"),
    ITT_Median = c(
      median(results_weak$hrs1[, "All Patients"], na.rm = TRUE),
      median(results_strong$hrs1[, "All Patients"], na.rm = TRUE)
    ),
    PostMeno_Median = c(
      median(results_weak$hrs1[, "Post-menopausal"], na.rm = TRUE),
      median(results_strong$hrs1[, "Post-menopausal"], na.rm = TRUE)
    )
  )
  
  cat("\nSensitivity Analysis Results:\n")
  print(sensitivity_comparison, row.names = FALSE)
  
  write.csv(sensitivity_comparison, "output/sensitivity_comparison.csv", row.names = FALSE)
  cat("  ✓ Sensitivity results saved\n")
  
} else {
  cat("  Skipping parallel execution (parallel package not available)\n")
}

# ------------------------------------------------------------------------------
# STEP 12: GENERATE FINAL REPORT
# ------------------------------------------------------------------------------

cat("\n=== Step 12: Generating Final Report ===\n")

# Create summary report
report <- c(
  "=====================================",
  "SIMULATION STUDY SUMMARY REPORT",
  "=====================================",
  "",
  "Study Design:",
  sprintf("  - Base dataset: %d patients", nrow(df.case)),
  sprintf("  - Number of simulations: %d", full_results$metadata$n_sims),
  sprintf("  - Number of subgroups: %d", full_results$metadata$n_subgroups),
  sprintf("  - Number of analyses: %d", full_results$metadata$n_analyses),
  "",
  "True Parameters:",
  "  - Uniform treatment effect: HR = 0.60",
  "  - Prognostic factor effect: HR = 5.0 (meno)",
  "",
  "Key Findings:",
  sprintf("  - ITT median HR: %.3f", median(full_results$hrs1[, "All Patients"], na.rm = TRUE)),
  sprintf("  - ITT IQR: %.3f", IQR(full_results$hrs1[, "All Patients"], na.rm = TRUE)),
  sprintf("  - Smallest subgroup analyzed: %s (n=%.0f)",
          subgroups_name[which.min(apply(full_results$ns, 2, mean))],
          min(apply(full_results$ns, 2, mean))),
  "",
  "Extreme Results:",
  sprintf("  - Subgroups with HR<0.4 in >10%% sims: %d",
          sum(extreme_results$Prop_Extreme_Low > 0.10)),
  sprintf("  - Subgroups with HR>1.0 in >10%% sims: %d",
          sum(extreme_results$Prop_Extreme_High > 0.10)),
  "",
  "Output Files Generated:",
  "  - output/ahr_profiles.pdf",
  "  - output/example_km_curves.pdf",
  "  - output/forest_plot_sR.pdf",
  "  - output/subgroup_summary.csv",
  "  - output/extreme_results_summary.csv",
  "  - output/analysis_comparison.pdf",
  "  - output/full_results_uniform.Rdata",
  "",
  "Analysis Complete!",
  "====================================="
)

# Write report
writeLines(report, "output/SIMULATION_REPORT.txt")

# Print report to console
cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n\n")

cat("✓ All analyses complete! Check the output/ directory for results.\n\n")

# ------------------------------------------------------------------------------
# BONUS: INTERACTIVE EXPLORATION FUNCTION
# ------------------------------------------------------------------------------

#' Interactive function to explore specific subgroup results
#' 
#' @param results Output from get_SGanalyses()
#' @param subgroup Name of subgroup to explore
#' @export
explore_subgroup <- function(results, subgroup) {
  
  if (!subgroup %in% results$metadata$subgroup_names) {
    cat("Available subgroups:\n")
    print(results$metadata$subgroup_names)
    stop("Subgroup not found")
  }
  
  sg_results <- extract_subgroup_results(results, subgroup)
  
  cat(sprintf("\n=== Subgroup: %s ===\n\n", subgroup))
  
  # Sample sizes
  cat("Sample Size:\n")
  cat(sprintf("  Mean: %.1f\n", mean(sg_results$sample_size)))
  cat(sprintf("  Range: [%d, %d]\n", 
              min(sg_results$sample_size), 
              max(sg_results$sample_size)))
  
  # HR estimates for each analysis
  cat("\nHazard Ratio Estimates:\n")
  for (i in seq_along(results$metadata$analysis_names)) {
    analysis_name <- results$metadata$analysis_names[i]
    hr_col <- paste0("hr_", analysis_name)
    hrs <- sg_results[[hr_col]]
    
    cat(sprintf("  %s:\n", analysis_name))
    cat(sprintf("    Median: %.3f [IQR: %.3f]\n", 
                median(hrs, na.rm = TRUE),
                IQR(hrs, na.rm = TRUE)))
    cat(sprintf("    95%% range: [%.3f, %.3f]\n",
                quantile(hrs, 0.025, na.rm = TRUE),
                quantile(hrs, 0.975, na.rm = TRUE)))
  }
  
  # Create visualization
  par(mfrow = c(2, 2))
  
  # Histogram of HRs
  hist(sg_results$hr_sR, breaks = 30, 
       main = paste("HR Distribution:", subgroup),
       xlab = "Hazard Ratio", col = "skyblue", border = "white")
  abline(v = 0.60, col = "red", lwd = 2, lty = 2)
  
  # Sample size over simulations
  plot(sg_results$sample_size, type = "l",
       main = "Sample Size Across Simulations",
       xlab = "Simulation", ylab = "N")
  
  # QQ plot
  qqnorm(sg_results$hr_sR, main = "QQ Plot: HR Estimates")
  qqline(sg_results$hr_sR, col = "red")
  
  # Comparison across analyses
  boxplot(
    list(
      sR = sg_results$hr_sR,
      none = sg_results$hr_none,
      sW = sg_results$hr_sW
    ),
    main = "HR Estimates by Analysis",
    ylab = "Hazard Ratio",
    col = rainbow(3, alpha = 0.5)
  )
  abline(h = 0.60, col = "red", lty = 2)
  
  invisible(sg_results)
}

# Example usage:
# explore_subgroup(full_results, "Post-menopausal")

cat("\n")
cat("=================================================================\n")
cat("WORKFLOW COMPLETE!\n")
cat("=================================================================\n")
cat("\n")
cat("To explore specific subgroups interactively, use:\n")
cat("  explore_subgroup(full_results, 'Subgroup Name')\n")
cat("\n")
cat("All output files are in the output/ directory.\n")
cat("\n")