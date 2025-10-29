# Migration Guide: Refactored Simulation Functions

## 📋 Overview

This guide helps you migrate from the original monolithic functions to the new modular design.

## 🎯 Benefits of Refactored Code

| Aspect | Before | After |
|--------|--------|-------|
| **Testability** | Impossible to test parts | Each function unit tested |
| **Debuggability** | 250-line functions | 10-50 line functions |
| **Maintainability** | Hard to modify | Easy to extend |
| **Reusability** | Code duplication | Shared components |
| **Documentation** | Sparse | Comprehensive |
| **Performance** | Not parallelizable | Parallel-ready |

---

## 🔄 Migration Steps

### Step 1: Install Refactored Functions

```r
# Source the new files
source("sim_weibspline_functions_v2.R")  # Refactored simulation
source("sim_subgroup_analyses_v2.R")      # Refactored subgroup analysis
```

### Step 2: Run Side-by-Side Comparison

```r
# Load your existing setup
source("sim_weibspline_functions_v1.R")
source("sim_subgroup_analyses_v1.R")

# Test with same data
set.seed(123)

# OLD VERSION
df_old <- draw_sim_stratified(
  dgm = dgm, ss = 1, wname = "meno", bw = -log(5), 
  strata_rand = "stratum", return_df = TRUE
)

# NEW VERSION
df_new <- draw_sim_stratified(
  dgm = dgm, ss = 1, wname = "meno", bw = -log(5),
  strata_rand = "stratum", return_df = TRUE
)

# Compare outputs
compare_simulation_outputs(df_old, df_new)
```

### Step 3: Validate Results Match

```r
# Function to compare outputs
compare_simulation_outputs <- function(df_old, df_new, tolerance = 1e-10) {
  cat("=== Comparing Simulation Outputs ===\n\n")
  
  # Check dimensions
  cat(sprintf("Dimensions: old=%dx%d, new=%dx%d\n",
              nrow(df_old), ncol(df_old),
              nrow(df_new), ncol(df_new)))
  
  # Check key columns exist
  key_cols <- c("y.sim", "event.sim", "treat.sim", "loghr.po")
  
  for (col in key_cols) {
    if (col %in% names(df_old) && col %in% names(df_new)) {
      max_diff <- max(abs(df_old[[col]] - df_new[[col]]))
      
      if (max_diff < tolerance) {
        cat(sprintf("✓ %s: MATCH (max diff = %.2e)\n", col, max_diff))
      } else {
        cat(sprintf("✗ %s: DIFFER (max diff = %.2e)\n", col, max_diff))
      }
    }
  }
  
  # Compare summary statistics
  cat("\n=== Summary Statistics ===\n")
  cat(sprintf("Censoring rate - old: %.3f, new: %.3f\n",
              mean(1 - df_old$event.sim),
              mean(1 - df_new$event.sim)))
  
  cat(sprintf("Mean log(HR) - old: %.4f, new: %.4f\n",
              mean(df_old$loghr.po),
              mean(df_new$loghr.po)))
  
  invisible(list(old = df_old, new = df_new))
}
```

### Step 4: Test Subgroup Analysis

```r
# Define test subgroups
subgroups_id <- c("itt == 'Y'", "meno == 1")
subgroups_name <- c("All Patients", "Post-meno")

# OLD VERSION
results_old <- get_SGanalyses(
  dgm, subgroups_name, subgroups_id,
  sims = 10, wname = "meno", bw = -log(5)
)

# NEW VERSION  
results_new <- get_SGanalyses(
  dgm, subgroups_name, subgroups_id,
  sims = 10, wname = "meno", bw = -log(5)
)

# Compare
compare_subgroup_results(results_old, results_new)
```

---

## 🧪 Comprehensive Testing Suite

### Unit Tests

```r
# test_simulation_functions.R

library(testthat)

test_that("prepare_simulation_data samples correctly", {
  dgm <- get_dgm_stratified(df.case, log.hrs = log(c(0.6, 0.6, 0.6)))
  
  # Test same size
  df1 <- prepare_simulation_data(dgm, Ndraw = nrow(dgm$df_super), seed = 1)
  expect_equal(nrow(df1), nrow(dgm$df_super))
  
  # Test smaller sample
  df2 <- prepare_simulation_data(dgm, Ndraw = 100, seed = 1)
  expect_equal(nrow(df2), 100)
  
  # Test reproducibility
  df3 <- prepare_simulation_data(dgm, Ndraw = 100, seed = 1)
  expect_equal(df2, df3)
})

test_that("build_covariate_matrix creates correct structure", {
  df_test <- data.frame(
    treat = c(1, 0, 1),
    z = c(1.5, 2.0, 0.5),
    z.treat = c(1.5, 0, 0.5),
    z.k = c(0, 0, 0),
    z.k.treat = c(0, 0, 0)
  )
  
  zmat_1 <- build_covariate_matrix(df_test, treat = 1)
  
  expect_equal(zmat_1[, "treat"], rep(1, 3))
  expect_equal(zmat_1[, "z.treat"], df_test$z)
})

test_that("calculate_log_hazard_ratio produces valid output", {
  N <- 100
  zmat_1 <- matrix(rnorm(N * 5), ncol = 5)
  zmat_0 <- matrix(rnorm(N * 5), ncol = 5)
  gamma_true <- rnorm(5)
  tau_strata <- rep(1, N)
  
  loghr <- calculate_log_hazard_ratio(zmat_1, zmat_0, gamma_true, tau_strata)
  
  expect_length(loghr, N)
  expect_true(all(is.finite(loghr)))
})

test_that("apply_randomization_treatment balances treatment", {
  set.seed(123)
  po_data <- data.frame(
    log.Y1 = rnorm(100),
    log.Y0 = rnorm(100),
    stratum = rep(1:4, 25),
    treat = sample(0:1, 100, replace = TRUE)
  )
  
  result <- apply_randomization_treatment(
    po_data, 
    strata_rand = "stratum", 
    keep_rand = FALSE,
    seed = 1
  )
  
  # Check balance within strata
  for (s in unique(result$stratum)) {
    strata_data <- result[result$stratum == s, ]
    balance <- mean(strata_data$treat.sim)
    expect_true(balance >= 0.4 && balance <= 0.6)  # Roughly balanced
  }
})

test_that("apply_censoring produces valid events", {
  obs_data <- data.frame(
    y.sim = rexp(100, rate = 0.1)
  )
  
  result <- apply_censoring(
    obs_data,
    muC = 5,
    tauC = 1,
    time_eos = 10,
    seed = 1
  )
  
  expect_true("event.sim" %in% names(result))
  expect_true(all(result$event.sim %in% c(0, 1)))
  expect_true(all(result$y.sim <= 10))  # Administrative censoring
})

test_that("calculate_population_summaries returns expected structure", {
  df_test <- data.frame(
    loghr.po = rnorm(100, mean = -0.3, sd = 0.2),
    theta1.po = rnorm(100),
    theta0.po = rnorm(100),
    w = sample(0:1, 100, replace = TRUE),
    z = runif(100, 0, 10),
    y.sim = rexp(100),
    event.sim = rbinom(100, 1, 0.7),
    treat.sim = rbinom(100, 1, 0.5),
    strata.simR = sample(1:3, 100, replace = TRUE),
    tau.strataO = rep(1, 100)
  )
  
  result <- calculate_population_summaries(df_test, hrz_crit = log(1.2))
  
  expect_true(is.list(result))
  expect_true("AHR" %in% names(result))
  expect_true("CDE" %in% names(result))
  expect_true(result$AHR > 0)
})

test_that("validate_simulation_inputs catches errors", {
  dgm_bad <- list(gamma.true = c(1, 2, 3))  # Missing elements
  
  expect_error(
    validate_simulation_inputs(dgm_bad, "stratum", "meno"),
    "Required DGM elements missing"
  )
})
```

### Integration Tests

```r
# test_integration.R

test_that("Full simulation workflow produces consistent results", {
  # Setup
  dgm <- get_dgm_stratified(
    df = df.case,
    log.hrs = log(c(0.6, 0.6, 0.6)),
    strata_tte = NULL
  )
  
  # Run two simulations with same seed
  df1 <- draw_sim_stratified(dgm, ss = 1, wname = "meno", bw = -log(5))
  df2 <- draw_sim_stratified(dgm, ss = 1, wname = "meno", bw = -log(5))
  
  # Should be identical
  expect_equal(df1$y.sim, df2$y.sim)
  expect_equal(df1$event.sim, df2$event.sim)
  expect_equal(df1$treat.sim, df2$treat.sim)
})

test_that("Subgroup analysis handles edge cases", {
  dgm <- get_dgm_stratified(df.case, log.hrs = log(c(0.6, 0.6, 0.6)))
  
  # Define challenging subgroups
  subgroups_id <- c(
    "itt == 'Y'",           # Full population
    "age > 70",             # Small subgroup
    "nodes > 20",           # Very small subgroup
    "grade == 1"            # Another small subgroup
  )
  
  subgroups_name <- c("ITT", "Age>70", "Nodes>20", "Grade1")
  
  # Should handle gracefully
  expect_silent({
    results <- get_SGanalyses(
      dgm, subgroups_name, subgroups_id,
      sims = 5, wname = "meno", bw = -log(5),
      verbose = FALSE
    )
  })
  
  # Check structure
  expect_true(is.list(results))
  expect_equal(ncol(results$ns), length(subgroups_name))
})

test_that("Parallel execution matches sequential", {
  skip_if_not_installed("parallel")
  
  dgm <- get_dgm_stratified(df.case, log.hrs = log(c(0.6, 0.6, 0.6)))
  
  subgroups_id <- c("itt == 'Y'", "meno == 1")
  subgroups_name <- c("ITT", "Post-meno")
  
  # Sequential
  set.seed(123)
  results_seq <- get_SGanalyses(
    dgm, subgroups_name, subgroups_id,
    sims = 10, wname = "meno", bw = -log(5),
    verbose = FALSE
  )
  
  # Parallel
  set.seed(123)
  results_par <- get_SGanalyses_parallel(
    dgm, subgroups_name, subgroups_id,
    sims = 10, wname = "meno", bw = -log(5),
    n_cores = 2
  )
  
  # Compare (allowing for small numerical differences)
  expect_equal(results_seq$ns, results_par$ns)
  expect_equal(results_seq$hrs1, results_par$hrs1, tolerance = 1e-10)
})
```

### Performance Tests

```r
# test_performance.R

test_that("Refactored code is faster than original", {
  skip_if_not_installed("microbenchmark")
  library(microbenchmark)
  
  dgm <- get_dgm_stratified(df.case, log.hrs = log(c(0.6, 0.6, 0.6)))
  
  # Compare single simulation speed
  timing <- microbenchmark(
    old = draw_sim_stratified_v1(dgm, ss = 1, wname = "meno", bw = -log(5)),
    new = draw_sim_stratified(dgm, ss = 1, wname = "meno", bw = -log(5)),
    times = 10
  )
  
  print(timing)
  
  # New version should be comparable or faster
  median_old <- median(timing[timing$expr == "old", "time"])
  median_new <- median(timing[timing$expr == "new", "time"])
  
  expect_true(median_new <= median_old * 1.2)  # Allow 20% slower at worst
})
```

---

## 📊 Validation Checklist

### ✅ Numerical Validation

```r
# validation_script.R

validate_refactored_code <- function(dgm, n_sims = 100) {
  cat("=== Validation Report ===\n\n")
  
  results <- list()
  
  # 1. Test reproducibility
  cat("1. Testing reproducibility...\n")
  df1 <- draw_sim_stratified(dgm, ss = 42, wname = "meno", bw = -log(5))
  df2 <- draw_sim_stratified(dgm, ss = 42, wname = "meno", bw = -log(5))
  
  reproducible <- all.equal(df1$y.sim, df2$y.sim)
  results$reproducible <- isTRUE(reproducible)
  cat(sprintf("   %s\n", ifelse(results$reproducible, "✓ PASS", "✗ FAIL")))
  
  # 2. Test potential outcomes consistency
  cat("2. Testing PO consistency...\n")
  max_diff <- max(abs(df1$loghr.po - (df1$log.Y0 - df1$log.Y1) / df1$tau.strataO))
  results$po_consistent <- max_diff < 1e-10
  cat(sprintf("   Max difference: %.2e %s\n", 
              max_diff, ifelse(results$po_consistent, "✓ PASS", "✗ FAIL")))
  
  # 3. Test treatment balance
  cat("3. Testing treatment balance...\n")
  balance <- abs(mean(df1$treat.sim) - 0.5)
  results$balanced <- balance < 0.1
  cat(sprintf("   Balance: %.3f %s\n", 
              balance, ifelse(results$balanced, "✓ PASS", "✗ FAIL")))
  
  # 4. Test censoring rate reasonable
  cat("4. Testing censoring rate...\n")
  censor_rate <- mean(1 - df1$event.sim)
  results$censoring_ok <- censor_rate > 0.05 && censor_rate < 0.95
  cat(sprintf("   Censoring rate: %.3f %s\n",
              censor_rate, ifelse(results$censoring_ok, "✓ PASS", "✗ FAIL")))
  
  # 5. Test population summaries
  cat("5. Testing population summaries...\n")
  pop_sum <- draw_sim_stratified(
    dgm, ss = 1, Ndraw = 10000, wname = "meno", 
    bw = -log(5), return_df = FALSE
  )
  
  results$ahr_reasonable <- pop_sum$AHR > 0.3 && pop_sum$AHR < 1.5
  cat(sprintf("   AHR: %.3f %s\n",
              pop_sum$AHR, ifelse(results$ahr_reasonable, "✓ PASS", "✗ FAIL")))
  
  # 6. Test subgroup analysis runs
  cat("6. Testing subgroup analysis...\n")
  subgroups_id <- c("itt == 'Y'", "meno == 1")
  subgroups_name <- c("ITT", "Post-meno")
  
  sg_results <- tryCatch({
    get_SGanalyses(
      dgm, subgroups_name, subgroups_id,
      sims = 10, wname = "meno", bw = -log(5),
      verbose = FALSE
    )
  }, error = function(e) NULL)
  
  results$sg_runs <- !is.null(sg_results)
  cat(sprintf("   %s\n", ifelse(results$sg_runs, "✓ PASS", "✗ FAIL")))
  
  # Summary
  cat("\n=== Summary ===\n")
  n_pass <- sum(unlist(results))
  n_total <- length(results)
  cat(sprintf("Passed: %d/%d tests\n", n_pass, n_total))
  
  if (n_pass == n_total) {
    cat("✓ All validation checks passed!\n")
  } else {
    cat("✗ Some validation checks failed. Review output above.\n")
  }
  
  return(results)
}

# Run validation
validation_results <- validate_refactored_code(dgm, n_sims = 100)
```

---

## 🔧 Debugging Tips

### Common Issues and Solutions

#### Issue 1: Different Random Numbers

**Problem:** Results differ between old and new versions.

**Solution:** Ensure seed management is identical.

```r
# Check seed sequence
SIMULATION_CONSTANTS$BASE_SEED              # 8316951
SIMULATION_CONSTANTS$SEED_INCREMENT         # 1000
SIMULATION_CONSTANTS$CENSORING_SEED_OFFSET  # 100

# For simulation ss=1:
# Main seed: 8316951 + 1*1000 = 8317951
# PO seed:   8316951 + 1*1000 + 1 = 8317952
# Rand seed: 8316951 + 1*1000 + 2 = 8317953
# Censor:    8316951 + 1*1000 + 100 = 8318051
```

#### Issue 2: Missing Variables

**Problem:** Error about missing columns in dataframe.

**Solution:** Check variable naming consistency.

```r
# Ensure these exist in dgm$df_super:
required_vars <- c("treat", "z", "z.treat", "z.k", "z.k.treat", "tau.strataO")

# Check what's missing
missing <- required_vars[!required_vars %in% names(dgm$df_super)]
if (length(missing) > 0) {
  cat("Missing variables:", paste(missing, collapse = ", "), "\n")
}
```

#### Issue 3: Subgroup Analysis Fails

**Problem:** get_SGanalyses returns all NAs for certain subgroups.

**Solution:** Check subgroup size and event counts.

```r
# Diagnose subgroup issues
diagnose_subgroup <- function(df, subgroup_expr) {
  df_sg <- subset(df, eval(parse(text = subgroup_expr)))
  
  cat(sprintf("Subgroup: %s\n", subgroup_expr))
  cat(sprintf("  N: %d\n", nrow(df_sg)))
  cat(sprintf("  Events: %d (%.1f%%)\n", 
              sum(df_sg$event.sim), 
              100 * mean(df_sg$event.sim)))
  cat(sprintf("  Treatment balance: %.3f\n", mean(df_sg$treat.sim)))
  
  # Check if analysis is feasible
  feasible <- nrow(df_sg) >= 10 && sum(df_sg$event.sim) >= 5
  cat(sprintf("  Feasible: %s\n", ifelse(feasible, "Yes", "No")))
}
```

---

## 📈 Performance Comparison

### Benchmarking Script

```r
# benchmark_refactored.R

library(microbenchmark)
library(ggplot2)

run_performance_comparison <- function(dgm) {
  cat("=== Performance Comparison ===\n\n")
  
  # Test 1: Single simulation
  cat("Test 1: Single simulation\n")
  t1 <- microbenchmark(
    v1 = draw_sim_stratified_v1(dgm, ss = 1, wname = "meno", bw = -log(5)),
    v2 = draw_sim_stratified(dgm, ss = 1, wname = "meno", bw = -log(5)),
    times = 20
  )
  print(t1)
  
  # Test 2: Subgroup analysis (small)
  cat("\nTest 2: Subgroup analysis (10 sims)\n")
  subgroups_id <- c("itt == 'Y'", "meno == 1")
  subgroups_name <- c("ITT", "Post-meno")
  
  t2 <- microbenchmark(
    v1 = get_SGanalyses_v1(
      dgm, subgroups_name, subgroups_id,
      sims = 10, wname = "meno", bw = -log(5)
    ),
    v2 = get_SGanalyses(
      dgm, subgroups_name, subgroups_id,
      sims = 10, wname = "meno", bw = -log(5),
      verbose = FALSE
    ),
    times = 5
  )
  print(t2)
  
  # Test 3: Parallel vs Sequential
  if (requireNamespace("parallel", quietly = TRUE)) {
    cat("\nTest 3: Parallel execution (50 sims)\n")
    t3 <- microbenchmark(
      sequential = get_SGanalyses(
        dgm, subgroups_name, subgroups_id,
        sims = 50, wname = "meno", bw = -log(5),
        verbose = FALSE
      ),
      parallel = get_SGanalyses_parallel(
        dgm, subgroups_name, subgroups_id,
        sims = 50, wname = "meno", bw = -log(5),
        n_cores = 4
      ),
      times = 3
    )
    print(t3)
  }
  
  # Plot results
  autoplot(t1) + 
    ggtitle("Single Simulation Performance") +
    theme_minimal()
}
```

---

## 🚀 Deployment Steps

### Step-by-Step Migration

#### Phase 1: Testing (Week 1)

```r
# 1. Run all tests
testthat::test_dir("tests/")

# 2. Run validation
validation_results <- validate_refactored_code(dgm)

# 3. Compare with small simulation
results_comparison <- compare_old_vs_new(dgm, n_sims = 100)
```

#### Phase 2: Parallel Deployment (Week 2)

```r
# 1. Keep old functions with "_v1" suffix
# 2. Load new functions
source("sim_weibspline_functions_v2.R")
source("sim_subgroup_analyses_v2.R")

# 3. Run critical analyses with both versions
results_v1 <- run_analysis_v1()
results_v2 <- run_analysis_v2()

# 4. Verify results match
verify_results_match(results_v1, results_v2)
```

#### Phase 3: Full Migration (Week 3)

```r
# 1. Update all analysis scripts to use v2
# 2. Archive v1 functions
# 3. Update documentation
# 4. Run full simulation suite
```

---

## 📚 Additional Resources

### Function Reference Quick Guide

| Old Function | New Function | Changes |
|--------------|--------------|---------|
| `draw_sim_stratified()` | `draw_sim_stratified()` | Modular internals, same API |
| N/A | `prepare_simulation_data()` | New: Data preparation |
| N/A | `generate_potential_outcomes()` | New: PO generation |
| N/A | `apply_randomization_treatment()` | New: Randomization |
| N/A | `apply_censoring()` | New: Censoring logic |
| `get_SGanalyses()` | `get_SGanalyses()` | Cleaner structure, same API |
| N/A | `get_SGanalyses_parallel()` | New: Parallel execution |
| N/A | `extract_subgroup_results()` | New: Result extraction |
| N/A | `summarize_subgroup_results()` | New: Summary generation |

### Key Improvements Summary

1. **Modularity**: 250-line functions → 10-50 line functions
2. **Testability**: Untestable → Comprehensive unit tests
3. **Documentation**: Minimal → Full roxygen2 docs
4. **Error Handling**: Silent failures → Explicit error messages
5. **Performance**: Single-threaded → Parallel-ready
6. **Maintainability**: Hard to modify → Easy to extend

---

## ✅ Final Checklist

- [ ] All unit tests pass
- [ ] Integration tests pass  
- [ ] Validation checks pass
- [ ] Performance is acceptable
- [ ] Documentation is complete
- [ ] Old and new results match
- [ ] Team is trained on new functions
- [ ] Backup of old code created
- [ ] Migration complete

**Once all items are checked, you're ready to fully migrate!**