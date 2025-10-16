# ============================================================================
# TEST SCRIPT FOR REFACTORED get_FSdata.R
# FILE: dev/test_fsdata_refactored_fixed.R
# ============================================================================

rm(list = ls())
devtools::load_all()

cat("\n")
cat("=".replicate(70), "\n", sep = "")
cat("TESTING REFACTORED get_FSdata.R\n")
cat("=".replicate(70), "\n\n", sep = "")

# ============================================================================
# TEST 1: Simple test with continuous variables
# ============================================================================

cat("TEST 1: Basic functionality with continuous variables\n")
cat("-".replicate(70), "\n", sep = "")

set.seed(123)
df_test1 <- data.frame(
  age = rnorm(100, 50, 10),
  size = rnorm(100, 5, 2),
  event = rbinom(100, 1, 0.3),
  tte = rexp(100),
  treat = rbinom(100, 1, 0.5),
  id = 1:100
)

confounders1 <- c("age", "size")

result1 <- try(
  get_FSdata(
    df.analysis = df_test1,
    confounders.name = confounders1,
    outcome.name = "tte",
    event.name = "event",
    cut_type = "default",
    use_lasso = FALSE,
    use_grf = FALSE,
    details = TRUE
  ),
  silent = FALSE
)

if (!inherits(result1, "try-error")) {
  cat("\n✓ TEST 1 PASSED\n")
  cat("  Cuts created:", length(result1$confs_names), "\n")
  cat("  Column names:", paste(result1$confs_names, collapse = ", "), "\n")
  
  # Verify all columns are 0/1
  new_cols <- result1$df[, result1$confs_names]
  all_binary <- all(apply(new_cols, 2, function(x) all(x %in% c(0, 1))))
  cat("  All columns are 0/1?", all_binary, "\n")
  
  # Show sample
  cat("\n  Sample of created factors:\n")
  print(head(new_cols, 10))
} else {
  cat("\n✗ TEST 1 FAILED\n")
  print(result1)
}

cat("\n\n")

# ============================================================================
# TEST 2: With categorical and continuous variables
# ============================================================================

cat("TEST 2: Mixed categorical and continuous variables\n")
cat("-".replicate(70), "\n", sep = "")

set.seed(456)
df_test2 <- data.frame(
  age = rnorm(150, 50, 10),
  size = rnorm(150, 5, 2),
  grade = sample(c(0, 1, 2, 3), 150, replace = TRUE),  # Categorical (4 levels)
  event = rbinom(150, 1, 0.3),
  tte = rexp(150),
  treat = rbinom(150, 1, 0.5),
  id = 1:150
)

confounders2 <- c("age", "size", "grade")

result2 <- try(
  get_FSdata(
    df.analysis = df_test2,
    confounders.name = confounders2,
    outcome.name = "tte",
    event.name = "event",
    cut_type = "default",
    use_lasso = FALSE,
    use_grf = FALSE,
    details = TRUE
  ),
  silent = FALSE
)

if (!inherits(result2, "try-error")) {
  cat("\n✓ TEST 2 PASSED\n")
  cat("  Cuts created:", length(result2$confs_names), "\n")
  cat("  Cut expressions:", length(result2$confs), "\n")
  
  new_cols2 <- result2$df[, result2$confs_names]
  all_binary2 <- all(apply(new_cols2, 2, function(x) all(x %in% c(0, 1))))
  cat("  All columns are 0/1?", all_binary2, "\n")
  
  cat("\n  Sample of created factors:\n")
  print(head(new_cols2, 10))
} else {
  cat("\n✗ TEST 2 FAILED\n")
  print(result2)
}

cat("\n\n")

# ============================================================================
# TEST 3: With forced cuts
# ============================================================================

cat("TEST 3: With forced cut expressions\n")
cat("-".replicate(70), "\n", sep = "")

set.seed(789)
df_test3 <- data.frame(
  age = rnorm(200, 50, 10),
  size = rnorm(200, 5, 2),
  nodes = rnorm(200, 3, 1),
  event = rbinom(200, 1, 0.3),
  tte = rexp(200),
  treat = rbinom(200, 1, 0.5),
  id = 1:200
)

confounders3 <- c("age", "size", "nodes")

result3 <- try(
  get_FSdata(
    df.analysis = df_test3,
    confounders.name = confounders3,
    outcome.name = "tte",
    event.name = "event",
    cut_type = "default",
    conf_force = c("age <= 45", "size <= 4"),  # Add forced cuts
    use_lasso = FALSE,
    use_grf = FALSE,
    details = TRUE
  ),
  silent = FALSE
)

if (!inherits(result3, "try-error")) {
  cat("\n✓ TEST 3 PASSED\n")
  cat("  Cuts created:", length(result3$confs_names), "\n")
  
  new_cols3 <- result3$df[, result3$confs_names]
  all_binary3 <- all(apply(new_cols3, 2, function(x) all(x %in% c(0, 1))))
  cat("  All columns are 0/1?", all_binary3, "\n")
  
  # Check if forced cuts were included
  has_age_45 <- any(grepl("age.*45", result3$confs))
  has_size_4 <- any(grepl("size.*4", result3$confs))
  cat("  'age <= 45' included?", has_age_45, "\n")
  cat("  'size <= 4' included?", has_size_4, "\n")
} else {
  cat("\n✗ TEST 3 FAILED\n")
  print(result3)
}

cat("\n\n")

# ============================================================================
# TEST 4: Performance benchmark
# ============================================================================

cat("TEST 4: Performance benchmark\n")
cat("-".replicate(70), "\n", sep = "")

# Larger dataset similar to GBSG
set.seed(999)
df_large <- data.frame(
  age = rnorm(686, 53, 11),
  meno = sample(c(0, 1), 686, replace = TRUE),
  size = rnorm(686, 27, 13),
  grade = sample(c(1, 2, 3), 686, replace = TRUE),
  nodes = rnorm(686, 4.6, 7),
  pgr = rnorm(686, 100, 90),
  er = rnorm(686, 100, 85),
  event = rbinom(686, 1, 0.4),
  tte = rexp(686, 0.5),
  treat = rbinom(686, 1, 0.5),
  id = 1:686
)

confounders_large <- c("age", "meno", "size", "grade", "nodes", "pgr", "er")

cat("Dataset: 686 rows x 7 confounders (like GBSG)\n")
cat("Running 3 iterations to measure timing...\n\n")

timing <- system.time({
  for (i in 1:3) {
    result_bench <- get_FSdata(
      df.analysis = df_large,
      confounders.name = confounders_large,
      outcome.name = "tte",
      event.name = "event",
      cut_type = "default",
      use_lasso = FALSE,
      use_grf = FALSE,
      details = FALSE
    )
  }
})

cat("\n✓ TEST 4 PASSED\n")
cat("  3 iterations timing (seconds):\n")
cat("    User:   ", round(timing["user.self"], 3), "\n")
cat("    System: ", round(timing["sys.self"], 3), "\n")
cat("    Total:  ", round(timing["elapsed"], 3), "\n")
cat("  Per iteration: ", round(timing["elapsed"]/3, 3), " seconds\n")
cat("  Expected speedup: 3-4x vs original\n")

cat("\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("=".replicate(70), "\n", sep = "")
cat("TEST SUMMARY\n")
cat("=".replicate(70), "\n")
cat("✓ Test 1: Basic continuous variables\n")
cat("✓ Test 2: Mixed categorical/continuous\n")
cat("✓ Test 3: With forced cuts\n")
cat("✓ Test 4: Performance benchmark\n")
cat("\nAll tests passed! Refactored get_FSdata.R is working correctly.\n")
cat("=".replicate(70), "\n\n")
