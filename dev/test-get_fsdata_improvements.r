
# ============================================================================
# TESTING YOUR CHANGES
# ============================================================================
#
# Add this temporary test code to verify correctness:
# FILE: dev/test_fsdata_improvement.R

rm(list = ls())
devtools::load_all()

# Create simple test data
set.seed(123)
df_test <- data.frame(
  age = rnorm(100, 50, 10),
  size = rnorm(100, 5, 2),
  grade = sample(c(0, 1, 2), 100, replace = TRUE),
  event = rbinom(100, 1, 0.3),
  tte = rexp(100),
  treat = rbinom(100, 1, 0.5),
  id = 1:100
)

# Define confounders
confounders <- c("age", "size", "grade")

# Test 1: Basic functionality
cat("TEST 1: Basic functionality\n")
result <- get_FSdata(
  df.analysis = df_test,
  confounders.name = confounders,
  outcome.name = "tte",
  event.name = "event",
  cut_type = "default",
  use_lasso = FALSE,
  use_grf = FALSE,
  details = TRUE
)

cat("Result contains:", names(result), "\n")
cat("Number of cuts created:", length(result$confs_names), "\n")

# Test 2: Verify column creation worked
cat("\nTEST 2: Verify columns created\n")
new_cols <- result$df[, result$confs_names]
print(head(new_cols))
cat("All columns contain only 0/1?",
    all(apply(new_cols, 2, function(x) all(x %in% c(0, 1)))), "\n")

# Test 3: Compare with old version (if you kept it)
# Uncomment after creating comparison version
# cat("\nTEST 3: Compare old vs new (timing)\n")
# system.time(result_old <- get_FSdata_old(...))
# system.time(result_new <- get_FSdata(...))
# Verify same results: all.equal(result_old, result_new)


# ============================================================================
# OPTIONAL: PERFORMANCE BENCHMARK
# ============================================================================
# Add to test file to measure actual speedup

cat("\n===== PERFORMANCE BENCHMARK =====\n")

# Larger test dataset
df_large <- data.frame(
  age = rnorm(700, 50, 10),
  size = rnorm(700, 5, 2),
  nodes = rnorm(700, 3, 1),
  grade = sample(c(0, 1, 2), 700, replace = TRUE),
  pgr = rnorm(700, 100, 50),
  er = rnorm(700, 150, 70),
  meno = sample(c(0, 1), 700, replace = TRUE),
  event = rbinom(700, 1, 0.3),
  tte = rexp(700),
  treat = rbinom(700, 1, 0.5),
  id = 1:700
)

confounders_large <- c("age", "size", "nodes", "grade", "pgr", "er", "meno")

cat("Dataset: 700 rows x 7 confounders\n")
cat("Running timing test (this will take ~30 seconds)...\n\n")

t1 <- system.time({
  for (i in 1:3) {
    result <- get_FSdata(
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

cat("\nImproved version (3 iterations):\n")
print(t1)
cat("\nPer iteration: ", round(t1[3]/3, 3), " seconds\n")
cat("Expected improvement: 3-4x faster\n")
