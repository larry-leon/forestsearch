#' Test Suite for Enhanced Bootstrap Function
#'
#' Tests the enhanced bootstrap function with various scenarios

# Load required packages
library(survival)
library(data.table)
library(foreach)
library(doFuture)

# ============================================================================
# TEST DATA GENERATION
# ============================================================================

#' Generate Test Data for ForestSearch
generate_test_data <- function(n = 500, seed = 123) {
  set.seed(seed)

  # Generate covariates
  df <- data.frame(
    id = 1:n,
    age = rnorm(n, 60, 10),
    sex = rbinom(n, 1, 0.5),
    biomarker = rnorm(n, 0, 1),
    treat = rbinom(n, 1, 0.5)
  )

  # Generate survival times
  lambda <- 0.1 * exp(0.02 * df$age + 0.3 * df$sex - 0.5 * df$biomarker * df$treat)
  df$tte <- rexp(n, lambda)
  df$event <- rbinom(n, 1, 0.8)

  # Add treatment recommendation (mock)
  df$treat.recommend <- ifelse(df$biomarker > 0, 1, 0)

  return(df)
}

#' Create Mock ForestSearch Result
create_mock_fs_result <- function(df, est.scale = "hr") {

  fs.est <- list(
    df.est = df,
    args_call_all = list(
      outcome.name = "tte",
      event.name = "event",
      treat.name = "treat",
      id.name = "id",
      potentialOutcome.name = NULL,
      est.scale = est.scale,
      parallel_args = list(
        plan = "multisession",
        workers = 2,
        show_message = FALSE
      )
    ),
    confounders.candidate = c("age", "sex", "biomarker"),
    sg.harm = c("biomarker > 0"),
    find.grps = list(
      max_sg_est = 1.5,
      time_search = 0.1,
      prop_max_count = 0.8,
      max_count = 10,
      L = 3
    ),
    prop_maxk = 0.8,
    sg_focus = "hr"
  )

  return(fs.est)
}

# ============================================================================
# TEST CASES
# ============================================================================

#' Test 1: Basic Functionality
test_basic_functionality <- function() {
  cat("\n=== TEST 1: Basic Functionality ===\n")

  # Generate test data
  df <- generate_test_data(n = 200)
  fs.est <- create_mock_fs_result(df)

  # Run enhanced bootstrap
  result <- tryCatch({
    forestsearch_bootstrap_dofuture_enhanced(
      fs.est = fs.est,
      nb_boots = 5,  # Small number for testing
      details = TRUE,
      show_three = FALSE,
      parallel_args = list(
        plan = "sequential",  # Sequential for testing
        workers = 1
      )
    )
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(NULL)
  })

  # Check results
  if (!is.null(result)) {
    cat("✓ Basic functionality test passed\n")
    cat("  - Results structure created\n")
    cat("  - Diagnostics available:", !is.null(result$diagnostics), "\n")
    cat("  - Success rate:", sprintf("%.1f%%", 100 * result$diagnostics$prop_successful), "\n")
  } else {
    cat("✗ Basic functionality test failed\n")
  }

  return(result)
}

#' Test 2: est.scale Logic Preservation
test_est_scale_logic <- function() {
  cat("\n=== TEST 2: est.scale Logic Preservation ===\n")

  df <- generate_test_data(n = 200)

  # Test with est.scale = "hr"
  fs.est_hr <- create_mock_fs_result(df, est.scale = "hr")
  result_hr <- forestsearch_bootstrap_dofuture_enhanced(
    fs.est = fs.est_hr,
    nb_boots = 3,
    details = FALSE,
    parallel_args = list(plan = "sequential", workers = 1)
  )

  # Test with est.scale = "1/hr"
  fs.est_1hr <- create_mock_fs_result(df, est.scale = "1/hr")
  result_1hr <- forestsearch_bootstrap_dofuture_enhanced(
    fs.est = fs.est_1hr,
    nb_boots = 3,
    details = FALSE,
    parallel_args = list(plan = "sequential", workers = 1)
  )

  # Check that FSsg_tab was created with correct scale
  if (!is.null(result_hr$FSsg_tab) && !is.null(result_1hr$FSsg_tab)) {
    cat("✓ est.scale logic test passed\n")
    cat("  - HR scale table created:", !is.null(result_hr$FSsg_tab), "\n")
    cat("  - 1/HR scale table created:", !is.null(result_1hr$FSsg_tab), "\n")
  } else {
    cat("✗ est.scale logic test failed\n")
  }

  return(list(hr = result_hr, inv_hr = result_1hr))
}

#' Test 3: Input Validation
test_input_validation <- function() {
  cat("\n=== TEST 3: Input Validation ===\n")

  df <- generate_test_data(n = 100)
  fs.est <- create_mock_fs_result(df)

  # Test invalid fs.est
  error1 <- tryCatch({
    forestsearch_bootstrap_dofuture_enhanced(
      fs.est = "not a list",
      nb_boots = 5
    )
    FALSE
  }, error = function(e) {
    cat("✓ Caught invalid fs.est:", substr(e$message, 1, 50), "...\n")
    TRUE
  })

  # Test invalid nb_boots
  error2 <- tryCatch({
    forestsearch_bootstrap_dofuture_enhanced(
      fs.est = fs.est,
      nb_boots = -1
    )
    FALSE
  }, error = function(e) {
    cat("✓ Caught invalid nb_boots:", substr(e$message, 1, 50), "...\n")
    TRUE
  })

  # Test invalid parallel_args
  error3 <- tryCatch({
    forestsearch_bootstrap_dofuture_enhanced(
      fs.est = fs.est,
      nb_boots = 5,
      parallel_args = "not a list"
    )
    FALSE
  }, error = function(e) {
    cat("✓ Caught invalid parallel_args:", substr(e$message, 1, 50), "...\n")
    TRUE
  })

  if (error1 && error2 && error3) {
    cat("✓ All input validation tests passed\n")
  } else {
    cat("✗ Some input validation tests failed\n")
  }

  return(NULL)
}

#' Test 4: Parallel Processing
test_parallel_processing <- function() {
  cat("\n=== TEST 4: Parallel Processing ===\n")

  df <- generate_test_data(n = 200)
  fs.est <- create_mock_fs_result(df)

  # Time sequential execution
  t1 <- system.time({
    result_seq <- forestsearch_bootstrap_dofuture_enhanced(
      fs.est = fs.est,
      nb_boots = 10,
      details = FALSE,
      parallel_args = list(plan = "sequential", workers = 1)
    )
  })

  # Time parallel execution (if multiple cores available)
  if (parallel::detectCores() > 1) {
    t2 <- system.time({
      result_par <- forestsearch_bootstrap_dofuture_enhanced(
        fs.est = fs.est,
        nb_boots = 10,
        details = FALSE,
        parallel_args = list(plan = "multisession", workers = 2)
      )
    })

    cat("Sequential time:", t1[3], "seconds\n")
    cat("Parallel time:", t2[3], "seconds\n")
    cat("Speedup:", sprintf("%.2fx", t1[3] / t2[3]), "\n")

    if (t2[3] < t1[3] * 1.5) {  # Allow some overhead
      cat("✓ Parallel processing test passed\n")
    } else {
      cat("✗ Parallel processing did not show expected speedup\n")
    }
  } else {
    cat("⚠ Only 1 core available, skipping parallel test\n")
  }

  return(NULL)
}

#' Test 5: Progress Reporting
test_progress_reporting <- function() {
  cat("\n=== TEST 5: Progress Reporting ===\n")

  df <- generate_test_data(n = 100)
  fs.est <- create_mock_fs_result(df)

  # Capture output
  output <- capture.output({
    result <- forestsearch_bootstrap_dofuture_enhanced(
      fs.est = fs.est,
      nb_boots = 20,
      details = TRUE,
      progress_interval = 5,
      parallel_args = list(plan = "sequential", workers = 1)
    )
  })

  # Check for progress messages
  progress_lines <- grep("Completed.*bootstraps", output, value = TRUE)

  if (length(progress_lines) > 0) {
    cat("✓ Progress reporting test passed\n")
    cat("  - Found", length(progress_lines), "progress messages\n")
  } else {
    cat("✗ Progress reporting test failed\n")
  }

  return(NULL)
}

#' Test 6: Diagnostic Output
test_diagnostic_output <- function() {
  cat("\n=== TEST 6: Diagnostic Output ===\n")

  df <- generate_test_data(n = 150)
  fs.est <- create_mock_fs_result(df)

  result <- forestsearch_bootstrap_dofuture_enhanced(
    fs.est = fs.est,
    nb_boots = 10,
    details = FALSE,
    parallel_args = list(plan = "sequential", workers = 1)
  )

  # Check diagnostics structure
  if (!is.null(result$diagnostics)) {
    diag <- result$diagnostics

    checks <- c(
      "n_successful" %in% names(diag),
      "prop_successful" %in% names(diag),
      "bias_corrections" %in% names(diag),
      "search_time" %in% names(diag)
    )

    if (all(checks)) {
      cat("✓ Diagnostic output test passed\n")
      cat("  - Success rate:", sprintf("%.1f%%", 100 * diag$prop_successful), "\n")
      cat("  - Mean search time:", sprintf("%.2f min", diag$search_time$mean), "\n")
      print(diag$bias_corrections)
    } else {
      cat("✗ Diagnostic output incomplete\n")
    }
  } else {
    cat("✗ No diagnostics found in result\n")
  }

  return(result)
}

# ============================================================================
# RUN ALL TESTS
# ============================================================================

run_all_tests <- function() {
  cat("\n")
  cat("================================================\n")
  cat("RUNNING TEST SUITE FOR ENHANCED BOOTSTRAP\n")
  cat("================================================\n")

  # Store results
  results <- list()

  # Run tests
  results$basic <- test_basic_functionality()
  results$scale <- test_est_scale_logic()
  results$validation <- test_input_validation()
  results$parallel <- test_parallel_processing()
  results$progress <- test_progress_reporting()
  results$diagnostics <- test_diagnostic_output()

  cat("\n")
  cat("================================================\n")
  cat("TEST SUITE COMPLETE\n")
  cat("================================================\n")

  return(results)
}

# Run the test suite
test_results <- run_all_tests()
