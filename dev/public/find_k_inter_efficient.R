# Efficient Method to Find k_inter for Target Harm Subgroup HR
# This provides optimized approaches to calibrate the interaction effect

# ================================================================================
# METHOD 1: Root-Finding with Uniroot (Most Efficient)
# ================================================================================

#' Find k_inter value to achieve target harm subgroup HR
#'
#' @param target_hr_harm Target hazard ratio for the harm subgroup
#' @param data Dataset to use
#' @param continuous_vars Character vector of continuous variables
#' @param factor_vars Character vector of factor variables
#' @param outcome_var Name of outcome variable
#' @param event_var Name of event variable
#' @param treatment_var Name of treatment variable
#' @param subgroup_vars Variables defining subgroup
#' @param subgroup_cuts Cutpoint specifications
#' @param k_treat Treatment effect modifier
#' @param k_inter_range Search range for k_inter (default c(-10, 10))
#' @param tol Tolerance for root finding (default 0.001)
#' @param n_super Size of super population for HR calculation
#' @param verbose Print progress information
#'
#' @return List with k_inter value and achieved HR

find_k_inter_for_target_hr <- function(target_hr_harm,
                                       data,
                                       continuous_vars,
                                       factor_vars,
                                       outcome_var,
                                       event_var,
                                       treatment_var,
                                       subgroup_vars,
                                       subgroup_cuts,
                                       k_treat = 1,
                                       k_inter_range = c(-10, 10),
                                       tol = 0.001,
                                       n_super = 5000,
                                       verbose = TRUE) {

  # Load the generate_aft_dgm_flex function if not already loaded
  if (!exists("generate_aft_dgm_flex")) {
    source("generate_aft_dgm_flexible.R")
  }

  # Objective function: returns difference between achieved and target HR
  objective_function <- function(k_inter_val) {

    # Generate DGM with current k_inter
    dgm_temp <- generate_aft_dgm_flex(
      data = data,
      continuous_vars = continuous_vars,
      factor_vars = factor_vars,
      outcome_var = outcome_var,
      event_var = event_var,
      treatment_var = treatment_var,
      subgroup_vars = subgroup_vars,
      subgroup_cuts = subgroup_cuts,
      model = "alt",
      k_treat = k_treat,
      k_inter = k_inter_val,
      n_super = n_super,
      verbose = FALSE  # Suppress output during search
    )

    # Get the harm subgroup HR
    hr_harm_achieved <- dgm_temp$hazard_ratios$harm_subgroup

    # Return difference from target
    return(hr_harm_achieved - target_hr_harm)
  }

  if (verbose) {
    cat("Searching for k_inter to achieve harm subgroup HR =", target_hr_harm, "\n")
    cat("Search range: [", k_inter_range[1], ",", k_inter_range[2], "]\n")
  }

  # Use uniroot to find k_inter
  tryCatch({
    result <- uniroot(
      f = objective_function,
      interval = k_inter_range,
      tol = tol
    )

    k_inter_optimal <- result$root

    # Verify the result
    dgm_final <- generate_aft_dgm_flex(
      data = data,
      continuous_vars = continuous_vars,
      factor_vars = factor_vars,
      outcome_var = outcome_var,
      event_var = event_var,
      treatment_var = treatment_var,
      subgroup_vars = subgroup_vars,
      subgroup_cuts = subgroup_cuts,
      model = "alt",
      k_treat = k_treat,
      k_inter = k_inter_optimal,
      n_super = n_super,
      verbose = FALSE
    )

    hr_harm_final <- dgm_final$hazard_ratios$harm_subgroup

    if (verbose) {
      cat("\n=== RESULTS ===\n")
      cat("Optimal k_inter:", round(k_inter_optimal, 4), "\n")
      cat("Achieved harm subgroup HR:", round(hr_harm_final, 4), "\n")
      cat("Target harm subgroup HR:", target_hr_harm, "\n")
      cat("Absolute error:", round(abs(hr_harm_final - target_hr_harm), 6), "\n")
      cat("Overall HR:", round(dgm_final$hazard_ratios$overall, 4), "\n")
      cat("No-harm subgroup HR:", round(dgm_final$hazard_ratios$no_harm_subgroup, 4), "\n")
    }

    return(list(
      k_inter = k_inter_optimal,
      achieved_hr_harm = hr_harm_final,
      target_hr_harm = target_hr_harm,
      error = abs(hr_harm_final - target_hr_harm),
      dgm = dgm_final,
      convergence = result$iter
    ))

  }, error = function(e) {
    if (verbose) {
      cat("\nError in root finding. Trying boundary search...\n")
      cat("Error message:", e$message, "\n")
    }

    # If uniroot fails, try boundary values
    obj_lower <- objective_function(k_inter_range[1])
    obj_upper <- objective_function(k_inter_range[2])

    if (verbose) {
      cat("HR at k_inter =", k_inter_range[1], ":", round(obj_lower + target_hr_harm, 4), "\n")
      cat("HR at k_inter =", k_inter_range[2], ":", round(obj_upper + target_hr_harm, 4), "\n")
      cat("\nSolution may be outside the search range.\n")
      cat("Try adjusting k_inter_range parameter.\n")
    }

    return(NULL)
  })
}

# ================================================================================
# METHOD 2: Grid Search (More Robust but Slower)
# ================================================================================

#' Grid search for k_inter to achieve target HR
#'
#' @param target_hr_harm Target hazard ratio for harm subgroup
#' @param k_inter_grid Grid of k_inter values to search
#' @param ... Other parameters passed to generate_aft_dgm_flex
#'
#' @return Best k_inter value and results

find_k_inter_grid_search <- function(target_hr_harm,
                                     k_inter_grid = seq(-5, 5, by = 0.5),
                                     data,
                                     continuous_vars,
                                     factor_vars,
                                     outcome_var,
                                     event_var,
                                     treatment_var,
                                     subgroup_vars,
                                     subgroup_cuts,
                                     k_treat = 1,
                                     n_super = 5000,
                                     verbose = TRUE) {

  if (!exists("generate_aft_dgm_flex")) {
    source("generate_aft_dgm_flexible.R")
  }

  results <- data.frame(
    k_inter = numeric(),
    hr_harm = numeric(),
    hr_overall = numeric(),
    error = numeric()
  )

  if (verbose) {
    cat("Grid search for k_inter...\n")
    pb <- txtProgressBar(min = 0, max = length(k_inter_grid), style = 3)
  }

  for (i in seq_along(k_inter_grid)) {
    k_val <- k_inter_grid[i]
    dgm_temp <- generate_aft_dgm_flex(
      data = data,
      continuous_vars = continuous_vars,
      factor_vars = factor_vars,
      outcome_var = outcome_var,
      event_var = event_var,
      treatment_var = treatment_var,
      subgroup_vars = subgroup_vars,
      subgroup_cuts = subgroup_cuts,
      model = "alt",
      k_treat = k_treat,
      k_inter = k_val,
      n_super = n_super,
      verbose = FALSE
    )

    results <- rbind(results, data.frame(
      k_inter = k_val,
      hr_harm = dgm_temp$hazard_ratios$harm_subgroup,
      hr_overall = dgm_temp$hazard_ratios$overall,
      error = abs(dgm_temp$hazard_ratios$harm_subgroup - target_hr_harm)
    ))

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    cat("\n")
  }

  # Find best k_inter
  best_idx <- which.min(results$error)
  best_result <- results[best_idx, ]

  if (verbose) {
    cat("\n=== GRID SEARCH RESULTS ===\n")
    cat("Target harm subgroup HR:", target_hr_harm, "\n")
    cat("Best k_inter:", round(best_result$k_inter, 4), "\n")
    cat("Achieved HR:", round(best_result$hr_harm, 4), "\n")
    cat("Error:", round(best_result$error, 6), "\n")

    # Show nearby values
    cat("\nNearby values:\n")
    nearby <- results[max(1, best_idx-2):min(nrow(results), best_idx+2), ]
    print(round(nearby, 4))
  }

  return(list(
    k_inter = best_result$k_inter,
    achieved_hr_harm = best_result$hr_harm,
    target_hr_harm = target_hr_harm,
    error = best_result$error,
    all_results = results
  ))
}

# ================================================================================
# METHOD 3: Optimization with optim (Alternative approach)
# ================================================================================

#' Optimization approach to find k_inter
#'
#' @param target_hr_harm Target hazard ratio
#' @param initial_k_inter Starting value for optimization
#' @param ... Other parameters

find_k_inter_optim <- function(target_hr_harm,
                              initial_k_inter = 0,
                              data,
                              continuous_vars,
                              factor_vars,
                              outcome_var,
                              event_var,
                              treatment_var,
                              subgroup_vars,
                              subgroup_cuts,
                              k_treat = 1,
                              n_super = 5000,
                              verbose = TRUE) {

  if (!exists("generate_aft_dgm_flex")) {
    source("generate_aft_dgm_flexible.R")
  }

  # Objective function to minimize (squared error)
  objective <- function(k_inter_val) {
    dgm_temp <- generate_aft_dgm_flex(
      data = data,
      continuous_vars = continuous_vars,
      factor_vars = factor_vars,
      outcome_var = outcome_var,
      event_var = event_var,
      treatment_var = treatment_var,
      subgroup_vars = subgroup_vars,
      subgroup_cuts = subgroup_cuts,
      model = "alt",
      k_treat = k_treat,
      k_inter = k_inter_val,
      n_super = n_super,
      verbose = FALSE
    )

    error <- (dgm_temp$hazard_ratios$harm_subgroup - target_hr_harm)^2
    return(error)
  }

  # Run optimization
  result <- optim(
    par = initial_k_inter,
    fn = objective,
    method = "Brent",
    lower = -10,
    upper = 10
  )

  k_inter_optimal <- result$par

  # Get final DGM
  dgm_final <- generate_aft_dgm_flex(
    data = data,
    continuous_vars = continuous_vars,
    factor_vars = factor_vars,
    outcome_var = outcome_var,
    event_var = event_var,
    treatment_var = treatment_var,
    subgroup_vars = subgroup_vars,
    subgroup_cuts = subgroup_cuts,
    model = "alt",
    k_treat = k_treat,
    k_inter = k_inter_optimal,
    n_super = n_super,
    verbose = FALSE
  )

  if (verbose) {
    cat("\n=== OPTIMIZATION RESULTS ===\n")
    cat("Optimal k_inter:", round(k_inter_optimal, 4), "\n")
    cat("Achieved harm subgroup HR:", round(dgm_final$hazard_ratios$harm_subgroup, 4), "\n")
    cat("Target harm subgroup HR:", target_hr_harm, "\n")
    cat("Convergence code:", result$convergence, "(0 = success)\n")
  }

  return(list(
    k_inter = k_inter_optimal,
    achieved_hr_harm = dgm_final$hazard_ratios$harm_subgroup,
    target_hr_harm = target_hr_harm,
    dgm = dgm_final,
    convergence = result$convergence
  ))
}

# ================================================================================
# EXAMPLE USAGE
# ================================================================================

if (FALSE) {  # Set to TRUE to run examples

  library(survival)
  data(cancer)

  # Define the model parameters
  params <- list(
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
    n_super = 5000
  )

  # Example 1: Find k_inter for harm HR = 2.0
  cat("Example 1: Finding k_inter for harm subgroup HR = 2.0\n")
  cat("=" * 60, "\n")

  result1 <- do.call(find_k_inter_for_target_hr, c(
    list(target_hr_harm = 2.0, verbose = TRUE),
    params
  ))

  # Example 2: Find k_inter for harm HR = 0.5 (protective effect)
  cat("\n\nExample 2: Finding k_inter for harm subgroup HR = 0.5\n")
  cat("=" * 60, "\n")

  result2 <- do.call(find_k_inter_for_target_hr, c(
    list(target_hr_harm = 0.5, verbose = TRUE),
    params
  ))

  # Example 3: Grid search approach
  cat("\n\nExample 3: Grid search for harm subgroup HR = 1.5\n")
  cat("=" * 60, "\n")

  result3 <- do.call(find_k_inter_grid_search, c(
    list(
      target_hr_harm = 1.5,
      k_inter_grid = seq(-3, 3, by = 0.2),
      verbose = TRUE
    ),
    params
  ))

  # Plot the grid search results
  if (!is.null(result3$all_results)) {
    plot(result3$all_results$k_inter, result3$all_results$hr_harm,
         type = "b", pch = 19,
         xlab = "k_inter", ylab = "Harm Subgroup HR",
         main = "Relationship between k_inter and Harm Subgroup HR")
    abline(h = 1.5, col = "red", lty = 2)
    abline(v = result3$k_inter, col = "blue", lty = 2)
    legend("topright",
           legend = c("Target HR", "Optimal k_inter"),
           col = c("red", "blue"), lty = 2)
  }
}

# ================================================================================
# BATCH PROCESSING: Multiple Target HRs
# ================================================================================

#' Find k_inter values for multiple target HRs
#'
#' @param target_hrs Vector of target hazard ratios
#' @param ... Parameters for DGM generation
#'
#' @return Data frame with results for each target HR

find_k_inter_batch <- function(target_hrs, ...) {

  results <- data.frame(
    target_hr = numeric(),
    k_inter = numeric(),
    achieved_hr = numeric(),
    error = numeric(),
    converged = logical()
  )

  cat("Processing", length(target_hrs), "target HRs...\n")
  pb <- txtProgressBar(min = 0, max = length(target_hrs), style = 3)

  for (i in seq_along(target_hrs)) {
    result <- find_k_inter_for_target_hr(
      target_hr_harm = target_hrs[i],
      verbose = FALSE,
      ...
    )

    if (!is.null(result)) {
      results <- rbind(results, data.frame(
        target_hr = target_hrs[i],
        k_inter = result$k_inter,
        achieved_hr = result$achieved_hr_harm,
        error = result$error,
        converged = TRUE
      ))
    } else {
      results <- rbind(results, data.frame(
        target_hr = target_hrs[i],
        k_inter = NA,
        achieved_hr = NA,
        error = NA,
        converged = FALSE
      ))
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)
  cat("\n\nBatch Results:\n")
  print(round(results, 4))

  return(results)
}

# ================================================================================
# SENSITIVITY ANALYSIS: How k_inter affects HRs
# ================================================================================

#' Analyze sensitivity of HRs to k_inter
#'
#' @param k_inter_range Range of k_inter values to analyze
#' @param n_points Number of points to evaluate
#' @param ... Parameters for DGM generation

sensitivity_analysis_k_inter <- function(k_inter_range = c(-5, 5),
                                        n_points = 21,
                                        ...) {

  k_inter_vals <- seq(k_inter_range[1], k_inter_range[2], length.out = n_points)

  results <- data.frame(
    k_inter = numeric(),
    hr_harm = numeric(),
    hr_no_harm = numeric(),
    hr_overall = numeric(),
    subgroup_size = numeric()
  )

  cat("Running sensitivity analysis...\n")
  pb <- txtProgressBar(min = 0, max = n_points, style = 3)

  for (i in seq_along(k_inter_vals)) {
    dgm <- generate_aft_dgm_flex(
      k_inter = k_inter_vals[i],
      verbose = FALSE,
      ...
    )

    results <- rbind(results, data.frame(
      k_inter = k_inter_vals[i],
      hr_harm = dgm$hazard_ratios$harm_subgroup,
      hr_no_harm = dgm$hazard_ratios$no_harm_subgroup,
      hr_overall = dgm$hazard_ratios$overall,
      subgroup_size = dgm$subgroup_info$size
    ))

    setTxtProgressBar(pb, i)
  }

  close(pb)

  # Create visualization
  par(mfrow = c(2, 2))

  # Plot 1: Harm subgroup HR vs k_inter
  plot(results$k_inter, results$hr_harm, type = "l", lwd = 2,
       xlab = "k_inter", ylab = "Harm Subgroup HR",
       main = "Harm Subgroup HR vs k_inter")
  abline(h = 1, lty = 2, col = "gray")
  grid()

  # Plot 2: All HRs vs k_inter
  plot(results$k_inter, results$hr_overall, type = "l", lwd = 2,
       xlab = "k_inter", ylab = "Hazard Ratio",
       main = "All HRs vs k_inter", ylim = range(c(results$hr_harm, results$hr_no_harm, results$hr_overall)))
  lines(results$k_inter, results$hr_harm, col = "red", lwd = 2)
  lines(results$k_inter, results$hr_no_harm, col = "blue", lwd = 2)
  abline(h = 1, lty = 2, col = "gray")
  legend("topright",
         legend = c("Overall", "Harm", "No-harm"),
         col = c("black", "red", "blue"), lty = 1, lwd = 2)
  grid()

  # Plot 3: Ratio of HRs
  plot(results$k_inter, results$hr_harm / results$hr_no_harm, type = "l", lwd = 2,
       xlab = "k_inter", ylab = "HR Ratio (Harm/No-harm)",
       main = "Relative Effect Modification")
  abline(h = 1, lty = 2, col = "gray")
  grid()

  # Plot 4: Table of key values
  plot.new()
  key_indices <- seq(1, nrow(results), length.out = min(10, nrow(results)))
  key_results <- results[key_indices, ]
  text(0.5, 0.9, "Key Values", cex = 1.2, font = 2)
  text(0.5, 0.8, paste("k_inter: HR_harm / HR_no_harm"), cex = 0.9)
  for (i in 1:nrow(key_results)) {
    text(0.5, 0.7 - i*0.07,
         sprintf("%.1f: %.2f / %.2f",
                 key_results$k_inter[i],
                 key_results$hr_harm[i],
                 key_results$hr_no_harm[i]),
         cex = 0.8)
  }

  par(mfrow = c(1, 1))

  return(results)
}
