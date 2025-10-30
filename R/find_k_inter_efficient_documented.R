#' @title Efficient Methods for Calibrating Interaction Effects in AFT Models
#' @description Functions to find the interaction parameter (k_inter) that achieves 
#'   target hazard ratios in subgroups for AFT data generating mechanisms.
#' @author Your Name
#' @keywords survival simulation AFT interaction calibration
#' @import survival
#' @importFrom stats uniroot optim quantile
#' @importFrom utils txtProgressBar setTxtProgressBar

# ================================================================================
# METHOD 1: Root-Finding with Uniroot (Most Efficient)
# ================================================================================

#' Find k_inter Value to Achieve Target Harm Subgroup Hazard Ratio
#' 
#' Uses numerical root-finding to determine the interaction parameter (k_inter) 
#' that achieves a specified target hazard ratio in the harm subgroup. This is 
#' the most efficient method for single target calibration.
#'
#' @param target_hr_harm Numeric value specifying the target hazard ratio for 
#'   the harm subgroup. Must be positive.
#' @param data A data.frame containing the dataset to use for model fitting.
#' @param continuous_vars Character vector of continuous variable names to be 
#'   standardized and included as covariates.
#' @param factor_vars Character vector of factor/categorical variable names to be 
#'   converted to dummy variables.
#' @param outcome_var Character string specifying the name of the outcome/time variable.
#' @param event_var Character string specifying the name of the event/status variable 
#'   (1 = event, 0 = censored).
#' @param treatment_var Character string specifying the name of the treatment variable.
#' @param subgroup_vars Character vector of variable names defining the subgroup.
#' @param subgroup_cuts Named list of cutpoint specifications for subgroup variables. 
#'   See \code{\link{generate_aft_dgm_flex}} for details on flexible specifications.
#' @param k_treat Numeric value for treatment effect modifier. Default is 1 
#'   (no modification).
#' @param k_inter_range Numeric vector of length 2 specifying the search range 
#'   for k_inter. Default is c(-10, 10).
#' @param tol Numeric value specifying tolerance for root finding convergence. 
#'   Default is 0.001.
#' @param n_super Integer specifying size of super population for hazard ratio 
#'   calculation. Default is 5000.
#' @param verbose Logical indicating whether to print progress information. 
#'   Default is TRUE.
#'
#' @return A list of class "k_inter_result" containing:
#' \describe{
#'   \item{k_inter}{Numeric value of optimal k_inter parameter}
#'   \item{achieved_hr_harm}{Numeric value of achieved hazard ratio in harm subgroup}
#'   \item{target_hr_harm}{Numeric value of target hazard ratio (for reference)}
#'   \item{error}{Numeric value of absolute error between achieved and target HR}
#'   \item{dgm}{Object of class "aft_dgm_flex" containing the final DGM}
#'   \item{convergence}{Integer number of iterations to convergence}
#'   \item{method}{Character string "root-finding" indicating method used}
#' }
#'
#' @details
#' This function uses the \code{uniroot} algorithm to solve the equation:
#' \deqn{HR_{harm}(k_{inter}) - HR_{target} = 0}
#' 
#' The algorithm typically converges within 5-10 iterations and achieves high 
#' precision (within the specified tolerance). If the root-finding fails, the 
#' function evaluates the boundaries and provides diagnostic information.
#'
#' @seealso 
#' \code{\link{find_k_inter_grid_search}} for grid search method
#' \code{\link{find_k_inter_batch}} for processing multiple targets
#' \code{\link{sensitivity_analysis_k_inter}} for sensitivity analysis
#' \code{\link{generate_aft_dgm_flex}} for DGM generation
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data(cancer)
#' 
#' # Find k_inter for target HR = 2.0 in harm subgroup
#' result <- find_k_inter_for_target_hr(
#'   target_hr_harm = 2.0,
#'   data = gbsg,
#'   continuous_vars = c("age", "er", "pgr"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(
#'     er = list(type = "quantile", value = 0.25),
#'     meno = 0
#'   ),
#'   k_treat = 1.0,
#'   verbose = TRUE
#' )
#' 
#' cat("Optimal k_inter:", result$k_inter, "\n")
#' cat("Achieved HR:", result$achieved_hr_harm, "\n")
#' }
#' 
#' @export
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
  
  # Input validation
  if (target_hr_harm <= 0) {
    stop("target_hr_harm must be positive")
  }
  if (length(k_inter_range) != 2 || k_inter_range[1] >= k_inter_range[2]) {
    stop("k_inter_range must be a vector of length 2 with min < max")
  }
  
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
    
    output <- list(
      k_inter = k_inter_optimal,
      achieved_hr_harm = hr_harm_final,
      target_hr_harm = target_hr_harm,
      error = abs(hr_harm_final - target_hr_harm),
      dgm = dgm_final,
      convergence = result$iter,
      method = "root-finding"
    )
    
    class(output) <- c("k_inter_result", "list")
    return(output)
    
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

#' Grid Search for k_inter to Achieve Target Hazard Ratio
#' 
#' Uses grid search to find the interaction parameter (k_inter) that achieves 
#' a specified target hazard ratio in the harm subgroup. More robust than 
#' root-finding but slower.
#'
#' @param target_hr_harm Numeric value specifying the target hazard ratio for 
#'   the harm subgroup. Must be positive.
#' @param k_inter_grid Numeric vector of k_inter values to evaluate. Default is 
#'   seq(-5, 5, by = 0.5).
#' @param data A data.frame containing the dataset to use for model fitting.
#' @param continuous_vars Character vector of continuous variable names.
#' @param factor_vars Character vector of factor/categorical variable names.
#' @param outcome_var Character string specifying the outcome/time variable name.
#' @param event_var Character string specifying the event/status variable name.
#' @param treatment_var Character string specifying the treatment variable name.
#' @param subgroup_vars Character vector of subgroup-defining variable names.
#' @param subgroup_cuts Named list of cutpoint specifications.
#' @param k_treat Numeric treatment effect modifier. Default is 1.
#' @param n_super Integer size of super population. Default is 5000.
#' @param verbose Logical indicating whether to show progress. Default is TRUE.
#'
#' @return A list of class "k_inter_grid_result" containing:
#' \describe{
#'   \item{k_inter}{Numeric value of best k_inter from grid}
#'   \item{achieved_hr_harm}{Numeric value of achieved hazard ratio}
#'   \item{target_hr_harm}{Numeric value of target hazard ratio}
#'   \item{error}{Numeric value of minimum absolute error}
#'   \item{all_results}{Data.frame with all grid search results}
#'   \item{method}{Character string "grid-search"}
#' }
#'
#' @details
#' This function evaluates the hazard ratio at each point in the grid and 
#' selects the k_inter value that minimizes the absolute error from the target. 
#' The complete results are returned for visualization of the relationship 
#' between k_inter and hazard ratios.
#'
#' @seealso 
#' \code{\link{find_k_inter_for_target_hr}} for root-finding method
#' \code{\link{plot.k_inter_grid_result}} for plotting results
#'
#' @examples
#' \dontrun{
#' # Grid search for k_inter
#' result <- find_k_inter_grid_search(
#'   target_hr_harm = 1.5,
#'   k_inter_grid = seq(-3, 3, by = 0.2),
#'   data = gbsg,
#'   continuous_vars = c("age", "er", "pgr"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(
#'     er = list(type = "quantile", value = 0.25),
#'     meno = 0
#'   )
#' )
#' 
#' # Plot results
#' plot(result$all_results$k_inter, result$all_results$hr_harm,
#'      type = "b", xlab = "k_inter", ylab = "Harm HR")
#' abline(h = 1.5, col = "red")
#' }
#' 
#' @export
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
  
  # Input validation
  if (target_hr_harm <= 0) {
    stop("target_hr_harm must be positive")
  }
  if (length(k_inter_grid) < 2) {
    stop("k_inter_grid must contain at least 2 values")
  }
  
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
  
  output <- list(
    k_inter = best_result$k_inter,
    achieved_hr_harm = best_result$hr_harm,
    target_hr_harm = target_hr_harm,
    error = best_result$error,
    all_results = results,
    method = "grid-search"
  )
  
  class(output) <- c("k_inter_grid_result", "list")
  return(output)
}

# ================================================================================
# METHOD 3: Optimization with optim (Alternative approach)
# ================================================================================

#' Optimization Method for Finding k_inter
#' 
#' Uses optimization (Brent's method) to find the interaction parameter that 
#' minimizes the squared error between achieved and target hazard ratio.
#'
#' @param target_hr_harm Numeric target hazard ratio for harm subgroup.
#' @param initial_k_inter Numeric starting value for optimization. Default is 0.
#' @param data A data.frame containing the dataset.
#' @param continuous_vars Character vector of continuous variable names.
#' @param factor_vars Character vector of factor variable names.
#' @param outcome_var Character string of outcome variable name.
#' @param event_var Character string of event variable name.
#' @param treatment_var Character string of treatment variable name.
#' @param subgroup_vars Character vector of subgroup variable names.
#' @param subgroup_cuts Named list of cutpoint specifications.
#' @param k_treat Numeric treatment effect modifier. Default is 1.
#' @param n_super Integer size of super population. Default is 5000.
#' @param verbose Logical for progress output. Default is TRUE.
#'
#' @return A list of class "k_inter_optim_result" containing:
#' \describe{
#'   \item{k_inter}{Numeric optimal k_inter value}
#'   \item{achieved_hr_harm}{Numeric achieved hazard ratio}
#'   \item{target_hr_harm}{Numeric target hazard ratio}
#'   \item{dgm}{Object of class "aft_dgm_flex"}
#'   \item{convergence}{Integer convergence code (0 = success)}
#'   \item{method}{Character string "optimization"}
#' }
#'
#' @details
#' This function minimizes the squared error:
#' \deqn{(HR_{harm}(k_{inter}) - HR_{target})^2}
#' 
#' using bounded optimization (Brent's method) which is suitable for 
#' one-dimensional optimization problems.
#'
#' @examples
#' \dontrun{
#' result <- find_k_inter_optim(
#'   target_hr_harm = 2.0,
#'   initial_k_inter = 0,
#'   data = gbsg,
#'   continuous_vars = c("age", "er", "pgr"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(er = 20, meno = 0)
#' )
#' }
#' 
#' @export
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
  
  if (target_hr_harm <= 0) {
    stop("target_hr_harm must be positive")
  }
  
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
  
  output <- list(
    k_inter = k_inter_optimal,
    achieved_hr_harm = dgm_final$hazard_ratios$harm_subgroup,
    target_hr_harm = target_hr_harm,
    dgm = dgm_final,
    convergence = result$convergence,
    method = "optimization"
  )
  
  class(output) <- c("k_inter_optim_result", "list")
  return(output)
}

# ================================================================================
# BATCH PROCESSING: Multiple Target HRs
# ================================================================================

#' Batch Processing of Multiple Target Hazard Ratios
#' 
#' Finds k_inter values for multiple target hazard ratios using the root-finding 
#' method. Useful for creating multiple simulation scenarios.
#'
#' @param target_hrs Numeric vector of target hazard ratios to process.
#' @param ... Additional arguments passed to \code{\link{find_k_inter_for_target_hr}}.
#'
#' @return A data.frame of class "k_inter_batch_result" with columns:
#' \describe{
#'   \item{target_hr}{Numeric target hazard ratio}
#'   \item{k_inter}{Numeric optimal k_inter value}
#'   \item{achieved_hr}{Numeric achieved hazard ratio}
#'   \item{error}{Numeric absolute error}
#'   \item{converged}{Logical indicating successful convergence}
#' }
#'
#' @details
#' This function processes multiple target hazard ratios sequentially, using the 
#' root-finding method for each. Progress is displayed via a text progress bar. 
#' Failed convergences are marked in the output.
#'
#' @examples
#' \dontrun{
#' # Process multiple scenarios
#' target_hrs <- c(0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
#' 
#' batch_results <- find_k_inter_batch(
#'   target_hrs = target_hrs,
#'   data = gbsg,
#'   continuous_vars = c("age", "er", "pgr"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(er = 20, meno = 0)
#' )
#' 
#' print(batch_results)
#' }
#' 
#' @export
find_k_inter_batch <- function(target_hrs, ...) {
  
  if (any(target_hrs <= 0)) {
    stop("All target_hrs must be positive")
  }
  
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
  
  class(results) <- c("k_inter_batch_result", "data.frame")
  return(results)
}

# ================================================================================
# SENSITIVITY ANALYSIS: How k_inter affects HRs
# ================================================================================

#' Sensitivity Analysis of Hazard Ratios to k_inter
#' 
#' Analyzes how the interaction parameter k_inter affects hazard ratios in 
#' different populations (overall, harm subgroup, no-harm subgroup).
#'
#' @param k_inter_range Numeric vector of length 2 specifying the range of 
#'   k_inter values to analyze. Default is c(-5, 5).
#' @param n_points Integer number of points to evaluate within the range. 
#'   Default is 21.
#' @param plot Logical indicating whether to create visualization plots. 
#'   Default is TRUE.
#' @param ... Additional arguments passed to \code{\link{generate_aft_dgm_flex}}.
#'
#' @return A data.frame of class "k_inter_sensitivity" with columns:
#' \describe{
#'   \item{k_inter}{Numeric k_inter value}
#'   \item{hr_harm}{Numeric hazard ratio in harm subgroup}
#'   \item{hr_no_harm}{Numeric hazard ratio in no-harm subgroup}
#'   \item{hr_overall}{Numeric overall hazard ratio}
#'   \item{subgroup_size}{Integer size of harm subgroup}
#' }
#'
#' @details
#' This function evaluates the hazard ratios at evenly spaced points across 
#' the k_inter range. If plot = TRUE, it creates a 4-panel visualization showing:
#' \enumerate{
#'   \item Harm subgroup HR vs k_inter
#'   \item All HRs (overall, harm, no-harm) vs k_inter
#'   \item Ratio of HRs (harm/no-harm) showing effect modification
#'   \item Table of key values
#' }
#'
#' @examples
#' \dontrun{
#' # Analyze sensitivity to k_inter
#' sensitivity_results <- sensitivity_analysis_k_inter(
#'   k_inter_range = c(-2, 2),
#'   n_points = 11,
#'   data = gbsg,
#'   continuous_vars = c("age", "er", "pgr"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(er = 20, meno = 0),
#'   model = "alt",
#'   plot = TRUE
#' )
#' 
#' # Results show relationship between k_inter and HRs
#' print(sensitivity_results)
#' }
#' 
#' @export
#' @importFrom graphics par plot lines abline legend grid text
sensitivity_analysis_k_inter <- function(k_inter_range = c(-5, 5),
                                        n_points = 21,
                                        plot = TRUE,
                                        ...) {
  
  if (length(k_inter_range) != 2 || k_inter_range[1] >= k_inter_range[2]) {
    stop("k_inter_range must be a vector of length 2 with min < max")
  }
  if (n_points < 2) {
    stop("n_points must be at least 2")
  }
  
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
  
  if (plot) {
    # Create visualization
    oldpar <- par(mfrow = c(2, 2))
    on.exit(par(oldpar))
    
    # Plot 1: Harm subgroup HR vs k_inter
    plot(results$k_inter, results$hr_harm, type = "l", lwd = 2,
         xlab = "k_inter", ylab = "Harm Subgroup HR",
         main = "Harm Subgroup HR vs k_inter")
    abline(h = 1, lty = 2, col = "gray")
    grid()
    
    # Plot 2: All HRs vs k_inter
    plot(results$k_inter, results$hr_overall, type = "l", lwd = 2,
         xlab = "k_inter", ylab = "Hazard Ratio",
         main = "All HRs vs k_inter", 
         ylim = range(c(results$hr_harm, results$hr_no_harm, results$hr_overall)))
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
  }
  
  class(results) <- c("k_inter_sensitivity", "data.frame")
  return(results)
}

# ================================================================================
# S3 Methods for Result Objects
# ================================================================================

#' Print Method for k_inter Result Objects
#' 
#' @param x An object of class "k_inter_result"
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisible NULL
#' @export
#' @method print k_inter_result
print.k_inter_result <- function(x, digits = 4, ...) {
  cat("\nk_inter Calibration Result (", x$method, ")\n", sep = "")
  cat("==================================\n")
  cat("Target HR (harm subgroup):", round(x$target_hr_harm, digits), "\n")
  cat("Optimal k_inter:", round(x$k_inter, digits), "\n")
  cat("Achieved HR:", round(x$achieved_hr_harm, digits), "\n")
  cat("Absolute error:", format(x$error, scientific = TRUE, digits = 2), "\n")
  if (!is.null(x$convergence)) {
    cat("Convergence iterations:", x$convergence, "\n")
  }
  invisible(NULL)
}

#' Plot Method for Grid Search Results
#' 
#' @param x An object of class "k_inter_grid_result"
#' @param ... Additional graphical parameters
#' @return Invisible NULL
#' @export
#' @method plot k_inter_grid_result
plot.k_inter_grid_result <- function(x, ...) {
  plot(x$all_results$k_inter, x$all_results$hr_harm,
       type = "b", pch = 19,
       xlab = "k_inter", ylab = "Harm Subgroup HR",
       main = "Grid Search: k_inter vs Harm Subgroup HR", ...)
  abline(h = x$target_hr_harm, col = "red", lty = 2)
  abline(v = x$k_inter, col = "blue", lty = 2)
  legend("topleft",
         legend = c(paste("Target HR =", x$target_hr_harm),
                   paste("Optimal k_inter =", round(x$k_inter, 3))),
         col = c("red", "blue"), lty = 2)
  grid()
  invisible(NULL)
}

#' Summary Method for Batch Results
#' 
#' @param object An object of class "k_inter_batch_result"
#' @param ... Additional arguments (ignored)
#' @return Summary statistics
#' @export
#' @method summary k_inter_batch_result
summary.k_inter_batch_result <- function(object, ...) {
  cat("\nBatch k_inter Calibration Summary\n")
  cat("==================================\n")
  cat("Number of targets:", nrow(object), "\n")
  cat("Successful convergence:", sum(object$converged), "/", nrow(object), "\n")
  cat("\nTarget HR range:", range(object$target_hr), "\n")
  cat("k_inter range:", range(object$k_inter, na.rm = TRUE), "\n")
  cat("Mean absolute error:", mean(object$error, na.rm = TRUE), "\n")
  
  if (any(!object$converged)) {
    cat("\nFailed targets:", object$target_hr[!object$converged], "\n")
  }
  
  invisible(object)
}
