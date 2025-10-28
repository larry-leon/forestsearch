# ==============================================================================
# REFACTORED SUBGROUP ANALYSIS FUNCTIONS
# Complete implementation with modular design
# ==============================================================================

# ------------------------------------------------------------------------------
# ANALYSIS CONFIGURATION
# ------------------------------------------------------------------------------

#' Define Cox model analysis specifications
#'
#' @param wname Name of prognostic covariate
#' @return List of analysis specifications
#' @keywords internal
define_cox_analyses <- function(wname) {
  list(
    # Analysis 1: Stratified by randomization factors
    list(
      id = 1,
      name = "sR",
      description = "Stratified by randomization factors",
      formula_template = "Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR)",
      requires_strata = TRUE,
      min_strata_size = NULL
    ),
    
    # Analysis 2: Unadjusted (no stratification)
    list(
      id = 2,
      name = "none",
      description = "Unadjusted Cox model",
      formula_template = "Surv(y.sim, event.sim) ~ treat.sim",
      requires_strata = FALSE,
      min_strata_size = NULL
    ),
    
    # Analysis 3: Stratified by prognostic factor W
    list(
      id = 3,
      name = "sW",
      description = "Stratified by prognostic factor",
      formula_template = sprintf("Surv(y.sim, event.sim) ~ treat.sim + strata(%s)", wname),
      requires_strata = TRUE,
      min_strata_size = NULL
    ),
    
    # Analysis 4: Stratified by biomarker cutpoint
    list(
      id = 4,
      name = "sBM",
      description = "Stratified by biomarker cutpoint",
      formula_template = "Surv(y.sim, event.sim) ~ treat.sim + strata(cut_opt)",
      requires_strata = TRUE,
      min_strata_size = 5,
      strata_var = "cut_opt"
    ),
    
    # Analysis 5: Stratified by both W and randomization
    list(
      id = 5,
      name = "sW+sR",
      description = "Stratified by W and randomization factors",
      formula_template = sprintf(
        "Surv(y.sim, event.sim) ~ treat.sim + strata(%s) + strata(strata.simR)",
        wname
      ),
      requires_strata = TRUE,
      min_strata_size = NULL
    ),
    
    # Analysis 6: Adjusted for W and age, stratified by randomization
    list(
      id = 6,
      name = "W+age+sR",
      description = "Adjusted for W and age, stratified by randomization",
      formula_template = sprintf(
        "Surv(y.sim, event.sim) ~ treat.sim + %s + age + strata(strata.simR)",
        wname
      ),
      requires_strata = TRUE,
      min_strata_size = NULL
    )
  )
}

# ------------------------------------------------------------------------------
# MAIN SUBGROUP ANALYSIS FUNCTION (REFACTORED)
# ------------------------------------------------------------------------------

#' Conduct Multiple Subgroup Analyses Across Simulations
#'
#' @param dgm Data generating mechanism from get_dgm_stratified()
#' @param subgroups_name Character vector of subgroup names
#' @param subgroups_id Character vector of subgroup definitions (R expressions)
#' @param sims Number of simulations to run
#' @param wname Name of prognostic covariate
#' @param bw Log hazard ratio for prognostic covariate
#' @param Ndraw Number of observations per simulation
#' @param bmcut Biomarker cutpoint for defining subgroups
#' @param outfile Optional path to save results
#' @param verbose Logical. Print progress?
#'
#' @return List containing HR and upper bound estimates for all analyses
#' @export
get_SGanalyses <- function(dgm, 
                           subgroups_name, 
                           subgroups_id, 
                           sims, 
                           wname, 
                           bw, 
                           Ndraw = nrow(dgm$df_super),
                           bmcut = log(1),
                           outfile = NULL,
                           verbose = TRUE) {
  
  # Validate inputs
  validate_subgroup_inputs(subgroups_name, subgroups_id)
  
  # Define analysis specifications
  analyses <- define_cox_analyses(wname)
  
  # Initialize results storage
  results <- initialize_results_storage(sims, subgroups_name, analyses)
  
  # Progress tracking
  if (verbose) {
    cat(sprintf("\nRunning %d simulations across %d subgroups...\n", 
                sims, length(subgroups_name)))
    cat(sprintf("Biomarker cutpoint: %.3f\n", bmcut))
  }
  
  # Main simulation loop
  for (ss in seq_len(sims)) {
    if (verbose && ss %% 100 == 0) {
      cat(sprintf("  Completed %d/%d simulations (%.1f%%)\n", 
                  ss, sims, 100 * ss / sims))
    }
    
    # Generate simulation dataset
    df_sim <- generate_simulation_dataset(dgm, ss, Ndraw, wname, bw, bmcut)
    
    # Analyze all subgroups
    for (gg in seq_along(subgroups_id)) {
      # Extract subgroup
      df_sg <- subset_by_expression(df_sim, subgroups_id[gg])
      
      # Store sample size
      results$ns[ss, gg] <- nrow(df_sg)
      
      # Fit all analysis models
      results <- fit_all_analyses_for_subgroup(
        df_sg = df_sg,
        analyses = analyses,
        results = results,
        sim_idx = ss,
        sg_idx = gg
      )
    }
    
    # Optional: Plot first few simulations
    if (ss <= 10) {
      plot_km_for_simulation(df_sim, ss)
    }
  }
  
  if (verbose) {
    cat(sprintf("\nCompleted all %d simulations!\n", sims))
    print_analysis_summary(results)
  }
  
  # Save results if requested
  if (!is.null(outfile)) {
    save_subgroup_results(results, dgm, subgroups_name, outfile)
  }
  
  return(results)
}

# ------------------------------------------------------------------------------
# VALIDATION
# ------------------------------------------------------------------------------

#' Validate subgroup inputs
#' @keywords internal
validate_subgroup_inputs <- function(subgroups_name, subgroups_id) {
  if (length(subgroups_name) != length(subgroups_id)) {
    stop(sprintf(
      "Subgroup names (%d) and IDs (%d) have different lengths",
      length(subgroups_name), length(subgroups_id)
    ), call. = FALSE)
  }
  
  # Check for duplicates
  if (any(duplicated(subgroups_name))) {
    stop("Duplicate subgroup names detected", call. = FALSE)
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# RESULTS STORAGE INITIALIZATION
# ------------------------------------------------------------------------------

#' Initialize results storage matrices
#' @keywords internal
initialize_results_storage <- function(sims, subgroups_name, analyses) {
  # Create template matrix
  template_matrix <- matrix(
    nrow = sims,
    ncol = length(subgroups_name),
    dimnames = list(NULL, subgroups_name)
  )
  
  results <- list(
    ns = template_matrix  # Sample sizes
  )
  
  # Create storage for each analysis
  for (i in seq_along(analyses)) {
    results[[paste0("hrs", i)]] <- template_matrix  # Hazard ratios
    results[[paste0("ubs", i)]] <- template_matrix  # Upper bounds
  }
  
  # Store metadata
  results$metadata <- list(
    n_sims = sims,
    n_subgroups = length(subgroups_name),
    subgroup_names = subgroups_name,
    n_analyses = length(analyses),
    analysis_names = sapply(analyses, function(x) x$name)
  )
  
  return(results)
}

# ------------------------------------------------------------------------------
# SIMULATION DATASET GENERATION
# ------------------------------------------------------------------------------

#' Generate single simulation dataset
#' @keywords internal
generate_simulation_dataset <- function(dgm, ss, Ndraw, wname, bw, bmcut) {
  # Draw simulated data
  df_sim <- draw_sim_stratified(
    dgm = dgm,
    ss = ss,
    Ndraw = Ndraw,
    wname = wname,
    bw = bw,
    strata_rand = "stratum",
    checking = FALSE,
    details = FALSE,
    return_df = TRUE
  )
  
  # Add ITT flag
  df_sim$itt <- "Y"
  
  # Add biomarker cutpoint indicator
  df_sim$cut_opt <- ifelse(df_sim$z < bmcut, 0, 1)
  
  # Add prognostic factor alias for formulas
  df_sim$w <- df_sim[[wname]]
  
  return(df_sim)
}

#' Subset data by subgroup expression
#' @keywords internal
subset_by_expression <- function(df, expression_string) {
  tryCatch({
    subset(df, eval(parse(text = expression_string)))
  }, error = function(e) {
    warning(sprintf(
      "Failed to subset using expression: %s\nError: %s",
      expression_string, e$message
    ))
    # Return empty dataframe with same structure
    df[0, ]
  })
}

# ------------------------------------------------------------------------------
# MODEL FITTING
# ------------------------------------------------------------------------------

#' Fit all analysis models for a single subgroup
#' @keywords internal
fit_all_analyses_for_subgroup <- function(df_sg, analyses, results, 
                                         sim_idx, sg_idx) {
  # Skip if subgroup is too small
  if (nrow(df_sg) < 10) {
    return(results)
  }
  
  for (analysis in analyses) {
    # Check if analysis is feasible
    if (!is_analysis_feasible(df_sg, analysis)) {
      next
    }
    
    # Fit Cox model
    fit_result <- fit_cox_model_safely(analysis$formula_template, df_sg)
    
    # Store results
    if (!is.null(fit_result)) {
      hr_col <- paste0("hrs", analysis$id)
      ub_col <- paste0("ubs", analysis$id)
      
      results[[hr_col]][sim_idx, sg_idx] <- fit_result$hr
      results[[ub_col]][sim_idx, sg_idx] <- fit_result$upper_bound
    }
  }
  
  return(results)
}

#' Check if analysis is feasible for given data
#' @keywords internal
is_analysis_feasible <- function(df, analysis) {
  # Check minimum sample size
  if (nrow(df) < 10) {
    return(FALSE)
  }
  
  # Check minimum events
  if (sum(df$event.sim) < 5) {
    return(FALSE)
  }
  
  # Check stratification requirements
  if (!is.null(analysis$min_strata_size) && !is.null(analysis$strata_var)) {
    strata_var <- analysis$strata_var
    
    if (!strata_var %in% names(df)) {
      return(FALSE)
    }
    
    strata_counts <- table(df[[strata_var]])
    
    if (any(strata_counts < analysis$min_strata_size)) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Fit Cox model with error handling
#' @keywords internal
fit_cox_model_safely <- function(formula_string, data) {
  tryCatch({
    # Fit Cox model
    fit <- survival::coxph(
      as.formula(formula_string),
      data = data
    )
    
    # Extract summary
    fit_summary <- summary(fit)
    
    # Return hazard ratio and confidence interval
    list(
      hr = fit_summary$conf.int[1, 1],
      lower_bound = fit_summary$conf.int[1, 3],
      upper_bound = fit_summary$conf.int[1, 4],
      pvalue = fit_summary$coefficients[1, 5]
    )
    
  }, error = function(e) {
    # Return NULL on failure
    NULL
  })
}

# ------------------------------------------------------------------------------
# PLOTTING
# ------------------------------------------------------------------------------

#' Plot Kaplan-Meier curves for simulation
#' @keywords internal
plot_km_for_simulation <- function(df_sim, sim_idx) {
  tryCatch({
    # Only plot if plotting function is available
    if (exists("KM.plot.2sample.weighted", mode = "function")) {
      KM.plot.2sample.weighted(
        df = df_sim,
        tte.name = "y.sim",
        event.name = "event.sim",
        treat.name = "treat.sim",
        strata.name = "strata.simO",
        stop.onerror = FALSE,
        risk.set = TRUE,
        by.risk = 12,
        risk.cex = 0.8,
        Xlab = "Months",
        Ylab = "Overall Survival",
        details = FALSE,
        show.ticks = TRUE,
        arms = c("Treat", "Control")
      )
      title(main = sprintf("Simulation %d", sim_idx))
    }
  }, error = function(e) {
    # Silently skip plotting on error
    invisible(NULL)
  })
}

# ------------------------------------------------------------------------------
# REPORTING AND SAVING
# ------------------------------------------------------------------------------

#' Print analysis summary
#' @keywords internal
print_analysis_summary <- function(results) {
  cat("\n=== Analysis Summary ===\n")
  
  # Count successful fits
  for (i in seq_len(results$metadata$n_analyses)) {
    hr_col <- paste0("hrs", i)
    analysis_name <- results$metadata$analysis_names[i]
    
    n_success <- sum(!is.na(results[[hr_col]]))
    n_total <- length(results[[hr_col]])
    success_rate <- 100 * n_success / n_total
    
    cat(sprintf(
      "Analysis %d (%s): %d/%d successful (%.1f%%)\n",
      i, analysis_name, n_success, n_total, success_rate
    ))
  }
  
  cat("\nSubgroup sample sizes (mean across simulations):\n")
  mean_ns <- colMeans(results$ns, na.rm = TRUE)
  for (i in seq_along(mean_ns)) {
    cat(sprintf("  %s: %.1f\n", names(mean_ns)[i], mean_ns[i]))
  }
}

#' Save results to file
#' @keywords internal
save_subgroup_results <- function(results, dgm, subgroups_name, outfile) {
  # Extract results components for saving
  subgroup_ns <- results$ns
  
  # Extract HR and UB matrices
  n_analyses <- results$metadata$n_analyses
  
  save_list <- list(
    dgm = dgm,
    subgroups_name = subgroups_name,
    subgroup_ns = subgroup_ns
  )
  
  for (i in seq_len(n_analyses)) {
    save_list[[paste0("subgroup_hrs", i)]] <- results[[paste0("hrs", i)]]
    save_list[[paste0("subgroup_ubs", i)]] <- results[[paste0("ubs", i)]]
  }
  
  # Save to file
  do.call(save, c(save_list, list(file = outfile)))
  
  cat(sprintf("\nResults saved to: %s\n", outfile))
}

# ------------------------------------------------------------------------------
# POST-PROCESSING UTILITIES
# ------------------------------------------------------------------------------

#' Extract results for specific subgroup
#'
#' @param results Output from get_SGanalyses()
#' @param subgroup_name Name of subgroup to extract
#' @return List with HR and UB for all analyses
#' @export
extract_subgroup_results <- function(results, subgroup_name) {
  if (!subgroup_name %in% results$metadata$subgroup_names) {
    stop(sprintf("Subgroup '%s' not found in results", subgroup_name))
  }
  
  sg_idx <- which(results$metadata$subgroup_names == subgroup_name)
  
  output <- list(
    subgroup = subgroup_name,
    sample_size = results$ns[, sg_idx]
  )
  
  for (i in seq_len(results$metadata$n_analyses)) {
    analysis_name <- results$metadata$analysis_names[i]
    output[[paste0("hr_", analysis_name)]] <- results[[paste0("hrs", i)]][, sg_idx]
    output[[paste0("ub_", analysis_name)]] <- results[[paste0("ubs", i)]][, sg_idx]
  }
  
  return(output)
}

#' Calculate summary statistics for subgroup analyses
#'
#' @param results Output from get_SGanalyses()
#' @return Data frame with summary statistics
#' @export
summarize_subgroup_results <- function(results) {
  summary_df <- data.frame()
  
  for (sg in results$metadata$subgroup_names) {
    sg_results <- extract_subgroup_results(results, sg)
    
    for (i in seq_len(results$metadata$n_analyses)) {
      analysis_name <- results$metadata$analysis_names[i]
      hr_col <- paste0("hr_", analysis_name)
      ub_col <- paste0("ub_", analysis_name)
      
      hrs <- sg_results[[hr_col]]
      ubs <- sg_results[[ub_col]]
      
      row <- data.frame(
        subgroup = sg,
        analysis = analysis_name,
        mean_n = mean(sg_results$sample_size, na.rm = TRUE),
        n_success = sum(!is.na(hrs)),
        median_hr = median(hrs, na.rm = TRUE),
        q25_hr = quantile(hrs, 0.25, na.rm = TRUE),
        q75_hr = quantile(hrs, 0.75, na.rm = TRUE),
        prop_ub_lt_1 = mean(ubs <= 1, na.rm = TRUE),
        median_ub = median(ubs, na.rm = TRUE)
      )
      
      summary_df <- rbind(summary_df, row)
    }
  }
  
  return(summary_df)
}

# ------------------------------------------------------------------------------
# PARALLEL EXECUTION VARIANT
# ------------------------------------------------------------------------------

#' Run subgroup analyses in parallel
#'
#' @param dgm Data generating mechanism
#' @param subgroups_name Character vector of subgroup names
#' @param subgroups_id Character vector of subgroup definitions
#' @param sims Number of simulations
#' @param wname Name of prognostic covariate
#' @param bw Log hazard ratio for prognostic covariate
#' @param Ndraw Number of observations per simulation
#' @param bmcut Biomarker cutpoint
#' @param n_cores Number of parallel cores to use
#' @param outfile Optional path to save results
#'
#' @return List containing results
#' @export
get_SGanalyses_parallel <- function(dgm, subgroups_name, subgroups_id, sims,
                                   wname, bw, Ndraw = nrow(dgm$df_super),
                                   bmcut = log(1), n_cores = 4,
                                   outfile = NULL) {
  
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' required for parallel execution")
  }
  
  # Validate inputs
  validate_subgroup_inputs(subgroups_name, subgroups_id)
  
  # Define analyses
  analyses <- define_cox_analyses(wname)
  
  cat(sprintf("Running %d simulations in parallel using %d cores...\n", 
              sims, n_cores))
  
  # Run simulations in parallel
  cl <- parallel::makeCluster(n_cores)
  
  # Export necessary objects
  parallel::clusterExport(
    cl,
    varlist = c("dgm", "subgroups_id", "subgroups_name", "wname", "bw",
                "Ndraw", "bmcut", "analyses"),
    envir = environment()
  )
  
  # Load required packages on each worker
  parallel::clusterEvalQ(cl, {
    library(survival)
    library(data.table)
  })
  
  # Run simulations
  sim_results <- parallel::parLapply(cl, seq_len(sims), function(ss) {
    # Generate dataset
    df_sim <- generate_simulation_dataset(dgm, ss, Ndraw, wname, bw, bmcut)
    
    # Analyze all subgroups
    sg_results <- list()
    
    for (gg in seq_along(subgroups_id)) {
      df_sg <- subset_by_expression(df_sim, subgroups_id[gg])
      
      sg_results[[gg]] <- list(
        n = nrow(df_sg),
        hrs = numeric(length(analyses)),
        ubs = numeric(length(analyses))
      )
      
      for (i in seq_along(analyses)) {
        if (is_analysis_feasible(df_sg, analyses[[i]])) {
          fit_result <- fit_cox_model_safely(analyses[[i]]$formula_template, df_sg)
          
          if (!is.null(fit_result)) {
            sg_results[[gg]]$hrs[i] <- fit_result$hr
            sg_results[[gg]]$ubs[i] <- fit_result$upper_bound
          } else {
            sg_results[[gg]]$hrs[i] <- NA
            sg_results[[gg]]$ubs[i] <- NA
          }
        } else {
          sg_results[[gg]]$hrs[i] <- NA
          sg_results[[gg]]$ubs[i] <- NA
        }
      }
    }
    
    sg_results
  })
  
  parallel::stopCluster(cl)
  
  # Combine results
  results <- combine_parallel_results(sim_results, subgroups_name, analyses)
  
  # Save if requested
  if (!is.null(outfile)) {
    save_subgroup_results(results, dgm, subgroups_name, outfile)
  }
  
  return(results)
}

#' Combine results from parallel execution
#' @keywords internal
combine_parallel_results <- function(sim_results, subgroups_name, analyses) {
  sims <- length(sim_results)
  n_sg <- length(subgroups_name)
  
  results <- initialize_results_storage(sims, subgroups_name, analyses)
  
  for (ss in seq_len(sims)) {
    for (gg in seq_len(n_sg)) {
      results$ns[ss, gg] <- sim_results[[ss]][[gg]]$n
      
      for (i in seq_along(analyses)) {
        results[[paste0("hrs", i)]][ss, gg] <- sim_results[[ss]][[gg]]$hrs[i]
        results[[paste0("ubs", i)]][ss, gg] <- sim_results[[ss]][[gg]]$ubs[i]
      }
    }
  }
  
  return(results)
}

# ------------------------------------------------------------------------------
# EXAMPLE USAGE
# ------------------------------------------------------------------------------

#' Example workflow for subgroup analyses
#' @examples
#' \dontrun{
#' # Define subgroups
#' subgroups_id <- c(
#'   "itt == 'Y'",
#'   "bm_low == 0",
#'   "bm_low == 1",
#'   "meno == 0",
#'   "meno == 1"
#' )
#' 
#' subgroups_name <- c(
#'   "All Patients",
#'   "BM non-low",
#'   "BM low",
#'   "Pre-meno",
#'   "Post-meno"
#' )
#' 
#' # Run analyses
#' results <- get_SGanalyses(
#'   dgm = dgm,
#'   subgroups_name = subgroups_name,
#'   subgroups_id = subgroups_id,
#'   sims = 1000,
#'   wname = "meno",
#'   bw = -log(5),
#'   bmcut = log(2),
#'   outfile = "results/subgroup_analysis.Rdata"
#' )
#' 
#' # Summarize results
#' summary_df <- summarize_subgroup_results(results)
#' print(summary_df)
#' 
#' # Extract specific subgroup
#' postmeno_results <- extract_subgroup_results(results, "Post-meno")
#' }
#' @export
example_subgroup_workflow <- function() {
  message("See examples in function documentation")
}