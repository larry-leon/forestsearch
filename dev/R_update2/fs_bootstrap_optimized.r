#' ForestSearch Bootstrap with doFuture Parallelization
#'
#' Orchestrates bootstrap analysis for ForestSearch using doFuture parallelization.
#' This optimized version maintains all original functionality while improving
#' performance and reliability.
#'
#' @param fs.est ForestSearch results object.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @param details Logical. Print details during execution.
#' @param show_three Logical. Show details for first three bootstraps.
#' @param parallel_args List. Parallelization arguments (plan, workers, show_message).
#'
#' @return List with bootstrap results, confidence intervals, summary table, Ystar matrix, and estimates.
#'
#' @importFrom future plan
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @importFrom data.table data.table
#' @export

forestsearch_bootstrap_dofuture <- function(fs.est, 
                                           nb_boots, 
                                           details = FALSE, 
                                           show_three = FALSE,
                                           parallel_args = list()) {
  
  # ============================================================================
  # STEP 1: INPUT VALIDATION
  # ============================================================================
  
  if (!is.list(fs.est) || is.null(fs.est$args_call_all)) {
    stop("fs.est must be a valid ForestSearch results object with args_call_all")
  }
  
  if (!is.numeric(nb_boots) || nb_boots < 1) {
    stop("nb_boots must be a positive integer")
  }
  
  if (!is.logical(details) || !is.logical(show_three)) {
    stop("details and show_three must be logical values")
  }
  
  # ============================================================================
  # STEP 2: SETUP PARALLELIZATION
  # ============================================================================
  
  args_forestsearch_call <- fs.est$args_call_all
  
  # Configure parallel arguments
  if (length(parallel_args) == 0) {
    if (details) message("Using parallel plan from observed data analysis")
    parallel_args <- as.list(args_forestsearch_call$parallel_args)
    
    # Optimize worker allocation
    max_cores <- parallel::detectCores()
    if (details) message("Available cores: ", max_cores)
    
    # Use sensible defaults if not specified
    if (is.null(parallel_args$workers)) {
      parallel_args$workers <- min(nb_boots, max(1, max_cores - 1))
    }
  }
  
  # ============================================================================
  # STEP 3: ENSURE PACKAGES
  # ============================================================================
  
  ensure_packages(c("data.table", "foreach", "doFuture", "doRNG", "survival"))
  
  # ============================================================================
  # STEP 4: BUILD FORMULA AND FIT INITIAL MODELS
  # ============================================================================
  
  # Build Cox formula using only needed arguments
  cox.formula.boot <- build_cox_formula(
    outcome.name = args_forestsearch_call$outcome.name,
    event.name = args_forestsearch_call$event.name,
    treat.name = args_forestsearch_call$treat.name
  )
  
  # Fit initial Cox models for observed estimates
  cox_fits <- fit_cox_models(fs.est$df.est, cox.formula.boot)
  H_obs <- cox_fits$H_obs
  seH_obs <- cox_fits$seH_obs
  Hc_obs <- cox_fits$Hc_obs
  seHc_obs <- cox_fits$seHc_obs
  
  # ============================================================================
  # STEP 5: SETUP PARALLEL BACKEND
  # ============================================================================
  
  # Save and restore original plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  # Setup parallel processing
  setup_parallel_SGcons(parallel_args)
  
  # ============================================================================
  # STEP 6: GENERATE YSTAR MATRIX
  # ============================================================================
  
  if (details) cat("Generating Ystar matrix...\n")
  
  Ystar_mat <- bootstrap_ystar(fs.est$df.est, nb_boots)
  
  # Validate Ystar dimensions
  if (nrow(Ystar_mat) != nb_boots || ncol(Ystar_mat) != nrow(fs.est$df.est)) {
    stop("Ystar_mat dimensions incorrect: expected (", nb_boots, " x ", 
         nrow(fs.est$df.est), "), got (", nrow(Ystar_mat), " x ", 
         ncol(Ystar_mat), ")")
  }
  
  if (details) cat("Ystar matrix complete\n")
  
  # ============================================================================
  # STEP 7: RUN BOOTSTRAP
  # ============================================================================
  
  if (details) {
    cat("\nRunning bootstrap with ", nb_boots, " iterations...\n")
    if (show_three) cat("Showing details for first 3 iterations\n")
  }
  
  # Ensure add_id_column is available (define locally if needed)
  if (!exists("add_id_column", mode = "function")) {
    add_id_column <- function(df.analysis, id.name = NULL) {
      if (is.null(id.name)) {
        df.analysis$id <- seq_len(nrow(df.analysis))
        id.name <- "id"
      } else if (!(id.name %in% names(df.analysis))) {
        df.analysis[[id.name]] <- seq_len(nrow(df.analysis))
      }
      return(df.analysis)
    }
  }
  
  results <- bootstrap_results_optimized(
    fs.est = fs.est,
    df_boot_analysis = fs.est$df.est,
    cox.formula.boot = cox.formula.boot,
    nb_boots = nb_boots,
    show_three = show_three,
    H_obs = H_obs,
    Hc_obs = Hc_obs
  )
  
  # ============================================================================
  # STEP 8: CALCULATE CONFIDENCE INTERVALS
  # ============================================================================
  
  if (details) cat("\nCalculating confidence intervals...\n")
  
  est.scale <- args_forestsearch_call$est.scale
  
  # Check how many bootstrap iterations found subgroups
  n_found <- sum(!is.na(results$H_biasadj_2))
  
  if (details) {
    cat("Bootstrap iterations with subgroups found: ", n_found, " of ", nb_boots, "\n")
  }
  
  # Warn if too few subgroups were found
  if (n_found < nb_boots * 0.5) {
    warning("Only ", n_found, " of ", nb_boots, 
            " bootstrap iterations found subgroups. Results may be unreliable.")
  }
  
  # Calculate H estimates with error handling and NA checking
  H_estimates <- NULL
  if (n_found > 0 && !all(is.na(results$H_biasadj_1))) {
    H_estimates <- tryCatch(
      get_dfRes(
        Hobs = H_obs,
        seHobs = seH_obs,
        H1_adj = results$H_biasadj_1,
        H2_adj = results$H_biasadj_2,
        ystar = Ystar_mat,
        cov_method = "standard",
        cov_trim = 0.0,
        est.scale = est.scale,
        est.loghr = TRUE
      ),
      error = function(e) {
        if (details) cat("H estimates error: ", e$message, "\n")
        NULL
      }
    )
  }
  
  # Calculate Hc estimates with error handling and NA checking
  Hc_estimates <- NULL
  if (n_found > 0 && !all(is.na(results$Hc_biasadj_1))) {
    Hc_estimates <- tryCatch(
      get_dfRes(
        Hobs = Hc_obs,
        seHobs = seHc_obs,
        H1_adj = results$Hc_biasadj_1,
        H2_adj = results$Hc_biasadj_2,
        ystar = Ystar_mat,
        cov_method = "standard",
        cov_trim = 0.0,
        est.scale = est.scale,
        est.loghr = TRUE
      ),
      error = function(e) {
        if (details) cat("Hc estimates error: ", e$message, "\n")
        NULL
      }
    )
  }
  
  # ============================================================================
  # STEP 9: FORMAT RESULTS
  # ============================================================================
  
  # Check if estimates were successful
  if (is.null(H_estimates) || is.null(Hc_estimates)) {
    warning("Bootstrap estimation failed - returning partial results")
    return(list(
      results = results,
      SG_CIs = NULL,
      FSsg_tab = NULL,
      Ystar_mat = Ystar_mat,
      H_estimates = H_estimates,
      Hc_estimates = Hc_estimates
    ))
  }
  
  # Format confidence intervals
  H_res1 <- format_CI(H_estimates, c("H0", "H0_lower", "H0_upper"))
  H_res2 <- format_CI(H_estimates, c("H2", "H2_lower", "H2_upper"))
  Hc_res1 <- format_CI(Hc_estimates, c("H0", "H0_lower", "H0_upper"))
  Hc_res2 <- format_CI(Hc_estimates, c("H2", "H2_lower", "H2_upper"))
  
  # Print summary if requested
  if (details) {
    prop_found <- sum(!is.na(results$H_biasadj_2)) / nb_boots
    cat("\n")
    cat("Bootstrap Summary:\n")
    cat("  Proportion subgroups found: ", round(prop_found, 3), "\n")
    cat("  H un-adjusted:    ", H_res1, "\n")
    cat("  H bias-corrected: ", H_res2, "\n")
    cat("  H^c un-adjusted:  ", Hc_res1, "\n")
    cat("  H^c bias-corrected:", Hc_res2, "\n")
  }
  
  SG_CIs <- list(
    H_raw = H_res1,
    H_bc = H_res2,
    Hc_raw = Hc_res1,
    Hc_bc = Hc_res2
  )
  
  # ============================================================================
  # STEP 10: CREATE SUMMARY TABLE
  # ============================================================================
  
  # Prepare arguments for summary table
  sg_tab_args <- list(
    df = fs.est$df.est,
    SG_flag = "treat.recommend",
    draws = 0,
    details = FALSE,
    outcome.name = args_forestsearch_call$outcome.name,
    event.name = args_forestsearch_call$event.name,
    treat.name = args_forestsearch_call$treat.name,
    strata.name = NULL,
    potentialOutcome.name = args_forestsearch_call$potentialOutcome.name,
    est.scale = est.scale
  )
  
  # Add adjusted CIs based on scale
  if (est.scale == "1/hr") {
    sg_tab_args$hr_1a <- SG_CIs$H_bc
    sg_tab_args$hr_0a <- SG_CIs$Hc_bc
    sg_tab_args$sg0_name <- "Questionable"
    sg_tab_args$sg1_name <- "Recommend"
  } else {
    sg_tab_args$hr_1a <- SG_CIs$Hc_bc
    sg_tab_args$hr_0a <- SG_CIs$H_bc
    sg_tab_args$sg0_name <- "Questionable"
    sg_tab_args$sg1_name <- "Recommend"
  }
  
  # Generate summary table
  FSsg_tab <- do.call(SG_tab_estimates, sg_tab_args)
  
  # ============================================================================
  # STEP 11: RETURN RESULTS
  # ============================================================================
  
  list(
    results = results,
    SG_CIs = SG_CIs,
    FSsg_tab = FSsg_tab,
    Ystar_mat = Ystar_mat,
    H_estimates = H_estimates,
    Hc_estimates = Hc_estimates
  )
}

#' Optimized Bootstrap Results for ForestSearch
#'
#' Runs bootstrap analysis with improved efficiency and error handling.
#' Key optimizations:
#' - Pre-allocated result structures
#' - Simplified argument passing
#' - Better memory management
#' - Cleaner error handling
#'
#' @inheritParams bootstrap_results
#' @return Data.table with bias-adjusted estimates and search metrics.
#' @importFrom foreach foreach
#' @importFrom data.table data.table
#' @importFrom doFuture %dofuture%
#' @keywords internal

bootstrap_results_optimized <- function(fs.est, df_boot_analysis, cox.formula.boot, 
                                       nb_boots, show_three, H_obs, Hc_obs) {
  
  NN <- nrow(df_boot_analysis)
  
  # Extract ForestSearch arguments once
  args_FS_base <- fs.est$args_call_all
  
  # Pre-process arguments for efficiency
  # Remove arguments that will be set per bootstrap
  args_to_remove <- c("df.analysis", "df.predict", "details", 
                      "showten_subgroups", "plot.sg", "plot.grf",
                      "grf_res", "grf_cuts", "parallel_args")
  
  args_FS_clean <- args_FS_base[!names(args_FS_base) %in% args_to_remove]
  
  # Set static arguments
  args_FS_clean$details <- FALSE
  args_FS_clean$showten_subgroups <- FALSE
  args_FS_clean$plot.sg <- FALSE
  args_FS_clean$plot.grf <- FALSE
  args_FS_clean$grf_res <- NULL
  args_FS_clean$grf_cuts <- NULL
  
  # Force sequential processing inside bootstrap
  args_FS_clean$parallel_args <- list(
    plan = "sequential",
    workers = 1,
    show_message = FALSE
  )
  
  # Prepare columns to drop
  drop.vars <- c(fs.est$confounders.candidate, "treat.recommend")
  keep.vars <- setdiff(names(df_boot_analysis), drop.vars)
  
  # Run bootstrap iterations
  foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(
      seed = TRUE,
      # Add all necessary functions to worker environment
      add = c(
        # Core forestsearch functions
        "forestsearch", "get_FSdata", "subgroup.search", 
        "subgroup.consistency", "get_Cox_sg",
        
        # Helper functions  
        "dummy", "dummy2", "acm.disjctif", "is.continuous",
        "get_conf_force", "lasso_selection", "FS_labels",
        "get_dfpred", "extract_subgroup", "sort_subgroups",
        
        # Bootstrap specific (note the comma above)
        "count.id", "calc_cov", "ci_est", "get_targetEst"
      )
    ),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {
    
    # Set seed for reproducibility
    set.seed(8316951 + boot * 100)
    
    # Show details for first 3 if requested
    show_details <- show_three && (boot <= 3)
    if (show_details) {
      cat("Bootstrap iteration", boot, "\n")
    }
    
    # Generate bootstrap sample
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df_boot_analysis[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))
    
    # Calculate bootstrap estimates at observed subgroups
    H_star <- tryCatch({
      fit <- get_Cox_sg(
        df_sg = subset(df_boot, treat.recommend == 0),
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      fit$est_obs
    }, error = function(e) NA_real_)
    
    Hc_star <- tryCatch({
      fit <- get_Cox_sg(
        df_sg = subset(df_boot, treat.recommend == 1),
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      fit$est_obs
    }, error = function(e) NA_real_)
    
    # Initialize results
    H_biasadj_1 <- H_biasadj_2 <- NA_real_
    Hc_biasadj_1 <- Hc_biasadj_2 <- NA_real_
    tmins_search <- max_sg_est <- prop_maxk <- L <- max_count <- NA_real_
    
    # Prepare data for ForestSearch
    dfnew_boot <- df_boot[, keep.vars, drop = FALSE]
    dfnew <- df_boot_analysis[, keep.vars, drop = FALSE]
    
    # Setup arguments for this bootstrap
    args_FS_boot <- args_FS_clean
    args_FS_boot$df.analysis <- dfnew_boot
    args_FS_boot$df.predict <- dfnew
    args_FS_boot$details <- show_details
    
    # Run ForestSearch with error handling
    run_bootstrap <- tryCatch(
      do.call(forestsearch, args_FS_boot),
      error = function(e) {
        if (show_details) {
          cat("Bootstrap", boot, "failed:", e$message, "\n")
        }
        NULL
      }
    )
    
    # Process results if successful
    if (!is.null(run_bootstrap) && !is.null(run_bootstrap$sg.harm)) {
      
      # Extract results
      df_PredBoot <- run_bootstrap$df.predict
      dfboot_PredBoot <- run_bootstrap$df.est
      
      # Store metrics
      max_sg_est <- as.numeric(run_bootstrap$find.grps$max_sg_est)
      tmins_search <- as.numeric(run_bootstrap$find.grps$time_search)
      prop_maxk <- as.numeric(run_bootstrap$prop_maxk)
      max_count <- run_bootstrap$find.grps$max_count
      L <- run_bootstrap$find.grps$L
      
      # Calculate bias adjustments
      # H subgroup (treat.recommend == 0)
      Hstar_obs <- tryCatch({
        fit <- get_Cox_sg(
          df_sg = subset(df_PredBoot, treat.recommend == 0),
          cox.formula = cox.formula.boot,
          est.loghr = TRUE
        )
        fit$est_obs
      }, error = function(e) NA_real_)
      
      Hstar_star <- tryCatch({
        fit <- get_Cox_sg(
          df_sg = subset(dfboot_PredBoot, treat.recommend == 0),
          cox.formula = cox.formula.boot,
          est.loghr = TRUE
        )
        fit$est_obs
      }, error = function(e) NA_real_)
      
      if (!is.na(Hstar_obs) && !is.na(Hstar_star)) {
        H_biasadj_1 <- H_obs - (Hstar_star - Hstar_obs)
        H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)
      }
      
      # Hc subgroup (treat.recommend == 1)
      Hcstar_obs <- tryCatch({
        fit <- get_Cox_sg(
          df_sg = subset(df_PredBoot, treat.recommend == 1),
          cox.formula = cox.formula.boot,
          est.loghr = TRUE
        )
        fit$est_obs
      }, error = function(e) NA_real_)
      
      Hcstar_star <- tryCatch({
        fit <- get_Cox_sg(
          df_sg = subset(dfboot_PredBoot, treat.recommend == 1),
          cox.formula = cox.formula.boot,
          est.loghr = TRUE
        )
        fit$est_obs
      }, error = function(e) NA_real_)
      
      if (!is.na(Hcstar_obs) && !is.na(Hcstar_star)) {
        Hc_biasadj_1 <- Hc_obs - (Hcstar_star - Hcstar_obs)
        Hc_biasadj_2 <- 2 * Hc_obs - (Hc_star + Hcstar_star - Hcstar_obs)
      }
      
      if (show_details) {
        cat("  Subgroup found in bootstrap", boot, "\n")
        cat("  H bias adjustments:", round(c(H_biasadj_1, H_biasadj_2), 4), "\n")
        cat("  Hc bias adjustments:", round(c(Hc_biasadj_1, Hc_biasadj_2), 4), "\n")
      }
    }
    
    # Return results as data.table row
    data.table::data.table(
      H_biasadj_1 = H_biasadj_1,
      H_biasadj_2 = H_biasadj_2,
      Hc_biasadj_1 = Hc_biasadj_1,
      Hc_biasadj_2 = Hc_biasadj_2,
      tmins_search = tmins_search,
      max_sg_est = max_sg_est,
      prop_maxk = prop_maxk,
      L = L,
      max_count = max_count
    )
  }
}

# ============================================================================
# BACKWARD COMPATIBILITY WRAPPER
# ============================================================================

#' Bootstrap Results for ForestSearch
#'
#' Runs bootstrap analysis for ForestSearch, fitting Cox models and performing
#' bias correction for subgroup treatment effect estimates. This function 
#' performs the core bootstrap iterations, re-running ForestSearch on each 
#' bootstrap sample to assess the stability of subgroup identification and 
#' calculate bias-corrected confidence intervals.
#'
#' @param fs.est ForestSearch results object from \code{\link{forestsearch}}.
#'   Must contain:
#'   \itemize{
#'     \item \code{df.est}: Data frame with treatment recommendations
#'     \item \code{args_call_all}: Original ForestSearch arguments
#'     \item \code{confounders.candidate}: Candidate confounder names
#'     \item \code{find.grps}: Subgroup search results
#'   }
#' @param df_boot_analysis Data frame for bootstrap analysis. Typically the same
#'   as \code{fs.est$df.est} but can be a different dataset for external validation.
#' @param cox.formula.boot Formula object for Cox regression in bootstrap.
#'   Usually created by \code{\link{build_cox_formula}} with outcome, event,
#'   and treatment variable names.
#' @param nb_boots Integer. Number of bootstrap samples to generate. Recommended
#'   minimum is 100, with 500-1000 for final analyses. Higher values give more
#'   stable estimates but increase computation time.
#' @param show_three Logical. If TRUE, displays detailed output for the first
#'   three bootstrap iterations to help with debugging and understanding the
#'   process. Default is FALSE.
#' @param H_obs Numeric. Observed hazard ratio (or log hazard ratio if on log scale)
#'   for subgroup H (treatment not recommended group) from the original analysis.
#' @param Hc_obs Numeric. Observed hazard ratio (or log hazard ratio if on log scale)
#'   for subgroup H^c (treatment recommended group) from the original analysis.
#'
#' @return A data.table with the following columns for each bootstrap iteration:
#' \describe{
#'   \item{H_biasadj_1}{First bias-adjusted estimate for subgroup H using 
#'     bootstrap bias correction method 1 (simple bias correction)}
#'   \item{H_biasadj_2}{Second bias-adjusted estimate for subgroup H using 
#'     bootstrap bias correction method 2 (double bootstrap correction)}
#'   \item{Hc_biasadj_1}{First bias-adjusted estimate for subgroup H^c}
#'   \item{Hc_biasadj_2}{Second bias-adjusted estimate for subgroup H^c}
#'   \item{tmins_search}{Time in minutes for subgroup search in this bootstrap}
#'   \item{max_sg_est}{Maximum subgroup effect estimate found}
#'   \item{prop_maxk}{Proportion of maximum k-factor combinations evaluated}
#'   \item{L}{Number of candidate variables in subgroup search}
#'   \item{max_count}{Maximum number of subgroup combinations considered}
#' }
#'
#' @details
#' The bootstrap procedure works as follows:
#' \enumerate{
#'   \item Generate bootstrap sample with replacement from original data
#'   \item Calculate HR estimates at original subgroups (H* and Hc*)
#'   \item Re-run ForestSearch on bootstrap sample to find new subgroups
#'   \item Calculate bias adjustments using the double bootstrap method
#'   \item Store results for confidence interval calculation
#' }
#'
#' Two bias correction methods are implemented:
#' \itemize{
#'   \item \strong{Method 1}: Simple bias correction: \code{H_obs - (H** - H*obs)}
#'   \item \strong{Method 2}: Double bootstrap: \code{2*H_obs - (H* + H** - H*obs)}
#' }
#'
#' Where:
#' \itemize{
#'   \item \code{H_obs}: Original estimate from full data
#'   \item \code{H*}: Bootstrap estimate at original subgroup
#'   \item \code{H*obs}: Original data evaluated at bootstrap subgroup
#'   \item \code{H**}: Bootstrap data evaluated at bootstrap subgroup
#' }
#'
#' @note 
#' \itemize{
#'   \item This function is called internally by \code{\link{forestsearch_bootstrap_dofuture}}
#'   \item Parallelization is handled by the parent function
#'   \item NA values are returned for bootstrap iterations that don't find subgroups
#'   \item At least 50\% of bootstrap iterations should find subgroups for reliable results
#' }
#'
#' @seealso 
#' \code{\link{forestsearch_bootstrap_dofuture}} for the main bootstrap interface,
#' \code{\link{forestsearch}} for the subgroup identification method,
#' \code{\link{get_dfRes}} for confidence interval calculation,
#' \code{\link{get_Cox_sg}} for Cox model fitting
#'
#' @references
#' Efron, B. and Tibshirani, R.J. (1993). An Introduction to the Bootstrap.
#' Chapman & Hall/CRC.
#'
#' @examples
#' \dontrun{
#' # This function is typically called by forestsearch_bootstrap_dofuture
#' # but can be used directly for custom bootstrap analyses
#' 
#' # Prepare Cox formula
#' cox.formula <- build_cox_formula(
#'   outcome.name = "time",
#'   event.name = "event", 
#'   treat.name = "treatment"
#' )
#' 
#' # Fit initial models
#' cox_fits <- fit_cox_models(fs.est$df.est, cox.formula)
#' 
#' # Run bootstrap (usually done in parallel by parent function)
#' results <- bootstrap_results(
#'   fs.est = fs_result,
#'   df_boot_analysis = fs_result$df.est,
#'   cox.formula.boot = cox.formula,
#'   nb_boots = 100,
#'   show_three = TRUE,
#'   H_obs = cox_fits$H_obs,
#'   Hc_obs = cox_fits$Hc_obs
#' )
#' 
#' # Check proportion of successful bootstraps
#' prop_found <- sum(!is.na(results$H_biasadj_2)) / nrow(results)
#' cat("Proportion finding subgroups:", prop_found, "\n")
#' }
#'
#' @importFrom data.table data.table
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @export

bootstrap_results <- function(fs.est, df_boot_analysis, cox.formula.boot, nb_boots,
                              show_three, H_obs, Hc_obs) {
  # Call the optimized version
  bootstrap_results_optimized(
    fs.est = fs.est,
    df_boot_analysis = df_boot_analysis,
    cox.formula.boot = cox.formula.boot,
    nb_boots = nb_boots,
    show_three = show_three,
    H_obs = H_obs,
    Hc_obs = Hc_obs
  )
}

# ============================================================================
# ORIGINAL HELPER FUNCTIONS (Maintained for compatibility)
# ============================================================================

# Note: The following helper functions are used by forestsearch_bootstrap_dofuture
# and are kept exactly as in the original implementation to ensure compatibility:
#
# - bootstrap_ystar: Generates bootstrap Ystar matrix
# - ensure_packages: Checks and loads required packages  
# - build_cox_formula: Builds Cox model formula
# - fit_cox_models: Fits Cox models for subgroups
# - format_CI: Formats confidence intervals
# - count.id: Counts ID occurrences in bootstrap sample
# - ci_est: Calculates confidence interval for estimate
# - get_Cox_sg: Fits Cox model for subgroup
#
# These functions are defined in forestsearch_bootstrap_helpers.R
# and should be loaded from that file to avoid duplication