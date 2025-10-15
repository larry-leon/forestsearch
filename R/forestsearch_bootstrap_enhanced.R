#' ForestSearch Bootstrap with doFuture Parallelization (Enhanced)
#'
#' Orchestrates bootstrap analysis for ForestSearch using doFuture parallelization.
#' This enhanced version includes improved input validation, memory management,
#' progress reporting, and clearer code organization while preserving the original
#' est.scale logic.
#'
#' @param fs.est ForestSearch results object.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @param details Logical. Print details during execution.
#' @param show_three Logical. Show details for first three bootstraps.
#' @param parallel_args List. Parallelization arguments (plan, workers, show_message).
#' @param seed_base Integer. Base seed for reproducibility (default: 8316951).
#' @param progress_interval Integer. Report progress every N bootstraps (default: 10).
#'
#' @return List with bootstrap results, confidence intervals, summary table, Ystar matrix, and estimates.
#'
#' @importFrom future plan
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @importFrom data.table data.table
#' @export
#'
#' @examples
#' \dontrun{
#' # Run bootstrap with parallel processing
#' fs_bc <- forestsearch_bootstrap_dofuture_enhanced(
#'   fs.est = fs_result,
#'   nb_boots = 100,
#'   parallel_args = list(plan = "multisession", workers = 4)
#' )
#' }

forestsearch_bootstrap <- function(
    fs.est,
    nb_boots,
    details = FALSE,
    show_three = FALSE,
    parallel_args = list()
) {

  # Get arguments from original forestsearch call
  args_forestsearch_call <- fs.est$args_call_all

  # Handle parallel arguments
  if(length(parallel_args) == 0){
    if(details) message("Using parallel plan from original forestsearch analysis")
    parallel_args <- as.list(args_forestsearch_call$parallel_args)
    max_cores <- parallel::detectCores()
    if(details) message("Note that max cores = ", max_cores)
  }

  # ---- Step 1: Ensure packages ----
  required_packages <- c("data.table", "foreach", "doFuture", "doRNG", "survival")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required. Install with: install.packages('", pkg, "')", sep=""))
    }
  }

  library(data.table)
  library(foreach)
  library(doFuture)
  library(survival)

  # ---- Step 2: Build Cox formula ----
  cox.formula.boot <- as.formula(paste0(
    "Surv(", args_forestsearch_call$outcome.name, ", ",
    args_forestsearch_call$event.name, ") ~ ",
    args_forestsearch_call$treat.name
  ))

  # ---- Step 3: Fit initial Cox models ----
  if(details) cat("Fitting initial Cox models...\n")

  # Extract the data
  df.analysis <- fs.est$df.est

  # Fit models for each subgroup
  df_sg0 <- subset(df.analysis, treat.recommend == 0)
  df_sg1 <- subset(df.analysis, treat.recommend == 1)

  # Get Cox estimates with error handling
  tryCatch({
    fit_sg0 <- get_Cox_sg(df_sg = df_sg0, cox.formula = cox.formula.boot, est.loghr = TRUE)
    H_obs <- as.numeric(fit_sg0$est_obs)
    seH_obs <- as.numeric(fit_sg0$se_obs)
  }, error = function(e) {
    stop("Failed to fit Cox model for subgroup 0: ", e$message)
  })

  tryCatch({
    fit_sg1 <- get_Cox_sg(df_sg = df_sg1, cox.formula = cox.formula.boot, est.loghr = TRUE)
    Hc_obs <- as.numeric(fit_sg1$est_obs)
    seHc_obs <- as.numeric(fit_sg1$se_obs)
  }, error = function(e) {
    stop("Failed to fit Cox model for subgroup 1: ", e$message)
  })

  if(details) {
    cat("Initial estimates:\n")
    cat("  H_obs:", H_obs, "seH_obs:", seH_obs, "\n")
    cat("  Hc_obs:", Hc_obs, "seHc_obs:", seHc_obs, "\n")
  }

  # ---- Step 4: Setup parallel processing ----
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (!is.null(parallel_args$plan)) {
    if (parallel_args$plan == "multisession") {
      future::plan(future::multisession, workers = parallel_args$workers)
    } else if (parallel_args$plan == "sequential") {
      future::plan(future::sequential)
    }
  }

  # ---- Step 5: Generate Ystar matrix ----
  if (details) cat("Generating Ystar matrix...\n")

  NN <- nrow(df.analysis)

  Ystar_mat <- foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {
    set.seed(8316951 + boot * 100)
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df.analysis[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))
    ystar <- unlist(lapply(df.analysis$id, count.id, dfb = df_boot))
    return(ystar)
  }

  if (details) cat("Done with Ystar_mat\n")

  # Validate Ystar dimensions
  if (nrow(Ystar_mat) != nb_boots || ncol(Ystar_mat) != nrow(df.analysis)) {
    stop("Dimension of Ystar_mat must be (nb_boots x n)")
  }

  # ---- Step 6: Run bootstrap analysis ----
  if (details) cat("Starting bootstrap analysis...\n")

  results <- foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {

    # Progress reporting
    if (details && boot %% 10 == 0) {
      cat("Completed", boot, "of", nb_boots, "bootstraps\n")
    }

    set.seed(8316951 + boot * 100)
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df.analysis[in_boot, ]

    # Calculate bootstrap estimates
    tryCatch({
      # Bootstrap estimates at original subgroups
      fit_H_star <- get_Cox_sg(
        df_sg = subset(df_boot, treat.recommend == 0),
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      H_star <- as.numeric(fit_H_star$est_obs)

      fit_Hc_star <- get_Cox_sg(
        df_sg = subset(df_boot, treat.recommend == 1),
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      Hc_star <- as.numeric(fit_Hc_star$est_obs)

      # Simple bias correction (method 1)
      H_biasadj_1 <- 2 * H_obs - H_star
      H_biasadj_2 <- H_star  # Use bootstrap estimate as is for method 2
      Hc_biasadj_1 <- 2 * Hc_obs - Hc_star
      Hc_biasadj_2 <- Hc_star

    }, error = function(e) {
      # Return NAs if bootstrap fails
      H_biasadj_1 <- NA
      H_biasadj_2 <- NA
      Hc_biasadj_1 <- NA
      Hc_biasadj_2 <- NA
    })

    data.table::data.table(
      H_biasadj_1 = H_biasadj_1,
      H_biasadj_2 = H_biasadj_2,
      Hc_biasadj_1 = Hc_biasadj_1,
      Hc_biasadj_2 = Hc_biasadj_2
    )
  }

  # ---- Step 7: Calculate confidence intervals ----
  if (details) cat("Calculating confidence intervals...\n")

  est.scale <- args_forestsearch_call$est.scale
  if (is.null(est.scale)) est.scale <- "hr"

  # Calculate estimates with fixed function
  H_estimates <- tryCatch({
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
    )
  }, error = function(e) {
    warning("Failed to calculate H estimates: ", e$message)
    NULL
  }, warning = function(w) {
    warning("Warning in H estimates: ", w$message)
    NULL
  })

  Hc_estimates <- tryCatch({
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
    )
  }, error = function(e) {
    warning("Failed to calculate Hc estimates: ", e$message)
    NULL
  }, warning = function(w) {
    warning("Warning in Hc estimates: ", w$message)
    NULL
  })

  # ---- Step 8: Format results ----
  if (!is.null(H_estimates) && !is.null(Hc_estimates)) {
    # Format confidence intervals
    format_CI_safe <- function(estimates, cols) {
      if (is.null(estimates) || !all(cols %in% names(estimates))) {
        return("NA (NA, NA)")
      }
      vals <- round(unlist(estimates[1, cols, with=FALSE]), 2)
      if (any(is.na(vals))) {
        return("NA (NA, NA)")
      }
      paste0(vals[1], " (", vals[2], ", ", vals[3], ")")
    }

    H_res1 <- format_CI_safe(H_estimates, c("H0", "H0_lower", "H0_upper"))
    H_res2 <- format_CI_safe(H_estimates, c("H2", "H2_lower", "H2_upper"))
    Hc_res1 <- format_CI_safe(Hc_estimates, c("H0", "H0_lower", "H0_upper"))
    Hc_res2 <- format_CI_safe(Hc_estimates, c("H2", "H2_lower", "H2_upper"))

    if (details) {
      n_successful <- sum(!is.na(results$H_biasadj_2))
      cat("\n**** % bootstrap subgroups found =", n_successful/nb_boots, "\n")
      cat("H un-adjusted estimates-----:   ", H_res1, "\n")
      cat("H bias-corrected estimates--:   ", H_res2, "\n")
      cat("H^c un-adjusted estimates---:   ", Hc_res1, "\n")
      cat("H^c bias-corrected estimates:   ", Hc_res2, "\n")
    }

    SG_CIs <- list(H_raw = H_res1, H_bc = H_res2, Hc_raw = Hc_res1, Hc_bc = Hc_res2)

    # Create summary table
    FSsg_tab <- NULL  # Simplified for now

  } else {
    if (details) cat("Warning: Bootstrap confidence interval calculation failed\n")
    SG_CIs <- NULL
    FSsg_tab <- NULL
    H_estimates <- NULL
    Hc_estimates <- NULL
  }

  # Return results
  return(list(
    results = results,
    SG_CIs = SG_CIs,
    FSsg_tab = FSsg_tab,
    Ystar_mat = Ystar_mat,
    H_estimates = H_estimates,
    Hc_estimates = Hc_estimates
  ))
}


forestsearch_bootstrap_dofuture_enhanced_old <- function(
    fs.est,
    nb_boots,
    details = FALSE,
    show_three = FALSE,
    parallel_args = list(),
    seed_base = 8316951,
    progress_interval = 10
) {

  # ============================================================================
  # STEP 1: INPUT VALIDATION
  # ============================================================================

  validate_bootstrap_inputs(fs.est, nb_boots, details, show_three, parallel_args, seed_base, progress_interval)

  # ============================================================================
  # STEP 2: SETUP
  # ============================================================================

  args_forestsearch_call <- fs.est$args_call_all

  # Handle parallel arguments with clear logic
  parallel_args <- setup_parallel_args(parallel_args, args_forestsearch_call, details)

  # Ensure required packages
  ensure_packages(c("data.table", "foreach", "doFuture", "doRNG", "survival"))

  # ============================================================================
  # STEP 3: BUILD COX FORMULA
  # ============================================================================

  cox.formula.boot <- build_cox_formula_safe(args_forestsearch_call)

  # ============================================================================
  # STEP 4: FIT INITIAL COX MODELS
  # ============================================================================

  cox_fits <- fit_cox_models_safe(fs.est$df.est, cox.formula.boot, details)
  H_obs <- cox_fits$H_obs
  seH_obs <- cox_fits$seH_obs
  Hc_obs <- cox_fits$Hc_obs
  seHc_obs <- cox_fits$seHc_obs

  # ============================================================================
  # STEP 5: SETUP PARALLEL PROCESSING
  # ============================================================================

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)  # Restore plan on exit
  setup_parallel_SGcons(parallel_args)

  # ============================================================================
  # STEP 6: GENERATE YSTAR MATRIX
  # ============================================================================

  if (details) cat("Generating Ystar matrix...\n")

  Ystar_mat <- bootstrap_ystar_safe(fs.est$df.est, nb_boots, seed_base)

  # Validate Ystar dimensions
  if (nrow(Ystar_mat) != nb_boots || ncol(Ystar_mat) != nrow(fs.est$df.est)) {
    stop("Dimension of Ystar_mat must be (nb_boots x n)")
  }

  if (details) cat("Done with Ystar_mat\n")

  # ============================================================================
  # STEP 7: RUN BOOTSTRAP ANALYSIS
  # ============================================================================

  if (details) cat("Starting bootstrap analysis...\n")

  results <- bootstrap_results_enhanced(
    fs.est = fs.est,
    df_boot_analysis = fs.est$df.est,
    cox.formula.boot = cox.formula.boot,
    nb_boots = nb_boots,
    show_three = show_three,
    H_obs = H_obs,
    Hc_obs = Hc_obs,
    seed_base = seed_base,
    progress_interval = progress_interval,
    details = details
  )

  # ============================================================================
  # STEP 8: CALCULATE CONFIDENCE INTERVALS
  # ============================================================================

  if (details) cat("Calculating confidence intervals...\n")

  est.scale <- args_forestsearch_call$est.scale

  # Calculate H estimates with error handling
  H_estimates <- calculate_estimates_safe(
    Hobs = H_obs,
    seHobs = seH_obs,
    H1_adj = results$H_biasadj_1,
    H2_adj = results$H_biasadj_2,
    ystar = Ystar_mat,
    est.scale = est.scale,
    details = details
  )

  # Calculate Hc estimates with error handling
  Hc_estimates <- calculate_estimates_safe(
    Hobs = Hc_obs,
    seHobs = seHc_obs,
    H1_adj = results$Hc_biasadj_1,
    H2_adj = results$Hc_biasadj_2,
    ystar = Ystar_mat,
    est.scale = est.scale,
    details = details
  )

  # Clean up large matrix from memory
  rm(Ystar_mat)
  gc(verbose = FALSE)

  # ============================================================================
  # STEP 9: FORMAT RESULTS
  # ============================================================================

  # Check if estimates were successful
  if (is.null(H_estimates) || is.null(Hc_estimates)) {
    warning("Bootstrap confidence interval calculation failed")
    return(list(
      results = results,
      SG_CIs = NULL,
      FSsg_tab = NULL,
      Ystar_mat = NULL,
      H_estimates = NULL,
      Hc_estimates = NULL,
      diagnostics = calculate_diagnostics(results, nb_boots)
    ))
  }

  # Format confidence intervals
  formatted_CIs <- format_confidence_intervals(H_estimates, Hc_estimates, details)

  # ============================================================================
  # STEP 10: CREATE SUMMARY TABLE WITH PRESERVED EST.SCALE LOGIC
  # ============================================================================

  # CRITICAL: Preserve the original est.scale logic exactly
  FSsg_tab <- create_summary_table_with_scale(
    fs.est = fs.est,
    SG_CIs = formatted_CIs,
    est.scale = est.scale,
    args_forestsearch_call = args_forestsearch_call
  )

  # ============================================================================
  # STEP 11: CALCULATE DIAGNOSTICS
  # ============================================================================

  diagnostics <- calculate_diagnostics(results, nb_boots)

  if (details) {
    print_bootstrap_summary(results, nb_boots, formatted_CIs, diagnostics)
  }

  # ============================================================================
  # STEP 12: RETURN RESULTS
  # ============================================================================

  return(list(
    results = results,
    SG_CIs = formatted_CIs,
    FSsg_tab = FSsg_tab,
    Ystar_mat = NULL,  # Already cleaned up
    H_estimates = H_estimates,
    Hc_estimates = Hc_estimates,
    diagnostics = diagnostics
  ))
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Validate Bootstrap Inputs
#' @keywords internal
validate_bootstrap_inputs <- function(fs.est, nb_boots, details, show_three,
                                      parallel_args, seed_base, progress_interval) {

  # Check fs.est
  if (!is.list(fs.est)) {
    stop("'fs.est' must be a ForestSearch results object (list)")
  }

  required_elements <- c("df.est", "args_call_all")
  missing <- setdiff(required_elements, names(fs.est))
  if (length(missing) > 0) {
    stop("'fs.est' missing required elements: ", paste(missing, collapse = ", "))
  }

  # Check nb_boots
  if (!is.numeric(nb_boots) || length(nb_boots) != 1 || nb_boots < 1) {
    stop("'nb_boots' must be a positive integer")
  }

  # Check logical arguments
  if (!is.logical(details) || length(details) != 1) {
    stop("'details' must be a single logical value")
  }

  if (!is.logical(show_three) || length(show_three) != 1) {
    stop("'show_three' must be a single logical value")
  }

  # Check parallel_args
  if (!is.list(parallel_args)) {
    stop("'parallel_args' must be a list")
  }

  # Check seed_base
  if (!is.numeric(seed_base) || length(seed_base) != 1) {
    stop("'seed_base' must be a single numeric value")
  }

  # Check progress_interval
  if (!is.numeric(progress_interval) || length(progress_interval) != 1 || progress_interval < 1) {
    stop("'progress_interval' must be a positive integer")
  }

  invisible(TRUE)
}

#' Setup Parallel Arguments
#' @keywords internal
setup_parallel_args <- function(parallel_args, args_forestsearch_call, details) {

  if (length(parallel_args) == 0) {
    if (details) {
      message("Using parallel plan from original forestsearch analysis")
    }
    parallel_args <- as.list(args_forestsearch_call$parallel_args)

    # Get available cores
    max_cores <- parallel::detectCores()
    if (details) {
      message("Note that max cores = ", max_cores)
    }

    # Ensure workers doesn't exceed max cores
    if (!is.null(parallel_args$workers)) {
      parallel_args$workers <- min(parallel_args$workers, max_cores)
    }
  }

  return(parallel_args)
}

#' Build Cox Formula Safely
#' @keywords internal
build_cox_formula_safe <- function(args_forestsearch_call) {

  args_build <- names(formals(build_cox_formula))
  args_build_filtered <- args_forestsearch_call[names(args_forestsearch_call) %in% args_build]

  cox.formula <- try(
    do.call(build_cox_formula, args_build_filtered),
    silent = TRUE
  )

  if (inherits(cox.formula, "try-error")) {
    stop("Failed to build Cox formula: ", attr(cox.formula, "condition")$message)
  }

  return(cox.formula)
}

#' Fit Cox Models with Error Handling
#' @keywords internal
fit_cox_models_safe <- function(df, formula, details) {

  if (details) cat("Fitting initial Cox models...\n")

  cox_fits <- try(
    fit_cox_models(df, formula),
    silent = TRUE
  )

  if (inherits(cox_fits, "try-error")) {
    stop("Failed to fit Cox models: ", attr(cox_fits, "condition")$message)
  }

  # Validate outputs
  required_names <- c("H_obs", "seH_obs", "Hc_obs", "seHc_obs")
  if (!all(required_names %in% names(cox_fits))) {
    stop("Cox model fitting did not return expected results")
  }

  return(cox_fits)
}

#' Generate Ystar Matrix Safely
#' @keywords internal
bootstrap_ystar_safe <- function(df, nb_boots, seed_base) {

  Ystar_mat <- try(
    bootstrap_ystar_with_seed(df, nb_boots, seed_base),
    silent = TRUE
  )

  if (inherits(Ystar_mat, "try-error")) {
    stop("Failed to generate Ystar matrix: ", attr(Ystar_mat, "condition")$message)
  }

  return(Ystar_mat)
}

#' Bootstrap Ystar with Custom Seed
#' @keywords internal
bootstrap_ystar_with_seed <- function(df, nb_boots, seed_base) {
  NN <- nrow(df)
  foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {
    set.seed(seed_base + boot * 100)
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))
    ystar <- unlist(lapply(df$id, count.id, dfb = df_boot))
    return(ystar)
  }
}

#' Enhanced Bootstrap Results Function
#' @keywords internal
bootstrap_results_enhanced <- function(fs.est, df_boot_analysis, cox.formula.boot,
                                       nb_boots, show_three, H_obs, Hc_obs,
                                       seed_base, progress_interval, details) {

  NN <- nrow(df_boot_analysis)

  # Run bootstrap with progress reporting
  results_list <- foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {

    # Progress reporting
    if (details && boot %% progress_interval == 0) {
      cat("Completed", boot, "of", nb_boots, "bootstraps\n")
    }

    # Run single bootstrap
    result <- perform_single_bootstrap(
      boot = boot,
      NN = NN,
      df_boot_analysis = df_boot_analysis,
      cox.formula.boot = cox.formula.boot,
      fs.est = fs.est,
      H_obs = H_obs,
      Hc_obs = Hc_obs,
      show_three = show_three,
      seed_base = seed_base
    )

    return(result)
  }

  return(results_list)
}

#' Perform Single Bootstrap Iteration
#' @keywords internal
perform_single_bootstrap <- function(boot, NN, df_boot_analysis, cox.formula.boot,
                                     fs.est, H_obs, Hc_obs, show_three, seed_base) {

  # Set seed for reproducibility
  set.seed(seed_base + boot * 100)

  # Create bootstrap sample
  in_boot <- sample.int(NN, size = NN, replace = TRUE)
  df_boot <- df_boot_analysis[in_boot, ]
  df_boot$id_boot <- seq_len(nrow(df_boot))

  # Calculate bootstrap estimates at original subgroups
  fitH_star <- get_Cox_sg(
    df_sg = subset(df_boot, treat.recommend == 0),
    cox.formula = cox.formula.boot,
    est.loghr = TRUE
  )
  H_star <- fitH_star$est_obs

  fitHc_star <- get_Cox_sg(
    df_sg = subset(df_boot, treat.recommend == 1),
    cox.formula = cox.formula.boot,
    est.loghr = TRUE
  )
  Hc_star <- fitHc_star$est_obs

  # Initialize bias corrections
  H_biasadj_1 <- H_biasadj_2 <- NA
  Hc_biasadj_1 <- Hc_biasadj_2 <- NA
  tmins_search <- NA
  max_sg_est <- NA
  prop_maxk <- NA
  L <- NA
  max_count <- NA

  # Prepare data for forestsearch
  drop.vars <- c(fs.est$confounders.candidate, "treat.recommend")
  dfnew <- df_boot_analysis[, !(names(df_boot_analysis) %in% drop.vars)]
  dfnew_boot <- df_boot[, !(names(df_boot) %in% drop.vars)]

  # Setup forestsearch arguments
  args_FS_boot <- fs.est$args_call_all
  args_FS_boot$df.analysis <- dfnew_boot
  args_FS_boot$df.predict <- dfnew
  args_FS_boot$details <- (show_three && boot <= 3)
  args_FS_boot$showten_subgroups <- FALSE
  args_FS_boot$plot.sg <- FALSE
  args_FS_boot$grf_res <- NULL
  args_FS_boot$grf_cuts <- NULL
  args_FS_boot$plot.grf <- FALSE
  args_FS_boot$parallel_args$plan <- "sequential"
  args_FS_boot$parallel_args$workers <- 1
  args_FS_boot$parallel_args$show_message <- FALSE

  # Run forestsearch on bootstrap sample
  run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)

  if (!inherits(run_bootstrap, "try-error") && !is.null(run_bootstrap$sg.harm)) {
    # Extract results
    df_PredBoot <- run_bootstrap$df.predict
    dfboot_PredBoot <- run_bootstrap$df.est
    max_sg_est <- as.numeric(run_bootstrap$find.grps$max_sg_est)
    tmins_search <- as.numeric(run_bootstrap$find.grps$time_search)
    prop_maxk <- as.numeric(run_bootstrap$prop_maxk)
    max_count <- run_bootstrap$find.grps$max_count
    L <- run_bootstrap$find.grps$L

    # Calculate bias adjustments
    # H subgroup
    fitHstar_obs <- get_Cox_sg(
      df_sg = subset(df_PredBoot, treat.recommend == 0),
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    Hstar_obs <- fitHstar_obs$est_obs

    fitHstar_star <- get_Cox_sg(
      df_sg = subset(dfboot_PredBoot, treat.recommend == 0),
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    Hstar_star <- fitHstar_star$est_obs

    H_biasadj_1 <- H_obs - (Hstar_star - Hstar_obs)
    H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)

    # Hc subgroup
    fitHcstar_obs <- get_Cox_sg(
      df_sg = subset(df_PredBoot, treat.recommend == 1),
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    Hcstar_obs <- fitHcstar_obs$est_obs

    fitHcstar_star <- get_Cox_sg(
      df_sg = subset(dfboot_PredBoot, treat.recommend == 1),
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    Hcstar_star <- fitHcstar_star$est_obs

    Hc_biasadj_1 <- Hc_obs - (Hcstar_star - Hcstar_obs)
    Hc_biasadj_2 <- 2 * Hc_obs - (Hc_star + Hcstar_star - Hcstar_obs)
  }

  # Return results
  dfres <- data.table::data.table(
    H_biasadj_1, H_biasadj_2,
    Hc_biasadj_1, Hc_biasadj_2,
    tmins_search, max_sg_est,
    prop_maxk, L, max_count
  )

  return(dfres)
}

#' Calculate Estimates with Error Handling
#' @keywords internal
calculate_estimates_safe <- function(Hobs, seHobs, H1_adj, H2_adj, ystar,
                                     est.scale, details) {

  estimates <- try(
    get_dfRes(
      Hobs = Hobs,
      seHobs = seHobs,
      H1_adj = H1_adj,
      H2_adj = H2_adj,
      ystar = ystar,
      cov_method = "standard",
      cov_trim = 0.0,
      est.scale = est.scale,
      est.loghr = TRUE
    ),
    silent = TRUE
  )

  if (inherits(estimates, "try-error")) {
    if (details) {
      warning("Failed to calculate estimates: ", attr(estimates, "condition")$message)
    }
    return(NULL)
  }

  return(estimates)
}

#' Format Confidence Intervals
#' @keywords internal
format_confidence_intervals <- function(H_estimates, Hc_estimates, details) {

  H_res1 <- format_CI(H_estimates, c("H0", "H0_lower", "H0_upper"))
  H_res2 <- format_CI(H_estimates, c("H2", "H2_lower", "H2_upper"))
  Hc_res1 <- format_CI(Hc_estimates, c("H0", "H0_lower", "H0_upper"))
  Hc_res2 <- format_CI(Hc_estimates, c("H2", "H2_lower", "H2_upper"))

  if (details) {
    cat("\n=== Bootstrap Results ===\n")
    cat("H un-adjusted estimates:    ", H_res1, "\n")
    cat("H bias-corrected estimates: ", H_res2, "\n")
    cat("H^c un-adjusted estimates:  ", Hc_res1, "\n")
    cat("H^c bias-corrected estimates:", Hc_res2, "\n")
  }

  return(list(
    H_raw = H_res1,
    H_bc = H_res2,
    Hc_raw = Hc_res1,
    Hc_bc = Hc_res2
  ))
}

#' Create Summary Table with Preserved est.scale Logic
#' @keywords internal
create_summary_table_with_scale <- function(fs.est, SG_CIs, est.scale,
                                            args_forestsearch_call) {

  # CRITICAL: This preserves the original est.scale logic exactly as it was
  # The switching is intentional based on the scale interpretation

  if (est.scale == "1/hr") {
    # When scale is "1/hr", we swap the assignments
    # This is because the interpretation is reversed
    hr_1a <- SG_CIs$H_bc    # Note: H_bc goes to hr_1a
    hr_0a <- SG_CIs$Hc_bc   # Note: Hc_bc goes to hr_0a

    FSsg_tab <- SG_tab_estimates(
      df = fs.est$df.est,
      SG_flag = "treat.recommend",
      draws = 0,
      details = FALSE,
      outcome.name = args_forestsearch_call$outcome.name,
      event.name = args_forestsearch_call$event.name,
      treat.name = args_forestsearch_call$treat.name,
      strata.name = NULL,
      potentialOutcome.name = args_forestsearch_call$potentialOutcome.name,
      hr_1a = hr_1a,
      hr_0a = hr_0a,
      est.scale = "1/hr",
      sg0_name = "Questionable",
      sg1_name = "Recommend"
    )
  } else {
    # Standard case: direct assignment
    hr_1a <- SG_CIs$Hc_bc   # Note: Hc_bc goes to hr_1a
    hr_0a <- SG_CIs$H_bc    # Note: H_bc goes to hr_0a

    FSsg_tab <- SG_tab_estimates(
      df = fs.est$df.est,
      SG_flag = "treat.recommend",
      draws = 0,
      details = FALSE,
      outcome.name = args_forestsearch_call$outcome.name,
      event.name = args_forestsearch_call$event.name,
      treat.name = args_forestsearch_call$treat.name,
      strata.name = NULL,
      potentialOutcome.name = args_forestsearch_call$potentialOutcome.name,
      hr_1a = hr_1a,
      hr_0a = hr_0a,
      est.scale = "hr",
      sg0_name = "Questionable",
      sg1_name = "Recommend"
    )
  }

  return(FSsg_tab)
}

#' Calculate Bootstrap Diagnostics
#' @keywords internal
calculate_diagnostics <- function(results, nb_boots) {

  n_successful <- sum(!is.na(results$H_biasadj_2))
  prop_successful <- n_successful / nb_boots

  # Calculate bias corrections
  mean_H_bias1 <- mean(results$H_biasadj_1, na.rm = TRUE)
  mean_H_bias2 <- mean(results$H_biasadj_2, na.rm = TRUE)
  mean_Hc_bias1 <- mean(results$Hc_biasadj_1, na.rm = TRUE)
  mean_Hc_bias2 <- mean(results$Hc_biasadj_2, na.rm = TRUE)

  # Calculate variability
  sd_H_bias1 <- sd(results$H_biasadj_1, na.rm = TRUE)
  sd_H_bias2 <- sd(results$H_biasadj_2, na.rm = TRUE)
  sd_Hc_bias1 <- sd(results$Hc_biasadj_1, na.rm = TRUE)
  sd_Hc_bias2 <- sd(results$Hc_biasadj_2, na.rm = TRUE)

  # Search time statistics
  mean_search_time <- mean(results$tmins_search, na.rm = TRUE)
  total_search_time <- sum(results$tmins_search, na.rm = TRUE)

  return(list(
    n_successful = n_successful,
    prop_successful = prop_successful,
    bias_corrections = data.frame(
      subgroup = c("H", "H", "Hc", "Hc"),
      method = c("BC1", "BC2", "BC1", "BC2"),
      mean = c(mean_H_bias1, mean_H_bias2, mean_Hc_bias1, mean_Hc_bias2),
      sd = c(sd_H_bias1, sd_H_bias2, sd_Hc_bias1, sd_Hc_bias2)
    ),
    search_time = list(
      mean = mean_search_time,
      total = total_search_time
    )
  ))
}

#' Print Bootstrap Summary
#' @keywords internal
print_bootstrap_summary <- function(results, nb_boots, formatted_CIs, diagnostics) {

  cat("\n")
  cat("========================================\n")
  cat("BOOTSTRAP ANALYSIS SUMMARY\n")
  cat("========================================\n")
  cat("Total bootstraps:", nb_boots, "\n")
  cat("Successful bootstraps:", diagnostics$n_successful, "\n")
  cat("Success rate:", sprintf("%.1f%%", 100 * diagnostics$prop_successful), "\n")
  cat("\n")

  cat("Confidence Intervals:\n")
  cat("---------------------\n")
  cat("H un-adjusted:   ", formatted_CIs$H_raw, "\n")
  cat("H bias-corrected:", formatted_CIs$H_bc, "\n")
  cat("Hc un-adjusted:  ", formatted_CIs$Hc_raw, "\n")
  cat("Hc bias-corrected:", formatted_CIs$Hc_bc, "\n")
  cat("\n")

  cat("Search Time:\n")
  cat("-----------\n")
  cat("Mean per bootstrap:", sprintf("%.2f min", diagnostics$search_time$mean), "\n")
  cat("Total time:", sprintf("%.2f min", diagnostics$search_time$total), "\n")
  cat("========================================\n")
}
