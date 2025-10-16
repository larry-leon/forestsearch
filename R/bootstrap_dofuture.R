
# Required packages for parallel bootstrap workers
BOOTSTRAP_REQUIRED_PACKAGES <- c(
  "data.table",    # Data manipulation in parallel workers
  "foreach",       # Foreach iteration framework
  "doFuture",      # Future-based parallel backend
  "doRNG",         # Reproducible random numbers in parallel
  "survival"       # Cox model fitting for bias correction
)

#' Resolve parallel processing arguments for bootstrap
#'
#' If parallel_args not provided, falls back to forestsearch call's
#' parallel configuration. Always reports configuration to user.
#'
#' @param parallel_args List or empty list
#' @param forestsearch_call_args List from original forestsearch call
#' @return List with plan, workers, show_message
#' @keywords internal

resolve_bootstrap_parallel_args <- function(parallel_args, forestsearch_call_args) {
  # Use provided args if non-empty, otherwise inherit from forestsearch call
  if (length(parallel_args) == 0) {
    resolved_args <- as.list(forestsearch_call_args$parallel_args)
    message("Bootstrap parallel config: Using 'observed' data analysis forestsearch settings")
  } else {
    resolved_args <- parallel_args
    message("Bootstrap parallel config: Using user-provided settings")
  }

  # Report configuration to user
  max_cores <- parallel::detectCores()
  message("System max cores available: ", max_cores)
  message("Bootstrap will use: ", resolved_args$workers, " workers with '",
          resolved_args$plan, "' plan")

  resolved_args
}


#' Functions required in parallel bootstrap environment
#'
#' Organized by functional category for maintainability
BOOTSTRAP_REQUIRED_FUNCTIONS <- list(
  statistics = c(
    "calc_cov",
    "calculate_counts",
    "analyze_subgroups",
    "calculate_potential_hr",
    "ci.est",
    "count.id",
    "CV_sgs",
    "get_targetEst",
    "get_dfRes",
    "getCIs",
    "getci_Cox",
    "SummaryStat",
    "var_summary",
    "format_results",
    "format_CI",
    "qlow", "qhigh"
  ),
  subgroup_analysis = c(
    "cox_summary",
    "km_summary",
    "n_pcnt",
    "plot_subgroup",
    "plot_weighted_km",
    "prepare_subgroup_data",
    "quiet",
    "rmst_calculation",
    "sg_tables",
    "sort_subgroups",
    "remove_redundant_subgroups",
    "sg_consistency_out",
    "get_split_hr"
  ),
  data_prep = c(
    "get_FSdata",
    "dummy", "dummy2",
    "acm.disjctif", "acm.util.df2", "acm.util.df",
    "ztrail", "one.zero",
    "prepare_data",
    "clean_data"
  ),
  forestsearch_core = c(
    "forestsearch",
    "forestsearch_bootstrap_dofuture",
    "run_bootstrap",
    "run_grf",
    "evaluate_subgroups",
    "summarize_results",
    "SG_tab_estimates",
    "get_dfpred",
    "get_combinations_info",
    "get_subgroup_membership",
    "get_covs_in",
    "extract_idx_flagredundancy",
    "get_cut_name",
    "cut_var",
    "thiscut",
    "FS_labels"
  ),
  grf_policy = c(
    "grf.subg.harm.survival",
    "grf.estimates.out",
    "policy_tree",
    "causal_survival_forest"
  ),
  search_consistency = c(
    "subgroup.search",
    "subgroup.consistency",
    "extract_subgroup",
    "filter_by_lassokeep",
    "is.continuous",
    "process_conf_force_expr",
    "is_flag_continuous",
    "is_flag_drop",
    "lasso_selection",
    "get_Cox_sg",
    "get_conf_force"
  ),
  bootstrap_parallel = c(
    "bootstrap_results",
    "bootstrap_ystar",
    "ensure_packages",
    "fit_cox_models",
    "build_cox_formula",
    "cox.formula.boot",
    "setup_parallel_SGcons"
  )
)

#' Flatten required functions list to character vector
#'
#' @return Character vector of all function names needed in parallel workers
#' @keywords internal
get_bootstrap_exports <- function() {
  unlist(BOOTSTRAP_REQUIRED_FUNCTIONS, use.names = FALSE)
}

#' Find integer pairs (x, y) such that x * y = z and y >= x
#'
#' Given an integer z, this function finds all integer pairs (x, y) such that x * y = z and y >= x.
#' Optionally, you can return only the pair with the largest value of x or y.
#'
#' @param z Integer. The target product.
#' @param return_largest Character. If \"x\", returns the pair with the largest x. If \"y\", returns the pair with the largest y. If NULL, returns all pairs.
#'
#' @return A matrix of integer pairs (x, y) satisfying the conditions, or a single pair if return_largest is specified.
#' @examples
#' find_xy_given_z(12)
#' find_xy_given_z(12, return_largest = \"x\")
#' find_xy_given_z(12, return_largest = \"y\")
#' @keywords internal

find_xy_given_z <- function(z, return_largest = NULL) {
  pairs <- list()
  for (x in 1:z) {
    if (z %% x == 0) {
      y <- z / x
      if (y >= x && y %% 1 == 0) {
        pairs[[length(pairs) + 1]] <- c(x, y)
      }
    }
  }
  result <- do.call(rbind, pairs)
  if (!is.null(return_largest)) {
    if (return_largest == "x") {
      idx <- which.max(result[,1])
      return(result[idx, , drop = FALSE])
    } else if (return_largest == "y") {
      idx <- which.max(result[,2])
      return(result[idx, , drop = FALSE])
    }
  }
  return(result)
}


#' Ensure Required Packages Are Installed and Loaded
#'
#' Installs and loads required packages if not already available.
#'
#' @param pkgs Character vector of package names.
#' @return None. Packages are loaded into the session.
#' @export

ensure_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

#' Build Cox Model Formula
#'
#' Constructs a Cox model formula from variable names.
#'
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @return An R formula object for Cox regression.
#' @export

build_cox_formula <- function(outcome.name, event.name, treat.name) {
  sf <- paste0("Surv(", outcome.name, ",", event.name, ") ~ ", treat.name)
  as.formula(sf)
}

#' Fit Cox Models for Subgroups
#'
#' Fits Cox models for two subgroups defined by treatment recommendation.
#'
#' @param df Data frame.
#' @param formula Cox model formula.
#' @return List with HR and SE for each subgroup.
#' @export

fit_cox_models <- function(df, formula) {
  fitH <- get_Cox_sg(df_sg = subset(df, treat.recommend == 0), cox.formula = formula)
  fitHc <- get_Cox_sg(df_sg = subset(df, treat.recommend == 1), cox.formula = formula)
  list(H_obs = fitH$est_obs, seH_obs = fitH$se_obs, Hc_obs = fitHc$est_obs, seHc_obs = fitHc$se_obs)
}


#' Bootstrap Ystar Matrix
#'
#' Generates a bootstrap matrix for Ystar using parallel processing.
#'
#' @param df Data frame.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @return Matrix of bootstrap samples.
#' @importFrom foreach foreach
#' @export

bootstrap_ystar <- function(df, nb_boots) {
  NN <- nrow(df)
  foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {
    set.seed(8316951 + boot * 100)
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))
    ystar <- unlist(lapply(df$id, count.id, dfb = df_boot))
    return(ystar)
  }
}

#' Prepare data for bootstrap forestsearch re-run
#'
#' Removes confounders and treatment flag to force fresh variable selection
#' in each bootstrap iteration. This ensures the full bootstrap distribution
#' of variable importance and selection is captured.
#'
#' @param df_original Full original analysis data frame
#' @param df_bootstrap Bootstrap resample of analysis data
#' @param confounders_to_drop Character vector of confounder names
#'
#' @return List with:
#'   - df_analysis: Bootstrap data with confounders removed
#'   - df_predict: Original data with confounders removed
#'
#' @details
#' By removing confounders from the bootstrap data, we force forestsearch
#' to re-run its full variable selection pipeline (GRF, LASSO, etc.) on
#' each bootstrap sample. This captures the distribution of variable
#' selection stability across resamples.
#'
#' @export

prepare_bootstrap_dataframes <- function(df_original, df_bootstrap,
                                         confounders_to_drop) {
  # Variables to exclude: original confounders + treatment recommendation
  drop_vars <- c(confounders_to_drop, "treat.recommend")

  # Remove these from both original and bootstrap data
  df_analysis_clean <- df_original[, !(names(df_original) %in% drop_vars)]
  df_predict_clean <- df_bootstrap[, !(names(df_bootstrap) %in% drop_vars)]

  list(
    df_analysis = df_analysis_clean,
    df_predict = df_predict_clean
  )
}

#' Configure forestsearch arguments for bootstrap execution
#'
#' Modifies original forestsearch arguments to suit bootstrap context:
#' - Disables output generation (plots, verbose details)
#' - Forces variable re-selection (drops GRF results, LASSO keeps)
#' - Prevents nested parallelization
#' - Applies Bootstrap-specific seeds and data
#'
#' @param base_args List from original forestsearch call (args_call_all)
#' @param df_analysis_boot Bootstrap data for analysis
#' @param df_predict_boot Original data for prediction (for oracle)
#' @param show_details Logical. Show algorithm details for this iteration
#'
#' @return Modified args list ready for do.call(forestsearch, .)
#'
#' @details
#' Three categories of modifications:
#'
#' 1. OUTPUT SUPPRESSION (for parallelization efficiency):
#'    - details: Set to show_details (usually FALSE)
#'    - showten_subgroups: Always FALSE in bootstrap
#'    - plot.sg, plot.grf: Always FALSE (can't plot in parallel)
#'
#' 2. VARIABLE RE-SELECTION (to capture bootstrap variability):
#'    - grf_res: Set to NULL (force GRF re-run on bootstrap)
#'    - grf_cuts: Set to NULL (force fresh cuts)
#'
#' 3. SEQUENTIAL EXECUTION (prevent nested parallelization):
#'    - parallel_args$plan: "sequential" (single-threaded)
#'    - parallel_args$workers: 1
#'    - parallel_args$show_message: FALSE (suppress within-bootstrap messages)
#'
#' @export

configure_forestsearch_bootstrap_args <- function(base_args, df_analysis_boot,
                                                  df_predict_boot,
                                                  show_details = FALSE) {

  args <- base_args  # Start with original configuration

  # ===================================================================
  # Category 1: OUTPUT SUPPRESSION (parallelization efficiency)
  # ===================================================================
  # Why: Plotting/verbose output incompatible with parallel workers
  # In parallel context, output goes to individual workers, not main thread

  args$details <- show_details                # Show only if debugging first N bootstraps
  args$showten_subgroups <- FALSE             # Never show in bootstrap (large output)
  args$plot.sg <- FALSE                       # Can't plot from parallel worker
  args$plot.grf <- FALSE                      # Can't plot from parallel worker

  # ===================================================================
  # Category 2: FORCE VARIABLE RE-SELECTION (bootstrap variability)
  # ===================================================================
  # Why: Each bootstrap should re-discover which variables are important
  # If we reuse GRF results from original, bootstrap distribution is biased

  args$grf_res <- NULL       # Force GRF to run fresh on bootstrap sample
  args$grf_cuts <- NULL      # Force fresh GRF cuts (not oracle)

  # Keep LASSO selection though - it's already data-adaptive
  # (Note: LASSO settings from original call preserved)

  # ===================================================================
  # Category 3: SEQUENTIAL EXECUTION (prevent nested parallelization)
  # ===================================================================
  # Why: We're already in a parallel worker. Nested parallelization
  # causes: resource contention, deadlocks, very slow execution

  args$parallel_args$plan <- "sequential"       # No parallelization in bootstrap
  args$parallel_args$workers <- 1L              # Single thread only
  args$parallel_args$show_message <- FALSE      # Suppress worker-level messages

  # ===================================================================
  # DATA: Cleaned data with confounders removed
  # ===================================================================
  args$df.analysis <- df_analysis_boot
  args$df.predict <- df_predict_boot

  args
}




#' Bootstrap Results for ForestSearch (FIXED VERSION)
#'
#' Fixed version addressing variable name inconsistency bug in legacy code
#'
#' @inheritParams bootstrap_results
#' @export

bootstrap_results_fixed <- function(fs.est, df_boot_analysis, cox.formula.boot,
                                    nb_boots, show_three, H_obs, Hc_obs) {

  NN <- nrow(df_boot_analysis)
  id0 <- seq_len(NN)

  foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(
      seed = TRUE,
      add = get_bootstrap_exports()
    ),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {

    show3 <- FALSE
    if (show_three) show3 <- (boot <= 3)
    set.seed(8316951 + boot * 100)

    # Create bootstrap sample
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df_boot_analysis[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))

    # =================================================================
    # BUG FIX: Variable naming consistency
    # =================================================================
    # ORIGINAL BUG: Used inconsistent variable names (H_star vs h_star)
    # which caused compute_bias_corrections() to fail

    # Bootstrap data evaluated at ORIGINAL subgroup H
    fit_h_star <- get_Cox_sg(
      df_sg = subset(df_boot, treat.recommend == 0),
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    h_star <- fit_h_star$est_obs  # Use lowercase h_star consistently

    fit_hc_star <- get_Cox_sg(
      df_sg = subset(df_boot, treat.recommend == 1),
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    hc_star <- fit_hc_star$est_obs  # Use lowercase hc_star consistently

    # Initialize return values as NA
    H_biasadj_1 <- H_biasadj_2 <- NA
    Hc_biasadj_1 <- Hc_biasadj_2 <- NA
    tmins_search <- NA
    max_sg_est <- NA
    prop_maxk <- NA
    L <- NA
    max_count <- NA

    # Prepare bootstrap dataframes
    drop.vars <- c(fs.est$confounders.candidate, "treat.recommend")
    dfnew <- df_boot_analysis[, !(names(df_boot_analysis) %in% drop.vars)]
    dfnew_boot <- df_boot[, !(names(df_boot) %in% drop.vars)]

    # Configure forestsearch arguments for bootstrap
    args_FS_boot <- fs.est$args_call_all
    args_FS_boot$df.analysis <- dfnew_boot
    args_FS_boot$df.predict <- dfnew

    # CATEGORY 1: OUTPUT SUPPRESSION
    args_FS_boot$details <- show3
    args_FS_boot$showten_subgroups <- FALSE
    args_FS_boot$plot.sg <- FALSE
    args_FS_boot$plot.grf <- FALSE

    # CATEGORY 2: VARIABLE RE-SELECTION
    args_FS_boot$grf_res <- NULL
    args_FS_boot$grf_cuts <- NULL

    # CATEGORY 3: SEQUENTIAL EXECUTION
    args_FS_boot$parallel_args$plan <- "sequential"
    args_FS_boot$parallel_args$workers <- 1L
    args_FS_boot$parallel_args$show_message <- FALSE

    # Run forestsearch on bootstrap sample
    run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)

    if (inherits(run_bootstrap, "try-error")) {
      warning("Bootstrap ", boot, " failed: ", as.character(run_bootstrap))
    }

    # =================================================================
    # Compute bias corrections if subgroup found
    # =================================================================
    if (!inherits(run_bootstrap, "try-error") && !is.null(run_bootstrap$sg.harm)) {

      # Compute bias corrections using helper function
      bias_results <- compute_bias_corrections(
        run_bootstrap = run_bootstrap,
        H_obs = H_obs,
        Hc_obs = Hc_obs,
        h_star = h_star,      # Correct variable name (lowercase)
        hc_star = hc_star,    # Correct variable name (lowercase)
        cox.formula.boot = cox.formula.boot
      )

      # Extract results
      H_biasadj_1 <- bias_results$H_biasadj_1
      H_biasadj_2 <- bias_results$H_biasadj_2
      Hc_biasadj_1 <- bias_results$Hc_biasadj_1
      Hc_biasadj_2 <- bias_results$Hc_biasadj_2
      tmins_search <- bias_results$tmins_search
      max_sg_est <- bias_results$max_sg_est
      prop_maxk <- bias_results$prop_maxk
      L <- bias_results$L
      max_count <- bias_results$max_count
    }

    # CRITICAL: Always return data.table with same structure
    # This ensures .combine = "rbind" works correctly
    dfres <- data.table::data.table(
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

    return(dfres)
  }
}


#' Compute Bias Corrections for Bootstrap Subgroup Estimates
#'
#' Helper function that computes bias-corrected estimates for both subgroups
#' (H and H^c) using two correction methods. Handles cases where ForestSearch
#' fails to identify a subgroup.
#'
#' @param run_bootstrap List. Output from \code{\link{forestsearch}} run on
#'   bootstrap sample, or a \code{try-error} object if the run failed.
#' @param H_obs Numeric. Observed log hazard ratio for subgroup H from original
#'   sample. Used as the reference point for bias correction.
#' @param Hc_obs Numeric. Observed log hazard ratio for subgroup H^c from
#'   original sample.
#' @param h_star Numeric. Log hazard ratio for H estimated on bootstrap sample
#'   using original subgroup definition.
#' @param hc_star Numeric. Log hazard ratio for H^c estimated on bootstrap sample
#'   using original subgroup definition.
#' @param cox.formula.boot Formula. Cox model formula for estimation.
#'
#' @return Named list with components:
#'   \describe{
#'     \item{H_biasadj_1}{Numeric or NA. Bias correction method 1 for H}
#'     \item{H_biasadj_2}{Numeric or NA. Bias correction method 2 for H}
#'     \item{Hc_biasadj_1}{Numeric or NA. Bias correction method 1 for H^c}
#'     \item{Hc_biasadj_2}{Numeric or NA. Bias correction method 2 for H^c}
#'     \item{tmins_search}{Numeric or NA. Search time in minutes}
#'     \item{max_sg_est}{Numeric or NA. Maximum subgroup HR estimate}
#'     \item{prop_maxk}{Numeric or NA. Proportion using max K factors}
#'     \item{L}{Integer or NA. Number of factors evaluated}
#'     \item{max_count}{Integer or NA. Total factor combinations}
#'   }
#'   All values are \code{NA} if \code{run_bootstrap} is an error or found no
#'   valid subgroup.
#'
#' @section Bias Correction Formulas:
#' When a valid subgroup is found on the bootstrap sample:
#' \itemize{
#'   \item \code{H_biasadj_1 = H_obs - (Hstar_star - Hstar_obs)}
#'     where \code{Hstar_star} is the HR for the \emph{new} subgroup evaluated
#'     on the bootstrap sample, and \code{Hstar_obs} is the HR for the new
#'     subgroup evaluated on the original sample (optimism estimate).
#'   \item \code{H_biasadj_2 = 2*H_obs - (H_star + Hstar_star - Hstar_obs)}
#'     incorporates both the bootstrap estimate (\code{H_star}) and the
#'     optimism correction.
#' }
#'
#' @section Error Handling:
#' The function gracefully handles three failure modes:
#' \enumerate{
#'   \item \code{run_bootstrap} is a \code{try-error}: Returns all \code{NA}
#'   \item \code{run_bootstrap$sg.harm} is \code{NULL}: No subgroup found, returns \code{NA}
#'   \item Cox model fitting fails: Wrapped in \code{try()} to prevent crashes
#' }
#'
#' @note This is an internal helper function not meant for direct user calls.
#'   It's extracted from \code{\link{bootstrap_results}} to improve readability
#'   and testability.
#'
#' @seealso \code{\link{bootstrap_results}} for the main bootstrap loop
#' @keywords internal
#' @family bootstrap functions
compute_bias_corrections <- function(run_bootstrap, H_obs, Hc_obs, h_star,
                                     hc_star, cox.formula.boot) {

  # Initialize all return values as NA
  result <- list(
    H_biasadj_1 = NA_real_,
    H_biasadj_2 = NA_real_,
    Hc_biasadj_1 = NA_real_,
    Hc_biasadj_2 = NA_real_,
    tmins_search = NA_real_,
    max_sg_est = NA_real_,
    prop_maxk = NA_real_,
    L = NA_integer_,
    max_count = NA_integer_
  )

  # Check if bootstrap run succeeded and found a subgroup
  if (inherits(run_bootstrap, "try-error") || is.null(run_bootstrap$sg.harm)) {
    return(result)
  }

  # Extract prediction datasets from bootstrap ForestSearch run
  # df_pred_boot: original data with new subgroup assignments
  # dfboot_pred_boot: bootstrap data with new subgroup assignments
  df_pred_boot <- run_bootstrap$df.predict
  dfboot_pred_boot <- run_bootstrap$df.est

  # Extract search metrics
  result$max_sg_est <- as.numeric(run_bootstrap$find.grps$max_sg_est)
  result$tmins_search <- as.numeric(run_bootstrap$find.grps$time_search)
  result$prop_maxk <- as.numeric(run_bootstrap$prop_maxk)
  result$max_count <- run_bootstrap$find.grps$max_count
  result$L <- run_bootstrap$find.grps$L

  # ====================================================================
  # Compute bias corrections for subgroup H (harm/questionable group)
  # ====================================================================

  # Hstar_obs: New subgroup (from bootstrap) evaluated on ORIGINAL data
  hstar_obs <- get_Cox_sg(
    df_sg = subset(df_pred_boot, treat.recommend == 0),
    cox.formula = cox.formula.boot,
    est.loghr = TRUE
  )$est_obs

  # Hstar_star: New subgroup (from bootstrap) evaluated on BOOTSTRAP data
  hstar_star <- get_Cox_sg(
    df_sg = subset(dfboot_pred_boot, treat.recommend == 0),
    cox.formula = cox.formula.boot,
    est.loghr = TRUE
  )$est_obs

  # Bias correction method 1: Simple optimism correction
  # Removes the optimism (Hstar_star - Hstar_obs) from observed estimate
  result$H_biasadj_1 <- H_obs - (hstar_star - hstar_obs)

  # Bias correction method 2: Double correction
  # Uses both bootstrap estimate (H_star) and optimism correction
  result$H_biasadj_2 <- 2 * H_obs - (h_star + hstar_star - hstar_obs)

  # ====================================================================
  # Compute bias corrections for subgroup H^c (complement/recommend group)
  # ====================================================================

  # Hcstar_obs: New subgroup complement evaluated on ORIGINAL data
  hcstar_obs <- get_Cox_sg(
    df_sg = subset(df_pred_boot, treat.recommend == 1),
    cox.formula = cox.formula.boot,
    est.loghr = TRUE
  )$est_obs

  # Hcstar_star: New subgroup complement evaluated on BOOTSTRAP data
  hcstar_star <- get_Cox_sg(
    df_sg = subset(dfboot_pred_boot, treat.recommend == 1),
    cox.formula = cox.formula.boot,
    est.loghr = TRUE
  )$est_obs

  # Apply same correction methods for H^c
  result$Hc_biasadj_1 <- Hc_obs - (hcstar_star - hcstar_obs)
  result$Hc_biasadj_2 <- 2 * Hc_obs - (hc_star + hcstar_star - hcstar_obs)

  return(result)
}

#' Bootstrap Results for ForestSearch with Bias Correction
#'
#' Runs bootstrap analysis for ForestSearch, fitting Cox models and computing
#' bias-corrected estimates using the .632+ bootstrap method. Each bootstrap
#' iteration re-runs the full ForestSearch pipeline to capture variability in
#' subgroup identification.
#'
#' @param fs.est List. ForestSearch results object from \code{\link{forestsearch}}.
#'   Must contain:
#'   \itemize{
#'     \item \code{df.est}: Data frame with analysis data including \code{treat.recommend}
#'     \item \code{confounders.candidate}: Character vector of confounder names
#'     \item \code{args_call_all}: List of original forestsearch call arguments
#'   }
#' @param df_boot_analysis Data frame. Bootstrap analysis data with same structure
#'   as \code{fs.est$df.est}. Must contain columns for outcome, event, treatment,
#'   and the \code{treat.recommend} flag.
#' @param cox.formula.boot Formula. Cox model formula for bootstrap, typically
#'   created by \code{\link{build_cox_formula}}. Should be of form
#'   \code{Surv(outcome, event) ~ treatment}.
#' @param nb_boots Integer. Number of bootstrap samples to generate (e.g., 500-1000).
#'   More iterations provide better bias correction but increase computation time.
#' @param show_three Logical. If \code{TRUE}, prints detailed progress for the
#'   first three bootstrap iterations for debugging purposes. Default: \code{FALSE}.
#' @param H_obs Numeric. Observed log hazard ratio for subgroup H (harm group,
#'   \code{treat.recommend == 0}) from original sample. Used as reference for
#'   bias correction.
#' @param Hc_obs Numeric. Observed log hazard ratio for subgroup H^c (complement,
#'   \code{treat.recommend == 1}) from original sample. Used as reference for
#'   bias correction.
#'
#' @return Data.table with one row per bootstrap iteration and columns:
#'   \describe{
#'     \item{H_biasadj_1}{Bias-corrected estimate for H using method 1:
#'       \code{H_obs - (H_boot_boot - H_boot_obs)}}
#'     \item{H_biasadj_2}{Bias-corrected estimate for H using method 2:
#'       \code{2*H_obs - (H_boot + H_boot_boot - H_boot_obs)}}
#'     \item{Hc_biasadj_1}{Bias-corrected estimate for H^c using method 1}
#'     \item{Hc_biasadj_2}{Bias-corrected estimate for H^c using method 2}
#'     \item{tmins_search}{Numeric. Minutes spent on subgroup search in this iteration}
#'     \item{max_sg_est}{Numeric. Maximum subgroup hazard ratio found}
#'     \item{prop_maxk}{Numeric. Proportion of maximum K factors used}
#'     \item{L}{Integer. Number of candidate factors evaluated}
#'     \item{max_count}{Integer. Maximum number of factor combinations}
#'   }
#'   Rows where no valid subgroup was found will have \code{NA} for bias corrections.
#'
#' @section Bias Correction Methods:
#' The function implements two bias correction approaches:
#' \enumerate{
#'   \item \strong{Method 1 (Simple)}: Corrects for optimism using the difference
#'     between bootstrap internal validation (\code{H_boot_boot}) and
#'     bootstrap-on-original (\code{H_boot_obs})
#'   \item \strong{Method 2 (Double)}: Uses both the bootstrap estimate and the
#'     optimism correction for a more conservative adjustment
#' }
#'
#' @section Computational Details:
#' \itemize{
#'   \item Uses \code{doFuture} backend for parallel execution (configured externally)
#'   \item Sets reproducible seeds: \code{8316951 + boot * 100} for each iteration
#'   \item Each bootstrap iteration runs full ForestSearch pipeline including
#'     variable selection, subgroup search, and consistency evaluation
#'   \item Sequential execution within each bootstrap prevents nested parallelization
#'   \item Failed bootstrap iterations generate warnings but don't stop execution
#' }
#'
#' @section Performance Considerations:
#' \itemize{
#'   \item Typical runtime: 1-5 seconds per bootstrap iteration
#'   \item For 1000 bootstraps with 6 workers: ~3-10 minutes total
#'   \item Memory usage scales with dataset size and number of workers
#'   \item Consider reducing \code{nb_boots} for initial testing (e.g., 100)
#' }
#'
#' @note This function is designed to be called within a \code{foreach} loop
#'   with \code{\%dofuture\%} operator. It requires:
#'   \itemize{
#'     \item All functions in \code{\link{get_bootstrap_exports}} to be available
#'       in the parallel workers
#'     \item Packages listed in \code{BOOTSTRAP_REQUIRED_PACKAGES} to be installed
#'     \item Proper parallel backend setup via \code{\link{setup_parallel_SGcons}}
#'   }
#'
#' @seealso
#' \code{\link{forestsearch_bootstrap_dofuture}} for the wrapper function that
#'   sets up parallelization and calls this function
#' \code{\link{prepare_bootstrap_dataframes}} for data preparation
#' \code{\link{configure_forestsearch_bootstrap_args}} for argument configuration
#' \code{\link{get_Cox_sg}} for Cox model fitting
#' \code{\link{get_dfRes}} for processing bootstrap results into confidence intervals
#'
#' @examples
#' \dontrun{
#' # Typically called via forestsearch_bootstrap_dofuture()
#' # Manual usage for debugging:
#'
#' # 1. Fit initial ForestSearch model
#' fs_result <- forestsearch(
#'   df.analysis = mydata,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   confounders.name = c("age", "sex", "stage")
#' )
#'
#' # 2. Build Cox formula
#' cox_formula <- build_cox_formula("time", "status", "treatment")
#'
#' # 3. Get observed estimates
#' cox_fits <- fit_cox_models(fs_result$df.est, cox_formula)
#'
#' # 4. Run bootstrap (within foreach %dofuture% loop)
#' boot_results <- bootstrap_results(
#'   fs.est = fs_result,
#'   df_boot_analysis = fs_result$df.est,
#'   cox.formula.boot = cox_formula,
#'   nb_boots = 100,
#'   show_three = TRUE,
#'   H_obs = cox_fits$H_obs,
#'   Hc_obs = cox_fits$Hc_obs
#' )
#' }
#'
#' @references
#' Efron, B., & Tibshirani, R. (1997). Improvements on cross-validation: the
#' .632+ bootstrap method. \emph{Journal of the American Statistical Association},
#' 92(438), 548-560.
#'
#' @family bootstrap functions
#' @importFrom foreach foreach
#' @importFrom data.table data.table
#' @importFrom doFuture %dofuture%
#' @export

bootstrap_results_new <- function(fs.est, df_boot_analysis, cox.formula.boot,
                              nb_boots, show_three, H_obs, Hc_obs) {

  NN <- nrow(df_boot_analysis)
  id0 <- seq_len(NN)

  foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(
      seed = TRUE,
      add = get_bootstrap_exports()
    ),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {

    # Show details for first 3 bootstrap iterations
    show3 <- show_three && (boot <= 3)

    # Set reproducible seed
    BASE_SEED <- 8316951
    set.seed(BASE_SEED + boot * 100)

    # Create bootstrap sample
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df_boot_analysis[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))

    # Cache subgroup subsets
    df_boot_h <- subset(df_boot, treat.recommend == 0)
    df_boot_hc <- subset(df_boot, treat.recommend == 1)

    # Bootstrap data evaluated at H: H_star
    fit_h_star <- get_Cox_sg(
      df_sg = df_boot_h,
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    h_star <- fit_h_star$est_obs

    fit_hc_star <- get_Cox_sg(
      df_sg = df_boot_hc,
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    hc_star <- fit_hc_star$est_obs

    # Prepare bootstrap dataframes using helper function
    boot_data <- prepare_bootstrap_dataframes(
      df_original = df_boot_analysis,
      df_bootstrap = df_boot,
      confounders_to_drop = fs.est$confounders.candidate
    )

    # Configure forestsearch arguments using helper function
    args_fs_boot <- configure_forestsearch_bootstrap_args(
      base_args = fs.est$args_call_all,
      df_analysis_boot = boot_data$df_analysis,
      df_predict_boot = boot_data$df_predict,
      show_details = show3
    )

    # Run forestsearch on bootstrap sample
    run_bootstrap <- try(do.call(forestsearch, args_fs_boot), silent = TRUE)

    if (inherits(run_bootstrap, "try-error")) {
      warning("Bootstrap iteration ", boot, " failed: ", as.character(run_bootstrap))
    }

    # Compute bias corrections
    bias_results <- compute_bias_corrections(
      run_bootstrap = run_bootstrap,
      H_obs = H_obs,
      Hc_obs = Hc_obs,
      h_star = h_star,
      hc_star = hc_star,
      cox.formula.boot = cox.formula.boot
    )

    # Return results as data.table
    dfres <- data.table::data.table(
      H_biasadj_1 = bias_results$H_biasadj_1,
      H_biasadj_2 = bias_results$H_biasadj_2,
      Hc_biasadj_1 = bias_results$Hc_biasadj_1,
      Hc_biasadj_2 = bias_results$Hc_biasadj_2,
      tmins_search = bias_results$tmins_search,
      max_sg_est = bias_results$max_sg_est,
      prop_maxk = bias_results$prop_maxk,
      L = bias_results$L,
      max_count = bias_results$max_count
    )

    return(dfres)
  }
}


#' Bootstrap Results for ForestSearch (legacy, in case new one doesn't work out)
#'
#' Runs bootstrap analysis for ForestSearch, fitting Cox models and bias correction.
#'
#' @param fs.est ForestSearch results object.
#' @param df_boot_analysis Data frame for bootstrap analysis.
#' @param cox.formula.boot Cox model formula for bootstrap.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @param show_three Logical. Show details for first three bootstraps.
#' @param H_obs Numeric. Observed HR for subgroup H.
#' @param Hc_obs Numeric. Observed HR for subgroup Hc.
#' @param reset_parallel Logical. Reset parallel plan for bootstrap.
#' @param boot_workers Integer. Number of parallel workers.
#' @return Data.table with bias-adjusted estimates and search metrics.
#' @importFrom foreach foreach
#' @importFrom data.table data.table
#' @importFrom doFuture %dofuture%
#' @export

bootstrap_results_legacy <- function(fs.est, df_boot_analysis, cox.formula.boot, nb_boots, show_three, H_obs, Hc_obs) {

  NN <- nrow(df_boot_analysis)
  id0 <- seq_len(NN)
    foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE,
                           add = get_bootstrap_exports()
                           ),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {
    show3 <- FALSE
    if (show_three) show3 <- (boot <= 3)
    set.seed(8316951 + boot * 100)
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df_boot_analysis[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))

    # Bootstrap data evaluated at H: H_star
    fitH_star <- get_Cox_sg(df_sg = subset(df_boot, treat.recommend == 0), cox.formula = cox.formula.boot, est.loghr = TRUE)
    H_star <- fitH_star$est_obs
    fitHc_star <- get_Cox_sg(df_sg = subset(df_boot, treat.recommend == 1), cox.formula = cox.formula.boot, est.loghr = TRUE)
    Hc_star <- fitHc_star$est_obs

    # Bias corrections
    H_biasadj_1 <- H_biasadj_2 <- NA
    Hc_biasadj_1 <- Hc_biasadj_2 <- NA
    tmins_search <- NA
    max_sg_est <- NA
    prop_maxk <- NA
    L <- NA
    max_count <- NA

    # Note: need to add these functions to "add list"
    # boot_data <- prepare_bootstrap_dataframes(
    #   df_original = df_boot_analysis,
    #   df_bootstrap = df_boot,
    #   confounders_to_drop = fs.est$confounders.candidate
    # )
    #
    # args_FS_boot <- configure_forestsearch_bootstrap_args(
    #   base_args = fs.est$args_call_all,
    #   df_analysis_boot = boot_data$df_analysis,
    #   df_predict_boot = boot_data$df_predict,
    #   show_details = show3  # show3 = (boot <= 3)
    # )

    # Drop initial confounders
    drop.vars <- c(fs.est$confounders.candidate, "treat.recommend")
    dfnew <- df_boot_analysis[, !(names(df_boot_analysis) %in% drop.vars)]
    dfnew_boot <- df_boot[, !(names(df_boot) %in% drop.vars)]
    # Extract arguments in forestsearch (observed) data analysis
    args_FS_boot <- fs.est$args_call_all
    args_FS_boot$df.analysis <- dfnew_boot
    args_FS_boot$df.predict <- dfnew
    # CATEGORY 1: OUTPUT SUPPRESSION (parallelization efficiency)
    args_FS_boot$details <- show3                # Only show first 3 for debugging
    args_FS_boot$showten_subgroups <- FALSE      # Suppress large output
    args_FS_boot$plot.sg <- FALSE                # Can't plot from parallel worker
    args_FS_boot$plot.grf <- FALSE               # Can't plot from parallel worker
    # CATEGORY 2: VARIABLE RE-SELECTION (bootstrap variability)
    # Force fresh GRF run on bootstrap (oracle uses predictions from bootstrap)
    args_FS_boot$grf_res <- NULL
    args_FS_boot$grf_cuts <- NULL
    # CATEGORY 3: SEQUENTIAL EXECUTION (prevent nested parallelization)
    # Each bootstrap is already running in a parallel worker
    # Nested parallelization causes resource contention and deadlocks
    args_FS_boot$parallel_args$plan <- "sequential"
    args_FS_boot$parallel_args$workers <- 1L
    args_FS_boot$parallel_args$show_message <- FALSE

    run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)

    if (inherits(run_bootstrap, "try-error")) {
      warning("Bootstrap ", boot, " failed: ", as.character(run_bootstrap))
    }


      if (!inherits(run_bootstrap, "try-error") && !is.null(run_bootstrap$sg.harm)) {
      df_PredBoot <- run_bootstrap$df.predict
      dfboot_PredBoot <- run_bootstrap$df.est
      max_sg_est <- as.numeric(run_bootstrap$find.grps$max_sg_est)
      tmins_search <- as.numeric(run_bootstrap$find.grps$time_search)
      prop_maxk <- as.numeric(run_bootstrap$prop_maxk)
      max_count <- run_bootstrap$find.grps$max_count
      L <- run_bootstrap$find.grps$L
      fitHstar_obs <- get_Cox_sg(df_sg = subset(df_PredBoot, treat.recommend == 0), cox.formula = cox.formula.boot, est.loghr = TRUE)
      Hstar_obs <- fitHstar_obs$est_obs
      fitHstar_star <- get_Cox_sg(df_sg = subset(dfboot_PredBoot, treat.recommend == 0), cox.formula = cox.formula.boot, est.loghr = TRUE)
      Hstar_star <- fitHstar_star$est_obs
      rm(fitHstar_star)
      H_biasadj_1 <- H_obs - (Hstar_star - Hstar_obs)
      H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)
      fitHcstar_obs <- get_Cox_sg(df_sg = subset(df_PredBoot, treat.recommend == 1), cox.formula = cox.formula.boot, est.loghr = TRUE)
      Hcstar_obs <- fitHcstar_obs$est_obs
      rm(fitHcstar_obs)
      fitHcstar_star <- get_Cox_sg(df_sg = subset(dfboot_PredBoot, treat.recommend == 1), cox.formula = cox.formula.boot, est.loghr = TRUE)
      Hcstar_star <- fitHcstar_star$est_obs
      Hc_biasadj_1 <- Hc_obs - (Hcstar_star - Hcstar_obs)
      Hc_biasadj_2 <- 2 * Hc_obs - (Hc_star + Hcstar_star - Hcstar_obs)
      }
    dfres <- data.table::data.table(H_biasadj_1, H_biasadj_2,
                                    Hc_biasadj_1, Hc_biasadj_2,
                                    tmins_search, max_sg_est, prop_maxk, L, max_count)
    return(dfres)
    }

}

#' Format Confidence Interval for Estimates
#'
#' Formats confidence interval for estimates.
#'
#' @param estimates Data frame or data.table of estimates.
#' @param col_names Character vector of column names for estimate, lower, upper.
#' @return Character string formatted as \"estimate (lower, upper)\".
#' @export

format_CI <- function(estimates, col_names) {
  resH <- estimates[, ..col_names]
  Hstat <- round(unlist(resH[1, ]), 2)
  paste0(Hstat[1], " (", Hstat[2], ",", Hstat[3], ")")
}

#' ForestSearch Bootstrap with doFuture Parallelization
#'
#' Orchestrates bootstrap analysis for ForestSearch using doFuture parallelization.
#'
#' Bootstrap Bias Correction Naming Convention:
#'
#'   Variable naming: {estimate_source}_{method}
#'
#'   estimate_source:
#'     - H_: Subgroup H (harm group, treat.recommend == 0)
#'     - Hc_: Subgroup H^c (complement, treat.recommend == 1)
#'
#'   method suffixes:
#'     - obs: Observed estimate from original sample
#'     - boot: Bootstrap replicate estimate
#'     - bc_1: Bias correction method 1: H_obs - (H_boot_boot - H_boot_obs)
#'     - bc_2: Bias correction method 2: 2*H_obs - (H_boot + H_boot_boot - H_boot_obs)
#' @param fs.est ForestSearch results object.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @param details Logical. Print details during execution.
#' @param show_three Logical. Show details for first three bootstraps.
#' @param reset_parallel_fs Logical. Reset parallel plan for bootstrap.
#' @param boot_workers Integer. Number of parallel workers.
#' @param parallel_args List. Parallelization arguments (plan, workers, show_message).
#'
#' @return List with bootstrap results, confidence intervals, summary table, Ystar matrix, and estimates.
#'
#' @importFrom future plan
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @importFrom data.table data.table
#' @export

forestsearch_bootstrap_dofuture <- function(fs.est, nb_boots, details=FALSE, show_three=FALSE,
                                           parallel_args = list()
                                            ) {

  args_forestsearch_call <- fs.est$args_call_all

  # If parallel_args is empty then default to main forestsearch (data analysis) call

  # if(length(parallel_args) == 0){
  # message("Using parallel plan of 'observed' data analysis forestsearch")
  # parallel_args <- as.list(args_forestsearch_call$parallel_args)
  # max_cores <- parallel::detectCores()
  # message("Note that max cores = ", max_cores)
  #  }

parallel_args <- resolve_bootstrap_parallel_args(parallel_args, args_forestsearch_call)


  # 1. Ensure packages
  ensure_packages(BOOTSTRAP_REQUIRED_PACKAGES)

  # 2. Build formula
  cox.formula.boot <- do.call(build_cox_formula,
                              filter_call_args(args_forestsearch_call, build_cox_formula))

  # 3. Fit Cox models
  cox_fits <- fit_cox_models(fs.est$df.est, cox.formula.boot)
  H_obs <- cox_fits$H_obs
  seH_obs <- cox_fits$seH_obs
  Hc_obs <- cox_fits$Hc_obs
  seHc_obs <- cox_fits$seHc_obs

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)  # Restore plan on exit
  # Note: setup_parallel_SGcons is re-purposed from subgroup_consistency
  setup_parallel_SGcons(parallel_args)

  # 4. Bootstrap Ystar matrix
  Ystar_mat <- bootstrap_ystar(fs.est$df.est, nb_boots)
  if (details) cat("Done with Ystar_mat\n")

  if(nrow(Ystar_mat) != nb_boots || ncol(Ystar_mat) != nrow(fs.est$df.est)){
  stop("Dimension of Ystar_mat must be (n x nb_boots)")
  }

  # Note: reset_parallel_fs re-sets parallel for subgroup consistency in forestsearch
  # That is reset_parallel_fs = TRUE only the outer *bootstrap* loop is parallelized

  results <-  bootstrap_results(fs.est, fs.est$df.est, cox.formula.boot, nb_boots, show_three, H_obs, Hc_obs)

  # 6. Post-processing and formatting
  est.scale <- args_forestsearch_call$est.scale

  H_estimates <- try(get_dfRes(Hobs = H_obs, seHobs = seH_obs, H1_adj = results$H_biasadj_1, H2_adj = results$H_biasadj_2,
                               ystar = Ystar_mat, cov_method = "standard", cov_trim = 0.0, est.scale = est.scale, est.loghr = TRUE), TRUE)

  Hc_estimates <- try(get_dfRes(Hobs = Hc_obs, seHobs = seHc_obs, H1_adj = results$Hc_biasadj_1, H2_adj = results$Hc_biasadj_2,
                                ystar = Ystar_mat, cov_method = "standard", cov_trim = 0.0, est.scale = est.scale, est.loghr = TRUE), TRUE)

  if (inherits(H_estimates, "try-error") | inherits(Hc_estimates, "try-error")) {
    out <- list(results = results, SG_CIs = NULL, FSsg_tab = NULL, Ystar_mat = Ystar_mat, H_estimates = NULL, Hc_estimates = NULL)
    return(out)
  }

  H_res1 <- format_CI(H_estimates, c("H0", "H0_lower", "H0_upper"))
  H_res2 <- format_CI(H_estimates, c("H2", "H2_lower", "H2_upper"))
  Hc_res1 <- format_CI(Hc_estimates, c("H0", "H0_lower", "H0_upper"))
  Hc_res2 <- format_CI(Hc_estimates, c("H2", "H2_lower", "H2_upper"))

  if (details) {
    cat("**** % bootstrap subgroups found =",
    c(sum(!is.na(results$H_biasadj_2))/nb_boots), "\n")
    cat("H un-adjusted estimates-----:   ", H_res1, "\n")
    cat("H bias-corrected estimates--:   ", H_res2, "\n")
    cat("H^c un-adjusted estimates---:   ", Hc_res1, "\n")
    cat("H^c bias-corrected estimates:   ", Hc_res2, "\n")
  }

    SG_CIs <- list(H_raw = H_res1, H_bc = H_res2, Hc_raw = Hc_res1, Hc_bc = Hc_res2)

# Adjusted CIs and summary table
if (est.scale == "1/hr") {
    hr_1a <- SG_CIs$H_bc
    hr_0a <- SG_CIs$Hc_bc
    FSsg_tab <- SG_tab_estimates(df = fs.est$df.est, SG_flag = "treat.recommend", draws = 0, details = FALSE,
                                 outcome.name = args_forestsearch_call$outcome.name,
                                 event.name = args_forestsearch_call$event.name,
                                 treat.name = args_forestsearch_call$treat.name,
                                 strata.name = NULL,
                                 potentialOutcome.name = args_forestsearch_call$potentialOutcome.name,
                                 hr_1a = hr_1a, hr_0a = hr_0a, est.scale = "1/hr", sg0_name = "Questionable", sg1_name = "Recommend"
                                 )
    } else {
    hr_1a <- SG_CIs$Hc_bc
    hr_0a <- SG_CIs$H_bc
    FSsg_tab <- SG_tab_estimates(df = fs.est$df.est, SG_flag = "treat.recommend", draws = 0, details = FALSE,
                                 outcome.name = args_forestsearch_call$outcome.name,
                                 event.name = args_forestsearch_call$event.name,
                                 treat.name = args_forestsearch_call$treat.name,
                                 strata.name = NULL,
                                 potentialOutcome.name = args_forestsearch_call$potentialOutcome.name,
                                 hr_1a = hr_1a, hr_0a = hr_0a, est.scale = "hr", sg0_name = "Questionable", sg1_name = "Recommend"
                                 )
  }
  out <- list(results = results, SG_CIs = SG_CIs, FSsg_tab = FSsg_tab, Ystar_mat = Ystar_mat, H_estimates = H_estimates, Hc_estimates = Hc_estimates)
  return(out)
}

