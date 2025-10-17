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
