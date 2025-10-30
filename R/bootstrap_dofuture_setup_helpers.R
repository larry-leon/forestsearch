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
    "count_boot_it",
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
    "get_conf_force",
    "evaluate_subgroup_consistency"
  ),
  bootstrap_parallel = c(
    "bootstrap_results",
    "bootstrap_ystar",
    "ensure_packages",
    "fit_cox_models",
    "build_cox_formula",
    "cox.formula.boot",
    "setup_parallel_SGcons",
    "filter_call_args",
    "summarize_bootstrap_results"
  )
)

#' Flatten required functions list to character vector
#'
#' @return Character vector of all function names needed in parallel workers
#' @keywords internal
get_bootstrap_exports <- function() {
  unlist(BOOTSTRAP_REQUIRED_FUNCTIONS, use.names = FALSE)
}

ensure_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Required package(s) not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install with: install.packages(c('",
      paste(missing, collapse = "', '"), "'))",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


