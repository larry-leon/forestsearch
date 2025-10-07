
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
#' @importFrom foreach foreach %dopar%
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

#' Bootstrap Results for ForestSearch
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
#' @export

bootstrap_results <- function(fs.est, df_boot_analysis, cox.formula.boot, nb_boots, show_three, H_obs, Hc_obs, reset_parallel, boot_workers) {
  NN <- nrow(df_boot_analysis)
  id0 <- seq_len(NN)
  foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE,
                           add = c("calc_cov",  "calculate_counts", "analyze_subgroups", "calculate_potential_hr","ci.est","count.id","CV_sgs",
                           "cox_summary","df_counting","double_robust_scores", "extract_subgroup","format_results", "get_targetEst","getci_Cox",
                           "getCIs","grf.estimates.out","hrCI_format","km_summary","n_pcnt","plot_subgroup","plot_weighted_km",
                           "prepare_subgroup_data","quiet","rmst_calculation","sg_tables","sort_subgroups","SummaryStat","var_summary",
                           "get_FSdata", "dummy","run_bootstrap",
                             "forestsearch", "forestsearch_bootstrap_dofuture","get_combinations_info",
                             "get_dfpred",
                             "grf.subg.harm.survival",
                             "subgroup.search",
                             "subgroup.consistency",
                             "lasso_selection",
                             "get_Cox_sg",
                             "get_conf_force",
                             "filter_by_lassokeep",
                             "is.continuous",
                             "process_conf_force_expr",
                             "is_flag_continuous",
                             "is_flag_drop", "acm.disjctif",  "acm.util.df2", "acm.util.df", "dummy2","ztrail","one.zero",
                             "get_dfRes", "get_subgroup_membership",
                             "SG_tab_estimates",
                             "prepare_data",
                             "run_grf",
                             "evaluate_subgroups",
                             "summarize_results",
                              "clean_data", "qlow", "qhigh","FS_labels","thiscut","get_cut_name",
                             "bootstrap_results", "remove_redundant_subgroups", "sg_consistency_out","get_split_hr","cut_var",
                             "bootstrap_ystar", "ensure_packages", "fit_cox_models", "build_cox_formula", "cox.formula.boot",
                             "format_CI","setup_parallel_SGcons", "get_covs_in", "extract_idx_flagredundancy"
                           )),
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
    args_FS_boot$details <- show3
    args_FS_boot$showten_subgroups <- FALSE
    args_FS_boot$plot.sg <- FALSE
    # Re-do grf for bootstrap data
    # NOTE: setting grf_res and grf_cuts to NULL
    # Induces evaluation to revert back to default
    args_FS_boot$grf_res <- NULL
    args_FS_boot$grf_cuts <- NULL
    # In bootstrap re-set parallel_args per specification in this call
    # For parallel_args we do NOT want to revert back to default
    # because default is not a null list
    if(reset_parallel){
    args_FS_boot[["parallel_args"]] <- list()
    } else {
    args_FS_boot$parallel_args$workers <- boot_workers
    args_FS_boot$parallel_args$show_message <- FALSE
    }

    #print(args_FS_boot)
    #cat("Length of parallel args",c(length(args_FS_boot$parallel_args)),"\n")

    run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)

    if (inherits(run_bootstrap, "try-error")) warning("Bootstrap failure")

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
      rm(fitHstar_obs)
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
      rm(fitHcstar_star)
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


# Do not export
# For checking bootstrap to initiate defaults
forchecking <- function(fs){
fs.est <- fs
nb_boots <- 3
details <- TRUE
show_three <- TRUE
 reset_parallel_fs <- TRUE
 boot_workers <- 6
 parallel_args <- list(plan = "multisession", workers = 6, show_message = TRUE)
}


#' ForestSearch Bootstrap with doFuture Parallelization
#'
#' Orchestrates bootstrap analysis for ForestSearch using doFuture parallelization.
#'
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
#' @importFrom data.table data.table
#' @export

forestsearch_bootstrap_dofuture <- function(fs.est, nb_boots, details=FALSE, show_three=FALSE, reset_parallel_fs = TRUE, boot_workers = 3,
                                            parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE)) {

  args_forestsearch_call <- fs.est$args_call_all

  # 1. Ensure packages
  ensure_packages(c("data.table", "foreach", "doFuture", "doRNG", "survival"))

  # 2. Build formula
  args_build <- base::names(formals(build_cox_formula))
  # align with args_call_all
  args_build_filtered <- args_forestsearch_call[names(args_forestsearch_call) %in% args_build]

  #cox.formula.boot <- build_cox_formula(fs.est$outcome.name, fs.est$event.name, fs.est$treat.name)

  cox.formula.boot <- base::do.call(build_cox_formula, args_build_filtered)

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

  # 5. Bootstrap results

  # Remove
  # Remove after testing
  # Suppress printing message
  #if(length(parallel_args) > 0) parallel_args$show_message <- FALSE

  # Note: reset_parallel_fs re-sets parallel for subgroup consistency in forestsearch
  # That is reset_parallel_fs = TRUE only the outer *bootstrap* loop is parallelized

  # Remove
  #cat("Running bootstraps now","\n")

  results <- bootstrap_results(fs.est, fs.est$df.est, cox.formula.boot, nb_boots, show_three, H_obs, Hc_obs, reset_parallel_fs, boot_workers)

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

