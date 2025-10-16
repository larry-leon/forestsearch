
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
    # Do NOT modify the above seed this needs to align with resampling
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))
    ystar <- unlist(lapply(df$id, count.id, dfb = df_boot))
    return(ystar)
  }
}

#' Bootstrap Results for ForestSearch with Bias Correction
#'
#' Runs bootstrap analysis for ForestSearch, fitting Cox models and computing
#' bias-corrected estimates and valid CIs (see vignette for references)
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
#' @param H_obs Numeric. Observed log hazard ratio for subgroup H (harm/questionable group,
#'   \code{treat.recommend == 0}) from original sample. Used as reference for
#'   bias correction.
#' @param Hc_obs Numeric. Observed log hazard ratio for subgroup H^c (complement/recommend,
#'   \code{treat.recommend == 1}) from original sample. Used as reference for
#'   bias correction.
#'
#' @return Data.table with one row per bootstrap iteration and columns:
#'   \describe{
#'     \item{H_biasadj_1}{Bias-corrected estimate for H using method 1:
#'       \code{H_obs - (Hstar_star - Hstar_obs)}}
#'     \item{H_biasadj_2}{Bias-corrected estimate for H using method 2:
#'       \code{2*H_obs - (H_star + Hstar_star - Hstar_obs)}}
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
#'   \item \strong{Method 1 (Simple Optimism)}: Corrects for optimism using the difference
#'     between bootstrap internal validation (\code{Hstar_star}) and
#'     bootstrap-on-original (\code{Hstar_obs}):
#'     \deqn{H_{adj1} = H_{obs} - (H*_{star} - H*_{obs})}
#'   \item \strong{Method 2 (Double Bootstrap)}: Uses both the bootstrap estimate and the
#'     optimism correction for a more conservative adjustment:
#'     \deqn{H_{adj2} = 2 \times H_{obs} - (H_{star} + H*_{star} - H*_{obs})}
#' }
#' where:
#' \itemize{
#'   \item \code{H_obs}: Original subgroup HR on original data
#'   \item \code{H_star}: Original subgroup HR on bootstrap data
#'   \item \code{Hstar_obs}: New subgroup (found in bootstrap) HR on original data
#'   \item \code{Hstar_star}: New subgroup (found in bootstrap) HR on bootstrap data
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
#'   \item Confounders are removed from bootstrap data to force fresh variable selection
#' }
#'
#' @section Bootstrap Configuration:
#' Each bootstrap iteration modifies ForestSearch arguments to:
#' \itemize{
#'   \item \strong{Suppress output}: \code{details}, \code{showten_subgroups},
#'     \code{plot.sg}, \code{plot.grf} all set to \code{FALSE}
#'   \item \strong{Force re-selection}: \code{grf_res} and \code{grf_cuts} set to \code{NULL}
#'   \item \strong{Prevent nested parallel}: \code{parallel_args$plan = "sequential"},
#'     \code{workers = 1}
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
#' @section Error Handling:
#' The function gracefully handles three failure modes:
#' \enumerate{
#'   \item Bootstrap sample creation fails: Returns row with all \code{NA}
#'   \item ForestSearch fails to run: Warns and returns row with all \code{NA}
#'   \item ForestSearch runs but finds no subgroup: Returns row with all \code{NA}
#' }
#' All three cases ensure the foreach loop can still combine results via \code{rbind}.
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
#' \code{\link{build_cox_formula}} for creating the Cox formula
#' \code{\link{fit_cox_models}} for initial Cox model fitting
#' \code{\link{get_Cox_sg}} for Cox model fitting on subgroups
#' \code{\link{get_dfRes}} for processing bootstrap results into confidence intervals
#' \code{\link{bootstrap_ystar}} for generating the Ystar matrix
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
#' # 4. Set up parallel backend
#' library(doFuture)
#' registerDoFuture()
#' plan(multisession, workers = 6)
#'
#' # 5. Run bootstrap (note: this is already parallelized internally)
#' boot_results <- bootstrap_results(
#'   fs.est = fs_result,
#'   df_boot_analysis = fs_result$df.est,
#'   cox.formula.boot = cox_formula,
#'   nb_boots = 100,
#'   show_three = TRUE,
#'   H_obs = cox_fits$H_obs,
#'   Hc_obs = cox_fits$Hc_obs
#' )
#'
#' # 6. Check results
#' summary(boot_results)
#'
#' # Proportion of bootstraps that found a subgroup
#' mean(!is.na(boot_results$H_biasadj_2))
#' }
#'
#' @family bootstrap functions
#' @importFrom foreach foreach
#' @importFrom data.table data.table
#' @importFrom progressr progressor handlers
#' @importFrom doFuture %dofuture%
#' @export

bootstrap_results <- function(fs.est, df_boot_analysis, cox.formula.boot,
                              nb_boots, show_three, H_obs, Hc_obs, show_progress = TRUE) {
  NN <- nrow(df_boot_analysis)
  id0 <- seq_len(NN)

  # Set up progress bar if requested and package available
  if (show_progress) {
    if (!requireNamespace("progressr", quietly = TRUE)) {
      warning("Package 'progressr' needed for progress bars. Install with: install.packages('progressr')")
      show_progress <- FALSE
    } else {
      progressr::handlers(global = TRUE)
      progressr::handlers("progress")
      p <- progressr::progressor(steps = nb_boots)
    }
  }

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
    # Do NOT modify the above seed this needs to align with ystar calculation
    # Create bootstrap sample
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df_boot_analysis[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))

    # =================================================================
    # Bootstrap data evaluated at ORIGINAL subgroup H
    # =================================================================
    # IMPORTANT: Use uppercase variable names to match legacy version
    fitH_star <- get_Cox_sg(
      df_sg = subset(df_boot, treat.recommend == 0),
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    H_star <- fitH_star$est_obs  # Use uppercase H_star

    fitHc_star <- get_Cox_sg(
      df_sg = subset(df_boot, treat.recommend == 1),
      cox.formula = cox.formula.boot,
      est.loghr = TRUE
    )
    Hc_star <- fitHc_star$est_obs  # Use uppercase Hc_star

    # =================================================================
    # Initialize bias corrections as NA (will remain NA if bootstrap fails)
    # =================================================================
    H_biasadj_1 <- H_biasadj_2 <- NA
    Hc_biasadj_1 <- Hc_biasadj_2 <- NA
    tmins_search <- NA
    max_sg_est <- NA
    prop_maxk <- NA
    L <- NA
    max_count <- NA

    # =================================================================
    # Prepare bootstrap dataframes - drop confounders and treat.recommend
    # =================================================================
    drop.vars <- c(fs.est$confounders.candidate, "treat.recommend")
    dfnew <- df_boot_analysis[, !(names(df_boot_analysis) %in% drop.vars)]
    dfnew_boot <- df_boot[, !(names(df_boot) %in% drop.vars)]

    # =================================================================
    # Configure forestsearch arguments for bootstrap
    # =================================================================
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

    # =================================================================
    # Run forestsearch on bootstrap sample
    # =================================================================
    run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)

    if (inherits(run_bootstrap, "try-error")) {
      warning("Bootstrap ", boot, " failed: ", as.character(run_bootstrap))
    }

    # =================================================================
    # CRITICAL FIX: Only compute bias corrections if bootstrap succeeded
    # AND found a valid subgroup
    # =================================================================
    if (!inherits(run_bootstrap, "try-error") && !is.null(run_bootstrap$sg.harm)) {

      # Extract prediction datasets from bootstrap ForestSearch run
      df_PredBoot <- run_bootstrap$df.predict
      dfboot_PredBoot <- run_bootstrap$df.est

      # Extract search metrics
      max_sg_est <- as.numeric(run_bootstrap$find.grps$max_sg_est)
      tmins_search <- as.numeric(run_bootstrap$find.grps$time_search)
      prop_maxk <- as.numeric(run_bootstrap$prop_maxk)
      max_count <- run_bootstrap$find.grps$max_count
      L <- run_bootstrap$find.grps$L

      # ==============================================================
      # Compute bias corrections for subgroup H (harm/questionable)
      # ==============================================================

      # Hstar_obs: New subgroup (from bootstrap) evaluated on ORIGINAL data
      fitHstar_obs <- get_Cox_sg(
        df_sg = subset(df_PredBoot, treat.recommend == 0),
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      Hstar_obs <- fitHstar_obs$est_obs

      # Hstar_star: New subgroup (from bootstrap) evaluated on BOOTSTRAP data
      fitHstar_star <- get_Cox_sg(
        df_sg = subset(dfboot_PredBoot, treat.recommend == 0),
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      Hstar_star <- fitHstar_star$est_obs

      # Bias correction method 1: Simple optimism correction
      H_biasadj_1 <- H_obs - (Hstar_star - Hstar_obs)

      # Bias correction method 2: Double correction
      H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)

      # ==============================================================
      # Compute bias corrections for subgroup H^c (complement/recommend)
      # ==============================================================

      # Hcstar_obs: New subgroup complement evaluated on ORIGINAL data
      fitHcstar_obs <- get_Cox_sg(
        df_sg = subset(df_PredBoot, treat.recommend == 1),
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      Hcstar_obs <- fitHcstar_obs$est_obs

      # Hcstar_star: New subgroup complement evaluated on BOOTSTRAP data
      fitHcstar_star <- get_Cox_sg(
        df_sg = subset(dfboot_PredBoot, treat.recommend == 1),
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
      Hcstar_star <- fitHcstar_star$est_obs

      # Apply same correction methods for H^c
      Hc_biasadj_1 <- Hc_obs - (Hcstar_star - Hcstar_obs)
      Hc_biasadj_2 <- 2 * Hc_obs - (Hc_star + Hcstar_star - Hcstar_obs)
    }

    # =================================================================
    # CRITICAL: Always return data.table with same structure
    # This ensures .combine = "rbind" works correctly
    # =================================================================
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

    # Update progress bar (only if available)
    if (show_progress) {
      p(message = sprintf("Bootstrap %d/%d complete", boot, nb_boots))
    }

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

