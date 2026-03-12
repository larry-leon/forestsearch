# ==============================================================================
# ForestSearch Cross-Validation Functions
# ==============================================================================
#
# This file contains K-fold cross-validation functions for ForestSearch:
#   - forestsearch_Kfold: Single K-fold cross-validation
#   - forestsearch_tenfold: Repeated K-fold simulations
#   - forestsearch_KfoldOut: Summarize CV results
#   - CV_sgs: Cross-validation subgroup match summary
#   - resolve_cv_parallel_args: Helper for parallel configuration
#
# ==============================================================================


#' ForestSearch K-Fold Cross-Validation
#'
#' Performs K-fold cross-validation for ForestSearch, evaluating subgroup
#' identification and agreement between training and test sets.
#'
#' @description
#' This function assesses the stability and reproducibility of ForestSearch
#' subgroup identification through cross-validation. For each fold:
#' \enumerate{
#'   \item Train ForestSearch on (K-1) folds
#'   \item Apply the identified subgroup to the held-out fold
#'   \item Compare predictions to the original full-data analysis
#' }
#'
#' @section Cross-Validation Types:
#' \itemize{
#'   \item \strong{Leave-One-Out (LOO)}: When \code{Kfolds = nrow(df)}, each
#'     observation is held out once. Most thorough but computationally intensive.
#'   \item \strong{K-Fold}: When \code{Kfolds < nrow(df)}, data is split into K
#'     roughly equal folds. Good balance of bias-variance tradeoff.
#' }
#'
#' @section Output Metrics:
#' The returned \code{resCV} data frame contains:
#' \itemize{
#'   \item \code{treat.recommend}: Prediction from CV model
#'   \item \code{treat.recommend.original}: Prediction from full-data model
#'   \item \code{cvindex}: Fold assignment
#'   \item \code{sg1}, \code{sg2}: Subgroup definitions found in each fold
#' }
#'
#' @param fs.est List. ForestSearch results object from \code{\link{forestsearch}}.
#'   Must contain \code{df.est} (data frame) and \code{args_call_all} (list of arguments).
#' @param Kfolds Integer. Number of folds (default: \code{nrow(fs.est$df.est)} for LOO).
#' @param seedit Integer. Random seed for fold assignment (default: 8316951).
#' @param parallel_args List. Parallelization configuration with elements:
#'   \itemize{
#'     \item \code{plan}: Character. One of "multisession", "multicore", "sequential"
#'     \item \code{workers}: Integer. Number of parallel workers
#'     \item \code{show_message}: Logical. Show parallel setup messages
#'   }
#' @param sg0.name Character. Label for subgroup 0 (default: "Not recommend").
#' @param sg1.name Character. Label for subgroup 1 (default: "Recommend").
#' @param details Logical. Print progress details (default: FALSE).
#'
#' @return List with components:
#' \describe{
#'   \item{resCV}{Data frame with CV predictions for each observation}
#'   \item{cv_args}{Arguments used for CV ForestSearch calls}
#'   \item{timing_minutes}{Execution time in minutes}
#'   \item{prop_SG_found}{Percentage of folds where a subgroup was found}
#'   \item{sg_analysis}{Original subgroup definition from full-data analysis}
#'   \item{sg0.name, sg1.name}{Subgroup labels}
#'   \item{Kfolds}{Number of folds used}
#'   \item{sens_summary}{Named vector of sensitivity metrics (sens_H, sens_Hc, ppv_H, ppv_Hc)}
#'   \item{find_summary}{Named vector of subgroup-finding metrics (Any, Exact, etc.)}
#' }
#'
#' @examples
#' \dontrun{
#' # Run initial ForestSearch
#' fs_result <- forestsearch(
#'   df.analysis = trial_data,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   confounders.name = c("age", "biomarker")
#' )
#'
#' # Run 10-fold cross-validation
#' cv_results <- forestsearch_Kfold(
#'   fs.est = fs_result,
#'   Kfolds = 10,
#'   parallel_args = list(plan = "multisession", workers = 4),
#'   details = TRUE
#' )
#'
#' # Summarize results
#' cv_summary <- forestsearch_KfoldOut(cv_results, outall = TRUE)
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for initial subgroup identification
#' \code{\link{forestsearch_KfoldOut}} for summarizing CV results
#' \code{\link{forestsearch_tenfold}} for repeated K-fold simulations
#'
#' @importFrom future plan
#' @importFrom data.table data.table copy
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @export

forestsearch_Kfold <- function(
    fs.est,
    Kfolds = nrow(fs.est$df.est),
    seedit = 8316951L,
    parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE),
    sg0.name = "Not recommend",
    sg1.name = "Recommend",
    details = FALSE
) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION
  # ===========================================================================

  if (!is.list(fs.est)) {
    stop("fs.est must be a list (ForestSearch results object)")
  }

  required_components <- c("df.est", "args_call_all", "sg.harm")
  missing <- setdiff(required_components, names(fs.est))
  if (length(missing) > 0) {
    stop("fs.est missing required components: ", paste(missing, collapse = ", "))
  }

  if (is.null(fs.est$df.est) || !is.data.frame(fs.est$df.est) || nrow(fs.est$df.est) == 0) {
    stop("fs.est$df.est must be a non-empty data frame")
  }

  if (!is.numeric(Kfolds) || length(Kfolds) != 1 || Kfolds < 2) {
    stop("Kfolds must be an integer >= 2")
  }

  if (Kfolds > nrow(fs.est$df.est)) {
    stop("Kfolds cannot exceed number of observations (", nrow(fs.est$df.est), ")")
  }

  # ===========================================================================
  # SECTION 2: PARALLEL SETUP
  # ===========================================================================

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  # Resolve parallel arguments
  resolved_args <- resolve_cv_parallel_args(parallel_args, fs.est$args_call_all, details)
  setup_parallel_SGcons(resolved_args)

  t_start <- proc.time()[3]

  # ===========================================================================
  # SECTION 3: DATA PREPARATION
  # ===========================================================================

  fs_args <- fs.est$args_call_all

  # Extract variable names
  outcome.name <- fs_args$outcome.name
  event.name <- fs_args$event.name
  treat.name <- fs_args$treat.name
  id.name <- fs_args$id.name
  est.scale <- fs_args$est.scale
  confounders.name <- fs_args$confounders.name

  get_names <- c(confounders.name, outcome.name, event.name, id.name, treat.name)

  # Prepare analysis data
  dfa <- fs.est$df.est[, c(get_names, "treat.recommend")]
  names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"
  dfnew <- as.data.frame(dfa)

  # ===========================================================================
  # SECTION 4: FOLD ASSIGNMENT
  # ===========================================================================

  df_scrambled <- data.table::copy(dfnew)

  # Scramble only if not doing LOO
  if (Kfolds < nrow(fs.est$df.est)) {
    set.seed(seedit)
    id_sample <- sample(seq_len(nrow(df_scrambled)), replace = FALSE)
    df_scrambled <- df_scrambled[id_sample, ]
  }

  folds <- cut(seq_len(nrow(df_scrambled)), breaks = Kfolds, labels = FALSE)

  if (details) {
    cat("Cross-validation setup:\n")
    cat("  - Observations:", nrow(df_scrambled), "\n")
    cat("  - Folds:", Kfolds, "\n")
    cat("  - Fold sizes (range):", paste(range(table(folds)), collapse = "-"), "\n")
  }

  # ===========================================================================
  # SECTION 5: CONFIGURE CV FORESTSEARCH ARGUMENTS
  # ===========================================================================

  cv_args <- fs_args
  cv_args$parallel_args <- list(plan = "sequential", workers = 1, show_message = FALSE)
  cv_args$details <- FALSE
  cv_args$plot.sg <- FALSE


  # ===========================================================================
  # SECTION 5.1: VERIFY CV FORESTSEARCH PARAMETERS (when details = TRUE)
  # ===========================================================================

  if(details) print_cv_params(cv_args)

  # ===========================================================================
  # SECTION 6: RUN CROSS-VALIDATION (PARALLEL)
  # ===========================================================================

  resCV <- suppressWarnings({foreach::foreach(
    cv_index = seq_len(Kfolds),
    .options.future = list(seed = TRUE),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {

    # Split data
    testIndexes <- which(folds == cv_index, arr.ind = TRUE)
    x.test <- df_scrambled[testIndexes, ]
    x.train <- df_scrambled[-testIndexes, ]

    # Update training data
    cv_args$df.analysis <- x.train

    # Run ForestSearch on training fold
    fs.train <- suppressWarnings(try(do.call(forestsearch, cv_args), TRUE))

    # Process results
    if (!inherits(fs.train, "try-error") && !is.null(fs.train$sg.harm)) {
      # Subgroup found - apply to test data
      sg1 <- fs.train$sg.harm[1]
      sg2 <- if (length(fs.train$sg.harm) > 1) fs.train$sg.harm[2] else NA

      df.test <- get_dfpred(df.predict = x.test, sg.harm = fs.train$sg.harm, version = 2)
      df.test$cvindex <- cv_index
      df.test$sg1 <- sg1
      df.test$sg2 <- sg2

    } else {
      # No subgroup found - default to ITT (recommend all)
      df.test <- x.test
      df.test$cvindex <- cv_index
      df.test$sg1 <- NA
      df.test$sg2 <- NA
      df.test$treat.recommend <- 1.0
    }

    # Return standardized columns
    data.table::data.table(
      df.test[, c(id.name, outcome.name, event.name, treat.name,
                  "treat.recommend", "treat.recommend.original",
                  "cvindex", "sg1", "sg2")]
    )
  }
})
  # ===========================================================================
  # SECTION 7: POST-PROCESSING AND VALIDATION
  # ===========================================================================

  t_end <- proc.time()[3]
  t_min <- (t_end - t_start) / 60

  # Validate results
  if (length(unique(resCV[[id.name]])) != nrow(df_scrambled)) {
    stop("K-fold cross-validation sample size does not equal full analysis dataset")
  }

  chk1 <- range(table(folds))
  chk2 <- range(table(resCV$cvindex))
  if (!all(chk1 == chk2)) {
    stop("Mismatch on range of fold sample sizes")
  }

  # Convert to data.frame
  resCV <- as.data.frame(resCV)

  # Handle est.scale for treatment direction
  if (est.scale == "1/hr") {
    resCV$treat2 <- 1 - resCV[, treat.name]
    treat.name <- "treat2"
    sg0.name <- "Recommend"
    sg1.name <- "Not recommend"
  }

  # Calculate proportion of folds where subgroup was found
  unique_folds_with_sg <- length(unique(resCV$cvindex[!is.na(resCV$sg1) | !is.na(resCV$sg2)]))
  prop_SG_found <- 100 * round(unique_folds_with_sg / Kfolds, 3)

  if (details) {
    cat("\nCross-validation complete:\n")
    cat("  - Time:", round(t_min, 2), "minutes\n")
    cat("  - Subgroup found in", prop_SG_found, "% of folds\n")
  }

  # ===========================================================================
  # SECTION 8: COMPUTE SUMMARY METRICS (align with forestsearch_tenfold output)
  # ===========================================================================

  # Build intermediate result for forestsearch_KfoldOut
  res_for_kfoldout <- list(
    resCV = resCV,
    Kfolds = Kfolds,
    cv_args = cv_args,
    sg_analysis = fs.est$sg.harm,
    sg0.name = sg0.name,
    sg1.name = sg1.name
  )

  # Get metrics from forestsearch_KfoldOut
  kfold_out <- tryCatch(
    forestsearch_KfoldOut(res = res_for_kfoldout, outall = FALSE, details = details),
    error = function(e) {
      if (details) warning("forestsearch_KfoldOut failed: ", e$message)
      NULL
    }
  )

  # Extract sens_summary and find_summary (single run, not median)
  if (!is.null(kfold_out)) {
    sens_summary <- kfold_out$sens_metrics_original
    find_summary <- kfold_out$find_metrics
  } else {
    sens_summary <- c(sens_H = NA, sens_Hc = NA, ppv_H = NA, ppv_Hc = NA)
    find_summary <- c(Any = NA, Exact = NA, `At least 1` = NA, Cov1 = NA,
                      Cov2 = NA, `Cov 1 & 2` = NA, `Cov1 exact` = NA, `Cov2 exact` = NA)
  }

  # ===========================================================================
  # SECTION 9: RETURN RESULTS
  # ===========================================================================

  result <- list(
    resCV = resCV,
    cv_args = cv_args,
    timing_minutes = t_min,
    prop_SG_found = prop_SG_found,
    sg_analysis = fs.est$sg.harm,
    sg0.name = sg0.name,
    sg1.name = sg1.name,
    Kfolds = Kfolds,
    # Add summary outputs to align with forestsearch_tenfold
    sens_summary = sens_summary,
    find_summary = find_summary
  )

  class(result) <- c("fs_kfold", "list")
  return(result)
}


#' ForestSearch Repeated K-Fold Cross-Validation
#'
#' Runs repeated K-fold cross-validation simulations for ForestSearch and
#' summarizes subgroup identification stability across repetitions.
#'
#' @description
#' This function performs multiple independent K-fold cross-validations to
#' assess the variability in subgroup identification. Each simulation:
#' \enumerate{
#'   \item Randomly shuffles the data
#'   \item Performs K-fold CV
#'   \item Records sensitivity and agreement metrics
#' }
#' Results are summarized across all simulations.
#'
#' @section Parallelization Strategy:
#' Unlike the single K-fold function which parallelizes across folds, this
#' function parallelizes across simulations for better efficiency when
#' running many repetitions. Each simulation runs its K-fold CV sequentially.
#'
#' @param fs.est List. ForestSearch results object from \code{\link{forestsearch}}.
#' @param sims Integer. Number of simulation repetitions.
#' @param seed Integer. Base random seed for fold shuffling. Default 8316951L.
#'   Each simulation uses seed + 1000 * ksim for reproducibility.
#' @param Kfolds Integer. Number of folds per simulation (default: 10).
#' @param details Logical. Print progress details (default: TRUE).
#' @param parallel_args List. Parallelization configuration.
#'
#' @return List with components:
#' \describe{
#'   \item{sens_summary}{Named vector of median sensitivity metrics across simulations}
#'   \item{find_summary}{Named vector of median subgroup-finding metrics}
#'   \item{sens_out}{Matrix of sensitivity metrics (sims x metrics)}
#'   \item{find_out}{Matrix of finding metrics (sims x metrics)}
#'   \item{timing_minutes}{Total execution time}
#'   \item{sims}{Number of simulations run}
#'   \item{Kfolds}{Number of folds per simulation}
#' }
#'
#' @examples
#' \dontrun{
#' # Run 100 repetitions of 10-fold CV
#' tenfold_results <- forestsearch_tenfold(
#'   fs.est = fs_result,
#'   sims = 100,
#'   Kfolds = 10,
#'   parallel_args = list(plan = "multisession", workers = 6),
#'   details = TRUE
#' )
#'
#' # View summary
#' print(tenfold_results$sens_summary)
#' print(tenfold_results$find_summary)
#' }
#'
#' @seealso
#' \code{\link{forestsearch_Kfold}} for single K-fold CV
#' \code{\link{forestsearch_KfoldOut}} for summarizing CV results
#'
#' @importFrom future plan
#' @importFrom data.table data.table copy
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @export

forestsearch_tenfold <- function(
    fs.est,
    sims,
    Kfolds = 10,
    details = TRUE,
    seed = 8316951L,
    parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE)
) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION
  # ===========================================================================

  if (!is.list(fs.est)) {
    stop("fs.est must be a list (ForestSearch results object)")
  }

  required_components <- c("df.est", "args_call_all", "sg.harm")
  missing <- setdiff(required_components, names(fs.est))
  if (length(missing) > 0) {
    stop("fs.est missing required components: ", paste(missing, collapse = ", "))
  }

  if (!is.numeric(sims) || length(sims) != 1 || sims < 1) {
    stop("sims must be a positive integer")
  }

  if (!is.numeric(Kfolds) || length(Kfolds) != 1 || Kfolds < 2) {
    stop("Kfolds must be an integer >= 2")
  }

  # ===========================================================================
  # SECTION 2: PARALLEL SETUP
  # ===========================================================================

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  resolved_args <- resolve_cv_parallel_args(parallel_args, fs.est$args_call_all, details)
  setup_parallel_SGcons(resolved_args)

  t_start <- proc.time()[3]

  if (details) {
    cat("Starting repeated K-fold cross-validation:\n")
    cat("  - Simulations:", sims, "\n")
    cat("  - Folds per simulation:", Kfolds, "\n")
    cat("  - Workers:", resolved_args$workers, "\n")
  }

  # ===========================================================================
  # SECTION 3: DATA PREPARATION
  # ===========================================================================

  fs_args <- fs.est$args_call_all

  # Extract variable names
  outcome.name <- fs_args$outcome.name
  event.name <- fs_args$event.name
  treat.name <- fs_args$treat.name
  id.name <- fs_args$id.name
  confounders.name <- fs_args$confounders.name

  get_names <- c(confounders.name, outcome.name, event.name, id.name, treat.name)

  # Prepare base data
  dfa <- fs.est$df.est[, c(get_names, "treat.recommend")]
  names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"
  dfnew <- as.data.frame(dfa)

  # Configure CV arguments (sequential within each simulation)
  cv_args <- fs_args
  cv_args$parallel_args <- list(plan = "sequential", workers = 1, show_message = FALSE)
  cv_args$details <- FALSE
  cv_args$plot.sg <- FALSE

  # ===========================================================================
  # SECTION 3.1: VERIFY CV FORESTSEARCH PARAMETERS (when details = TRUE)
  # ===========================================================================

  if (details) print_cv_params(cv_args)

  # ===========================================================================
  # SECTION 4: RUN SIMULATIONS (PARALLEL ACROSS SIMS)
  # ===========================================================================

  simulation_results <- suppressWarnings({foreach::foreach(
    ksim = seq_len(sims),
    .options.future = list(seed = TRUE),
    .errorhandling = "pass"
  ) %dofuture% {

    # Shuffle data for this simulation
    df_scrambled <- data.table::copy(dfnew)
    set.seed(seed + 1000 * ksim)
    id_sample <- sample(seq_len(nrow(df_scrambled)), replace = FALSE)
    df_scrambled <- df_scrambled[id_sample, ]

    # Create fold assignments
    folds <- cut(seq_len(nrow(df_scrambled)), breaks = Kfolds, labels = FALSE)

    # Run K-fold CV sequentially within this simulation
    resCV_list <- vector("list", Kfolds)

    for (cv_index in seq_len(Kfolds)) {
      testIndexes <- which(folds == cv_index, arr.ind = TRUE)
      x.test <- df_scrambled[testIndexes, ]
      x.train <- df_scrambled[-testIndexes, ]

      cv_args$df.analysis <- x.train

      fs.train <- suppressWarnings(try(do.call(forestsearch, cv_args), TRUE))

      if (!inherits(fs.train, "try-error") && !is.null(fs.train$sg.harm)) {
        sg1 <- fs.train$sg.harm[1]
        sg2 <- if (length(fs.train$sg.harm) > 1) fs.train$sg.harm[2] else NA

        df.test <- get_dfpred(df.predict = x.test, sg.harm = fs.train$sg.harm, version = 2)
        df.test$cvindex <- cv_index
        df.test$sg1 <- sg1
        df.test$sg2 <- sg2

      } else {
        df.test <- x.test
        df.test$cvindex <- cv_index
        df.test$sg1 <- NA
        df.test$sg2 <- NA
        df.test$treat.recommend <- 1.0
      }

      resCV_list[[cv_index]] <- df.test[, c(id.name, outcome.name, event.name, treat.name,
                                            "treat.recommend", "treat.recommend.original",
                                            "cvindex", "sg1", "sg2")]
    }

    resCV <- do.call(rbind, resCV_list)
    resCV <- as.data.frame(resCV)

    # Summarize this simulation
    res <- list(
      resCV = resCV,
      Kfolds = Kfolds,
      confounders.name = confounders.name,
      cv_args = cv_args,
      outcome.name = outcome.name,
      event.name = event.name,
      id.name = id.name,
      treat.name = treat.name,
      sg_analysis = fs.est$sg.harm,
      sg0.name = "Not recommend",
      sg1.name = "Recommend"
    )

    out <- forestsearch_KfoldOut(res = res, outall = FALSE, details = FALSE)

    list(
      sens_metrics_original = out$sens_metrics_original,
      find_metrics = out$find_metrics,
      sim_id = ksim
    )
  }
})

  # ===========================================================================
  # SECTION 5: COMBINE AND SUMMARIZE RESULTS
  # ===========================================================================

  # Extract metrics from successful simulations
  valid_results <- Filter(function(x) !inherits(x, "error") && is.list(x), simulation_results)

  if (length(valid_results) == 0) {
    stop("All simulations failed. Check ForestSearch configuration.")
  }

  if (length(valid_results) < sims && details) {
    warning(sims - length(valid_results), " simulations failed and were excluded.")
  }

  sens_out <- do.call(rbind, lapply(valid_results, `[[`, "sens_metrics_original"))
  find_out <- do.call(rbind, lapply(valid_results, `[[`, "find_metrics"))

  # Compute summaries
  sens_summary <- apply(sens_out, 2, median, na.rm = TRUE)
  find_summary <- apply(find_out, 2, median, na.rm = TRUE)

  t_end <- proc.time()[3]
  t_min <- (t_end - t_start) / 60

  if (details) {
    cat("\nRepeated K-fold CV complete:\n")
    cat("  - Time:", round(t_min, 2), "minutes\n")
    cat("  - Successful simulations:", length(valid_results), "/", sims, "\n")
    cat("  - Projected hours per 100 sims:", round((t_min / 60) * (100 / sims), 2), "\n")
  }

  # ===========================================================================
  # SECTION 6: RETURN RESULTS
  # ===========================================================================

  result <- list(
    sens_summary = sens_summary,
    find_summary = find_summary,
    sens_out = sens_out,
    find_out = find_out,
    timing_minutes = t_min,
    sims = length(valid_results),
    Kfolds = Kfolds
  )

  class(result) <- c("fs_tenfold", "list")
  return(result)
}


#' ForestSearch K-Fold Cross-Validation Output Summary
#'
#' Summarizes cross-validation results for ForestSearch, including subgroup
#' agreement and performance metrics.
#'
#' @param res List. Result object from ForestSearch cross-validation, must contain
#'   elements: \code{cv_args}, \code{sg_analysis}, \code{sg0.name}, \code{sg1.name},
#'   \code{Kfolds}, \code{resCV}.
#' @param details Logical. Print details during execution (default: FALSE).
#' @param outall Logical. If TRUE, returns all summary tables; if FALSE, returns
#'   only metrics (default: FALSE).
#'
#' @return If \code{outall=FALSE}, a list with \code{sens_metrics_original} and
#'   \code{find_metrics}. If \code{outall=TRUE}, a list with summary tables and metrics.
#'
#' @importFrom stringr str_sub str_length
#' @export

forestsearch_KfoldOut <- function(res, details = FALSE, outall = FALSE) {

  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================

  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' required: If issues loading consider substituting
         'stringr::str_sub(x,start,end)' with 'substr(x,start,end)' and
         'stringr::str_length(x)' with 'nchar(x)'.")
  }

  # List of required elements in res
  required_elements <- c("cv_args", "sg_analysis", "sg0.name", "sg1.name", "Kfolds", "resCV")
  missing_elements <- required_elements[!sapply(required_elements, function(x) !is.null(res[[x]]))]
  if (length(missing_elements) > 0) {
    stop("The following elements are missing in 'res': ", paste(missing_elements, collapse = ", "))
  }

  # For nested elements in cv_args
  cv_args_elements <- c("confounders.name", "outcome.name", "event.name", "treat.name")
  missing_cv_args <- cv_args_elements[!sapply(cv_args_elements, function(x) !is.null(res$cv_args[[x]]))]
  if (length(missing_cv_args) > 0) {
    stop("The following elements are missing in 'res$cv_args': ", paste(missing_cv_args, collapse = ", "))
  }

  # ===========================================================================
  # EXTRACT CONFIGURATION
  # ===========================================================================

  confounders.name <- res$cv_args$confounders.name
  outcome.name <- res$cv_args$outcome.name
  event.name <- res$cv_args$event.name
  treat.name <- res$cv_args$treat.name
  sg_analysis <- res$sg_analysis
  sg0.name <- res$sg0.name
  sg1.name <- res$sg1.name
  est.scale <- if (!is.null(res$cv_args$est.scale)) res$cv_args$est.scale else "hr"

  Kfolds <- res$Kfolds
  df_CV <- res$resCV

  # ===========================================================================
  # EXTRACT SUBGROUPS FROM EACH FOLD
  # ===========================================================================

  sg1 <- sg2 <- rep(NA, Kfolds)
  for (ks in seq_len(Kfolds)) {
    fold_data <- subset(df_CV, cvindex == ks)
    if (nrow(fold_data) > 0) {
      sg1[ks] <- fold_data$sg1[1]
      sg2[ks] <- fold_data$sg2[1]
    }
  }
  SGs_found <- cbind(sg1, sg2)

  # ===========================================================================
  # COMPUTE CV SUMMARY METRICS
  # ===========================================================================

  CV_summary <- CV_sgs(
    sg1 = sg1,
    sg2 = sg2,
    confs = confounders.name,
    sg_analysis = sg_analysis
  )

  if (details) {
    cat("Any found:", mean(CV_summary$any_found), "\n")
    cat("Exact match:", mean(CV_summary$exact_match), "\n")
    cat("At least 1 match:", mean(CV_summary$one_match), "\n")
    cat("Cov 1 any:", mean(CV_summary$cov1_any), "\n")
    cat("Cov 2 any:", mean(CV_summary$cov2_any, na.rm = TRUE), "\n")
    cat("Cov 1 and 2 any:", mean(ifelse(CV_summary$cov1_any & CV_summary$cov2_any, 1, 0), na.rm = TRUE), "\n")
    cat("Cov 1 exact:", mean(CV_summary$cov1_exact), "\n")
    cat("Cov 2 exact:", mean(CV_summary$cov2_exact, na.rm = TRUE), "\n")
  }

  find_metrics <- c(
    mean(CV_summary$any_found),
    mean(CV_summary$exact_match),
    mean(CV_summary$one_match),
    mean(CV_summary$cov1_any),
    mean(CV_summary$cov2_any, na.rm = TRUE),
    mean(ifelse(CV_summary$cov1_any & CV_summary$cov2_any, 1, 0), na.rm = TRUE),
    mean(CV_summary$cov1_exact),
    mean(CV_summary$cov2_exact, na.rm = TRUE)
  )
  names(find_metrics) <- c("Any", "Exact", "At least 1", "Cov1", "Cov2",
                           "Cov 1 & 2", "Cov1 exact", "Cov2 exact")

  # ===========================================================================
  # COMPUTE SENSITIVITY/PPV METRICS
  # ===========================================================================

  n_sgfound <- length(unique(df_CV$treat.recommend))

  if (n_sgfound == 2) {
    tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
    sensH <- tabit[1, 1] / sum(tabit[, 1])
    sensHc <- tabit[2, 2] / sum(tabit[, 2])
    ppvH <- tabit[1, 1] / sum(tabit[1, ])
    ppvHc <- tabit[2, 2] / sum(tabit[2, ])
  } else {
    tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
    sensH <- if (nrow(tabit) > 0 && ncol(tabit) > 0) tabit[1, 1] / sum(tabit[, 1]) else NA
    sensHc <- NA
    ppvH <- if (nrow(tabit) > 0) tabit[1, 1] / sum(tabit[1, ]) else NA
    ppvHc <- NA
  }

  sens_metrics_original <- c(sensH, sensHc, ppvH, ppvHc)
  names(sens_metrics_original) <- c("sens_H", "sens_Hc", "ppv_H", "ppv_Hc")

  if (details) {
    cat("Agreement (sens, ppv) in H and Hc:", c(sensH, sensHc, ppvH, ppvHc), "\n")
  }

  # ===========================================================================
  # RETURN RESULTS
  # ===========================================================================

  if (outall) {
    # Generate detailed tables
    itt_tab <- SG_tab_estimates(
      df = as.data.frame(df_CV),
      SG_flag = "ITT",
      draws = 0,
      details = FALSE,
      outcome.name = outcome.name,
      event.name = event.name,
      treat.name = treat.name,
      strata.name = NULL,
      potentialOutcome.name = NULL,
      est.scale = est.scale,
      sg0_name = "Questionable",
      sg1_name = "Recommend"
    )

    if (n_sgfound == 2) {
      SG_tab_Kfold <- SG_tab_estimates(
        df = as.data.frame(df_CV),
        SG_flag = "treat.recommend",
        sg1_name = sg1.name,
        sg0_name = sg0.name,
        outcome.name = outcome.name,
        event.name = event.name,
        treat.name = treat.name,
        strata.name = NULL,
        draws = 0,
        details = FALSE
      )
    } else {
      SG_tab_Kfold <- itt_tab
    }

    SG_tab_original <- SG_tab_estimates(
      df = as.data.frame(df_CV),
      SG_flag = "treat.recommend.original",
      sg1_name = sg1.name,
      sg0_name = sg0.name,
      outcome.name = outcome.name,
      event.name = event.name,
      treat.name = treat.name,
      strata.name = NULL,
      draws = 0,
      details = FALSE
    )

    # Combine tables
    estimates_toget <- c("Subgroup", "n", "n1", "m1", "m0", "RMST", "HR (95% CI)")

    if (n_sgfound == 2) {
      temp1 <- itt_tab[estimates_toget]
      temp2a <- SG_tab_original[1, estimates_toget]
      temp2b <- SG_tab_original[2, estimates_toget]
      temp3a <- SG_tab_Kfold[1, estimates_toget]
      temp3b <- SG_tab_Kfold[2, estimates_toget]
      tab_all <- rbind(temp1, temp2a, temp3a, temp2b, temp3b)
      colnames(tab_all) <- c("Subgroup", "n", "n1", "m1", "m0", "RMST", "Hazard ratio")
      rownames(tab_all) <- c("Overall", "FA_0", "KfA_0", "FA_1", "KfA_1")
    } else {
      temp1 <- itt_tab[estimates_toget]
      temp2a <- SG_tab_original[1, estimates_toget]
      temp2b <- if (nrow(SG_tab_original) > 1) SG_tab_original[2, estimates_toget] else NULL
      tab_all <- rbind(temp1, temp2a, temp2b)
      colnames(tab_all) <- c("Subgroup", "n", "n1", "m1", "m0", "RMST", "Hazard ratio")
      rownames(tab_all) <- if (is.null(temp2b)) c("Overall", "FA_0") else c("Overall", "FA_0", "FA_1")
    }

    if (details) print(tab_all)

    out <- list(
      itt_tab = itt_tab,
      SG_tab_original = SG_tab_original,
      SG_tab_Kfold = SG_tab_Kfold,
      CV_summary = CV_summary,
      sens_metrics_original = sens_metrics_original,
      find_metrics = find_metrics,
      SGs_found = SGs_found,
      tab_all = tab_all
    )
  } else {
    out <- list(
      sens_metrics_original = sens_metrics_original,
      find_metrics = find_metrics
    )
  }

  return(out)
}


#' Cross-Validation Subgroup Match Summary
#'
#' Summarizes the match between cross-validation subgroups and analysis subgroups.
#'
#' @param sg1 Character vector. Subgroup 1 labels for each fold.
#' @param sg2 Character vector. Subgroup 2 labels for each fold.
#' @param confs Character vector. Confounder names.
#' @param sg_analysis Character vector. Subgroup analysis labels.
#'
#' @return List with indicators for any match, exact match, one match, and
#'   covariate-specific matches.
#'
#' @importFrom stringr str_sub str_length
#' @export

CV_sgs <- function(sg1, sg2, confs, sg_analysis) {

  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' required.")
  }

  any_found <- ifelse(!is.na(sg1) | !is.na(sg2), 1, 0)
  sg_depth <- length(sg_analysis)

  if (sg_depth == 2) {
    sg1a <- sg_analysis[1]
    sg2a <- sg_analysis[2]

    # Exact match on both to analysis data
    exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a) & (sg1 == sg2a | sg2 == sg2a), 1, 0)
    exact_match[is.na(exact_match)] <- 0.0

    # At least 1 exact match
    one_match <- ifelse((sg1 == sg1a | sg2 == sg1a) | (sg1 == sg2a | sg2 == sg2a), 1, 0)
    one_match[is.na(one_match)] <- 0.0

    # Cov 1 exact
    cov1_match <- ifelse((sg1 == sg1a | sg2 == sg1a), 1, 0)
    cov1_match[is.na(cov1_match)] <- 0.0

    cov2_match <- ifelse((sg1 == sg2a | sg2 == sg2a), 1, 0)
    cov2_match[is.na(cov2_match)] <- 0.0

    # Find confounder names involved in sg1a and sg2a
    cov1_any <- find_covariate_any_match(sg1a, sg1, sg2, confs)
    cov2_any <- find_covariate_any_match(sg2a, sg1, sg2, confs)

  } else if (sg_depth == 1) {
    sg1a <- sg_analysis[1]

    # Exact match
    exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a), 1, 0)
    exact_match[is.na(exact_match)] <- 0.0

    one_match <- exact_match
    cov1_match <- exact_match
    cov2_match <- NA

    # Find covariate any match
    cov1_any <- find_covariate_any_match(sg1a, sg1, sg2, confs)
    cov2_any <- NA

  } else {
    # No subgroup in analysis
    exact_match <- rep(0, length(sg1))
    one_match <- rep(0, length(sg1))
    cov1_match <- rep(0, length(sg1))
    cov2_match <- NA
    cov1_any <- rep(0, length(sg1))
    cov2_any <- NA
  }

  return(list(
    any_found = any_found,
    exact_match = exact_match,
    one_match = one_match,
    cov1_any = cov1_any,
    cov2_any = cov2_any,
    cov1_exact = cov1_match,
    cov2_exact = cov2_match
  ))
}


#' Find Covariate Any Match
#'
#' Helper function to determine if any CV fold found a subgroup involving
#' the same covariate (not necessarily same cut).
#'
#' @param sg_target Character. Target subgroup definition to match.
#' @param sg1 Character vector. Subgroup 1 labels for each fold.
#' @param sg2 Character vector. Subgroup 2 labels for each fold.
#' @param confs Character vector. Confounder names.
#'
#' @return Numeric vector (0/1) indicating match for each fold.
#'
#' @keywords internal

find_covariate_any_match <- function(sg_target, sg1, sg2, confs) {

  if (is.na(sg_target) || is.null(sg_target)) {
    return(rep(0, length(sg1)))
  }

  # Determine prefix ({ or !{)
  dda <- charmatch("{", sg_target, nomatch = 0)
  ddb <- charmatch("!{", sg_target, nomatch = 0)

  if (dda == 1) {
    aa <- rep("{", length(confs))
  } else if (ddb == 1) {
    aa <- rep("!{", length(confs))
  } else {
    return(rep(0, length(sg1)))
  }

  # Build search patterns
  temp <- paste0(aa, confs)
  temp <- paste0(temp, "}")

  # Find which confounder is in sg_target
  loc_name <- charmatch(temp, sg_target)
  index_name <- which(loc_name == 1)

  if (length(index_name) == 0) {
    # Try without closing brace
    temp <- paste0(aa, confs)
    loc_name <- charmatch(temp, sg_target)
    index_name <- which(loc_name == 1)
  }

  if (length(index_name) == 0) {
    return(rep(0, length(sg1)))
  }

  # Handle case where multiple confounders match (e.g., "z1" and "z11")
  if (length(index_name) > 1) {
    confs2 <- confs[index_name]
    lc <- stringr::str_length(confs2)
    if (dda == 1) {
      ctoget <- stringr::str_sub(sg_target, 2, max(lc))
    } else {
      ctoget <- stringr::str_sub(sg_target, 3, max(lc))
    }
    itoget <- which(confs == ctoget)
    cfs <- confs[itoget]
  } else {
    cfs <- confs[which(loc_name == 1)]
  }

  if (length(cfs) == 0) {
    return(rep(0, length(sg1)))
  }

  # Check if covariate appears in any CV fold subgroups
  bb1 <- grepl(cfs, sg1)
  bb2 <- grepl(cfs, sg2)
  cov_any <- ifelse(bb1 | bb2, 1, 0)

  return(cov_any)
}


#' Resolve Parallel Arguments for Cross-Validation
#'
#' Helper function to resolve and validate parallel processing arguments,
#' similar to bootstrap's \code{resolve_bootstrap_parallel_args}.
#'
#' @param parallel_args List. User-provided parallel arguments.
#' @param fs_args List. Original ForestSearch call arguments.
#' @param details Logical. Print configuration messages.
#'
#' @return List with resolved plan, workers, show_message.
#'
#' @keywords internal

resolve_cv_parallel_args <- function(parallel_args, fs_args, details = FALSE) {

  # Default values
  default_args <- list(
    plan = "multisession",
    workers = max(1, parallel::detectCores() - 1),
    show_message = TRUE
  )

  # Merge user args with defaults
  resolved_args <- default_args

  if (length(parallel_args) > 0) {
    for (nm in names(parallel_args)) {
      resolved_args[[nm]] <- parallel_args[[nm]]
    }
  } else if (!is.null(fs_args$parallel_args) && length(fs_args$parallel_args) > 0) {
    # Fall back to ForestSearch settings
    for (nm in names(fs_args$parallel_args)) {
      resolved_args[[nm]] <- fs_args$parallel_args[[nm]]
    }
    if (details) {
      message("CV parallel config: Using ForestSearch analysis settings")
    }
  }

  # Validate workers
  max_cores <- parallel::detectCores()
  if (resolved_args$workers > max_cores) {
    warning("Requested workers (", resolved_args$workers, ") exceeds available cores (",
            max_cores, "). Using ", max_cores - 1, " workers.")
    resolved_args$workers <- max(1, max_cores - 1)
  }

  if (details && resolved_args$show_message) {
    message("CV will use: ", resolved_args$workers, " workers with '",
            resolved_args$plan, "' plan")
  }

  resolved_args
}


# ==============================================================================
# PRINT METHODS
# ==============================================================================

#' Print Method for K-Fold CV Results
#'
#' @param x An fs_kfold object
#' @param ... Additional arguments (ignored)
#' @export

print.fs_kfold <- function(x, ...) {
  cat("ForestSearch K-Fold Cross-Validation Results\n")
  cat("=============================================\n")
  cat("Folds:", x$Kfolds, "\n")
  cat("Observations:", nrow(x$resCV), "\n")
  cat("Subgroup found in:", x$prop_SG_found, "% of folds\n")
  cat("Original subgroup:", paste(x$sg_analysis, collapse = " & "), "\n")
  cat("Time:", round(x$timing_minutes, 2), "minutes\n")
  cat("\nUse forestsearch_KfoldOut() for detailed metrics.\n")
  invisible(x)
}


#' Print Method for Repeated K-Fold CV Results
#'
#' @param x An fs_tenfold object
#' @param ... Additional arguments (ignored)
#' @export

print.fs_tenfold <- function(x, ...) {
  cat("ForestSearch Repeated K-Fold Cross-Validation Results\n")
  cat("======================================================\n")
  cat("Simulations:", x$sims, "\n")
  cat("Folds per simulation:", x$Kfolds, "\n")
  cat("Time:", round(x$timing_minutes, 2), "minutes\n")
  cat("\nSensitivity Summary (median across simulations):\n")
  print(round(x$sens_summary, 3))
  cat("\nSubgroup Finding Summary (median across simulations):\n")
  print(round(x$find_summary, 3))
  invisible(x)
}

#' Print CV ForestSearch Parameters
#' @keywords internal
print_cv_params <- function(cv_args) {
  cat("\nForestSearch parameters for CV folds:\n")
  cat("  - sg_focus:", cv_args$sg_focus, "\n")
  cat("  - maxk:", cv_args$maxk, "\n")
  cat("  - fs.splits:", cv_args$fs.splits, "\n")
  cat("  - max_subgroups_search:", cv_args$max_subgroups_search, "\n")
  cat("  - hr.threshold:", cv_args$hr.threshold, "\n")
  cat("  - hr.consistency:", cv_args$hr.consistency, "\n")
  cat("  - pconsistency.threshold:", cv_args$pconsistency.threshold, "\n")
  cat("  - n.min:", cv_args$n.min, "\n")
  cat("  - use_twostage:", cv_args$use_twostage, "\n")
  if (isTRUE(cv_args$use_twostage) && length(cv_args$twostage_args) > 0) {
    cat("  - twostage_args:\n")
    cat("      n.splits.screen:", cv_args$twostage_args$n.splits.screen, "\n")
    cat("      batch.size:", cv_args$twostage_args$batch.size, "\n")
  }
  cat("  - use_lasso:", cv_args$use_lasso, "\n")
  cat("  - use_grf:", cv_args$use_grf, "\n")
  cat("  - (per-fold parallel: sequential)\n")
  cat("  - (per-fold details: FALSE)\n")
  cat("  - (per-fold plot.sg: FALSE)\n")
}
