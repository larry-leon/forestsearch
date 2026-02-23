# =============================================================================
# MRCT Simulation Functions for ForestSearch Package
# =============================================================================
#
# Multi-Regional Clinical Trial (MRCT) simulation framework for evaluating
# subgroup identification methods using ForestSearch. Supports scenarios with
# heterogeneous treatment effects across regions.
#
# Key functions:
#   - mrct_region_sims(): Main simulation function
#   - create_dgm_for_mrct(): DGM creation wrapper for MRCT scenarios
#   - summaryout_mrct(): Summary table generation
#   - SGplot_estimates(): Visualization of HR estimates
#
# =============================================================================


# =============================================================================
# Main Simulation Function: mrct_region_sims
# =============================================================================

#' MRCT Regional Subgroup Simulation
#'
#' Simulates multi-regional clinical trials and evaluates ForestSearch subgroup
#' identification. Splits data by region into training and testing populations,
#' identifies subgroups using ForestSearch on training data, and evaluates
#' performance on the testing region.
#'
#' @param dgm Data generating mechanism object from \code{\link{generate_aft_dgm_flex}}
#' @param n_sims Integer. Number of simulations to run
#' @param n_sample Integer. Sample size per simulation. If NULL (default), uses
#'   the entire super-population from dgm
#' @param region_var Character. Name of the region indicator variable used to split
#'   data into training (region_var == 0) and testing (region_var == 1) populations.
#'   Default: "z_regA"
#' @param sg_focus Character. Subgroup selection criterion passed to
#'   \code{\link{forestsearch}}: "minSG", "hr", or "maxSG". Default: "minSG"
#' @param maxk Integer. Maximum number of factors in subgroup combinations (1 or 2).
#'   Default: 1
#' @param hr.threshold Numeric. Hazard ratio threshold for subgroup identification.
#'   Default: 0.90
#' @param hr.consistency Numeric. Consistency threshold for hazard ratio.
#'   Default: 0.80
#' @param pconsistency.threshold Numeric. Probability threshold for consistency.
#'   Default: 0.90
#' @param confounders.name Character vector. Confounder variable names for ForestSearch.
#'   If NULL, automatically extracted from dgm
#' @param conf_force Character vector. Forced cuts to consider in ForestSearch.
#'   Default: c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5")
#' @param analysis_time Numeric. Time of analysis for administrative censoring.
#'   Default: 60
#' @param cens_adjust Numeric. Adjustment factor for censoring rate on log scale.
#'   Default: 0
#' @param parallel_args List. Parallel processing configuration with components:
#'   \itemize{
#'     \item plan: "multisession", "multicore", "callr", or "sequential"
#'     \item workers: Number of workers (NULL for auto-detect)
#'     \item show_message: Logical for progress messages
#'   }
#' @param details Logical. Print detailed progress information. Default: FALSE
#' @param seed Integer. Base random seed for reproducibility. Default: NULL
#'
#' @return A data.table with simulation results containing:
#' \describe{
#'   \item{sim}{Simulation index}
#'   \item{n_itt}{ITT sample size}
#'   \item{hr_itt}{ITT hazard ratio (stratified if strat variable present)}
#'   \item{hr_ittX}{ITT hazard ratio stratified by region}
#'   \item{n_train}{Training (non-region A) sample size}
#'   \item{hr_train}{Training population hazard ratio}
#'   \item{n_test}{Testing (region A) sample size}
#'   \item{hr_test}{Testing population hazard ratio}
#'   \item{any_found}{Indicator: 1 if subgroup identified, 0 otherwise}
#'   \item{sg_found}{Character description of identified subgroup}
#'   \item{n_sg}{Subgroup sample size}
#'   \item{hr_sg}{Subgroup hazard ratio}
#'   \item{POhr_sg}{Potential outcome hazard ratio in subgroup}
#'   \item{prev_sg}{Subgroup prevalence (proportion of testing population)}
#'   \item{hr_sg_null}{Subgroup HR when found, NA otherwise}
#' }
#'
#' @details
#' ## Simulation Process
#'
#' For each simulation:
#' 1. Sample from super-population using \code{\link{simulate_from_dgm}}
#' 2. Split by region_var into training and testing populations
#' 3. Estimate HRs in ITT, training, and testing populations
#' 4. Run \code{\link{forestsearch}} on training population
#' 5. Apply identified subgroup to testing population
#' 6. Calculate subgroup-specific estimates
#'
#' ## Region Variable
#'
#' The region_var parameter is used ONLY for splitting data into training/testing
#' populations. It does not imply any prognostic effect. To include prognostic
#' confounder effects, specify them when creating the DGM using
#' \code{\link{create_dgm_for_mrct}} or \code{\link{generate_aft_dgm_flex}}.
#'
#' @examples
#' \dontrun{
#' # Create DGM for alternative hypothesis
#' dgm_alt <- create_dgm_for_mrct(
#'   df_case = df_case,
#'   model_type = "alt",
#'   log_hrs = log(c(3, 1.25, 0.50)),
#'   verbose = TRUE
#' )
#'
#' # Run simulations
#' results <- mrct_region_sims(
#'   dgm = dgm_alt,
#'   n_sims = 100,
#'   region_var = "z_regA",
#'   sg_focus = "minSG",
#'   parallel_args = list(plan = "multisession", workers = 4),
#'   details = TRUE
#' )
#'
#' # Summarize results
#' cat("Subgroup identification rate:", mean(results$any_found) * 100, "%\n")
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for subgroup identification algorithm
#' \code{\link{generate_aft_dgm_flex}} for DGM creation
#' \code{\link{simulate_from_dgm}} for data simulation
#' \code{\link{create_dgm_for_mrct}} for MRCT-specific DGM wrapper
#' \code{\link{summaryout_mrct}} for summarizing simulation results
#'
#' @importFrom survival coxph Surv
#' @importFrom data.table data.table as.data.table fifelse copy
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @importFrom future plan sequential multisession multicore nbrOfWorkers
#' @importFrom progressr progressor handlers with_progress
#' @export

mrct_region_sims <- function(
    dgm,
    n_sims,
    n_sample = NULL,
    region_var = "z_regA",
    sg_focus = "minSG",
    maxk = 1,
    hr.threshold = 0.90,
    hr.consistency = 0.80,
    pconsistency.threshold = 0.90,
    confounders.name = NULL,
    conf_force = NULL,
    analysis_time = 60,
    cens_adjust = 0,
    parallel_args = list(plan = "multisession", workers = NULL, show_message = TRUE),
    details = FALSE,
    seed = NULL
) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION
  # ===========================================================================

  if (!inherits(dgm, c("aft_dgm_flex", "aft_dgm"))) {
    stop("dgm must be an object created by generate_aft_dgm_flex()")
  }

  if (!is.numeric(n_sims) || n_sims < 1) {
    stop("n_sims must be a positive integer")
  }

  if (!region_var %in% names(dgm$df_super)) {
    stop("region_var '", region_var, "' not found in dgm$df_super")
  }

  t_start <- proc.time()[3]

  # ===========================================================================
  # SECTION 2: SET DEFAULTS
  # ===========================================================================

  # Default sample size from dgm

  if (is.null(n_sample)) {
    n_sample <- nrow(dgm$df_super)
  }

  # Default confounders from dgm
  if (is.null(confounders.name)) {
    all_vars <- names(dgm$df_super)
    exclude_vars <- c("id", "treat", "tte", "event", "y_sim", "treat_sim",
                      "event_sim", "t_true", "c_time", "lin_pred_0", "lin_pred_1",
                      "loghr_po", "entrytime", "flag_harm")
    confounders.name <- setdiff(
      all_vars[grepl("^z_|^ecog|^strat", all_vars)],
      exclude_vars
    )
  }

  # Default forced cuts
  if (is.null(conf_force)) {
    conf_force <- c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5")
  }

  # ===========================================================================
  # SECTION 3: SETUP PARALLEL PROCESSING
  # ===========================================================================

  plan_type <- null_or(parallel_args$plan, "multisession")
  workers <- parallel_args$workers
  show_message <- null_or(parallel_args$show_message, TRUE)

  if (plan_type == "sequential") {
    future::plan(future::sequential)
  } else {
    if (is.null(workers)) {
      workers <- max(1, parallel::detectCores() - 1)
    }

    if (plan_type == "multisession") {
      future::plan(future::multisession, workers = workers)
    } else if (plan_type == "callr") {
      if (requireNamespace("future.callr", quietly = TRUE)) {
        future::plan(future.callr::callr, workers = workers)
      } else {
        warning("future.callr not available, using multisession")
        future::plan(future::multisession, workers = workers)
      }
    } else if (plan_type == "multicore") {
      if (.Platform$OS.type == "unix") {
        future::plan(future::multicore, workers = workers)
      } else {
        future::plan(future::multisession, workers = workers)
      }
    }

    if (show_message) {
      message(sprintf("Parallel processing: %d workers using %s backend",
                      workers, class(future::plan())[1]))
    }
  }

  # ===========================================================================
  # SECTION 4: RUN SIMULATIONS
  # ===========================================================================

  progressr::handlers(global = TRUE)
  progressr::handlers("progress")

  results <- progressr::with_progress({
    p <- progressr::progressor(along = seq_len(n_sims))

    foreach::foreach(
      sim = seq_len(n_sims),
      .options.future = list(seed = TRUE),
      .combine = "rbind",
      .errorhandling = "pass",
      .packages = c("survival", "data.table", "forestsearch")
    ) %dofuture% {

      p(sprintf("Simulation %d/%d", sim, n_sims))

      # -----------------------------------------------------------------------
      # Simulate data from DGM
      # -----------------------------------------------------------------------
      dfs <- simulate_from_dgm(
        dgm = dgm,
        n = n_sample,
        rand_ratio = 1,
        draw_treatment = TRUE,
        entry_var = if ("entrytime" %in% names(dgm$df_super)) "entrytime" else NULL,
        analysis_time = analysis_time,
        cens_adjust = cens_adjust,
        seed = if (!is.null(seed)) seed + sim else NULL
      )

      dfs <- data.table::as.data.table(dfs)

      # -----------------------------------------------------------------------
      # Split by region
      # -----------------------------------------------------------------------
      df_nonRegA <- dfs[get(region_var) == 0]
      df_regA <- dfs[get(region_var) == 1]

      n_test <- nrow(df_regA)
      n_train <- nrow(df_nonRegA)
      n_itt <- nrow(dfs)

      # -----------------------------------------------------------------------
      # Safe Cox model fitting
      # -----------------------------------------------------------------------
      safe_coxph <- function(formula, data) {
        fit <- tryCatch(
          summary(survival::coxph(formula, data = data))$conf.int,
          error = function(e) NULL
        )
        if (!is.null(fit) && is.numeric(fit[1])) return(fit[1])
        return(NA_real_)
      }

      # -----------------------------------------------------------------------
      # Estimate HRs
      # -----------------------------------------------------------------------
      hr_train <- safe_coxph(
        survival::Surv(y_sim, event_sim) ~ treat_sim,
        df_nonRegA
      )

      if ("strat" %in% names(dfs)) {
        hr_itt <- safe_coxph(
          survival::Surv(y_sim, event_sim) ~ treat_sim + strata(strat),
          dfs
        )
      } else {
        hr_itt <- safe_coxph(
          survival::Surv(y_sim, event_sim) ~ treat_sim,
          dfs
        )
      }

      hr_ittX <- safe_coxph(
        stats::as.formula(
          paste0("survival::Surv(y_sim, event_sim) ~ treat_sim + strata(", region_var, ")")
        ),
        dfs
      )

      # -----------------------------------------------------------------------
      # Initialize results
      # -----------------------------------------------------------------------
      hr_test <- NA_real_
      any_found <- 0
      sg_found <- "none"
      n_sg <- n_test
      hr_sg <- NA_real_
      prev_sg <- 1.0
      POhr_sg <- NA_real_

      # Testing region HR
      hr_test <- safe_coxph(
        survival::Surv(y_sim, event_sim) ~ treat_sim,
        df_regA
      )

      # -----------------------------------------------------------------------
      # Run ForestSearch if testing HR is valid
      # -----------------------------------------------------------------------
      if (!is.na(hr_test)) {

        if ("loghr_po" %in% names(df_regA)) {
          POhr_sg <- exp(mean(df_regA$loghr_po, na.rm = TRUE))
        }

        fs_result <- tryCatch({
          forestsearch(
            df.analysis = as.data.frame(df_nonRegA),
            df.test = as.data.frame(df_regA),
            confounders.name = confounders.name,
            outcome.name = "y_sim",
            treat.name = "treat_sim",
            event.name = "event_sim",
            id.name = if ("id" %in% names(df_nonRegA)) "id" else NULL,
            potentialOutcome.name = if ("loghr_po" %in% names(df_nonRegA)) "loghr_po" else NULL,
            hr.threshold = hr.threshold,
            hr.consistency = hr.consistency,
            pconsistency.threshold = pconsistency.threshold,
            sg_focus = sg_focus,
            conf_force = conf_force,
            maxk = maxk,
            n.min = 60,
            d0.min = 10,
            d1.min = 10,
            showten_subgroups = FALSE,
            details = FALSE,
            plot.sg = FALSE,
            parallel_args = list(plan = "sequential")
          )
        }, error = function(e) NULL)

        # ---------------------------------------------------------------------
        # Process ForestSearch results
        # ---------------------------------------------------------------------
        if (!is.null(fs_result)) {
          if (is.null(fs_result$sg.harm)) {
            any_found <- 0
            sg_found <- "none"
            n_sg <- n_test
            hr_sg <- hr_test
            prev_sg <- 1.0
          } else {
            any_found <- 1
            sg_found <- paste(fs_result$sg.harm, collapse = " & ")

            df_test <- fs_result$df.test
            if (!is.null(df_test) && "treat.recommend" %in% names(df_test)) {
              df_sg <- df_test[df_test$treat.recommend == 1, ]
              n_sg <- nrow(df_sg)
              prev_sg <- n_sg / n_test
              hr_sg <- safe_coxph(
                survival::Surv(y_sim, event_sim) ~ treat_sim,
                df_sg
              )

              if ("loghr_po" %in% names(df_sg)) {
                POhr_sg <- exp(mean(df_sg$loghr_po, na.rm = TRUE))
              }
            }
          }
        }
      }

      # -----------------------------------------------------------------------
      # Return results
      # -----------------------------------------------------------------------
      data.table::data.table(
        sim = sim,
        n_itt = n_itt,
        hr_itt = hr_itt,
        hr_ittX = hr_ittX,
        n_train = n_train,
        hr_train = hr_train,
        n_test = n_test,
        hr_test = hr_test,
        any_found = any_found,
        sg_found = sg_found,
        n_sg = n_sg,
        hr_sg = hr_sg,
        POhr_sg = POhr_sg,
        prev_sg = prev_sg
      )
    }
  })

  # ===========================================================================
  # SECTION 5: POST-PROCESSING
  # ===========================================================================

  future::plan(future::sequential)

  results <- data.table::as.data.table(results)
  results[, hr_sg_null := data.table::fifelse(any_found == 0, NA_real_, hr_sg)]

  # ===========================================================================
  # SECTION 6: REPORT SUMMARY
  # ===========================================================================

  t_now <- proc.time()[3]
  t_min <- (t_now - t_start) / 60

  if (details) {
    message(sprintf("\n=== Simulation Complete ==="))
    message(sprintf("Total simulations: %d", n_sims))
    message(sprintf("Time elapsed: %.2f minutes", t_min))
    message(sprintf("Projection per 1000 sims: %.2f minutes", t_min * (1000 / n_sims)))
    message(sprintf("Proportion subgroups found: %.3f", mean(results$any_found, na.rm = TRUE)))
    message(sprintf("Mean HR (ITT): %.3f", mean(results$hr_itt, na.rm = TRUE)))
    message(sprintf("Mean HR (Test/Region A): %.3f", mean(results$hr_test, na.rm = TRUE)))
    message(sprintf("Mean HR (Subgroup): %.3f", mean(results$hr_sg_null, na.rm = TRUE)))
  }

  return(results)
}


# =============================================================================
# Data Validation Function
# =============================================================================

#' Validate Dataset for MRCT Simulations
#'
#' Checks that a dataset contains all required variables for MRCT simulation
#' functions and reports any issues. Required variables include outcome (tte, event),
#' treatment (treat), continuous covariates (age, bm), and factor covariates
#' (male, histology, prior_treat, regA).
#'
#' @param df.case Data frame to validate
#' @param verbose Logical. Print detailed validation results. Default: TRUE
#'
#' @return Logical. TRUE if all requirements met, FALSE otherwise (invisibly)
#'
#' @details
#' ## Required Variables
#'
#' The function checks for the following variables:
#' \itemize{
#'   \item \strong{Outcome}: tte (time-to-event), event (0/1 indicator)
#'   \item \strong{Treatment}: treat (0/1 indicator)
#'   \item \strong{Continuous}: age, bm (biomarker)
#'   \item \strong{Factor}: male (0/1), histology, prior_treat (0/1), regA (0/1)
#' }
#'
#' The function also validates variable types and value ranges.
#'
#' @examples
#' \dontrun{
#' # Check if dataset is ready for MRCT simulations
#' is_valid <- validate_mrct_data(df.case)
#'
#' if (!is_valid) {
#'   stop("Please fix data issues before running simulations")
#' }
#' }
#'
#' @seealso \code{\link{create_dgm_for_mrct}} for creating DGM from validated data
#'
#' @export

validate_mrct_data <- function(df.case, verbose = TRUE) {

  # Define required variables by category
  required <- list(
    outcome = c("tte", "event"),
    treatment = "treat",
    continuous = c("age", "bm"),
    factor = c("male", "histology", "prior_treat", "regA")
  )

  all_required <- unlist(required)

  # Check for missing variables
  present <- names(df.case)
  missing <- setdiff(all_required, present)
  found <- intersect(all_required, present)

  # Initialize results
  is_valid <- TRUE
  issues <- list()

  # ---------------------------------------------------------------------------
  # Check presence of variables
  # ---------------------------------------------------------------------------
  if (length(missing) > 0) {
    is_valid <- FALSE
    issues$missing <- missing
  }

  # ---------------------------------------------------------------------------
  # Check variable types and values for found variables
  # ---------------------------------------------------------------------------
  type_issues <- character()
  value_issues <- character()

  if ("tte" %in% found) {
    if (!is.numeric(df.case$tte)) {
      type_issues <- c(type_issues, "tte should be numeric")
    } else if (any(df.case$tte < 0, na.rm = TRUE)) {
      value_issues <- c(value_issues, "tte contains negative values")
    }
  }

  if ("event" %in% found) {
    if (!all(df.case$event %in% c(0, 1, NA))) {
      value_issues <- c(value_issues, "event should be binary (0/1)")
    }
  }

  if ("treat" %in% found) {
    if (!all(df.case$treat %in% c(0, 1, NA))) {
      value_issues <- c(value_issues, "treat should be binary (0/1)")
    }
  }

  if ("age" %in% found) {
    if (!is.numeric(df.case$age)) {
      type_issues <- c(type_issues, "age should be numeric")
    }
  }

  if ("bm" %in% found) {
    if (!is.numeric(df.case$bm)) {
      type_issues <- c(type_issues, "bm (biomarker) should be numeric")
    }
  }

  if ("male" %in% found) {
    if (!all(df.case$male %in% c(0, 1, NA))) {
      value_issues <- c(value_issues, "male should be binary (0/1)")
    }
  }

  if ("regA" %in% found) {
    if (!all(df.case$regA %in% c(0, 1, NA))) {
      value_issues <- c(value_issues, "regA should be binary (0/1)")
    }
  }

  if (length(type_issues) > 0) {
    is_valid <- FALSE
    issues$type_issues <- type_issues
  }

  if (length(value_issues) > 0) {
    is_valid <- FALSE
    issues$value_issues <- value_issues
  }

  # ---------------------------------------------------------------------------
  # Report results
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("=== MRCT Data Validation ===\n\n")
    cat("Dataset:", deparse(substitute(df.case)), "\n")
    cat("Dimensions:", nrow(df.case), "rows x", ncol(df.case), "columns\n\n")

    # Report by category
    for (cat_name in names(required)) {
      cat_vars <- required[[cat_name]]
      cat_found <- intersect(cat_vars, found)
      cat_missing <- setdiff(cat_vars, found)

      cat(sprintf("%-12s: ", tools::toTitleCase(cat_name)))

      if (length(cat_missing) == 0) {
        cat("All present (", paste(cat_found, collapse = ", "), ")\n", sep = "")
      } else {
        cat("Missing: ", paste(cat_missing, collapse = ", "), "\n", sep = "")
      }
    }

    cat("\n")

    # Report issues
    if (length(type_issues) > 0) {
      cat("Type issues:\n")
      for (issue in type_issues) {
        cat("  - ", issue, "\n", sep = "")
      }
    }

    if (length(value_issues) > 0) {
      cat("Value issues:\n")
      for (issue in value_issues) {
        cat("  - ", issue, "\n", sep = "")
      }
    }

    # Summary
    cat("\n")
    if (is_valid) {
      cat("Dataset is valid for MRCT simulations\n")
    } else {
      cat("Dataset requires fixes before use\n")
    }

    # Additional info
    if (is_valid && "regA" %in% found) {
      cat("\nRegion distribution:\n")
      cat("  Region A (test):     ", sum(df.case$regA == 1, na.rm = TRUE),
          " (", round(100 * mean(df.case$regA == 1, na.rm = TRUE), 1), "%)\n", sep = "")
      cat("  Non-Region A (train): ", sum(df.case$regA == 0, na.rm = TRUE),
          " (", round(100 * mean(df.case$regA == 0, na.rm = TRUE), 1), "%)\n", sep = "")
    }
  }

  invisible(is_valid)
}


# =============================================================================
# DGM Creation Wrapper for MRCT
# =============================================================================

#' Create Data Generating Mechanism for MRCT Simulations
#'
#' Wrapper function to create a data generating mechanism (DGM) for MRCT
#' simulation scenarios using \code{\link{generate_aft_dgm_flex}}.
#'
#' @param df_case Data frame containing case study data
#' @param model_type Character. Either "alt" (alternative hypothesis with
#'   heterogeneous treatment effects) or "null" (uniform treatment effect)
#' @param log_hrs Numeric vector. Log hazard ratios for spline specification.
#'   If NULL, defaults are used based on model_type
#' @param confounder_var Character. Name of a confounder variable to include with
#'   a forced prognostic effect. Default: NULL (no forced effect)
#' @param confounder_effect Numeric. Log hazard ratio for confounder_var effect.
#'   Only used if confounder_var is specified
#' @param include_regA Logical. Include regA as a factor in the model.
#'   Default: TRUE
#' @param verbose Logical. Print detailed output. Default: FALSE
#'
#' @return An object of class "aft_dgm_flex" for use with
#'   \code{\link{simulate_from_dgm}} and \code{\link{mrct_region_sims}}
#'
#' @details
#' ## Model Types
#'
#' \describe{
#'   \item{alt}{Alternative hypothesis: Treatment effect varies by biomarker
#'     level (heterogeneous treatment effect). Default log_hrs create HR
#'     ranging from 2.0 (harm) to 0.5 (benefit) across biomarker range}
#'   \item{null}{Null hypothesis: Uniform treatment effect regardless of
#'     biomarker level. Default log_hrs = log(0.7) uniformly}
#' }
#'
#' ## Confounder Effects
#'
#' By default, NO prognostic confounder effect is forced. The confounder_var
#' and confounder_effect parameters allow optionally specifying ANY baseline
#' covariate to have a fixed prognostic effect in the outcome model.
#'
#' The regA variable (region indicator) is included as a factor by default
#' but without a forced effect - its coefficient is estimated from data.
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' df_case <- read.csv("data/dfsynthetic.csv")
#' df_case$regA <- df_case$region_asia
#'
#' # Alternative hypothesis (heterogeneous treatment effect)
#' dgm_alt <- create_dgm_for_mrct(
#'   df_case = df_case,
#'   model_type = "alt",
#'   log_hrs = log(c(3, 1.25, 0.50)),
#'   verbose = TRUE
#' )
#'
#' # Null hypothesis (uniform effect)
#' dgm_null <- create_dgm_for_mrct(
#'   df_case = df_case,
#'   model_type = "null",
#'   verbose = TRUE
#' )
#'
#' # With forced confounder effect
#' dgm_conf <- create_dgm_for_mrct(
#'   df_case = df_case,
#'   model_type = "alt",
#'   confounder_var = "prior_treat",
#'   confounder_effect = log(1.5),
#'   verbose = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{generate_aft_dgm_flex}} for underlying DGM creation
#' \code{\link{mrct_region_sims}} for running simulations with the DGM
#'
#' @export

create_dgm_for_mrct <- function(
    df_case,
    model_type = c("alt", "null"),
    log_hrs = NULL,
    confounder_var = NULL,
    confounder_effect = NULL,
    include_regA = TRUE,
    verbose = FALSE
) {

  model_type <- match.arg(model_type)

  # ---------------------------------------------------------------------------
  # Set default log hazard ratios based on model type
  # ---------------------------------------------------------------------------
  if (is.null(log_hrs)) {
    if (model_type == "alt") {
      # HTE: HR varies from 2 (harmful) to 0.5 (protective) across biomarker
      log_hrs <- log(c(2, 1.25, 0.5))
    } else {
      # Null: uniform HR = 0.7
      log_hrs <- log(c(0.7, 0.7, 0.7))
    }
  }

  # ---------------------------------------------------------------------------
  # Build factor_vars list
  # ---------------------------------------------------------------------------
  factor_vars <- c("male", "histology", "prior_treat")
  if (include_regA) {
    factor_vars <- c(factor_vars, "regA")
  }

  # ---------------------------------------------------------------------------
  # Build set_beta_spec for optional confounder effect
  # ---------------------------------------------------------------------------
  set_beta_spec <- list(set_var = NULL, beta_var = NULL)

  if (!is.null(confounder_var) && !is.null(confounder_effect)) {
    # Add z_ prefix if not already present
    z_confounder_var <- if (grepl("^z_", confounder_var)) {
      confounder_var
    } else {
      paste0("z_", confounder_var)
    }

    set_beta_spec <- list(
      set_var = z_confounder_var,
      beta_var = confounder_effect
    )

    if (verbose) {
      message(sprintf("Forcing confounder effect: %s = %.4f (HR = %.3f)",
                      z_confounder_var, confounder_effect, exp(confounder_effect)))
    }
  }

  # ---------------------------------------------------------------------------
  # Generate DGM
  # ---------------------------------------------------------------------------
  dgm <- generate_aft_dgm_flex(
    data = df_case,
    continuous_vars = c("age", "bm"),
    factor_vars = factor_vars,
    set_beta_spec = set_beta_spec,
    continuous_vars_cens = c("age"),
    factor_vars_cens = c("prior_treat"),
    cens_type = "weibull",
    outcome_var = "tte",
    event_var = "event",
    treatment_var = "treat",
    subgroup_vars = NULL,
    subgroup_cuts = NULL,
    model = model_type,
    spline_spec = list(
      var = "z_bm",
      knot = 5,
      zeta = 10,
      log_hrs = log_hrs
    ),
    k_inter = 0.0,
    verbose = verbose,
    standardize = FALSE
  )

  return(dgm)
}


# =============================================================================
# Summary Output Function
# =============================================================================

#' Summary Tables for MRCT Simulation Results
#'
#' Creates summary tables from MRCT simulation results using the gt package.
#' Summarizes hazard ratio estimates, subgroup identification rates, and
#' classification of identified subgroups.
#'
#' @param pop_summary List. Population summary from large sample approximation
#'   (optional). Default: NULL
#' @param mrct_sims data.table. Simulation results from \code{\link{mrct_region_sims}}
#' @param sg_type Integer. Type of subgroup summary:
#'   \itemize{
#'     \item 1: Basic summary (found, biomarker, age)
#'     \item 2: Extended summary (all subgroup types)
#'   }
#'   Default: 1
#' @param tab_caption Character. Caption for the output table.
#'   Default: "Identified subgroups and estimation summaries"
#' @param digits Integer. Number of decimal places for numeric summaries.
#'   Default: 3
#' @param showtable Logical. Print the table. Default: TRUE
#'
#' @return List with components:
#' \describe{
#'   \item{res}{List of summary statistics from population (if provided)}
#'   \item{out_table}{Formatted gt table object (or data.frame if gt unavailable)}
#'   \item{data}{Processed mrct_sims data.table with derived variables}
#'   \item{summary_df}{Data frame of computed summary statistics}
#' }
#'
#' @examples
#' \dontrun{
#' # After running simulations
#' summary_results <- summaryout_mrct(
#'   pop_summary = NULL,
#'   mrct_sims = results_alt,
#'   sg_type = 1,
#'   tab_caption = "Subgroups under strong biomarker effect"
#' )
#' }
#'
#' @seealso \code{\link{mrct_region_sims}} for generating simulation results
#'
#' @importFrom data.table as.data.table copy fifelse `:=`
#' @importFrom stats median sd quantile
#' @importFrom gt gt tab_header tab_style cell_text cell_borders cells_body
#'   tab_source_note cols_label cols_width tab_options px
#' @export

summaryout_mrct <- function(
    pop_summary = NULL,
    mrct_sims,
    sg_type = 1,
    tab_caption = "Identified subgroups and estimation summaries",
    digits = 3,
    showtable = TRUE
) {

 # Ensure data.table
  mrct_sims <- data.table::as.data.table(data.table::copy(mrct_sims))

  res <- list()

  # ---------------------------------------------------------------------------
  # Process population summary if provided
  # ---------------------------------------------------------------------------
  if (!is.null(pop_summary)) {
    res$ahr_true <- round(null_or(pop_summary$AHR, NA), 4)
    res$ahr_unadj <- round(null_or(pop_summary$ITT_unadj, NA), 4)
    res$ahr_sR <- round(null_or(pop_summary$ITT_sR, NA), 4)
    res$ahr_sRw <- round(null_or(pop_summary$ITT_sRw, NA), 4)

    if (!is.na(res$ahr_true) && res$ahr_true != 0) {
      res$bias_unadj <- round(100 * (res$ahr_unadj - res$ahr_true) / res$ahr_true, 1)
      res$bias_sR <- round(100 * (res$ahr_sR - res$ahr_true) / res$ahr_true, 1)
      res$bias_sRw <- round(100 * (res$ahr_sRw - res$ahr_true) / res$ahr_true, 1)
    }

    res$ahr_w1 <- round(null_or(pop_summary$W_1, NA), 4)
    if (!is.na(res$ahr_true) && res$ahr_true != 0 && !is.na(res$ahr_w1)) {
      res$bias_w1 <- round(100 * (res$ahr_w1 - res$ahr_true) / res$ahr_true, 1)
    }
  }

  # ---------------------------------------------------------------------------
  # Create derived variables
  # ---------------------------------------------------------------------------
  mrct_sims[, `:=`(
    regAflag = data.table::fifelse(hr_test > 0.9, "RegA > 0.9", "RegA <= 0.9"),
    sg_le85 = data.table::fifelse(hr_sg <= 0.85, "RegA(sg) <= 0.85", "RegA(sg) > 0.85"),
    regAflag2 = data.table::fifelse(
      hr_test > 0.9 & hr_sg <= 0.85,
      "RegA > 0.9 & RegA(sg) <= 0.85", "Not"
    ),
    regAflag3 = data.table::fifelse(
      hr_test > 0.9 & hr_sg <= 0.80,
      "RegA > 0.9 & RegA(sg) <= 0.80", "Not"
    ),
    found = as.factor(any_found)
  )]

  # ---------------------------------------------------------------------------
  # Classify subgroups by factor type
  # ---------------------------------------------------------------------------
  mrct_sims[, `:=`(
    sg_biomarker = data.table::fifelse(
      grepl("bm|biomarker", sg_found, ignore.case = TRUE),
      "biomarker",
      data.table::fifelse(sg_found != "none", "other", "none")
    ),
    sg_age = data.table::fifelse(
      grepl("age", sg_found, ignore.case = TRUE),
      "age",
      data.table::fifelse(sg_found != "none", "other", "none")
    ),
    sg_male = data.table::fifelse(
      grepl("male|sex|gender", sg_found, ignore.case = TRUE),
      "sex",
      data.table::fifelse(sg_found != "none", "other", "none")
    ),
    sg_ecog = data.table::fifelse(
      grepl("ecog", sg_found, ignore.case = TRUE),
      "ecog",
      data.table::fifelse(sg_found != "none", "other", "none")
    ),
    sg_histology = data.table::fifelse(
      grepl("histology", sg_found, ignore.case = TRUE),
      "histology",
      data.table::fifelse(sg_found != "none", "other", "none")
    ),
    sg_CTregimen = data.table::fifelse(
      grepl("CTregimen|chemo", sg_found, ignore.case = TRUE),
      "CT_regimen",
      data.table::fifelse(sg_found != "none", "other", "none")
    ),
    sg_region = data.table::fifelse(
      grepl("region|EU", sg_found, ignore.case = TRUE),
      "region",
      data.table::fifelse(sg_found != "none", "other", "none")
    ),
    sg_surgery = data.table::fifelse(
      grepl("surgery", sg_found, ignore.case = TRUE),
      "surgery",
      data.table::fifelse(sg_found != "none", "other", "none")
    ),
    sg_prior_treat = data.table::fifelse(
      grepl("prior_treat", sg_found, ignore.case = TRUE),
      "prior_treat",
      data.table::fifelse(sg_found != "none", "other", "none")
    )
  )]

  # ---------------------------------------------------------------------------
  # Helper function to format mean (SD)
  # ---------------------------------------------------------------------------
  fmt_mean_sd <- function(x, d = digits) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return("--")
    sprintf("%.*f (%.*f)", d, mean(x), d, stats::sd(x))
  }

  # Helper function to format median [IQR]
  fmt_median_iqr <- function(x, d = digits) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return("--")
    q <- stats::quantile(x, c(0.25, 0.5, 0.75))
    sprintf("%.*f [%.*f, %.*f]", d, q[2], d, q[1], d, q[3])
  }

  # Helper function for proportion
  fmt_prop <- function(x) {
    sprintf("%.1f%%", 100 * mean(x, na.rm = TRUE))
  }

  # Helper for category counts
  fmt_cat <- function(x, level) {
    n <- sum(x == level, na.rm = TRUE)
    pct <- 100 * n / length(x)
    sprintf("%d (%.1f%%)", n, pct)
  }

  # ---------------------------------------------------------------------------
  # Build summary data frame
  # ---------------------------------------------------------------------------
  n_sims <- nrow(mrct_sims)

  summary_rows <- list(
    # HR Estimates
    data.frame(Category = "HR Estimates", Metric = "HR ITT",
               Value = fmt_mean_sd(mrct_sims$hr_itt), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "HR ITT (stratified by region)",
               Value = fmt_mean_sd(mrct_sims$hr_ittX), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "HR Training (non-Region A)",
               Value = fmt_mean_sd(mrct_sims$hr_train), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "HR Testing (Region A)",
               Value = fmt_mean_sd(mrct_sims$hr_test), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "HR Subgroup",
               Value = fmt_mean_sd(mrct_sims$hr_sg), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "PO HR Subgroup",
               Value = fmt_mean_sd(mrct_sims$POhr_sg), stringsAsFactors = FALSE),

    # Subgroup Identification
    data.frame(Category = "Subgroup Identification", Metric = "Subgroup Found",
               Value = fmt_prop(mrct_sims$any_found), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "Prevalence (when found)",
               Value = fmt_mean_sd(mrct_sims$prev_sg[mrct_sims$any_found == 1]),
               stringsAsFactors = FALSE),

    # Region A Flags
    data.frame(Category = "Region A Classification", Metric = "RegA HR > 0.9",
               Value = fmt_cat(mrct_sims$regAflag, "RegA > 0.9"), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "RegA(sg) HR <= 0.85",
               Value = fmt_cat(mrct_sims$sg_le85, "RegA(sg) <= 0.85"), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "RegA > 0.9 & RegA(sg) <= 0.85",
               Value = fmt_cat(mrct_sims$regAflag2, "RegA > 0.9 & RegA(sg) <= 0.85"),
               stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "RegA > 0.9 & RegA(sg) <= 0.80",
               Value = fmt_cat(mrct_sims$regAflag3, "RegA > 0.9 & RegA(sg) <= 0.80"),
               stringsAsFactors = FALSE)
  )

  # Subgroup classification rows
  sg_class_rows <- list(
    data.frame(Category = "Subgroup Classification", Metric = "Biomarker",
               Value = fmt_cat(mrct_sims$sg_biomarker, "biomarker"), stringsAsFactors = FALSE),
    data.frame(Category = "", Metric = "Age",
               Value = fmt_cat(mrct_sims$sg_age, "age"), stringsAsFactors = FALSE)
  )

  # Extended classification for sg_type == 2
  if (sg_type == 2) {
    sg_class_rows <- c(sg_class_rows, list(
      data.frame(Category = "", Metric = "Sex",
                 Value = fmt_cat(mrct_sims$sg_male, "sex"), stringsAsFactors = FALSE),
      data.frame(Category = "", Metric = "ECOG",
                 Value = fmt_cat(mrct_sims$sg_ecog, "ecog"), stringsAsFactors = FALSE),
      data.frame(Category = "", Metric = "Histology",
                 Value = fmt_cat(mrct_sims$sg_histology, "histology"), stringsAsFactors = FALSE),
      data.frame(Category = "", Metric = "CT Regimen",
                 Value = fmt_cat(mrct_sims$sg_CTregimen, "CT_regimen"), stringsAsFactors = FALSE),
      data.frame(Category = "", Metric = "Region",
                 Value = fmt_cat(mrct_sims$sg_region, "region"), stringsAsFactors = FALSE),
      data.frame(Category = "", Metric = "Surgery",
                 Value = fmt_cat(mrct_sims$sg_surgery, "surgery"), stringsAsFactors = FALSE),
      data.frame(Category = "", Metric = "Prior Treatment",
                 Value = fmt_cat(mrct_sims$sg_prior_treat, "prior_treat"), stringsAsFactors = FALSE)
    ))
  }

  summary_rows <- c(summary_rows, sg_class_rows)
  summary_df <- do.call(rbind, summary_rows)

  # ---------------------------------------------------------------------------
  # Create gt table
  # ---------------------------------------------------------------------------
  if (!requireNamespace("gt", quietly = TRUE)) {
    warning("Package 'gt' not available. Returning data frame instead of formatted table.")
    if (showtable) print(summary_df)
    return(list(res = res, out_table = summary_df, data = mrct_sims, summary_df = summary_df))
  }

  # Build gt table
  out_table <- gt::gt(summary_df) |>
    gt::tab_header(
      title = tab_caption,
      subtitle = sprintf("Based on %d simulations", n_sims)
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        columns = "Category",
        rows = Category != ""
      )
    ) |>
    gt::tab_style(
      style = gt::cell_borders(
        sides = "bottom",
        color = "#D3D3D3",
        weight = gt::px(1)
      ),
      locations = gt::cells_body()
    ) |>
    gt::cols_label(
      Category = "Category",
      Metric = "Metric",
      Value = "Mean (SD) / N (%)"
    ) |>
    gt::cols_width(
      Category ~ gt::px(180),
      Metric ~ gt::px(250),
      Value ~ gt::px(180)
    ) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(14),
      heading.subtitle.font.size = gt::px(12),
      column_labels.font.weight = "bold"
    )

  # Add population summary if available
  if (!is.null(pop_summary) && length(res) > 0) {
    pop_note <- sprintf("Population AHR: %.4f", null_or(res$ahr_true, NA))
    out_table <- gt::tab_source_note(out_table, source_note = pop_note)
  }

  if (showtable) print(out_table)

  return(list(res = res, out_table = out_table, data = mrct_sims, summary_df = summary_df))
}


# =============================================================================
# Visualization Function
# =============================================================================

#' Violin/Boxplot Visualization of HR Estimates
#'
#' Creates violin plots with embedded boxplots showing the distribution of
#' hazard ratio estimates across simulations for different analysis populations.
#'
#' @param df data.frame or data.table. Simulation results from
#'   \code{\link{mrct_region_sims}}
#' @param label_training Character. Label for training data estimates.
#'   Default: "Training"
#' @param label_testing Character. Label for testing data estimates.
#'   Default: "Testing"
#' @param label_itt Character. Label for ITT estimates.
#'   Default: "ITT (stratified)"
#' @param label_sg Character. Label for subgroup estimates.
#'   Default: "Testing (subgroup)"
#'
#' @return List with components:
#' \describe{
#'   \item{dfPlot_estimates}{data.table formatted for plotting}
#'   \item{plot_estimates}{ggplot2 object}
#' }
#'
#' @examples
#' \dontrun{
#' # After running simulations
#' plot_results <- SGplot_estimates(
#'   results_alt,
#'   label_training = "Non-Region A, ITT",
#'   label_itt = "Overall, ITT",
#'   label_testing = "Region A, ITT",
#'   label_sg = "Region A, identified subgroup"
#' )
#' print(plot_results$plot_estimates)
#' }
#'
#' @seealso \code{\link{mrct_region_sims}} for generating simulation results
#'
#' @importFrom data.table data.table
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_hline
#'   scale_fill_brewer labs theme_minimal theme element_text
#' @export

SGplot_estimates <- function(
    df,
    label_training = "Training",
    label_testing = "Testing",
    label_itt = "ITT (stratified)",
    label_sg = "Testing (subgroup)"
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting. Please install it.")
  }

  # ---------------------------------------------------------------------------
  # Prepare data for plotting
  # ---------------------------------------------------------------------------
  df_itt <- data.table::data.table(est = df$hr_itt, analysis = label_itt)
  df_training <- data.table::data.table(est = df$hr_train, analysis = label_training)
  df_testing <- data.table::data.table(est = df$hr_test, analysis = label_testing)
  df_sg <- data.table::data.table(est = df$hr_sg, analysis = label_sg)

  hr_estimates <- rbind(df_itt, df_training, df_testing, df_sg)

  # Set factor order
  est_order <- c(label_itt, label_training, label_testing, label_sg)
  hr_estimates[, analysis := factor(analysis, levels = est_order)]

  # ---------------------------------------------------------------------------
  # Create plot
  # ---------------------------------------------------------------------------
  p <- ggplot2::ggplot(hr_estimates, ggplot2::aes(x = analysis, y = est, fill = analysis)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", alpha = 0.7) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      x = "Analysis Population",
      y = "Hazard Ratio Estimate",
      title = "Distribution of HR Estimates Across Simulations"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 15, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  return(list(dfPlot_estimates = hr_estimates, plot_estimates = p))
}


# =============================================================================
# Helper Function: Null Coalescing Operator
# =============================================================================

# Null Coalescing Operator (internal helper, no documentation generated)
# Returns the left-hand side if not NULL, otherwise returns the right-hand side.
# @noRd
null_or <- function(x, y) {
  if (is.null(x)) y else x
}
