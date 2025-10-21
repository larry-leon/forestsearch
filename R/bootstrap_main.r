#' ForestSearch Bootstrap with doFuture Parallelization
#'
#' Orchestrates bootstrap analysis for ForestSearch using doFuture parallelization.
#' Implements bias correction methods to adjust for optimism in subgroup selection.
#'
#' @section Bias Correction Methods:
#' Two bias correction approaches are implemented:
#' \enumerate{
#'   \item \strong{Method 1 (Simple Optimism)}:
#'     \deqn{H_{adj1} = H_{obs} - (H^{*}_{*} - H^{*}_{obs})}
#'     where H^{*}_{*} is the new subgroup HR on bootstrap data and
#'     H^{*}_{obs} is the new subgroup HR on original data.
#'
#'   \item \strong{Method 2 (Double Bootstrap)}:
#'     \deqn{H_{adj2} = 2 \times H_{obs} - (H_{*} + H^{*}_{*} - H^{*}_{obs})}
#'     where H_{*} is the original subgroup HR on bootstrap data.
#' }
#'
#' @section Variable Naming Convention:
#' - \code{H}: Original subgroup (harm/questionable, treat.recommend == 0)
#' - \code{Hc}: Complement subgroup (recommend, treat.recommend == 1)
#' - \code{_obs}: Estimate from original data
#' - \code{_star}: Estimate from bootstrap data
#' - \code{_biasadj_1}: Bias correction method 1
#' - \code{_biasadj_2}: Bias correction method 2
#'
#' @param fs.est List. ForestSearch results object from \code{\link{forestsearch}}.
#'   Must contain \code{df.est} (data frame) and \code{args_call_all} (list of arguments).
#' @param nb_boots Integer. Number of bootstrap samples (recommend 500-1000).
#' @param details Logical. If \code{TRUE}, prints detailed progress information.
#'   Default: \code{FALSE}.
#' @param show_three Logical. If \code{TRUE}, shows verbose output for first 3
#'   bootstrap iterations for debugging. Default: \code{FALSE}.
#' @param parallel_args List. Parallelization configuration with elements:
#'   \itemize{
#'     \item \code{plan}: Character. One of "multisession", "multicore", "callr",
#'       or "sequential"
#'     \item \code{workers}: Integer. Number of parallel workers
#'     \item \code{show_message}: Logical. Show parallel setup messages
#'   }
#'   If empty list, inherits settings from original forestsearch call.
#' @param create_summary Logical. If \code{TRUE}, generates enhanced summary
#'   output with formatted tables and diagnostics. Default: \code{TRUE}.
#' @param create_plots Logical. If \code{TRUE}, generates diagnostic plots
#'   of bootstrap distributions. Requires \code{ggplot2}. Default: \code{FALSE}.
#'
#' @return List with the following components:
#' \describe{
#'   \item{results}{Data.table with bias-corrected estimates for each bootstrap iteration}
#'   \item{SG_CIs}{List of confidence intervals for H and Hc (raw and bias-corrected)}
#'   \item{FSsg_tab}{Formatted table of subgroup estimates}
#'   \item{Ystar_mat}{Matrix (nb_boots x n) of bootstrap sample indicators}
#'   \item{H_estimates}{Detailed estimates for subgroup H}
#'   \item{Hc_estimates}{Detailed estimates for subgroup Hc}
#'   \item{summary}{(If create_summary=TRUE) Enhanced summary with tables and diagnostics}
#' }
#'
#' @section Performance:
#' Typical runtime: 1-5 seconds per bootstrap iteration. For 1000 bootstraps with
#' 6 workers, expect 3-10 minutes total. Memory usage scales with dataset size
#' and number of workers.
#'
#' @section Requirements:
#' \itemize{
#'   \item Original \code{fs.est} must have identified a valid subgroup
#'   \item Requires packages: \code{data.table}, \code{foreach}, \code{doFuture},
#'     \code{survival}
#'   \item For plots: requires \code{ggplot2}
#' }
#'
#' @examples
#' \dontrun{
#' # Run ForestSearch
#' fs_result <- forestsearch(
#'   df.analysis = mydata,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   confounders.name = c("age", "sex", "stage")
#' )
#'
#' # Run bootstrap with bias correction
#' boot_results <- forestsearch_bootstrap_dofuture(
#'   fs.est = fs_result,
#'   nb_boots = 1000,
#'   parallel_args = list(
#'     plan = "multisession",
#'     workers = 6,
#'     show_message = TRUE
#'   ),
#'   create_summary = TRUE,
#'   create_plots = TRUE
#' )
#'
#' # View results
#' print(boot_results$FSsg_tab)
#' print(boot_results$summary$table)
#'
#' # Check success rate
#' mean(!is.na(boot_results$results$H_biasadj_2))
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for initial subgroup identification
#' \code{\link{bootstrap_results}} for the core bootstrap worker function
#' \code{\link{build_cox_formula}} for Cox formula construction
#' \code{\link{fit_cox_models}} for Cox model fitting
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
                                            parallel_args = list(),
                                            create_summary = TRUE,
                                            create_plots = FALSE) {

  # =======================================================================
  # SECTION 1: INPUT VALIDATION
  # =======================================================================

  # Validate fs.est structure
  if (!is.list(fs.est)) {
    stop("fs.est must be a list (ForestSearch results object)")
  }

  required_components <- c("df.est", "args_call_all")
  missing <- setdiff(required_components, names(fs.est))
  if (length(missing) > 0) {
    stop("fs.est missing required components: ", paste(missing, collapse = ", "))
  }

  if (is.null(fs.est$df.est) || !is.data.frame(fs.est$df.est) || nrow(fs.est$df.est) == 0) {
    stop("fs.est$df.est must be a non-empty data frame")
  }

  # Validate nb_boots
  if (!is.numeric(nb_boots) || length(nb_boots) != 1 || nb_boots < 1) {
    stop("nb_boots must be a positive integer")
  }

  if (nb_boots < 100) {
    warning("nb_boots < 100 may produce unreliable confidence intervals. ",
            "Recommend at least 500 for production use.")
  }

  # =======================================================================
  # SECTION 2: EXTRACT ARGUMENTS AND RESOLVE PARALLEL CONFIGURATION
  # =======================================================================

  args_forestsearch_call <- fs.est$args_call_all
  parallel_args <- resolve_bootstrap_parallel_args(parallel_args, args_forestsearch_call)

  # =======================================================================
  # SECTION 3: ENSURE REQUIRED PACKAGES
  # =======================================================================

  ensure_packages(BOOTSTRAP_REQUIRED_PACKAGES)

  # =======================================================================
  # SECTION 4: BUILD COX FORMULA
  # =======================================================================

  cox.formula.boot <- do.call(
    build_cox_formula,
    filter_call_args(args_forestsearch_call, build_cox_formula)
  )

  # =======================================================================
  # SECTION 5: FIT COX MODELS ON ORIGINAL DATA
  # =======================================================================

  # Note: Identified subgroups meet minimum size/event requirements
  cox_fits <- fit_cox_models(fs.est$df.est, cox.formula.boot)
  H_obs <- cox_fits$H_obs
  seH_obs <- cox_fits$seH_obs
  Hc_obs <- cox_fits$Hc_obs
  seHc_obs <- cox_fits$seHc_obs

  # =======================================================================
  # SECTION 6: SETUP PARALLEL PROCESSING
  # =======================================================================

  old_plan <- future::plan()
  on.exit({
    future::plan(old_plan)
    # Optional: Force garbage collection after parallel work
    # Uncomment if experiencing memory issues with many workers:
    # gc()
  }, add = TRUE)

  setup_parallel_SGcons(parallel_args)

  # =======================================================================
  # SECTION 7: GENERATE YSTAR MATRIX (BOOTSTRAP INDICATORS)
  # =======================================================================

  Ystar_mat <- bootstrap_ystar(fs.est$df.est, nb_boots)

  if (details) {
    cat("Ystar matrix generated should be 'boots x N': ", nrow(Ystar_mat), " x ",
        ncol(Ystar_mat), "\n", sep = "")
  }

  # Validate Ystar dimensions
  if (nrow(Ystar_mat) != nb_boots || ncol(Ystar_mat) != nrow(fs.est$df.est)) {
    stop("Ystar_mat dimension mismatch. Expected (", nb_boots, " x ",
         nrow(fs.est$df.est), "), got (", nrow(Ystar_mat), " x ",
         ncol(Ystar_mat), ")")
  }

  # =======================================================================
  # SECTION 8: RUN BOOTSTRAP ANALYSIS
  # =======================================================================

  results <- bootstrap_results(
    fs.est = fs.est,
    df_boot_analysis = fs.est$df.est,
    cox.formula.boot = cox.formula.boot,
    nb_boots = nb_boots,
    show_three = show_three,
    H_obs = H_obs,
    Hc_obs = Hc_obs
  )

  # =======================================================================
  # SECTION 9: POST-PROCESSING AND CONFIDENCE INTERVALS
  # =======================================================================

  est.scale <- args_forestsearch_call$est.scale

  # Compute H estimates with error handling
  H_estimates <- try(
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
    silent = TRUE
  )

  # Handle H estimates error
  if (inherits(H_estimates, "try-error")) {
    warning("Failed to compute H estimates: ", as.character(H_estimates))
    H_estimates <- NULL
  }

  # Compute Hc estimates with error handling
  Hc_estimates <- try(
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
    silent = TRUE
  )

  # Handle Hc estimates error
  if (inherits(Hc_estimates, "try-error")) {
    warning("Failed to compute Hc estimates: ", as.character(Hc_estimates))
    Hc_estimates <- NULL
  }

  # Return early if BOTH failed
  if (is.null(H_estimates) && is.null(Hc_estimates)) {
    warning("Both H and Hc estimate calculations failed. ",
            "Returning partial results without confidence intervals.")
    out <- list(
      results = results,
      SG_CIs = NULL,
      FSsg_tab = NULL,
      Ystar_mat = Ystar_mat,
      H_estimates = NULL,
      Hc_estimates = NULL
    )
    return(out)
  }

  # =======================================================================
  # SECTION 10: FORMAT CONFIDENCE INTERVALS
  # =======================================================================

  # Format CIs (handling NULL estimates gracefully)
  H_res1 <- if (!is.null(H_estimates)) {
    format_CI(H_estimates, c("H0", "H0_lower", "H0_upper"))
  } else {
    "NA (NA, NA)"
  }

  H_res2 <- if (!is.null(H_estimates)) {
    format_CI(H_estimates, c("H2", "H2_lower", "H2_upper"))
  } else {
    "NA (NA, NA)"
  }

  Hc_res1 <- if (!is.null(Hc_estimates)) {
    format_CI(Hc_estimates, c("H0", "H0_lower", "H0_upper"))
  } else {
    "NA (NA, NA)"
  }

  Hc_res2 <- if (!is.null(Hc_estimates)) {
    format_CI(Hc_estimates, c("H2", "H2_lower", "H2_upper"))
  } else {
    "NA (NA, NA)"
  }

  # Print details if requested
  if (details) {
    boot_success_rate <- sum(!is.na(results$H_biasadj_2)) / nb_boots
    cat("\n=== Bootstrap Analysis Complete ===\n")
    cat("Success rate: ", sprintf("%.1f%%", boot_success_rate * 100),
        " (", sum(!is.na(results$H_biasadj_2)), "/", nb_boots, ")\n", sep = "")
    cat("\nH (Questionable) Estimates:\n")
    cat("  Unadjusted:      ", H_res1, "\n")
    cat("  Bias-corrected: ", H_res2, "\n")
    cat("\nHc (Recommend) Estimates:\n")
    cat("  Unadjusted:      ", Hc_res1, "\n")
    cat("  Bias-corrected: ", Hc_res2, "\n")
    cat("===================================\n\n")
  }

  SG_CIs <- list(
    H_raw = H_res1,
    H_bc = H_res2,
    Hc_raw = Hc_res1,
    Hc_bc = Hc_res2
  )

  # =======================================================================
  # SECTION 11: CREATE SUMMARY TABLE
  # =======================================================================

  # Determine HR parameters based on scale
  if (est.scale == "1/hr") {
    hr_1a <- SG_CIs$H_bc
    hr_0a <- SG_CIs$Hc_bc
    sg0_name <- "Questionable"
    sg1_name <- "Recommend"
  } else {
    hr_1a <- SG_CIs$Hc_bc
    hr_0a <- SG_CIs$H_bc
    sg0_name <- "Questionable"
    sg1_name <- "Recommend"
  }

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
    est.scale = est.scale,
    sg0_name = sg0_name,
    sg1_name = sg1_name
  )

  # =======================================================================
  # SECTION 12: COMPILE OUTPUT
  # =======================================================================

  out <- list(
    results = results,
    SG_CIs = SG_CIs,
    FSsg_tab = FSsg_tab,
    Ystar_mat = Ystar_mat,
    H_estimates = H_estimates,
    Hc_estimates = Hc_estimates
  )

  # =======================================================================
  # SECTION 13: ENHANCED SUMMARY (OPTIONAL)
  # =======================================================================

  if (create_summary) {
    summary_output <- summarize_bootstrap_results(
      boot_results = out,
      create_plots = create_plots,
      est.scale = est.scale
    )
    out$summary <- summary_output
  }

  # =======================================================================
  # SECTION 14: CLEANUP AND RETURN
  # =======================================================================

  # Force cleanup of parallel resources (handled by on.exit)
  invisible(gc())

  return(out)
}
