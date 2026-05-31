#' ForestSearch Bootstrap with doFuture Parallelization
#'
#' Orchestrates bootstrap analysis for ForestSearch using doFuture parallelization.
#' Implements bias correction methods to adjust for optimism in subgroup selection.
#'
#' @section Bias Correction Methods:
#' Two bias correction approaches are implemented:
#' \enumerate{
#'   \item \strong{Method 1 (Simple Optimism)}:
#'     \deqn{H_{adj1} = H_{obs} - (H^*_{*} - H^*_{obs})}
#'     where \eqn{H^*_{*}} is the new subgroup HR on bootstrap data and
#'     \eqn{H^*_{obs}} is the new subgroup HR on original data.
#'
#'   \item \strong{Method 2 (Double Bootstrap)}:
#'     \deqn{H_{adj2} = 2 \times H_{obs} - (H_{*} + H^*_{*} - H^*_{obs})}
#'     where \eqn{H_{*}} is the original subgroup HR on bootstrap data.
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
#' @param seed Integer. Random seed for reproducibility of bootstrap sample
#'   generation. Default \code{8316951L}. The value is passed to both
#'   \code{\link{bootstrap_ystar}} (which constructs the \eqn{n \times B}
#'   bootstrap index matrix) and \code{\link{bootstrap_results}} (which
#'   re-runs the ForestSearch algorithm on each replicate); both calls must
#'   use the same seed to ensure bootstrap index alignment. Set to
#'   \code{NULL} for a non-reproducible run.
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
#' @param digits Integer. Decimal places for construction-time numeric
#'   formatting in \code{FSsg_tab} (arm summaries, Diff, effect-estimate
#'   CI, bias-corrected CI).  Default 4 to provide round-down headroom
#'   for downstream display reformatting (e.g.,
#'   \code{\link{summarize_bootstrap_results}(digits = ...)}, which can
#'   round down but cannot recover precision lost at construction).
#'   Threaded through \code{\link{format_CI}} and
#'   \code{\link{SG_tab_estimates_glm}}.
#'
#' @return An \code{fs_bootstrap} object (a list with class
#'   \code{c("fs_bootstrap", "list")}) containing:
#' \describe{
#'   \item{results}{Data.table with bias-corrected estimates for each bootstrap iteration}
#'   \item{SG_CIs}{List of confidence intervals for H and Hc (raw and bias-corrected)}
#'   \item{FSsg_tab}{Formatted table of subgroup estimates}
#'   \item{Ystar_mat}{Matrix (nb_boots x n) of bootstrap sample indicators}
#'   \item{H_estimates}{Detailed estimates for subgroup H}
#'   \item{Hc_estimates}{Detailed estimates for subgroup Hc}
#'   \item{summary}{(If create_summary=TRUE) Enhanced summary with tables and diagnostics}
#'   \item{nb_boots}{Integer. Number of bootstrap iterations requested.
#'     Used by \code{\link{print.fs_bootstrap}} to compute identification
#'     percentages without re-inspecting \code{results}.  Added in v0.2.0.}
#'   \item{original_sg}{Character vector. The subgroup identified by the
#'     primary \code{\link{forestsearch}} call (\code{fs.est$sg.harm}),
#'     carried forward so that \code{\link{print.fs_bootstrap}} can
#'     report exact- and partial-match rates without retaining a
#'     reference to \code{fs.est}.  Added in v0.2.0.}
#'   \item{outcome_type}{Character. Outcome type of the underlying
#'     analysis: one of \code{"survival"}, \code{"binary"},
#'     \code{"continuous"}, or \code{"count"}.  Added in v0.2.0.}
#'   \item{effect_measure}{Character. Effect measure label
#'     (\code{"HR"} for survival; \code{"OR"}, \code{"RR"}, \code{"RD"},
#'     \code{"IRR"}, \code{"IRD"}, or \code{"MD"} for GLM outcomes).
#'     Added in v0.2.0.}
#'   \item{est.scale}{Character. Estimation scale for confidence
#'     intervals (\code{"hr"} or \code{"1/hr"}).  Added in v0.2.0.}
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
                                            seed      = 8316951L,
                                            details = FALSE,
                                            show_three = FALSE,
                                            parallel_args = list(),
                                            digits = 4
){

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

  # No identified subgroup => nothing to bias-correct.  Fail loudly here
  # rather than deeper in fit_cox_models()/fit_effect_models(), where the
  # absent 'treat.recommend' column surfaces as a cryptic
  # "object 'treat.recommend' not found".
  if (is.null(fs.est$sg.harm) || length(fs.est$sg.harm) == 0L ||
      !"treat.recommend" %in% names(fs.est$df.est)) {
    stop("forestsearch_bootstrap_dofuture(): no subgroup was identified ",
         "(fs.est$sg.harm is empty, or fs.est$df.est has no 'treat.recommend' ",
         "column), so there is nothing to bias-correct.  The bootstrap ",
         "requires an identified subgroup -- check fs.est$sg.harm before ",
         "calling.", call. = FALSE)
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
  # SECTION 4: BUILD FORMULA / ESTIMATOR
  # =======================================================================

  # Detect GLM outcome type from the original forestsearch call
  outcome_type <- args_forestsearch_call$outcome_type
  if (is.null(outcome_type)) outcome_type <- "survival"
  if (length(outcome_type) > 1L) outcome_type <- outcome_type[1L]
  is_glm <- outcome_type != "survival"

  estimator_fn_boot <- NULL
  cox.formula.boot  <- NULL

  if (is_glm) {
    # GLM path: build the estimator closure from saved arguments
    estimator_fn_boot <- make_effect_estimator(
      outcome_type   = outcome_type,
      treat.name     = args_forestsearch_call$treat.name,
      outcome.name   = args_forestsearch_call$outcome.name,
      event.name     = NULL,
      offset.name    = args_forestsearch_call$offset.name,
      effect_measure = args_forestsearch_call$effect_measure
    )
  } else {
    # Survival path: build Cox formula (unchanged)
    cox.formula.boot <- do.call(
      build_cox_formula,
      filter_call_args(args_forestsearch_call, build_cox_formula)
    )
  }

  # =======================================================================
  # SECTION 5: FIT MODELS ON ORIGINAL DATA
  # =======================================================================

  if (is_glm) {
    # GLM: use estimator closure
    effect_fits <- fit_effect_models(fs.est$df.est, estimator_fn_boot)
  } else {
    # Survival: use Cox (unchanged)
    # Note: Identified subgroups meet minimum size/event requirements
    effect_fits <- fit_cox_models(fs.est$df.est, cox.formula.boot)
  }
  H_obs    <- effect_fits$H_obs
  seH_obs  <- effect_fits$seH_obs
  Hc_obs   <- effect_fits$Hc_obs
  seHc_obs <- effect_fits$seHc_obs

  # =======================================================================
  # SECTION 6: SETUP PARALLEL PROCESSING
  # =======================================================================


  # Suppress warnings during entire parallel section
  old_warn <- getOption("warn")
  options(warn = -1)

  on.exit({
    options(warn = old_warn)  # Restore first

    if (exists(".Last.future.plan")) {
      future::plan(.Last.future.plan)
    } else {
      future::plan("sequential")
    }
    big_objs <- c("Ystar_mat", "boot_index_mat")
    if (any(vapply(big_objs, function(nm) {
      exists(nm, inherits = FALSE) &&
        as.numeric(object.size(get(nm, inherits = FALSE))) > 1e9
    }, logical(1)))) {
      gc(verbose = FALSE, reset = TRUE)
    }
  }, add = TRUE)

  setup_parallel_SGcons(parallel_args)

  # =======================================================================
  # SECTION 7: GENERATE BOOTSTRAP INDEX MATRIX AND DERIVE YSTAR
  # =======================================================================
  # The B x N integer matrix `boot_index_mat` is the single source of
  # truth for bootstrap row indices: it is consumed both by
  # bootstrap_results() (each iteration uses row b for resampling) and by
  # the count matrix Ystar_mat (passed to get_dfRes() for the
  # infinitesimal-jackknife variance estimator).  Generating it once on
  # the main process replaces a prior contract where two independent
  # parallel passes (bootstrap_ystar() and bootstrap_results()) had to
  # produce identical sample.int() outputs by virtue of doFuture's
  # L'Ecuyer-CMRG seed-derivation - an alignment that was easy to break
  # silently.

  boot_seed <- seed
  boot_index_mat <- .make_boot_index_matrix(
    n        = nrow(fs.est$df.est),
    nb_boots = nb_boots,
    seed     = boot_seed
  )
  Ystar_mat <- .boot_indices_to_ystar(boot_index_mat)

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
  # SECTION 7.1: VERIFY BOOTSTRAP FORESTSEARCH PARAMETERS (when details = TRUE)
  # =======================================================================

  if (details) {
    cat("\nForestSearch parameters for bootstrap iterations:\n")
    sm <- args_forestsearch_call$subgroup_method
    if (is.null(sm)) sm <- "consistency"
    cat("  - subgroup_method:", sm, "\n")
    cat("  - sg_focus:", args_forestsearch_call$sg_focus, "\n")
    if (identical(sm, "dina")) {
      # DINA-selection mode bypasses GRF / LASSO / the consistency search;
      # only the selection knobs below are operative.
      cat("  - hr.threshold:", args_forestsearch_call$hr.threshold, "\n")
      cat("  - n.min:", args_forestsearch_call$n.min, "\n")
      cat("  - selection_rule:", args_forestsearch_call$selection_rule, "\n")
      cat("  - effect_neighborhood:", args_forestsearch_call$effect_neighborhood, "\n")
      cat("  Bootstrap-specific overrides:\n")
      cat("  - dina_res / dina_cuts: NULL (forces DINA re-fit + re-selection)\n")
      cat("  - GRF / LASSO / consistency search: bypassed in DINA mode\n")
    } else {
      cat("  - maxk:", args_forestsearch_call$maxk, "\n")
      cat("  - fs.splits:", args_forestsearch_call$fs.splits, "\n")
      cat("  - max_subgroups_search:", args_forestsearch_call$max_subgroups_search, "\n")
      cat("  - hr.threshold:", args_forestsearch_call$hr.threshold, "\n")
      cat("  - hr.consistency:", args_forestsearch_call$hr.consistency, "\n")
      cat("  - pconsistency.threshold:", args_forestsearch_call$pconsistency.threshold, "\n")
      cat("  - n.min:", args_forestsearch_call$n.min, "\n")
      cat("  - use_twostage:", args_forestsearch_call$use_twostage, "\n")
      if (isTRUE(args_forestsearch_call$use_twostage) &&
          length(args_forestsearch_call$twostage_args) > 0) {
        cat("  - twostage_args:\n")
        cat("      n.splits.screen:", args_forestsearch_call$twostage_args$n.splits.screen, "\n")
        cat("      batch.size:", args_forestsearch_call$twostage_args$batch.size, "\n")
      }
      cat("  - use_lasso:", args_forestsearch_call$use_lasso, "\n")
      cat("  - use_grf:", args_forestsearch_call$use_grf, "\n")
      cat("  Bootstrap-specific overrides:\n")
      cat("  - grf_res: NULL (forces re-selection)\n")
      cat("  - grf_cuts: NULL (forces re-selection)\n")
    }
    cat("  - parallel_args: sequential (prevents nested parallelism)\n")
    cat("  - details: per-iteration (show_three -> first 3 only)\n")
    cat("  - plot.sg: FALSE\n")
    cat("  - plot.grf: FALSE\n")
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
    Hc_obs = Hc_obs,
    seed = boot_seed,
    estimator_fn = estimator_fn_boot,
    boot_index_mat = boot_index_mat
  )

  options(warn = old_warn)

  # =======================================================================
  # SECTION 9: POST-PROCESSING AND CONFIDENCE INTERVALS
  # =======================================================================

  # Zero-success guard: bias correction is only computed for resamples that
  # re-identified a subgroup (run_bootstrap$sg.harm non-NULL).  If NO resample
  # did, every *_biasadj_* is NA and the variance/CI machinery has nothing to
  # work with (get_targetEst() would return NaN).  Report this explicitly as a
  # partial result rather than emitting a NaN-filled table.  This is most
  # likely a configuration signal (e.g. an effect floor the resamples never
  # clear), not a numerical failure.
  n_success <- sum(!is.na(results$H_biasadj_2))
  if (n_success == 0L) {
    warning("No bootstrap resample re-identified a subgroup (0/", nb_boots,
            " succeeded), so no bias-corrected estimates could be computed. ",
            "This usually means the selection criteria (e.g. the effect/HR ",
            "floor) are not reached on resampled data -- inspect details = TRUE ",
            "output and consider relaxing the floor. Returning partial results ",
            "without confidence intervals.", call. = FALSE)
    out <- list(
      results = results,
      SG_CIs = NULL,
      FSsg_tab = NULL,
      Ystar_mat = Ystar_mat,
      H_estimates = NULL,
      Hc_estimates = NULL,
      nb_boots       = nb_boots,
      boot_success_rate = 0,
      original_sg    = fs.est$sg.harm,
      outcome_type   = if (is_glm) outcome_type else "survival",
      effect_measure = if (is_glm) args_forestsearch_call$effect_measure
                       else "HR",
      est.scale      = args_forestsearch_call$est.scale
    )
    class(out) <- c("fs_bootstrap", "list")
    return(out)
  }

  est.scale <- args_forestsearch_call$est.scale

  # For GLM outcomes, determine whether estimates are on the log scale.
  # Log-scale: HR, OR, RR, IRR → est.loghr = TRUE
  # Identity-scale: RD, IRD, MD → est.loghr = FALSE
  if (is_glm) {
    effect_measure <- args_forestsearch_call$effect_measure
    est_loghr <- effect_measure %in% c("OR", "RR", "IRR")
    # Identity-scale measures use "hr" est.scale for CI formatting
    if (!est_loghr) est.scale <- "hr"
  } else {
    est_loghr <- TRUE
  }

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
      est.loghr = est_loghr
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
      est.loghr = est_loghr
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
      Hc_estimates = NULL,
      # Metadata mirrors the normal-return path so print.fs_bootstrap can
      # still produce a partial summary on the early-return object.
      nb_boots       = nb_boots,
      original_sg    = fs.est$sg.harm,
      outcome_type   = if (is_glm) outcome_type else "survival",
      effect_measure = if (is_glm) args_forestsearch_call$effect_measure
                       else "HR",
      est.scale      = est.scale
    )
    class(out) <- c("fs_bootstrap", "list")
    return(out)
  }

  # =======================================================================
  # SECTION 10: FORMAT CONFIDENCE INTERVALS
  # =======================================================================

  # Format CIs (handling NULL estimates gracefully).  digits controls
  # construction-time precision; downstream display can round down via
  # summarize_bootstrap_results(digits = ...) but cannot recover precision
  # rounded away here.
  H_res1 <- if (!is.null(H_estimates)) {
    format_CI(H_estimates, c("H0", "H0_lower", "H0_upper"), digits = digits)
  } else {
    "NA (NA, NA)"
  }

  H_res2 <- if (!is.null(H_estimates)) {
    format_CI(H_estimates, c("H2", "H2_lower", "H2_upper"), digits = digits)
  } else {
    "NA (NA, NA)"
  }

  Hc_res1 <- if (!is.null(Hc_estimates)) {
    format_CI(Hc_estimates, c("H0", "H0_lower", "H0_upper"), digits = digits)
  } else {
    "NA (NA, NA)"
  }

  Hc_res2 <- if (!is.null(Hc_estimates)) {
    format_CI(Hc_estimates, c("H2", "H2_lower", "H2_upper"), digits = digits)
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

  FSsg_tab <- NULL

  if (is_glm) {
    # GLM summary table using SG_tab_estimates_glm.
    # Note: outcome_type MUST be passed explicitly.  The function defaults
    # to "binary", which mislabels columns ("Rate(C)/(T)" vs "Mean(C)/(T)")
    # and percent-formats values for continuous/count outcomes.
    FSsg_tab <- tryCatch(
      SG_tab_estimates_glm(
        df             = fs.est$df.est,
        SG_flag        = "treat.recommend",
        outcome.name   = args_forestsearch_call$outcome.name,
        treat.name     = args_forestsearch_call$treat.name,
        estimator_fn   = estimator_fn_boot,
        effect_measure = args_forestsearch_call$effect_measure,
        outcome_type   = outcome_type,
        effect_a_1     = SG_CIs$Hc_bc,
        effect_a_0     = SG_CIs$H_bc,
        sg1_name       = "Recmnd",
        sg0_name       = "Qstnbl",
        est.scale      = est.scale,
        digits         = digits
      ),
      error = function(e) {
        warning("GLM SG_tab_estimates failed: ", e$message)
        list(
          outcome_type   = outcome_type,
          effect_measure = args_forestsearch_call$effect_measure,
          H_raw  = H_res1, H_bc  = H_res2,
          Hc_raw = Hc_res1, Hc_bc = Hc_res2
        )
      }
    )
  } else {
    # Survival path: existing SG_tab_estimates (unchanged)
    # Determine HR parameters based on scale
    if (est.scale == "1/hr") {
      hr_1a <- SG_CIs$H_bc
      hr_0a <- SG_CIs$Hc_bc
      sg0_name <- "Qstnbl"
      sg1_name <- "Recmnd"
    } else {
      hr_1a <- SG_CIs$Hc_bc
      hr_0a <- SG_CIs$H_bc
      sg0_name <- "Qstnbl"
      sg1_name <- "Recmnd"
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
  }

  # =======================================================================
  # SECTION 12: COMPILE OUTPUT
  # =======================================================================

  out <- list(
    results = results,
    SG_CIs = SG_CIs,
    FSsg_tab = FSsg_tab,
    Ystar_mat = Ystar_mat,
    H_estimates = H_estimates,
    Hc_estimates = Hc_estimates,
    # Metadata needed by print.fs_bootstrap so that the print method does
    # not need to retain a reference to fs.est.  Added in v0.2.0 alongside
    # the fs_bootstrap class tag.
    nb_boots       = nb_boots,
    original_sg    = fs.est$sg.harm,
    outcome_type   = if (is_glm) outcome_type else "survival",
    effect_measure = if (is_glm) args_forestsearch_call$effect_measure
                     else "HR",
    est.scale      = est.scale
  )

  class(out) <- c("fs_bootstrap", "list")
  return(out)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for ForestSearch Bootstrap Results
#'
#' Prints a concise but informative summary of an \code{fs_bootstrap}
#' object returned by \code{\link{forestsearch_bootstrap_dofuture}}.
#' Covers identification rate, agreement with the primary subgroup,
#' top identified subgroups, and bias-corrected effect estimates.
#'
#' @details
#' The print method is deliberately richer than the corresponding
#' \code{print.fs_kfold} / \code{print.fs_tenfold} methods because the
#' bootstrap is the primary inferential machinery of the package, while
#' the CV routines are complementary diagnostics.  For the full
#' per-factor breakdown, consistency distribution, size distribution,
#' and GRF cut tabulation, call \code{\link{summarize_bootstrap_subgroups}}.
#'
#' @param x An \code{fs_bootstrap} object.
#' @param top_n Integer.  Maximum number of top identified subgroups
#'   to print.  Default 3.
#' @param ... Additional arguments (ignored; present for S3
#'   consistency).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso
#' \code{\link{forestsearch_bootstrap_dofuture}} to produce the object.
#' \code{\link{summarize_bootstrap_subgroups}} for full tabulations.
#'
#' @examples
#' \dontrun{
#' fs <- forestsearch(
#'   df.analysis = mydata,
#'   confounders.name = c("age", "sex", "biomarker"),
#'   outcome.name = "time", event.name = "status", treat.name = "treat"
#' )
#' fs_bc <- forestsearch_bootstrap_dofuture(fs, nb_boots = 500)
#' print(fs_bc)          # rich default summary
#' print(fs_bc, top_n = 5)
#' }
#'
#' @export
print.fs_bootstrap <- function(x, top_n = 3L, ...) {

  # ---------------------------------------------------------------------------
  # Header
  # ---------------------------------------------------------------------------
  cat("ForestSearch Bootstrap Results\n")
  cat("==============================\n")

  # ---------------------------------------------------------------------------
  # Iteration count, timing
  # ---------------------------------------------------------------------------
  nb <- if (!is.null(x$nb_boots)) x$nb_boots else nrow(x$results)
  cat(sprintf("Iterations: %d\n", nb))

  timing <- attr(x$results, "timing")
  if (!is.null(timing) && !is.null(timing$total_minutes)) {
    cat(sprintf("Total time: %.2f minutes  (%.2f sec/iter)\n",
                timing$total_minutes, timing$avg_seconds_per_boot))
  }

  # ---------------------------------------------------------------------------
  # Outcome context (helps when printing a loaded/cached object)
  # ---------------------------------------------------------------------------
  if (!is.null(x$outcome_type)) {
    cat(sprintf("Outcome type: %s  (effect measure: %s)\n",
                x$outcome_type,
                x$effect_measure %||% "HR"))
  }

  # ---------------------------------------------------------------------------
  # Identification rate
  # ---------------------------------------------------------------------------
  if (is.null(x$results) || nrow(x$results) == 0L) {
    cat("\nNo bootstrap results to summarize.\n")
    return(invisible(x))
  }

  # Prefer any_found when present (v0.2.0+), fall back to Pcons for
  # backward compatibility with cached objects from older versions.
  any_found_vec <- if ("any_found" %in% names(x$results)) {
    x$results$any_found == 1L
  } else {
    !is.na(x$results$Pcons)
  }
  n_found <- sum(any_found_vec, na.rm = TRUE)
  pct_found <- 100 * n_found / nb

  cat(sprintf(
    "\nIterations identifying a subgroup: %d / %d (%.1f%%)\n",
    n_found, nb, pct_found
  ))

  # ---------------------------------------------------------------------------
  # Bias-corrected effect estimates (if available)
  # ---------------------------------------------------------------------------
  if (!is.null(x$SG_CIs)) {
    effect_label <- x$effect_measure %||% "HR"
    cat(sprintf("\nSubgroup effect estimates (%s, 95%% CI):\n", effect_label))
    cat(sprintf("  H  (Questionable)  unadjusted:     %s\n",
                x$SG_CIs$H_raw  %||% "NA"))
    cat(sprintf("                     bias-corrected: %s\n",
                x$SG_CIs$H_bc   %||% "NA"))
    cat(sprintf("  Hc (Recommended)   unadjusted:     %s\n",
                x$SG_CIs$Hc_raw %||% "NA"))
    cat(sprintf("                     bias-corrected: %s\n",
                x$SG_CIs$Hc_bc  %||% "NA"))
  } else {
    cat("\nSubgroup effect estimates: not available (both H and Hc estimates failed).\n")
  }

  # ---------------------------------------------------------------------------
  # Top identified subgroups
  # ---------------------------------------------------------------------------
  if (n_found > 0L) {

    # Build canonical-form subgroup strings from M.1..M.7
    m_cols <- intersect(paste0("M.", 1:7), names(x$results))
    if (length(m_cols) > 0L) {
      # Coerce to data.frame BEFORE subsetting so standard `[` semantics
      # apply.  If `x$results` is a data.table (which is what
      # `bootstrap_results()` returns), writing the inner subset first
      # would trigger NSE and interpret `m_cols` as a column name.
      res_ok <- as.data.frame(x$results)[any_found_vec, m_cols, drop = FALSE]

      sg_strs <- vapply(seq_len(nrow(res_ok)), function(i) {
        vals <- unlist(res_ok[i, , drop = TRUE])
        vals <- vals[!is.na(vals) & vals != ""]
        if (length(vals) == 0L) NA_character_
        else paste(sort(vals), collapse = " & ")
      }, character(1))

      sg_strs <- sg_strs[!is.na(sg_strs)]
      if (length(sg_strs) > 0L) {
        tt <- sort(table(sg_strs), decreasing = TRUE)
        k <- min(as.integer(top_n), length(tt))
        cat(sprintf("\nTop %d identified subgroup(s):\n", k))
        for (i in seq_len(k)) {
          cat(sprintf("  %d. %s  (%d / %d, %.1f%%)\n",
                      i, names(tt)[i], as.integer(tt[i]),
                      n_found, 100 * as.integer(tt[i]) / n_found))
        }
      }
    }

    # -------------------------------------------------------------------------
    # Agreement with original (primary-analysis) subgroup
    # -------------------------------------------------------------------------
    if (!is.null(x$original_sg) && length(x$original_sg) > 0L) {
      orig_chars <- as.character(x$original_sg)
      orig_chars <- orig_chars[!is.na(orig_chars) & orig_chars != ""]

      if (length(orig_chars) > 0L && length(m_cols) > 0L) {
        orig_canonical <- paste(sort(orig_chars), collapse = " & ")

        # sg_strs already computed above if m_cols existed
        exact_n <- sum(sg_strs == orig_canonical, na.rm = TRUE)

        # Partial match: any shared factor
        partial_n <- sum(vapply(sg_strs, function(s) {
          if (is.na(s)) return(FALSE)
          any(strsplit(s, " & ", fixed = TRUE)[[1]] %in% orig_chars)
        }, logical(1)), na.rm = TRUE)

        cat(sprintf("\nAgreement with primary subgroup %s:\n", orig_canonical))
        cat(sprintf("  Exact match:   %d / %d (%.1f%%)\n",
                    exact_n, n_found, 100 * exact_n / n_found))
        cat(sprintf("  Partial match (shared factor): %d / %d (%.1f%%)\n",
                    partial_n, n_found, 100 * partial_n / n_found))
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Footer
  # ---------------------------------------------------------------------------
  cat("\nUse summarize_bootstrap_subgroups() for full diagnostics.\n")

  invisible(x)
}

# Helper: %||% (NULL-coalescing operator), kept local to avoid a package-
# wide dependency on rlang.  Used only in print.fs_bootstrap above.
`%||%` <- function(a, b) if (is.null(a)) b else a
