#' @title ForestSearch: Exploratory Subgroup Identification
#'
#' @description
#' Identifies subgroups with differential treatment effects in clinical trials
#' using a combination of Generalized Random Forests (GRF), LASSO variable
#' selection, and exhaustive combinatorial search with split-sample validation.
#'
#' @param df.analysis Data frame. Analysis dataset with required columns.
#' @param outcome.name Character. Name of time-to-event outcome variable. Default "tte".
#' @param event.name Character. Name of event indicator (1=event, 0=censored). Default "event".
#' @param treat.name Character. Name of treatment variable (1=treatment, 0=control). Default "treat".
#' @param id.name Character. Name of subject ID variable. Default "id".
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param flag_harm.name Character. Name of true harm flag for simulations (optional).
#' @param confounders.name Character vector. Names of candidate subgroup-defining variables.
#' @param parallel_args List. Parallel processing configuration:
#'   \describe{
#'     \item{plan}{Character. One of "multisession", "multicore", "callr", "sequential"}
#'     \item{workers}{Integer. Number of parallel workers}
#'     \item{show_message}{Logical. Show parallel setup messages}
#'   }
#' @param df.predict Data frame. Prediction dataset (optional).
#' @param df.test Data frame. Test dataset (optional).
#' @param is.RCT Logical. Is this a randomized controlled trial? Default TRUE.
#' @param seedit Integer. Random seed. Default 8316951.
#' @param est.scale Character. Estimation scale ("hr" or "rmst"). Default "hr".
#' @param use_lasso Logical. Use LASSO for variable selection. Default TRUE.
#' @param use_grf Logical. Use GRF for variable importance. Default TRUE.
#' @param grf_res GRF results object (optional, for reuse).
#' @param grf_cuts List. Custom GRF cut points (optional).
#' @param max_n_confounders Integer. Maximum confounders to consider. Default 1000.
#' @param grf_depth Integer. GRF tree depth. Default 2.
#' @param dmin.grf Integer. Minimum events for GRF. Default 12.
#' @param frac.tau Numeric. Fraction of tau for RMST. Default 0.6.
#' @param return_selected_cuts_only Logical. If TRUE (default), GRF returns only cuts from the
#'   tree depth that identified the selected subgroup meeting `dmin.grf`. If FALSE
#'   returns all cuts from all fitted trees (depths 1 through `grf_depth`).
#'   See \code{\link{grf.subg.harm.survival}} for details.
#' @param conf_force Character vector. Variables to force include (optional).
#' @param defaultcut_names Character vector. Default cut variable names (optional).
#' @param cut_type Character. Cut type ("default" or "custom"). Default "default".
#' @param exclude_cuts Character vector. Variables to exclude from cutting (optional).
#' @param replace_med_grf Logical. Replace median with GRF cuts. Default FALSE.
#' @param cont.cutoff Integer. Cutoff for continuous vs categorical. Default 4.
#' @param conf.cont_medians Named numeric vector. Median values for continuous variables (optional).
#' @param conf.cont_medians_force Named numeric vector. Forced median values (optional).
#' @param n.min Integer. Minimum subgroup size. Default 60.
#' @param hr.threshold Numeric. Minimum HR for candidate subgroups. Default 1.25.
#' @param hr.consistency Numeric. Minimum HR for consistency validation. Default 1.0.
#' @param sg_focus Character. Subgroup selection focus. One of "hr", "hrMaxSG", "maxSG",
#'   "hrMinSG", "minSG". Default "hr".
#' @param stop_threshold Numeric in \code{[0, 1]} or \code{NULL}.
#'   Early stopping threshold for consistency evaluation. When a candidate
#'   subgroup's consistency probability (\code{Pcons}) meets or exceeds this
#'   threshold, evaluation stops early -- remaining candidates are skipped.
#'   Set to \code{NULL} to disable early stopping and force evaluation of
#'   all candidates up to \code{max_subgroups_search}. Default \code{0.95}.
#'
#'   \strong{Note:} Values > 1.0 are not permitted. To disable early
#'   stopping, use \code{stop_threshold = NULL}, not a value above 1.
#'
#'   Automatically reset to \code{NULL} (with a warning) when
#'   \code{sg_focus} is \code{"hrMaxSG"} or \code{"hrMinSG"}, as these
#'   compound criteria require comparing HR \emph{and} size across all
#'   candidates.
#' @param fs.splits Integer. Number of splits for consistency evaluation (or maximum
#'   splits when \code{use_twostage = TRUE}). Default 1000.
#' @param m1.threshold Numeric. Maximum median survival threshold. Default Inf.
#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.
#' @param showten_subgroups Logical. Show top 10 subgroups. Default FALSE.
#' @param d0.min Integer. Minimum control arm events. Default 12.
#' @param d1.min Integer. Minimum treatment arm events. Default 12.
#' @param max.minutes Numeric. Maximum search time in minutes. Default 3.
#' @param minp Numeric. Minimum prevalence threshold. Default 0.025.
#' @param details Logical. Print progress details. Default FALSE.
#' @param maxk Integer. Maximum number of factors per subgroup. Default 2.
#' @param by.risk Integer. Risk table interval. Default 12.
#' @param plot.sg Logical. Plot subgroup survival curves. Default FALSE.
#' @param plot.grf Logical. Plot GRF results. Default FALSE.
#' @param max_subgroups_search Integer. Maximum subgroups to evaluate. Default 10.
#' @param vi.grf.min Numeric. Minimum GRF variable importance. Default -0.2.
#' @param use_twostage Logical. Use two-stage sequential consistency algorithm for
#'   improved performance. Default FALSE for backward compatibility. When TRUE,
#'   \code{fs.splits} becomes the maximum number of splits, and early stopping
#'   is enabled. See Details.
#' @param twostage_args List. Parameters for two-stage algorithm (only used when
#'   \code{use_twostage = TRUE}):
#'   \describe{
#'     \item{n.splits.screen}{Integer. Splits for Stage 1 screening. Default 30.}
#'     \item{screen.threshold}{Numeric. Consistency threshold for Stage 1. Default
#'       is automatically calculated to provide ~2.5 SE margin.}
#'     \item{batch.size}{Integer. Splits per batch in Stage 2. Default 20.}
#'     \item{conf.level}{Numeric. Confidence level for early stopping. Default 0.95.}
#'     \item{min.valid.screen}{Integer. Minimum valid Stage 1 splits. Default 10.}
#'   }
#' @param outcome_type Character. One of \code{"survival"} (default),
#'   \code{"binary"}, \code{"continuous"}, or \code{"count"}.
#' @param effect_measure Character or \code{NULL}. Effect measure for GLM
#'   outcomes.  For binary: \code{"RD"}, \code{"OR"}, \code{"RR"},
#'   \code{"IRR"}, \code{"IRD"}.  For continuous: \code{"MD"}.
#'   For count: \code{"IRR"} (default), \code{"IRD"}.
#'   Default depends on \code{outcome_type}.
#' @param offset.name Character or \code{NULL}. Name of the follow-up time
#'   column for rate-based measures (IRR, IRD).
#' @param adverse_outcome Logical or \code{NULL}. If \code{TRUE}, higher
#'   outcome values indicate harm (default for binary and count).  If \code{FALSE},
#'   higher values indicate benefit (default for continuous).
#' @param overdispersion Character. Overdispersion correction for
#'   \code{outcome_type = "count"}.  One of \code{"none"} (default,
#'   standard Poisson), \code{"quasi"} (quasi-Poisson), or \code{"negbin"}
#'   (negative binomial via \code{MASS::glm.nb}).  Ignored for other
#'   outcome types.
#' @param grf_count_transform Character. Transformation applied to the
#'   count outcome before passing to \code{grf::causal_forest()} for
#'   factor screening.  \code{"log"} (default) applies
#'   \eqn{\log(Y + 0.5)}; \code{"identity"} passes raw counts.
#'   Ignored for non-count outcomes.
#' @param ps_method Character or \code{NULL}. Propensity score estimation
#'   method: \code{"grf"}, \code{"lasso"}, \code{"logistic"}, or
#'   \code{"none"}.  Default: \code{"none"} for RCT, \code{"grf"} for
#'   observational.
#' @param ps_adjust_method Character. PS adjustment method: \code{"none"}
#'   (default), \code{"iptw"} (stabilized inverse probability of treatment
#'   weights), or \code{"dr_gcomp"} (IPS covariate).
#' @param ps_hat Numeric vector or \code{NULL}. User-supplied propensity
#'   scores.  Must have length equal to \code{nrow(df.analysis)}.
#'
#' @return A list of class "forestsearch" containing:
#'   \describe{
#'     \item{grp.consistency}{Consistency evaluation results including:
#'       \itemize{
#'         \item out_sg: Selected subgroup based on sg_focus
#'         \item sg_focus: Focus criterion used
#'         \item df_flag: Treatment recommendations
#'         \item algorithm: "twostage" or "fixed"
#'         \item n_candidates_evaluated: Number evaluated
#'         \item n_passed: Number passing threshold
#'       }}
#'     \item{find.grps}{Subgroup search results}
#'     \item{confounders.candidate}{Candidate confounders considered}
#'     \item{confounders.evaluated}{Confounders after variable selection}
#'     \item{df.est}{Analysis data with treatment recommendations}
#'     \item{df.predict}{Prediction data with recommendations (if provided)}
#'     \item{df.test}{Test data with recommendations (if provided)}
#'     \item{minutes_all}{Total computation time}
#'     \item{grf_res}{GRF results object}
#'     \item{sg_focus}{Subgroup focus criterion used}
#'     \item{sg.harm}{Selected subgroup definition}
#'     \item{grf_cuts}{GRF cut points used}
#'     \item{prop_maxk}{Proportion of max combinations searched}
#'     \item{max_sg_est}{Maximum subgroup HR estimate}
#'     \item{grf_plot}{GRF plot object (if plot.grf = TRUE)}
#'     \item{args_call_all}{All arguments for reproducibility}
#'   }
#'
#' @details
#' \strong{Algorithm Overview:}
#' \enumerate{
#'   \item \strong{Variable Selection}: GRF identifies variables with treatment
#'     effect heterogeneity; LASSO selects most predictive
#'   \item \strong{Subgroup Discovery}: Exhaustive search over factor combinations
#'     up to \code{maxk}
#'   \item \strong{Consistency Validation}: Split-sample validation ensures
#'     reproducibility
#'   \item \strong{Selection}: Choose subgroup based on \code{sg_focus} criterion
#' }
#'
#' \strong{Two-Stage Consistency Algorithm:}
#' When \code{use_twostage = TRUE}, the consistency evaluation uses an optimized
#' algorithm that can provide 3-10x speedup:
#' \itemize{
#'   \item \strong{Stage 1}: Quick screening with \code{n.splits.screen} splits
#'     eliminates clearly non-viable candidates
#'   \item \strong{Stage 2}: Sequential batched evaluation with early stopping
#'     for candidates passing Stage 1
#' }
#'
#' The two-stage algorithm is recommended for:
#' \itemize{
#'   \item Exploratory analyses with many candidate subgroups
#'   \item Large \code{fs.splits} values (>200)
#'   \item Iterative model development
#' }
#'
#' For final regulatory submissions, \code{use_twostage = FALSE} may be preferred
#' for exact reproducibility.
#'
#' @examples
#' \dontrun{
#' # Example 1: Standard analysis (backward compatible)
#' result <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",
#'   hr.threshold = 1.25,
#'   pconsistency.threshold = 0.90,
#'   fs.splits = 400,
#'   details = TRUE
#' )
#'
#' # Example 2: Fast exploratory analysis with two-stage
#' result_fast <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "maxSG",
#'   hr.threshold = 1.15,
#'   pconsistency.threshold = 0.85,
#'   fs.splits = 500,
#'   use_twostage = TRUE,
#'   details = TRUE
#' )
#'
#' # Example 3: Two-stage with custom parameters
#' result_custom <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",
#'   hr.threshold = 1.3,
#'   pconsistency.threshold = 0.95,
#'   fs.splits = 600,
#'   use_twostage = TRUE,
#'   twostage_args = list(
#'     n.splits.screen = 50,
#'     batch.size = 25,
#'     conf.level = 0.99
#'   ),
#'   parallel_args = list(plan = "multisession", workers = 4),
#'   details = TRUE
#' )
#' }
#'
#' @references
#' \itemize{
#'   \item FDA Guidance for Industry: Enrichment Strategies for Clinical Trials
#'   \item Athey & Imbens (2016). Recursive partitioning for heterogeneous
#'     causal effects. PNAS.
#'   \item Wager & Athey (2018). Estimation and inference of heterogeneous
#'     treatment effects using random forests. JASA.
#' }
#'
#' @seealso
#' \code{\link{subgroup.consistency}} for consistency evaluation details
#' \code{\link{forestsearch_bootstrap_dofuture}} for bootstrap inference
#' \code{\link{forestsearch_Kfold}} for cross-validation
#'
#' Package website: \url{https://larry-leon.github.io/forestsearch/}
#'
#' Source code: \url{https://github.com/larry-leon/forestsearch}
#'
#' @importFrom survival coxph Surv
#' @importFrom grf causal_survival_forest variable_importance
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table setorder
#' @importFrom stats quantile sd median
#' @importFrom stats complete.cases
#' @importFrom future.apply future_lapply
#' @importFrom randomForest randomForest
#' @importFrom weightedsurv df_counting
#' @export

forestsearch <- function(df.analysis,
                         outcome.name = "tte",
                         event.name = "event",
                         treat.name = "treat",
                         id.name = "id",
                         potentialOutcome.name = NULL,
                         flag_harm.name = NULL,
                         confounders.name = NULL,
                         parallel_args = list(plan = "callr", workers = 6, show_message = TRUE),
                         df.predict = NULL,
                         df.test = NULL,
                         is.RCT = TRUE,
                         seedit = 8316951,
                         est.scale = "hr",
                         use_lasso = TRUE,
                         use_grf = TRUE,
                         grf_res = NULL,
                         grf_cuts = NULL,
                         max_n_confounders = 1000,
                         grf_depth = 2,
                         dmin.grf = 12,
                         frac.tau = 0.6,
                         return_selected_cuts_only = TRUE,
                         conf_force = NULL,
                         defaultcut_names = NULL,
                         cut_type = "default",
                         exclude_cuts = NULL,
                         replace_med_grf = FALSE,
                         cont.cutoff = 4,
                         conf.cont_medians = NULL,
                         conf.cont_medians_force = NULL,
                         n.min = 60,
                         hr.threshold = 1.25,
                         hr.consistency = 1.0,
                         sg_focus = "hr",
                         fs.splits = 1000,
                         m1.threshold = Inf,
                         pconsistency.threshold = 0.90,
                         stop_threshold = 0.95,
                         showten_subgroups = FALSE,
                         d0.min = 12,
                         d1.min = 12,
                         max.minutes = 3,
                         minp = 0.025,
                         details = FALSE,
                         maxk = 2,
                         by.risk = 12,
                         plot.sg = FALSE,
                         plot.grf = FALSE,
                         max_subgroups_search = 10,
                         vi.grf.min = -0.2,
                         # NEW: Two-stage consistency parameters
                         use_twostage = TRUE,
                         twostage_args = list(),
                         # NEW: GLM outcome support
                         outcome_type = c("survival", "binary", "continuous",
                                          "count"),
                         effect_measure = NULL,
                         offset.name = NULL,
                         adverse_outcome = NULL,
                         overdispersion = c("none", "quasi", "negbin"),
                         grf_count_transform = c("log", "identity"),
                         # Propensity score adjustment
                         ps_method = NULL,
                         ps_adjust_method = c("none", "iptw", "dr_gcomp"),
                         ps_hat = NULL) {

  # ===========================================================================
  # SECTION 1: CAPTURE ALL ARGUMENTS FOR REPRODUCIBILITY
  # ===========================================================================

  args_names <- names(formals())
  args_call_all <- mget(args_names, envir = environment())

  # ===========================================================================
  # SECTION 2: VALIDATE INPUTS
  # ===========================================================================

  # Validate parallel arguments
  if (length(parallel_args) > 0) {
    allowed_plans <- c("multisession", "multicore", "callr", "sequential")
    plan_type <- parallel_args$plan
    n_workers <- parallel_args$workers
    max_cores <- parallel::detectCores()

    if (is.null(plan_type)) {
      stop("parallel_args$plan must be specified.")
    }

    if (!plan_type %in% allowed_plans) {
      stop("plan_type must be one of: ", paste(allowed_plans, collapse = ", "))
    }

    if (is.null(n_workers) || !is.numeric(n_workers) || n_workers < 1) {
      n_workers <- 1
    } else {
      n_workers <- min(n_workers, max_cores)
    }
  }

  # Early stopping is incompatible with compound sg_focus criteria that
  # require comparing HR *and* size across all candidates.
  if (!is.null(stop_threshold) && sg_focus %in% c("hrMaxSG", "hrMinSG")) {

    # Detect whether user explicitly set stop_threshold (vs. inheriting default)
    user_explicit <- !missing(stop_threshold)

    if (user_explicit) {
      warning(
        "stop_threshold = ", stop_threshold,
        " reset to NULL for sg_focus = '", sg_focus, "'.\n",
        "Compound selection criteria require evaluating all candidates.\n",
        "To suppress this warning, pass stop_threshold = NULL explicitly.",
        call. = FALSE
      )
    }

    stop_threshold <- NULL
    args_call_all$stop_threshold <- NULL
  }

  # Validate two-stage parameters
 if (!is.logical(use_twostage) || length(use_twostage) != 1) {
    stop("'use_twostage' must be a single logical value (TRUE or FALSE)")
  }

  if (use_twostage && length(twostage_args) > 0) {
    valid_ts_args <- c("n.splits.screen", "screen.threshold", "batch.size",
                       "conf.level", "min.valid.screen")
    invalid_args <- setdiff(names(twostage_args), valid_ts_args)
    if (length(invalid_args) > 0) {
      warning("Unknown twostage_args parameters ignored: ",
              paste(invalid_args, collapse = ", "))
    }
  }

  # ===========================================================================
  # SECTION 2B: GLM OUTCOME TYPE SETUP
  # ===========================================================================

  outcome_type <- match.arg(outcome_type)
  overdispersion <- match.arg(overdispersion)
  grf_count_transform <- match.arg(grf_count_transform)

  # Store resolved value so bootstrap/CV pick up the scalar, not the default vector
  args_call_all$outcome_type <- outcome_type
  args_call_all$overdispersion <- overdispersion
  args_call_all$grf_count_transform <- grf_count_transform

  # Resolve adverse_outcome default:
  #   Binary outcomes: TRUE (event = adverse, the clinical trial convention)
  #   Count outcomes:  TRUE (higher count = more events = worse)
  #   Continuous outcomes: FALSE (user must set TRUE if higher Y = worse)
  #   Survival: not used (causal_survival_forest handles sign correctly)
  if (is.null(adverse_outcome)) {
    adverse_outcome <- (outcome_type %in% c("binary", "count"))
  }
  args_call_all$adverse_outcome <- adverse_outcome

  # Build effect estimator closure for non-survival outcomes
  estimator_fn <- NULL
  consistency_threshold <- NULL
  effect_threshold <- NULL


  if (outcome_type != "survival") {

    # Resolve default effect measure
    if (is.null(effect_measure)) {
      effect_measure <- switch(outcome_type,
        binary     = "RD",
        continuous = "MD",
        count      = "IRR"
      )
    }
    # Store resolved effect_measure
    args_call_all$effect_measure <- effect_measure

    # Validate offset requirement for rate-based measures
    if (effect_measure %in% c("IRR", "IRD") && is.null(offset.name)) {
      stop(
        "effect_measure = '", effect_measure,
        "' requires `offset.name` (follow-up time column).",
        call. = FALSE
      )
    }

    # Build the estimator closure
    estimator_fn <- make_effect_estimator(
      outcome_type   = outcome_type,
      treat.name     = treat.name,
      outcome.name   = outcome.name,
      event.name     = if (outcome_type == "survival") event.name else NULL,
      offset.name    = offset.name,
      effect_measure = effect_measure
    )

    # Resolve screening and consistency thresholds from effect_measure
    # The user can override via hr.threshold and hr.consistency, which
    # serve as the generic threshold params regardless of outcome type.
    # Defaults differ by measure:
    effect_threshold <- hr.threshold
    consistency_threshold <- hr.consistency

    # For identity-scale measures, consistency_threshold = 0 means
    # "both halves show any treatment harm" (analogous to HR > 1.0)
    if (effect_measure %in% c("RD", "IRD", "MD")) {
      if (missing(hr.consistency) || hr.consistency == 1.0) {
        consistency_threshold <- 0.0
      }
      if (missing(hr.threshold) || hr.threshold == 1.25) {
        effect_threshold <- switch(effect_measure,
          RD  = 0.05,
          IRD = 0.01,
          MD  = 0.0
        )
      }
    } else {
      # Log-scale measures (OR, RR, IRR): convert to log scale
      if (effect_threshold > 0) {
        effect_threshold <- log(effect_threshold)
      }
      if (consistency_threshold > 0) {
        consistency_threshold <- log(consistency_threshold)
      }
    }

    # Resolve dmin.grf for GLM outcomes.
    # The survival default (dmin.grf = 12) is on the RMST scale (months)
    # and would filter out all candidates for binary/continuous outcomes
    # where DR score differences are on the probability/mean scale (0--1).
    # Default to 0.0 (any positive harm qualifies) unless the user
    # explicitly specified a value.
    if (missing(dmin.grf) || dmin.grf == 12) {
      dmin.grf <- 0.0
    }
    args_call_all$dmin.grf <- dmin.grf

    if (details) {
      cat("GLM mode: outcome_type =", outcome_type,
          ", effect_measure =", effect_measure, "\n")
      cat("  Screening threshold:", effect_threshold, "\n")
      cat("  Consistency threshold:", consistency_threshold, "\n")
      cat("  dmin.grf:", dmin.grf, "\n")
    }
  }

  # ===========================================================================
  # SECTION 2C: PROPENSITY SCORE ESTIMATION (optional)
  # ===========================================================================

  is_glm <- outcome_type != "survival"
  #
  # When ps_method != "none", estimate P(W=1|X) and use either:
  #   - IPTW: stabilized weights in glm (coefficient IS the effect)
  #   - DR G-comp: Bang & Robins (2005) IPS covariate + G-computation
  #
  # Defaults:
  #   RCT = TRUE  -> ps_method = "none"  (randomization ensures balance)
  #   RCT = FALSE -> ps_method = "grf"   (non-parametric, cross-fitted)
  #   ps_adjust_method = "iptw" (default)
  # ===========================================================================

  # Resolve ps_method default
  if (is.null(ps_method)) {
    ps_method <- if (is.RCT) "none" else "grf"
  }
  args_call_all$ps_method <- ps_method

  # Resolve ps_adjust_method
  ps_adjust_method <- match.arg(ps_adjust_method)
  ps_adjust_resolved <- if (ps_method == "none") "none" else ps_adjust_method
  args_call_all$ps_adjust_method <- ps_adjust_resolved

  adjust_covariates <- NULL

  if (ps_method != "none") {
    if (!is.null(ps_hat)) {
      # User-supplied propensity scores
      if (length(ps_hat) != nrow(df.analysis)) {
        stop("ps_hat must have length equal to nrow(df.analysis)", call. = FALSE)
      }
      df.analysis$ps_hat <- ps_hat
      W_vec <- df.analysis[[treat.name]]
      p_treat <- mean(W_vec, na.rm = TRUE)
      df.analysis$sw <- ifelse(W_vec == 1,
                               p_treat / ps_hat,
                               (1 - p_treat) / (1 - ps_hat))
      df.analysis$ips_covar <- ifelse(W_vec == 1,
                                      1 / ps_hat,
                                      1 / (1 - ps_hat))
    } else {
      # Estimate PS using the requested method
      ps_result <- estimate_propensity_scores(
        data             = df.analysis,
        treat.name       = treat.name,
        confounders.name = confounders.name,
        method           = ps_method,
        seed             = if (exists("seedit")) seedit else 8316951L
      )
      df.analysis <- ps_result$data
      if (details) {
        cat("  PS estimation: method =", ps_result$method,
            ", adjust =", ps_adjust_resolved,
            ", trimmed =", ps_result$trimmed, "\n")
      }
    }

    # Rebuild estimator closure with PS adjustment (GLM outcomes only)
    if (is_glm) {
      estimator_fn <- make_effect_estimator(
        outcome_type     = outcome_type,
        treat.name       = treat.name,
        outcome.name     = outcome.name,
        offset.name      = offset.name,
        effect_measure   = effect_measure,
        ps_adjust_method = ps_adjust_resolved
      )
    }
  }

  args_call_all$adjust_covariates <- adjust_covariates

  # ===========================================================================
  # SECTION 3: INITIALIZE TIMING AND DATA
  # ===========================================================================

  t.start_all <- proc.time()[3]

  df <- df.analysis
  grf_plot <- NULL

  if (details && !is.null(conf_force)) {
    cat("Forced confounders:", paste(conf_force, collapse = ", "), "\n")
  }

  # ===========================================================================
  # SECTION 3A: GRF CUT GENERATION (if use_grf = TRUE)
  # ===========================================================================

  # If using grf and cuts not already populated, run GRF subgroup identification
  if (use_grf && (is.null(grf_res) || is.null(grf_res$tree.cuts))) {

    grf_res <- tryCatch({
      if (outcome_type == "survival") {
        # Survival path: causal_survival_forest (unchanged)
        grf.subg.harm.survival(
          data = df.analysis,
          confounders.name = confounders.name,
          outcome.name = outcome.name,
          event.name = event.name,
          id.name = id.name,
          treat.name = treat.name,
          frac.tau = frac.tau,
          n.min = n.min,
          dmin.grf = dmin.grf,
          RCT = is.RCT,
          details = FALSE,
          maxdepth = grf_depth,
          seedit = seedit,
          return_selected_cuts_only = return_selected_cuts_only
        )
      } else {
        # GLM path: causal_forest (no event/horizon)
        # Map forestsearch outcome_type to grf.subg.harm.glm outcome_type
        grf_glm_args <- list(
          data = df.analysis,
          confounders.name = confounders.name,
          outcome.name = outcome.name,
          treat.name = treat.name,
          id.name = id.name,
          outcome_type = if (outcome_type == "count") "count" else
                         if (outcome_type == "binary") "binary" else "continuous",
          n.min = n.min,
          dmin.grf = dmin.grf,
          RCT = is.RCT,
          details = FALSE,
          maxdepth = grf_depth,
          seedit = seedit,
          return_selected_cuts_only = return_selected_cuts_only,
          adverse_outcome = adverse_outcome
        )
        # Count-specific parameters
        if (outcome_type == "count") {
          grf_glm_args$offset.name         <- offset.name
          grf_glm_args$overdispersion      <- overdispersion
          grf_glm_args$grf_count_transform <- grf_count_transform
        }
        do.call(grf.subg.harm.glm, grf_glm_args)
      }
    },
      error = function(e) {
        warning("GRF analysis failed: ", e$message)
        return(NULL)
      }
    )

    # Extract GRF cuts if successful
    if (!is.null(grf_res) && !inherits(grf_res, "try-error")) {
      grf_cuts <- grf_res$tree.cuts

      if (details) {
        # Concise GRF summary: subgroup found + cuts
        grf_sg <- grf_res$sg.harm.id
        if (!is.null(grf_sg) && length(grf_sg) > 0 &&
            !all(is.na(grf_sg)) && !all(grf_sg == "")) {
          cat("GRF subgroup:", paste(grf_sg, collapse = " & "), "\n")
        } else {
          cat("GRF: no subgroup identified\n")
        }
        cat("GRF cuts identified:", length(grf_cuts), "\n")
        if (length(grf_cuts) > 0) {
          cat("  Cuts:", paste(grf_cuts, collapse = ", "), "\n")
        }
      }

      # Generate GRF plot if requested
      if (plot.grf && !is.null(grf_res$tree)) {
        grf_plot <- tryCatch({
          plot(grf_res$tree)
        }, error = function(e) {
          warning("GRF plot failed: ", e$message)
          NULL
        })
      }
    } else {
      grf_cuts <- NULL
      if (details) {
        cat("GRF analysis did not produce cuts\n")
      }
    }
  }

  # ===========================================================================
  # SECTION 4: DATA PREPARATION (get_FSdata)
  # ===========================================================================

  # For GLM outcomes, get_FSdata validates that event.name is numeric.
  # Map it to something valid so the factor-construction machinery works.
  if (outcome_type == "binary") {
    # Binary outcome IS the event indicator
    if (!event.name %in% names(df.analysis)) {
      event.name <- outcome.name
    }
  } else if (outcome_type %in% c("continuous", "count")) {
    # Continuous / count outcomes have no censoring event; create a placeholder
    if (!event.name %in% names(df.analysis)) {
      df.analysis[[".event_placeholder"]] <- 1L
      event.name <- ".event_placeholder"
      df <- df.analysis
    }
  }

  # Update args_call_all with resolved event.name for downstream calls
  args_call_all$event.name <- event.name

  FSdata <- tryCatch(
    do.call(
      get_FSdata,
      filter_call_args(
        args_call_all,
        get_FSdata,
        list(df.analysis = df.analysis, grf_cuts = grf_cuts,
             details = FALSE)
      )
    ),
    error = function(e) {
      warning("Error in get_FSdata: ", e$message)
      return(NULL)
    }
  )

  if (inherits(FSdata, "try-error") || is.null(FSdata)) {
    warning("FSdata failure - returning NULL result")
    return(list(sg.harm = NULL))
  }

  # Extract FSdata components
  lassoomit <- FSdata$lassoomit
  lassokeep <- FSdata$lassokeep
  df <- FSdata$df

  # Concise LASSO summary
  if (details && use_lasso) {
    n_kept <- length(lassokeep)
    n_total <- n_kept + length(lassoomit)
    cat("Cox-LASSO selected:", n_kept, "of", n_total, "candidate factors\n")
    if (length(lassoomit) > 0) {
      cat("  Omitted:", paste(lassoomit, collapse = ", "), "\n")
    }
  }

  Y <- df[, outcome.name]
  if (outcome_type == "survival") {
    Event <- df[, event.name]
  } else {
    # For GLM outcomes, Event is not meaningful but some downstream
    # code references it.  Use a vector of 1s as placeholder.
    Event <- rep(1L, nrow(df))
  }
  Treat <- df[, treat.name]

  FSconfounders.name <- FSdata$confs_names
  confs_labels <- FSdata$confs

  if (details) {
    cat("Candidate factors:", length(confs_labels), "\n")
    print(confs_labels)
  }

  if (is.null(df.predict)) df.predict <- df

  # ===========================================================================
  # SECTION 5: GRF VARIABLE IMPORTANCE SCREENING
  # ===========================================================================

  if (!is.null(vi.grf.min)) {
    X <- as.matrix(df[, FSconfounders.name])
    X <- apply(X, 2, as.numeric)

    if (outcome_type == "survival") {
      # Survival path: causal_survival_forest (unchanged)
      tau.rmst <- min(c(max(Y[Treat == 1 & Event == 1]), max(Y[Treat == 0 & Event == 1])))

      if (!is.RCT) {
        cs.forest <- try(suppressWarnings(
          grf::causal_survival_forest(X, Y, Treat, Event,
                                       horizon = 0.9 * tau.rmst, seed = 8316951)
        ), TRUE)
      } else {
        cs.forest <- try(suppressWarnings(
          grf::causal_survival_forest(X, Y, Treat, Event, W.hat = 0.5,
                                       horizon = 0.9 * tau.rmst, seed = 8316951)
        ), TRUE)
      }
    } else {
      # GLM path: causal_forest (no event/horizon)
      if (outcome_type == "count") {
        # Count data: variance-stabilising log transform (same as
        # grf.subg.harm.glm with grf_count_transform = "log")
        Y_vi <- if (grf_count_transform == "log") log(Y + 0.5) else Y
      } else {
        # Binary / continuous: flip for adverse outcomes
        Y_vi <- if (adverse_outcome) 1L - Y else Y
      }
      cs.forest <- try(suppressWarnings(
        fit_causal_forest_glm(X, Y_vi, Treat, is.RCT, seedit = 8316951)
      ), TRUE)
    }

    vi.cs <- round(grf::variable_importance(cs.forest), 4)
    vi.cs2 <- data.frame(confs_labels, FSconfounders.name, vi.cs)
    vi.order <- order(vi.cs, decreasing = TRUE)
    vi.cs2 <- vi.cs2[vi.order, ]

    conf.screen <- vi.cs2[, 2]
    vi_ratio <- vi.cs2[, 3] / max(vi.cs2[, 3])
    selected.vars <- which(vi_ratio > vi.grf.min)
    conf.screen <- conf.screen[selected.vars]
    conf.screen <- conf.screen[seq_len(min(length(conf.screen), max_n_confounders))]
  } else {
    conf.screen <- FSconfounders.name
  }

  # ===========================================================================
  # SECTION 6: SUBGROUP SEARCH
  # ===========================================================================

  # Prepare data for subgroup search
  # NOTE: df[, conf.screen] columns are factors with levels {0, 1}.
  # dummy() expands each 2-level factor into two indicator columns
  # (e.g., q1 -> q1.0, q1.1), so the search explores BOTH directions
  # of each cut (subgroup and complement).

  df.confounders <- df[, conf.screen, drop = FALSE]
  df.confounders <- dummy(df.confounders)


  id <- df[, c(id.name)]
  df.fs <- data.frame(Y, Event, Treat, id, df.confounders)
  Z <- as.matrix(df.confounders)
  colnames(Z) <- names(df.confounders)

  # For GLM outcomes: the estimator closure expects columns with the original
  # names (e.g., "treat", "response"), not the renamed Y/Event/Treat.
  # Attach the original columns so the closure can find them.
  if (!is.null(estimator_fn)) {
    if (!outcome.name %in% names(df.fs) && outcome.name %in% names(df)) {
      df.fs[[outcome.name]] <- df[[outcome.name]]
    }
    if (!treat.name %in% names(df.fs) && treat.name %in% names(df)) {
      df.fs[[treat.name]] <- df[[treat.name]]
    }
    if (!is.null(offset.name) &&
        !offset.name %in% names(df.fs) && offset.name %in% names(df)) {
      df.fs[[offset.name]] <- df[[offset.name]]
    }
  }


  search_overrides <- list(
        Y = Y,
        Event = Event,
        Treat = Treat,
        Z = Z,
        d0.min = d0.min,
        d1.min = d1.min,
        n.min = n.min,
        hr.threshold = if (!is.null(effect_threshold)) effect_threshold else hr.threshold,
        max.minutes = max.minutes,
        minp = minp,
        maxk = maxk,
        parallel_workers = 1,  # Ensure sequential execution if inside parallel context
        details = details       # Optionally suppress details
  )

  # Only pass GLM params when active (same rationale as consistency_overrides)
  if (!is.null(estimator_fn)) {
    search_overrides$estimator_fn <- estimator_fn
    search_overrides$df_analysis <- df
    search_overrides$effect_threshold <- effect_threshold
  }

  # Merge and filter arguments
  search_args <- modifyList(args_call_all, search_overrides)
  # Optionally, filter only valid arguments for subgroup.search
  valid_args <- names(formals(subgroup.search))
  search_args <- search_args[names(search_args) %in% valid_args]

  find.grps <- tryCatch(
    do.call(subgroup.search, search_args),
    error = function(e) {
      warning("Error in subgroup.search: ", e$message)
      return(NULL)
    }
  )
  if (is.null(find.grps) || inherits(find.grps, "try-error")) {
    warning("Subgroup search failed")
    return(list(sg.harm = NULL, args_call_all = args_call_all))
  }


  # ===========================================================================
  # SECTION 7: INITIALIZE OUTPUT VARIABLES
  # ===========================================================================

  sg.harm <- NULL
  df.est_out <- NULL
  df.predict_out <- NULL
  df.test_out <- NULL
  grp.consistency <- NULL

  max_sg_est <- find.grps$max_sg_est
  prop_maxk <- find.grps$prop_max_count

  t.end_all <- proc.time()[3]
  t.min_all <- (t.end_all - t.start_all) / 60

  # ===========================================================================
  # SECTION 8: CHECK FOR VALID SUBGROUPS
  # ===========================================================================

  has_subgroups <- FALSE

  if (!is.null(find.grps) &&
      !inherits(find.grps, "try-error") &&
      !is.null(find.grps$out.found) &&
      !is.null(find.grps$out.found$hr.subgroups)) {

    hr_values <- find.grps$out.found$hr.subgroups$HR
    # For survival, check HR > hr.consistency (natural scale)
    # For GLM, check effect > consistency_threshold (already on correct scale)
    check_threshold <- if (!is.null(consistency_threshold)) {
      consistency_threshold
    } else {
      hr.consistency
    }
    has_subgroups <- any(hr_values > check_threshold, na.rm = TRUE)
  }

  # ===========================================================================
  # SECTION 9: CONSISTENCY EVALUATION (WITH TWO-STAGE OPTION)
  # ===========================================================================

  if (has_subgroups) {
    # Set plotting parameter if needed
    if (plot.sg && is.null(by.risk)) {
      by.risk <- round(max(Y) / 12, 0)
    }

    if (details) {
      n_candidates <- nrow(find.grps$out.found$hr.subgroups)
      cat("# of candidate subgroups (meeting all criteria) =", n_candidates, "\n")
    }

    # -------------------------------------------------------------------------
    # Build arguments for subgroup.consistency
    # -------------------------------------------------------------------------

    # Base override arguments
    consistency_overrides <- list(
      df = df.fs,
      hr.subgroups = find.grps$out.found$hr.subgroups,
      Lsg = find.grps$L,
      confs_labels = confs_labels,
      n.splits = fs.splits,
      stop_Kgroups = max_subgroups_search,
      # NEW: Pass two-stage parameters
      use_twostage = use_twostage,
      twostage_args = twostage_args
    )

    # Only pass GLM closure params when active (non-NULL).
    # This ensures the survival path never sends unknown args to
    # subgroup.consistency -- critical for compatibility with
    # the CRAN-installed version.
    if (!is.null(estimator_fn)) {
      consistency_overrides$estimator_fn <- estimator_fn
      consistency_overrides$consistency_threshold <- consistency_threshold
      # The HR column in hr.subgroups contains log-scale values for GLM
      # outcomes; hr.threshold must be on the same scale so the re-filtering
      # inside subgroup.consistency (line: found.hrs$HR >= hr.threshold) works.
      consistency_overrides$hr.threshold <- effect_threshold
    }

    # Run subgroup consistency analysis with error handling
    sc_filtered_args <- filter_call_args(
      args_call_all,
      subgroup.consistency,
      consistency_overrides
    )

    # Diagnostic: detect zero-length arguments that would cause
    # "argument is of length zero" errors in if() statements
    # Exclude twostage_args and parallel_args: empty list() is valid (use defaults)
    sc_zero_len <- names(sc_filtered_args)[
      vapply(sc_filtered_args, function(x) length(x) == 0 && !is.null(x), logical(1))
    ]
    sc_zero_len <- setdiff(sc_zero_len, c("twostage_args", "parallel_args",
                                          "estimator_fn", "consistency_threshold"))
    if (length(sc_zero_len) > 0) {
      warning("Zero-length arguments passed to subgroup.consistency: ",
              paste(sc_zero_len, collapse = ", "))
      if (details) {
        cat("*** DIAGNOSTIC: Zero-length args for subgroup.consistency:\n")
        for (zarg in sc_zero_len) {
          cat("    ", zarg, " = ", deparse(sc_filtered_args[[zarg]]), "\n")
        }
      }
    }

    grp.consistency <- tryCatch({
      do.call(subgroup.consistency, sc_filtered_args)
    }, error = function(e) {
      warning("Error in subgroup.consistency: ", e$message)
      if (details) {
        cat("*** subgroup.consistency error traceback:\n")
        cat("    Message: ", e$message, "\n")
        cat("    Arg names: ", paste(names(sc_filtered_args), collapse = ", "), "\n")
        # Report lengths of key arguments
        sc_key_args <- c("df", "hr.subgroups", "Lsg", "confs_labels",
                         "sg_focus", "stop_Kgroups", "parallel_args")
        for (ka in intersect(sc_key_args, names(sc_filtered_args))) {
          val <- sc_filtered_args[[ka]]
          if (is.data.frame(val)) {
            cat("    ", ka, ": data.frame [", nrow(val), " x ", ncol(val), "]\n")
          } else {
            cat("    ", ka, ": length=", length(val), " class=", class(val)[1], "\n")
          }
        }
      }
      return(NULL)
    })

    # Handle errors gracefully
    if (is.null(grp.consistency) || inherits(grp.consistency, "try-error")) {
      if (details) {
        cat("Consistency analysis failed - proceeding without results\n")
      }
      grp.consistency <- NULL
    }

    # Update timing
    t.end_all <- proc.time()[3]
    t.min_all <- (t.end_all - t.start_all) / 60

    if (details) {
      cat("Seconds and minutes forestsearch overall =", round(c(t.min_all*60,t.min_all), 4), "\n")
      if (!is.null(grp.consistency$algorithm)) {
        cat("Consistency algorithm used:", grp.consistency$algorithm, "\n")
      }
    }

    # -------------------------------------------------------------------------
    # Process results if consistency analysis succeeded
    # -------------------------------------------------------------------------

    if (!is.null(grp.consistency) && !is.null(grp.consistency$sg.harm)) {
      sg.harm <- grp.consistency$sg.harm

      if (details) {
        cat("Subgroup identified:", paste(sg.harm, collapse = " & "), "\n")
      }

      # Extract prediction datasets
      temp <- grp.consistency$df_flag

      # Merge to analysis data
      df.est_out <- merge(df, temp, by.x = id.name, by.y = "id", all.x = TRUE)

      # Return df.predict
      if (!is.null(df.predict)) {
        df.predict_out <- merge(df.predict, temp, by.x = id.name, by.y = "id", all.x = TRUE)
      }

      # Return df.test
      if (!is.null(df.test)) {
        df.test_out <- get_dfpred(
          df.predict = df.test,
          sg.harm = grp.consistency$sg.harm,
          version = 2
        )
      }
    }
  } # End has_subgroups

  # ===========================================================================
  # SECTION 10: COMPILE AND RETURN OUTPUT
  # ===========================================================================

  out <- list(
    grp.consistency = grp.consistency,
    find.grps = find.grps,
    confounders.candidate = FSconfounders.name,
    confounders.evaluated = confs_labels,
    df.est = df.est_out,
    df.predict = df.predict_out,
    df.test = df.test_out,
    minutes_all = t.min_all,
    grf_res = grf_res,
    sg_focus = sg_focus,
    sg.harm = sg.harm,
    grf_cuts = grf_cuts,
    prop_maxk = prop_maxk,
    max_sg_est = max_sg_est,
    grf_plot = grf_plot,
    args_call_all = args_call_all,
    # NEW: Include algorithm information
    consistency_algorithm = if (!is.null(grp.consistency)) grp.consistency$algorithm else NA,
    # NEW: GLM outcome information
    outcome_type = outcome_type,
    effect_measure = effect_measure
  )

  # Return early if FSdata or find.grps failed
  if (inherits(FSdata, "try-error") || inherits(find.grps, "try-error")) {
    out <- list(sg.harm = NULL, args_call_all = args_call_all)
  }

  class(out) <- c("forestsearch", "list")
  return(out)
}


