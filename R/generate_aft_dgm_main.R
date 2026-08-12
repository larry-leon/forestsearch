# Enhanced AFT Data Generating Mechanism with Flexible Subgroup Definitions
# Refactored version with improved modularity
# Supports multiple cutpoint specifications: fixed values, quantiles, functions, etc.

#' Generate Synthetic Survival Data using AFT Model with Flexible Subgroups
#'
#' Creates a data generating mechanism (DGM) for survival data using an Accelerated
#' Failure Time (AFT) model with Weibull distribution. Supports flexible subgroup
#' definitions and treatment-subgroup interactions.
#'
#' @param data A data.frame containing the input dataset to base the simulation on
#' @param continuous_vars Character vector of continuous variable names to be
#'   standardized and included as covariates
#' @param factor_vars Character vector of factor/categorical variable names to be
#'   converted to dummy variables (largest value as reference)
#' @param outcome_var Character string specifying the name of the outcome/time variable
#' @param event_var Character string specifying the name of the event/status variable
#'   (1 = event, 0 = censored)
#' @param treatment_var Character string specifying the name of the treatment variable.
#'   If NULL, treatment will be randomly simulated with 50/50 allocation
#' @param subgroup_vars Character vector of variable names defining the subgroup.
#'   Default is NULL (no subgroups)
#' @param subgroup_cuts Named list of cutpoint specifications for subgroup variables.
#'   See Details section for flexible specification options
#' @param draw_treatment Logical indicating whether to redraw treatment assignment
#'   in simulation. Default is FALSE (use original assignments)
#' @param model Character string: "alt" for alternative model with subgroup effects,
#'   "null" for null model without subgroup effects. Default is "alt"
#' @param k_treat Numeric treatment effect modifier. Values >1 increase treatment
#'   effect, <1 decrease it. Default is 1 (no modification)
#' @param k_inter Numeric interaction effect modifier for treatment-subgroup interaction.
#'   Default is 1 (no modification)
#' @param n_super Integer specifying size of super population to generate.
#'   Default is 5000
#' @param select_censoring Logical. If \code{TRUE} (default), fits the censoring
#'   distribution to the observed censoring times in \code{data} using
#'   \code{survreg} with AIC-based selection among Weibull and log-normal models
#'   (with and without covariates). If \code{FALSE}, no model is fitted; the
#'   censoring distribution is specified entirely by \code{cens_params}. Default
#'   \code{TRUE}.
#' @param cens_type Character string specifying censoring distribution type:
#'   \code{"weibull"} or \code{"uniform"}. Controls which parametric family is
#'   considered when \code{select_censoring = TRUE}, and determines the required
#'   structure of \code{cens_params} when \code{select_censoring = FALSE}.
#'   Default \code{"weibull"}.
#' @param cens_params Named list of censoring distribution parameters.
#'   Interpretation depends on \code{select_censoring} and \code{cens_type}:
#'   \describe{
#'     \item{\code{select_censoring = TRUE}}{Ignored; all parameters are estimated
#'       from data.}
#'     \item{\code{select_censoring = FALSE, cens_type = "uniform"}}{Must supply
#'       \code{min} and \code{max}. If either is absent, defaults to
#'       \code{0.5 * min(y)} and \code{1.5 * max(y)} with a message.}
#'     \item{\code{select_censoring = FALSE, cens_type = "weibull"}}{Must supply
#'       \code{mu} (log-scale location) and \code{tau} (scale). Optionally supply
#'       \code{type} (\code{"weibull"} or \code{"lognormal"}); defaults to
#'       \code{"weibull"}. Censoring is treated as intercept-only (no covariate
#'       or treatment dependence): \code{lin_pred_cens_0 = lin_pred_cens_1 = mu}.}
#'   }
#'   Default \code{list()}.
#' @param cens_intercept_only Logical. Honored only in force-fit mode
#'   (\code{select_censoring = FALSE} with empty \code{cens_params}).  If
#'   \code{TRUE}, fits \code{Surv(y, 1 - event) ~ 1} with no treat or
#'   covariate dependence.  If \code{FALSE} (default), the censoring
#'   formula is determined by the censoring covariate vectors:
#'   \itemize{
#'     \item Both \code{*_cens = NULL} (default): inherit from outcome
#'           model; formula is \code{~ treat + zcens_<outcome covariates>}.
#'     \item Both \code{*_cens = character(0)}: treat-only; formula is
#'           \code{~ treat}.
#'     \item One or both \code{*_cens} non-empty: formula is
#'           \code{~ treat + zcens_<supplied covariates>}.
#'   }
#'   Setting \code{cens_intercept_only = TRUE} with
#'   \code{select_censoring = TRUE} or with analytical \code{cens_params}
#'   (mu, tau supplied) raises an error.
#' @param seed Integer random seed for reproducibility. Default is 8316951
#' @param verbose Logical indicating whether to print diagnostic information during
#'   execution. Default is TRUE
#' @param standardize Logical indicating whether to standardize continuous variables.
#'   Default is FALSE
#' @param continuous_vars_cens Character vector of continuous variable names to be
#'   used for censoring model. If NULL, uses same as continuous_vars. Default NULL
#' @param factor_vars_cens Character vector of factor variable names to be used
#'   for censoring model. If NULL, uses same as factor_vars. Default NULL
#' @param set_beta_spec List with elements 'set_var' and 'beta_var' for manually
#'   setting specific beta coefficients. Default list(set_var = NULL, beta_var = NULL)
#' @param spline_spec List specifying spline configuration for treatment effect.
#'   Must include 'var' (variable name), 'knot', 'zeta', and 'log_hrs' (vector of
#'   length 3). Default NULL (no spline)
#'
#' @details
#' ## Subgroup Cutpoint Specifications
#'
#' The `subgroup_cuts` parameter accepts multiple flexible specifications:
#'
#' ### Fixed Value
#' ```r
#' subgroup_cuts = list(er = 20)  # er <= 20
#' ```
#'
#' ### Quantile-based
#' ```r
#' subgroup_cuts = list(
#'   er = list(type = "quantile", value = 0.25)  # er <= 25th percentile
#' )
#' ```
#'
#' ### Function-based
#' ```r
#' subgroup_cuts = list(
#'   er = list(type = "function", fun = median)  # er <= median
#' )
#' ```
#'
#' ### Range
#' ```r
#' subgroup_cuts = list(
#'   age = list(type = "range", min = 40, max = 60)  # 40 <= age <= 60
#' )
#' ```
#'
#' ### Greater than
#' ```r
#' subgroup_cuts = list(
#'   nodes = list(type = "greater", quantile = 0.75)  # nodes > 75th percentile
#' )
#' ```
#'
#' ### Multiple values (for categorical)
#' ```r
#' subgroup_cuts = list(
#'   grade = list(type = "multiple", values = c(2, 3))  # grade in (2, 3)
#' )
#' ```
#'
#' ### Custom function
#' ```r
#' subgroup_cuts = list(
#'   er = list(
#'     type = "custom",
#'     fun = function(x) x <= quantile(x, 0.3) | x >= quantile(x, 0.9)
#'   )
#' )
#' ```
#'
#' ## Model Structure
#'
#' The AFT model with Weibull distribution is specified as:
#' \deqn{\log(T) = \mu + \gamma' X + \sigma \epsilon}
#'
#' Where:
#' - \eqn{T} is the survival time
#' - \eqn{\mu} is the intercept
#' - \eqn{\gamma} contains the covariate effects
#' - \eqn{X} includes treatment, covariates, and treatment x subgroup interaction
#' - \eqn{\sigma} is the scale parameter
#' - \eqn{\epsilon} follows an extreme value distribution
#'
#' ## Interaction Term
#'
#' The model creates a SINGLE interaction term representing the treatment effect
#' modification when ALL subgroup conditions are simultaneously satisfied. This
#' is not multiple separate interactions but one combined indicator.
#'
#' @references
#' Leon, L.F., et al. (2024). Statistics in Medicine.
#'
#' Kalbfleisch, J.D. and Prentice, R.L. (2002). The Statistical Analysis of
#'   Failure Time Data (2nd ed.). Wiley.
#'
#' @author Your Name
#' @examples
#' \dontrun{
#' df <- survival::gbsg
#' dgm <- generate_aft_dgm_flex(
#'   data            = df,
#'   outcome_var     = "rfstime",
#'   event_var       = "status",
#'   treatment_var   = "hormon",
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars     = "meno",
#'   model           = "null",
#'   verbose         = FALSE
#' )
#' str(dgm)
#' }
#' @export
#' @importFrom survival survreg coxph Surv
#' @importFrom stats quantile median uniroot rexp runif rnorm rbinom model.matrix coef reformulate
#'
#' @return An object of class \code{"aft_dgm_flex"} (a list) with components:
#'   \describe{
#'     \item{\code{df_super}}{Data frame containing the super-population,
#'       including covariates, treatment, counterfactual linear predictors,
#'       and subgroup indicator \code{flag_harm}.}
#'     \item{\code{df_source}}{Data frame with every source patient exactly
#'       once (the prepared trial data), carrying the same counterfactual
#'       outcome and censoring linear predictors as \code{df_super}.
#'       Consumed by \code{simulate_from_dgm(baseline = "fixed")} to hold
#'       baseline covariates -- and hence subgroup composition -- fixed
#'       across simulated trials.}
#'     \item{\code{model_params}}{List with AFT parameters \code{mu},
#'       \code{tau}, \code{gamma}, \code{b0}, the fitted \code{censoring}
#'       model, and optional \code{spline_info}.}
#'     \item{\code{subgroup_info}}{List describing the true subgroup:
#'       \code{vars}, \code{cuts}, \code{definitions}, \code{size},
#'       \code{proportion}.}
#'     \item{\code{hazard_ratios}}{List of true HR/AHR/CDE values on the
#'       super-population (see \code{\link{compute_dgm_cde}}).}
#'     \item{\code{analysis_vars}}{Named list of column roles
#'       (continuous, factor, covariates, treatment, outcome, event).}
#'     \item{\code{model_type}}{\code{"null"} or \code{"alt"}.}
#'     \item{\code{n_super}}{Super-population size.}
#'     \item{\code{seed}}{Seed used for super-population generation.}
#'   }

generate_aft_dgm_flex <- function(data,
                                  continuous_vars,
                                  factor_vars,
                                  continuous_vars_cens = NULL,
                                  factor_vars_cens = NULL,
                                  set_beta_spec = list(set_var = NULL, beta_var = NULL),
                                  outcome_var,
                                  event_var,
                                  treatment_var = NULL,
                                  subgroup_vars = NULL,
                                  subgroup_cuts = NULL,
                                  draw_treatment = FALSE,
                                  model = "alt",
                                  k_treat = 1,
                                  k_inter = 1,
                                  n_super = 5000,
                                  select_censoring = TRUE,
                                  cens_type = "weibull",
                                  cens_params = list(),
                                  cens_intercept_only = FALSE,
                                  seed = 8316951,
                                  verbose = TRUE,
                                  standardize = FALSE,
                                  spline_spec = NULL) {

  if (standardize)
    message("Standardizing continuous covariates, take care interpreting ",
            "predictive effects (especially spline)")

  # Set seed for reproducibility
  set.seed(seed)

  # When the censoring model will be fitted from data, ensure censoring
  # covariate vectors default to the outcome model vectors if not supplied.
  # The censoring model is fitted in two cases:
  #   (i)  AIC selection                : select_censoring = TRUE
  #   (ii) Force-fit                    : select_censoring = FALSE AND
  #                                       cens_params lacks analytical specs
  # The analytical path (FALSE + cens_params with mu/tau) does not fit
  # anything; covariate processing is irrelevant there and is skipped.
  #
  # cens_intercept_only = TRUE forcibly clears any covariate vectors so
  # the censoring formula reduces to ~ 1.  This must run *before* the
  # inheritance step.
  if (cens_intercept_only) {
    if (select_censoring) {
      stop("cens_intercept_only = TRUE requires select_censoring = FALSE ",
           "(force-fit mode). For an unconditional model in AIC mode, the ",
           "comparison already includes Weibull0 / LogNormal0 candidates.",
           call. = FALSE)
    }
    if (!is.null(cens_params$mu) || !is.null(cens_params$tau)) {
      stop("cens_intercept_only = TRUE is incompatible with analytical ",
           "cens_params (mu/tau supplied).  Choose one: force-fit ",
           "intercept-only (cens_params empty), or analytical with ",
           "explicit mu/tau (cens_intercept_only = FALSE).",
           call. = FALSE)
    }
    factor_vars_cens     <- character(0)
    continuous_vars_cens <- character(0)
  }

  # Inheritance (default *_cens to outcome covariates) when fitting.
  will_fit_cens <- select_censoring ||
    (cens_type != "uniform" &&
     is.null(cens_params$mu) && is.null(cens_params$tau))

  if (will_fit_cens && !cens_intercept_only) {
    if (is.null(factor_vars_cens) && is.null(continuous_vars_cens))
      factor_vars_cens <- factor_vars
    if (is.null(continuous_vars_cens))
      continuous_vars_cens <- continuous_vars
  }

  # ============================================================================
  # Step 1: Input Validation
  # ============================================================================
  validate_inputs(data, model, cens_type, outcome_var, event_var,
                  treatment_var, continuous_vars, factor_vars)

  # ============================================================================
  # Step 2: Data Preparation
  # ============================================================================
  df_work <- prepare_working_dataset(data, outcome_var, event_var,
                                     treatment_var, continuous_vars, factor_vars,
                                     standardize, continuous_vars_cens,
                                     factor_vars_cens, verbose)

  # ============================================================================
  # Step 3: Define Subgroups with Flexible Cutpoints
  # ============================================================================
  subgroup_result <- define_subgroups(df_work, data, subgroup_vars,
                                      subgroup_cuts, continuous_vars,
                                      model, verbose)
  df_work$flag_harm       <- subgroup_result$flag_harm
  interaction_term        <- subgroup_result$interaction_term
  subgroup_definitions    <- subgroup_result$definitions

  # ============================================================================
  # Step 4: Fit AFT Model (Weibull) - with optional spline
  # ============================================================================
  set_var  <- set_beta_spec$set_var
  beta_var <- set_beta_spec$beta_var

  aft_params  <- fit_aft_model(df_work, interaction_term, k_treat,
                               k_inter, verbose, spline_spec, set_var, beta_var)
  mu          <- aft_params$mu
  tau         <- aft_params$tau
  gamma       <- aft_params$gamma
  b0          <- aft_params$b0
  spline_info <- aft_params$spline_info

  # If spline was used, update df_work with spline variables
  if (!is.null(spline_info)) {
    df_work <- spline_info$df_work
  }

  # ============================================================================
  # Step 5: Generate Super Population
  # ============================================================================
  df_super <- generate_super_population(df_work, n_super, draw_treatment,
                                        gamma, b0, mu, tau, verbose,
                                        spline_info)

  # ============================================================================
  # Step 6: Calculate Hazard Ratios
  # ============================================================================
  hr_results <- calculate_hazard_ratios(df_super, n_super, mu, tau,
                                        model, verbose)

  # ============================================================================
  # Step 7: Prepare Censoring Parameters
  # ============================================================================
  # Three modes (see prepare_censoring_model docstring for details):
  #   (1) AIC selection : select_censoring = TRUE
  #   (2) Analytical    : select_censoring = FALSE; cens_params supplies
  #                       mu+tau (Weibull/LogNormal) or min+max (uniform)
  #   (3) Force-fit     : select_censoring = FALSE; cens_params empty;
  #                       formula determined by cens_intercept_only and
  #                       the *_cens vectors
  cens_result <- prepare_censoring_model(
    df_work             = df_work,
    cens_type           = cens_type,
    cens_params         = cens_params,
    df_super            = df_super,
    select_censoring    = select_censoring,
    cens_intercept_only = cens_intercept_only,
    verbose             = verbose
  )
  cens_model <- cens_result$cens_model
  df_super   <- cens_result$df_super

  # ============================================================================
  # Step 7b: Fixed-Baseline Source Frame (df_source)
  # ============================================================================
  # df_work holds every source patient exactly once (post covariate prep,
  # flag_harm, and any spline columns).  Attaching the same potential-outcome
  # and censoring linear predictors used for df_super makes it directly
  # consumable by simulate_from_dgm(baseline = "fixed"), which holds baseline
  # covariates -- and hence subgroup memberships and sizes -- fixed at the
  # observed trial composition while re-drawing outcomes each simulation.
  #
  # The outcome linear predictors are computed with the same
  # calculate_linear_predictors() call used for df_super.  The censoring
  # linear predictors are attached via a second prepare_censoring_model()
  # call: the refit is deterministic on identical df_work (AIC selection
  # included), so lin_pred_cens_* are exactly consistent with cens_model;
  # its returned model object is discarded.  This trades one redundant
  # survreg fit at DGM-build time for zero duplication of the per-mode
  # (weibull / lognormal / uniform / analytical) branch logic.
  df_source    <- df_work
  df_source$id <- seq_len(nrow(df_source))
  covariate_cols_src <- grep("^z_", names(df_source), value = TRUE)
  df_source <- calculate_linear_predictors(df_source, covariate_cols_src,
                                           gamma, b0, spline_info)
  cens_result_src <- prepare_censoring_model(
    df_work             = df_work,
    cens_type           = cens_type,
    cens_params         = cens_params,
    df_super            = df_source,      # attaches lin_pred_cens_0/1
    select_censoring    = select_censoring,
    cens_intercept_only = cens_intercept_only,
    verbose             = FALSE           # already reported for df_super
  )
  df_source <- cens_result_src$df_super

  # ============================================================================
  # Step 8: Assemble and Return Results
  # ============================================================================
  results <- assemble_results(
    df_super             = df_super,
    mu                   = mu,
    tau                  = tau,
    gamma                = gamma,
    b0                   = b0,
    cens_model           = cens_model,
    subgroup_vars        = subgroup_vars,
    subgroup_cuts        = subgroup_cuts,
    subgroup_definitions = subgroup_definitions,
    hr_results           = hr_results,
    continuous_vars      = continuous_vars,
    factor_vars          = factor_vars,
    model                = model,
    n_super              = n_super,
    seed                 = seed,
    spline_info          = spline_info,
    df_source            = df_source
  )

  return(results)
}
