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
#' @export
#' @importFrom survival survreg coxph Surv
#' @importFrom stats quantile median uniroot rexp runif rnorm rbinom model.matrix coef

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
  # Skipped when select_censoring = FALSE: no covariate processing is needed
  # for a censoring distribution that will be specified analytically.
  if (select_censoring) {
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
  # select_censoring = TRUE  : fits censoring distribution from observed data
  #                            via AIC comparison across Weibull / lognormal
  #                            models (with and without covariates)
  # select_censoring = FALSE : uses caller-supplied cens_params directly;
  #                            no model fitting is performed
  cens_result <- prepare_censoring_model(
    df_work          = df_work,
    cens_type        = cens_type,
    cens_params      = cens_params,
    df_super         = df_super,
    select_censoring = select_censoring,
    verbose          = verbose
  )
  cens_model <- cens_result$cens_model
  df_super   <- cens_result$df_super

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
    spline_info          = spline_info
  )

  return(results)
}
