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
#' @param cens_type Character string specifying censoring distribution: "weibull"
#'   or "uniform". Default is "weibull"
#' @param cens_params List of parameters for censoring distribution. For uniform:
#'   list(min = value, max = value). For Weibull: fitted from data
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
#' @param select_censoring Logical indicating whether to fit censoring model from
#'   data. If FALSE, censoring is not modeled. Default FALSE
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
#' - \eqn{X} includes treatment, covariates, and treatment×subgroup interaction
#' - \eqn{\sigma} is the scale parameter
#' - \eqn{\epsilon} follows an extreme value distribution
#'
#' ## Interaction Term
#'
#' The model creates a SINGLE interaction term representing the treatment effect
#' modification when ALL subgroup conditions are simultaneously satisfied. This
#' is not multiple separate interactions but one combined indicator.
#'
#' @return An object of class `c("aft_dgm_flex", "list")` containing:
#' \describe{
#'   \item{df_super}{Data frame with the super population including all covariates,
#'     linear predictors, and potential outcomes}
#'   \item{model_params}{List containing model parameters:
#'     \describe{
#'       \item{mu}{Intercept from AFT model}
#'       \item{sigma}{Scale parameter}
#'       \item{gamma}{Vector of regression coefficients}
#'       \item{b_weibull}{Weibull parameterization coefficients}
#'       \item{b0_weibull}{Weibull baseline hazard coefficients}
#'       \item{censoring}{Censoring distribution parameters}
#'     }}
#'   \item{subgroup_info}{List with subgroup information:
#'     \describe{
#'       \item{vars}{Variables used to define subgroup}
#'       \item{cuts}{Cutpoint specifications used}
#'       \item{definitions}{Human-readable subgroup definitions}
#'       \item{size}{Number of observations in subgroup}
#'       \item{proportion}{Proportion of observations in subgroup}
#'     }}
#'   \item{hazard_ratios}{List of true hazard ratios:
#'     \describe{
#'       \item{overall}{Overall treatment HR}
#'       \item{harm_subgroup}{HR within subgroup (if model="alt")}
#'       \item{no_harm_subgroup}{HR outside subgroup (if model="alt")}
#'     }}
#'   \item{analysis_vars}{List of variable classifications for analysis}
#'   \item{model_type}{Character: "alt" or "null"}
#'   \item{n_super}{Size of super population}
#'   \item{seed}{Random seed used}
#' }
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data(cancer)
#'
#' # Example 1: Simple fixed cutpoints
#' dgm1 <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(
#'     er = 20,         # Fixed value
#'     meno = 0         # Factor level
#'   ),
#'   model = "alt",
#'   verbose = TRUE
#' )
#'
#' # Example 2: Quantile-based cutpoints
#' dgm2 <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "pgr", "age"),
#'   subgroup_cuts = list(
#'     er = list(type = "quantile", value = 0.25),
#'     pgr = list(type = "function", fun = median),
#'     age = list(type = "range", min = 40, max = 60)
#'   ),
#'   model = "alt",
#'   k_inter = 2,  # Double the interaction effect
#'   verbose = TRUE
#' )
#'
#' # Print summary
#' print(dgm2)
#' }
#'
#' @seealso
#' \code{\link{simulate_from_dgm}} for generating simulated data from the DGM
#' \code{\link{find_quantile_for_proportion}} for finding quantiles that achieve
#'   target subgroup proportions
#'
#' @references
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
                                  select_censoring = FALSE,
                                  cens_type = "weibull",
                                  cens_params = list(),
                                  seed = 8316951,
                                  verbose = TRUE,
                                  standardize = FALSE,
                                  spline_spec = NULL) {
 if(standardize) message("Standardizing continuous covariates, take care interpreting predictive effects (especially spline)")
  # Set seed for reproducibility
  set.seed(seed)
  if(cens_type %in% c("weibull","uniform")) select_censoring = TRUE
  if(select_censoring){
  if(is.null(factor_vars_cens) & is.null(continuous_vars_cens)) factor_vars_cens = factor_vars
  if(is.null(continuous_vars_cens)) continuous_vars_cens = continuous_vars
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
  treatment_var, continuous_vars, factor_vars, standardize,
  continuous_vars_cens, factor_vars_cens, verbose)

  # ============================================================================
  # Step 3: Define Subgroups with Flexible Cutpoints
  # ============================================================================
  subgroup_result <- define_subgroups(df_work, data, subgroup_vars,
                                      subgroup_cuts, continuous_vars,
                                      model, verbose)
  df_work$flag_harm <- subgroup_result$flag_harm
  interaction_term <- subgroup_result$interaction_term
  subgroup_definitions <- subgroup_result$definitions

  # ============================================================================
  # Step 4: Fit AFT Model (Weibull) - with optional spline
  # ============================================================================
  set_var <- set_beta_spec$set_var
  beta_var <- set_beta_spec$beta_var

  aft_params <- fit_aft_model(df_work, interaction_term, k_treat,
                              k_inter, verbose, spline_spec, set_var, beta_var)
  mu <- aft_params$mu
  tau <- aft_params$tau
  gamma <- aft_params$gamma
  b0 <- aft_params$b0
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
                                        spline_info)  # PASS spline_info
  # ============================================================================
  # Step 6: Calculate Hazard Ratios
  # ============================================================================
  hr_results <- calculate_hazard_ratios(df_super, n_super, mu, tau,
                                        model, verbose)

  # ============================================================================
  # Step 7: Prepare Censoring Parameters
  # ============================================================================
  cens_model <- NULL

  if(select_censoring){
  cens_result <- prepare_censoring_model(df_work, cens_type, cens_params,
                                         df_super)
  cens_model <- cens_result$cens_model
  df_super <- cens_result$df_super
  }

  # ============================================================================
  # Step 8: Assemble and Return Results
  # ============================================================================
  results <- assemble_results(
    df_super = df_super,
    mu = mu,
    tau = tau,
    gamma = gamma,
    b0 = b0,
    cens_model = cens_model,
    subgroup_vars = subgroup_vars,
    subgroup_cuts = subgroup_cuts,
    subgroup_definitions = subgroup_definitions,
    hr_results = hr_results,
    continuous_vars = continuous_vars,
    factor_vars = factor_vars,
    model = model,
    n_super = n_super,
    seed = seed,
    spline_info = spline_info  # This is now always passed (can be NULL)
  )

  return(results)
}

