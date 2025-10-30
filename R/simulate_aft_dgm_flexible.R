# Enhanced AFT Data Generating Mechanism with Flexible Subgroup Definitions
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
#' - T is the survival time
#' - μ is the intercept
#' - γ contains the covariate effects
#' - X includes treatment, covariates, and treatment×subgroup interaction
#' - σ is the scale parameter
#' - ε follows an extreme value distribution
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
                                 cens_type = "weibull",
                                 cens_params = list(),
                                 seed = 8316951,
                                 verbose = TRUE, standardize = FALSE) {

  # ============================================================================
  # Helper function to process cutpoint specifications
  # ============================================================================

  process_cutpoint <- function(var_data, cut_spec, var_name = "") {
    # If cut_spec is a simple numeric value, treat as fixed cutpoint
    if (is.numeric(cut_spec) && length(cut_spec) == 1) {
      return(var_data <= cut_spec)
    }

    # If it's a list, process based on type
    if (is.list(cut_spec)) {
      cut_type <- cut_spec$type

      if (cut_type == "quantile") {
        # Quantile-based cutpoint
        cutpoint <- quantile(var_data, probs = cut_spec$value, na.rm = TRUE)
        return(var_data <= cutpoint)

      } else if (cut_type == "function") {
        # Function-based cutpoint (e.g., median, mean)
        cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
        return(var_data <= cutpoint)

      } else if (cut_type == "range") {
        # Range-based (between min and max)
        return(var_data >= cut_spec$min & var_data <= cut_spec$max)

      } else if (cut_type == "greater") {
        # Greater than cutpoint
        if (!is.null(cut_spec$value)) {
          cutpoint <- cut_spec$value
        } else if (!is.null(cut_spec$quantile)) {
          cutpoint <- quantile(var_data, probs = cut_spec$quantile, na.rm = TRUE)
        } else if (!is.null(cut_spec$fun)) {
          cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
        }
        return(var_data > cutpoint)

      } else if (cut_type == "multiple") {
        # Multiple cutpoints (var in specified values)
        return(var_data %in% cut_spec$values)

      } else if (cut_type == "custom") {
        # Custom function that returns logical vector
        return(cut_spec$fun(var_data))

      } else {
        stop(paste("Unknown cutpoint type:", cut_type, "for variable:", var_name))
      }
    }

    # Default: if no specification, use median
    if (is.null(cut_spec)) {
      if (verbose) cat("  Using median as default cutpoint for", var_name, "\n")
      return(var_data <= median(var_data, na.rm = TRUE))
    }

    stop(paste("Invalid cutpoint specification for variable:", var_name))
  }

  # ============================================================================
  # Input Validation
  # ============================================================================

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if (!model %in% c("alt", "null")) {
    stop("'model' must be either 'alt' or 'null'")
  }

  if (!cens_type %in% c("weibull", "uniform")) {
    stop("'cens_type' must be either 'weibull' or 'uniform'")
  }

  # Check that required variables exist
  required_vars <- c(outcome_var, event_var, treatment_var)

  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Required variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Check continuous and factor variables
  all_covars <- c(continuous_vars, factor_vars)
  missing_covars <- setdiff(all_covars, names(data))
  if (length(missing_covars) > 0) {
    stop("Covariate variables not found in data: ", paste(missing_covars, collapse = ", "))
  }

  # ============================================================================
  # Data Preparation
  # ============================================================================

  set.seed(seed)

  # Create working dataset
  df_work <- data.frame(
    id = 1:nrow(data),
    y = data[[outcome_var]],
    treat = data[[treatment_var]],
    event = ifelse(data[[event_var]] == 1, 1, 0)
  )


  # ============================================================================
  # Process Covariates
  # ============================================================================

  # Process continuous variables
  for (var in continuous_vars) {
 if(standardize){
   df_work[[paste0("z_", var)]] <- scale(data[[var]])[, 1]  # Standardize
 } else {
   df_work[[paste0("z_", var)]] <- data[[var]]
 }
   }

  # Process factor variables with LARGEST value as reference (UNLESS already binary, then retain)
  for (var in factor_vars) {
    # Get unique values
    if (is.factor(data[[var]])) {
      all_levels <- levels(data[[var]])
    } else {
      all_levels <- sort(unique(data[[var]]))
    }

    n_levels <- length(all_levels)

    if (n_levels == 1) {
      # Skip variables with only one level
      next

    } else if (n_levels == 2) {
      # NOTE: revising to keep as-is (switching max/min roles)
      # Binary variable
      if (is.numeric(all_levels)) {
        ref_level <- min(all_levels)
        other_level <- max(all_levels)
      } else {
        # For character, last alphabetically
        ref_level <- sort(all_levels, decreasing = FALSE)[1]
        other_level <- setdiff(all_levels, ref_level)
      }

      # Create single indicator for non-reference level
      df_work[[paste0("z_", var)]] <- as.numeric(data[[var]] == other_level)

    } else {
      # Multiple levels - create dummies with largest as reference
      if (is.numeric(all_levels)) {
        ref_level <- max(all_levels)
      } else {
        ref_level <- sort(all_levels, decreasing = TRUE)[1]
      }

      other_levels <- setdiff(all_levels, ref_level)
      other_levels <- sort(other_levels)

      # Create dummy for each non-reference level
      for (level in other_levels) {
        # Use the actual level value in the variable name for clarity
        dummy_name <- paste0("z_", var, "_", level)
        df_work[[dummy_name]] <- as.numeric(data[[var]] == level)
      }
    }
  }

  # Add any baseline variables in original data that are not specified in the DGM outcome process
  # Identify processed variables
  processed <- c(continuous_vars, factor_vars, outcome_var, event_var, treatment_var)

  # Find unprocessed variables
  unprocessed <- setdiff(names(data), processed)

  if (length(unprocessed) > 0) {
    # Add them to df_work
    for (var in unprocessed) {
      if (!var %in% names(df_work)) {  # Avoid duplicates
        df_work[[var]] <- data[[var]]
      }
    }

    if (verbose) {
      cat("\nAdded", length(unprocessed), "unprocessed variables:",
          paste(unprocessed, collapse = ", "), "\n")
    }
  }

  # ============================================================================
  # Define Subgroups with Flexible Cutpoints
  # ============================================================================

  df_work$flag_harm <- 0
  interaction_term <- NULL
  subgroup_definitions <- list()

  if (model == "alt" && !is.null(subgroup_vars)) {
    # Create subgroup indicators
    subgroup_indicators <- list()

    if (verbose) {
      cat("\n=== Subgroup Definitions ===\n")
    }

    for (var in subgroup_vars) {
      # Get the cutpoint specification for this variable
      cut_spec <- subgroup_cuts[[var]]

      if (var %in% continuous_vars || is.numeric(data[[var]])) {
        # Process continuous variable with flexible cutpoint
        subgroup_indicators[[var]] <- process_cutpoint(
          var_data = data[[var]],
          cut_spec = cut_spec,
          var_name = var
        )

        # Store the actual cutpoint used for reporting
        if (is.numeric(cut_spec) && length(cut_spec) == 1) {
          actual_cutpoint <- cut_spec
          subgroup_definitions[[var]] <- paste(var, "<=", actual_cutpoint)
        } else if (is.list(cut_spec)) {
          if (cut_spec$type == "quantile") {
            actual_cutpoint <- quantile(data[[var]], probs = cut_spec$value, na.rm = TRUE)
            subgroup_definitions[[var]] <- paste(var, "<=", round(actual_cutpoint, 2),
                                                "(", cut_spec$value*100, "th percentile)")
          } else if (cut_spec$type == "function") {
            actual_cutpoint <- cut_spec$fun(data[[var]], na.rm = TRUE)
            fun_name <- deparse(substitute(cut_spec$fun))
            subgroup_definitions[[var]] <- paste(var, "<=", round(actual_cutpoint, 2),
                                                "(", fun_name, ")")
          } else if (cut_spec$type == "range") {
            subgroup_definitions[[var]] <- paste(cut_spec$min, "<=", var, "<=", cut_spec$max)
          } else if (cut_spec$type == "greater") {
            if (!is.null(cut_spec$value)) {
              actual_cutpoint <- cut_spec$value
            } else if (!is.null(cut_spec$quantile)) {
              actual_cutpoint <- quantile(data[[var]], probs = cut_spec$quantile, na.rm = TRUE)
            }
            subgroup_definitions[[var]] <- paste(var, ">", round(actual_cutpoint, 2))
          } else if (cut_spec$type == "custom") {
            subgroup_definitions[[var]] <- paste(var, "(custom function)")
          }
        }

      } else {
        # For factor variables
        if (!is.null(cut_spec)) {
          if (is.character(cut_spec) || is.factor(cut_spec)) {
            subgroup_indicators[[var]] <- data[[var]] == cut_spec
            subgroup_definitions[[var]] <- paste(var, "==", cut_spec)
          } else if (is.list(cut_spec) && cut_spec$type == "multiple") {
            subgroup_indicators[[var]] <- data[[var]] %in% cut_spec$values
            subgroup_definitions[[var]] <- paste(var, "in",
                                                paste(cut_spec$values, collapse = ", "))
          }
        } else {
          # Default: use first level
          first_level <- levels(as.factor(data[[var]]))[1]
          subgroup_indicators[[var]] <- data[[var]] == first_level
          subgroup_definitions[[var]] <- paste(var, "==", first_level)
        }
      }

      if (verbose) {
        cat("  ", subgroup_definitions[[var]], "\n")
        cat("    Proportion in subgroup:",
            round(mean(subgroup_indicators[[var]], na.rm = TRUE), 3), "\n")
      }
    }

    # Create harm flag (all subgroup conditions met)
    df_work$flag_harm <- as.numeric(Reduce("&", subgroup_indicators))

    # Create interaction term
    if (length(subgroup_indicators) > 0) {
      interaction_term <- df_work$treat * df_work$flag_harm
    }

    if (verbose) {
      cat("\nOverall subgroup (all conditions met):\n")
      cat("  Size:", sum(df_work$flag_harm), "out of", nrow(df_work), "\n")
      cat("  Proportion:", round(mean(df_work$flag_harm), 3), "\n")
    }
  }

  # ============================================================================
  # Fit AFT Model (Weibull)
  # ============================================================================

  # Prepare model matrix
  covariate_cols <- grep("^z_", names(df_work), value = TRUE)
  X <- as.matrix(df_work[, c("treat", covariate_cols)])

  # Add interaction term if needed
  if (!is.null(interaction_term)) {
    X <- cbind(X, treat_harm = interaction_term)
  }

  # Fit Weibull AFT model
  formula_str <- paste("Surv(y, event) ~ ", paste(c("treat", covariate_cols), collapse = " + "))
  if (!is.null(interaction_term)) {
    df_work$treat_harm <- interaction_term
    formula_str <- paste(formula_str, "+ treat_harm")
  }

  fit_aft <- survreg(as.formula(formula_str),
                     data = df_work,
                     dist = "weibull")

  # Extract parameters
  mu <- coef(fit_aft)[1]  # Intercept
  tau <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]  # Coefficients (excluding intercept)

  # Apply effect modifiers AFT parameterization
  # gamma["treat"] <- k_treat * gamma["treat"]
  # if ("treat_harm" %in% names(gamma)) {
  #   gamma["treat_harm"] <- k_inter * gamma["treat_harm"]
  # }

  # Weibull parameterization
  b0 <- -gamma / tau

  # Apply effect modifiers Weibull log(hazard-ratio) parameterization
  b0["treat"] <- k_treat * b0["treat"]
  if ("treat_harm" %in% names(b0)) {
    b0["treat_harm"] <- k_inter * b0["treat_harm"]
  }

# Transform to corresponding revised gamma

    gamma <- -b0 * tau

  # ============================================================================
  # Generate Super Population
  # ============================================================================

  # Sample with replacement to create super population

  # Sample indices
  idx_sample <- sample(1:nrow(df_work), size = n_super, replace = TRUE)
  df_super <- df_work[idx_sample, ]


  if (verbose) {
    cat("\nOverall subgroup in super-population (all conditions met):\n")
    cat("  Proportion:", round(mean(df_super$flag_harm), 3), "\n")
  }

  if (verbose) {
    cat("\n=== Model Parameters (AFT, log(T)) ===\n")
    cat("Intercept (mu):", round(mu, 3), "\n")
    cat("Scale (tau):", round(tau, 3), "\n")
    cat("Treatment effect:", round(gamma["treat"], 3), "\n")
    if ("treat_harm" %in% names(gamma)) {
      cat("Interaction effect:", round(gamma["treat_harm"], 3), "\n")
    }
  }

  # Now, in the super population we can either retain the original treatment
  # pattern from the observed "data" or re-draw
  # Default is to leave as-is to retain original randomization pattern/scheme
  if(draw_treatment){
  # Assign treatment simple 1/2
  # May extend to allow for block randomization and unequal
  n_treat <- round(n_super / 2)
  n_control <- n_super - n_treat
  df_super$treat[1:n_treat] <- 1
  df_super$treat[(n_treat + 1):n_super] <- 0
  }

  # Reset IDs
  df_super$id <- 1:n_super

  # Calculate linear predictors for potential outcomes

  X_super <- as.matrix(df_super[, c("treat", covariate_cols)])
  if (!is.null(interaction_term)) {
    # Recalculate interaction for super population
    df_super$treat_harm <- df_super$treat * df_super$flag_harm
    X_super <- cbind(X_super, treat_harm = df_super$treat_harm)
  }
  # Potential outcomes under treatment
  X_treat <- X_super
  X_treat[, "treat"] <- 1
  if ("treat_harm" %in% colnames(X_treat)) {
    X_treat[, "treat_harm"] <- df_super$flag_harm
  }

  # Potential outcomes under control
  X_control <- X_super
  X_control[, "treat"] <- 0
  if ("treat_harm" %in% colnames(X_control)) {
    X_control[, "treat_harm"] <- 0
  }

  # Linear predictors per AFT
  df_super$lin_pred_1 <- X_treat %*% gamma
  df_super$lin_pred_0 <- X_control %*% gamma
  df_super$lin_pred_obs <- X_super %*% gamma

  ## PO log-hazards excluding baseline
  #  Used for calculating empirical version of CDEs (controlled direct effects)
  #  theta0 = L0'beta
  #  theta1 - L1'beta
  # If additional confounder "w" is included
  # theta0 <- -c(zmat.0%*%gamma.true + w*gamma_w)/tau.strataO
  # Here "no special w":

  df_super$theta_0 <- X_control %*% b0
  df_super$theta_1 <- X_treat %*% b0

  df_super$loghr_po <- with(df_super, theta_1 - theta_0)

  # Notes:
  # linear predictor (setting treat=1)
  # z.1 <- zmat.1[,"z"]
  # # log(Y) AFT parameterization
  # eta1 <- mu + c(zmat.1%*%gamma.true)+w*gamma_w
  #  Weibull log(hazard ratio) setting treat=1
  # # and excluding mu and w*bw (since taking difference below)
  # phi1 <- (-1)*c(zmat.1%*%gamma.true)/tau.strataO
  # log.Y1 <- eta1 + tau.strataO*epsilon
  #
  #  Setting treat=0
  # z.0 <- zmat.0[,"z"]
  # eta0 <- mu + c(zmat.0%*%gamma.true)+w*gamma_w
  # log.Y0 <- eta0 + tau.strataO*epsilon
  # phi0 <- (-1)*c(zmat.0%*%gamma.true)/tau.strataO
  #
  # Potential outcome log(hr) difference
  #loghr.po <- phi1-phi0
  #ahr_causal <- with(dfs,exp(mean(loghr.po)))
  #res$AHR <- ahr_causal
  # "Empirical causal" in sense that we are averaging across any covariates in dgm

  # Generate potential outcomes for HR calculation
  epsilon <- log(rexp(n_super))  # Extreme value distribution

  # Under treatment
  logT_1 <- mu + tau * epsilon + df_super$lin_pred_1
  T_1 <- exp(logT_1)

  # Under control
  logT_0 <- mu + tau * epsilon + df_super$lin_pred_0
  T_0 <- exp(logT_0)

  # Calculate empirical hazard ratios
  # based on potential outcomes T_1 = T(treat==1) as-if all receive treatment;
  # T_0 = T(treat=0) under control; and no censoring

    df_temp <- data.frame(
    time = c(T_1, T_0),
    event = 1,
    treat = c(rep(1, n_super), rep(0, n_super)),
    flag_harm = rep(df_super$flag_harm, 2)
  )

  hr_overall <- exp(coxph(Surv(time, event) ~ treat, data = df_temp)$coefficients)

  hr_results <- list(overall = hr_overall)

  AHR <- with(df_super, exp(mean(loghr_po)))
  hr_results$AHR <- AHR

  AHR_harm <- with(subset(df_super, flag_harm == 1), exp(mean(loghr_po)))
  hr_results$AHR_harm <- AHR_harm


  AHR_no_harm <- with(subset(df_super, flag_harm == 0), exp(mean(loghr_po)))
  hr_results$AHR_no_harm <- AHR_no_harm

  if (model == "alt" && sum(df_super$flag_harm) > 0) {
    hr_harm <- exp(coxph(Surv(time, event) ~ treat,
                        data = subset(df_temp, flag_harm == 1))$coefficients)
    hr_no_harm <- exp(coxph(Surv(time, event) ~ treat,
                           data = subset(df_temp, flag_harm == 0))$coefficients)

    hr_results$harm_subgroup <- hr_harm
    hr_results$no_harm_subgroup <- hr_no_harm

    if (verbose) {
      cat("\n=== Hazard Ratios (super popln)===\n")
      cat("Overall HR:", round(hr_overall, 3), "\n")
      cat("Causal AHR:", round(AHR, 3), "\n")
      cat("Harm subgroup HR:", round(hr_harm, 3), "\n")
      cat("Harm subgroup AHR:", round(AHR_harm, 3), "\n")
      cat("No-harm subgroup HR:", round(hr_no_harm, 3), "\n")
      cat("No-harm subgroup AHR:", round(AHR_no_harm, 3), "\n")
    }
  } else if (verbose) {
    cat("\n=== Hazard Ratios (super popln) ===\n")
    cat("Overall HR:", round(hr_overall, 3), "\n")
    cat("Causal AHR:", round(AHR, 3), "\n")
  }

  # ============================================================================
  # Prepare Censoring Parameters
  # ============================================================================

  cens_model <- NULL

  if (cens_type == "weibull") {
    # Fit censoring model
    X_cens <- as.matrix(df_work[, c("treat", covariate_cols)])

    fit_cens <- survreg(Surv(y, 1 - event) ~ X_cens,
                       data = df_work,
                       dist = "weibull")

    mu_cens <- coef(fit_cens)[1]
    tau_cens <- fit_cens$scale
    gamma_cens <- coef(fit_cens)[-1]

    # Store censoring parameters
    cens_model <- list(
      mu = mu_cens,
      tau = tau_cens,
      gamma = gamma_cens,
      type = "weibull"
    )

    # Calculate censoring linear predictors for super population
    df_super$lin_pred_cens_1 <- X_treat[, 1:ncol(X_cens)] %*% gamma_cens
    df_super$lin_pred_cens_0 <- X_control[, 1:ncol(X_cens)] %*% gamma_cens

  } else if (cens_type == "uniform") {
    # Use provided or default uniform censoring parameters
    if (is.null(cens_params$min) || is.null(cens_params$max)) {
      # Default: use range of observed times
      cens_params$min <- min(df_work$y) * 0.5
      cens_params$max <- max(df_work$y) * 1.5
    }

    cens_model <- list(
      min = cens_params$min,
      max = cens_params$max,
      type = "uniform"
    )
  }

  # ============================================================================
  # Prepare Output
  # ============================================================================

  # Model parameters
  model_params <- list(
    mu = mu,
    tau = tau,
    gamma = gamma,
    b0 = b0,
    censoring = cens_model
  )

  # Subgroup information
  subgroup_info <- list(
    vars = subgroup_vars,
    cuts = subgroup_cuts,
    definitions = subgroup_definitions,
    size = sum(df_super$flag_harm),
    proportion = mean(df_super$flag_harm)
  )

  # Analysis variables (for downstream use)
  analysis_vars <- list(
    continuous = continuous_vars,
    factor = factor_vars,
    covariates = covariate_cols,
    treatment = "treat",
    outcome = "y_sim",
    event = "event_sim"
  )

  # Return comprehensive results
  results <- list(
    df_super = df_super,
    model_params = model_params,
    subgroup_info = subgroup_info,
    hazard_ratios = hr_results,
    analysis_vars = analysis_vars,
    model_type = model,
    n_super = n_super,
    seed = seed
  )

  class(results) <- c("aft_dgm_flex", "list")

  return(results)
}

#' Find Quantile for Target Subgroup Proportion
#'
#' Determines the quantile cutpoint that achieves a target proportion of
#' observations in a subgroup. Useful for calibrating subgroup sizes.
#'
#' @param data A data.frame containing the variable of interest
#' @param var_name Character string specifying the variable name to analyze
#' @param target_prop Numeric value between 0 and 1 specifying the target
#'   proportion of observations to be included in the subgroup
#' @param direction Character string: "less" for values <= cutpoint (default),
#'   "greater" for values > cutpoint
#' @param tol Numeric tolerance for root finding algorithm. Default is 0.0001
#'
#' @return A list containing:
#' \describe{
#'   \item{quantile}{The quantile value (between 0 and 1) that achieves the
#'     target proportion}
#'   \item{cutpoint}{The actual data value corresponding to this quantile}
#'   \item{actual_proportion}{The achieved proportion (should equal target_prop
#'     within tolerance)}
#' }
#'
#' @details
#' This function uses root finding (\code{uniroot}) to determine the quantile
#' that results in exactly the target proportion of observations being classified
#' into the subgroup. This is particularly useful when you want to ensure a
#' specific subgroup size regardless of the data distribution.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data(cancer)
#'
#' # Find ER cutpoint for 12.5% subgroup
#' result <- find_quantile_for_proportion(
#'   data = gbsg,
#'   var_name = "er",
#'   target_prop = 0.125,
#'   direction = "less"
#' )
#'
#' print(result)
#' # Use in subgroup definition
#' subgroup_cuts = list(
#'   er = list(type = "quantile", value = result$quantile)
#' )
#' }
#'
#' @seealso \code{\link{generate_aft_dgm_flex}}
#' @export
#' @importFrom stats quantile uniroot

find_quantile_for_proportion <- function(data, var_name, target_prop,
                                        direction = "less", tol = 0.0001) {

  var_data <- data[[var_name]]

  # Objective function
  obj_fun <- function(q) {
    cutpoint <- quantile(var_data, probs = q, na.rm = TRUE)
    if (direction == "less") {
      actual_prop <- mean(var_data <= cutpoint, na.rm = TRUE)
    } else {
      actual_prop <- mean(var_data > cutpoint, na.rm = TRUE)
    }
    return(actual_prop - target_prop)
  }

  # Find root
  result <- uniroot(obj_fun, interval = c(0, 1), tol = tol)

  return(list(
    quantile = result$root,
    cutpoint = quantile(var_data, probs = result$root, na.rm = TRUE),
    actual_proportion = target_prop
  ))
}

#' Simulate Survival Data from AFT Data Generating Mechanism
#'
#' Generates simulated survival data from a previously created AFT data generating
#' mechanism (DGM). Samples from the super population and generates survival times
#' with specified censoring.
#'
#' @param dgm An object of class "aft_dgm_flex" created by
#'   \code{\link{generate_aft_dgm_flex}}
#' @param n Integer specifying the sample size. If NULL (default), uses the
#'   entire super population
#' @param rand_ratio Numeric specifying the randomization ratio (treatment:control).
#'   Default is 1 (1:1 allocation)
#' @param max_follow Numeric specifying maximum follow-up time for administrative
#'   censoring. Default is Inf (no administrative censoring)
#' @param cens_adjust Numeric adjustment to censoring distribution on log scale.
#'   Positive values increase censoring, negative values decrease it. Default is 0
#' @param draw_treatment Logical indicating whether to redraw treatment assignment.
#'   If TRUE (default), reassigns treatment according to rand_ratio. If FALSE,
#'   keeps original treatment assignments from super population
#' @param seed Integer random seed for reproducibility. Default is NULL (no seed set)
#'
#' @return A data.frame containing simulated survival data with columns:
#' \describe{
#'   \item{id}{Subject identifier}
#'   \item{treat}{Treatment assignment (0 or 1)}
#'   \item{treat_sim}{Simulated treatment assignment (may differ from treat if
#'     draw_treatment = TRUE)}
#'   \item{flag_harm}{Subgroup indicator (1 if all subgroup conditions met, 0 otherwise)}
#'   \item{z_*}{Standardized covariate values}
#'   \item{y_sim}{Observed survival time (minimum of true time and censoring time)}
#'   \item{event_sim}{Event indicator (1 = event observed, 0 = censored)}
#'   \item{t_true}{True underlying survival time (before censoring)}
#'   \item{c_time}{Censoring time}
#' }
#'
#' @details
#' ## Simulation Process
#'
#' 1. **Sampling**: Draws n observations with replacement from the super population
#' 2. **Treatment Assignment**:
#'    - If `draw_treatment = TRUE`: Reassigns treatment based on `rand_ratio`
#'    - If `draw_treatment = FALSE`: Keeps original treatment assignments
#' 3. **Survival Times**: Generates from Weibull AFT model:
#'    \deqn{\log(T) = \mu + \sigma \epsilon + X'\gamma}
#'    where ε ~ extreme value distribution
#' 4. **Censoring**: Applies specified censoring distribution (Weibull or uniform)
#' 5. **Administrative Censoring**: Applies max_follow cutoff if specified
#'
#' ## Censoring Adjustment
#'
#' The `cens_adjust` parameter modifies the censoring distribution:
#' - `cens_adjust = log(2)` doubles expected censoring times
#' - `cens_adjust = log(0.5)` halves expected censoring times
#'
#' @examples
#' \dontrun{
#' # Create DGM first
#' dgm <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(er = 20, meno = 0),
#'   model = "alt"
#' )
#'
#' # Simulate data with 1:1 randomization
#' sim_data <- simulate_from_dgm(
#'   dgm = dgm,
#'   n = 1000,
#'   rand_ratio = 1,
#'   max_follow = 84,
#'   cens_adjust = log(1.5),
#'   seed = 123
#' )
#'
#' # Check results
#' table(sim_data$treat_sim)
#' mean(sim_data$event_sim)
#'
#' # Simulate with 2:1 randomization
#' sim_data_2to1 <- simulate_from_dgm(
#'   dgm = dgm,
#'   n = 900,
#'   rand_ratio = 2,  # 2:1 treatment:control
#'   seed = 456
#' )
#' }
#'
#' @seealso
#' \code{\link{generate_aft_dgm_flex}} for creating the DGM
#'
#' @export
#' @importFrom stats rexp runif pmin

simulate_from_dgm <- function(dgm,
                             n = NULL,
                             rand_ratio = 1,
                             max_follow = Inf,
                             cens_adjust = 0,
                             draw_treatment = TRUE,
                             seed = NULL) {

  if (!inherits(dgm, c("aft_dgm_flex", "aft_dgm"))) {
    stop("dgm must be an object created by generate_aft_dgm_flex()")
  }

  if (!is.null(seed)) set.seed(seed)

  df_super <- dgm$df_super
  params <- dgm$model_params

  # Determine sample size
  if (is.null(n)) {
    df_sim <- df_super
    n <- nrow(df_sim)
  } else {
    # Sample from super population

    idx_sample <- sample(1:nrow(df_super), size = n, replace = TRUE)
    df_sim <- df_super[idx_sample, ]

    # Reassign treatment if draw_treatment, otherwise retain per super_population
    # Set to original assignent in df_super
    df_sim$treat_sim <- df_sim$treat

    if(draw_treatment){
    n_treat <- round(n * rand_ratio / (1 + rand_ratio))
    n_control <- n - n_treat
    df_sim$treat_sim[1:n_treat] <- 1
    df_sim$treat_sim[(n_treat + 1):n] <- 0
    }

    # Update linear predictors based on assigned treatment
    df_sim$lin_pred_obs <- ifelse(df_sim$treat_sim == 1,
                                  df_sim$lin_pred_1,
                                  df_sim$lin_pred_0)

    # Reset IDs
    df_sim$id <- 1:n
  }

  # Generate survival times
  epsilon <- log(rexp(n))  # Extreme value distribution
  logT_sim <- params$mu + params$tau * epsilon + df_sim$lin_pred_obs
  T_sim <- exp(logT_sim)

  # Generate censoring times
  if (params$censoring$type == "weibull") {
    # Weibull censoring
    lin_pred_cens <- ifelse(df_sim$treat_sim == 1,
                           df_sim$lin_pred_cens_1,
                           df_sim$lin_pred_cens_0)

    epsilon_cens <- log(rexp(n))
    logC_sim <- params$censoring$mu + cens_adjust +
                params$censoring$tau * epsilon_cens + lin_pred_cens
    C_sim <- exp(logC_sim)

  } else if (params$censoring$type == "uniform") {
    # Uniform censoring
    C_sim <- runif(n,
                  min = params$censoring$min,
                  max = params$censoring$max)
  }

  # Apply administrative censoring
  C_sim <- pmin(C_sim, max_follow)

  # Observed times and events
  df_sim$y_sim <- pmin(T_sim, C_sim)
  df_sim$event_sim <- ifelse(T_sim <= C_sim, 1, 0)
  df_sim$t_true <- T_sim
  df_sim$c_time <- C_sim

  return(df_sim)
}

#' Print Method for AFT DGM Objects
#'
#' Provides a formatted summary of an AFT data generating mechanism object.
#'
#' @param x An object of class "aft_dgm_flex"
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#'
#' @examples
#' \dontrun{
#' dgm <- generate_aft_dgm_flex(data = gbsg, ...)
#' print(dgm)
#' }
#'
#' @export
#' @method print aft_dgm_flex

print.aft_dgm_flex <- function(x, ...) {
  cat("AFT Data Generating Mechanism\n")
  cat("=============================\n")
  cat("Model type:", x$model_type, "\n")
  cat("Super population size:", x$n_super, "\n")
  cat("Number of covariates:", length(x$analysis_vars$covariates), "\n")

  if (x$model_type == "alt" && !is.null(x$subgroup_info$vars)) {
    cat("\nSubgroup Information:\n")
    cat("  Variables:", paste(x$subgroup_info$vars, collapse = ", "), "\n")
    cat("  Size:", x$subgroup_info$size, "\n")
    cat("  Proportion:", round(x$subgroup_info$proportion, 3), "\n")

    if (!is.null(x$subgroup_info$definitions)) {
      cat("  Definitions:\n")
      for (def in x$subgroup_info$definitions) {
        cat("    -", def, "\n")
      }
    }
  }

  cat("\nHazard Ratios:\n")
  for (name in names(x$hazard_ratios)) {
    cat("  ", gsub("_", " ", name), ":", round(x$hazard_ratios[[name]], 3), "\n")
  }

  cat("\nModel Parameters:\n")
  cat("  Intercept (mu):", round(x$model_params$mu, 3), "\n")
  cat("  Scale (sigma):", round(x$model_params$sigma, 3), "\n")
  if ("treat" %in% names(x$model_params$gamma)) {
    cat("  Treatment effect:", round(x$model_params$gamma["treat"], 3), "\n")
  }
  if ("treat_harm" %in% names(x$model_params$gamma)) {
    cat("  Interaction effect:", round(x$model_params$gamma["treat_harm"], 3), "\n")
  }

  invisible(x)
}

#' Summary Method for AFT DGM Objects
#'
#' Provides a comprehensive summary of an AFT data generating mechanism object.
#'
#' @param object An object of class "aft_dgm_flex"
#' @param ... Additional arguments (currently unused)
#'
#' @return A list of class "summary.aft_dgm_flex" containing summary information
#'
#' @examples
#' \dontrun{
#' dgm <- generate_aft_dgm_flex(data = gbsg, ...)
#' summary(dgm)
#' }
#'
#' @export
#' @method summary aft_dgm_flex

summary.aft_dgm_flex <- function(object, ...) {
  summary_list <- list(
    model_type = object$model_type,
    n_super = object$n_super,
    n_covariates = length(object$analysis_vars$covariates),
    continuous_vars = object$analysis_vars$continuous,
    factor_vars = object$analysis_vars$factor,
    subgroup_info = object$subgroup_info,
    hazard_ratios = object$hazard_ratios,
    model_params = object$model_params
  )

  class(summary_list) <- "summary.aft_dgm_flex"
  return(summary_list)
}

