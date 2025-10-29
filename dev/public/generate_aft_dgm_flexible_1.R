# Enhanced AFT Data Generating Mechanism with Flexible Subgroup Definitions
# Supports multiple cutpoint specifications: fixed values, quantiles, functions, etc.

#' Generate Synthetic Survival Data using AFT Model with Flexible Subgroups
#'
#' @param data The input dataset (data.frame)
#' @param continuous_vars Character vector of continuous variable names
#' @param factor_vars Character vector of factor/categorical variable names
#' @param outcome_var Name of the outcome/time variable
#' @param event_var Name of the event/status variable
#' @param treatment_var Name of the treatment variable (if NULL, will be simulated)
#' @param subgroup_vars Character vector of variables defining the subgroup (optional)
#' @param subgroup_cuts Named list of cutpoint specifications (see details)
#' @param model Character: "alt" (with subgroup effects) or "null" (no subgroup effects)
#' @param k_treat Numeric: treatment effect modifier (default = 1)
#' @param k_inter Numeric: interaction effect modifier (default = 1)
#' @param n_super Integer: size of super population (default = 5000)
#' @param cens_type Character: "weibull" or "uniform" censoring
#' @param cens_params List: parameters for censoring distribution
#' @param seed Integer: random seed for reproducibility
#' @param verbose Logical: print diagnostic information
#'
#' @details
#' The subgroup_cuts parameter accepts flexible specifications:
#'   - Numeric value: subgroup_cuts = list(er = 20) means er <= 20
#'   - Quantile specification: subgroup_cuts = list(er = list(type = "quantile", value = 0.25))
#'   - Function: subgroup_cuts = list(er = list(type = "function", fun = median))
#'   - Range: subgroup_cuts = list(er = list(type = "range", min = 10, max = 50))
#'   - Multiple conditions: subgroup_cuts = list(er = list(type = "multiple", values = c(10, 20, 30)))
#'   - Custom function: subgroup_cuts = list(er = list(type = "custom", fun = function(x) x <= quantile(x, 0.3)))
#'
#' @return List containing:
#'   - df_super: Super population dataset
#'   - model_params: Model parameters (coefficients, scale, etc.)
#'   - subgroup_info: Information about subgroups
#'   - hazard_ratios: True hazard ratios
#'
#' @examples
#' # Example with various cutpoint specifications
#' dgm <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "pgr", "age"),
#'   subgroup_cuts = list(
#'     er = list(type = "quantile", value = 0.25),        # er <= 25th percentile
#'     pgr = list(type = "function", fun = median),       # pgr <= median
#'     age = list(type = "range", min = 40, max = 60)     # age between 40 and 60
#'   ),
#'   model = "alt",
#'   verbose = TRUE
#' )

generate_aft_dgm_flex <- function(data,
                                 continuous_vars,
                                 factor_vars,
                                 outcome_var,
                                 event_var,
                                 treatment_var = NULL,
                                 subgroup_vars = NULL,
                                 subgroup_cuts = NULL,
                                 model = "alt",
                                 k_treat = 1,
                                 k_inter = 1,
                                 n_super = 5000,
                                 cens_type = "weibull",
                                 cens_params = list(),
                                 seed = 8316951,
                                 verbose = TRUE) {

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
  required_vars <- c(outcome_var, event_var)
  if (!is.null(treatment_var)) required_vars <- c(required_vars, treatment_var)

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
    event = ifelse(data[[event_var]] == 1, 1, 0)
  )

  # Add treatment (or simulate if not provided)
  if (!is.null(treatment_var)) {
    df_work$treat <- data[[treatment_var]]
  } else {
    df_work$treat <- rbinom(nrow(data), size = 1, prob = 0.5)
    if (verbose) cat("Treatment variable simulated (50/50 randomization)\n")
  }

  # ============================================================================
  # Process Covariates
  # ============================================================================

  # Process continuous variables
  for (var in continuous_vars) {
    df_work[[paste0("z_", var)]] <- scale(data[[var]])[, 1]  # Standardize
  }

  # Process factor variables
  for (var in factor_vars) {
    if (is.factor(data[[var]])) {
      # Create dummy variables for factors
      dummies <- model.matrix(~ data[[var]] - 1)
      # Keep all but first level (reference)
      if (ncol(dummies) > 1) {
        for (j in 2:ncol(dummies)) {
          dummy_name <- paste0("z_", var, "_", j-1)
          df_work[[dummy_name]] <- dummies[, j]
        }
      }
    } else {
      # Treat as binary
      df_work[[paste0("z_", var)]] <- as.numeric(data[[var]])
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
  sigma <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]  # Coefficients (excluding intercept)

  # Apply effect modifiers
  gamma["treat"] <- k_treat * gamma["treat"]
  if ("treat_harm" %in% names(gamma)) {
    gamma["treat_harm"] <- k_inter * gamma["treat_harm"]
  }

  # Weibull parameterization
  b_weibull <- gamma
  b0_weibull <- -gamma / sigma

  if (verbose) {
    cat("\n=== Model Parameters ===\n")
    cat("Intercept (mu):", round(mu, 3), "\n")
    cat("Scale (sigma):", round(sigma, 3), "\n")
    cat("Treatment effect:", round(gamma["treat"], 3), "\n")
    if ("treat_harm" %in% names(gamma)) {
      cat("Interaction effect:", round(gamma["treat_harm"], 3), "\n")
    }
  }

  # ============================================================================
  # Generate Super Population
  # ============================================================================

  # Sample with replacement to create super population
  n_treat <- round(n_super / 2)
  n_control <- n_super - n_treat

  # Sample indices
  idx_sample <- sample(1:nrow(df_work), size = n_super, replace = TRUE)
  df_super <- df_work[idx_sample, ]

  # Assign treatment
  df_super$treat[1:n_treat] <- 1
  df_super$treat[(n_treat + 1):n_super] <- 0

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

  # Linear predictors
  df_super$lin_pred_1 <- X_treat %*% b_weibull
  df_super$lin_pred_0 <- X_control %*% b_weibull
  df_super$lin_pred_obs <- X_super %*% b_weibull

  # Hazard ratios (individual level)
  df_super$hr_individual <- exp((df_super$lin_pred_1 - df_super$lin_pred_0) / sigma)


  #   # linear predictor (setting treat=1)
  # z.1 <- zmat.1[,"z"]
  # # log(Y) AFT parameterization
  # # psi.true are AFT parameters
  # eta1 <- mu + c(zmat.1%*%gamma.true)+w*gamma_w
  # # Weibull log(hazard ratio) setting treat=1
  # # and excluding mu and w*bw (since taking difference below)
  # phi1 <- (-1)*c(zmat.1%*%gamma.true)/tau.strataO
  # log.Y1 <- eta1 + tau.strataO*epsilon
  #
  # # Setting treat=0
  # z.0 <- zmat.0[,"z"]
  # eta0 <- mu + c(zmat.0%*%gamma.true)+w*gamma_w
  # log.Y0 <- eta0 + tau.strataO*epsilon
  # phi0 <- (-1)*c(zmat.0%*%gamma.true)/tau.strataO
  #
  # # PO hazards excluding baseline
  # # Used for calculating empirical version of CDEs (controlled direct effects)
  # # theta0 = exp(L0'beta)
  # theta0 <- -c(zmat.0%*%gamma.true + w*gamma_w)/tau.strataO
  # # theta1 = exp(L1'beta)
  # theta1 <- -c(zmat.1%*%gamma.true + w*gamma_w)/tau.strataO

  # Potential outcome log(hr) difference
  #loghr.po <- phi1-phi0

  #ahr_empirical <- with(dfs,exp(mean(loghr.po)))
  #res$AHR <- ahr_empirical

  # ============================================================================
  # Calculate True Hazard Ratios
  # ============================================================================

  # Generate potential outcomes for HR calculation
  epsilon <- log(rexp(n_super))  # Extreme value distribution

  # Under treatment
  logT_1 <- mu + sigma * epsilon + df_super$lin_pred_1
  T_1 <- exp(logT_1)

  # Under control
  logT_0 <- mu + sigma * epsilon + df_super$lin_pred_0
  T_0 <- exp(logT_0)

  # Calculate empirical hazard ratios
  df_temp <- data.frame(
    time = c(T_1, T_0),
    event = 1,
    treat = c(rep(1, n_super), rep(0, n_super)),
    flag_harm = rep(df_super$flag_harm, 2)
  )

  hr_overall <- exp(coxph(Surv(time, event) ~ treat, data = df_temp)$coefficients)

  hr_results <- list(overall = hr_overall)

  if (model == "alt" && sum(df_super$flag_harm) > 0) {
    hr_harm <- exp(coxph(Surv(time, event) ~ treat,
                        data = subset(df_temp, flag_harm == 1))$coefficients)
    hr_no_harm <- exp(coxph(Surv(time, event) ~ treat,
                           data = subset(df_temp, flag_harm == 0))$coefficients)

    hr_results$harm_subgroup <- hr_harm
    hr_results$no_harm_subgroup <- hr_no_harm

    if (verbose) {
      cat("\n=== Hazard Ratios ===\n")
      cat("Overall HR:", round(hr_overall, 3), "\n")
      cat("Harm subgroup HR:", round(hr_harm, 3), "\n")
      cat("No-harm subgroup HR:", round(hr_no_harm, 3), "\n")
    }
  } else if (verbose) {
    cat("\n=== Hazard Ratios ===\n")
    cat("Overall HR:", round(hr_overall, 3), "\n")
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
    sigma_cens <- fit_cens$scale
    gamma_cens <- coef(fit_cens)[-1]

    # Store censoring parameters
    cens_model <- list(
      mu = mu_cens,
      sigma = sigma_cens,
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
    sigma = sigma,
    gamma = gamma,
    b_weibull = b_weibull,
    b0_weibull = b0_weibull,
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

# ================================================================================
# Convenience function to find quantile for target subgroup proportion
# ================================================================================

#' Find Quantile for Target Subgroup Proportion
#'
#' @param data Data containing the variable
#' @param var_name Name of the variable
#' @param target_prop Target proportion for the subgroup
#' @param direction "less" for <=, "greater" for >
#' @param tol Tolerance for root finding
#'
#' @return The quantile value that achieves the target proportion

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

# ================================================================================
# The simulation function remains the same as before
# ================================================================================

simulate_from_dgm <- function(dgm,
                             n = NULL,
                             rand_ratio = 1,
                             max_follow = Inf,
                             cens_adjust = 0,
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
    n_treat <- round(n * rand_ratio / (1 + rand_ratio))
    n_control <- n - n_treat

    idx_sample <- sample(1:nrow(df_super), size = n, replace = TRUE)
    df_sim <- df_super[idx_sample, ]

    # Reassign treatment
    df_sim$treat[1:n_treat] <- 1
    df_sim$treat[(n_treat + 1):n] <- 0

    # Update linear predictors based on assigned treatment
    df_sim$lin_pred_obs <- ifelse(df_sim$treat == 1,
                                  df_sim$lin_pred_1,
                                  df_sim$lin_pred_0)

    # Reset IDs
    df_sim$id <- 1:n
  }

  # Generate survival times
  epsilon <- log(rexp(n))  # Extreme value distribution
  logT_sim <- params$mu + params$sigma * epsilon + df_sim$lin_pred_obs
  T_sim <- exp(logT_sim)

  # Generate censoring times
  if (params$censoring$type == "weibull") {
    # Weibull censoring
    lin_pred_cens <- ifelse(df_sim$treat == 1,
                           df_sim$lin_pred_cens_1,
                           df_sim$lin_pred_cens_0)

    epsilon_cens <- log(rexp(n))
    logC_sim <- params$censoring$mu + cens_adjust +
                params$censoring$sigma * epsilon_cens + lin_pred_cens
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

  # Add analysis variables
  analysis_cols <- c("id", "treat", "flag_harm",
                    dgm$analysis_vars$covariates,
                    "y_sim", "event_sim", "t_true")

  # Keep only necessary columns
  df_sim <- df_sim[, intersect(analysis_cols, names(df_sim))]

  return(df_sim)
}



# Testing

# Your original request - using quantiles
dgm <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "pgr"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),   # er <= 25th percentile
    pgr = list(type = "quantile", value = 0.50)   # pgr <= median
  ),
  model = "alt"
)



# Your original request - using quantiles
dgm <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),   # er <= 25th percentile
    meno = 0
  ),
  model = "alt",
  k_inter = 1.5
)


result <- find_k_inter_for_target_hr(
  target_hr_harm = 2.0,
  data = gbsg,
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
    continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),
    meno = 0
  ),
  k_treat = 1.0
)

# Result: k_inter = 1.00 achieves HR_harm = 2.0



base_params <- list(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),
    meno = 0
  ),
  k_treat = 1.0,
  n_super = 50000  # Using smaller for faster demonstration
)


sensitivity_results <- do.call(sensitivity_analysis_k_inter, c(
  list(
    k_inter_range = c(-3, 3),
    n_points = 11,
    model = "alt"
  ),
  base_params
))

cat("\nSensitivity results:\n")
print(round(sensitivity_results, 3))




