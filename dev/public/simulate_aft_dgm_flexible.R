
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


