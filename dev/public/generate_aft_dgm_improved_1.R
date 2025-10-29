# Improved and Generalized AFT Data Generating Mechanism
# This is a refactored version of dgm_aftm4_gbsg that can work with any dataset

#' Generate Synthetic Survival Data using AFT Model
#'
#' @param data The input dataset (data.frame)
#' @param continuous_vars Character vector of continuous variable names
#' @param factor_vars Character vector of factor/categorical variable names
#' @param outcome_var Name of the outcome/time variable
#' @param event_var Name of the event/status variable
#' @param treatment_var Name of the treatment variable (if NULL, will be simulated)
#' @param subgroup_vars Character vector of variables defining the subgroup (optional)
#' @param subgroup_cuts Named list of cutpoints for subgroup variables (optional)
#' @param model Character: "alt" (with subgroup effects) or "null" (no subgroup effects)
#' @param k_treat Numeric: treatment effect modifier (default = 1)
#' @param k_inter Numeric: interaction effect modifier (default = 1)
#' @param n_super Integer: size of super population (default = 5000)
#' @param cens_type Character: "weibull" or "uniform" censoring
#' @param cens_params List: parameters for censoring distribution
#' @param seed Integer: random seed for reproducibility
#' @param verbose Logical: print diagnostic information
#'
#' @return List containing:
#'   - df_super: Super population dataset
#'   - model_params: Model parameters (coefficients, scale, etc.)
#'   - subgroup_info: Information about subgroups
#'   - hazard_ratios: True hazard ratios
#'
#' @examples
#' # Using with GBSG data
#' library(survival)
#' data(cancer)
#'
#' dgm <- generate_aft_dgm(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   model = "alt",
#'   verbose = TRUE
#' )

generate_aft_dgm <- function(data,
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
  # Define Subgroups (if specified)
  # ============================================================================

  df_work$flag_harm <- 0
  interaction_term <- NULL

  if (model == "alt" && !is.null(subgroup_vars)) {
    # Create subgroup indicators
    subgroup_indicators <- list()

    for (var in subgroup_vars) {
      if (var %in% continuous_vars) {
        # Use median or provided cutpoint
        cutpoint <- ifelse(!is.null(subgroup_cuts[[var]]),
                          subgroup_cuts[[var]],
                          median(data[[var]], na.rm = TRUE))
        subgroup_indicators[[var]] <- data[[var]] <= cutpoint
      } else {
        # For factor variables, use first level or provided specification
        if (!is.null(subgroup_cuts[[var]])) {
          subgroup_indicators[[var]] <- data[[var]] == subgroup_cuts[[var]]
        } else {
          subgroup_indicators[[var]] <- data[[var]] == levels(as.factor(data[[var]]))[1]
        }
      }
    }

    # Create harm flag (all subgroup conditions met)
    df_work$flag_harm <- as.numeric(Reduce("&", subgroup_indicators))

    # Create interaction term
    if (length(subgroup_indicators) > 0) {
      interaction_term <- df_work$treat * df_work$flag_harm
    }

    if (verbose) {
      cat("Subgroup defined by:", paste(subgroup_vars, collapse = " AND "), "\n")
      cat("Subgroup size:", sum(df_work$flag_harm), "out of", nrow(df_work), "\n")
      cat("Subgroup proportion:", mean(df_work$flag_harm), "\n")
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
    cat("Intercept (mu):", mu, "\n")
    cat("Scale (sigma):", sigma, "\n")
    cat("Treatment effect:", gamma["treat"], "\n")
    if ("treat_harm" %in% names(gamma)) {
      cat("Interaction effect:", gamma["treat_harm"], "\n")
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
      cat("Overall HR:", hr_overall, "\n")
      cat("Harm subgroup HR:", hr_harm, "\n")
      cat("No-harm subgroup HR:", hr_no_harm, "\n")
    }
  } else if (verbose) {
    cat("\n=== Hazard Ratios ===\n")
    cat("Overall HR:", hr_overall, "\n")
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

  class(results) <- c("aft_dgm", "list")

  return(results)
}

# ================================================================================
# Simulation Function - Generate Data from DGM
# ================================================================================

#' Simulate Data from AFT DGM
#'
#' @param dgm Object returned by generate_aft_dgm
#' @param n Sample size (if NULL, uses super population size)
#' @param rand_ratio Randomization ratio (treatment:control)
#' @param max_follow Maximum follow-up time
#' @param cens_adjust Adjustment to censoring (on log scale)
#' @param seed Random seed
#'
#' @return Data frame with simulated survival data

simulate_from_dgm <- function(dgm,
                             n = NULL,
                             rand_ratio = 1,
                             max_follow = Inf,
                             cens_adjust = 0,
                             seed = NULL) {

  if (!inherits(dgm, "aft_dgm")) {
    stop("dgm must be an object created by generate_aft_dgm()")
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

# ================================================================================
# Helper Functions
# ================================================================================

#' Print summary of AFT DGM
#'
#' @param x An aft_dgm object
#' @param ... Additional arguments

print.aft_dgm <- function(x, ...) {
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
  }

  cat("\nHazard Ratios:\n")
  for (name in names(x$hazard_ratios)) {
    cat("  ", name, ":", round(x$hazard_ratios[[name]], 3), "\n")
  }

  cat("\nModel Parameters:\n")
  cat("  Intercept (mu):", round(x$model_params$mu, 3), "\n")
  cat("  Scale (sigma):", round(x$model_params$sigma, 3), "\n")
  cat("  Treatment effect:", round(x$model_params$gamma["treat"], 3), "\n")
}

#' Summary of simulated data
#'
#' @param data Simulated dataset
#' @param dgm Optional: original DGM for comparison

summarize_simulation <- function(data, dgm = NULL) {
  cat("Simulated Data Summary\n")
  cat("======================\n")
  cat("Sample size:", nrow(data), "\n")
  cat("Treatment allocation:", table(data$treat), "\n")
  cat("Event rate:", mean(data$event_sim), "\n")
  cat("Median follow-up:", median(data$y_sim), "\n")

  if ("flag_harm" %in% names(data)) {
    cat("\nSubgroup sizes:\n")
    cat("  Harm subgroup:", sum(data$flag_harm), "\n")
    cat("  No-harm subgroup:", sum(1 - data$flag_harm), "\n")
  }

  # Calculate observed hazard ratio
  hr_obs <- exp(coxph(Surv(y_sim, event_sim) ~ treat, data = data)$coefficients)
  cat("\nObserved HR (overall):", round(hr_obs, 3), "\n")

  if ("flag_harm" %in% names(data) && sum(data$flag_harm) > 0) {
    hr_harm_obs <- exp(coxph(Surv(y_sim, event_sim) ~ treat,
                            data = subset(data, flag_harm == 1))$coefficients)
    hr_no_harm_obs <- exp(coxph(Surv(y_sim, event_sim) ~ treat,
                               data = subset(data, flag_harm == 0))$coefficients)

    cat("Observed HR (harm subgroup):", round(hr_harm_obs, 3), "\n")
    cat("Observed HR (no-harm subgroup):", round(hr_no_harm_obs, 3), "\n")
  }

  if (!is.null(dgm)) {
    cat("\nTrue hazard ratios from DGM:\n")
    for (name in names(dgm$hazard_ratios)) {
      cat("  ", name, ":", round(dgm$hazard_ratios[[name]], 3), "\n")
    }
  }
}

# ================================================================================
# Example Usage
# ================================================================================

if (FALSE) {  # Set to TRUE to run example

  library(survival)

  # Example 1: Using with GBSG data
  data(cancer)

  # Define the DGM
  dgm <- generate_aft_dgm(
    data = gbsg,
    continuous_vars = c("age", "size", "nodes", "pgr", "er"),
    factor_vars = c("meno", "grade"),
    outcome_var = "rfstime",
    event_var = "status",
    treatment_var = "hormon",
    subgroup_vars = c("er", "meno"),
    subgroup_cuts = list(er = 20, meno = 0),  # er <= 20 and meno == 0
    model = "alt",
    k_treat = 0.9,
    k_inter = 2,
    n_super = 5000,
    verbose = TRUE
  )

  # Print DGM summary
  print(dgm)

  # Simulate data
  sim_data <- simulate_from_dgm(
    dgm = dgm,
    n = 700,
    rand_ratio = 1,
    max_follow = 84,
    cens_adjust = log(1.5),
    seed = 123
  )

  # Summarize simulation
  summarize_simulation(sim_data, dgm)

  # Example 2: Using with custom dataset
  set.seed(42)
  custom_data <- data.frame(
    time = rexp(500, rate = 0.01),
    status = rbinom(500, 1, 0.7),
    age = rnorm(500, 50, 10),
    biomarker = rgamma(500, 2, 1),
    stage = sample(1:4, 500, replace = TRUE),
    sex = sample(c("M", "F"), 500, replace = TRUE),
    treatment = rbinom(500, 1, 0.5)
  )

  dgm_custom <- generate_aft_dgm(
    data = custom_data,
    continuous_vars = c("age", "biomarker"),
    factor_vars = c("stage", "sex"),
    outcome_var = "time",
    event_var = "status",
    treatment_var = "treatment",
    subgroup_vars = c("biomarker", "sex"),
    subgroup_cuts = list(biomarker = 2, sex = "F"),
    model = "alt",
    verbose = TRUE
  )

  sim_custom <- simulate_from_dgm(dgm_custom, n = 1000)
  summarize_simulation(sim_custom, dgm_custom)
}
