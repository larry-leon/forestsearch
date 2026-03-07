#' Prepare Censoring Model Parameters
#'
#' Constructs the censoring model object and appends per-subject counterfactual
#' censoring linear predictors (\code{lin_pred_cens_0}, \code{lin_pred_cens_1})
#' to the super-population data frame.
#'
#' ## Linear predictor convention
#'
#' \code{lin_pred_cens_0} and \code{lin_pred_cens_1} store the
#' \strong{covariate contribution only} — i.e. \eqn{\gamma_c' X}, with the
#' intercept \eqn{\mu_c} excluded.  This matches the convention used for the
#' outcome model (\code{lin_pred_0}, \code{lin_pred_1} = \eqn{\gamma' X},
#' no intercept) computed in \code{calculate_linear_predictors()}.
#'
#' \code{simulate_from_dgm()} reconstructs the full log-censoring time as:
#' \deqn{\log C = \mu_c + \delta + \tau_c \epsilon + \gamma_c' X}
#' where \eqn{\mu_c} = \code{params$censoring$mu},
#' \eqn{\delta} = \code{cens_adjust},
#' \eqn{\tau_c} = \code{params$censoring$tau}, and
#' \eqn{\gamma_c' X} = \code{lin_pred_cens_{0|1}}.
#'
#' When \code{select_censoring = TRUE}, \code{predict(survreg, type = "linear")}
#' returns the full linear predictor \eqn{\mu_c + \gamma_c' X}.  The stored
#' intercept \eqn{\mu_c} is therefore subtracted before writing
#' \code{lin_pred_cens_*}, so that \code{simulate_from_dgm()} can add
#' \code{params$censoring$mu} exactly once.  Omitting this subtraction causes
#' \eqn{\mu_c} to be counted twice, producing astronomically large censoring
#' times and universal censoring.
#'
#' When \code{select_censoring = FALSE} with a Weibull/lognormal
#' \code{cens_type}, the intercept-only model has zero covariate contribution,
#' so \code{lin_pred_cens_0 = lin_pred_cens_1 = 0}.  Storing \code{mu} instead
#' of 0 causes the same double-counting.
#'
#' @param df_work Working data frame (output of \code{prepare_working_dataset}).
#' @param cens_type Character. \code{"weibull"} or \code{"uniform"}.
#' @param cens_params Named list of user-supplied censoring parameters.
#' @param df_super Super-population data frame; receives
#'   \code{lin_pred_cens_0} and \code{lin_pred_cens_1} columns.
#' @param select_censoring Logical. If \code{TRUE} (default), fits the
#'   censoring distribution from observed data using AIC-based \code{survreg}
#'   model comparison. If \code{FALSE}, uses \code{cens_params} directly with
#'   no model fitting. See \code{\link{generate_aft_dgm_flex}} for the required
#'   \code{cens_params} structure under each combination of
#'   \code{select_censoring} and \code{cens_type}.
#' @param verbose Logical. If \code{TRUE} (default), prints the censoring
#'   model comparison table and recommendation. Set to \code{FALSE} to
#'   suppress all censoring model selection output.
#'
#' @return A named list:
#' \describe{
#'   \item{cens_model}{List of censoring distribution parameters stored in
#'     \code{dgm$model_params$censoring}.}
#'   \item{df_super}{Updated super-population data frame with
#'     \code{lin_pred_cens_0} and \code{lin_pred_cens_1} appended.  These
#'     hold covariate contributions only (\eqn{\gamma_c' X}); the intercept
#'     is excluded.}
#' }
#'
#' @keywords internal
prepare_censoring_model <- function(df_work,
                                    cens_type,
                                    cens_params,
                                    df_super,
                                    select_censoring = TRUE,
                                    verbose = TRUE) {

  cens_model <- NULL

  # ===========================================================================
  # PATH A: select_censoring = FALSE
  # Use caller-supplied cens_params directly; no model fitting.
  # ===========================================================================

  if (!select_censoring) {

    if (cens_type == "uniform") {
      # -----------------------------------------------------------------------
      # Uniform: use supplied min/max, defaulting to data range with a message.
      # simulate_from_dgm() does not read lin_pred_cens for the uniform branch,
      # but columns are added for structural consistency.
      # -----------------------------------------------------------------------
      if (is.null(cens_params$min)) {
        cens_params$min <- min(df_work$y, na.rm = TRUE) * 0.5
        message("select_censoring = FALSE: cens_params$min not supplied; ",
                "defaulting to 0.5 * min(y) = ", round(cens_params$min, 3))
      }
      if (is.null(cens_params$max)) {
        cens_params$max <- max(df_work$y, na.rm = TRUE) * 1.5
        message("select_censoring = FALSE: cens_params$max not supplied; ",
                "defaulting to 1.5 * max(y) = ", round(cens_params$max, 3))
      }

      cens_model <- list(
        min  = cens_params$min,
        max  = cens_params$max,
        type = "uniform"
      )

      df_super$lin_pred_cens_0 <- 0
      df_super$lin_pred_cens_1 <- 0

    } else {
      # -----------------------------------------------------------------------
      # Weibull / lognormal: require mu + tau from cens_params.
      # Intercept-only model => covariate contribution gamma'X = 0.
      # simulate_from_dgm() will compute:
      #   logC = mu_cens + cens_adjust + tau*epsilon + 0    (CORRECT)
      # Storing mu here instead of 0 would double-count mu_cens.
      # -----------------------------------------------------------------------
      if (is.null(cens_params$mu) || is.null(cens_params$tau)) {
        stop(
          "When select_censoring = FALSE and cens_type = '", cens_type, "', ",
          "cens_params must supply both 'mu' (log-scale location) and ",
          "'tau' (scale). Elements received: ",
          paste(names(cens_params), collapse = ", ")
        )
      }

      cens_dist <- if (!is.null(cens_params$type)) cens_params$type else "weibull"

      cens_model <- list(
        mu    = cens_params$mu,
        tau   = cens_params$tau,
        gamma = numeric(0),   # intercept-only: no covariate effects
        type  = cens_dist
      )

      # Intercept-only => covariate contribution is zero for all subjects.
      # DO NOT store mu here: simulate_from_dgm() adds params$censoring$mu
      # separately, so lin_pred_cens must hold gamma'X only.
      df_super$lin_pred_cens_0 <- 0
      df_super$lin_pred_cens_1 <- 0
    }

    return(list(cens_model = cens_model, df_super = df_super))
  }

  # ===========================================================================
  # PATH B: select_censoring = TRUE
  # Fit censoring distribution from observed data via AIC comparison.
  # ===========================================================================

  if (cens_type != "uniform") {

    covariate_cols <- grep("^zcens_", names(df_work), value = TRUE)

    if (length(covariate_cols) >= 1) {
      X_cens <- as.matrix(df_work[, c("treat", covariate_cols)])
    } else {
      X_cens <- as.matrix(df_work[, c("treat")])
    }

    fit_cens1 <- survreg(Surv(y, 1 - event) ~ X_cens,
                         data = df_work, dist = "weibull")
    fit_cens2 <- survreg(Surv(y, 1 - event) ~ X_cens,
                         data = df_work, dist = "lognormal")
    fit_cens3 <- survreg(Surv(y, 1 - event) ~ 1,
                         data = df_work, dist = "weibull")
    fit_cens4 <- survreg(Surv(y, 1 - event) ~ 1,
                         data = df_work, dist = "lognormal")

    # Compare all 4 models; select best by AIC
    comparison <- compare_multiple_survreg(
      fit_cens1, fit_cens2, fit_cens3, fit_cens4,
      model_names = c("Weibull", "LogNormal", "Weibull0", "LogNormal0"),
      verbose = verbose
    )

    fit_cens   <- get_best_survreg(comparison)
    mu_cens    <- coef(fit_cens)[1]    # intercept on log-time scale
    tau_cens   <- fit_cens$scale
    gamma_cens <- coef(fit_cens)[-1]   # covariate coefficients (no intercept)

    cens_model <- list(
      mu    = mu_cens,
      tau   = tau_cens,
      gamma = gamma_cens,
      type  = c(fit_cens$dist)
    )

    # -------------------------------------------------------------------------
    # KEY FIX: store COVARIATE-ONLY contribution (gamma'X), intercept excluded.
    #
    # predict(survreg, type = "linear") returns the FULL linear predictor:
    #   mu_cens + gamma_cens' X
    #
    # Subtracting mu_cens leaves gamma_cens' X only, matching the contract
    # expected by simulate_from_dgm():
    #   logC = params$censoring$mu   <- added by simulate_from_dgm()
    #        + cens_adjust
    #        + params$censoring$tau * epsilon_cens
    #        + lin_pred_cens          <- gamma'X only (stored here)
    #
    # Without the subtraction, mu_cens appears twice and censoring times
    # become exp(2 * mu_cens + ...), which are astronomically large and
    # collapse all subjects to administrative censoring.
    # -------------------------------------------------------------------------

    # Counterfactual: control arm (treat := 0)
    newdata       <- df_super
    newdata$treat <- 0
    if (!all(newdata$treat == 0))
      stop("Error in creating counterfactual: could not set treat := 0")
    df_super$lin_pred_cens_0 <- predict(fit_cens, newdata = newdata,
                                        type = "linear") - mu_cens

    # Counterfactual: treatment arm (treat := 1)
    newdata       <- df_super
    newdata$treat <- 1
    if (!all(newdata$treat == 1))
      stop("Error in creating counterfactual: could not set treat := 1")
    df_super$lin_pred_cens_1 <- predict(fit_cens, newdata = newdata,
                                        type = "linear") - mu_cens

  } else {
    # -------------------------------------------------------------------------
    # Uniform censoring (select_censoring = TRUE)
    # Use supplied cens_params or default to observed data range.
    # -------------------------------------------------------------------------
    if (is.null(cens_params$min) || is.null(cens_params$max)) {
      cens_params$min <- min(df_work$y, na.rm = TRUE) * 0.5
      cens_params$max <- max(df_work$y, na.rm = TRUE) * 1.5
    }

    cens_model <- list(
      min  = cens_params$min,
      max  = cens_params$max,
      type = "uniform"
    )

    # Placeholder columns for structural consistency across censoring types
    df_super$lin_pred_cens_0 <- 0
    df_super$lin_pred_cens_1 <- 0
  }

  return(list(cens_model = cens_model, df_super = df_super))
}




#' Validate Input Parameters
#' @keywords internal
validate_inputs <- function(data, model, cens_type, outcome_var, event_var,
                            treatment_var, continuous_vars, factor_vars) {

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
    stop("Required variables not found in data: ",
         paste(missing_vars, collapse = ", "))
  }

  # Check continuous and factor variables
  all_covars <- c(continuous_vars, factor_vars)
  missing_covars <- setdiff(all_covars, names(data))
  if (length(missing_covars) > 0) {
    stop("Covariate variables not found in data: ",
         paste(missing_covars, collapse = ", "))
  }
}


#' Prepare Working Dataset with Processed Covariates
#' @keywords internal
prepare_working_dataset <- function(data, outcome_var, event_var, treatment_var,
                                    continuous_vars, factor_vars, standardize,
                                    continuous_vars_cens, factor_vars_cens,
                                    verbose) {

  # Create base working dataset
  df_work <- data.frame(
    id = 1:nrow(data),
    y = data[[outcome_var]],
    treat = data[[treatment_var]],
    event = ifelse(data[[event_var]] == 1, 1, 0)
  )

  # Process continuous variables
  df_work <- process_continuous_vars(df_work, data, continuous_vars,
                                     standardize)


  # Process factor variables
  df_work <- process_factor_vars(df_work, data, factor_vars)

  # Add unprocessed variables
  df_work <- add_unprocessed_vars(df_work, data, outcome_var, event_var,
                                  treatment_var, continuous_vars, factor_vars,
                                  verbose)

  if(!is.null(continuous_vars_cens)){
    df_work <- process_continuous_vars(df_work, data, continuous_vars_cens,
                                       standardize, marker ="zcens_")
  }

  if(!is.null(factor_vars_cens)){
    df_work <- process_factor_vars(df_work, data, factor_vars_cens, marker = "zcens_")
  }

  return(df_work)
}


#' Process Continuous Variables
#' @keywords internal
process_continuous_vars <- function(df_work, data, continuous_vars,
                                    standardize, marker = "z_") {

  for (var in continuous_vars) {
    if (standardize) {
      df_work[[paste0(marker, var)]] <- scale(data[[var]])[, 1]  # Standardize
    } else {
      df_work[[paste0(marker, var)]] <- data[[var]]
    }
  }

  return(df_work)
}


#' Process Factor Variables with Largest Value as Reference
#' @keywords internal
process_factor_vars <- function(df_work, data, factor_vars, marker = "z_") {

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
      # Binary variable - keep as-is (min as reference)
      if (is.numeric(all_levels)) {
        ref_level <- min(all_levels)
        other_level <- max(all_levels)
      } else {
        # For character, first alphabetically
        ref_level <- sort(all_levels, decreasing = FALSE)[1]
        other_level <- setdiff(all_levels, ref_level)
      }

      # Create single indicator for non-reference level
      df_work[[paste0(marker, var)]] <- as.numeric(data[[var]] == other_level)

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
        dummy_name <- paste0(marker, var, "_", level)
        df_work[[dummy_name]] <- as.numeric(data[[var]] == level)
      }
    }
  }

  return(df_work)
}


#' Add Unprocessed Variables from Original Data
#' @keywords internal
add_unprocessed_vars <- function(df_work, data, outcome_var, event_var,
                                 treatment_var, continuous_vars, factor_vars,
                                 verbose) {

  # Retain original columns for ALL covariates (both continuous and factor),
  # alongside their z_* / dummy counterparts created by the process_* helpers.
  # This ensures that downstream Cox strata() calls, subgroup subset
  # expressions such as "age <= 50" or "grade == 3", and any user-defined
  # binary cuts work correctly on simulated datasets returned by
  # simulate_from_dgm().  Only the three structural columns (outcome, event,
  # treatment) are excluded because they are replaced by simulation-specific
  # counterparts (y_sim, event_sim, treat_sim).
  processed <- c(outcome_var, event_var, treatment_var)

  # Find columns in source data not already present in df_work
  unprocessed <- setdiff(names(data), processed)

  if (length(unprocessed) > 0) {
    for (var in unprocessed) {
      if (!var %in% names(df_work)) {   # avoid duplicates
        df_work[[var]] <- data[[var]]
      }
    }
  }

  return(df_work)
}


#' Define Subgroups with Flexible Cutpoints
#' @keywords internal
define_subgroups <- function(df_work, data, subgroup_vars, subgroup_cuts,
                             continuous_vars, model, verbose) {

  flag_harm <- rep(0, nrow(df_work))
  interaction_term <- NULL
  definitions <- list()

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
        # Process continuous variable
        result <- process_continuous_subgroup(data[[var]], cut_spec, var,
                                              verbose)
        subgroup_indicators[[var]] <- result$indicator
        definitions[[var]] <- result$definition

      } else {
        # Process factor variable
        result <- process_factor_subgroup(data[[var]], cut_spec, var,
                                          verbose)
        subgroup_indicators[[var]] <- result$indicator
        definitions[[var]] <- result$definition
      }
    }

    # Validate that all indicators are non-NULL
    null_indicators <- sapply(subgroup_indicators, is.null)
    if (any(null_indicators)) {
      null_vars <- names(subgroup_indicators)[null_indicators]
      stop(paste("Failed to create subgroup indicators for variables:",
                 paste(null_vars, collapse = ", ")))
    }

    # Create harm flag (all subgroup conditions met)
    flag_harm <- as.numeric(Reduce("&", subgroup_indicators))

    # Create interaction term
    if (length(subgroup_indicators) > 0) {
      interaction_term <- df_work$treat * flag_harm
    }

    if (verbose) {
      cat("\nOverall subgroup (all conditions met):\n")
      cat("  Size:", sum(flag_harm), "out of", nrow(df_work), "\n")
      cat("  Proportion:", round(mean(flag_harm), 3), "\n")
    }
  }

  return(list(
    flag_harm = flag_harm,
    interaction_term = interaction_term,
    definitions = definitions
  ))
}


#' Process Continuous Variable for Subgroup Definition
#' @keywords internal
process_continuous_subgroup <- function(var_data, cut_spec, var_name, verbose) {

  indicator <- process_cutpoint(var_data, cut_spec, var_name, verbose)
  definition <- format_continuous_definition(var_data, cut_spec, var_name)

  if (verbose) {
    cat("  ", definition, "\n")
    cat("    Proportion in subgroup:",
        round(mean(indicator, na.rm = TRUE), 3), "\n")
  }

  return(list(indicator = indicator, definition = definition))
}


#' Process Factor Variable for Subgroup Definition
#' @keywords internal
process_factor_subgroup <- function(var_data, cut_spec, var_name, verbose) {

  # Initialize to avoid NULL issues
  indicator <- NULL
  definition <- NULL

  if (!is.null(cut_spec)) {
    # Handle simple numeric or character value (e.g., meno = 0 or meno = "1")
    if ((is.numeric(cut_spec) || is.character(cut_spec) || is.factor(cut_spec)) &&
        length(cut_spec) == 1 && !is.list(cut_spec)) {
      indicator <- var_data == cut_spec
      definition <- paste(var_name, "==", cut_spec)

    } else if (is.list(cut_spec) && !is.null(cut_spec$type)) {
      # Handle list specifications
      if (cut_spec$type == "multiple") {
        indicator <- var_data %in% cut_spec$values
        definition <- paste(var_name, "in",
                            paste(cut_spec$values, collapse = ", "))
      } else {
        # Unknown list type, use default
        first_level <- levels(as.factor(var_data))[1]
        indicator <- var_data == first_level
        definition <- paste(var_name, "==", first_level)
      }

    } else {
      # If cut_spec doesn't match expected patterns, use default
      first_level <- levels(as.factor(var_data))[1]
      indicator <- var_data == first_level
      definition <- paste(var_name, "==", first_level)
    }
  } else {
    # Default: use first level
    first_level <- levels(as.factor(var_data))[1]
    indicator <- var_data == first_level
    definition <- paste(var_name, "==", first_level)
  }

  if (verbose) {
    cat("  ", definition, "\n")
    cat("    Proportion in subgroup:",
        round(mean(indicator, na.rm = TRUE), 3), "\n")
  }

  return(list(indicator = indicator, definition = definition))
}


#' Process Cutpoint Specification for Subgroup Definition
#' @keywords internal
process_cutpoint <- function(var_data, cut_spec, var_name = "",
                             verbose = FALSE) {

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
      cutpoint <- NULL
      if (!is.null(cut_spec$value)) {
        cutpoint <- cut_spec$value
      } else if (!is.null(cut_spec$quantile)) {
        cutpoint <- quantile(var_data, probs = cut_spec$quantile,
                             na.rm = TRUE)
      } else if (!is.null(cut_spec$fun)) {
        cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
      }

      if (is.null(cutpoint)) {
        stop(paste("For 'greater' type, must provide 'value', 'quantile', or 'fun' for variable:",
                   var_name))
      }

      return(var_data > cutpoint)

    } else if (cut_type == "multiple") {
      # Multiple cutpoints (var in specified values)
      return(var_data %in% cut_spec$values)

    } else if (cut_type == "custom") {
      # Custom function that returns logical vector
      return(cut_spec$fun(var_data))

    } else {
      stop(paste("Unknown cutpoint type:", cut_type, "for variable:",
                 var_name))
    }
  }

  # Default: if no specification, use median
  if (is.null(cut_spec)) {
    if (verbose) cat("  Using median as default cutpoint for", var_name, "\n")
    return(var_data <= median(var_data, na.rm = TRUE))
  }

  stop(paste("Invalid cutpoint specification for variable:", var_name))
}


#' Format Continuous Variable Definition for Display
#' @keywords internal
format_continuous_definition <- function(var_data, cut_spec, var_name) {

  if (is.numeric(cut_spec) && length(cut_spec) == 1) {
    actual_cutpoint <- cut_spec
    return(paste(var_name, "<=", actual_cutpoint))

  } else if (is.list(cut_spec)) {
    if (cut_spec$type == "quantile") {
      actual_cutpoint <- quantile(var_data, probs = cut_spec$value,
                                  na.rm = TRUE)
      return(paste(var_name, "<=", round(actual_cutpoint, 2),
                   "(", cut_spec$value * 100, "th percentile)"))

    } else if (cut_spec$type == "function") {
      actual_cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
      fun_name <- deparse(substitute(cut_spec$fun))
      return(paste(var_name, "<=", round(actual_cutpoint, 2),
                   "(", fun_name, ")"))

    } else if (cut_spec$type == "range") {
      return(paste(cut_spec$min, "<=", var_name, "<=", cut_spec$max))

    } else if (cut_spec$type == "greater") {
      actual_cutpoint <- NULL
      if (!is.null(cut_spec$value)) {
        actual_cutpoint <- cut_spec$value
      } else if (!is.null(cut_spec$quantile)) {
        actual_cutpoint <- quantile(var_data, probs = cut_spec$quantile,
                                    na.rm = TRUE)
      } else if (!is.null(cut_spec$fun)) {
        actual_cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
      }

      if (!is.null(actual_cutpoint)) {
        return(paste(var_name, ">", round(actual_cutpoint, 2)))
      } else {
        return(paste(var_name, "> (unspecified cutpoint)"))
      }

    } else if (cut_spec$type == "custom") {
      return(paste(var_name, "(custom function)"))
    }
  }

  return(paste(var_name, "(unknown specification)"))
}




#' Generate Super Population and Calculate Linear Predictors
#' @keywords internal
generate_super_population <- function(df_work, n_super, draw_treatment,
                                      gamma, b0, mu, tau, verbose,
                                      spline_info = NULL) {  # ADD spline_info

  # Sample with replacement to create super population
  idx_sample <- sample(1:nrow(df_work), size = n_super, replace = TRUE)
  df_super <- df_work[idx_sample, ]

  if (verbose) {
    cat("\nOverall subgroup in super-population (all conditions met):\n")
    cat("  Proportion:", round(mean(df_super$flag_harm), 3), "\n")
  }

  # Handle treatment assignment
  if (draw_treatment) {
    # Assign treatment simple 1/2
    n_treat <- round(n_super / 2)
    n_control <- n_super - n_treat
    df_super$treat[1:n_treat] <- 1
    df_super$treat[(n_treat + 1):n_super] <- 0

    # If spline is present, need to recalculate spline interaction terms
    if (!is.null(spline_info)) {
      spline_var <- spline_info$var
      knot <- spline_info$knot

      # Recreate spline terms with new treatment assignment
      df_super[[paste0(spline_var, "_treat")]] <-
        df_super[[spline_var]] * df_super$treat

      df_super[[paste0(spline_var, "_k_treat")]] <-
        df_super[[paste0(spline_var, "_k")]] * df_super$treat
    }
  }

  # Reset IDs
  df_super$id <- 1:n_super

  # Calculate linear predictors for potential outcomes
  covariate_cols <- grep("^z_", names(df_super), value = TRUE)
  df_super <- calculate_linear_predictors(df_super, covariate_cols, gamma, b0,
                                          spline_info)  # PASS spline_info

  return(df_super)
}



#' Calculate Linear Predictors for Potential Outcomes
#' @keywords internal
calculate_linear_predictors <- function(df_super, covariate_cols, gamma, b0,
                                        spline_info = NULL) {

  # Separate spline interaction terms from regular covariates
  if (!is.null(spline_info)) {
    spline_var <- spline_info$var
    spline_interaction_terms <- c(
      paste0(spline_var, "_treat"),
      paste0(spline_var, "_k"),
      paste0(spline_var, "_k_treat")
    )
    # Base covariates (excluding spline interactions)
    base_covariate_cols <- setdiff(covariate_cols, spline_interaction_terms)
  } else {
    base_covariate_cols <- covariate_cols
    spline_interaction_terms <- character(0)
  }

  # Check if subgroup interaction term exists
  has_interaction <- "treat_harm" %in% names(gamma)

  if (has_interaction) {
    # Recalculate interaction for super population
    df_super$treat_harm <- df_super$treat * df_super$flag_harm
  }

  # ==========================================================================
  # Build Design Matrices for Observed Data
  # ==========================================================================

  # Start with base covariates and treatment
  X_super <- as.matrix(df_super[, c("treat", base_covariate_cols)])

  # Add spline terms if present (these were created in df_super already)
  if (!is.null(spline_info)) {
    for (spline_term in spline_interaction_terms) {
      if (spline_term %in% names(df_super)) {
        X_super <- cbind(X_super, df_super[[spline_term]])
        colnames(X_super)[ncol(X_super)] <- spline_term
      }
    }
  }

  # Add subgroup interaction if present
  if (has_interaction) {
    X_super <- cbind(X_super, treat_harm = df_super$treat_harm)
  }

  # ==========================================================================
  # Build Design Matrices for Potential Outcomes Under Treatment (A=1)
  # ==========================================================================

  X_treat <- as.matrix(df_super[, c("treat", base_covariate_cols)])
  X_treat[, "treat"] <- 1  # Set treatment to 1

  # Add spline terms for treatment arm
  if (!is.null(spline_info)) {
    spline_var <- spline_info$var

    # For A=1: spline_var_treat = spline_var * 1 = spline_var
    if (paste0(spline_var, "_treat") %in% names(df_super)) {
      X_treat <- cbind(X_treat,
                       df_super[[spline_var]])  # spline_var * treat where treat=1
      colnames(X_treat)[ncol(X_treat)] <- paste0(spline_var, "_treat")
    }

    # spline_var_k stays the same (doesn't depend on treatment)
    if (paste0(spline_var, "_k") %in% names(df_super)) {
      X_treat <- cbind(X_treat,
                       df_super[[paste0(spline_var, "_k")]])
      colnames(X_treat)[ncol(X_treat)] <- paste0(spline_var, "_k")
    }

    # For A=1: spline_var_k_treat = spline_var_k * 1 = spline_var_k
    if (paste0(spline_var, "_k_treat") %in% names(df_super)) {
      X_treat <- cbind(X_treat,
                       df_super[[paste0(spline_var, "_k")]])  # spline_var_k * treat where treat=1
      colnames(X_treat)[ncol(X_treat)] <- paste0(spline_var, "_k_treat")
    }
  }

  # Add subgroup interaction for treatment arm
  if (has_interaction) {
    X_treat <- cbind(X_treat, treat_harm = df_super$flag_harm)  # treat=1 * flag_harm
  }

  # ==========================================================================
  # Build Design Matrices for Potential Outcomes Under Control (A=0)
  # ==========================================================================

  X_control <- as.matrix(df_super[, c("treat", base_covariate_cols)])
  X_control[, "treat"] <- 0  # Set treatment to 0

  # Add spline terms for control arm
  if (!is.null(spline_info)) {
    spline_var <- spline_info$var

    # For A=0: spline_var_treat = spline_var * 0 = 0
    if (paste0(spline_var, "_treat") %in% names(df_super)) {
      X_control <- cbind(X_control,
                         rep(0, nrow(df_super)))  # spline_var * treat where treat=0
      colnames(X_control)[ncol(X_control)] <- paste0(spline_var, "_treat")
    }

    # spline_var_k stays the same (doesn't depend on treatment)
    if (paste0(spline_var, "_k") %in% names(df_super)) {
      X_control <- cbind(X_control,
                         df_super[[paste0(spline_var, "_k")]])
      colnames(X_control)[ncol(X_control)] <- paste0(spline_var, "_k")
    }

    # For A=0: spline_var_k_treat = spline_var_k * 0 = 0
    if (paste0(spline_var, "_k_treat") %in% names(df_super)) {
      X_control <- cbind(X_control,
                         rep(0, nrow(df_super)))  # spline_var_k * treat where treat=0
      colnames(X_control)[ncol(X_control)] <- paste0(spline_var, "_k_treat")
    }
  }

  # Add subgroup interaction for control arm
  if (has_interaction) {
    X_control <- cbind(X_control, treat_harm = 0)  # treat=0 * flag_harm = 0
  }

  # ==========================================================================
  # Verify Column Alignment
  # ==========================================================================

  # Ensure all matrices have the same columns in the same order
  expected_cols <- names(gamma)

  if (!all(colnames(X_super) == expected_cols)) {
    stop("Column mismatch in X_super")
  }
  if (!all(colnames(X_treat) == expected_cols)) {
    stop("Column mismatch in X_treat")
  }
  if (!all(colnames(X_control) == expected_cols)) {
    stop("Column mismatch in X_control")
  }

  # ==========================================================================
  # Calculate Linear Predictors
  # ==========================================================================

  # Linear predictors per AFT (log-time scale)
  df_super$lin_pred_1 <- as.vector(X_treat %*% gamma)
  df_super$lin_pred_0 <- as.vector(X_control %*% gamma)
  df_super$lin_pred_obs <- as.vector(X_super %*% gamma)

  # PO log-hazards (excluding baseline hazard)
  # Used for calculating empirical version of CDEs (controlled direct effects)
  df_super$theta_0 <- as.vector(X_control %*% b0)
  df_super$theta_1 <- as.vector(X_treat %*% b0)
  df_super$loghr_po <- df_super$theta_1 - df_super$theta_0

  return(df_super)
}

#' Calculate Hazard Ratios from Potential Outcomes
#'
#' @param df_super Data frame with super population
#' @param n_super Size of super population
#' @param mu Intercept parameter
#' @param tau Scale parameter
#' @param model Model type ("alt" or "null")
#' @param verbose Logical for verbose output
#'
#' @return List of hazard ratios
#'
#' @importFrom survival coxph Surv
#' @importFrom stats rexp
#' @keywords internal
calculate_hazard_ratios <- function(df_super, n_super, mu, tau, model,
                                    verbose) {
  flag_harm <- NULL
  # Generate potential outcomes for HR calculation
  epsilon <- log(rexp(n_super))  # Extreme value distribution

  # Under treatment
  logT_1 <- mu + tau * epsilon + df_super$lin_pred_1
  T_1 <- exp(logT_1)

  # Under control
  logT_0 <- mu + tau * epsilon + df_super$lin_pred_0
  T_0 <- exp(logT_0)

  # Create data frame for Cox models
  df_temp <- data.frame(
    time = c(T_1, T_0),
    event = 1,
    treat = c(rep(1, n_super), rep(0, n_super)),
    flag_harm = rep(df_super$flag_harm, 2)
  )

  # Calculate overall HR
  hr_overall <- exp(coxph(Surv(time, event) ~ treat,
                          data = df_temp)$coefficients)

  # Calculate average hazard ratios (AHR)
  AHR <- with(df_super, exp(mean(loghr_po)))
  AHR_harm <- with(subset(df_super, flag_harm == 1), exp(mean(loghr_po)))
  AHR_no_harm <- with(subset(df_super, flag_harm == 0), exp(mean(loghr_po)))

  # Calculate controlled direct effects (CDE)
  # CDE(S) = mean(exp(theta_1[S])) / mean(exp(theta_0[S]))
  # Differs from AHR due to Jensen's inequality (theta-ddagger in paper)
  CDE <- mean(exp(df_super$theta_1)) / mean(exp(df_super$theta_0))

  if (sum(df_super$flag_harm) > 0) {
    in_harm <- df_super$flag_harm == 1
    CDE_harm <- mean(exp(df_super$theta_1[in_harm])) /
      mean(exp(df_super$theta_0[in_harm]))
    CDE_no_harm <- mean(exp(df_super$theta_1[!in_harm])) /
      mean(exp(df_super$theta_0[!in_harm]))
  } else {
    CDE_harm <- NA_real_
    CDE_no_harm <- CDE
  }

  hr_results <- list(
    overall     = hr_overall,
    AHR         = AHR,
    AHR_harm    = AHR_harm,
    AHR_no_harm = AHR_no_harm,
    CDE         = CDE,
    CDE_harm    = CDE_harm,
    CDE_no_harm = CDE_no_harm
  )

  # Calculate subgroup-specific HRs if applicable
  if (model == "alt" && sum(df_super$flag_harm) > 0) {
    hr_harm <- exp(coxph(Surv(time, event) ~ treat,
                         data = subset(df_temp, flag_harm == 1))$coefficients)
    hr_no_harm <- exp(coxph(Surv(time, event) ~ treat,
                            data = subset(df_temp, flag_harm == 0))$coefficients)

    hr_results$harm_subgroup <- hr_harm
    hr_results$no_harm_subgroup <- hr_no_harm

    if (verbose) {
      cat("\n=== Hazard Ratios (super popln) ===\n")
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

  return(hr_results)
}


#' Prepare Censoring Model Parameters
#' @keywords internal
prepare_censoring_model_legacy <- function(df_work, cens_type, cens_params,
                                    df_super) {

  cens_model <- NULL

  if(cens_type != "uniform"){
    covariate_cols <- grep("^zcens_", names(df_work), value = TRUE)

    if(length(covariate_cols) >= 1){
      X_cens <- as.matrix(df_work[, c("treat", covariate_cols)])
    } else {
      X_cens <- as.matrix(df_work[, c("treat")])
    }

    fit_cens1 <- survreg(Surv(y, 1 - event) ~ X_cens,
                         data = df_work, dist = "weibull")


    fit_cens2 <- survreg(Surv(y, 1 - event) ~ X_cens,
                         data = df_work, dist = "lognormal")

    fit_cens3 <- survreg(Surv(y, 1 - event) ~ 1,
                         data = df_work, dist = "weibull")


    fit_cens4 <- survreg(Surv(y, 1 - event) ~ 1,
                         data = df_work, dist = "lognormal")


    # Compare all 4 models
    comparison <- compare_multiple_survreg(
      fit_cens1, fit_cens2, fit_cens3, fit_cens4,
      model_names = c("Weibull", "LogNormal", "Weibull0", "LogNormal0"),
      verbose = TRUE
    )

    # Get the best model
    fit_cens <- get_best_survreg(comparison)

    mu_cens <- coef(fit_cens)[1]
    tau_cens <- fit_cens$scale
    gamma_cens <- coef(fit_cens)[-1]

    # Store censoring parameters
    cens_model <- list(
      mu = mu_cens,
      tau = tau_cens,
      gamma = gamma_cens,
      type = c(fit_cens$dist)
    )

    newdata <- df_super
    newdata$treat <- 0
    if(!all(newdata$treat == 0)) stop("Error in creating counterfactual setting treat := 0")
    df_super$lin_pred_cens_0  <- predict(fit_cens, newdata = newdata, type = "linear")

    newdata <- df_super
    newdata$treat <- 1
    if(!all(newdata$treat == 1)) stop("Error in creating counterfactual setting treat := 1")
    df_super$lin_pred_cens_1  <- predict(fit_cens, newdata = newdata, type = "linear")


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

  return(list(
    cens_model = cens_model,
    df_super = df_super
  ))
}


#' Assemble Final Results Object
#' @keywords internal
assemble_results <- function(df_super,
                             mu,
                             tau,
                             gamma,
                             b0,
                             cens_model,
                             subgroup_vars,
                             subgroup_cuts,
                             subgroup_definitions,
                             hr_results,
                             continuous_vars,
                             factor_vars,
                             model,
                             n_super,
                             seed,
                             spline_info = NULL) {  # DEFAULT = NULL

  # Model parameters
  model_params <- list(
    mu = mu,
    tau = tau,
    gamma = gamma,
    b0 = b0,
    censoring = cens_model,
    spline_info = spline_info  # Can be NULL
  )

  # Subgroup information
  subgroup_info <- list(
    vars = subgroup_vars,
    cuts = subgroup_cuts,
    definitions = subgroup_definitions,
    size = sum(df_super$flag_harm),
    proportion = mean(df_super$flag_harm)
  )

  # Get covariate column names
  covariate_cols <- grep("^z_", names(df_super), value = TRUE)

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


#' Compare Multiple Survival Regression Models
#'
#' Performs comprehensive comparison of multiple survreg models including
#' convergence checking, information criteria comparison, and model selection.
#'
#' @param ... survreg model objects to compare
#' @param model_names Optional character vector of model names
#' @param verbose Logical, whether to print detailed output (default: TRUE)
#' @param criteria Character vector of criteria to use ("AIC", "BIC", or both)
#'
#' @return A list of class "multi_survreg_comparison" containing:
#' \describe{
#'   \item{models}{Named list of input models}
#'   \item{convergence}{Convergence status for each model}
#'   \item{comparison}{Model comparison statistics}
#'   \item{rankings}{Model rankings by different criteria}
#'   \item{best_model}{Name of the best model}
#'   \item{recommendation}{Text recommendation}
#' }
#'
#' @importFrom stats vcov AIC BIC model.frame coef
#' @export
#' @examples
#' \dontrun{
#' fit1 <- survreg(Surv(time, status) ~ x, dist = "weibull")
#' fit2 <- survreg(Surv(time, status) ~ x, dist = "lognormal")
#' comparison <- compare_multiple_survreg(fit1, fit2)
#' }
#'
compare_multiple_survreg <- function(...,
                                     model_names = NULL,
                                     verbose = TRUE,
                                     criteria = c("AIC", "BIC")) {

  # Collect models
  models <- list(...)
  n_models <- length(models)

  if(n_models < 2) {
    stop("At least 2 models required for comparison")
  }

  # Check that all inputs are survreg objects
  if (!all(vapply(models, inherits, logical(1), what = "survreg"))) {
    stop("All inputs must be survreg model objects")
  }

  # Auto-generate names if not provided
  if(is.null(model_names)) {
    model_names <- vapply(models, function(x) {
      dist_name <- x$dist
      if (is.list(dist_name)) dist_name <- dist_name$name
      paste0(dist_name, "_model")
    }, character(1))
  } else if(length(model_names) != n_models) {
    stop("Number of model names must match number of models")
  }

  # Initialize results
  results <- list(
    models = setNames(models, model_names),
    model_names = model_names,
    n_models = n_models,
    convergence = list(),
    comparison = list(),
    rankings = list(),
    best_model = NULL,
    recommendation = NULL
  )

  # ---- CONVERGENCE CHECKING ----
  check_convergence_single <- function(fit) {
    conv_status <- list(
      converged = TRUE,
      issues = character(),
      iterations = NA,
      log_likelihood = NA,
      scale = NA
    )

    # Check for NA coefficients
    if(any(is.na(coef(fit)))) {
      conv_status$converged <- FALSE
      conv_status$issues <- c(conv_status$issues, "NA coefficients")
    }

    # Check variance-covariance matrix
    tryCatch({
      vc <- vcov(fit)
      se <- sqrt(diag(vc))
      if(any(is.na(se) | is.infinite(se))) {
        conv_status$converged <- FALSE
        conv_status$issues <- c(conv_status$issues, "Invalid standard errors")
      }
      # Check for negative variances
      if(any(diag(vc) < 0)) {
        conv_status$converged <- FALSE
        conv_status$issues <- c(conv_status$issues, "Negative variance estimates")
      }
    }, error = function(e) {
      conv_status$converged <- FALSE
      conv_status$issues <- c(conv_status$issues, "Cannot compute variance matrix")
    })

    # Check iterations
    if(!is.null(fit$iter)) {
      conv_status$iterations <- fit$iter[1]
      if(!is.null(fit$maxiter) && fit$iter[1] >= fit$maxiter) {
        conv_status$converged <- FALSE
        conv_status$issues <- c(conv_status$issues,
                                paste("Max iterations reached:", fit$iter[1]))
      }
    }

    # Check log-likelihood
    if(!is.null(fit$loglik)) {
      ll_final <- fit$loglik[length(fit$loglik)]
      conv_status$log_likelihood <- ll_final
      if(is.na(ll_final) || is.infinite(ll_final)) {
        conv_status$converged <- FALSE
        conv_status$issues <- c(conv_status$issues, "Invalid log-likelihood")
      }
    }

    # Check scale parameter
    if(!is.null(fit$scale)) {
      conv_status$scale <- fit$scale
      if(is.na(fit$scale) || fit$scale <= 0 || is.infinite(fit$scale)) {
        conv_status$converged <- FALSE
        conv_status$issues <- c(conv_status$issues, "Invalid scale parameter")
      }
    }

    return(conv_status)
  }

  # Check all models
  convergence_summary <- data.frame(
    Model = model_names,
    Converged = logical(n_models),
    Issues = character(n_models),
    Iterations = integer(n_models),
    LogLik = numeric(n_models),
    Scale = numeric(n_models),
    stringsAsFactors = FALSE
  )

  for(i in 1:n_models) {
    conv_check <- check_convergence_single(models[[i]])
    results$convergence[[model_names[i]]] <- conv_check

    convergence_summary$Converged[i] <- conv_check$converged
    convergence_summary$Issues[i] <- paste(conv_check$issues, collapse = "; ")
    convergence_summary$Iterations[i] <- ifelse(is.na(conv_check$iterations), 0, conv_check$iterations)
    convergence_summary$LogLik[i] <- conv_check$log_likelihood
    convergence_summary$Scale[i] <- conv_check$scale
  }

  results$convergence_summary <- convergence_summary
  n_converged <- sum(convergence_summary$Converged)

  # ---- MODEL COMPARISON ----

  if(n_converged >= 2) {
    # Only compare converged models
    converged_idx <- which(convergence_summary$Converged)
    converged_models <- models[converged_idx]
    converged_names <- model_names[converged_idx]

    # Initialize comparison data frame
    comparison_df <- data.frame(
      Model = converged_names,
      Distribution = character(length(converged_models)),
      Converged = TRUE,
      LogLik = numeric(length(converged_models)),
      AIC = numeric(length(converged_models)),
      BIC = numeric(length(converged_models)),
      df = integer(length(converged_models)),
      Scale = numeric(length(converged_models)),
      n_obs = integer(length(converged_models)),
      stringsAsFactors = FALSE
    )

    # Fill in comparison metrics
    for(i in seq_along(converged_models)) {
      fit <- converged_models[[i]]
      comparison_df$Distribution[i] <- if(is.list(fit$dist)) fit$dist$name else fit$dist
      comparison_df$LogLik[i] <- fit$loglik[length(fit$loglik)]
      comparison_df$AIC[i] <- AIC(fit)
      comparison_df$BIC[i] <- BIC(fit)
      comparison_df$df[i] <- length(coef(fit)) + 1  # +1 for scale
      comparison_df$Scale[i] <- fit$scale
      comparison_df$n_obs[i] <- nrow(model.frame(fit))
    }

    # Add delta values for each criterion
    comparison_df$delta_AIC <- comparison_df$AIC - min(comparison_df$AIC)
    comparison_df$delta_BIC <- comparison_df$BIC - min(comparison_df$BIC)

    # Add AIC weights (Burnham & Anderson)
    comparison_df$AIC_weight <- exp(-0.5 * comparison_df$delta_AIC) /
      sum(exp(-0.5 * comparison_df$delta_AIC))

    # Sort by AIC
    comparison_df <- comparison_df[order(comparison_df$AIC), ]

    # Create rankings
    rankings <- data.frame(
      Rank = 1:nrow(comparison_df),
      By_AIC = comparison_df$Model[order(comparison_df$AIC)],
      AIC_value = sort(comparison_df$AIC),
      By_BIC = comparison_df$Model[order(comparison_df$BIC)],
      BIC_value = sort(comparison_df$BIC),
      By_LogLik = comparison_df$Model[order(comparison_df$LogLik, decreasing = TRUE)],
      LogLik_value = sort(comparison_df$LogLik, decreasing = TRUE)
    )

    results$comparison$table <- comparison_df
    results$rankings <- rankings

    # Determine best models
    best_aic <- comparison_df$Model[which.min(comparison_df$AIC)]
    best_bic <- comparison_df$Model[which.min(comparison_df$BIC)]

    results$comparison$best_aic <- best_aic
    results$comparison$best_bic <- best_bic

    # Evidence strength for top model
    if(nrow(comparison_df) > 1) {
      delta_second_best <- sort(comparison_df$delta_AIC)[2]
      weight_top_model <- comparison_df$AIC_weight[1]

      if(delta_second_best < 2) {
        evidence_strength <- "No clear winner"
        evidence_detail <- "Multiple models have similar support"
      } else if(delta_second_best < 4) {
        evidence_strength <- "Weak evidence"
        evidence_detail <- paste0("Top model has ", round(weight_top_model * 100, 1), "% of AIC weight")
      } else if(delta_second_best < 7) {
        evidence_strength <- "Moderate evidence"
        evidence_detail <- paste0("Top model has ", round(weight_top_model * 100, 1), "% of AIC weight")
      } else if(delta_second_best < 10) {
        evidence_strength <- "Strong evidence"
        evidence_detail <- paste0("Top model has ", round(weight_top_model * 100, 1), "% of AIC weight")
      } else {
        evidence_strength <- "Very strong evidence"
        evidence_detail <- paste0("Top model has ", round(weight_top_model * 100, 1), "% of AIC weight")
      }

      results$comparison$evidence_strength <- evidence_strength
      results$comparison$evidence_detail <- evidence_detail
    }

    # Model selection recommendation
    if(best_aic == best_bic) {
      results$best_model <- best_aic
      results$recommendation <- paste0(
        "Clear winner: ", best_aic, " (best by both AIC and BIC)\n",
        "  Evidence: ", evidence_strength, " - ", evidence_detail
      )
    } else {
      # Check if AIC winner has strong support
      aic_winner_weight <- comparison_df$AIC_weight[comparison_df$Model == best_aic]
      if(aic_winner_weight > 0.9) {
        results$best_model <- best_aic
        results$recommendation <- paste0(
          "Recommended: ", best_aic, " (overwhelming AIC support: ",
          round(aic_winner_weight * 100, 1), "%)\n",
          "  Note: BIC prefers ", best_bic, " but AIC evidence is very strong"
        )
      } else if(aic_winner_weight > 0.7) {
        results$best_model <- best_aic
        results$recommendation <- paste0(
          "Recommended: ", best_aic, " (strong AIC support: ",
          round(aic_winner_weight * 100, 1), "%)\n",
          "  Note: BIC prefers ", best_bic
        )
      } else {
        results$best_model <- best_aic
        results$recommendation <- paste0(
          "Mixed results: AIC prefers ", best_aic,
          " (", round(aic_winner_weight * 100, 1), "% weight)",
          ", BIC prefers ", best_bic, "\n",
          "  Consider model assumptions and purpose of analysis"
        )
      }
    }

  } else if(n_converged == 1) {
    results$best_model <- model_names[which(convergence_summary$Converged)]
    results$recommendation <- paste0(results$best_model, " (only converged model)")
  } else {
    results$best_model <- NA
    results$recommendation <- "No models converged successfully"
  }

  # ---- PRINT OUTPUT ----
  if(verbose) {
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat(" CENSORING MODEL SELECTION:\n")
    cat(" MULTIPLE SURVIVAL REGRESSION MODEL COMPARISONS\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")

    # Convergence summary
    cat("\nCONVERGENCE SUMMARY (", n_converged, "/", n_models, " converged):\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")

    for(i in 1:nrow(convergence_summary)) {
      cat("\n", convergence_summary$Model[i],
          " (", models[[i]]$dist, "):\n", sep = "")
      cat("  Status: ", ifelse(convergence_summary$Converged[i],
                               "\u2713 Converged", "\u2717 Failed"), "\n")
      if(!convergence_summary$Converged[i] && convergence_summary$Issues[i] != "") {
        cat("  Issues: ", convergence_summary$Issues[i], "\n")
      }
      if(convergence_summary$Converged[i]) {
        cat("  Iterations: ", convergence_summary$Iterations[i], "\n")
        cat("  Log-likelihood: ", round(convergence_summary$LogLik[i], 3), "\n")
        cat("  Scale parameter: ", round(convergence_summary$Scale[i], 4), "\n")
      }
    }

    # Model comparison
    if(n_converged >= 2) {
      cat("\nMODEL COMPARISON TABLE:\n")
      cat(paste(rep("-", 50), collapse = ""), "\n")
      # Print key columns with nice formatting
      print_df <- results$comparison$table[, c("Model", "AIC", "delta_AIC",
                                               "AIC_weight", "BIC", "delta_BIC")]
      print_df$AIC <- round(print_df$AIC, 2)
      print_df$delta_AIC <- round(print_df$delta_AIC, 2)
      print_df$AIC_weight <- round(print_df$AIC_weight, 3)
      print_df$BIC <- round(print_df$BIC, 2)
      print_df$delta_BIC <- round(print_df$delta_BIC, 2)
      print(print_df, row.names = FALSE)

      cat("\nMODEL RANKINGS:\n")
      cat(paste(rep("-", 50), collapse = ""), "\n")
      cat("By AIC: ", paste(rankings$By_AIC, collapse = " > "), "\n")
      cat("By BIC: ", paste(rankings$By_BIC, collapse = " > "), "\n")

      if(exists("evidence_strength")) {
        cat("\nEVIDENCE ASSESSMENT:\n")
        cat(paste(rep("-", 50), collapse = ""), "\n")
        cat("  Strength: ", evidence_strength, "\n")
        cat("  Details: ", evidence_detail, "\n")
      }
    }

    cat("\nFINAL RECOMMENDATION:\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat(" ", results$recommendation, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n\n")
  }

  class(results) <- c("multi_survreg_comparison", "list")
  invisible(results)
}




#' Get Best Model from Comparison
#'
#' Extracts the best fitting model object from a comparison result.
#' If no single best model can be determined, returns the Weibull model if
#' selected by either AIC or BIC. Defaults to Weibull0 model if no model
#' can be determined.
#'
#' @param comparison_result Output from compare_survreg_models or compare_multiple_survreg
#'
#' @return A survreg model object (defaults to Weibull0 model)
#' @keywords internal
#'
get_best_survreg <- function(comparison_result) {
  if(!inherits(comparison_result, c("multi_survreg_comparison", "survreg_comparison"))) {
    stop("Input must be output from model comparison function")
  }

  # Determine which model name to use
  best_name <- NULL

  # First, check if there's a clear best model
  if(!is.na(comparison_result$best_model)) {
    best_name <- comparison_result$best_model
  } else {
    # No best model determined - check if Weibull selected by AIC or BIC
    best_aic <- if("best_by_aic" %in% names(comparison_result)) {
      comparison_result$best_by_aic
    } else {
      NA
    }

    best_bic <- if("best_by_bic" %in% names(comparison_result)) {
      comparison_result$best_by_bic
    } else {
      NA
    }

    # Use Weibull if selected by either criterion
    if(!is.na(best_aic) && tolower(best_aic) == "weibull") {
      message("No single best model, but Weibull selected by AIC")
      best_name <- best_aic
    } else if(!is.na(best_bic) && tolower(best_bic) == "weibull") {
      message("No single best model, but Weibull selected by BIC")
      best_name <- best_bic
    } else {
      # Default to Weibull0 if no best model could be determined
      message("No best model could be determined - defaulting to Weibull0")
      best_name <- "Weibull0"
    }
  }

  # Return the actual model object
  if(!is.null(best_name) && best_name %in% names(comparison_result$models)) {
    return(comparison_result$models[[best_name]])
  } else {
    stop("Model '", best_name, "' not found in comparison_result$models")
  }
}

#' Print method for survreg_comparison objects
#'
#' @param x A survreg_comparison object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @exportS3Method
#'
print.multi_survreg_comparison <- function(x, ...) {
  cat("\nMultiple Survival Model Comparison (", x$n_models, " models)\n", sep = "")
  cat("Models compared: ", paste(x$model_names, collapse = ", "), "\n")
  cat("Models converged: ", sum(x$convergence_summary$Converged), "/", x$n_models, "\n")
  if(!is.na(x$best_model)) {
    cat("Best model: ", x$best_model, "\n")
  }
  invisible(x)
}
