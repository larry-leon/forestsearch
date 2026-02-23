#' Fit Cox Model for Subgroup
#'
#' Fits a Cox model for a subgroup and returns estimate and standard error.
#'
#' Function is utilized throughout codebase
#'
#' @param df_sg Data frame for subgroup.
#' @param cox_initial Optional pre-fitted Cox model object to use instead of
#'   fitting a new model. Default NULL
#' @param cox.formula Cox model formula.
#' @param est.loghr Logical. Is estimate on log(HR) scale?
#' @return List with estimate and standard error.
#' @importFrom survival coxph
#' @export
get_Cox_sg <- function(df_sg, cox.formula, est.loghr = TRUE, cox_initial = log(1)) {

  # Validate inputs (keep your existing validation)
  names_tocheck <- all.vars(cox.formula)
  check <- unlist(lapply(names_tocheck, grep, names(df_sg), value = TRUE))
  check2 <- match(names_tocheck, check)
  if (sum(!is.na(check2)) != length(names_tocheck)) {
    stop("df_sg dataset NOT contain cox.formula variables")
  }

  # Fit model
  fit <- suppressWarnings(
    coxph(cox.formula, data = df_sg, model = FALSE, x = FALSE, y = FALSE, robust = TRUE, init = cox_initial)
  )

  # OPTIMIZATION: Call summary() once, cache result
  fit_sum <- summary(fit)
  coef_matrix <- fit_sum$coefficients

  # Extract coefficients
  bhat <- coef_matrix[, "coef"]

  # Extract appropriate SE
  if (est.loghr) {
    est_obs <- bhat
    se_obs <- coef_matrix[, "robust se"]
  } else {
    est_obs <- exp(bhat)
    se_obs <- exp(bhat) * coef_matrix[, "robust se"]
  }

  return(list(est_obs = est_obs, se_obs = se_obs))
}


#' Build Cox Model Formula
#'
#' Constructs a Cox model formula from variable names.
#'
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @return An R formula object for Cox regression.
#' @export

build_cox_formula <- function(outcome.name, event.name, treat.name) {
  sf <- paste0("Surv(", outcome.name, ",", event.name, ") ~ ", treat.name)
  as.formula(sf)
}

#' Fit Cox Models for Subgroups
#'
#' Fits Cox models for two subgroups defined by treatment recommendation.
#'
#' @param df Data frame.
#' @param formula Cox model formula.
#' @return List with HR and SE for each subgroup.
#' @export

fit_cox_models <- function(df, formula) {
  fitH <- get_Cox_sg(df_sg = subset(df, treat.recommend == 0), cox.formula = formula)
  fitHc <- get_Cox_sg(df_sg = subset(df, treat.recommend == 1), cox.formula = formula)
  list(H_obs = fitH$est_obs, seH_obs = fitH$se_obs, Hc_obs = fitHc$est_obs, seHc_obs = fitHc$se_obs)
}


#' Cox model summary for subgroup
#'
#' Called in analyze_subgroup() <-- SG_tab_estimates
#'
#' Calculates hazard ratio and confidence interval for a subgroup using Cox regression.
#'
#' @param Y Numeric vector of outcome.
#' @param E Numeric vector of event indicators.
#' @param Treat Numeric vector of treatment indicators.
#' @param Strata Vector of strata (optional).
#' @return Character string with formatted HR and CI.
#' @importFrom survival coxph Surv
#' @keywords internal

cox_summary_legacy <- function(Y, E, Treat, Strata) {
  fit <- survival::coxph(survival::Surv(Y, E) ~ Treat + strata(Strata), robust = TRUE)
  hr <- summary(fit)$conf.int[c(1, 3, 4)]
  hrCI_format(hr)
}


#' Cox model summary for subgroup (OPTIMIZED)
#'
#' Called in analyze_subgroup() <-- SG_tab_estimates
#'
#' Calculates hazard ratio and confidence interval for a subgroup using Cox regression.
#' Optimized version with reduced overhead and better error handling.
#'
#' @param Y Numeric vector of outcome.
#' @param E Numeric vector of event indicators.
#' @param Treat Numeric vector of treatment indicators.
#' @param Strata Vector of strata (optional).
#' @param use_strata Logical. Whether to use strata in the model (default: TRUE if Strata provided).
#' @param return_format Character. "formatted" (default) or "numeric" for downstream use.
#' @return Character string with formatted HR and CI (or numeric vector if return_format="numeric").
#' @importFrom survival coxph Surv
#' @export

cox_summary <- function(Y, E, Treat, Strata = NULL,
                        use_strata = !is.null(Strata),
                        return_format = c("formatted", "numeric")) {

  return_format <- match.arg(return_format)

  # =========================================================================
  # OPTIMIZATION 1: Input validation (fail fast)
  # =========================================================================

  n <- length(Y)
  if (length(E) != n || length(Treat) != n) {
    stop("Y, E, and Treat must have the same length")
  }

  if (use_strata && !is.null(Strata) && length(Strata) != n) {
    stop("Strata must have the same length as Y")
  }

  # Quick check for sufficient events
  n_events <- sum(E)
  if (n_events < 2) {
    warning("Fewer than 2 events; returning NA")
    if (return_format == "formatted") {
      return("NA (NA, NA)")
    } else {
      return(c(HR = NA_real_, Lower = NA_real_, Upper = NA_real_))
    }
  }

  # Check treatment variation
  if (length(unique(Treat)) < 2) {
    warning("No variation in treatment; returning NA")
    if (return_format == "formatted") {
      return("NA (NA, NA)")
    } else {
      return(c(HR = NA_real_, Lower = NA_real_, Upper = NA_real_))
    }
  }

  # =========================================================================
  # OPTIMIZATION 2: Efficient formula construction
  # =========================================================================

  # Pre-build the formula based on whether strata is used
  if (use_strata && !is.null(Strata)) {
    # Build formula with strata only once
    fit <- tryCatch({
      survival::coxph(
        survival::Surv(Y, E) ~ Treat + strata(Strata),
        robust = TRUE,
        model = FALSE,   # Don't store model frame (saves memory)
        x = FALSE,       # Don't store design matrix
        y = FALSE        # Don't store response
      )
    }, error = function(e) {
      warning("Cox model failed: ", e$message)
      return(NULL)
    })
  } else {
    fit <- tryCatch({
      survival::coxph(
        survival::Surv(Y, E) ~ Treat,
        robust = TRUE,
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
    }, error = function(e) {
      warning("Cox model failed: ", e$message)
      return(NULL)
    })
  }

  # Handle fitting errors
  if (is.null(fit)) {
    if (return_format == "formatted") {
      return("NA (NA, NA)")
    } else {
      return(c(HR = NA_real_, Lower = NA_real_, Upper = NA_real_))
    }
  }

  # =========================================================================
  # OPTIMIZATION 3: Call summary() once and extract efficiently
  # =========================================================================

  # Call summary once and cache the result
  fit_summary <- summary(fit)

  # Extract confidence interval directly (no need for intermediate variables)
  conf_int <- fit_summary$conf.int

  # Handle edge case where conf.int might not exist
  if (is.null(conf_int) || nrow(conf_int) == 0) {
    warning("No confidence interval available from Cox model")
    if (return_format == "formatted") {
      return("NA (NA, NA)")
    } else {
      return(c(HR = NA_real_, Lower = NA_real_, Upper = NA_real_))
    }
  }

  # Extract HR and CI bounds (columns 1, 3, 4 of conf.int)
  hr <- conf_int[1, 1]    # HR
  lower <- conf_int[1, 3]  # Lower CI
  upper <- conf_int[1, 4]  # Upper CI

  # =========================================================================
  # OPTIMIZATION 4: Conditional formatting (only if needed)
  # =========================================================================

  if (return_format == "formatted") {
    return(hrCI_format(c(hr, lower, upper)))
  } else {
    return(c(HR = hr, Lower = lower, Upper = upper))
  }
}

#' Cox model summary for subgroup - vectorized version
#'
#' Efficiently processes multiple subgroups at once.
#' Useful when analyzing many subgroups (e.g., in cross-validation).
#'
#' @param data Data frame with columns for Y, E, Treat, and optionally Strata.
#' @param outcome_col Character. Name of outcome column.
#' @param event_col Character. Name of event column.
#' @param treat_col Character. Name of treatment column.
#' @param strata_col Character. Name of strata column (optional).
#' @param subgroup_col Character. Name of subgroup indicator column.
#' @param return_format Character. "formatted" or "numeric".
#'
#' @return Data frame with one row per subgroup and HR results.
#' @importFrom survival coxph Surv
#' @keywords internal

cox_summary_vectorized <- function(data,
                                   outcome_col,
                                   event_col,
                                   treat_col,
                                   strata_col = NULL,
                                   subgroup_col = "subgroup",
                                   return_format = c("formatted", "numeric")) {

  return_format <- match.arg(return_format)

  # Get unique subgroups
  subgroups <- unique(data[[subgroup_col]])
  n_subgroups <- length(subgroups)

  # Pre-allocate results
  results <- vector("list", n_subgroups)

  for (i in seq_along(subgroups)) {
    sg <- subgroups[i]

    # Subset once
    data_sg <- data[data[[subgroup_col]] == sg, ]

    # Extract vectors
    Y <- data_sg[[outcome_col]]
    E <- data_sg[[event_col]]
    Treat <- data_sg[[treat_col]]
    Strata <- if (!is.null(strata_col)) data_sg[[strata_col]] else NULL

    # Call optimized cox_summary
    result <- cox_summary(Y, E, Treat, Strata,
                          use_strata = !is.null(strata_col),
                          return_format = return_format)

    if (return_format == "formatted") {
      results[[i]] <- data.frame(
        subgroup = sg,
        n = nrow(data_sg),
        events = sum(E),
        HR_CI = result,
        stringsAsFactors = FALSE
      )
    } else {
      results[[i]] <- data.frame(
        subgroup = sg,
        n = nrow(data_sg),
        events = sum(E),
        HR = result["HR"],
        Lower = result["Lower"],
        Upper = result["Upper"],
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine results
  do.call(rbind, results)
}

#' Batch Cox summaries with caching
#'
#' For repeated calls with the same data structure but different subsets,
#' this version pre-processes the data structure once.
#'
#' @param Y Numeric vector of outcome (full dataset).
#' @param E Numeric vector of event indicators (full dataset).
#' @param Treat Numeric vector of treatment indicators (full dataset).
#' @param Strata Vector of strata (optional, full dataset).
#' @param subset_indices List of integer vectors, each defining a subset to analyze.
#' @param return_format Character. "formatted" or "numeric".
#'
#' @return List of results, one per subset.
#' @importFrom survival coxph Surv
#' @keywords internal

cox_summary_batch <- function(Y, E, Treat, Strata = NULL,
                              subset_indices,
                              return_format = c("formatted", "numeric")) {

  return_format <- match.arg(return_format)

  n_subsets <- length(subset_indices)
  results <- vector("list", n_subsets)

  # Pre-validate full dataset
  n <- length(Y)
  if (length(E) != n || length(Treat) != n) {
    stop("Y, E, and Treat must have the same length")
  }

  use_strata <- !is.null(Strata)
  if (use_strata && length(Strata) != n) {
    stop("Strata must have the same length as Y")
  }

  for (i in seq_along(subset_indices)) {
    idx <- subset_indices[[i]]

    # Extract subset
    Y_sub <- Y[idx]
    E_sub <- E[idx]
    Treat_sub <- Treat[idx]
    Strata_sub <- if (use_strata) Strata[idx] else NULL

    # Call optimized function
    results[[i]] <- cox_summary(
      Y_sub, E_sub, Treat_sub, Strata_sub,
      use_strata = use_strata,
      return_format = return_format
    )
  }

  return(results)
}


