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
#' @param treat.name Character or `NULL`. Name of the treatment term whose
#'   coefficient is the effect of interest.  When the formula is adjusted
#'   for covariates the fit has multiple coefficients; supplying
#'   `treat.name` extracts only that term's estimate and SE.  `NULL`
#'   (default) preserves the historical behaviour of returning every
#'   coefficient (correct for the treatment-only formula).
#' @return List with estimate and standard error.
#' @examples
#' \dontrun{
#' library(survival)
#' df <- data.frame(
#'   tte   = gbsg$rfstime / 30.4375,
#'   event = gbsg$status,
#'   treat = gbsg$hormon
#' )
#' formula <- build_cox_formula("tte", "event", "treat")
#' result  <- get_Cox_sg(df, cox.formula = formula)
#' exp(result$est_obs)  # hazard ratio
#' }
#' @importFrom survival coxph
#' @export
get_Cox_sg <- function(df_sg, cox.formula, est.loghr = TRUE, cox_initial = log(1),
                       treat.name = NULL) {

  # Validate inputs (keep your existing validation)
  names_tocheck <- all.vars(cox.formula)
  check <- unlist(lapply(names_tocheck, grep, names(df_sg), value = TRUE))
  check2 <- match(names_tocheck, check)
  if (sum(!is.na(check2)) != length(names_tocheck)) {
    stop("df_sg dataset NOT contain cox.formula variables")
  }

  # Fit model.  When the formula is adjusted (more than the treatment term),
  # a scalar `init` would be the wrong length for coxph; since `init = NULL`
  # is itself an error, the argument is omitted in that case so coxph uses
  # its default starting values.
  n_terms <- length(attr(stats::terms(cox.formula), "term.labels"))
  fit_args <- list(
    formula = cox.formula, data = df_sg,
    model = FALSE, x = FALSE, y = FALSE, robust = TRUE
  )
  if (n_terms <= 1L) fit_args$init <- cox_initial
  fit <- suppressWarnings(do.call(coxph, fit_args))

  # OPTIMIZATION: Call summary() once, cache result
  fit_sum <- summary(fit)
  coef_matrix <- fit_sum$coefficients

  # Restrict to the treatment coefficient when requested.  Adjusted models
  # carry multiple rows; without this the returned est/se would be vectors.
  if (!is.null(treat.name)) {
    if (!treat.name %in% rownames(coef_matrix)) {
      stop("treat.name '", treat.name, "' not among fitted Cox coefficients: ",
           paste(rownames(coef_matrix), collapse = ", "))
    }
    coef_matrix <- coef_matrix[treat.name, , drop = FALSE]
  }

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
#' Constructs a Cox model formula from variable names, optionally adjusted
#' for additional covariates.
#'
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @param adjust_covariates Character vector or `NULL`. Additional model terms
#'   appended to the right-hand side after the treatment term.  Terms are
#'   pasted verbatim, so survival modelling functions such as `strata()` and
#'   `pspline()` may be used directly, e.g. `adjust_covariates = "strata(x1)"`
#'   produces a stratified Cox model rather than including `x1` as a linear
#'   covariate.  `NULL` (default) reproduces the unadjusted, treatment-only
#'   formula.
#' @return An R formula object for Cox regression.  When `adjust_covariates`
#'   is supplied, the treatment term is always placed first on the
#'   right-hand side.
#' @examples
#' build_cox_formula("time_months", "event", "treat")
#' build_cox_formula("time_months", "event", "treat",
#'                   adjust_covariates = c("strata(site)", "age"))
#' @export

build_cox_formula <- function(outcome.name, event.name, treat.name,
                              adjust_covariates = NULL) {
  rhs <- treat.name
  adj <- .fs_adjust_terms(adjust_covariates)
  if (length(adj) > 0L) {
    rhs <- paste(c(treat.name, adj), collapse = " + ")
  }
  sf <- paste0("Surv(", outcome.name, ",", event.name, ") ~ ", rhs)
  as.formula(sf)
}


#' Normalise Adjustment Terms
#'
#' Coerces an `adjust_covariates` argument (character vector or list of
#' single-element character strings) into a plain character vector of model
#' terms, dropping `NULL`/empty entries.  Used wherever an adjusted Cox
#' right-hand side is assembled.
#'
#' @param adjust_covariates Character vector, list, or `NULL`.
#' @return Character vector of terms (possibly length zero).
#' @keywords internal
#' @noRd
.fs_adjust_terms <- function(adjust_covariates) {
  if (is.null(adjust_covariates)) return(character(0))
  terms <- unlist(adjust_covariates, use.names = FALSE)
  terms <- terms[!is.na(terms)]
  terms <- trimws(as.character(terms))
  terms[nzchar(terms)]
}


#' Extract Bare Variable Names From Adjustment Terms
#'
#' Returns the underlying column names referenced by a set of adjustment
#' terms, unwrapping survival modelling helpers.  For example,
#' `"strata(x1)"` yields `"x1"` and `"pspline(x3, df = 4)"` yields `"x3"`.
#' Used to determine which raw columns must be carried into the internal
#' scoring frame so that the adjusted formula resolves.
#'
#' @param adjust_covariates Character vector, list, or `NULL`.
#' @return Character vector of unique variable names (possibly length zero).
#' @keywords internal
#' @noRd
.fs_adjust_vars <- function(adjust_covariates) {
  terms <- .fs_adjust_terms(adjust_covariates)
  if (length(terms) == 0L) return(character(0))
  vars <- tryCatch(
    all.vars(stats::reformulate(terms)),
    error = function(e) character(0)
  )
  unique(vars)
}

#' Fit Cox Models for Subgroups
#'
#' Fits Cox models for two subgroups defined by treatment recommendation.
#'
#' @param df Data frame.
#' @param formula Cox model formula.
#' @param treat.name Character or `NULL`. Passed to \code{\link{get_Cox_sg}}
#'   to extract the treatment coefficient by name when \code{formula} is
#'   adjusted for covariates.  `NULL` (default) preserves the treatment-only
#'   behaviour.
#' @return List with HR and SE for each subgroup.
#' @examples
#' \dontrun{
#' library(survival)
#' df <- data.frame(
#'   tte            = gbsg$rfstime / 30.4375,
#'   event          = gbsg$status,
#'   treat          = gbsg$hormon,
#'   treat.recommend = as.integer(gbsg$er > 0)
#' )
#' formula <- build_cox_formula("tte", "event", "treat")
#' fit_cox_models(df, formula)
#' }
#' @export

fit_cox_models <- function(df, formula, treat.name = NULL) {
  fitH <- get_Cox_sg(df_sg = subset(df, treat.recommend == 0), cox.formula = formula,
                     treat.name = treat.name)
  fitHc <- get_Cox_sg(df_sg = subset(df, treat.recommend == 1), cox.formula = formula,
                      treat.name = treat.name)
  list(H_obs = fitH$est_obs, seH_obs = fitH$se_obs, Hc_obs = fitHc$est_obs, seHc_obs = fitHc$se_obs)
}


#' Estimate Subgroup Effect via Estimator Closure
#'
#' GLM counterpart to \code{\link{get_Cox_sg}}.  Calls the pre-built estimator closure
#' on a data slice and returns `est_obs` and `se_obs` in the same format
#' that `get_Cox_sg()` returns.
#'
#' @param df_sg Data frame.  The subgroup data slice.
#' @param estimator_fn Closure from `make_effect_estimator()`.
#'
#' @return List with `est_obs` (effect estimate) and `se_obs` (standard error).
#'
#' @keywords internal
fit_subgroup_effect <- function(df_sg, estimator_fn) {
  res <- tryCatch(
    estimator_fn(df_sg),
    error = function(e) NULL
  )

  if (is.null(res) || is.na(res$estimate)) {
    return(list(est_obs = NA_real_, se_obs = NA_real_))
  }

  list(est_obs = res$estimate, se_obs = res$se)
}


#' Fit Effect Models for Both Subgroups (GLM path)
#'
#' GLM counterpart to \code{\link{fit_cox_models}}.  Calls the estimator closure on
#' H (treat.recommend == 0) and Hc (treat.recommend == 1) subgroups.
#'
#' @param df Data frame with a `treat.recommend` column.
#' @param estimator_fn Closure from `make_effect_estimator()`.
#'
#' @return List with `H_obs`, `seH_obs`, `Hc_obs`, `seHc_obs`.
#'
#' @keywords internal
fit_effect_models <- function(df, estimator_fn) {
  fitH  <- fit_subgroup_effect(subset(df, treat.recommend == 0), estimator_fn)
  fitHc <- fit_subgroup_effect(subset(df, treat.recommend == 1), estimator_fn)
  list(
    H_obs    = fitH$est_obs,
    seH_obs  = fitH$se_obs,
    Hc_obs   = fitHc$est_obs,
    seHc_obs = fitHc$se_obs
  )
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
#' @examples
#' \dontrun{
#' library(survival)
#' cox_summary(
#'   Y     = gbsg$rfstime / 30.4375,
#'   E     = gbsg$status,
#'   Treat = gbsg$hormon
#' )
#' }
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


