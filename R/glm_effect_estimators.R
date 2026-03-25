# =============================================================================
# glm_effect_estimators.R
#
# Internal effect-estimator closures for all outcome types.
# Each factory returns a function of the form:
#
#   function(data_slice) -> list(estimate, se, converged, n0, n1,
#                                measure, method_used)
#
# where `estimate` is on the log scale (log-OR, log-RR, log-IRR) or identity
# scale (RD, IRD, MD) consistent with the chosen `effect_measure`.
#
# All functions are internal (@noRd) -- not exported.
# =============================================================================


# ---------------------------------------------------------------------------
# Top-level dispatcher
# ---------------------------------------------------------------------------

#' Make an effect-estimator closure
#'
#' Returns a closure `function(data_slice)` that estimates the within-subgroup
#' treatment effect for the requested outcome type and effect measure.  The
#' closure is called repeatedly inside the splitting-consistency loop and the
#' bootstrap worker, so it must be lightweight and handle small-sample edge
#' cases gracefully.
#'
#' @param outcome_type Character. One of `"survival"`, `"binary"`,
#'   `"continuous"`.
#' @param treat.name  Character. Name of the binary treatment column (0/1).
#' @param outcome.name Character. Name of the outcome column.
#' @param event.name  Character or `NULL`. Name of the event indicator
#'   (survival only).
#' @param offset.name Character or `NULL`. Name of the follow-up time column
#'   for rate-based measures (IRR, IRD).  The log of this column is used as
#'   the Poisson offset.
#' @param effect_measure Character or `NULL`.
#'   For binary: `"RD"` (default), `"OR"`, `"RR"`, `"IRR"`, `"IRD"`.
#'   For continuous: `"MD"` (default).
#'   For survival: `"HR"` (default, existing pipeline).
#'   `NULL` uses the default for the outcome type.
#' @param robust_se Logical. For rate-based measures (IRR, IRD) and modified
#'   Poisson (RR), use sandwich robust SE when the `sandwich` package is
#'   available.
#'   Default `TRUE`.
#' @param ... Additional arguments (currently unused; reserved for future
#'   covariate-adjusted estimators).
#'
#' @return A closure `function(data_slice)` returning a list with components
#'   `estimate`, `se`, `converged`, `n0`, `n1`, `measure`, `method_used`.
#'
#' @noRd
make_effect_estimator <- function(
    outcome_type,
    treat.name,
    outcome.name,
    event.name        = NULL,
    offset.name       = NULL,
    effect_measure    = NULL,
    robust_se         = TRUE,
    adjust_covariates = NULL,
    ...
) {
  outcome_type <- match.arg(
    outcome_type,
    choices = c("survival", "binary", "continuous")
  )

  # Resolve default effect measure
  if (is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      survival   = "HR",
      binary     = "RD",
      continuous = "MD"
    )
  }

  # Validate offset requirement for rate-based measures
  if (effect_measure %in% c("IRR", "IRD") && is.null(offset.name)) {
    stop(
      "effect_measure = '", effect_measure,
      "' requires `offset.name` (follow-up time column).",
      call. = FALSE
    )
  }

  switch(outcome_type,

    survival = {
      .make_cox_estimator(
        treat.name        = treat.name,
        outcome.name      = outcome.name,
        event.name        = event.name,
        adjust_covariates = adjust_covariates,
        ...
      )
    },

    binary = {
      effect_measure <- match.arg(
        effect_measure,
        choices = c("RD", "OR", "RR", "IRR", "IRD")
      )
      if (effect_measure %in% c("IRR", "IRD")) {
        .make_poisson_rate_estimator(
          treat.name        = treat.name,
          outcome.name      = outcome.name,
          offset.name       = offset.name,
          effect_measure    = effect_measure,
          robust_se         = robust_se,
          adjust_covariates = adjust_covariates,
          ...
        )
      } else {
        .make_glm_binary_estimator(
          treat.name        = treat.name,
          outcome.name      = outcome.name,
          effect_measure    = effect_measure,
          robust_se         = robust_se,
          adjust_covariates = adjust_covariates,
          ...
        )
      }
    },

    continuous = {
      .make_lm_estimator(
        treat.name        = treat.name,
        outcome.name      = outcome.name,
        adjust_covariates = adjust_covariates,
        ...
      )
    }
  )
}


# ---------------------------------------------------------------------------
# Survival estimator (Cox) -- wraps existing inline logic into closure form
# ---------------------------------------------------------------------------

#' @noRd
.make_cox_estimator <- function(treat.name, outcome.name, event.name,
                                adjust_covariates = NULL, ...) {

  force(treat.name)
  force(outcome.name)
  force(event.name)
  force(adjust_covariates)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    result <- tryCatch({
      surv_obj <- survival::Surv(
        time  = data_slice[[outcome.name]],
        event = data_slice[[event.name]]
      )
      rhs <- treat.name
      if (!is.null(adjust_covariates) && length(adjust_covariates) > 0) {
        rhs <- paste(c(treat.name, adjust_covariates), collapse = " + ")
      }
      fmla <- stats::as.formula(paste0("surv_obj ~ ", rhs))
      fit <- survival::coxph(fmla, data = data_slice)
      coef_val <- stats::coef(fit)[[treat.name]]
      se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
      list(
        estimate    = coef_val,
        se          = se_val,
        converged   = TRUE,
        n0          = n0,
        n1          = n1,
        measure     = "log-HR",
        method_used = "coxph"
      )
    },
    error = function(e) {
      list(
        estimate = NA_real_, se = NA_real_, converged = FALSE,
        n0 = n0, n1 = n1, measure = "log-HR", method_used = "coxph_failed"
      )
    })
    result
  }
}


# ---------------------------------------------------------------------------
# Binary GLM estimator (RD, OR, RR)
# ---------------------------------------------------------------------------

#' Make a binary-outcome GLM estimator closure
#'
#' Supports three non-rate effect measures:
#'   - `"RD"` : risk difference (identity-link binomial, with fallbacks)
#'   - `"OR"` : odds ratio (logistic regression); returns log-OR
#'   - `"RR"` : risk ratio (log-binomial with modified-Poisson fallback)
#'
#' For RD the estimate is on the identity scale.  For OR and RR the estimate
#' is on the log scale, consistent with the consistency-criterion comparisons.
#'
#' @noRd
.make_glm_binary_estimator <- function(
    treat.name,
    outcome.name,
    effect_measure    = "RD",
    robust_se         = TRUE,
    adjust_covariates = NULL,
    ...
) {
  force(treat.name)
  force(outcome.name)
  force(effect_measure)
  force(robust_se)
  force(adjust_covariates)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    fmla <- .build_adjusted_formula(outcome.name, treat.name,
                                    adjust_covariates)

    result <- switch(effect_measure,

      # ---- Risk Difference (primary) ----------------------------------------
      RD = {
        .estimate_rd(data_slice, fmla, treat.name, outcome.name, n0, n1,
                     adjust_covariates)
      },

      # ---- Odds Ratio -------------------------------------------------------
      OR = {
        tryCatch({
          fit <- stats::glm(
            fmla,
            data   = data_slice,
            family = stats::binomial(link = "logit")
          )
          coef_val <- stats::coef(fit)[[treat.name]]
          se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
          list(
            estimate    = coef_val,
            se          = se_val,
            converged   = fit$converged,
            n0          = n0,
            n1          = n1,
            measure     = "log-OR",
            method_used = "logistic"
          )
        },
        error = function(e) {
          list(
            estimate = NA_real_, se = NA_real_, converged = FALSE,
            n0 = n0, n1 = n1, measure = "log-OR",
            method_used = "logistic_failed"
          )
        })
      },

      # ---- Risk Ratio -------------------------------------------------------
      RR = {
        .estimate_rr(data_slice, fmla, treat.name, outcome.name,
                     n0, n1, robust_se)
      }
    )

    result
  }
}


# ---------------------------------------------------------------------------
# RD estimation with three-tier fallback
# ---------------------------------------------------------------------------

#' @noRd
.estimate_rd <- function(data_slice, fmla, treat.name, outcome.name, n0, n1,
                         adjust_covariates = NULL) {
  # Arm-level proportions (always available for fallback)
  y0 <- data_slice[[outcome.name]][data_slice[[treat.name]] == 0]
  y1 <- data_slice[[outcome.name]][data_slice[[treat.name]] == 1]
  p0 <- mean(y0, na.rm = TRUE)
  p1 <- mean(y1, na.rm = TRUE)

  n_adj <- if (!is.null(adjust_covariates)) length(adjust_covariates) else 0L

  # Tier 1: Identity-link binomial GLM
  tier1 <- tryCatch({
    start_vals <- c(p0, p1 - p0, rep(0, n_adj))
    fit <- stats::glm(
      fmla,
      data   = data_slice,
      family = stats::binomial(link = "identity"),
      start  = start_vals
    )
    if (!fit$converged) stop("identity-link binomial did not converge")
    coef_val <- stats::coef(fit)[[treat.name]]
    se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
    list(
      estimate    = coef_val,
      se          = se_val,
      converged   = TRUE,
      n0          = n0,
      n1          = n1,
      measure     = "RD",
      method_used = if (n_adj > 0) "identity_binomial_adj" else "identity_binomial"
    )
  },
  error   = function(e) NULL,
  warning = function(w) {
    if (grepl("converge", conditionMessage(w), ignore.case = TRUE)) {
      return(NULL)
    }
    invokeRestart("muffleWarning")
  })

  if (!is.null(tier1)) return(tier1)

  # Tier 2: G-computation from logistic regression
  tier2 <- tryCatch({
    fit_logit <- stats::glm(
      fmla,
      data   = data_slice,
      family = stats::binomial(link = "logit")
    )
    if (!fit_logit$converged) stop("logistic did not converge")

    # G-computation: predict under each treatment arm, average over subjects
    nd0 <- nd1 <- data_slice
    nd0[[treat.name]] <- 0L
    nd1[[treat.name]] <- 1L
    pred0 <- stats::predict(fit_logit, newdata = nd0, type = "response")
    pred1 <- stats::predict(fit_logit, newdata = nd1, type = "response")
    rd_est <- mean(pred1) - mean(pred0)

    # SE via coefficient SE as approximation (exact delta method is complex
    # with covariates; this is conservative for subgroup estimation)
    se_treat <- sqrt(diag(stats::vcov(fit_logit)))[[treat.name]]
    # Convert logistic SE to RD scale using average marginal effect
    avg_deriv <- mean(pred1 * (1 - pred1))
    se_rd <- abs(avg_deriv) * se_treat

    list(
      estimate    = rd_est,
      se          = se_rd,
      converged   = TRUE,
      n0          = n0,
      n1          = n1,
      measure     = "RD",
      method_used = if (n_adj > 0) "logistic_gcomp_adj" else "logistic_margins"
    )
  },
  error = function(e) NULL)

  if (!is.null(tier2)) return(tier2)

  # Tier 3: Raw proportions
  rd_raw <- p1 - p0
  se_raw <- sqrt(
    p1 * (1 - p1) / max(n1, 1L) + p0 * (1 - p0) / max(n0, 1L)
  )
  list(
    estimate    = rd_raw,
    se          = se_raw,
    converged   = FALSE,
    n0          = n0,
    n1          = n1,
    measure     = "RD",
    method_used = "raw_means"
  )
}


# ---------------------------------------------------------------------------
# RR estimation: log-binomial -> modified Poisson fallback
# ---------------------------------------------------------------------------

#' @noRd
.estimate_rr <- function(data_slice, fmla, treat.name, outcome.name,
                         n0, n1, robust_se) {

  # Attempt log-binomial first
  fit_rr <- tryCatch({
    p_start <- mean(data_slice[[outcome.name]], na.rm = TRUE)
    fit <- stats::glm(
      fmla,
      data   = data_slice,
      family = stats::binomial(link = "log"),
      start  = c(log(max(p_start, 0.01)), 0)
    )
    if (!fit$converged) stop("log-binomial did not converge")
    coef_val <- stats::coef(fit)[[treat.name]]
    se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
    list(
      estimate    = coef_val,
      se          = se_val,
      converged   = TRUE,
      n0          = n0,
      n1          = n1,
      measure     = "log-RR",
      method_used = "log_binomial"
    )
  },
  error   = function(e) NULL,
  warning = function(w) {
    if (grepl("converge", conditionMessage(w), ignore.case = TRUE)) {
      return(NULL)
    }
    invokeRestart("muffleWarning")
  })

  if (!is.null(fit_rr)) return(fit_rr)

  # Fallback: modified Poisson (Zou 2004)
  tryCatch({
    fit <- stats::glm(
      fmla,
      data   = data_slice,
      family = stats::poisson(link = "log")
    )
    coef_val <- stats::coef(fit)[[treat.name]]
    se_val   <- .robust_or_model_se(fit, treat.name, robust_se)
    list(
      estimate    = coef_val,
      se          = se_val,
      converged   = TRUE,
      n0          = n0,
      n1          = n1,
      measure     = "log-RR",
      method_used = "modified_poisson"
    )
  },
  error = function(e) {
    list(
      estimate = NA_real_, se = NA_real_, converged = FALSE,
      n0 = n0, n1 = n1, measure = "log-RR",
      method_used = "rr_failed"
    )
  })
}


# ---------------------------------------------------------------------------
# Poisson rate estimator (IRR and IRD)
# ---------------------------------------------------------------------------

#' Make a Poisson rate-model estimator closure
#'
#' Fits `glm(event ~ treat, family = poisson(link = "log"),
#'           offset = log(time))` and returns:
#'
#'   - For `"IRR"`: the log incidence rate ratio (treatment coefficient)
#'   - For `"IRD"`: the incidence rate difference computed from predicted rates
#'
#' Under constant baseline hazard (exponential survival), `exp(log-IRR)`
#' equals the Cox hazard ratio, providing a sanity check against the
#' existing survival pipeline.
#'
#' @noRd
.make_poisson_rate_estimator <- function(
    treat.name,
    outcome.name,
    offset.name,
    effect_measure    = "IRR",
    robust_se         = TRUE,
    adjust_covariates = NULL,
    ...
) {
  force(treat.name)
  force(outcome.name)
  force(offset.name)
  force(effect_measure)
  force(robust_se)
  force(adjust_covariates)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    # Validate offset column
    time_vec <- data_slice[[offset.name]]
    if (any(time_vec <= 0, na.rm = TRUE)) {
      time_vec <- pmax(time_vec, .Machine$double.eps)
    }

    fmla <- .build_adjusted_formula(outcome.name, treat.name,
                                    adjust_covariates)

    result <- tryCatch({
      # Add log-offset as a column so glm can resolve it in model.frame
      data_slice$.log_offset <- log(time_vec)
      fmla_off <- stats::update(fmla, . ~ . + offset(.log_offset))

      fit <- stats::glm(
        fmla_off,
        data   = data_slice,
        family = stats::poisson(link = "log")
      )

      if (effect_measure == "IRR") {
        # ---- Incidence Rate Ratio (log scale) ----
        coef_val <- stats::coef(fit)[[treat.name]]
        se_val   <- .robust_or_model_se(fit, treat.name, robust_se)
        list(
          estimate    = coef_val,
          se          = se_val,
          converged   = fit$converged,
          n0          = n0,
          n1          = n1,
          measure     = "log-IRR",
          method_used = "poisson_offset"
        )
      } else {
        # ---- Incidence Rate Difference (identity scale) ----
        # Predict rates for each arm at mean offset
        beta <- stats::coef(fit)
        V    <- stats::vcov(fit)
        mean_log_t <- mean(log(time_vec), na.rm = TRUE)

        # Rate in control:   lambda_0 = exp(beta_0 + mean_log_t)
        # Rate in treatment:  lambda_1 = exp(beta_0 + beta_1 + mean_log_t)
        lam0 <- exp(beta[[1L]] + mean_log_t)
        lam1 <- exp(beta[[1L]] + beta[[treat.name]] + mean_log_t)
        ird  <- lam1 - lam0

        # Delta-method SE: d(lam1 - lam0)/d(beta)
        grad <- c(lam1 - lam0, lam1)
        se_ird <- sqrt(as.numeric(t(grad) %*% V %*% grad))

        list(
          estimate    = ird,
          se          = se_ird,
          converged   = fit$converged,
          n0          = n0,
          n1          = n1,
          measure     = "IRD",
          method_used = "poisson_offset_delta"
        )
      }
    },
    error = function(e) {
      measure_label <- if (effect_measure == "IRR") "log-IRR" else "IRD"
      list(
        estimate = NA_real_, se = NA_real_, converged = FALSE,
        n0 = n0, n1 = n1, measure = measure_label,
        method_used = "poisson_failed"
      )
    })

    result
  }
}


# ---------------------------------------------------------------------------
# Continuous outcome estimator (linear regression / mean difference)
# ---------------------------------------------------------------------------

#' Make a continuous-outcome estimator closure
#'
#' Returns the within-subgroup OLS treatment coefficient (mean difference,
#' experimental minus control) and its standard error.
#'
#' @noRd
.make_lm_estimator <- function(treat.name, outcome.name,
                               adjust_covariates = NULL, ...) {

  force(treat.name)
  force(outcome.name)
  force(adjust_covariates)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    fmla <- .build_adjusted_formula(outcome.name, treat.name,
                                    adjust_covariates)

    result <- tryCatch({
      fit      <- stats::lm(fmla, data = data_slice)
      coef_val <- stats::coef(fit)[[treat.name]]
      se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
      list(
        estimate    = coef_val,
        se          = se_val,
        converged   = TRUE,
        n0          = n0,
        n1          = n1,
        measure     = "MD",
        method_used = "ols"
      )
    },
    error = function(e) {
      list(
        estimate = NA_real_, se = NA_real_, converged = FALSE,
        n0 = n0, n1 = n1, measure = "MD",
        method_used = "ols_failed"
      )
    })

    result
  }
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

#' Build model formula with optional covariate adjustment
#'
#' Constructs a formula string of the form `outcome ~ treat + cov1 + cov2`.
#' When `adjust_covariates` is NULL or empty, returns `outcome ~ treat`.
#'
#' @noRd
.build_adjusted_formula <- function(outcome.name, treat.name,
                                    adjust_covariates = NULL) {
  rhs <- treat.name
  if (!is.null(adjust_covariates) && length(adjust_covariates) > 0) {
    rhs <- paste(c(treat.name, adjust_covariates), collapse = " + ")
  }
  stats::as.formula(paste0(outcome.name, " ~ ", rhs))
}


#' Use sandwich robust SE if available; else model-based
#' @noRd
.robust_or_model_se <- function(fit, treat.name, robust_se = TRUE) {
  if (robust_se) {
    se_val <- tryCatch({
      if (!requireNamespace("sandwich", quietly = TRUE)) stop("no sandwich")
      sqrt(sandwich::sandwich(fit)[treat.name, treat.name])
    },
    error = function(e) {
      sqrt(diag(stats::vcov(fit)))[[treat.name]]
    })
  } else {
    se_val <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
  }
  se_val
}


#' Exponentiate an estimate for log-scale measures; pass through identity-scale
#'
#' @param estimate Numeric scalar.
#' @param measure  Character. One of `"log-HR"`, `"log-OR"`, `"log-RR"`,
#'   `"log-IRR"`, `"RD"`, `"IRD"`, `"MD"`.
#'
#' @return Numeric scalar on the natural (exponentiated or identity) scale.
#'
#' @noRd
back_transform_estimate <- function(estimate, measure) {
  if (grepl("^log-", measure)) {
    exp(estimate)
  } else {
    estimate
  }
}


#' Is this effect measure on the log scale?
#' @noRd
is_log_scale <- function(measure) {
  measure %in% c("log-HR", "log-OR", "log-RR", "log-IRR")
}


#' Map user-facing effect_measure to internal measure label
#' @noRd
effect_measure_to_label <- function(effect_measure) {
  switch(effect_measure,
    HR  = "log-HR",
    OR  = "log-OR",
    RR  = "log-RR",
    IRR = "log-IRR",
    RD  = "RD",
    IRD = "IRD",
    MD  = "MD",
    stop("Unknown effect_measure: ", effect_measure, call. = FALSE)
  )
}


#' Map user-facing effect_measure to consistency comparison direction
#'
#' Returns a function `f(estimate, threshold)` that evaluates `TRUE` when
#' the estimate exceeds the threshold in the direction consistent with harm.
#'
#' For log-scale measures: `estimate >= log(threshold)`
#' For identity-scale measures: `estimate >= threshold`
#'
#' @noRd
make_consistency_comparator <- function(effect_measure) {
  if (effect_measure %in% c("HR", "OR", "RR", "IRR")) {
    # Log-scale: compare log(estimate) >= log(threshold)
    function(estimate, threshold) estimate >= log(threshold)
  } else {
    # Identity scale (RD, IRD, MD): compare estimate >= threshold
    function(estimate, threshold) estimate >= threshold
  }
}


# =============================================================================
# PROPENSITY SCORE ESTIMATION MODULE
# =============================================================================

#' Estimate Propensity Scores
#'
#' Computes P(W = 1 | X) using the requested method, then constructs
#' the inverse propensity score (IPS) covariate following the doubly-robust
#' methodology of Bang & Robins (2005).  The IPS covariate takes the value
#' `1/e(x)` for treated subjects and `1/(1 - e(x))` for controls.
#' Including this as a covariate in the outcome model yields a
#' doubly-robust treatment effect estimate: consistent if either the
#' propensity model or the outcome model is correctly specified.
#'
#' @param data Data frame.
#' @param treat.name Character.  Name of binary treatment column (0/1).
#' @param confounders.name Character vector.  Covariate names.
#' @param method Character.  One of `"grf"` (default), `"lasso"`,
#'   `"logistic"`, or `"none"`.
#' @param grf_forest Optional.  A fitted `causal_forest` or
#'   `causal_survival_forest` object from which to extract `W.hat`.
#'   If provided and `method = "grf"`, the PS is taken directly from
#'   the forest rather than re-estimating.
#' @param seed Integer.  Random seed for GRF.
#' @param trim Numeric vector of length 2.  Quantile bounds for PS trimming
#'   (e.g., `c(0.01, 0.99)`).  Default: `c(0.025, 0.975)`.
#'
#' @return A list with:
#'   \item{data}{Data frame with `ips_covar` column appended (the Bang &
#'     Robins inverse PS covariate) and `ps_hat` column (the raw PS).}
#'   \item{ps_hat}{Numeric vector of estimated propensity scores.}
#'   \item{ips_covar}{Numeric vector of inverse PS covariates.}
#'   \item{method}{Character.  Method actually used.}
#'   \item{trimmed}{Integer.  Number of observations trimmed.}
#'
#' @references
#' Bang H, Robins JM (2005). "Doubly Robust Estimation in Missing Data
#' and Causal Inference Models." *Biometrics* 61(4): 962--973.
#'
#' @noRd
estimate_propensity_scores <- function(
    data,
    treat.name,
    confounders.name,
    method     = c("grf", "lasso", "logistic", "none"),
    grf_forest = NULL,
    seed       = 8316951L,
    trim       = c(0.025, 0.975)
) {
  method <- match.arg(method)

  if (method == "none") {
    return(list(
      data      = data,
      ps_hat    = NULL,
      ips_covar = NULL,
      method    = "none",
      trimmed   = 0L
    ))
  }

  W <- data[[treat.name]]
  X <- as.matrix(data[, confounders.name, drop = FALSE])

  ps_hat <- switch(method,
    grf      = .ps_grf(X, W, grf_forest, seed),
    lasso    = .ps_lasso(X, W),
    logistic = .ps_logistic(X, W, confounders.name, treat.name, data)
  )

  # Trim to avoid extreme weights / positivity violations
  n_trimmed <- 0L
  if (!is.null(trim) && length(trim) == 2L) {
    lo <- trim[1]
    hi <- trim[2]
    too_low  <- ps_hat < lo
    too_high <- ps_hat > hi
    n_trimmed <- sum(too_low | too_high, na.rm = TRUE)
    ps_hat <- pmax(pmin(ps_hat, hi), lo)
  }

  # -----------------------------------------------------------------------
  # Bang & Robins (2005) inverse propensity score covariate
  #
  # IPS_i = 1/e(X_i)         if W_i = 1 (treated)
  #       = 1/(1 - e(X_i))   if W_i = 0 (control)
  #
  # Including this as a covariate in the outcome model gives the
  # doubly-robust property: consistency requires only ONE of the
  # propensity or outcome model to be correctly specified.
  # -----------------------------------------------------------------------
  ips_covar <- ifelse(W == 1, 1 / ps_hat, 1 / (1 - ps_hat))

  data$ps_hat    <- ps_hat
  data$ips_covar <- ips_covar

  list(
    data      = data,
    ps_hat    = ps_hat,
    ips_covar = ips_covar,
    method    = method,
    trimmed   = n_trimmed
  )
}


#' GRF-based propensity score estimation
#'
#' Extracts W.hat from an existing causal forest, or fits a new
#' regression_forest on the treatment indicator.
#'
#' @noRd
.ps_grf <- function(X, W, grf_forest = NULL, seed = 8316951L) {

  # If a forest is already available, extract W.hat directly
  if (!is.null(grf_forest) && !is.null(grf_forest$W.hat)) {
    return(as.numeric(grf_forest$W.hat))
  }

  # Otherwise, fit a regression forest on W
  rf_w <- grf::regression_forest(X, W, seed = seed)
  ps <- predict(rf_w)$predictions
  as.numeric(ps)
}


#' LASSO logistic propensity score estimation
#'
#' Uses `glmnet` with `family = "binomial"` and `lambda.min` from CV.
#'
#' @noRd
.ps_lasso <- function(X, W) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' required for ps_method = 'lasso'.", call. = FALSE)
  }

  cv_fit <- glmnet::cv.glmnet(
    x      = X,
    y      = W,
    family  = "binomial",
    alpha   = 1,     # LASSO
    nfolds  = 5
  )
  ps <- as.numeric(stats::predict(
    cv_fit, newx = X, s = "lambda.min", type = "response"
  ))
  ps
}


#' Logistic regression propensity score estimation
#'
#' Simple logistic regression with main effects only.
#'
#' @noRd
.ps_logistic <- function(X, W, confounders.name, treat.name, data) {
  fmla <- stats::as.formula(
    paste0(treat.name, " ~ ", paste(confounders.name, collapse = " + "))
  )
  fit <- stats::glm(fmla, data = data, family = stats::binomial(link = "logit"))
  ps <- stats::predict(fit, type = "response")
  as.numeric(ps)
}
