# =============================================================================
# glm_adjusted_estimators.R
#
# Covariate-adjusted, propensity-score, and doubly-robust effect estimators.
# These are the adjusted counterparts to the unadjusted closures in
# glm_effect_estimators.R.
#
# Three tiers:
#   "regression"     — Outcome regression with g-computation (no PS model)
#   "ps_covariate"   — Bang & Robins (2005) DR: inv(PS) as covariate
#   "aipw"           — Funk & Davidian (2011) augmented IPW
#
# All functions are internal (@noRd) — not exported.
#
# References:
#   Bang H, Robins JM (2005). Biometrics 61:962–972.
#   Funk MJ et al. (2011). Am J Epidemiol 173:761–767.
# =============================================================================


# ---------------------------------------------------------------------------
# Outcome-regression estimator with g-computation
# Corresponds to delta.OR in BR.TE.fit()
# ---------------------------------------------------------------------------

#' Make a covariate-adjusted (regression) estimator closure
#'
#' Fits `outcome ~ treat + confounders` and uses g-computation
#' (standardisation) to estimate the risk difference, odds ratio, or other
#' effect measure.  The treatment effect is obtained by predicting under
#' each arm and averaging, NOT from the model coefficient — this is the
#' correct nonparametric-style estimand even when the model is nonlinear.
#'
#' @param treat.name Character. Treatment column name (binary 0/1).
#' @param outcome.name Character. Outcome column name.
#' @param confounders.name Character vector. Covariate column names.
#' @param effect_measure Character. `"RD"` (default), `"OR"`, `"RR"`.
#'   For `"RD"`, returns the g-computed risk difference.
#'   For `"OR"`, returns the g-computed log-odds ratio (logistic model).
#'   For `"RR"`, returns the g-computed log-risk ratio.
#' @param ... Additional arguments (currently unused).
#'
#' @return A closure `function(data_slice)`.
#'
#' @noRd
.make_regression_estimator <- function(
    treat.name,
    outcome.name,
    confounders.name,
    effect_measure = "RD",
    ...
) {
  force(treat.name)
  force(outcome.name)
  force(confounders.name)
  force(effect_measure)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    # Build outcome model formula: Y ~ treat + confounders
    rhs <- paste(c(treat.name, confounders.name), collapse = " + ")
    fmla <- stats::as.formula(paste0(outcome.name, " ~ ", rhs))

    result <- tryCatch({
      fit <- stats::glm(
        fmla,
        data   = data_slice,
        family = stats::binomial(link = "logit")
      )

      # G-computation: predict under each arm for ALL subjects
      nd1 <- nd0 <- data_slice
      nd1[[treat.name]] <- 1L
      nd0[[treat.name]] <- 0L
      mu1 <- stats::predict(fit, newdata = nd1, type = "response")
      mu0 <- stats::predict(fit, newdata = nd0, type = "response")

      .gcomp_result(mu1, mu0, effect_measure, n0, n1,
                    converged = fit$converged, method = "regression")
    },
    error = function(e) {
      .gcomp_fallback(data_slice, treat.name, outcome.name,
                      effect_measure, n0, n1, method = "regression_failed")
    })

    result
  }
}


# ---------------------------------------------------------------------------
# Bang & Robins PS-covariate DR estimator
# Corresponds to delta.DR (v1) in BR.TE.fit()
# ---------------------------------------------------------------------------

#' Make a Bang & Robins DR estimator closure
#'
#' Augments the outcome regression model with the inverse of the estimated
#' propensity score as an additional covariate:
#'
#'   `Y ~ treat + confounders + inv_ps`
#'
#' where `inv_ps = 1/pihat` for treated subjects and `inv_ps = 1/(1-pihat)`
#' for control subjects.  The estimate is doubly robust: consistent if either
#' the PS model or the outcome model is correctly specified.
#'
#' The `inv_ps` column (and optionally `.pihat`) must already exist in the
#' data slice — they are pre-computed in `forestsearch()` before entering
#' the consistency loop.
#'
#' @param treat.name Character. Treatment column name (binary 0/1).
#' @param outcome.name Character. Outcome column name.
#' @param confounders.name Character vector. Covariate column names.
#' @param ps.name Character. Column name of the pre-computed inverse
#'   propensity score.  Default `"inv_ps"`.
#' @param effect_measure Character. `"RD"` (default), `"OR"`, `"RR"`.
#' @param pihat.name Character. Column name of the raw propensity score
#'   (needed for counterfactual prediction).  Default `".pihat"`.
#' @param ... Additional arguments.
#'
#' @return A closure `function(data_slice)`.
#'
#' @noRd
.make_ps_covariate_estimator <- function(
    treat.name,
    outcome.name,
    confounders.name,
    ps.name        = "inv_ps",
    pihat.name     = ".pihat",
    effect_measure = "RD",
    ...
) {
  force(treat.name)
  force(outcome.name)
  force(confounders.name)
  force(ps.name)
  force(pihat.name)
  force(effect_measure)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    # Validate that PS columns exist
    if (!ps.name %in% names(data_slice)) {
      return(list(
        estimate = NA_real_, se = NA_real_, converged = FALSE,
        n0 = n0, n1 = n1,
        measure = effect_measure_to_label(effect_measure),
        method_used = "ps_covariate_missing_column"
      ))
    }

    # Build DR model formula: Y ~ treat + confounders + inv_ps
    rhs <- paste(c(treat.name, confounders.name, ps.name), collapse = " + ")
    fmla <- stats::as.formula(paste0(outcome.name, " ~ ", rhs))

    result <- tryCatch({
      fit <- stats::glm(
        fmla,
        data   = data_slice,
        family = stats::binomial(link = "logit")
      )

      # G-computation under each arm:
      # Set treat = 1, inv_ps = 1/pihat
      # Set treat = 0, inv_ps = 1/(1 - pihat)
      pihat <- data_slice[[pihat.name]]
      if (is.null(pihat)) {
        # Reconstruct from inv_ps if pihat column not available
        pihat <- ifelse(
          data_slice[[treat.name]] == 1L,
          1 / data_slice[[ps.name]],
          1 - 1 / data_slice[[ps.name]]
        )
      }

      nd1 <- data_slice
      nd1[[treat.name]] <- 1L
      nd1[[ps.name]]    <- 1 / pihat

      nd0 <- data_slice
      nd0[[treat.name]] <- 0L
      nd0[[ps.name]]    <- 1 / (1 - pihat)

      mu1 <- stats::predict(fit, newdata = nd1, type = "response")
      mu0 <- stats::predict(fit, newdata = nd0, type = "response")

      .gcomp_result(mu1, mu0, effect_measure, n0, n1,
                    converged = fit$converged, method = "ps_covariate")
    },
    error = function(e) {
      .gcomp_fallback(data_slice, treat.name, outcome.name,
                      effect_measure, n0, n1,
                      method = "ps_covariate_failed")
    })

    result
  }
}


# ---------------------------------------------------------------------------
# Funk & Davidian AIPW estimator
# Corresponds to delta.DR3 (v3) in BR.TE.fit()
# ---------------------------------------------------------------------------

#' Make an AIPW (augmented inverse probability weighted) estimator closure
#'
#' Implements the Funk & Davidian (2011) doubly-robust estimator using the
#' efficient influence function:
#'
#'   phi1_i = Z_i * Y_i / pihat_i  -  (Z_i - pihat_i) / pihat_i * muhat_1(L_i)
#'   phi0_i = (1-Z_i) * Y_i / (1-pihat_i)  +  (Z_i - pihat_i) / (1-pihat_i) * muhat_0(L_i)
#'   RD_AIPW = mean(phi1) - mean(phi0)
#'
#' This estimator achieves the semiparametric efficiency bound when both
#' models are correct.
#'
#' @param treat.name Character. Treatment column name.
#' @param outcome.name Character. Outcome column name.
#' @param confounders.name Character vector. Covariate column names.
#' @param pihat.name Character. Column name of pre-computed propensity score.
#'   Default `".pihat"`.
#' @param effect_measure Character. Currently only `"RD"` is implemented
#'   for AIPW.
#' @param ... Additional arguments.
#'
#' @return A closure `function(data_slice)`.
#'
#' @noRd
.make_aipw_estimator <- function(
    treat.name,
    outcome.name,
    confounders.name,
    pihat.name     = ".pihat",
    effect_measure = "RD",
    ...
) {
  force(treat.name)
  force(outcome.name)
  force(confounders.name)
  force(pihat.name)
  force(effect_measure)

  if (effect_measure != "RD") {
    warning(
      "AIPW estimator currently returns the risk difference (RD) ",
      "regardless of effect_measure = '", effect_measure, "'. ",
      "The RD is the natural scale for the AIPW influence function.",
      call. = FALSE
    )
  }

  function(data_slice) {
    n  <- nrow(data_slice)
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    if (!pihat.name %in% names(data_slice)) {
      return(list(
        estimate = NA_real_, se = NA_real_, converged = FALSE,
        n0 = n0, n1 = n1, measure = "RD",
        method_used = "aipw_missing_pihat"
      ))
    }

    result <- tryCatch({
      Z     <- data_slice[[treat.name]]
      Y     <- data_slice[[outcome.name]]
      pihat <- data_slice[[pihat.name]]

      # Fit outcome regression: separate models per arm
      # (corresponds to v4 in BR.TE.fit — recommended by Funk et al.)
      rhs_or <- paste(confounders.name, collapse = " + ")
      fmla_or <- stats::as.formula(paste0(outcome.name, " ~ ", rhs_or))

      data1 <- data_slice[Z == 1L, , drop = FALSE]
      data0 <- data_slice[Z == 0L, , drop = FALSE]

      fit1 <- stats::glm(fmla_or, data = data1,
                          family = stats::binomial(link = "logit"))
      fit0 <- stats::glm(fmla_or, data = data0,
                          family = stats::binomial(link = "logit"))

      # Predict potential outcomes for ALL subjects
      muhat_1 <- stats::predict(fit1, newdata = data_slice, type = "response")
      muhat_0 <- stats::predict(fit0, newdata = data_slice, type = "response")

      # AIPW influence function components
      # Treated arm: phi1_i = Z*Y/pihat - (Z - pihat)/pihat * muhat_1
      ipw_1  <- 1 / pihat
      term1a <- Z * Y * ipw_1
      term1b <- (Z - pihat) * ipw_1
      phi1   <- term1a - term1b * muhat_1

      # Control arm: phi0_i = (1-Z)*Y/(1-pihat) + (Z - pihat)/(1-pihat) * muhat_0
      ipw_0  <- 1 / (1 - pihat)
      term0a <- (1 - Z) * Y * ipw_0
      term0b <- (Z - pihat) * ipw_0
      phi0   <- term0a + term0b * muhat_0

      # AIPW risk difference
      rd_aipw <- mean(phi1) - mean(phi0)

      # Influence-function-based SE:
      # Var(RD) ≈ (1/n) * var(phi1 - phi0)
      psi     <- phi1 - phi0
      se_aipw <- sqrt(stats::var(psi) / n)

      list(
        estimate    = rd_aipw,
        se          = se_aipw,
        converged   = fit1$converged && fit0$converged,
        n0          = n0,
        n1          = n1,
        measure     = "RD",
        method_used = "aipw"
      )
    },
    error = function(e) {
      # Fall back to unadjusted RD on failure
      .gcomp_fallback(data_slice, treat.name, outcome.name,
                      "RD", n0, n1, method = "aipw_failed")
    })

    result
  }
}


# ---------------------------------------------------------------------------
# Pre-computation: fit PS model and augment the data frame
# Called once in forestsearch() before the consistency loop
# ---------------------------------------------------------------------------

#' Pre-compute propensity scores and inverse-PS weights
#'
#' Fits a logistic regression PS model on the full analysis dataset and
#' appends `.pihat` and `inv_ps` columns.  For randomized trials where
#' the true propensity is known, `ps_known` bypasses the model fit.
#'
#' @param df A `data.frame`.
#' @param treat.name Character. Treatment column name.
#' @param confounders.name Character vector. Covariate names for the PS model.
#' @param ps_known Numeric scalar or vector, or `NULL`.  If not `NULL`,
#'   the known propensity score (e.g., 0.5 for 1:1 randomisation) is used
#'   directly; no PS model is fitted.
#' @param ps_truncate Numeric. Lower bound for pihat truncation (and
#'   `1 - ps_truncate` for the upper bound).  Default 0.01.
#' @param verbose Logical. Print PS model summary. Default `FALSE`.
#'
#' @return The input `data.frame` with two additional columns:
#'   `.pihat` (estimated or known propensity score) and `inv_ps`
#'   (inverse PS weight: `1/pihat` for treated, `1/(1-pihat)` for control).
#'
#' @noRd
precompute_ps <- function(
    df,
    treat.name,
    confounders.name,
    ps_known     = NULL,
    ps_truncate  = 0.01,
    verbose      = FALSE
) {
  treat <- df[[treat.name]]

  if (!is.null(ps_known)) {
    # Known propensity score (RCT with known allocation ratio)
    pihat <- rep(ps_known, nrow(df))
    if (verbose) {
      message("Using known propensity score: ", ps_known)
    }
  } else {
    # Estimate PS via logistic regression
    ps_fmla <- stats::as.formula(
      paste0(treat.name, " ~ ", paste(confounders.name, collapse = " + "))
    )
    ps_fit <- stats::glm(ps_fmla, data = df, family = stats::binomial())
    pihat  <- stats::fitted(ps_fit)

    if (verbose) {
      message("PS model summary:")
      print(summary(ps_fit))
      message("pihat range: [", round(min(pihat), 4), ", ",
              round(max(pihat), 4), "]")
    }
  }

  # Truncate extreme probabilities
  pihat <- pmax(pmin(pihat, 1 - ps_truncate), ps_truncate)

  # Compute inverse PS weights (Bang & Robins convention)
  inv_ps <- ifelse(treat == 1L, 1 / pihat, 1 / (1 - pihat))

  df[[".pihat"]] <- pihat
  df[["inv_ps"]] <- inv_ps

  df
}


# ---------------------------------------------------------------------------
# Shared utilities for g-computation
# ---------------------------------------------------------------------------

#' Compute g-computation result from predicted potential outcomes
#' @noRd
.gcomp_result <- function(mu1, mu0, effect_measure, n0, n1,
                          converged, method) {
  switch(effect_measure,
    RD = {
      rd  <- mean(mu1, na.rm = TRUE) - mean(mu0, na.rm = TRUE)
      # Delta-method SE from individual-level predictions
      psi <- mu1 - mu0
      se  <- sqrt(stats::var(psi, na.rm = TRUE) / length(psi))
      list(
        estimate    = rd,
        se          = se,
        converged   = converged,
        n0          = n0,
        n1          = n1,
        measure     = "RD",
        method_used = method
      )
    },
    OR = {
      # G-computed log-OR: log( mean(mu1)/(1-mean(mu1)) ) - same for mu0
      p1 <- mean(mu1, na.rm = TRUE)
      p0 <- mean(mu0, na.rm = TRUE)
      log_or <- log(p1 / (1 - p1)) - log(p0 / (1 - p0))
      # Approximate SE via delta method on the probabilities
      n_eff  <- length(mu1)
      se_p1  <- sqrt(p1 * (1 - p1) / n_eff)
      se_p0  <- sqrt(p0 * (1 - p0) / n_eff)
      se_lor <- sqrt(1 / (p1 * (1 - p1) * n_eff) +
                     1 / (p0 * (1 - p0) * n_eff))
      list(
        estimate    = log_or,
        se          = se_lor,
        converged   = converged,
        n0          = n0,
        n1          = n1,
        measure     = "log-OR",
        method_used = method
      )
    },
    RR = {
      p1 <- mean(mu1, na.rm = TRUE)
      p0 <- mean(mu0, na.rm = TRUE)
      log_rr <- log(p1) - log(p0)
      n_eff  <- length(mu1)
      se_lrr <- sqrt((1 - p1) / (p1 * n_eff) + (1 - p0) / (p0 * n_eff))
      list(
        estimate    = log_rr,
        se          = se_lrr,
        converged   = converged,
        n0          = n0,
        n1          = n1,
        measure     = "log-RR",
        method_used = method
      )
    },
    stop("Unsupported effect_measure for g-computation: ", effect_measure,
         call. = FALSE)
  )
}


#' Fallback: raw group-level estimate when the model fails
#' @noRd
.gcomp_fallback <- function(data_slice, treat.name, outcome.name,
                            effect_measure, n0, n1, method) {
  y0 <- data_slice[[outcome.name]][data_slice[[treat.name]] == 0]
  y1 <- data_slice[[outcome.name]][data_slice[[treat.name]] == 1]
  p0 <- mean(y0, na.rm = TRUE)
  p1 <- mean(y1, na.rm = TRUE)

  switch(effect_measure,
    RD = {
      rd <- p1 - p0
      se <- sqrt(p1 * (1 - p1) / max(n1, 1L) + p0 * (1 - p0) / max(n0, 1L))
      list(estimate = rd, se = se, converged = FALSE,
           n0 = n0, n1 = n1, measure = "RD", method_used = method)
    },
    OR = {
      log_or <- log(p1 / (1 - p1)) - log(p0 / (1 - p0))
      se <- sqrt(1 / (p1 * (1 - p1) * max(n1, 1L)) +
                 1 / (p0 * (1 - p0) * max(n0, 1L)))
      list(estimate = log_or, se = se, converged = FALSE,
           n0 = n0, n1 = n1, measure = "log-OR", method_used = method)
    },
    RR = {
      log_rr <- log(p1) - log(p0)
      se <- sqrt((1 - p1) / (p1 * max(n1, 1L)) +
                 (1 - p0) / (p0 * max(n0, 1L)))
      list(estimate = log_rr, se = se, converged = FALSE,
           n0 = n0, n1 = n1, measure = "log-RR", method_used = method)
    },
    {
      # Default: raw RD
      rd <- p1 - p0
      se <- sqrt(p1 * (1 - p1) / max(n1, 1L) + p0 * (1 - p0) / max(n0, 1L))
      list(estimate = rd, se = se, converged = FALSE,
           n0 = n0, n1 = n1, measure = "RD", method_used = method)
    }
  )
}
