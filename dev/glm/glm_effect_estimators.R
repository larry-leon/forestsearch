# =============================================================================
# glm_effect_estimators.R
#
# Internal effect-estimator closures for GLM outcomes (binary, continuous).
# These are the counterparts to the inline Cox calls in subgroup_consistency.R
# and subgroup_search.R.  Each factory returns a function of the form:
#
#   function(data_slice) -> list(estimate, se, converged, n0, n1)
#
# where `estimate` is on the log scale (log-OR, log-RR) or identity scale (RD,
# MD) consistent with the chosen `effect_measure`.
#
# All functions are internal (@noRd) — not exported.
# =============================================================================


# -----------------------------------------------------------------------------
# Top-level dispatcher
# -----------------------------------------------------------------------------

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
#' @param event.name  Character or `NULL`. Name of the event indicator (survival
#'   only).
#' @param effect_measure Character or `NULL`.  For binary: `"OR"`, `"RR"`,
#'   `"RD"`.  For continuous: `"MD"`.  `NULL` uses the default for the outcome
#'   type (`"HR"` for survival, `"OR"` for binary, `"MD"` for continuous).
#' @param ... Additional arguments passed to the specific estimator factory.
#'
#' @return A closure `function(data_slice)`.
#'
#' @noRd
make_effect_estimator <- function(
    outcome_type,
    treat.name,
    outcome.name,
    event.name    = NULL,
    effect_measure = NULL,
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
      binary     = "OR",
      continuous = "MD"
    )
  }

  switch(outcome_type,
    survival = {
      .make_cox_estimator(
        treat.name   = treat.name,
        outcome.name = outcome.name,
        event.name   = event.name,
        ...
      )
    },
    binary = {
      effect_measure <- match.arg(effect_measure, choices = c("OR", "RR", "RD"))
      .make_glm_binary_estimator(
        treat.name     = treat.name,
        outcome.name   = outcome.name,
        effect_measure = effect_measure,
        ...
      )
    },
    continuous = {
      .make_lm_estimator(
        treat.name   = treat.name,
        outcome.name = outcome.name,
        ...
      )
    }
  )
}


# -----------------------------------------------------------------------------
# Survival estimator (Cox) — wraps existing inline logic into the closure form
# -----------------------------------------------------------------------------

#' @noRd
.make_cox_estimator <- function(treat.name, outcome.name, event.name, ...) {

  force(treat.name)
  force(outcome.name)
  force(event.name)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    result <- tryCatch({
      surv_obj <- survival::Surv(
        time  = data_slice[[outcome.name]],
        event = data_slice[[event.name]]
      )
      fmla <- stats::as.formula(
        paste0("surv_obj ~ ", treat.name)
      )
      fit <- survival::coxph(fmla, data = data_slice)
      coef_val <- stats::coef(fit)[[treat.name]]
      se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
      list(estimate = coef_val, se = se_val, converged = TRUE, n0 = n0, n1 = n1)
    },
    error = function(e) {
      list(estimate = NA_real_, se = NA_real_, converged = FALSE, n0 = n0, n1 = n1)
    })
    result
  }
}


# -----------------------------------------------------------------------------
# Binary GLM estimator
# -----------------------------------------------------------------------------

#' Make a binary-outcome GLM estimator closure
#'
#' Supports three effect measures:
#'   - `"OR"` : logistic regression; returns log-OR
#'   - `"RR"` : log-binomial (with modified-Poisson fallback); returns log-RR
#'   - `"RD"` : identity-link binomial; returns risk difference
#'
#' All estimates are returned on the natural (non-exponentiated) scale to match
#' the consistency-criterion comparisons (i.e., `log(OR) >= log(or.threshold)`).
#' For `"RD"` no transformation is applied.
#'
#' @noRd
.make_glm_binary_estimator <- function(
    treat.name,
    outcome.name,
    effect_measure = "OR",
    fallback_or    = TRUE,   # fall back to logistic if log-binomial fails
    ...
) {
  force(treat.name)
  force(outcome.name)
  force(effect_measure)
  force(fallback_or)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    fmla <- stats::as.formula(
      paste0(outcome.name, " ~ ", treat.name)
    )

    result <- switch(effect_measure,

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
            estimate  = coef_val,   # log-OR
            se        = se_val,
            converged = fit$converged,
            n0 = n0, n1 = n1,
            measure   = "log-OR"
          )
        },
        error = function(e) {
          list(
            estimate = NA_real_, se = NA_real_, converged = FALSE,
            n0 = n0, n1 = n1, measure = "log-OR"
          )
        })
      },

      RR = {
        # Attempt log-binomial; fall back to modified Poisson on non-convergence
        fit_rr <- tryCatch({
          fit <- stats::glm(
            fmla,
            data   = data_slice,
            family = stats::binomial(link = "log"),
            # start values often help log-binomial converge
            start  = c(
              log(mean(data_slice[[outcome.name]], na.rm = TRUE)),
              0
            )
          )
          if (!fit$converged) stop("log-binomial did not converge")
          coef_val <- stats::coef(fit)[[treat.name]]
          se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
          list(
            estimate  = coef_val,
            se        = se_val,
            converged = TRUE,
            link_used = "log-binomial",
            n0 = n0, n1 = n1, measure = "log-RR"
          )
        },
        error = function(e) NULL)

        if (!is.null(fit_rr)) {
          fit_rr
        } else if (fallback_or) {
          # Modified Poisson (Zou 2004): Poisson with robust sandwich SE
          tryCatch({
            fit <- stats::glm(
              fmla,
              data   = data_slice,
              family = stats::poisson(link = "log")
            )
            coef_val <- stats::coef(fit)[[treat.name]]
            # Use sandwich robust SE if available; else model-based
            se_val <- tryCatch({
              sqrt(
                sandwich::sandwich(fit)[treat.name, treat.name]
              )
            }, error = function(e2) {
              sqrt(diag(stats::vcov(fit)))[[treat.name]]
            })
            list(
              estimate  = coef_val,
              se        = se_val,
              converged = TRUE,
              link_used = "modified-poisson",
              n0 = n0, n1 = n1, measure = "log-RR"
            )
          },
          error = function(e) {
            list(
              estimate = NA_real_, se = NA_real_, converged = FALSE,
              link_used = "failed", n0 = n0, n1 = n1, measure = "log-RR"
            )
          })
        } else {
          list(
            estimate = NA_real_, se = NA_real_, converged = FALSE,
            link_used = "failed", n0 = n0, n1 = n1, measure = "log-RR"
          )
        }
      },

      RD = {
        tryCatch({
          # Identity-link binomial; start at group means
          p0    <- mean(data_slice[[outcome.name]][
            data_slice[[treat.name]] == 0], na.rm = TRUE)
          p1    <- mean(data_slice[[outcome.name]][
            data_slice[[treat.name]] == 1], na.rm = TRUE)
          start <- c(p0, p1 - p0)

          fit <- stats::glm(
            fmla,
            data   = data_slice,
            family = stats::binomial(link = "identity"),
            start  = start
          )
          coef_val <- stats::coef(fit)[[treat.name]]
          se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
          list(
            estimate  = coef_val,   # risk difference (no log transform)
            se        = se_val,
            converged = fit$converged,
            n0 = n0, n1 = n1, measure = "RD"
          )
        },
        error = function(e) {
          # Last-resort: raw risk difference from group means
          p0 <- mean(data_slice[[outcome.name]][
            data_slice[[treat.name]] == 0], na.rm = TRUE)
          p1 <- mean(data_slice[[outcome.name]][
            data_slice[[treat.name]] == 1], na.rm = TRUE)
          rd   <- p1 - p0
          # Delta-method SE: sqrt(p1(1-p1)/n1 + p0(1-p0)/n0)
          se   <- sqrt(p1 * (1 - p1) / max(n1, 1) + p0 * (1 - p0) / max(n0, 1))
          list(
            estimate  = rd,
            se        = se,
            converged = FALSE,   # flag as fallback
            n0 = n0, n1 = n1, measure = "RD"
          )
        })
      }
    )

    result
  }
}


# -----------------------------------------------------------------------------
# Continuous outcome estimator (linear regression / mean difference)
# -----------------------------------------------------------------------------

#' Make a continuous-outcome estimator closure
#'
#' Returns the within-subgroup OLS treatment coefficient (mean difference,
#' experimental minus control) and its standard error.
#'
#' @noRd
.make_lm_estimator <- function(treat.name, outcome.name, ...) {

  force(treat.name)
  force(outcome.name)

  function(data_slice) {
    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    fmla <- stats::as.formula(
      paste0(outcome.name, " ~ ", treat.name)
    )

    result <- tryCatch({
      fit      <- stats::lm(fmla, data = data_slice)
      coef_val <- stats::coef(fit)[[treat.name]]
      se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
      list(
        estimate  = coef_val,   # mean difference (identity scale)
        se        = se_val,
        converged = TRUE,
        n0 = n0, n1 = n1, measure = "MD"
      )
    },
    error = function(e) {
      list(
        estimate = NA_real_, se = NA_real_, converged = FALSE,
        n0 = n0, n1 = n1, measure = "MD"
      )
    })

    result
  }
}


# -----------------------------------------------------------------------------
# Utility: back-transform for reporting
# -----------------------------------------------------------------------------

#' Exponentiate an estimate for a given effect measure
#'
#' Applies `exp()` for log-scale measures (HR, OR, RR); returns the estimate
#' unchanged for identity-scale measures (RD, MD).
#'
#' @param estimate Numeric scalar.
#' @param measure  Character. One of `"HR"`, `"OR"`, `"RR"`, `"RD"`, `"MD"`.
#'
#' @return Numeric scalar on the natural scale.
#'
#' @noRd
back_transform_estimate <- function(estimate, measure) {
  if (measure %in% c("HR", "OR", "RR")) {
    exp(estimate)
  } else {
    estimate
  }
}
