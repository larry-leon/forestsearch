# File: R/dina_subgroup_refit.R
# Part of the forestsearch package.
#
# Standard within-subgroup treatment-effect model for a subgroup
# discovered by dina_subgroup().
#
# dina_subgroup() reports the DINA / HTE quantity for the identified
# subgroup S*: the subgroup-mean of the per-patient effect estimate,
# a*^T beta-hat, a linear projection of the cross-fitted CATE surface
# (a "best linear projection" analog).  That is NOT the same as the
# treatment contrast a clinician usually wants: the *standard*
# regression model -- a Cox model for survival, or a GLM otherwise --
# fit on the patients in S*, comparing the treatment arms directly.
#
# This file provides that standard within-subgroup estimate.  The two
# quantities coincide only under correct linearity; in general they
# answer related but distinct questions, so both are worth reporting.
#
# IMPORTANT INFERENCE CAVEAT.  The Wald CI returned here treats the
# subgroup definition (the discovered signature) as *pre-specified*.
# It does NOT adjust for the data-driven selection of that signature.
# Removing selection ("winner's curse") bias and obtaining valid
# selection-adjusted coverage requires the more involved bias-corrected
# bootstrap of Leon et al. (2024, Stat. Med., DOI 10.1002/sim.10163,
# Section 3.2) -- a two-term bias correction with infinitesimal-jackknife
# variance -- which simple refitting does not deliver.  Use this estimate
# as the clinically interpretable within-signature contrast, not as a
# selection-adjusted inference.


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Standard within-subgroup treatment-effect model for a DINA subgroup
#'
#' Given a subgroup discovered by \code{\link{dina_subgroup}}, fit the
#' \emph{standard} treatment-effect model on the patients in that
#' subgroup: a Cox proportional-hazards model when
#' \code{sg$family == "cox"}, or a generalized linear model otherwise.
#' The reported effect is the treatment coefficient on the
#' natural-parameter scale (log-HR for Cox, log-OR for binomial,
#' log-IRR for Poisson, mean difference for Gaussian), directly
#' comparable to the DINA effect returned by \code{dina_subgroup()} but
#' obtained from a conventional within-subgroup regression rather than
#' from the cross-fitted CATE surface.
#'
#' The model is unadjusted by default (\code{~ treatment}); optional
#' covariate adjustment and stratification are available via
#' \code{confounders} and \code{strata}.
#'
#' @section Adjustment set (the confounders argument):
#' Three modes:
#' \itemize{
#'   \item \code{"none"} (default) -- unadjusted: \code{~ treatment}.
#'   \item \code{NULL} -- automatic: the DINA covariates supplied in
#'     \code{covariates}, with the subgroup-defining covariate
#'     (\code{sg$covariate}) omitted.  Omitting it avoids adjusting on
#'     a variable whose range is restricted by the subgroup definition.
#'   \item a character vector -- exactly those columns.  A warning is
#'     issued if it includes \code{sg$covariate}.
#' }
#'
#' @section Inference caveat:
#' The Wald confidence interval treats the subgroup definition as
#' \emph{pre-specified}.  It does \strong{not} adjust for the
#' data-driven selection of the signature; it is a within-signature
#' (conditional) interval, not a selection-adjusted one.  A
#' selection-adjusted interval requires a bias-corrected bootstrap of
#' the full search procedure (Leon et al. 2024, Section 3.2) and is not
#' provided here.
#'
#' @param sg a \code{"dina_subgroup"} object from
#'   \code{\link{dina_subgroup}} with \code{sg$found == TRUE}.
#' @param df data frame the subgroup was identified on; must contain the
#'   outcome, treatment, covariate, and any \code{strata} columns.
#' @param treatment character(1); name of the binary treatment column.
#' @param outcome character(1); name of the response column.  For
#'   \code{family = "cox"} this is the event/censoring time column.
#' @param covariates character vector of the DINA covariate names.  Used
#'   only to construct the automatic adjustment set when
#'   \code{confounders = NULL}; ignored otherwise.
#' @param status character(1) or \code{NULL}; for Cox only, the event
#'   indicator column name.  Required when \code{sg$family == "cox"}.
#' @param strata \code{NULL} (default) or a character vector of column
#'   names entered as Cox \code{strata()} terms.  Only valid for
#'   \code{family = "cox"}.
#' @param confounders \code{"none"} (default, unadjusted), \code{NULL}
#'   (automatic; see Details), or a character vector of column names.
#' @param alpha confidence level for the Wald interval.  Default
#'   \code{0.05} for a 95% CI.
#'
#' @return An object of class \code{"dina_subgroup_refit"}, a list with
#'   components:
#'   \describe{
#'     \item{effect}{treatment coefficient on the natural-parameter
#'       scale (the point estimate).}
#'     \item{se}{Wald standard error of \code{effect}.}
#'     \item{ci}{named length-2 numeric vector (\code{lower},
#'       \code{upper}); Wald \code{1 - alpha} CI on the
#'       natural-parameter scale.}
#'     \item{effect_scale}{character label for the scale: one of
#'       \code{"log-HR"}, \code{"log-OR"}, \code{"log-IRR"},
#'       \code{"MD"}.}
#'     \item{ratio_scale}{logical; \code{TRUE} when \code{effect} is a
#'       log-ratio (Cox/binomial/Poisson) so that \code{exp()} maps it
#'       to HR/OR/IRR.}
#'     \item{signature}{character; the subgroup rule, e.g.
#'       \code{"nodes >= 12"}.}
#'     \item{formula}{the fitted model formula.}
#'     \item{confounders_used, strata_used}{the resolved adjustment set
#'       and stratification columns actually used.}
#'     \item{n_subgroup, n_treated, n_control}{subgroup size and per-arm
#'       counts.}
#'     \item{n_events, n_events_treated, n_events_control}{event counts
#'       within the subgroup (Cox only; \code{NA} otherwise).}
#'     \item{family, alpha}{echoed inputs.}
#'     \item{model}{the fitted \code{coxph} / \code{glm} object.}
#'     \item{call}{the matched call.}
#'   }
#'
#' @seealso \code{\link{dina_subgroup}} for the subgroup search and its
#'   DINA (BLP-analog) effect; \code{\link{dina_subgroup_bootstrap}} for
#'   bootstrap inference that reports this within-subgroup model
#'   alongside the DINA effect.
#'
#' @examples
#' \dontrun{
#' fit <- dina(df, outcome = "time", treatment = "trt",
#'             covariates = c("age", "nodes", "er"),
#'             family = "cox", status = "status", seed = 1L)
#' sg  <- dina_subgroup(fit, df, covariates = c("age", "nodes", "er"),
#'                      m_diff = 0, n_min = 60L)
#'
#' # Unadjusted within-subgroup Cox model (default):
#' dina_subgroup_refit(sg, df, treatment = "trt", outcome = "time",
#'                     covariates = c("age", "nodes", "er"),
#'                     status = "status")
#'
#' # Adjusted for DINA covariates except the subgroup-defining one:
#' dina_subgroup_refit(sg, df, treatment = "trt", outcome = "time",
#'                     covariates = c("age", "nodes", "er"),
#'                     status = "status", confounders = NULL)
#' }
#'
#' @importFrom stats coef vcov qnorm as.formula
#' @export
dina_subgroup_refit <- function(sg,
                                df,
                                treatment,
                                outcome,
                                covariates,
                                status = NULL,
                                strata = NULL,
                                confounders = "none",
                                alpha = 0.05) {

  call <- match.call()

  if (!inherits(sg, "dina_subgroup")) {
    stop("`sg` must be a \"dina_subgroup\" object from dina_subgroup().")
  }
  if (!isTRUE(sg$found)) {
    stop("`sg` did not identify a subgroup (sg$found is not TRUE); ",
         "there is nothing to refit.")
  }
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  }
  if (!is.null(sg$n_total) && !identical(as.integer(sg$n_total),
                                         as.integer(nrow(df)))) {
    stop("`sg` was computed on ", sg$n_total, " rows but `df` has ",
         nrow(df), " rows; `sg` must be computed on the same `df`.")
  }
  if (length(sg$mask) != nrow(df)) {
    stop("`sg$mask` length (", length(sg$mask), ") does not match ",
         "nrow(df) (", nrow(df), "); `sg` must be computed on the same `df`.")
  }
  if (length(alpha) != 1L || !is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single numeric in (0, 1).")
  }

  family <- sg$family

  # ---- Resolve the adjustment set ---------------------------------------
  confounders_used <- .dina_resolve_confounders(
    confounders = confounders,
    covariates  = covariates,
    sg_covariate = sg$covariate
  )

  # ---- Subset to the discovered subgroup --------------------------------
  df_sub <- df[sg$mask, , drop = FALSE]

  # ---- Fit the standard within-subgroup model ---------------------------
  fitres <- .dina_fit_subgroup_model(
    df_sub      = df_sub,
    treatment   = treatment,
    outcome     = outcome,
    status      = status,
    family      = family,
    confounders = confounders_used,
    strata      = strata
  )

  # ---- Wald CI on the natural-parameter scale ---------------------------
  z   <- stats::qnorm(1 - alpha / 2)
  ci  <- c(lower = fitres$effect - z * fitres$se,
           upper = fitres$effect + z * fitres$se)

  scale_info  <- .dina_effect_scale(family)
  signature   <- .dina_signature_string(sg)

  out <- list(
    effect            = fitres$effect,
    se                = fitres$se,
    ci                = ci,
    effect_scale      = scale_info$label,
    ratio_scale       = scale_info$ratio,
    signature         = signature,
    formula           = fitres$formula,
    confounders_used  = confounders_used,
    strata_used       = if (is.null(strata)) character(0) else strata,
    n_subgroup        = fitres$n_subgroup,
    n_treated         = fitres$n_treated,
    n_control         = fitres$n_control,
    n_events          = fitres$n_events,
    n_events_treated  = fitres$n_events_treated,
    n_events_control  = fitres$n_events_control,
    family            = family,
    alpha             = alpha,
    model             = fitres$model,
    call              = call
  )
  class(out) <- "dina_subgroup_refit"
  out
}


# ---------------------------------------------------------------------------
# Internal: resolve the confounder adjustment set
# ---------------------------------------------------------------------------

#' Resolve the `confounders` argument into an explicit character vector.
#'
#' `"none"` -> character(0) (unadjusted); `NULL` -> DINA covariates with
#' the subgroup-defining covariate omitted; a character vector -> itself
#' (with a warning if it contains the subgroup covariate).
#'
#' @noRd
.dina_resolve_confounders <- function(confounders, covariates, sg_covariate) {

  # Sentinel "none" -> unadjusted.
  if (is.character(confounders) && length(confounders) == 1L &&
      identical(confounders, "none")) {
    return(character(0))
  }

  # NULL -> automatic: DINA covariates minus the subgroup-defining factor.
  if (is.null(confounders)) {
    if (!is.character(covariates) || length(covariates) < 1L) {
      stop("`confounders = NULL` requests the automatic adjustment set, ",
           "but `covariates` is empty; supply the DINA covariate names.")
    }
    auto <- setdiff(covariates, sg_covariate)
    return(auto)
  }

  # Explicit character vector.
  if (!is.character(confounders) || anyNA(confounders)) {
    stop("`confounders` must be \"none\", NULL, or a character vector ",
         "of column names.")
  }
  if (sg_covariate %in% confounders) {
    warning("`confounders` includes the subgroup-defining covariate \"",
            sg_covariate, "\"; its range is restricted within the subgroup, ",
            "making it a poor adjustment term.  Consider removing it.",
            call. = FALSE)
  }
  confounders
}


# ---------------------------------------------------------------------------
# Internal: fit the standard model on a given subset
# ---------------------------------------------------------------------------

#' Fit the standard within-subgroup model and extract the treatment effect.
#'
#' Shared by `dina_subgroup_refit()` (original-data point estimate) and
#' the per-iteration bootstrap.  Builds the model formula, fits a Cox
#' model (`family == "cox"`) or a GLM (otherwise) on `df_sub`, and
#' returns the treatment coefficient, its Wald SE, arm/event counts, the
#' fitted model, and the formula.
#'
#' Loud failures: errors if a treatment arm is empty in the subset, if
#' the treatment coefficient is absent from the fit, or if the
#' underlying model fit errors.
#'
#' @noRd
.dina_fit_subgroup_model <- function(df_sub, treatment, outcome, status,
                                     family, confounders, strata) {

  if (!all(c(treatment, outcome) %in% names(df_sub))) {
    stop("`df` must contain the treatment (\"", treatment,
         "\") and outcome (\"", outcome, "\") columns.")
  }
  miss_conf <- setdiff(confounders, names(df_sub))
  if (length(miss_conf) > 0L) {
    stop("Confounder column(s) not found in `df`: ",
         paste(miss_conf, collapse = ", "), ".")
  }
  if (!is.null(strata)) {
    if (family != "cox") {
      stop("`strata` is only supported for family = \"cox\".")
    }
    miss_strata <- setdiff(strata, names(df_sub))
    if (length(miss_strata) > 0L) {
      stop("Strata column(s) not found in `df`: ",
           paste(miss_strata, collapse = ", "), ".")
    }
  }

  trt <- df_sub[[treatment]]
  n_treated <- sum(trt == 1, na.rm = TRUE)
  n_control <- sum(trt == 0, na.rm = TRUE)
  if (n_treated == 0L || n_control == 0L) {
    stop("Within the subgroup, the treatment arm sizes are ",
         "treated = ", n_treated, ", control = ", n_control,
         "; both arms must be non-empty to estimate a contrast.")
  }

  # ---- Build the model formula ------------------------------------------
  rhs_terms <- c(treatment, confounders)
  if (family == "cox" && !is.null(strata) && length(strata) > 0L) {
    rhs_terms <- c(rhs_terms, sprintf("strata(%s)", strata))
  }
  rhs <- paste(rhs_terms, collapse = " + ")

  if (family == "cox") {
    if (is.null(status) || !nzchar(status)) {
      stop("`status` (event indicator column) is required for ",
           "family = \"cox\".")
    }
    if (!status %in% names(df_sub)) {
      stop("Status column \"", status, "\" not found in `df`.")
    }
    lhs  <- sprintf("Surv(%s, %s)", outcome, status)
    form <- stats::as.formula(paste(lhs, "~", rhs))

    model <- tryCatch(
      survival::coxph(form, data = df_sub),
      error = function(e) {
        stop("Within-subgroup coxph() failed: ", conditionMessage(e),
             call. = FALSE)
      }
    )

    ev <- df_sub[[status]]
    n_events         <- sum(ev == 1, na.rm = TRUE)
    n_events_treated <- sum(ev == 1 & trt == 1, na.rm = TRUE)
    n_events_control <- sum(ev == 1 & trt == 0, na.rm = TRUE)

  } else {
    lhs  <- outcome
    form <- stats::as.formula(paste(lhs, "~", rhs))

    model <- tryCatch(
      stats::glm(form, data = df_sub, family = family),
      error = function(e) {
        stop("Within-subgroup glm() failed: ", conditionMessage(e),
             call. = FALSE)
      }
    )

    n_events         <- NA_integer_
    n_events_treated <- NA_integer_
    n_events_control <- NA_integer_
  }

  # ---- Extract the treatment coefficient and Wald SE --------------------
  cf <- stats::coef(model)
  if (!treatment %in% names(cf)) {
    stop("Treatment coefficient \"", treatment, "\" not found in the ",
         "fitted model; is the treatment column binary 0/1 numeric?")
  }
  effect <- unname(cf[[treatment]])

  V <- stats::vcov(model)
  if (!treatment %in% rownames(V)) {
    stop("Treatment term \"", treatment, "\" not found in the model ",
         "variance-covariance matrix.")
  }
  se <- sqrt(V[treatment, treatment])

  list(
    effect           = effect,
    se               = se,
    n_subgroup       = nrow(df_sub),
    n_treated        = n_treated,
    n_control        = n_control,
    n_events         = n_events,
    n_events_treated = n_events_treated,
    n_events_control = n_events_control,
    formula          = form,
    model            = model
  )
}


# ---------------------------------------------------------------------------
# Internal: signature helpers
# ---------------------------------------------------------------------------

#' Logical mask for a subgroup rule applied to an arbitrary data frame.
#'
#' Re-applies the discovered (covariate, direction, threshold) rule to
#' `df` -- e.g. to a bootstrap resample -- returning a logical vector of
#' length `nrow(df)`.  `direction == "left"` is `x <= q`, otherwise
#' `x >= q`, matching `dina_subgroup()`.
#'
#' @noRd
.dina_signature_mask <- function(df, covariate, direction, threshold) {
  if (!covariate %in% names(df)) {
    stop("Subgroup covariate \"", covariate, "\" not found in `df`.")
  }
  x <- df[[covariate]]
  if (direction == "left") x <= threshold else x >= threshold
}

#' Human-readable signature string, e.g. "nodes >= 12".
#'
#' @noRd
.dina_signature_string <- function(sg) {
  cmp <- if (sg$direction == "left") "<=" else ">="
  sprintf("%s %s %s", sg$covariate, cmp, format(sg$threshold))
}

#' Effect-scale label and ratio flag for a family.
#'
#' @noRd
.dina_effect_scale <- function(family) {
  switch(
    family,
    cox      = list(label = "log-HR",  ratio = TRUE),
    binomial = list(label = "log-OR",  ratio = TRUE),
    poisson  = list(label = "log-IRR", ratio = TRUE),
    gaussian = list(label = "MD",      ratio = FALSE),
    stop("Unsupported family: ", family)
  )
}


# ---------------------------------------------------------------------------
# S3 method
# ---------------------------------------------------------------------------

#' Print a `dina_subgroup_refit` result.
#'
#' @param x a \code{"dina_subgroup_refit"} object.
#' @param digits number of digits for the numeric summary.
#' @param ... unused.
#' @return invisibly returns \code{x}.
#' @export
print.dina_subgroup_refit <- function(x,
                                      digits = max(3L, getOption("digits") - 3L),
                                      ...) {
  cat("Standard within-subgroup model (DINA-identified subgroup)\n")
  cat("  Family:       ", x$family, "\n", sep = "")
  cat("  Signature:    ", x$signature, "\n", sep = "")
  cat("  Subgroup n:   ", x$n_subgroup,
      " (treated ", x$n_treated, ", control ", x$n_control, ")\n", sep = "")
  if (!is.na(x$n_events)) {
    cat("  Events:       ", x$n_events,
        " (treated ", x$n_events_treated,
        ", control ", x$n_events_control, ")\n", sep = "")
  }

  adj <- if (length(x$confounders_used) == 0L) "unadjusted"
         else paste(x$confounders_used, collapse = ", ")
  cat("  Adjustment:   ", adj, "\n", sep = "")
  if (length(x$strata_used) > 0L) {
    cat("  Strata:       ", paste(x$strata_used, collapse = ", "), "\n",
        sep = "")
  }

  ci_pct <- format(100 * (1 - x$alpha), digits = 3L)
  cat("\n  Treatment effect (", x$effect_scale, "):  ",
      format(x$effect, digits = digits),
      "  (SE ", format(x$se, digits = digits), ")\n", sep = "")
  cat("  Wald ", ci_pct, "% CI (", x$effect_scale, "):  [",
      format(x$ci[["lower"]], digits = digits), ", ",
      format(x$ci[["upper"]], digits = digits), "]\n", sep = "")

  if (isTRUE(x$ratio_scale)) {
    ratio_lab <- sub("^log-", "", x$effect_scale)
    cat("  ", ratio_lab, ":  ", format(exp(x$effect), digits = digits),
        "   95% CI (", ratio_lab, "):  [",
        format(exp(x$ci[["lower"]]), digits = digits), ", ",
        format(exp(x$ci[["upper"]]), digits = digits), "]\n", sep = "")
  }

  cat("\n  Note: CI is conditional on the discovered signature treated as\n")
  cat("  pre-specified; it does NOT adjust for signature selection.\n")
  invisible(x)
}
