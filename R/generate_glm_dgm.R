# =============================================================================
# generate_glm_dgm() -- GLM Data Generating Mechanism
# =============================================================================
#
# Mirrors generate_aft_dgm_flex() for non-survival outcomes.
# Supports binary (logistic), continuous (Gaussian), and count (Poisson)
# endpoints with optional offset for rate-based measures.
#
# Key design decisions:
#   - Reuses define_subgroups() from generate_aft_dgm_helpers.R for subgroup
#     flag creation and flexible cutpoint specifications.
#   - Returns hazard_ratios list with harm_subgroup / no_harm_subgroup /
#     overall for compatibility with get_dgm_hr() and reporting functions.
#   - Stores individual-level potential outcomes (p0, p1) in df_super for
#     simulate_from_glm_dgm() to draw from.
# =============================================================================


#' Generate GLM-Based Data Generating Mechanism
#'
#' Creates a data generating mechanism (DGM) for non-survival outcomes using
#' generalised linear models.
#'
#' Fits a GLM to the source data, creates a super-population with
#' individual-level potential outcomes, and stores true subgroup effects.
#' The returned object can be passed to \code{\link{simulate_from_glm_dgm}}
#' for simulation studies, and is compatible with
#' \code{\link{run_simulation_analysis}} (which dispatches on class).
#'
#' @param data A data frame containing the source dataset.
#' @param factor_vars Character vector of factor/categorical variable names
#'   to include as covariates. These are the candidate confounders for
#'   ForestSearch (e.g., \code{paste0("z", 1:12)}).
#' @param outcome_var Character string naming the outcome variable.
#' @param treatment_var Character string naming the treatment variable
#'   (must be 0/1 coded).
#' @param outcome_type Character. One of \code{"binary"},
#'   \code{"continuous"}, or \code{"count"}.
#' @param effect_measure Character. Effect measure for the GLM.
#'   Default is \code{NULL}, which selects \code{"OR"} for binary,
#'   \code{"MD"} for continuous, \code{"IRR"} for count outcomes.
#' @param offset_var Character or \code{NULL}. Name of the exposure /
#'   follow-up time column for count outcomes with an offset.
#'   Required when \code{effect_measure} is \code{"IRR"} or
#'   \code{"IRD"}.  The offset enters the Poisson GLM as
#'   \code{offset(log(offset_var))}.  Default \code{NULL}.
#' @param subgroup_vars Character vector of variable names defining the
#'   subgroup. Default \code{NULL} (no subgroup).
#' @param subgroup_cuts Named list of cutpoint specifications for
#'   subgroup variables. Uses the same flexible format as
#'   \code{\link{generate_aft_dgm_flex}}: fixed values, quantiles, or
#'   functions.
#' @param model Character. \code{"alt"} for alternative hypothesis with
#'   treatment--subgroup interaction, \code{"null"} for uniform treatment
#'   effect. Default \code{"alt"}.
#' @param k_treat Numeric. Scaling factor for the fitted treatment
#'   coefficient on the linear predictor scale.  \code{k_treat = 1}
#'   (default) preserves the fitted treatment effect as-is.
#'   \code{k_treat > 1} amplifies the treatment main effect, useful
#'   when the source dataset has a weak treatment effect and a
#'   stronger ITT is desired for simulation.  \code{k_treat < 1}
#'   attenuates the effect.  Analogous to the \code{k_treat}
#'   parameter in \code{\link{generate_aft_dgm_flex}}.
#'   Note: the null DGM (\code{model = "null"}) represents a
#'   \emph{homogeneous} treatment effect (no HTE), not a zero
#'   treatment effect.  \code{k_treat} scales the magnitude of
#'   that homogeneous effect.
#' @param k_inter Numeric. Direct additive shift on the linear predictor
#'   for Q members under treatment.  When \code{adverse_outcome = FALSE}
#'   (the default), positive values increase \eqn{P(Y=1)} for Q under
#'   treatment.  When \code{adverse_outcome = TRUE} for binary outcomes,
#'   the shift is \strong{negated internally} so that positive
#'   \code{k_inter} consistently means "amplify the treatment--subgroup
#'   contrast on the \emph{beneficial} (1-Y) scale."
#'
#'   This ensures that the \emph{same} \code{k_inter} produces the
#'   \emph{same} biological heterogeneity regardless of whether Y is
#'   coded as an adverse or beneficial event.  Without this adjustment,
#'   flipping Y from beneficial to adverse changes the baseline
#'   treatment-coefficient sign, and the same \code{target_effect} in
#'   \code{\link{calibrate_glm_interaction}} can produce vastly
#'   different interaction strengths (see Details).
#'
#'   Default \code{0} (no interaction, equivalent to \code{model = "null"}).
#' @param adverse_outcome Logical.  If \code{TRUE}, the outcome
#'   variable represents an \emph{adverse} event (e.g., failure,
#'   non-improvement) where Y = 1 is clinically \emph{bad}.  For
#'   binary outcomes this causes \code{k_inter} to be negated
#'   internally so that positive \code{k_inter} amplifies the
#'   treatment contrast on the beneficial (1-Y) scale.  Ignored
#'   for continuous and count outcomes.  Default \code{FALSE}
#'   (backward compatible: Y = 1 is beneficial or direction is
#'   irrelevant).
#' @param n_super Integer. Size of the super-population. Default
#'   \code{5000L}.
#' @param seed Integer. Random seed for super-population sampling.
#'   Default \code{8316951L}.
#' @param verbose Logical. Print diagnostic information. Default
#'   \code{FALSE}.
#'
#' @return An object of class \code{c("glm_dgm", "list")} with:
#'   \describe{
#'     \item{\code{df_super}}{Super-population data frame with potential
#'       outcome columns: \code{p0}, \code{p1} (binary probabilities),
#'       or \code{mu0}, \code{mu1} (conditional means for continuous /
#'       count outcomes), and \code{flag_harm}.}
#'     \item{\code{hazard_ratios}}{Named list with \code{harm_subgroup},
#'       \code{no_harm_subgroup}, and \code{overall} -- effect estimates
#'       on the scale determined by \code{effect_measure}.  Field name
#'       retained for compatibility with \code{get_dgm_hr()} and
#'       reporting functions.}
#'     \item{\code{outcome_type}}{Character: \code{"binary"},
#'       \code{"continuous"}, or \code{"count"}.}
#'     \item{\code{effect_measure}}{Character: the effect measure used.}
#'     \item{\code{model_params}}{List with fitted coefficients, family,
#'       \code{k_inter}, residual SD (continuous), and offset variable
#'       name (count).}
#'     \item{\code{subgroup_info}}{List with subgroup definition,
#'       size, and proportion.}
#'     \item{\code{model_type}}{Character: \code{"alt"} or \code{"null"}.}
#'   }
#'
#' @seealso \code{\link{simulate_from_glm_dgm}},
#'   \code{\link{calibrate_glm_interaction}},
#'   \code{\link{generate_aft_dgm_flex}}
#'
#' @examples
#' \dontrun{
#' library(speff2trial)
#' actg <- subset(ACTG175, arms %in% c(1, 3))
#' actg$treat  <- 1L - as.integer(actg$arms == 1)
#' actg$y      <- as.integer(actg$cd420 > actg$cd40)
#' actg$z1     <- as.factor(ifelse(actg$age > 34, 1L, 0L))
#' actg$z2     <- as.factor(ifelse(actg$preanti <= 745, 1L, 0L))
#'
#' dgm <- generate_glm_dgm(
#'   data          = actg,
#'   factor_vars   = c("z1", "z2"),
#'   outcome_var   = "y",
#'   treatment_var = "treat",
#'   outcome_type  = "binary",
#'   subgroup_vars = c("z1", "z2"),
#'   subgroup_cuts = list(z1 = 1L, z2 = 1L),
#'   model         = "alt",
#'   k_inter       = 1.0,   # additive log-odds shift for Q under treatment
#'   verbose       = TRUE
#' )
#' }
#'
#' @export
#' @importFrom stats glm binomial coef predict qlogis plogis as.formula
generate_glm_dgm <- function(
    data,
    factor_vars,
    outcome_var,
    treatment_var,
    outcome_type    = c("binary", "continuous", "count"),
    effect_measure  = NULL,
    offset_var      = NULL,
    subgroup_vars   = NULL,
    subgroup_cuts   = NULL,
    model           = c("alt", "null"),
    k_treat         = 1,
    k_inter         = 0,
    adverse_outcome = FALSE,
    n_super         = 5000L,
    seed            = 8316951L,
    verbose         = FALSE
) {

  # -- Argument matching -----------------------------------------------------
  outcome_type <- match.arg(outcome_type)
  model        <- match.arg(model)

  if (is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      binary     = "OR",
      continuous = "MD",
      count      = "IRR"
    )
  }

  # -- Validate inputs -------------------------------------------------------
  stopifnot(
    is.data.frame(data),
    outcome_var %in% names(data),
    treatment_var %in% names(data),
    all(factor_vars %in% names(data)),
    is.numeric(k_inter), length(k_inter) == 1
  )

  if (outcome_type == "binary") {
    y_vals <- unique(data[[outcome_var]])
    if (!all(y_vals %in% c(0L, 1L, 0, 1))) {
      stop("For binary outcomes, '", outcome_var,
           "' must contain only 0/1 values.")
    }
  }

  if (outcome_type == "count") {
    y_vals <- data[[outcome_var]]
    if (any(y_vals < 0, na.rm = TRUE) || any(y_vals != floor(y_vals), na.rm = TRUE)) {
      stop("For count outcomes, '", outcome_var,
           "' must contain non-negative integers.")
    }
    if (!is.null(offset_var)) {
      if (!offset_var %in% names(data)) {
        stop("offset_var '", offset_var, "' not found in data.")
      }
      if (any(data[[offset_var]] <= 0, na.rm = TRUE)) {
        stop("offset_var '", offset_var,
             "' must contain strictly positive values (exposure/person-time).")
      }
    }
  }

  if (outcome_type == "continuous") {
    if (!is.numeric(data[[outcome_var]])) {
      stop("For continuous outcomes, '", outcome_var, "' must be numeric.")
    }
  }

  if (verbose) {
    cat("=== generate_glm_dgm() ===\n")
    cat(sprintf("  outcome_type:    %s\n", outcome_type))
    cat(sprintf("  effect_measure:  %s\n", effect_measure))
    cat(sprintf("  adverse_outcome: %s\n", adverse_outcome))
    cat(sprintf("  model:           %s\n", model))
    cat(sprintf("  k_inter:         %.3f\n", k_inter))
    cat(sprintf("  k_treat:         %.3f\n", k_treat))
    cat(sprintf("  n_super:         %d\n", n_super))
    cat(sprintf("  n observations: %d\n", nrow(data)))
  }

  # -- Ensure factor types ---------------------------------------------------
  for (v in factor_vars) {
    if (!is.factor(data[[v]])) {
      data[[v]] <- as.factor(data[[v]])
    }
  }

  # -- Drop single-level factors ---------------------------------------------
  n_levels <- vapply(factor_vars, function(v) {
    length(unique(data[[v]]))
  }, integer(1))

  dropped <- factor_vars[n_levels < 2]
  if (length(dropped) > 0) {
    if (verbose) {
      cat(sprintf("  Dropped single-level factors: %s\n",
          paste(dropped, collapse = ", ")))
    }
    factor_vars <- factor_vars[n_levels >= 2]
  }

  if (length(factor_vars) == 0) {
    stop("No factor variables with 2+ levels remain after filtering.")
  }

  # -- Define subgroup (reuse generate_aft_dgm_helpers.R) --------------------
  # Always create flag_harm with model = "alt" so the subgroup indicator
  # is available for the coefficient fitting pipeline (interaction_term).
  # Under the null, flag_harm is zeroed out after subgroup definition
  # because Q = emptyset (León et al. 2024, Section 3).
  df_work <- data
  df_work$treat <- data[[treatment_var]]

  sg_result <- define_subgroups(
    df_work        = df_work,
    data           = data,
    subgroup_vars  = subgroup_vars,
    subgroup_cuts  = subgroup_cuts,
    continuous_vars = character(0),
    model          = if (!is.null(subgroup_vars)) "alt" else "null",
    verbose        = verbose
  )

  data$flag_harm        <- sg_result$flag_harm
  data$interaction_term <- sg_result$interaction_term

  # Under null: Q = emptyset. Zero out flag_harm so classification
  # metrics are correctly undefined (sens=NA, ppv=0).
  if (model == "null") {
    data$flag_harm <- 0L
  }

  # -- Fit baseline GLM ------------------------------------------------------
  glm_family <- switch(outcome_type,
    binary     = stats::binomial(),
    continuous = stats::gaussian(),
    count      = stats::poisson()
  )

  # -- Fit baseline model WITHOUT interaction ----------------------------------
  # The interaction is added as a direct additive shift to the linear
  # predictor AFTER fitting.  This avoids the non-collapsibility problem
  # in logistic regression where including/excluding covariates changes
  # other coefficient estimates.
  rhs <- c(treatment_var, factor_vars)
  # For count outcomes with offset, include offset in the model.
  # We use the formula approach (offset(log(var))) rather than the
  # offset= argument so that predict() with newdata handles offset
  # correctly without recycling the original offset vector.
  if (outcome_type == "count" && !is.null(offset_var)) {
    rhs <- c(rhs, sprintf("offset(log(%s))", offset_var))
  }

  fml <- stats::as.formula(paste(outcome_var, "~", paste(rhs, collapse = " + ")))
  fit_base <- stats::glm(fml, data = data, family = glm_family)

  if (verbose) {
    cat("\nBaseline GLM fit (no interaction):\n")
    cat(sprintf("  AIC: %.1f\n", stats::AIC(fit_base)))
    coef_label <- switch(outcome_type,
      binary = sprintf("(OR = %.3f)", exp(stats::coef(fit_base)[treatment_var])),
      count  = sprintf("(IRR = %.3f)", exp(stats::coef(fit_base)[treatment_var])),
      continuous = sprintf("(MD = %.3f)", stats::coef(fit_base)[treatment_var])
    )
    cat(sprintf("  Treatment coefficient: %.4f %s\n",
        stats::coef(fit_base)[treatment_var], coef_label))
    if (outcome_type == "count" && !is.null(offset_var)) {
      cat(sprintf("  Offset variable: %s\n", offset_var))
    }
  }

  # -- Determine interaction shift ---------------------------------------------
  # k_inter is a DIRECT additive shift on the linear predictor (log-odds
  # for binary) for Q members under treatment.  This matches the inline
  # approach where beta_inter is added to qlogis(p1) for Q members.
  #   k_inter = 0: no interaction (null model)
  #   k_inter > 0: treatment increases the outcome for Q
  #   k_inter < 0: treatment decreases the outcome for Q
  beta_inter <- if (model == "null") 0 else k_inter

  # Adverse-outcome direction correction (binary only).
  # When Y is adverse (e.g., non-improvement), the fitted beta_treat has the

  # OPPOSITE sign vs a beneficial coding of the same endpoint.  A naive
  # positive k_inter would then amplify the adverse direction, producing a
  # much weaker contrast than the same k_inter on the beneficial scale.
  # Negating beta_inter ensures that positive k_inter consistently amplifies
  # the treatment contrast on the beneficial (1-Y) scale, making
  # calibrated DGMs invariant to outcome coding.
  if (isTRUE(adverse_outcome) && outcome_type == "binary" && beta_inter != 0) {
    beta_inter <- -beta_inter
    if (verbose) {
      cat(sprintf(
        "  adverse_outcome = TRUE: beta_inter negated to %.4f\n",
        beta_inter
      ))
      cat("  (positive k_inter now amplifies contrast on the beneficial 1-Y scale)\n")
    }
  }

  if (verbose) {
    cat(sprintf("  Interaction shift (k_inter): %.4f", k_inter))
    if (model == "null") cat(" (null model: forced to 0)")
    if (isTRUE(adverse_outcome) && outcome_type == "binary" && k_inter != 0) {
      cat(sprintf(" -> beta_inter applied: %.4f", beta_inter))
    }
    cat("\n")
    if (k_treat != 1) {
      cat(sprintf("  Treatment scaling (k_treat): %.4f\n", k_treat))
    }
  }

  # -- Create super-population -----------------------------------------------
  set.seed(seed)
  idx <- sample(nrow(data), n_super, replace = TRUE)
  df_super <- data[idx, , drop = FALSE]
  df_super$id <- seq_len(n_super)
  rownames(df_super) <- NULL

  # -- Compute potential outcomes --------------------------------------------
  # All computation uses the linear predictor (link scale):
  #   binary:     logit scale
  #   continuous: identity scale (link = response)
  #   count:      log scale (includes offset if in formula)
  #
  # Treatment effect scaling:
  #   eta_0 = predict(treatment=0, type="link")
  #   eta_1 = predict(treatment=1, type="link")
  #   treatment_effect = eta_1 - eta_0  (= fitted beta_treat, constant)
  #   eta_1_scaled = eta_0 + k_treat * treatment_effect
  #
  # Then add interaction (k_inter) for Q members under treatment.

  df_0 <- df_super; df_0[[treatment_var]] <- 0L
  df_1 <- df_super; df_1[[treatment_var]] <- 1L

  in_Q <- df_super$flag_harm == 1

  if (outcome_type == "binary") {
    # Link scale: logit
    eta_0 <- stats::predict(fit_base, newdata = df_0, type = "link")
    eta_1 <- stats::predict(fit_base, newdata = df_1, type = "link")

    # Scale treatment effect
    treat_effect <- eta_1 - eta_0   # = beta_treat (constant)
    eta_1 <- eta_0 + k_treat * treat_effect

    # Add interaction for Q under treatment
    eta_1[in_Q] <- eta_1[in_Q] + beta_inter

    df_super$p0 <- stats::plogis(eta_0)
    df_super$p1 <- stats::plogis(eta_1)

  } else if (outcome_type == "continuous") {
    # Identity link: response = link
    mu0 <- stats::predict(fit_base, newdata = df_0, type = "response")
    mu1 <- stats::predict(fit_base, newdata = df_1, type = "response")

    # Scale treatment effect
    treat_effect <- mu1 - mu0       # = beta_treat (constant)
    mu1 <- mu0 + k_treat * treat_effect

    # Add interaction for Q under treatment
    mu1[in_Q] <- mu1[in_Q] + beta_inter

    df_super$mu0 <- mu0
    df_super$mu1 <- mu1

  } else if (outcome_type == "count") {
    # Log link (includes offset when in formula)
    eta_0 <- stats::predict(fit_base, newdata = df_0, type = "link")
    eta_1 <- stats::predict(fit_base, newdata = df_1, type = "link")

    # Scale treatment effect on the log scale
    treat_effect <- eta_1 - eta_0   # = beta_treat (constant)
    eta_1 <- eta_0 + k_treat * treat_effect

    # Add interaction for Q under treatment
    eta_1[in_Q] <- eta_1[in_Q] + beta_inter

    df_super$mu0 <- exp(eta_0)
    df_super$mu1 <- exp(eta_1)
  }

  # -- Compute true effects --------------------------------------------------
  effects <- .compute_glm_effects(
    df_super, outcome_type, effect_measure, verbose
  )

  # -- Store residual SD (for future continuous support) ---------------------
  sigma_resid <- if (outcome_type == "continuous") {
    stats::sigma(fit_base)
  } else {
    NULL
  }

  # -- Assemble result -------------------------------------------------------
  result <- list(
    # Super-population
    df_super = df_super,

    # True effects (named for reporting compatibility)
    hazard_ratios = list(
      harm_subgroup    = effects$effect_Q,
      no_harm_subgroup = effects$effect_Qc,
      overall          = effects$effect_ITT
    ),

    # GLM metadata
    outcome_type   = outcome_type,
    effect_measure = effect_measure,

    # Model parameters
    model_params = list(
      coefficients = stats::coef(fit_base),
      family       = glm_family,
      k_treat      = k_treat,
      k_inter      = k_inter,
      beta_inter   = beta_inter,
      sigma        = sigma_resid,
      formula      = fml,
      offset_var   = offset_var
    ),

    # Subgroup information
    subgroup_info = list(
      vars        = subgroup_vars,
      cuts        = subgroup_cuts,
      definitions = sg_result$definitions,
      size        = sum(df_super$flag_harm),
      proportion  = mean(df_super$flag_harm)
    ),

    # Covariate info (factors that survived filtering)
    factor_vars = factor_vars,

    model_type      = model,
    model           = model,       # For downstream null detection
    adverse_outcome = adverse_outcome,
    n_super         = n_super,
    seed            = seed
  )

  class(result) <- c("glm_dgm", "list")

  if (verbose) {
    cat("\n=== DGM Summary ===\n")
    cat(sprintf("  Effect(Q):   %.4f  (%s)\n",
        effects$effect_Q, effect_measure))
    cat(sprintf("  Effect(Qc):  %.4f  (%s)\n",
        effects$effect_Qc, effect_measure))
    cat(sprintf("  Effect(ITT): %.4f  (%s)\n",
        effects$effect_ITT, effect_measure))
    cat(sprintf("  Subgroup prevalence: %.1f%%\n",
        100 * mean(df_super$flag_harm)))
    if (isTRUE(adverse_outcome) && outcome_type == "binary" &&
        effect_measure %in% c("OR", "RR")) {
      cat(sprintf(
        "  Beneficial-scale equivalent: Effect(Q) = %.4f, Effect(Qc) = %.4f\n",
        1 / effects$effect_Q, 1 / effects$effect_Qc
      ))
    }
  }


  result
}


# =============================================================================
# Internal: compute true effects from super-population potential outcomes
# =============================================================================

#' @keywords internal
.compute_glm_effects <- function(df_super, outcome_type, effect_measure,
                                  verbose = FALSE) {

  .effect_one <- function(df) {
    if (outcome_type == "binary") {
      p1 <- mean(df$p1)
      p0 <- mean(df$p0)
      switch(effect_measure,
        "OR" = (p1 / (1 - p1)) / (p0 / (1 - p0)),
        "RD" = p1 - p0,
        "RR" = p1 / p0,
        stop("Unknown binary effect_measure: ", effect_measure)
      )
    } else if (outcome_type == "continuous") {
      mu1 <- mean(df$mu1)
      mu0 <- mean(df$mu0)
      switch(effect_measure,
        "MD" = mu1 - mu0,
        stop("Unknown continuous effect_measure: ", effect_measure)
      )
    } else if (outcome_type == "count") {
      mu1 <- mean(df$mu1)
      mu0 <- mean(df$mu0)
      switch(effect_measure,
        "IRR" = mu1 / mu0,
        "IRD" = mu1 - mu0,
        stop("Unknown count effect_measure: ", effect_measure)
      )
    }
  }

  in_Q <- df_super$flag_harm == 1

  effect_Q   <- if (sum(in_Q) > 0)  .effect_one(df_super[in_Q, ])  else NA_real_
  effect_Qc  <- if (sum(!in_Q) > 0) .effect_one(df_super[!in_Q, ]) else NA_real_
  effect_ITT <- .effect_one(df_super)

  list(
    effect_Q   = effect_Q,
    effect_Qc  = effect_Qc,
    effect_ITT = effect_ITT
  )
}


# =============================================================================
# Print method
# =============================================================================

#' Print method for \code{glm_dgm} objects
#' @param x A \code{glm_dgm} object from \code{\link{generate_glm_dgm}}.
#' @param ... Unused; present for S3 compatibility.
#' @return The input \code{x}, invisibly.
#' @export
print.glm_dgm <- function(x, ...) {
  cat("GLM Data Generating Mechanism\n")
  cat(sprintf("  outcome_type:   %s\n", x$outcome_type))
  cat(sprintf("  effect_measure: %s\n", x$effect_measure))
  cat(sprintf("  model:          %s\n", x$model_type))
  cat(sprintf("  n_super:        %d\n", x$n_super))
  cat(sprintf("  Effect(Q):      %.4f\n", x$hazard_ratios$harm_subgroup))
  cat(sprintf("  Effect(Qc):     %.4f\n", x$hazard_ratios$no_harm_subgroup))
  cat(sprintf("  Effect(ITT):    %.4f\n", x$hazard_ratios$overall))
  if (!is.null(x$subgroup_info)) {
    cat(sprintf("  Subgroup:       %d / %d (%.1f%%)\n",
        x$subgroup_info$size, x$n_super,
        100 * x$subgroup_info$proportion))
  }
  invisible(x)
}
