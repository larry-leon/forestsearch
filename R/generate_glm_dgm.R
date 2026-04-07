# =============================================================================
# generate_glm_dgm() -- GLM Data Generating Mechanism
# =============================================================================
#
# Mirrors generate_aft_dgm_flex() for non-survival outcomes.
# Phase 1: binary endpoints (logistic regression).
# Phase 2/3: continuous and count endpoints (planned).
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
#' @param outcome_type Character. Currently only \code{"binary"} is
#'   supported (Phase 1). Future phases will add \code{"continuous"} and
#'   \code{"count"}.
#' @param effect_measure Character. Effect measure for the GLM.
#'   Default is \code{NULL}, which selects \code{"OR"} for binary outcomes.
#' @param subgroup_vars Character vector of variable names defining the
#'   subgroup. Default \code{NULL} (no subgroup).
#' @param subgroup_cuts Named list of cutpoint specifications for
#'   subgroup variables. Uses the same flexible format as
#'   \code{\link{generate_aft_dgm_flex}}: fixed values, quantiles, or
#'   functions.
#' @param model Character. \code{"alt"} for alternative hypothesis with
#'   treatment--subgroup interaction, \code{"null"} for uniform treatment
#'   effect. Default \code{"alt"}.
#' @param k_inter Numeric interaction-effect modifier applied to the
#'   treatment x subgroup coefficient. Values > 1 strengthen the
#'   interaction; 0 removes it. Default \code{1}.
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
#'       outcome columns \code{p0}, \code{p1} (binary), and
#'       \code{flag_harm}.}
#'     \item{\code{hazard_ratios}}{Named list with \code{harm_subgroup},
#'       \code{no_harm_subgroup}, and \code{overall} -- effect estimates
#'       on the scale determined by \code{effect_measure}.  Field name
#'       retained for compatibility with \code{get_dgm_hr()} and
#'       reporting functions.}
#'     \item{\code{outcome_type}}{Character: \code{"binary"}.}
#'     \item{\code{effect_measure}}{Character: the effect measure used.}
#'     \item{\code{model_params}}{List with fitted coefficients, family,
#'       \code{k_inter}, and residual SD (for continuous outcomes).}
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
#'   k_inter       = 1.5,
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
    outcome_type    = c("binary"),
    effect_measure  = NULL,
    subgroup_vars   = NULL,
    subgroup_cuts   = NULL,
    model           = c("alt", "null"),
    k_inter         = 1,
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

  if (verbose) {
    cat("=== generate_glm_dgm() ===\n")
    cat(sprintf("  outcome_type:   %s\n", outcome_type))
    cat(sprintf("  effect_measure: %s\n", effect_measure))
    cat(sprintf("  model:          %s\n", model))
    cat(sprintf("  k_inter:        %.3f\n", k_inter))
    cat(sprintf("  n_super:        %d\n", n_super))
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
  # is available for classification metrics under both hypotheses.
  # The null model zeros out the interaction coefficient later.
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

  # -- Fit baseline GLM ------------------------------------------------------
  glm_family <- switch(outcome_type,
    binary     = stats::binomial(),
    continuous = stats::gaussian(),
    count      = stats::poisson()
  )

  # Build formula: outcome ~ treatment + factors [+ interaction]
  rhs <- c(treatment_var, factor_vars)
  # Always include the interaction when a subgroup is defined, so the
  # coefficient exists in both alt (k_inter > 0) and null (k_inter = 0).
  has_subgroup <- !is.null(subgroup_vars) && sum(data$flag_harm) > 0
  if (has_subgroup) {
    data$treat_x_sg <- data[[treatment_var]] * data$flag_harm
    rhs <- c(rhs, "treat_x_sg")
  }

  fml <- stats::as.formula(paste(outcome_var, "~", paste(rhs, collapse = " + ")))
  fit_base <- stats::glm(fml, data = data, family = glm_family)

  if (verbose) {
    cat("\nBaseline GLM fit:\n")
    cat(sprintf("  AIC: %.1f\n", stats::AIC(fit_base)))
    cat(sprintf("  Treatment coefficient: %.4f\n",
        stats::coef(fit_base)[treatment_var]))
    if ("treat_x_sg" %in% names(stats::coef(fit_base))) {
      cat(sprintf("  Interaction coefficient: %.4f\n",
          stats::coef(fit_base)["treat_x_sg"]))
    }
  }

  # -- Modify interaction coefficient ----------------------------------------
  beta <- stats::coef(fit_base)
  if ("treat_x_sg" %in% names(beta)) {
    if (model == "null") {
      beta["treat_x_sg"] <- 0
      if (verbose) cat("  Null model: interaction set to 0\n")
    } else {
      beta["treat_x_sg"] <- beta["treat_x_sg"] * k_inter
      if (verbose) {
        cat(sprintf("  Modified interaction (k_inter = %.3f): %.4f\n",
            k_inter, beta["treat_x_sg"]))
      }
    }
  }

  # -- Create super-population -----------------------------------------------
  set.seed(seed)
  idx <- sample(nrow(data), n_super, replace = TRUE)
  df_super <- data[idx, , drop = FALSE]
  df_super$id <- seq_len(n_super)
  rownames(df_super) <- NULL

  # -- Compute potential outcomes --------------------------------------------
  # Build model matrices for control and treatment
  df_0 <- df_super
  df_0[[treatment_var]] <- 0L
  if ("treat_x_sg" %in% names(df_0)) df_0$treat_x_sg <- 0

  df_1 <- df_super
  df_1[[treatment_var]] <- 1L
  if ("treat_x_sg" %in% names(df_1)) {
    df_1$treat_x_sg <- 1L * df_1$flag_harm
  }

  # Compute linear predictors using modified coefficients
  X_0 <- stats::model.matrix(fml, data = df_0)
  X_1 <- stats::model.matrix(fml, data = df_1)

  # Ensure coefficient vector aligns with model matrix columns
  beta_aligned <- beta[colnames(X_0)]

  eta_0 <- as.numeric(X_0 %*% beta_aligned)
  eta_1 <- as.numeric(X_1 %*% beta_aligned)

  if (outcome_type == "binary") {
    df_super$p0 <- stats::plogis(eta_0)
    df_super$p1 <- stats::plogis(eta_1)
  }
  # Future: continuous -> mu0/mu1; count -> mu0/mu1 with offset

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
      coefficients = beta,
      family       = glm_family,
      k_inter      = k_inter,
      sigma        = sigma_resid,
      formula      = fml
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

    model_type = model,
    n_super    = n_super,
    seed       = seed
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
        stop("Unknown effect_measure: ", effect_measure)
      )
    }
    # Future: continuous -> mean(mu1) - mean(mu0); count -> exp(...)
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
