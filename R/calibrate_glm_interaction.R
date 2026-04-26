# =============================================================================
# calibrate_glm_interaction() -- Find k_inter for a target effect size
# =============================================================================
#
# Replaces the ad-hoc grid search + add_interaction() pattern used in Quarto
# simulation documents.
# =============================================================================


#' Calibrate GLM Interaction for a Target Subgroup Effect Size
#'
#' Searches over a grid of interaction multipliers (\code{k_inter}) to find
#' the value that produces a target effect size in the subgroup \eqn{Q}.
#' Returns a modified \code{"glm_dgm"} object with the calibrated
#' interaction and updated true effects.
#'
#' @param data The source data frame (same as passed to
#'   \code{\link{generate_glm_dgm}}).
#' @param factor_vars Character vector of factor variable names.
#' @param outcome_var Character string naming the outcome variable.
#' @param treatment_var Character string naming the treatment variable.
#' @param target_effect Numeric. Target effect size in Q.  When
#'   \code{adverse_outcome = FALSE}, this is the target on the
#'   raw Y scale (e.g., OR = 2.0 means treatment doubles the
#'   odds of Y = 1 in Q).  When \code{adverse_outcome = TRUE}
#'   and \code{effect_measure} is a ratio (OR, RR, IRR), the
#'   target is interpreted on the \emph{beneficial} (1-Y) scale
#'   and automatically inverted (\code{1 / target_effect}) for
#'   the grid search, so that \code{target_effect = 2.0} produces
#'   the same biological heterogeneity regardless of outcome coding.
#' @param outcome_type Character. One of \code{"binary"},
#'   \code{"continuous"}, or \code{"count"}.
#' @param effect_measure Character. Effect measure. Default \code{NULL}
#'   (selects default per \code{outcome_type}).
#' @param offset_var Character or \code{NULL}. Name of the exposure /
#'   follow-up time column for count outcomes with an offset (IRR, IRD).
#'   Default \code{NULL}.
#' @param subgroup_vars Character vector of subgroup-defining variables.
#' @param subgroup_cuts Named list of cutpoint specifications.
#' @param k_treat Numeric. Treatment effect scaling factor passed to
#'   \code{\link{generate_glm_dgm}}.  Default \code{1} (preserve fitted
#'   treatment effect).  Set to \code{0} to zero out the ITT effect.
#' @param adverse_outcome Logical.  Passed through to
#'   \code{\link{generate_glm_dgm}}.  When \code{TRUE} for binary
#'   outcomes with ratio measures, \code{target_effect} is interpreted
#'   on the beneficial (1-Y) scale and inverted before searching.
#'   Default \code{FALSE}.
#' @param k_inter_range Numeric vector of length 2. Search range for
#'   \code{k_inter}. Default \code{c(-10, 10)}.
#' @param grid_step Numeric. Grid resolution. Default \code{0.05}.
#' @param n_super Integer. Super-population size. Default \code{5000L}.
#' @param seed Integer. Random seed. Default \code{8316951L}.
#' @param verbose Logical. Print calibration progress. Default \code{FALSE}.
#'
#' @return An object of class \code{"glm_dgm"} with the calibrated
#'   \code{k_inter} and updated \code{hazard_ratios}.
#'
#' @details
#' The function constructs a grid of \code{k_inter} values, calls
#' \code{\link{generate_glm_dgm}} for each, and selects the value
#' whose subgroup effect is closest to \code{target_effect}.
#'
#' When \code{adverse_outcome = TRUE} for binary outcomes with a ratio
#' measure (OR, RR), the target is automatically inverted
#' (\code{1 / target_effect}) for the grid search.  This ensures that
#' specifying \code{target_effect = 2.0} always means "the treatment
#' contrast in Q is equivalent to OR = 2.0 on the beneficial scale,"
#' regardless of whether Y is coded as adverse or beneficial.  The
#' returned \code{hazard_ratios} are always on the \emph{actual Y
#' scale} (i.e., on the adverse scale when Y is adverse).
#'
#' @seealso \code{\link{generate_glm_dgm}},
#'   \code{\link{simulate_from_glm_dgm}}
#'
#' @examples
#' \dontrun{
#' # Beneficial outcome (Y = improvement) -- target OR(Q) = 2.0 on Y scale
#' dgm_benefit <- calibrate_glm_interaction(
#'   data = actg_df, factor_vars = paste0("z", 1:12),
#'   outcome_var = "y_improve", treatment_var = "treat",
#'   target_effect = 2.0, outcome_type = "binary",
#'   subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
#'   verbose = TRUE
#' )
#'
#' # Adverse outcome (Y = non-improvement) -- same biological target
#' # target_effect = 2.0 is automatically inverted to 0.5 on the adverse scale
#' dgm_adverse <- calibrate_glm_interaction(
#'   data = actg_df, factor_vars = paste0("z", 1:12),
#'   outcome_var = "y_adverse", treatment_var = "treat",
#'   target_effect = 2.0, outcome_type = "binary",
#'   adverse_outcome = TRUE,
#'   subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
#'   verbose = TRUE
#' )
#' # Both DGMs produce the same heterogeneity (same k_inter)
#' }
#'
#' @export
calibrate_glm_interaction <- function(
    data,
    factor_vars,
    outcome_var,
    treatment_var,
    target_effect,
    outcome_type    = c("binary", "continuous", "count"),
    effect_measure  = NULL,
    offset_var      = NULL,
    subgroup_vars   = NULL,
    subgroup_cuts   = NULL,
    k_treat         = 1,
    adverse_outcome = FALSE,
    k_inter_range   = c(-10, 10),
    grid_step       = 0.05,
    n_super         = 5000L,
    seed            = 8316951L,
    verbose         = FALSE
) {

  outcome_type <- match.arg(outcome_type)

  if (is.null(subgroup_vars) || is.null(subgroup_cuts)) {
    stop("subgroup_vars and subgroup_cuts are required for calibration.")
  }

  # -- Resolve effect measure for target-inversion logic ---------------------
  if (is.null(effect_measure)) {
    em_resolved <- switch(outcome_type,
      binary     = "OR",
      continuous = "MD",
      count      = "IRR"
    )
  } else {
    em_resolved <- effect_measure
  }

  # -- Adverse-outcome target inversion (ratio measures only) ----------------
  # When Y is adverse and the measure is a ratio (OR, RR, IRR), the user's
  # target_effect is on the *beneficial* (1-Y) scale.  Since
  # generate_glm_dgm() with adverse_outcome = TRUE negates beta_inter,
  # the calibrated DGM will have OR(Q) = 1/target_effect on the actual
  # (adverse) Y scale.  We search for that inverted target.
  is_ratio <- em_resolved %in% c("OR", "RR", "IRR")
  if (isTRUE(adverse_outcome) && outcome_type == "binary" && is_ratio) {
    search_target <- 1 / target_effect
    if (verbose) {
      cat(sprintf(
        "adverse_outcome = TRUE: user target %s(Q) = %.3f on beneficial scale\n",
        em_resolved, target_effect
      ))
      cat(sprintf(
        "  -> searching for %s(Q) = %.4f on the adverse (actual Y) scale\n",
        em_resolved, search_target
      ))
    }
  } else if (isTRUE(adverse_outcome) && outcome_type == "binary" &&
             em_resolved == "RD") {
    search_target <- -target_effect
    if (verbose) {
      cat(sprintf(
        "adverse_outcome = TRUE: user target RD(Q) = %.3f -> searching %.4f\n",
        target_effect, search_target
      ))
    }
  } else {
    search_target <- target_effect
  }

  # -- Build grid ------------------------------------------------------------
  k_grid <- seq(k_inter_range[1], k_inter_range[2], by = grid_step)

  if (verbose) {
    cat(sprintf("Calibrating: target %s = %.3f",
        em_resolved, target_effect))
    if (search_target != target_effect) {
      cat(sprintf(" (search target on actual Y scale: %.4f)", search_target))
    }
    cat("\n")
    cat(sprintf("  Grid: %.2f to %.2f by %.3f (%d values)\n",
        k_inter_range[1], k_inter_range[2], grid_step, length(k_grid)))
  }

  # -- Evaluate each k_inter -------------------------------------------------
  effects_Q <- vapply(k_grid, function(k) {
    dgm_k <- generate_glm_dgm(
      data            = data,
      factor_vars     = factor_vars,
      outcome_var     = outcome_var,
      treatment_var   = treatment_var,
      outcome_type    = outcome_type,
      effect_measure  = effect_measure,
      offset_var      = offset_var,
      subgroup_vars   = subgroup_vars,
      subgroup_cuts   = subgroup_cuts,
      model           = "alt",
      k_treat         = k_treat,
      k_inter         = k,
      adverse_outcome = adverse_outcome,
      n_super         = n_super,
      seed            = seed,
      verbose         = FALSE
    )
    dgm_k$hazard_ratios$harm_subgroup
  }, numeric(1))

  # -- Find best -------------------------------------------------------------
  best_idx <- which.min(abs(effects_Q - search_target))
  best_k   <- k_grid[best_idx]
  best_eff <- effects_Q[best_idx]

  if (verbose) {
    cat(sprintf(
      "  Best: k_inter = %.3f -> Effect(Q) = %.4f (search target: %.4f)\n",
      best_k, best_eff, search_target
    ))
    if (search_target != target_effect && is_ratio) {
      cat(sprintf(
        "  Beneficial-scale equivalent: Effect(Q) = %.4f (user target: %.3f)\n",
        1 / best_eff, target_effect
      ))
    }
  }

  # -- Build final DGM with the calibrated k_inter --------------------------
  dgm_cal <- generate_glm_dgm(
    data            = data,
    factor_vars     = factor_vars,
    outcome_var     = outcome_var,
    treatment_var   = treatment_var,
    outcome_type    = outcome_type,
    effect_measure  = effect_measure,
    offset_var      = offset_var,
    subgroup_vars   = subgroup_vars,
    subgroup_cuts   = subgroup_cuts,
    model           = "alt",
    k_treat         = k_treat,
    k_inter         = best_k,
    adverse_outcome = adverse_outcome,
    n_super         = n_super,
    seed            = seed,
    verbose         = verbose
  )

  if (verbose) {
    cat(sprintf("\nCalibration complete: k_inter = %.3f\n", best_k))
    cat(sprintf("  Effect(Q):   %.4f  (actual Y scale)\n",
        dgm_cal$hazard_ratios$harm_subgroup))
    cat(sprintf("  Effect(Qc):  %.4f\n",
        dgm_cal$hazard_ratios$no_harm_subgroup))
    cat(sprintf("  Effect(ITT): %.4f\n",
        dgm_cal$hazard_ratios$overall))
    if (isTRUE(adverse_outcome) && outcome_type == "binary" && is_ratio) {
      cat(sprintf(
        "  Beneficial-scale: Effect(Q) = %.4f, Effect(Qc) = %.4f\n",
        1 / dgm_cal$hazard_ratios$harm_subgroup,
        1 / dgm_cal$hazard_ratios$no_harm_subgroup
      ))
    }
  }

  dgm_cal
}


# Null-coalescing helper (local)
`%||%` <- function(a, b) if (!is.null(a)) a else b
