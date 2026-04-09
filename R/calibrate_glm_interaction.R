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
#' @param target_effect Numeric. Target effect size in Q on the scale
#'   determined by \code{effect_measure} (e.g., OR = 2.0 for a binary
#'   outcome).
#' @param outcome_type Character. One of \code{"binary"},
#'   \code{"continuous"}, or \code{"count"}.
#' @param effect_measure Character. Effect measure. Default \code{NULL}
#'   (selects default per \code{outcome_type}).
#' @param offset_var Character or \code{NULL}. Name of the exposure /
#'   follow-up time column for count outcomes with an offset (IRR, IRD).
#'   Default \code{NULL}.
#' @param subgroup_vars Character vector of subgroup-defining variables.
#' @param subgroup_cuts Named list of cutpoint specifications.
#' @param k_inter_range Numeric vector of length 2. Search range for
#'   \code{k_inter}. Default \code{c(0, 10)}.
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
#' For binary outcomes with \code{effect_measure = "OR"}, the target
#' is on the switched-treatment OR scale (OR > 1 = treatment increases
#' the outcome).
#'
#' @seealso \code{\link{generate_glm_dgm}},
#'   \code{\link{simulate_from_glm_dgm}}
#'
#' @examples
#' \dontrun{
#' dgm <- calibrate_glm_interaction(
#'   data          = actg_df,
#'   factor_vars   = paste0("z", 1:12),
#'   outcome_var   = "y_binary",
#'   treatment_var = "treat",
#'   target_effect = 2.0,
#'   outcome_type  = "binary",
#'   subgroup_vars = c("z1", "z2"),
#'   subgroup_cuts = list(z1 = 1L, z2 = 1L),
#'   verbose       = TRUE
#' )
#' print(dgm)
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

  # -- Build grid ------------------------------------------------------------
  k_grid <- seq(k_inter_range[1], k_inter_range[2], by = grid_step)

  if (verbose) {
    cat(sprintf("Calibrating: target %s = %.3f\n",
        effect_measure %||% "OR", target_effect))
    cat(sprintf("  Grid: %.2f to %.2f by %.3f (%d values)\n",
        k_inter_range[1], k_inter_range[2], grid_step, length(k_grid)))
  }

  # -- Evaluate each k_inter -------------------------------------------------
  effects_Q <- vapply(k_grid, function(k) {
    dgm_k <- generate_glm_dgm(
      data          = data,
      factor_vars   = factor_vars,
      outcome_var   = outcome_var,
      treatment_var = treatment_var,
      outcome_type  = outcome_type,
      effect_measure = effect_measure,
      offset_var    = offset_var,
      subgroup_vars = subgroup_vars,
      subgroup_cuts = subgroup_cuts,
      model         = "alt",
      k_inter       = k,
      n_super       = n_super,
      seed          = seed,
      verbose       = FALSE
    )
    dgm_k$hazard_ratios$harm_subgroup
  }, numeric(1))

  # -- Find best -------------------------------------------------------------
  best_idx <- which.min(abs(effects_Q - target_effect))
  best_k   <- k_grid[best_idx]
  best_eff <- effects_Q[best_idx]

  if (verbose) {
    cat(sprintf("  Best: k_inter = %.3f -> Effect(Q) = %.4f (target: %.3f)\n",
        best_k, best_eff, target_effect))
  }

  # -- Build final DGM with the calibrated k_inter --------------------------
  dgm_cal <- generate_glm_dgm(
    data          = data,
    factor_vars   = factor_vars,
    outcome_var   = outcome_var,
    treatment_var = treatment_var,
    outcome_type  = outcome_type,
    effect_measure = effect_measure,
    offset_var    = offset_var,
    subgroup_vars = subgroup_vars,
    subgroup_cuts = subgroup_cuts,
    model         = "alt",
    k_inter       = best_k,
    n_super       = n_super,
    seed          = seed,
    verbose       = verbose
  )

  if (verbose) {
    cat(sprintf("\nCalibration complete: k_inter = %.3f\n", best_k))
    cat(sprintf("  Effect(Q):   %.4f\n", dgm_cal$hazard_ratios$harm_subgroup))
    cat(sprintf("  Effect(Qc):  %.4f\n", dgm_cal$hazard_ratios$no_harm_subgroup))
    cat(sprintf("  Effect(ITT): %.4f\n", dgm_cal$hazard_ratios$overall))
  }

  dgm_cal
}


# Null-coalescing helper (local)
`%||%` <- function(a, b) if (!is.null(a)) a else b
