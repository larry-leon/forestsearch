# =============================================================================
# fs1_params.R
#
# Canonical FS parameter bundle "fs1" for ACTG175 binary harm simulations.
#
# Sourced by per-config qmds at:
#   quarto/simulations/actg175/actg175_binary_m1_harm_<focus>_fs1.qmd
#
# Provides one function: get_fs1_params(sg_focus, ...)
#   Returns a list of FS parameters with sg_focus injected.
#   Note: stop_threshold is NOT set here — the calling qmd must set it
#   based on its own sg_focus (NULL for hrMaxSG/hrMinSG, 0.90 otherwise).
#   This is "Approach B": stop_threshold visible in the qmd.
#
# -----------------------------------------------------------------------------
# Regime: fs1 — permissive search, broad candidate space, full cut retention.
#
# Designed for: exploratory analyses where the goal is high sensitivity to
# harm subgroups; tolerates a higher false positive rate as the cost.
#
# Key choices:
#   max_subgroups_search = 10  → broad candidate space (vs. fs2's 5)
#   return_selected_cuts_only = FALSE  → retain all cut data per call
#   pconsistency.threshold = 0.80  → moderate consistency requirement
#   c1 (hr.threshold) = 1.25  → OR scale; modest signal threshold
#
# To create a sibling bundle "fs2" with different parameter choices:
#   1. Copy this file to fs2_params.R
#   2. Rename the function to get_fs2_params()
#   3. Edit the values that differ
#   4. Update the regime description above
#   5. Use clone_config.sh fs1 fs2 to spawn the parallel sg_focus qmds
# =============================================================================

#' Get FS parameters for the "fs1" bundle
#'
#' Returns the FS parameter list used by ACTG175 binary harm simulations
#' under the "fs1" regime. The calling qmd is responsible for setting
#' \code{stop_threshold} based on its own \code{sg_focus} (NULL for
#' \code{hrMaxSG}/\code{hrMinSG}, 0.90 otherwise) and for injecting
#' \code{conf_force} after data prep defines the necessary factor
#' variables.
#'
#' @param sg_focus Character. The sg_focus value (e.g., "hr", "hrMaxSG",
#'   "maxSG"). Injected into the parameter list as-is.
#' @param outcome_name Character. Outcome variable name in the simulated
#'   data. Default \code{"y_sim"} matches \code{simulate_from_glm_dgm()}.
#' @param treat_name Character. Treatment variable name. Default
#'   \code{"treat_sim"}.
#' @param id_name Character. Subject identifier name. Default \code{"id"}.
#' @param seedit Integer. Master seed. Default \code{8316951L}.
#'
#' @return A named list of FS parameters suitable for passing to
#'   \code{forestsearch()} (single-trial) or
#'   \code{run_simulation_analysis(fs_params = ...)} (simulations).
get_fs1_params <- function(sg_focus,
                           outcome_name = "y_sim",
                           treat_name   = "treat_sim",
                           id_name      = "id",
                           seedit       = 8316951L) {
  list(
    # -- Identifiers (data column names) -------------------------------------
    outcome.name              = outcome_name,
    event.name                = outcome_name,
    treat.name                = treat_name,
    id.name                   = id_name,

    # -- Outcome typing ------------------------------------------------------
    outcome_type              = "binary",
    effect_measure            = "OR",

    # -- Search controls -----------------------------------------------------
    use_lasso                 = FALSE,
    use_grf                   = TRUE,
    return_selected_cuts_only = FALSE,
    max_subgroups_search      = 10,        # fs1: broad
    use_twostage              = TRUE,

    # -- Effect-size thresholds ----------------------------------------------
    hr.threshold              = 1.25,      # OR scale
    hr.consistency            = 1.0,
    pconsistency.threshold    = 0.80,

    # -- Search dispatch -----------------------------------------------------
    sg_focus                  = sg_focus,  # injected from caller
    # NOTE: stop_threshold is intentionally NOT set here. Caller sets it
    # based on its own sg_focus value.

    # -- Sample-size constraints ---------------------------------------------
    fs.splits                 = 1000L,
    n.min                     = 60L,
    d0.min                    = 10L,
    d1.min                    = 0L,
    maxk                      = 2L,
    dmin.grf                  = 0.0,

    # -- RCT / seed ----------------------------------------------------------
    is.RCT                    = TRUE,
    seedit                    = seedit,

    # -- Output / verbosity --------------------------------------------------
    showten_subgroups         = TRUE,
    details                   = TRUE,
    quiet                     = FALSE,

    # -- Inner-loop parallel plan (forced sequential) -----------------------
    parallel_args             = list(plan         = "sequential",
                                     workers      = 1,
                                     show_message = FALSE)
  )
}
