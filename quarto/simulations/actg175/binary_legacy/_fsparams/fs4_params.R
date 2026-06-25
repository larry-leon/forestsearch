# =============================================================================
# fs4_params.R
#
# Canonical FS parameter bundle "fs4" for ACTG175 binary harm simulations.
#
# Sourced by per-config qmds at:
#   quarto/simulations/actg175/actg175_binary_m1_harm_<focus>_fs4.qmd
#
# Provides one function: get_fs4_params(sg_focus, ...)
#   Returns a list of FS parameters with sg_focus injected and
#   stop_threshold set to 0.90 (or NULL when sg_focus is "hrMaxSG"
#   or "hrMinSG", a package-level constraint). The calling qmd uses
#   the returned list directly with no overrides.
#
# -----------------------------------------------------------------------------
# Regime: fs4 — permissive search, broad candidate space, full cut retention.
#
# Designed for: exploratory analyses where the goal is high sensitivity to
# harm subgroups; tolerates a higher false positive rate as the cost.
#
# Key choices:
#   max_subgroups_search = 10  → broad candidate space (vs. fs4's 5)
#   return_selected_cuts_only = FALSE  → retain all cut data per call
#   pconsistency.threshold = 0.80  → moderate consistency requirement
#   c1 (hr.threshold) = 1.25  → OR scale; modest signal threshold
#   stop_threshold = 0.90  → early-stopping gate (returned as NULL when
#                            sg_focus is hrMaxSG or hrMinSG; rule
#                            applied inside get_fs4_params())
#
# To create a sibling bundle "fs4" with different parameter choices:
#   1. Copy this file to fs4_params.R
#   2. Rename the function to get_fs4_params()
#   3. Edit the values that differ
#   4. Update the regime description above
#   5. Use clone_config.sh fs4 fs4 to spawn the parallel sg_focus qmds
# =============================================================================

#' Get FS parameters for the "fs4" bundle
#'
#' Returns the FS parameter list used by ACTG175 binary harm simulations
#' under the "fs4" regime. The bundle's \code{stop_threshold} is set to
#' 0.90, with one exception: when \code{sg_focus} is \code{"hrMaxSG"} or
#' \code{"hrMinSG"}, \code{stop_threshold} is returned as \code{NULL}
#' (a package-level constraint). The calling qmd uses the returned list
#' directly. The qmd is still responsible for injecting \code{conf_force}
#' after data prep defines the necessary factor variables.
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
get_fs4_params <- function(sg_focus,
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
    max_subgroups_search      = 50,        # fs4: broad
    use_twostage              = TRUE,

    # -- Effect-size thresholds ----------------------------------------------
    hr.threshold              = 1.5,      # OR scale
    hr.consistency            = 1.0,
    pconsistency.threshold    = 0.80,
    # stop_threshold: 0.90 is the fs4 bundle's value. Package requires
    # NULL when sg_focus uses neighborhood-based selection (hrMaxSG /
    # hrMinSG); applied here so qmds can call get_fs4_params() and use
    # the result directly with no overrides.
    stop_threshold            = if (sg_focus %in% c("hrMaxSG", "hrMinSG"))
                                  NULL else 0.95,

    # -- Search dispatch -----------------------------------------------------
    sg_focus                  = sg_focus,
    conf_force             = c("ar_naive == 1","prior_6mo == 1"),
    conf.cont_jcuts        = list(cd40 = 12, wtkg = 12, preanti =12),

    # -- Sample-size constraints ---------------------------------------------
    fs.splits                 = 1000L,
    # Assuming potential subgroup is at least 10% of the sample size
    n.min                     = 100L,
    d0.min                    = 10L,
    d1.min                    = 0L,
    maxk                      = 2L,
    dmin.grf                  = 0.0,

    # -- RCT / seed ----------------------------------------------------------
    is.RCT                    = TRUE,
    seedit                    = seedit,

    # -- Output / verbosity --------------------------------------------------
    details                   = TRUE,
    quiet                     = FALSE,

    # -- Inner-loop parallel plan (forced sequential) -----------------------
    parallel_args             = list(plan         = "sequential",
                                     workers      = 1,
                                     show_message = FALSE)
  )
}
