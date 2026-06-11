# ============================================================================
# DINA single-trial diagnostic
# ----------------------------------------------------------------------------
# Paste this as a NEW chunk immediately AFTER the setup chunk of
# dina-tier2_gate_standalone_simulation.qmd (so dgm, confounders_base, the
# column-name objects, and all knobs are in scope), then render or run once.
#
# It runs ONE forestsearch() DINA call WITHOUT the harness's error-swallowing
# tryCatch and with details = TRUE / quiet = FALSE, so the real error message
# (or the DINA selection trace) is printed instead of being collapsed to an
# all-NA record.  Try it with dina_select_statistic = "effect" first, then
# flip the knob to "dina" to localise whether the effect re-rank is involved.
# ============================================================================

.dx_df    <- simulate_from_dgm(dgm, n = n_sample, analysis_time = analysis_time,
                               cens_adjust = cens_adjust, seed = seed_base + 1L)
.dx_confs <- intersect(confounders_base, names(.dx_df))
cat(sprintf("n = %d | events = %d | arms = %d | confs = %s\n\n",
            nrow(.dx_df), sum(.dx_df[[event_name]]),
            length(unique(.dx_df[[treat_name]])),
            paste(.dx_confs, collapse = ", ")))

.dx <- tryCatch(
  forestsearch(
    df.analysis            = .dx_df,
    outcome.name           = outcome_name, event.name = event_name,
    treat.name             = treat_name,   id.name    = id_name,
    flag_harm.name         = harm_col,     confounders.name = .dx_confs,
    is.RCT                 = TRUE,         seedit = seed_base + 1L,
    sg_focus               = sg_focus,    subgroup_method = "dina",
    hr.threshold           = hr_threshold, hr.consistency = hr_consistency,
    pconsistency.threshold = pconsistency, n.min = n_min,
    selection_rule         = selection_rule,
    parallel_args          = list(plan = "sequential"),
    debias_gate            = TRUE,
    debias_gate_args       = list(ci_method = "ij", draws = gate_draws,
                                  include_complement = TRUE),
    dina_args = modifyList(dina_args, list(select_statistic = dina_select_statistic)),
    details = TRUE, quiet = FALSE),                # <- expose DINA's decisions
  error = function(e) {
    message("\n>>> forestsearch() ERROR: ", conditionMessage(e))
    message(">>> call: ", deparse(conditionCall(e))[1])
    NULL
  })

if (!is.null(.dx)) {
  cat("\nsg.harm :", if (is.null(.dx$sg.harm)) "NULL (no subgroup)"
                     else paste(.dx$sg.harm, collapse = " & "), "\n")
  cat("gate    :", if (is.null(.dx$debias_gate)) "NULL" else "present", "\n")
}
