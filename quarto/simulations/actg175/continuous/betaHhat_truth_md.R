# betaHhat_truth_md.R -- SHIM
# ---------------------------------------------------------------------------
# This module is now a thin wrapper over the package implementation in
# R/betaHhat_truth.R.  It exists only so that sweep .qmd files sourcing it keep
# working unchanged; every name and signature below is exactly what it was
# before consolidation.
#
# THE LOGIC LIVES IN THE PACKAGE.  Do not add arithmetic here.  Four
# near-identical copies of this file were the reason four defects (D1-D4) each
# had to be found and fixed up to four times; re-introducing logic in a shim
# recreates that.  See
# dev/betaHhat-consolidation/SPEC_betaHhat_package_function.md and
# dev/betaHhat-consolidation/T10_GATE_RESULT.md.
#
# Migration status: step 3 of 4.  Step 4 re-points the sweeps at the package
# functions directly and deletes this file.
#
# T10 gate: every value below reproduces the pre-consolidation module exactly
# (all six continuous rules bitwise identical, diff == 0).
# ---------------------------------------------------------------------------

stopifnot(requireNamespace("forestsearch", quietly = TRUE))


# Evaluation population for beta(Hhat): the fixed super-population itself.  The
# mean difference is collapsible, so beta(Hhat) is an exact finite mean over
# df_super -- no evaluation frame, no eval_seed, no model fit, and no residual
# Monte-Carlo noise.  There is deliberately no build_eval_frame_md().


# --- delegating wrappers ----------------------------------------------------

.check_focus_md <- function(focus) forestsearch:::.fs_check_focus(focus)

.beta_region_md <- function(df_super, idx) {
  forestsearch:::.fs_region_effect(df_super, idx,
                                   outcome_type = "continuous",
                                   effect_measure = "MD")
}

betaHhat_one_md <- function(rule, df_super, focus) {
  .check_focus_md(focus)
  r <- forestsearch::fs_betaHhat_one(
    rule, df_super, focus = focus,
    outcome_type = "continuous", effect_measure = "MD")
  c(betaHhat_H  = r$betaHhat_H,  betaHhat_Hc = r$betaHhat_Hc,
    nH_eval     = r$nH_eval,     nHc_eval    = r$nHc_eval)
}

betaHhat_table_md <- function(sg_defs, df_super, focus) {
  .check_focus_md(focus)
  tb <- forestsearch::fs_betaHhat_table(
    sg_defs, df_super, focus = focus,
    outcome_type = "continuous", effect_measure = "MD")
  data.frame(sg_def = tb$sg_def, betaHhat_H = tb$betaHhat_H,
             betaHhat_Hc = tb$betaHhat_Hc,
             nH_eval = as.integer(tb$nH_eval),
             nHc_eval = as.integer(tb$nHc_eval),
             stringsAsFactors = FALSE)
}

attach_betaHhat_md <- function(results, df_super, focus) {
  .check_focus_md(focus)
  out <- forestsearch::fs_attach_betaHhat(
    results, df_super, focus = focus,
    outcome_type = "continuous", effect_measure = "MD")
  # Pre-consolidation shape: only the two target columns were added.
  out$betaHhat_status <- NULL
  out$nH_eval <- NULL
  out$nHc_eval <- NULL
  out
}

betaHhat_theta_dagger_check_md <- function(df_super, harm.name = "flag_harm") {
  if (is.null(df_super[[harm.name]]))
    stop("`df_super` has no column \"", harm.name, "\".", call. = FALSE)
  inQ <- df_super[[harm.name]] == 1L
  c(thetaDagger_H  = .beta_region_md(df_super,  inQ),
    thetaDagger_Hc = .beta_region_md(df_super, !inQ))
}
