# betaHhat_truth_glm.R -- SHIM
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
# (all five binary rules bitwise identical, diff == 0).
# ---------------------------------------------------------------------------

stopifnot(requireNamespace("forestsearch", quietly = TRUE))


# Evaluation population for beta(Hhat): the ENTIRE fixed super-population
# dgm$df_super -- every subject exactly once (replace = FALSE, n = nrow) --
# under the SAME randomized analysis the trials run.  beta(Hhat) is
# n-invariant.  The fixed eval_seed pins one treatment/outcome realization on
# the full pool and sits far from the trial seed band seed_base + sim_id.
build_eval_frame_glm <- function(dgm, eval_seed = 20260628L) {
  forestsearch::simulate_from_glm_dgm(dgm, n = nrow(dgm$df_super),
                                      replace = FALSE, seed = eval_seed)
}


# --- delegating wrappers ----------------------------------------------------

.beta_region_or <- function(eval_df, idx, outcome.name, treat.name) {
  forestsearch:::.fs_region_effect(eval_df, idx,
                                   outcome_type = "binary",
                                   effect_measure = NULL,
                                   outcome.name = outcome.name,
                                   treat.name = treat.name)
}

# Membership convention matches get_dfpred(): treat.recommend == 0L <=> in
# Hhat.  The binary pathway has never taken a `focus` argument, so this pins
# focus = "harm", which is what the hard-coded comparison meant.
betaHhat_one_or <- function(rule, eval_df,
                            outcome.name = "y_sim", treat.name = "treat_sim") {
  r <- forestsearch::fs_betaHhat_one(
    rule, eval_df, focus = "harm", outcome_type = "binary",
    outcome.name = outcome.name, treat.name = treat.name)
  c(betaHhat_H  = r$betaHhat_H,  betaHhat_Hc = r$betaHhat_Hc,
    nH_eval     = r$nH_eval,     nHc_eval    = r$nHc_eval)
}

betaHhat_table_or <- function(sg_defs, eval_df,
                              outcome.name = "y_sim", treat.name = "treat_sim") {
  tb <- forestsearch::fs_betaHhat_table(
    sg_defs, eval_df, focus = "harm", outcome_type = "binary",
    outcome.name = outcome.name, treat.name = treat.name)
  data.frame(sg_def = tb$sg_def, betaHhat_H = tb$betaHhat_H,
             betaHhat_Hc = tb$betaHhat_Hc,
             nH_eval = as.integer(tb$nH_eval),
             nHc_eval = as.integer(tb$nHc_eval),
             stringsAsFactors = FALSE)
}

attach_betaHhat_or <- function(results, eval_df,
                               outcome.name = "y_sim", treat.name = "treat_sim") {
  out <- forestsearch::fs_attach_betaHhat(
    results, eval_df, focus = "harm", outcome_type = "binary",
    outcome.name = outcome.name, treat.name = treat.name)
  out$betaHhat_status <- NULL
  out$nH_eval <- NULL
  out$nHc_eval <- NULL
  out
}

betaHhat_theta_dagger_check_or <- function(eval_df, harm.name = "flag_harm",
                                           outcome.name = "y_sim",
                                           treat.name = "treat_sim") {
  inH <- eval_df[[harm.name]] == 1L
  c(thetaDagger_H  = .beta_region_or(eval_df,  inH, outcome.name, treat.name),
    thetaDagger_Hc = .beta_region_or(eval_df, !inH, outcome.name, treat.name))
}

betaHhat_neff_parity_or <- function(cov_df, strict = TRUE,
                                    target_beta = "C_betaHhat",
                                    target_ref = "C_dagger") {
  forestsearch::fs_betaHhat_neff_parity(cov_df, strict = strict,
                                        target_beta = target_beta,
                                        target_ref = target_ref)
}

betaHhat_neff_report_or <- function(bundle, block = c("H", "Hc")) {
  forestsearch::fs_betaHhat_neff_report(bundle, block = block)
}
