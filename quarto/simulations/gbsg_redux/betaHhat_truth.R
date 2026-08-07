# betaHhat_truth.R -- SHIM
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
# functions directly and deletes this file, together with the byte-identical
# copy at dev/identifier-alignment/rerun/betaHhat_truth.R.
#
# BEHAVIOURAL CHANGE ON THIS PATHWAY -- read before comparing to old output.
# D1 (the " & " split ahead of disjunction dispatch) was still live here when
# the binary and continuous modules were fixed at e6f6024.  A GRF disjunction
# such as "(er > 125 & size > 20) | (nodes > 5)" previously shredded into
# fragments naming no real column, and the target came back NA.  It now
# resolves.  T10 records this as the single sanctioned movement:
#   old  betaHhat_H = NA               new  betaHhat_H = 0.688255552881
#   old  nH_eval    = NA               new  nH_eval    = 42793  (+ 57207 = 100000)
# All four non-disjunction survival rules are bitwise identical (diff == 0).
# Any consumer that previously received NA for a disjunction now receives a
# value.
# ---------------------------------------------------------------------------

stopifnot(requireNamespace("forestsearch", quietly = TRUE))


# Evaluation population for beta(Hhat): the ENTIRE fixed super-population,
# every subject exactly once (replace = FALSE, n = nrow), under the SAME
# randomized/censored analysis the trials run.
#
# NOTE: this file is sourced verbatim by BOTH engines, so build_eval_frame()
# does not take a draw size -- it always evaluates the full pool.  A caller
# still passing the legacy `n_eval` triggers a hard error rather than silently
# scoring a different target, so the MR grid and the FB batches stay identical.
#
# The `cens_adjust = log(1.5)` default and the `n_eval` trap are carried over
# UNCHANGED from the pre-shim module.  All 106 call sites pass cens_adjust
# explicitly and none passes n_eval, so both are inert today -- but the shim's
# contract is identical signatures, and the trap is a guard, not decoration.
build_eval_frame <- function(dgm, analysis_time = 84,
                             cens_adjust = log(1.5), eval_seed = 20260628L,
                             n_eval = NULL) {
  if (!is.null(n_eval))
    stop("build_eval_frame() now evaluates the FULL super-population ",
         "(each subject once); the legacy `n_eval` argument has been removed. ",
         "Drop n_eval from the call so both engines score the identical ",
         "full-pool beta(Hhat).")
  forestsearch::simulate_from_dgm(dgm, n = nrow(dgm$df_super), replace = FALSE,
                                  analysis_time = analysis_time,
                                  cens_adjust = cens_adjust, seed = eval_seed)
}


# --- delegating wrappers ----------------------------------------------------

.beta_region <- function(eval_df, idx, outcome.name, event.name, treat.name) {
  forestsearch:::.fs_region_effect(eval_df, idx,
                                   outcome_type = "survival",
                                   effect_measure = NULL,
                                   outcome.name = outcome.name,
                                   event.name = event.name,
                                   treat.name = treat.name)
}

# Membership convention matches get_dfpred(): treat.recommend == 0L <=> in
# Hhat.  The survival pathway has never taken a `focus` argument, so this pins
# focus = "harm", which is what the hard-coded comparison meant.
betaHhat_one <- function(rule, eval_df,
                         outcome.name = "y_sim", event.name = "event_sim",
                         treat.name = "treat_sim") {
  r <- forestsearch::fs_betaHhat_one(
    rule, eval_df, focus = "harm", outcome_type = "survival",
    outcome.name = outcome.name, event.name = event.name,
    treat.name = treat.name)
  c(betaHhat_H  = r$betaHhat_H,  betaHhat_Hc = r$betaHhat_Hc,
    nH_eval     = r$nH_eval,     nHc_eval    = r$nHc_eval)
}

betaHhat_table <- function(sg_defs, eval_df,
                           outcome.name = "y_sim", event.name = "event_sim",
                           treat.name = "treat_sim") {
  tb <- forestsearch::fs_betaHhat_table(
    sg_defs, eval_df, focus = "harm", outcome_type = "survival",
    outcome.name = outcome.name, event.name = event.name,
    treat.name = treat.name)
  data.frame(sg_def = tb$sg_def, betaHhat_H = tb$betaHhat_H,
             betaHhat_Hc = tb$betaHhat_Hc,
             nH_eval = as.integer(tb$nH_eval),
             nHc_eval = as.integer(tb$nHc_eval),
             stringsAsFactors = FALSE)
}

attach_betaHhat <- function(results, eval_df,
                            outcome.name = "y_sim", event.name = "event_sim",
                            treat.name = "treat_sim") {
  out <- forestsearch::fs_attach_betaHhat(
    results, eval_df, focus = "harm", outcome_type = "survival",
    outcome.name = outcome.name, event.name = event.name,
    treat.name = treat.name)
  out$betaHhat_status <- NULL
  out$nH_eval <- NULL
  out$nHc_eval <- NULL
  out
}

betaHhat_theta_dagger_check <- function(eval_df, harm.name = "flag_harm",
                                        outcome.name = "y_sim",
                                        event.name = "event_sim",
                                        treat.name = "treat_sim") {
  inH <- eval_df[[harm.name]] == 1L
  c(thetaDagger_H  = .beta_region(eval_df,  inH, outcome.name, event.name,
                                  treat.name),
    thetaDagger_Hc = .beta_region(eval_df, !inH, outcome.name, event.name,
                                  treat.name))
}

# Covariates referenced by a rule.  Parsing mirrors get_dfpred(): strip
# "!"/braces, split conjunctions and GRF disjunctions, take the leading
# variable token of each comparison.
rule_covs <- function(sg.harm) forestsearch:::.fs_rule_columns(sg.harm)
