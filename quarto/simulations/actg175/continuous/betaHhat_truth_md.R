# betaHhat_truth_md.R
# ---------------------------------------------------------------------------
# Super-population conditional targets  beta(Hhat)  and  beta(Hhat^c) on the
# CONTINUOUS / MEAN-DIFFERENCE scale.  This is the MD analogue of
# betaHhat_truth_glm.R (odds ratio) and betaHhat_truth.R (hazard ratio).
#
# beta(Hhat) is the POPULATION value of the STANDARD within-subgroup analysis --
# the marginal subgroup mean difference on the observed (randomized) data --
# evaluated at the REALIZED selection region Hhat.  It is the estimand
# Theorem 2 certifies coverage for; it is NOT the marginal true-subgroup target
# theta-dagger(Q).  The complement beta(Hhat^c) is the same quantity on Hhat^c.
#
# WHY THIS MODULE HAS NO EVALUATION FRAME.  The hazard ratio and the odds ratio
# are non-collapsible, so their modules must fit a model on a realized
# evaluation frame.  The mean difference is COLLAPSIBLE: under randomization
# (A independent of X), for ANY region R,
#
#     beta(R) = E[Y | A=1, R] - E[Y | A=0, R] = E[tau(X) | R],
#
# because the prognostic part E[m(X) | R] is common to both arms and cancels.
# No independence between covariates, no rectangle shape, no step-shaped CATE
# is required.  df_super is a fixed finite population, so beta(Hhat) is an
# EXACT finite mean over it: no evaluation frame, no eval_seed, no model fit,
# and no residual Monte-Carlo noise -- not "negligible", zero.  Repeated
# evaluation is bit-identical.  See
# dev/glm-continuous-sims/NOTE_target_is_collapsibility.md, section 2.
#
# WHAT THE SUBJECT-LEVEL EFFECT LOOKS LIKE HERE.  Read from
# R/generate_glm_dgm.R:414-427 (and confirmed against the additive base
# formula at R/generate_glm_dgm.R:310-323, which carries no
# treatment-by-covariate interactions), the continuous DGM is
#
#     mu1_i = mu0_i + k_treat * beta_treat + beta_inter * 1{i in Q},
#
# so the subject-level effect is TWO-VALUED: delta := k_treat * beta_treat
# outside Q, and delta + beta_inter inside.  Hence
#
#     beta(Hhat) = delta + beta_inter * P(Q | Hhat).
#
# The originating handoff's tau * P(Q | Hhat) sets the complement effect to
# zero; it is not zero, and that form is low by delta on EVERY region --
# including a region disjoint from Q, where it gives 0 against a truth of
# delta.  See NOTE_target_is_collapsibility.md, section 3.  This module never
# forms either expression: it averages the potential outcomes directly, which
# is correct whatever the complement effect happens to be.
#
# WHY THE _md SUFFIX.  Every exported name here carries an _md suffix so this
# module can never be confused with -- or silently mask -- the survival
# betaHhat_truth.R (whose .beta_region() fits a coxph()) or the binary
# betaHhat_truth_glm.R (whose .beta_region_or() fits a glm()).  The three files
# are never meant to be sourced into a single session; the suffix makes a
# mistaken co-source fail loudly (object-not-found) rather than silently score
# an MD target with Cox or GLM code.  The COLUMN names it produces
# (betaHhat_H / betaHhat_Hc / nH_eval / nHc_eval) are kept scale-agnostic and
# identical to the other two modules so the downstream scoring code
# (paste0("betaHhat_", suffix)) is shared verbatim.
#
# Requires (from forestsearch, already loaded by the engine):
#   compute_aor(), get_dfpred().
# Note there is NO simulate_from_glm_dgm() dependency: the target is read off
# dgm$df_super directly, never off a simulated trial.
# ---------------------------------------------------------------------------


# --- focus ------------------------------------------------------------------
# The harm/benefit fork is OPEN for the continuous pathway and this module does
# not settle it.  betaHhat_one_or() reads treat.recommend == 0L as in-region,
# which is the HARM convention.  The continuous template searches for BENEFIT
# via treatment switching, so its realized regions are Ghat/Ghat^c and the
# sense inverts.  `focus` is therefore a required argument with NO default: a
# missing focus is an error, not a guess.  Making it explicit keeps the choice
# visible at the call site and in the bundle.
.check_focus_md <- function(focus) {
  if (missing(focus) || is.null(focus)) {
    stop("`focus` is required and has no default: supply \"harm\" or ",
         "\"benefit\".  The harm/benefit convention for the continuous ",
         "pathway is an open decision; it must be stated at the call site.",
         call. = FALSE)
  }
  if (!(is.character(focus) && length(focus) == 1L &&
        focus %in% c("harm", "benefit"))) {
    stop("`focus` must be exactly one of \"harm\" or \"benefit\"; got: ",
         paste(utils::capture.output(utils::str(focus)), collapse = " "),
         call. = FALSE)
  }
  focus
}


# --- the region computation -------------------------------------------------
# Population MD within an index set of df_super.  DELEGATED to the package:
# compute_aor() is exported, accepts an arbitrary subset_indicator, and
# dispatches "MD" to mean(mu1) - mean(mu0).
#
# Do NOT restate mean(mu1) - mean(mu0) inline here.  Under the code-theory
# alignment standard a quantity reconstructed at a second site is a defect, and
# the package already carries two sites for this arithmetic (compute_aor() and
# the local .effect_one() closure inside .compute_glm_effects()).  This module
# adds no third site and refactors neither of the existing two.
#
# The only degenerate case is an empty region -- there is no model fit, so
# separation and minimum-cell guards have no analogue.
.beta_region_md <- function(df_super, idx) {
  if (!length(idx) || !any(idx)) {
    return(NA_real_)
  }
  compute_aor(df_super,
              subset_indicator = as.integer(idx),
              effect_measure   = "MD")
}


# --- one rule ---------------------------------------------------------------
# beta(Hhat) and beta(Hhat^c) for ONE rule.  `rule` may be the named sg.harm
# character vector OR the " & "-joined sg_def string (the inverse of
# paste(sg.harm, collapse = " & ")); GRF "|" disjunctions live inside a single
# element and so survive the split.  Membership is resolved through
# get_dfpred(), exactly as betaHhat_one_or() does; only the sense of
# treat.recommend depends on `focus`.
betaHhat_one_md <- function(rule, df_super, focus) {
  .check_focus_md(focus)
  # Dispatch order matters and mirrors get_dfpred() itself
  # (R/forestsearch_helpers.R:101): the DISJUNCTION form is tested FIRST, on
  # the UNSPLIT string.  A GRF disjunction contains " & " inside each
  # conjunct, so splitting first shreds it --
  #   "(age > 34 & preanti <= 744.5) | (wtkg <= 60)"
  #     -> c("(age > 34", "preanti <= 744.5) | (wtkg <= 60)")
  # which get_dfpred() then fails to resolve, yielding NA targets rather than
  # a loud error.  The binary sibling betaHhat_truth_glm.R:76-77 splits first
  # and carries exactly that defect; do not copy it back in.
  sg <- if (length(rule) == 1L && grepl("|", rule, fixed = TRUE)) {
          rule
        } else if (length(rule) == 1L && grepl(" & ", rule, fixed = TRUE)) {
          strsplit(rule, " & ", fixed = TRUE)[[1L]]
        } else {
          rule
        }
  pred <- tryCatch(get_dfpred(df_super, sg), error = function(e) NULL)
  if (is.null(pred) || is.null(pred$treat.recommend))
    return(c(betaHhat_H = NA_real_, betaHhat_Hc = NA_real_,
             nH_eval = NA_real_, nHc_eval = NA_real_))
  in_region <- if (focus == "harm") {
    pred$treat.recommend == 0L
  } else {
    pred$treat.recommend == 1L
  }
  c(betaHhat_H  = .beta_region_md(df_super,  in_region),
    betaHhat_Hc = .beta_region_md(df_super, !in_region),
    nH_eval = sum(in_region), nHc_eval = sum(!in_region))
}


# --- deduplicated table -----------------------------------------------------
# Deduplicated beta(Hhat)/beta(Hhat^c) over the DISTINCT rules in `sg_defs`.
# Deduplication is by distinct sg_def, NOT by replicate: two replicates that
# land on the same rule score the same target, and the target does not depend
# on the trial.  Returns a data.frame keyed by sg_def for a match()-merge back
# onto the per-replicate records.
betaHhat_table_md <- function(sg_defs, df_super, focus) {
  .check_focus_md(focus)
  u <- unique(sg_defs[!is.na(sg_defs) & nzchar(sg_defs)])
  empty <- data.frame(sg_def = character(0), betaHhat_H = numeric(0),
                      betaHhat_Hc = numeric(0), nH_eval = integer(0),
                      nHc_eval = integer(0), stringsAsFactors = FALSE)
  if (!length(u)) return(empty)
  do.call(rbind, lapply(u, function(g) {
    v <- betaHhat_one_md(g, df_super, focus)
    data.frame(sg_def      = g,
               betaHhat_H  = unname(v["betaHhat_H"]),
               betaHhat_Hc = unname(v["betaHhat_Hc"]),
               nH_eval     = as.integer(unname(v["nH_eval"])),
               nHc_eval    = as.integer(unname(v["nHc_eval"])),
               stringsAsFactors = FALSE)
  }))
}


# --- the one call each engine adds ------------------------------------------
# Attach betaHhat_H / betaHhat_Hc to a results data.frame carrying an sg_def
# column (NA for undetected reps -> NA targets).  Unchanged from the engines'
# point of view: this remains the one call each engine adds to run_cell(),
# right before assembling the bundle.
attach_betaHhat_md <- function(results, df_super, focus) {
  .check_focus_md(focus)
  if (is.null(results$sg_def)) {
    results$betaHhat_H <- NA_real_; results$betaHhat_Hc <- NA_real_
    return(results)
  }
  bt <- betaHhat_table_md(results$sg_def, df_super, focus)
  j  <- match(results$sg_def, bt$sg_def)
  results$betaHhat_H  <- bt$betaHhat_H[j]
  results$betaHhat_Hc <- bt$betaHhat_Hc[j]
  results
}


# --- theta-dagger identity gate ---------------------------------------------
# theta-dagger on the SAME population at the TRUE subgroup flag
# (flag_harm == 1 -> Q).  Because the target is an exact finite mean over
# df_super and the DGM's own effect_Q / effect_Qc are computed from the same
# mu0/mu1 columns by the same compute_aor() dispatch, this is an EXACT
# IDENTITY, not agreement to Monte-Carlo error.  Any nonzero difference is a
# real defect; do not meet it with a tolerance.
#
# No `focus` argument: flag_harm == 1 IS Q by construction, whichever way the
# realized regions are read.
betaHhat_theta_dagger_check_md <- function(df_super, harm.name = "flag_harm") {
  if (is.null(df_super[[harm.name]]))
    stop("`df_super` has no column \"", harm.name, "\".", call. = FALSE)
  inQ <- df_super[[harm.name]] == 1L
  c(thetaDagger_H  = .beta_region_md(df_super,  inQ),
    thetaDagger_Hc = .beta_region_md(df_super, !inQ))
}
