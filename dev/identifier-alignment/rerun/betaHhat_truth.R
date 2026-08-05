# betaHhat_truth.R
# ---------------------------------------------------------------------------
# Super-population conditional targets  beta(Hhat)  and  beta(Hhat^c).
#
# beta(Hhat) is the POPULATION value of the STANDARD within-subgroup analysis --
#   coxph(Surv(time, event) ~ treat)  on observed (randomized, censored) data --
# evaluated at the REALIZED selection region Hhat.  It is the estimand Theorem 2
# certifies coverage for; it is NOT the true-subgroup potential-outcome target
# theta-dagger(H).  The complement beta(Hhat^c) is the same quantity on Hhat^c.
#
# These travel WITH the simulation bundle: the PRODUCER (each simulation engine)
# computes them once per distinct realized rule on a large, fixed evaluation
# population; the manuscript fragments are pure CONSUMERS that only read them.
# Sourced verbatim by BOTH engines so the scored target is identical across the
# MR grid (methods_estimation_sweep) and the FB batches (fs_t1_t2).
#
# Requires (from forestsearch, already loaded by the engines):
#   simulate_from_dgm(), get_dfpred(), survival::coxph().
# ---------------------------------------------------------------------------

# Evaluation population for beta(Hhat): the ENTIRE fixed super-population
# dgm$df_super -- every subject exactly once (replace = FALSE, n = nrow) -- under
# the SAME randomized, censored analysis the trials run (matching analysis_time /
# cens_adjust).  This makes beta(Hhat) literally "run the standard within-subgroup
# analysis on ALL Hhat subjects of the universe"; because the simulated trials are
# sub-samples of this same df_super, the resampling population and the beta(Hhat)
# evaluation population are ONE fixed universe.  beta(Hhat) is n-invariant (a
# property of the realized rule and the population, not the trial size).  The fixed
# eval_seed only fixes the single outcome/censoring realization on the full pool
# (reproducible; residual Monte-Carlo noise is negligible at this size) and sits
# far from the trial seed band seed_base + sim_id.
#
# NOTE: this file is sourced verbatim by BOTH engines, so build_eval_frame() no
# longer takes a draw size -- it always evaluates the full pool.  A caller still
# passing the legacy `n_eval` triggers a hard error (below) rather than silently
# scoring a different target, so the MR grid and the FB batches stay identical.
build_eval_frame <- function(dgm, analysis_time = 84,
                             cens_adjust = log(1.5), eval_seed = 20260628L,
                             n_eval = NULL) {
  if (!is.null(n_eval))
    stop("build_eval_frame() now evaluates the FULL super-population ",
         "(each subject once); the legacy `n_eval` argument has been removed. ",
         "Drop n_eval from the call so both engines score the identical ",
         "full-pool beta(Hhat).")
  simulate_from_dgm(dgm, n = nrow(dgm$df_super), replace = FALSE,
                    analysis_time = analysis_time,
                    cens_adjust = cens_adjust, seed = eval_seed)
}

# Population HR of treat within an index set of the evaluation frame.  Mirrors the
# engines' .cox_hr_ci() point estimate (the coef is robust-invariant), guarding
# separation / too-few-events with NA rather than an Inf.
.beta_region <- function(eval_df, idx, outcome.name, event.name, treat.name) {
  d <- eval_df[idx, , drop = FALSE]
  if (sum(d[[event.name]], na.rm = TRUE) < 5L ||
      length(unique(d[[treat.name]])) < 2L) return(NA_real_)
  fml <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s",
                                   outcome.name, event.name, treat.name))
  fit <- tryCatch(survival::coxph(fml, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  b <- unname(stats::coef(fit)[1L])
  if (is.finite(b)) exp(b) else NA_real_
}

# beta(Hhat) and beta(Hhat^c) for ONE rule.  `rule` may be the named sg.harm
# character vector OR the " & "-joined sg_def string (the inverse of
# paste(sg.harm, collapse = " & ")); GRF "|" disjunctions live inside a single
# element and so survive the split.  Membership convention matches get_dfpred():
# treat.recommend == 0L  <=>  in Hhat (the harm region).
betaHhat_one <- function(rule, eval_df,
                         outcome.name = "y_sim", event.name = "event_sim",
                         treat.name = "treat_sim") {
  sg <- if (length(rule) == 1L && grepl(" & ", rule, fixed = TRUE))
          strsplit(rule, " & ", fixed = TRUE)[[1L]] else rule
  pred <- tryCatch(get_dfpred(eval_df, sg), error = function(e) NULL)
  if (is.null(pred) || is.null(pred$treat.recommend))
    return(c(betaHhat_H = NA_real_, betaHhat_Hc = NA_real_,
             nH_eval = NA_real_, nHc_eval = NA_real_))
  inH <- pred$treat.recommend == 0L
  c(betaHhat_H  = .beta_region(eval_df,  inH, outcome.name, event.name, treat.name),
    betaHhat_Hc = .beta_region(eval_df, !inH, outcome.name, event.name, treat.name),
    nH_eval = sum(inH), nHc_eval = sum(!inH))
}

# Deduplicated beta(Hhat)/beta(Hhat^c) over the DISTINCT rules in `sg_defs`.
# Runs in the MAIN process after the replicate loop (NOT inside the parallel
# workers -- the big eval frame stays in one place), computing each distinct rule
# once.  Returns a data.frame keyed by sg_def for a match()-merge back onto the
# per-replicate records.
betaHhat_table <- function(sg_defs, eval_df,
                           outcome.name = "y_sim", event.name = "event_sim",
                           treat.name = "treat_sim") {
  u <- unique(sg_defs[!is.na(sg_defs) & nzchar(sg_defs)])
  empty <- data.frame(sg_def = character(0), betaHhat_H = numeric(0),
                      betaHhat_Hc = numeric(0), nH_eval = integer(0),
                      nHc_eval = integer(0), stringsAsFactors = FALSE)
  if (!length(u)) return(empty)
  do.call(rbind, lapply(u, function(g) {
    v <- betaHhat_one(g, eval_df, outcome.name, event.name, treat.name)
    data.frame(sg_def      = g,
               betaHhat_H  = unname(v["betaHhat_H"]),
               betaHhat_Hc = unname(v["betaHhat_Hc"]),
               nH_eval     = as.integer(unname(v["nH_eval"])),
               nHc_eval    = as.integer(unname(v["nHc_eval"])),
               stringsAsFactors = FALSE)
  }))
}

# Attach betaHhat_H / betaHhat_Hc to a results data.frame carrying an sg_def
# column (NA for undetected reps -> NA targets).  This is the one call each engine
# adds to run_cell(), right before assembling the bundle.
attach_betaHhat <- function(results, eval_df,
                            outcome.name = "y_sim", event.name = "event_sim",
                            treat.name = "treat_sim") {
  if (is.null(results$sg_def)) {
    results$betaHhat_H <- NA_real_; results$betaHhat_Hc <- NA_real_
    return(results)
  }
  bt <- betaHhat_table(results$sg_def, eval_df, outcome.name, event.name, treat.name)
  j  <- match(results$sg_def, bt$sg_def)
  results$betaHhat_H  <- bt$betaHhat_H[j]
  results$betaHhat_Hc <- bt$betaHhat_Hc[j]
  results
}

# Optional sanity gate: theta-dagger on the SAME eval frame at the TRUE harm flag
# (flag_harm == 1 -> H).  Should reproduce dgm$hr_H_true / dgm$hr_Hc_true to MC
# error, confirming the eval frame is large enough to read as "population" and
# that beta(Hhat) and theta-dagger sit on one common population realization.
betaHhat_theta_dagger_check <- function(eval_df, harm.name = "flag_harm",
                                        outcome.name = "y_sim",
                                        event.name = "event_sim",
                                        treat.name = "treat_sim") {
  inH <- eval_df[[harm.name]] == 1L
  c(thetaDagger_H  = .beta_region(eval_df,  inH, outcome.name, event.name, treat.name),
    thetaDagger_Hc = .beta_region(eval_df, !inH, outcome.name, event.name, treat.name))
}

# ---------------------------------------------------------------------------
# Identification-structure enrichment: covariate NAMES of a realized rule.
#
# fs_identification_figures.qmd drives its involvement / structure / accuracy
# summaries off a "+"-joined covariate-name string `covs` (e.g. "er+meno"), which
# it splits at strsplit(covs, "+").  rule_covs() projects that out of the SAME
# fs.est$sg.harm used for sg_def -- sg_def is the full rule ("er <= 8 & meno == 0"),
# covs is just its variable names ("er+meno") -- so the single t1_t2 bundle drives
# the identification figures with no separate id-sweep.  Stored at record time:
#
#   rec$covs <- paste(rule_covs(fs.est$sg.harm), collapse = "+")
#
# Parsing mirrors get_dfpred(): strip "!"/braces, split conjunctions ("&") and GRF
# disjunctions ("|"), then take the leading variable token of each comparison.
rule_covs <- function(sg.harm) {
  if (is.null(sg.harm) || !length(sg.harm)) return(character(0))
  comps <- unlist(strsplit(unlist(strsplit(sg.harm, "\\s*\\|\\s*")), "\\s*&\\s*"))
  comps <- gsub("^\\s*\\(|\\)\\s*$", "", comps)            # strip parens
  comps <- gsub("^!?\\{(.*)\\}$", "\\1", trimws(comps))    # strip negation + braces
  v <- sub("^\\s*([A-Za-z.][A-Za-z0-9._]*).*$", "\\1", comps)   # leading variable token
  unique(v[nzchar(v)])
}
