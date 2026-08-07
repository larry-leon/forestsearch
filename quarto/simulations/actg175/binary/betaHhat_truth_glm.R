# betaHhat_truth_glm.R
# ---------------------------------------------------------------------------
# Super-population conditional targets  beta(Hhat)  and  beta(Hhat^c) on the
# BINARY / ODDS-RATIO scale.  This is the GLM analogue of betaHhat_truth.R
# (which is on the Cox / hazard-ratio scale).
#
# beta(Hhat) is the POPULATION value of the STANDARD within-subgroup analysis --
#   glm(y ~ treat, family = binomial("logit"))  (the marginal subgroup odds
# ratio) on the observed (randomized) data -- evaluated at the REALIZED
# selection region Hhat.  It is the estimand Theorem 2 certifies coverage for;
# it is NOT the marginal true-subgroup target theta-dagger(H) (truth$marg_H).
# The complement beta(Hhat^c) is the same quantity on Hhat^c.
#
# These travel WITH the simulation bundle: the PRODUCER (the sweep engine)
# computes them once per distinct realized rule on the fixed evaluation
# population; the manuscript fragment is a pure CONSUMER that only reads them.
#
# WHY THE _or / _glm SUFFIXES.  Every exported name here carries an _or (or
# _glm) suffix so this module can never be confused with -- or silently mask --
# the survival betaHhat_truth.R, whose .beta_region() fits a coxph().  The two
# files are never meant to be sourced into a single session; the suffix makes a
# mistaken co-source fail loudly (object-not-found) rather than silently score a
# Cox target with the GLM code (or vice versa).  The COLUMN names it produces
# (betaHhat_H / betaHhat_Hc) are kept scale-agnostic and identical to the
# survival module so the downstream scoring code (paste0("betaHhat_", suffix))
# is shared verbatim.
#
# Requires (from forestsearch, already loaded by the engine):
#   simulate_from_glm_dgm()  [with the `replace` argument], get_dfpred(),
#   stats::glm().
# ---------------------------------------------------------------------------

# Evaluation population for beta(Hhat): the ENTIRE fixed super-population
# dgm$df_super -- every subject exactly once (replace = FALSE, n = nrow) -- under
# the SAME randomized analysis the trials run.  Because the simulated trials are
# sub-samples of this same df_super, the resampling population and the beta(Hhat)
# evaluation population are ONE fixed universe.  beta(Hhat) is n-invariant (a
# property of the realized rule and the population, not the trial size).  The
# fixed eval_seed only fixes the single treatment/outcome realization on the full
# pool (reproducible; residual Monte-Carlo noise is negligible at this size) and
# sits far from the trial seed band seed_base + sim_id.
build_eval_frame_glm <- function(dgm, eval_seed = 20260628L) {
  simulate_from_glm_dgm(dgm, n = nrow(dgm$df_super), replace = FALSE,
                        seed = eval_seed)
}

# Population OR of treat within an index set of the evaluation frame.  Mirrors the
# engine's .logit_or_ci() point estimate (glm logit, exp(coef)), guarding
# separation / a near-constant outcome / too-few subjects in either outcome cell
# with NA rather than an Inf.  (Binary: the outcome IS the event, so there is no
# separate event column.)
.beta_region_or <- function(eval_df, idx, outcome.name, treat.name) {
  d <- eval_df[idx, , drop = FALSE]
  y <- d[[outcome.name]]
  if (length(unique(d[[treat.name]])) < 2L ||
      sum(y == 1L, na.rm = TRUE) < 5L || sum(y == 0L, na.rm = TRUE) < 5L)
    return(NA_real_)
  fml <- stats::as.formula(sprintf("%s ~ %s", outcome.name, treat.name))
  fit <- tryCatch(suppressWarnings(
           stats::glm(fml, data = d, family = stats::binomial("logit"))),
         error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  cf <- stats::coef(fit)
  if (!(treat.name %in% names(cf))) return(NA_real_)
  b <- unname(cf[[treat.name]])
  if (is.finite(b)) exp(b) else NA_real_
}

# beta(Hhat) and beta(Hhat^c) for ONE rule.  `rule` may be the named sg.harm
# character vector OR the " & "-joined sg_def string (the inverse of
# paste(sg.harm, collapse = " & ")).  Membership convention matches
# get_dfpred(): treat.recommend == 0L  <=>  in Hhat (the harm region).
betaHhat_one_or <- function(rule, eval_df,
                            outcome.name = "y_sim", treat.name = "treat_sim") {
  # Dispatch order matters and mirrors get_dfpred() itself
  # (R/forestsearch_helpers.R:101): the DISJUNCTION form is tested FIRST, on
  # the UNSPLIT string.  A GRF disjunction contains " & " inside each
  # conjunct, so splitting first shreds it --
  #   "(age > 34 & preanti <= 744.5) | (wtkg <= 60)"
  #     -> c("(age > 34", "preanti <= 744.5) | (wtkg <= 60)")
  # which get_dfpred() then fails to resolve, yielding NA targets rather than
  # a loud error.  An earlier version of this comment claimed disjunctions
  # "live inside a single element and so survive the split"; they do not.
  # The defect was latent, not live: .grf_build_subgroup_definition() emits the
  # " | " form only for a multi-leaf policy tree, and every GRF selection in
  # the committed mr_sweep bundles collapsed to a single leaf (0 of 250,025
  # non-empty sg_def rows under quarto/ contain "|").  Multi-leaf selection is
  # reachable at maxdepth >= 2, so this is a correction, not a repair.
  sg <- if (length(rule) == 1L && grepl("|", rule, fixed = TRUE)) {
          rule
        } else if (length(rule) == 1L && grepl(" & ", rule, fixed = TRUE)) {
          strsplit(rule, " & ", fixed = TRUE)[[1L]]
        } else {
          rule
        }
  pred <- tryCatch(get_dfpred(eval_df, sg), error = function(e) NULL)
  if (is.null(pred) || is.null(pred$treat.recommend))
    return(c(betaHhat_H = NA_real_, betaHhat_Hc = NA_real_,
             nH_eval = NA_real_, nHc_eval = NA_real_))
  inH <- pred$treat.recommend == 0L
  c(betaHhat_H  = .beta_region_or(eval_df,  inH, outcome.name, treat.name),
    betaHhat_Hc = .beta_region_or(eval_df, !inH, outcome.name, treat.name),
    nH_eval = sum(inH), nHc_eval = sum(!inH))
}

# Deduplicated beta(Hhat)/beta(Hhat^c) over the DISTINCT rules in `sg_defs`.
# Runs in the MAIN process after the replicate loop (NOT inside the parallel
# workers -- the big eval frame stays in one place), computing each distinct rule
# once.  Returns a data.frame keyed by sg_def for a match()-merge back onto the
# per-replicate records.
betaHhat_table_or <- function(sg_defs, eval_df,
                              outcome.name = "y_sim", treat.name = "treat_sim") {
  u <- unique(sg_defs[!is.na(sg_defs) & nzchar(sg_defs)])
  empty <- data.frame(sg_def = character(0), betaHhat_H = numeric(0),
                      betaHhat_Hc = numeric(0), nH_eval = integer(0),
                      nHc_eval = integer(0), stringsAsFactors = FALSE)
  if (!length(u)) return(empty)
  do.call(rbind, lapply(u, function(g) {
    v <- betaHhat_one_or(g, eval_df, outcome.name, treat.name)
    data.frame(sg_def      = g,
               betaHhat_H  = unname(v["betaHhat_H"]),
               betaHhat_Hc = unname(v["betaHhat_Hc"]),
               nH_eval     = as.integer(unname(v["nH_eval"])),
               nHc_eval    = as.integer(unname(v["nHc_eval"])),
               stringsAsFactors = FALSE)
  }))
}

# Attach betaHhat_H / betaHhat_Hc to a results data.frame carrying an sg_def
# column (NA for undetected reps -> NA targets).  This is the one call the engine
# adds to run_cell(), right before assembling the bundle.
attach_betaHhat_or <- function(results, eval_df,
                               outcome.name = "y_sim", treat.name = "treat_sim") {
  if (is.null(results$sg_def)) {
    results$betaHhat_H <- NA_real_; results$betaHhat_Hc <- NA_real_
    return(results)
  }
  bt <- betaHhat_table_or(results$sg_def, eval_df, outcome.name, treat.name)
  j  <- match(results$sg_def, bt$sg_def)
  results$betaHhat_H  <- bt$betaHhat_H[j]
  results$betaHhat_Hc <- bt$betaHhat_Hc[j]
  results
}

# Optional sanity gate: theta-dagger on the SAME eval frame at the TRUE harm flag
# (flag_harm == 1 -> H).  Should reproduce the DGM marginal ORs (truth$marg_H /
# truth$marg_Hc) to MC error, confirming the eval frame is large enough to read as
# "population" and that beta(Hhat) and theta-dagger sit on one common population
# realization.
betaHhat_theta_dagger_check_or <- function(eval_df, harm.name = "flag_harm",
                                           outcome.name = "y_sim",
                                           treat.name = "treat_sim") {
  inH <- eval_df[[harm.name]] == 1L
  c(thetaDagger_H  = .beta_region_or(eval_df,  inH, outcome.name, treat.name),
    thetaDagger_Hc = .beta_region_or(eval_df, !inH, outcome.name, treat.name))
}


# ---------------------------------------------------------------------------
# n_eff PARITY GUARD
#
# betaHhat_H/_Hc enter the coverage tables as the C_betaHhat target.
# .coverage_meta() computes  ok <- is.finite(target) & is.finite(lo) &
# is.finite(hi);  n_eff <- sum(ok).  A non-finite target is therefore dropped
# SILENTLY: no error, no warning.  Coverage is then computed over the subset of
# replicates whose targets happened to score -- a non-random subset, since
# whatever made the target NA is a property of the realized rule.
#
# C_dagger's target (truth$marg_H) is a finite scalar, so within any
# (src, block, estimator) cell n_eff(C_dagger) counts exactly the replicates
# with a finite interval, and n_eff(C_betaHhat) <= that.  Equality holds if and
# only if no betaHhat target was dropped.  Parity is thus the exact check, and
# it is free -- n_eff is already in the table.
#
# TWO causes break parity, and the guard deliberately does not distinguish
# them, because both make C_betaHhat and C_dagger incomparable as printed:
#   (1) a rule get_dfpred() could not resolve -> NA.  This is the class the
#       disjunction defect above belonged to.
#   (2) a genuinely degenerate region -> .beta_region_or() returns NA BY
#       DESIGN (one treatment arm, or <5 in an outcome cell).
# Cause (2) is legitimate behaviour, so a mismatch is not automatically a bug
# -- but it is always a reason the printed coverage is computed on a different
# denominator than the row beside it.  Use betaHhat_neff_report_or() on the
# bundle to tell the two apart before deciding.
#
# `strict = TRUE` stops; set FALSE to warn and carry on once the cause is known.
betaHhat_neff_parity_or <- function(cov_df, strict = TRUE,
                                    target_beta = "C_betaHhat",
                                    target_ref  = "C_dagger") {
  need <- c("target", "n_eff", "block", "estimator")
  if (!all(need %in% names(cov_df))) return(invisible(cov_df))
  keycols <- intersect(c("src", "subgroup_method", "n_sample", "block", "estimator"),
                       names(cov_df))
  key <- function(d) do.call(paste, c(d[keycols], sep = "|"))

  ref  <- cov_df[cov_df$target == target_ref,  , drop = FALSE]
  beta <- cov_df[cov_df$target == target_beta, , drop = FALSE]
  if (!nrow(ref) || !nrow(beta)) return(invisible(cov_df))

  ref$.k <- key(ref); beta$.k <- key(beta)
  j <- match(ref$.k, beta$.k)
  cmp <- data.frame(ref[keycols],
                    n_eff_ref  = ref$n_eff,
                    n_eff_beta = ifelse(is.na(j), NA_integer_, beta$n_eff[j]),
                    stringsAsFactors = FALSE)
  cmp$dropped <- cmp$n_eff_ref - cmp$n_eff_beta
  bad <- cmp[is.na(cmp$n_eff_beta) | cmp$dropped != 0L, , drop = FALSE]

  if (nrow(bad)) {
    msg <- paste0(
      sprintf("n_eff parity FAILED in %d of %d cells: %s is computed on a ",
              nrow(bad), nrow(cmp), target_beta),
      sprintf("smaller denominator than %s.\n", target_ref),
      "  Its coverage is therefore NOT comparable to the row beside it.\n",
      paste(utils::capture.output(print(bad, row.names = FALSE)), collapse = "\n"),
      "\n  Run betaHhat_neff_report_or(bundle) to see which rules produced the ",
      "NA targets\n  and whether they are unresolvable rules or degenerate regions.")
    if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(cov_df)
}


# Diagnostic companion: for ONE bundle, list the realized rules whose betaHhat
# target is NA, and flag whether each is a disjunction (cause 1 above) or a
# region small enough for .beta_region_or()'s guards to fire (cause 2).
betaHhat_neff_report_or <- function(bundle, block = c("H", "Hc")) {
  if (is.character(bundle)) bundle <- readRDS(bundle)
  r <- if (is.data.frame(bundle)) bundle else bundle$results
  block <- match.arg(block)
  col <- paste0("betaHhat_", block)
  if (is.null(r[[col]]) || is.null(r$sg_def)) return(invisible(NULL))
  bad <- r[is.na(r[[col]]) & !is.na(r$sg_def) & nzchar(r$sg_def), , drop = FALSE]
  if (!nrow(bad)) { cat(sprintf("no NA %s among rows with a realized rule\n", col))
                    return(invisible(NULL)) }
  out <- data.frame(
    sg_def         = bad$sg_def,
    is_disjunction = grepl("|", bad$sg_def, fixed = TRUE),
    nH_eval        = if (!is.null(bad$nH_eval)) bad$nH_eval else NA_integer_,
    stringsAsFactors = FALSE)
  out <- out[!duplicated(out$sg_def), , drop = FALSE]
  cat(sprintf("%d distinct rule(s) with NA %s; %d are disjunctions\n",
              nrow(out), col, sum(out$is_disjunction)))
  out
}
