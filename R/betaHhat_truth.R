# betaHhat_truth.R
# ---------------------------------------------------------------------------
# Package-level implementation of the conditional estimand beta(Hhat), shared
# by the survival, binary and continuous simulation pathways.
#
# beta(Hhat) is the POPULATION value of the standard within-subgroup analysis,
# evaluated at the REALIZED selection region Hhat.  It is the estimand
# Theorem 2 certifies coverage for; it is NOT the marginal true-subgroup target
# theta-dagger.  The complement beta(Hhat^c) is the same quantity on Hhat^c.
#
# This file consolidates logic previously duplicated across four simulation
# modules.  It changes NO estimand: every family computes exactly the target
# its module computes today.  See
# dev/betaHhat-consolidation/SPEC_betaHhat_package_function.md.
# ---------------------------------------------------------------------------


# --- rule column extraction -------------------------------------------------
# Column names referenced by a rule, used to report which ones a frame lacks.
# Handles the disjunctive form "(a & b) | (c & d)", the AND-composed brace
# vector c("{x <= 1}", "!{y > 2}"), and bare column names.
#
# @noRd
.fs_rule_columns <- function(rule) {
  if (is.null(rule) || !length(rule)) return(character(0))
  txt <- unlist(strsplit(as.character(rule), "|", fixed = TRUE))
  txt <- unlist(strsplit(txt, "&", fixed = TRUE))
  txt <- gsub("[()]", "", txt)
  txt <- gsub("^\\s*!?\\s*\\{?", "", txt)
  txt <- gsub("\\}?\\s*$", "", txt)
  txt <- trimws(txt)
  ops <- c("<=", ">=", "!=", "==", "<", ">")
  vars <- vapply(txt, function(e) {
    for (op in ops) {
      if (grepl(op, e, fixed = TRUE)) {
        return(trimws(strsplit(e, op, fixed = TRUE)[[1L]][1L]))
      }
    }
    e
  }, character(1), USE.NAMES = FALSE)
  unique(vars[nzchar(vars)])
}


# --- LAYER 1: membership resolution -----------------------------------------
# The single place a rule is interpreted.  All four historical defects lived
# here.
#
# Returns a list, never a bare vector:
#   in_region : logical(nrow(df)) under the HARM convention
#               (treat.recommend == 0L), or NULL when unresolved
#   status    : "ok" | "unresolved" | "empty"
#   missing   : character(0), or the referenced columns the frame lacks
#
# Contract:
#   1. STRUCTURED FIRST.  With `sg_def_struct`, membership comes from
#      .grf_evaluate_subgroup() -- the package's own answer for GRF membership,
#      chosen precisely because get_dfpred() cannot evaluate the disjunction
#      string (R/forestsearch_helpers.R:1772).  String parsing is the fallback.
#   2. The string path mirrors get_dfpred()'s dispatch order EXACTLY: the
#      disjunction form is tested BEFORE any " & " split, on the unsplit
#      string (R/forestsearch_helpers.R:101).  Splitting first makes
#      length > 1, skips that branch, and shreds the rule -- defect D1.
#   3. REJECT PARTIAL NA.  evaluate_comparison() warns and yields NA for a
#      missing column, and get_dfpred() propagates it asymmetrically:
#      conjunction gives FALSE & NA = FALSE but TRUE & NA = NA; disjunction
#      gives TRUE | NA = TRUE but FALSE | NA = NA.  Either way NA reaches
#      treat.recommend.  A rule the frame cannot express has NO membership,
#      and inventing one is worse than reporting none -- so this returns
#      `unresolved` and never a vector containing NA.  Defect D2.
#   4. EMPTY IS NOT UNRESOLVED.  A rule resolving to zero members is valid:
#      an all-FALSE vector, status "empty", target NA, partition intact.
#
# @noRd
.fs_resolve_membership <- function(df, rule, sg_def_struct = NULL) {
  n <- nrow(df)
  none <- list(in_region = NULL, status = "unresolved", missing = character(0))

  # ---- structured path (preferred) ----
  if (!is.null(sg_def_struct)) {
    cj <- sg_def_struct$conjunctions
    vars <- unique(unlist(lapply(cj, function(x) x$variable)))
    miss <- setdiff(vars, names(df))
    if (length(miss)) {
      none$missing <- miss
      return(none)
    }
    tr <- tryCatch(.grf_evaluate_subgroup(sg_def_struct, df),
                   error = function(e) NULL)
    if (is.null(tr) || length(tr) != n) return(none)
    in_region <- tr == 0L
    if (anyNA(in_region)) {
      none$missing <- miss
      return(none)
    }
    return(list(in_region = in_region,
                status = if (any(in_region)) "ok" else "empty",
                missing = character(0)))
  }

  # ---- string path ----
  if (is.null(rule) || !length(rule) || all(is.na(rule)) ||
      (length(rule) == 1L && !nzchar(rule))) {
    # No rule: not an error.  Caller treats this as "no subgroup identified".
    return(list(in_region = rep(FALSE, n), status = "empty",
                missing = character(0)))
  }

  # Dispatch order copied from get_dfpred(): disjunction FIRST, unsplit.
  sg <- if (length(rule) == 1L && grepl("|", rule, fixed = TRUE)) {
          rule
        } else if (length(rule) == 1L && grepl(" & ", rule, fixed = TRUE)) {
          strsplit(rule, " & ", fixed = TRUE)[[1L]]
        } else {
          rule
        }

  # Muffle only the two warnings that signal an unresolvable reference; any
  # other warning is a real signal and is allowed to propagate.
  pred <- withCallingHandlers(
    tryCatch(get_dfpred(df, sg), error = function(e) NULL),
    warning = function(w) {
      if (grepl("not found in data frame|Could not parse expression",
                conditionMessage(w))) invokeRestart("muffleWarning")
    })

  if (is.null(pred) || is.null(pred$treat.recommend)) {
    none$missing <- setdiff(.fs_rule_columns(rule), names(df))
    return(none)
  }

  in_region <- pred$treat.recommend == 0L
  if (length(in_region) != n || anyNA(in_region)) {
    none$missing <- setdiff(.fs_rule_columns(rule), names(df))
    return(none)
  }

  list(in_region = in_region,
       status = if (any(in_region)) "ok" else "empty",
       missing = character(0))
}


# --- LAYER 2: per-family target on a region ---------------------------------
# Each branch reproduces EXACTLY what its module computes today.  No branch
# restates arithmetic available from compute_aor(); the package already carries
# two sites for the GLM marginal effect (compute_aor() and the .effect_one()
# closure inside .compute_glm_effects()) and this adds no third.
#
# @noRd
.fs_region_effect <- function(df, idx, outcome_type, effect_measure,
                              outcome.name = NULL, event.name = NULL,
                              treat.name = NULL) {
  if (is.null(idx) || !length(idx) || !any(idx)) return(NA_real_)

  if (identical(outcome_type, "survival")) {
    d <- df[idx, , drop = FALSE]
    if (sum(d[[event.name]], na.rm = TRUE) < 5L ||
        length(unique(d[[treat.name]])) < 2L) return(NA_real_)
    fml <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s",
                                     outcome.name, event.name, treat.name))
    fit <- tryCatch(survival::coxph(fml, data = d), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    b <- unname(stats::coef(fit)[1L])
    return(if (is.finite(b)) exp(b) else NA_real_)
  }

  if (identical(outcome_type, "binary")) {
    d <- df[idx, , drop = FALSE]
    y <- d[[outcome.name]]
    if (length(unique(d[[treat.name]])) < 2L ||
        sum(y == 1L, na.rm = TRUE) < 5L ||
        sum(y == 0L, na.rm = TRUE) < 5L) return(NA_real_)
    fml <- stats::as.formula(sprintf("%s ~ %s", outcome.name, treat.name))
    fit <- tryCatch(suppressWarnings(
             stats::glm(fml, data = d, family = stats::binomial("logit"))),
           error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    cf <- stats::coef(fit)
    if (!(treat.name %in% names(cf))) return(NA_real_)
    b <- unname(cf[[treat.name]])
    return(if (is.finite(b)) exp(b) else NA_real_)
  }

  # continuous / count: potential-outcome average, delegated to the package
  compute_aor(df, subset_indicator = as.integer(idx),
              effect_measure = effect_measure)
}


# --- focus ------------------------------------------------------------------
# @noRd
.fs_check_focus <- function(focus) {
  if (missing(focus) || is.null(focus)) {
    stop("`focus` is required and has no default: supply \"harm\" or ",
         "\"benefit\".  The harm/benefit convention is an open decision; it ",
         "must be stated at the call site.", call. = FALSE)
  }
  if (!(is.character(focus) && length(focus) == 1L &&
        focus %in% c("harm", "benefit"))) {
    stop("`focus` must be exactly one of \"harm\" or \"benefit\".",
         call. = FALSE)
  }
  focus
}


#' Conditional estimand `beta(Hhat)` for one realized rule
#'
#' Computes the population value of the standard within-subgroup analysis on a
#' realized selection region and its complement, on a fixed evaluation frame.
#' This is the estimand Theorem 2 certifies coverage for, evaluated at the
#' rule a detector actually returned -- not the marginal true-subgroup target.
#'
#' Membership is resolved once, by the internal `.fs_resolve_membership()`, and the same
#' resolution serves every outcome family: membership is a property of the rule
#' and the frame, not of the effect measure.
#'
#' @section The partition invariant:
#'
#' `nH_eval + nHc_eval` must equal `nrow(frame)` on every call where the rule
#' resolves. Hhat and its complement partition the frame by construction, so a
#' sum below `nrow(frame)` means some rows fell into neither side -- which
#' happens when a membership vector carries `NA`. That is a hard error, not a
#' warning: if membership is incoherent, every number downstream of it is
#' meaningless.
#'
#' @section Unresolvable rules:
#'
#' A rule naming a column the frame lacks has no membership on that frame. The
#' record is all-`NA` with `status = "unresolved"` and `missing_cols` naming
#' the columns, rather than a partially-resolved region. Silently dropping the
#' offending clause returns a finite value for the wrong region.
#'
#' @param rule Character. The realized rule: a named `sg.harm` vector, the
#'   `" & "`-joined `sg_def` string, or a GRF disjunction
#'   `"(a & b) | (c & d)"`. `NULL` or `NA` means no subgroup was identified.
#' @param frame Data frame. The fixed evaluation population.
#' @param focus Character, required, no default. `"harm"` reads
#'   `treat.recommend == 0` as in-region; `"benefit"` inverts.
#' @param outcome_type Character. One of `"survival"`, `"binary"`,
#'   `"continuous"`, `"count"`.
#' @param effect_measure Character. Passed to [compute_aor()] for the
#'   continuous and count families; ignored for survival and binary, which keep
#'   their fitted targets.
#' @param outcome.name,event.name,treat.name Character. Column names in
#'   `frame`. `event.name` is used only by the survival family.
#' @param sg_def_struct Optional structured GRF subgroup definition. When
#'   supplied it takes precedence over `rule` and membership is computed by
#'   `.grf_evaluate_subgroup()`, which is the package's own answer for GRF
#'   membership.
#'
#' @return A one-row data frame with columns `betaHhat_H`, `betaHhat_Hc`,
#'   `nH_eval`, `nHc_eval`, `status` and `missing_cols`. The column names are
#'   scale-agnostic and identical across families so downstream
#'   `paste0("betaHhat_", suffix)` scoring is shared verbatim.
#'
#' @keywords internal
#' @export
fs_betaHhat_one <- function(rule, frame, focus,
                            outcome_type = c("survival", "binary",
                                             "continuous", "count"),
                            effect_measure = NULL,
                            outcome.name = "y_sim",
                            event.name = "event_sim",
                            treat.name = "treat_sim",
                            sg_def_struct = NULL) {
  .fs_check_focus(focus)
  outcome_type <- match.arg(outcome_type)
  stopifnot(is.data.frame(frame))
  N <- nrow(frame)

  rec <- function(H, Hc, nH, nHc, status, missing) {
    data.frame(betaHhat_H = H, betaHhat_Hc = Hc,
               nH_eval = nH, nHc_eval = nHc,
               status = status,
               missing_cols = paste(missing, collapse = ","),
               stringsAsFactors = FALSE)
  }

  res <- .fs_resolve_membership(frame, rule, sg_def_struct)

  if (identical(res$status, "unresolved")) {
    return(rec(NA_real_, NA_real_, NA_integer_, NA_integer_,
               "unresolved", res$missing))
  }

  in_region <- res$in_region
  if (identical(focus, "benefit")) in_region <- !in_region

  nH  <- sum(in_region)
  nHc <- sum(!in_region)

  # HARD INVARIANT.  Checked whenever the rule resolved; skipped only on the
  # unresolved branch above, where both counts are NA by construction.
  if (!identical(nH + nHc, N)) {
    stop(sprintf(paste0("partition invariant violated: nH_eval + nHc_eval = ",
                        "%d + %d = %d, but nrow(frame) = %d.  Membership is ",
                        "incoherent; refusing to return a target."),
                 nH, nHc, nH + nHc, N), call. = FALSE)
  }

  eff <- function(idx) .fs_region_effect(frame, idx, outcome_type,
                                         effect_measure, outcome.name,
                                         event.name, treat.name)
  rec(eff(in_region), eff(!in_region), nH, nHc, "ok", character(0))
}


# --- LAYER 4: deduplicated table and the engine-facing attach ---------------

#' Deduplicated `beta(Hhat)` targets over distinct realized rules
#'
#' Scores each **distinct** realized rule once. Deduplication is by rule, not
#' by replicate: two replicates landing on the same rule score the same target,
#' and the target does not depend on the trial.
#'
#' @section Resolution accounting:
#'
#' Six counters are attached to the result and are the fix for targets being
#' dropped silently. `n_eff` reported beside a target is what makes a reduced
#' denominator visible; a coverage figure computed on an unknown fraction of
#' replicates is not interpretable, and the fraction has to travel with it.
#'
#' A seventh, `n_reps_undetected`, counts replicates with no rule at all, so
#' that `n_reps_resolved + n_reps_unresolved + n_reps_undetected` closes
#' against `n_reps_total`.
#'
#' @param sg_defs Character vector of realized rules, one entry per replicate.
#'   `NA` or `""` marks a replicate where no subgroup was identified.
#' @param frame Data frame. The fixed evaluation population.
#' @param focus,outcome_type,effect_measure,outcome.name,event.name,treat.name
#'   Passed to [fs_betaHhat_one()].
#' @param sg_def_structs Optional named list of structured GRF definitions,
#'   keyed by the rule string, used in preference to string parsing.
#'
#' @return A data frame keyed by `sg_def`, one row per distinct rule, with the
#'   [fs_betaHhat_one()] schema. Resolution counters are attached as
#'   attributes and retrievable with [fs_betaHhat_counts()].
#'
#' @keywords internal
#' @export
fs_betaHhat_table <- function(sg_defs, frame, focus,
                              outcome_type = c("survival", "binary",
                                               "continuous", "count"),
                              effect_measure = NULL,
                              outcome.name = "y_sim",
                              event.name = "event_sim",
                              treat.name = "treat_sim",
                              sg_def_structs = NULL) {
  .fs_check_focus(focus)
  outcome_type <- match.arg(outcome_type)

  empty <- data.frame(sg_def = character(0), betaHhat_H = numeric(0),
                      betaHhat_Hc = numeric(0), nH_eval = integer(0),
                      nHc_eval = integer(0), status = character(0),
                      missing_cols = character(0), stringsAsFactors = FALSE)

  has_rule <- !is.na(sg_defs) & nzchar(sg_defs)
  u <- unique(sg_defs[has_rule])

  out <- if (!length(u)) empty else do.call(rbind, lapply(u, function(g) {
    r <- fs_betaHhat_one(g, frame, focus = focus, outcome_type = outcome_type,
                         effect_measure = effect_measure,
                         outcome.name = outcome.name, event.name = event.name,
                         treat.name = treat.name,
                         sg_def_struct = if (!is.null(sg_def_structs))
                                           sg_def_structs[[g]] else NULL)
    cbind(data.frame(sg_def = g, stringsAsFactors = FALSE), r)
  }))

  unres <- out$sg_def[out$status == "unresolved"]
  attr(out, "betaHhat_counts") <- list(
    n_rules_total      = length(u),
    n_rules_resolved   = sum(out$status != "unresolved"),
    n_rules_unresolved = sum(out$status == "unresolved"),
    n_reps_total       = length(sg_defs),
    n_reps_resolved    = sum(has_rule & !(sg_defs %in% unres)),
    n_reps_unresolved  = sum(sg_defs %in% unres),
    n_reps_undetected  = sum(!has_rule))
  out
}


#' Resolution counters from a `beta(Hhat)` table or attached results
#'
#' @param x Output of [fs_betaHhat_table()] or [fs_attach_betaHhat()].
#' @return A named list of counters, or `NULL` if none are attached.
#' @keywords internal
#' @export
fs_betaHhat_counts <- function(x) attr(x, "betaHhat_counts")


#' Attach `beta(Hhat)` targets to a results frame
#'
#' The one call each engine adds to its `run_cell()`, immediately before the
#' bundle is assembled. Runs in the main process, after the replicate loop, so
#' the evaluation frame stays in one place.
#'
#' @section Replicates with no realized rule:
#'
#' An undetected replicate is not a missing measurement -- it is a run in which
#' the whole population is the complement. It therefore receives the
#' no-subgroup record [fs_betaHhat_one()] already returns: `nH_eval = 0`,
#' `nHc_eval = nrow(frame)`, `betaHhat_Hc` the ITT effect, `betaHhat_H` `NA`
#' (an empty region has no target), and `status = "ok"`.
#'
#' This is the partition invariant reaching this layer: `nH_eval + nHc_eval`
#' now equals `nrow(frame)` on **every** row that is not `"unresolved"`, not
#' only on rows carrying a rule. Before this change such rows were all-`NA`,
#' which silently excluded them from any consumer keying on a finite target.
#' Their count is still reported separately as `n_reps_undetected`.
#'
#' @param results Data frame carrying an `sg_def` column.
#' @param frame Data frame. The fixed evaluation population.
#' @param focus,outcome_type,effect_measure,outcome.name,event.name,treat.name
#'   Passed to [fs_betaHhat_one()].
#' @param sg_def_structs Optional named list of structured GRF definitions.
#'
#' @return `results` with `betaHhat_H`, `betaHhat_Hc`, `betaHhat_status` and
#'   `nH_eval` / `nHc_eval` columns added, and counters attached.
#'
#' @keywords internal
#' @export
fs_attach_betaHhat <- function(results, frame, focus,
                               outcome_type = c("survival", "binary",
                                                "continuous", "count"),
                               effect_measure = NULL,
                               outcome.name = "y_sim",
                               event.name = "event_sim",
                               treat.name = "treat_sim",
                               sg_def_structs = NULL) {
  .fs_check_focus(focus)
  outcome_type <- match.arg(outcome_type)

  # The no-subgroup record: Hhat is empty, so Hhat^c IS the ITT population.
  # Computed once -- it depends only on the frame, not on any replicate.
  none <- fs_betaHhat_one(NULL, frame, focus = focus,
                          outcome_type = outcome_type,
                          effect_measure = effect_measure,
                          outcome.name = outcome.name,
                          event.name = event.name, treat.name = treat.name)

  if (is.null(results$sg_def)) {
    results$betaHhat_H      <- none$betaHhat_H
    results$betaHhat_Hc     <- none$betaHhat_Hc
    results$betaHhat_status <- none$status
    results$nH_eval         <- none$nH_eval
    results$nHc_eval        <- none$nHc_eval
    attr(results, "betaHhat_counts") <- list(
      n_rules_total = 0L, n_rules_resolved = 0L, n_rules_unresolved = 0L,
      n_reps_total = nrow(results), n_reps_resolved = 0L,
      n_reps_unresolved = 0L, n_reps_undetected = nrow(results))
    return(results)
  }

  bt <- fs_betaHhat_table(results$sg_def, frame, focus = focus,
                          outcome_type = outcome_type,
                          effect_measure = effect_measure,
                          outcome.name = outcome.name,
                          event.name = event.name, treat.name = treat.name,
                          sg_def_structs = sg_def_structs)
  j <- match(results$sg_def, bt$sg_def)
  results$betaHhat_H      <- bt$betaHhat_H[j]
  results$betaHhat_Hc     <- bt$betaHhat_Hc[j]
  results$betaHhat_status <- bt$status[j]
  results$nH_eval         <- bt$nH_eval[j]
  results$nHc_eval        <- bt$nHc_eval[j]

  # Undetected replicates (no realized rule) take the no-subgroup record, not
  # all-NA: the whole frame is the complement, so the target exists.
  und <- is.na(results$sg_def) | !nzchar(results$sg_def)
  if (any(und)) {
    results$betaHhat_H[und]      <- none$betaHhat_H
    results$betaHhat_Hc[und]     <- none$betaHhat_Hc
    results$betaHhat_status[und] <- none$status
    results$nH_eval[und]         <- none$nH_eval
    results$nHc_eval[und]        <- none$nHc_eval
  }
  attr(results, "betaHhat_counts") <- attr(bt, "betaHhat_counts")
  results
}


# --- LAYER 5: n_eff parity ---------------------------------------------------

#' Refuse to print coverage targets computed on different denominators
#'
#' `.coverage_meta()`-style assembly computes
#' `ok <- is.finite(target) & is.finite(lo) & is.finite(hi)` and
#' `n_eff <- sum(ok)`, so a non-finite target is dropped with no error and no
#' warning. Coverage is then computed over the replicates whose targets
#' happened to score, which is not a random subset: whatever made a target
#' non-finite is a property of the realized rule.
#'
#' Because a reference target that is a finite scalar is dropped only when the
#' interval itself is non-finite, `n_eff` parity between the two is the exact
#' check for a target-specific loss, and it is free -- `n_eff` is already in
#' the table.
#'
#' @section Two causes, deliberately not distinguished:
#'
#' A parity break can mean an unresolvable rule, or a genuinely degenerate
#' region whose per-family guard returned `NA` by design. Both make the two
#' coverage figures incomparable as printed, so both stop. Use
#' [fs_betaHhat_neff_report()] on the bundle to tell them apart before
#' deciding.
#'
#' @param cov_df Coverage table with `target`, `n_eff`, `block` and
#'   `estimator` columns.
#' @param strict Logical. `TRUE` stops; `FALSE` warns, for a run whose cause is
#'   known.
#' @param target_beta,target_ref Character. The target to check and the
#'   reference to check it against.
#'
#' @return `cov_df`, invisibly.
#'
#' @keywords internal
#' @export
fs_betaHhat_neff_parity <- function(cov_df, strict = TRUE,
                                    target_beta = "C_betaHhat",
                                    target_ref = "C_dagger") {
  need <- c("target", "n_eff", "block", "estimator")
  if (!all(need %in% names(cov_df))) return(invisible(cov_df))
  keycols <- intersect(c("src", "subgroup_method", "n_sample", "block",
                         "estimator"), names(cov_df))
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
      "  Its coverage is NOT comparable to the row beside it.\n",
      paste(utils::capture.output(print(bad, row.names = FALSE)),
            collapse = "\n"),
      "\n  Run fs_betaHhat_neff_report() on the bundle to see which rules ",
      "produced the\n  NA targets and whether they are unresolvable rules or ",
      "degenerate regions.")
    if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(cov_df)
}


#' Which rules produced NA `beta(Hhat)` targets, and why
#'
#' Diagnostic companion to [fs_betaHhat_neff_parity()]. Separates the two
#' causes of a parity break: a rule the evaluation frame cannot express, and a
#' region small enough for a per-family guard to fire.
#'
#' @param bundle A results data frame, a bundle list with a `results` element,
#'   or a path to an `.rds` holding either.
#' @param block Character. `"H"` or `"Hc"`.
#'
#' @return A data frame of distinct offending rules with `is_disjunction`,
#'   `nH_eval` and `status` where available, or `NULL` when there are none.
#'
#' @keywords internal
#' @export
fs_betaHhat_neff_report <- function(bundle, block = c("H", "Hc")) {
  if (is.character(bundle)) bundle <- readRDS(bundle)
  r <- if (is.data.frame(bundle)) bundle else bundle$results
  block <- match.arg(block)
  col <- paste0("betaHhat_", block)
  if (is.null(r[[col]]) || is.null(r$sg_def)) return(invisible(NULL))
  bad <- r[is.na(r[[col]]) & !is.na(r$sg_def) & nzchar(r$sg_def), ,
           drop = FALSE]
  if (!nrow(bad)) {
    message(sprintf("no NA %s among rows with a realized rule", col))
    return(invisible(NULL))
  }
  out <- data.frame(
    sg_def = bad$sg_def,
    is_disjunction = grepl("|", bad$sg_def, fixed = TRUE),
    nH_eval = if (!is.null(bad$nH_eval)) bad$nH_eval else NA_integer_,
    status = if (!is.null(bad$betaHhat_status)) bad$betaHhat_status
             else NA_character_,
    stringsAsFactors = FALSE)
  out[!duplicated(out$sg_def), , drop = FALSE]
}


# --- evaluation frame and theta-dagger, per family ---------------------------

#' The frame `beta(Hhat)` is scored on
#'
#' One entry point per outcome family, dispatching on `outcome_type` rather
#' than exposing a separate function per family. Defaults reproduce the
#' simulation modules this replaces exactly.
#'
#' @section What each family returns:
#'
#' * **survival** -- the entire fixed super-population, every subject exactly
#'   once (`replace = FALSE`), under the same randomized/censored analysis the
#'   trials run.
#' * **binary** -- the same construction on the GLM simulator: every subject
#'   once, under one fixed treatment/outcome realization.
#' * **continuous / count** -- `dgm$df_super` **unchanged**. No simulation
#'   occurs. The mean difference is collapsible, so the target is an exact
#'   finite mean over the super-population: the scoring frame *is* the
#'   population. `eval_seed` is accepted and ignored here so that generic
#'   harness code need not branch, and the target carries **zero Monte Carlo
#'   error** -- uniformity of the call surface is not sameness of the object.
#'
#' @section Rejected arguments:
#'
#' `analysis_time`, `cens_adjust` and `n_eval` are survival-only. Supplying any
#' of them on a non-survival path is an error rather than a silent no-op:
#' quietly ignoring an argument that means something on another path is how
#' conventions drift apart.
#'
#' @param dgm A DGM object carrying `df_super`.
#' @param outcome_type Character. One of `"survival"`, `"binary"`,
#'   `"continuous"`, `"count"`.
#' @param eval_seed Integer. Fixes the single realization on the full pool.
#'   Ignored for `"continuous"` and `"count"`.
#' @param analysis_time,cens_adjust Survival only.
#' @param n_eval Defunct. Non-`NULL` is a hard error: this always evaluates the
#'   full super-population, so both engines score the identical target.
#'
#' @return A data frame: the evaluation population.
#'
#' @keywords internal
#' @export
fs_build_eval_frame <- function(dgm,
                                outcome_type = c("survival", "binary",
                                                 "continuous", "count"),
                                eval_seed = 20260628L,
                                analysis_time = 84,
                                cens_adjust = log(1.5),
                                n_eval = NULL) {
  outcome_type <- match.arg(outcome_type)
  if (is.null(dgm$df_super))
    stop("`dgm` has no `df_super`.", call. = FALSE)

  if (!is.null(n_eval))
    stop("fs_build_eval_frame() evaluates the FULL super-population ",
         "(each subject once); the legacy `n_eval` argument has been removed. ",
         "Drop n_eval from the call so both engines score the identical ",
         "full-pool beta(Hhat).", call. = FALSE)

  if (!identical(outcome_type, "survival")) {
    surv_only <- c(analysis_time = !missing(analysis_time),
                   cens_adjust   = !missing(cens_adjust))
    if (any(surv_only))
      stop(sprintf(paste0("`%s` is survival-only and was supplied with ",
                          "outcome_type = \"%s\".  Ignoring it silently is how ",
                          "conventions drift; drop it from the call."),
                   paste(names(surv_only)[surv_only], collapse = "`, `"),
                   outcome_type), call. = FALSE)
  }

  if (identical(outcome_type, "survival")) {
    return(simulate_from_dgm(dgm, n = nrow(dgm$df_super), replace = FALSE,
                             analysis_time = analysis_time,
                             cens_adjust = cens_adjust, seed = eval_seed))
  }
  if (identical(outcome_type, "binary")) {
    return(simulate_from_glm_dgm(dgm, n = nrow(dgm$df_super),
                                 replace = FALSE, seed = eval_seed))
  }
  # continuous / count: the super-population IS the scoring frame.
  dgm$df_super
}


#' `theta-dagger` at the true subgroup flag, on the scoring frame
#'
#' The marginal target at the DGM's own harm flag, computed on the same frame
#' `beta(Hhat)` is scored on. Use it as a sanity gate: it should reproduce the
#' DGM's own subgroup effects.
#'
#' For `"continuous"` and `"count"` this is an **exact identity**, not
#' agreement to Monte Carlo error -- the frame is the super-population and the
#' arithmetic is the same [compute_aor()] dispatch the DGM used. A tolerance
#' there would hide a real defect.
#'
#' This is a thin dispatch onto the same per-family effect used for every
#' region; it introduces no arithmetic of its own.
#'
#' @param frame Data frame. The evaluation population, normally from
#'   [fs_build_eval_frame()].
#' @param outcome_type Character. One of `"survival"`, `"binary"`,
#'   `"continuous"`, `"count"`.
#' @param harm.name Character. The true-subgroup flag column; `== 1L` is in.
#' @param outcome.name,event.name,treat.name Character. Column names in
#'   `frame`. `event.name` is survival only.
#' @param effect_measure Character. Required for `"continuous"` and
#'   `"count"`; ignored for `"survival"` and `"binary"`, which keep their
#'   fitted targets.
#'
#' @return A named numeric vector, `thetaDagger_H` and `thetaDagger_Hc`.
#'
#' @keywords internal
#' @export
fs_betaHhat_theta_dagger_check <- function(frame,
                                           outcome_type = c("survival",
                                                            "binary",
                                                            "continuous",
                                                            "count"),
                                           harm.name = "flag_harm",
                                           outcome.name = "y_sim",
                                           event.name = "event_sim",
                                           treat.name = "treat_sim",
                                           effect_measure = NULL) {
  outcome_type <- match.arg(outcome_type)
  stopifnot(is.data.frame(frame))
  if (is.null(frame[[harm.name]]))
    stop("`frame` has no column \"", harm.name, "\".", call. = FALSE)
  if (outcome_type %in% c("continuous", "count") && is.null(effect_measure))
    stop("`effect_measure` is required for outcome_type = \"", outcome_type,
         "\"; there is no default to guess at.", call. = FALSE)

  inH <- frame[[harm.name]] == 1L
  c(thetaDagger_H  = .fs_region_effect(frame,  inH, outcome_type,
                                       effect_measure, outcome.name,
                                       event.name, treat.name),
    thetaDagger_Hc = .fs_region_effect(frame, !inH, outcome_type,
                                       effect_measure, outcome.name,
                                       event.name, treat.name))
}
