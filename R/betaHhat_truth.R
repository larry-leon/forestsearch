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
#' Membership is resolved once, by [.fs_resolve_membership()], and the same
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
