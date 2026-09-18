# =============================================================================
# helper-threshold-sync.R
#
# Replicate-equality probe for TASK_threshold_sync_2026-09-18.md, Steps 2 & 5,
# and the fixture behind test-threshold-sync.R.  dev/tasks/probe_threshold_
# sync_2026-09-18.R sources this file; there is no second copy.
#
# NO COMPUTE.  Nothing here fits a model, resamples, or touches data.  The
# probe constructs forestsearch() argument lists and resolves the two effect
# thresholds from them, twice:
#
#   parent     -- the argument list as a user writes it, so an argument the
#                 user did not write is genuinely absent and missing() is TRUE;
#   replicate  -- the argument list as bootstrap_analysis_dofuture.R builds
#                 args_FS_boot from fs.est$args_call_all (and as
#                 forestsearch_cross_validation.R builds cv_args): EVERY
#                 formal is present, so missing() is FALSE for all of them.
#
# The resolution is NOT transcribed.  .probe_resolver() lifts the relevant
# statements out of body(forestsearch) and rebuilds them as a small function
# whose formals are forestsearch()'s own threshold formals -- so missing(),
# the alias merge and every remap branch are the package's code, evaluated
# under the package namespace, and the probe cannot drift from the source.
#
#   forestsearch_main.R  effect_measure default per outcome_type
#                        user_set_* detection + alias merge
#                        threshold resolution (GLM branches)
#                        SECTION 2B-ii threshold sync (present only after the
#                          Step 3 edit; the probe detects its absence)
#   bootstrap_analysis_dofuture.R:406, 558, 614       the bootstrap replay
#   forestsearch_cross_validation.R:345, 479; 859, 1001   the CV replay
# =============================================================================

# ---------------------------------------------------------------------------
# Statement lifting
# ---------------------------------------------------------------------------

.probe_deparse1 <- function(e) paste(deparse(e), collapse = " ")

.probe_pick <- function(stmts, pattern, what, n = 1L) {
  d  <- vapply(stmts, .probe_deparse1, character(1))
  ix <- which(grepl(pattern, d, fixed = TRUE))
  if (length(ix) != n)
    stop(sprintf("probe: expected %d statement(s) matching '%s' (%s), found %d",
                 n, pattern, what, length(ix)), call. = FALSE)
  ix
}

#' Build the resolver: forestsearch()'s own threshold code, as a function.
#'
#' @param sync TRUE to include the SECTION 2B-ii sync statements.
#' @return a function(effect.threshold, consistency.threshold, hr.threshold,
#'   hr.consistency, outcome_type, effect_measure) returning the resolved
#'   record plus the args_call_all entries a replay would receive.
.probe_resolver <- function(sync = FALSE) {
  fm   <- formals(forestsearch)
  top  <- as.list(body(forestsearch))[-1L]   # drop the `{`

  # -- effect_measure default (top level) ----------------------------------
  i_em <- .probe_pick(top, "is.null(effect_measure)", "effect_measure default")
  s_em <- top[[i_em]]

  # -- user_set_* detection + alias merge (top level, four statements) ------
  i_ust <- .probe_pick(top, "user_set_threshold <- ",   "user_set_threshold")
  i_usc <- .probe_pick(top, "user_set_consistency <- ", "user_set_consistency")
  i_a1  <- .probe_pick(top, "hr.threshold <- effect.threshold",
                       "alias merge c1")
  i_a2  <- .probe_pick(top, "hr.consistency <- consistency.threshold",
                       "alias merge c2")
  s_flags <- top[c(i_ust, i_usc, i_a1, i_a2)]

  # -- the GLM resolution block (nested in `if (outcome_type != "survival")`)
  i_glm <- .probe_pick(top, 'if (outcome_type != "survival") {',
                       "GLM branch", n = 1L)
  inner <- as.list(top[[i_glm]][[3L]])[-1L]
  j_e   <- .probe_pick(inner, "effect_threshold <- hr.threshold",
                       "effect_threshold seed")
  j_c   <- .probe_pick(inner, "consistency_threshold <- hr.consistency",
                       "consistency_threshold seed")
  j_br  <- .probe_pick(inner, 'if (effect_measure %in% c("RD", "IRD")) {',
                       "measure branch")
  s_res <- inner[c(j_e, j_c, j_br)]

  # -- the sync block (top level; absent before the Step 3 edit) ------------
  s_sync <- list()
  if (isTRUE(sync)) {
    i_id <- .probe_pick(top, ".fs_identity_scale <- ", "identity-scale flag")
    i_s1 <- .probe_pick(top, "effect.threshold <- if (.fs_identity_scale)",
                        "sync c1")
    i_s2 <- .probe_pick(top, "consistency.threshold <- if (.fs_identity_scale)",
                        "sync c2")
    i_sy <- .probe_pick(top, '"effect.threshold", "consistency.threshold"',
                        "sync call")
    s_sync <- top[c(i_id, i_s1, i_s2, i_sy)]
  }

  body_exprs <- c(
    list(quote(outcome_type <- match.arg(
      outcome_type, c("survival", "binary", "continuous", "count")))),
    list(s_em),
    s_flags,
    list(quote(args_call_all <- list(
      effect.threshold      = effect.threshold,
      consistency.threshold = consistency.threshold,
      hr.threshold          = hr.threshold,
      hr.consistency        = hr.consistency,
      outcome_type          = outcome_type,
      effect_measure        = effect_measure))),
    list(quote(effect_threshold <- NULL),
         quote(consistency_threshold <- NULL)),
    list(as.call(c(quote(`if`), quote(outcome_type != "survival"),
                   as.call(c(quote(`{`), s_res))))),
    s_sync,
    list(quote(.probe_record(
      outcome_type, effect_measure, hr.threshold, hr.consistency,
      effect_threshold, consistency_threshold,
      user_set_threshold, user_set_consistency, args_call_all)))
  )

  f <- as.function(c(
    fm[c("effect.threshold", "consistency.threshold",
         "hr.threshold", "hr.consistency")],
    alist(outcome_type = "survival", effect_measure = NULL),
    as.call(c(quote(`{`), body_exprs))))
  environment(f) <- asNamespace("forestsearch")
  f
}

# ---------------------------------------------------------------------------
# The record compared between parent and replicate.
#
# `screening` / `consistency` are on the comparison scale the search uses;
# the naturals and `scale` mirror threshold_config (forestsearch_main.R,
# SECTION 2B "Store resolved config for the return object").
# ---------------------------------------------------------------------------
.probe_record <- function(outcome_type, effect_measure,
                          hr.threshold, hr.consistency,
                          effect_threshold, consistency_threshold,
                          user_set_threshold, user_set_consistency,
                          args_call_all) {
  if (outcome_type == "survival") {
    out <- list(screening           = log(hr.threshold),
                consistency         = log(max(hr.consistency, 0.001)),
                screening_natural   = hr.threshold,
                consistency_natural = hr.consistency,
                scale               = "log")
  } else {
    is_identity <- effect_measure %in% c("RD", "IRD", "MD")
    out <- list(
      screening           = effect_threshold,
      consistency         = consistency_threshold,
      screening_natural   = if (is_identity) effect_threshold
                            else exp(effect_threshold),
      consistency_natural = if (is_identity) consistency_threshold
                            else exp(consistency_threshold),
      scale               = if (is_identity) "identity" else "log")
  }
  out$user_set_threshold   <- user_set_threshold
  out$user_set_consistency <- user_set_consistency
  out$args_call_all        <- args_call_all
  out
}

.probe_eval <- function(resolver, supplied) {
  warns <- character(0)
  res <- withCallingHandlers(
    tryCatch(do.call(resolver, supplied, quote = TRUE),
             error = function(e) structure(list(error = conditionMessage(e)),
                                           class = "probe_error")),
    warning = function(w) {
      warns <<- c(warns, .probe_warn_tag(conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  if (inherits(res, "probe_error")) {
    return(list(screening = NA_real_, consistency = NA_real_,
                screening_natural = NA_real_, consistency_natural = NA_real_,
                scale = NA_character_, args_call_all = NULL,
                warnings = "", error = res$error))
  }
  res$warnings <- if (length(warns)) paste(unique(warns), collapse = "|") else ""
  res$error    <- NA_character_
  res
}

.probe_warn_tag <- function(msg) {
  if (grepl("^effect.threshold", msg))      "c1-ratio-scale-remap"
  else if (grepl("^consistency.threshold", msg)) "c2-ratio-scale-remap"
  else "other"
}

# ---------------------------------------------------------------------------
# The matrix
# ---------------------------------------------------------------------------
.PROBE_ESTIMANDS <- list(
  list(key = "survival",       outcome_type = "survival",   effect_measure = NULL),
  list(key = "binary-unset",   outcome_type = "binary",     effect_measure = NULL),
  list(key = "binary-RD",      outcome_type = "binary",     effect_measure = "RD"),
  list(key = "binary-OR",      outcome_type = "binary",     effect_measure = "OR"),
  list(key = "continuous-MD",  outcome_type = "continuous", effect_measure = "MD"),
  list(key = "count-IRR",      outcome_type = "count",      effect_measure = "IRR"),
  list(key = "count-IRD",      outcome_type = "count",      effect_measure = "IRD")
)

.probe_custom_values <- function(key) {
  switch(key,
    "binary-RD"     = c(c1 = 0.07, c2 = 0.03),
    "count-IRD"     = c(c1 = 0.02, c2 = 0.01),
    "continuous-MD" = c(c1 = 30,   c2 = 10),
    c(c1 = 1.5, c2 = 1.2))
}

.probe_supplied_list <- function(est, spelling, subset, values = "default") {
  s <- list(outcome_type = est$outcome_type)
  if (!is.null(est$effect_measure)) s$effect_measure <- est$effect_measure
  if (identical(spelling, "none") || identical(subset, "neither")) return(s)

  fm <- formals(forestsearch)
  v  <- if (identical(values, "default"))
          c(c1 = eval(fm[["hr.threshold"]]), c2 = eval(fm[["hr.consistency"]]))
        else .probe_custom_values(est$key)

  nm_c1 <- if (identical(spelling, "legacy")) "hr.threshold"   else "effect.threshold"
  nm_c2 <- if (identical(spelling, "legacy")) "hr.consistency" else "consistency.threshold"
  s[[nm_c1]] <- unname(v[["c1"]])
  if (identical(subset, "both")) s[[nm_c2]] <- unname(v[["c2"]])
  s
}

.probe_fmt <- function(x) {
  if (is.null(x) || length(x) != 1L) return(NA_character_)
  if (is.character(x)) return(x)
  if (is.na(x)) return(NA_character_)
  formatC(as.numeric(x), format = "g", digits = 17)
}

.probe_grid <- function() {
  g <- expand.grid(est = seq_along(.PROBE_ESTIMANDS),
                   spelling = c("none", "legacy", "new"),
                   subset   = c("neither", "c1", "both"),
                   values   = c("default", "custom"),
                   stringsAsFactors = FALSE)
  keep <- !(g$spelling == "none" & (g$subset != "neither" | g$values != "default"))
  keep <- keep & !(g$spelling != "none" & g$subset == "neither")
  g[keep, , drop = FALSE]
}

#' Run the replicate-equality probe.
#'
#' @param sync TRUE to include the SECTION 2B-ii sync statements when building
#'   the resolver -- i.e. to probe the edited tree.
#' @return a data.frame, one row per cell of the matrix.
probe_threshold_sync <- function(sync = FALSE) {
  resolver <- .probe_resolver(sync = sync)
  g <- .probe_grid()

  rows <- lapply(seq_len(nrow(g)), function(i) {
    est <- .PROBE_ESTIMANDS[[g$est[i]]]
    sup <- .probe_supplied_list(est, g$spelling[i], g$subset[i], g$values[i])

    par <- .probe_eval(resolver, sup)

    # The replay: every formal supplied, exactly as args_FS_boot is built.
    rep <- if (is.null(par$args_call_all)) par
           else .probe_eval(resolver, par$args_call_all)

    flds <- c("screening", "consistency", "screening_natural",
              "consistency_natural", "scale")
    pv <- vapply(flds, function(f) .probe_fmt(par[[f]]), character(1))
    rv <- vapply(flds, function(f) .probe_fmt(rep[[f]]), character(1))

    a <- par$args_call_all
    data.frame(
      estimand = est$key,
      spelling = g$spelling[i],
      subset   = if (g$spelling[i] == "none") "neither" else g$subset[i],
      values   = if (g$spelling[i] == "none") "-" else g$values[i],
      parent_screening   = pv[["screening"]],
      parent_consistency = pv[["consistency"]],
      parent_scr_natural = pv[["screening_natural"]],
      parent_con_natural = pv[["consistency_natural"]],
      parent_scale       = pv[["scale"]],
      parent_error       = par$error,
      parent_warnings    = par$warnings,
      rep_screening      = rv[["screening"]],
      rep_consistency    = rv[["consistency"]],
      rep_scr_natural    = rv[["screening_natural"]],
      rep_con_natural    = rv[["consistency_natural"]],
      rep_scale          = rv[["scale"]],
      rep_error          = rep$error,
      rep_warnings       = rep$warnings,
      acall_effect.threshold =
        .probe_fmt(if (is.null(a)) NULL else a$effect.threshold),
      acall_consistency.threshold =
        .probe_fmt(if (is.null(a)) NULL else a$consistency.threshold),
      violation = !identical(unname(pv), unname(rv)) ||
                  !identical(par$error, rep$error),
      stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
