# =============================================================================
# probe_threshold_sync_2026-09-18.R
#
# Replicate-equality probe for TASK_threshold_sync_2026-09-18.md, Step 2.
#
# NO COMPUTE.  Nothing here fits a model, resamples, or touches data.  The
# probe constructs forestsearch() argument lists and resolves the two effect
# thresholds from them, twice:
#
#   parent     -- the argument list as a user writes it (an argument the user
#                 did not write is genuinely absent, so missing() is TRUE);
#   replicate  -- the argument list as bootstrap_analysis_dofuture.R builds
#                 args_FS_boot from fs.est$args_call_all (and as
#                 forestsearch_cross_validation.R builds cv_args): EVERY
#                 formal is present, so missing() is FALSE for all of them.
#
# The resolution below is a transcription of forestsearch() at the cited
# lines.  The four threshold DEFAULTS are read from formals(forestsearch) so
# the probe cannot drift from the source on those values; the branch structure
# is transcribed and is checked by the acceptance tests.
#
#   forestsearch_main.R:1446-1452   effect_measure default per outcome_type
#   forestsearch_main.R:1462-1465   user_set_* detection + alias merge
#   forestsearch_main.R:1512-1513   args_call_all <- mget(names(formals()))
#   forestsearch_main.R:1896-1990   threshold resolution (GLM branches)
#   forestsearch_main.R:2053-2067   threshold_config (GLM)
#   forestsearch_main.R:2103-2117   threshold_config (survival)
#   bootstrap_analysis_dofuture.R:406,558,614   the replay
#   forestsearch_cross_validation.R:345,388-421,479   the CV replay
# =============================================================================

.probe_fs_defaults <- function() {
  fm <- formals(forestsearch)
  list(
    effect.threshold      = eval(fm[["effect.threshold"]]),
    consistency.threshold = eval(fm[["consistency.threshold"]]),
    hr.threshold          = eval(fm[["hr.threshold"]]),
    hr.consistency        = eval(fm[["hr.consistency"]]),
    outcome_type          = eval(fm[["outcome_type"]])[[1L]]
  )
}

# ---------------------------------------------------------------------------
# Resolve the thresholds from one argument list.
#
# `supplied` is a named list.  A name PRESENT in the list counts as supplied
# (so missing() is FALSE), even when its value is NULL -- that is exactly the
# state of args_call_all's NULL-defaulted entries after mget().
# ---------------------------------------------------------------------------
.probe_resolve <- function(supplied) {
  d <- .probe_fs_defaults()
  has  <- function(n) n %in% names(supplied)
  get1 <- function(n, dflt) if (has(n)) supplied[[n]] else dflt

  outcome_type   <- get1("outcome_type",   d$outcome_type)
  effect_measure <- get1("effect_measure", NULL)

  effect.threshold      <- get1("effect.threshold",      d$effect.threshold)
  consistency.threshold <- get1("consistency.threshold", d$consistency.threshold)
  hr.threshold          <- get1("hr.threshold",          d$hr.threshold)
  hr.consistency        <- get1("hr.consistency",        d$hr.consistency)

  # :1446-1452
  if (outcome_type != "survival" && is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      binary = "OR", continuous = "MD", count = "IRR")
  }

  # :1462-1465 -- detection, then merge.  Detection for the legacy spellings
  # is missing()-based, which is what the replay defeats.
  user_set_threshold   <- !is.null(effect.threshold)      || has("hr.threshold")
  user_set_consistency <- !is.null(consistency.threshold) || has("hr.consistency")
  if (!is.null(effect.threshold))      hr.threshold   <- effect.threshold
  if (!is.null(consistency.threshold)) hr.consistency <- consistency.threshold

  warns <- character(0)
  note  <- NA_character_

  if (outcome_type != "survival") {

    # :1896-1897
    effect_threshold      <- hr.threshold
    consistency_threshold <- hr.consistency

    if (effect_measure %in% c("RD", "IRD")) {
      # :1902-1907
      if (!user_set_threshold && hr.threshold == 1.25)
        effect_threshold <- switch(effect_measure, RD = 0.05, IRD = 0.01)
      if (!user_set_consistency && hr.consistency == 1.0)
        consistency_threshold <- 0.0
      # :1912-1928
      if (user_set_threshold && effect_threshold > 1.0) {
        warns <- c(warns, "effect.threshold-ratio-scale-remap")
        effect_threshold <- switch(effect_measure, RD = 0.05, IRD = 0.01)
      }
      # :1930-1945
      if (consistency_threshold >= 1.0 && !user_set_consistency) {
        consistency_threshold <- 0.0
      } else if (consistency_threshold > 1.0) {
        warns <- c(warns, "consistency.threshold-ratio-scale-remap")
        consistency_threshold <- 0.0
      }

    } else if (effect_measure == "MD") {
      # :1953-1959
      if (!user_set_threshold   && hr.threshold   == 1.25) effect_threshold      <- 0.0
      if (!user_set_consistency && hr.consistency == 1.0)  consistency_threshold <- 0.0

    } else {
      # :1968-1988 -- ratio measures: reject non-positive, then log()
      if (effect_threshold <= 0)
        return(.probe_row_error("effect.threshold <= 0 for ratio measure"))
      if (consistency_threshold <= 0)
        return(.probe_row_error("consistency.threshold <= 0 for ratio measure"))
      effect_threshold      <- log(effect_threshold)
      consistency_threshold <- log(consistency_threshold)
    }

    # :2053-2067
    is_identity <- effect_measure %in% c("RD", "IRD", "MD")
    out <- list(
      screening           = effect_threshold,
      consistency         = consistency_threshold,
      screening_natural   = if (is_identity) effect_threshold   else exp(effect_threshold),
      consistency_natural = if (is_identity) consistency_threshold else exp(consistency_threshold),
      scale               = if (is_identity) "identity" else "log")

  } else {
    # :2103-2117
    out <- list(
      screening           = log(hr.threshold),
      consistency         = log(max(hr.consistency, 0.001)),
      screening_natural   = hr.threshold,
      consistency_natural = hr.consistency,
      scale               = "log")
  }

  out$user_set_threshold   <- user_set_threshold
  out$user_set_consistency <- user_set_consistency
  out$warnings <- if (length(warns)) paste(warns, collapse = "|") else ""
  out$error    <- NA_character_
  out$note     <- note
  out
}

.probe_row_error <- function(msg) {
  list(screening = NA_real_, consistency = NA_real_,
       screening_natural = NA_real_, consistency_natural = NA_real_,
       scale = NA_character_, user_set_threshold = NA,
       user_set_consistency = NA, warnings = "", error = msg, note = NA_character_)
}

# ---------------------------------------------------------------------------
# Build the replay argument list the way the bootstrap / CV build it:
# args_call_all = mget(names(formals())) taken AFTER the alias merge, so every
# threshold name is present.  `sync = TRUE` additionally writes the parent's
# resolved NATURAL-scale values into the two NULL-defaulted spellings, which is
# the Step 3 mechanism under test.
# ---------------------------------------------------------------------------
.probe_args_call_all <- function(supplied, sync = FALSE) {
  d <- .probe_fs_defaults()
  has  <- function(n) n %in% names(supplied)
  get1 <- function(n, dflt) if (has(n)) supplied[[n]] else dflt

  outcome_type   <- get1("outcome_type",   d$outcome_type)
  effect_measure <- get1("effect_measure", NULL)
  if (outcome_type != "survival" && is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      binary = "OR", continuous = "MD", count = "IRR")
  }

  eff <- get1("effect.threshold",      d$effect.threshold)
  con <- get1("consistency.threshold", d$consistency.threshold)
  hrt <- get1("hr.threshold",          d$hr.threshold)
  hrc <- get1("hr.consistency",        d$hr.consistency)
  if (!is.null(eff)) hrt <- eff
  if (!is.null(con)) hrc <- con

  a <- list()
  a["effect.threshold"]      <- list(eff)
  a["consistency.threshold"] <- list(con)
  a$hr.threshold             <- hrt
  a$hr.consistency           <- hrc
  a$outcome_type             <- outcome_type
  a["effect_measure"]        <- list(effect_measure)

  if (isTRUE(sync)) {
    p <- .probe_resolve(supplied)
    if (is.na(p$error)) {
      a["effect.threshold"]      <- list(p$screening_natural)
      a["consistency.threshold"] <- list(p$consistency_natural)
    }
  }
  a
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

# Custom (non-default) supply values, measure-appropriate.
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

  d <- .probe_fs_defaults()
  v <- if (identical(values, "default"))
         c(c1 = d$hr.threshold, c2 = d$hr.consistency)
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

#' Run the full replicate-equality probe.
#'
#' @param sync TRUE to build the replay list the Step 3 way (resolved naturals
#'   written into effect.threshold / consistency.threshold).
#' @return a data.frame, one row per cell.
probe_threshold_sync <- function(sync = FALSE) {
  grid <- expand.grid(
    est      = seq_along(.PROBE_ESTIMANDS),
    spelling = c("none", "legacy", "new"),
    subset   = c("neither", "c1", "both"),
    values   = c("default", "custom"),
    stringsAsFactors = FALSE)

  # 'none' carries no value variant and no subset variant; keep one row each.
  keep <- !(grid$spelling == "none" & (grid$subset != "neither" | grid$values != "default"))
  # a supplied spelling with subset 'neither' is the same cell as 'none'
  keep <- keep & !(grid$spelling != "none" & grid$subset == "neither")
  grid <- grid[keep, , drop = FALSE]

  rows <- lapply(seq_len(nrow(grid)), function(i) {
    est <- .PROBE_ESTIMANDS[[grid$est[i]]]
    sup <- .probe_supplied_list(est, grid$spelling[i], grid$subset[i], grid$values[i])
    par <- .probe_resolve(sup)
    rep <- .probe_resolve(.probe_args_call_all(sup, sync = sync))

    flds <- c("screening", "consistency", "screening_natural",
              "consistency_natural", "scale")
    pv <- vapply(flds, function(f) .probe_fmt(par[[f]]), character(1))
    rv <- vapply(flds, function(f) .probe_fmt(rep[[f]]), character(1))

    data.frame(
      estimand = est$key,
      spelling = grid$spelling[i],
      subset   = if (grid$spelling[i] == "none") "neither" else grid$subset[i],
      values   = if (grid$spelling[i] == "none") "-" else grid$values[i],
      parent_screening      = pv[["screening"]],
      parent_consistency    = pv[["consistency"]],
      parent_scr_natural    = pv[["screening_natural"]],
      parent_con_natural    = pv[["consistency_natural"]],
      parent_scale          = pv[["scale"]],
      parent_error          = par$error,
      parent_warnings       = par$warnings,
      rep_screening         = rv[["screening"]],
      rep_consistency       = rv[["consistency"]],
      rep_scr_natural       = rv[["screening_natural"]],
      rep_con_natural       = rv[["consistency_natural"]],
      rep_scale             = rv[["scale"]],
      rep_error             = rep$error,
      rep_warnings          = rep$warnings,
      violation = !identical(unname(pv), unname(rv)) ||
                  !identical(par$error, rep$error),
      stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
