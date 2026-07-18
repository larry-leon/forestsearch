# guohe_from_forestsearch.R
#
# Bridge: apply the Guo & He (2021) de-biasing comparator to a fitted
# `forestsearch` object. Requires `guohe_algorithm3.R` (and, for adaptive r,
# `guohe_adaptive_r.R`) to be sourced first.
#
# ---------------------------------------------------------------------------
# SCOPE -- read this before using the output.
#
# Guo & He is argmax-specific. Their estimand is beta_s = beta at
# argmax_i beta_hat_i; their offsets d_i = (1 - n^(r-1/2))(beta_hat_max -
# beta_hat_i) are anchored to the maximum; and their limiting law is
# max_{i in H} T_i over H = {i : beta_i = beta_max}. Section 3 generalises to
# gamma_max = sup_c beta(c) -- still a supremum.
#
# `forestsearch()`'s `sg_focus` therefore matters:
#
#   * "maxeff"   -- which.max(beta) among passers. THIS is the argmax primitive,
#                   and the configuration Guo & He covers.
#   * "effMaxSG" -- max SIZE within an effect band (the neighbourhood rule), and
#                   "pareto" / "both" bands. These are NOT suprema of the
#                   effect, so none of the four pillars above transfers.
#                   Extending Guo & He to them would be new theory, not new code.
#
# This bridge therefore refuses non-`maxeff` focus unless `force = TRUE`, and
# labels the output as out of scope when forced.
#
# WHICH FAMILY. The bridge reconstructs the FULL enumerated <= maxk candidate
# space filtered by `n.min` -- the same space `forestsearch_main.R` hands the
# Tier-2 gate ("Full candidate space: enumerate all <= maxk combinations of the
# cut matrix Z ... so the family is identical to the space subgroup.search()
# ranked over"). It does NOT apply the consistency screen, because a plain
# argmax over the full family is exactly the primitive Guo & He analyse:
# removing the replicability screen reduces forest search's selection to their
# argmax primitive. The screen is itself a data-dependent selection step that
# their theory does not cover, so conditioning on it would misrepresent the
# comparator rather than improve it.
#
# SELF-CHECK. If the fit carries a Tier-2 gate, its `n_family` is compared with
# the reconstructed family size. A mismatch means the reconstruction diverged
# from the search's own enumeration and the result should not be used.

#' Apply the Guo & He comparator to a fitted forestsearch object
#'
#' @param fs A fitted `forestsearch` object.
#' @param r Shrinkage tuning parameter in (0, 0.5).
#' @param B Bootstrap resamples.
#' @param level One-sided level; the two-sided interval has coverage `1 - level`.
#' @param seed Optional integer for reproducibility.
#' @param orient `+1` to select the most harmful effect (the forest-search harm
#'   convention), `-1` for the most protective. Default `+1`.
#' @param adaptive Logical; if `TRUE`, choose `r` by Guo & He's Algorithm 2
#'   (requires `guohe_adaptive_r.R`).
#' @param r_grid Grid for `adaptive = TRUE`.
#' @param force Logical; proceed even when `sg_focus` is outside Guo & He's
#'   scope. The result is then labelled out of scope.
#' @param verbose Logical; report the reconstruction and self-check.
#' @param ... Passed to `guohe_algorithm3()` (e.g. `min_events`, `parallel`).
#'
#' @return The `guohe_algorithm3()` (or `guohe_adaptive_r()`) result, with a
#'   `bridge` element recording family size, self-check, and scope.
#' @export
guohe_from_forestsearch <- function(fs,
                                    r = 0.03,
                                    B = 2000L,
                                    level = 0.05,
                                    seed = NULL,
                                    orient = 1,
                                    adaptive = FALSE,
                                    r_grid = c(0.03, 0.10, 0.20, 0.30, 0.45),
                                    force = FALSE,
                                    verbose = TRUE,
                                    ...) {
  if (!exists("guohe_algorithm3", mode = "function")) {
    stop("guohe_algorithm3() not found; source 'guohe_algorithm3.R' first.")
  }
  if (!inherits(fs, "forestsearch")) {
    warning("`fs` is not of class 'forestsearch'; proceeding on structure alone.")
  }
  args <- fs$args_call_all
  if (is.null(args)) stop("`fs$args_call_all` is missing; cannot read the search configuration.")

  # ---- scope gate --------------------------------------------------------
  focus <- fs$sg_focus %||% args$sg_focus
  in_scope <- identical(focus, "maxeff")
  if (!in_scope) {
    msg <- paste0(
      "sg_focus = '", focus, "' is outside Guo & He's scope. Their correction is ",
      "built on the argmax functional; '", focus, "' selects by size within an ",
      "effect band, which is not a supremum of the effect. Use sg_focus = ",
      "'maxeff', or set force = TRUE and report the result as out of scope."
    )
    if (!isTRUE(force)) stop(msg) else warning(msg)
  }

  # ---- data and cut matrix ----------------------------------------------
  df <- fs$df.est
  if (is.null(df)) stop("`fs$df.est` is missing; refit with the estimation frame retained.")
  cuts <- fs$confounders.candidate
  if (!length(cuts)) stop("`fs$confounders.candidate` is empty; no cut columns to enumerate.")
  cuts <- intersect(cuts, names(df))
  if (!length(cuts)) stop("None of `fs$confounders.candidate` are columns of `fs$df.est`.")

  maxk  <- args$maxk  %||% 2L
  n.min <- args$n.min %||% 60L

  # column-name resolution: internal (Y/Event/Treat) or user-facing
  pick <- function(internal, user) {
    if (!is.null(internal) && internal %in% names(df)) return(internal)
    if (!is.null(user) && user %in% names(df)) return(user)
    NA_character_
  }
  treat_col <- pick("Treat", args$treat.name)
  otype <- fs$outcome_type %||% "survival"
  if (identical(otype, "survival")) {
    time_col  <- pick("Y", args$outcome.name)
    event_col <- pick("Event", args$event.name)
    y_col <- NULL
    if (anyNA(c(treat_col, time_col, event_col))) {
      stop("Could not resolve treatment/time/event columns in `fs$df.est`.")
    }
  } else {
    y_col <- pick("Y", args$outcome.name); time_col <- event_col <- NULL
    if (anyNA(c(treat_col, y_col))) {
      stop("Could not resolve treatment/outcome columns in `fs$df.est`.")
    }
  }

  # ---- reconstruct the full <= maxk family, filtered by n.min -----------
  Z <- as.data.frame(lapply(df[cuts], function(z) as.integer(z == 1)))
  names(Z) <- cuts
  L <- ncol(Z)
  memb <- list()
  for (k in seq_len(min(maxk, L))) {
    combos <- utils::combn(L, k, simplify = FALSE)
    for (cb in combos) {
      v <- Reduce(`&`, lapply(cb, function(j) Z[[j]] == 1L))
      if (sum(v) >= n.min) {
        memb[[paste(cuts[cb], collapse = " & ")]] <- as.integer(v)
      }
    }
  }
  if (!length(memb)) stop("No candidate met `n.min`; nothing to search over.")
  fam <- as.data.frame(memb, check.names = FALSE)

  # ---- self-check against the gate's own family size ---------------------
  gate_n <- tryCatch(as.integer(fs$debias_gate$n_family), error = function(e) NA_integer_)
  if (!length(gate_n)) gate_n <- NA_integer_          # no gate on this fit
  check <- if (is.na(gate_n)) "not available (no Tier-2 gate on this fit)" else
    if (identical(gate_n, ncol(fam))) sprintf("PASS (%d == %d)", ncol(fam), gate_n) else
      sprintf("MISMATCH (reconstructed %d vs gate %d)", ncol(fam), gate_n)
  if (grepl("^MISMATCH", check)) {
    warning("Family reconstruction does not match the gate's n_family: ", check,
            ". Do not use this result until reconciled.")
  }
  if (verbose) {
    cat(sprintf("Reconstructed family: %d candidates from %d cuts (maxk = %d, n.min = %d)\n",
                ncol(fam), L, maxk, n.min))
    cat(sprintf("Self-check vs gate n_family: %s\n", check))
    cat(sprintf("sg_focus = %s  (%s)\n", focus,
                if (in_scope) "in scope" else "OUT OF SCOPE -- forced"))
  }

  dat <- cbind(df[, c(treat_col, time_col, event_col, y_col), drop = FALSE], fam)

  common <- list(
    data = dat, outcome = otype, treatment = treat_col,
    candidates = names(fam), time = time_col, event = event_col, y = y_col,
    orient = orient, B = B, level = level, seed = seed
  )
  out <- if (isTRUE(adaptive)) {
    if (!exists("guohe_adaptive_r", mode = "function")) {
      stop("adaptive = TRUE requires 'guohe_adaptive_r.R'.")
    }
    do.call(guohe_adaptive_r, c(common, list(r_grid = r_grid), list(...)))
  } else {
    do.call(guohe_algorithm3, c(common, list(r = r), list(...)))
  }

  out$bridge <- list(n_family = ncol(fam), n_cuts = L, maxk = maxk, n.min = n.min,
                     sg_focus = focus, in_scope = in_scope,
                     gate_n_family = gate_n, self_check = check,
                     screen_applied = FALSE)
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
