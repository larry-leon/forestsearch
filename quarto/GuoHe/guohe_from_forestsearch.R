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
#                   and the ONLY configuration Guo & He covers.
#   * "hr"/"eff"  -- lexicographic (-Pcons, -hr): consistency rate first, effect
#                   only as a tiebreaker; and under early stopping the search
#                   returns the FIRST passer in effect order, not the effect
#                   maximiser. Not a supremum of the effect -- out of scope.
#   * "effMaxSG" -- max SIZE within an effect band (the neighbourhood rule), and
#                   "pareto" / "both" bands. These are NOT suprema of the
#                   effect, so none of the four pillars above transfers.
#                   Extending Guo & He to them would be new theory, not new code.
#
# This bridge therefore accepts ONLY `sg_focus = "maxeff"`; every other focus
# requires `force = TRUE` and is labelled out of scope when forced.
#
# PRECONDITION. Even for "maxeff" the bridge re-derives the argmax over its OWN
# reconstructed family and asserts it equals `fs$sg.harm`. The two families
# differ by construction -- the bridge filters on `n.min` only, whereas the
# search now evaluates a smaller space (per-arm minima, rmin = 0 containment
# pruning, exact-membership dedup). They agree on GBSG, but if they diverge the
# bridge would de-bias a subgroup the analysis does not report, so it STOPs.
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
  # Guo & He cover the argmax-of-effect primitive, which is forestsearch()'s
  # user-facing sg_focus = "maxeff": which.max(beta) among passers, with no
  # auxiliary selection condition.  "hr"/"eff" are demoted to force-only: they
  # rank lexicographically by (-Pcons, -hr) -- consistency first, effect only as
  # a tiebreaker -- and under early stopping return the first passer in effect
  # order rather than the maximiser, so they are NOT an effect argmax.
  in_scope <- identical(focus, "maxeff")
  if (!in_scope) {
    msg <- paste0(
      "sg_focus = '", focus, "' is outside Guo & He's scope. Their correction is ",
      "built on the argmax-of-effect functional (which.max(beta) among passers). ",
      "'hr'/'eff' select lexicographically by (-Pcons, -hr) -- consistency first, ",
      "effect only as a tiebreaker -- and under early stopping return the first ",
      "passer in effect order, not the maximiser; other foci select by size ",
      "within an effect band. None is a supremum of the effect. Use sg_focus = ",
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
  # Enumerate with the SEARCH'S OWN helpers over the dummy-expanded cut matrix
  # Z -- i.e. BOTH directions of each cut (x.0 and x.1), exactly as
  # forestsearch_main.R builds the space it hands the Tier-2 gate
  # (Z <- as.matrix(dummy(df[, conf.screen]))).  Coercing the single-direction
  # `cuts` with combn (the previous approach) omitted every complement column
  # and could never reproduce the gate's family.  `as.factor` before
  # dummy_encode guarantees each 0/1 cut expands to two indicators even if the
  # estimation frame stored it numeric.
  Zdf <- forestsearch::dummy_encode(
    as.data.frame(lapply(df[cuts], as.factor), check.names = FALSE))
  Z   <- as.matrix(Zdf)
  colnames(Z) <- names(Zdf)
  L   <- ncol(Z)
  combo <- forestsearch:::generate_combination_indices(L, maxk)
  tot   <- forestsearch:::calculate_max_combinations(L, maxk)
  memb  <- list()
  for (kk in seq_len(tot)) {
    covs.in <- forestsearch::get_covs_in(
      kk, maxk, L,
      combo$counts_1, combo$indices_1,
      combo$counts_2, combo$indices_2,
      combo$counts_3, combo$indices_3)
    k_sel <- sum(covs.in)
    if (k_sel < 1L || k_sel > maxk) next
    v <- as.integer(forestsearch::get_subgroup_membership(Z, covs.in))
    if (sum(v) >= n.min) {
      memb[[paste(colnames(Z)[covs.in == 1], collapse = " & ")]] <- v
    }
  }
  if (!length(memb)) stop("No candidate met `n.min`; nothing to search over.")
  fam <- as.data.frame(memb, check.names = FALSE)

  # ---- self-check against the gate's own family size ---------------------
  mr_n <- tryCatch(as.integer(fs$mr_inference$n_family), error = function(e) NA_integer_)
  if (!length(mr_n)) mr_n <- NA_integer_          # no gate on this fit
  check <- if (is.na(mr_n)) "not available (no Tier-2 gate on this fit)" else
    if (identical(mr_n, ncol(fam))) sprintf("PASS (%d == %d)", ncol(fam), mr_n) else
      sprintf("MISMATCH (reconstructed %d vs gate %d)", ncol(fam), mr_n)
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

  # ---- precondition: the bridge's OWN argmax must equal fs$sg.harm -------
  # The reconstructed family (n.min only) and the space forestsearch actually
  # evaluated (per-arm minima, rmin = 0 containment pruning, membership dedup)
  # differ by construction.  They agree on GBSG, but nothing guarantees it
  # elsewhere.  If the bridge's effect-argmax over ITS family is not the
  # subgroup the analysis reports, Guo & He would de-bias a subgroup that is
  # never selected -- so STOP rather than mislabel the output.  The argmax is
  # computed with the SAME within-subgroup estimator guohe_algorithm3() uses
  # (.g3_coef), oriented by `orient`.
  sel_sg <- fs$sg.harm
  if (is.null(sel_sg) || !length(sel_sg)) {
    stop("`fs$sg.harm` is empty; the fit selected no subgroup to de-bias.")
  }
  .dots <- list(...)
  .me   <- if (!is.null(.dots$min_events)) .dots$min_events else 5L
  .adjc <- .dots$adjust_covariates                       # NULL unless supplied
  beta_hat <- vapply(names(fam), function(nm) {
    sub <- df[fam[[nm]] == 1L, , drop = FALSE]
    orient * .g3_coef(sub, otype, treat_col, time_col, event_col, y_col, .me, .adjc)
  }, numeric(1))
  if (!any(is.finite(beta_hat))) {
    stop("No estimable candidate in the reconstructed family; cannot check the argmax.")
  }
  argmax_label  <- names(fam)[which.max(replace(beta_hat, !is.finite(beta_hat), -Inf))]
  argmax_memb   <- as.integer(fam[[argmax_label]] == 1L)
  sg_harm_label <- paste(sel_sg, collapse = " & ")

  # fs$sg.harm membership over the SAME estimation frame (row-aligned with fam;
  # both are built from fs$df.est).  Compare memberships, not labels, so a
  # different but equivalent cut-encoding of the same subgroup is not a spurious
  # mismatch.
  sel_memb <- fs$grp.consistency$sg.harm.id
  if (is.null(sel_memb) || length(sel_memb) != nrow(df)) {
    stop("Cannot recover fs$sg.harm membership (grp.consistency$sg.harm.id) ",
         "aligned to fs$df.est; precondition unverifiable.")
  }
  sel_memb        <- as.integer(sel_memb == 1L)
  precondition_ok <- identical(argmax_memb, sel_memb)
  if (verbose) {
    cat(sprintf("Precondition (bridge argmax == fs$sg.harm): %s\n",
                if (precondition_ok) "PASS" else "FAIL"))
    cat(sprintf("  bridge argmax : %s\n", argmax_label))
    cat(sprintf("  fs$sg.harm    : %s\n", sg_harm_label))
  }
  if (!precondition_ok) {
    stop(sprintf(paste0(
      "Precondition FAILED: the bridge's argmax over its reconstructed family ",
      "('%s') is not the subgroup forestsearch reports (fs$sg.harm = '%s'). Guo & ",
      "He would de-bias a subgroup the analysis does not select. The reconstructed ",
      "family (n.min only) and the evaluated space differ by construction; ",
      "reconcile the two before using this comparator."),
      argmax_label, sg_harm_label))
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
                     mr_n_family = mr_n, self_check = check,
                     screen_applied = FALSE,
                     precondition_ok = precondition_ok,
                     argmax_label = argmax_label,
                     sg_harm_label = sg_harm_label)
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
