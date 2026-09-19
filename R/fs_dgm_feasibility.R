# ---------------------------------------------------------------------------
# fs_dgm_feasibility(): can the planted region be recovered at this n?
# ---------------------------------------------------------------------------
# A design-time question, asked of a DGM before a campaign is launched: over
# repeated draws, how often is the planted region too small for the search to
# declare at all, how often does it fall under the per-arm events floor, and
# how often does the estimand fail to exist on it.
#
# This is a diagnostic.  It imposes nothing, changes no search behaviour, and
# is never consulted by forestsearch().

#' Design-time feasibility of a DGM's planted region
#'
#' Draws replicates through the DGM's own generator and reports, per sample
#' size, how often the planted harm region \eqn{Q} would be **undeclarable**
#' (too small for the search's own size test), how often it falls under the
#' per-arm events floor, and how often the analysis estimand does not exist
#' on it.
#'
#' The motivating case: at \code{n = 500} a region averaging 48.5 subjects is
#' at or below \code{n.min = 60} in 96\% of replicates, so the search is
#' structurally unable to recover a region of the planted size -- what it
#' declares is a larger overlapping region. That is a property of the design,
#' knowable before any replicate is run, and this function is how to know it.
#'
#' @section What is counted:
#' Per \code{n}, over \code{n_rep} draws:
#' \describe{
#'   \item{\code{share_undeclarable}}{share with \eqn{|Q| \le} \code{n.min}.
#'     The comparison is \code{<=}, exactly as \code{subgroup.search()} tests
#'     it (\code{nx <= n.min} rejects), so a region must *exceed*
#'     \code{n.min} to be declarable.}
#'   \item{\code{share_under_events}}{share with fewer than \code{d0.min}
#'     events in the control arm or \code{d1.min} in the treated arm. Binary
#'     and survival only -- the floor is skipped entirely for continuous and
#'     count outcomes.}
#'   \item{\code{share_nonestimable}}{share on which the estimand does not
#'     exist, under the same per-estimand condition the estimator boundary
#'     applies: OR needs all four cells \eqn{\ge 1}; RR, IRR and HR need an
#'     event in each arm; RD, IRD and MD have no condition and this share is
#'     always 0.}
#' }
#' Together with per-arm cell summaries (mean, 5th and 95th percentile, and
#' minimum of each cell across replicates).
#'
#' @section Drawing through the DGM's own generator:
#' Replicates are drawn with \code{\link{simulate_from_glm_dgm}()} -- the same
#' function the campaign templates call -- never a re-implementation, so the
#' population this reports on is the population a campaign would run on.
#'
#' This function **does not change the RNG kind**. It seeds each replicate
#' with \code{set.seed()} under whatever generator is current and restores the
#' caller's RNG state on exit. The caller's own obligation is the mirror of
#' that one: **build and calibrate the DGM before any replicate switches the
#' RNG kind.** Calibrating a DGM after a switch to \code{"L'Ecuyer-CMRG"}
#' yields a *different* super-population, and every count here would then
#' describe a population no campaign ever used.
#'
#' @param dgm An object of class \code{"glm_dgm"} -- the interface the
#'   operating-characteristics family accepts (\code{\link{fs_oc_predict}},
#'   \code{\link{fs_oc_grid}}, which validate \code{dgm$df_super} and its
#'   \code{flag_harm} column). Survival DGMs are not accepted; see Note.
#' @param n Numeric vector of sample sizes to evaluate.
#' @param n.min,d0.min,d1.min The search's own floors, passed so the report
#'   describes the campaign that will actually run. Defaults are the package
#'   defaults (60, 10, 10). Nothing is imposed: these are read, not applied.
#' @param n_rep Integer. Replicates per sample size (default 200).
#' @param tolerance Numeric. The largest undeclarable share that still counts
#'   as feasible (default 0.05).
#' @param effect_measure Character or \code{NULL}. The estimand whose
#'   existence condition is checked. \code{NULL} (default) reads
#'   \code{dgm$effect_measure}.
#' @param seed Integer or \code{NULL}. Base seed; replicate \code{i} of
#'   sample size \code{n} is drawn at \code{seed + i}.
#' @param rand_ratio Numeric. Passed to \code{simulate_from_glm_dgm()}.
#'
#' @return An object of class \code{"fs_dgm_feasibility"}: a list with
#'   \code{$table} (one row per \code{n}), \code{$feasible} (a single logical
#'   -- \code{TRUE} when every undeclarable share is at or below
#'   \code{tolerance}), and \code{$args}.
#'
#' @note **Survival DGMs are not supported, and what they would need is
#'   this:** the OC family's contract is \code{inherits(dgm, "glm_dgm")} with
#'   a \code{df_super} carrying \code{flag_harm}. \code{setup_gbsg_dgm()}
#'   produces a different object whose own generator is
#'   \code{simulate_from_dgm(dgm, n, analysis_time, cens_adjust, seed)} --
#'   a different signature, and it takes two censoring arguments that have no
#'   GLM counterpart and that change the per-arm event counts this function
#'   reports. Supporting it means a generator-dispatch layer plus passing
#'   \code{analysis_time} and \code{cens_adjust} through; that is a separate
#'   change and is deliberately not guessed at here.
#'
#' @examples
#' \donttest{
#' dgm <- generate_glm_dgm(n_super = 2000, outcome_type = "binary",
#'                         effect_measure = "OR", seed = 8316951)
#' feas <- fs_dgm_feasibility(dgm, n = c(500, 1000), n_rep = 50)
#' feas
#' feas$feasible
#' }
#'
#' @seealso \code{\link{simulate_from_glm_dgm}}, \code{\link{fs_oc_predict}}
#' @export
fs_dgm_feasibility <- function(dgm,
                               n              = c(500, 750, 1000, 2000),
                               n.min          = 60,
                               d0.min         = 10,
                               d1.min         = 10,
                               n_rep          = 200L,
                               tolerance      = 0.05,
                               effect_measure = NULL,
                               seed           = NULL,
                               rand_ratio     = 1) {

  # ---- validation: the OC family's contract, established from source -------
  if (!inherits(dgm, "glm_dgm")) {
    stop("'dgm' must be an object of class 'glm_dgm' -- the interface ",
         "fs_oc_predict() and fs_oc_grid() accept. Survival DGMs from ",
         "setup_gbsg_dgm() are not supported; see ?fs_dgm_feasibility (Note).",
         call. = FALSE)
  }
  if (!is.data.frame(dgm$df_super)) {
    stop("'dgm$df_super' must be a data frame.", call. = FALSE)
  }
  if (!"flag_harm" %in% names(dgm$df_super)) {
    stop("'dgm$df_super' has no 'flag_harm' column.", call. = FALSE)
  }
  if (!is.numeric(n) || !length(n) || anyNA(n) || any(n <= 0)) {
    stop("'n' must be a vector of positive sample sizes.", call. = FALSE)
  }
  if (!is.numeric(n_rep) || length(n_rep) != 1L || is.na(n_rep) || n_rep < 1) {
    stop("'n_rep' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) ||
      tolerance < 0 || tolerance > 1) {
    stop("'tolerance' must be a single number in [0, 1].", call. = FALSE)
  }

  outcome_type <- dgm[["outcome_type"]] %||% "binary"
  if (is.null(effect_measure)) {
    effect_measure <- dgm[["effect_measure"]] %||%
      switch(outcome_type, binary = "RD", continuous = "MD",
             count = "IRR", "RD")
  }
  effect_measure <- as.character(effect_measure)[1]
  # The events floor is skipped entirely for continuous and count outcomes
  # (subgroup.search() Status 3), so reporting it there would be misleading.
  events_apply <- identical(outcome_type, "binary")

  # ---- RNG: seed, never switch kind; restore the caller's stream -----------
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (has_seed) assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)

  n_rep <- as.integer(n_rep)
  q05 <- function(x) unname(stats::quantile(x, 0.05, names = FALSE))
  q95 <- function(x) unname(stats::quantile(x, 0.95, names = FALSE))

  rows <- vector("list", length(n))

  for (j in seq_along(n)) {
    nj <- n[[j]]
    size <- integer(n_rep)
    e0 <- ne0 <- e1 <- ne1 <- integer(n_rep)

    for (i in seq_len(n_rep)) {
      s_i <- if (is.null(seed)) NULL else as.integer(seed) + i
      if (!is.null(s_i)) set.seed(s_i)          # kind unchanged: seed only
      d <- simulate_from_glm_dgm(dgm, n = nj, rand_ratio = rand_ratio,
                                 seed = s_i)
      inQ <- d$flag_harm == 1
      size[i] <- sum(inQ)
      tt <- d$treat_sim[inQ]
      yy <- d$y_sim[inQ]
      if (events_apply) {
        e0[i]  <- sum(tt == 0L & yy == 1L, na.rm = TRUE)
        ne0[i] <- sum(tt == 0L & yy == 0L, na.rm = TRUE)
        e1[i]  <- sum(tt == 1L & yy == 1L, na.rm = TRUE)
        ne1[i] <- sum(tt == 1L & yy == 0L, na.rm = TRUE)
      } else {
        # subjects per arm; no events floor and no existence condition
        e0[i] <- sum(tt == 0L); e1[i] <- sum(tt == 1L)
        ne0[i] <- NA_integer_;  ne1[i] <- NA_integer_
      }
    }

    # Undeclarable: the search's own test is `nx <= n.min` (status 4), so a
    # region must EXCEED n.min.  Strict, as the search is.
    undecl <- size <= n.min
    under_ev <- if (events_apply) (e0 < d0.min | e1 < d1.min)
                else rep(NA, n_rep)
    nonest <- if (events_apply) {
      vapply(seq_len(n_rep), function(i) {
        !is.null(.fs_existence_reason(
          effect_measure,
          c(e0 = e0[i], ne0 = ne0[i], e1 = e1[i], ne1 = ne1[i])))
      }, logical(1))
    } else rep(FALSE, n_rep)

    rows[[j]] <- data.frame(
      n                  = nj,
      n_rep              = n_rep,
      size_mean          = mean(size),
      size_q05           = q05(size),
      size_q95           = q95(size),
      size_min           = min(size),
      share_undeclarable = mean(undecl),
      share_under_events = if (events_apply) mean(under_ev) else NA_real_,
      share_nonestimable = mean(nonest),
      e0_mean = mean(e0), e0_min = min(e0),
      e1_mean = mean(e1), e1_min = min(e1),
      ne0_mean = if (events_apply) mean(ne0) else NA_real_,
      ne0_min  = if (events_apply) min(ne0)  else NA_integer_,
      ne1_mean = if (events_apply) mean(ne1) else NA_real_,
      ne1_min  = if (events_apply) min(ne1)  else NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  tab <- do.call(rbind, rows)

  structure(
    list(
      table    = tab,
      feasible = all(tab$share_undeclarable <= tolerance),
      args     = list(n = n, n.min = n.min, d0.min = d0.min, d1.min = d1.min,
                      n_rep = n_rep, tolerance = tolerance,
                      outcome_type = outcome_type,
                      effect_measure = effect_measure,
                      events_apply = events_apply, seed = seed)
    ),
    class = "fs_dgm_feasibility"
  )
}

#' @param x An \code{fs_dgm_feasibility} object.
#' @param ... Ignored.
#' @rdname fs_dgm_feasibility
#' @export
print.fs_dgm_feasibility <- function(x, ...) {
  a <- x$args
  cat("DGM feasibility -- can the planted region be declared at this n?\n")
  cat(sprintf("  outcome_type = %s; effect_measure = %s; %d replicates per n\n",
              a$outcome_type, a$effect_measure, a$n_rep))
  cat(sprintf("  floors read (not imposed): n.min = %s, d0.min = %s, d1.min = %s\n",
              format(a$n.min), format(a$d0.min), format(a$d1.min)))
  if (!a$events_apply) {
    cat(sprintf(paste0("  NOTE: the per-arm events floor is skipped entirely ",
                       "for outcome_type = \"%s\";\n        e0/e1 are ",
                       "SUBJECTS per arm and no existence condition applies.\n"),
                a$outcome_type))
  }
  cat("\n")
  tb <- x$table
  cat(sprintf("  %6s %9s %8s %8s %14s %14s %14s\n",
              "n", "size_mean", "size_q05", "size_min",
              "undeclarable", "under_events", "non-estimable"))
  for (i in seq_len(nrow(tb))) {
    cat(sprintf("  %6s %9.1f %8.1f %8d %13.1f%% %13s %13.1f%%\n",
                format(tb$n[i]), tb$size_mean[i], tb$size_q05[i],
                as.integer(tb$size_min[i]),
                100 * tb$share_undeclarable[i],
                if (is.na(tb$share_under_events[i])) "--"
                else sprintf("%.1f%%", 100 * tb$share_under_events[i]),
                100 * tb$share_nonestimable[i]))
  }
  cat("\n")
  bad <- tb$n[tb$share_undeclarable > a$tolerance]
  if (x$feasible) {
    cat(sprintf("  FEASIBLE: every undeclarable share is at or below the %.0f%% tolerance.\n",
                100 * a$tolerance))
  } else {
    cat(sprintf("  NOT FEASIBLE at n = %s: the planted region is at or below n.min = %s\n",
                paste(format(bad), collapse = ", "), format(a$n.min)))
    cat(sprintf("    more often than the %.0f%% tolerance allows. At those n the search\n",
                100 * a$tolerance))
    cat("    cannot declare a region of the planted size; what it declares is a\n")
    cat("    LARGER overlapping region, and recovery metrics read accordingly.\n")
  }
  invisible(x)
}
