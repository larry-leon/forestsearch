# =============================================================================
# fs_oc_grid() / fs_oc_invert() -- sweep (c1, c2) on one draw set per n, and
#                                  invert the declaration rate for a target
# -----------------------------------------------------------------------------
# At fixed n the family, its means, ses and covariance are fixed, and BOTH
# gates are thresholds on draws that do not depend on c1 or c2:
#     resample:  (Bhat >= c1) & (Bhat - c2 >= z_p * se_g)
#     split:     (W1 + W2 >= 2*c1) & (W1 >= c2) & (W2 >= c2)
# so one draw set serves the whole (c1, c2) grid and any inversion at that n.
# Only n forces re-enumeration and re-drawing (the size floor is n.min / n).
#
# The draw, gate and functional code is fs_oc_predict()'s own (.fs_oc_draw,
# .fs_oc_gate, .fs_oc_functionals, .fs_oc_result in R/fs_oc_predict.R), so a
# single grid point evaluated in one block is identical() to fs_oc_predict().
#
# Blocking: with block < draws the draws are generated in row blocks and the
# functionals accumulated (declaration count, per-member pass and selection
# counts, winner-noise sum).  Blocked draws lay the SAME RNG stream out
# differently across the Rdraw x M matrix, so blocked results agree with the
# one-block results to Monte-Carlo precision, not bit-for-bit; and the
# accumulated proportions are count / draws rather than colMeans()/mean(),
# which can differ in the last bit.  The one-block path is exact.
# =============================================================================


#' Sweep the declaration thresholds on one draw set per trial size
#'
#' Evaluates every \code{fs_oc_predict()} quantity over the full crossing of
#' \code{n}, \code{c1} and \code{c2}, enumerating the family and drawing the
#' candidates' fitted effects \strong{once per \code{n}} (and per gate) and
#' sweeping the thresholds against those draws.  Thresholds enter only the
#' gate, so the sweep costs arithmetic, not draws.
#'
#' @details
#' For each \code{n} the family is enumerated with
#' \code{\link{fs_oc_family_enumerate}} (or the supplied \code{family} is used
#' when its \code{n} matches), the draw set is generated under \code{seed}
#' (the same seed for every \code{n}: common random numbers across trial
#' sizes), and every \code{(c1, c2)} is evaluated on it.  With
#' \code{block = Inf} the whole draw set is held and the computation is exactly
#' \code{\link{fs_oc_predict}}'s -- one grid point is \code{identical()} to a
#' call of \code{fs_oc_predict()} at the same settings and seed.  With a finite
#' \code{block} the draws are generated in row blocks of that size and the
#' quantities accumulated, so memory is O(block x M) instead of O(draws x M);
#' blocked results agree with the one-block results to Monte-Carlo precision
#' (the RNG stream is laid out differently and proportions are formed as
#' count / draws), not bit-for-bit.
#'
#' @inheritParams fs_oc_predict
#' @param n Numeric vector of trial sizes.
#' @param c1,c2 Numeric vectors of screening and consistency floors; the grid
#'   is their full crossing with \code{n}.
#' @param family \code{NULL} (enumerate per \code{n}) or an \code{fs_oc_family}
#'   object, accepted only when its \code{n} equals the single \code{n}
#'   requested -- a family built at one \code{n} has the wrong prevalence
#'   floor for another.
#' @param consistency_method \code{"resample"}, \code{"split"}, or both.
#' @param block Integer or \code{Inf}.  Draw-block size; \code{Inf} means one
#'   block (exact).  Default \code{5e4}.
#' @param verbose Logical.  Report per-\code{n} timings.
#'
#' @return An object of class \code{c("fs_oc_grid", "list")}:
#'   \describe{
#'     \item{\code{table}}{Data frame, one row per \code{(n, gate, c1, c2)}:
#'       \code{n, consistency_method, pconsistency, c1, c2, M, draws, block,
#'       seed, det_rate, det_rate_se, EnH, Esens, Espec, Eppv, Enpv, EbetaH,
#'       Enaive_bias, mass_below}.}
#'     \item{\code{results}}{List, parallel to the rows of \code{table}, of the
#'       full per-point objects (class \code{fs_oc_predict}, with \code{P1},
#'       \code{p_sel}, \code{sel_c} and their MC SEs).}
#'     \item{\code{families}}{Per \code{n}: \code{M}, the floor, the stage
#'       counts and the family.}
#'     \item{\code{timing}}{Per \code{n} and gate: seconds for enumeration,
#'       drawing and sweeping.}
#'     \item{\code{settings}}{The call's settings.}
#'   }
#'
#' @seealso \code{\link{fs_oc_predict}}, \code{\link{fs_oc_invert}}
#'
#' @examples
#' piQ <- 0.34
#' fam <- structure(list(
#'   lab = c("Q", "P", "D"), Pg = c(piQ, 0.45, 0.31),
#'   PQg = c(1, 0.28 / 0.45, 1), sens_g = c(1, 0.28 / piQ, 0.31 / piQ),
#'   spec_g = c(1, 1 - 0.17 / (1 - piQ), 1),
#'   ovl = matrix(c(piQ, 0.28, 0.31, 0.28, 0.45, 0.28 * 0.31 / piQ,
#'                  0.31, 0.28 * 0.31 / piQ, 0.31), 3, 3),
#'   M = 3L, PQ = piQ, n = 500), class = c("fs_oc_family", "list"))
#' fam$beta_g <- 26 + 14 * fam$PQg
#' fam$se_g   <- 13.7 * sqrt(2) * sqrt(piQ / fam$Pg)
#' g <- fs_oc_grid(family = fam, n = 500, c1 = c(20, 30, 40), c2 = 10,
#'                 consistency_method = "resample", draws = 2e4, seed = 1)
#' g
#'
#' @export
fs_oc_grid <- function(dgm = NULL, forestsearch_args = list(), n, c1, c2,
                       family = NULL,
                       consistency_method = c("resample", "split"),
                       pconsistency = NULL, draws = 2e5, block = 5e4,
                       seed = NULL, verbose = FALSE, ...) {

  methods <- match.arg(consistency_method, several.ok = TRUE)
  if (missing(n) || !is.numeric(n) || length(n) < 1L || anyNA(n) || any(n <= 0)) {
    stop("'n' must be a vector of positive numbers.", call. = FALSE)
  }
  if (missing(c1) || is.null(c1)) c1 <- forestsearch_args$effect.threshold
  if (missing(c2) || is.null(c2)) c2 <- forestsearch_args$consistency.threshold
  if (is.null(c1) || !is.numeric(c1) || anyNA(c1)) {
    stop("'c1' is required (a numeric vector).", call. = FALSE)
  }
  if (is.null(c2) || !is.numeric(c2) || anyNA(c2)) {
    stop("'c2' is required (a numeric vector).", call. = FALSE)
  }
  if (!is.numeric(draws) || length(draws) != 1L || is.na(draws) || draws < 1) {
    stop("'draws' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(block) || length(block) != 1L || is.na(block) || block < 1) {
    stop("'block' must be a single positive number or Inf.", call. = FALSE)
  }
  pcons <- .fs_oc_resolve_pcons(pconsistency, forestsearch_args)

  if (!is.null(family)) {
    if (!inherits(family, "fs_oc_family")) {
      stop("'family' must be an object of class 'fs_oc_family'.", call. = FALSE)
    }
    if (length(n) != 1L || is.null(family$n) || !isTRUE(all.equal(family$n, n))) {
      stop(sprintf(paste0(
        "the supplied family was built at n = %s but the grid asks for n = %s; ",
        "a family carries the prevalence floor n.min/n of its own n, so pass ",
        "family = NULL to enumerate per n, or match n."),
        if (is.null(family$n)) "<unrecorded>" else format(family$n),
        paste(format(n), collapse = ", ")), call. = FALSE)
    }
  } else if (is.null(dgm)) {
    stop("either 'family' or 'dgm' (with 'forestsearch_args') is required.",
         call. = FALSE)
  }

  results <- list(); rows <- list(); families <- list(); timing <- list()
  for (nn in n) {
    key <- as.character(nn)
    t0 <- proc.time()[["elapsed"]]
    fam <- if (!is.null(family)) family else
      fs_oc_family_enumerate(dgm, forestsearch_args, nn, ...)
    .fs_oc_check_family(fam)
    t_enum <- proc.time()[["elapsed"]] - t0
    families[[key]] <- list(n = nn, M = fam$M, floor = fam$args_used$n.min / nn,
                            counts = fam$counts, family = fam)
    for (gate in methods) {
      pc <- if (identical(gate, "resample")) pcons else NA_real_
      t0 <- proc.time()[["elapsed"]]
      ev <- .fs_oc_sweep(fam, nn, c1, c2, gate, pc, draws, block, seed)
      timing[[paste(key, gate)]] <- data.frame(
        n = nn, consistency_method = gate, M = fam$M,
        enumerate_secs = t_enum, draw_secs = ev$draw_secs,
        sweep_secs = ev$sweep_secs, stringsAsFactors = FALSE)
      if (verbose) {
        cat(sprintf("fs_oc_grid: n = %s  gate = %-8s  M = %d  enumerate %.1fs  draw %.1fs  sweep %.1fs (%d points)\n",
                    format(nn), gate, fam$M, t_enum, ev$draw_secs, ev$sweep_secs,
                    length(ev$results)))
      }
      for (r in ev$results) {
        results[[length(results) + 1L]] <- r
        rows[[length(rows) + 1L]] <- data.frame(
          n = nn, consistency_method = gate, pconsistency = pc,
          c1 = r$settings$c1, c2 = r$settings$c2, M = r$M,
          draws = draws, block = block,
          seed = if (is.null(seed)) NA_real_ else seed,
          det_rate = r$det_rate, det_rate_se = r$det_rate_se,
          EnH = r$EnH, Esens = r$Esens, Espec = r$Espec, Eppv = r$Eppv,
          Enpv = r$Enpv, EbetaH = r$EbetaH, Enaive_bias = r$Enaive_bias,
          mass_below = r$mass_below, stringsAsFactors = FALSE)
      }
    }
  }
  out <- list(
    table    = do.call(rbind, rows),
    results  = results,
    families = families,
    timing   = do.call(rbind, timing),
    settings = list(n = n, c1 = c1, c2 = c2, consistency_method = methods,
                    pconsistency = pcons, draws = draws, block = block,
                    seed = seed))
  rownames(out$table) <- NULL; rownames(out$timing) <- NULL
  class(out) <- c("fs_oc_grid", "list")
  out
}


#' Invert the family declaration rate for a target
#'
#' Finds the threshold -- \code{c1} at fixed \code{c2}, or \code{c2} at fixed
#' \code{c1} -- at which the family declaration rate equals \code{target}, on
#' one fixed draw set: "the \code{c1} giving 80\% power", the same quantity
#' read as a type-I error under a null DGM.
#'
#' @details
#' On fixed draws the rate is a monotone non-increasing step function of the
#' threshold.  The function first evaluates the \strong{ceiling}: the rate as
#' the solved-for threshold goes to \code{-Inf}, where only the other threshold
#' binds (\code{c2 + z_p * se_g} for \code{"resample"}; both halves clearing
#' \code{c2} for \code{"split"}; or the \code{c1} screen when solving for
#' \code{c2}).  A \code{target} above the ceiling is \strong{unattainable}: the
#' result is \code{NA} with the ceiling and the binding threshold named -- no
#' extrapolation.  Otherwise a bracket is established from the draws' range and
#' bisected to \code{tol} on the threshold scale; the returned value is the
#' largest threshold whose rate is still at least \code{target} (the step the
#' target falls on), reported with the achieved rate and its MC SE so the
#' resolution the draw count supports is visible.
#'
#' @inheritParams fs_oc_grid
#' @param n Single trial size.
#' @param target Numeric in (0, 1).  Target family declaration rate.
#' @param solve_for \code{"c1"} (at fixed \code{c2}) or \code{"c2"} (at fixed
#'   \code{c1}).
#' @param c1,c2 The fixed threshold (the one not solved for; the solved-for
#'   argument is ignored).  Defaults from \code{forestsearch_args}.
#' @param consistency_method One gate.
#' @param tol Bisection tolerance on the threshold scale.
#' @param ... Passed to \code{\link{fs_oc_family_enumerate}} when
#'   \code{family = NULL}.
#'
#' @return An object of class \code{c("fs_oc_invert", "list")}: \code{value}
#'   (the threshold, or \code{NA}), \code{achieved} (rate at \code{value}),
#'   \code{achieved_se}, \code{target}, \code{ceiling}, \code{ceiling_se},
#'   \code{binding} (which threshold sets the ceiling), \code{attainable},
#'   \code{solve_for}, \code{fixed} (the other threshold), \code{iterations},
#'   \code{bracket}, \code{settings}.
#'
#' @seealso \code{\link{fs_oc_grid}}
#'
#' @examples
#' piQ <- 0.34
#' fam <- structure(list(
#'   lab = c("Q", "P", "D"), Pg = c(piQ, 0.45, 0.31),
#'   PQg = c(1, 0.28 / 0.45, 1), sens_g = c(1, 0.28 / piQ, 0.31 / piQ),
#'   spec_g = c(1, 1 - 0.17 / (1 - piQ), 1),
#'   ovl = matrix(c(piQ, 0.28, 0.31, 0.28, 0.45, 0.28 * 0.31 / piQ,
#'                  0.31, 0.28 * 0.31 / piQ, 0.31), 3, 3),
#'   M = 3L, PQ = piQ, n = 500), class = c("fs_oc_family", "list"))
#' fam$beta_g <- 26 + 14 * fam$PQg
#' fam$se_g   <- 13.7 * sqrt(2) * sqrt(piQ / fam$Pg)
#' fs_oc_invert(family = fam, n = 500, target = 0.5, solve_for = "c1", c2 = 10,
#'              consistency_method = "resample", draws = 2e4, seed = 1)
#'
#' @export
fs_oc_invert <- function(dgm = NULL, forestsearch_args = list(), n, target,
                         solve_for = c("c1", "c2"), c1 = NULL, c2 = NULL,
                         family = NULL,
                         consistency_method = c("resample", "split"),
                         pconsistency = NULL, draws = 2e5, seed = NULL,
                         tol = 1e-3, ...) {

  solve_for <- match.arg(solve_for)
  gate <- match.arg(consistency_method)
  if (missing(n) || !is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0) {
    stop("'n' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(target) || length(target) != 1L || is.na(target) ||
      target <= 0 || target >= 1) {
    stop("'target' must be a single number in (0, 1).", call. = FALSE)
  }
  if (is.null(c1)) c1 <- forestsearch_args$effect.threshold
  if (is.null(c2)) c2 <- forestsearch_args$consistency.threshold
  fixed <- if (solve_for == "c1") c2 else c1
  if (is.null(fixed) || !is.numeric(fixed) || length(fixed) != 1L || is.na(fixed)) {
    stop(sprintf("the fixed threshold '%s' is required.",
                 if (solve_for == "c1") "c2" else "c1"), call. = FALSE)
  }
  pc <- if (identical(gate, "resample")) .fs_oc_resolve_pcons(pconsistency, forestsearch_args) else NA_real_

  if (!is.null(family)) {
    if (!inherits(family, "fs_oc_family")) {
      stop("'family' must be an object of class 'fs_oc_family'.", call. = FALSE)
    }
    if (is.null(family$n) || !isTRUE(all.equal(family$n, n))) {
      stop(sprintf("the supplied family was built at n = %s, not n = %s.",
                   if (is.null(family$n)) "<unrecorded>" else format(family$n),
                   format(n)), call. = FALSE)
    }
  } else {
    if (is.null(dgm)) stop("either 'family' or 'dgm' is required.", call. = FALSE)
    family <- fs_oc_family_enumerate(dgm, forestsearch_args, n, ...)
  }
  .fs_oc_check_family(family)

  # ---- one draw set ----------------------------------------------------------
  Sg <- (family$ovl / sqrt(outer(family$Pg, family$Pg))) * outer(family$se_g, family$se_g)
  if (!is.null(seed)) set.seed(seed)
  dr <- .fs_oc_draw(Sg, family$beta_g, family$M, draws, gate)

  rate_at <- function(x) {
    pass <- if (solve_for == "c1") .fs_oc_gate(dr, x, fixed, family$se_g, pc, gate)
            else                   .fs_oc_gate(dr, fixed, x, family$se_g, pc, gate)
    mean(rowSums(pass) > 0)
  }
  se_of <- function(p) .fs_oc_mc_se_prop(p, draws)

  # ---- ceiling: the solved-for threshold at -Inf ------------------------------
  ceiling <- rate_at(-Inf)
  binding <- if (solve_for == "c1") {
    if (identical(gate, "resample")) "c2 + z_p * se_g (resample consistency screen)"
    else "both halves >= c2 (split consistency screen)"
  } else {
    if (identical(gate, "resample")) "Bhat >= c1 (effect screen)"
    else "W1 + W2 >= 2 * c1 (effect screen)"
  }
  settings <- list(n = n, solve_for = solve_for, fixed = fixed, target = target,
                   consistency_method = gate, pconsistency = pc, draws = draws,
                   seed = seed, tol = tol, M = family$M)
  if (target > ceiling) {
    out <- list(value = NA_real_, achieved = NA_real_, achieved_se = NA_real_,
                target = target, ceiling = ceiling, ceiling_se = se_of(ceiling),
                binding = binding, attainable = FALSE, solve_for = solve_for,
                fixed = fixed, iterations = 0L, bracket = c(NA_real_, NA_real_),
                settings = settings)
    class(out) <- c("fs_oc_invert", "list")
    return(out)
  }

  # ---- bracket from the draws' range, then bisect ----------------------------
  # Rate is non-increasing in x.  lo: rate >= target; hi: rate < target.
  lo <- min(dr$Bhat) - 1; hi <- max(dr$Bhat) + 1
  if (identical(gate, "split")) { lo <- min(dr$W1, dr$W2) - 1; hi <- max(dr$W1, dr$W2) + 1 }
  if (solve_for == "c2" && identical(gate, "resample")) {
    zp <- stats::qnorm((1 + pc) / 2); lo <- lo - zp * max(family$se_g)
  }
  r_lo <- rate_at(lo); r_hi <- rate_at(hi)
  if (r_lo < target) { lo <- -Inf; r_lo <- ceiling }
  if (r_hi >= target) stop("bisection bracket failed: rate at the upper end still >= target.", call. = FALSE)
  it <- 0L
  while (is.finite(lo) && is.finite(hi) && (hi - lo) > tol && it < 200L) {
    mid <- (lo + hi) / 2; r_mid <- rate_at(mid); it <- it + 1L
    if (r_mid >= target) { lo <- mid; r_lo <- r_mid } else { hi <- mid; r_hi <- r_mid }
  }
  if (!is.finite(lo)) {
    # target equals the ceiling exactly at -Inf and nowhere finite: report the
    # lowest finite bracket end that still attains it
    lo <- hi; r_lo <- rate_at(lo)
    while (r_lo < target && it < 400L) { lo <- lo - 1; r_lo <- rate_at(lo); it <- it + 1L }
  }
  out <- list(value = lo, achieved = r_lo, achieved_se = se_of(r_lo),
              target = target, ceiling = ceiling, ceiling_se = se_of(ceiling),
              binding = binding, attainable = TRUE, solve_for = solve_for,
              fixed = fixed, iterations = it, bracket = c(lo, hi),
              next_step_rate = r_hi, settings = settings)
  class(out) <- c("fs_oc_invert", "list")
  out
}


# -----------------------------------------------------------------------------
# methods
# -----------------------------------------------------------------------------

#' @export
print.fs_oc_grid <- function(x, ...) {
  s <- x$settings
  cat("Threshold sweep on one draw set per n (fs_oc_grid)\n")
  cat(sprintf("  n: %s | gates: %s | c1: %d values | c2: %d values | draws = %s, block = %s%s\n",
              paste(format(s$n), collapse = ", "),
              paste(s$consistency_method, collapse = ", "),
              length(s$c1), length(s$c2), format(s$draws, big.mark = ","),
              format(s$block), if (is.null(s$seed)) "" else sprintf(", seed = %s", format(s$seed))))
  for (k in names(x$families)) {
    f <- x$families[[k]]
    cat(sprintf("  n = %s: M = %d (floor Pg >= %.4f)\n", k, f$M, f$floor))
  }
  cat(sprintf("  %d grid rows; table columns: %s\n", nrow(x$table),
              paste(names(x$table), collapse = ", ")))
  invisible(x)
}

#' @export
summary.fs_oc_grid <- function(object, ...) {
  tb <- object$table
  out <- tb[, c("n", "consistency_method", "c1", "c2", "M", "det_rate",
                "det_rate_se", "EnH", "Esens", "Espec", "Eppv", "Enpv",
                "EbetaH", "Enaive_bias", "mass_below")]
  class(out) <- c("summary.fs_oc_grid", "data.frame")
  out
}

#' @export
print.summary.fs_oc_grid <- function(x, digits = 3, ...) {
  print.data.frame(x, digits = digits, row.names = FALSE, ...)
  invisible(x)
}

#' @export
print.fs_oc_invert <- function(x, ...) {
  s <- x$settings
  cat("Declaration-rate inversion (fs_oc_invert)\n")
  cat(sprintf("  solve for %s at fixed %s = %s; gate = %s; n = %s; M = %d; draws = %s\n",
              x$solve_for, if (x$solve_for == "c1") "c2" else "c1", format(x$fixed),
              s$consistency_method, format(s$n), s$M, format(s$draws, big.mark = ",")))
  cat(sprintf("  ceiling (threshold -> -Inf): %.4f (MC SE %.4f), set by %s\n",
              x$ceiling, x$ceiling_se, x$binding))
  if (!x$attainable) {
    cat(sprintf("  target %.3f: NA -- infeasible (ceiling %.4f)\n", x$target, x$ceiling))
  } else {
    cat(sprintf("  target %.3f: %s = %.4f, achieved rate %.4f (MC SE %.4f), next step %.4f; %d bisections\n",
                x$target, x$solve_for, x$value, x$achieved, x$achieved_se,
                x$next_step_rate, x$iterations))
  }
  invisible(x)
}


# -----------------------------------------------------------------------------
# internals
# -----------------------------------------------------------------------------

.fs_oc_resolve_pcons <- function(pconsistency, forestsearch_args) {
  if (is.null(pconsistency)) pconsistency <- forestsearch_args$pconsistency.threshold
  if (is.null(pconsistency)) pconsistency <- eval(formals(forestsearch)$pconsistency.threshold)
  if (!is.numeric(pconsistency) || length(pconsistency) != 1L ||
      is.na(pconsistency) || pconsistency <= 0 || pconsistency >= 1) {
    stop("'pconsistency' must be a single number in (0, 1).", call. = FALSE)
  }
  pconsistency
}

# One family, one gate, one draw set (possibly blocked), every (c1, c2).
.fs_oc_sweep <- function(fam, n, c1, c2, gate, pc, draws, block, seed) {
  Sg <- (fam$ovl / sqrt(outer(fam$Pg, fam$Pg))) * outer(fam$se_g, fam$se_g)
  grid <- expand.grid(c1 = c1, c2 = c2, KEEP.OUT.ATTRS = FALSE)
  if (!is.null(seed)) set.seed(seed)
  results <- vector("list", nrow(grid))

  if (!is.finite(block) || block >= draws) {
    # ---- one block: exactly fs_oc_predict()'s computation ----------------
    t0 <- proc.time()[["elapsed"]]
    dr <- .fs_oc_draw(Sg, fam$beta_g, fam$M, draws, gate)
    draw_secs <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(nrow(grid))) {
      pass <- .fs_oc_gate(dr, grid$c1[i], grid$c2[i], fam$se_g, pc, gate)
      fx   <- .fs_oc_functionals(pass, dr$Bhat, fam, n, grid$c1[i], draws)
      results[[i]] <- .fs_oc_result(fx, fam, n, grid$c1[i], grid$c2[i], gate,
                                    pc, draws, seed)
    }
    return(list(results = results, draw_secs = draw_secs,
                sweep_secs = proc.time()[["elapsed"]] - t0))
  }

  # ---- blocked: accumulate ---------------------------------------------------
  M <- fam$M; G <- nrow(grid)
  det_n  <- numeric(G)               # declaring draws
  pass_n <- matrix(0, G, M)          # per-member declarations
  sel_n  <- matrix(0, G, M)          # per-member wins
  noise_s <- numeric(G)              # sum of winner noise over declaring draws
  draw_secs <- 0; sweep_secs <- 0
  done <- 0
  while (done < draws) {
    b <- min(block, draws - done)
    t0 <- proc.time()[["elapsed"]]
    dr <- .fs_oc_draw(Sg, fam$beta_g, M, b, gate)
    draw_secs <- draw_secs + proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    noise <- dr$Bhat - rep(fam$beta_g, each = b)
    for (i in seq_len(G)) {
      pass <- .fs_oc_gate(dr, grid$c1[i], grid$c2[i], fam$se_g, pc, gate)
      det  <- rowSums(pass) > 0
      Bmask <- dr$Bhat; Bmask[!pass] <- -Inf
      winner <- max.col(Bmask, ties.method = "first")
      det_n[i]    <- det_n[i] + sum(det)
      pass_n[i, ] <- pass_n[i, ] + colSums(pass)
      sel_n[i, ]  <- sel_n[i, ] + tabulate(winner[det], M)
      noise_s[i]  <- noise_s[i] + sum(noise[cbind(which(det), winner[det])])
    }
    sweep_secs <- sweep_secs + proc.time()[["elapsed"]] - t0
    done <- done + b
  }
  for (i in seq_len(G)) {
    det_rate <- det_n[i] / draws
    p_sel <- sel_n[i, ] / draws
    sel_c <- p_sel / det_rate
    Pg <- fam$Pg; PQg <- fam$PQg; PQ <- fam$PQ
    fx <- list(
      det_rate = det_rate, P1 = pass_n[i, ] / draws, p_sel = p_sel, sel_c = sel_c,
      EnH = n * sum(sel_c * Pg),
      Esens = sum(sel_c * fam$sens_g), Espec = sum(sel_c * fam$spec_g),
      Eppv = sum(sel_c * PQg),
      Enpv = sum(sel_c * vapply(seq_len(M), function(m) {
        pHc <- 1 - Pg[m]; (pHc - (PQ - PQg[m] * Pg[m])) / pHc }, numeric(1))),
      EbetaH = sum(sel_c * fam$beta_g),
      Enaive_bias = noise_s[i] / det_n[i],
      mass_below = sum(sel_c[fam$beta_g < grid$c1[i]]),
      Rdraw = draws)
    results[[i]] <- .fs_oc_result(fx, fam, n, grid$c1[i], grid$c2[i], gate,
                                  pc, draws, seed)
  }
  list(results = results, draw_secs = draw_secs, sweep_secs = sweep_secs)
}
