# =============================================================================
# STAGE M -- UNIT-COST MICROBENCHMARKS
# =============================================================================
# Single-core, deterministic. These confirm on your BLAS / R build what a
# sandbox already suggested; they are NOT the headline (see stage P and E).

.mb_row <- function(id, what, n, cur_ms, new_ms, note = "") {
  data.frame(stringsAsFactors = FALSE,
             id = id, what = what, n = n,
             current_ms = cur_ms, candidate_ms = new_ms,
             speedup = cur_ms / new_ms, note = note)
}

#' Median elapsed time in milliseconds
#'
#' Two things this must get right, both of which a naive implementation gets
#' wrong:
#'
#' 1. `expr` arrives as a promise. `force(expr)` evaluates it ONCE and caches
#'    the value, so `replicate(times, system.time(force(expr)))` times the
#'    first iteration and nothing afterwards -- silently reporting ~0 ms.
#'    Capturing with `substitute()` and re-`eval()`ing in the caller's frame
#'    forces genuine re-evaluation each iteration.
#' 2. `system.time()` elapsed resolution is far too coarse for sub-millisecond
#'    work (glm.fit, coxph.fit). The inner loop is auto-calibrated so each
#'    timed block runs ~20 ms, then divided back out.
.med_ms <- function(expr, times) {
  e  <- substitute(expr)
  pf <- parent.frame()

  one <- system.time(eval(e, pf))[["elapsed"]]
  inner <- if (one < 0.02) {
    min(2000L, max(1L, as.integer(ceiling(0.02 / max(one, 1e-6)))))
  } else 1L

  t <- numeric(times)
  for (i in seq_len(times)) {
    t[i] <- system.time(
      for (j in seq_len(inner)) eval(e, pf)
    )[["elapsed"]] / inner
  }
  stats::median(t) * 1000
}


# M1 -- per-fit cost, coxph() vs coxph.fit(), across subgroup sizes
micro_cox_fit <- function(times = 50L, seed = 11) {
  set.seed(seed)
  out <- list()
  for (N in c(60L, 120L, 250L, 500L)) {
    tt <- rep(0:1, length.out = N)
    Y  <- stats::rexp(N, 0.002); E <- stats::rbinom(N, 1, 0.6)
    d  <- data.frame(Y = Y, Event = E, Treat = tt)
    o  <- order(Y)
    yy <- as.double(Y)[o]; ee <- as.double(E)[o]; xx <- as.double(tt)[o]
    ctl <- survival::coxph.control()

    cur <- .med_ms(suppressWarnings(survival::coxph(
      survival::Surv(Y, Event) ~ Treat, data = d, robust = FALSE,
      model = FALSE, x = FALSE, y = FALSE, init = 0)), times)
    new <- .med_ms(cand_cox_beta(yy, ee, xx, 0, ctl), times)

    out[[length(out) + 1L]] <- .mb_row(
      "M1", "single Cox fit (treatment-only)", N, cur, new,
      "dispatch overhead removed; numerics unchanged")
  }
  do.call(rbind, out)
}


# M2 -- search scorer, with vs without summary.survfit() medians
micro_search_scorer <- function(times = 50L, seed = 12) {
  set.seed(seed)
  out <- list()
  for (N in c(300L, 686L, 1500L)) {
    tt <- rep(0:1, length.out = N)
    Y  <- stats::rexp(N, 0.002); E <- stats::rbinom(N, 1, 0.55)
    idx <- stats::rbinom(N, 1, 0.45)
    ctl <- survival::coxph.control()

    cur <- .med_ms(forestsearch:::fit_cox_for_subgroup(Y, E, tt, idx), times)
    new <- .med_ms(cand_fit_cox_subgroup(Y, E, tt, idx, ctl), times)

    out[[length(out) + 1L]] <- .mb_row(
      "M2", "fit_cox_for_subgroup (per candidate)", N, cur, new,
      "summary.survfit medians replaced by direct $time/$surv read")
  }
  do.call(rbind, out)
}


# M3 -- combination filtering, per-combination vs precomputed
micro_prescreen <- function(times = 10L, seed = 13) {
  set.seed(seed)
  out <- list()
  N <- 686; minp <- 0.05; rmin <- 5
  for (L in c(20L, 40L)) {
    prev <- c(stats::runif(round(L * 0.75), 0.10, 0.60),
              stats::runif(L - round(L * 0.75), 0.005, 0.045))
    zz <- vapply(prev, function(p) stats::rbinom(N, 1, p), numeric(N))
    colnames(zz) <- paste0("z", seq_len(L))
    L_eff <- sum(colMeans(zz) >= minp)

    for (k in c(2L, 3L)) {
      combos <- t(utils::combn(L, k))
      cur <- .med_ms({
        for (i in seq_len(nrow(combos))) {
          x <- zz[, combos[i, ], drop = FALSE]
          if (!forestsearch:::has_positive_variance(x)) next
          if (!forestsearch:::meets_prevalence_threshold(x, minp)) next
          forestsearch:::extract_idx_flagredundancy(x, rmin)
        }
      }, times)
      new <- .med_ms(cand_screen_combos(zz, combos, minp, rmin), times)

      ntot <- function(l) if (k == 2L) l + choose(l, 2) else
        l + choose(l, 2) + choose(l, 3)
      out[[length(out) + 1L]] <- .mb_row(
        sprintf("M3.k%d", k), sprintf("screen %d-factor combinations", k),
        nrow(combos), cur, new,
        sprintf("L=%d -> L_eff=%d; search space %d -> %d (%.0f%% fewer)",
                L, L_eff, ntot(L), ntot(L_eff),
                100 * (1 - ntot(L_eff) / ntot(L))))
    }
  }
  do.call(rbind, out)
}


# M4 -- GLM estimator
micro_glm <- function(times = 50L, seed = 14) {
  set.seed(seed)
  out <- list()
  for (N in c(120L, 300L, 800L)) {
    d <- data.frame(treat = rep(0:1, length.out = N))
    d$Y <- stats::rbinom(N, 1, stats::plogis(-0.4 + 0.6 * d$treat))
    X <- cbind(`(Intercept)` = 1, treat = d$treat)

    cur <- .med_ms({
      f <- stats::glm(Y ~ treat, data = d, family = stats::binomial())
      c(stats::coef(f)[["treat"]], sqrt(diag(stats::vcov(f)))[["treat"]])
    }, times)
    new <- .med_ms(cand_glm_effect(d$Y, X), times)

    out[[length(out) + 1L]] <- .mb_row(
      "M4", "GLM effect + SE (binary/OR)", N, cur, new,
      "glm.fit on prebuilt design; SE from QR")
  }
  do.call(rbind, out)
}


run_micro <- function(cfg) {
  rows <- list(micro_cox_fit(cfg$micro_times),
               micro_search_scorer(cfg$micro_times),
               micro_prescreen(max(5L, cfg$micro_times %/% 5L)))
  if (isTRUE(cfg$run_glm)) rows <- c(rows, list(micro_glm(cfg$micro_times)))
  do.call(rbind, rows)
}
