# guohe_adaptive_r.R
#
# Cross-validated choice of the shrinkage tuning parameter r -- Algorithm 2 of
# Guo & He (2021), Section 2.5 -- for use with `guohe_algorithm3()`.
#
# WHY THIS MATTERS FOR THE COMPARISON. The offsets leave a residual gap of
# n^(r - 1/2) times the observed gap between a candidate and the winner. With a
# very small r that residual gap can fall BELOW bootstrap noise, so genuinely
# inferior candidates re-enter the resampled competition and the correction
# over-shoots. Guo & He note in Section 2.5 that a small r preserves coverage in
# finite samples at the cost of possibly conservative bounds. A head-to-head
# against a TUNING-FREE method (MR) with r pinned at one arbitrary value would
# therefore be a weak comparison; using the authors' own adaptive rule removes
# that objection.
#
# ALGORITHM 2, as published:
#   1. Partition the data into v approximately equal subsamples.
#   2. For each candidate r_l and each fold j:
#        - training data = all folds except j; reference data = fold j;
#        - compute the bias-reduced estimator beta_hat_max,reduced,j(r_l) on the
#          TRAINING data;
#        - on the REFERENCE data estimate each subgroup effect beta_hat_{i,j}
#          and its standard error sigma_hat_{i,j};
#        - h_{i,j}(r_l) = (beta_hat_max,reduced,j(r_l) - beta_hat_{i,j})^2
#                         - sigma_hat_{i,j}^2 .
#   3. Choose r = argmin_l { min_i [ sum_j h_{i,j}(r_l) / v ] }.
#
# The sigma_hat^2 subtraction removes the sampling variance of the reference
# estimate, so h targets the squared distance to the true beta_i rather than to a
# noisy estimate of it. Theorem 3 justifies the inner min over i.
#
# TWO EFFICIENCIES over a naive transcription:
#   (a) The reference-fold quantities beta_hat_{i,j} and sigma_hat_{i,j} do NOT
#       depend on r, so they are computed ONCE per fold and reused across the
#       whole r grid -- a length(r_grid)-fold saving.
#   (b) For UNADJUSTED binary/continuous outcomes both the coefficient and its
#       standard error have closed forms, so all candidates in a fold are scored
#       with a single matrix product instead of k model fits. These are EXACT:
#         binary     : log-OR = log(ad/bc), SE = sqrt(1/a + 1/b + 1/c + 1/d)
#                      (Woolf; the saturated-model vcov that glm returns)
#         continuous : MD = m1 - m0, SE = sigma_hat * sqrt(1/n1 + 1/n0) with
#                      sigma_hat^2 the pooled residual variance (the lm SE)
#       The gate matches `guohe_algorithm3()`: adjustment or a survival outcome
#       disables it. This duplicates a little of that file's fast machinery
#       deliberately, to avoid refactoring an already-validated engine.
#
# All arithmetic is on the ORIENTED SCORE scale (orient * coefficient).

suppressWarnings(suppressMessages(library(survival)))

# ---- internal: coefficient AND standard error for one subgroup (slow path) --
.ar_coef_se <- function(sub, outcome, treatment, time, event, y, min_events,
                        adjust_covariates = NULL) {
  na2 <- c(NA_real_, NA_real_)
  tr <- sub[[treatment]]
  if (length(unique(tr)) < 2L) return(na2)
  if (min(table(tr)) < min_events) return(na2)
  tryCatch(
    suppressWarnings({
      if (outcome == "survival") {
        ev <- sub[[event]]
        if (sum(ev == 1) < min_events) return(na2)
        if (is.null(adjust_covariates)) {
          fit <- survival::coxph(survival::Surv(sub[[time]], ev) ~ tr)
          c(unname(fit$coefficients[[1L]]),
            unname(sqrt(diag(stats::vcov(fit)))[1L]))
        } else {
          md <- data.frame(.tr = tr, .time = sub[[time]], .ev = ev)
          md[adjust_covariates] <- sub[adjust_covariates]
          fml <- stats::as.formula(paste("survival::Surv(.time, .ev) ~ .tr +",
                                         paste(adjust_covariates, collapse = " + ")))
          fit <- survival::coxph(fml, data = md)
          c(unname(fit$coefficients[[".tr"]]),
            unname(sqrt(diag(stats::vcov(fit)))[".tr"]))
        }
      } else {
        yy <- sub[[y]]
        if (outcome == "binary" && length(unique(yy)) < 2L) return(na2)
        fam <- if (outcome == "binary") stats::binomial() else stats::gaussian()
        if (is.null(adjust_covariates)) {
          sm <- summary(stats::glm(yy ~ tr, family = fam))$coefficients
          if (!("tr" %in% rownames(sm))) return(na2)
          c(unname(sm["tr", 1L]), unname(sm["tr", 2L]))
        } else {
          md <- data.frame(.y = yy, .tr = tr)
          md[adjust_covariates] <- sub[adjust_covariates]
          fml <- stats::as.formula(paste(".y ~ .tr +",
                                         paste(adjust_covariates, collapse = " + ")))
          sm <- summary(stats::glm(fml, family = fam, data = md))$coefficients
          if (!(".tr" %in% rownames(sm))) return(na2)
          c(unname(sm[".tr", 1L]), unname(sm[".tr", 2L]))
        }
      }
    }),
    error = function(e) na2
  )
}

# ---- internal: closed-form coefficient AND SE for every candidate at once ----
# `w` is a 0/1 weight vector selecting the reference fold.
.ar_ref_fast <- function(Mmat, cellmat, w, outcome, min_events, orient) {
  cnt <- crossprod(Mmat, w * cellmat)
  if (outcome == "binary") {
    a <- cnt[, 1]; b <- cnt[, 2]; cc <- cnt[, 3]; dd <- cnt[, 4]
    est <- log(a) + log(dd) - log(b) - log(cc)
    se  <- sqrt(1 / a + 1 / b + 1 / cc + 1 / dd)
    bad <- (a + b) < min_events | (cc + dd) < min_events |
      a <= 0 | b <= 0 | cc <= 0 | dd <= 0
  } else {
    n1 <- cnt[, 1]; n0 <- cnt[, 2]; s1 <- cnt[, 3]
    s0 <- cnt[, 4]; q1 <- cnt[, 5]; q0 <- cnt[, 6]
    est <- s1 / n1 - s0 / n0
    rss <- pmax((q1 - s1^2 / n1) + (q0 - s0^2 / n0), 0)
    dfr <- n1 + n0 - 2
    se  <- sqrt((rss / dfr) * (1 / n1 + 1 / n0))
    bad <- n1 < min_events | n0 < min_events | dfr <= 0
  }
  est[bad] <- NA_real_; se[bad] <- NA_real_
  list(est = unname(orient * est), se = unname(se))   # SE invariant to sign flip
}

# report scale -> oriented score scale (exact inverse of guohe_algorithm3's map)
.ar_to_score <- function(report, orient, outcome) {
  if (outcome == "continuous") orient * report else orient * log(report)
}

#' Cross-validated choice of the shrinkage tuning parameter r
#'
#' Implements Algorithm 2 of Guo and He (2021). Requires `guohe_algorithm3()`
#' to be available (source `guohe_algorithm3.R` first).
#'
#' @inheritParams guohe_algorithm3
#' @param r_grid Numeric vector of candidate values in (0, 0.5).
#' @param v Number of cross-validation folds.
#' @param refit Logical; if `TRUE`, also return a full-data fit at the selected
#'   `r`.
#'
#' @return An object of class `"guohe_ar"`: the selected `r`, the CV objective
#'   over the grid, the per-candidate objective at the selected `r`, and
#'   optionally the refitted full-data result.
#' @export
guohe_adaptive_r <- function(data,
                             outcome = c("survival", "binary", "continuous"),
                             treatment,
                             candidates,
                             time = NULL,
                             event = NULL,
                             y = NULL,
                             orient = -1,
                             r_grid = c(0.03, 0.10, 0.20, 0.30, 0.40, 0.45),
                             v = 5L,
                             B = 200L,
                             level = 0.05,
                             seed = NULL,
                             min_events = 5L,
                             refit = TRUE,
                             adjust_covariates = NULL,
                             fast = NULL) {
  outcome <- match.arg(outcome)
  if (!exists("guohe_algorithm3", mode = "function")) {
    stop("guohe_algorithm3() not found; source 'guohe_algorithm3.R' first.")
  }
  stopifnot(
    is.data.frame(data),
    all(r_grid > 0), all(r_grid < 0.5),
    v >= 2L, B >= 2L,
    orient %in% c(-1, 1)
  )
  n <- nrow(data)
  k <- length(candidates)
  if (!is.null(seed)) set.seed(seed)

  fast_ok <- outcome %in% c("binary", "continuous") && is.null(adjust_covariates)
  if (is.null(fast)) {
    fast <- fast_ok
  } else if (isTRUE(fast) && !fast_ok) {
    stop("fast = TRUE is invalid here: closed forms exist only for unadjusted ",
         "binary or continuous outcomes.")
  }

  folds <- sample(rep_len(seq_len(v), n))
  memb <- as.data.frame(lapply(candidates, function(cn) as.integer(data[[cn]] == 1)))
  names(memb) <- candidates

  if (isTRUE(fast)) {
    Mmat <- as.matrix(memb) * 1.0
    tr_v <- as.numeric(data[[treatment]] == 1)
    y_v  <- as.numeric(data[[y]])
    cellmat <- if (outcome == "binary") {
      cbind(tr_v * y_v, tr_v * (1 - y_v), (1 - tr_v) * y_v, (1 - tr_v) * (1 - y_v))
    } else {
      cbind(tr_v, 1 - tr_v, tr_v * y_v, (1 - tr_v) * y_v,
            tr_v * y_v^2, (1 - tr_v) * y_v^2)
    }
  }

  # ---- (a) reference-fold statistics: r-INDEPENDENT, so computed once -------
  ref <- lapply(seq_len(v), function(j) {
    if (isTRUE(fast)) {
      .ar_ref_fast(Mmat, cellmat, as.numeric(folds == j), outcome, min_events,
                   orient)
    } else {
      rf <- data[folds == j, , drop = FALSE]
      mr <- memb[folds == j, , drop = FALSE]
      cs <- vapply(seq_len(k), function(i) {
        .ar_coef_se(rf[mr[[i]] == 1L, , drop = FALSE], outcome, treatment,
                    time, event, y, min_events, adjust_covariates)
      }, numeric(2))
      list(est = orient * cs[1, ], se = cs[2, ])
    }
  })

  # ---- (b) training-fold bias-reduced estimator, per r ----------------------
  h <- array(NA_real_, dim = c(k, v, length(r_grid)))
  for (l in seq_along(r_grid)) {
    for (j in seq_len(v)) {
      fit_tr <- try(suppressWarnings(
        guohe_algorithm3(data[folds != j, , drop = FALSE], outcome = outcome,
                         treatment = treatment, candidates = candidates,
                         time = time, event = event, y = y, orient = orient,
                         B = B, r = r_grid[l], level = level,
                         min_events = min_events, diagnostics = FALSE,
                         adjust_covariates = adjust_covariates, fast = fast)
      ), silent = TRUE)
      if (inherits(fit_tr, "try-error")) next
      d_score <- .ar_to_score(fit_tr$debiased, orient, outcome)
      h[, j, l] <- (d_score - ref[[j]]$est)^2 - ref[[j]]$se^2
    }
  }

  # ---- (c) average over folds, minimise over candidates, then over r --------
  H <- apply(h, c(1, 3),
             function(z) if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE))
  obj <- apply(H, 2, function(z) if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE))
  if (all(is.na(obj))) stop("Cross-validation objective undefined for every r.")
  r_hat <- r_grid[which.min(obj)]

  structure(list(
    r_hat         = r_hat,
    r_grid        = r_grid,
    objective     = obj,
    per_candidate = stats::setNames(H[, which.min(obj)], candidates),
    v = v, B = B, n = n, n_candidates = k, fast = isTRUE(fast),
    outcome = outcome, orient = orient, level = level,
    adjust_covariates = adjust_covariates,
    fit = if (refit) {
      suppressWarnings(
        guohe_algorithm3(data, outcome = outcome, treatment = treatment,
                         candidates = candidates, time = time, event = event,
                         y = y, orient = orient, B = B, r = r_hat,
                         level = level, seed = seed, min_events = min_events,
                         adjust_covariates = adjust_covariates, fast = fast)
      )
    } else NULL
  ), class = "guohe_ar")
}

#' @export
print.guohe_ar <- function(x, ...) {
  cat("Guo-He Algorithm 2: cross-validated tuning parameter\n")
  cat(sprintf("  n = %d, candidates = %d, folds = %d, B = %d%s\n",
              x$n, x$n_candidates, x$v, x$B,
              if (isTRUE(x$fast)) "  [closed-form]" else ""))
  cat("  CV objective over the grid (smaller is better):\n")
  for (i in seq_along(x$r_grid)) {
    cat(sprintf("    r = %.2f  ->  %s%s\n", x$r_grid[i],
                if (is.na(x$objective[i])) "NA" else sprintf("%.5f", x$objective[i]),
                if (isTRUE(x$r_grid[i] == x$r_hat)) "   <-- selected" else ""))
  }
  if (!is.null(x$fit)) {
    cat("\n  Full-data fit at the selected r:\n")
    print(x$fit)
  }
  invisible(x)
}
