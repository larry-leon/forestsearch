# =============================================================================
# CANDIDATE IMPLEMENTATIONS
# =============================================================================
# Defined entirely inside the harness. Nothing here is written back to the
# package. Each candidate is paired with the package function it would
# replace, so the gates in 02_gates.R can compare them directly.
#
# SCOPE: unadjusted (treatment-only) models only. When `adjust_covariates` is
# supplied -- especially with `strata()` -- these candidates DO NOT APPLY and
# the package's existing formula path must be retained. The harness does not
# test that path and makes no claim about it.
# =============================================================================


# -----------------------------------------------------------------------------
# C1. Direct Cox fitter  (replaces the coxph() call inside get_split_hr_fast())
# -----------------------------------------------------------------------------
# ~91% of a coxph() call on subgroup-sized data is formula / model.frame /
# model.matrix dispatch, not partial-likelihood iteration. coxph.fit() skips it.
#
# Requirements: `y` ascending; all three vectors double (an integer column
# raises "REAL() can only be applied to a 'numeric', not a 'integer'").

cand_cox_beta <- function(y, e, x, init = 0, ctl = survival::coxph.control()) {
  srv <- structure(cbind(time = y, status = e), class = "Surv", type = "right")
  fit <- tryCatch(
    survival::coxph.fit(
      x = matrix(x, ncol = 1L), y = srv, strata = NULL, offset = NULL,
      init = init, control = ctl, weights = NULL, method = "efron",
      rownames = NULL, resid = FALSE, nocenter = c(0, 1)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  unname(fit$coefficients[1L])
}


# -----------------------------------------------------------------------------
# C2. Consistency split loop
# -----------------------------------------------------------------------------
# Two INDEPENDENT changes, separately switchable so the gates can attribute
# any difference to the right one:
#   engine  = "fast"  -> C1 instead of coxph()
#   futility = TRUE   -> exact early rejection
#
# Futility bound. Pcons = round(mean(flags, na.rm = TRUE), digits) and the
# caller rejects when Pcons < threshold. After i splits with v valid and f
# failures, r = n.splits - i remaining, the largest attainable final mean is
# reached by making every remaining split a valid success:
#     max = (v - f + r) / (v + r)
# Rejection is certain once round(max, digits) < threshold.
#
# Reproducibility. subgroup.consistency() pre-generates one L'Ecuyer stream
# PER CANDIDATE (.make_candidate_rng_seeds), so truncating candidate m cannot
# perturb any other candidate. The sample() call keeps prob = c(0.5, 0.5)
# verbatim: dropping it selects a different algorithm in R and changes draws.

cand_consistency_splits <- function(df.x, N.x, n.splits, hr.consistency,
                                    cox_init = 0, pconsistency.digits = 2,
                                    pconsistency.threshold = NULL,
                                    engine = c("fast", "coxph"),
                                    futility = TRUE) {
  engine <- match.arg(engine)

  o  <- order(df.x$Y)
  yy <- as.double(df.x$Y)[o]
  ee <- as.double(df.x$Event)[o]
  tt <- as.double(df.x$Treat)[o]

  ctl <- survival::coxph.control()
  if (is.null(cox_init) || !is.finite(cox_init)) cox_init <- 0

  fit_half <- if (engine == "fast") {
    function(ii) cand_cox_beta(yy[ii], ee[ii], tt[ii], cox_init, ctl)
  } else {
    function(ii) {
      d <- data.frame(Y = yy[ii], Event = ee[ii], Treat = tt[ii])
      f <- tryCatch(suppressWarnings(survival::coxph(
        survival::Surv(Y, Event) ~ Treat, data = d, robust = FALSE,
        model = FALSE, x = FALSE, y = FALSE, init = cox_init)),
        error = function(e) NULL)
      if (is.null(f)) NA_real_ else unname(f$coefficients[1L])
    }
  }

  flags <- rep(NA_real_, n.splits)
  v <- 0L; f <- 0L

  for (i in seq_len(n.splits)) {
    s1 <- sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
    s1 <- s1[o]
    i1 <- which(s1); i2 <- which(!s1)

    ok <- length(i1) >= 5L && length(i2) >= 5L &&
      sum(ee[i1]) >= 2 && sum(ee[i2]) >= 2

    if (ok) {
      b1 <- fit_half(i1); b2 <- fit_half(i2)
      if (!is.na(b1) && !is.na(b2)) {
        # natural scale, exactly as the package compares; testing
        # b > log(hr.consistency) is equivalent in exact arithmetic but not
        # guaranteed bit-identical at the boundary when hr.consistency != 1
        flags[i] <- as.numeric(exp(b1) > hr.consistency &&
                                 exp(b2) > hr.consistency)
        v <- v + 1L
        if (flags[i] == 0) f <- f + 1L
      }
    }

    if (isTRUE(futility) && !is.null(pconsistency.threshold)) {
      r <- n.splits - i
      if ((v + r) >= 10L) {
        mx <- (v - f + r) / (v + r)
        if (round(mx, pconsistency.digits) < pconsistency.threshold) {
          return(round(mx, pconsistency.digits))   # rejection certain
        }
      }
    }
  }

  if (v < 10) return(NA_real_)
  round(mean(flags, na.rm = TRUE), pconsistency.digits)
}


# -----------------------------------------------------------------------------
# C3. Search scorer without summary.survfit()
# -----------------------------------------------------------------------------
# fit_cox_for_subgroup() spends ~1/3 of its time in
# summary(survfit(...))$table[, "median"]. Under the default
# m1.threshold = Inf those medians never enter selection
# (subgroup_consistency_main.R:514) -- they are display columns only.
# Reading survfit$time / $surv directly reproduces them.
#
# NOTE: summary.coxph builds its CI with qnorm(0.975) = 1.959964, not 1.96.
# The candidate uses qnorm to stay comparable to the package.

cand_fit_cox_subgroup <- function(yy, dd, tt, id.x,
                                  ctl = survival::coxph.control(),
                                  z = stats::qnorm(0.975)) {
  keep <- which(id.x == 1)
  y <- as.double(yy[keep]); e <- as.double(dd[keep]); x <- as.double(tt[keep])
  o <- order(y); y <- y[o]; e <- e[o]; x <- x[o]

  srv <- structure(cbind(time = y, status = e), class = "Surv", type = "right")
  f <- tryCatch(survival::coxph.fit(
    x = matrix(x, ncol = 1L), y = srv, strata = NULL, offset = NULL,
    init = 0, control = ctl, weights = NULL, method = "efron",
    rownames = NULL, resid = FALSE, nocenter = c(0, 1)),
    error = function(e) NULL)
  if (is.null(f)) return(NULL)

  b  <- f$coefficients[1L]
  se <- sqrt(f$var[1, 1])

  meds <- tryCatch({
    sf  <- survival::survfit(srv ~ x)
    st  <- sf$strata
    idx <- rep(seq_along(st), st)
    vapply(seq_along(st), function(g) {
      tg <- sf$time[idx == g]; sg <- sf$surv[idx == g]
      w  <- which(sg <= 0.5)
      if (!length(w)) NA_real_ else tg[w[1L]]
    }, numeric(1))
  }, error = function(e) NULL)
  if (is.null(meds)) return(NULL)

  list(hr = unname(exp(b)), lower = unname(exp(b - z * se)),
       upper = unname(exp(b + z * se)), med0 = meds[1], med1 = meds[2])
}


# -----------------------------------------------------------------------------
# C4. Vectorised combination pre-screen
# -----------------------------------------------------------------------------
# Both reformulations are EXACT, not approximations:
#   * prevalence is per-column -- meets_prevalence_threshold() demands
#     all(colMeans(x) >= minp), so a column under minp fails in EVERY
#     combination containing it and can be dropped before combinations are
#     generated, shrinking L (and the count as O(L^maxk));
#   * has_positive_variance() computes t(x) %*% x per combination; every such
#     sub-block is a slice of crossprod(zz), one BLAS call.

cand_prescreen <- function(zz, minp) {
  list(keep_col = colMeans(zz) >= minp, CP = crossprod(zz))
}

cand_screen_combos <- function(zz, combos, minp, rmin, pre = NULL) {
  if (is.null(pre)) pre <- cand_prescreen(zz, minp)
  keep <- pre$keep_col; CP <- pre$CP
  out <- logical(nrow(combos))
  for (i in seq_len(nrow(combos))) {
    cid <- combos[i, ]
    if (!all(keep[cid])) next
    if (!all(CP[cid, cid, drop = FALSE] > 0)) next
    if (forestsearch:::extract_idx_flagredundancy(
      zz[, cid, drop = FALSE], rmin)$flag.redundant) next
    out[i] <- TRUE
  }
  out
}


# -----------------------------------------------------------------------------
# C5. GLM estimator via glm.fit()
# -----------------------------------------------------------------------------
# The package path is stats::glm(formula, data) followed by stats::vcov(fit),
# which routes through summary.glm. glm.fit() on a prebuilt design matrix with
# the SE read from the QR reproduces both to machine precision.

cand_glm_effect <- function(y, X, family = stats::binomial(link = "logit"),
                            col = 2L) {
  fit <- tryCatch(stats::glm.fit(x = X, y = y, family = family),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) {
    return(list(estimate = NA_real_, se = NA_real_, converged = FALSE))
  }
  V <- chol2inv(qr.R(fit$qr))
  list(estimate = unname(fit$coefficients[col]),
       se = sqrt(V[col, col]), converged = TRUE)
}
