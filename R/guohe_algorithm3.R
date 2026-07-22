# guohe_algorithm3.R
#
# Guo & He (2021) de-biasing for a LARGE enumerated candidate family: the finite
# realization of their Section 3 ("post-hoc identified subgroups") Algorithm 3,
# where the family {S(c) : c in D} is supplied as materialized membership
# columns. This is the FS-side comparator for `fs-glms-interpretable`, built to
# be compared against MR and FB on the same selected subgroup.
#
# Relationship to `guohe_comparator.R`: that file implements the Section 2
# (small predefined family) engine and is validated to machine precision against
# the authors' `synthetic code.R` and to ~1e-3 against the published Table 2.
# The arithmetic here is identical -- sup over the family, shrinkage offsets
# against the GLOBAL maximum, pair bootstrap -- so the one-sided bound returned
# by this file reproduces that file's `bound` exactly on the same inputs. What
# Section 3 adds is (i) theoretical licence for a large/continuum family under
# Donsker-type conditions, and (ii) the interval machinery below.
#
# INTERVALS. Three objects are returned, deliberately distinguished:
#
#   * `bound_one_sided` -- Guo & He's CERTIFIED object. Inverts the bootstrap law
#     of T* = sqrt(n) * (sup_c(beta*(c) + d(c)) - gamma_max_hat) at q_{1-level}.
#     This is what their Theorems 1/4/5 and their simulations cover.
#
#   * `ci_two_sided` -- quantile inversion of the SAME T* at (level/2, 1-level/2).
#     A defensible extension: Theorems 1/4/5 are Kolmogorov-distance statements
#     (uniform in x), and the limiting sup-of-Gaussian cdf is atomless by
#     anti-concentration, so every quantile is consistently estimated. NOTE the
#     conservative endpoint uses q_{1-level/2}, so it is STRICTLY WIDER than
#     `bound_one_sided`; it equals the one-sided bound at level/2. Use this for
#     the MR head-to-head.
#
#   * `ci_wald` -- Wald interval (debiased point estimate +/- z * bootstrap SD),
#     the Zhao, Ivanova & Fine (2023) style construction. HEURISTIC: it is not an
#     inversion of T*, it is not covered by Guo & He's theory, and because T* is
#     skewed (it is a maximum) it need not agree with `ci_two_sided`. Provided
#     only for comparability with that literature; do not present as certified.
#
# Depends on {survival}, {stats}, and (optionally) {future.apply}.
#
# References:
#   Guo, X. and He, X. (2021). Inference on selected subgroups in clinical
#     trials. JASA 116(535):1498-1506.  [Algorithm 3, Theorems 4-5]
#   Zhao, B., Ivanova, A. and Fine, J. (2023). Inference on subgroups identified
#     based on a heterogeneous treatment effect in a post hoc analysis of a
#     clinical trial. Clinical Trials 20(4):370-379.

# ---- internal: one within-subgroup treatment coefficient ------------------
# Returns the treatment coefficient on the model's natural scale
# (log-HR / log-OR / mean difference), or NA_real_ if not estimable.
#' @noRd
.g3_coef <- function(sub, outcome, treatment, time, event, y, min_events,
                     adjust_covariates = NULL) {
  tr <- sub[[treatment]]
  if (length(unique(tr)) < 2L) {
    return(NA_real_)
  }
  if (min(table(tr)) < min_events) {
    return(NA_real_)
  }
  tryCatch(
    suppressWarnings({
      if (outcome == "survival") {
        ev <- sub[[event]]
        if (sum(ev == 1) < min_events) {
          return(NA_real_)
        }
        if (is.null(adjust_covariates)) {
          fit <- survival::coxph(survival::Surv(sub[[time]], ev) ~ tr)
          unname(fit$coefficients[[1L]])
        } else {
          md <- data.frame(.tr = tr, .time = sub[[time]], .ev = ev)
          md[adjust_covariates] <- sub[adjust_covariates]
          fml <- stats::as.formula(paste("survival::Surv(.time, .ev) ~ .tr +",
                                         paste(adjust_covariates, collapse = " + ")))
          unname(survival::coxph(fml, data = md)$coefficients[[".tr"]])
        }
      } else {
        yy <- sub[[y]]
        if (outcome == "binary" && length(unique(yy)) < 2L) {
          return(NA_real_)
        }
        fam <- if (outcome == "binary") stats::binomial() else stats::gaussian()
        if (is.null(adjust_covariates)) {
          unname(stats::coef(stats::glm(yy ~ tr, family = fam))[["tr"]])
        } else {
          md <- data.frame(.y = yy, .tr = tr)
          md[adjust_covariates] <- sub[adjust_covariates]
          fml <- stats::as.formula(paste(".y ~ .tr +", paste(adjust_covariates, collapse = " + ")))
          unname(stats::coef(stats::glm(fml, family = fam, data = md))[[".tr"]])
        }
      }
    }),
    error = function(e) NA_real_
  )
}

# ---- internal: lean unadjusted two-group Cox ------------------------------
# `.g3_coef()` receives a data-frame subset of the FULL frame -- which at
# forest-search scale carries ~1700 membership columns -- and then uses exactly
# three of them. This variant takes plain vectors and calls
# `survival::coxph.fit()` directly, skipping formula parsing, model.frame
# construction, and assembly of the large coxph return object.
#
# It is the SAME estimator, not an approximation: `method = "efron"` matches
# coxph()'s default tie handling, and the two agree to ~1e-15 (asserted in the
# smoke test). Guards are replicated from `.g3_coef()` exactly so the two paths
# declare the same candidates non-estimable.
#' @noRd
.g3_cox_lean <- function(time_v, ev_v, tr_v, min_events, ctl) {
  if (length(unique(tr_v)) < 2L) return(NA_real_)
  if (min(table(tr_v)) < min_events) return(NA_real_)
  if (sum(ev_v == 1) < min_events) return(NA_real_)
  tryCatch(
    suppressWarnings(unname(
      survival::coxph.fit(
        x = matrix(as.double(tr_v), ncol = 1L),
        y = survival::Surv(time_v, ev_v),
        strata = NULL, offset = NULL, init = 0, control = ctl,
        weights = NULL, method = "efron", rownames = NULL
      )$coefficients[[1L]]
    )),
    error = function(e) NA_real_
  )
}

#' De-biased inference for the best subgroup over a large enumerated family
#'
#' Finite realization of Algorithm 3 of Guo and He (2021): the candidate family
#' is supplied as materialized 0/1 membership columns, the most extreme oriented
#' effect is selected, and its selection optimism is removed by a pair (case)
#' bootstrap with deterministic shrinkage offsets computed against the global
#' maximum. Survival, binary, and continuous outcomes are supported.
#'
#' Membership is computed once on the full data and carried through the
#' bootstrap; cut points are NOT re-derived per resample. This is required by
#' the method (the index set is held fixed) and is what makes the comparison
#' against multiplier resampling like-for-like; regenerating the family per
#' resample is the full-bootstrap object instead.
#'
#' The oriented score is `orient * b`, with `b` the within-subgroup treatment
#' coefficient on the natural scale. `orient = -1` treats a protective ratio
#' (HR or OR below 1) as the effect of interest, matching Guo and He's MONET1
#' convention; `orient = +1` treats a harmful ratio as the effect of interest,
#' which is the setting for a forest-search harm subgroup.
#'
#' This is an independent implementation of the published method, validated by
#' reproducing the simulation study of Section 5 of Guo and He (2021).
#'
#' @param data A data frame.
#' @param outcome One of `"survival"`, `"binary"`, `"continuous"`.
#' @param treatment Name of the 0/1 treatment column.
#' @param candidates Character vector of names of 0/1 membership columns forming
#'   the enumerated family.
#' @param time,event Survival only: names of the time and event columns.
#' @param y Binary or continuous only: name of the outcome column.
#' @param orient Either `-1` or `+1`; sign orienting the effect of interest.
#' @param B Number of bootstrap resamples.
#' @param r Shrinkage tuning parameter, strictly between 0 and 0.5. Guo and He's
#'   Algorithm 2 gives a cross-validated choice; a fixed small value is used here.
#' @param level One-sided level; the two-sided interval has coverage `1 - level`.
#' @param seed Optional integer for reproducible resampling.
#' @param min_events Minimum events (survival) or minimum per-arm count for a
#'   candidate to be estimable.
#' @param diagnostics Logical; if `TRUE`, record per-resample estimability.
#' @param adjust_covariates Optional character vector of covariate names for the
#'   within-subgroup model, matching the argument of the same name in
#'   `forestsearch()`'s `.estimate_or()` / `.make_lm_estimator()`. Default `NULL`
#'   is the unadjusted model. Supplying it DISABLES the closed-form fast path,
#'   because an adjusted coefficient is a partial regression coefficient, not a
#'   2x2 log-odds-ratio or a difference of arm means.
#'
#'   NOTE: `forestsearch()` has a SECOND adjustment axis, `ps_adjust_method`
#'   (`"iptw"` weighted fitting, `"dr_gcomp"` Bang-Robins G-computation). Those
#'   modes are not implemented in this comparator, and the closed forms would be
#'   invalid under them too -- IPTW fits a weighted logistic, and DR reports a
#'   marginal log-OR from G-computed probabilities; neither equals the 2x2
#'   statistic. If the comparator is extended to those modes, the `fast_ok` gate
#'   must be extended with them.
#' @param fast Logical or `NULL`. Closed-form evaluation of the within-subgroup
#'   effect, available only for UNADJUSTED binary and continuous outcomes. These
#'   are EXACT, not approximations: with a single binary covariate the logistic
#'   model is saturated, so its MLE is the 2x2 log-odds-ratio, and the Gaussian
#'   coefficient is the difference of arm means. `NULL` (default) auto-enables
#'   it whenever it is valid. `TRUE` errors if the configuration makes it
#'   invalid (survival, or any `adjust_covariates`), rather than silently falling back.
#' @param parallel Logical; if `TRUE`, distribute the bootstrap model fitting
#'   with `future.apply`. Requires a `future::plan()` to have been set by the
#'   caller. Results are identical to `parallel = FALSE` for the same `seed`,
#'   because resample indices are drawn serially before any distribution.
#'   Do NOT enable this inside an already-parallel outer loop.
#'
#' @return An object of class `"guohe_a3"`: a list with the selected candidate
#'   (`selected`), the naive and de-biased effects on the reporting scale
#'   (`naive`, `debiased`), the certified one-sided bound (`bound_one_sided`,
#'   with its conservative side in `bound_side`), the two-sided
#'   quantile-inversion interval (`ci_two_sided`), the heuristic Wald interval
#'   (`ci_wald`), and bookkeeping fields (`n`, `B`, `B_requested`, `r`,
#'   `level`, `n_candidates`, `n_estimable`, `fast`, `lean_cox`,
#'   `adjust_covariates`, and, when `diagnostics = TRUE`, the per-candidate
#'   resample `drop_rate`).
#'
#' @references
#' Guo, X. and He, X. (2021). Inference on Selected Subgroups in Clinical
#' Trials. \emph{Journal of the American Statistical Association},
#' \strong{116}(535), 1498--1506. \doi{10.1080/01621459.2020.1740096}
#'
#' @examples
#' set.seed(1)
#' n <- 160
#' df <- data.frame(
#'   treat = rbinom(n, 1, 0.5),
#'   y     = rbinom(n, 1, 0.3),
#'   g1    = rbinom(n, 1, 0.5),
#'   g2    = rbinom(n, 1, 0.5),
#'   g3    = rbinom(n, 1, 0.5)
#' )
#' fit <- guohe_algorithm3(
#'   df, outcome = "binary", treatment = "treat",
#'   candidates = c("g1", "g2", "g3"), y = "y",
#'   orient = +1, B = 50, seed = 2
#' )
#' fit
#'
#' @importFrom survival coxph Surv coxph.fit coxph.control
#' @importFrom stats as.formula coef glm binomial gaussian quantile qnorm sd setNames
#' @importFrom future.apply future_lapply
#' @export
guohe_algorithm3 <- function(data,
                             outcome = c("survival", "binary", "continuous"),
                             treatment,
                             candidates,
                             time = NULL,
                             event = NULL,
                             y = NULL,
                             orient = -1,
                             B = 2000L,
                             r = 0.03,
                             level = 0.05,
                             seed = NULL,
                             min_events = 5L,
                             diagnostics = TRUE,
                             parallel = FALSE,
                             adjust_covariates = NULL,
                             fast = NULL) {
  outcome <- match.arg(outcome)
  stopifnot(
    is.data.frame(data),
    orient %in% c(-1, 1),
    r > 0, r < 0.5,
    level > 0, level < 1,
    B >= 2L,
    length(candidates) >= 1L
  )
  need <- c(treatment, candidates)
  need <- if (outcome == "survival") c(need, time, event) else c(need, y)
  miss <- setdiff(need, names(data))
  if (length(miss)) {
    stop("Columns not found in `data`: ", paste(miss, collapse = ", "))
  }
  n <- nrow(data)

  membership <- as.data.frame(
    lapply(candidates, function(cn) as.integer(data[[cn]] == 1))
  )
  names(membership) <- candidates
  n_cand <- ncol(membership)

  # ---- closed-form eligibility --------------------------------------------
  # EXACT, not approximate: with a single binary covariate the logistic model is
  # saturated (its MLE is the 2x2 log-odds-ratio) and the Gaussian coefficient
  # is the difference of arm means. Neither identity survives adjustment, and
  # Cox has no closed form -- hence the gate.
  fast_ok <- outcome %in% c("binary", "continuous") && is.null(adjust_covariates)
  if (is.null(fast)) {
    fast <- fast_ok
  } else if (isTRUE(fast) && !fast_ok) {
    stop("fast = TRUE is invalid here: closed forms exist only for unadjusted ",
         "binary or continuous outcomes (got outcome = '", outcome, "'",
         if (!is.null(adjust_covariates)) ", with adjustment covariates" else "", ").")
  }
  if (!is.null(adjust_covariates)) {
    miss_adj <- setdiff(adjust_covariates, names(data))
    if (length(miss_adj)) {
      stop("Adjustment covariates not found: ", paste(miss_adj, collapse = ", "))
    }
  }

  # Closed-form scorer. `w` is the per-subject resample multiplicity, so ONE
  # crossprod scores the whole family: t(M) %*% (w * cells) returns, for every
  # candidate, the sufficient statistics of its resampled subgroup. This
  # replaces n_cand model fits per draw with a single matrix product.
  if (isTRUE(fast)) {
    Mmat <- as.matrix(membership) * 1.0
    tr_v <- as.numeric(data[[treatment]] == 1)
    y_v  <- as.numeric(data[[y]])
    cellmat <- if (outcome == "binary") {
      cbind(tr_v * y_v, tr_v * (1 - y_v), (1 - tr_v) * y_v, (1 - tr_v) * (1 - y_v))
    } else {
      cbind(tr_v, 1 - tr_v, tr_v * y_v, (1 - tr_v) * y_v)
    }
    score_fast <- function(w) {
      cnt <- crossprod(Mmat, w * cellmat)
      if (outcome == "binary") {
        a <- cnt[, 1]; b <- cnt[, 2]; cc <- cnt[, 3]; dd <- cnt[, 4]
        val <- log(a) + log(dd) - log(b) - log(cc)
        bad <- (a + b) < min_events | (cc + dd) < min_events |
          a <= 0 | b <= 0 | cc <= 0 | dd <= 0
      } else {
        n1 <- cnt[, 1]; n0 <- cnt[, 2]; s1 <- cnt[, 3]; s0 <- cnt[, 4]
        val <- s1 / n1 - s0 / n0
        bad <- n1 < min_events | n0 < min_events
      }
      val[bad] <- NA_real_
      unname(orient * val)
    }
  }

  # ---- lean refit path ----------------------------------------------------
  # Engages exactly where the closed forms do NOT: unadjusted survival, which is
  # where the remaining cost lives. Pre-extracting the three columns coxph needs
  # avoids materialising a wide data-frame subset once per candidate per draw.
  lean_cox <- identical(outcome, "survival") && is.null(adjust_covariates) &&
    !isTRUE(fast)
  ctl_lean <- NULL
  if (isTRUE(lean_cox)) {
    ctl_lean <- survival::coxph.control()
    # coxph.fit() is exported, but probe once anyway: a signature change across
    # survival versions should degrade to the data-frame path, not error out.
    lean_cox <- isTRUE(tryCatch(
      is.finite(.g3_cox_lean(seq_len(10), rep(c(1L, 0L), 5),
                             rep(c(0, 1), each = 5), 1L, ctl_lean)),
      error = function(e) FALSE))
  }
  if (isTRUE(lean_cox)) {
    Mslow   <- as.matrix(membership)
    tr_slow <- data[[treatment]]
    tm_slow <- data[[time]]
    ev_slow <- data[[event]]
  }

  # oriented score for every candidate on a given row-index set
  score_vec <- function(idx) {
    if (isTRUE(fast)) return(score_fast(tabulate(idx, nbins = n)))
    if (isTRUE(lean_cox)) {
      tb <- tr_slow[idx]; yb <- tm_slow[idx]; eb <- ev_slow[idx]
      return(vapply(seq_len(n_cand), function(k) {
        s <- Mslow[idx, k] == 1L
        orient * .g3_cox_lean(yb[s], eb[s], tb[s], min_events, ctl_lean)
      }, numeric(1)))
    }
    d <- data[idx, , drop = FALSE]
    m <- membership[idx, , drop = FALSE]
    vapply(
      seq_len(n_cand),
      function(k) {
        sub <- d[m[[k]] == 1L, , drop = FALSE]
        orient * .g3_coef(sub, outcome, treatment, time, event, y, min_events,
                          adjust_covariates)
      },
      numeric(1)
    )
  }

  # ---- observed competition: gamma_max_hat = sup_c beta_hat(c) -------------
  obs <- score_vec(seq_len(n))
  ok <- is.finite(obs)
  if (!any(ok)) {
    stop("No estimable candidate on the full sample.")
  }
  if (!all(ok)) {
    warning(sum(!ok), " of ", n_cand,
            " candidates not estimable on the full sample; dropped.")
  }
  gamma_max <- max(obs[ok])
  sel <- candidates[which.max(replace(obs, !ok, -Inf))]

  # d(c) = (1 - n^(r - 0.5)) * (gamma_max_hat - beta_hat(c))
  offset <- (1 - n^(r - 0.5)) * (gamma_max - obs)

  # ---- pair bootstrap ------------------------------------------------------
  # U_b = sup_c (beta*_b(c) + d(c)) - gamma_max_hat   (T*_b = sqrt(n) * U_b;
  # the sqrt(n) cancels on inversion, so everything below is on the U scale)
  #
  # PARALLELISM. All B resample index vectors are drawn SERIALLY under `seed`
  # first, and only the (expensive) model fitting is distributed. Because the
  # RNG draws therefore occur in the same order regardless of backend, the
  # parallel path reproduces the serial path EXACTLY -- the property that makes
  # the speed-up verifiable rather than merely plausible. (Index pre-generation
  # costs O(B * n) integers, negligible at trial sizes; for very large n switch
  # to pre-generating per-draw seeds instead.)
  if (!is.null(seed)) set.seed(seed)
  idx_list <- lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))

  one_draw <- function(idx) {
    val <- score_vec(idx) + offset
    fin <- is.finite(val)
    list(u = if (any(fin)) max(val[fin]) - gamma_max else NA_real_,
         bad = if (diagnostics) !fin else NULL)
  }

  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel = TRUE requires the 'future.apply' package.")
    }
    # future.seed = FALSE is correct here: no RNG is used inside one_draw --
    # the indices were already drawn above. `.g3_coef` qualifies survival and
    # stats with `::`, so workers need no attached packages.
    res <- future.apply::future_lapply(
      idx_list, one_draw,
      future.seed = FALSE,
      future.globals = c("score_vec", "offset", "gamma_max", "diagnostics",
                         ".g3_coef", ".g3_cox_lean", "fast", "n")
    )
  } else {
    res <- lapply(idx_list, one_draw)
  }

  U <- vapply(res, function(z) z$u, numeric(1))
  drop_ct <- if (diagnostics) {
    Reduce(`+`, lapply(res, function(z) as.integer(z$bad)))
  } else {
    integer(n_cand)
  }
  n_bad <- sum(!is.finite(U))
  if (n_bad > 0) {
    warning(n_bad, " of ", B, " resamples had no estimable candidate; dropped.")
    U <- U[is.finite(U)]
  }

  q <- function(p) as.numeric(stats::quantile(U, p, names = FALSE))

  # ---- point estimate and the three interval objects ----------------------
  disc_point <- mean(U)                       # Guo & He bias-reduced correction
  s_debiased <- gamma_max - disc_point

  s_one      <- gamma_max - q(1 - level)      # CERTIFIED one-sided limit
  s_two_lo   <- gamma_max - q(1 - level / 2)  # two-sided, conservative endpoint
  s_two_hi   <- gamma_max - q(level / 2)      # two-sided, optimistic endpoint

  z <- stats::qnorm(1 - level / 2)            # HEURISTIC Wald (Zhao et al. style)
  se <- stats::sd(U)
  s_wald_lo <- s_debiased - z * se
  s_wald_hi <- s_debiased + z * se

  # ---- back to the reporting scale ----------------------------------------
  # score -> coefficient -> report. For ratio scales with orient = -1 the map is
  # decreasing, so endpoints swap; sort and label by role rather than by sign.
  to_report <- function(s) {
    b <- orient * s
    if (outcome == "continuous") b else exp(b)
  }
  pair <- function(lo_s, hi_s) sort(c(to_report(lo_s), to_report(hi_s)))

  ci2 <- pair(s_two_lo, s_two_hi)
  ciw <- pair(s_wald_lo, s_wald_hi)

  scale_lab <- switch(outcome,
    survival   = "hazard ratio",
    binary     = "odds ratio",
    continuous = "mean difference"
  )
  # the conservative (toward-null) side of the one-sided limit
  bound_side <- if (outcome == "continuous") {
    "toward zero"
  } else if (orient == -1) "upper" else "lower"

  structure(
    list(
      outcome         = outcome,
      scale           = scale_lab,
      orient          = orient,
      selected        = sel,
      naive           = to_report(gamma_max),
      debiased        = to_report(s_debiased),
      bound_one_sided = to_report(s_one),
      bound_side      = bound_side,
      ci_two_sided    = ci2,
      ci_wald         = ciw,
      n               = n,
      B               = length(U),
      B_requested     = B,
      r               = r,
      level           = level,
      n_candidates    = n_cand,
      n_estimable     = sum(ok),
      fast            = isTRUE(fast),
      lean_cox        = isTRUE(lean_cox),
      adjust_covariates          = adjust_covariates,
      drop_rate       = if (diagnostics) {
        stats::setNames(drop_ct / max(1L, B), candidates)
      } else NULL
    ),
    class = "guohe_a3"
  )
}

#' Print method for guohe_a3 objects
#'
#' Prints the family size, the selected candidate, the naive and de-biased
#' effects, and the three interval objects returned by [guohe_algorithm3()],
#' labelled by their inferential status (certified / extension / heuristic).
#'
#' @param x An object of class `"guohe_a3"`.
#' @param ... Ignored; present for method consistency.
#' @return `x`, invisibly.
#' @export
print.guohe_a3 <- function(x, ...) {
  cat(sprintf("Guo-He Algorithm 3 (enumerated family)  [%s, orient = %+d]\n",
              x$outcome, x$orient))
  cat(sprintf("  family           : %d candidates (%d estimable), n = %d\n",
              x$n_candidates, x$n_estimable, x$n))
  cat(sprintf("  selected         : %s\n", x$selected))
  cat(sprintf("  naive %-16s: %.4f\n", x$scale, x$naive))
  cat(sprintf("  de-biased %-12s: %.4f\n", x$scale, x$debiased))
  cat(sprintf("  %s bound (one-sided, %.3f) : %.4f   [CERTIFIED]\n",
              x$bound_side, 1 - x$level, x$bound_one_sided))
  cat(sprintf("  two-sided CI (%.2f, quantile inversion): (%.4f, %.4f)   [EXTENSION]\n",
              1 - x$level, x$ci_two_sided[1], x$ci_two_sided[2]))
  cat(sprintf("  Wald CI (%.2f, Zhao-style)             : (%.4f, %.4f)   [HEURISTIC]\n",
              1 - x$level, x$ci_wald[1], x$ci_wald[2]))
  cat(sprintf("  B = %d, r = %.3f\n", x$B, x$r))
  if (!is.null(x$drop_rate)) {
    worst <- sort(x$drop_rate, decreasing = TRUE)
    nz <- worst[worst > 0]
    if (length(nz)) {
      cat(sprintf("  resample drop-out: %d candidate(s) non-zero; worst %s = %.3f\n",
                  length(nz), names(nz)[1], nz[[1]]))
    } else {
      cat("  resample drop-out: none\n")
    }
  }
  invisible(x)
}
