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
# Depends on {survival} and base R only. Not part of the package surface.
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
.g3_coef <- function(sub, outcome, treatment, time, event, y, min_events) {
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
        fit <- survival::coxph(survival::Surv(sub[[time]], ev) ~ tr)
        unname(fit$coefficients[[1L]])
      } else {
        yy <- sub[[y]]
        if (outcome == "binary" && length(unique(yy)) < 2L) {
          return(NA_real_)
        }
        fam <- if (outcome == "binary") stats::binomial() else stats::gaussian()
        unname(stats::coef(stats::glm(yy ~ tr, family = fam))[["tr"]])
      }
    }),
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
#'
#' @return An object of class `"guohe_a3"`.
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
                             diagnostics = TRUE) {
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

  # oriented score for every candidate on a given row-index set
  score_vec <- function(idx) {
    d <- data[idx, , drop = FALSE]
    m <- membership[idx, , drop = FALSE]
    vapply(
      seq_len(n_cand),
      function(k) {
        sub <- d[m[[k]] == 1L, , drop = FALSE]
        orient * .g3_coef(sub, outcome, treatment, time, event, y, min_events)
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
  if (!is.null(seed)) set.seed(seed)
  U <- numeric(B)
  drop_ct <- integer(n_cand)
  for (b_i in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    val <- score_vec(idx) + offset
    fin <- is.finite(val)
    if (diagnostics) drop_ct <- drop_ct + as.integer(!fin)
    if (!any(fin)) {
      U[b_i] <- NA_real_
      next
    }
    U[b_i] <- max(val[fin]) - gamma_max
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
      drop_rate       = if (diagnostics) {
        stats::setNames(drop_ct / max(1L, B), candidates)
      } else NULL
    ),
    class = "guohe_a3"
  )
}

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
