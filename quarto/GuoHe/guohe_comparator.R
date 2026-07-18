# guohe_comparator.R
#
# De-biased inference for the best selected subgroup (Guo & He, 2021),
# generalized to survival (hazard ratio), binary (odds ratio), and continuous
# (mean difference) outcomes.
#
# Standalone FS-side comparator for `forestsearch` / `fs-glms-interpretable`.
# Depends only on {survival} and base R; it is NOT part of the package surface
# and imports nothing from `forestsearch`. It takes a supplied fixed family of
# candidate subgroups (as membership-indicator columns) and returns the naive
# best-subgroup effect, the bias-reduced point estimate, and a one-sided
# confidence bound toward the null.
#
# With `outcome = "survival"`, `orient = -1`, `symmetric = TRUE` it reproduces
# the survival reference `synthetic code.R` (MONET1) exactly.
#
# Reference: Guo, X. and He, X. (2021). Inference on selected subgroups in
# clinical trials. Journal of the American Statistical Association,
# 116(535):1498-1506.

# ---- internal: one within-subgroup treatment coefficient ------------------
# Returns the treatment coefficient b on the model's natural scale
# (log-HR / log-OR / mean difference), or NA_real_ if not estimable.
.guohe_coef <- function(sub, outcome, treatment, time, event, y, min_events) {
  tr <- sub[[treatment]]
  if (length(unique(tr)) < 2L) {
    return(NA_real_)
  }
  if (min(table(tr)) < min_events) {
    return(NA_real_)
  }
  out <- tryCatch(
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
        fit <- stats::glm(yy ~ tr, family = fam)
        unname(stats::coef(fit)[["tr"]])
      }
    }),
    error = function(e) NA_real_
  )
  out
}

#' De-biased inference for the best selected subgroup, generalized
#'
#' Model-free, fixed-family effect-debiasing comparator for the winner's-curse
#' bias of a data-driven "best" subgroup, after Guo and He (2021). Given a
#' supplied candidate family, the subgroup with the most extreme oriented effect
#' is selected and its selection optimism removed by the pair (case) bootstrap
#' with deterministic shrinkage offsets. Survival, binary, and continuous
#' outcomes are supported through a pluggable within-subgroup coefficient
#' (Cox log-hazard-ratio, logistic log-odds-ratio, Gaussian mean difference).
#'
#' The procedure maximizes an oriented benefit score `score = orient * b`, where
#' `b` is the within-subgroup treatment coefficient on the model's natural
#' (log or identity) scale. `orient = -1` makes a protective ratio (HR or OR
#' below 1) the effect of interest (the Guo and He / MONET1 convention);
#' `orient = +1` makes a harmful ratio (above 1) the effect of interest. For a
#' continuous outcome `orient = +1` maximizes the mean difference and
#' `orient = -1` minimizes it. Reported effects are `exp(b)` for the ratio
#' scales (survival, binary) and `b` for the mean-difference scale (continuous).
#'
#' The returned bound is the conservative one-sided limit toward the null at
#' level `level` (an upper HR/OR bound when `orient = -1`, a lower HR/OR bound
#' when `orient = +1`, the toward-zero limit for a mean difference); it is not a
#' two-sided interval. `symmetric = TRUE` expands each supplied indicator into
#' the subgroup `{col == 1}` and its complement `{col == 0}` (the
#' `synthetic code.R` behaviour); `symmetric = FALSE` treats each supplied
#' column as a single candidate, the natural choice for an enumerated
#' forest-search family.
#'
#' @param data A data frame.
#' @param outcome One of `"survival"`, `"binary"`, `"continuous"`.
#' @param treatment Name of the 0/1 treatment column.
#' @param candidates Character vector of names of 0/1 subgroup-membership
#'   columns forming the fixed candidate family.
#' @param time,event Survival only: names of the time and event columns
#'   (event = 1 for the event, 0 for censoring).
#' @param y Binary or continuous only: name of the outcome column.
#' @param orient Either `-1` or `+1`; sign orienting the effect of interest.
#' @param B Number of bootstrap resamples.
#' @param r Shrinkage tuning parameter, strictly between 0 and 0.5.
#' @param level One-sided level of the bound.
#' @param seed Optional integer for reproducible resampling.
#' @param min_events Minimum events (survival) or minimum per-arm count
#'   (binary, continuous) for a candidate to be fit; smaller candidates are
#'   dropped from the competition.
#' @param symmetric Logical; expand each candidate into subgroup and complement.
#'
#' @return An object of class `"guohe_debias"`: a list with the selected
#'   candidate, the `naive`, `debiased`, and one-sided `bound` effects on the
#'   reporting scale, plus run metadata.
#' @export
guohe_debias <- function(data,
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
                         symmetric = FALSE) {
  outcome <- match.arg(outcome)
  stopifnot(
    is.data.frame(data),
    orient %in% c(-1, 1),
    r > 0, r < 0.5,
    level > 0, level < 1,
    B >= 1L
  )
  need <- c(treatment, candidates)
  if (outcome == "survival") need <- c(need, time, event) else need <- c(need, y)
  missing_cols <- setdiff(need, names(data))
  if (length(missing_cols)) {
    stop("Columns not found in `data`: ", paste(missing_cols, collapse = ", "))
  }
  n <- nrow(data)

  # ---- build the fixed candidate membership matrix ------------------------
  memb <- lapply(candidates, function(cn) as.integer(data[[cn]] == 1))
  memb <- as.data.frame(memb)
  if (symmetric) {
    comp <- lapply(candidates, function(cn) as.integer(data[[cn]] == 0))
    comp <- as.data.frame(comp)
    names(memb) <- paste0(candidates, "{=1}")
    names(comp) <- paste0(candidates, "{=0}")
    membership <- cbind(memb, comp)
  } else {
    names(memb) <- candidates
    membership <- memb
  }
  cand_names <- names(membership)
  n_cand <- ncol(membership)

  # ---- oriented scores for a given row-index set --------------------------
  score_vec <- function(idx) {
    d <- data[idx, , drop = FALSE]
    m <- membership[idx, , drop = FALSE]
    vapply(
      seq_len(n_cand),
      function(k) {
        sub <- d[m[[k]] == 1L, , drop = FALSE]
        orient * .guohe_coef(sub, outcome, treatment, time, event, y, min_events)
      },
      numeric(1)
    )
  }

  # ---- observed competition ----------------------------------------------
  obs <- score_vec(seq_len(n))
  ok <- is.finite(obs)
  if (!any(ok)) {
    stop("No estimable candidate on the full sample.")
  }
  if (!all(ok)) {
    warning(sum(!ok), " of ", n_cand,
            " candidates not estimable on the full sample; dropped.")
  }
  mm <- max(obs[ok])
  offset <- (sqrt(n) - n^r) * (mm - obs) / sqrt(n)

  # ---- bootstrap ----------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  # bs holds the *unscaled* bootstrap discrepancy max(score + offset) - mm.
  # Guo and He's point correction (CB) uses its mean; the one-sided bound (BC)
  # uses its (1 - level) quantile. The sqrt(n) scaling in BC's original code
  # cancels against the sqrt(n) in its quantile, so both live on this scale.
  bs <- numeric(B)
  for (b_i in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    val <- score_vec(idx) + offset
    bs[b_i] <- max(val[is.finite(val)]) - mm
  }

  disc_point <- mean(bs)
  disc_bound <- as.numeric(stats::quantile(bs, 1 - level))

  to_report <- function(score) {
    b <- orient * score
    if (outcome == "continuous") b else exp(b)
  }

  scale_lab <- switch(outcome,
    survival   = "hazard ratio",
    binary     = "odds ratio",
    continuous = "mean difference"
  )
  bound_side <- if (outcome == "continuous") {
    "toward zero"
  } else if (orient == -1) {
    "upper"
  } else {
    "lower"
  }

  structure(
    list(
      outcome      = outcome,
      scale        = scale_lab,
      orient       = orient,
      selected     = cand_names[which.max(replace(obs, !ok, -Inf))],
      naive        = to_report(mm),
      debiased     = to_report(mm - disc_point),
      bound        = to_report(mm - disc_bound),
      bound_side   = bound_side,
      n            = n,
      B            = B,
      r            = r,
      level        = level,
      n_candidates = n_cand,
      n_estimable  = sum(ok)
    ),
    class = "guohe_debias"
  )
}

#' @export
print.guohe_debias <- function(x, ...) {
  cat(sprintf("Guo-He de-biased best subgroup  (%s, orient = %+d)\n",
              x$outcome, x$orient))
  cat(sprintf("  candidates       : %d (%d estimable)\n",
              x$n_candidates, x$n_estimable))
  cat(sprintf("  selected         : %s\n", x$selected))
  cat(sprintf("  naive %-11s: %.4f\n", x$scale, x$naive))
  cat(sprintf("  de-biased %-7s: %.4f\n", x$scale, x$debiased))
  cat(sprintf("  %s bound (1-%.2f) : %.4f\n", x$bound_side, x$level, x$bound))
  cat(sprintf("  B = %d, r = %.3f, n = %d\n", x$B, x$r, x$n))
  invisible(x)
}
