# File: R/dina_subgroup.R
# Part of the forestsearch package.
#
# Bridge between DINA's continuous per-patient HTE estimates and
# forestsearch's interpretable subgroup signature framework.  Given a
# fitted DINA object (single-pass or bagged) and a data frame, search
# over (covariate, threshold, direction) triples for the LARGEST
# subgroup whose mean tau-hat reaches a harm threshold m_diff.  This
# is the univariate analog of GRF's leaf-level subgroup identification
# and the simplest analog of forestsearch's "largest subgroup with
# consistency" rule.
#
# Variance for the subgroup-mean is a quadratic form in the DINA
# coefficient covariance.  For a fixed subgroup S with covariate-mean
# vector x_bar_S, the subgroup-mean tau-hat is a linear functional of
# the DINA coefficient vector,
#
#   mean_S(tau_hat) = a^T beta_hat,  with  a = (1, x_bar_S),
#
# so Var(mean_S(tau_hat)) = a^T Var(beta_hat) a.  This automatically
# inherits whatever variance is stored in fit$vcov: sandwich for
# `dina`/`dina_fit`, infinitesimal jackknife for `dina_bagged`/
# `dina_fit_bagged`.
#
# IMPORTANT INFERENCE CAVEAT.  The Wald CI returned here is conditional
# on the selected (covariate, threshold) being treated as pre-specified.
# It does NOT adjust for the selection across (j, q) candidates and
# will undercover when the search space is large.  A selection-adjusted
# CI requires bootstrap of the full search procedure, which is deferred.


# ---------------------------------------------------------------------------
# Main user-facing function
# ---------------------------------------------------------------------------

#' Identify a harm subgroup from a DINA fit via univariate threshold search
#'
#' Given a fitted DINA object (`"dina"` or `"dina_bagged"`) and a data
#' frame containing the covariates the model was fit on, search over
#' (covariate, threshold, direction) candidates for the *largest*
#' subgroup whose mean per-patient HTE estimate exceeds the harm
#' threshold `m_diff`.  This is the simplest bridge between DINA's
#' continuous estimator \eqn{\hat\tau(x) = \hat\beta_0 + x'\hat\beta}
#' and forestsearch's interpretable signature framework.
#'
#' For each covariate `x_j` in `covariates`, the search visits the
#' sorted unique values of `df[[x_j]]` as candidate thresholds `q` and
#' for each threshold considers the left-tail subgroup of patients with
#' \eqn{x_{j,i} \le q} (and optionally also the right-tail subgroup
#' with \eqn{x_{j,i} \ge q}).  A candidate
#' subgroup must satisfy:
#'   * subgroup size `|S| >= n_min`, AND
#'   * mean per-patient tau-hat over `S` `>= m_diff`.
#'
#' Among all qualifying candidates, the largest subgroup is returned;
#' ties on size are broken by the most extreme mean tau-hat.
#'
#' The variance of the subgroup-mean tau-hat is computed as
#' \deqn{\widehat{\mathrm{Var}}(\bar{\hat\tau}_S) =
#'        a_S^\top \widehat{\mathrm{Var}}(\hat\beta) a_S,}
#' with \eqn{a_S = (1, \bar{x}_S^\top)} and
#' \eqn{\widehat{\mathrm{Var}}(\hat\beta) = \mathrm{vcov}(\mathtt{fit})}.
#' The variance source therefore tracks the fit class: sandwich for
#' `"dina"` and infinitesimal jackknife for `"dina_bagged"`.
#'
#' @section Inference caveat:
#' The returned Wald confidence interval is *conditional on the
#' selected (covariate, threshold) being treated as pre-specified*.
#' It does not adjust for the selection across the candidate set and
#' will undercover when the search space is large.  Use this CI for
#' descriptive reporting; for hypothesis testing against `m_diff`,
#' a selection-adjusted interval via bootstrap of the full search
#' procedure is required.
#'
#' @param fit a fitted DINA object of class `"dina"` or
#'   `"dina_bagged"`, as returned by [dina_fit()], [dina_fit_bagged()],
#'   [dina()], or [dina_bagged()].
#' @param df data frame containing the covariate columns referenced by
#'   `covariates`.  Typically the same data frame used to fit `fit`.
#' @param covariates character vector of column names in `df` to search
#'   over.  All must be numeric.
#' @param m_diff scalar harm threshold on the natural-parameter scale
#'   of `fit$family`.  Subgroups are kept only if their mean tau-hat
#'   meets or exceeds this value.
#' @param n_min positive integer; minimum subgroup size.
#'   Default `60L`, matching the convention in GRF and `forestsearch()`.
#' @param direction one of `"both"` (default), `"left"`, or `"right"`.
#'   For each covariate `x_j` controls whether to search subgroups of
#'   the form `x_j <= q` (left), `x_j >= q` (right), or both.
#' @param alpha confidence level for the Wald interval on the
#'   subgroup-mean tau-hat.  Default `0.05`.
#'
#' @return An object of class `"dina_subgroup"`, a list with components:
#'   \describe{
#'     \item{found}{logical; `TRUE` if any candidate satisfied
#'       `m_diff` and `n_min`, otherwise `FALSE`.  When `FALSE`, the
#'       subgroup-specific components below are `NULL`.}
#'     \item{covariate, direction, threshold}{the chosen
#'       `(x_j, dir, q)` triple.}
#'     \item{n_subgroup}{integer size of the chosen subgroup.}
#'     \item{mean_tau_hat}{scalar subgroup-mean tau-hat.}
#'     \item{se_mean_tau_hat}{Wald standard error, computed as
#'       `sqrt(a_S^T vcov(fit) a_S)` (CONDITIONAL on the chosen
#'       subgroup -- not selection-adjusted).}
#'     \item{ci}{named length-2 vector `(lower, upper)` giving the
#'       Wald 1 - alpha CI.}
#'     \item{mask}{logical vector of length `nrow(df)` marking which
#'       rows are in the chosen subgroup.}
#'     \item{m_diff, n_min, alpha, family}{the inputs / fit family,
#'       echoed for reproducibility.}
#'     \item{n_total}{`nrow(df)`.}
#'     \item{n_candidates_searched}{total number of `(j, dir, q)`
#'       triples evaluated.}
#'     \item{n_candidates_qualifying}{number of those triples that met
#'       both the size and `m_diff` constraints.}
#'     \item{call}{the matched call.}
#'   }
#'
#' @seealso [dina_fit()] / [dina()] for the underlying estimator;
#'   [dina_fit_bagged()] / [dina_bagged()] for the bagged variant whose
#'   IJ variance propagates through the subgroup-mean SE here.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' df_demo <- data.frame(
#'   w  = stats::rbinom(n, 1, 0.5),
#'   x1 = stats::runif(n, -1, 1),
#'   x2 = stats::runif(n, -1, 1),
#'   x3 = stats::runif(n, -1, 1)
#' )
#' tau_x      <- 0.3 + 1.2 * df_demo$x1 - 0.4 * df_demo$x2
#' df_demo$y  <- 0.5 * df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)
#'
#' fit <- dina(df_demo, outcome = "y", treatment = "w",
#'             covariates = c("x1", "x2", "x3"),
#'             family = "gaussian", seed = 1L)
#'
#' sg <- dina_subgroup(fit, df_demo,
#'                     covariates = c("x1", "x2", "x3"),
#'                     m_diff = 0.5, n_min = 60L)
#' sg
#'
#' @export
dina_subgroup <- function(fit, df, covariates,
                          m_diff,
                          n_min = 60L,
                          direction = c("both", "left", "right"),
                          alpha = 0.05) {

  if (!inherits(fit, "dina")) {
    stop("`fit` must be a DINA object (class \"dina\" or \"dina_bagged\").")
  }
  if (missing(m_diff) || length(m_diff) != 1L || !is.numeric(m_diff) ||
      !is.finite(m_diff)) {
    stop("`m_diff` must be a single finite numeric value.")
  }
  direction <- match.arg(direction)
  n_min <- as.integer(n_min)
  if (length(n_min) != 1L || is.na(n_min) || n_min < 1L) {
    stop("`n_min` must be a single positive integer.")
  }
  if (length(alpha) != 1L || !is.numeric(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single numeric in (0, 1).")
  }

  X <- .dina_extract_X_from_df(df, covariates)
  n <- nrow(X)
  d <- ncol(X)

  if (n_min > n) {
    stop("`n_min` (", n_min, ") exceeds nrow(df) (", n, ").")
  }

  beta <- stats::coef(fit)
  if (length(beta) != d + 1L) {
    stop("Dimension mismatch: fit has ", length(beta) - 1L,
         " covariate coefficients but `covariates` has length ", d, ".")
  }
  V <- stats::vcov(fit)

  # Per-patient tau-hat
  tau_hat <- as.numeric(beta[1L] + X %*% beta[-1L])

  dirs_to_try <- if (direction == "both") c("left", "right") else direction

  # Cumulative tracking
  best <- NULL
  n_searched   <- 0L
  n_qualifying <- 0L

  for (j in seq_len(d)) {
    x_j <- X[, j]
    qs <- sort(unique(x_j))

    for (dir in dirs_to_try) {
      for (q in qs) {
        n_searched <- n_searched + 1L
        mask <- if (dir == "left") x_j <= q else x_j >= q
        n_S <- sum(mask)

        if (n_S < n_min) next

        mean_tau <- mean(tau_hat[mask])
        if (mean_tau < m_diff) next

        n_qualifying <- n_qualifying + 1L

        # Better than current best?  Largest subgroup; ties broken by
        # most extreme mean tau-hat.
        if (is.null(best) ||
            n_S > best$n_subgroup ||
            (n_S == best$n_subgroup && mean_tau > best$mean_tau)) {
          best <- list(
            covariate    = covariates[j],
            direction    = dir,
            threshold    = q,
            n_subgroup   = n_S,
            mean_tau     = mean_tau,
            mask         = mask
          )
        }
      }
    }
  }

  if (is.null(best)) {
    out <- list(
      found                   = FALSE,
      covariate               = NULL,
      direction               = NULL,
      threshold               = NULL,
      n_subgroup              = NULL,
      mean_tau_hat            = NULL,
      se_mean_tau_hat         = NULL,
      ci                      = NULL,
      mask                    = NULL,
      m_diff                  = m_diff,
      n_min                   = n_min,
      alpha                   = alpha,
      family                  = fit$family,
      n_total                 = n,
      n_candidates_searched   = n_searched,
      n_candidates_qualifying = n_qualifying,
      call                    = match.call()
    )
    class(out) <- "dina_subgroup"
    return(out)
  }

  # Variance of the subgroup-mean tau-hat
  x_bar_S <- colMeans(X[best$mask, , drop = FALSE])
  a       <- c(1, x_bar_S)
  var_mt  <- as.numeric(t(a) %*% V %*% a)
  se_mt   <- sqrt(max(var_mt, 0))

  z_crit <- stats::qnorm(1 - alpha / 2)
  ci_lo  <- best$mean_tau - z_crit * se_mt
  ci_hi  <- best$mean_tau + z_crit * se_mt

  out <- list(
    found                   = TRUE,
    covariate               = best$covariate,
    direction               = best$direction,
    threshold               = best$threshold,
    n_subgroup              = best$n_subgroup,
    mean_tau_hat            = best$mean_tau,
    se_mean_tau_hat         = se_mt,
    ci                      = c(lower = ci_lo, upper = ci_hi),
    mask                    = best$mask,
    m_diff                  = m_diff,
    n_min                   = n_min,
    alpha                   = alpha,
    family                  = fit$family,
    n_total                 = n,
    n_candidates_searched   = n_searched,
    n_candidates_qualifying = n_qualifying,
    call                    = match.call()
  )
  class(out) <- "dina_subgroup"
  out
}


# ---------------------------------------------------------------------------
# Internal helper: covariate-only extraction
# ---------------------------------------------------------------------------

#' Extract a numeric covariate matrix from a data frame by column name.
#'
#' Lighter-weight than `.dina_extract_from_df()`: validates `df` and
#' `covariates` and returns the covariate matrix as a numeric `matrix`
#' with `storage.mode == "double"`.  No outcome, treatment, or family
#' handling.
#'
#' @noRd
.dina_extract_X_from_df <- function(df, covariates) {

  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  }
  if (!is.character(covariates) || length(covariates) < 1L ||
      anyNA(covariates)) {
    stop("`covariates` must be a character vector of length >= 1 with no NAs.")
  }
  missing_cols <- setdiff(covariates, names(df))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `df`: ",
         paste(missing_cols, collapse = ", "), ".")
  }
  if (anyNA(df[, covariates])) {
    stop("`df` contains NA values in the covariate columns.  ",
         "Use `na.omit()` or impute before calling.")
  }

  cov_df <- df[, covariates, drop = FALSE]
  is_num <- vapply(cov_df, is.numeric, logical(1L))
  if (!all(is_num)) {
    stop("Covariate column(s) must be numeric: ",
         paste(covariates[!is_num], collapse = ", "), ".")
  }
  X <- as.matrix(cov_df)
  storage.mode(X) <- "double"
  X
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print a `dina_subgroup` result.
#'
#' @param x a `"dina_subgroup"` object.
#' @param digits number of digits for numeric summary.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @export
print.dina_subgroup <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("DINA subgroup search\n")
  cat("  Family:      ", x$family, "\n", sep = "")
  cat("  m_diff:      ", format(x$m_diff, digits = digits), "\n", sep = "")
  cat("  n_min:       ", x$n_min, "\n", sep = "")
  cat("  Candidates:  ", x$n_candidates_searched,
      " searched, ", x$n_candidates_qualifying,
      " satisfied m_diff and size\n", sep = "")

  if (!isTRUE(x$found)) {
    cat("\n  No subgroup satisfied m_diff = ",
        format(x$m_diff, digits = digits),
        " under n_min = ", x$n_min, ".\n", sep = "")
    return(invisible(x))
  }

  cmp <- if (x$direction == "left") "<=" else ">="
  cat("\n  Signature:   ", x$covariate, " ", cmp, " ",
      format(x$threshold, digits = digits), "\n", sep = "")
  cat("  Subgroup n:  ", x$n_subgroup, " / ", x$n_total,
      " (", format(100 * x$n_subgroup / x$n_total, digits = 3L),
      "%)\n", sep = "")
  cat("  mean tau-hat:", format(x$mean_tau_hat, digits = digits),
      " (SE ", format(x$se_mean_tau_hat, digits = digits), ")\n", sep = "")
  cat("  Wald ", format(100 * (1 - x$alpha), digits = 3L),
      "% CI:  [", format(x$ci[["lower"]], digits = digits), ", ",
      format(x$ci[["upper"]], digits = digits), "]\n", sep = "")
  cat("\n",
      "  Note: CI is conditional on the selected (covariate, threshold)\n",
      "  and does NOT adjust for the selection across candidates.\n",
      sep = "")
  invisible(x)
}
