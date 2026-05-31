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
#' (covariate, threshold, direction) candidates for subgroups whose mean
#' per-patient HTE estimate meets the harm threshold `m_diff`, then select
#' one according to `sg_focus` (default `"maxSG"`: the largest qualifying
#' subgroup).  This is the simplest bridge between DINA's continuous
#' estimator \eqn{\hat\tau(x) = \hat\beta_0 + x'\hat\beta} and
#' forestsearch's interpretable signature framework.
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
#' Among the qualifying candidates, `sg_focus` selects the returned
#' subgroup, mirroring [forestsearch()]'s selection vocabulary:
#'   * `"maxSG"` (default) -- largest qualifying subgroup (legacy
#'     behaviour; size ties broken by the most extreme mean tau-hat).
#'   * `"minSG"` -- smallest qualifying subgroup.
#'   * `"eff"` -- most extreme effect, ignoring size.
#'   * `"effMaxSG"` -- largest subgroup in the candidate inclusion band
#'     (see `selection_rule`).
#'   * `"effMinSG"` -- smallest subgroup in the inclusion band.
#' The `"eff*"` band foci reuse the inclusion-band logic shared with
#' [sort_subgroups()]; for ratio families (binomial, poisson, cox) the
#' band is applied on the natural (exponentiated) effect scale.
#' `selection_rule` defines the band, exactly as in [forestsearch()]:
#' `"neighborhood"` (1-D effect band within `effect_neighborhood` of the
#' maximum qualifying effect), `"pareto"` (the 2-D non-dominated frontier
#' in (effect, N) space), or `"both"` (their intersection).  DINA's
#' univariate (covariate, threshold, direction) candidates populate the
#' (effect, N) plane just as forestsearch's signatures do, so the
#' frontier is defined analogously.
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
#' @param sg_focus character; subgroup selection criterion applied to the
#'   qualifying candidates.  One of `"maxSG"` (default), `"minSG"`,
#'   `"eff"`, `"effMaxSG"`, `"effMinSG"`.  The canonical forms `"hr"`,
#'   `"hrMaxSG"`, `"hrMinSG"` (used internally by [forestsearch()]) are
#'   also accepted.  See Description for the semantics of each.
#' @param selection_rule character; rule defining the candidate inclusion
#'   band for `"effMaxSG"` / `"effMinSG"`.  One of `"neighborhood"`
#'   (default; 1-D effect band), `"pareto"` (2-D non-dominated frontier
#'   in (effect, N)), or `"both"` (intersection).  Must be
#'   `"neighborhood"` for the other foci.  Forwarded unchanged to the
#'   shared band computation, so semantics match [forestsearch()].
#' @param effect_neighborhood numeric in `[0, 1)`; relative tolerance for
#'   the effect band.  A candidate is in the neighborhood iff its
#'   (natural-scale) effect is at least
#'   `(1 - effect_neighborhood) * max(effect)` over the qualifying set.
#'   Default `0.10` (within 10% of the strongest effect).  Consulted by
#'   `selection_rule = "neighborhood"` and `"both"`; ignored by
#'   `"pareto"` (and must be left at its default there) and for
#'   `"maxSG"`, `"minSG"`, `"eff"`.
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
#'     \item{m_diff, n_min, sg_focus, selection_rule,
#'       effect_neighborhood, alpha, family}{the inputs / fit family,
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
#'   IJ variance propagates through the subgroup-mean SE here;
#'   [sort_subgroups()] and [forestsearch()] for the selection
#'   vocabulary mirrored by `sg_focus`.
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
#' # Tighter subgroup via the effect-neighborhood band:
#' dina_subgroup(fit, df_demo, covariates = c("x1", "x2", "x3"),
#'               m_diff = 0.5, n_min = 60L, sg_focus = "effMaxSG")
#'
#' # Frontier-aware selection (2-D non-dominated set in (effect, N)):
#' dina_subgroup(fit, df_demo, covariates = c("x1", "x2", "x3"),
#'               m_diff = 0.5, n_min = 60L, sg_focus = "effMaxSG",
#'               selection_rule = "both")
#'
#' @export
dina_subgroup <- function(fit, df, covariates,
                          m_diff,
                          n_min = 60L,
                          direction = c("both", "left", "right"),
                          sg_focus = "maxSG",
                          selection_rule = "neighborhood",
                          effect_neighborhood = 0.10,
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

  # Normalize the GLM-natural vocabulary ("eff", "effMaxSG", "effMinSG")
  # to the canonical internal form ("hr", "hrMaxSG", "hrMinSG") shared
  # with forestsearch(), then whitelist.
  sg_focus <- .normalize_sg_focus(sg_focus)
  valid_sg_focus <- c("hr", "maxSG", "minSG", "hrMaxSG", "hrMinSG")
  if (!is.character(sg_focus) || length(sg_focus) != 1L ||
      !sg_focus %in% valid_sg_focus) {
    stop("`sg_focus` must be one of \"maxSG\", \"minSG\", \"eff\", ",
         "\"effMaxSG\", \"effMinSG\" (canonical forms \"hr\", ",
         "\"hrMaxSG\", \"hrMinSG\" are also accepted).")
  }
  # Always range-check; only consulted for the band foci below.
  .validate_effect_neighborhood(effect_neighborhood)
  # selection_rule must be one of "neighborhood"/"pareto"/"both", is only
  # meaningful for the band foci, and "pareto" forbids a non-default
  # effect_neighborhood -- same rules forestsearch() enforces.
  .validate_selection_rule(selection_rule, sg_focus, effect_neighborhood)

  # DINA's tau-hat is the link-scale linear predictor: identity (MD) for
  # gaussian, log-OR for binomial, log-IRR for poisson, log-HR for cox.
  # Mirror forestsearch by applying the neighborhood band on the natural
  # scale for ratio families (exponentiate before the inclusion test).
  effect_log_scale <- isTRUE(fit$family %in% c("binomial", "poisson", "cox"))

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

  # Collect the qualifying candidate cloud (size >= n_min AND mean
  # tau-hat >= m_diff).  Shared with dina_frontier() so both define
  # "qualifying candidate" identically.  The full set is required so the
  # effect band (sg_focus %in% c("hrMaxSG", "hrMinSG")) can be anchored
  # to the maximum effect over all qualifying candidates.
  cl <- .dina_collect_candidates(X = X, tau_hat = tau_hat,
                                 covariates = covariates, m_diff = m_diff,
                                 n_min = n_min, direction = direction)
  n_searched   <- cl$n_searched
  n_qualifying <- cl$n_qualifying

  if (n_qualifying == 0L) {
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
      sg_focus                = sg_focus,
      selection_rule          = selection_rule,
      effect_neighborhood     = effect_neighborhood,
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

  # Trim preallocated vectors to the qualifying candidates.
  cand_j   <- cl$j
  cand_dir <- cl$direction
  cand_q   <- cl$threshold
  cand_n   <- cl$n_subgroup
  cand_tau <- cl$mean_tau

  # Effect on the natural scale for the band / ordering (ratio families
  # are exponentiated; the transform is monotone, so the ordering of the
  # non-band foci is unaffected).  `idx` is insertion order (covariate,
  # then left/right, then ascending threshold), used as the final
  # tiebreak so the default "maxSG" reproduces the legacy running-best
  # selection exactly.
  eff <- if (effect_log_scale) exp(cand_tau) else cand_tau
  idx <- seq_along(cand_tau)

  ord <- switch(
    sg_focus,
    maxSG   = order(-cand_n, -eff, idx),
    minSG   = order( cand_n, -eff, idx),
    hr      = order(-eff, idx),
    hrMaxSG = {
      in_band <- .compute_inclusion_band(
        hr_vec              = eff,
        n_vec               = cand_n,
        selection_rule      = selection_rule,
        effect_neighborhood = effect_neighborhood
      )
      order(-in_band, -cand_n, -eff, idx)
    },
    hrMinSG = {
      in_band <- .compute_inclusion_band(
        hr_vec              = eff,
        n_vec               = cand_n,
        selection_rule      = selection_rule,
        effect_neighborhood = effect_neighborhood
      )
      order(-in_band, cand_n, -eff, idx)
    }
  )

  w <- ord[1L]
  best <- list(
    covariate  = covariates[cand_j[w]],
    direction  = cand_dir[w],
    threshold  = cand_q[w],
    n_subgroup = cand_n[w],
    mean_tau   = cand_tau[w]
  )
  x_win     <- X[, cand_j[w]]
  best$mask <- if (best$direction == "left") {
    x_win <= best$threshold
  } else {
    x_win >= best$threshold
  }

  # Variance of the subgroup-mean tau-hat (link scale, as reported).
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
    sg_focus                = sg_focus,
    selection_rule          = selection_rule,
    effect_neighborhood     = effect_neighborhood,
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
# Internal helper: collect the qualifying candidate cloud
# ---------------------------------------------------------------------------

#' Collect qualifying (covariate, direction, threshold) candidates.
#'
#' Shared by \code{dina_subgroup()} (which then ranks and selects one)
#' and \code{dina_frontier()} (which then extracts the non-dominated
#' set), so both functions define "qualifying candidate" identically.
#' Visits the sorted unique values of each covariate as candidate
#' thresholds and keeps every (covariate, direction, threshold) whose
#' subgroup has size \code{>= n_min} and mean tau-hat \code{>= m_diff}.
#' Pass \code{m_diff = -Inf} to keep all size-qualifying candidates
#' (no effect floor).
#'
#' @param X numeric covariate matrix (rows = subjects).
#' @param tau_hat numeric per-subject tau-hat (length \code{nrow(X)}).
#' @param covariates character column names, aligned with \code{X}.
#' @param m_diff scalar effect floor on the link (tau) scale; may be
#'   \code{-Inf}.
#' @param n_min positive integer minimum subgroup size.
#' @param direction one of \code{"both"}, \code{"left"}, \code{"right"}.
#'
#' @return A list with vectors \code{covariate}, \code{j} (column index),
#'   \code{direction}, \code{threshold}, \code{n_subgroup},
#'   \code{mean_tau} (link scale), each of length \code{n_qualifying},
#'   plus scalars \code{n_searched} and \code{n_qualifying}.  Candidates
#'   are in insertion order (covariate, then left/right, then ascending
#'   threshold).
#'
#' @noRd
.dina_collect_candidates <- function(X, tau_hat, covariates, m_diff, n_min,
                                     direction) {
  d <- ncol(X)
  dirs_to_try <- if (direction == "both") c("left", "right") else direction

  # Upper bound on the number of (j, dir, q) candidates, for preallocation.
  n_cand_max <- sum(vapply(seq_len(d),
                           function(j) length(unique(X[, j])),
                           integer(1L))) * length(dirs_to_try)

  cand_j   <- integer(n_cand_max)
  cand_dir <- character(n_cand_max)
  cand_q   <- numeric(n_cand_max)
  cand_n   <- integer(n_cand_max)
  cand_tau <- numeric(n_cand_max)

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
        cand_j[n_qualifying]   <- j
        cand_dir[n_qualifying] <- dir
        cand_q[n_qualifying]   <- q
        cand_n[n_qualifying]   <- n_S
        cand_tau[n_qualifying] <- mean_tau
      }
    }
  }

  keep <- seq_len(n_qualifying)
  list(
    covariate    = covariates[cand_j[keep]],
    j            = cand_j[keep],
    direction    = cand_dir[keep],
    threshold    = cand_q[keep],
    n_subgroup   = cand_n[keep],
    mean_tau     = cand_tau[keep],
    n_searched   = n_searched,
    n_qualifying = n_qualifying
  )
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
  # Map canonical sg_focus to a readable label (guarded for objects
  # created before sg_focus was added).
  sgf <- if (is.null(x$sg_focus)) "maxSG" else x$sg_focus
  focus_label <- switch(
    sgf,
    hr      = "eff (most extreme effect)",
    maxSG   = "maxSG (largest qualifying subgroup)",
    minSG   = "minSG (smallest qualifying subgroup)",
    hrMaxSG = "effMaxSG (largest within effect band)",
    hrMinSG = "effMinSG (smallest within effect band)",
    sgf
  )
  uses_band <- isTRUE(sgf %in% c("hrMaxSG", "hrMinSG"))

  cat("DINA subgroup search\n")
  cat("  Family:      ", x$family, "\n", sep = "")
  cat("  Focus:       ", focus_label, "\n", sep = "")
  if (uses_band) {
    en <- if (is.null(x$effect_neighborhood)) 0.10 else x$effect_neighborhood
    sr <- if (is.null(x$selection_rule)) "neighborhood" else x$selection_rule
    band_desc <- switch(
      sr,
      neighborhood = sprintf("within %s%% of max effect",
                             format(100 * en, digits = digits)),
      pareto       = "Pareto frontier on (effect, N)",
      both         = sprintf("within %s%% of max effect AND Pareto frontier",
                             format(100 * en, digits = digits)),
      sr
    )
    cat("  Inclusion:   ", band_desc, "\n", sep = "")
  }
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


# ---------------------------------------------------------------------------
# Pareto-frontier extractor (forestsearch candidate generation)
# ---------------------------------------------------------------------------

#' Validate a cap argument that may be a positive integer or Inf.
#'
#' @param x candidate value.
#' @param nm argument name, for the error message.
#' @return `Inf` unchanged, or `x` coerced to a positive integer.
#' @noRd
.dina_check_cap <- function(x, nm) {
  ok <- length(x) == 1L && is.numeric(x) && !is.na(x) &&
    (is.infinite(x) && x > 0 ||
       (is.finite(x) && x >= 1 && x == round(x)))
  if (!ok) {
    stop("`", nm, "` must be a single positive integer or Inf.")
  }
  if (is.infinite(x)) Inf else as.integer(x)
}

#' Extract DINA per-covariate Pareto frontiers as forestsearch cuts
#'
#' Runs the same univariate (covariate, direction, threshold) search as
#' [dina_subgroup()], but instead of selecting one subgroup it returns a
#' pool of **per-covariate Pareto-non-dominated** cuts.  Each cut is
#' rendered as a forestsearch cut expression (`"x1 <= 0.5"`), so the pool
#' can be fed to [forestsearch()] as screening-stage candidates in the
#' same way GRF tree cuts are.  This is the DINA analog of GRF candidate
#' discovery.
#'
#' @details
#' **Why per-covariate and not a single global frontier.**  The global
#' Pareto frontier in (effect, N) space collapses onto whichever single
#' covariate dominates the trade-off boundary, returning many
#' micro-stepped cuts on that one covariate.  That is faithful to the
#' definition but the opposite of what forestsearch needs: its value is
#' in *composing* cuts drawn from several covariates.  So this function
#' computes a separate frontier for each covariate, then pools them,
#' mirroring how a GRF tree splits on multiple variables.
#'
#' **Two caps, each optional.**  Diversity and total burden are
#' controlled by two independent caps, either of which may be `Inf` to
#' disable it:
#' \itemize{
#'   \item `max_per_covariate` bounds how many cuts any single covariate
#'     contributes (its top cuts by effect), limiting within-covariate
#'     redundancy.
#'   \item `max_subgroups` bounds the total pool size across covariates,
#'     the DINA analog of forestsearch's `max_subgroups_search`.
#' }
#' When the pool exceeds `max_subgroups`, the global trim is round-robin
#' on within-covariate effect rank -- every covariate's best cut first
#' (ordered by effect), then every covariate's second-best, and so on --
#' so a finite budget preserves cross-covariate spread instead of
#' refilling with one covariate's micro-steps.  Setting
#' `max_per_covariate = Inf` recovers a pure global budget; setting
#' `max_subgroups = Inf` recovers per-covariate-only limits.
#'
#' **Cut form.**  The expression is always the canonical `"<covariate>
#' <= <threshold>"`.  forestsearch's factor machinery exposes both the
#' cut and its complement, so one expression covers both of DINA's
#' left/right directions -- the `direction` column is retained for
#' reference only.  Thresholds are rounded to `digits` significant
#' figures so they deduplicate cleanly against existing GRF / median cuts
#' and read well; raw values remain in `threshold`.
#'
#' These cuts are **candidates**, not forced selections: forestsearch
#' composes (up to `maxk`) and consistency-gatekeeps them, exactly as it
#' does GRF cuts.  To use them without disturbing any user-supplied
#' forced cuts, append rather than overwrite, e.g.
#' `forestsearch(..., conf_force = c(my_forced_cuts, fr$cut_expr))`.
#'
#' @param fit a fitted DINA object (class `"dina"` or `"dina_bagged"`).
#' @param df data frame containing the covariate columns.
#' @param covariates character vector of numeric covariate columns to
#'   search over.
#' @param scope one of `"wide"` (default) or `"harm"`.  `"wide"` casts a
#'   broad net: every candidate with size `>= n_min` is eligible, and the
#'   per-covariate frontiers span the whole effect-size boundary
#'   (mirroring GRF's wide-net role, leaving forestsearch's consistency
#'   stage as the gatekeeper).  `"harm"` restricts the eligible set to
#'   candidates with mean tau-hat `>= m_diff` before taking the frontiers.
#' @param m_diff scalar harm floor on the link (natural-parameter) scale.
#'   Required when `scope = "harm"`; ignored when `scope = "wide"`.
#'   Default `NULL`.
#' @param n_min positive integer minimum subgroup size.  Default `60L`.
#' @param direction one of `"both"` (default), `"left"`, `"right"`.
#' @param max_per_covariate positive integer, or `Inf` for no limit: the
#'   most cuts any single covariate may contribute (its top cuts by
#'   effect).  Default `3L`.
#' @param max_subgroups positive integer, or `Inf` for no limit: the cap
#'   on the total pooled cuts returned -- the DINA analog of
#'   forestsearch's `max_subgroups_search`.  When the pool is larger it is
#'   trimmed round-robin by within-covariate effect rank.  Default `10L`.
#' @param digits significant figures for rounding the emitted threshold
#'   in `cut_expr`.  Default `3L`.
#'
#' @return A data frame (one row per cut, after the caps and
#'   threshold-rounding deduplication) with columns:
#'   \describe{
#'     \item{covariate, direction, threshold}{the subgroup definition
#'       (raw threshold).}
#'     \item{n_subgroup}{subgroup size.}
#'     \item{mean_tau_hat}{subgroup-mean tau-hat on the link scale.}
#'     \item{effect}{the same on the natural scale (exponentiated for
#'       ratio families binomial/poisson/cox), used for the frontiers and
#'       the effect ranking.}
#'     \item{cut_expr}{the canonical forestsearch cut expression.}
#'   }
#'   The result carries attributes `n_frontier` (total per-covariate
#'   frontier size, deduped, before the caps), `n_qualifying`, and
#'   `scope`.  Returns a 0-row data frame (same columns) when no candidate
#'   qualifies.
#'
#' @seealso [dina_subgroup()] for the single-subgroup selector;
#'   [forestsearch()] for the screening stage that consumes the cuts;
#'   [compute_pareto_frontier()] for the post-hoc frontier reporting on a
#'   forestsearch result.
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
#' fr <- dina_frontier(fit, df_demo, covariates = c("x1", "x2", "x3"),
#'                     n_min = 60L, max_per_covariate = 3L,
#'                     max_subgroups = 10L)
#' fr
#' fr$cut_expr   # ready to append to forestsearch(conf_force = ...)
#'
#' @export
dina_frontier <- function(fit, df, covariates,
                          scope = c("wide", "harm"),
                          m_diff = NULL,
                          n_min = 60L,
                          direction = c("both", "left", "right"),
                          max_per_covariate = 3L,
                          max_subgroups = 10L,
                          digits = 3L) {

  if (!inherits(fit, "dina")) {
    stop("`fit` must be a DINA object (class \"dina\" or \"dina_bagged\").")
  }
  scope     <- match.arg(scope)
  direction <- match.arg(direction)
  n_min <- as.integer(n_min)
  if (length(n_min) != 1L || is.na(n_min) || n_min < 1L) {
    stop("`n_min` must be a single positive integer.")
  }
  # Both caps accept a positive integer OR Inf (Inf = no limit).  as.integer
  # would turn Inf into NA, so validate without coercing.
  max_per_covariate <- .dina_check_cap(max_per_covariate, "max_per_covariate")
  max_subgroups     <- .dina_check_cap(max_subgroups, "max_subgroups")
  if (length(digits) != 1L || is.na(digits) || digits < 1L) {
    stop("`digits` must be a single positive integer.")
  }

  # Resolve the effect floor from scope.  "wide" => no floor (-Inf).
  if (scope == "harm") {
    if (is.null(m_diff) || length(m_diff) != 1L || !is.numeric(m_diff) ||
        !is.finite(m_diff)) {
      stop("`scope = \"harm\"` requires a single finite `m_diff` ",
           "(the harm floor on the link scale).")
    }
    floor <- m_diff
  } else {
    floor <- -Inf
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
  tau_hat <- as.numeric(beta[1L] + X %*% beta[-1L])
  effect_log_scale <- isTRUE(fit$family %in% c("binomial", "poisson", "cox"))

  cl <- .dina_collect_candidates(X = X, tau_hat = tau_hat,
                                 covariates = covariates, m_diff = floor,
                                 n_min = n_min, direction = direction)

  empty <- data.frame(
    covariate = character(0), direction = character(0),
    threshold = numeric(0), n_subgroup = integer(0),
    mean_tau_hat = numeric(0), effect = numeric(0),
    cut_expr = character(0), stringsAsFactors = FALSE
  )

  if (cl$n_qualifying == 0L) {
    attr(empty, "n_frontier")   <- 0L
    attr(empty, "n_qualifying") <- 0L
    attr(empty, "scope")        <- scope
    return(empty)
  }

  # Natural-scale effect (ratio families exponentiated) for the frontier.
  eff <- if (effect_log_scale) exp(cl$mean_tau) else cl$mean_tau

  cand <- data.frame(
    covariate    = cl$covariate,
    direction    = cl$direction,
    threshold    = cl$threshold,
    n_subgroup   = as.integer(cl$n_subgroup),
    mean_tau_hat = cl$mean_tau,
    effect       = eff,
    stringsAsFactors = FALSE
  )

  # Per-covariate Pareto frontier.  The GLOBAL frontier collapses onto
  # whichever covariate dominates the (effect, N) boundary, so instead we
  # take each covariate's own non-dominated set; this yields the
  # cross-covariate diversity forestsearch composes over.  Within each
  # covariate: rank by effect (strongest first), render the canonical
  # "var <= t" cut, collapse rounding-equal micro-steps, then apply the
  # per-covariate cap.  `cov_rank` records the within-covariate effect
  # rank for the round-robin global trim below.
  per_cov <- lapply(split(cand, cand$covariate), function(cc) {
    dom <- .pareto_dominated_xy(cc$effect, cc$n_subgroup)
    ff  <- cc[!dom, , drop = FALSE]
    ff  <- ff[order(-ff$effect, -ff$n_subgroup), , drop = FALSE]
    # Render the cut from each row's OWN direction: "left" => x <= q (low
    # values are the subgroup), "right" => x >= q (high values).  Hardcoding
    # "<=" here mislabels "right" candidates and makes cut_expr contradict
    # both the direction column and the label dina_subgroup() selects.
    op <- ifelse(ff$direction == "left", "<=", ">=")
    ff$cut_expr <- paste0(ff$covariate, " ", op, " ", signif(ff$threshold, digits))
    ff[!duplicated(ff$cut_expr), , drop = FALSE]
  })

  # Raw per-covariate frontier size (deduped) before any cap.
  n_frontier <- sum(vapply(per_cov, nrow, integer(1L)))

  per_cov <- lapply(per_cov, function(ff) {
    if (nrow(ff) > max_per_covariate) {
      ff <- ff[seq_len(max_per_covariate), , drop = FALSE]
    }
    ff$cov_rank <- seq_len(nrow(ff))
    ff
  })
  fr <- do.call(rbind, per_cov)

  # Global cap via round-robin on within-covariate effect rank: every
  # covariate's best cut first (ordered by effect), then every covariate's
  # second-best, and so on -- so a finite budget preserves diversity
  # rather than refilling with one covariate's micro-steps.
  fr <- fr[order(fr$cov_rank, -fr$effect), , drop = FALSE]
  if (nrow(fr) > max_subgroups) {
    fr <- fr[seq_len(max_subgroups), , drop = FALSE]
  }
  fr$cov_rank <- NULL
  rownames(fr) <- NULL

  attr(fr, "n_frontier")   <- n_frontier
  attr(fr, "n_qualifying") <- cl$n_qualifying
  attr(fr, "scope")        <- scope
  fr
}
