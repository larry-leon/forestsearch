# =============================================================================
# fs_dgm_scale() -- finite-population sampling scale of the difference-in-means
#                   estimator, computed from a DGM's potential outcomes
# =============================================================================
#
# Replaces the practice of anchoring the analytic scale on a single MEASURED
# standard error.  Every quantity here is enumerable on the super-population, so
# the scale becomes a computed constant of the fixture rather than an estimate.
#
# The failure this prevents: a measured anchor of 13.6786 at n = 1000 implied an
# effective per-subject variance of 16,119 for the ACTG175 CD4-change fixture,
# which is BELOW that fixture's residual variance of 16,256.  Since the bracket
# on the true region equals sigma^2 + V_Q[mu_0] with V_Q[mu_0] >= 0, the
# anchored value was not merely imprecise but structurally impossible.  A
# computed scale cannot express that state.
# =============================================================================


#' Sampling Scale of the Difference-in-Means Estimator
#'
#' Computes the finite-population sampling scale of the within-region
#' difference-in-means estimator directly from a data generating mechanism's
#' individual-level potential outcomes.
#'
#' For a region \eqn{g} the estimator is the within-region difference in arm
#' means.  Conditioning on the realized arm counts,
#' \deqn{\mathrm{Var}[\hat\beta(g)] = V_{\mathrm{eff}}(g) / (n P(g)),}
#' with the region- and \eqn{n}-free constant
#' \deqn{V_{\mathrm{eff}}(g) = V_1(g)/p + V_0(g)/(1-p),}
#' where \eqn{p} is the randomisation probability and \eqn{V_w(g)} is the total
#' within-arm outcome variance in \eqn{g} under arm \eqn{w}:
#' \deqn{V_w(g) = m_g[v_w] + V_g[\mu_w].}
#' Here \eqn{\mu_w} is the arm mean surface, \eqn{v_w(x) = \mathrm{Var}(Y \mid
#' X = x, W = w)} the conditional outcome variance, and \eqn{m_g[\cdot]},
#' \eqn{V_g[\cdot]} finite-population moments over the rows of the
#' super-population lying in \eqn{g}.  The conditional variance is supplied by
#' the outcome family:
#' \describe{
#'   \item{continuous}{\eqn{v_w = \sigma^2}, the residual variance (constant).}
#'   \item{binary}{\eqn{v_w = p_w (1 - p_w)}.}
#'   \item{count}{\eqn{v_w = \mu_w} (Poisson).}
#' }
#'
#' @section Balanced-arm reading:
#' At 1:1 allocation \eqn{V_{\mathrm{eff}}(g) = 2\{V_0(g) + V_1(g)\}}, and
#' \code{bracket} \eqn{= V_{\mathrm{eff}}(g)/4} admits the decomposition
#' \deqn{\sigma^2 + V_g[\mu_0] + C_g[\mu_0, \tau] + \tfrac12 V_g[\tau],}
#' with \eqn{\tau = \mu_1 - \mu_0} the conditional average treatment effect.
#' The component columns are returned for every region and every outcome type;
#' the decomposition above is exact for continuous outcomes at 1:1 allocation.
#' Away from 1:1, \code{bracket} remains \eqn{V_{\mathrm{eff}}/4} but no longer
#' equals that sum, so prefer \code{V_eff} for any downstream calculation.
#'
#' @section Idealisations:
#' \code{V_eff} treats the region size and the arm split as fixed at their
#' expectations.  A trial draws \eqn{n_g \sim \mathrm{Bin}(n, P(g))} and
#' \eqn{n_1 \mid n_g \sim \mathrm{Bin}(n_g, p)}, and by Jensen's inequality
#' \eqn{E[1/n_1 + 1/n_0]} and \eqn{E[1/n_g]} exceed their fixed-count
#' counterparts.  Use \code{\link{fs_scale_se}} with \code{jensen = TRUE} for
#' the unconditional standard deviation.
#'
#' @section Scope:
#' Exact for identity-scale effect measures (\code{"MD"}, \code{"RD"},
#' \code{"IRD"}), where the estimator is a difference in means.  Ratio measures
#' (\code{"OR"}, \code{"RR"}, \code{"IRR"}) require a delta-method layer and are
#' rejected rather than silently approximated.
#'
#' @param dgm An object of class \code{"glm_dgm"} from
#'   \code{\link{generate_glm_dgm}}.
#' @param regions Optional named list of region specifications.  Each element is
#'   either a logical vector of length \code{nrow(dgm$df_super)} or an integer
#'   vector of row indices.  Defaults to \code{NULL}, which uses the true harm
#'   region, its complement, and the whole super-population.
#' @param harm_col Character.  Name of the true-region indicator column in
#'   \code{df_super}.  Default \code{"flag_harm"}.
#' @param rand_ratio Numeric.  Treatment:control randomisation ratio, so the
#'   randomisation probability is \code{rand_ratio / (1 + rand_ratio)}.  Default
#'   \code{1} (1:1), matching \code{\link{simulate_from_glm_dgm}}.
#' @param labels Optional character vector of length 3 renaming the default
#'   regions.  Ignored when \code{regions} is supplied.  Default
#'   \code{c("Q", "Qc", "S")}.
#'
#' @return An object of class \code{c("fs_dgm_scale", "list")}:
#'   \describe{
#'     \item{\code{regions}}{Data frame, one row per region, with columns
#'       \code{region}, \code{n_g}, \code{P_g}, \code{m_mu0}, \code{m_mu1},
#'       \code{m_tau}, \code{V_mu0}, \code{V_mu1}, \code{V_tau},
#'       \code{C_mu0_tau}, \code{v_cond0}, \code{v_cond1}, \code{V_arm0},
#'       \code{V_arm1}, \code{bracket}, \code{V_eff}.}
#'     \item{\code{sigma}}{Residual standard deviation for continuous outcomes,
#'       otherwise \code{NA_real_}.}
#'     \item{\code{outcome_type}, \code{effect_measure}}{Copied from \code{dgm}.}
#'     \item{\code{rand_ratio}, \code{p_treat}}{Allocation used.}
#'     \item{\code{n_super}}{Super-population size.}
#'   }
#'
#' @seealso \code{\link{fs_scale_se}}, \code{\link{generate_glm_dgm}}
#'
#' @examples
#' \dontrun{
#' dgm <- generate_glm_dgm(..., outcome_type = "continuous",
#'                         effect_measure = "MD")
#' sc <- fs_dgm_scale(dgm)
#' sc$regions[, c("region", "P_g", "bracket", "V_eff")]
#'
#' # Standard deviation of the estimator on the true region at n = 500
#' fs_scale_se(sc, n = 500, region = "Q")
#'
#' # Arbitrary regions, e.g. a candidate family
#' fs_dgm_scale(dgm, regions = list(big = dgm$df_super$age > 40))
#' }
#'
#' @export
fs_dgm_scale <- function(dgm,
                         regions = NULL,
                         harm_col = "flag_harm",
                         rand_ratio = 1,
                         labels = c("Q", "Qc", "S")) {

  if (!inherits(dgm, "glm_dgm")) {
    stop("'dgm' must be an object of class 'glm_dgm'.", call. = FALSE)
  }
  if (!is.numeric(rand_ratio) || length(rand_ratio) != 1L || rand_ratio <= 0) {
    stop("'rand_ratio' must be a single positive number.", call. = FALSE)
  }

  outcome_type   <- dgm$outcome_type
  effect_measure <- dgm$effect_measure

  if (!is_identity_scale(effect_measure)) {
    stop(
      "fs_dgm_scale() is exact only for identity-scale effect measures ",
      "(MD, RD, IRD); got '", effect_measure %||% "NULL", "'. Ratio measures ",
      "need a delta-method layer and are not approximated here.",
      call. = FALSE
    )
  }

  df   <- dgm$df_super
  po   <- .fs_scale_po(df, outcome_type)
  mu0  <- po$mu0
  mu1  <- po$mu1
  tau  <- mu1 - mu0
  n_super <- nrow(df)

  # Conditional outcome variance v_w(x) = Var(Y | X = x, W = w), by family.
  vc <- .fs_scale_condvar(dgm, mu0, mu1, outcome_type)

  regions <- .fs_scale_regions(regions, df, harm_col, n_super, labels)

  p_treat <- rand_ratio / (1 + rand_ratio)

  rows <- lapply(names(regions), function(nm) {
    g <- regions[[nm]]
    if (!any(g)) {
      stop("region '", nm, "' is empty.", call. = FALSE)
    }
    V_arm0 <- mean(vc$v0[g]) + .fs_var_fp(mu0[g])
    V_arm1 <- mean(vc$v1[g]) + .fs_var_fp(mu1[g])
    V_eff  <- V_arm1 / p_treat + V_arm0 / (1 - p_treat)

    data.frame(
      region    = nm,
      n_g       = sum(g),
      P_g       = mean(g),
      m_mu0     = mean(mu0[g]),
      m_mu1     = mean(mu1[g]),
      m_tau     = mean(tau[g]),
      V_mu0     = .fs_var_fp(mu0[g]),
      V_mu1     = .fs_var_fp(mu1[g]),
      V_tau     = .fs_var_fp(tau[g]),
      C_mu0_tau = .fs_cov_fp(mu0[g], tau[g]),
      v_cond0   = mean(vc$v0[g]),
      v_cond1   = mean(vc$v1[g]),
      V_arm0    = V_arm0,
      V_arm1    = V_arm1,
      bracket   = V_eff / 4,
      V_eff     = V_eff,
      stringsAsFactors = FALSE
    )
  })

  out <- list(
    regions        = do.call(rbind, rows),
    sigma          = if (identical(outcome_type, "continuous")) {
      dgm$model_params$sigma
    } else {
      NA_real_
    },
    outcome_type   = outcome_type,
    effect_measure = effect_measure,
    rand_ratio     = rand_ratio,
    p_treat        = p_treat,
    n_super        = n_super
  )
  rownames(out$regions) <- NULL
  class(out) <- c("fs_dgm_scale", "list")
  out
}


#' Estimator Standard Deviation at a Given Trial Size
#'
#' Converts a \code{\link{fs_dgm_scale}} object into the sampling standard
#' deviation of the within-region difference-in-means estimator at trial size
#' \code{n}.
#'
#' @param scale An object of class \code{"fs_dgm_scale"}.
#' @param n Integer.  Trial size.
#' @param region Character.  Region name, matching a row of
#'   \code{scale$regions}.  Default \code{NULL} returns every region.
#' @param jensen Logical.  If \code{TRUE}, inflate by the exact factor
#'   accounting for the random region size \eqn{n_g \sim \mathrm{Bin}(n, P(g))}
#'   and random arm split \eqn{n_1 \mid n_g \sim \mathrm{Bin}(n_g, p)},
#'   returning the unconditional standard deviation.  Default \code{FALSE},
#'   which returns the idealised (fixed-count) value.
#'
#' @return Named numeric vector of standard deviations.
#'
#' @seealso \code{\link{fs_dgm_scale}}
#' @export
fs_scale_se <- function(scale, n, region = NULL, jensen = FALSE) {
  if (!inherits(scale, "fs_dgm_scale")) {
    stop("'scale' must be an object of class 'fs_dgm_scale'.", call. = FALSE)
  }
  reg <- scale$regions
  if (!is.null(region)) {
    if (!all(region %in% reg$region)) {
      stop("unknown region(s): ",
           paste(setdiff(region, reg$region), collapse = ", "), call. = FALSE)
    }
    reg <- reg[reg$region %in% region, , drop = FALSE]
  }
  se <- sqrt(reg$V_eff / (n * reg$P_g))
  if (isTRUE(jensen)) {
    jf <- vapply(seq_len(nrow(reg)), function(i) {
      .fs_jensen_factor(reg$P_g[i], n = n, p_treat = scale$p_treat,
                        V0 = reg$V_arm0[i], V1 = reg$V_arm1[i])
    }, numeric(1))
    se <- se * sqrt(jf)
  }
  stats::setNames(se, reg$region)
}


#' Print method for \code{fs_dgm_scale} objects
#'
#' @param x An object of class \code{"fs_dgm_scale"} from
#'   \code{\link{fs_dgm_scale}}.
#' @param ... Unused; present for S3 compatibility.
#'
#' @return The input \code{x}, invisibly.
#'
#' @export
print.fs_dgm_scale <- function(x, ...) {
  cat("Difference-in-means sampling scale\n")
  cat(sprintf("  outcome_type   : %s (%s)\n", x$outcome_type, x$effect_measure))
  cat(sprintf("  allocation     : %g:1  (p = %.4f)\n", x$rand_ratio, x$p_treat))
  if (!is.na(x$sigma)) {
    cat(sprintf("  sigma          : %.6f   (sigma^2 = %.4f)\n",
                x$sigma, x$sigma^2))
  }
  cat(sprintf("  n_super        : %d\n\n", x$n_super))
  print(x$regions[, c("region", "n_g", "P_g", "V_mu0", "C_mu0_tau",
                      "V_tau", "bracket", "V_eff")], digits = 7)
  invisible(x)
}


# -- internals ----------------------------------------------------------------

#' Finite-population variance (divisor N)
#' @keywords internal
#' @noRd
.fs_var_fp <- function(a) mean((a - mean(a))^2)

#' Finite-population covariance (divisor N)
#' @keywords internal
#' @noRd
.fs_cov_fp <- function(a, b) mean((a - mean(a)) * (b - mean(b)))

#' Locate potential-outcome columns for the outcome family
#' @keywords internal
#' @noRd
.fs_scale_po <- function(df, outcome_type) {
  cols <- if (identical(outcome_type, "binary")) c("p0", "p1") else c("mu0", "mu1")
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols)) {
    stop("df_super is missing potential-outcome column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  list(mu0 = df[[cols[1]]], mu1 = df[[cols[2]]])
}

#' Conditional outcome variance by family
#' @keywords internal
#' @noRd
.fs_scale_condvar <- function(dgm, mu0, mu1, outcome_type) {
  switch(
    outcome_type,
    continuous = {
      s <- dgm$model_params$sigma
      if (is.null(s)) {
        stop("continuous DGM has no model_params$sigma.", call. = FALSE)
      }
      list(v0 = rep(s^2, length(mu0)), v1 = rep(s^2, length(mu1)))
    },
    binary = list(v0 = mu0 * (1 - mu0), v1 = mu1 * (1 - mu1)),
    count  = list(v0 = mu0, v1 = mu1),
    stop("unsupported outcome_type: '", outcome_type, "'.", call. = FALSE)
  )
}

#' Resolve and validate the region specification
#' @keywords internal
#' @noRd
.fs_scale_regions <- function(regions, df, harm_col, n_super, labels) {
  if (is.null(regions)) {
    if (!harm_col %in% names(df)) {
      stop("df_super has no column '", harm_col, "'.", call. = FALSE)
    }
    if (length(labels) != 3L) {
      stop("'labels' must have length 3.", call. = FALSE)
    }
    in_q <- df[[harm_col]] == 1
    out <- list(in_q, !in_q, rep(TRUE, n_super))
    names(out) <- labels
    return(out)
  }

  if (!is.list(regions) || is.null(names(regions)) || any(names(regions) == "")) {
    stop("'regions' must be a fully named list.", call. = FALSE)
  }
  lapply(regions, function(g) {
    if (is.logical(g)) {
      if (length(g) != n_super) {
        stop("logical region vectors must have length ", n_super, ".",
             call. = FALSE)
      }
      if (anyNA(g)) stop("region vectors must not contain NA.", call. = FALSE)
      return(g)
    }
    if (is.numeric(g)) {
      if (anyNA(g) || any(g < 1) || any(g > n_super)) {
        stop("region indices must lie in 1:", n_super, ".", call. = FALSE)
      }
      idx <- logical(n_super)
      idx[as.integer(g)] <- TRUE
      return(idx)
    }
    stop("each region must be a logical vector or integer indices.",
         call. = FALSE)
  })
}

#' Exact Jensen inflation factor for random region size and arm split
#'
#' Returns `E[V_1/n_1 + V_0/n_0]` divided by its fixed-count counterpart
#' `(V_1/p + V_0/(1-p)) / (n P_g)`, over `n_g ~ Bin(n, P_g)` and
#' `n_1 | n_g ~ Bin(n_g, p)`, conditioning on both arms being non-empty.  The
#' factor depends on `V_0` and `V_1` separately, not only on their combination,
#' so both are required; it reduces to a function of the counts alone when
#' `V_0 = V_1`, which holds exactly on a region where tau is constant, since
#' then `V_g[mu_1] = V_g[mu_0]`.
#' @keywords internal
#' @noRd
.fs_jensen_factor <- function(P_g, n, p_treat, V0, V1) {
  k  <- seq_len(n)
  wg <- stats::dbinom(k, n, P_g)
  s  <- sum(wg)
  if (s <= 0) return(NA_real_)
  wg <- wg / s
  contrib <- vapply(k, function(m) {
    if (m < 2L) return(0)
    j  <- seq_len(m - 1L)
    wj <- stats::dbinom(j, m, p_treat)
    sj <- sum(wj)
    if (sj <= 0) return(0)
    sum((wj / sj) * (V1 / j + V0 / (m - j)))
  }, numeric(1))
  denom <- (V1 / p_treat + V0 / (1 - p_treat)) / (n * P_g)
  sum(wg * contrib) / denom
}
