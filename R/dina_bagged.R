# File: R/dina_bagged.R
# Part of the forestsearch package.
#
# Bagged DINA estimator with infinitesimal-jackknife (IJ) variance, per
#
#   Wager, S., Hastie, T. and Efron, B. (2014).  Confidence intervals for
#   random forests: the jackknife and the infinitesimal jackknife.
#   Journal of Machine Learning Research 15: 1625-1651.
#
# applied to the DINA two-step estimator of
#
#   Gao, Z. and Hastie, T. (2025).  Estimating heterogeneous treatment
#   effects for general responses.  Biometrics 81(4): ujaf162.
#
# Builds on `dina_fit()` defined in R/dina.R: each bootstrap bag calls
# `dina_fit()` with cross-fitting (optionally) and we aggregate the
# per-bag coefficients into a bagged estimator with IJ variance.


# ---------------------------------------------------------------------------
# Main user-facing function
# ---------------------------------------------------------------------------

#' Bagged DINA with infinitesimal-jackknife variance
#'
#' `dina_fit_bagged()` fits the DINA estimator on `n_bags` bootstrap
#' replicates of the data and aggregates the resulting per-bag coefficient
#' vectors into a bagged estimator.  The variance is estimated by the
#' infinitesimal jackknife (IJ) of Wager, Hastie and Efron (2014), with
#' the finite-sample-corrected form
#' \deqn{\widehat{V}_{IJ-U} = \widehat{V}_{IJ} - (n / B) \widehat{V}_{boot},}
#' where \eqn{\widehat{V}_{boot}} is the empirical covariance of the
#' per-bag coefficients.  Negative eigenvalues introduced by the
#' Monte-Carlo correction (an artifact of small `n_bags`) are optionally
#' clipped to zero via PSD projection so that downstream `confint()` is
#' always well-defined.
#'
#' The implementation calls [dina_fit()] once per bag, with cross-fitting
#' optionally enabled within each bag.  All other arguments
#' (`propensity_method`, `baseline_method`, `cens_type`, `cens_params`,
#' `n_grid`, `eps`) are passed through unchanged.
#'
#' @inheritParams dina_fit
#' @param n_bags positive integer giving the number of bootstrap bags
#'   (default `100L`).  For stable IJ variance estimates the literature
#'   recommends `n_bags` on the order of the sample size `n`; smaller
#'   values are fine for diagnostic use but produce noisier variance.
#' @param cross_fitting_per_bag logical; if `TRUE` (the default) each bag
#'   uses `n_folds_per_bag`-fold cross-fitting internally, matching the
#'   paper.  If `FALSE`, each bag uses a single in-bag fit (faster but
#'   loses the Chernozhukov debiasing step).
#' @param n_folds_per_bag positive integer giving the number of
#'   cross-fitting folds within each bag (default `2L`).  Ignored when
#'   `cross_fitting_per_bag = FALSE`.
#' @param ij_finite_sample_correction logical; if `TRUE` (the default)
#'   the IJ variance includes the Monte-Carlo correction
#'   \eqn{(n / B) \widehat{V}_{boot}}.  If `FALSE` the raw IJ formula
#'   is returned, which tends to overstate variance for small `n_bags`.
#' @param project_psd logical; if `TRUE` (the default) the corrected IJ
#'   variance matrix is projected to its nearest positive
#'   semi-definite matrix by clipping negative eigenvalues to zero.
#'   Strongly recommended unless you want to inspect raw variance
#'   values yourself.
#' @param verbose logical; if `TRUE`, print a one-line progress message
#'   every 10 bags.  Default `FALSE`.  Per-bag progress messages are
#'   only emitted in `parallel = "none"` mode; in `parallel = "bags"`
#'   mode a single start-of-loop summary is printed if `verbose = TRUE`.
#' @param parallel character, one of `"none"` (default) or `"bags"`.
#'   When `"bags"`, the per-bag loop is dispatched via
#'   `foreach::foreach() %dofuture% { ... }`, using whatever
#'   `future::plan()` is currently active.  Call e.g.
#'   `future::plan(future::multisession, workers = N)` *before* this
#'   function to enable actual parallel execution; otherwise the
#'   parallel branch will execute serially with a warning.
#'
#' @importFrom foreach foreach
#' @importFrom future nbrOfWorkers
#'
#' @return An object of class `c("dina_bagged", "dina")`, a list
#'   containing all the components of a `"dina"` fit (so that the standard
#'   S3 methods work via inheritance) plus:
#'   \describe{
#'     \item{vcov}{the corrected, PSD-projected IJ variance matrix
#'       (or the raw IJ if `project_psd = FALSE`).}
#'     \item{vcov_ij_raw}{the IJ variance without the finite-sample
#'       correction.}
#'     \item{vcov_ij_corrected}{the IJ variance with the finite-sample
#'       correction, before PSD projection.}
#'     \item{vcov_boot}{the empirical bootstrap covariance of the per-bag
#'       coefficients, \eqn{\widehat{V}_{boot}}.}
#'     \item{bag_coefficients}{`n_bags_used` by `(d + 1)` matrix of
#'       per-bag coefficient vectors.}
#'     \item{bag_inclusion}{`n` by `n_bags_used` integer matrix of bag
#'       inclusion counts, \eqn{N_{bi}^*}.}
#'     \item{n_bags, n_bags_used}{requested and successfully-fit bag
#'       counts.}
#'     \item{n_folds_per_bag, cross_fitting_per_bag}{bag-level
#'       cross-fitting settings.}
#'   }
#'
#' @references
#' Wager, S., Hastie, T. and Efron, B. (2014).  Confidence intervals for
#' random forests: the jackknife and the infinitesimal jackknife.
#' *Journal of Machine Learning Research* 15: 1625--1651.
#'
#' Gao, Z. and Hastie, T. (2025).  Estimating heterogeneous treatment
#' effects for general responses.  *Biometrics* 81(4): ujaf162.
#' \doi{10.1093/biomtc/ujaf162}.
#'
#' @seealso [dina_fit()] for the single-pass DINA estimator with
#'   sandwich variance.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- 3
#' X <- matrix(stats::runif(n * d, -1, 1), n, d)
#' colnames(X) <- paste0("x", seq_len(d))
#' W <- stats::rbinom(n, 1, stats::plogis(0.5 * X[, 1]))
#' tau_x <- 0.5 + 0.8 * X[, 1] - 0.3 * X[, 2]
#' Y <- 1 + X[, 1] + W * tau_x + stats::rnorm(n)
#' fit <- dina_fit_bagged(X, Y, W, family = "gaussian",
#'                        n_bags = 30L, seed = 1L)
#' coef(fit)
#' sqrt(diag(vcov(fit)))
#'
#' @export
dina_fit_bagged <- function(X, Y, W,
                            family = c("gaussian", "binomial", "poisson", "cox"),
                            propensity_method = "logistic",
                            baseline_method = NULL,
                            n_bags = 100L,
                            cross_fitting_per_bag = TRUE,
                            n_folds_per_bag = 2L,
                            cens_type = NULL,
                            cens_params = list(),
                            n_grid = 500L,
                            ij_finite_sample_correction = TRUE,
                            project_psd = TRUE,
                            eps = 1e-6,
                            seed = NULL,
                            parallel = c("none", "bags"),
                            verbose = FALSE) {

  family <- match.arg(family)
  parallel <- match.arg(parallel)
  call <- match.call()

  # ---- Light input validation; dina_fit() will do the rest -------------
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("`X` must be a matrix or data frame.")
  }
  colnames_X <- colnames(X)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  d <- ncol(X)
  if (d < 1L) stop("`X` must have at least one column.")

  n_bags <- as.integer(n_bags)
  if (length(n_bags) != 1L || is.na(n_bags) || n_bags < 2L) {
    stop("`n_bags` must be a single integer >= 2.")
  }
  if (n_bags < 30L) {
    warning("`n_bags` < 30; IJ variance estimates will be noisy. ",
            "Consider increasing `n_bags`.", call. = FALSE)
  }

  n_folds_per_bag <- as.integer(n_folds_per_bag)
  if (length(n_folds_per_bag) != 1L ||
      is.na(n_folds_per_bag) || n_folds_per_bag < 1L) {
    stop("`n_folds_per_bag` must be a single positive integer.")
  }

  if (parallel == "bags") {
    n_workers <- future::nbrOfWorkers()
    if (n_workers <= 1L) {
      warning(
        "`parallel = \"bags\"` requested but `future::plan()` reports ",
        "only 1 worker; bags will run serially. ",
        "Call e.g. `future::plan(future::multisession, workers = N)` ",
        "before this function to enable actual parallelism.",
        call. = FALSE
      )
    } else if (verbose) {
      message(sprintf(
        "dina_fit_bagged: running %d bags across %d workers",
        n_bags, n_workers
      ))
    }
  }

  # ---- Bootstrap sampling ----------------------------------------------
  boot <- .dina_make_bootstrap_indices(n, n_bags, seed = seed)
  bag_inclusion <- boot$inclusion

  # ---- Per-bag fits ----------------------------------------------------
  # Each bag returns a list(coef, failed, message) so that parallel and
  # sequential code paths share the same aggregation step downstream.
  if (parallel == "bags") {
    b <- NULL  # silence R CMD check "no visible binding" for iterator
    bag_results <- foreach::foreach(
      b = seq_len(n_bags),
      .options.future = list(
        packages = "forestsearch",
        seed = TRUE
      )
    ) %dofuture% {
      .dina_safe_one_bag(
        boot_idx_b = boot$indices[, b],
        X = X, Y = Y, W = W,
        family = family, d = d,
        propensity_method = propensity_method,
        baseline_method = baseline_method,
        cross_fitting = cross_fitting_per_bag,
        n_folds = n_folds_per_bag,
        cens_type = cens_type, cens_params = cens_params,
        n_grid = n_grid, eps = eps
      )
    }
  } else {
    bag_results <- lapply(seq_len(n_bags), function(b) {
      if (verbose && b %% 10L == 0L) {
        message(sprintf("dina_fit_bagged: bag %d of %d", b, n_bags))
      }
      .dina_safe_one_bag(
        boot_idx_b = boot$indices[, b],
        X = X, Y = Y, W = W,
        family = family, d = d,
        propensity_method = propensity_method,
        baseline_method = baseline_method,
        cross_fitting = cross_fitting_per_bag,
        n_folds = n_folds_per_bag,
        cens_type = cens_type, cens_params = cens_params,
        n_grid = n_grid, eps = eps
      )
    })
  }

  # Aggregate: build matrices/vectors from the list of per-bag results
  bag_coefs    <- do.call(rbind, lapply(bag_results, `[[`, "coef"))
  bag_failed   <- vapply(bag_results, `[[`, logical(1L),   "failed")
  bag_messages <- vapply(bag_results, `[[`, character(1L), "message")

  n_bags_failed <- sum(bag_failed)
  if (n_bags_failed > 0L) {
    failure_summary <- if (any(nzchar(bag_messages))) {
      paste0(": ", unique(bag_messages[nzchar(bag_messages)])[1L])
    } else {
      ""
    }
    warning(sprintf("%d of %d bags failed and were excluded%s",
                    n_bags_failed, n_bags, failure_summary),
            call. = FALSE)
  }

  ok <- !bag_failed
  bag_coefs <- bag_coefs[ok, , drop = FALSE]
  bag_inclusion <- bag_inclusion[, ok, drop = FALSE]
  n_bags_used <- nrow(bag_coefs)

  if (n_bags_used < 2L) {
    stop("Fewer than 2 bags succeeded; cannot compute IJ variance.")
  }

  # ---- Aggregate -------------------------------------------------------
  beta_hat <- colMeans(bag_coefs)

  vcov_boot <- stats::cov(bag_coefs)

  vcov_ij_raw <- .dina_ij_variance(
    bag_coefs, bag_inclusion, finite_sample_correction = FALSE
  )
  vcov_ij_corrected <- if (ij_finite_sample_correction) {
    .dina_ij_variance(
      bag_coefs, bag_inclusion, finite_sample_correction = TRUE
    )
  } else {
    vcov_ij_raw
  }

  # PSD projection: clip negative eigenvalues to zero
  if (project_psd) {
    vcov_out <- .dina_psd_project(vcov_ij_corrected, tol = 0)
    neg_eig <- attr(vcov_out, "n_negative_eigenvalues")
    if (!is.null(neg_eig) && neg_eig > 0L) {
      warning(sprintf(
        "IJ variance had %d negative eigenvalue(s); clipped to zero. ",
        neg_eig
      ), "Consider increasing `n_bags`.", call. = FALSE)
    }
  } else {
    vcov_out <- vcov_ij_corrected
  }

  # ---- Naming ----------------------------------------------------------
  out_names <- c(
    "(Intercept)",
    if (!is.null(colnames_X)) colnames_X else paste0("X", seq_len(d))
  )
  names(beta_hat) <- out_names
  dimnames(vcov_out) <- list(out_names, out_names)
  dimnames(vcov_ij_raw) <- list(out_names, out_names)
  dimnames(vcov_ij_corrected) <- list(out_names, out_names)
  dimnames(vcov_boot) <- list(out_names, out_names)
  colnames(bag_coefs) <- out_names

  # ---- Assemble result -------------------------------------------------
  out <- list(
    coefficients          = beta_hat,
    vcov                  = vcov_out,
    vcov_ij_raw           = vcov_ij_raw,
    vcov_ij_corrected     = vcov_ij_corrected,
    vcov_boot             = vcov_boot,
    bag_coefficients      = bag_coefs,
    bag_inclusion         = bag_inclusion,
    n_bags                = n_bags,
    n_bags_used           = n_bags_used,
    n_folds_per_bag       = n_folds_per_bag,
    cross_fitting_per_bag = cross_fitting_per_bag,
    family                = family,
    n                     = n,
    d                     = d,
    cross_fitting         = cross_fitting_per_bag,
    n_folds               = n_folds_per_bag,
    call                  = call
  )
  class(out) <- c("dina_bagged", "dina")
  out
}


# ---------------------------------------------------------------------------
# Internal: fit DINA on one bag
# ---------------------------------------------------------------------------

#' Fit DINA on a single bootstrap bag and return only the coefficient vector.
#'
#' Thin wrapper around `dina_fit()` that strips down to just the coefficient
#' vector and avoids paying for `summary` / `confint` materialization on
#' every bag.  Errors propagate to the caller for the `tryCatch` upstream.
#'
#' @noRd
.dina_fit_one_bag <- function(X_b, Y_b, W_b, family,
                              propensity_method, baseline_method,
                              cross_fitting, n_folds,
                              cens_type, cens_params, n_grid, eps) {
  fit_b <- dina_fit(
    X = X_b, Y = Y_b, W = W_b,
    family = family,
    propensity_method = propensity_method,
    baseline_method = baseline_method,
    cross_fitting = cross_fitting,
    n_folds = n_folds,
    cens_type = cens_type, cens_params = cens_params,
    n_grid = n_grid, eps = eps,
    seed = NULL
  )
  unname(fit_b$coefficients)
}

#' Subset a DINA response object by row indices.
#' @noRd
.dina_subset_Y <- function(Y, idx) {
  if (is.data.frame(Y)) Y[idx, , drop = FALSE] else Y[idx]
}

#' Fit DINA on a single bag with error handling, returning a list result.
#'
#' Wraps `.dina_fit_one_bag()` in `tryCatch()` so that bag-level failures
#' surface as `list(coef = NA, failed = TRUE, message = ...)` instead of
#' as R conditions.  This uniform list contract is what the parallel
#' (`foreach %dofuture%`) and sequential (`lapply`) code paths in
#' `dina_fit_bagged()` both consume.
#'
#' @noRd
.dina_safe_one_bag <- function(boot_idx_b, X, Y, W, family, d,
                               propensity_method, baseline_method,
                               cross_fitting, n_folds,
                               cens_type, cens_params, n_grid, eps) {
  tryCatch(
    {
      cf <- .dina_fit_one_bag(
        X_b = X[boot_idx_b, , drop = FALSE],
        Y_b = .dina_subset_Y(Y, boot_idx_b),
        W_b = W[boot_idx_b],
        family = family,
        propensity_method = propensity_method,
        baseline_method = baseline_method,
        cross_fitting = cross_fitting,
        n_folds = n_folds,
        cens_type = cens_type, cens_params = cens_params,
        n_grid = n_grid, eps = eps
      )
      list(coef = unname(cf), failed = FALSE, message = "")
    },
    error = function(e) {
      list(
        coef    = rep(NA_real_, d + 1L),
        failed  = TRUE,
        message = conditionMessage(e)
      )
    }
  )
}


# ---------------------------------------------------------------------------
# Internal: bootstrap indices and bag-inclusion counts
# ---------------------------------------------------------------------------

#' Build `n_bags` bootstrap samples and the associated inclusion-count matrix.
#'
#' Returns:
#'   * `indices`   --  an `n` x `n_bags` integer matrix whose b-th column
#'                     is a sample-with-replacement of `seq_len(n)`.
#'   * `inclusion` --  an `n` x `n_bags` integer matrix whose `(i, b)`
#'                     entry is the number of times row `i` of the original
#'                     data appears in bag `b` (i.e., \eqn{N_{bi}^*}).
#'
#' @noRd
.dina_make_bootstrap_indices <- function(n, n_bags, seed = NULL) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      on.exit(
        assign(".Random.seed", old_seed, envir = .GlobalEnv),
        add = TRUE
      )
    } else {
      on.exit(
        rm(list = ".Random.seed", envir = .GlobalEnv, inherits = FALSE),
        add = TRUE
      )
    }
    set.seed(seed)
  }

  indices <- matrix(
    sample.int(n, n * n_bags, replace = TRUE),
    nrow = n, ncol = n_bags
  )

  inclusion <- matrix(0L, nrow = n, ncol = n_bags)
  for (b in seq_len(n_bags)) {
    inclusion[, b] <- tabulate(indices[, b], nbins = n)
  }

  list(indices = indices, inclusion = inclusion)
}


# ---------------------------------------------------------------------------
# Internal: infinitesimal jackknife variance
# ---------------------------------------------------------------------------

#' Compute the infinitesimal-jackknife variance of a bagged estimator.
#'
#' Implements the IJ-for-bagging variance estimator of Wager, Hastie and
#' Efron (2014), with optional Monte-Carlo correction for finite `B`.
#' Specifically, with
#' \deqn{C_i = \frac{1}{B} \sum_{b=1}^B (N_{bi}^* - \bar N_i)
#'            (\hat\theta_b - \bar\theta),}
#' the raw IJ variance matrix is
#' \deqn{\widehat V_{IJ} = \sum_{i=1}^n C_i C_i^\top}
#' and the corrected form is
#' \deqn{\widehat V_{IJ-U} = \widehat V_{IJ} - (n / B) \widehat V_{boot},}
#' where \eqn{\widehat V_{boot}} is the empirical covariance of the
#' per-bag coefficient vectors (with `B - 1` denominator, matching
#' `stats::cov`).
#'
#' @param bag_coefs `B` x `p` matrix of per-bag coefficient vectors.
#' @param bag_inclusion `n` x `B` matrix of bag-inclusion counts
#'   \eqn{N_{bi}^*}.
#' @param finite_sample_correction logical; if `TRUE` subtract
#'   \eqn{(n / B) \widehat V_{boot}}.
#'
#' @return a `p` x `p` symmetric matrix.
#'
#' @noRd
.dina_ij_variance <- function(bag_coefs, bag_inclusion,
                              finite_sample_correction = TRUE) {
  B <- nrow(bag_coefs)
  p <- ncol(bag_coefs)
  n <- nrow(bag_inclusion)

  beta_bar <- colMeans(bag_coefs)
  N_bar    <- rowMeans(bag_inclusion)

  bag_coefs_c <- sweep(bag_coefs, 2L, beta_bar, FUN = "-")  # B x p
  N_c <- sweep(bag_inclusion, 1L, N_bar, FUN = "-")         # n x B

  # C is n x p:  C[i, j] = (1/B) sum_b N_c[i, b] * bag_coefs_c[b, j]
  C <- (N_c %*% bag_coefs_c) / B

  V_ij <- crossprod(C)  # p x p

  if (finite_sample_correction) {
    V_boot <- stats::cov(bag_coefs)
    V_ij <- V_ij - (n / B) * V_boot
  }

  # Symmetrize to suppress numerical asymmetry from floating-point arithmetic
  (V_ij + t(V_ij)) / 2
}


# ---------------------------------------------------------------------------
# Internal: PSD projection
# ---------------------------------------------------------------------------

#' Project a symmetric matrix onto the cone of positive semi-definite
#' matrices in the Frobenius norm by clipping negative eigenvalues to
#' `tol` (default `0`).
#'
#' Returns the projected matrix with an attribute
#' `n_negative_eigenvalues` recording how many eigenvalues were clipped.
#'
#' @noRd
.dina_psd_project <- function(M, tol = 0) {
  M_sym <- (M + t(M)) / 2
  eig <- eigen(M_sym, symmetric = TRUE)
  neg <- eig$values < tol
  n_neg <- sum(neg)
  eig$values[neg] <- tol
  out <- eig$vectors %*% (eig$values * t(eig$vectors))
  out <- (out + t(out)) / 2
  attr(out, "n_negative_eigenvalues") <- n_neg
  out
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print a bagged DINA fit.
#'
#' Extends `print.dina()` with bag count and per-bag cross-fitting info.
#'
#' @param x an object of class `"dina_bagged"`.
#' @param digits number of digits to show for coefficients.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @export
print.dina_bagged <- function(x,
                              digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("Bagged DINA fit (family = \"", x$family, "\", n = ", x$n,
      ")\n", sep = "")
  cat("Bags: ", x$n_bags_used, "/", x$n_bags, " successful",
      sep = "")
  if (x$cross_fitting_per_bag) {
    cat(", ", x$n_folds_per_bag, "-fold cross-fitting per bag\n",
        sep = "")
  } else {
    cat(", no cross-fitting per bag\n")
  }
  cat("Variance: infinitesimal jackknife (IJ)",
      if (isTRUE(attr(x$vcov, "n_negative_eigenvalues") > 0L))
        ", PSD-projected" else "",
      "\n", sep = "")
  cat("\nCoefficients:\n")
  print.default(
    format(x$coefficients, digits = digits),
    print.gap = 2L, quote = FALSE
  )
  invisible(x)
}
