# File: R/dina.R
# Part of the forestsearch package.
#
# Clean-room reimplementation of the DINA (DIfference in NAtural parameters)
# estimator following:
#
#   Gao, Z. and Hastie, T. (2025).  Estimating heterogeneous treatment effects
#   for general responses.  Biometrics 81(4): ujaf162.
#   doi:10.1093/biomtc/ujaf162.
#
# Implements Algorithm 1 (exponential family; Section 3.2) and the
# partial-likelihood Cox extension (Section 4) of the paper.  Variance
# estimation uses the standard sandwich estimator for the GLM step 2 and the
# Cox model's information-based variance for the Cox step 2.
#
# This is an independent implementation written from the mathematics in the
# Biometrics article.  No code has been taken from the published
# supplementary materials, which are distributed under CC BY-NC 4.0 and
# therefore incompatible with CRAN inclusion.  Methodology is attributed via
# the `@references` field below.


# ---------------------------------------------------------------------------
# Main user-facing function
# ---------------------------------------------------------------------------

#' Estimate heterogeneous treatment effects via DINA
#'
#' `dina_fit()` implements the two-step DINA (DIfference in NAtural parameters)
#' estimator of Gao and Hastie (2025) for heterogeneous treatment effects
#' (HTE) with general response types: Gaussian (mean difference), binomial
#' (log odds ratio), Poisson (log mean ratio), and Cox (log hazard ratio).
#' The treatment effect is assumed to be linear in the covariates,
#' \deqn{\tau(x) = \beta_0 + x^\top \beta_{1:d},}
#' following Proposition 1 of the paper.  The procedure has two steps:
#' (1) estimate nuisance functions \eqn{a(x)} and \eqn{\nu(x)} on a training
#' fold, and (2) plug those estimators into a GLM (or Cox) likelihood on the
#' held-out fold to estimate \eqn{\beta}.  Cross-fitting averages the two
#' estimators obtained by swapping the folds.
#'
#' The Gaussian family uses the simplification of Section 3.2: \eqn{a(x) =
#' e(x)} and \eqn{\nu(x) = E[Y \mid X = x]}, so the conditional outcome means
#' \eqn{\eta_0(x)}, \eqn{\eta_1(x)} are not required separately.  For the
#' binomial and Poisson families, separate within-arm fits to \eqn{\eta_0(x)}
#' and \eqn{\eta_1(x)} are used to construct \eqn{a(x)} via the modified
#' propensity in Eq. (12) of the paper.  For the Cox family, \eqn{a(x)} uses
#' the not-censoring probabilities of Eq. (18); how those are estimated is
#' controlled by `cens_type`.
#'
#' Variance estimation follows the standard sandwich form:
#' \eqn{\widehat{V}(\hat\beta) = \widehat\Sigma_1^{-1} \widehat\Sigma_2
#' \widehat\Sigma_1^{-1} / n}, where \eqn{\widehat\Sigma_1} and
#' \eqn{\widehat\Sigma_2} are the bread (negative score derivative) and meat
#' (outer product of scores) matrices.  Under cross-fitting these are
#' averaged across folds before the sandwich is formed.  For the Cox step 2
#' we use the partial-likelihood information-matrix variance returned by
#' [survival::coxph()].
#'
#' @param X a numeric covariate matrix or data frame with `n` rows and `d`
#'   columns.  All columns are treated as numeric predictors.  Column names,
#'   when present, are propagated to the coefficient names of the fitted
#'   object.
#' @param Y the response.  For `family = "cox"`, a data frame (or matrix
#'   coercible to one) with columns `time` (event/censoring time) and
#'   `status` (1 = event, 0 = censored).  For all other families, a numeric
#'   vector of length `n`.
#' @param W binary treatment indicator (0 = control, 1 = treatment), length
#'   `n`.
#' @param family character; one of `"gaussian"`, `"binomial"`, `"poisson"`,
#'   `"cox"`.
#' @param propensity_method how to estimate the propensity score
#'   \eqn{e(x) = P(W = 1 \mid X = x)}.  One of:
#'   * the string `"logistic"` (the default): a logistic GLM on `X`;
#'   * a numeric scalar in `(0, 1)`: a known constant propensity (e.g.,
#'     `0.5` for a balanced randomized trial);
#'   * a function taking a covariate matrix and returning a numeric vector
#'     of fitted propensities, allowing the user to plug in random forests,
#'     boosting, etc.
#' @param baseline_method how to estimate the natural-parameter functions
#'   \eqn{\eta_0(x)}, \eqn{\eta_1(x)} (or the marginal mean for the Gaussian
#'   family).  One of:
#'   * `NULL` (the default): use `"glm"` for non-Cox families and `"cox"` for
#'     the Cox family;
#'   * `"glm"`: a GLM with the natural link for the family;
#'   * `"cox"`: [survival::coxph()] (Cox family only);
#'   * a list with components `eta0` and `eta1`, each a function taking a
#'     covariate matrix and returning a vector of linear predictors.  For
#'     Gaussian, a list with a single component `nu` may be supplied
#'     instead.
#' @param cross_fitting logical; if `TRUE` (the default) use `n_folds`-fold
#'   cross-fitting to separate nuisance estimation from target-parameter
#'   estimation.  If `FALSE`, the same data are used for both steps.
#' @param n_folds positive integer giving the number of cross-fitting folds
#'   (default `2`, matching the paper).
#' @param cens_type for `family = "cox"`, the censoring mechanism.  One of
#'   `NULL` (no censoring), `"fixed"`, `"uniform"`, `"exponential"`, or
#'   `"unknown"` (estimate the probability of *not* being censored from `X`
#'   by logistic regression).  Ignored for non-Cox families.
#' @param cens_params named list of censoring parameters.  Required entries
#'   depend on `cens_type`:
#'   * `"fixed"`, `"uniform"`: `T0` (the administrative censoring time).
#'   * `"exponential"`: `rate_c` (the exponential rate).
#'   * `"unknown"`, `NULL`: ignored.
#' @param n_grid positive integer giving the number of time-grid points used
#'   for numerical integration of the not-censoring probability for
#'   `cens_type %in% c("uniform", "exponential")` (default `500`).
#' @param eps small positive constant added to denominators for numerical
#'   stability (default `1e-6`).
#' @param seed optional integer seed for the random cross-fitting fold
#'   assignment.  When supplied, the global RNG state is restored on exit.
#'
#' @return An object of class `"dina"`, a list with components:
#' \describe{
#'   \item{coefficients}{named numeric vector of estimated DINA coefficients
#'     \eqn{\hat\beta}, with `(Intercept)` first followed by one entry per
#'     covariate.}
#'   \item{vcov}{estimated variance-covariance matrix of \eqn{\hat\beta}.}
#'   \item{family}{the response family, as supplied.}
#'   \item{n, d}{sample size and number of covariates.}
#'   \item{cross_fitting, n_folds}{cross-fitting settings.}
#'   \item{call}{the matched call.}
#' }
#'
#' @references
#' Gao, Z. and Hastie, T. (2025).  Estimating heterogeneous treatment
#' effects for general responses.  *Biometrics* 81(4): ujaf162.
#' \doi{10.1093/biomtc/ujaf162}.
#'
#' Robinson, P. M. (1988).  Root-n-consistent semiparametric regression.
#' *Econometrica* 56(4): 931--954.
#'
#' Nie, X. and Wager, S. (2021).  Quasi-oracle estimation of heterogeneous
#' treatment effects.  *Biometrika* 108(2): 299--319.
#'
#' Dandl, S., Bender, A. and Hothorn, T. (2024).  Heterogeneous treatment
#' effect estimation for observational data using model-based forests.
#' *Statistical Methods in Medical Research* 33(3): 392--413.
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
#' fit <- dina_fit(X, Y, W, family = "gaussian", seed = 1)
#' coef(fit)
#' confint(fit)
#'
#' @importFrom stats glm glm.fit binomial poisson gaussian predict residuals
#' @importFrom stats coef qnorm pnorm qlogis plogis approx printCoefmat
#' @importFrom stats setNames model.frame as.formula
#' @importFrom survival coxph Surv
#' @export
dina_fit <- function(X, Y, W,
                     family = c("gaussian", "binomial", "poisson", "cox"),
                     propensity_method = "logistic",
                     baseline_method = NULL,
                     cross_fitting = TRUE,
                     n_folds = 2L,
                     cens_type = NULL,
                     cens_params = list(),
                     n_grid = 500L,
                     eps = 1e-6,
                     seed = NULL) {

  family <- match.arg(family)
  call <- match.call()

  # ---- Input validation -------------------------------------------------
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("`X` must be a matrix or data frame.")
  }
  colnames_X <- colnames(X)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  d <- ncol(X)
  if (d < 1L) stop("`X` must have at least one column.")

  if (family == "cox") {
    if (is.matrix(Y)) Y <- as.data.frame(Y)
    if (!is.data.frame(Y) || !all(c("time", "status") %in% names(Y))) {
      stop("For `family = \"cox\"`, `Y` must be a data frame with columns ",
           "'time' and 'status'.")
    }
    if (nrow(Y) != n) stop("`nrow(Y)` must equal `nrow(X)`.")
    if (any(Y$time <= 0)) stop("All `Y$time` values must be positive.")
    if (!all(Y$status %in% c(0, 1))) {
      stop("`Y$status` must contain only 0 and 1.")
    }
  } else {
    Y <- as.numeric(Y)
    if (length(Y) != n) stop("`length(Y)` must equal `nrow(X)`.")
    if (family == "binomial" && !all(Y %in% c(0, 1))) {
      stop("For `family = \"binomial\"`, `Y` must contain only 0 and 1.")
    }
    if (family == "poisson" && (any(Y < 0) || any(Y != round(Y)))) {
      stop("For `family = \"poisson\"`, `Y` must be non-negative integers.")
    }
  }

  W <- as.numeric(W)
  if (length(W) != n) stop("`length(W)` must equal `nrow(X)`.")
  if (!all(W %in% c(0, 1))) stop("`W` must contain only 0 and 1.")

  if (!is.logical(cross_fitting) || length(cross_fitting) != 1L) {
    stop("`cross_fitting` must be a single logical value.")
  }

  n_folds <- as.integer(n_folds)
  if (length(n_folds) != 1L || is.na(n_folds) || n_folds < 1L) {
    stop("`n_folds` must be a single positive integer.")
  }

  if (is.null(baseline_method)) {
    baseline_method <- if (family == "cox") "cox" else "glm"
  }

  # Standardize internal column names; preserve user names for output
  colnames(X) <- paste0("X", seq_len(d))

  # ---- Run the procedure ------------------------------------------------
  if (!cross_fitting || n_folds == 1L) {
    res <- .dina_one_pass(
      X_train = X, Y_train = Y, W_train = W,
      X_test  = X, Y_test  = Y, W_test  = W,
      family = family,
      propensity_method = propensity_method,
      baseline_method = baseline_method,
      cens_type = cens_type, cens_params = cens_params,
      n_grid = n_grid, eps = eps
    )
    beta_hat <- res$beta
    sigma1 <- res$sigma1
    sigma2 <- res$sigma2
    n_eff <- n
  } else {
    folds <- .dina_make_folds(n, n_folds, seed)
    beta_list <- vector("list", n_folds)
    sigma1_list <- vector("list", n_folds)
    sigma2_list <- vector("list", n_folds)

    for (k in seq_len(n_folds)) {
      train_idx <- folds$train[[k]]
      test_idx  <- folds$test[[k]]
      Y_train <- if (is.data.frame(Y)) Y[train_idx, , drop = FALSE] else Y[train_idx]
      Y_test  <- if (is.data.frame(Y)) Y[test_idx,  , drop = FALSE] else Y[test_idx]
      res_k <- .dina_one_pass(
        X_train = X[train_idx, , drop = FALSE],
        Y_train = Y_train,
        W_train = W[train_idx],
        X_test  = X[test_idx, , drop = FALSE],
        Y_test  = Y_test,
        W_test  = W[test_idx],
        family = family,
        propensity_method = propensity_method,
        baseline_method = baseline_method,
        cens_type = cens_type, cens_params = cens_params,
        n_grid = n_grid, eps = eps
      )
      beta_list[[k]] <- res_k$beta
      sigma1_list[[k]] <- res_k$sigma1
      sigma2_list[[k]] <- res_k$sigma2
    }

    beta_mat <- do.call(rbind, beta_list)
    beta_hat <- colMeans(beta_mat, na.rm = TRUE)
    sigma1 <- Reduce(`+`, sigma1_list) / n_folds
    sigma2 <- Reduce(`+`, sigma2_list) / n_folds
    n_eff <- n
  }

  # ---- Sandwich variance ------------------------------------------------
  vcov_beta <- tryCatch(
    solve(sigma1, tol = .Machine$double.eps) %*% sigma2 %*%
      solve(sigma1, tol = .Machine$double.eps) / n_eff,
    error = function(e) {
      warning("Variance matrix computation failed: ",
              conditionMessage(e),
              ".  Returning matrix of NA.",
              call. = FALSE)
      matrix(NA_real_, length(beta_hat), length(beta_hat))
    }
  )

  out_names <- c(
    "(Intercept)",
    if (!is.null(colnames_X)) colnames_X else paste0("X", seq_len(d))
  )
  names(beta_hat) <- out_names
  dimnames(vcov_beta) <- list(out_names, out_names)

  out <- list(
    coefficients   = beta_hat,
    vcov           = vcov_beta,
    family         = family,
    n              = n,
    d              = d,
    cross_fitting  = cross_fitting,
    n_folds        = n_folds,
    call           = call
  )
  class(out) <- "dina"
  out
}


# ---------------------------------------------------------------------------
# Internal: one cross-fitting pass (train -> test)
# ---------------------------------------------------------------------------

#' One cross-fitting pass: fit nuisance functions on the training data, then
#' estimate the target parameter on the test data.
#'
#' @return list with `beta`, `sigma1`, `sigma2` (the bread / meat matrices).
#'   For the Cox family, `sigma1 == sigma2 == solve(coxph$var) / n_test`
#'   so that the sandwich reduces to the partial-likelihood variance.
#' @noRd
.dina_one_pass <- function(X_train, Y_train, W_train,
                           X_test, Y_test, W_test,
                           family,
                           propensity_method, baseline_method,
                           cens_type, cens_params,
                           n_grid, eps) {

  prop_fn <- .dina_fit_propensity(X_train, W_train, propensity_method)

  if (family == "cox") {
    base_fns <- .dina_fit_baselines_cox(
      X_train, Y_train, W_train, baseline_method
    )
    censor_fns <- .dina_fit_censoring_cox(
      X_train, Y_train, W_train,
      eta0_train = base_fns$eta0_train,
      eta1_train = base_fns$eta1_train,
      cens_type = cens_type, cens_params = cens_params, n_grid = n_grid
    )
    nuisance <- .dina_combine_nuisance_cox(
      prop_fn = prop_fn,
      eta0_fn = base_fns$eta0_fn, eta1_fn = base_fns$eta1_fn,
      censor0_fn = censor_fns$censor0_fn,
      censor1_fn = censor_fns$censor1_fn,
      eps = eps
    )
    .dina_step2_cox(X_test, Y_test, W_test, nuisance)
  } else {
    lv <- .dina_link_and_variance(family)
    base_fns <- .dina_fit_baselines_glm(
      X_train, Y_train, W_train, family, baseline_method
    )
    nuisance <- .dina_combine_nuisance_glm(
      family = family,
      prop_fn = prop_fn,
      eta0_fn = base_fns$eta0_fn, eta1_fn = base_fns$eta1_fn,
      nu_fn   = base_fns$nu_fn,
      link_fn = lv$link, var_fn = lv$var,
      eps = eps
    )
    .dina_step2_glm(X_test, Y_test, W_test, family, nuisance)
  }
}


# ---------------------------------------------------------------------------
# Internal: link / variance functions per family
# ---------------------------------------------------------------------------

#' Family link, inverse-link, and variance functions.
#'
#' For the exponential family with canonical link, the variance function
#' \eqn{V(\mu) = d\mu/d\eta} appears in the modified propensity
#' \eqn{a(x)} of Eq. (12) in Gao and Hastie (2025).
#'
#' @noRd
.dina_link_and_variance <- function(family) {
  switch(
    family,
    gaussian = list(
      link = function(mu) mu,
      inv  = function(eta) eta,
      var  = function(mu) rep_len(1, length(mu))
    ),
    binomial = list(
      link = function(mu) stats::qlogis(.dina_clip01(mu)),
      inv  = stats::plogis,
      var  = function(mu) {
        mu <- .dina_clip01(mu)
        mu * (1 - mu)
      }
    ),
    poisson = list(
      link = function(mu) log(pmax(mu, .Machine$double.eps)),
      inv  = exp,
      var  = function(mu) pmax(mu, .Machine$double.eps)
    ),
    cox = list(link = NULL, inv = NULL, var = NULL)
  )
}

#' Clip probabilities to (eps_p, 1 - eps_p) for numerical stability.
#' @noRd
.dina_clip01 <- function(p, eps_p = 1e-8) {
  pmin(pmax(p, eps_p), 1 - eps_p)
}


# ---------------------------------------------------------------------------
# Internal: cross-fitting folds
# ---------------------------------------------------------------------------

#' Build k cross-fitting folds.
#'
#' @return list with `train` and `test` (each a list of length `k` of
#'   integer indices into `seq_len(n)`).
#' @noRd
.dina_make_folds <- function(n, k, seed = NULL) {
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
  perm <- sample.int(n)
  fold_id <- ((seq_len(n) - 1L) %% k) + 1L
  test_idx <- split(perm, fold_id)
  train_idx <- lapply(test_idx, function(i) setdiff(seq_len(n), i))
  names(train_idx) <- names(test_idx) <- NULL
  list(train = train_idx, test = test_idx)
}


# ---------------------------------------------------------------------------
# Internal: propensity score estimation
# ---------------------------------------------------------------------------

#' Fit propensity on training data; return a function evaluating on new X.
#' @noRd
.dina_fit_propensity <- function(X_train, W_train, propensity_method) {
  if (is.function(propensity_method)) {
    return(function(X_new) as.numeric(propensity_method(X_new)))
  }
  if (is.numeric(propensity_method) && length(propensity_method) == 1L) {
    p_const <- as.numeric(propensity_method)
    if (is.na(p_const) || p_const <= 0 || p_const >= 1) {
      stop("Scalar `propensity_method` must be in (0, 1).")
    }
    return(function(X_new) rep_len(p_const, nrow(X_new)))
  }
  if (identical(propensity_method, "logistic")) {
    train_df <- as.data.frame(X_train)
    train_df$.W <- W_train
    fit <- stats::glm(
      .W ~ ., data = train_df, family = stats::binomial(link = "logit")
    )
    return(function(X_new) {
      new_df <- as.data.frame(X_new)
      names(new_df) <- setdiff(names(train_df), ".W")
      as.numeric(stats::predict(fit, newdata = new_df, type = "response"))
    })
  }
  stop(
    "`propensity_method` must be \"logistic\", a function, or a scalar in (0, 1)."
  )
}


# ---------------------------------------------------------------------------
# Internal: baseline natural-parameter estimation (GLM families)
# ---------------------------------------------------------------------------

#' Fit baseline natural-parameter functions for GLM families.
#'
#' Returns predictors on the linear-predictor scale (`eta0_fn`, `eta1_fn`)
#' and, for Gaussian, also `nu_fn` --- the marginal mean on the original
#' scale, which equals \eqn{\nu(x)} directly per Section 3.2.
#'
#' @noRd
.dina_fit_baselines_glm <- function(X_train, Y_train, W_train, family,
                                    baseline_method) {

  # User-supplied closures
  if (is.list(baseline_method)) {
    if (family == "gaussian" && !is.null(baseline_method$nu)) {
      return(list(
        eta0_fn = NULL, eta1_fn = NULL,
        nu_fn = function(X_new) as.numeric(baseline_method$nu(X_new))
      ))
    }
    if (!is.null(baseline_method$eta0) && !is.null(baseline_method$eta1)) {
      return(list(
        eta0_fn = function(X_new) as.numeric(baseline_method$eta0(X_new)),
        eta1_fn = function(X_new) as.numeric(baseline_method$eta1(X_new)),
        nu_fn   = NULL
      ))
    }
    stop("`baseline_method` list must supply either `nu` (Gaussian only) ",
         "or both `eta0` and `eta1`.")
  }

  if (!identical(baseline_method, "glm")) {
    stop("`baseline_method` must be \"glm\" or a list of functions ",
         "for non-Cox families.")
  }

  glm_family <- switch(
    family,
    gaussian = stats::gaussian(link = "identity"),
    binomial = stats::binomial(link = "logit"),
    poisson  = stats::poisson(link = "log")
  )

  if (family == "gaussian") {
    # Direct estimation of nu(x) = E[Y | X] using all training data
    train_df <- as.data.frame(X_train)
    train_df$.Y <- Y_train
    fit_nu <- stats::glm(.Y ~ ., data = train_df, family = glm_family)
    nu_fn <- function(X_new) {
      new_df <- as.data.frame(X_new)
      names(new_df) <- setdiff(names(train_df), ".Y")
      as.numeric(stats::predict(fit_nu, newdata = new_df, type = "response"))
    }
    return(list(eta0_fn = NULL, eta1_fn = NULL, nu_fn = nu_fn))
  }

  # Binomial / Poisson: separate within-arm fits
  fit_arm <- function(idx) {
    df_arm <- as.data.frame(X_train[idx, , drop = FALSE])
    df_arm$.Y <- Y_train[idx]
    stats::glm(.Y ~ ., data = df_arm, family = glm_family)
  }
  fit_0 <- fit_arm(W_train == 0)
  fit_1 <- fit_arm(W_train == 1)

  make_eta_fn <- function(fit) {
    function(X_new) {
      new_df <- as.data.frame(X_new)
      names(new_df) <- setdiff(names(stats::model.frame(fit)), ".Y")
      as.numeric(stats::predict(fit, newdata = new_df, type = "link"))
    }
  }
  list(
    eta0_fn = make_eta_fn(fit_0),
    eta1_fn = make_eta_fn(fit_1),
    nu_fn   = NULL
  )
}


# ---------------------------------------------------------------------------
# Internal: baseline natural-parameter estimation (Cox family)
# ---------------------------------------------------------------------------

#' Fit Cox baselines on training data.
#'
#' Returns prediction functions for the linear predictors and the training
#' linear predictors themselves (the latter are needed to estimate the
#' baseline cumulative hazard for the not-censoring probability).
#'
#' @noRd
.dina_fit_baselines_cox <- function(X_train, Y_train, W_train,
                                    baseline_method) {

  if (is.list(baseline_method)) {
    if (is.null(baseline_method$eta0) || is.null(baseline_method$eta1)) {
      stop("`baseline_method` list must supply both `eta0` and `eta1`.")
    }
    eta0_fn <- function(X_new) as.numeric(baseline_method$eta0(X_new))
    eta1_fn <- function(X_new) as.numeric(baseline_method$eta1(X_new))
    return(list(
      eta0_fn = eta0_fn, eta1_fn = eta1_fn,
      eta0_train = eta0_fn(X_train[W_train == 0, , drop = FALSE]),
      eta1_train = eta1_fn(X_train[W_train == 1, , drop = FALSE])
    ))
  }

  if (!identical(baseline_method, "cox")) {
    stop("`baseline_method` must be \"cox\" or a list of functions for ",
         "the Cox family.")
  }

  fit_arm <- function(idx) {
    df_arm <- as.data.frame(X_train[idx, , drop = FALSE])
    df_arm$.time <- Y_train$time[idx]
    df_arm$.status <- Y_train$status[idx]
    rhs <- paste(setdiff(names(df_arm), c(".time", ".status")), collapse = " + ")
    fmla <- stats::as.formula(paste(
      "survival::Surv(.time, .status) ~", rhs
    ))
    survival::coxph(fmla, data = df_arm, x = FALSE, y = FALSE)
  }
  fit_0 <- fit_arm(W_train == 0)
  fit_1 <- fit_arm(W_train == 1)

  make_eta_fn <- function(fit) {
    cov_names <- attr(fit$terms, "term.labels")
    function(X_new) {
      new_df <- as.data.frame(X_new)
      names(new_df) <- cov_names
      as.numeric(stats::predict(fit, newdata = new_df, type = "lp"))
    }
  }
  eta0_fn <- make_eta_fn(fit_0)
  eta1_fn <- make_eta_fn(fit_1)

  list(
    eta0_fn = eta0_fn, eta1_fn = eta1_fn,
    eta0_train = eta0_fn(X_train[W_train == 0, , drop = FALSE]),
    eta1_train = eta1_fn(X_train[W_train == 1, , drop = FALSE])
  )
}


# ---------------------------------------------------------------------------
# Internal: not-censoring probability estimation (Cox family)
# ---------------------------------------------------------------------------

#' Build per-arm not-censoring probability functions P(C >= Y | W = w, X).
#'
#' Implements the four mechanisms described after Eq. (18) of
#' Gao and Hastie (2025).  Uses the integration-by-parts identity
#' \eqn{E[S_C(Y)] = \int F_Y(t) f_C(t) dt} for "uniform" and "exponential"
#' censoring; closed-form expressions for "fixed" and "unknown".
#'
#' @noRd
.dina_fit_censoring_cox <- function(X_train, Y_train, W_train,
                                    eta0_train, eta1_train,
                                    cens_type, cens_params, n_grid) {

  if (is.null(cens_type)) {
    return(list(
      censor0_fn = function(X_new, eta_new) rep_len(1, nrow(X_new)),
      censor1_fn = function(X_new, eta_new) rep_len(1, nrow(X_new))
    ))
  }

  if (identical(cens_type, "unknown")) {
    # Estimate P(NOT censored | X) directly by logistic regression
    train_df <- as.data.frame(X_train)
    train_df$.S <- Y_train$status
    fit <- stats::glm(
      .S ~ ., data = train_df, family = stats::binomial(link = "logit")
    )
    pred_fn <- function(X_new, eta_new) {
      new_df <- as.data.frame(X_new)
      names(new_df) <- setdiff(names(train_df), ".S")
      as.numeric(stats::predict(fit, newdata = new_df, type = "response"))
    }
    return(list(censor0_fn = pred_fn, censor1_fn = pred_fn))
  }

  # For "fixed" / "uniform" / "exponential" we use per-arm Breslow Lambda
  build_arm <- function(eta_arm, time_arm, status_arm) {

    if (identical(cens_type, "fixed")) {
      if (is.null(cens_params$T0)) {
        stop("`cens_params$T0` is required for cens_type = \"fixed\".")
      }
      T0 <- cens_params$T0
      Lambda_T0 <- .dina_base_cumhaz(eta_arm, time_arm, status_arm, T0)
      # P(C >= Y) = F_Y(T0) = 1 - exp(-Lambda(T0) * exp(eta_new))
      function(X_new, eta_new) {
        1 - exp(-Lambda_T0 * exp(eta_new))
      }

    } else if (identical(cens_type, "uniform")) {
      if (is.null(cens_params$T0)) {
        stop("`cens_params$T0` is required for cens_type = \"uniform\".")
      }
      T0 <- cens_params$T0
      t_grid <- seq.int(0, T0, length.out = n_grid + 1L)[-1L]
      dt <- diff(c(0, t_grid))
      f_C <- rep_len(1 / T0, length(t_grid))
      Lambda_grid <- .dina_base_cumhaz(eta_arm, time_arm, status_arm, t_grid)
      # P(C >= Y | x) = sum_k F_Y(t_k | x) f_C(t_k) dt_k
      function(X_new, eta_new) {
        # F_Y at grid for each new observation: (n_new x n_grid)
        exp_eta <- exp(eta_new)
        F_Y <- 1 - exp(-outer(exp_eta, Lambda_grid))
        as.numeric(F_Y %*% (f_C * dt))
      }

    } else if (identical(cens_type, "exponential")) {
      if (is.null(cens_params$rate_c)) {
        stop("`cens_params$rate_c` is required for cens_type = \"exponential\".")
      }
      rate_c <- cens_params$rate_c
      t_max <- max(time_arm) * 2
      t_grid <- seq.int(0, t_max, length.out = n_grid + 1L)[-1L]
      dt <- diff(c(0, t_grid))
      f_C <- rate_c * exp(-rate_c * t_grid)
      Lambda_grid <- .dina_base_cumhaz(eta_arm, time_arm, status_arm, t_grid)
      function(X_new, eta_new) {
        exp_eta <- exp(eta_new)
        F_Y <- 1 - exp(-outer(exp_eta, Lambda_grid))
        as.numeric(F_Y %*% (f_C * dt))
      }

    } else {
      stop("Unknown `cens_type`: ", cens_type, ".")
    }
  }

  list(
    censor0_fn = build_arm(eta0_train,
                           Y_train$time[W_train == 0],
                           Y_train$status[W_train == 0]),
    censor1_fn = build_arm(eta1_train,
                           Y_train$time[W_train == 1],
                           Y_train$status[W_train == 1])
  )
}

#' Breslow estimator of the baseline cumulative hazard, evaluated at `t_eval`.
#'
#' Step-function (constant-between-events) Breslow estimator with linear
#' interpolation between event times for smoothness; values beyond the
#' observed range are held constant at the nearest extreme.
#'
#' @noRd
.dina_base_cumhaz <- function(eta_train, time_train, status_train, t_eval) {
  ord <- order(time_train)
  t_sorted <- time_train[ord]
  s_sorted <- status_train[ord]
  exp_eta <- exp(eta_train[ord])
  # Risk-set sums sum_{l: t_l >= t_i} exp(eta_l), descending then reversed
  risk_set_sum <- rev(cumsum(rev(exp_eta)))
  # Increments at sorted times: 0 at censored, status / risk-set-sum at events
  dLambda <- s_sorted / pmax(risk_set_sum, .Machine$double.eps)
  Lambda_sorted <- cumsum(dLambda)
  # Linear interpolation, hold at extremes outside the observed range
  stats::approx(
    x = t_sorted, y = Lambda_sorted,
    xout = t_eval, method = "linear", rule = 2, ties = "mean"
  )$y
}


# ---------------------------------------------------------------------------
# Internal: combine nuisance estimators into (a(x), nu(x)) for test data
# ---------------------------------------------------------------------------

#' Build a closure that evaluates (a(x), nu(x)) on a new covariate matrix
#' for GLM families.
#'
#' @noRd
.dina_combine_nuisance_glm <- function(family,
                                       prop_fn, eta0_fn, eta1_fn,
                                       nu_fn, link_fn, var_fn, eps) {
  function(X_new) {
    e_x <- .dina_clip01(prop_fn(X_new))

    if (family == "gaussian") {
      a_x  <- e_x
      nu_x <- nu_fn(X_new)
    } else {
      eta0 <- eta0_fn(X_new)
      eta1 <- eta1_fn(X_new)
      mu0 <- switch(family,
                    binomial = stats::plogis(eta0),
                    poisson  = exp(eta0))
      mu1 <- switch(family,
                    binomial = stats::plogis(eta1),
                    poisson  = exp(eta1))
      V0 <- var_fn(mu0)
      V1 <- var_fn(mu1)
      num  <- e_x * V1
      den  <- num + (1 - e_x) * V0 + eps
      a_x  <- num / den
      nu_x <- a_x * eta1 + (1 - a_x) * eta0
    }
    list(a = a_x, nu = nu_x)
  }
}

#' Build a closure that evaluates (a(x), nu(x)) on a new covariate matrix
#' for the Cox family.
#'
#' Following Eq. (18) of Gao and Hastie (2025), a(x) uses the not-censoring
#' probabilities in place of the variance functions of the exponential
#' family case.
#'
#' @noRd
.dina_combine_nuisance_cox <- function(prop_fn, eta0_fn, eta1_fn,
                                       censor0_fn, censor1_fn, eps) {
  function(X_new) {
    e_x  <- .dina_clip01(prop_fn(X_new))
    eta0 <- eta0_fn(X_new)
    eta1 <- eta1_fn(X_new)
    p0_nc <- censor0_fn(X_new, eta0)
    p1_nc <- censor1_fn(X_new, eta1)
    num <- e_x * p1_nc
    den <- num + (1 - e_x) * p0_nc + eps
    a_x  <- num / den
    nu_x <- a_x * eta1 + (1 - a_x) * eta0
    list(a = a_x, nu = nu_x)
  }
}


# ---------------------------------------------------------------------------
# Internal: step 2 (target parameter) for GLM families
# ---------------------------------------------------------------------------

#' GLM step 2: fit
#'   Y ~ offset(nu) + resW + X:resW - 1
#' (with family-appropriate link), where resW = W - a(X).
#'
#' Returns beta plus sandwich components Sigma_1 (sum of (resW^2 V(mu))
#' x x') and Sigma_2 (sum of (resW (Y - mu))^2 x x'), each divided by
#' n_test.
#'
#' @noRd
.dina_step2_glm <- function(X_test, Y_test, W_test, family, nuisance_fn) {

  nu_a <- nuisance_fn(X_test)
  res_W <- W_test - nu_a$a
  n_t <- nrow(X_test)
  d <- ncol(X_test)

  X_aug <- cbind(`(Intercept)` = 1, X_test)            # n_t x (d+1)
  X_resW <- X_aug * res_W                              # n_t x (d+1)

  glm_family <- switch(
    family,
    gaussian = stats::gaussian(link = "identity"),
    binomial = stats::binomial(link = "logit"),
    poisson  = stats::poisson(link = "log")
  )

  # Fit: outcome Y, predictor matrix X_resW (no intercept), offset nu
  fit <- stats::glm.fit(
    x = X_resW, y = Y_test, offset = nu_a$nu,
    family = glm_family,
    intercept = FALSE
  )
  if (isFALSE(fit$converged)) {
    warning("GLM step 2 did not converge for family = \"", family,
            "\".  Results may be unreliable.", call. = FALSE)
  }
  beta_hat <- stats::setNames(fit$coefficients, colnames(X_aug))

  # Residuals on response and pearson scales
  mu_hat <- fit$fitted.values
  resp_res <- Y_test - mu_hat
  V_mu <- glm_family$variance(mu_hat)
  V_mu <- pmax(V_mu, .Machine$double.eps)

  # Sigma_1 = X' diag(resW^2 V(mu)) X / n  (bread; sample mean of score derivative)
  # Sigma_2 = X' diag(resW^2 (Y - mu)^2) X / n  (meat; sample mean of score outer product)
  W_bread <- res_W^2 * V_mu
  W_meat  <- (res_W * resp_res)^2
  sigma1 <- crossprod(X_aug, W_bread * X_aug) / n_t
  sigma2 <- crossprod(X_aug, W_meat  * X_aug) / n_t

  list(beta = beta_hat, sigma1 = sigma1, sigma2 = sigma2)
}


# ---------------------------------------------------------------------------
# Internal: step 2 (target parameter) for the Cox family
# ---------------------------------------------------------------------------

#' Cox step 2: fit
#'   coxph(Surv(time, status) ~ resW + X:resW + offset(nu) - 1)
#'
#' Returns beta plus sandwich components Sigma_1 = Sigma_2 = solve(coxph$var)
#' / n_test, so that the sandwich algebra in `dina_fit()` collapses to
#' coxph$var (the model-based information variance).
#'
#' @noRd
.dina_step2_cox <- function(X_test, Y_test, W_test, nuisance_fn) {

  nu_a <- nuisance_fn(X_test)
  res_W <- W_test - nu_a$a
  n_t <- nrow(X_test)
  d <- ncol(X_test)

  # Use formula-safe internal names; recover user-facing names via setNames()
  pretty_names <- c("(Intercept)", colnames(X_test))
  safe_names <- c("b0", paste0("b", seq_len(d)))
  X_aug <- cbind(1, X_test)
  colnames(X_aug) <- safe_names
  X_resW <- X_aug * res_W
  colnames(X_resW) <- paste0("resW_", safe_names)

  df <- data.frame(
    .time = Y_test$time,
    .status = Y_test$status,
    X_resW,
    .nu = nu_a$nu,
    check.names = FALSE
  )
  rhs <- paste(colnames(X_resW), collapse = " + ")
  fmla <- stats::as.formula(paste(
    "survival::Surv(.time, .status) ~", rhs, "+ offset(.nu) - 1"
  ))

  fit <- survival::coxph(fmla, data = df)
  beta_hat <- stats::setNames(stats::coef(fit), pretty_names)

  # Information-matrix variance from coxph; ensure the algebra in dina_fit()
  # recovers it: Sigma_1 = Sigma_2 = solve(coxph$var) / n  =>
  # Sigma_1^-1 Sigma_2 Sigma_1^-1 / n = coxph$var.
  info_over_n <- tryCatch(
    solve(fit$var) / n_t,
    error = function(e) {
      warning("Cox information matrix is singular: ",
              conditionMessage(e), call. = FALSE)
      matrix(NA_real_, length(beta_hat), length(beta_hat))
    }
  )
  list(beta = beta_hat, sigma1 = info_over_n, sigma2 = info_over_n)
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @export
coef.dina <- function(object, ...) {
  if (!inherits(object, "dina")) {
    stop("`object` must be of class \"dina\".")
  }
  object$coefficients
}

#' @export
vcov.dina <- function(object, ...) {
  if (!inherits(object, "dina")) {
    stop("`object` must be of class \"dina\".")
  }
  object$vcov
}

#' Predicted DINA values for new covariate data.
#'
#' Returns the estimated treatment-effect contrast
#' \eqn{\hat\tau(x) = \hat\beta_0 + x^\top \hat\beta_{1:d}}, on the natural
#' parameter scale of the fitted family (mean difference for Gaussian, log
#' odds ratio for binomial, log mean ratio for Poisson, log hazard ratio for
#' Cox).
#'
#' @param object an object of class `"dina"`.
#' @param newdata a numeric matrix or data frame with the same number of
#'   columns as the training `X` and columns in the same order.
#' @param ... unused.
#' @return numeric vector of predicted treatment-effect contrasts, one per
#'   row of `newdata`.
#' @export
predict.dina <- function(object, newdata, ...) {
  if (!inherits(object, "dina")) {
    stop("`object` must be of class \"dina\".")
  }
  if (missing(newdata)) {
    stop("`newdata` must be supplied.")
  }
  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
    stop("`newdata` must be a matrix or data frame.")
  }
  newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"
  if (ncol(newdata) != object$d) {
    stop("`newdata` must have ", object$d, " columns, got ",
         ncol(newdata), ".")
  }
  as.numeric(cbind(1, newdata) %*% object$coefficients)
}

#' Wald confidence intervals for DINA coefficients.
#'
#' Computes normal-approximation confidence intervals
#' \eqn{\hat\beta_j \pm z_{1 - \alpha/2} \sqrt{\widehat{V}_{jj}}} using the
#' sandwich (or Cox model-based) variance returned by `dina_fit()`.
#'
#' @param object an object of class `"dina"`.
#' @param parm optional character or integer vector selecting which
#'   coefficients to return intervals for; defaults to all.
#' @param level confidence level (default `0.95`).
#' @param ... unused.
#' @return a matrix with two columns giving lower and upper limits, named
#'   in standard `confint()` style (e.g. `"2.5 %"`, `"97.5 %"`).
#' @export
confint.dina <- function(object, parm, level = 0.95, ...) {
  if (!inherits(object, "dina")) {
    stop("`object` must be of class \"dina\".")
  }
  if (length(level) != 1L || is.na(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).")
  }
  est <- object$coefficients
  se  <- sqrt(diag(object$vcov))
  z   <- stats::qnorm((1 + level) / 2)
  result <- cbind(est - z * se, est + z * se)
  rownames(result) <- names(est)
  alpha <- (1 - level) / 2
  colnames(result) <- paste0(
    format(100 * c(alpha, 1 - alpha), trim = TRUE, nsmall = 1L), " %"
  )
  if (!missing(parm) && !is.null(parm)) {
    result <- result[parm, , drop = FALSE]
  }
  result
}

#' Summary method for a DINA fit.
#'
#' Builds a coefficient table with estimate, standard error, Wald z-statistic,
#' and two-sided p-value for each coefficient.
#'
#' @param object an object of class `"dina"`.
#' @param ... unused.
#' @return an object of class `"summary.dina"`, printable via
#'   `print()`.
#' @export
summary.dina <- function(object, ...) {
  if (!inherits(object, "dina")) {
    stop("`object` must be of class \"dina\".")
  }
  est <- object$coefficients
  se  <- sqrt(diag(object$vcov))
  z   <- est / se
  p   <- 2 * stats::pnorm(-abs(z))
  coef_tab <- cbind(
    Estimate    = est,
    `Std. Error` = se,
    `z value`   = z,
    `Pr(>|z|)`  = p
  )
  out <- list(
    call          = object$call,
    family        = object$family,
    coefficients  = coef_tab,
    n             = object$n,
    cross_fitting = object$cross_fitting,
    n_folds       = object$n_folds
  )
  class(out) <- "summary.dina"
  out
}

#' @export
print.dina <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("DINA fit (family = \"", x$family, "\", n = ", x$n, ")\n", sep = "")
  if (x$cross_fitting) {
    cat("Cross-fitting: ", x$n_folds, " folds\n", sep = "")
  }
  cat("\nCoefficients:\n")
  print.default(
    format(x$coefficients, digits = digits),
    print.gap = 2L, quote = FALSE
  )
  invisible(x)
}

#' @export
print.summary.dina <- function(x,
                               digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  cat("DINA fit (family = \"", x$family, "\", n = ", x$n, ")\n", sep = "")
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }
  cat("\nCoefficients:\n")
  stats::printCoefmat(
    x$coefficients,
    digits = digits,
    signif.stars = signif.stars,
    has.Pvalue = TRUE
  )
  if (isTRUE(x$cross_fitting)) {
    cat("\nCross-fitting: ", x$n_folds, " folds\n", sep = "")
  }
  invisible(x)
}
