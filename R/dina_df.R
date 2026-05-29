# File: R/dina_df.R
# Part of the forestsearch package.
#
# Data-frame + column-name interface wrappers for DINA estimation,
# matching the conventions of forestsearch().  These are the high-level
# user-facing API; the low-level matrix-based functions dina_fit() and
# dina_fit_bagged() are kept for power-user / internal use.
#
# The pattern mirrors base R's glm() (high-level, formula/data) versus
# glm.fit() (low-level, matrix x and vector y).


# ---------------------------------------------------------------------------
# dina() -- data-frame wrapper around dina_fit()
# ---------------------------------------------------------------------------

#' Estimate heterogeneous treatment effects via DINA (data-frame interface)
#'
#' `dina()` is the high-level user-facing wrapper around [dina_fit()].
#' It accepts a data frame plus column-name arguments, matching the
#' conventions used elsewhere in the `forestsearch` package, and
#' extracts the covariate matrix, response, and treatment indicator
#' before dispatching to [dina_fit()].
#'
#' All other arguments are forwarded to [dina_fit()] via `...`,
#' including `propensity_method`, `baseline_method`, `cross_fitting`,
#' `n_folds`, `cens_type`, `cens_params`, `n_grid`, `eps`, and `seed`.
#' See [dina_fit()] for their full documentation.
#'
#' Covariate columns must be numeric (or integer).  Convert factors to
#' dummy variables (e.g., with [stats::model.matrix()]) before calling
#' if needed.  Rows with any `NA` value in the referenced columns are
#' rejected with an error; use [stats::na.omit()] or impute beforehand.
#'
#' @param df data frame (or tibble / `data.table`) containing the
#'   outcome, treatment, and covariate columns.
#' @param outcome character(1); the name of the response column in
#'   `df`.  For `family = "cox"` this is the event/censoring time
#'   column.
#' @param treatment character(1); the name of the binary treatment
#'   indicator column in `df`.  Values must be 0/1 or logical.
#' @param covariates character vector of length `>= 1`; the names of
#'   the covariate columns in `df`.  All must be numeric.
#' @param family one of `"gaussian"`, `"binomial"`, `"poisson"`, `"cox"`.
#' @param status character(1) or `NULL`; for `family = "cox"` only, the
#'   name of the event indicator column in `df` (`1` = event, `0` =
#'   censored).  Ignored for non-Cox families.
#' @param ... further arguments forwarded to [dina_fit()].
#'
#' @return An object of class `"dina"`, as returned by [dina_fit()].
#'   The `$call` component is overwritten with the original `dina()`
#'   call so that `print()` and `summary()` show the high-level
#'   invocation rather than the internal matrix-based dispatch.
#'
#' @seealso [dina_fit()] for the underlying matrix-based estimator;
#'   [dina_bagged()] for the bagged version with infinitesimal-jackknife
#'   variance.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' df_demo <- data.frame(
#'   y         = NA,
#'   w         = stats::rbinom(n, 1, 0.5),
#'   x1        = stats::runif(n, -1, 1),
#'   x2        = stats::runif(n, -1, 1),
#'   x3        = stats::runif(n, -1, 1)
#' )
#' tau_x      <- 0.5 + 0.8 * df_demo$x1 - 0.3 * df_demo$x2
#' df_demo$y  <- 1 + df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)
#'
#' fit <- dina(
#'   df         = df_demo,
#'   outcome    = "y",
#'   treatment  = "w",
#'   covariates = c("x1", "x2", "x3"),
#'   family     = "gaussian",
#'   seed       = 1L
#' )
#' coef(fit)
#' confint(fit)
#'
#' @export
dina <- function(df,
                 outcome,
                 treatment,
                 covariates,
                 family = c("gaussian", "binomial", "poisson", "cox"),
                 status = NULL,
                 ...) {

  family <- match.arg(family)
  parts <- .dina_extract_from_df(
    df = df, outcome = outcome, treatment = treatment,
    covariates = covariates, family = family, status = status
  )

  fit <- dina_fit(
    X = parts$X, Y = parts$Y, W = parts$W,
    family = family,
    ...
  )
  fit$call <- match.call()
  fit
}


# ---------------------------------------------------------------------------
# dina_bagged() -- data-frame wrapper around dina_fit_bagged()
# ---------------------------------------------------------------------------

#' Bagged DINA with IJ variance (data-frame interface)
#'
#' `dina_bagged()` is the high-level user-facing wrapper around
#' [dina_fit_bagged()].  It accepts a data frame plus column-name
#' arguments and dispatches to [dina_fit_bagged()], which fits the
#' DINA estimator on `n_bags` bootstrap replicates and returns the
#' infinitesimal-jackknife variance per Wager, Hastie and Efron (2014).
#'
#' All bagging-related arguments (`n_bags`, `cross_fitting_per_bag`,
#' `n_folds_per_bag`, `parallel`, `ij_finite_sample_correction`,
#' `project_psd`, `verbose`, ...) are forwarded via `...`.  See
#' [dina_fit_bagged()] for full documentation.
#'
#' For `parallel = "bags"`, set `future::plan()` to a parallel backend
#' *before* calling `dina_bagged()`.  Restore with
#' `future::plan("sequential"); gc()` afterwards.
#'
#' @inheritParams dina
#' @param ... further arguments forwarded to [dina_fit_bagged()].
#'
#' @return An object of class `c("dina_bagged", "dina")`, as returned
#'   by [dina_fit_bagged()].  The `$call` component is overwritten with
#'   the original `dina_bagged()` call.
#'
#' @seealso [dina_fit_bagged()] for the underlying matrix-based bagged
#'   estimator; [dina()] for the single-pass DINA with sandwich variance.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' df_demo <- data.frame(
#'   y         = NA,
#'   w         = stats::rbinom(n, 1, 0.5),
#'   x1        = stats::runif(n, -1, 1),
#'   x2        = stats::runif(n, -1, 1),
#'   x3        = stats::runif(n, -1, 1)
#' )
#' tau_x      <- 0.5 + 0.8 * df_demo$x1 - 0.3 * df_demo$x2
#' df_demo$y  <- 1 + df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)
#'
#' fit <- dina_bagged(
#'   df         = df_demo,
#'   outcome    = "y",
#'   treatment  = "w",
#'   covariates = c("x1", "x2", "x3"),
#'   family     = "gaussian",
#'   n_bags     = 30L,
#'   seed       = 1L
#' )
#' coef(fit)
#' sqrt(diag(vcov(fit)))
#'
#' @export
dina_bagged <- function(df,
                        outcome,
                        treatment,
                        covariates,
                        family = c("gaussian", "binomial", "poisson", "cox"),
                        status = NULL,
                        ...) {

  family <- match.arg(family)
  parts <- .dina_extract_from_df(
    df = df, outcome = outcome, treatment = treatment,
    covariates = covariates, family = family, status = status
  )

  fit <- dina_fit_bagged(
    X = parts$X, Y = parts$Y, W = parts$W,
    family = family,
    ...
  )
  fit$call <- match.call()
  fit
}


# ---------------------------------------------------------------------------
# Internal helper: data-frame extraction and validation
# ---------------------------------------------------------------------------

#' Extract `X`, `Y`, `W` from a data frame given column names.
#'
#' Centralized validation for the `dina()` and `dina_bagged()`
#' wrappers.  Returns a list with components `X` (numeric matrix),
#' `Y` (numeric vector for non-Cox families; data frame with `time`
#' and `status` for Cox), and `W` (numeric 0/1 vector).
#'
#' Treats logical `treatment` as 0/1; anything else with non-0/1
#' values triggers a downstream error in the called `dina_fit()` or
#' `dina_fit_bagged()`.  Covariates must be numeric; factors are
#' rejected with a clear message.
#'
#' @noRd
.dina_extract_from_df <- function(df, outcome, treatment, covariates,
                                  family, status) {

  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  }

  if (!is.character(outcome) || length(outcome) != 1L || is.na(outcome)) {
    stop("`outcome` must be a single non-NA character string.")
  }
  if (!is.character(treatment) || length(treatment) != 1L || is.na(treatment)) {
    stop("`treatment` must be a single non-NA character string.")
  }
  if (!is.character(covariates) || length(covariates) < 1L ||
      anyNA(covariates)) {
    stop("`covariates` must be a character vector of length >= 1 with no NAs.")
  }
  if (family == "cox") {
    if (is.null(status)) {
      stop("`status` (event indicator column name) is required for ",
           "family = \"cox\".")
    }
    if (!is.character(status) || length(status) != 1L || is.na(status)) {
      stop("`status` must be a single non-NA character string.")
    }
  } else if (!is.null(status)) {
    warning("`status` is ignored for family = \"", family,
            "\"; only used when family = \"cox\".", call. = FALSE)
  }

  required_cols <- c(outcome, treatment, covariates)
  if (family == "cox") required_cols <- c(required_cols, status)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `df`: ",
         paste(missing_cols, collapse = ", "), ".")
  }

  if (anyNA(df[, required_cols])) {
    stop("`df` contains NA values in the columns referenced by ",
         "`outcome`, `treatment`, `covariates`",
         if (family == "cox") ", or `status`" else "",
         ".  Use `na.omit()` or impute before calling.")
  }

  # ---- Treatment ---------------------------------------------------------
  W <- df[[treatment]]
  if (is.logical(W)) W <- as.integer(W)

  # ---- Covariates --------------------------------------------------------
  cov_df <- df[, covariates, drop = FALSE]
  is_num <- vapply(cov_df, is.numeric, logical(1L))
  if (!all(is_num)) {
    non_num <- covariates[!is_num]
    stop("Covariate column(s) must be numeric: ",
         paste(non_num, collapse = ", "),
         ".  Convert factors to dummy variables (e.g., via ",
         "`stats::model.matrix()`) before calling.")
  }
  X <- as.matrix(cov_df)
  storage.mode(X) <- "double"

  # ---- Response ----------------------------------------------------------
  Y <- if (family == "cox") {
    data.frame(
      time   = as.numeric(df[[outcome]]),
      status = as.integer(df[[status]])
    )
  } else {
    as.numeric(df[[outcome]])
  }

  list(X = X, Y = Y, W = W)
}
