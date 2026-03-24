# =============================================================================
# glm_simulate.R
#
# Data-generating mechanisms for binary and continuous GLM outcomes.
# Used in testing, simulation studies, and vignettes.
# =============================================================================


#' Simulate a Two-Arm Randomized Trial with a Binary Outcome
#'
#' Generates a randomized trial dataset where a pre-specified subgroup has a
#' detrimental (or differential) treatment effect on a binary endpoint.
#'
#' The data-generating model is a logistic regression with an interaction
#' between treatment and subgroup membership:
#'
#' \deqn{\text{logit}(P(Y=1)) = \mu_0 + \alpha \cdot A + \beta_Z \cdot Z +
#'       \gamma \cdot A \cdot \mathbf{1}(Z \in H)}
#'
#' where \eqn{H} is the true harmful subgroup, \eqn{\gamma > 0} represents
#' additional log-odds of harm in \eqn{H}, and \eqn{A} is the treatment
#' indicator.
#'
#' @param n               Integer. Total sample size (default 500).
#' @param treat_prob      Numeric in (0, 1). Probability of treatment allocation
#'   (default 0.5 for 1:1 randomisation).
#' @param n_covariates    Integer. Number of binary baseline covariates
#'   (default 5).
#' @param subgroup_defn   A named list defining the true harmful subgroup,
#'   e.g., `list(x1 = 1, x2 = 0)` means the subgroup `{x1 == 1} & {x2 == 0}`.
#'   `NULL` for no subgroup effect.
#' @param baseline_risk   Numeric. Baseline event probability in the control arm
#'   for subjects not in the subgroup (default 0.20).
#' @param itt_log_or      Numeric. Log-odds ratio for treatment in the ITT
#'   population (default 0 — no overall treatment effect).
#' @param subgroup_log_or Numeric. **Additional** log-odds ratio for treatment
#'   in the subgroup (default `log(2)` — OR of 2 for harm in subgroup).
#' @param prog_log_or     Numeric vector of length `n_covariates`.  Log-odds
#'   ratio for each covariate on the outcome (prognostic effects).  Default is
#'   0.5 for the first two covariates, 0 for the rest.
#' @param seed            Integer or `NULL`. Random seed.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{`id`}{Subject identifier.}
#'   \item{`treat`}{Binary treatment indicator (0 = control, 1 = experimental).}
#'   \item{`x1`, ..., `xK`}{Binary baseline covariates.}
#'   \item{`response`}{Binary outcome (0/1).}
#'   \item{`in_subgroup`}{Logical. `TRUE` if the subject is in the true subgroup.}
#' }
#'
#' @examples
#' \dontrun{
#'   df <- simulate_binary_trial(
#'     n             = 600,
#'     subgroup_defn = list(x1 = 1),
#'     baseline_risk  = 0.25,
#'     subgroup_log_or = log(3)
#'   )
#'   head(df)
#'   table(df$treat, df$in_subgroup, df$response)
#' }
#'
#' @export
simulate_binary_trial <- function(
    n               = 500L,
    treat_prob      = 0.5,
    n_covariates    = 5L,
    subgroup_defn   = list(x1 = 1),
    baseline_risk   = 0.20,
    itt_log_or      = 0,
    subgroup_log_or = log(2),
    prog_log_or     = NULL,
    seed            = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  n             <- as.integer(n)
  n_covariates  <- as.integer(n_covariates)

  # Default prognostic effects
  if (is.null(prog_log_or)) {
    prog_log_or <- c(0.5, 0.5, rep(0, max(n_covariates - 2L, 0L)))
  }
  if (length(prog_log_or) < n_covariates) {
    prog_log_or <- c(
      prog_log_or,
      rep(0, n_covariates - length(prog_log_or))
    )
  }

  # Covariate names
  cov_names <- paste0("x", seq_len(n_covariates))

  # Simulate covariates (binary, prevalence 0.5)
  X <- vapply(
    seq_len(n_covariates),
    function(j) stats::rbinom(n, 1L, 0.5),
    integer(n)
  )
  colnames(X) <- cov_names
  df <- as.data.frame(X)

  # Treatment assignment
  df$treat <- stats::rbinom(n, 1L, treat_prob)

  # Subgroup membership
  if (!is.null(subgroup_defn) && length(subgroup_defn) > 0L) {
    in_sg <- Reduce(
      `&`,
      lapply(names(subgroup_defn), function(v) df[[v]] == subgroup_defn[[v]])
    )
  } else {
    in_sg <- rep(FALSE, n)
  }
  df$in_subgroup <- in_sg

  # Linear predictor
  mu0 <- log(baseline_risk / (1 - baseline_risk))   # intercept
  lp  <- mu0 +
    itt_log_or        * df$treat +
    subgroup_log_or   * df$treat * as.integer(in_sg) +
    as.matrix(df[, cov_names, drop = FALSE]) %*% prog_log_or

  df$response <- stats::rbinom(n, 1L, plogis(lp))
  df$id       <- seq_len(n)

  # Reorder columns
  df[, c("id", "treat", cov_names, "response", "in_subgroup")]
}


#' Simulate a Two-Arm Randomized Trial with a Continuous Outcome
#'
#' Generates a randomized trial dataset where a pre-specified subgroup has a
#' detrimental treatment effect on a continuous (normally distributed) endpoint.
#'
#' The data-generating model is:
#'
#' \deqn{Y_i = \mu_0 + \alpha \cdot A_i + \beta_Z^T Z_i +
#'       \gamma \cdot A_i \cdot \mathbf{1}(Z_i \in H) + \varepsilon_i,
#'       \quad \varepsilon_i \sim N(0, \sigma^2)}
#'
#' @param n               Integer. Total sample size (default 500).
#' @param treat_prob      Numeric. Probability of treatment allocation (default 0.5).
#' @param n_covariates    Integer. Number of continuous baseline covariates
#'   (default 5, simulated from N(0,1)).
#' @param subgroup_defn   A named list defining the true harmful subgroup using
#'   covariate cut rules, e.g., `list(x1 = list(op = "<=", val = 0))` for
#'   `{x1 <= 0}`.  `NULL` for no subgroup effect.
#' @param mu0             Numeric. Intercept (mean outcome in control arm with
#'   all covariates at 0; default 50).
#' @param itt_md          Numeric. Overall (ITT) mean difference (experimental
#'   minus control; default 0).
#' @param subgroup_md     Numeric. **Additional** mean difference in the subgroup
#'   (harm: positive values indicate worse outcome in experimental arm;
#'   default 5).
#' @param prog_coef       Numeric vector of prognostic coefficients for
#'   covariates (default 1 for first two, 0 for the rest).
#' @param sigma           Numeric. Residual standard deviation (default 10).
#' @param seed            Integer or `NULL`. Random seed.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{`id`}{Subject identifier.}
#'   \item{`treat`}{Binary treatment indicator.}
#'   \item{`x1`, ..., `xK`}{Continuous baseline covariates (N(0,1)).}
#'   \item{`outcome`}{Continuous numeric outcome.}
#'   \item{`in_subgroup`}{Logical. `TRUE` if the subject is in the true subgroup.}
#' }
#'
#' @examples
#' \dontrun{
#'   df <- simulate_continuous_trial(
#'     n            = 500,
#'     subgroup_defn = list(x1 = list(op = "<=", val = 0)),
#'     subgroup_md   = 8,
#'     sigma         = 12
#'   )
#'   head(df)
#'   with(df, tapply(outcome, list(treat, in_subgroup), mean))
#' }
#'
#' @export
simulate_continuous_trial <- function(
    n             = 500L,
    treat_prob    = 0.5,
    n_covariates  = 5L,
    subgroup_defn = list(x1 = list(op = "<=", val = 0)),
    mu0           = 50,
    itt_md        = 0,
    subgroup_md   = 5,
    prog_coef     = NULL,
    sigma         = 10,
    seed          = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  n            <- as.integer(n)
  n_covariates <- as.integer(n_covariates)

  if (is.null(prog_coef)) {
    prog_coef <- c(1, 1, rep(0, max(n_covariates - 2L, 0L)))
  }
  if (length(prog_coef) < n_covariates) {
    prog_coef <- c(prog_coef, rep(0, n_covariates - length(prog_coef)))
  }

  cov_names <- paste0("x", seq_len(n_covariates))
  X <- matrix(
    stats::rnorm(n * n_covariates),
    nrow = n, ncol = n_covariates,
    dimnames = list(NULL, cov_names)
  )
  df        <- as.data.frame(X)
  df$treat  <- stats::rbinom(n, 1L, treat_prob)

  # Subgroup membership via cut rules
  if (!is.null(subgroup_defn) && length(subgroup_defn) > 0L) {
    in_sg <- Reduce(
      `&`,
      lapply(names(subgroup_defn), function(v) {
        rule <- subgroup_defn[[v]]
        switch(rule$op,
          "<="  = df[[v]] <= rule$val,
          ">="  = df[[v]] >= rule$val,
          "<"   = df[[v]] <  rule$val,
          ">"   = df[[v]] >  rule$val,
          "=="  = df[[v]] == rule$val,
          stop("Unknown op '", rule$op, "' in subgroup_defn.", call. = FALSE)
        )
      })
    )
  } else {
    in_sg <- rep(FALSE, n)
  }
  df$in_subgroup <- in_sg

  # Outcome
  df$outcome <- mu0 +
    itt_md      * df$treat +
    subgroup_md * df$treat * as.integer(in_sg) +
    X %*% prog_coef +
    stats::rnorm(n, 0, sigma)

  df$id <- seq_len(n)
  df[, c("id", "treat", cov_names, "outcome", "in_subgroup")]
}
