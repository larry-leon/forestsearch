# =============================================================================
# glm_simulate.R
#
# Data-generating mechanisms for binary, continuous, and rate-based GLM
# outcomes.  Used in testing, simulation studies, and vignettes.
# =============================================================================


#' Simulate a Two-Arm Randomized Trial with a Binary Outcome
#'
#' Generates a randomized trial dataset where a pre-specified subgroup has a
#' detrimental (or differential) treatment effect on a binary endpoint.  The
#' primary estimand is the **risk difference** (RD) in each subgroup.
#'
#' The data-generating model is a logistic regression with an interaction
#' between treatment and subgroup membership:
#'
#' \deqn{\mathrm{logit}(P(Y=1)) = \mu_0 + \alpha A + \beta_Z Z +
#'       \gamma A \mathbf{1}(Z \in H)}
#'
#' where \eqn{H} is the true harmful subgroup, \eqn{\gamma > 0} represents
#' additional log-odds of harm in \eqn{H}, and \eqn{A} is the treatment
#' indicator.
#'
#' @param n               Integer. Total sample size (default 500).
#' @param treat_prob      Numeric in (0, 1). Probability of treatment
#'   allocation (default 0.5 for 1:1 randomisation).
#' @param n_covariates    Integer. Number of binary baseline covariates
#'   (default 5).
#' @param subgroup_defn   A named list defining the true harmful subgroup,
#'   e.g., `list(x1 = 1, x2 = 0)` means the subgroup
#'   \eqn{\{x_1 = 1\} \cap \{x_2 = 0\}}.
#'   `NULL` for no subgroup effect (null DGM).
#' @param baseline_risk   Numeric. Baseline event probability in the control
#'   arm for subjects not in the subgroup (default 0.20).
#' @param itt_log_or      Numeric. Log-odds ratio for treatment in the ITT
#'   population (default 0 --- no overall treatment effect).
#' @param subgroup_log_or Numeric. **Additional** log-odds ratio for treatment
#'   in the subgroup (default `log(2)` --- OR of 2 for harm in subgroup).
#' @param prog_log_or     Numeric vector of length `n_covariates`.  Log-odds
#'   ratios for each covariate on the outcome (prognostic effects).  Default
#'   is 0.5 for the first two covariates, 0 for the rest.
#' @param seed            Integer or `NULL`. Random seed.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{`id`}{Subject identifier.}
#'   \item{`treat`}{Binary treatment indicator (0 = control, 1 = experimental).}
#'   \item{`x1`, ..., `xK`}{Binary baseline covariates.}
#'   \item{`response`}{Binary outcome (0/1).}
#'   \item{`in_subgroup`}{Logical. `TRUE` if the subject is in the true
#'     subgroup.}
#'   \item{`true_prob`}{Numeric.  True event probability for this subject
#'     (useful for verifying DGM).}
#' }
#'
#' @examples
#' \dontrun{
#'   df <- simulate_binary_trial(
#'     n             = 600,
#'     subgroup_defn = list(x1 = 1),
#'     baseline_risk  = 0.25,
#'     subgroup_log_or = log(3),
#'     seed           = 42
#'   )
#'   # Verify risk differences
#'   with(df[df$in_subgroup, ],
#'        tapply(response, treat, mean) |> diff())
#'   with(df[!df$in_subgroup, ],
#'        tapply(response, treat, mean) |> diff())
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
  prog_log_or <- .pad_coef(prog_log_or, n_covariates)

  cov_names <- paste0("x", seq_len(n_covariates))

  # Simulate binary covariates (prevalence 0.5 each)
  X <- vapply(
    seq_len(n_covariates),
    function(j) stats::rbinom(n, 1L, 0.5),
    integer(n)
  )
  colnames(X) <- cov_names
  df <- as.data.frame(X)

  df$treat <- stats::rbinom(n, 1L, treat_prob)

  # Subgroup membership
  in_sg <- .eval_binary_subgroup(df, subgroup_defn)
  df$in_subgroup <- in_sg

  # Linear predictor on logit scale
  mu0 <- stats::qlogis(baseline_risk)
  lp  <- mu0 +
    itt_log_or        * df$treat +
    subgroup_log_or   * df$treat * as.integer(in_sg) +
    as.matrix(df[, cov_names, drop = FALSE]) %*% prog_log_or

  prob <- stats::plogis(as.numeric(lp))
  df$response  <- stats::rbinom(n, 1L, prob)
  df$true_prob <- prob
  df$id        <- seq_len(n)

  df[, c("id", "treat", cov_names, "response", "in_subgroup", "true_prob")]
}


#' Simulate a Two-Arm Randomized Trial with a Continuous Outcome
#'
#' Generates a randomized trial dataset where a pre-specified subgroup has a
#' differential treatment effect on a continuous (normally distributed)
#' endpoint.  The primary estimand is the **mean difference** (MD).
#'
#' @param n               Integer. Total sample size (default 500).
#' @param treat_prob      Numeric. Probability of treatment allocation
#'   (default 0.5).
#' @param n_covariates    Integer. Number of continuous baseline covariates
#'   (default 5, simulated from N(0,1)).
#' @param subgroup_defn   A named list defining the true harmful subgroup
#'   using covariate cut rules, e.g.,
#'   `list(x1 = list(op = "<=", val = 0))` for \eqn{\{x_1 \le 0\}}.
#'   `NULL` for no subgroup effect.
#' @param mu0             Numeric. Intercept (mean outcome in control arm
#'   with all covariates at 0; default 50).
#' @param itt_md          Numeric. Overall (ITT) mean difference
#'   (experimental minus control; default 0).
#' @param subgroup_md     Numeric. **Additional** mean difference in the
#'   subgroup (positive = worse outcome in experimental arm; default 5).
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
#'   \item{`in_subgroup`}{Logical.}
#' }
#'
#' @examples
#' \dontrun{
#'   df <- simulate_continuous_trial(
#'     n             = 500,
#'     subgroup_defn = list(x1 = list(op = "<=", val = 0)),
#'     subgroup_md   = 8,
#'     sigma         = 12,
#'     seed          = 42
#'   )
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
  prog_coef <- .pad_coef(prog_coef, n_covariates)

  cov_names <- paste0("x", seq_len(n_covariates))
  X <- matrix(
    stats::rnorm(n * n_covariates),
    nrow = n, ncol = n_covariates,
    dimnames = list(NULL, cov_names)
  )
  df        <- as.data.frame(X)
  df$treat  <- stats::rbinom(n, 1L, treat_prob)

  # Subgroup membership via cut rules
  in_sg <- .eval_continuous_subgroup(df, subgroup_defn)
  df$in_subgroup <- in_sg

  df$outcome <- mu0 +
    itt_md      * df$treat +
    subgroup_md * df$treat * as.integer(in_sg) +
    X %*% prog_coef +
    stats::rnorm(n, 0, sigma)

  df$id <- seq_len(n)
  df[, c("id", "treat", cov_names, "outcome", "in_subgroup")]
}


#' Simulate a Two-Arm Randomized Trial with Survival-Like Rate Data
#'
#' Generates a dataset with both a **binary event indicator** and
#' **follow-up time**, suitable for both Cox regression and Poisson
#' regression with a log-time offset.  Under the exponential baseline
#' hazard used here, the Poisson-offset log-IRR equals the Cox log-HR
#' exactly, providing a clean sanity check.
#'
#' The data-generating mechanism is an exponential survival model:
#'
#' \deqn{h_i(t) = \lambda_0 \exp(\alpha A_i + \gamma A_i
#'       \mathbf{1}(i \in H) + Z_i^\top \beta_Z)}
#'
#' with administrative censoring \eqn{C_i \sim \text{Unif}(c_{lo},
#' c_{hi})}.  The columns `event` (binary) and `tte` (time) are provided
#' for both analysis paths.
#'
#' @param n               Integer. Total sample size (default 700).
#' @param treat_prob      Numeric. Treatment allocation probability
#'   (default 0.5).
#' @param n_covariates    Integer. Number of binary baseline covariates
#'   (default 5).
#' @param subgroup_defn   Named list defining the true subgroup.  Same
#'   syntax as [simulate_binary_trial()].
#' @param baseline_rate   Numeric. Baseline event rate (hazard) in the
#'   control arm, non-subgroup subjects (default 0.10).
#' @param itt_log_irr     Numeric. Log incidence-rate ratio for treatment
#'   in the ITT population (default `log(0.80)` --- 20\% rate reduction).
#' @param subgroup_log_irr Numeric. **Additional** log-IRR for treatment in
#'   the subgroup (default `log(2.5)` --- substantially elevated rate in
#'   subgroup, yielding overall harm).
#' @param prog_log_hr     Numeric vector. Prognostic log-hazard-ratios for
#'   covariates (default 0.3 for first two, 0 for rest).
#' @param cens_range      Numeric vector of length 2. Uniform censoring
#'   support `c(c_lo, c_hi)` (default `c(2, 14)` years, matching the
#'   Leon et al. simulation setup).
#' @param seed            Integer or `NULL`. Random seed.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{`id`}{Subject identifier.}
#'   \item{`treat`}{Binary treatment indicator.}
#'   \item{`x1`, ..., `xK`}{Binary baseline covariates.}
#'   \item{`event`}{Binary event indicator (1 = event observed, 0 = censored).}
#'   \item{`tte`}{Numeric follow-up time (min of event time and censoring time).}
#'   \item{`in_subgroup`}{Logical.}
#'   \item{`true_log_irr`}{Numeric. True individual log-IRR (= log-HR under
#'     the exponential model).}
#' }
#'
#' @details
#' The exponential baseline hazard ensures that the Poisson regression
#' `glm(event ~ treat, family = poisson, offset = log(tte))` yields
#' a log-IRR that converges to the Cox log-HR as sample size increases.
#' This makes `simulate_rate_trial()` the primary tool for validating
#' the Poisson-offset path against the existing survival pipeline.
#'
#' @examples
#' \dontrun{
#'   df <- simulate_rate_trial(n = 700, seed = 42)
#'
#'   # Cox log-HR
#'   fit_cox <- survival::coxph(
#'     survival::Surv(tte, event) ~ treat, data = df
#'   )
#'   coef(fit_cox)["treat"]
#'
#'   # Poisson-offset log-IRR (should approximate the above)
#'   fit_pois <- glm(
#'     event ~ treat, family = poisson(link = "log"),
#'     offset = log(tte), data = df
#'   )
#'   coef(fit_pois)["treat"]
#' }
#'
#' @export
simulate_rate_trial <- function(
    n                  = 700L,
    treat_prob         = 0.5,
    n_covariates       = 5L,
    subgroup_defn      = list(x1 = 1, x3 = 1),
    baseline_rate      = 0.10,
    itt_log_irr        = log(0.80),
    subgroup_log_irr   = log(2.5),
    prog_log_hr        = NULL,
    cens_range         = c(2, 14),
    seed               = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  n            <- as.integer(n)
  n_covariates <- as.integer(n_covariates)

  if (is.null(prog_log_hr)) {
    prog_log_hr <- c(0.3, 0.3, rep(0, max(n_covariates - 2L, 0L)))
  }
  prog_log_hr <- .pad_coef(prog_log_hr, n_covariates)

  cov_names <- paste0("x", seq_len(n_covariates))

  # Binary covariates (prevalence 0.5)
  X <- vapply(
    seq_len(n_covariates),
    function(j) stats::rbinom(n, 1L, 0.5),
    integer(n)
  )
  colnames(X) <- cov_names
  df <- as.data.frame(X)

  df$treat <- stats::rbinom(n, 1L, treat_prob)

  # Subgroup membership
  in_sg <- .eval_binary_subgroup(df, subgroup_defn)
  df$in_subgroup <- in_sg

  # Individual log-IRR (= log-HR under exponential model)
  log_irr_i <- itt_log_irr + subgroup_log_irr * as.integer(in_sg)
  df$true_log_irr <- log_irr_i

  # Hazard for each subject
  lp_prog <- as.numeric(as.matrix(df[, cov_names, drop = FALSE]) %*%
                           prog_log_hr)
  lambda_i <- baseline_rate * exp(log_irr_i * df$treat + lp_prog)

  # Exponential event times (rate = lambda_i)
  event_time <- stats::rexp(n, rate = lambda_i)

  # Administrative censoring
  cens_time <- stats::runif(n, min = cens_range[1L], max = cens_range[2L])

  df$tte   <- pmin(event_time, cens_time)
  df$event <- as.integer(event_time <= cens_time)
  df$id    <- seq_len(n)

  df[, c("id", "treat", cov_names, "event", "tte",
         "in_subgroup", "true_log_irr")]
}


# ---------------------------------------------------------------------------
# Internal helpers for subgroup evaluation
# ---------------------------------------------------------------------------

#' Evaluate binary-covariate subgroup membership
#' @noRd
.eval_binary_subgroup <- function(df, subgroup_defn) {
  if (is.null(subgroup_defn) || length(subgroup_defn) == 0L) {
    return(rep(FALSE, nrow(df)))
  }
  Reduce(
    `&`,
    lapply(names(subgroup_defn), function(v) df[[v]] == subgroup_defn[[v]])
  )
}

#' Evaluate continuous-covariate subgroup membership (op/val rules)
#' @noRd
.eval_continuous_subgroup <- function(df, subgroup_defn) {
  if (is.null(subgroup_defn) || length(subgroup_defn) == 0L) {
    return(rep(FALSE, nrow(df)))
  }
  Reduce(
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
}

#' Zero-pad a coefficient vector to length k
#' @noRd
.pad_coef <- function(coef_vec, k) {
  if (length(coef_vec) < k) {
    c(coef_vec, rep(0, k - length(coef_vec)))
  } else {
    coef_vec[seq_len(k)]
  }
}
