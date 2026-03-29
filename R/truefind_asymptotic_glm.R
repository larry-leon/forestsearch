# =============================================================================
# Asymptotic Detection Probability for GLM Outcomes
# =============================================================================
#
# Extends the Section 2.1 approximation (Leon et al., 2024, eq. 3)
# to arbitrary GLM outcomes via the effective information d_eff.
#
# The key identity: for any canonical GLM with balanced arms,
#   Var(beta_hat_full) = 4 / d_eff
#   Var(beta_hat_half) = 8 / d_eff
#
# where d_eff depends on the family:
#   Survival (Cox):    d_eff = n_H * (1 - prop_cens)       [events]
#   Poisson + offset:  d_eff = D                            [events]
#   Binomial (logit):  d_eff = n_H * p_bar * (1 - p_bar)   [Bernoulli info]
#   Gaussian (id):     d_eff = n_H / sigma_y^2             [precision-weighted n]
#   Quasi-Poisson:     d_eff = D / phi                     [dispersion-adjusted]
#
# Reference: Leon et al. (2024) "Exploratory subgroup identification
# in the heterogeneous Cox model", Statistics in Medicine, 43, 4253-4275.
#
# =============================================================================

#' @importFrom stats dnorm rnorm
NULL


# =============================================================================
# Effective Information Helpers
# =============================================================================

#' Effective Information for Survival (Cox PH)
#'
#' Computes \eqn{d_{\text{eff}} = n_H (1 - p_c)} where \eqn{p_c} is the
#' censoring proportion.
#'
#' @param n_sg Integer. Subgroup sample size.
#' @param prop_cens Numeric. Proportion censored (0--1).
#'
#' @return Numeric. Effective information (expected events).
#'
#' @examples
#' d_eff_survival(n_sg = 60, prop_cens = 0.45)
#' # 33
#'
#' @export
d_eff_survival <- function(n_sg, prop_cens) {
  stopifnot(n_sg > 0, prop_cens >= 0, prop_cens < 1)
  n_sg * (1 - prop_cens)
}


#' Effective Information for Binary Outcomes (Logistic)
#'
#' Computes \eqn{d_{\text{eff}} = n_H \bar{p} (1 - \bar{p})} from the
#' Fisher information of the logistic regression treatment coefficient.
#'
#' @param n_sg Integer. Subgroup sample size.
#' @param p_event Numeric. Event probability (0--1), averaged across arms.
#'
#' @return Numeric. Effective information.
#'
#' @details
#' For the model \eqn{\text{logit}(p_i) = \beta_0 + \beta V_i} with
#' balanced arms (\eqn{n_1 = n_0 = n_H/2}):
#' \deqn{\text{Var}(\hat\beta) = \frac{4}{n_H \bar{p}(1-\bar{p})}}
#'
#' @examples
#' d_eff_binary(n_sg = 100, p_event = 0.30)
#' # 21
#'
#' @export
d_eff_binary <- function(n_sg, p_event) {
  stopifnot(n_sg > 0, p_event > 0, p_event < 1)
  n_sg * p_event * (1 - p_event)
}


#' Effective Information for Count Outcomes (Poisson)
#'
#' Computes \eqn{d_{\text{eff}} = D / \phi} where \eqn{D} is the total
#' expected events and \eqn{\phi} is the dispersion parameter.
#'
#' @param n_sg Integer. Subgroup sample size.  Required if
#'   \code{total_events} is not provided.
#' @param event_rate Numeric. Mean events per patient.  Required if
#'   \code{total_events} is not provided.
#' @param total_events Numeric. Total observed events \eqn{D}.
#'   Overrides \code{n_sg * event_rate} if provided.
#' @param dispersion Numeric. Overdispersion parameter (1.0 for Poisson,
#'   greater than 1 for quasi-Poisson).  Default: 1.0.
#'
#' @return Numeric. Effective information.
#'
#' @details
#' For Poisson regression with log-offset:
#' \deqn{\text{Var}(\hat\beta) = \phi \left(\frac{1}{D_0} + \frac{1}{D_1}\right)
#'   \approx \frac{4\phi}{D}}
#' Under proportional hazards, \eqn{D = d} (Cox events), so
#' \eqn{d_{\text{eff}}^{\text{Poisson}} = d_{\text{eff}}^{\text{Cox}}}
#' (Laird and Olivier, 1981).
#'
#' @examples
#' # From total events
#' d_eff_count(total_events = 55)
#'
#' # From sample size and rate
#' d_eff_count(n_sg = 100, event_rate = 0.8)
#'
#' # Quasi-Poisson with overdispersion
#' d_eff_count(total_events = 55, dispersion = 1.5)
#'
#' @export
d_eff_count <- function(n_sg = NULL, event_rate = NULL,
                        total_events = NULL, dispersion = 1.0) {
  if (!is.null(total_events)) {
    D <- total_events
  } else if (!is.null(n_sg) && !is.null(event_rate)) {
    D <- n_sg * event_rate
  } else {
    stop("Provide either total_events or both n_sg and event_rate")
  }
  stopifnot(D > 0, dispersion > 0)
  D / dispersion
}


#' Effective Information for Continuous Outcomes (Gaussian)
#'
#' Computes \eqn{d_{\text{eff}} = n_H / \sigma_Y^2} from the
#' Fisher information of the Gaussian linear model treatment coefficient.
#'
#' @param n_sg Integer. Subgroup sample size.
#' @param sigma_y Numeric. Residual standard deviation of the outcome.
#'
#' @return Numeric. Effective information.
#'
#' @details
#' For the model \eqn{Y_i = \beta_0 + \beta V_i + \epsilon_i},
#' \eqn{\epsilon_i \sim N(0, \sigma^2)}, with balanced arms:
#' \deqn{\text{Var}(\hat\beta) = \frac{4\sigma^2}{n_H}}
#'
#' @examples
#' d_eff_continuous(n_sg = 100, sigma_y = 1.5)
#' # 44.4
#'
#' @export
d_eff_continuous <- function(n_sg, sigma_y) {
  stopifnot(n_sg > 0, sigma_y > 0)
  n_sg / sigma_y^2
}


# =============================================================================
# Main Detection Probability Function
# =============================================================================

#' Compute Detection Probability for GLM Outcomes
#'
#' Extends the Section 2.1 approximation (Leon et al., 2024, eq. 3) to
#' arbitrary GLM outcomes via the effective information
#' \eqn{d_{\text{eff}}}.  This is the GLM generalization of
#' \code{\link{compute_detection_probability}}.
#'
#' @param theta Numeric (scalar or vector).  True treatment effect in the
#'   subgroup on the natural scale:
#'   \describe{
#'     \item{Ratio measures (HR, OR, IRR)}{The ratio itself, e.g., 2.0
#'       for doubling.}
#'     \item{Difference measures (MD)}{The mean difference, e.g., 0.5.}
#'   }
#' @param d_eff Numeric.  Effective information in the subgroup.
#'   Use the \code{d_eff_*} helpers:
#'   \code{\link{d_eff_survival}}, \code{\link{d_eff_binary}},
#'   \code{\link{d_eff_count}}, \code{\link{d_eff_continuous}}.
#' @param c1 Numeric. Screening threshold on the natural scale.
#'   For ratio measures: the ratio (e.g., 1.25).
#'   For difference measures: the difference (e.g., 0.2).
#' @param c2 Numeric. Consistency threshold on the natural scale.
#'   Default: 1.0 for ratio, 0.0 for difference.
#' @param effect_scale Character. \code{"ratio"} (default) for log-scale
#'   effects (HR, OR, IRR); \code{"difference"} for identity-scale (MD).
#' @param method Character. \code{"cubature"} or \code{"monte_carlo"}.
#' @param n_mc Integer. Monte Carlo samples if method = "monte_carlo".
#' @param seed Integer. RNG seed for Monte Carlo.
#' @param verbose Logical. Progress messages for vector theta.
#'
#' @return If \code{theta} is scalar, a single probability.
#'   If vector, a data.frame with columns \code{theta}, \code{probability}.
#'
#' @details
#' The detection criterion (Leon et al., 2024, Section 2.1):
#' \enumerate{
#'   \item Screening: \eqn{W_1 + W_2 \ge 2k_1}
#'   \item Consistency: \eqn{\min(W_1, W_2) \ge k_2}
#' }
#' with \eqn{W_1, W_2 \sim N(\beta, 8/d_{\text{eff}})} independently.
#'
#' For ratio-scale effects, \eqn{\beta = \log(\theta)},
#' \eqn{k_j = \log(c_j)}.
#' For difference-scale effects, \eqn{\beta = \theta}, \eqn{k_j = c_j}.
#'
#' @examples
#' \donttest{
#' # Survival: reproduce Figure 2
#' d <- d_eff_survival(n_sg = 60, prop_cens = 0.45)
#' compute_detection_probability_glm(theta = 2.0, d_eff = d)
#'
#' # Binary: OR = 2.0, event rate 30%
#' d <- d_eff_binary(n_sg = 100, p_event = 0.30)
#' compute_detection_probability_glm(theta = 2.0, d_eff = d)
#'
#' # Count: IRR = 2.0, 80 total events
#' d <- d_eff_count(total_events = 80)
#' compute_detection_probability_glm(theta = 2.0, d_eff = d)
#'
#' # Continuous: mean difference = 0.5, sigma = 1.5
#' d <- d_eff_continuous(n_sg = 100, sigma_y = 1.5)
#' compute_detection_probability_glm(
#'   theta = 0.5, d_eff = d,
#'   c1 = 0.2, c2 = 0.0,
#'   effect_scale = "difference"
#' )
#'
#' # Power curve
#' d <- d_eff_binary(100, 0.30)
#' result <- compute_detection_probability_glm(
#'   theta = seq(0.5, 3.0, by = 0.1), d_eff = d
#' )
#' plot(result$theta, result$probability, type = "l")
#' }
#'
#' @seealso \code{\link{compute_detection_probability}} for the
#'   survival-specific version.
#'
#' @export
compute_detection_probability_glm <- function(
    theta,
    d_eff,
    c1 = 1.25,
    c2 = 1.0,
    effect_scale = c("ratio", "difference"),
    method = c("cubature", "monte_carlo"),
    n_mc = 500000L,
    seed = 42L,
    verbose = FALSE
) {

  # --- Argument matching ---

  effect_scale <- match.arg(effect_scale)
  method <- match.arg(method)

  # --- Input validation ---
  if (d_eff <= 0) {
    stop("d_eff must be positive")
  }
  if (effect_scale == "ratio" && any(theta <= 0)) {
    stop("theta must be positive for ratio-scale effects")
  }
  if (effect_scale == "ratio" && (c1 <= 0 || c2 <= 0)) {
    stop("c1 and c2 must be positive for ratio-scale effects")
  }
  if (d_eff < 5) {
    warning("d_eff < 5: normal approximation may be unreliable")
  }

  # --- Transform to natural parameter scale ---
  if (effect_scale == "ratio") {
    beta <- log(theta)
    k1 <- log(c1)
    k2 <- log(c2)
  } else {
    beta <- theta
    k1 <- c1
    k2 <- c2
  }

  # --- Split-half variance ---
  sigma2_s <- 8 / d_eff

  # --- Handle vector theta ---
  if (length(beta) > 1) {
    probs <- vapply(beta, function(b) {
      .detect_prob_glm_single(b, sigma2_s, k1, k2, method, n_mc, seed)
    }, numeric(1))
    return(data.frame(theta = theta, probability = probs))
  }

  # --- Single theta ---
  .detect_prob_glm_single(beta, sigma2_s, k1, k2, method, n_mc, seed)
}


#' Internal Workhorse for Single-Theta Detection Probability
#'
#' @param beta Numeric. True effect on natural parameter scale.
#' @param sigma2_s Numeric. Split-half variance (8 / d_eff).
#' @param k1 Numeric. Screening threshold (natural scale).
#' @param k2 Numeric. Consistency threshold (natural scale).
#' @param method Character. Integration method.
#' @param n_mc Integer. Monte Carlo samples.
#' @param seed Integer. RNG seed.
#'
#' @return Numeric probability.
#'
#' @keywords internal
.detect_prob_glm_single <- function(beta, sigma2_s, k1, k2,
                                    method, n_mc, seed) {
  sd_s <- sqrt(sigma2_s)

  if (method == "cubature") {
    if (!requireNamespace("cubature", quietly = TRUE)) {
      stop("Package 'cubature' required for method = 'cubature'.",
           call. = FALSE)
    }

    integrand <- function(x) {
      meets <- (x[1] + x[2] >= 2 * k1) &&
               (x[1] >= k2) && (x[2] >= k2)
      if (!meets) return(0)
      dnorm(x[1], beta, sd_s) * dnorm(x[2], beta, sd_s)
    }

    # Integration bounds: the indicator restricts to the region
    # where both w1 >= k2 and w2 >= k2.  Use generous upper bounds.
    upper <- beta + 6 * sd_s
    result <- cubature::adaptIntegrate(
      f = integrand,
      lowerLimit = c(k2, k2),
      upperLimit = c(upper, upper),
      tol = 1e-4
    )
    result$integral

  } else {
    # Monte Carlo integration
    set.seed(seed)
    w1 <- rnorm(n_mc, mean = beta, sd = sd_s)
    w2 <- rnorm(n_mc, mean = beta, sd = sd_s)
    mean((w1 + w2 >= 2 * k1) & (w1 >= k2) & (w2 >= k2))
  }
}
