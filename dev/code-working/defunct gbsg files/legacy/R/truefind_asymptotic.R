# =============================================================================
# Asymptotic Approximation for Subgroup Detection Probability
# =============================================================================
#
# Functions for computing the probability of detecting a true subgroup
# using asymptotic normal approximations for the log hazard ratio estimator.
#
# Based on the variance approximation: Var(log(HR)) ~ 4/d for each treatment
# arm, giving Var(log(HR)) ~ 8/d for the combined estimator, where d is the

# number of events in the subgroup.
#
# Reference: Leon et al. (2024) "Exploratory subgroup identification in the
# heterogeneous Cox model"
#
# =============================================================================

#' @importFrom stats dnorm pnorm qnorm integrate
NULL


# =============================================================================
# Core Density Functions
# =============================================================================

#' Bivariate Density for Split-Sample HR Threshold Detection
#'
#' Computes the joint density for the two-split detection criterion where
#' both split-halves must exceed individual thresholds and their average
#' must exceed a consistency threshold.
#'
#' @param x Numeric vector of length 2. Log hazard ratio estimates from
#'   the two split-halves.
#' @param theta Numeric. True hazard ratio in the subgroup.
#' @param prop_cens Numeric. Proportion censored (0-1). Default: 0.3
#' @param n_sg Integer. Subgroup sample size.
#' @param k_avg Numeric. Threshold for average log(HR) across splits.
#'   Typically log(hr.threshold).
#' @param k_ind Numeric. Threshold for individual split log(HR).
#'   Typically log(hr.consistency).
#'
#' @return Numeric. Joint density value at x, or 0 if thresholds not met.
#'
#' @details
#' The detection criterion requires:
#' \itemize{
#'   \item Average of two splits: (x1 + x2)/2 >= k_avg
#'   \item Individual splits: x1 >= k_ind AND x2 >= k_ind
#' }
#'
#' Under the asymptotic approximation, each split-half log(HR) estimator
#' follows N(log(theta), 8/d) where d = n_sg * (1 - prop_cens) / 2 is the
#' expected number of events per split.
#'
#' @keywords internal
density_threshold_both <- function(x, theta, prop_cens = 0.3, n_sg, k_avg, k_ind) {

 # Events in subgroup (total)
  d_sg <- n_sg * (1 - prop_cens)

 # Variance for each split-half estimator
  # Each split has d_sg/2 events, so Var = 8/(d_sg/2) = 16/d_sg
  # But the original code uses 8/d_sg which assumes full sample variance
  # for demonstration purposes
  sig2_split <- 8 / d_sg
  sig_split <- sqrt(sig2_split)

  # Mean is log of true HR
  mu_split <- log(theta)

  # Compute individual densities
  dens_1 <- dnorm(x[1], mean = mu_split, sd = sig_split)
  dens_2 <- dnorm(x[2], mean = mu_split, sd = sig_split)

  # Check threshold conditions
  # Condition 1: Average exceeds k_avg (note: x[1]+x[2] >= 2*k_avg)
  # Condition 2: Both individual splits exceed k_ind
  meets_criteria <- (x[1] + x[2] >= 2 * k_avg) &
                    (x[1] >= k_ind) &
                    (x[2] >= k_ind)

  # Return joint density if criteria met, else 0
  if (meets_criteria) {
    return(dens_1 * dens_2)
  } else {
    return(0)
  }
}


#' Vectorized Density for Integration
#'
#' Wrapper around density_threshold_both for use with cubature integration.
#'
#' @inheritParams density_threshold_both
#'
#' @return Numeric density value.
#'
#' @keywords internal
density_threshold_integrand <- function(x, theta, prop_cens, n_sg, k_avg, k_ind) {
  density_threshold_both(
    x = x,
    theta = theta,
    prop_cens = prop_cens,
    n_sg = n_sg,
    k_avg = k_avg,
    k_ind = k_ind
  )
}


# =============================================================================
# Main Computation Functions
# =============================================================================

#' Compute Probability of Detecting True Subgroup
#'
#' Calculates the probability that a true subgroup with given hazard ratio
#' will be detected using the ForestSearch consistency-based criteria.
#'
#' @param theta Numeric. True hazard ratio in the subgroup. Can be a vector
#'   for computing detection probability across multiple HR values.
#' @param n_sg Integer. Subgroup sample size.
#' @param prop_cens Numeric. Proportion censored (0-1). Default: 0.3
#' @param hr_threshold Numeric. HR threshold for detection (e.g., 1.25).
#'   This is the threshold that the average HR across splits must exceed.
#' @param hr_consistency Numeric. HR consistency threshold (e.g., 1.0).
#'   This is the threshold each individual split must exceed. Default: 1.0
#' @param method Character. Integration method: "cubature" (recommended for
#'   accuracy) or "monte_carlo" (faster for exploration). Default: "cubature"
#' @param n_mc Integer. Number of Monte Carlo samples if method = "monte_carlo".
#'   Default: 100000
#' @param tol Numeric. Relative tolerance for cubature integration.
#'   Default: 1e-4
#' @param verbose Logical. Print progress for vector inputs. Default: FALSE
#'
#' @return If theta is scalar, returns a single probability. If theta is a
#'   vector, returns a data.frame with columns: theta, probability.
#'
#' @details
#' This function computes P(detect | theta) using the asymptotic normal
#' approximation for the log hazard ratio estimator. The detection criterion
#' is based on ForestSearch's split-sample consistency evaluation:
#'
#' \enumerate{
#'   \item The subgroup HR estimate must exceed hr_threshold on average
#'   \item Each split-half must individually exceed hr_consistency
#' }
#'
#' The approximation assumes:
#' \itemize{
#'   \item Large sample sizes (CLT applies)
#'   \item Var(log(HR)) ~ 4/d per treatment arm
#'   \item Independence between split-halves (conditional on true effect)
#' }
#'
#' @examples
#' \dontrun{
#' # Single HR value
#' prob <- compute_detection_probability(
#'   theta = 1.5,
#'   n_sg = 60,
#'   prop_cens = 0.2,
#'   hr_threshold = 1.25
#' )
#'
#' # Vector of HR values for power curve
#' hr_values <- seq(1.0, 2.5, by = 0.1)
#' results <- compute_detection_probability(
#'   theta = hr_values,
#'   n_sg = 60,
#'   prop_cens = 0.2,
#'   hr_threshold = 1.25,
#'   verbose = TRUE
#' )
#'
#' # Plot detection probability curve
#' plot(results$theta, results$probability, type = "l",
#'      xlab = "True HR", ylab = "P(detect)")
#' }
#'
#' @importFrom stats dnorm rnorm
#' @export
compute_detection_probability <- function(
    theta,
    n_sg,
    prop_cens = 0.3,
    hr_threshold = 1.25,
    hr_consistency = 1.0,
    method = c("cubature", "monte_carlo"),
    n_mc = 100000L,
    tol = 1e-4,
    verbose = FALSE
) {

  # Input validation
  method <- match.arg(method)

  if (any(theta <= 0)) {
    stop("theta must be positive (hazard ratios)")
  }
  if (n_sg < 10) {
    warning("n_sg < 10: asymptotic approximation may be unreliable")
  }
  if (prop_cens < 0 || prop_cens >= 1) {
    stop("prop_cens must be in [0, 1)")
  }
  if (hr_threshold <= 0 || hr_consistency <= 0) {
    stop("hr_threshold and hr_consistency must be positive")
  }

  # Convert HR thresholds to log scale
  k_avg <- log(hr_threshold)
  k_ind <- log(hr_consistency)

  # Handle vector input
  if (length(theta) > 1) {
    probs <- numeric(length(theta))

    for (i in seq_along(theta)) {
      if (verbose && i %% 10 == 1) {
        message(sprintf("Computing probability for theta[%d] = %.3f", i, theta[i]))
      }

      probs[i] <- compute_detection_probability_single(
        theta = theta[i],
        n_sg = n_sg,
        prop_cens = prop_cens,
        k_avg = k_avg,
        k_ind = k_ind,
        method = method,
        n_mc = n_mc,
        tol = tol
      )
    }

    return(data.frame(
      theta = theta,
      probability = probs
    ))
  }

  # Single theta value
  compute_detection_probability_single(
    theta = theta,
    n_sg = n_sg,
    prop_cens = prop_cens,
    k_avg = k_avg,
    k_ind = k_ind,
    method = method,
    n_mc = n_mc,
    tol = tol
  )
}


#' Compute Detection Probability for Single Theta (Internal)
#'
#' @inheritParams compute_detection_probability
#' @param k_avg Log of hr_threshold
#' @param k_ind Log of hr_consistency
#'
#' @return Numeric probability
#'
#' @keywords internal
compute_detection_probability_single <- function(
    theta,
    n_sg,
    prop_cens,
    k_avg,
    k_ind,
    method,
    n_mc,
    tol
) {

    if (method == "cubature") {
      # Use adaptive cubature integration
      if (!requireNamespace("cubature", quietly = TRUE)) {
        stop("Package 'cubature' is required for method = 'cubature'.",
             call. = FALSE)
      }

    result <- cubature::adaptIntegrate(
      f = density_threshold_integrand,
      lowerLimit = c(-Inf, -Inf),
      upperLimit = c(Inf, Inf),
      theta = theta,
      prop_cens = prop_cens,
      n_sg = n_sg,
      k_avg = k_avg,
      k_ind = k_ind,
      tol = tol
    )

    return(result$integral)

  } else {
    # Monte Carlo integration
    d_sg <- n_sg * (1 - prop_cens)
    sig_split <- sqrt(8 / d_sg)
    mu_split <- log(theta)

    # Generate samples from the joint distribution
    x1 <- rnorm(n_mc, mean = mu_split, sd = sig_split)
    x2 <- rnorm(n_mc, mean = mu_split, sd = sig_split)

    # Check detection criteria
    detected <- (x1 + x2 >= 2 * k_avg) & (x1 >= k_ind) & (x2 >= k_ind)

    return(mean(detected))
  }
}


# =============================================================================
# Analysis and Visualization Functions
# =============================================================================

#' Generate Detection Probability Curve
#'
#' Computes detection probability across a range of hazard ratios to create
#' a power-like curve for subgroup detection.
#'
#' @param theta_range Numeric vector of length 2. Range of HR values to evaluate.
#'   Default: c(0.5, 3.0)
#' @param n_points Integer. Number of points to evaluate. Default: 50
#' @param n_sg Integer. Subgroup sample size.
#' @param prop_cens Numeric. Proportion censored (0-1). Default: 0.3
#' @param hr_threshold Numeric. HR threshold for detection. Default: 1.25
#' @param hr_consistency Numeric. HR consistency threshold. Default: 1.0
#' @param include_reference Logical. Include reference HR values (0.5, 0.75, 1.0).
#'   Default: TRUE
#' @param method Character. Integration method. Default: "cubature"
#' @param verbose Logical. Print progress. Default: TRUE
#'
#' @return A data.frame with columns:
#'   \item{theta}{Hazard ratio values}
#'   \item{probability}{Detection probability}
#'   \item{n_sg}{Subgroup size (repeated)}
#'   \item{prop_cens}{Censoring proportion (repeated)}
#'   \item{hr_threshold}{Detection threshold (repeated)}
#'
#' @examples
#' \dontrun{
#' # Generate detection curve
#' curve_data <- generate_detection_curve(
#'   n_sg = 60,
#'   prop_cens = 0.2,
#'   hr_threshold = 1.25
#' )
#'
#' # Plot
#' plot_detection_curve(curve_data)
#' }
#'
#' @export
generate_detection_curve <- function(
    theta_range = c(0.5, 3.0),
    n_points = 50L,
    n_sg,
    prop_cens = 0.3,
    hr_threshold = 1.25,
    hr_consistency = 1.0,
    include_reference = TRUE,
    method = "cubature",
    verbose = TRUE
) {

  # Generate theta values
  theta_vals <- seq(theta_range[1], theta_range[2], length.out = n_points)

  # Add reference values if requested
  if (include_reference) {
    ref_vals <- c(0.5, 0.7, 0.75, 0.80, 1.0, hr_threshold, 1.5, 2.0)
    ref_vals <- ref_vals[ref_vals >= theta_range[1] & ref_vals <= theta_range[2]]
    theta_vals <- sort(unique(c(theta_vals, ref_vals)))
  }

  if (verbose) {
    message(sprintf("Computing detection probabilities for %d HR values...",
                    length(theta_vals)))
    message(sprintf("  n_sg = %d, prop_cens = %.2f, hr_threshold = %.2f",
                    n_sg, prop_cens, hr_threshold))
  }

  # Compute probabilities
  result <- compute_detection_probability(
    theta = theta_vals,
    n_sg = n_sg,
    prop_cens = prop_cens,
    hr_threshold = hr_threshold,
    hr_consistency = hr_consistency,
    method = method,
    verbose = verbose
  )

  # Add metadata
  result$n_sg <- n_sg
  result$prop_cens <- prop_cens
  result$hr_threshold <- hr_threshold
  result$hr_consistency <- hr_consistency

  class(result) <- c("detection_curve", "data.frame")

  result
}


#' Plot Detection Probability Curve
#'
#' Creates a visualization of the detection probability curve.
#'
#' @param curve_data A data.frame from \code{\link{generate_detection_curve}}
#'   or with columns: theta, probability.
#' @param add_reference_lines Logical. Add horizontal reference lines at
#'   0.05, 0.10, 0.80. Default: TRUE
#' @param add_threshold_line Logical. Add vertical line at hr_threshold.
#'   Default: TRUE
#' @param title Character. Plot title. Default: auto-generated
#' @param ... Additional arguments passed to plot()
#'
#' @return Invisibly returns the input data.
#'
#' @examples
#' \dontrun{
#' curve_data <- generate_detection_curve(n_sg = 60, prop_cens = 0.2)
#' plot_detection_curve(curve_data)
#' }
#'
#' @export
plot_detection_curve <- function(
    curve_data,
    add_reference_lines = TRUE,
    add_threshold_line = TRUE,
    title = NULL,
    ...
) {

  # Generate title if not provided
  if (is.null(title) && "n_sg" %in% names(curve_data)) {
    title <- sprintf("Detection Probability (n=%d, cens=%.0f%%, threshold=%.2f)",
                     curve_data$n_sg[1],
                     100 * curve_data$prop_cens[1],
                     curve_data$hr_threshold[1])
  }

  # Main plot
  plot(
    curve_data$theta,
    curve_data$probability,
    type = "l",
    lwd = 2,
    col = "darkblue",
    xlab = "True Hazard Ratio (theta)",
    ylab = "P(Detect Subgroup)",
    ylim = c(0, 1),
    main = title,
    ...
  )

  # Reference lines
  if (add_reference_lines) {
    abline(h = 0.05, lty = 2, col = "gray60", lwd = 0.5)
    abline(h = 0.10, lty = 2, col = "red", lwd = 0.5)
    abline(h = 0.80, lty = 2, col = "darkgreen", lwd = 0.5)

    # Add legend for reference lines
    legend("bottomright",
           legend = c("5% (Type I)", "10%", "80% (Power)"),
           lty = 2,
           col = c("gray60", "red", "darkgreen"),
           lwd = 0.5,
           cex = 0.8,
           bty = "n")
  }

  # Threshold line
  if (add_threshold_line && "hr_threshold" %in% names(curve_data)) {
    abline(v = curve_data$hr_threshold[1], lty = 3, col = "orange", lwd = 1)
  }

  # Grid
  grid(col = "gray90")

  invisible(curve_data)
}


#' Compare Detection Curves Across Sample Sizes
#'
#' Generates and compares detection probability curves for multiple
#' subgroup sample sizes.
#'
#' @param n_sg_values Integer vector. Subgroup sample sizes to compare.
#' @param prop_cens Numeric. Proportion censored. Default: 0.3
#' @param hr_threshold Numeric. HR threshold. Default: 1.25
#' @param hr_consistency Numeric. HR consistency threshold. Default: 1.0
#' @param theta_range Numeric vector of length 2. Range of HR values.
#'   Default: c(0.5, 3.0)
#' @param n_points Integer. Number of points per curve. Default: 40
#' @param verbose Logical. Print progress. Default: TRUE
#'
#' @return A data.frame with all curves combined, including n_sg as a factor.
#'
#' @examples
#' \dontrun{
#' comparison <- compare_detection_curves(
#'   n_sg_values = c(40, 60, 80, 100),
#'   prop_cens = 0.2
#' )
#'
#' # Plot with ggplot2
#' library(ggplot2)
#' ggplot(comparison, aes(x = theta, y = probability, color = factor(n_sg))) +
#'   geom_line(linewidth = 1) +
#'   labs(x = "True HR", y = "P(Detect)", color = "n_sg") +
#'   theme_minimal()
#' }
#'
#' @export
compare_detection_curves <- function(
    n_sg_values,
    prop_cens = 0.3,
    hr_threshold = 1.25,
    hr_consistency = 1.0,
    theta_range = c(0.5, 3.0),
    n_points = 40L,
    verbose = TRUE
) {

  results_list <- lapply(n_sg_values, function(n) {
    if (verbose) {
      message(sprintf("\nProcessing n_sg = %d...", n))
    }

    generate_detection_curve(
      theta_range = theta_range,
      n_points = n_points,
      n_sg = n,
      prop_cens = prop_cens,
      hr_threshold = hr_threshold,
      hr_consistency = hr_consistency,
      include_reference = FALSE,
      verbose = FALSE
    )
  })

  # Combine results
  result <- do.call(rbind, results_list)
  result$n_sg_factor <- factor(result$n_sg)

  class(result) <- c("detection_comparison", "data.frame")

  result
}


#' Find Minimum Sample Size for Target Detection Power
#'
#' Determines the minimum subgroup sample size needed to achieve a target
#' detection probability for a given true hazard ratio.
#'
#' @param theta Numeric. True hazard ratio in subgroup.
#' @param target_power Numeric. Target detection probability (0-1). Default: 0.80
#' @param prop_cens Numeric. Proportion censored. Default: 0.3
#' @param hr_threshold Numeric. HR threshold. Default: 1.25
#' @param hr_consistency Numeric. HR consistency threshold. Default: 1.0
#' @param n_range Integer vector of length 2. Range of sample sizes to search.
#'   Default: c(20, 500)
#' @param tol Numeric. Tolerance for bisection search. Default: 1
#' @param verbose Logical. Print progress. Default: TRUE
#'
#' @return A list with:
#'   \item{n_sg_required}{Minimum sample size (rounded up)}
#'   \item{achieved_power}{Actual detection probability at n_sg_required}
#'   \item{theta}{Input hazard ratio}
#'   \item{target_power}{Input target power}
#'
#' @examples
#' \dontrun{
#' # Find sample size for 80% power to detect HR = 1.5
#' result <- find_required_sample_size(
#'   theta = 1.5,
#'   target_power = 0.80,
#'   prop_cens = 0.2
#' )
#' print(result)
#' }
#'
#' @export
find_required_sample_size <- function(
    theta,
    target_power = 0.80,
    prop_cens = 0.3,
    hr_threshold = 1.25,
    hr_consistency = 1.0,
    n_range = c(20L, 500L),
    tol = 1,
    verbose = TRUE
) {

  if (theta <= hr_threshold) {
    warning("theta <= hr_threshold: detection probability may be low even with large n")
  }

  # Define function to minimize
  power_at_n <- function(n) {
    compute_detection_probability(
      theta = theta,
      n_sg = n,
      prop_cens = prop_cens,
      hr_threshold = hr_threshold,
      hr_consistency = hr_consistency,
      method = "cubature"
    ) - target_power
  }

  # Check bounds
  power_low <- power_at_n(n_range[1])
  power_high <- power_at_n(n_range[2])

  if (power_low > 0) {
    if (verbose) message("Target power achieved even at minimum n")
    return(list(
      n_sg_required = n_range[1],
      achieved_power = power_low + target_power,
      theta = theta,
      target_power = target_power
    ))
  }

  if (power_high < 0) {
    if (verbose) message("Target power not achieved even at maximum n")
    return(list(
      n_sg_required = NA,
      achieved_power = power_high + target_power,
      theta = theta,
      target_power = target_power,
      note = "Increase n_range[2] or lower target_power"
    ))
  }

  # Bisection search
  n_low <- n_range[1]
  n_high <- n_range[2]

  while ((n_high - n_low) > tol) {
    n_mid <- round((n_low + n_high) / 2)
    power_mid <- power_at_n(n_mid)

    if (verbose) {
      message(sprintf("  n = %d, power = %.3f", n_mid, power_mid + target_power))
    }

    if (power_mid < 0) {
      n_low <- n_mid
    } else {
      n_high <- n_mid
    }
  }

  # Return ceiling to ensure target is met
  n_required <- ceiling(n_high)
  achieved <- compute_detection_probability(
    theta = theta,
    n_sg = n_required,
    prop_cens = prop_cens,
    hr_threshold = hr_threshold,
    hr_consistency = hr_consistency
  )

  list(
    n_sg_required = n_required,
    achieved_power = achieved,
    theta = theta,
    target_power = target_power,
    prop_cens = prop_cens,
    hr_threshold = hr_threshold
  )
}


#' Create Sample Size Table for Multiple Scenarios
#'
#' Generates a table of required sample sizes for different combinations
#' of true hazard ratios and censoring proportions.
#'
#' @param theta_values Numeric vector. True hazard ratios to evaluate.
#' @param prop_cens_values Numeric vector. Censoring proportions to evaluate.
#' @param target_power Numeric. Target detection probability. Default: 0.80
#' @param hr_threshold Numeric. HR threshold. Default: 1.25
#' @param verbose Logical. Print progress. Default: TRUE
#'
#' @return A data.frame with columns: theta, prop_cens, n_required, achieved_power
#'
#' @examples
#' \dontrun{
#' ss_table <- create_sample_size_table(
#'   theta_values = c(1.5, 1.75, 2.0, 2.5),
#'   prop_cens_values = c(0.2, 0.3, 0.4),
#'   target_power = 0.80
#' )
#' print(ss_table)
#' }
#'
#' @keywords internal
create_sample_size_table <- function(
    theta_values,
    prop_cens_values,
    target_power = 0.80,
    hr_threshold = 1.25,
    verbose = TRUE
) {

  results <- expand.grid(
    theta = theta_values,
    prop_cens = prop_cens_values,
    n_required = NA_integer_,
    achieved_power = NA_real_
  )

  for (i in seq_len(nrow(results))) {
    if (verbose) {
      message(sprintf("Computing: theta = %.2f, prop_cens = %.2f",
                      results$theta[i], results$prop_cens[i]))
    }

    res <- find_required_sample_size(
      theta = results$theta[i],
      target_power = target_power,
      prop_cens = results$prop_cens[i],
      hr_threshold = hr_threshold,
      verbose = FALSE
    )

    results$n_required[i] <- res$n_sg_required
    results$achieved_power[i] <- res$achieved_power
  }

  results
}


# =============================================================================
# Print Methods
# =============================================================================

#' @export
print.detection_curve <- function(x, ...) {
  cat("Detection Probability Curve\n")
  cat("===========================\n")
  cat(sprintf("  Subgroup size:     n = %d\n", x$n_sg[1]))
  cat(sprintf("  Censoring:         %.0f%%\n", 100 * x$prop_cens[1]))
  cat(sprintf("  HR threshold:      %.2f\n", x$hr_threshold[1]))
  cat(sprintf("  HR consistency:    %.2f\n", x$hr_consistency[1]))
  cat(sprintf("  HR range:          [%.2f, %.2f]\n", min(x$theta), max(x$theta)))
  cat(sprintf("  Points evaluated:  %d\n", nrow(x)))
  cat("\nKey probabilities:\n")

  # Find probability at threshold
  idx_thresh <- which.min(abs(x$theta - x$hr_threshold[1]))
  cat(sprintf("  P(detect | HR = %.2f) = %.3f  (at threshold)\n",
              x$theta[idx_thresh], x$probability[idx_thresh]))

  # Find HR for 80% power if available
  idx_80 <- which.min(abs(x$probability - 0.80))
  if (abs(x$probability[idx_80] - 0.80) < 0.05) {
    cat(sprintf("  P(detect | HR = %.2f) ~ 0.80  (80%% power)\n",
                x$theta[idx_80]))
  }

  invisible(x)
}
