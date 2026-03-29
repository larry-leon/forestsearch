# =============================================================================
# Null FPR Correction: Calibrating L_eff
# =============================================================================
#
# The per-subgroup approximation P1 (Section 2.1) gives the probability
# that a single pre-specified subgroup passes the screening and
# consistency criteria.  The procedure-level FPR is higher because
# ForestSearch evaluates multiple candidate subgroups:
#
#   FPR_corrected = 1 - (1 - P1)^L_eff
#
# where L_eff is the effective number of approximately independent
# candidates.  L_eff depends on the candidate generation process
# (confounders, GRF, LASSO, maxk, n.min) and grows with N.
#
# These functions calibrate L_eff from a small reference simulation.
#
# =============================================================================

#' @importFrom stats lm coef predict
NULL


# =============================================================================
# Layer 1: Fit L_eff from Empirical Data
# =============================================================================

#' Calibrate L_eff from Simulation Results
#'
#' Given per-subgroup approximation probabilities and simulated FPR
#' values at multiple sample sizes, fits the power-law model
#' \eqn{L_{\text{eff}}(N) = C \cdot (N / n_{\min})^\alpha}.
#'
#' @param N Integer vector. Sample sizes at which simulations were run.
#' @param P1 Numeric vector. Per-subgroup detection probability from
#'   \code{\link{compute_detection_probability_glm}} at each N.
#' @param sim_fpr Numeric vector. Simulated procedure-level FPR at each N.
#' @param n_min Integer. Minimum subgroup size used in ForestSearch.
#'   Default: 60.
#'
#' @return A list of class \code{"leff_calibration"} with components:
#'   \describe{
#'     \item{C}{Numeric. Fitted constant.}
#'     \item{alpha}{Numeric. Fitted power-law exponent.}
#'     \item{n_min}{Integer. Minimum subgroup size.}
#'     \item{r_squared}{Numeric. R-squared of the log-linear fit.}
#'     \item{data}{Data frame with N, P1, sim_fpr, L_eff, and fitted values.}
#'     \item{fit}{The \code{lm} object.}
#'   }
#'
#' @details
#' At each sample size, the implied \eqn{L_{\text{eff}}} is:
#' \deqn{L_{\text{eff}} = \frac{\log(1 - \text{FPR}_{\text{sim}})}
#'   {\log(1 - P_1)}}
#'
#' A log-linear model is then fit:
#' \eqn{\log(L_{\text{eff}}) = \log(C) + \alpha \log(N / n_{\min})}.
#'
#' @examples
#' \donttest{
#' # From binary threshold calibration results
#' cal <- calibrate_L_eff(
#'   N       = c(200, 500, 700, 1000),
#'   P1      = c(0.152, 0.109, 0.091, 0.071),
#'   sim_fpr = c(0.185, 0.220, 0.270, 0.280),
#'   n_min   = 60
#' )
#' print(cal)
#' }
#'
#' @export
calibrate_L_eff <- function(N, P1, sim_fpr, n_min = 60) {

  # --- Validation ---
  stopifnot(length(N) >= 2, length(N) == length(P1),
            length(N) == length(sim_fpr))
  stopifnot(all(P1 > 0), all(P1 < 1))
  stopifnot(all(sim_fpr > 0), all(sim_fpr < 1))
  stopifnot(n_min > 0)

  # --- Implied L_eff ---
  L_eff <- log(1 - sim_fpr) / log(1 - P1)
  ratio <- N / n_min

  # --- Fit power law ---
  fit <- lm(log(L_eff) ~ log(ratio))
  C <- exp(coef(fit)[1])
  alpha <- unname(coef(fit)[2])
  r2 <- summary(fit)$r.squared

  # --- Build result ---
  cal_data <- data.frame(
    N = N, P1 = P1, sim_fpr = sim_fpr,
    L_eff = L_eff,
    L_eff_fitted = C * ratio^alpha,
    fpr_corrected = 1 - (1 - P1)^(C * ratio^alpha)
  )

  result <- list(
    C = unname(C),
    alpha = alpha,
    n_min = n_min,
    r_squared = r2,
    data = cal_data,
    fit = fit
  )
  class(result) <- "leff_calibration"
  result
}


#' @export
print.leff_calibration <- function(x, ...) {
  cat("L_eff Calibration\n")
  cat(sprintf("  Model: L_eff = %.3f * (N / %d)^%.3f\n",
              x$C, x$n_min, x$alpha))
  cat(sprintf("  R-squared: %.3f\n", x$r_squared))
  cat(sprintf("  Calibrated at N = %s\n",
              paste(x$data$N, collapse = ", ")))
  cat("\n  Implied L_eff:\n")
  for (i in seq_len(nrow(x$data))) {
    cat(sprintf("    N = %4d: L_eff = %.2f (fitted: %.2f)\n",
                x$data$N[i], x$data$L_eff[i], x$data$L_eff_fitted[i]))
  }
  invisible(x)
}


# =============================================================================
# Layer 2: Predict Corrected FPR
# =============================================================================

#' Predict Corrected Null FPR
#'
#' Uses a calibrated \eqn{L_{\text{eff}}} model to predict the
#' procedure-level FPR at arbitrary sample sizes and thresholds.
#'
#' @param calibration An \code{"leff_calibration"} object from
#'   \code{\link{calibrate_L_eff}}.
#' @param theta Numeric. True treatment effect (1.0 for null).
#' @param d_eff Numeric. Effective information in the subgroup.
#' @param c1 Numeric. Screening threshold.
#' @param c2 Numeric. Consistency threshold.
#' @param N Integer. Total sample size.
#' @param effect_scale Character. \code{"ratio"} or \code{"difference"}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{P1}{Per-subgroup detection probability.}
#'     \item{L_eff}{Effective number of candidates at this N.}
#'     \item{fpr_corrected}{Corrected procedure-level FPR.}
#'   }
#'
#' @examples
#' \donttest{
#' cal <- calibrate_L_eff(
#'   N = c(200, 500, 700, 1000),
#'   P1 = c(0.152, 0.109, 0.091, 0.071),
#'   sim_fpr = c(0.185, 0.220, 0.270, 0.280)
#' )
#'
#' # Predict FPR at N = 800
#' d <- d_eff_binary(n_sg = round(0.30 * 800), p_event = 0.30)
#' predict_fpr_corrected(cal, theta = 1.0, d_eff = d,
#'                       c1 = 1.25, c2 = 1.25, N = 800)
#' }
#'
#' @export
predict_fpr_corrected <- function(calibration, theta, d_eff, c1, c2, N,
                                   effect_scale = "ratio") {

  stopifnot(inherits(calibration, "leff_calibration"))
  stopifnot(d_eff > 0, N > 0)

  # Per-subgroup probability
  P1 <- compute_detection_probability_glm(
    theta = theta, d_eff = d_eff, c1 = c1, c2 = c2,
    effect_scale = effect_scale, method = "monte_carlo"
  )

  # Effective candidates
  L_eff <- calibration$C * (N / calibration$n_min)^calibration$alpha

  # Corrected FPR
  fpr <- 1 - (1 - P1)^L_eff

  list(P1 = P1, L_eff = L_eff, fpr_corrected = fpr)
}


# =============================================================================
# Layer 3: Full Calibration Wrapper
# =============================================================================

#' Run Null Calibration Simulation
#'
#' Runs ForestSearch under the null hypothesis at multiple sample
#' sizes to calibrate \eqn{L_{\text{eff}}}, then returns a
#' calibration object for predicting corrected FPR.
#'
#' @param sim_null_fn Function. Generates one H0 dataset.
#'   Signature: \code{sim_null_fn(n, seed)} returning a data.frame
#'   suitable for \code{forestsearch()}.
#' @param fs_args Named list. Arguments passed to \code{forestsearch()}
#'   (excluding \code{df.analysis}, \code{hr.threshold},
#'   \code{hr.consistency}, and \code{seedit}, which are set by the
#'   wrapper).
#' @param d_eff_fn Function. Computes \eqn{d_{\text{eff}}} for a
#'   subgroup.  Signature: \code{d_eff_fn(n_sg)} returning numeric.
#'   Example: \code{function(n) d_eff_binary(n, 0.30)}.
#' @param N_values Integer vector. Sample sizes to simulate.
#'   At least 2 required.  Default: \code{c(200, 500, 1000)}.
#' @param c1 Numeric. Screening threshold for calibration.
#'   Default: 1.25.
#' @param c2 Numeric. Consistency threshold for calibration.
#'   Default: 1.0.
#' @param n_sims Integer. Simulations per sample size.
#'   Default: 100.
#' @param prop_harm Numeric. Expected proportion in harm subgroup
#'   (for computing d_eff). Default: 0.30.
#' @param n_min Integer. Minimum subgroup size. Default: 60.
#' @param seed_base Integer. Base seed. Default: 42.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An \code{"leff_calibration"} object.
#'
#' @details
#' For each N in \code{N_values}, the function:
#' \enumerate{
#'   \item Generates \code{n_sims} datasets via \code{sim_null_fn()}
#'   \item Runs \code{forestsearch()} on each
#'   \item Computes the simulated FPR = proportion with any H found
#'   \item Computes \eqn{P_1} via the GLM approximation
#'   \item Calls \code{calibrate_L_eff()} to fit the power-law model
#' }
#'
#' ForestSearch is run with
#' \code{parallel_args = list(plan = "sequential")} to avoid
#' nested parallelism.  The outer simulation loop is sequential.
#' For faster execution, wrap in \code{foreach \%dofuture\%}
#' externally.
#'
#' @examples
#' \donttest{
#' # Binary DGM under H0
#' sim_null <- function(n, seed) {
#'   set.seed(seed)
#'   data.frame(
#'     id = seq_len(n),
#'     treat = rbinom(n, 1, 0.5),
#'     bm1 = as.factor(rbinom(n, 1, 0.70)),
#'     bm2 = as.factor(rbinom(n, 1, 0.50)),
#'     age = round(rnorm(n, 55, 10)),
#'     ecog = as.factor(sample(0:1, n, TRUE, c(0.6, 0.4))),
#'     progressed = rbinom(n, 1, 0.30)
#'   )
#' }
#'
#' fs_args <- list(
#'   confounders.name = c("bm1", "bm2", "age", "ecog"),
#'   outcome.name = "progressed",
#'   event.name = "progressed",
#'   treat.name = "treat",
#'   id.name = "id",
#'   outcome_type = "binary",
#'   effect_measure = "OR",
#'   adverse_outcome = TRUE,
#'   pconsistency.threshold = 0.90,
#'   fs.splits = 200,
#'   n.min = 60,
#'   maxk = 2,
#'   use_lasso = TRUE,
#'   use_grf = TRUE,
#'   use_twostage = TRUE,
#'   is.RCT = TRUE,
#'   details = FALSE,
#'   plot.sg = FALSE,
#'   parallel_args = list(plan = "sequential",
#'                        workers = 1,
#'                        show_message = FALSE)
#' )
#'
#' cal <- run_null_calibration(
#'   sim_null_fn = sim_null,
#'   fs_args = fs_args,
#'   d_eff_fn = function(n) d_eff_binary(n, 0.30),
#'   N_values = c(200, 500, 1000),
#'   n_sims = 50,
#'   verbose = TRUE
#' )
#' print(cal)
#' }
#'
#' @export
run_null_calibration <- function(
    sim_null_fn,
    fs_args,
    d_eff_fn,
    N_values = c(200, 500, 1000),
    c1 = 1.25,
    c2 = 1.0,
    n_sims = 100L,
    prop_harm = 0.30,
    n_min = 60L,
    seed_base = 42L,
    verbose = TRUE
) {

  # --- Validation ---
  stopifnot(is.function(sim_null_fn))
  stopifnot(is.list(fs_args))
  stopifnot(is.function(d_eff_fn))
  stopifnot(length(N_values) >= 2)
  stopifnot(n_sims >= 10)

  # --- Ensure sequential inner plan ---
  fs_args$parallel_args <- list(
    plan = "sequential", workers = 1, show_message = FALSE
  )

  if (verbose) {
    cat(sprintf("Null calibration: %d N values x %d sims = %d tasks\n",
                length(N_values), n_sims, length(N_values) * n_sims))
  }

  # --- Run simulations at each N ---
  results <- data.frame(
    N = integer(0), fpr = numeric(0), P1 = numeric(0)
  )

  for (N_val in N_values) {
    t0 <- proc.time()[3]
    n_found <- 0L

    for (s in seq_len(n_sims)) {
      seed_s <- seed_base + (N_val * 1000L) + s
      df_null <- sim_null_fn(N_val, seed_s)

      fs_call_args <- c(
        list(df.analysis = df_null),
        fs_args,
        list(
          hr.threshold = c1,
          hr.consistency = c2,
          seedit = seed_s
        )
      )

      fs_result <- tryCatch(
        do.call(forestsearch, fs_call_args),
        error = function(e) NULL
      )

      found <- !is.null(fs_result) &&
               !is.null(fs_result$sg.harm) &&
               length(fs_result$sg.harm) > 0
      if (found) n_found <- n_found + 1L
    }

    sim_fpr <- n_found / n_sims

    # Per-subgroup approximation
    n_sg <- round(prop_harm * N_val)
    d <- d_eff_fn(n_sg)
    P1_val <- compute_detection_probability_glm(
      theta = 1.0, d_eff = d, c1 = c1, c2 = c2,
      effect_scale = "ratio", method = "monte_carlo"
    )

    elapsed <- proc.time()[3] - t0

    if (verbose) {
      cat(sprintf("  N = %4d: FPR = %.3f (%d/%d) | P1 = %.3f | %.1f sec\n",
                  N_val, sim_fpr, n_found, n_sims, P1_val, elapsed))
    }

    results <- rbind(results, data.frame(
      N = N_val, fpr = sim_fpr, P1 = P1_val
    ))
  }

  # --- Handle edge cases ---
  # Clamp FPR away from 0 and 1 for log transform
  results$fpr <- pmax(results$fpr, 0.5 / n_sims)
  results$fpr <- pmin(results$fpr, 1 - 0.5 / n_sims)

  # --- Fit L_eff ---
  cal <- calibrate_L_eff(
    N = results$N,
    P1 = results$P1,
    sim_fpr = results$fpr,
    n_min = n_min
  )

  if (verbose) {
    cat("\n")
    print(cal)
  }

  cal
}
