# =============================================================================
# Threshold Suggestion Helper
# =============================================================================
#
# Provides recommended (c1, c2) thresholds for ForestSearch based on
# the corrected approximation.  Does NOT change any defaults.
#
# The paper defaults (c1=1.25, c2=1.0) remain the function signature
# defaults for backward compatibility.  This helper provides guidance
# for non-survival outcomes where different thresholds may be needed.
#
# =============================================================================


#' Suggest Screening and Consistency Thresholds
#'
#' Recommends \code{(c1, c2)} threshold pairs for ForestSearch based on
#' the corrected null approximation.  Uses the per-subgroup detection
#' probability and an \eqn{L_{\text{eff}}} multiplicity correction to
#' find the smallest \code{c1} (screening threshold) such that the
#' corrected FPR is at most \code{fpr_target}.
#'
#' @param d_eff Numeric. Effective information in the expected harm
#'   subgroup.  Use \code{\link{d_eff_survival}},
#'   \code{\link{d_eff_binary}}, \code{\link{d_eff_count}}, or
#'   \code{\link{d_eff_continuous}}.
#' @param N Integer. Total sample size.
#' @param fpr_target Numeric. Maximum acceptable procedure-level FPR.
#'   Default: 0.10.
#' @param n_min Integer. Minimum subgroup size. Default: 60.
#' @param L_eff_C Numeric. Calibration constant for
#'   \eqn{L_{\text{eff}} = C (N/n_{\min})^\alpha}.  Default: 1.0
#'   (conservative; see Details).
#' @param L_eff_alpha Numeric. Calibration exponent. Default: 0.0
#'   (no N-dependence; see Details).
#' @param c1_range Numeric vector of length 2. Search range for c1.
#'   Default: \code{c(1.05, 2.50)}.
#' @param c1_step Numeric. Grid step for c1. Default: 0.05.
#' @param c2_candidates Numeric vector. Candidate c2 values to evaluate
#'   for each c1. Default: \code{seq(0.90, 2.0, by = 0.05)}.
#' @param effect_scale Character. \code{"ratio"} (default) or
#'   \code{"difference"}.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{c1}{Screening threshold.}
#'     \item{min_c2}{Minimum consistency threshold achieving the target.}
#'     \item{FPR_corr}{Corrected FPR at this (c1, min_c2).}
#'     \item{P1}{Per-subgroup detection probability.}
#'     \item{L_eff}{Effective number of candidates at this N.}
#'   }
#'   Sorted by c1 ascending.  The first row is the smallest c1 that
#'   achieves the target FPR.
#'
#' @details
#' The corrected FPR is:
#' \deqn{\text{FPR}_{\text{corr}} = 1 - (1 - P_1)^{L_{\text{eff}}}}
#'
#' where \eqn{P_1} = \code{compute_detection_probability_glm(1.0, d_eff, c1, c2)}
#' and \eqn{L_{\text{eff}} = C (N / n_{\min})^\alpha}.
#'
#' **Default L_eff values (C=1, alpha=0)** assume \eqn{L_{\text{eff}} = 1},
#' which gives the per-subgroup approximation with no multiplicity
#' correction.  This is appropriate when:
#' \itemize{
#'   \item No calibration data is available
#'   \item The analysis has few confounders and small N
#' }
#'
#' For more accurate suggestions, calibrate \eqn{L_{\text{eff}}} from
#' an H0 simulation using \code{\link{calibrate_L_eff}} or
#' \code{\link{run_null_calibration}}.
#'
#' **Reference calibrations** (from the forestsearch GLM extension
#' analysis documents):
#' \describe{
#'   \item{Binary, 4 confounders}{C = 0.220, alpha = 1.298
#'     (binary_threshold_calibration.qmd, 5000 sims)}
#'   \item{GBSG survival, 7 confounders}{C = 0.029, alpha = 0.882
#'     (cox_vs_glm_approximation.qmd, 500 sims per N)}
#' }
#'
#' @examples
#' \donttest{
#' # Survival: paper defaults should be recovered
#' d <- d_eff_survival(n_sg = 60, prop_cens = 0.45)
#' suggest_thresholds(d_eff = d, N = 700)
#'
#' # Binary: needs higher thresholds
#' d <- d_eff_binary(n_sg = 150, p_event = 0.30)
#' suggest_thresholds(d_eff = d, N = 500,
#'                    L_eff_C = 0.220, L_eff_alpha = 1.298)
#'
#' # Use calibration object from run_null_calibration()
#' # cal <- run_null_calibration(...)
#' # suggest_thresholds(d_eff = d, N = 700,
#' #                    L_eff_C = cal$C, L_eff_alpha = cal$alpha)
#' }
#'
#' @seealso \code{\link{compute_detection_probability_glm}},
#'   \code{\link{calibrate_L_eff}},
#'   \code{\link{predict_fpr_corrected}}
#'
#' @export
suggest_thresholds <- function(
    d_eff,
    N,
    fpr_target = 0.10,
    n_min = 60,
    L_eff_C = 1.0,
    L_eff_alpha = 0.0,
    c1_range = c(1.05, 2.50),
    c1_step = 0.05,
    c2_candidates = seq(0.90, 2.0, by = 0.05),
    effect_scale = c("ratio", "difference")
) {

  effect_scale <- match.arg(effect_scale)

  # --- Validation ---
  stopifnot(d_eff > 0, N > 0, fpr_target > 0, fpr_target < 1)
  stopifnot(n_min > 0, L_eff_C > 0, L_eff_alpha >= 0)

  # --- L_eff ---
  L_eff <- L_eff_C * (N / n_min)^L_eff_alpha

  # --- Search grid ---
  c1_seq <- seq(c1_range[1], c1_range[2], by = c1_step)

  # --- For each c1, find minimum c2 meeting the target ---
  frontier <- data.frame(
    c1 = numeric(0), min_c2 = numeric(0),
    FPR_corr = numeric(0), P1 = numeric(0),
    L_eff = numeric(0)
  )

  for (c1_val in c1_seq) {
    # Only consider c2 <= c1
    c2_valid <- c2_candidates[c2_candidates <= c1_val]
    if (length(c2_valid) == 0) next

    # Start from smallest c2 and find first that meets target
    for (c2_val in sort(c2_valid)) {
      P1 <- compute_detection_probability_glm(
        theta = if (effect_scale == "ratio") 1.0 else 0.0,
        d_eff = d_eff, c1 = c1_val, c2 = c2_val,
        effect_scale = effect_scale, method = "monte_carlo"
      )
      fpr_corr <- 1 - (1 - P1)^L_eff

      if (fpr_corr <= fpr_target) {
        frontier <- rbind(frontier, data.frame(
          c1 = c1_val, min_c2 = c2_val,
          FPR_corr = fpr_corr, P1 = P1,
          L_eff = L_eff
        ))
        break
      }
    }
  }

  if (nrow(frontier) == 0) {
    warning(sprintf(
      paste0("No (c1, c2) pair achieves FPR <= %.2f with ",
             "d_eff=%.1f, L_eff=%.1f. Consider larger c1_range."),
      fpr_target, d_eff, L_eff
    ))
  }

  frontier
}
