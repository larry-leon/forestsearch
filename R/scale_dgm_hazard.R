## =============================================================================
## scale_dgm_hazard.R
##
## Multiply every hazard of a Weibull AFT DGM by a constant m, leaving every
## hazard ratio unchanged: mu <- mu_fitted - tau * log(m).
## =============================================================================

#' Multiply every hazard of a Weibull AFT DGM by a constant
#'
#' Sets the AFT intercept of a DGM to `mu_fitted - tau * log(m)`, where
#' `mu_fitted` is the intercept as fitted when the DGM was built. Under the
#' Weibull AFT model the baseline cumulative hazard is
#' \eqn{\Lambda_0(s) = (s e^{-\mu})^{1/\tau}}, so every patient's hazard
#' under either arm is multiplied by `m` and every hazard ratio (individual,
#' AHR, CDE) is unchanged; only event probabilities in a finite window move.
#' The scaling is not cumulative: it is always applied to the fitted
#' intercept, so `scale_dgm_hazard(scale_dgm_hazard(dgm, 2), 3)` equals
#' `scale_dgm_hazard(dgm, 3)`.
#'
#' The fitted intercept is taken, in order, from `dgm$mrct$mu_fitted` (stored
#' by an earlier call), `attr(dgm, "mu_fitted")` (the attribute the
#' mrct-consistency builders set), or, for a DGM never scaled,
#' `dgm$model_params$mu`. A DGM that records a multiplier but no fitted
#' intercept is refused.
#'
#' @param dgm An `aft_dgm_flex` object, typically from
#'   [inject_mrct_structure()].
#' @param m Positive number: the hazard multiplier. `m = 1` restores the
#'   fitted intercept (for an unscaled DGM, the DGM is returned unchanged
#'   except for the record).
#' @return The DGM with `model_params$mu` set, the multiplier recorded in
#'   `dgm$mrct$hazard_multiplier` and the fitted intercept in
#'   `dgm$mrct$mu_fitted`, and the attribute `"m"` set to `m` (as the
#'   mrct-consistency helper `ideal_set_m()` sets it, so scripts reading
#'   `attr(dgm, "m")` keep working).
#' @seealso [inject_mrct_structure()], [ahr_population()]
#' @examples
#' \dontrun{
#' dgm3 <- scale_dgm_hazard(dgm, 3)   # every hazard x 3, HRs unchanged
#' }
#' @export
scale_dgm_hazard <- function(dgm, m) {
  if (!inherits(dgm, c("aft_dgm_flex", "aft_dgm")))
    stop("dgm must be an object created by generate_aft_dgm_flex() or ",
         "inject_mrct_structure().", call. = FALSE)
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m <= 0)
    stop("m must be one finite positive number.", call. = FALSE)
  mp <- dgm$model_params
  if (is.null(mp$tau) || !is.finite(mp$tau) || mp$tau <= 0)
    stop("dgm$model_params$tau must be a positive number.", call. = FALSE)
  mu_fitted <- dgm$mrct$mu_fitted
  if (is.null(mu_fitted)) mu_fitted <- attr(dgm, "mu_fitted")
  if (is.null(mu_fitted)) {
    if (!is.null(dgm$mrct$hazard_multiplier) || !is.null(attr(dgm, "m")))
      stop("dgm records a hazard multiplier but not its fitted intercept ",
           "(dgm$mrct$mu_fitted or attr(dgm, \"mu_fitted\")).", call. = FALSE)
    mu_fitted <- mp$mu
  }
  if (is.null(mu_fitted) || length(mu_fitted) != 1L || !is.finite(mu_fitted))
    stop("The fitted intercept mu is missing from dgm.", call. = FALSE)
  dgm$model_params$mu <- mu_fitted - mp$tau * log(m)
  dgm$mrct$mu_fitted <- mu_fitted
  dgm$mrct$hazard_multiplier <- m
  attr(dgm, "m") <- m
  dgm
}
