# calibrate_helpers.R
# ---------------------------------------------------------------------------
# Shared root-finding backbone for the DGM effect-calibrators
# (calibrate_k_inter(), calibrate_k_treat(), calibrate_glm_interaction()).
#
# Every effect-calibrator solves the same one-dimensional problem: choose the
# scalar shift (k_inter or k_treat) so that a monotone realised effect (the
# harm-subgroup HR, the overall HR, or the subgroup OR) equals a target.  This
# helper centralises that solve so the Cox and GLM paths behave identically:
#
#   * uniroot() root-finding -- exact, not snapped to a grid;
#   * extendInt = "yes" -- the search interval is widened automatically until it
#     brackets the root, so a too-narrow range self-corrects instead of silently
#     saturating at an edge;
#   * a HARD relative-tolerance gate -- if the achieved effect is not within
#     `tol_rel` percent of the target, the call STOPS (fails loudly) rather than
#     returning a mis-calibrated DGM.
# ---------------------------------------------------------------------------

#' Solve a monotone effect-vs-shift calibration by bracketed root-finding
#'
#' Internal backbone shared by [calibrate_k_inter()], [calibrate_k_treat()], and
#' [calibrate_glm_interaction()].  Finds the scalar `k` for which
#' `effect_fun(k)` equals `target`, extending the search interval automatically
#' (via [stats::uniroot()]'s `extendInt`) and failing loudly if the target
#' cannot be met to within `tol_rel` percent.
#'
#' @param effect_fun Function of one numeric argument returning the realised
#'   effect (e.g. the harm-subgroup HR or the subgroup OR) at that shift.  It is
#'   assumed monotone in its argument over the region of interest.
#' @param target Numeric target effect (finite and non-zero).
#' @param interval Numeric length-2 initial search interval for `k`.
#' @param tol Numeric tolerance on `k` passed to [stats::uniroot()].
#'   Default `1e-6`.
#' @param tol_rel Numeric. Maximum tolerated relative error, in percent: the
#'   call stops if `100 * abs(achieved - target) / abs(target) >= tol_rel`.
#'   Default `2.5`.
#' @param label,param Character labels used in messages -- the effect name and
#'   the shift name (e.g. `"HR(H)"` and `"k_inter"`).
#' @param verbose Logical. Print the solved value and relative error.
#'
#' @return A list with elements `root` (the calibrated `k`), `achieved`,
#'   `rel_error` (percent), and `iter`.
#'
#' @keywords internal
#' @noRd
.calibrate_by_root <- function(effect_fun, target, interval,
                               tol = 1e-6, tol_rel = 2.5,
                               label = "effect", param = "k",
                               verbose = FALSE) {
  stopifnot(
    is.function(effect_fun),
    is.numeric(target), length(target) == 1L, is.finite(target), target != 0,
    length(interval) == 2L, interval[1L] < interval[2L],
    is.numeric(tol_rel), length(tol_rel) == 1L, tol_rel > 0
  )

  objective <- function(k) effect_fun(k) - target

  root <- tryCatch(
    stats::uniroot(objective, interval = interval, extendInt = "yes", tol = tol),
    error = function(e)
      stop(sprintf(
        paste0("calibration for %s could not bracket a root over [%.4g, %.4g] ",
               "even with automatic interval extension (%s). The target %s = ",
               "%.4g is likely unreachable for any %s; check the DGM/subgroup ",
               "specification."),
        label, interval[1L], interval[2L], conditionMessage(e),
        label, target, param),
        call. = FALSE)
  )

  k_hat    <- root$root
  achieved <- effect_fun(k_hat)
  rel_err  <- 100 * abs(achieved - target) / abs(target)

  if (!is.finite(rel_err) || rel_err >= tol_rel) {
    stop(sprintf(
      paste0("calibration for %s FAILED the tolerance gate: achieved %.4f vs ",
             "target %.4f (%.2f%% relative error >= %.1f%% allowed) at %s = ",
             "%.4f. Widen the search range or revise the target."),
      label, achieved, target, rel_err, tol_rel, param, k_hat),
      call. = FALSE)
  }

  if (verbose) {
    message(sprintf(
      "  %s calibrated: %s = %.6f  ->  %s = %.4f  (target %.4f; %.3f%% rel. error; %d iter)",
      label, param, k_hat, label, achieved, target, rel_err, root$iter))
  }

  list(root = k_hat, achieved = achieved, rel_error = rel_err, iter = root$iter)
}
