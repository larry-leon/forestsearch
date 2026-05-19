#' Calibrate k_treat for Target Overall Hazard Ratio
#'
#' Finds the treatment effect multiplier (\code{k_treat}) that achieves a
#' target overall hazard ratio in a DGM built by
#' \code{\link{generate_aft_dgm_flex}}.  The root is found by \code{uniroot}
#' over \code{k_treat} with every other argument to
#' \code{generate_aft_dgm_flex()} held fixed.
#'
#' This is the natural sibling of \code{\link{calibrate_k_inter}} (which
#' targets the harm-subgroup HR under \code{model = "alt"}) and
#' \code{\link{calibrate_cens_adjust}} (which targets the censoring rate).
#' It is typically used with \code{model = "null"} to fix the overall HR
#' for a uniform-effect DGM, but also works with \code{model = "alt"}
#' when the marginal HR is the calibration target rather than the
#' subgroup-specific HR.
#'
#' @param target_hr_overall Numeric. Target overall hazard ratio
#'   (must be positive).
#' @param base_args Named list of arguments to
#'   \code{generate_aft_dgm_flex()} excluding \code{k_treat}.  Every other
#'   parameter (data, covariates, model, subgroup_vars, seed, n_super,
#'   etc.) is taken from this list.  See Details below for reproducibility
#'   considerations.
#' @param k_treat_range Numeric vector of length 2.  Search range for
#'   \code{k_treat} passed to \code{uniroot}.  Default \code{c(-5, 5)} is
#'   typically wide enough; widen if the boundary HRs do not bracket the
#'   target (uniroot will error otherwise).
#' @param tol Numeric.  Tolerance for root finding.  Default \code{1e-6}.
#' @param use_ahr Logical.  If \code{TRUE}, calibrate to the average
#'   hazard ratio (\code{dgm$hazard_ratios$AHR}) instead of the Cox-based
#'   overall HR (\code{dgm$hazard_ratios$overall}).  Default \code{FALSE}.
#' @param verbose Logical.  Print diagnostic information.  Default
#'   \code{FALSE}.
#'
#' @return Numeric scalar.  The calibrated \code{k_treat} value.  Returns
#'   \code{NA_real_} with a warning if root finding fails — most commonly
#'   when the target is not bracketed by \code{k_treat_range}.
#'
#' @details
#' Reproducibility requires the \code{uniroot} objective to be deterministic
#' in \code{k_treat}.  This is achieved by including an integer
#' \code{seed} field in \code{base_args} so that every call to
#' \code{generate_aft_dgm_flex()} during the search uses the same
#' random-number stream.  If \code{base_args$seed} is \code{NULL} the
#' objective becomes stochastic and \code{uniroot} may not converge; the
#' function emits a warning in this case.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data(gbsg)
#' gbsg$time_months <- gbsg$rfstime / 30.4375
#'
#' base_args <- list(
#'   data            = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars     = c("meno", "grade"),
#'   outcome_var     = "time_months",
#'   event_var       = "status",
#'   treatment_var   = "hormon",
#'   model           = "null",
#'   n_super         = 5000,
#'   seed            = 99,
#'   verbose         = FALSE
#' )
#'
#' # Calibrate k_treat so that the overall HR is exactly 0.70
#' k <- calibrate_k_treat(target_hr_overall = 0.70,
#'                        base_args         = base_args,
#'                        verbose           = TRUE)
#'
#' # Verify
#' dgm <- do.call(generate_aft_dgm_flex,
#'                c(base_args, list(k_treat = k)))
#' print(dgm$hazard_ratios$overall)   # ~ 0.70
#'
#' # Calibrate to AHR instead
#' k_ahr <- calibrate_k_treat(target_hr_overall = 0.70,
#'                            base_args         = base_args,
#'                            use_ahr           = TRUE)
#' dgm_ahr <- do.call(generate_aft_dgm_flex,
#'                    c(base_args, list(k_treat = k_ahr)))
#' print(dgm_ahr$hazard_ratios$AHR)   # ~ 0.70
#' }
#'
#' @seealso \code{\link{calibrate_k_inter}},
#'   \code{\link{calibrate_cens_adjust}},
#'   \code{\link{generate_aft_dgm_flex}}
#' @importFrom stats uniroot
#' @importFrom utils modifyList
#' @export
calibrate_k_treat <- function(target_hr_overall,
                              base_args,
                              k_treat_range = c(-5, 5),
                              tol           = 1e-6,
                              use_ahr       = FALSE,
                              verbose       = FALSE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  stopifnot(
    "target_hr_overall must be positive" = is.numeric(target_hr_overall) &&
                                            length(target_hr_overall) == 1L &&
                                            target_hr_overall > 0,
    "base_args must be a list"           = is.list(base_args),
    "base_args must not contain k_treat" = is.null(base_args$k_treat),
    "k_treat_range must be length 2"     = length(k_treat_range) == 2L,
    "k_treat_range must be increasing"   = k_treat_range[1] < k_treat_range[2],
    "tol must be positive"               = is.numeric(tol) && tol > 0,
    "use_ahr must be logical"            = is.logical(use_ahr) &&
                                            length(use_ahr) == 1L
  )

  if (is.null(base_args$seed)) {
    warning("base_args$seed is NULL; uniroot objective may be noisy in ",
            "k_treat. Set an integer seed in base_args for a reproducible ",
            "calibration.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Objective:  realised HR (or AHR) minus target
  #
  # NB: modifyList() (not c()) — base_args may already carry `verbose`,
  # and we want to *override* it during calibration rather than create a
  # duplicate argument that do.call() would reject.
  # ---------------------------------------------------------------------------
  objective <- function(k_val) {
    args_call <- utils::modifyList(
      base_args,
      list(k_treat = k_val, verbose = FALSE)
    )
    dgm <- do.call(generate_aft_dgm_flex, args_call)
    if (use_ahr) {
      dgm$hazard_ratios$AHR     - target_hr_overall
    } else {
      dgm$hazard_ratios$overall - target_hr_overall
    }
  }

  hr_type <- if (use_ahr) "AHR" else "HR(overall)"

  if (verbose) {
    message(sprintf("Calibrating k_treat to achieve %s = %.4f",
                    hr_type, target_hr_overall))
    message(sprintf("Search range: [%.2f, %.2f]",
                    k_treat_range[1], k_treat_range[2]))

    # Evaluate at boundaries so the user sees the bracket interval
    lo_val <- objective(k_treat_range[1]) + target_hr_overall
    hi_val <- objective(k_treat_range[2]) + target_hr_overall
    message(sprintf("%s at k_treat = %.2f: %.4f",
                    hr_type, k_treat_range[1], lo_val))
    message(sprintf("%s at k_treat = %.2f: %.4f",
                    hr_type, k_treat_range[2], hi_val))
  }

  # ---------------------------------------------------------------------------
  # Root finding
  # ---------------------------------------------------------------------------
  result <- tryCatch(
    stats::uniroot(objective, interval = k_treat_range, tol = tol),
    error = function(e) {
      warning("Root finding failed: ", e$message,
              "\nWiden k_treat_range if the boundary HRs do not bracket the target.",
              call. = FALSE)
      NULL
    }
  )

  if (is.null(result)) return(NA_real_)

  k_treat <- result$root

  if (verbose) {
    args_verify <- utils::modifyList(
      base_args,
      list(k_treat = k_treat, verbose = FALSE)
    )
    dgm_verify <- do.call(generate_aft_dgm_flex, args_verify)
    achieved <- if (use_ahr) {
      dgm_verify$hazard_ratios$AHR
    } else {
      dgm_verify$hazard_ratios$overall
    }
    message("")
    message("=== Calibration Result ===")
    message(sprintf("Found k_treat = %.6f", k_treat))
    message(sprintf("Achieved %s = %.4f (target: %.4f)",
                    hr_type, achieved, target_hr_overall))
    message(sprintf("Error: %.6f", abs(achieved - target_hr_overall)))
    message(sprintf("Iterations: %d", result$iter))
  }

  k_treat
}
