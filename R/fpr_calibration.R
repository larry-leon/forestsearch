# =============================================================================
# fpr_calibration.R
#
# Direct simulation-based FPR calibration for a specific (c1, c2) threshold
# pair and a specific ForestSearch analysis configuration.
#
# Rationale
# ---------
# The L_eff power-law correction assumes that the implied effective candidate
# count is approximately constant across threshold pairs at a given N
# (threshold-independence).  In practice this assumption fails — particularly
# for binary outcomes with small d_eff — so a universal L_eff overestimates
# the correction at the operational thresholds used for selection.
#
# fpr_calibration() bypasses the power-law entirely.  It estimates the
# procedure-level FPR directly at the analyst's chosen (c1, c2) by running
# ForestSearch under H0 a small number of times using either:
#   (a) permutation of treatment within the analyst's actual dataset, or
#   (b) a user-supplied null data-generating function.
#
# It also computes the pair-specific implied L_eff as a diagnostic, and
# returns the per-subgroup approximation P1 for comparison.
#
# Supports: binary (logistic, OR / RR / RD), continuous (Gaussian, MD),
#           count (Poisson + offset, IRR), survival (Cox PH, HR).
#
# Author: Larry F. Leon  |  forestsearch >= 0.2.0
# =============================================================================

#' Direct Simulation-Based FPR Calibration at a Specific Threshold Pair
#'
#' Estimates the procedure-level false-positive rate (FPR) for a
#' \code{\link{forestsearch}} analysis at a **specific** \eqn{(c_1, c_2)}
#' threshold pair by running \code{nsims} null simulations.  No power-law
#' model or threshold-independence assumption is required.
#'
#' Two null-generation modes are supported:
#' \describe{
#'   \item{\code{null_fn} provided}{Each replicate draws a fresh null dataset
#'     from \code{null_fn(seed)}.  Use this when you have a parametric data
#'     generating mechanism (DGM) or want to mimic a specific study design.}
#'   \item{\code{df.analysis} only (permutation mode)}{Each replicate
#'     permutes the treatment column within the analyst's actual dataset.
#'     This preserves the observed covariate and outcome distributions and
#'     is the natural H0 reference for a real analysis.  Requires the
#'     \code{treat.name} column to be present in \code{df.analysis}.}
#' }
#'
#' @param c1 Numeric scalar. Screening threshold (hr.threshold).
#' @param c2 Numeric scalar. Consistency threshold (hr.consistency).
#'   Must satisfy \code{c2 <= c1}.
#' @param df.analysis Data frame. The analyst's dataset.  Used directly in
#'   permutation mode; also used to extract N for the per-subgroup
#'   approximation.  Can be \code{NULL} if \code{null_fn} is supplied and
#'   \code{n_sample} is given explicitly.
#' @param fs_params Named list.  All forestsearch() arguments except
#'   \code{hr.threshold} and \code{hr.consistency} (which are set to
#'   \code{c1}/\code{c2}) and \code{df.analysis} (supplied separately).
#'   Must include at minimum \code{outcome.name}, \code{treat.name},
#'   \code{id.name}, and — for binary/continuous/count —
#'   \code{outcome_type} and \code{effect_measure}.
#' @param nsims Integer. Number of null replicates. Default 200.
#'   At FPR ~ 10%, \code{nsims = 200} gives a 95% CI half-width of
#'   ~±2.1 pp; \code{nsims = 500} gives ~±1.3 pp.
#' @param nsims_verify Integer. Number of independent verification
#'   replicates.  Uses a separate seed range from the calibration batch
#'   to provide an honest check of the corrected FPR.  Default 0 (skip
#'   verification).  Set to \code{nsims} or larger for production use.
#' @param null_fn Function or \code{NULL}.  If provided, called as
#'   \code{null_fn(seed)} and must return a data frame with the same
#'   column structure as \code{df.analysis}.  Overrides permutation mode.
#'   Default \code{NULL} (permutation mode).
#' @param n_sample Integer or \code{NULL}.  Sample size used when
#'   \code{null_fn} is provided.  Ignored in permutation mode (N is taken
#'   from \code{nrow(df.analysis)}).
#' @param n_min Integer. Minimum subgroup size for the per-subgroup
#'   approximation \eqn{P_1}.  Should match \code{fs_params$n.min}.
#'   Default 60.
#' @param outcome_type_approx Character. One of \code{"binary"},
#'   \code{"survival"}, \code{"count"}, \code{"continuous"}.  Used only
#'   for computing the asymptotic \eqn{P_1} approximation.  If \code{NULL}
#'   (default), inferred from \code{fs_params$outcome_type}.
#' @param d_eff_approx Numeric or \code{NULL}.  Effective information for
#'   \eqn{P_1}.  If \code{NULL} (default), estimated from \code{df.analysis}
#'   (or \code{n_sample}) using \code{d_eff_*()}.  Supply explicitly if you
#'   want to override (e.g. use an expected event count rather than observed).
#' @param prop_cens Numeric. Censoring proportion for survival outcomes.
#'   Used only when \code{d_eff_approx} is not supplied and
#'   \code{outcome_type_approx = "survival"}.  Default 0.
#' @param p_event Numeric. Event probability for binary outcomes.
#'   Used only when \code{d_eff_approx} is not supplied and
#'   \code{outcome_type_approx = "binary"}.  Default 0.5.
#' @param sigma_y Numeric. Residual SD for continuous outcomes.
#'   Default 1.
#' @param sim_seed Integer. Base random seed.  Each replicate uses
#'   \code{sim_seed + i}.  Default 42.
#' @param n_workers Integer. Parallel workers for the foreach loop.
#'   Default 1 (sequential).  Set to \code{parallel::detectCores() - 1}
#'   for production runs.
#' @param quiet Logical. Suppress progress messages. Default \code{FALSE}.
#'
#' @return An object of class \code{"fpr_calibration"} with components:
#' \describe{
#'   \item{\code{fpr}}{Numeric. Estimated procedure-level FPR.}
#'   \item{\code{fpr_ci}}{Numeric vector of length 2. 95% Wilson CI.}
#'   \item{\code{n_found}}{Integer. Number of replicates where a subgroup
#'     was found.}
#'   \item{\code{nsims}}{Integer. Total replicates attempted.}
#'   \item{\code{n_failed}}{Integer. Replicates that errored.}
#'   \item{\code{P1}}{Numeric. Per-subgroup approximation \eqn{P_1} at
#'     \code{(c1, c2)}.}
#'   \item{\code{d_eff}}{Numeric. Effective information used for \eqn{P_1}.}
#'   \item{\code{L_eff_implied}}{Numeric. Implied \eqn{L_{\text{eff}}} from
#'     the simulation result: \eqn{\log(1 - \text{FPR}) / \log(1 - P_1)}.}
#'   \item{\code{c1, c2}}{Numeric. The threshold pair.}
#'   \item{\code{N}}{Integer. Sample size.}
#'   \item{\code{mode}}{Character. \code{"permutation"} or \code{"dgm"}.}
#'   \item{\code{call_args}}{List. Key call arguments for reproducibility.}
#' }
#'
#' @details
#' ## Permutation mode
#' Under the null hypothesis of no subgroup heterogeneity, the treatment
#' assignment is exchangeable with the observed outcomes.  Permuting
#' \code{treat.name} within each replicate breaks any true treatment-by-
#' subgroup interaction while preserving the marginal distributions of all
#' other variables.  This is a model-free H0 reference that adapts
#' automatically to the actual dataset structure.
#'
#' ## P1 approximation
#' The asymptotic per-subgroup \eqn{P_1} at \code{(c1, c2)} is computed
#' via \code{compute_detection_probability_glm()} using the effective
#' information \eqn{d_{\text{eff}}} at \code{n.min}.  The ratio
#' \eqn{L_{\text{eff}} = \log(1 - \hat{\text{FPR}}) / \log(1 - P_1)}
#' is the implied multiplicity factor specific to this threshold pair and
#' dataset — it replaces the universal power-law model and does not assume
#' threshold-independence.
#'
#' @examples
#' \dontrun{
#' # ---- Binary example (permutation mode) ----
#' set.seed(1)
#' df <- data.frame(
#'   id = 1:500, treat = rbinom(500, 1, 0.5),
#'   bm1 = as.factor(rbinom(500, 1, 0.7)),
#'   bm2 = as.factor(rbinom(500, 1, 0.5)),
#'   age = round(rnorm(500, 55, 10)),
#'   y   = rbinom(500, 1, 0.30)   # null: no treatment effect
#' )
#'
#' fs_params <- list(
#'   confounders.name = c("bm1", "bm2", "age"),
#'   outcome.name = "y", event.name = "y",
#'   treat.name = "treat", id.name = "id",
#'   outcome_type = "binary", effect_measure = "OR",
#'   adverse_outcome = TRUE,
#'   pconsistency.threshold = 0.90, fs.splits = 400,
#'   n.min = 60, maxk = 2,
#'   use_lasso = TRUE, use_grf = TRUE,
#'   is.RCT = TRUE, details = FALSE, quiet = TRUE,
#'   parallel_args = list(plan = "sequential", workers = 1)
#' )
#'
#' res <- fpr_calibration(
#'   c1 = 1.70, c2 = 1.70,
#'   df.analysis = df,
#'   fs_params = fs_params,
#'   nsims = 200,
#'   n_min = 60, outcome_type_approx = "binary", p_event = 0.30,
#'   sim_seed = 42, n_workers = 4
#' )
#' print(res)
#' }
#'
#' @seealso \code{\link{forestsearch}},
#'   \code{\link{compute_detection_probability_glm}},
#'   \code{\link{d_eff_binary}}, \code{\link{d_eff_survival}}
#'
#' @export
fpr_calibration <- function(c1,
                             c2,
                             df.analysis   = NULL,
                             fs_params     = list(),
                             nsims         = 200L,
                             nsims_verify  = 0L,
                             null_fn       = NULL,
                             n_sample      = NULL,
                             n_min         = 60L,
                             outcome_type_approx = NULL,
                             d_eff_approx  = NULL,
                             prop_cens     = 0,
                             p_event       = 0.5,
                             sigma_y       = 1,
                             sim_seed      = 42L,
                             n_workers     = 1L,
                             quiet         = FALSE) {

  # ── Input validation ────────────────────────────────────────────
  stopifnot(
    "c1 must be a positive scalar"   = is.numeric(c1) && length(c1) == 1 && c1 > 0,
    "c2 must be a positive scalar"   = is.numeric(c2) && length(c2) == 1 && c2 > 0,
    "c2 must be <= c1"               = c2 <= c1,
    "nsims must be a positive integer" = is.numeric(nsims) && nsims >= 1,
    "nsims_verify must be non-negative" = is.numeric(nsims_verify) && nsims_verify >= 0,
    "fs_params must be a list"       = is.list(fs_params)
  )

  mode <- if (!is.null(null_fn)) "dgm" else "permutation"

  if (mode == "permutation") {
    if (is.null(df.analysis))
      stop("df.analysis must be supplied when null_fn = NULL (permutation mode).")
    treat_col <- fs_params$treat.name %||% "treat"
    if (!treat_col %in% names(df.analysis))
      stop(sprintf("treat column '%s' not found in df.analysis.", treat_col))
    N <- nrow(df.analysis)
  } else {
    if (is.null(n_sample) && !is.null(df.analysis))
      n_sample <- nrow(df.analysis)
    if (is.null(n_sample))
      stop("n_sample must be supplied when null_fn is provided and df.analysis is NULL.")
    N <- n_sample
  }

  # ── Infer outcome type for approximation ───────────────────────
  otype <- outcome_type_approx %||% fs_params$outcome_type %||% "binary"

  # ── Compute d_eff and P1 ────────────────────────────────────────
  if (is.null(d_eff_approx)) {
    d_eff_approx <- switch(otype,
      "survival"   = d_eff_survival(n_min, prop_cens),
      "binary"     = d_eff_binary(n_min, p_event),
      "count"      = d_eff_binary(n_min, p_event),   # approximate
      "continuous" = d_eff_continuous(n_min, sigma_y),
      d_eff_binary(n_min, p_event)  # fallback
    )
  }

  P1 <- tryCatch(
    compute_detection_probability_glm(
      theta = 1.0, d_eff = d_eff_approx,
      c1 = c1, c2 = c2,
      effect_scale = "ratio", method = "monte_carlo"
    ),
    error = function(e) NA_real_
  )

  # ── Build base forestsearch call args ───────────────────────────
  fs_args_base <- modifyList(fs_params, list(
    hr.threshold  = c1,
    hr.consistency = c2,
    details = FALSE,
    quiet   = TRUE,
    plot.sg = FALSE
  ))
  if (is.null(fs_args_base$parallel_args))
    fs_args_base$parallel_args <- list(
      plan = "sequential", workers = 1, show_message = FALSE)

  # ── Set up parallel workers ──────────────────────────────────────
  if (n_workers > 1L) {
    future::plan("multisession", workers = n_workers)
    on.exit({
      future::plan("sequential")
      gc(verbose = FALSE)
    }, add = TRUE)
  }

  if (!quiet)
    message(sprintf(
      "fpr_calibration: c1=%.2f, c2=%.2f | mode=%s | N=%d | nsims=%d%s | workers=%d",
      c1, c2, mode, N, nsims,
      if (nsims_verify > 0) sprintf(" + %d verify", nsims_verify) else "",
      n_workers))

  treat_col <- fs_params$treat.name %||% "treat"

  # ── Internal helper: run a batch of null simulations ─────────────
  .run_batch <- function(n_batch, seed_offset) {
    i <- NULL  # suppress R CMD check NOTE: foreach iterator variable
    raw <- foreach::foreach(
      i = seq_len(n_batch),
      .errorhandling = "pass",
      .options.future = list(packages = c("forestsearch"), seed = TRUE)
    ) %dofuture% {
      future::plan("sequential")
      seed_i <- sim_seed + seed_offset + i
      null_df <- if (mode == "dgm") {
        null_fn(seed_i)
      } else {
        set.seed(seed_i)
        df_perm <- df.analysis
        df_perm[[treat_col]] <- sample(df_perm[[treat_col]])
        df_perm
      }
      fs_call <- c(list(df.analysis = null_df, seedit = seed_i), fs_args_base)
      fs_out  <- tryCatch(do.call(forestsearch, fs_call), error = function(e) e)
      if (inherits(fs_out, "error")) fs_out
      else list(found = as.integer(
        !is.null(fs_out$sg.harm) && length(fs_out$sg.harm) > 0))
    }
    is_err    <- vapply(raw, inherits, logical(1), "error")
    good      <- raw[!is_err]
    found_vec <- vapply(good, function(x) x$found, integer(1))
    list(n_found = sum(found_vec), n_success = length(found_vec),
         n_failed = sum(is_err))
  }

  # ── Batch 1: calibration sims (derive P1, L_eff, corrected FPR) ──
  b1 <- .run_batch(nsims, seed_offset = 0L)

  if (b1$n_failed > 0 && !quiet)
    warning(sprintf(
      "fpr_calibration (calibration batch): %d of %d replicates failed.",
      b1$n_failed, nsims))

  fpr_hat <- if (b1$n_success > 0) b1$n_found / b1$n_success else NA_real_
  fpr_ci  <- .wilson_ci(b1$n_found, b1$n_success)

  L_eff_implied <- if (!is.na(fpr_hat) && !is.na(P1) &&
                       fpr_hat > 0 && fpr_hat < 1 && P1 > 0 && P1 < 1)
    log(1 - fpr_hat) / log(1 - P1) else NA_real_

  fpr_corrected <- if (!is.na(P1) && !is.na(L_eff_implied))
    1 - (1 - P1)^L_eff_implied else NA_real_

  # ── Batch 2: verification sims (independent check of corrected FPR) ──
  # Seeds are offset by nsims + 10000 to ensure independence from batch 1.
  fpr_verify    <- NA_real_
  fpr_verify_ci <- c(NA_real_, NA_real_)
  n_verify_success <- 0L
  n_verify_failed  <- 0L

  if (nsims_verify > 0L) {
    if (!quiet)
      message(sprintf(
        "fpr_calibration: running %d verification sims (independent seed)...",
        nsims_verify))
    b2 <- .run_batch(nsims_verify, seed_offset = nsims + 10000L)
    if (b2$n_failed > 0 && !quiet)
      warning(sprintf(
        "fpr_calibration (verification batch): %d of %d replicates failed.",
        b2$n_failed, nsims_verify))
    fpr_verify       <- if (b2$n_success > 0) b2$n_found / b2$n_success else NA_real_
    fpr_verify_ci    <- .wilson_ci(b2$n_found, b2$n_success)
    n_verify_success <- b2$n_success
    n_verify_failed  <- b2$n_failed
  }

  # ── Return structured result ─────────────────────────────────────
  structure(
    list(
      # Calibration batch (source of L_eff and corrected FPR)
      fpr           = fpr_hat,
      fpr_ci        = fpr_ci,
      n_found       = b1$n_found,
      nsims         = b1$n_success,
      n_failed      = b1$n_failed,
      P1            = P1,
      d_eff         = d_eff_approx,
      L_eff_implied = L_eff_implied,
      fpr_corrected = fpr_corrected,    # 1-(1-P1)^L_eff; tautological at this pair
      # Verification batch (independent; validates corrected FPR)
      fpr_verify       = fpr_verify,
      fpr_verify_ci    = fpr_verify_ci,
      nsims_verify     = n_verify_success,
      n_verify_failed  = n_verify_failed,
      # Metadata
      c1           = c1,
      c2           = c2,
      N            = N,
      mode         = mode,
      call_args    = list(
        nsims        = nsims,
        nsims_verify = nsims_verify,
        n_min        = n_min,
        otype        = otype,
        n_workers    = n_workers,
        sim_seed     = sim_seed
      )
    ),
    class = "fpr_calibration"
  )
}


# ── S3 methods ───────────────────────────────────────────────────────────────

#' Print method for \code{fpr_calibration} objects
#' @param x An \code{fpr_calibration} object.
#' @param ... Unused; present for S3 compatibility.
#' @return The input \code{x}, invisibly.
#' @export
print.fpr_calibration <- function(x, ...) {
  cat("\n=== fpr_calibration result ===\n")
  cat(sprintf("  Thresholds:     c1 = %.2f, c2 = %.2f\n", x$c1, x$c2))
  cat(sprintf("  Mode:           %s  (N = %d)\n", x$mode, x$N))
  cat(sprintf("  Cal sims:       %d completed, %d failed\n",
              x$nsims, x$n_failed))

  cat(sprintf("\n  Sim FPR (cal):     %.1f%%  [95%% CI: %.1f%% -- %.1f%%]\n",
              100 * x$fpr, 100 * x$fpr_ci[1], 100 * x$fpr_ci[2]))
  cat(sprintf("  P1 (per-subgroup): %.1f%%  (d_eff = %.1f, uncorrected)\n",
              100 * x$P1, x$d_eff))
  if (!is.na(x$L_eff_implied)) {
    cat(sprintf("  Implied L_eff:     %.3f\n", x$L_eff_implied))
    cat(sprintf("  Corrected FPR:     %.1f%%  [= 1-(1-P1)^L_eff]\n",
                100 * x$fpr_corrected))
  }
  if (!is.null(x$fpr_verify) && !is.na(x$fpr_verify)) {
    cat(sprintf("\n  Sim FPR (verify):  %.1f%%  [95%% CI: %.1f%% -- %.1f%%]  (%d independent sims)\n",
                100 * x$fpr_verify,
                100 * x$fpr_verify_ci[1], 100 * x$fpr_verify_ci[2],
                x$nsims_verify))
    in_ci <- x$fpr_verify >= x$fpr_verify_ci[1] &&
             x$fpr_verify <= x$fpr_verify_ci[2]
    cat(sprintf("  Corrected vs verify: %.1f%% vs %.1f%%\n",
                100 * x$fpr_corrected, 100 * x$fpr_verify))
  }
  cat("\n")
  invisible(x)
}


#' Summary method for \code{fpr_calibration} objects
#'
#' Delegates to \code{\link{print.fpr_calibration}}.
#'
#' @param object An \code{fpr_calibration} object.
#' @param ... Passed to the print method.
#' @return The input \code{object}, invisibly.
#' @export
summary.fpr_calibration <- function(object, ...) {
  print(object, ...)
}


# =============================================================================
# fpr_calibration_transfer(): predict corrected FPR at (c1+delta, c2+delta)
# using L_eff_implied from a calibrated fpr_calibration object.
# =============================================================================

#' Transfer Corrected FPR Prediction to a Nearby Threshold Pair
#'
#' Uses the implied \eqn{L_{\text{eff}}} from a calibrated
#' \code{fpr_calibration} object to predict the procedure-level FPR at
#' \eqn{(c_1 + \delta_1, c_2 + \delta_2)}, then optionally runs a direct
#' simulation at the new pair to validate the prediction.
#'
#' @param cal_result An object of class \code{"fpr_calibration"} produced
#'   by \code{\link{fpr_calibration}}.  Must contain a finite
#'   \code{L_eff_implied}.
#' @param delta Numeric scalar.  Added to both \eqn{c_1} and \eqn{c_2}.
#'   Positive values tighten both thresholds; negative values relax them.
#'   Ignored if \code{delta_c1} or \code{delta_c2} are supplied.
#' @param delta_c1 Numeric or \code{NULL}.  Increment for \eqn{c_1} only.
#' @param delta_c2 Numeric or \code{NULL}.  Increment for \eqn{c_2} only.
#' @param nsims_check Integer.  Simulations to run at the new pair for
#'   empirical validation.  Default 200.  Set to 0 to skip.
#' @param ... Additional arguments passed to \code{fpr_calibration} for the
#'   validation run (e.g. \code{n_workers}, \code{quiet}).
#'
#' @return A data frame with one row per delta value containing:
#'   \code{c1_new}, \code{c2_new}, \code{P1_new}, \code{fpr_pred}
#'   (predicted corrected FPR using transferred \code{L_eff_implied}),
#'   \code{fpr_sim} (direct simulation, if \code{nsims_check > 0}),
#'   \code{fpr_sim_lo}, \code{fpr_sim_hi} (95\% CI), \code{in_ci}
#'   (whether \code{fpr_pred} falls within the simulation CI).
#'
#' @examples
#' \dontrun{
#' # Calibrate at base pair
#' cal <- fpr_calibration(c1 = 1.70, c2 = 1.70, null_fn = my_dgm,
#'                         n_sample = 500, fs_params = fs_binary,
#'                         nsims = 200, n_min = 60,
#'                         outcome_type_approx = "binary", p_event = 0.30)
#'
#' # Predict at base + 0.05, +0.10, +0.20 and validate
#' transfer_results <- fpr_calibration_transfer(
#'   cal, delta = c(0.05, 0.10, 0.20),
#'   nsims_check = 200, n_workers = 4)
#'
#' print(transfer_results)
#' }
#' @export
fpr_calibration_transfer <- function(cal_result,
                                      delta    = c(0.05, 0.10, 0.20),
                                      delta_c1 = NULL,
                                      delta_c2 = NULL,
                                      nsims_check = 200L,
                                      ...) {

  stopifnot(inherits(cal_result, "fpr_calibration"))
  if (is.na(cal_result$L_eff_implied))
    stop("cal_result has no finite L_eff_implied; cannot transfer.")

  # Resolve per-threshold deltas
  if (!is.null(delta_c1) || !is.null(delta_c2)) {
    dc1 <- if (!is.null(delta_c1)) delta_c1 else delta
    dc2 <- if (!is.null(delta_c2)) delta_c2 else delta
    if (length(dc1) != length(dc2))
      stop("delta_c1 and delta_c2 must have the same length.")
  } else {
    dc1 <- delta
    dc2 <- delta
  }

  L_eff <- cal_result$L_eff_implied
  d_eff <- cal_result$d_eff

  results <- lapply(seq_along(dc1), function(k) {
    c1_new <- cal_result$c1 + dc1[k]
    c2_new <- cal_result$c2 + dc2[k]

    # New P1 at the shifted pair
    P1_new <- tryCatch(
      compute_detection_probability_glm(
        theta = 1.0, d_eff = d_eff,
        c1 = c1_new, c2 = c2_new,
        effect_scale = "ratio", method = "monte_carlo"),
      error = function(e) NA_real_)

    # Predicted corrected FPR: transfer L_eff from base pair
    fpr_pred <- if (!is.na(P1_new)) 1 - (1 - P1_new)^L_eff else NA_real_

    # Optional direct simulation at the new pair
    fpr_sim <- NA_real_
    fpr_sim_lo <- NA_real_
    fpr_sim_hi <- NA_real_
    n_sim_actual <- 0L

    if (nsims_check > 0L && c2_new <= c1_new) {
      sim_args <- c(
        list(c1 = c1_new, c2 = c2_new,
             df.analysis = cal_result$call_args$df.analysis %||%
               eval(cal_result$call_args$null_fn_expr),
             fs_params    = list(),
             nsims        = nsims_check,
             n_min        = cal_result$call_args$n_min,
             d_eff_approx = d_eff,
             sim_seed     = cal_result$call_args$sim_seed + 20000L +
               round(dc1[k] * 1000)),
        list(...))
      # Re-run using same null_fn / df.analysis context from parent call
      # (these must be supplied via ...)
      chk <- tryCatch(do.call(fpr_calibration, sim_args), error = function(e) e)
      if (!inherits(chk, "error")) {
        fpr_sim    <- chk$fpr
        fpr_sim_lo <- chk$fpr_ci[1]
        fpr_sim_hi <- chk$fpr_ci[2]
        n_sim_actual <- chk$nsims
      }
    }

    in_ci <- if (!is.na(fpr_pred) && !is.na(fpr_sim_lo))
      fpr_pred >= fpr_sim_lo & fpr_pred <= fpr_sim_hi
    else NA

    data.frame(
      delta_c1   = dc1[k],
      delta_c2   = dc2[k],
      c1_new     = c1_new,
      c2_new     = c2_new,
      P1_base    = cal_result$P1,
      P1_new     = P1_new,
      L_eff      = L_eff,
      fpr_pred   = fpr_pred,
      fpr_sim    = fpr_sim,
      fpr_sim_lo = fpr_sim_lo,
      fpr_sim_hi = fpr_sim_hi,
      n_sim      = n_sim_actual,
      in_ci      = in_ci,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}


# ── Helpers ───────────────────────────────────────────────────────────────────

# Null-coalescing operator (unexported, used internally)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Wilson (score) confidence interval for a proportion
.wilson_ci <- function(x, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p_hat <- x / n
  denom <- 1 + z^2 / n
  centre <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - margin), min(1, centre + margin))
}
