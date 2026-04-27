#' Approximate Procedure-Level False Positive Rate
#'
#' Computes an analytical approximation to the ForestSearch false positive
#' rate under the null hypothesis of homogeneous treatment effect.  The
#' calculation accounts for (i) the sampling variability of the effect
#' estimate in a subgroup of size \code{n_min}, (ii) the number of
#' candidate subgroups searched (driven by the covariate space and
#' \code{maxk}), and (iii) the screening threshold \code{c1}.
#'
#' The approximation assumes:
#' \itemize{
#'   \item Normal approximation to the log-OR (or log-HR) distribution.
#'   \item Independent candidates (upper bound; true FPR is lower due
#'         to overlap).  A correlation-adjusted estimate is also returned.
#'   \item Equal allocation to treatment arms within each subgroup.
#'   \item No consistency check adjustment (the reported FPR is the
#'         \emph{screening-stage} FPR before \code{pconsistency.threshold}
#'         filtering).
#' }
#'
#' @param effect_true Numeric.  True effect under H0 on the ratio scale
#'   (e.g., OR or HR).  Default 0.65 (protective treatment).
#' @param p_event Numeric in (0, 1).  Overall event rate (control arm
#'   baseline rate for binary; used to compute SE).  Default 0.40.
#' @param n_min Integer.  Minimum subgroup size (\code{n.min} in
#'   \code{forestsearch()}).  Default 60.
#' @param c1 Numeric.  Screening threshold (\code{hr.threshold}).
#'   Candidates with effect > \code{c1} pass screening.  Default 1.25.
#' @param n_confounders_cont Integer.  Number of continuous confounders.
#'   Default 6.
#' @param n_confounders_bin Integer.  Number of binary confounders.
#'   Default 6.
#' @param cont_cutoff Integer.  Number of cut points per continuous
#'   covariate (\code{cont.cutoff} in \code{forestsearch()}).  Each
#'   continuous variable produces \code{cont_cutoff} candidate
#'   subgroups.  Default 4.
#' @param maxk Integer.  Maximum interaction depth.  1 = single-variable
#'   subgroups only; 2 = up to two-variable interactions.  Default 2.
#' @param effect_type Character.  \code{"ratio"} (OR, HR, IRR: log-scale
#'   SE) or \code{"difference"} (RD, MD: identity-scale SE).
#'   Default \code{"ratio"}.
#' @param treat_alloc Numeric in (0, 1).  Treatment allocation fraction.
#'   Default 0.5 (1:1 randomisation).
#' @param rho Numeric in [0, 1).  Assumed average pairwise correlation
#'   between candidate test statistics, used for the correlation-adjusted
#'   FPR estimate.  Default 0.30.  Set to 0 for the independence
#'   (upper-bound) estimate only.
#' @param thresholds Numeric vector.  Additional thresholds to evaluate
#'   besides \code{c1}.  Default \code{c(1.0, 1.5, 2.0)}.
#' @param verbose Logical.  Print a formatted summary.  Default
#'   \code{TRUE}.
#'
#' @return A list (class \code{"fpr_approximation"}) with components:
#'   \describe{
#'     \item{\code{se_log_effect}}{SE of log(effect) in a subgroup of
#'       size \code{n_min}.}
#'     \item{\code{n_candidates}}{Total number of candidate subgroups
#'       searched.}
#'     \item{\code{n_single}}{Number of single-variable candidates.}
#'     \item{\code{n_pairs}}{Number of two-variable candidates
#'       (0 if \code{maxk == 1}).}
#'     \item{\code{per_candidate}}{Named numeric vector: P(effect > thr)
#'       for each threshold evaluated.}
#'     \item{\code{fpr_indep}}{Named numeric vector: procedure-level FPR
#'       assuming independent candidates, for each threshold.}
#'     \item{\code{fpr_adjusted}}{Named numeric vector: correlation-adjusted
#'       FPR, for each threshold.}
#'     \item{\code{params}}{List of input parameters.}
#'   }
#'
#' @examples
#' # ACTG175 binary harm search
#' fpr_approx(effect_true = 0.65, p_event = 0.40, n_min = 60, c1 = 1.25)
#'
#' # Tighter threshold
#' fpr_approx(effect_true = 0.65, p_event = 0.40, n_min = 60, c1 = 2.0)
#'
#' # Larger minimum subgroup
#' fpr_approx(effect_true = 0.65, p_event = 0.40, n_min = 100, c1 = 1.25)
#'
#' @export
fpr_approx <- function(
    effect_true        = 0.65,
    p_event            = 0.40,
    n_min              = 60L,
    c1                 = 1.25,
    n_confounders_cont = 6L,
    n_confounders_bin  = 6L,
    cont_cutoff        = 4L,
    maxk               = 2L,
    effect_type        = c("ratio", "difference"),
    treat_alloc        = 0.5,
    rho                = 0.30,
    thresholds         = NULL,
    verbose            = TRUE
) {

  effect_type <- match.arg(effect_type)

  # ── Validate inputs ────────────────────────────────────────────────────────
  stopifnot(
    is.numeric(effect_true), length(effect_true) == 1, effect_true > 0,
    is.numeric(p_event), length(p_event) == 1, p_event > 0, p_event < 1,
    is.numeric(n_min), length(n_min) == 1, n_min >= 10,
    is.numeric(c1), length(c1) == 1, c1 > 0,
    is.numeric(treat_alloc), treat_alloc > 0, treat_alloc < 1,
    is.numeric(rho), rho >= 0, rho < 1,
    maxk %in% c(1L, 2L)
  )

  # ── Per-arm sample sizes in the smallest subgroup ──────────────────────────
  n1 <- n_min * treat_alloc
  n0 <- n_min * (1 - treat_alloc)

  # ── SE of the effect estimate ──────────────────────────────────────────────
  if (effect_type == "ratio") {
    # Log-OR SE from expected 2×2 cell counts
    # Control arm: p_event is the baseline rate
    p_c <- p_event
    # Treatment arm rate implied by the true OR
    p_t <- p_c * effect_true / (1 - p_c + p_c * effect_true)

    a <- n1 * p_t          # treat + event
    b <- n1 * (1 - p_t)    # treat + no event
    c <- n0 * p_c           # control + event
    d <- n0 * (1 - p_c)     # control + no event

    # Guard against degenerate cells
    a <- max(a, 0.5); b <- max(b, 0.5)
    c <- max(c, 0.5); d <- max(d, 0.5)

    se_log <- sqrt(1/a + 1/b + 1/c + 1/d)
    log_true <- log(effect_true)

  } else {
    # Difference scale (RD, MD): SE ≈ sqrt(p(1-p)/n0 + p(1-p)/n1)
    # Simplified for binary; for continuous, user would need to supply SD
    se_log <- sqrt(p_event * (1 - p_event) * (1/n0 + 1/n1))
    log_true <- effect_true  # already on identity scale
  }

  # ── Number of candidate subgroups ──────────────────────────────────────────
  L_cont <- n_confounders_cont * cont_cutoff
  L_bin  <- n_confounders_bin * 1L  # each binary factor → 1 split
  L <- L_cont + L_bin

  if (maxk >= 2L) {
    n_pairs <- L * (L - 1L) / 2L
  } else {
    n_pairs <- 0L
  }
  n_candidates <- L + n_pairs

  # ── Per-candidate exceedance probability ───────────────────────────────────
  all_thresholds <- sort(unique(c(c1, thresholds)))

  per_candidate <- vapply(all_thresholds, function(thr) {
    if (effect_type == "ratio") {
      z <- (log(thr) - log_true) / se_log
    } else {
      z <- (thr - log_true) / se_log
    }
    stats::pnorm(z, lower.tail = FALSE)
  }, numeric(1))
  names(per_candidate) <- sprintf("%.2f", all_thresholds)

  # ── Procedure-level FPR ────────────────────────────────────────────────────
  # Independence assumption (upper bound)
  fpr_indep <- 1 - (1 - per_candidate)^n_candidates

  # Correlation-adjusted (Šidák-like with effective number of tests)
  # n_eff ≈ n_candidates / (1 + (n_candidates - 1) * rho)
  if (rho > 0 && n_candidates > 1) {
    n_eff <- n_candidates / (1 + (n_candidates - 1) * rho)
    fpr_adjusted <- 1 - (1 - per_candidate)^n_eff
  } else {
    n_eff <- n_candidates
    fpr_adjusted <- fpr_indep
  }

  # ── Output ─────────────────────────────────────────────────────────────────
  result <- structure(
    list(
      se_log_effect = se_log,
      n_candidates  = n_candidates,
      n_single      = L,
      n_pairs       = n_pairs,
      n_eff         = n_eff,
      per_candidate = per_candidate,
      fpr_indep     = fpr_indep,
      fpr_adjusted  = fpr_adjusted,
      params        = list(
        effect_true        = effect_true,
        p_event            = p_event,
        n_min              = n_min,
        c1                 = c1,
        n_confounders_cont = n_confounders_cont,
        n_confounders_bin  = n_confounders_bin,
        cont_cutoff        = cont_cutoff,
        maxk               = maxk,
        effect_type        = effect_type,
        treat_alloc        = treat_alloc,
        rho                = rho,
        p_treat            = if (effect_type == "ratio") p_t else NA,
        thresholds         = all_thresholds
      )
    ),
    class = "fpr_approximation"
  )

  if (verbose) print(result)

  invisible(result)
}


#' @export
print.fpr_approximation <- function(x, ...) {

  p <- x$params
  cat("============================================================\n")
  cat("  ForestSearch FPR Approximation (screening stage)\n")
  cat("============================================================\n\n")

  cat(sprintf("  True effect (H0):   %.3f", p$effect_true))
  if (p$effect_type == "ratio") {
    cat(sprintf("  (log = %+.3f)\n", log(p$effect_true)))
  } else {
    cat("\n")
  }
  cat(sprintf("  Event rate (ctrl):  %.1f%%\n", 100 * p$p_event))
  if (p$effect_type == "ratio" && !is.na(p$p_treat)) {
    cat(sprintf("  Event rate (trt):   %.1f%%  (implied by OR = %.2f)\n",
        100 * p$p_treat, p$effect_true))
  }
  cat(sprintf("  Min subgroup (n):   %d  (~%d per arm)\n",
      p$n_min, round(p$n_min * p$treat_alloc)))
  cat(sprintf("  SE(log effect):     %.3f\n\n", x$se_log_effect))

  cat(sprintf("  Confounders:        %d continuous (×%d cuts) + %d binary\n",
      p$n_confounders_cont, p$cont_cutoff, p$n_confounders_bin))
  cat(sprintf("  Single candidates:  %d\n", x$n_single))
  cat(sprintf("  Pair candidates:    %d  (maxk = %d)\n", x$n_pairs, p$maxk))
  cat(sprintf("  Total candidates:   %d\n", x$n_candidates))
  if (p$rho > 0) {
    cat(sprintf("  Effective tests:    %.0f  (rho = %.2f)\n",
        x$n_eff, p$rho))
  }

  cat("\n  Threshold   P(per candidate)   FPR (indep)   FPR (adj)\n")
  cat(  "  ---------   ----------------   -----------   ---------\n")

  for (i in seq_along(p$thresholds)) {
    thr <- p$thresholds[i]
    nm  <- sprintf("%.2f", thr)
    marker <- if (abs(thr - p$c1) < 1e-6) " <-- c1" else ""
    cat(sprintf("    %.2f        %5.1f%%             %5.1f%%        %5.1f%%%s\n",
        thr,
        100 * x$per_candidate[nm],
        100 * x$fpr_indep[nm],
        100 * x$fpr_adjusted[nm],
        marker))
  }

  cat("\n  Note: Screening-stage only.  Consistency check (pconsistency)\n")
  cat("  and bootstrap correction further reduce the operational FPR.\n")
  cat("  Independence FPR is an upper bound; adjusted uses Sidak with\n")
  cat(sprintf("  effective tests (rho = %.2f).\n", p$rho))
  cat("============================================================\n")
}
