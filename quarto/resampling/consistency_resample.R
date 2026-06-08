# =============================================================================
# Resampling approximation to the consistency (splitting) calculation
# -----------------------------------------------------------------------------
# Un-adjusted Cox path: Surv(Y, Event) ~ Treat (the forestsearch default).
#
# The repeated 50/50 split + two Cox refits in run_single_consistency_split()
# is replaced by a SINGLE subgroup fit plus a score-residual multiplier
# representation of the split pair.
#
# Theory (Leon et al. 2024, Sec 2.1; Dobler-Beyersmann-Pauly 2017):
#   The half-split log-HRs are, to first order, hatbeta_s +/- D where
#       D = sum_i G_i * dfbeta_i,
#   dfbeta_i = residuals(fit, type = "dfbeta") are the per-subject Cox
#   influence contributions (c_i / I_s), and G_i is any mean-zero,
#   unit-variance multiplier:
#       Rademacher (G_i = 2*B_i - 1)  <=> the actual 50/50 Bernoulli split
#       Normal                        <=> Lin (1997) Gaussian multiplier
#       Centred unit Poisson          <=> the DBP "weird"-type multiplier
#   The split perturbation SD is sigma_D = sqrt(sum_i dfbeta_i^2), i.e. the
#   robust (sandwich) SE of the subgroup log-HR. A split is "consistent with
#   harm" iff both halves clear the threshold, i.e. |D| < hatbeta_s - log(c2),
#   giving the closed form 2 * Phi((hatbeta_s - log(c2)) / sigma_D) - 1.
#
# Dependency note: this unadjusted path is intentionally survival-only and does
# NOT use weightedsurv; it mirrors get_split_hr_fast()'s toolchain exactly.
# =============================================================================


#' Per-subject Cox pieces for the consistency approximation (unadjusted)
#'
#' Fits the unadjusted treatment-only Cox model once on a subgroup and returns
#' the quantities needed to approximate the splitting consistency rate without
#' refitting: the log hazard ratio, the per-subject `dfbeta` influence
#' contributions, and the implied split-perturbation standard deviation.
#'
#' @param df data.frame or data.table with the survival columns.
#' @param tte.name,event.name,treat.name Character; column names for the
#'   time-to-event, event indicator, and (0/1) treatment. Defaults `"Y"`,
#'   `"Event"`, `"Treat"` match the forestsearch subgroup data contract.
#' @param cox_init Numeric; warm-start for the single Cox coefficient
#'   (default 0), matching `get_split_hr_fast()`.
#'
#' @return A list with elements `beta_hat` (log HR), `dfbeta` (numeric vector of
#'   per-subject influence values), `sigma_D` (robust SE of `beta_hat`, equal to
#'   the split-perturbation SD), `n`, and `d` (events); or `NULL` if the fit
#'   fails or is degenerate.
#'
#' @keywords internal
#' @importFrom survival coxph Surv
#' @importFrom stats residuals coef
.consistency_cox_pieces <- function(df, tte.name = "Y", event.name = "Event",
                                    treat.name = "Treat", cox_init = 0) {
  df <- as.data.frame(df)
  y <- df[[tte.name]]
  d <- df[[event.name]]
  z <- df[[treat.name]]
  if (is.null(y) || is.null(d) || is.null(z)) {
    stop("df must contain columns: ", tte.name, ", ", event.name, ", ", treat.name)
  }
  if (length(y) < 2L || sum(d) < 2L) {
    return(NULL)
  }

  fit <- tryCatch(
    suppressWarnings(survival::coxph(
      survival::Surv(y, d) ~ z,
      init = cox_init, robust = FALSE, model = FALSE, x = FALSE, y = FALSE
    )),
    error = function(e) NULL
  )
  if (is.null(fit) || length(stats::coef(fit)) == 0L || any(is.na(stats::coef(fit)))) {
    return(NULL)
  }

  beta_hat <- unname(stats::coef(fit)[1L])
  # dfbeta_i = c_i / I_s : per-subject influence on the log-HR.
  dfbeta <- as.numeric(stats::residuals(fit, type = "dfbeta"))
  sigma_D <- sqrt(sum(dfbeta^2))           # robust (sandwich) SE of beta_hat
  if (!is.finite(beta_hat) || !is.finite(sigma_D) || sigma_D <= 0) {
    return(NULL)
  }

  list(beta_hat = beta_hat, dfbeta = dfbeta, sigma_D = sigma_D,
       n = length(y), d = sum(d))
}


#' Resampling approximation to the splitting consistency rate (unadjusted Cox)
#'
#' Approximates the forestsearch consistency rate for a candidate subgroup
#' without performing repeated sample splits and Cox refits. The subgroup is fit
#' once; the random-split pair of half-sample log-HRs is represented as
#' `beta_hat +/- D`, where `D` is a multiplier sum of the per-subject Cox
#' influence (`dfbeta`) contributions. See the file header for the underlying
#' representation.
#'
#' The default `multiplier = "rademacher"` reproduces the exact 50/50 Bernoulli
#' split mechanism used by [run_single_consistency_split()]. `"normal"` gives the
#' Lin (1997) Gaussian multiplier and `"poisson"` a centred unit-Poisson
#' multiplier (a Dobler-Beyersmann-Pauly "weird"-type choice); both are offered
#' for experimentation with small-sample multiplier accuracy.
#'
#' @param df data.frame or data.table with the survival columns.
#' @param hr.consistency Numeric; HR consistency threshold `c2`
#'   (default 1.0). A split counts as consistent when both half log-HRs exceed
#'   `log(hr.consistency)`, matching the strict comparison in
#'   [run_single_consistency_split()].
#' @param method Character; `"closed"` (normal-limit closed form), `"mc"`
#'   (multiplier Monte Carlo), or `"both"` (default).
#' @param multiplier Character; multiplier distribution for the Monte Carlo
#'   path: `"rademacher"` (default), `"normal"`, or `"poisson"`.
#' @param draws Integer; number of multiplier draws for the Monte Carlo path
#'   (default 1000).
#' @param tte.name,event.name,treat.name Character; survival column names
#'   (defaults `"Y"`, `"Event"`, `"Treat"`).
#' @param cox_init Numeric; warm-start for the Cox coefficient (default 0).
#' @param seed Integer or `NULL`; if non-`NULL`, sets the RNG seed for the
#'   Monte Carlo draws (default `NULL`).
#'
#' @return A list with `rate_closed`, `rate_mc` (whichever are requested; the
#'   other is `NA`), and the fit summaries `beta_hat`, `sigma_D`, `n`, `d`, and
#'   `delta` (`beta_hat - log(hr.consistency)`). Returns all-`NA` rates with the
#'   available summaries when the subgroup fit fails.
#'
#' @seealso [run_single_consistency_split()], [consistency_resample_compare()]
#' @examples
#' \dontrun{
#' library(survival)
#' df <- within(gbsg, {
#'   Y <- rfstime / 30.4375; Event <- status; Treat <- hormon
#' })
#' sub <- df[df$er <= 0, c("Y", "Event", "Treat")]   # paper's H-hat subgroup
#' consistency_resample(sub, hr.consistency = 1.0)
#' }
#' @importFrom stats pnorm rbinom rnorm rpois
#' @export
consistency_resample <- function(df, hr.consistency = 1.0,
                                 method = c("both", "closed", "mc"),
                                 multiplier = c("rademacher", "normal", "poisson"),
                                 draws = 1000L,
                                 tte.name = "Y", event.name = "Event",
                                 treat.name = "Treat", cox_init = 0,
                                 seed = NULL) {
  method <- match.arg(method)
  multiplier <- match.arg(multiplier)

  pieces <- .consistency_cox_pieces(df, tte.name, event.name, treat.name, cox_init)
  out <- list(rate_closed = NA_real_, rate_mc = NA_real_,
              beta_hat = NA_real_, sigma_D = NA_real_,
              n = NA_integer_, d = NA_integer_, delta = NA_real_)
  if (is.null(pieces)) {
    return(out)
  }

  thr <- log(hr.consistency)
  delta <- pieces$beta_hat - thr
  out$beta_hat <- pieces$beta_hat
  out$sigma_D <- pieces$sigma_D
  out$n <- pieces$n
  out$d <- pieces$d
  out$delta <- delta

  if (method %in% c("closed", "both")) {
    out$rate_closed <- max(0, 2 * stats::pnorm(delta / pieces$sigma_D) - 1)
  }

  if (method %in% c("mc", "both")) {
    if (!is.null(seed)) set.seed(seed)
    n <- pieces$n
    G <- switch(multiplier,
      rademacher = matrix(2 * stats::rbinom(draws * n, 1, 0.5) - 1, nrow = draws),
      normal     = matrix(stats::rnorm(draws * n), nrow = draws),
      poisson    = matrix(stats::rpois(draws * n, 1) - 1, nrow = draws)
    )
    # D_r = sum_i G_ri * dfbeta_i ; halves are beta_hat +/- D_r
    D <- as.numeric(G %*% pieces$dfbeta)
    out$rate_mc <- mean((pieces$beta_hat - abs(D)) > thr)
  }

  out
}


#' Validate the consistency approximation against literal splitting
#'
#' Runs the true repeated-split consistency calculation (a faithful inline copy
#' of [run_single_consistency_split()]: independent Bernoulli(1/2) split, each
#' half requiring at least `min_rows` rows and `min_events` events, both half
#' hazard ratios strictly exceeding `hr.consistency`) and compares the resulting
#' Monte Carlo rate to the closed-form and multiplier approximations from
#' [consistency_resample()].
#'
#' @param df data.frame or data.table with the survival columns.
#' @param hr.consistency Numeric; HR consistency threshold (default 1.0).
#' @param R_true Integer; number of true random splits (default 400).
#' @param draws Integer; multiplier draws for the approximation (default 1000).
#' @param multiplier Character; multiplier for the approximation
#'   (default `"rademacher"`, matching the true split mechanism).
#' @param min_rows,min_events Integer; per-half validity guards mirroring
#'   [run_single_consistency_split()] (defaults 5 and 2).
#' @param tte.name,event.name,treat.name Character; survival column names.
#' @param cox_init Numeric; Cox warm-start (default 0).
#' @param seed Integer or `NULL`; RNG seed applied before the true splits and
#'   reused for the approximation draws (default `NULL`).
#'
#' @return A one-row data.frame: `n`, `d`, `beta_hat`, `sigma_D`,
#'   `rate_true`, `valid_true` (number of non-discarded splits), `rate_closed`,
#'   `rate_mc`, and the signed errors `err_closed`, `err_mc`.
#'
#' @seealso [consistency_resample()], [run_single_consistency_split()]
#' @examples
#' \dontrun{
#' library(survival)
#' df <- within(gbsg, {
#'   Y <- rfstime / 30.4375; Event <- status; Treat <- hormon
#' })
#' sub <- df[df$er <= 0, c("Y", "Event", "Treat")]
#' consistency_resample_compare(sub, hr.consistency = 1.0, R_true = 400, seed = 1)
#' }
#' @importFrom survival coxph Surv
#' @importFrom stats rbinom coef
#' @export
consistency_resample_compare <- function(df, hr.consistency = 1.0,
                                         R_true = 400L, draws = 1000L,
                                         multiplier = "rademacher",
                                         min_rows = 5L, min_events = 2L,
                                         tte.name = "Y", event.name = "Event",
                                         treat.name = "Treat", cox_init = 0,
                                         seed = NULL) {
  df <- as.data.frame(df)
  y <- df[[tte.name]]; d <- df[[event.name]]; z <- df[[treat.name]]
  n <- length(y)
  thr <- log(hr.consistency)

  # --- fast unadjusted half-HR, faithful to get_split_hr_fast() ---
  half_hr <- function(idx) {
    yy <- y[idx]; dd <- d[idx]; zz <- z[idx]
    if (length(yy) < min_rows || sum(dd) < min_events) return(NA_real_)
    fit <- tryCatch(
      suppressWarnings(survival::coxph(
        survival::Surv(yy, dd) ~ zz,
        init = cox_init, robust = FALSE, model = FALSE, x = FALSE, y = FALSE
      )),
      error = function(e) NULL
    )
    if (is.null(fit) || length(stats::coef(fit)) == 0L || any(is.na(stats::coef(fit)))) {
      return(NA_real_)
    }
    exp(unname(stats::coef(fit)[1L]))
  }

  # --- true splitting loop (mirrors run_single_consistency_split) ---
  if (!is.null(seed)) set.seed(seed)
  n_success <- 0L
  n_valid <- 0L
  for (r in seq_len(R_true)) {
    inA <- stats::rbinom(n, 1, 0.5) == 1L
    if (sum(inA) < min_rows || sum(!inA) < min_rows) next
    hrA <- half_hr(which(inA))
    hrB <- half_hr(which(!inA))
    if (is.na(hrA) || is.na(hrB)) next
    n_valid <- n_valid + 1L
    if (hrA > hr.consistency && hrB > hr.consistency) n_success <- n_success + 1L
  }
  rate_true <- if (n_valid > 0L) n_success / n_valid else NA_real_

  # --- approximation (reuse the same seed stream for the multiplier draws) ---
  approx <- consistency_resample(
    df, hr.consistency = hr.consistency, method = "both",
    multiplier = multiplier, draws = draws,
    tte.name = tte.name, event.name = event.name, treat.name = treat.name,
    cox_init = cox_init, seed = seed
  )

  data.frame(
    n          = approx$n,
    d          = approx$d,
    beta_hat   = approx$beta_hat,
    sigma_D    = approx$sigma_D,
    rate_true  = rate_true,
    valid_true = n_valid,
    rate_closed = approx$rate_closed,
    rate_mc     = approx$rate_mc,
    err_closed  = approx$rate_closed - rate_true,
    err_mc      = approx$rate_mc - rate_true
  )
}
