#' Generate an Easy-Case Synthetic Benchmark Dataset for ForestSearch CV
#'
#' Creates a survival dataset with a known harm subgroup for benchmarking
#' ForestSearch subgroup identification and its cross-validation
#' reproducibility under \strong{ideal conditions}: three binary
#' covariates (one true harm driver, two prognostic-only factors),
#' large sample size, and no noise covariates.  This is a deliberately
#' easy test --- the algorithm should not need any pre-screening
#' (LASSO, GRF) to identify the correct subgroup, and cross-validation
#' should recover it on nearly every fold.
#'
#' @section Data-Generating Mechanism:
#' Let \eqn{W_i} be the treatment indicator (RCT, 1:1 randomization) and
#' \eqn{z_{1i}} the true harm-driver indicator.  Event times follow a
#' Weibull distribution with hazard
#'   \deqn{\lambda(t \mid W, z_1, z_2, z_3) = \lambda_0(t) \cdot
#'         \exp(\beta_W W + \beta_{W z_1} W z_1
#'              + \gamma_2 z_2 + \gamma_3 z_3),}
#' where
#' \itemize{
#'   \item \eqn{\beta_W = \log(\mathrm{hr\_benefit})} --- treatment
#'         log-HR in the benefit group (\eqn{z_1 = 0}).
#'   \item \eqn{\beta_{W z_1} = \log(\mathrm{hr\_harm}) -
#'         \log(\mathrm{hr\_benefit})} --- treatment-by-subgroup
#'         interaction that flips the sign of the treatment effect
#'         inside the harm subgroup (\eqn{z_1 = 1}).
#'   \item \eqn{\gamma_2 = \log(\mathrm{prog\_hr\_z2})} and
#'         \eqn{\gamma_3 = \log(\mathrm{prog\_hr\_z3})} ---
#'         prognostic effects.  These shift the baseline hazard but
#'         \strong{do not} interact with treatment.
#' }
#' With defaults \code{hr_benefit = 0.65} and \code{hr_harm = 2.20}:
#' \itemize{
#'   \item Subjects with \code{z1 = 0} (~80\% of population): treatment
#'         HR \eqn{\approx} 0.65.
#'   \item Subjects with \code{z1 = 1} (~20\% of population): treatment
#'         HR \eqn{\approx} 2.20.
#' }
#' Censoring is administrative at \code{admin_censor_time}.  Defaults
#' yield roughly 60\% overall event rate.
#'
#' @section Expected ForestSearch Performance:
#' This is an \strong{easy} benchmark by design:
#' \itemize{
#'   \item Three binary candidate covariates, no cut grid, no
#'         compound false-positive strips to resist.
#'   \item \code{z2}, \code{z3} are prognostic-only --- they shift
#'         baseline hazard but do not interact with treatment, so
#'         their subgroup HRs should track the ITT-level HR, not
#'         \code{z1}'s subgroup HR.
#'   \item N = 700 yields ~140 subjects and ~100+ events in the
#'         harm subgroup --- well above Peduzzi's \eqn{\geq 10}
#'         events-per-covariate heuristic.
#' }
#' With \code{use_lasso = FALSE}, \code{use_grf = FALSE}, \code{maxk = 2}
#' a correctly-functioning implementation should give:
#' \itemize{
#'   \item Full-data: identify \code{"\{z1\}"} with ~100\% reliability.
#'   \item 10-fold CV: \code{Exact Match} \eqn{\geq} 0.95;
#'         \code{sens_H}, \code{PPV_H} \eqn{\geq} 0.95.
#'   \item LOO CV: \code{Exact Match} \eqn{\geq} 0.98.
#'   \item KfA_0 HR matches FA_0 HR within ratio [0.85, 1.18].
#' }
#' Performance substantially below these thresholds indicates a real
#' algorithmic problem, not a dataset-difficulty issue.
#'
#' @section Expected event counts per arm (defaults):
#' Harm subgroup (\eqn{n \approx 140}):
#' \itemize{
#'   \item Treatment (\eqn{n \approx 70}): P(event) \eqn{\approx} 0.90
#'         \eqn{\Rightarrow} ~63 events.
#'   \item Control (\eqn{n \approx 70}): P(event) \eqn{\approx} 0.65
#'         \eqn{\Rightarrow} ~45 events.
#' }
#' Both exceed the Peduzzi \eqn{\geq 10}-event heuristic for stable
#' Cox estimation.
#'
#' @param n Integer. Sample size.  Default \code{700L}.
#' @param prev_z1 Numeric in (0, 1).  Prevalence of the true harm
#'   driver \code{z1 = 1}.  Default \code{0.20}.
#' @param prev_z2 Numeric in (0, 1).  Prevalence of the first
#'   prognostic factor.  Default \code{0.40}.
#' @param prev_z3 Numeric in (0, 1).  Prevalence of the second
#'   prognostic factor.  Default \code{0.50}.
#' @param hr_benefit Numeric > 0.  Treatment hazard ratio in the
#'   benefit population (\code{z1 = 0}).  Default \code{0.65}.
#' @param hr_harm Numeric > 0.  Treatment hazard ratio in the harm
#'   subgroup (\code{z1 = 1}).  Default \code{2.20}.
#' @param prog_hr_z2 Numeric > 0.  Prognostic hazard ratio for
#'   \code{z2} (baseline-hazard multiplier; NOT treatment-interacting).
#'   Default \code{1.4}.
#' @param prog_hr_z3 Numeric > 0.  Prognostic hazard ratio for
#'   \code{z3}.  Default \code{0.8}.
#' @param shape Numeric > 0.  Weibull shape parameter; \code{shape = 1}
#'   gives exponential times.  Default \code{1.0}.
#' @param baseline_scale Numeric > 0.  Weibull scale for untreated
#'   (\code{W = 0}, all \code{z = 0}) hazard.  Default \code{60}.
#' @param admin_censor_time Numeric > 0.  Time of administrative
#'   censoring.  Default \code{60}.
#' @param seed Integer.  RNG seed for reproducibility.  Default
#'   \code{20250421L}.
#'
#' @return A data frame with \code{n} rows and the columns:
#' \describe{
#'   \item{\code{id}}{Subject identifier (1..n).}
#'   \item{\code{time}}{Observed time = min(event time, admin censor).}
#'   \item{\code{status}}{Event indicator (1 = event, 0 = censored).}
#'   \item{\code{treat}}{Treatment indicator (1 = treated, 0 = control).}
#'   \item{\code{z1}}{Binary.  The TRUE harm driver (ground truth).}
#'   \item{\code{z2}, \code{z3}}{Binary prognostic-only factors.
#'         Shift baseline hazard but do not interact with treatment.}
#'   \item{\code{true_harm}}{Integer indicator equal to
#'         \code{as.integer(z1 == 1)}.  Validation only ---
#'         \strong{do not} pass to \code{forestsearch()}.}
#'   \item{\code{true_log_hr}}{Realized treatment-arm log-HR for each
#'         subject.  Validation only.}
#' }
#'
#' @examples
#' \dontrun{
#' df <- make_benchmark_dataset(n = 700L)
#'
#' candidate_names <- c("z1", "z2", "z3")
#'
#' fs <- forestsearch::forestsearch(
#'   df.analysis      = df,
#'   confounders.name = candidate_names,
#'   outcome.name     = "time",
#'   event.name       = "status",
#'   treat.name       = "treat",
#'   id.name          = "id",
#'   is.RCT           = TRUE,
#'   use_lasso        = FALSE,
#'   use_grf          = FALSE,
#'   sg_focus         = "maxSG",
#'   use_twostage     = FALSE,
#'   quiet            = TRUE
#' )
#'
#' # Expect: paste(fs$sg.harm, collapse = " & ") == "{z1}"
#' }
#'
#' @export
make_benchmark_dataset <- function(
    n                 = 700L,
    prev_z1           = 0.20,
    prev_z2           = 0.40,
    prev_z3           = 0.50,
    hr_benefit        = 0.65,
    hr_harm           = 2.20,
    prog_hr_z2        = 1.4,
    prog_hr_z3        = 0.8,
    shape             = 1.0,
    baseline_scale    = 60,
    admin_censor_time = 60,
    seed              = 20250421L
) {

  # ---- Input validation -------------------------------------------------
  stopifnot(
    is.numeric(n),                    length(n) == 1L,
    n >= 100L,                        n == as.integer(n),
    is.numeric(prev_z1),              length(prev_z1) == 1L,
    prev_z1 > 0,                      prev_z1 < 1,
    is.numeric(prev_z2),              length(prev_z2) == 1L,
    prev_z2 > 0,                      prev_z2 < 1,
    is.numeric(prev_z3),              length(prev_z3) == 1L,
    prev_z3 > 0,                      prev_z3 < 1,
    is.numeric(hr_benefit),           length(hr_benefit) == 1L,
    hr_benefit > 0,
    is.numeric(hr_harm),              length(hr_harm) == 1L,
    hr_harm > 0,
    is.numeric(prog_hr_z2),           length(prog_hr_z2) == 1L,
    prog_hr_z2 > 0,
    is.numeric(prog_hr_z3),           length(prog_hr_z3) == 1L,
    prog_hr_z3 > 0,
    is.numeric(shape),                length(shape) == 1L,
    shape > 0,
    is.numeric(baseline_scale),       length(baseline_scale) == 1L,
    baseline_scale > 0,
    is.numeric(admin_censor_time),    length(admin_censor_time) == 1L,
    admin_censor_time > 0,
    is.numeric(seed),                 length(seed) == 1L
  )

  set.seed(as.integer(seed))
  n <- as.integer(n)

  # ---- Covariates (all binary; z2, z3 are prognostic-only noise) --------
  z1 <- rbinom(n, 1L, prev_z1)   # TRUE harm driver (treatment interaction)
  z2 <- rbinom(n, 1L, prev_z2)   # Prognostic-only (no treatment interaction)
  z3 <- rbinom(n, 1L, prev_z3)   # Prognostic-only (no treatment interaction)

  # ---- Treatment (RCT, 1:1) --------------------------------------------
  treat <- rbinom(n, 1L, 0.5)

  # ---- Weibull event-time generation ------------------------------------
  # Hazard multiplier =
  #   exp( beta_W * W + beta_WZ1 * W * z1 + gamma_2 * z2 + gamma_3 * z3 )
  # Prognostic terms depend on z2, z3 only (no treatment interaction).
  beta_W   <- log(hr_benefit)
  beta_WZ1 <- log(hr_harm) - log(hr_benefit)
  gamma_2  <- log(prog_hr_z2)
  gamma_3  <- log(prog_hr_z3)

  eta <- beta_W * treat +
         beta_WZ1 * treat * z1 +
         gamma_2 * z2 +
         gamma_3 * z3

  # Inverse-CDF sampling: T = scale * (-log U / exp(eta))^(1/shape)
  U       <- runif(n)
  T_event <- baseline_scale * (-log(U) / exp(eta))^(1 / shape)

  # ---- Administrative censoring ----------------------------------------
  time   <- pmin(T_event, admin_censor_time)
  status <- as.integer(T_event <= admin_censor_time)

  data.frame(
    id          = seq_len(n),
    time        = time,
    status      = status,
    treat       = treat,
    z1          = z1,
    z2          = z2,
    z3          = z3,
    true_harm   = as.integer(z1 == 1L),
    true_log_hr = eta
  )
}


#' Summarize an Easy-Case Synthetic Benchmark Dataset
#'
#' Computes empirical event counts, event rates, and per-subgroup Cox
#' hazard-ratio estimates from a dataset produced by
#' \code{\link{make_benchmark_dataset}}.  Useful as a quick sanity
#' check before running ForestSearch: verifies that the realized DGM
#' matches its nominal specification.
#'
#' @param df Data frame from \code{make_benchmark_dataset()}.
#' @return A list with:
#' \describe{
#'   \item{\code{n}}{Total sample size.}
#'   \item{\code{harm_prevalence}}{Empirical \code{mean(z1)}.}
#'   \item{\code{event_rate}}{Overall fraction of events.}
#'   \item{\code{events_by_arm}}{Matrix of event counts by treatment
#'         arm \eqn{\times} harm subgroup.}
#'   \item{\code{hr_overall}}{Cox HR, full data.}
#'   \item{\code{hr_harm}}{Cox HR within the true harm subgroup
#'         (\code{z1 = 1}).}
#'   \item{\code{hr_benefit}}{Cox HR within the true benefit group
#'         (\code{z1 = 0}).}
#' }
#' @export
summarize_benchmark_dataset <- function(df) {
  stopifnot(
    is.data.frame(df),
    all(c("time", "status", "treat", "z1") %in% names(df))
  )

  events_by_arm <- with(
    df,
    tapply(
      status, list(treat = treat, harm = z1), sum
    )
  )

  cox_safe <- function(subset_df) {
    if (nrow(subset_df) < 10L || sum(subset_df$status) < 5L) {
      return(NA_real_)
    }
    fit <- tryCatch(
      survival::coxph(
        survival::Surv(time, status) ~ treat,
        data = subset_df, robust = FALSE,
        model = FALSE, x = FALSE, y = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) NA_real_ else exp(unname(fit$coefficients[1]))
  }

  list(
    n               = nrow(df),
    harm_prevalence = mean(df$z1),
    event_rate      = mean(df$status),
    events_by_arm   = events_by_arm,
    hr_overall      = cox_safe(df),
    hr_harm         = cox_safe(df[df$z1 == 1L, , drop = FALSE]),
    hr_benefit      = cox_safe(df[df$z1 == 0L, , drop = FALSE])
  )
}
