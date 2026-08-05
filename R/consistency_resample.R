# =============================================================================
# Resampling approximation to the consistency (splitting) calculation
# -----------------------------------------------------------------------------
# Package-integration version for forestsearch/R/. Covers the Cox (survival)
# path AND the GLM paths (binary OR/RR/RD, continuous MD, count IRR), unadjusted
# or covariate-adjusted.
#
# The repeated 50/50 split + two refits in run_single_consistency_split() is
# replaced by a SINGLE subgroup fit plus a multiplier representation of the
# split pair. The per-subject treatment influence contributions are the
# treatment column of the model's dfbeta; the random-split half effects are
#       beta_hat +/- D,   D = sum_i G_i * dfbeta[i, treat],
# with G_i any mean-zero, unit-variance multiplier (Rademacher = the 50/50
# split, normal = Lin 1997, centred Poisson = Dobler-Beyersmann-Pauly). The
# split-perturbation SD is sigma_D = sqrt(sum_i dfbeta[i, treat]^2), the robust
# (sandwich) SE of the treatment coefficient. A split is consistent iff both
# halves clear the threshold on the COMPARISON scale c:
#       rate ~ 2 * Phi((beta_hat - c) / sigma_D) - 1,
#   c = log(threshold) for ratio measures (HR, OR, RR, IRR),
#   c = threshold      for identity measures (RD, IRD, MD).
# beta_hat is the treatment coefficient of the same model the effect estimator
# fits: coxph -> log-HR; logistic -> log-OR; log-binomial/Poisson -> log-RR;
# identity-binomial -> RD; lm -> MD; Poisson + offset(log time) -> log-IRR.
#
# Dependency note: survival is used only on the Cox path; the GLM paths use base
# stats (glm/lm/poisson + dfbeta). Adjusted paths call the forestsearch internal
# .fs_adjust_terms() when available (standalone fallback otherwise).
# =============================================================================


# -----------------------------------------------------------------------------
# Cox (survival) pieces
# -----------------------------------------------------------------------------

#' Per-subject Cox pieces for the consistency approximation
#'
#' Fits the treatment-only (or covariate-adjusted) Cox model once on a subgroup
#' and returns the treatment log hazard ratio, the per-subject treatment
#' `dfbeta` influence contributions, and the split-perturbation SD.
#'
#' @param df data.frame/data.table with the survival columns (and any columns
#'   referenced by `adjust_covariates`).
#' @param tte.name,event.name,treat.name Character; time, event, and (0/1)
#'   treatment column names (defaults `"Y"`, `"Event"`, `"Treat"`).
#' @param cox_init Numeric; warm-start for the unadjusted Cox coefficient.
#' @param adjust_covariates Character vector or `NULL`; terms appended after
#'   `treat.name`.
#'
#' @return List with `beta_hat`, `dfbeta`, `sigma_D`, `n`, `d`, `log_scale`
#'   (TRUE), `measure` (`"HR"`); or `NULL` if degenerate.
#'
#' @keywords internal
#' @importFrom survival coxph Surv
#' @importFrom stats residuals coef reformulate
.consistency_cox_pieces <- function(df, tte.name = "Y", event.name = "Event",
                                    treat.name = "Treat", cox_init = 0,
                                    adjust_covariates = NULL) {
  df <- as.data.frame(df)
  if (!all(c(tte.name, event.name, treat.name) %in% names(df))) {
    stop("df must contain columns: ", tte.name, ", ", event.name, ", ", treat.name)
  }
  if (nrow(df) < 2L || sum(df[[event.name]]) < 2L) return(NULL)

  adj <- .consistency_adj_terms(adjust_covariates)
  resp <- sprintf("survival::Surv(%s, %s)", tte.name, event.name)
  fmla <- stats::reformulate(c(treat.name, adj), response = resp)

  fit_args <- list(formula = fmla, data = df,
                   robust = FALSE, model = FALSE, x = FALSE, y = FALSE)
  if (length(adj) == 0L) fit_args$init <- cox_init

  fit <- tryCatch(suppressWarnings(do.call(survival::coxph, fit_args)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- stats::coef(fit)
  if (length(cf) == 0L || !(treat.name %in% names(cf)) || any(is.na(cf))) return(NULL)

  beta_hat <- unname(cf[[treat.name]])
  dfbeta <- stats::residuals(fit, type = "dfbeta")
  if (is.matrix(dfbeta)) {
    j <- match(treat.name, colnames(dfbeta))
    if (is.na(j)) return(NULL)
    dfbeta <- dfbeta[, j]
  }
  dfbeta <- as.numeric(dfbeta)
  sigma_D <- sqrt(sum(dfbeta^2))
  if (!is.finite(beta_hat) || !is.finite(sigma_D) || sigma_D <= 0) return(NULL)

  list(beta_hat = beta_hat, dfbeta = dfbeta, sigma_D = sigma_D,
       n = nrow(df), d = sum(df[[event.name]]), log_scale = TRUE, measure = "HR")
}


# -----------------------------------------------------------------------------
# GLM pieces (binary OR/RR/RD, continuous MD, count IRR)
#' One-step treatment `dfbeta` for a GLM fit
#'
#' Replaces \code{stats::dfbeta(fit)} on the GLM influence path. The
#' influence is built directly from the one-step formula rather than taken
#' from \code{lm.influence()}:
#'
#' \deqn{\mathrm{dfbeta}_i = (X'WX)^{-1} x_i w_i r_i / (1 - h_i)}
#'
#' with \eqn{w} the IRLS working weights, \eqn{r} the working residuals
#' (\code{fit$residuals}, which is the working residual for every GLM family,
#' so this is general rather than family-specific) and \eqn{h} the hat values.
#'
#' \strong{Why not \code{stats::dfbeta()}.} R 4.6.0 corrected
#' \code{stats::dfbeta()} for \code{glm} objects:
#' \code{weighted.residuals()} previously returned \emph{deviance} residuals,
#' which -- per R's own NEWS -- "do not give the one-step approximated dfbeta
#' values". \code{sigma_D} computed under R < 4.6.0 therefore differs from the
#' same quantity computed under R >= 4.6.0. Building the influence here makes
#' it correct on every R version rather than documenting the dependency.
#' Validated against exact leave-one-out refits in
#' \code{dev/identifier-alignment/dfbeta_glm_validation.qmd}.
#'
#' \strong{Behaviour at leverage 1, and why it differs from
#' \code{stats::dfbeta()}.} When \eqn{h_i = 1} -- one observation entirely
#' determines the coefficient, e.g. a single subject in one arm -- the
#' one-step is \strong{undefined}: \eqn{0/0}, returned as \code{NaN}. That is
#' the honest answer, because the leave-one-out estimate it approximates does
#' not exist either; removing that observation leaves the coefficient
#' unidentifiable.
#'
#' \code{stats::dfbeta()} instead returns \strong{zero} influence for that
#' observation. Zero is wrong, and wrong in the direction that matters: the
#' observation which alone determines the estimate is recorded as having no
#' influence on it, so \code{sigma_D} comes out \emph{smaller} than it would
#' with two such observations. Measured on a probit fit: one subject in the
#' treated arm gives \code{sigma_D} 0.447, two give 1.529. A subgroup whose
#' effect rests on a single observation is thereby scored with a falsely
#' narrow interval and a falsely low \code{t_g} -- anticonservative exactly
#' where caution is warranted.
#'
#' No epsilon guard is applied. Clamping \eqn{1 - h_i} would make
#' \code{sigma_D} scale as \eqn{1/\epsilon}, letting an arbitrary constant set
#' the magnitude of a reported quantity. \code{NaN} propagates to
#' \code{sigma_D}, the caller's \code{is.finite()} check rejects it, and the
#' candidate is dropped -- the same treatment every other non-estimable
#' candidate receives.
#'
#' @param fit A fitted \code{glm} (or \code{lm}) object.
#' @return Numeric matrix, one row per observation and one column per
#'   coefficient.
#' @noRd
.dfbeta_glm <- function(fit) {
  w <- fit$weights
  h <- stats::lm.influence(fit, do.coef = FALSE)$hat
  X <- stats::model.matrix(fit)
  A <- solve(crossprod(X, w * X))
  r <- sqrt(w) * fit$residuals          # weighted working residual
  t(A %*% t(X * (sqrt(w) * r / (1 - h))))
}


# -----------------------------------------------------------------------------

#' Normalise adjustment terms (mirrors forestsearch:::.fs_adjust_terms)
#' @keywords internal
.consistency_adj_terms <- function(adjust_covariates) {
  if (exists(".fs_adjust_terms", mode = "function")) {
    return(.fs_adjust_terms(adjust_covariates))
  }
  if (is.null(adjust_covariates)) return(character(0))
  tt <- trimws(as.character(adjust_covariates))
  tt[nzchar(tt)]
}

#' Per-subject GLM pieces for the consistency approximation
#'
#' Fits the same treatment-effect model the corresponding
#' [make_effect_estimator()] closure fits (matching `effect_measure` and
#' `adverse_outcome` conventions), and returns the treatment coefficient, its
#' per-subject `dfbeta` influence contributions, and the split-perturbation SD.
#'
#' Supported: `"OR"` (logistic), `"RR"` (log-binomial with Poisson fallback),
#' `"RD"` (identity-link binomial), `"MD"` (OLS), `"IRR"` (Poisson with
#' `offset(log(time))`). `"IRD"` and propensity-adjusted (`iptw`/`dr_gcomp`)
#' effects are not single coefficients, so they return `NULL` (the caller falls
#' back to splitting).
#'
#' @param df data.frame/data.table with `outcome.name`, `treat.name`, any
#'   `adjust_covariates`, and `offset.name` (rate measures).
#' @param outcome_type Character; `"binary"`, `"continuous"`, or `"count"`.
#' @param effect_measure Character or `NULL`; resolved to the outcome-type
#'   default when `NULL` (binary `"RD"`, continuous `"MD"`, count `"IRR"`).
#' @param treat.name,outcome.name,offset.name Character; column names.
#' @param adjust_covariates Character vector or `NULL`.
#' @param adverse_outcome Logical; if `FALSE`, binary `Y` is flipped to `1 - Y`
#'   and continuous `Y` negated, matching the estimator so that a larger
#'   coefficient consistently means treatment harm.
#'
#' @return List with `beta_hat`, `dfbeta`, `sigma_D`, `n`, `d`, `log_scale`,
#'   `measure`; or `NULL` if unsupported/degenerate/non-convergent.
#'
#' @keywords internal
#' @importFrom stats glm lm binomial poisson coef dfbeta reformulate
.consistency_glm_pieces <- function(df, outcome_type, effect_measure = NULL,
                                    treat.name = "Treat", outcome.name = "Y",
                                    offset.name = NULL, adjust_covariates = NULL,
                                    adverse_outcome = TRUE) {
  df <- as.data.frame(df)
  if (is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      binary = "RD", continuous = "MD", count = "IRR")
  }
  if (effect_measure == "IRD") return(NULL)   # rate difference: not a coefficient
  if (!all(c(outcome.name, treat.name) %in% names(df))) {
    stop("df must contain columns: ", outcome.name, ", ", treat.name)
  }

  # adverse-outcome flip (match make_effect_estimator)
  if (!adverse_outcome) {
    if (outcome_type == "continuous") {
      df[[outcome.name]] <- -df[[outcome.name]]
    } else if (outcome_type == "binary") {
      df[[outcome.name]] <- 1L - df[[outcome.name]]
    }
  }

  z <- df[[treat.name]]
  n0 <- sum(z == 0L, na.rm = TRUE); n1 <- sum(z == 1L, na.rm = TRUE)
  if (n0 < 3L || n1 < 3L) return(NULL)

  adj <- .consistency_adj_terms(adjust_covariates)
  rhs <- c(treat.name, adj)

  fit <- tryCatch(suppressWarnings(switch(effect_measure,
    OR = stats::glm(stats::reformulate(rhs, outcome.name), data = df,
                    family = stats::binomial("logit")),
    RR = {
      f <- stats::reformulate(rhs, outcome.name)
      p0 <- mean(df[[outcome.name]], na.rm = TRUE)
      g <- tryCatch(stats::glm(f, data = df, family = stats::binomial("log"),
                               start = c(log(max(p0, 0.01)), rep(0, length(rhs)))),
                    error = function(e) NULL, warning = function(w) NULL)
      if (is.null(g) || isFALSE(g$converged)) {
        g <- stats::glm(f, data = df, family = stats::poisson("log"))  # modified Poisson
      }
      g
    },
    RD = {
      y0 <- df[[outcome.name]][z == 0]; y1 <- df[[outcome.name]][z == 1]
      stats::glm(stats::reformulate(rhs, outcome.name), data = df,
                 family = stats::binomial("identity"),
                 start = c(mean(y0, na.rm = TRUE),
                           mean(y1, na.rm = TRUE) - mean(y0, na.rm = TRUE),
                           rep(0, length(adj))))
    },
    MD = stats::lm(stats::reformulate(rhs, outcome.name), data = df),
    IRR = {
      if (is.null(offset.name)) stop("IRR requires offset.name")
      df$.log_offset <- log(pmax(df[[offset.name]], .Machine$double.eps))
      stats::glm(stats::reformulate(c(rhs, "offset(.log_offset)"), outcome.name),
                 data = df, family = stats::poisson("log"))
    },
    return(NULL)
  )), error = function(e) NULL)

  if (is.null(fit)) return(NULL)
  cf <- stats::coef(fit)
  if (!(treat.name %in% names(cf)) || any(is.na(cf))) return(NULL)
  if (!is.null(fit$converged) && isFALSE(fit$converged)) return(NULL)

  beta_hat <- unname(cf[[treat.name]])
  # .dfbeta_glm(), not stats::dfbeta(): the latter changed meaning at R 4.6.0
  # and returns a wrong (zero) influence at leverage 1.  See .dfbeta_glm().
  dfb <- tryCatch(.dfbeta_glm(fit), error = function(e) NULL)
  if (is.null(dfb)) return(NULL)
  if (is.matrix(dfb)) {
    j <- match(treat.name, colnames(dfb))
    if (is.na(j)) return(NULL)
    dfb <- dfb[, j]
  }
  dfb <- as.numeric(dfb)
  sigma_D <- sqrt(sum(dfb^2))
  if (!is.finite(beta_hat) || !is.finite(sigma_D) || sigma_D <= 0) return(NULL)

  list(beta_hat = beta_hat, dfbeta = dfb, sigma_D = sigma_D,
       n = nrow(df), d = sum(df[[outcome.name]], na.rm = TRUE),
       log_scale = effect_measure %in% c("OR", "RR", "IRR"),
       measure = effect_measure)
}


# -----------------------------------------------------------------------------
# Public entry point (Cox + GLM)
# -----------------------------------------------------------------------------

#' Resampling approximation to the splitting consistency rate
#'
#' Approximates the forestsearch consistency rate for a candidate subgroup
#' without repeated sample splits and refits. The subgroup is fit once; the
#' random-split pair of half-sample treatment effects is represented as
#' `beta_hat +/- D`, a multiplier sum of the per-subject treatment `dfbeta`
#' contributions (see the file header). Handles the Cox (survival) path and the
#' GLM paths (binary OR/RR/RD, continuous MD, count IRR).
#'
#' `multiplier = "rademacher"` (default) reproduces the 50/50 Bernoulli split of
#' [run_single_consistency_split()]; `"normal"` is the Lin (1997) Gaussian
#' multiplier and `"poisson"` a centred unit-Poisson (DBP) multiplier.
#'
#' @param df data.frame/data.table with the relevant outcome columns.
#' @param hr.consistency Numeric; survival HR consistency threshold (also the
#'   default GLM threshold when `consistency_threshold` is `NULL`). Default 1.0.
#' @param consistency_threshold Numeric or `NULL`; GLM consistency threshold on
#'   the natural scale (OR/RR/IRR, or RD/MD/IRD). `NULL` uses `hr.consistency`.
#' @param comparison_threshold Numeric or `NULL`; the consistency threshold
#'   already resolved to the comparison scale `c` (`log(threshold)` for ratio
#'   measures, identity for identity measures). When supplied it is used
#'   directly and the internal transform is bypassed. This is the value
#'   [forestsearch()] threads through for GLM outcomes (it resolves the scale
#'   upstream), so callers avoid a double log-transform.
#' @param outcome_type Character; `"survival"` (default), `"binary"`,
#'   `"continuous"`, or `"count"`.
#' @param effect_measure Character or `NULL`; GLM effect measure (resolved to
#'   the outcome-type default when `NULL`). Ignored for survival.
#' @param method Character; `"closed"`, `"mc"`, or `"both"` (default).
#' @param multiplier Character; `"rademacher"` (default), `"normal"`,
#'   `"poisson"`.
#' @param draws Integer; multiplier draws for the Monte Carlo path.
#' @param tte.name,event.name,treat.name Character; survival/treatment column
#'   names.
#' @param outcome.name,offset.name Character; GLM outcome and (rate) offset
#'   column names.
#' @param cox_init Numeric; Cox warm-start.
#' @param adjust_covariates Character vector or `NULL`; adjustment terms.
#' @param adverse_outcome Logical; GLM adverse-outcome convention (see
#'   [.consistency_glm_pieces()]).
#' @param seed Integer or `NULL`; RNG seed for the Monte Carlo draws.
#'
#' @return List with `rate_closed`, `rate_mc`, `beta_hat`, `sigma_D`, `n`, `d`,
#'   `delta`, `measure`; all-`NA` rates (with available summaries) if the fit
#'   fails or the configuration is unsupported.
#'
#' @seealso [run_single_consistency_split()], [consistency_resample_compare()]
#' @examples
#' \dontrun{
#' # survival (Cox)
#' consistency_resample(sub, hr.consistency = 1.0)
#' # binary odds ratio
#' consistency_resample(sub, outcome_type = "binary", effect_measure = "OR",
#'                      consistency_threshold = 1.0,
#'                      outcome.name = "Y", treat.name = "Treat")
#' # count incidence-rate ratio
#' consistency_resample(sub, outcome_type = "count", effect_measure = "IRR",
#'                      consistency_threshold = 1.0,
#'                      outcome.name = "events", offset.name = "fu_time")
#' }
#' @importFrom stats pnorm rbinom rnorm rpois
#' @export
consistency_resample <- function(df, hr.consistency = 1.0,
                                 consistency_threshold = NULL,
                                 comparison_threshold = NULL,
                                 outcome_type = c("survival", "binary",
                                                  "continuous", "count"),
                                 effect_measure = NULL,
                                 method = c("both", "closed", "mc"),
                                 multiplier = c("rademacher", "normal", "poisson"),
                                 draws = 1000L,
                                 tte.name = "Y", event.name = "Event",
                                 treat.name = "Treat",
                                 outcome.name = "Y", offset.name = NULL,
                                 cox_init = 0, adjust_covariates = NULL,
                                 adverse_outcome = TRUE, seed = NULL) {
  outcome_type <- match.arg(outcome_type)
  method <- match.arg(method)
  multiplier <- match.arg(multiplier)

  out <- list(rate_closed = NA_real_, rate_mc = NA_real_,
              beta_hat = NA_real_, sigma_D = NA_real_,
              n = NA_integer_, d = NA_integer_, delta = NA_real_,
              measure = NA_character_)

  if (outcome_type == "survival") {
    pieces  <- .consistency_cox_pieces(df, tte.name, event.name, treat.name,
                                       cox_init, adjust_covariates)
    thr_nat <- hr.consistency
  } else {
    pieces  <- .consistency_glm_pieces(df, outcome_type, effect_measure,
                                       treat.name, outcome.name, offset.name,
                                       adjust_covariates, adverse_outcome)
    thr_nat <- if (!is.null(consistency_threshold)) consistency_threshold
               else hr.consistency
  }
  if (is.null(pieces)) return(out)

  c_cmp <- if (!is.null(comparison_threshold)) {
    comparison_threshold                 # caller-supplied, already on comparison scale
  } else if (pieces$log_scale) {
    log(thr_nat)                         # ratio measure: natural -> log
  } else {
    thr_nat                              # identity measure
  }
  delta <- pieces$beta_hat - c_cmp
  out$beta_hat <- pieces$beta_hat
  out$sigma_D  <- pieces$sigma_D
  out$n        <- pieces$n
  out$d        <- pieces$d
  out$delta    <- delta
  out$measure  <- pieces$measure

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
    D <- as.numeric(G %*% pieces$dfbeta)
    out$rate_mc <- mean((pieces$beta_hat - abs(D)) > c_cmp)
  }

  out
}


#' Validate the consistency approximation against literal splitting (Cox)
#'
#' Calls the package's [run_single_consistency_split()] `R_true` times and
#' compares the Monte Carlo rate to the closed-form and multiplier
#' approximations from [consistency_resample()]. Survival path only; for the
#' GLM paths use the standalone GLM validator.
#'
#' @inheritParams consistency_resample
#' @param R_true Integer; number of literal random splits used for the truth.
#'
#' @return One-row data.frame with `n`, `d`, `beta_hat`, `sigma_D`,
#'   `rate_true`, `valid_true`, `rate_closed`, `rate_mc`, `err_closed`,
#'   `err_mc`.
#' @seealso [consistency_resample()], [run_single_consistency_split()]
#' @export
consistency_resample_compare <- function(df, hr.consistency = 1.0,
                                         R_true = 400L, draws = 1000L,
                                         multiplier = "rademacher",
                                         tte.name = "Y", event.name = "Event",
                                         treat.name = "Treat", cox_init = 0,
                                         adjust_covariates = NULL, seed = NULL) {
  if (!exists("run_single_consistency_split", mode = "function")) {
    stop("run_single_consistency_split() not found; load forestsearch first.")
  }

  df <- as.data.frame(df)
  if (tte.name   != "Y")     df[["Y"]]     <- df[[tte.name]]
  if (event.name != "Event") df[["Event"]] <- df[[event.name]]
  if (treat.name != "Treat") df[["Treat"]] <- df[[treat.name]]
  N.x <- nrow(df)

  if (!is.null(seed)) set.seed(seed)
  res <- vapply(seq_len(R_true), function(r) {
    run_single_consistency_split(
      df.x = df, N.x = N.x, hr.consistency = hr.consistency,
      cox_init = cox_init, adjust_covariates = adjust_covariates)
  }, numeric(1))
  valid <- sum(!is.na(res))
  rate_true <- if (valid > 0L) mean(res, na.rm = TRUE) else NA_real_

  approx <- consistency_resample(
    df, hr.consistency = hr.consistency, method = "both",
    multiplier = multiplier, draws = draws,
    tte.name = "Y", event.name = "Event", treat.name = "Treat",
    cox_init = cox_init, adjust_covariates = adjust_covariates, seed = seed)

  data.frame(
    n = approx$n, d = approx$d, beta_hat = approx$beta_hat,
    sigma_D = approx$sigma_D, rate_true = rate_true, valid_true = valid,
    rate_closed = approx$rate_closed, rate_mc = approx$rate_mc,
    err_closed = approx$rate_closed - rate_true,
    err_mc = approx$rate_mc - rate_true)
}
