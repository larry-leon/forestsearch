#' Fit Cox Model for Subgroup
#'
#' Fits a Cox model for a subgroup and returns estimate and standard error.
#'
#' @param df_sg Data frame for subgroup.
#' @param cox.formula Cox model formula.
#' @param est.loghr Logical. Is estimate on log(HR) scale?
#' @return List with estimate and standard error.
#' @importFrom survival coxph
#' @export
get_Cox_sg <- function(df_sg, cox.formula, est.loghr = TRUE, cox_initial = log(1)) {

  # Validate inputs (keep your existing validation)
  names_tocheck <- all.vars(cox.formula)
  check <- unlist(lapply(names_tocheck, grep, names(df_sg), value = TRUE))
  check2 <- match(names_tocheck, check)
  if (sum(!is.na(check2)) != length(names_tocheck)) {
    stop("df_sg dataset NOT contain cox.formula variables")
  }

  # Fit model
  fit <- suppressWarnings(
    coxph(cox.formula, data = df_sg, model = FALSE, x = FALSE, y = FALSE, robust = TRUE, init = cox_initial)
  )

  # OPTIMIZATION: Call summary() once, cache result
  fit_sum <- summary(fit)
  coef_matrix <- fit_sum$coefficients

  # Extract coefficients
  bhat <- coef_matrix[, "coef"]

  # Extract appropriate SE
  if (est.loghr) {
    est_obs <- bhat
    se_obs <- coef_matrix[, "robust se"]
  } else {
    est_obs <- exp(bhat)
    se_obs <- exp(bhat) * coef_matrix[, "robust se"]
  }

  return(list(est_obs = est_obs, se_obs = se_obs))
}




get_Cox_sg_legacy <- function(df_sg, cox.formula, est.loghr = TRUE, cox_initial = log(1)) {
  names_tocheck <- all.vars(cox.formula)
  check <- unlist(lapply(names_tocheck, grep, names(df_sg), value = TRUE))
  check2 <- match(names_tocheck, check)
  if (sum(!is.na(check2)) != length(names_tocheck)) stop("df_sg dataset NOT contain cox.formula variables")
  # Fit Cox model with robust standard errors

  fit <- summary(coxph(cox.formula, data = df_sg, robust = TRUE, init = cox_initial))$coefficients

  fit <- summary(fit)$coefficients

  # log(hr) parameters
  if (est.loghr) {
    bhat <- c(fit[, "coef"])
    est_obs <- bhat
    se_obs <- c(fit[, "robust se"])
  }
  # Otherwise, hr
  if (!est.loghr) {
    bhat <- c(fit[, "coef"])
    est_obs <- exp(bhat)
    sebhat <- c(fit[, "robust se"])
    se_obs <- est_obs * sebhat
  }
  return(list(est_obs = est_obs, se_obs = se_obs))
}



#' Count ID Occurrences in Bootstrap Sample
#'
#' Counts the number of times an ID appears in a bootstrap sample.
#'
#' @param x ID value.
#' @param dfb Data frame of bootstrap sample.
#' @return Integer count of occurrences.
#' @export

count.id <- function(x,dfb){
  sum(dfb$id == x)
  }

#' Calculate Covariance for Bootstrap Estimates
#'
#' Calculates the covariance between a vector and bootstrap estimates.
#'
#' @param x Numeric vector.
#' @param Est Numeric vector of bootstrap estimates.
#' @return Numeric value of covariance.
#' @export

calc_cov <- function(x,Est){
  mean(c((x-mean(x,na.rm=TRUE)) * Est),na.rm=TRUE)
}

#' Confidence Interval for Estimate
#'
#' Calculates confidence interval for an estimate, optionally on log(HR) scale.
#'
#' @param x Numeric estimate.
#' @param sd Numeric standard deviation.
#' @param alpha Numeric significance level (default: 0.025).
#' @param scale Character. "hr" or "1/hr".
#' @param est.loghr Logical. Is estimate on log(HR) scale?
#' @return List with length, lower, upper, sd, and estimate.
#' @importFrom stats qnorm
#' @export

ci_est <- function(x, sd, alpha = 0.025, scale = "hr", est.loghr = TRUE) {
  # Input validation
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || is.infinite(x)) {
    stop("'x' must be a single, finite numeric value.")
  }
  if (!is.numeric(sd) || length(sd) != 1 || is.na(sd) || is.infinite(sd) || sd < 0) {
    stop("'sd' must be a single, non-negative, finite numeric value.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 0.5) {
    stop("'alpha' must be a single numeric value between 0 and 0.5.")
  }
  if (!scale %in% c("hr", "1/hr")) {
    stop("'scale' must be either 'hr' or '1/hr'.")
  }
  if (!is.logical(est.loghr) || length(est.loghr) != 1 || is.na(est.loghr)) {
    stop("'est.loghr' must be a single logical value.")
  }

  c_alpha <- qnorm(1 - alpha)
  c_low <- x - c_alpha * sd
  c_up <- x + c_alpha * sd
  est <- x
  new_low <- c_low
  new_up <- c_up
  out_sd <- sd

  if (scale == "hr" && est.loghr) {
    est <- exp(x)
    out_sd <- exp(x) * sd
    new_low <- exp(c_low)
    new_up <- exp(c_up)
  }
  if (scale == "1/hr" && est.loghr) {
    est <- exp(-x)
    out_sd <- exp(-x) * sd
    new_low <- exp(-c_up)
    new_up <- exp(-c_low)
  }
  length <- new_up - new_low

  # Return as named list
  return(list(
    length = length,
    lower = new_low,
    upper = new_up,
    sd = out_sd,
    est = est
  ))
}

#' Bootstrap Confidence Interval and Bias Correction Results
#'
#' Calculates confidence intervals and bias-corrected estimates for bootstrap results.
#'
#' @param Hobs Numeric. Observed estimate.
#' @param seHobs Numeric. Standard error of observed estimate.
#' @param H1_adj Numeric. Bias-corrected estimate 1.
#' @param H2_adj Numeric. Bias-corrected estimate 2 (optional).
#' @param ystar Matrix of bootstrap samples.
#' @param cov_method Character. Covariance method ("standard" or "nocorrect").
#' @param cov_trim Numeric. Trimming proportion for covariance (default: 0.0).
#' @param est.scale Character. "hr" or "1/hr".
#' @param est.loghr Logical. Is estimate on log(HR) scale?
#' @return Data.table with confidence intervals and estimates.
#' @importFrom data.table data.table
#' @export

get_dfRes <- function(Hobs, seHobs, H1_adj, H2_adj = NULL, ystar, cov_method = "standard", cov_trim = 0.0, est.scale = "hr", est.loghr = TRUE) {
  # Check required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  # Check required function
  if (!exists("get_targetEst")) {
    stop("Function 'get_targetEst' must be defined in the environment.")
  }
  # Input validation
  if (!is.numeric(Hobs) || length(Hobs) != 1 || is.na(Hobs) || is.infinite(Hobs)) {
    stop("'Hobs' must be a single, finite numeric value.")
  }
  if (!is.numeric(seHobs) || length(seHobs) != 1 || is.na(seHobs) || is.infinite(seHobs) || seHobs < 0) {
    stop("'seHobs' must be a single, non-negative, finite numeric value.")
  }
  if (missing(ystar)) {
    stop("'ystar' must be provided.")
  }
  if (!is.character(cov_method) || length(cov_method) != 1) {
    stop("'cov_method' must be a single character string.")
  }
  if (!is.numeric(cov_trim) || length(cov_trim) != 1 || is.na(cov_trim) || cov_trim < 0 || cov_trim > 1) {
    stop("'cov_trim' must be a single numeric value between 0 and 1.")
  }
  if (!est.scale %in% c("hr", "1/hr")) {
    stop("'est.scale' must be either 'hr' or '1/hr'.")
  }
  if (!is.logical(est.loghr) || length(est.loghr) != 1 || is.na(est.loghr)) {
    stop("'est.loghr' must be a single logical value.")
  }

  # Un-adjusted
  cest <- ci_est(x = Hobs, sd = seHobs, scale = est.scale, est.loghr = est.loghr)
  H0_lower <- cest$lower
  H0_upper <- cest$upper
  sdH0 <- cest$sd
  H0 <- cest$est

  est <- get_targetEst(x = H1_adj, ystar = ystar, cov_method = cov_method, cov_trim = cov_trim)
  q <- est$target_est
  se_new <- est$sehat_new

  cest <- ci_est(x = q, sd = se_new, scale = est.scale, est.loghr = est.loghr)
  H1_lower <- cest$lower
  H1_upper <- cest$upper
  sdH1 <- cest$sd
  H1 <- cest$est

  outres <- data.table::data.table(H0, sdH0, H0_lower, H0_upper, H1, sdH1, H1_lower, H1_upper)

  if (!is.null(H2_adj)) {
    est <- get_targetEst(x = H2_adj, ystar = ystar, cov_method = cov_method, cov_trim = cov_trim)
    q <- est$target_est
    se_new <- est$sehat_new
    cest <- ci_est(x = q, sd = se_new, scale = est.scale, est.loghr = est.loghr)
    H2_lower <- cest$lower
    H2_upper <- cest$upper
    sdH2 <- cest$sd
    H2 <- cest$est
    outres <- data.table::data.table(H0, sdH0, H0_lower, H0_upper,
                                     H1, sdH1, H1_lower, H1_upper,
                                     H2, sdH2, H2_lower, H2_upper)
  }
  return(outres)
}

#' Target Estimate and Standard Error for Bootstrap
#'
#' Calculates target estimate and standard error for bootstrap samples.
#'
#' @param x Numeric vector of estimates.
#' @param ystar Matrix of bootstrap samples.
#' @param cov_method Character. Covariance method ("standard" or "nocorrect").
#' @param cov_trim Numeric. Trimming proportion for covariance (default: 0.0).
#' @return List with target estimate, standard errors, and correction term.
#' @export

get_targetEst <- function(x, ystar, cov_method = "standard", cov_trim = 0.0) {
  mx <- mean(x, na.rm = TRUE)
  xc <- c(x - mx)
  N <- ncol(ystar)
  # Ystar_mat is (B x N) matrix
  # Across columns of Ystar
  if (cov_method == "standard" | cov_method == "nocorrect") {
    cov_i <- apply(ystar, 2, calc_cov, Est = xc)
  }
  # Implement correction
  if (cov_method != "nocorrect") {
    varhat <- N * mean(cov_i^2)
    seH <- sqrt(varhat)
    # Denominator in cov_i is B.eval:
    B.eval <- sum(!is.na(xc))
    # correction
    nb_ratio <- N / (B.eval^2)
    termc <- nb_ratio * B.eval * mean(xc^2, na.rm = TRUE, trim = cov_trim)
    if (varhat <= termc) {
      varhat_new <- varhat
      seH_new <- sqrt(varhat)
    }
    if (varhat > termc) {
      varhat_new <- varhat - termc
      seH_new <- sqrt(varhat_new)
    }
  }
  if (cov_method == "nocorrect") {
    varhat <- N * mean(cov_i^2, trim = cov_trim)
    seH <- sqrt(varhat)
    termc <- 0.0
    varhat_new <- varhat
    seH_new <- seH
  }
  out <- (list(target_est = mx, sehat = seH, sehat_new = seH_new, term_correct = termc, varhat = varhat_new))
  return(out)
}



