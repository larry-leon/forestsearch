#' Count ID Occurrences in Bootstrap Sample
#'
#' Counts the number of times an ID appears in a bootstrap sample.
#'
#' @param x ID value.
#' @param dfb Data frame of bootstrap sample.
#' @return Integer count of occurrences.
#' @export
count.id <- function(x, dfb) {
  # Input validation
  if (is.null(x)) stop("'x' cannot be NULL")
  if (!is.data.frame(dfb)) stop("'dfb' must be a data.frame")
  if (!"id" %in% names(dfb)) stop("'dfb' must contain column 'id'")
  
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
calc_cov <- function(x, Est) {
  # Input validation
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (!is.numeric(Est)) stop("'Est' must be numeric")
  if (length(x) != length(Est)) stop("'x' and 'Est' must have same length")
  
  mean(c((x - mean(x, na.rm = TRUE)) * Est), na.rm = TRUE)
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
get_dfRes <- function(Hobs, seHobs, H1_adj, H2_adj = NULL, ystar, 
                      cov_method = "standard", cov_trim = 0.0, 
                      est.scale = "hr", est.loghr = TRUE) {
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
  if (!is.numeric(seHobs) || length(seHobs) != 1 || is.na(seHobs) || 
      is.infinite(seHobs) || seHobs < 0) {
    stop("'seHobs' must be a single, non-negative, finite numeric value.")
  }
  if (missing(ystar) || !is.matrix(ystar)) {
    stop("'ystar' must be provided as a matrix.")
  }
  if (!is.character(cov_method) || length(cov_method) != 1) {
    stop("'cov_method' must be a single character string.")
  }
  if (!is.numeric(cov_trim) || length(cov_trim) != 1 || is.na(cov_trim) || 
      cov_trim < 0 || cov_trim > 1) {
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
  # Input validation
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (!is.matrix(ystar)) stop("'ystar' must be a matrix")
  if (!is.character(cov_method)) stop("'cov_method' must be character")
  if (!is.numeric(cov_trim) || cov_trim < 0 || cov_trim > 1) {
    stop("'cov_trim' must be numeric between 0 and 1")
  }
  
  mx <- mean(x, na.rm = TRUE)
  xc <- c(x - mx)
  N <- ncol(ystar)
  
  # Ystar_mat is (B x N) matrix
  # Across columns of Ystar
  if (cov_method == "standard" || cov_method == "nocorrect") {
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
  out <- (list(target_est = mx, sehat = seH, sehat_new = seH_new, 
               term_correct = termc, varhat = varhat_new))
  return(out)
}

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
get_Cox_sg <- function(df_sg, cox.formula, est.loghr = TRUE) {
  # Input validation
  if (!is.data.frame(df_sg)) stop("'df_sg' must be a data.frame")
  if (!inherits(cox.formula, "formula")) stop("'cox.formula' must be a formula")
  if (!is.logical(est.loghr)) stop("'est.loghr' must be logical")
  
  names_tocheck <- all.vars(cox.formula)
  check <- unlist(lapply(names_tocheck, grep, names(df_sg), value = TRUE))
  check2 <- match(names_tocheck, check)
  if (sum(!is.na(check2)) != length(names_tocheck)) {
    stop("df_sg dataset does NOT contain cox.formula variables")
  }
  
  # Fit Cox model with robust standard errors
  fit <- summary(coxph(cox.formula, data = df_sg, robust = TRUE))$coefficients
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

#' Coverage Indicator for Confidence Interval
#'
#' Checks if a target value is covered by a confidence interval.
#'
#' @param lower Numeric lower bound.
#' @param upper Numeric upper bound.
#' @param target Numeric target value (default: 0).
#' @return List with coverage indicator and interval length.
#' @export
ci_cover <- function(lower, upper, target = 0) {
  # Input validation
  if (!is.numeric(lower)) stop("'lower' must be numeric")
  if (!is.numeric(upper)) stop("'upper' must be numeric")
  if (!is.numeric(target)) stop("'target' must be numeric")
  
  if (length(target) == 1) {
    cover <- ifelse(lower <= target & upper >= target, 1, 0)
    LC <- upper - lower
    return(list(cover = cover, LC = LC))
  }
  if (length(target) > 1) {
    LC <- upper - lower
    covers <- NULL
    for (tt in 1:length(target)) {
      cover <- ifelse(lower <= target[tt] & upper >= target[tt], 1, 0)
      covers <- c(covers, cover)
    }
    return(list(cover = covers, LC = LC))
  }
}

#' Bootstrap Confidence Intervals for Two Targets
#'
#' Calculates confidence intervals for two bootstrap targets.
#'
#' @param Q1 Numeric vector of first target estimates.
#' @param Q2 Numeric vector of second target estimates.
#' @param ystar Matrix of bootstrap samples.
#' @return List with estimates, standard errors, and confidence intervals.
#' @export
getCIs <- function(Q1, Q2, ystar) {
  # Input validation
  if (!is.numeric(Q1)) stop("'Q1' must be numeric")
  if (!is.numeric(Q2)) stop("'Q2' must be numeric")
  if (!is.matrix(ystar)) stop("'ystar' must be a matrix")
  
  # Target 1
  est <- get_targetEst(x = Q1, ystar = ystar)
  # If est.loghr=TRUE then on log(HR) scale
  q1 <- est$target_est
  se1_new <- est$sehat_new
  se1 <- est$sehat
  rm("est")
  # use SE new
  cest <- ci_est(x = q1, sd = se1_new)
  H1_lower <- cest$lower
  H1_upper <- cest$upper
  rm("cest")

  # Target 2
  est <- get_targetEst(x = Q2, ystar = ystar)
  # Call this H.bc
  q2 <- est$target_est
  se2 <- est$sehat
  se2_new <- est$sehat_new
  rm("est")
  cest <- ci_est(x = q2, sd = se2_new)
  H2_lower <- cest$lower
  H2_upper <- cest$upper
  rm("cest")
  return(list(q1 = q1, se1 = se1, se1_new = se1_new, H1_lower = H1_lower, H1_upper = H1_upper,
              q2 = q2, se2 = se2, se2_new = se2_new, H2_lower = H2_lower, H2_upper = H2_upper))
}

#' Prepare Data for Bias Plot
#'
#' Prepares a data frame for plotting bias in bootstrap estimates.
#'
#' @param res Data frame of results.
#' @param dgm Data-generating mechanism (truth) for simulation.
#' @return Data frame for bias plot.
#' @importFrom data.table data.table
#' @export
get_dfPlot <- function(res, dgm) {
  # Input validation
  if (!is.data.frame(res)) stop("'res' must be a data.frame")
  if (!is.list(dgm)) stop("'dgm' must be a list")
  if (!all(c("hr.H.true", "hr.Hc.true") %in% names(dgm))) {
    stop("'dgm' must contain 'hr.H.true' and 'hr.Hc.true'")
  }
  
  hrH_true <- dgm$hr.H.true
  hrHc_true <- dgm$hr.Hc.true
  if (est.loghr && est.scale == "loghr") {
    hrH_true <- log(dgm$hr.H.true)
    hrHc_true <- log(dgm$hr.Hc.true)
  }
  res_new <- within(res, {
    b1H_1 <- 100 * (H1.bc - H_true) / H_true
    b1H_2 <- 100 * (H1.bc - hatH_causal) / hatH_causal
    b1H_3 <- 100 * (H1.bc - hrH_true) / hrH_true

    b2H_1 <- 100 * (H2.bc - H_true) / H_true
    b2H_2 <- 100 * (H2.bc - hatH_causal) / hatH_causal
    b2H_3 <- 100 * (H2.bc - hrH_true) / hrH_true

    b0H_1 <- 100 * (H_obs - H_true) / H_true
    b0H_2 <- 100 * (H_obs - hatH_causal) / hatH_causal
    b0H_3 <- 100 * (H_obs - hrH_true) / hrH_true
  })

  df_bc <- NULL
  # hr(H_true) target
  res_new <- data.table(res_new)
  df_res <- res_new[, c("b1H_1")]
  df_res$est <- "BC_1"
  df_res$target <- "H_true"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)

  df_res <- res_new[, c("b2H_1")]
  df_res$est <- "BC_2"
  df_res$target <- "H_true"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)

  df_res <- res_new[, c("b0H_1")]
  df_res$est <- "Obs"
  df_res$target <- "H_true"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)

  # hatH_causal target
  df_res <- res_new[, c("b1H_2")]
  df_res$est <- "BC_1"
  df_res$target <- "H_po"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)

  df_res <- res_new[, c("b2H_2")]
  df_res$est <- "BC_2"
  df_res$target <- "H_po"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)

  df_res <- res_new[, c("b0H_2")]
  df_res$est <- "Obs"
  df_res$target <- "H_po"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)

  # H_causal fixed
  df_res <- res_new[, c("b1H_3")]
  df_res$est <- "BC_1"
  df_res$target <- "H_fixed"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)

  df_res <- res_new[, c("b2H_3")]
  df_res$est <- "BC_2"
  df_res$target <- "H_fixed"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)

  df_res <- res_new[, c("b0H_3")]
  df_res$est <- "Obs"
  df_res$target <- "H_fixed"
  names(df_res)[1] <- "bias"
  df_bc <- rbind(df_bc, df_res)
  return(df_bc)
}

#' Variance Summary for Bootstrap Results
#'
#' Summarizes variance and coverage for bootstrap results.
#'
#' @param res Data frame of results.
#' @return Data frame with summary statistics.
#' @export
var_summary <- function(res) {
  # Input validation
  if (!is.data.frame(res)) stop("'res' must be a data.frame")
  
  df_var <- NULL
  # BC1 SD and Avg(est(sd))
  aa <- with(res, sqrt(var(H1.bc, na.rm = TRUE)))
  bb <- with(res, mean(seH1.bc, na.rm = TRUE))
  cc <- with(res, mean(H1.bc, na.rm = TRUE))
  dd <- with(res, mean(H1_cover1, na.rm = TRUE))
  ee <- with(res, mean(H1_cover2, na.rm = TRUE))
  ff <- with(res, mean(H1_cover3, na.rm = TRUE))
  gg <- with(res, mean(H1_length, na.rm = TRUE))
  df_var <- rbind(df_var, c(cc, aa, bb, dd, ee, ff, gg))

  aa <- with(res, sqrt(var(H2.bc, na.rm = TRUE)))
  bb <- with(res, mean(seH2.bc, na.rm = TRUE))
  cc <- with(res, mean(H2.bc, na.rm = TRUE))
  dd <- with(res, mean(H2_cover1, na.rm = TRUE))
  ee <- with(res, mean(H2_cover2, na.rm = TRUE))
  ff <- with(res, mean(H2_cover3, na.rm = TRUE))
  gg <- with(res, mean(H2_length, na.rm = TRUE))
  df_var <- rbind(df_var, c(cc, aa, bb, dd, ee, ff, gg))

  aa <- with(res, sqrt(var(Hc1.bc, na.rm = TRUE)))
  bb <- with(res, mean(seHc1.bc, na.rm = TRUE))
  cc <- with(res, mean(Hc1.bc, na.rm = TRUE))
  dd <- with(res, mean(Hc1_cover1, na.rm = TRUE))
  ee <- with(res, mean(Hc1_cover2, na.rm = TRUE))
  ff <- with(res, mean(Hc1_cover3, na.rm = TRUE))
  gg <- with(res, mean(Hc1_length, na.rm = TRUE))
  df_var <- rbind(df_var, c(cc, aa, bb, dd, ee, ff, gg))

  aa <- with(res, sqrt(var(Hc2.bc, na.rm = TRUE)))
  bb <- with(res, mean(seHc2.bc, na.rm = TRUE))
  cc <- with(res, mean(Hc2.bc, na.rm = TRUE))
  dd <- with(res, mean(Hc2_cover1, na.rm = TRUE))
  ee <- with(res, mean(Hc2_cover2, na.rm = TRUE))
  ff <- with(res, mean(Hc2_cover3, na.rm = TRUE))
  gg <- with(res, mean(Hc2_length, na.rm = TRUE))
  df_var <- rbind(df_var, c(cc, aa, bb, dd, ee, ff, gg))

  colnames(df_var) <- c("Est", "SD", "Est(SD)", "C1", "C2", "C3", "L")
  rownames(df_var) <- c("H1_bc", "H2_bc", "Hc1_bc", "Hc2_bc")
  return(round(df_var, digits = 3))
}

#' Confidence Interval for Cox Model Estimate
#'
#' Calculates confidence interval and coverage for Cox model estimate.
#'
#' @param df Data frame.
#' @param est Character. Name of estimate column.
#' @param se Character. Name of standard error column.
#' @param target Numeric or character. Target value or column name.
#' @param alpha Numeric significance level (default: 0.025).
#' @param digits Integer. Number of digits for rounding (default: 3).
#' @return List with lower bound, upper bound, and coverage indicator.
#' @export
getci_Cox <- function(df, est, se, target, alpha = 0.025, digits = 3) {
  # Input validation
  if (!is.data.frame(df)) stop("'df' must be a data.frame")
  if (!is.character(est) || !est %in% names(df)) {
    stop("'est' must be a column name in 'df'")
  }
  if (!is.character(se) || !se %in% names(df)) {
    stop("'se' must be a column name in 'df'")
  }
  
  a <- df[est]
  b <- df[se]
  if (!is.numeric(target)) {
    if (!target %in% names(df)) stop("'target' must be numeric or a column name in 'df'")
    c <- df[target]
  }
  if (is.numeric(target)) {
    c <- target
  }
  log_lb <- log(a) - qnorm(0.975) * (b / a)
  log_ub <- log(a) + qnorm(0.975) * (b / a)
  lb <- exp(log_lb)
  ub <- exp(log_ub)
  cov <- ifelse(c >= lb & c <= ub, 1, 0)
  return(list(lb = round(lb, 3), ub = round(ub, 3), cover = cov))
}

#' Summary Statistic for Data Frame Column
#'
#' Returns mean and standard deviation for a column, formatted as a string.
#'
#' @param df Data frame.
#' @param name Character. Column name.
#' @param sigdig Integer. Number of significant digits (default: 2).
#' @param includeSD Logical. Include "SD=" in output (default: FALSE).
#' @param showSD Logical. Show SD in output (default: TRUE).
#' @return Character string with mean and SD.
#' @export
SummaryStat <- function(df, name, sigdig = 2, includeSD = FALSE, showSD = TRUE) {
  # Input validation
  if (!is.data.frame(df)) stop("'df' must be a data.frame")
  if (!is.character(name) || !name %in% names(df)) {
    stop("'name' must be a column name in 'df'")
  }
  
  if (is.data.table(df)) {
    df <- as.data.frame(df)
  }
  # include "SD=" in output
  if (!includeSD) {
    temp <- na.omit(df[, c(name)])
    m <- round(mean(temp), sigdig)
    sig <- round(sqrt(var(temp)), sigdig)
    # Check if binary
    if (all(as.numeric(temp) %in% c(0, 1))) {
      sig <- round(sqrt(m * (1 - m) / length(temp)), sigdig)
    }
    if (sigdig <= 4 && m <= 0.001) m <- "<0.001"
    if (sigdig <= 4 && sig <= 0.001) sig <- "<0.001"
    out <- paste0(m, " (")
    out <- paste0(out, sig)
    out <- paste0(out, ")")
  }
  if (includeSD) {
    temp <- na.omit(df[, c(name)])
    m <- round(mean(temp), sigdig)
    sig <- round(sqrt(var(temp)), sigdig)
    # Check if binary
    if (all(as.numeric(temp) %in% c(0, 1))) {
      sig <- round(sqrt(m * (1 - m) / length(m)), sigdig)
    }
    if (sigdig <= 4 && m <= 0.001) m <- "<0.001"
    if (sigdig <= 4 && sig <= 0.001) sig <- "<0.001"
    out <- paste0(m, " (SD=")
    out <- paste0(out, sig)
    out <- paste0(out, ")")
  }
  if (!showSD) out <- c(m)
  return(out)
}

#' Difference in Rate Between Two Data Frames
#'
#' Calculates the difference in mean rate for a column between two data frames.
#'
#' @param df1 First data frame.
#' @param df2 Second data frame.
#' @param name Character. Column name.
#' @param sigdig Integer. Number of significant digits (default: 2).
#' @return Numeric value of difference in rate (percent).
#' @export
DiffRate <- function(df1, df2, name, sigdig = 2) {
  # Input validation
  if (!is.data.frame(df1)) stop("'df1' must be a data.frame")
  if (!is.data.frame(df2)) stop("'df2' must be a data.frame")
  if (!is.character(name)) stop("'name' must be character")
  if (!name %in% names(df1)) stop("'name' not found in 'df1'")
  if (!name %in% names(df2)) stop("'name' not found in 'df2'")
  
  if (is.data.table(df1)) {
    df1 <- as.data.frame(df1)
    df2 <- as.data.frame(df2)
  }
  temp <- na.omit(df1[, c(name)])
  m1 <- mean(temp)
  temp <- na.omit(df2[, c(name)])
  m2 <- mean(temp)
  out <- round(100 * (m2 - m1), 1)
  return(out)
}

#' Power Calculation for Subgroup Size
#'
#' Calculates power (rejection rate) for subgroups above a minimum size.
#'
#' @param df Data frame.
#' @param minsize Integer. Minimum subgroup size (default: 0).
#' @param sigdig Integer. Number of significant digits (default: 3).
#' @return Numeric value of power (rejection rate).
#' @export
pow_size <- function(df, minsize = 0, sigdig = 3) {
  # Input validation
  if (!is.data.frame(df)) stop("'df' must be a data.frame")
  if (!all(c("rej12", "size.Hc") %in% names(df))) {
    stop("'df' must contain 'rej12' and 'size.Hc' columns")
  }
  
  rej <- with(df, mean(rej12 & size.Hc >= minsize))
  return(round(rej, sigdig))
}
