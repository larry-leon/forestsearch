#' Format Hazard Ratio and Confidence Interval
#'
#' Formats a hazard ratio and confidence interval for display.
#'
#' @param hrest Numeric vector with HR, lower, and upper confidence limits.
#' @return Character string formatted as \"HR (lower, upper)\".
#' @export
hrCI_format <- function(hrest) {
  # Input validation
  if (!is.numeric(hrest)) stop("'hrest' must be numeric")
  if (length(hrest) < 3) stop("'hrest' must have at least 3 elements (HR, lower, upper)")

  # hrest: vector with HR, lower, upper
  sprintf("%.2f (%.2f, %.2f)", hrest[1], hrest[2], hrest[3])
}

#' Calculate n and percent
#'
#' Returns count and percent for a vector relative to a denominator.
#'
#' @param x Vector of values.
#' @param denom Denominator for percent calculation.
#' @return Character string formatted as \"n (percent%)\".
#' @export
n_pcnt <- function(x, denom) {
  # Input validation
  if (is.null(x)) stop("'x' cannot be NULL")
  if (!is.numeric(denom) || length(denom) != 1 || denom <= 0) {
    stop("'denom' must be a positive numeric value")
  }

  n <- length(x)

  # Handle division by zero
  if (denom == 0) {
    return("0 (NA%)")
  }

  sprintf("%d (%.1f%%)", n, 100 * n / denom)
}

#' Prepare subgroup data for analysis
#'
#' Splits a data frame into two subgroups based on a flag and treatment scale.
#'
#' @param df Data frame.
#' @param SG_flag Character. Name of subgroup flag variable.
#' @param est.scale Character. Effect scale ("hr" or "1/hr").
#' @param treat.name Character. Name of treatment variable.
#' @return List with subgroup data frames and treatment variable name.
#' @export
prepare_subgroup_data <- function(df, SG_flag, est.scale, treat.name) {
  # Input validation
  if (!is.data.frame(df)) stop("'df' must be a data.frame")
  if (!is.character(SG_flag) || !SG_flag %in% names(df)) {
    stop("'SG_flag' must be a column name in 'df'")
  }
  if (!est.scale %in% c("hr", "1/hr")) {
    stop("'est.scale' must be either 'hr' or '1/hr'")
  }
  if (!is.character(treat.name) || !treat.name %in% names(df)) {
    stop("'treat.name' must be a column name in 'df'")
  }

  if (est.scale == "1/hr") {
    df$treat2 <- 1 - df[, treat.name]
    treat.name <- "treat2"
  }
  df_0 <- subset(df, df[, SG_flag] == 0)
  df_1 <- subset(df, df[, SG_flag] == 1)
  list(df_0 = df_0, df_1 = df_1, treat.name = treat.name)
}

#' Cox model summary for subgroup
#'
#' Calculates hazard ratio and confidence interval for a subgroup using Cox regression.
#'
#' @param Y Numeric vector of outcome.
#' @param E Numeric vector of event indicators.
#' @param Treat Numeric vector of treatment indicators.
#' @param Strata Vector of strata (optional).
#' @return Character string with formatted HR and CI.
#' @importFrom survival coxph Surv
#' @export
cox_summary <- function(Y, E, Treat, Strata) {
  # Input validation
  if (!is.numeric(Y)) stop("'Y' must be numeric")
  if (!is.numeric(E)) stop("'E' must be numeric")
  if (!is.numeric(Treat)) stop("'Treat' must be numeric")
  if (length(Y) != length(E) || length(Y) != length(Treat)) {
    stop("'Y', 'E', and 'Treat' must have the same length")
  }
  if (!is.null(Strata) && length(Y) != length(Strata)) {
    stop("'Strata' must have the same length as 'Y'")
  }

  fit <- survival::coxph(survival::Surv(Y, E) ~ Treat + strata(Strata), robust = TRUE)
  hr <- summary(fit)$conf.int[c(1, 3, 4)]
  hrCI_format(hr)
}

#' KM median summary for subgroup
#'
#' Calculates median survival for each treatment group using Kaplan-Meier.
#'
#' @param Y Numeric vector of outcome.
#' @param E Numeric vector of event indicators.
#' @param Treat Numeric vector of treatment indicators.
#' @return Numeric vector of medians.
#' @importFrom survival survfit
#' @export
km_summary <- function(Y, E, Treat) {
  # Input validation
  if (!is.numeric(Y)) stop("'Y' must be numeric")
  if (!is.numeric(E)) stop("'E' must be numeric")
  if (!is.numeric(Treat)) stop("'Treat' must be numeric")
  if (length(Y) != length(E) || length(Y) != length(Treat)) {
    stop("'Y', 'E', and 'Treat' must have the same length")
  }

  fit <- summary(survival::survfit(survival::Surv(Y, E) ~ Treat))
  as.numeric(round(fit$table[, "median"], 1))
}

#' Calculate counts for subgroup summary
#'
#' Calculates sample size, treated count, and event count for a subgroup.
#'
#' @param Y Numeric vector of outcome.
#' @param E Numeric vector of event indicators.
#' @param Treat Numeric vector of treatment indicators.
#' @param N Integer. Total sample size.
#' @return List with formatted counts.
#' @export
calculate_counts <- function(Y, E, Treat, N) {
  # Input validation
  if (!is.numeric(Y)) stop("'Y' must be numeric")
  if (!is.numeric(E)) stop("'E' must be numeric")
  if (!is.numeric(Treat)) stop("'Treat' must be numeric")
  if (!is.numeric(N) || length(N) != 1 || N <= 0) {
    stop("'N' must be a positive numeric value")
  }

  n <- n_pcnt(Y, N)
  n_treat <- n_pcnt(Y[Treat == 1], length(Y))
  d <- n_pcnt(Y[E == 1], length(Y))
  list(n = n, n_treat = n_treat, d = d)
}

#' Calculate potential outcome hazard ratio
#'
#' Calculates the average hazard ratio from a potential outcome variable.
#'
#' @param df Data frame.
#' @param potentialOutcome.name Character. Name of potential outcome variable.
#' @return Numeric value of average hazard ratio.
#' @export
calculate_potential_hr <- function(df, potentialOutcome.name) {
  # Input validation
  if (!is.data.frame(df)) stop("'df' must be a data.frame")
  if (!is.character(potentialOutcome.name)) {
    stop("'potentialOutcome.name' must be character")
  }
  if (!potentialOutcome.name %in% names(df)) {
    stop("'potentialOutcome.name' must be a column name in 'df'")
  }

  loghr.po <- df[, potentialOutcome.name]
  round(exp(mean(loghr.po)), 2)
}

#' Format results for subgroup summary
#'
#' Formats results for subgroup summary table.
#'
#' @param subgroup.name Character. Subgroup name.
#' @param n Character. Sample size.
#' @param n.treat Character. Treated count.
#' @param d Character. Event count.
#' @param m1 Numeric. Median or RMST for treatment.
#' @param m0 Numeric. Median or RMST for control.
#' @param drmst Numeric. RMST difference.
#' @param hr Character. Hazard ratio (formatted).
#' @param hr.a Character. Adjusted hazard ratio (optional).
#' @param hr.po Numeric. Potential outcome hazard ratio (optional).
#' @param return.medians Logical. Use medians or RMST.
#' @return Character vector of results.
#' @export
format_results <- function(subgroup.name, n, n.treat, d, m1, m0, drmst, hr,
                          hr.a = NA, hr.po = NA, return.medians = TRUE) {
  # Input validation
  if (!is.character(subgroup.name)) stop("'subgroup.name' must be character")
  if (!is.character(n)) stop("'n' must be character")
  if (!is.character(n.treat)) stop("'n.treat' must be character")
  if (!is.character(d)) stop("'d' must be character")
  if (!is.numeric(m1)) stop("'m1' must be numeric")
  if (!is.numeric(m0)) stop("'m0' must be numeric")
  if (!is.numeric(drmst)) stop("'drmst' must be numeric")
  if (!is.character(hr)) stop("'hr' must be character")

  if (is.na(hr.po)) {
    if (is.na(hr.a)) {
      res <- c(subgroup.name, n, n.treat, d, m1, m0, drmst, hr)
    } else {
      res <- c(subgroup.name, n, n.treat, d, m1, m0, drmst, hr, hr.a)
    }
  } else {
    if (is.na(hr.a)) {
      res <- c(subgroup.name, n, n.treat, d, m1, m0, drmst, hr, hr.po)
    } else {
      res <- c(subgroup.name, n, n.treat, d, m1, m0, drmst, hr, hr.a, hr.po)
    }
  }
  res
}

#' RMST calculation for subgroup
#'
#' Calculates restricted mean survival time (RMST) for a subgroup.
#'
#' @param df Data frame.
#' @param tte.name Character. Name of time-to-event variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @return List with tau, RMST, RMST for treatment, RMST for control.
#' @importFrom weightedSurv df_counting
#' @export
rmst_calculation <- function(df, tte.name = "tte", event.name = "event", treat.name = "treat") {
  # Input validation
  if (!is.data.frame(df)) stop("'df' must be a data.frame")
  if (!is.character(tte.name) || !tte.name %in% names(df)) {
    stop("'tte.name' must be a column name in 'df'")
  }
  if (!is.character(event.name) || !event.name %in% names(df)) {
    stop("'event.name' must be a column name in 'df'")
  }
  if (!is.character(treat.name) || !treat.name %in% names(df)) {
    stop("'treat.name' must be a column name in 'df'")
  }

  if (!requireNamespace("weightedSurv", quietly = TRUE)) {
    stop("Package 'weightedSurv' needed for this function to work. Please install it install_github('larry-leon/weightedSurv').")
  }

  dfcount <- df_counting(df, tte.name = tte.name, event.name = event.name,
                        treat.name = treat.name, arms = c("treat", "control"), by.risk = 1)
  taumax <- with(dfcount, max(at_points[ybar1 > 0 & ybar0 > 0]))
  at_points <- dfcount$at_points
  tau_horizon <- which(at_points <= taumax)
  surv1 <- dfcount$surv1
  surv0 <- dfcount$surv0
  dhat <- c(surv1 - surv0)
  tpoints <- at_points[tau_horizon]
  dhat <- dhat[tau_horizon]
  surv1 <- surv1[tau_horizon]
  surv0 <- surv0[tau_horizon]
  dt <- diff(tpoints)
  mid_dhat <- (head(dhat, -1) + tail(dhat, -1)) / 2
  cumulative_rmst <- c(0, cumsum(mid_dhat * dt))
  rmst <- cumulative_rmst[length(dhat)]
  # m1
  mid_dhat <- (head(surv1, -1) + tail(surv1, -1)) / 2
  cumulative_rmst1 <- c(0, cumsum(mid_dhat * dt))
  rmst1 <- cumulative_rmst1[length(surv1)]
  # m0
  mid_dhat <- (head(surv0, -1) + tail(surv0, -1)) / 2
  cumulative_rmst0 <- c(0, cumsum(mid_dhat * dt))
  rmst0 <- cumulative_rmst0[length(surv0)]

  return(list(tau = taumax, rmst = signif(rmst, 2),
              rmst1 = signif(rmst1, 2), rmst0 = signif(rmst0, 2)))
}

#' Analyze subgroup for summary table
#'
#' Analyzes a subgroup and returns formatted results for summary table.
#'
#' @param df.sub Data frame for subgroup.
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @param strata.name Character. Name of strata variable (optional).
#' @param subgroup.name Character. Subgroup name.
#' @param hr.a Character. Adjusted hazard ratio (optional).
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param return.medians Logical. Use medians or RMST.
#' @param N Integer. Total sample size.
#' @return Character vector of results.
#' @export

analyze_subgroup <- function(df.sub, outcome.name, event.name, treat.name, strata.name,
                                   subgroup.name, hr.a, potentialOutcome.name, return.medians, N) {
  # Input validation
  if (!is.data.frame(df.sub)) stop("'df.sub' must be a data.frame")
  if (!all(c(outcome.name, event.name, treat.name) %in% names(df.sub))) {
    stop("Required columns not found in 'df.sub'")
  }

  Y <- df.sub[, outcome.name]
  E <- df.sub[, event.name]
  Treat <- df.sub[, treat.name]
  Strata <- if (is.null(strata.name)) rep("All", length(Y)) else df.sub[, strata.name]

  # Cox summary
  hr <- cox_summary(Y, E, Treat, Strata)

  # KM medians
  meds <- km_summary(Y, E, Treat)

  # RMST calculation
  rmst <- rmst_calculation(df.sub, outcome.name, event.name, treat.name)

  # Calculate counts - ENSURE THESE RETURN CHARACTER STRINGS
  n <- n_pcnt(Y, N)
  n_treat <- n_pcnt(Y[Treat == 1], length(Y))
  d <- n_pcnt(Y[E == 1], length(Y))

  # Ensure all are character
  n <- as.character(n)
  n_treat <- as.character(n_treat)
  d <- as.character(d)

  # Calculate potential outcome HR if provided
  hr.po <- if (!is.null(potentialOutcome.name)) {
    calculate_potential_hr(df.sub, potentialOutcome.name)
  } else NA

  # Choose median or RMST
  m1 <- if (return.medians) meds[2] else rmst$rmst1
  m0 <- if (return.medians) meds[1] else rmst$rmst0

  # Format results with explicit character conversion
  format_results(
    subgroup.name = as.character(subgroup.name),
    n = as.character(n),
    n.treat = as.character(n_treat),
    d = as.character(d),
    m1 = m1,
    m0 = m0,
    drmst = rmst$rmst,
    hr = as.character(hr),
    hr.a = if (!is.na(hr.a)) as.character(hr.a) else NA,
    hr.po = hr.po,
    return.medians = return.medians
  )
}





# Remove
analyze_subgroup_old <- function(df.sub, outcome.name, event.name, treat.name, strata.name,
                            subgroup.name, hr.a, potentialOutcome.name, return.medians, N) {
  # Input validation
  if (!is.data.frame(df.sub)) stop("'df.sub' must be a data.frame")
  if (!all(c(outcome.name, event.name, treat.name) %in% names(df.sub))) {
    stop("Required columns not found in 'df.sub'")
  }

  Y <- df.sub[, outcome.name]
  E <- df.sub[, event.name]
  Treat <- df.sub[, treat.name]
  Strata <- if (is.null(strata.name)) rep("All", length(Y)) else df.sub[, strata.name]

  hr <- cox_summary(Y, E, Treat, Strata)
  meds <- km_summary(Y, E, Treat)
  rmst <- rmst_calculation(df.sub, outcome.name, event.name, treat.name)
  counts <- calculate_counts(Y, E, Treat, N)
  hr.po <- if (!is.null(potentialOutcome.name)) {
    calculate_potential_hr(df.sub, potentialOutcome.name)
  } else NA

  m1 <- if (return.medians) meds[2] else rmst$rmst1
  m0 <- if (return.medians) meds[1] else rmst$rmst0

  format_results(subgroup.name, counts$n, counts$n.treat, counts$d, m1, m0,
                rmst$rmst, hr, hr.a, hr.po, return.medians)
}

#' Subgroup summary table estimates
#'
#' Returns a summary table of subgroup estimates (HR, RMST, medians, etc.).
#'
#' @param df Data frame.
#' @param SG_flag Character. Subgroup flag variable.
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @param strata.name Character. Name of strata variable (optional).
#' @param hr_1a Character. Adjusted HR for subgroup 1 (optional).
#' @param hr_0a Character. Adjusted HR for subgroup 0 (optional).
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param sg1_name Character. Name for subgroup 1.
#' @param sg0_name Character. Name for subgroup 0.
#' @param draws Integer. Number of draws for resampling (optional).
#' @param details Logical. Print details.
#' @param return.medians Logical. Use medians or RMST.
#' @param est.scale Character. Effect scale ("hr" or "1/hr").
#' @return Data frame of subgroup summary estimates.
#' @export

SG_tab_estimates <- function(df, SG_flag, outcome.name = "tte", event.name = "event",
                                   treat.name = "treat", strata.name = NULL,
                                   hr_1a = NA, hr_0a = NA, potentialOutcome.name = NULL,
                                   sg1_name = NULL, sg0_name = NULL, draws = 0,
                                   details = FALSE, return.medians = TRUE, est.scale = "hr") {

  # Input validation
  if (!is.data.frame(df)) stop("'df' must be a data.frame")
  if (!is.character(SG_flag)) stop("'SG_flag' must be character")
  if (SG_flag != "ITT" && !SG_flag %in% names(df)) {
    stop("'SG_flag' must be 'ITT' or a column name in 'df'")
  }
  if (!all(c(outcome.name, event.name, treat.name) %in% names(df))) {
    stop("Required columns not found in 'df'")
  }

  N <- nrow(df)

  if (SG_flag != "ITT") {
    # Prepare subgroup data
    subgroups <- prepare_subgroup_data(df, SG_flag, est.scale, treat.name)
    df_0 <- subgroups$df_0
    df_1 <- subgroups$df_1
    treat.name <- subgroups$treat.name

    # Analyze each subgroup
    res_0 <- analyze_subgroup(
      df.sub = df_0,
      outcome.name = outcome.name,
      event.name = event.name,
      treat.name = treat.name,
      strata.name = strata.name,
      subgroup.name = sg0_name,
      hr.a = hr_0a,
      potentialOutcome.name = potentialOutcome.name,
      return.medians = return.medians,
      N = N
    )

    res_1 <- analyze_subgroup(
      df.sub = df_1,
      outcome.name = outcome.name,
      event.name = event.name,
      treat.name = treat.name,
      strata.name = strata.name,
      subgroup.name = sg1_name,
      hr.a = hr_1a,
      potentialOutcome.name = potentialOutcome.name,
      return.medians = return.medians,
      N = N
    )

    res <- rbind(res_0, res_1)

    # Set column names based on what's available
    if (is.na(hr_1a)) {
      if (is.null(potentialOutcome.name)) {
        colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)")
      } else {
        colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                           "HR (95% CI)", "AHR(po)")
      }
    } else {
      if (is.null(potentialOutcome.name)) {
        colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                           "HR (95% CI)", "HR*")
      } else {
        colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                           "HR (95% CI)", "HR*", "AHR(po)")
      }
    }
    return(as.data.frame(res))

  } else {
    # ITT analysis
    res <- analyze_subgroup(
      df.sub = df,
      outcome.name = outcome.name,
      event.name = event.name,
      treat.name = treat.name,
      strata.name = strata.name,
      subgroup.name = "ITT",
      hr.a = NA,
      potentialOutcome.name = potentialOutcome.name,
      return.medians = return.medians,
      N = N
    )

    if (is.na(hr_1a)) {
      if (is.null(potentialOutcome.name)) {
        names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)")
      } else {
        names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                        "HR (95% CI)", "AHR(po)")
      }
    } else {
      if (is.null(potentialOutcome.name)) {
        names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                        "HR (95% CI)", "HR*")
      } else {
        names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                        "HR (95% CI)", "HR*", "AHR(po)")
      }
    }
    return(res)
  }
}



SG_tab_estimates_old <- function(df, SG_flag, outcome.name = "tte", event.name = "event",
                             treat.name = "treat", strata.name = NULL,
                             hr_1a = NA, hr_0a = NA, potentialOutcome.name = NULL,
                             sg1_name = NULL, sg0_name = NULL, draws = 0,
                             details = FALSE, return.medians = TRUE, est.scale = "hr") {
  # Input validation
  if (!is.data.frame(df)) stop("'df' must be a data.frame")
  if (!is.character(SG_flag)) stop("'SG_flag' must be character")
  if (SG_flag != "ITT" && !SG_flag %in% names(df)) {
    stop("'SG_flag' must be 'ITT' or a column name in 'df'")
  }
  if (!all(c(outcome.name, event.name, treat.name) %in% names(df))) {
    stop("Required columns not found in 'df'")
  }

  N <- nrow(df)
  if (SG_flag != "ITT") {
    subgroups <- prepare_subgroup_data(df, SG_flag, est.scale, treat.name)
    df_0 <- subgroups$df_0
    df_1 <- subgroups$df_1
    treat.name <- subgroups$treat.name

    res_0 <- analyze_subgroup(df_0, outcome.name, event.name, treat.name, strata.name,
                              sg0_name, hr_0a, potentialOutcome.name, return.medians, N)
    res_1 <- analyze_subgroup(df_1, outcome.name, event.name, treat.name, strata.name,
                              sg1_name, hr_1a, potentialOutcome.name, return.medians, N)
    res <- rbind(res_0, res_1)

    # Set column names
    if (is.na(hr_1a)) {
      if (is.null(potentialOutcome.name)) {
        colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)")
      } else {
        colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                          "HR (95% CI)", "AHR(po)")
      }
    } else {
      if (is.null(potentialOutcome.name)) {
        colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                          "HR (95% CI)", "HR*")
      } else {
        colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                          "HR (95% CI)", "HR*", "AHR(po)")
      }
    }
    return(res)
  } else {
    res <- analyze_subgroup(df, outcome.name, event.name, treat.name, strata.name,
                           "ITT", NA, potentialOutcome.name, return.medians, N)
    if (is.na(hr_1a)) {
      if (is.null(potentialOutcome.name)) {
        names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)")
      } else {
        names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                       "HR (95% CI)", "AHR(po)")
      }
    } else {
      if (is.null(potentialOutcome.name)) {
        names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                       "HR (95% CI)", "HR*")
      } else {
        names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                       "HR (95% CI)", "HR*", "AHR(po)")
      }
    }
    return(res)
  }
}

#' Subgroup summary tables (gt output)
#'
#' Returns formatted summary tables for subgroups using the gt package, with customizable decimal precision.
#'
#' @param fs ForestSearch results object.
#' @param which_df Character. Which data frame to use ("est" or "testing").
#' @param est_caption Character. Caption for estimates table.
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param hr_1a Character. Adjusted HR for subgroup 1 (optional).
#' @param hr_0a Character. Adjusted HR for subgroup 0 (optional).
#' @param ndecimals Integer. Number of decimals for formatted numbers (default: 3).
#'
#' @return List with gt tables for estimates and subgroups.
#'
#' @importFrom gt gt fmt_number tab_header
#' @export

sg_tables <- function(fs, which_df = "est", est_caption = "Training data estimates",
                            potentialOutcome.name = NULL, hr_1a = NA, hr_0a = NA, ndecimals = 3) {

  # Input validation
  if (!is.list(fs)) stop("'fs' must be a list (ForestSearch results object)")
  if (!which_df %in% c("est", "testing")) {
    stop("'which_df' must be either 'est' or 'testing'")
  }

  # Check for gt package
  use_gt <- requireNamespace("gt", quietly = TRUE)

  if (!use_gt) {
    message("Note: gt package not available. Returning standard data frames instead of gt tables.")
  }

  # Select data frame
  if (which_df == "est") {
    df <- fs$df.est
  } else if (which_df == "testing") {
    df <- fs$df.test
  } else {
    stop("which_df must be 'est' or 'testing'")
  }

  # Get arguments from ForestSearch call
  args_fs <- fs$args_call_all

  # ITT estimates
  itt_est <- SG_tab_estimates(
    df = df,
    SG_flag = "ITT",
    outcome.name = args_fs$outcome.name,
    event.name = args_fs$event.name,
    treat.name = args_fs$treat.name,
    strata.name = NULL,
    hr_1a = NA,
    hr_0a = NA,
    potentialOutcome.name = potentialOutcome.name,
    sg1_name = NULL,
    sg0_name = NULL,
    est.scale = args_fs$est.scale,
    return.medians = TRUE
  )

  # Subgroup estimates
  sg_est <- SG_tab_estimates(
    df = df,
    SG_flag = "treat.recommend",
    outcome.name = args_fs$outcome.name,
    event.name = args_fs$event.name,
    treat.name = args_fs$treat.name,
    strata.name = NULL,
    hr_1a = hr_1a,
    hr_0a = hr_0a,
    potentialOutcome.name = potentialOutcome.name,
    sg1_name = "Recommend",
    sg0_name = "Questionable",
    est.scale = args_fs$est.scale,
    return.medians = TRUE
  )

  # Combine tables
  tab_est <- as.data.frame(rbind(itt_est, sg_est))

  # Create gt table if package is available
  if (use_gt) {
    tab_estimates <- gt::gt(tab_est, caption = est_caption, auto_align = TRUE)
  } else {
    tab_estimates <- tab_est
  }

  # Get top subgroups based on sg_focus
  if (!is.null(fs$grp.consistency)) {
    if (fs$sg_focus == "hr") {
      sg10 <- as.data.frame(fs$grp.consistency$out_hr$result)
    } else if (fs$sg_focus %in% c("minSG", "hrMinSG")) {
      sg10 <- as.data.frame(fs$grp.consistency$out_minSG$result)
    } else if (fs$sg_focus %in% c("maxSG", "hrMaxSG")) {
      sg10 <- as.data.frame(fs$grp.consistency$out_maxSG$result)
    } else {
      sg10 <- as.data.frame(fs$grp.consistency$out_hr$result)  # Default to hr
    }

    # Format based on maxk
    maxk <- args_fs$maxk
    if (is.null(maxk)) maxk <- 2  # Default

    if (maxk == 1 && "M.1" %in% names(sg10)) {
      sg10_display <- sg10[, c("M.1", "N", "E", "hr", "Pcons"), drop = FALSE]

      if (args_fs$est.scale == "1/hr") {
        sg10_display$hr <- 1 / sg10_display$hr
      }

      if (use_gt) {
        sg10_out <- sg10_display |>
          gt::gt() |>
          gt::fmt_number(columns = c(4, 5), decimals = ndecimals) |>
          gt::fmt_number(columns = c(2, 3), decimals = 0) |>
          gt::tab_header(title = "Subgroups formed by single-factors",
                         subtitle = "maxk=1")
      } else {
        sg10_out <- sg10_display
      }

    } else if (maxk == 2 && all(c("M.1", "M.2") %in% names(sg10))) {
      sg10_display <- sg10[, c("M.1", "M.2", "N", "E", "hr", "Pcons"), drop = FALSE]

      if (args_fs$est.scale == "1/hr") {
        sg10_display$hr <- 1 / sg10_display$hr
      }

      if (use_gt) {
        sg10_out <- sg10_display |>
          gt::gt() |>
          gt::fmt_number(columns = c(5, 6), decimals = ndecimals) |>
          gt::fmt_number(columns = c(3, 4), decimals = 0) |>
          gt::tab_header(title = "Subgroups formed by two-factors",
                         subtitle = "maxk=2")
      } else {
        sg10_out <- sg10_display
      }

    } else {
      # Fallback for other cases
      sg10_out <- if (use_gt) gt::gt(sg10) else sg10
    }
  } else {
    # No consistency results
    sg10_out <- if (use_gt) {
      gt::gt(data.frame(Message = "No subgroup consistency results available"))
    } else {
      data.frame(Message = "No subgroup consistency results available")
    }
  }

  return(list(
    tab_estimates = tab_estimates,
    sg10_out = sg10_out
  ))
}



sg_tables_old <- function(fs, which_df = "est", est_caption = "Training data estimates",
                     potentialOutcome.name = NULL, hr_1a = NA, hr_0a = NA, ndecimals = 3) {
  # Input validation
  if (!is.list(fs)) stop("'fs' must be a list (ForestSearch results object)")
  if (!which_df %in% c("est", "testing")) {
    stop("'which_df' must be either 'est' or 'testing'")
  }

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' required.")
  }

  if (which_df == "est") df <- fs$df.est
  if (which_df == "testing") df <- fs$df.test

  args_fs <- fs$args_call_all
  args_tab <- names(formals(SG_tab_estimates))
  # align with args_call_all
  args_tab_filtered <- args_fs[args_fs %in% names(args_tab)]
  args_tab_filtered$df <- df
  args_tab_filtered$sg0_name <- "Questionable"
  args_tab_filtered$sg1_name <- "Recommend"
  args_tab_filtered$hr_1a <- hr_1a
  args_tab_filtered$hr_0a <- hr_0a

  # ITT estimates
  args_tab_filtered$SG_flag <- "ITT"

  aa <- do.call(SG_tab_estimates, args_tab_filtered)

  args_tab_filtered$SG_flag <- "treat.recommend"

  bb <- do.call(SG_tab_estimates, args_tab_filtered)

  tab_est <- as.data.frame(rbind(aa, bb))

  tab_estimates <- gt::gt(tab_est, caption = est_caption, auto_align = TRUE)

  if (fs$sg_focus == "hr") sg10 <- as.data.frame(fs$grp.consistency$out_hr$result)
  if (fs$sg_focus %in% c("minSG", "hrMinSG")) {
    sg10 <- as.data.frame(fs$grp.consistency$out_minSG$result)
  }
  if (fs$sg_focus %in% c("maxSG", "hrMaxSG")) {
    sg10 <- as.data.frame(fs$grp.consistency$out_maxSG$result)
  }

  if (args_fs$maxk == 1) {
    # Subgroup identification ="M.1"
    # Consistency rate = "Pcons"
    # hr = "hr"
    # sample size = "N"
    # Total events = "E"
    sg10 <- sg10[, c("M.1", "N", "E", "hr", "Pcons")]

    if (args_fs$est.scale == "1/hr") {
      sg10$hr <- c(1 / sg10$hr)
    }

    sg10_out <- sg10 |> gt::gt() |>
      gt::fmt_number(columns = c(4, 5), decimals = ndecimals) |>
      gt::fmt_number(columns = c(2, 3), decimals = 0) |>
      gt::tab_header(title = "Subgroups formed by single-factors",
                    subtitle = "maxk=1")
  }
  if (args_fs$maxk == 2) {
    sg10 <- sg10[, c("M.1", "M.2", "N", "E", "hr", "Pcons")]

    if (args_fs$est.scale == "1/hr") {
      sg10$hr <- c(1 / sg10$hr)
    }

    sg10_out <- sg10 |> gt::gt() |>
      gt::fmt_number(columns = c(5, 6), decimals = ndecimals) |>
      gt::fmt_number(columns = c(3, 4), decimals = 0) |>
      gt::tab_header(title = "Subgroups formed by two-factors",
                    subtitle = "maxk=2")
  }

  return(list(tab_estimates = tab_estimates, sg10_out = sg10_out))
}
