#' Format Hazard Ratio and Confidence Interval
#'
#' Formats a hazard ratio and confidence interval for display.
#'
#' @param hrest Numeric vector with HR, lower, and upper confidence limits.
#' @return Character string formatted as \"HR (lower, upper)\".
#' @examples
#' hrCI_format(c(1.45, 0.98, 2.14))
#' @export

hrCI_format <- function(hrest) {
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
#' @examples
#' n_pcnt(1:30, 100)
#' @export

n_pcnt <- function(x, denom) {
  n <- length(x)
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
#' @examples
#' df <- data.frame(
#'   treat = c(0, 1, 0, 1, 0),
#'   sg_flag = c(1, 1, 0, 0, 1)
#' )
#' result <- prepare_subgroup_data(df, SG_flag = "sg_flag",
#'                                  est.scale = "hr", treat.name = "treat")
#' nrow(result$df_1)
#' @export

prepare_subgroup_data <- function(df, SG_flag, est.scale, treat.name) {
  if (est.scale == "1/hr") {
    df$treat2 <- 1 - df[, treat.name]
    treat.name <- "treat2"
  }
  df_0 <- subset(df, df[, SG_flag] == 0)
  df_1 <- subset(df, df[, SG_flag] == 1)
  list(df_0 = df_0, df_1 = df_1, treat.name = treat.name)
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
#' @keywords internal

km_summary <- function(Y, E, Treat) {
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
#' @keywords internal

calculate_counts <- function(Y, E, Treat, N) {
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
#' @keywords internal

calculate_potential_hr <- function(df, potentialOutcome.name) {
  loghr.po <- df[, potentialOutcome.name]
  round(exp(mean(loghr.po)), 2)
}

#' Format results for subgroup summary
#'
#' Formats results for subgroup summary table.
#'
#' @param subgroup_name Character. Subgroup name.
#' @param n Character. Sample size.
#' @param n_treat Character. Treated count.
#' @param d Character. Event count.
#' @param m1 Numeric. Median or RMST for treatment.
#' @param m0 Numeric. Median or RMST for control.
#' @param drmst Numeric. RMST difference.
#' @param hr Character. Hazard ratio (formatted).
#' @param hr_a Character. Adjusted hazard ratio (optional).
#' @param hr_po Numeric. Potential outcome hazard ratio (optional).
#' @param return_medians Logical. Use medians or RMST.
#' @return Character vector of results.
#' @keywords internal

format_results <- function(subgroup_name, n, n_treat, d, m1, m0, drmst, hr, hr_a = NA, hr_po = NA, return_medians = TRUE) {
  if (is.na(hr_po)) {
    if (is.na(hr_a)) {
      res <- c(subgroup_name, n, n_treat, d, m1, m0, drmst, hr)
    } else {
      res <- c(subgroup_name, n, n_treat, d, m1, m0, drmst, hr, hr_a)
    }
  } else {
    if (is.na(hr_a)) {
      res <- c(subgroup_name, n, n_treat, d, m1, m0, drmst, hr, hr_po)
    } else {
      res <- c(subgroup_name, n, n_treat, d, m1, m0, drmst, hr, hr_a, hr_po)
    }
  }
  res
}


# =========================================================================
# GLM SUBGROUP SUMMARY FUNCTIONS
# =========================================================================

#' Analyze subgroup for GLM outcomes
#'
#' GLM counterpart to \code{\link{analyze_subgroup}}.  Computes per-arm outcome
#' rates, sample sizes, and the treatment effect estimate (RD, OR, etc.)
#' using the estimator closure.
#'
#' @param df_sub Data frame for the subgroup.
#' @param outcome.name Character. Name of outcome variable.
#' @param treat.name Character. Name of treatment variable.
#' @param subgroup_name Character. Label for this subgroup.
#' @param effect_a Character. Adjusted effect estimate string (optional).
#' @param estimator_fn Closure from `make_effect_estimator()`, or NULL.
#' @param effect_measure Character. Effect measure label (e.g., "RD", "OR").
#' @param outcome_type Character. One of \code{"binary"}, \code{"continuous"},
#'   or \code{"count"}.  Controls formatting: binary displays per-arm means
#'   as percentages (\code{"Rate"}); continuous and count display raw means
#'   (\code{"Mean"}).  Default \code{"binary"}.
#' @param N Integer. Total sample size (for percentage calculation).
#' @param digits Integer. Number of decimal places for arm-mean / Diff /
#'   effect-estimate / CI numeric formatting.  Default 4 to provide
#'   headroom for downstream display reformatting (e.g., via
#'   \code{summarize_bootstrap_results(digits = ...)}, which can round
#'   down but cannot recover precision lost at construction).  Note this
#'   does not affect the count-with-percent-of-total formatting in the
#'   \code{n} and \code{n1} columns, which remain at 1 decimal.
#'
#' @return One-row data frame of formatted subgroup results.  Column names
#'   are placeholders (\code{Subgroup}, \code{n}, \code{n1}, \code{rate0},
#'   \code{rate1}, \code{Diff}, \code{effect}, optionally \code{effect_a})
#'   and are overwritten by the calling function
#'   (\code{\link{SG_tab_estimates_glm}}) with the user-facing labels
#'   (e.g., \code{"Rate(C)"} vs \code{"Mean(C)"}, the effect-measure
#'   header).  Returning a data frame -- rather than a named character
#'   vector via \code{c(...)} -- ensures type stability across the ITT
#'   and subgroup call paths and prevents downstream rendering failures
#'   in consumers that go through \code{as.data.frame()} (e.g., gt).
#' @keywords internal
analyze_subgroup_glm <- function(df_sub, outcome.name, treat.name,
                                 subgroup_name, effect_a = NA,
                                 estimator_fn = NULL,
                                 effect_measure = "RD",
                                 outcome_type = "binary", N,
                                 digits = 4) {

  Y     <- df_sub[[outcome.name]]
  Treat <- df_sub[[treat.name]]
  n_tot <- length(Y)

  # Per-arm counts
  n0 <- sum(Treat == 0)
  n1 <- sum(Treat == 1)

  # Per-arm outcome summaries
  rate0 <- mean(Y[Treat == 0])
  rate1 <- mean(Y[Treat == 1])

  # Sprintf format strings parameterized by `digits`.  These produce a
  # uniformly-precise base representation; downstream consumers (e.g.,
  # summarize_bootstrap_results) can reformat to lower precision via
  # parse-and-reformat.  Construction-time precision is the ceiling for
  # what display-time can show — values rounded here cannot be recovered.
  fmt_val <- paste0("%.", digits, "f")
  fmt_signed <- paste0("%+.", digits, "f")
  fmt_pct <- paste0("%.", digits, "f%%")
  fmt_signed_pct <- paste0("%+.", digits, "f%%")
  fmt_ci <- paste0(fmt_val, " (", fmt_val, ", ", fmt_val, ")")

  # Effect estimate via closure (if provided)
  if (!is.null(estimator_fn)) {
    res_est <- tryCatch(estimator_fn(df_sub), error = function(e) NULL)
    if (!is.null(res_est) && !is.na(res_est$estimate)) {
      est   <- res_est$estimate
      se    <- res_est$se
      z     <- 1.96
      lower <- est - z * se
      upper <- est + z * se

      # Format depends on scale
      if (effect_measure %in% c("OR", "RR", "IRR")) {
        # Log scale: exponentiate for display
        effect_str <- sprintf(fmt_ci, exp(est), exp(lower), exp(upper))
      } else {
        # Identity scale (RD, IRD, MD)
        effect_str <- sprintf(fmt_ci, est, lower, upper)
      }
    } else {
      effect_str <- "NA"
    }
  } else {
    effect_str <- "NA"
  }

  # Formatted counts
  n_fmt   <- sprintf("%d (%.1f%%)", n_tot, 100 * n_tot / N)
  n1_fmt  <- sprintf("%d (%.1f%%)", n1, 100 * n1 / n_tot)

  # Per-arm summaries: percentage for binary, raw mean for continuous/count
  is_proportion <- outcome_type %in% c("binary")
  if (is_proportion) {
    r0_fmt <- sprintf(fmt_pct, 100 * rate0)
    r1_fmt <- sprintf(fmt_pct, 100 * rate1)
    rd_fmt <- sprintf(fmt_signed_pct, 100 * (rate1 - rate0))
  } else {
    r0_fmt <- sprintf(fmt_val, rate0)
    r1_fmt <- sprintf(fmt_val, rate1)
    rd_fmt <- sprintf(fmt_signed, rate1 - rate0)
  }

  # Return a 1-row data frame for type stability across the ITT and
  # subgroup call paths.  Using c(...) here would produce a named
  # character vector whose colnames() is NULL, breaking downstream
  # consumers that rely on column-name introspection (e.g.,
  # forestsearch_KfoldOut) and producing 1-column "scattered text" when
  # passed through as.data.frame() (e.g., cv_summary_tables -> gt::gt).
  # Placeholder column names are overwritten by SG_tab_estimates_glm().
  if (is.na(effect_a)) {
    data.frame(
      Subgroup = subgroup_name,
      n        = n_fmt,
      n1       = n1_fmt,
      rate0    = r0_fmt,
      rate1    = r1_fmt,
      Diff     = rd_fmt,
      effect   = effect_str,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Subgroup = subgroup_name,
      n        = n_fmt,
      n1       = n1_fmt,
      rate0    = r0_fmt,
      rate1    = r1_fmt,
      Diff     = rd_fmt,
      effect   = effect_str,
      effect_a = as.character(effect_a),
      stringsAsFactors = FALSE
    )
  }
}


#' Subgroup Summary Table for GLM Outcomes
#'
#' GLM counterpart to \code{\link{SG_tab_estimates}}.  Produces a summary data frame
#' with per-arm outcome rates and the GLM effect estimate for each subgroup.
#'
#' @param df Data frame with a `treat.recommend` column (or ITT).
#' @param SG_flag Character. "ITT" for the full population, or the name of
#'   the subgroup flag column (e.g., "treat.recommend").
#' @param outcome.name Character. Name of outcome variable.
#' @param treat.name Character. Name of treatment variable.
#' @param estimator_fn Closure from `make_effect_estimator()`, or NULL.
#' @param effect_measure Character. Effect measure label (e.g., "RD", "OR").
#' @param outcome_type Character. One of \code{"binary"}, \code{"continuous"},
#'   or \code{"count"}.  Controls column labels and formatting: binary uses
#'   \code{"Rate(C)"}/\code{"Rate(T)"} with percentages; continuous and count
#'   use \code{"Mean(C)"}/\code{"Mean(T)"} with raw values.
#'   Default \code{"binary"}.
#' @param effect_a_1 Character. Adjusted effect for subgroup 1 (optional).
#' @param effect_a_0 Character. Adjusted effect for subgroup 0 (optional).
#' @param sg1_name Character. Label for subgroup 1 (treat.recommend == 1).
#' @param sg0_name Character. Label for subgroup 0 (treat.recommend == 0).
#' @param est.scale Character. Effect scale ("hr" or "1/hr").
#' @param digits Integer. Decimal places for numeric formatting; passed
#'   through to \code{\link{analyze_subgroup_glm}}.  Default 4 to provide
#'   round-down headroom for downstream display reformatting (e.g.,
#'   \code{summarize_bootstrap_results(digits = ...)}).
#'
#' @return Data frame.  One row when \code{SG_flag = "ITT"}; two rows
#'   (sg0 then sg1) when \code{SG_flag = "treat.recommend"}.  Columns:
#'   \code{Subgroup}, \code{n}, \code{n1}, \code{Rate(C)} or
#'   \code{Mean(C)} (binary vs continuous/count), \code{Rate(T)} or
#'   \code{Mean(T)}, \code{Diff}, and the effect-measure header
#'   \code{"\{effect_measure\} (95% CI)"}.  When \code{effect_a_*} is
#'   supplied, an additional \code{"\{effect_measure\}*"} column carries
#'   the bias-corrected estimate.
#' @keywords internal
SG_tab_estimates_glm <- function(df, SG_flag,
                                 outcome.name, treat.name,
                                 estimator_fn = NULL,
                                 effect_measure = "RD",
                                 outcome_type = "binary",
                                 effect_a_1 = NA, effect_a_0 = NA,
                                 sg1_name = "Recommend",
                                 sg0_name = "Questionable",
                                 est.scale = "hr",
                                 digits = 4) {
  N <- nrow(df)

  # Effect label for column header
  eff_label <- paste0(effect_measure, " (95% CI)")

  # Column labels: Rate/% for binary, Mean for continuous/count
  is_proportion <- outcome_type %in% c("binary")
  c_label <- if (is_proportion) "Rate(C)" else "Mean(C)"
  t_label <- if (is_proportion) "Rate(T)" else "Mean(T)"

  if (SG_flag == "ITT") {
    res <- analyze_subgroup_glm(
      df, outcome.name, treat.name,
      subgroup_name  = "ITT",
      effect_a       = NA,
      estimator_fn   = estimator_fn,
      effect_measure = effect_measure,
      outcome_type   = outcome_type,
      N = N,
      digits = digits
    )
    col_names <- c("Subgroup", "n", "n1", c_label, t_label, "Diff",
                    eff_label)
    names(res) <- col_names
    return(res)
  }

  # Subgroup split
  subgroups <- prepare_subgroup_data(df, SG_flag, est.scale, treat.name)
  treat_name_use <- subgroups$treat.name

  res_0 <- analyze_subgroup_glm(
    subgroups$df_0, outcome.name, treat_name_use,
    subgroup_name  = sg0_name,
    effect_a       = effect_a_0,
    estimator_fn   = estimator_fn,
    effect_measure = effect_measure,
    outcome_type   = outcome_type,
    N = N,
    digits = digits
  )

  res_1 <- analyze_subgroup_glm(
    subgroups$df_1, outcome.name, treat_name_use,
    subgroup_name  = sg1_name,
    effect_a       = effect_a_1,
    estimator_fn   = estimator_fn,
    effect_measure = effect_measure,
    outcome_type   = outcome_type,
    N = N,
    digits = digits
  )

  res <- rbind(res_0, res_1)

  if (is.na(effect_a_1)) {
    colnames(res) <- c("Subgroup", "n", "n1", c_label, t_label, "Diff",
                        eff_label)
  } else {
    colnames(res) <- c("Subgroup", "n", "n1", c_label, t_label, "Diff",
                        eff_label, paste0(effect_measure, "*"))
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
#' @examples
#' \dontrun{
#' library(survival)
#' df <- data.frame(
#'   tte   = gbsg$rfstime / 30.4375,
#'   event = gbsg$status,
#'   treat = gbsg$hormon
#' )
#' rmst_calculation(df, tte.name = "tte", event.name = "event",
#'                  treat.name = "treat")
#' }
#' @export

rmst_calculation <- function(df,tte.name = "tte",event.name = "event",treat.name = "treat"){
  if (!requireNamespace("weightedsurv", quietly = TRUE)) {
    stop("Package 'weightedsurv' needed for this function to work. Please install it install_github('larry-leon/weightedsurv').")
  }
  dfcount <- weightedsurv::df_counting(df, tte.name = tte.name, event.name = event.name, treat.name = treat.name, arms = c("treat","control"), by.risk = 1)
  taumax <- with(dfcount, max(at_points[ybar1 > 0 & ybar0 >0]))
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
  mid_dhat <- (head(dhat, -1) + tail(dhat, -1))/2
  cumulative_rmst <- c(0, cumsum(mid_dhat * dt))
  rmst <- cumulative_rmst[length(dhat)]
  # m1
  mid_dhat <- (head(surv1, -1) + tail(surv1, -1))/2
  cumulative_rmst1 <- c(0, cumsum(mid_dhat * dt))
  rmst1 <- cumulative_rmst1[length(surv1)]
  # m0
  mid_dhat <- (head(surv0, -1) + tail(surv0, -1))/2
  cumulative_rmst0 <- c(0, cumsum(mid_dhat * dt))
  rmst0 <- cumulative_rmst0[length(surv0)]

  return(list(tau=taumax, rmst = signif(rmst,2), rmst1 = signif(rmst1,2), rmst0 = signif(rmst0,2)))
}

#' Analyze subgroup for summary table (OPTIMIZED)
#'
#' Analyzes a subgroup and returns formatted results for summary table.
#' Uses optimized cox_summary() and reduces redundant calculations.
#'
#' @param df_sub Data frame for subgroup.
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @param strata.name Character. Name of strata variable (optional).
#' @param subgroup_name Character. Subgroup name.
#' @param hr_a Character. Adjusted hazard ratio (optional).
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param return_medians Logical. Use medians or RMST.
#' @param N Integer. Total sample size.
#' @return Character vector of results.
#' @examples
#' \dontrun{
#' library(survival)
#' df <- data.frame(
#'   tte   = gbsg$rfstime / 30.4375,
#'   event = gbsg$status,
#'   treat = gbsg$hormon
#' )
#' analyze_subgroup(
#'   df_sub        = df,
#'   outcome.name  = "tte",
#'   event.name    = "event",
#'   treat.name    = "treat",
#'   strata.name   = NULL,
#'   subgroup_name = "All",
#'   hr_a          = NA,
#'   potentialOutcome.name = NULL,
#'   return_medians = TRUE,
#'   N             = nrow(df)
#' )
#' }
#' @export

analyze_subgroup <- function(df_sub, outcome.name, event.name, treat.name,
                             strata.name, subgroup_name, hr_a,
                             potentialOutcome.name, return_medians, N) {

  # =========================================================================
  # OPTIMIZATION 1: Extract vectors once (avoid repeated $ operations)
  # =========================================================================

  Y <- df_sub[[outcome.name]]
  E <- df_sub[[event.name]]
  Treat <- df_sub[[treat.name]]
  Strata <- if (is.null(strata.name)) {
    rep("All", length(Y))
  } else {
    df_sub[[strata.name]]
  }

  # =========================================================================
  # OPTIMIZATION 2: Use optimized cox_summary()
  # =========================================================================

  hr <- cox_summary(Y, E, Treat, Strata,
                    use_strata = !is.null(strata.name),
                    return_format = "formatted")

  # =========================================================================
  # OPTIMIZATION 3: Parallel independent calculations
  # =========================================================================

  # KM medians (unchanged, already efficient)
  meds <- km_summary(Y, E, Treat)

  # RMST calculation (unchanged)
  rmst <- rmst_calculation(df_sub, outcome.name, event.name, treat.name)

  # Counts (unchanged)
  counts <- calculate_counts(Y, E, Treat, N)

  # Potential outcome HR (only if needed)
  hr_po <- if (!is.null(potentialOutcome.name)) {
    calculate_potential_hr(df_sub, potentialOutcome.name)
  } else {
    NA
  }

  # =========================================================================
  # OPTIMIZATION 4: Efficient assignment based on return_medians
  # =========================================================================

  m1 <- if (return_medians) meds[2] else rmst$rmst1
  m0 <- if (return_medians) meds[1] else rmst$rmst0

  # Return formatted results
  format_results(subgroup_name, counts$n, counts$n_treat, counts$d,
                 m1, m0, rmst$rmst, hr, hr_a, hr_po, return_medians)
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
#' @param return_medians Logical. Use medians or RMST.
#' @param est.scale Character. Effect scale ("hr" or "1/hr").
#' @return Data frame of subgroup summary estimates.
#' @examples
#' \dontrun{
#' library(survival)
#' df <- data.frame(
#'   tte          = gbsg$rfstime / 30.4375,
#'   event        = gbsg$status,
#'   treat        = gbsg$hormon,
#'   treat.recommend = as.integer(gbsg$er > 0)
#' )
#' SG_tab_estimates(df, SG_flag = "ITT",
#'                  outcome.name = "tte",
#'                  event.name   = "event",
#'                  treat.name   = "treat")
#' }
#' @export

SG_tab_estimates <- function(df, SG_flag, outcome.name = "tte", event.name = "event", treat.name = "treat", strata.name = NULL,
                             hr_1a = NA, hr_0a = NA, potentialOutcome.name = NULL,
                             sg1_name = NULL, sg0_name = NULL, draws = 0,
                             details = FALSE, return_medians = TRUE, est.scale = "hr") {
  N <- nrow(df)
  if (SG_flag != "ITT") {
    subgroups <- prepare_subgroup_data(df, SG_flag, est.scale, treat.name)
    df_0 <- subgroups$df_0
    df_1 <- subgroups$df_1
    treat.name <- subgroups$treat.name

    res_0 <- analyze_subgroup(df_0, outcome.name, event.name, treat.name, strata.name, sg0_name, hr_0a, potentialOutcome.name, return_medians, N)
    res_1 <- analyze_subgroup(df_1, outcome.name, event.name, treat.name, strata.name, sg1_name, hr_1a, potentialOutcome.name, return_medians, N)
    res <- rbind(res_0, res_1)

    # Set column names
    if (is.na(hr_1a)) {
      if (is.null(potentialOutcome.name)) colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)")
      else colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)", "AHR(po)")
    } else {
      if (is.null(potentialOutcome.name)) colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)", "HR*")
      else colnames(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)", "HR*", "AHR(po)")
    }
    return(res)
  } else {
    res <- analyze_subgroup(df, outcome.name, event.name, treat.name, strata.name, "ITT", NA, potentialOutcome.name, return_medians, N)
    if (is.na(hr_1a)) {
      if (is.null(potentialOutcome.name)) names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)")
      else names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)", "AHR(po)")
    } else {
      if (is.null(potentialOutcome.name)) names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)", "HR*")
      else names(res) <- c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST", "HR (95% CI)", "HR*", "AHR(po)")
    }
    return(res)
  }
}



#' Filter and merge arguments for function calls
#'
#' Simplifies the common pattern of filtering arguments from a source list
#' to match a target function's formal parameters, then adding/overriding specific arguments.
#'
#' @param source_args List of all arguments (typically from `mget()` or a stored args list).
#' @param target_func Function whose formals define the filter criteria.
#' @param override_args List of arguments to add or override (optional).
#'
#' @return List of filtered arguments ready for `do.call()`.
#'
#' @details
#' This function:
#' 1. Extracts formal parameter names from `target_func`
#' 2. Keeps only arguments from `source_args` that match those names
#' 3. Adds or overrides with any `override_args` provided
#'
#' Reduces boilerplate and improves readability across the codebase.
#'
#' @examples
#' \dontrun{
#' # Instead of:
#' args_FS <- names(formals(get_FSdata))
#' args_FS_filtered <- args_call_all[names(args_call_all) %in% args_FS]
#' args_FS_filtered$df.analysis <- df.analysis
#' args_FS_filtered$grf_cuts <- grf_cuts
#' FSdata <- do.call(get_FSdata, args_FS_filtered)
#'
#' # You now write:
#' FSdata <- do.call(get_FSdata,
#'   filter_call_args(args_call_all, get_FSdata,
#'                    list(df.analysis = df.analysis, grf_cuts = grf_cuts)))
#' }
#'
#' @export

filter_call_args <- function(source_args, target_func, override_args = NULL) {
  target_params <- names(formals(target_func))
  filtered <- source_args[names(source_args) %in% target_params]

  if (!is.null(override_args)) {
    filtered[names(override_args)] <- override_args
  }

  filtered
}



#' Enhanced Subgroup Summary Tables (gt output)
#'
#' Returns formatted summary tables for subgroups using the gt package,
#' with search metadata and customizable decimal precision. Produces two
#' tables: a treatment effect estimates table and an identified subgroups
#' table, each with fully customizable titles and subtitles.
#'
#' @param fs ForestSearch results object, or a \code{grf_glm_result} object
#'   from \code{\link{grf.subg.harm.glm}}.  When a GRF object is supplied,
#'   only Table 1 (treatment effect estimates) is produced; Table 2
#'   (identified subgroups with consistency) returns \code{NULL}.
#' @param which_df Character. Which data frame to use ("est" or "testing").
#' @param est_title Character or NULL. Main title for the estimates table
#'   (default: "Treatment Effect Estimates"). Rendered as bold markdown.
#'   Set to NULL to suppress the title and display only `est_caption`.
#' @param est_caption Character. Subtitle for the estimates table
#'   (default: "Training data estimates").
#' @param sg_title Character or NULL. Main title for the identified subgroups
#'   table (default: "Identified Subgroups"). Rendered as bold markdown.
#'   Set to NULL to suppress the title and display only `sg_subtitle`.
#' @param sg_subtitle Character or NULL. Subtitle for the identified subgroups
#'   table. When NULL (default), an informative subtitle is auto-generated
#'   from `maxk` (e.g., "Two-factor subgroups (maxk=2)").
#' @param potentialOutcome.name Character. Name of potential outcome variable
#'   (optional).
#' @param hr_1a Character. Adjusted HR for subgroup 1 (optional).
#' @param hr_0a Character. Adjusted HR for subgroup 0 (optional).
#' @param ndecimals Integer. Number of decimals for formatted numbers
#'   (default: 3).  Controls precision in both \code{sg10_out}
#'   (formatted via \code{gt::fmt_number}) and \code{tab_estimates}
#'   (threaded as \code{digits} to \code{\link{SG_tab_estimates_glm}}
#'   for arm summaries, Diff, effect-estimate CI, and bias-corrected
#'   CI).  GLM and GRF-GLM paths only; the survival path uses its own
#'   formatting in \code{\link{SG_tab_estimates}}.
#' @param include_search_info Logical. Include search metadata table
#'   (default: TRUE).
#' @param subgroup_notation Character or \code{NULL}.  When set to
#'   \code{"harm"} or \code{"benefit"}, overrides the default
#'   "Recommend"/"Questionable" labels with the Unicode labels from
#'   \code{\link{get_sg_labels}} (e.g., \enc{Ĥ}{H-hat}/\enc{Ĥᶜ}{H-hat-c}
#'   for harm).  Default \code{NULL} (use "Recommend"/"Questionable").
#' @param font_size Numeric. Font size in pixels for table text (default: 12).
#'
#' @return List with gt tables for estimates, subgroups, and optionally
#'   search info.
#'
#' @importFrom gt gt fmt_number tab_header tab_spanner tab_source_note tab_options md px
#' @examples
#' \dontrun{
#' library(survival)
#' fs <- forestsearch(
#'   gbsg,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name     = "rfstime",
#'   treat.name       = "hormon",
#'   event.name       = "status",
#'   fs.splits        = 50,
#'   use_lasso        = FALSE,
#'   use_grf          = FALSE
#' )
#' tabs <- sg_tables(fs)
#' tabs$tab_estimates
#' }
#' @export
sg_tables <- function(fs,
                      which_df = "est",
                      est_title = "Treatment Effect Estimates",
                      est_caption = "Training data estimates",
                      sg_title = "Identified Subgroups",
                      sg_subtitle = NULL,
                      potentialOutcome.name = NULL,
                      hr_1a = NA,
                      hr_0a = NA,
                      ndecimals = 3,
                      include_search_info = TRUE,
                      subgroup_notation = NULL,
                      font_size = 12) {

  # -- Subgroup labels from notation ------------------------------------------
  # Default: "Recommend" / "Questionable".  When subgroup_notation is set,
  # use get_sg_labels() for consistent Unicode labeling across the package.
  if (!is.null(subgroup_notation)) {
    subgroup_notation <- match.arg(subgroup_notation, c("harm", "benefit"))
    L <- get_sg_labels(subgroup_notation)
    if (subgroup_notation == "harm") {
      # treat.recommend == 0 IS the harm subgroup
      sg0_label <- L$sg_hat     # Ĥ
      sg1_label <- L$sg_hat_c   # Ĥᶜ
    } else {
      # treat.recommend == 1 IS the benefit subgroup
      sg1_label <- L$sg_hat     # Ĝ
      sg0_label <- L$sg_hat_c   # Ĝᶜ
    }
  } else {
    sg0_label <- "Questionable"
    sg1_label <- "Recommend"
  }

  # Footnote target row: the "identified" subgroup
  # harm: sg0 (Ĥ);  benefit: sg1 (Ĝ);  default: sg0 (Questionable)
  fn_target_label <- if (!is.null(subgroup_notation) &&
                         subgroup_notation == "benefit") {
    sg1_label
  } else {
    sg0_label
  }

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' required.")
  }

  # =========================================================================
  # GRF STANDALONE DISPATCH
  # =========================================================================
  # When fs is a grf_glm_result object (from grf.subg.harm.glm), produce

  # Table 1 (ITT + subgroup estimates) only.  Table 2 (candidate
  # enumeration with Pcons) is FS-specific and returns NULL.
  if (inherits(fs, "grf_glm_result")) {
    df <- fs$data

    outcome.name   <- fs$outcome.name
    treat.name_grf <- fs$treat.name
    effect_measure <- fs$effect_measure %||% "OR"
    outcome_type   <- fs$outcome_type %||% "binary"
    offset.name    <- fs$offset.name

    # GRF stores effect_measure on the link scale (e.g., "log_OR");
    # SG_tab_estimates_glm and make_effect_estimator expect the display
    # scale (e.g., "OR").  Convert:
    em_map <- c(log_OR = "OR", log_RR = "RR", log_IRR = "IRR")
    if (effect_measure %in% names(em_map)) {
      effect_measure <- em_map[[effect_measure]]
    }

    # Build estimator closure
    estimator_fn_grf <- tryCatch(
      make_effect_estimator(
        outcome_type      = outcome_type,
        treat.name        = treat.name_grf,
        outcome.name      = outcome.name,
        offset.name       = offset.name,
        effect_measure    = effect_measure,
        adjust_covariates = fs$args_call_all$adjust_covariates
      ),
      error = function(e) NULL
    )

    # ITT row
    aa_grf <- SG_tab_estimates_glm(
      df, SG_flag = "ITT",
      outcome.name   = outcome.name,
      treat.name     = treat.name_grf,
      estimator_fn   = estimator_fn_grf,
      effect_measure = effect_measure,
      outcome_type   = outcome_type,
      est.scale      = "hr",
      digits         = ndecimals
    )

    # Subgroup rows
    bb_grf <- SG_tab_estimates_glm(
      df, SG_flag = "treat.recommend",
      outcome.name   = outcome.name,
      treat.name     = treat.name_grf,
      estimator_fn   = estimator_fn_grf,
      effect_measure = effect_measure,
      outcome_type   = outcome_type,
      effect_a_1     = hr_1a,
      effect_a_0     = hr_0a,
      sg1_name       = sg1_label,
      sg0_name       = sg0_label,
      est.scale      = "hr",
      digits         = ndecimals
    )

    tab_est_grf <- as.data.frame(rbind(aa_grf, bb_grf))

    tab_estimates_grf <- gt::gt(tab_est_grf, auto_align = TRUE) |>
      gt::tab_header(
        title = if (!is.null(est_title)) {
          gt::md(paste0("**", est_title, "**"))
        } else {
          est_caption
        },
        subtitle = if (!is.null(est_title)) est_caption
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size),
        heading.title.font.size = gt::px(font_size + 2),
        heading.subtitle.font.size = gt::px(font_size),
        column_labels.font.size = gt::px(font_size)
      )

    # Subgroup definition footnote
    sg_def_grf <- if (!is.null(fs$sg.harm.id) && length(fs$sg.harm.id) > 0) {
      paste(fs$sg.harm.id, collapse = " & ")
    } else {
      NULL
    }

    if (!is.null(sg_def_grf) && nzchar(sg_def_grf)) {
      fn_word <- if (!is.null(subgroup_notation)) {
        paste0("Identified ", subgroup_notation, " subgroup (GRF)")
      } else {
        "Identified subgroup (GRF)"
      }
      tab_estimates_grf <- tab_estimates_grf |>
        gt::tab_footnote(
          footnote = gt::md(
            paste0("**", fn_word, ":** ", sg_def_grf)
          ),
          locations = gt::cells_body(
            columns = 1,
            rows = tab_est_grf[[1]] == fn_target_label
          )
        )
    }

    # CATE footnote (GRF-specific context)
    if (!is.null(fs$cate_sg)) {
      tab_estimates_grf <- tab_estimates_grf |>
        gt::tab_source_note(
          source_note = gt::md(sprintf(
            "**GRF CATE:** H&#x302; = %.4f, H&#x302;<sup>c</sup> = %.4f (risk-difference scale)",
            fs$cate_sg, fs$cate_sgc
          ))
        )
    }

    return(list(
      tab_estimates = tab_estimates_grf,
      sg10_out      = NULL
    ))
  }

  # =========================================================================
  # FORESTSEARCH PATH (original code below)
  # =========================================================================

  # Select appropriate dataframe
  if (which_df == "est") df <- fs$df.est
  if (which_df == "testing") df <- fs$df.test

  args_fs <- fs$args_call_all

  # Detect GLM outcome
  outcome_type <- args_fs$outcome_type
  if (is.null(outcome_type) || length(outcome_type) > 1L) outcome_type <- "survival"
  is_glm <- outcome_type != "survival"

  # =========================================================================
  # TABLE 1: ITT AND SUBGROUP ESTIMATES
  # =========================================================================

  if (is_glm) {
    # GLM path: use SG_tab_estimates_glm
    effect_measure <- args_fs$effect_measure
    if (is.null(effect_measure)) effect_measure <- "RD"

    # Build estimator closure for the table
    estimator_fn_tab <- tryCatch(
      make_effect_estimator(
        outcome_type      = outcome_type,
        treat.name        = args_fs$treat.name,
        outcome.name      = args_fs$outcome.name,
        offset.name       = args_fs$offset.name,
        effect_measure    = effect_measure,
        adjust_covariates = args_fs$adjust_covariates
      ),
      error = function(e) NULL
    )

    # ITT estimates
    aa <- SG_tab_estimates_glm(
      df, SG_flag = "ITT",
      outcome.name   = args_fs$outcome.name,
      treat.name     = args_fs$treat.name,
      estimator_fn   = estimator_fn_tab,
      effect_measure = effect_measure,
      outcome_type   = outcome_type,
      est.scale      = args_fs$est.scale,
      digits         = ndecimals
    )

    # Subgroup estimates
    bb <- SG_tab_estimates_glm(
      df, SG_flag = "treat.recommend",
      outcome.name   = args_fs$outcome.name,
      treat.name     = args_fs$treat.name,
      estimator_fn   = estimator_fn_tab,
      effect_measure = effect_measure,
      outcome_type   = outcome_type,
      effect_a_1     = hr_1a,
      effect_a_0     = hr_0a,
      sg1_name       = sg1_label,
      sg0_name       = sg0_label,
      est.scale      = args_fs$est.scale,
      digits         = ndecimals
    )

  } else {
    # Survival path: use SG_tab_estimates (unchanged)
    args_tab_filtered <- filter_call_args(
      args_fs,
      SG_tab_estimates,
      list(
        df = df,
        sg0_name = sg0_label,
        sg1_name = sg1_label,
        hr_1a = hr_1a,
        hr_0a = hr_0a
      )
    )

    args_tab_filtered$SG_flag <- "ITT"
    aa <- do.call(SG_tab_estimates, args_tab_filtered)

    args_tab_filtered$SG_flag <- "treat.recommend"
    bb <- do.call(SG_tab_estimates, args_tab_filtered)
  }

  tab_est <- as.data.frame(rbind(aa, bb))

  tab_estimates <- gt::gt(tab_est, auto_align = TRUE) |>
    gt::tab_header(
      title = if (!is.null(est_title)) {
        gt::md(paste0("**", est_title, "**"))
      } else {
        est_caption
      },
      subtitle = if (!is.null(est_title)) est_caption
    ) |>
    gt::tab_options(
      table.font.size = gt::px(font_size),
      heading.title.font.size = gt::px(font_size + 2),
      heading.subtitle.font.size = gt::px(font_size),
      column_labels.font.size = gt::px(font_size)
    )

  # =========================================================================
  # ADD SUBGROUP DEFINITION FOOTNOTE TO tab_estimates
  # =========================================================================

  sg_definition <- NULL
  if (!is.null(fs$grp.consistency$out_sg$sg.harm_label)) {
    sg_definition <- paste(fs$grp.consistency$out_sg$sg.harm_label, collapse = " & ")
  } else if (!is.null(fs$sg.harm)) {
    sg_definition <- paste(fs$sg.harm, collapse = " & ")
  }

  # Add footnote to specific row
  if (!is.null(sg_definition) && nzchar(sg_definition)) {
    fn_word_fs <- if (!is.null(subgroup_notation)) {
      paste0("Identified ", subgroup_notation, " subgroup")
    } else {
      "Identified subgroup"
    }
    tab_estimates <- tab_estimates |>
      gt::tab_footnote(
        footnote = gt::md(paste0("**", fn_word_fs, ":** ", sg_definition)),
        locations = gt::cells_body(
          columns = 1,
          rows = tab_est[[1]] == fn_target_label
        )
      )
  }


  # =========================================================================
  # TABLE 2: IDENTIFIED SUBGROUPS (WITH EXPERIMENTAL ARM EVENTS)
  # =========================================================================

  # Extract from single out_sg object based on sg_focus
  if (!is.null(fs$grp.consistency) && !is.null(fs$grp.consistency$out_sg)) {
    sg10 <- as.data.frame(fs$grp.consistency$out_sg$result)
  } else {
    warning("No subgroup results found in fs$grp.consistency$out_sg")
    return(list(tab_estimates = tab_estimates, sg10_out = NULL))
  }

  # NEW: Calculate experimental arm events (d1) for each subgroup
  outcome.name <- args_fs$outcome.name
  event.name <- args_fs$event.name
  treat.name <- args_fs$treat.name

  # Get d1 for each subgroup in sg10
  sg10$d1 <- NA

  for (i in seq_len(nrow(sg10))) {
    # Get the subgroup ID
    sg_id <- sg10[i, "g"]

    # Look up d1 from the original search results
    if (!is.null(fs$find.grps$out.found$hr.subgroups)) {
      hr_sub <- fs$find.grps$out.found$hr.subgroups
      matching_row <- which(hr_sub$grp == sg_id)
      if (length(matching_row) > 0) {
        sg10$d1[i] <- hr_sub$d1[matching_row[1]]
      }
    }
  }

  # Format based on maxk
  # Resolve sg_subtitle: use user-supplied value or auto-generate from maxk
  if (is.null(sg_subtitle)) {
    sg_subtitle <- switch(
      as.character(args_fs$maxk),
      "1" = "Single-factor subgroups (maxk=1)",
      "2" = "Two-factor subgroups (maxk=2)",
      "3" = "Three-factor subgroups (maxk=3)",
      paste0(args_fs$maxk, "-factor subgroups (maxk=", args_fs$maxk, ")")
    )
  }

  # For GLM, transform the hr column for display and set appropriate labels
  if (is_glm) {
    effect_measure <- args_fs$effect_measure
    if (is.null(effect_measure)) effect_measure <- "RD"

    # For log-scale measures, exponentiate for display
    if (effect_measure %in% c("OR", "RR", "IRR")) {
      sg10$hr <- as.numeric(sg10$hr)
      sg10$hr <- exp(sg10$hr)
    }

    # Column labels for GLM
    hr_label  <- effect_measure
    e_label   <- "n(C+T)"
    d1_label  <- gt::md("n<sub>1</sub>")
  } else {
    # Survival column labels (unchanged)
    hr_label  <- "HR"
    e_label   <- "Events"
    d1_label  <- gt::md("E<sub>1</sub>")

    if (args_fs$est.scale == "1/hr") {
      sg10$hr <- 1 / as.numeric(sg10$hr)
    }
  }

  if (args_fs$maxk == 1) {
    sg10 <- sg10[, c("M.1", "N", "E", "d1", "hr", "Pcons")]

    sg10_out <- sg10 |>
      gt::gt() |>
      gt::fmt_number(columns = c("hr", "Pcons"), decimals = ndecimals) |>
      gt::fmt_number(columns = c("N", "E", "d1"), decimals = 0) |>
      gt::tab_header(
        title = if (!is.null(sg_title)) {
          gt::md(paste0("**", sg_title, "**"))
        } else {
          sg_subtitle
        },
        subtitle = if (!is.null(sg_title)) sg_subtitle
      ) |>
      gt::cols_label(
        M.1 = "Factor 1",
        N = "N",
        E = e_label,
        d1 = d1_label,
        hr = hr_label,
        Pcons = gt::md("P<sub>cons</sub>")
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size),
        heading.title.font.size = gt::px(font_size + 2),
        heading.subtitle.font.size = gt::px(font_size),
        column_labels.font.size = gt::px(font_size),
        source_notes.font.size = gt::px(font_size - 1)
      )

  } else if (args_fs$maxk == 2) {
    sg10 <- sg10[, c("M.1", "M.2", "N", "E", "d1", "hr", "Pcons")]

    sg10_out <- sg10 |>
      gt::gt() |>
      gt::fmt_number(columns = c("hr", "Pcons"), decimals = ndecimals) |>
      gt::fmt_number(columns = c("N", "E", "d1"), decimals = 0) |>
      gt::tab_header(
        title = if (!is.null(sg_title)) {
          gt::md(paste0("**", sg_title, "**"))
        } else {
          sg_subtitle
        },
        subtitle = if (!is.null(sg_title)) sg_subtitle
      ) |>
      gt::cols_label(
        M.1 = "Factor 1",
        M.2 = "Factor 2",
        N = "N",
        E = e_label,
        d1 = d1_label,
        hr = hr_label,
        Pcons = gt::md("P<sub>cons</sub>")
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size),
        heading.title.font.size = gt::px(font_size + 2),
        heading.subtitle.font.size = gt::px(font_size),
        column_labels.font.size = gt::px(font_size),
        source_notes.font.size = gt::px(font_size - 1)
      )

  } else if (args_fs$maxk == 3) {
    sg10 <- sg10[, c("M.1", "M.2", "M.3", "N", "E", "d1", "hr", "Pcons")]

    sg10_out <- sg10 |>
      gt::gt() |>
      gt::fmt_number(columns = c("hr", "Pcons"), decimals = ndecimals) |>
      gt::fmt_number(columns = c("N", "E", "d1"), decimals = 0) |>
      gt::tab_header(
        title = if (!is.null(sg_title)) {
          gt::md(paste0("**", sg_title, "**"))
        } else {
          sg_subtitle
        },
        subtitle = if (!is.null(sg_title)) sg_subtitle
      ) |>
      gt::cols_label(
        M.1 = "Factor 1",
        M.2 = "Factor 2",
        M.3 = "Factor 3",
        N = "N",
        E = e_label,
        d1 = d1_label,
        hr = hr_label,
        Pcons = gt::md("P<sub>cons</sub>")
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size),
        heading.title.font.size = gt::px(font_size + 2),
        heading.subtitle.font.size = gt::px(font_size),
        column_labels.font.size = gt::px(font_size),
        source_notes.font.size = gt::px(font_size - 1)
      )
  }

  # =========================================================================
  # ADD SEARCH METADATA AS FOOTNOTES TO sg10_out
  # =========================================================================

  if (include_search_info && !is.null(fs$find.grps)) {
    # Extract search metadata
    L <- fs$find.grps$L
    max_count <- fs$find.grps$max_count
    maxk <- args_fs$maxk
    n_candidates <- nrow(fs$find.grps$out.found$hr.subgroups)
    max_sg_est <- fs$find.grps$max_sg_est

    # Resolve display label for effect measure
    effect_display <- if (is_glm) {
      em <- args_fs$effect_measure
      if (is.null(em)) "Effect" else em
    } else {
      "HR"
    }
    # For GLM ratio measures, max_sg_est is on log scale; exponentiate
    max_display <- if (is_glm && !is.null(args_fs$effect_measure) &&
                       args_fs$effect_measure %in% c("OR", "RR", "IRR")) {
      round(exp(max_sg_est), 2)
    } else {
      round(max_sg_est, 2)
    }

    # Add search metadata as source notes
    sg10_out <- sg10_out |>
      gt::tab_source_note(
        source_note = gt::md(
          paste0(
            "**Search Configuration:** ",
            "Single-factor candidates (L) = ", L, "; ",
            "Maximum combinations evaluated = ", format(max_count, big.mark = ","), "; ",
            "Search depth (maxk) = ", maxk
          )
        )
      ) |>
      gt::tab_source_note(
        source_note = gt::md(
          paste0(
            "**Search Results:** ",
            "Candidate subgroups found = ", n_candidates, "; ",
            "Maximum ", effect_display, " estimate = ", max_display
          )
        )
      ) |>
      gt::tab_source_note(
        source_note = gt::md(
          "**Note:** E<sub>1</sub> = events in treatment arm; P<sub>cons</sub> = consistency proportion"
        )
      )
  }

  # =========================================================================
  # RETURN ALL TABLES
  # =========================================================================

  result <- list(
    tab_estimates = tab_estimates,
    sg10_out = sg10_out
  )

  return(result)
}


