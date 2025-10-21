#' Format Hazard Ratio and Confidence Interval
#'
#' Formats a hazard ratio and confidence interval for display.
#'
#' @param hrest Numeric vector with HR, lower, and upper confidence limits.
#' @return Character string formatted as \"HR (lower, upper)\".
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
#' @export

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
#' @export

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
#' @export

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
#' @export

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

#' RMST calculation for subgroup
#'
#' Calculates restricted mean survival time (RMST) for a subgroup.
#'
#' @param df Data frame.
#' @param tte.name Character. Name of time-to-event variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @return List with tau, RMST, RMST for treatment, RMST for control.
#' @importFrom weightedsurv df_counting
#' @export

rmst_calculation <- function(df,tte.name = "tte",event.name = "event",treat.name = "treat"){
  if (!requireNamespace("weightedsurv", quietly = TRUE)) {
    stop("Package 'weightedsurv' needed for this function to work. Please install it install_github('larry-leon/weightedSurv').")
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




#' Analyze subgroup for summary table (Legacy, to-be-removed)
#'
#' Analyzes a subgroup and returns formatted results for summary table.
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
#' @export

analyze_subgroup_legacy <- function(df_sub, outcome.name, event.name, treat.name, strata.name, subgroup_name, hr_a, potentialOutcome.name, return_medians, N) {
  Y <- df_sub[, outcome.name]
  E <- df_sub[, event.name]
  Treat <- df_sub[, treat.name]
  Strata <- if (is.null(strata.name)) rep("All", length(Y)) else df_sub[, strata.name]

  hr <- cox_summary(Y, E, Treat, Strata)
  meds <- km_summary(Y, E, Treat)
  rmst <- rmst_calculation(df_sub, outcome.name, event.name, treat.name)
  counts <- calculate_counts(Y, E, Treat, N)
  hr_po <- if (!is.null(potentialOutcome.name)) calculate_potential_hr(df_sub, potentialOutcome.name) else NA

  m1 <- if (return_medians) meds[2] else rmst$rmst1
  m0 <- if (return_medians) meds[1] else rmst$rmst0

  format_results(subgroup_name, counts$n, counts$n_treat, counts$d, m1, m0, rmst$rmst, hr, hr_a, hr_po, return_medians)
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
#' with search metadata and customizable decimal precision.
#'
#' @param fs ForestSearch results object.
#' @param which_df Character. Which data frame to use ("est" or "testing").
#' @param est_caption Character. Caption for estimates table.
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param hr_1a Character. Adjusted HR for subgroup 1 (optional).
#' @param hr_0a Character. Adjusted HR for subgroup 0 (optional).
#' @param ndecimals Integer. Number of decimals for formatted numbers (default: 3).
#' @param include_search_info Logical. Include search metadata table (default: TRUE).
#'
#' @return List with gt tables for estimates, subgroups, and optionally search info.
#'
#' @importFrom gt gt fmt_number tab_header tab_spanner tab_source_note md
#' @export
sg_tables <- function(fs,
                      which_df = "est",
                      est_caption = "Training data estimates",
                      potentialOutcome.name = NULL,
                      hr_1a = NA,
                      hr_0a = NA,
                      ndecimals = 3,
                      include_search_info = TRUE) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' required.")
  }

  # Select appropriate dataframe
  if (which_df == "est") df <- fs$df.est
  if (which_df == "testing") df <- fs$df.test

  args_fs <- fs$args_call_all

  # Prepare arguments for SG_tab_estimates
  args_tab_filtered <- filter_call_args(
    args_fs,
    SG_tab_estimates,
    list(
      df = df,
      sg0_name = "Questionable",
      sg1_name = "Recommend",
      hr_1a = hr_1a,
      hr_0a = hr_0a
    )
  )

  # =========================================================================
  # TABLE 1: ITT AND SUBGROUP ESTIMATES
  # =========================================================================

  # ITT estimates
  args_tab_filtered$SG_flag <- "ITT"
  aa <- do.call(SG_tab_estimates, args_tab_filtered)

  # Subgroup estimates
  args_tab_filtered$SG_flag <- "treat.recommend"
  bb <- do.call(SG_tab_estimates, args_tab_filtered)

  tab_est <- as.data.frame(rbind(aa, bb))

  tab_estimates <- gt::gt(tab_est, auto_align = TRUE) |>
    gt::tab_header(
      title = gt::md("**Treatment Effect Estimates**"),
      subtitle = est_caption
    )

  # =========================================================================
  # TABLE 2: IDENTIFIED SUBGROUPS (WITH EXPERIMENTAL ARM EVENTS)
  # =========================================================================

  # Select appropriate result based on sg_focus
  # if (fs$sg_focus == "hr") {
  #   sg10 <- as.data.frame(fs$grp.consistency$out_hr$result)
  # } else if (fs$sg_focus %in% c("minSG", "hrMinSG")) {
  #   sg10 <- as.data.frame(fs$grp.consistency$out_minSG$result)
  # } else if (fs$sg_focus %in% c("maxSG", "hrMaxSG")) {
  #   sg10 <- as.data.frame(fs$grp.consistency$out_maxSG$result)
  # }

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

  for (i in 1:nrow(sg10)) {
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
  if (args_fs$maxk == 1) {
    sg10 <- sg10[, c("M.1", "N", "E", "d1", "hr", "Pcons")]

    if (args_fs$est.scale == "1/hr") {
      sg10$hr <- 1 / sg10$hr
    }

    sg10_out <- sg10 |>
      gt::gt() |>
      gt::fmt_number(columns = c("hr", "Pcons"), decimals = ndecimals) |>
      gt::fmt_number(columns = c("N", "E", "d1"), decimals = 0) |>
      gt::tab_header(
        title = gt::md("**Identified Subgroups**"),
        subtitle = "Single-factor subgroups (maxk=1)"
      ) |>
      gt::cols_label(
        M.1 = "Factor",
        N = "N",
        E = "Events",
        d1 = "E₁",
        hr = "HR",
        Pcons = gt::md("P<sub>cons</sub>")
      )

  } else if (args_fs$maxk == 2) {
    sg10 <- sg10[, c("M.1", "M.2", "N", "E", "d1", "hr", "Pcons")]

    if (args_fs$est.scale == "1/hr") {
      sg10$hr <- 1 / sg10$hr
    }

    sg10_out <- sg10 |>
      gt::gt() |>
      gt::fmt_number(columns = c("hr", "Pcons"), decimals = ndecimals) |>
      gt::fmt_number(columns = c("N", "E", "d1"), decimals = 0) |>
      gt::tab_header(
        title = gt::md("**Identified Subgroups**"),
        subtitle = "Two-factor subgroups (maxk=2)"
      ) |>
      gt::cols_label(
        M.1 = "Factor 1",
        M.2 = "Factor 2",
        N = "N",
        E = "Events",
        d1 = "E₁",
        hr = "HR",
        Pcons = gt::md("P<sub>cons</sub>")
      )

  } else if (args_fs$maxk == 3) {
    sg10 <- sg10[, c("M.1", "M.2", "M.3", "N", "E", "d1", "hr", "Pcons")]

    if (args_fs$est.scale == "1/hr") {
      sg10$hr <- 1 / sg10$hr
    }

    sg10_out <- sg10 |>
      gt::gt() |>
      gt::fmt_number(columns = c("hr", "Pcons"), decimals = ndecimals) |>
      gt::fmt_number(columns = c("N", "E", "d1"), decimals = 0) |>
      gt::tab_header(
        title = gt::md("**Identified Subgroups**"),
        subtitle = "Three-factor subgroups (maxk=3)"
      ) |>
      gt::cols_label(
        M.1 = "Factor 1",
        M.2 = "Factor 2",
        M.3 = "Factor 3",
        N = "N",
        E = "Events",
        d1 = "E₁",
        hr = "HR",
        Pcons = gt::md("P<sub>cons</sub>")
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
            "Maximum HR estimate = ", round(max_sg_est, 2)
          )
        )
      ) |>
      gt::tab_source_note(
        source_note = gt::md(
          "**Note:** E₁ = events in treatment arm; P<sub>cons</sub> = consistency proportion"
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



