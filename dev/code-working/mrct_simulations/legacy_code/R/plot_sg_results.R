# =============================================================================
# plot_sg_results.R - Visualization Functions for ForestSearch Subgroup Results
# =============================================================================
#
# This file provides comprehensive plotting functions for visualizing subgroup
# results from ForestSearch, specifically working with fs$df.est output.
#
# Author: ForestSearch Development Team
# License: GPL-3
# =============================================================================


#' Plot ForestSearch Subgroup Results
#'
#' Creates comprehensive visualizations of subgroup results from ForestSearch,
#' including Kaplan-Meier survival curves, hazard ratio comparisons, and
#' summary statistics. This function is designed to work with the output
#' from \code{\link{forestsearch}}, specifically the \code{df.est} component.
#'
#' @param fs.est A forestsearch object or list containing at minimum:
#'   \describe{
#'     \item{df.est}{Data frame with analysis data including \code{treat.recommend}}
#'     \item{sg.harm}{Character vector of subgroup-defining variable names}
#'     \item{grp.consistency}{Optional. Consistency results from \code{sg_consistency_out}}
#'   }
#' @param outcome.name Character. Name of time-to-event outcome column.
#'   Default: "Y"
#' @param event.name Character. Name of event indicator column (1=event, 0=censored).
#'   Default: "Event"
#' @param treat.name Character. Name of treatment column (1=treatment, 0=control).
#'   Default: "Treat"
#' @param plot_type Character. Type of plot to create. One of:
#'   \describe{
#'     \item{"km"}{Kaplan-Meier survival curves}
#'     \item{"forest"}{Forest plot of hazard ratios}
#'     \item{"summary"}{Summary statistics panel}
#'     \item{"combined"}{All plots combined (default)}
#'   }
#' @param by.risk Numeric. Risk interval for KM survival curves.
#'   Default: NULL (auto-calculated)
#' @param conf.level Numeric. Confidence level for intervals. Default: 0.95
#' @param est.scale Character. Effect scale: "hr" (hazard ratio) or
#'   "1/hr" (inverse). Default: "hr"
#' @param sg0_name Character. Label for subgroup 0 (harm/questionable).
#'   Default: "Questionable (H)"
#' @param sg1_name Character. Label for subgroup 1 (recommend/complement).
#'   Default: "Recommend (H^c)"
#' @param treat_labels Named character vector. Labels for treatment arms.
#'   Default: c("0" = "Control", "1" = "Treatment")
#' @param colors Named character vector. Colors for plot elements.
#'   Default: uses package defaults
#' @param title Character. Main plot title. Default: auto-generated
#' @param show_events Logical. Show event counts on KM curves. Default: TRUE
#' @param show_ci Logical. Show confidence intervals. Default: TRUE
#' @param show_logrank Logical. Show log-rank p-value. Default: TRUE
#' @param show_hr Logical. Show hazard ratio annotation. Default: TRUE
#' @param verbose Logical. Print diagnostic messages. Default: FALSE
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return An object of class \code{fs_sg_plot} containing:
#'   \describe{
#'     \item{plots}{List of ggplot2 or base R plot objects}
#'     \item{summary}{Data frame of subgroup summary statistics}
#'     \item{hr_estimates}{Data frame of hazard ratio estimates}
#'     \item{call}{The matched call}
#'   }
#'
#' @details
#' The function extracts subgroup membership from \code{fs.est$df.est$treat.recommend}:
#' \itemize{
#'   \item \code{treat.recommend == 0}: Harm/questionable subgroup (H)
#'   \item \code{treat.recommend == 1}: Recommend/complement subgroup (H^c)
#' }
#'
#' For \code{est.scale = "1/hr"}, treatment labels and subgroup interpretation
#' are reversed to maintain clinical interpretability.
#'
#' @section Kaplan-Meier Plots:
#' When \code{plot_type = "km"}, creates side-by-side survival curves for:
#' \enumerate{
#'   \item The identified subgroup (H) with treatment vs control

#'   \item The complement subgroup (H^c) with treatment vs control
#' }
#'
#' @section Forest Plot:
#' When \code{plot_type = "forest"}, creates a forest plot showing hazard
#' ratios with confidence intervals for: ITT population, H subgroup,
#' and H^c complement.
#'
#' @examples
#' \dontrun{
#' # Run ForestSearch analysis
#' fs <- forestsearch(
#'   df.analysis = my_data,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "treatment",
#'   confounders.name = c("age_cat", "stage", "biomarker")
#' )
#'
#' # Create combined visualization
#' result <- plot_sg_results(
#'   fs.est = fs,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "treatment"
#' )
#'
#' # View the Kaplan-Meier plots only
#' plot_sg_results(fs, plot_type = "km")
#'
#' # Customize labels
#' plot_sg_results(
#'   fs,
#'   sg0_name = "High Risk",
#'   sg1_name = "Standard Risk",
#'   treat_labels = c("0" = "Placebo", "1" = "Active Drug")
#' )
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for running the subgroup analysis
#' \code{\link{sg_consistency_out}} for consistency evaluation
#' \code{\link{plot_subgroup_results_forestplot}} for publication-ready forest plots
#'
#' @importFrom survival survfit Surv coxph survdiff
#' @importFrom stats quantile confint pchisq
#' @importFrom graphics par plot lines legend axis box abline text mtext
#' @importFrom grDevices adjustcolor
#' @export
plot_sg_results <- function(
    fs.est,
    outcome.name = "Y",
    event.name = "Event",
    treat.name = "Treat",
    plot_type = c("combined", "km", "forest", "summary"),
    by.risk = NULL,
    conf.level = 0.95,
    est.scale = c("hr", "1/hr"),
    sg0_name = "Questionable (H)",
    sg1_name = "Recommend (H^c)",
    treat_labels = c("0" = "Control", "1" = "Treatment"),
    colors = NULL,
    title = NULL,
    show_events = TRUE,
    show_ci = TRUE,
    show_logrank = TRUE,
    show_hr = TRUE,
    verbose = FALSE,
    ...
) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION
  # ===========================================================================

  plot_type <- match.arg(plot_type)
  est.scale <- match.arg(est.scale)

  # Validate fs.est structure

  if (!inherits(fs.est, "forestsearch") && !is.list(fs.est)) {
    stop("fs.est must be a forestsearch object or list with df.est component")
  }

  if (is.null(fs.est$df.est)) {
    stop("fs.est$df.est is NULL. No subgroup was identified.")
  }

  df <- fs.est$df.est

  # Check required columns
  required_cols <- c(outcome.name, event.name, treat.name, "treat.recommend")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in df.est: ",
         paste(missing_cols, collapse = ", "))
  }

  # Validate treat.recommend values
  if (!all(df$treat.recommend %in% c(0, 1, NA))) {
    stop("treat.recommend must contain only 0, 1, or NA values")
  }

  if (verbose) {
    message("Input validation passed")
    message(sprintf("  N total: %d", nrow(df)))
    message(sprintf("  treat.recommend == 0: %d", sum(df$treat.recommend == 0, na.rm = TRUE)))
    message(sprintf("  treat.recommend == 1: %d", sum(df$treat.recommend == 1, na.rm = TRUE)))
  }

  # ===========================================================================
  # SECTION 2: SET UP DEFAULT COLORS AND PARAMETERS
  # ===========================================================================

  if (is.null(colors)) {
    colors <- list(
      treat = "#E41A1C",     # Red for treatment
      control = "#377EB8",   # Blue for control
      subgroup_h = "#FF7F00", # Orange for H subgroup
      subgroup_hc = "#4DAF4A", # Green for Hc complement
      ci_fill = "grey80"
    )
  }

  # Auto-calculate by.risk if not provided
  if (is.null(by.risk)) {
    by.risk <- round(max(df[[outcome.name]], na.rm = TRUE) / 12, 0)
    if (by.risk < 1) by.risk <- 1
  }

  # Handle est.scale for treatment reversal
  if (est.scale == "1/hr") {
    df$treat_plot <- 1 - df[[treat.name]]
    treat_labels <- rev(treat_labels)
    temp <- sg0_name
    sg0_name <- sg1_name
    sg1_name <- temp
  } else {
    df$treat_plot <- df[[treat.name]]
  }

  # ===========================================================================
  # SECTION 3: CREATE SUBGROUP DATA FRAMES
  # ===========================================================================

  # H subgroup (treat.recommend == 0): harm/questionable
  df_H <- df[df$treat.recommend == 0 & !is.na(df$treat.recommend), ]

  # Hc subgroup (treat.recommend == 1): complement/recommend
  df_Hc <- df[df$treat.recommend == 1 & !is.na(df$treat.recommend), ]

  if (verbose) {
    message(sprintf("  H subgroup (questionable): n = %d, events = %d",
                    nrow(df_H), sum(df_H[[event.name]])))
    message(sprintf("  Hc subgroup (recommend): n = %d, events = %d",
                    nrow(df_Hc), sum(df_Hc[[event.name]])))
  }

  # ===========================================================================
  # SECTION 4: COMPUTE HAZARD RATIO ESTIMATES
  # ===========================================================================

  hr_estimates <- compute_sg_hr_estimates(
    df = df,
    df_H = df_H,
    df_Hc = df_Hc,
    outcome.name = outcome.name,
    event.name = event.name,
    treat.name = treat.name,
    conf.level = conf.level,
    verbose = verbose
  )

  # ===========================================================================
  # SECTION 5: COMPUTE SUMMARY STATISTICS
  # ===========================================================================

  summary_stats <- compute_sg_summary(
    df = df,
    df_H = df_H,
    df_Hc = df_Hc,
    outcome.name = outcome.name,
    event.name = event.name,
    treat.name = treat.name,
    sg0_name = sg0_name,
    sg1_name = sg1_name
  )

  # ===========================================================================
  # SECTION 6: CREATE PLOTS BASED ON PLOT_TYPE
  # ===========================================================================

  plots <- list()

  if (plot_type %in% c("combined", "km")) {
    plots$km <- plot_sg_km(
      df_H = df_H,
      df_Hc = df_Hc,
      outcome.name = outcome.name,
      event.name = event.name,
      treat.name = treat.name,
      by.risk = by.risk,
      sg0_name = sg0_name,
      sg1_name = sg1_name,
      treat_labels = treat_labels,
      colors = colors,
      show_ci = show_ci,
      show_logrank = show_logrank,
      show_hr = show_hr,
      hr_estimates = hr_estimates,
      conf.level = conf.level,
      title = title,
      ...
    )
  }

  if (plot_type %in% c("combined", "forest")) {
    plots$forest <- plot_sg_forest(
      hr_estimates = hr_estimates,
      sg0_name = sg0_name,
      sg1_name = sg1_name,
      colors = colors,
      title = title,
      ...
    )
  }

  if (plot_type %in% c("combined", "summary")) {
    plots$summary <- plot_sg_summary_panel(
      summary_stats = summary_stats,
      hr_estimates = hr_estimates,
      sg0_name = sg0_name,
      sg1_name = sg1_name,
      colors = colors,
      ...
    )
  }

  # ===========================================================================
  # SECTION 7: RETURN OUTPUT
  # ===========================================================================

  # Get subgroup definition if available
  sg_definition <- NULL
  if (!is.null(fs.est$sg.harm)) {
    sg_definition <- fs.est$sg.harm
  }

  result <- list(
    plots = plots,
    summary = summary_stats,
    hr_estimates = hr_estimates,
    sg_definition = sg_definition,
    call = match.call()
  )

  class(result) <- c("fs_sg_plot", "list")

  return(result)
}


# =============================================================================
# HELPER FUNCTION: Compute Hazard Ratio Estimates
# =============================================================================

#' Compute Hazard Ratio Estimates for Subgroups
#'
#' Internal function to compute Cox model hazard ratio estimates with
#' confidence intervals for ITT, H, and Hc subgroups.
#'
#' @param df Full analysis data frame
#' @param df_H Data frame for H subgroup
#' @param df_Hc Data frame for Hc subgroup
#' @param outcome.name Character. Outcome variable name
#' @param event.name Character. Event indicator name
#' @param treat.name Character. Treatment variable name
#' @param conf.level Numeric. Confidence level
#' @param verbose Logical. Print messages
#'
#' @return Data frame with HR estimates
#' @keywords internal
compute_sg_hr_estimates <- function(
    df,
    df_H,
    df_Hc,
    outcome.name,
    event.name,
    treat.name,
    conf.level = 0.95,
    verbose = FALSE
) {

  # Helper function for safe Cox model fitting
  safe_cox <- function(data, outcome.name, event.name, treat.name) {
    if (nrow(data) < 10 || sum(data[[event.name]]) < 5) {
      return(list(hr = NA, lower = NA, upper = NA, n = nrow(data),
                  events = sum(data[[event.name]])))
    }

    tryCatch({
      formula_str <- sprintf("survival::Surv(%s, %s) ~ %s",
                             outcome.name, event.name, treat.name)
      fit <- survival::coxph(as.formula(formula_str), data = data)
      ci <- confint(fit, level = conf.level)
      list(
        hr = exp(coef(fit)[1]),
        lower = exp(ci[1]),
        upper = exp(ci[2]),
        n = nrow(data),
        events = sum(data[[event.name]])
      )
    }, error = function(e) {
      if (verbose) message("Cox model error: ", e$message)
      list(hr = NA, lower = NA, upper = NA, n = nrow(data),
           events = sum(data[[event.name]]))
    })
  }

  # Compute for each population
  itt_est <- safe_cox(df, outcome.name, event.name, treat.name)
  h_est <- safe_cox(df_H, outcome.name, event.name, treat.name)
  hc_est <- safe_cox(df_Hc, outcome.name, event.name, treat.name)

  # Create output data frame
  data.frame(
    Subgroup = c("ITT", "H", "Hc"),
    n = c(itt_est$n, h_est$n, hc_est$n),
    events = c(itt_est$events, h_est$events, hc_est$events),
    HR = c(itt_est$hr, h_est$hr, hc_est$hr),
    HR_lower = c(itt_est$lower, h_est$lower, hc_est$lower),
    HR_upper = c(itt_est$upper, h_est$upper, hc_est$upper),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# HELPER FUNCTION: Compute Summary Statistics
# =============================================================================

#' Compute Summary Statistics for Subgroups
#'
#' Internal function to compute summary statistics for each subgroup.
#'
#' @param df Full analysis data frame
#' @param df_H Data frame for H subgroup
#' @param df_Hc Data frame for Hc subgroup
#' @param outcome.name Character. Outcome variable name
#' @param event.name Character. Event indicator name
#' @param treat.name Character. Treatment variable name
#' @param sg0_name Character. Label for H subgroup
#' @param sg1_name Character. Label for Hc subgroup
#'
#' @return Data frame with summary statistics
#' @keywords internal
compute_sg_summary <- function(
    df,
    df_H,
    df_Hc,
    outcome.name,
    event.name,
    treat.name,
    sg0_name,
    sg1_name
) {

  # Helper for median survival
  get_median_surv <- function(data, outcome.name, event.name, treat.name) {
    if (nrow(data) < 5) return(c(NA, NA))
    tryCatch({
      formula_str <- sprintf("survival::Surv(%s, %s) ~ %s",
                             outcome.name, event.name, treat.name)
      fit <- survival::survfit(as.formula(formula_str), data = data)
      smry <- summary(fit)$table
      if (is.matrix(smry)) {
        round(smry[, "median"], 1)
      } else {
        round(smry["median"], 1)
      }
    }, error = function(e) c(NA, NA))
  }

  # Compute summaries for each group
  summary_list <- list(
    ITT = list(
      n = nrow(df),
      n_treat = sum(df[[treat.name]] == 1),
      n_control = sum(df[[treat.name]] == 0),
      events = sum(df[[event.name]]),
      events_treat = sum(df[[event.name]][df[[treat.name]] == 1]),
      events_control = sum(df[[event.name]][df[[treat.name]] == 0]),
      median_surv = get_median_surv(df, outcome.name, event.name, treat.name)
    ),
    H = list(
      n = nrow(df_H),
      n_treat = sum(df_H[[treat.name]] == 1),
      n_control = sum(df_H[[treat.name]] == 0),
      events = sum(df_H[[event.name]]),
      events_treat = sum(df_H[[event.name]][df_H[[treat.name]] == 1]),
      events_control = sum(df_H[[event.name]][df_H[[treat.name]] == 0]),
      median_surv = get_median_surv(df_H, outcome.name, event.name, treat.name)
    ),
    Hc = list(
      n = nrow(df_Hc),
      n_treat = sum(df_Hc[[treat.name]] == 1),
      n_control = sum(df_Hc[[treat.name]] == 0),
      events = sum(df_Hc[[event.name]]),
      events_treat = sum(df_Hc[[event.name]][df_Hc[[treat.name]] == 1]),
      events_control = sum(df_Hc[[event.name]][df_Hc[[treat.name]] == 0]),
      median_surv = get_median_surv(df_Hc, outcome.name, event.name, treat.name)
    )
  )

  # Convert to data frame
  data.frame(
    Subgroup = c("ITT", sg0_name, sg1_name),
    n = vapply(summary_list, `[[`, integer(1), "n"),
    n_percent = round(100 * vapply(summary_list, `[[`, integer(1), "n") / nrow(df), 1),
    n_treatment = vapply(summary_list, `[[`, integer(1), "n_treat"),
    n_control = vapply(summary_list, `[[`, integer(1), "n_control"),
    events = vapply(summary_list, `[[`, integer(1), "events"),
    event_rate = round(100 * vapply(summary_list, `[[`, integer(1), "events") /
                         vapply(summary_list, `[[`, integer(1), "n"), 1),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


# =============================================================================
# HELPER FUNCTION: Plot Kaplan-Meier Curves
# =============================================================================

#' Plot Kaplan-Meier Survival Curves for Subgroups
#'
#' Creates side-by-side Kaplan-Meier survival curves for the H and Hc subgroups.
#'
#' @param df_H Data frame for H subgroup
#' @param df_Hc Data frame for Hc subgroup
#' @param outcome.name Character. Outcome variable name
#' @param event.name Character. Event indicator name
#' @param treat.name Character. Treatment variable name
#' @param by.risk Numeric. Risk table interval
#' @param sg0_name Character. Label for H subgroup
#' @param sg1_name Character. Label for Hc subgroup
#' @param treat_labels Named character vector. Treatment labels
#' @param colors List. Color specifications
#' @param show_ci Logical. Show confidence intervals
#' @param show_logrank Logical. Show log-rank p-value
#' @param show_hr Logical. Show HR annotation
#' @param hr_estimates Data frame with HR estimates
#' @param conf.level Numeric. Confidence level
#' @param title Character. Plot title
#' @param ... Additional arguments
#'
#' @return Invisible NULL (creates plot as side effect)
#' @keywords internal
plot_sg_km <- function(
    df_H,
    df_Hc,
    outcome.name,
    event.name,
    treat.name,
    by.risk,
    sg0_name,
    sg1_name,
    treat_labels,
    colors,
    show_ci = TRUE,
    show_logrank = TRUE,
    show_hr = TRUE,
    hr_estimates = NULL,
    conf.level = 0.95,
    title = NULL,
    ...
) {

  # Set up plotting layout
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

  # Helper function to plot a single KM curve
  plot_single_km <- function(data, subtitle, hr_row = NULL) {

    if (nrow(data) < 5) {
      plot.new()
      text(0.5, 0.5, paste("Insufficient data for", subtitle),
           cex = 1.2, col = "gray50")
      return(invisible(NULL))
    }

    # Create survival formula
    formula_str <- sprintf("survival::Surv(%s, %s) ~ %s",
                           outcome.name, event.name, treat.name)
    surv_formula <- as.formula(formula_str)

    # Fit survival curves
    fit <- tryCatch({
      survival::survfit(surv_formula, data = data)
    }, error = function(e) NULL)

    if (is.null(fit)) {
      plot.new()
      text(0.5, 0.5, paste("Could not fit survival curve for", subtitle),
           cex = 1.2, col = "gray50")
      return(invisible(NULL))
    }

    # Extract time range
    max_time <- max(data[[outcome.name]], na.rm = TRUE)

    # Plot
    plot(fit,
         col = c(colors$control, colors$treat),
         lwd = 2,
         xlab = "Time",
         ylab = "Survival Probability",
         main = subtitle,
         conf.int = show_ci,
         mark.time = TRUE,
         ...)

    # Add legend
    legend("bottomleft",
           legend = treat_labels,
           col = c(colors$control, colors$treat),
           lwd = 2,
           bty = "n",
           cex = 0.9)

    # Add log-rank test
    if (show_logrank) {
      lr_test <- tryCatch({
        survival::survdiff(surv_formula, data = data)
      }, error = function(e) NULL)

      if (!is.null(lr_test)) {
        pval <- 1 - pchisq(lr_test$chisq, df = length(lr_test$n) - 1)
        pval_text <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
        mtext(paste("Log-rank:", pval_text), side = 3, line = 0, adj = 1, cex = 0.8)
      }
    }

    # Add HR annotation
    if (show_hr && !is.null(hr_row) && !is.na(hr_row$HR)) {
      hr_text <- sprintf("HR = %.2f (%.2f, %.2f)",
                         hr_row$HR, hr_row$HR_lower, hr_row$HR_upper)
      mtext(hr_text, side = 3, line = -1, adj = 1, cex = 0.8, col = "darkgray")
    }

    invisible(NULL)
  }

  # Get HR estimates for annotation
  hr_H <- if (!is.null(hr_estimates)) hr_estimates[hr_estimates$Subgroup == "H", ] else NULL
  hr_Hc <- if (!is.null(hr_estimates)) hr_estimates[hr_estimates$Subgroup == "Hc", ] else NULL

  # Plot H subgroup
  plot_single_km(
    data = df_H,
    subtitle = paste0(sg0_name, " (n=", nrow(df_H), ")"),
    hr_row = hr_H
  )

  # Plot Hc subgroup
  plot_single_km(
    data = df_Hc,
    subtitle = paste0(sg1_name, " (n=", nrow(df_Hc), ")"),
    hr_row = hr_Hc
  )

  # Add overall title
  if (!is.null(title)) {
    mtext(title, side = 3, outer = TRUE, line = -1.5, cex = 1.2, font = 2)
  }

  invisible(NULL)
}


# =============================================================================
# HELPER FUNCTION: Plot Forest Plot
# =============================================================================

#' Plot Forest Plot of Hazard Ratios
#'
#' Creates a forest plot showing hazard ratios with confidence intervals.
#'
#' @param hr_estimates Data frame with HR estimates
#' @param sg0_name Character. Label for H subgroup
#' @param sg1_name Character. Label for Hc subgroup
#' @param colors List. Color specifications
#' @param title Character. Plot title
#' @param ... Additional arguments
#'
#' @return Invisible NULL (creates plot as side effect)
#' @keywords internal
plot_sg_forest <- function(
    hr_estimates,
    sg0_name,
    sg1_name,
    colors,
    title = NULL,
    ...
) {

  # Set up plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 10, 4, 2) + 0.1)

  n_rows <- nrow(hr_estimates)

  # Define y positions
  y_pos <- n_rows:1

  # Define x limits
  all_vals <- c(hr_estimates$HR, hr_estimates$HR_lower, hr_estimates$HR_upper)
  all_vals <- all_vals[!is.na(all_vals)]
  xlim <- c(min(0.25, min(all_vals) * 0.9), max(2, max(all_vals) * 1.1))

  # Create labels
  labels <- c("ITT Population", sg0_name, sg1_name)

  # Create empty plot
  plot(NULL,
       xlim = xlim,
       ylim = c(0.5, n_rows + 0.5),
       xlab = "Hazard Ratio",
       ylab = "",
       yaxt = "n",
       log = "x",
       main = if (is.null(title)) "Forest Plot: Treatment Effect by Subgroup" else title)

  # Add reference line at HR = 1
  abline(v = 1, lty = 2, col = "gray50")

  # Add points and CIs
  for (i in 1:n_rows) {
    hr <- hr_estimates$HR[i]
    lower <- hr_estimates$HR_lower[i]
    upper <- hr_estimates$HR_upper[i]

    if (!is.na(hr)) {
      # Color by subgroup
      pt_col <- if (i == 1) "black" else if (i == 2) colors$subgroup_h else colors$subgroup_hc

      # Draw CI line
      if (!is.na(lower) && !is.na(upper)) {
        lines(c(lower, upper), c(y_pos[i], y_pos[i]), lwd = 2, col = pt_col)
      }

      # Draw point
      points(hr, y_pos[i], pch = 18, cex = 2, col = pt_col)
    }
  }

  # Add y-axis labels
  axis(2, at = y_pos, labels = labels, las = 2, tick = FALSE)

  # Add HR text annotations on right side
  for (i in 1:n_rows) {
    hr <- hr_estimates$HR[i]
    lower <- hr_estimates$HR_lower[i]
    upper <- hr_estimates$HR_upper[i]
    n <- hr_estimates$n[i]
    events <- hr_estimates$events[i]

    if (!is.na(hr)) {
      hr_text <- sprintf("%.2f (%.2f-%.2f)", hr, lower, upper)
      n_text <- sprintf("n=%d, E=%d", n, events)

      # Add to right of plot
      text(xlim[2] * 0.95, y_pos[i] + 0.15, hr_text, pos = 2, cex = 0.8)
      text(xlim[2] * 0.95, y_pos[i] - 0.15, n_text, pos = 2, cex = 0.7, col = "gray50")
    }
  }

  # Add legend
  legend("topright",
         legend = c("Favors Treatment", "Favors Control"),
         fill = c(adjustcolor("green", 0.3), adjustcolor("red", 0.3)),
         border = NA,
         bty = "n",
         cex = 0.8)

  # Shade regions
  rect(xlim[1], 0.5, 1, n_rows + 0.5, col = adjustcolor("green", 0.1), border = NA)
  rect(1, 0.5, xlim[2], n_rows + 0.5, col = adjustcolor("red", 0.1), border = NA)

  invisible(NULL)
}


# =============================================================================
# HELPER FUNCTION: Plot Summary Panel
# =============================================================================

#' Plot Summary Statistics Panel
#'
#' Creates a summary panel with subgroup characteristics.
#'
#' @param summary_stats Data frame with summary statistics
#' @param hr_estimates Data frame with HR estimates
#' @param sg0_name Character. Label for H subgroup
#' @param sg1_name Character. Label for Hc subgroup
#' @param colors List. Color specifications
#' @param ... Additional arguments
#'
#' @return Invisible NULL (creates plot as side effect)
#' @keywords internal
plot_sg_summary_panel <- function(
    summary_stats,
    hr_estimates,
    sg0_name,
    sg1_name,
    colors,
    ...
) {

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(2, 2, 3, 2))

  # Create empty plot
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  # Title
  text(0.5, 0.95, "Subgroup Summary Statistics", cex = 1.3, font = 2)

  # Create text for each row
  y_start <- 0.85
  y_step <- 0.12

  for (i in 1:nrow(summary_stats)) {
    y_pos <- y_start - (i - 1) * y_step

    # Subgroup name
    text(0.02, y_pos, summary_stats$Subgroup[i], adj = 0, font = 2, cex = 1.0)

    # Sample size
    n_text <- sprintf("n = %d (%.1f%%)", summary_stats$n[i], summary_stats$n_percent[i])
    text(0.35, y_pos, n_text, adj = 0, cex = 0.9)

    # Events
    event_text <- sprintf("Events = %d (%.1f%%)", summary_stats$events[i],
                          summary_stats$event_rate[i])
    text(0.65, y_pos, event_text, adj = 0, cex = 0.9)

    # Add HR if available
    if (i <= nrow(hr_estimates) && !is.na(hr_estimates$HR[i])) {
      hr_text <- sprintf("HR = %.2f", hr_estimates$HR[i])
      text(0.35, y_pos - 0.04, hr_text, adj = 0, cex = 0.8, col = "gray40")
    }
  }

  # Add horizontal lines
  for (i in 1:nrow(summary_stats)) {
    y_pos <- y_start - (i - 0.5) * y_step
    lines(c(0.02, 0.98), c(y_pos, y_pos), col = "gray80")
  }

  box()
  invisible(NULL)
}


# =============================================================================
# PRINT AND PLOT METHODS
# =============================================================================

#' Print Method for fs_sg_plot Objects
#'
#' @param x An fs_sg_plot object
#' @param ... Additional arguments (unused)
#'
#' @export
print.fs_sg_plot <- function(x, ...) {
  cat("ForestSearch Subgroup Visualization\n")
  cat("====================================\n\n")

  if (!is.null(x$sg_definition)) {
    cat("Subgroup Definition:\n")
    cat("  ", paste(x$sg_definition, collapse = " & "), "\n\n")
  }

  cat("Summary Statistics:\n")
  print(x$summary, row.names = FALSE)
  cat("\n")

  cat("Hazard Ratio Estimates:\n")
  hr_display <- x$hr_estimates
  hr_display$HR <- sprintf("%.3f", hr_display$HR)
  hr_display$HR_lower <- sprintf("%.3f", hr_display$HR_lower)
  hr_display$HR_upper <- sprintf("%.3f", hr_display$HR_upper)
  print(hr_display, row.names = FALSE)
  cat("\n")

  cat("Available plots:", paste(names(x$plots), collapse = ", "), "\n")

  invisible(x)
}


#' Plot Method for fs_sg_plot Objects
#'
#' @param x An fs_sg_plot object
#' @param which Character or integer. Which plot to display.
#'   Default: 1 (first available)
#' @param ... Additional arguments passed to plot functions
#'
#' @export
plot.fs_sg_plot <- function(x, which = 1, ...) {
  if (length(x$plots) == 0) {
    message("No plots available")
    return(invisible(x))
  }

  if (is.character(which)) {
    if (!which %in% names(x$plots)) {
      stop("Plot '", which, "' not found. Available: ",
           paste(names(x$plots), collapse = ", "))
    }
    # Re-create the plot
    message("To view plot '", which, "', call plot_sg_results() with plot_type = '", which, "'")
  } else {
    which <- min(which, length(x$plots))
    message("To view plot, call plot_sg_results() with appropriate plot_type")
  }

  invisible(x)
}



