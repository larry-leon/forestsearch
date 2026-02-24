#' Format Bootstrap Results Table with gt
#'
#' Creates a publication-ready table from ForestSearch bootstrap results,
#' with bias-corrected confidence intervals, informative formatting, and
#' optional subgroup definition footnote.
#'
#' @param FSsg_tab Data frame or matrix from forestsearch_bootstrap_dofuture()$FSsg_tab
#' @param nb_boots Integer. Number of bootstrap iterations performed
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#' @param boot_success_rate Numeric. Proportion of bootstraps that found subgroups
#' @param sg_definition Character. Subgroup definition string to display as footnote
#'   (e.g., "\{age>=50\} & \{nodes>=3\}"). If NULL, no subgroup footnote is added.
#' @param title Character. Custom title (optional)
#' @param subtitle Character. Custom subtitle (optional)
#'
#' @return A gt table object
#'
#' @importFrom gt gt tab_header tab_spanner tab_footnote tab_source_note md
#'   cols_label tab_style cell_fill cell_text cells_body cells_column_labels
#' @importFrom dplyr all_of
#' @export
format_bootstrap_table <- function(FSsg_tab,
                                   nb_boots,
                                   est.scale = "hr",
                                   boot_success_rate = NULL,
                                   sg_definition = NULL,
                                   title = NULL,
                                   subtitle = NULL) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table formatting.", call. = FALSE)
  }

  # Convert matrix to data frame if needed
  if (is.matrix(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  # Default title and subtitle
  if (is.null(title)) {
    title <- "Treatment Effect by Subgroup"
  }

  if (is.null(subtitle)) {
    subtitle <- sprintf(
      "Bootstrap bias-corrected estimates (%d iterations)",
      nb_boots
    )
  }

  # Column labels configuration
  col_names <- colnames(FSsg_tab)
  labels_list <- list(Subgroup = "Subgroup")

  if ("n" %in% col_names) labels_list$n <- "N"
  if ("n1" %in% col_names) labels_list$n1 <- gt::md("N<sub>T</sub>")
  if ("events" %in% col_names) labels_list$events <- "Events"
  if ("m1" %in% col_names) labels_list$m1 <- gt::md("Med<sub>T</sub>")
  if ("m0" %in% col_names) labels_list$m0 <- gt::md("Med<sub>C</sub>")
  if ("RMST" %in% col_names) labels_list$RMST <- gt::md("RMST<sub>d</sub>")

  # Handle HR column
  hr_col <- grep("HR.*CI", col_names, value = TRUE)[1]
  if (!is.na(hr_col) && length(hr_col) > 0) {
    labels_list[[hr_col]] <- gt::md("HR<br/>(95% CI)<sup>\u2020</sup>")
  }

  # Handle adjusted HR column
  hr_adj_col <- grep("HR\\*", col_names, value = TRUE)[1]
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    labels_list[[hr_adj_col]] <- gt::md("HR<sup>\u2021</sup><br/>(95% CI)")
  }

  # Create base gt table
  tbl <- FSsg_tab |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md(paste0("**", title, "**")),
      subtitle = subtitle
    ) |>
    gt::cols_label(.list = labels_list)

  # Add column spanners
  sample_size_cols <- intersect(c("n", "n1", "events"), col_names)
  if (length(sample_size_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(label = "Sample Size", columns = sample_size_cols)
  }

  survival_cols <- intersect(c("m1", "m0", "RMST"), col_names)
  if (length(survival_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(label = "Survival", columns = survival_cols)
  }

  hr_cols <- grep("HR", col_names, value = TRUE)
  if (length(hr_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(label = "Treatment Effect", columns = gt::starts_with("HR"))
  }

  # Add methodology footnotes
  if (!is.na(hr_col) && length(hr_col) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md(
          "**Unadjusted HR**: Standard Cox regression hazard ratio with robust standard errors"
        ),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_col))
      )
  }

  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md(sprintf(
          "**Bias-corrected HR**: Bootstrap-adjusted estimate using infinitesimal jackknife method (%d iterations). Corrects for optimism in subgroup selection.",
          nb_boots
        )),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_adj_col))
      )
  }

  # Add subgroup definition footnote (analogous to sg_tables tab_estimates)
  if (!is.null(sg_definition) && nzchar(sg_definition)) {
    # Match both abbreviated and full names used in various contexts:
    # - "Qstnbl" (abbreviated, from bootstrap)
    # - "Questionable" (full, from sg_tables)
    # - "H" (shorthand)
    # - "H (Harm)" or similar variants
    h_row_pattern <- "Qstnbl|Questionable|^H$|^H\\s|Harm"

    # Find matching rows in first column
    matching_rows <- grepl(h_row_pattern, FSsg_tab[[1]], ignore.case = TRUE)

    if (any(matching_rows)) {
      tbl <- tbl |>
        gt::tab_footnote(
          footnote = gt::md(paste0("**Identified subgroup**: ", sg_definition)),
          locations = gt::cells_body(
            columns = 1,
            rows = matching_rows
          )
        )
    }
  }

  # Add source note
  source_note_text <- paste0(
    "*Note*: Med = Median survival time (months). ",
    "RMST<sub>d</sub> = Restricted mean survival time difference."
  )

  if (!is.null(boot_success_rate)) {
    source_note_text <- paste0(
      source_note_text,
      sprintf(" Subgroup identified in %.1f%% of bootstrap samples.", boot_success_rate * 100)
    )
  }

  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(source_note_text))

  # Final styling
  tbl <- tbl |>
    gt::tab_options(
      table.border.top.style = "solid",
      table.border.top.width = gt::px(3),
      table.border.top.color = "#333333",
      heading.border.bottom.style = "solid",
      heading.border.bottom.width = gt::px(2),
      heading.border.bottom.color = "#333333",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = gt::px(2),
      column_labels.border.bottom.color = "#333333",
      table.font.size = gt::px(13),
      heading.align = "left"
    )

  return(tbl)
}


#' Create Timing Summary Table
#'
#' Creates a data frame summarizing bootstrap timing information.
#'
#' @param overall_timing List. Overall timing statistics
#' @param iteration_stats List. Per-iteration timing statistics
#' @param fs_stats List. ForestSearch-specific timing statistics
#' @param overhead_stats List. Overhead timing statistics
#' @param nb_boots Integer. Number of bootstrap iterations
#' @param boot_success_rate Numeric. Proportion of successful bootstraps
#'
#' @return Data frame with timing summary
#' @keywords internal
create_timing_summary_table <- function(overall_timing,
                                        iteration_stats,
                                        fs_stats,
                                        overhead_stats,
                                        nb_boots,
                                        boot_success_rate) {

  timing_rows <- list()

  # Overall timing
  if (!is.null(overall_timing) &&
      !is.null(overall_timing$total_minutes) &&
      !is.null(overall_timing$total_hours) &&
      !is.null(overall_timing$avg_minutes_per_boot) &&
      !is.null(overall_timing$avg_seconds_per_boot)) {

    projected_value <- if (nb_boots < 1000) {
      projected_min <- overall_timing$avg_minutes_per_boot * 1000
      sprintf("%.2f min (%.2f hrs)", projected_min, projected_min / 60)
    } else {
      "--"
    }

    timing_rows$overall <- data.frame(
      Category = c("Overall", "", "", ""),
      Metric = c("Total time", "Average per boot", "Successful boots", "Projected for 1000 boots"),
      Value = c(
        sprintf("%.2f min (%.2f hrs)", overall_timing$total_minutes, overall_timing$total_hours),
        sprintf("%.2f min (%.1f sec)", overall_timing$avg_minutes_per_boot, overall_timing$avg_seconds_per_boot),
        sprintf("%d (%.1f%%)", round(boot_success_rate * nb_boots), boot_success_rate * 100),
        projected_value
      ),
      stringsAsFactors = FALSE
    )
  }

  # Per-iteration statistics
  if (!is.null(iteration_stats) &&
      !is.null(iteration_stats$mean) &&
      !is.null(iteration_stats$median) &&
      !is.null(iteration_stats$sd) &&
      !is.null(iteration_stats$min) &&
      !is.null(iteration_stats$max) &&
      !is.null(iteration_stats$q25) &&
      !is.null(iteration_stats$q75)) {

    timing_rows$iteration <- data.frame(
      Category = c("Per-Iteration", "", "", "", "", "", ""),
      Metric = c("Mean", "Median", "Std Dev", "Minimum", "Maximum", "25th percentile", "75th percentile"),
      Value = c(
        sprintf("%.2f min (%.1f sec)", iteration_stats$mean, iteration_stats$mean * 60),
        sprintf("%.2f min (%.1f sec)", iteration_stats$median, iteration_stats$median * 60),
        sprintf("%.2f min", iteration_stats$sd),
        sprintf("%.2f min", iteration_stats$min),
        sprintf("%.2f min", iteration_stats$max),
        sprintf("%.2f min", iteration_stats$q25),
        sprintf("%.2f min", iteration_stats$q75)
      ),
      stringsAsFactors = FALSE
    )
  }

  # ForestSearch statistics
  if (!is.null(fs_stats) &&
      !is.null(fs_stats$n_runs) &&
      !is.null(fs_stats$pct_runs) &&
      !is.null(fs_stats$mean) &&
      !is.null(fs_stats$median) &&
      !is.null(fs_stats$total) &&
      !is.null(overall_timing) &&
      !is.null(overall_timing$total_minutes)) {

    timing_rows$fs <- data.frame(
      Category = c("ForestSearch", "", "", "", ""),
      Metric = c("Runs", "Mean (min)", "Median (min)", "Total (min)", "% of total time"),
      Value = c(
        sprintf("%d (%.1f%%)", fs_stats$n_runs, fs_stats$pct_runs),
        sprintf("%.2f", fs_stats$mean),
        sprintf("%.2f", fs_stats$median),
        sprintf("%.2f", fs_stats$total),
        sprintf("%.1f%%", 100 * fs_stats$total / overall_timing$total_minutes)
      ),
      stringsAsFactors = FALSE
    )
  }

  # Overhead statistics
  if (!is.null(overhead_stats) &&
      !is.null(overhead_stats$mean) &&
      !is.null(overhead_stats$median) &&
      !is.null(overhead_stats$total) &&
      !is.null(overhead_stats$pct_of_total)) {

    timing_rows$overhead <- data.frame(
      Category = c("Overhead", "", "", ""),
      Metric = c("Mean (min)", "Median (min)", "Total (min)", "% of total time"),
      Value = c(
        sprintf("%.2f", overhead_stats$mean),
        sprintf("%.2f", overhead_stats$median),
        sprintf("%.2f", overhead_stats$total),
        sprintf("%.1f%%", overhead_stats$pct_of_total)
      ),
      stringsAsFactors = FALSE
    )
  }

  # Combine all timing rows
  if (length(timing_rows) > 0) {
    timing_summary_table <- do.call(rbind, timing_rows)
    rownames(timing_summary_table) <- NULL
    return(timing_summary_table)
  }

  return(NULL)
}


#' Format Bootstrap Timing Table with gt
#'
#' Creates a publication-ready timing summary table from bootstrap results.
#'
#' @param timing_list List. Timing information from summarize_bootstrap_results()$timing
#' @param nb_boots Integer. Number of bootstrap iterations
#' @param boot_success_rate Numeric. Proportion of successful bootstraps
#'
#' @return A gt table object
#' @importFrom gt gt tab_header tab_spanner tab_footnote tab_source_note md
#'   cols_label tab_style cell_fill cell_text
#' @keywords internal
format_bootstrap_timing_table <- function(timing_list, nb_boots, boot_success_rate) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table formatting.")
  }

  if (is.null(timing_list)) {
    stop("timing_list cannot be NULL")
  }

  # Extract components
  overall <- timing_list$overall
  iteration_stats <- timing_list$iteration_stats
  fs_stats <- timing_list$fs_stats
  overhead_stats <- timing_list$overhead_stats

  # Build the data frame using helper
  timing_df <- create_timing_summary_table(
    overall_timing = overall,
    iteration_stats = iteration_stats,
    fs_stats = fs_stats,
    overhead_stats = overhead_stats,
    nb_boots = nb_boots,
    boot_success_rate = boot_success_rate
  )

  if (is.null(timing_df)) {
    return(NULL)
  }

  # Performance assessment
  if (!is.null(iteration_stats) && !is.null(overall)) {
    avg_sec <- overall$avg_seconds_per_boot

    if (avg_sec < 5) {
      performance <- "Excellent"
      perf_color <- "#d4edda"
    } else if (avg_sec < 15) {
      performance <- "Good"
      perf_color <- "#d1ecf1"
    } else if (avg_sec < 30) {
      performance <- "Acceptable"
      perf_color <- "#fff3cd"
    } else {
      performance <- "Slow"
      perf_color <- "#f8d7da"
    }

    perf_row <- data.frame(
      Category = c("Performance", ""),
      Metric = c("Rating", "Speed"),
      Value = c(performance, sprintf("%.1f sec/iteration", avg_sec)),
      stringsAsFactors = FALSE
    )

    timing_df <- rbind(timing_df, perf_row)
  }

  # Create the gt table
  tbl <- timing_df |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Bootstrap Timing Analysis**"),
      subtitle = sprintf("%d iterations (%.1f%% successful)", nb_boots, boot_success_rate * 100)
    ) |>
    gt::cols_label(
      Category = "Category",
      Metric = "Metric",
      Value = "Value"
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f8f9fa"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = Category,
        rows = Category %in% c("Overall", "Per-Iteration", "ForestSearch", "Overhead", "Performance")
      )
    ) |>
    gt::tab_options(
      table.border.top.style = "solid",
      table.border.top.width = gt::px(3),
      table.border.top.color = "#333333",
      table.font.size = gt::px(13),
      heading.align = "left"
    )

  return(tbl)
}


#' Format Bootstrap Diagnostics Table with gt
#'
#' Creates a publication-ready diagnostics table from bootstrap results.
#'
#' @param diagnostics List. Diagnostics information from summarize_bootstrap_results()
#' @param nb_boots Integer. Number of bootstrap iterations
#' @param results Data.table. Bootstrap results with bias-corrected estimates
#' @param H_estimates List. H subgroup estimates
#' @param Hc_estimates List. Hc subgroup estimates
#'
#' @return A gt table object
#' @importFrom gt gt tab_header tab_footnote tab_source_note md cols_label
#'   tab_style cell_fill cell_text cells_body cells_column_labels
#' @keywords internal
format_bootstrap_diagnostics_table <- function(diagnostics,
                                               nb_boots,
                                               results,
                                               H_estimates = NULL,
                                               Hc_estimates = NULL) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table formatting.")
  }

  # Determine success rating color
  success_rate <- diagnostics$success_rate
  if (success_rate >= 0.9) {
    success_color <- "#d4edda"
    success_rating <- "Excellent"
  } else if (success_rate >= 0.75) {
    success_color <- "#d1ecf1"
    success_rating <- "Good"
  } else if (success_rate >= 0.5) {
    success_color <- "#fff3cd"
    success_rating <- "Moderate"
  } else {
    success_color <- "#f8d7da"
    success_rating <- "Poor"
  }

  # Build diagnostics rows
  diagnostics_rows <- list()

  # Success rate section
  diagnostics_rows$success <- data.frame(
    Category = c("Success Rate", "", "", ""),
    Metric = c("Total iterations", "Successful", "Failed", "Success rating"),
    Value = c(
      as.character(diagnostics$n_boots),
      sprintf("%d (%.1f%%)", diagnostics$n_successful, diagnostics$success_rate * 100),
      sprintf("%d (%.1f%%)", diagnostics$n_failed, (1 - diagnostics$success_rate) * 100),
      success_rating
    ),
    stringsAsFactors = FALSE
  )

  # H estimates section
  if (!is.null(H_estimates)) {
    H_valid <- results$H_biasadj_2[!is.na(results$H_biasadj_2)]
    H_cv <- if (length(H_valid) > 1) sd(H_valid) / abs(mean(H_valid)) * 100 else NA

    diagnostics_rows$H <- data.frame(
      Category = c("Subgroup H (Questionable)", "", "", ""),
      Metric = c("Observed HR", "Bias-corrected HR", "Bootstrap CV (%)", "N estimates"),
      Value = c(
        sprintf("%.3f", H_estimates$H0),
        sprintf("%.3f", H_estimates$H2),
        if (!is.na(H_cv)) sprintf("%.1f%%", H_cv) else "--",
        as.character(length(H_valid))
      ),
      stringsAsFactors = FALSE
    )
  }

  # Hc estimates section
  if (!is.null(Hc_estimates)) {
    Hc_valid <- results$Hc_biasadj_2[!is.na(results$Hc_biasadj_2)]
    Hc_cv <- if (length(Hc_valid) > 1) sd(Hc_valid) / abs(mean(Hc_valid)) * 100 else NA

    diagnostics_rows$Hc <- data.frame(
      Category = c("Subgroup Hc (Recommend)", "", "", ""),
      Metric = c("Observed HR", "Bias-corrected HR", "Bootstrap CV (%)", "N estimates"),
      Value = c(
        sprintf("%.3f", Hc_estimates$H0),
        sprintf("%.3f", Hc_estimates$H2),
        if (!is.na(Hc_cv)) sprintf("%.1f%%", Hc_cv) else "--",
        as.character(length(Hc_valid))
      ),
      stringsAsFactors = FALSE
    )
  }

  # Combine all rows
  diagnostics_df <- do.call(rbind, diagnostics_rows)
  rownames(diagnostics_df) <- NULL

  # Create the gt table
  tbl <- diagnostics_df |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Bootstrap Diagnostics Summary**"),
      subtitle = sprintf("Analysis of %d bootstrap iterations", nb_boots)
    ) |>
    gt::cols_label(
      Category = "Category",
      Metric = "Metric",
      Value = "Value"
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f8f9fa"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = Category,
        rows = Category %in% c("Success Rate", "Subgroup H (Questionable)",
                               "Subgroup Hc (Recommend)")
      )
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = success_color),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c(Metric, Value),
        rows = Metric == "Success rating"
      )
    ) |>
    gt::tab_options(
      table.border.top.style = "solid",
      table.border.top.width = gt::px(3),
      table.border.top.color = "#333333",
      table.font.size = gt::px(13),
      heading.align = "left"
    )

  return(tbl)
}


#' Create Bootstrap Diagnostic Plots
#'
#' Generates diagnostic visualization plots for bootstrap analysis.
#'
#' @param results Data frame with bootstrap results
#' @param H_estimates List with H subgroup estimates
#' @param Hc_estimates List with Hc subgroup estimates
#' @param overall_timing List with overall timing information (optional)
#'
#' @return List of ggplot2 objects
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_minimal
#'   theme element_text
#' @importFrom rlang .data
#' @keywords internal
create_bootstrap_diagnostic_plots <- function(results,
                                              H_estimates,
                                              Hc_estimates,
                                              overall_timing = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }

  plots <- list()

  # Plot 1: Bootstrap distribution of bias-corrected estimates H
  if ("H_biasadj_2" %in% names(results) && !is.null(H_estimates)) {
    valid_H <- results[!is.na(results$H_biasadj_2), ]

    if (nrow(valid_H) > 0) {
      p1 <- ggplot2::ggplot(valid_H, ggplot2::aes(x = .data$H_biasadj_2)) +
        ggplot2::geom_histogram(bins = 30, fill = "#4472C4", alpha = 0.7, color = "white") +
        ggplot2::geom_vline(
          xintercept = c(log(H_estimates$H2), log(H_estimates$H0)),
          color = c("green", "red"),
          linetype = "dashed",
          linewidth = 1
        ) +
        ggplot2::labs(
          title = "Bootstrap Distribution: Subgroup H",
          subtitle = "Green = bias-corrected, Red = observed",
          x = "Log Hazard Ratio (bias-corrected)",
          y = "Frequency"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
        )

      plots$H_distribution <- p1
    }
  }

  # Plot 2: Bootstrap distribution of bias-corrected estimates Hc
  if ("Hc_biasadj_2" %in% names(results) && !is.null(Hc_estimates)) {
    valid_Hc <- results[!is.na(results$Hc_biasadj_2), ]

    if (nrow(valid_Hc) > 0) {
      p2 <- ggplot2::ggplot(valid_Hc, ggplot2::aes(x = .data$Hc_biasadj_2)) +
        ggplot2::geom_histogram(bins = 30, fill = "#70AD47", alpha = 0.7, color = "white") +
        ggplot2::geom_vline(
          xintercept = c(log(Hc_estimates$H2), log(Hc_estimates$H0)),
          color = c("green", "red"),
          linetype = "dashed",
          linewidth = 1
        ) +
        ggplot2::labs(
          title = "Bootstrap Distribution: Subgroup Hc",
          subtitle = "Green = bias-corrected, Red = observed",
          x = "Log Hazard Ratio (bias-corrected)",
          y = "Frequency"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
        )

      plots$Hc_distribution <- p2
    }
  }

  # Plot 3: Timing distribution (if available)
  if ("tmins_iteration" %in% names(results)) {
    valid_timing <- results[!is.na(results$tmins_iteration), ]

    if (nrow(valid_timing) > 0) {
      p3 <- ggplot2::ggplot(valid_timing, ggplot2::aes(x = .data$tmins_iteration)) +
        ggplot2::geom_histogram(bins = 30, fill = "#2E86AB", alpha = 0.7, color = "white") +
        ggplot2::geom_vline(
          xintercept = median(valid_timing$tmins_iteration),
          color = "red",
          linetype = "dashed",
          linewidth = 1
        ) +
        ggplot2::labs(
          title = "Bootstrap Iteration Timing Distribution",
          subtitle = sprintf(
            "Median: %.2f min, Mean: %.2f min",
            median(valid_timing$tmins_iteration),
            mean(valid_timing$tmins_iteration)
          ),
          x = "Minutes per iteration",
          y = "Frequency"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
        )

      plots$timing_distribution <- p3
    }
  }

  return(plots)
}


#' Calculate Bootstrap Table Caption
#'
#' Generates an interpretive caption for bootstrap results table.
#'
#' @param est.scale Character. "hr" or "1/hr"
#' @param nb_boots Integer. Number of bootstrap iterations
#' @param boot_success_rate Numeric. Proportion successful
#'
#' @return Character string with caption
#' @keywords internal
create_bootstrap_caption <- function(est.scale, nb_boots, boot_success_rate) {
  effect_name <- ifelse(est.scale == "hr", "Hazard Ratios", "Inverse Hazard Ratios")

  sprintf(
    paste0(
      "%s less than 1.0 indicate benefit from treatment. ",
      "Bias correction typically shifts estimates toward the null, ",
      "reflecting the optimism inherent in data-driven subgroup identification. ",
      "Bootstrap analysis based on %d iterations (%.1f%% successful)."
    ),
    effect_name,
    nb_boots,
    boot_success_rate * 100
  )
}


#' Calculate Skewness
#'
#' Helper function to calculate sample skewness.
#'
#' @param x Numeric vector
#' @return Numeric skewness value
#' @importFrom stats sd
#' @keywords internal
calculate_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA)

  x_bar <- mean(x)
  s <- sd(x)

  if (s == 0) return(0)

  skew <- (n / ((n - 1) * (n - 2))) * sum(((x - x_bar) / s)^3)
  return(skew)
}


#' Summarize Factor Presence Across Bootstrap Subgroups
#'
#' Analyzes how often each individual factor appears in identified subgroups,
#' extracting base factor names from full definitions and identifying common
#' specific definitions.
#'
#' @param results Data.table or data.frame. Bootstrap results with M.1, M.2, etc. columns
#' @param maxk Integer. Maximum number of factors allowed
#' @param threshold Numeric. Percentage threshold for including specific definitions (default: 10)
#' @param as_gt Logical. Return gt tables (TRUE) or data.frames (FALSE)
#'
#' @return List with base_factors and specific_factors data.frames or gt tables
#' @keywords internal
summarize_factor_presence_robust <- function(results,
                                             maxk = 2,
                                             threshold = 10,
                                             as_gt = TRUE) {

  # Input validation
 if (!is.data.frame(results) && !inherits(results, "data.table")) {
    stop("Input 'results' must be a data.frame or data.table")
  }

  results <- as.data.frame(results, stringsAsFactors = FALSE)

  # Filter to successful iterations
  if ("Pcons" %in% names(results)) {
    sg_found <- results[!is.na(results$Pcons), , drop = FALSE]
  } else {
    sg_found <- results
  }

  n_found <- nrow(sg_found)

  # Handle empty results
  if (n_found == 0) {
    base_factor_summary <- data.frame(
      Rank = integer(),
      Factor = character(),
      Count = integer(),
      Percent = numeric(),
      stringsAsFactors = FALSE
    )
    specific_factor_summary <- data.frame(
      Rank = integer(),
      Base_Factor = character(),
      Factor_Definition = character(),
      Count = integer(),
      Percent = numeric(),
      stringsAsFactors = FALSE
    )

    if (as_gt && requireNamespace("gt", quietly = TRUE)) {
      base_factor_summary <- gt::gt(base_factor_summary)
      specific_factor_summary <- gt::gt(specific_factor_summary)
    }

    return(list(
      base_factors = base_factor_summary,
      specific_factors = specific_factor_summary,
      threshold = threshold
    ))
  }

  # Collect all factors from M.1, M.2, etc. columns
  all_factors <- character()
  for (i in seq_len(maxk)) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      vals <- sg_found[[col]]
      if (is.factor(vals)) vals <- as.character(vals)
      vals <- vals[!is.na(vals) & vals != ""]
      all_factors <- c(all_factors, vals)
    }
  }

  # Handle no factors found
  if (length(all_factors) == 0) {
    base_factor_summary <- data.frame(
      Rank = integer(),
      Factor = character(),
      Count = integer(),
      Percent = numeric(),
      stringsAsFactors = FALSE
    )
    specific_factor_summary <- data.frame(
      Rank = integer(),
      Base_Factor = character(),
      Factor_Definition = character(),
      Count = integer(),
      Percent = numeric(),
      stringsAsFactors = FALSE
    )

    if (as_gt && requireNamespace("gt", quietly = TRUE)) {
      base_factor_summary <- gt::gt(base_factor_summary)
      specific_factor_summary <- gt::gt(specific_factor_summary)
    }

    return(list(
      base_factors = base_factor_summary,
      specific_factors = specific_factor_summary,
      threshold = threshold
    ))
  }

  # Helper function to extract base factor name
  extract_base_factor <- function(factor_string) {
    factor_string <- trimws(factor_string)
    # Remove leading negation and braces
    factor_string <- gsub("^!\\{", "{", factor_string)
    factor_string <- gsub("[{}]", "", factor_string)
    # Extract base name (before any operator like >=, <=, ==)
    base_name <- sub("^([a-zA-Z0-9_.]+).*", "\\1", factor_string)
    return(base_name)
  }

  # Extract base factors and count
  base_factors <- vapply(all_factors, extract_base_factor, character(1),
                         USE.NAMES = FALSE)

  base_factor_counts <- table(base_factors)

  base_factor_summary <- data.frame(
    Factor = names(base_factor_counts),
    Count = as.integer(base_factor_counts),
    Percent = 100 * as.integer(base_factor_counts) / n_found,
    stringsAsFactors = FALSE
  )
  base_factor_summary <- base_factor_summary[order(-base_factor_summary$Count), ]
  base_factor_summary$Rank <- seq_len(nrow(base_factor_summary))
  base_factor_summary <- base_factor_summary[, c("Rank", "Factor", "Count", "Percent")]

  # Count specific factor definitions
  specific_factor_counts <- table(all_factors)

  specific_factor_summary <- data.frame(
    Factor_Definition = names(specific_factor_counts),
    Count = as.integer(specific_factor_counts),
    Percent = 100 * as.integer(specific_factor_counts) / n_found,
    stringsAsFactors = FALSE
  )

  # Filter by threshold
  specific_factor_summary <- specific_factor_summary[
    specific_factor_summary$Percent >= threshold,
  ]

  if (nrow(specific_factor_summary) > 0) {
    specific_factor_summary$Base_Factor <- vapply(
      specific_factor_summary$Factor_Definition,
      extract_base_factor,
      character(1)
    )
    specific_factor_summary <- specific_factor_summary[
      order(specific_factor_summary$Base_Factor, -specific_factor_summary$Count),
    ]
    specific_factor_summary$Rank <- seq_len(nrow(specific_factor_summary))
    specific_factor_summary <- specific_factor_summary[
      , c("Rank", "Base_Factor", "Factor_Definition", "Count", "Percent")
    ]
  } else {
    specific_factor_summary <- data.frame(
      Rank = integer(),
      Base_Factor = character(),
      Factor_Definition = character(),
      Count = integer(),
      Percent = numeric(),
      stringsAsFactors = FALSE
    )
  }

  # Convert to gt if requested
  if (as_gt && requireNamespace("gt", quietly = TRUE)) {
    base_factor_summary <- gt::gt(base_factor_summary)
    specific_factor_summary <- gt::gt(specific_factor_summary)
  }

  return(list(
    base_factors = base_factor_summary,
    specific_factors = specific_factor_summary,
    threshold = threshold
  ))
}


#' Summarize Bootstrap Event Counts
#'
#' Provides summary statistics for event counts across bootstrap iterations,
#' helping assess the reliability of HR estimates when events are sparse.
#'
#' @param boot_results List. Output from forestsearch_bootstrap_dofuture()
#' @param threshold Integer. Minimum event threshold for flagging low counts (default: 5)
#'
#' @return Invisibly returns a list with summary statistics:
#'   \describe{
#'     \item{threshold}{The event threshold used}
#'     \item{nb_boots}{Total number of bootstrap iterations}
#'     \item{n_successful}{Number of iterations that found a new subgroup}
#'     \item{original_H}{List with low event counts for original H on bootstrap samples}
#'     \item{original_Hc}{List with low event counts for original Hc on bootstrap samples}
#'     \item{new_Hstar}{List with low event counts for new H* on original data}
#'     \item{new_Hcstar}{List with low event counts for new Hc* on original data}
#'   }
#'
#' @details
#' This function summarizes event counts in four scenarios:
#' \enumerate{
#'   \item ORIGINAL subgroup H evaluated on BOOTSTRAP samples
#'   \item ORIGINAL subgroup Hc evaluated on BOOTSTRAP samples
#'   \item NEW subgroup H* (found in bootstrap) evaluated on ORIGINAL data
#'   \item NEW subgroup Hc* (found in bootstrap) evaluated on ORIGINAL data
#' }
#'
#' Low event counts (below threshold) can lead to unstable HR estimates.
#' This summary helps identify potential issues with sparse events.
#'
#' @examples
#' \dontrun{
#' # After running bootstrap analysis
#' summarize_bootstrap_events(fs_bc, threshold = 10)
#' }
#'
#' @export
summarize_bootstrap_events <- function(boot_results, threshold = 5) {

  results <- boot_results$results
  nb_boots <- nrow(results)

  # Check low events for ORIGINAL subgroup on BOOTSTRAP samples
  low_H_0 <- sum(results$events_H_0 < threshold, na.rm = TRUE)
  low_H_1 <- sum(results$events_H_1 < threshold, na.rm = TRUE)
  low_H_either <- sum(
    results$events_H_0 < threshold | results$events_H_1 < threshold,
    na.rm = TRUE
  )

  low_Hc_0 <- sum(results$events_Hc_0 < threshold, na.rm = TRUE)
  low_Hc_1 <- sum(results$events_Hc_1 < threshold, na.rm = TRUE)
  low_Hc_either <- sum(
    results$events_Hc_0 < threshold | results$events_Hc_1 < threshold,
    na.rm = TRUE
  )

  # Check low events for NEW subgroup on ORIGINAL data
  low_Hstar_0 <- sum(results$events_Hstar_0 < threshold, na.rm = TRUE)
  low_Hstar_1 <- sum(results$events_Hstar_1 < threshold, na.rm = TRUE)
  low_Hstar_either <- sum(
    results$events_Hstar_0 < threshold | results$events_Hstar_1 < threshold,
    na.rm = TRUE
  )

  low_Hcstar_0 <- sum(results$events_Hcstar_0 < threshold, na.rm = TRUE)
  low_Hcstar_1 <- sum(results$events_Hcstar_1 < threshold, na.rm = TRUE)
  low_Hcstar_either <- sum(
    results$events_Hcstar_0 < threshold | results$events_Hcstar_1 < threshold,
    na.rm = TRUE
  )

  # Print summary statistics
  cat("\n=== Bootstrap Event Count Summary ===\n")
  cat(sprintf("Total bootstrap iterations: %d\n", nb_boots))
  cat(sprintf("Event threshold: <%d events\n\n", threshold))

  cat("ORIGINAL Subgroup H on BOOTSTRAP samples:\n")
  cat(sprintf("  Control arm <%d events: %d (%.1f%%)\n",
              threshold, low_H_0, 100 * low_H_0 / nb_boots))
  cat(sprintf("  Treatment arm <%d events: %d (%.1f%%)\n",
              threshold, low_H_1, 100 * low_H_1 / nb_boots))
  cat(sprintf("  Either arm <%d events: %d (%.1f%%)\n\n",
              threshold, low_H_either, 100 * low_H_either / nb_boots))

  cat("ORIGINAL Subgroup Hc on BOOTSTRAP samples:\n")
  cat(sprintf("  Control arm <%d events: %d (%.1f%%)\n",
              threshold, low_Hc_0, 100 * low_Hc_0 / nb_boots))
  cat(sprintf("  Treatment arm <%d events: %d (%.1f%%)\n",
              threshold, low_Hc_1, 100 * low_Hc_1 / nb_boots))
  cat(sprintf("  Either arm <%d events: %d (%.1f%%)\n\n",
              threshold, low_Hc_either, 100 * low_Hc_either / nb_boots))

  # Count successful bootstraps (found new subgroup)
  n_successful <- sum(!is.na(results$events_Hstar_0))

  if (n_successful > 0) {
    cat(sprintf("NEW Subgroups found: %d (%.1f%%)\n\n",
                n_successful, 100 * n_successful / nb_boots))

    cat("NEW Subgroup H* on ORIGINAL data:\n")
    cat(sprintf("  Control arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hstar_0, 100 * low_Hstar_0 / n_successful))
    cat(sprintf("  Treatment arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hstar_1, 100 * low_Hstar_1 / n_successful))
    cat(sprintf("  Either arm <%d events: %d (%.1f%% of successful)\n\n",
                threshold, low_Hstar_either, 100 * low_Hstar_either / n_successful))

    cat("NEW Subgroup Hc* on ORIGINAL data:\n")
    cat(sprintf("  Control arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_0, 100 * low_Hcstar_0 / n_successful))
    cat(sprintf("  Treatment arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_1, 100 * low_Hcstar_1 / n_successful))
    cat(sprintf("  Either arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_either, 100 * low_Hcstar_either / n_successful))
  } else {
    cat("No new subgroups found in bootstrap samples\n")
  }

  cat("\n")

  # Return invisible list with all counts
  invisible(list(
    threshold = threshold,
    nb_boots = nb_boots,
    n_successful = n_successful,
    original_H = list(
      low_control = low_H_0,
      low_treat = low_H_1,
      low_either = low_H_either
    ),
    original_Hc = list(
      low_control = low_Hc_0,
      low_treat = low_Hc_1,
      low_either = low_Hc_either
    ),
    new_Hstar = list(
      low_control = low_Hstar_0,
      low_treat = low_Hstar_1,
      low_either = low_Hstar_either
    ),
    new_Hcstar = list(
      low_control = low_Hcstar_0,
      low_treat = low_Hcstar_1,
      low_either = low_Hcstar_either
    )
  ))
}
