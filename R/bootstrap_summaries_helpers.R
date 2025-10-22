#' Format Bootstrap Results Table with gt
#'
#' Creates a publication-ready table from ForestSearch bootstrap results,
#' with bias-corrected confidence intervals and informative formatting.
#'
#' @param FSsg_tab Data frame or matrix from forestsearch_bootstrap_dofuture()$FSsg_tab
#' @param nb_boots Integer. Number of bootstrap iterations performed
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#' @param boot_success_rate Numeric. Proportion of bootstraps that found subgroups
#' @param title Character. Custom title (optional)
#' @param subtitle Character. Custom subtitle (optional)
#'
#' @return A gt table object
#' @importFrom gt gt tab_header tab_spanner tab_footnote md cols_label
#' @export

format_bootstrap_table <- function(FSsg_tab, nb_boots, est.scale = "hr",
                                   boot_success_rate = NULL,
                                   title = NULL, subtitle = NULL) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table formatting. Install with: install.packages('gt')")
  }

  # CRITICAL FIX: Convert matrix to data frame if needed
  if (is.matrix(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  # Ensure it's a data frame
  if (!is.data.frame(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  # Default title and subtitle
  if (is.null(title)) {
    title <- "Treatment Effect by Subgroup"
  }

  if (is.null(subtitle)) {
    effect_label <- ifelse(est.scale == "hr", "Hazard Ratio", "Inverse Hazard Ratio (1/HR)")
    subtitle <- sprintf("Bootstrap bias-corrected estimates (%d iterations)", nb_boots)
  }

  # Get column names (handle different possible names)
  col_names <- colnames(FSsg_tab)

  # Create base labels list
  labels_list <- list(
    Subgroup = "Subgroup"
  )

  # Add labels for columns that exist
  if ("n" %in% col_names) labels_list$n <- "N"
  if ("n1" %in% col_names) labels_list$n1 <- gt::md("N<sub>treat</sub>")
  if ("events" %in% col_names) labels_list$events <- "Events"
  if ("m1" %in% col_names) labels_list$m1 <- gt::md("Med<sub>treat</sub>")
  if ("m0" %in% col_names) labels_list$m0 <- gt::md("Med<sub>ctrl</sub>")
  if ("RMST" %in% col_names) labels_list$RMST <- gt::md("RMST<sub>diff</sub>")

  # Handle HR column (might be "HR (95% CI)" or similar)
  hr_col <- grep("HR.*CI", col_names, value = TRUE)[1]
  if (!is.na(hr_col) && length(hr_col) > 0) {
    labels_list[[hr_col]] <- gt::md("HR<br/>(95% CI)<sup>†</sup>")
  }

  # Handle adjusted HR column (might be "HR*" or similar)
  hr_adj_col <- grep("HR\\*", col_names, value = TRUE)[1]
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    labels_list[[hr_adj_col]] <- gt::md("HR<sup>‡</sup><br/>(95% CI)")
  }

  # Create the gt table
  tbl <- FSsg_tab |>
    gt::gt() |>

    # Title and subtitle
    gt::tab_header(
      title = gt::md(paste0("**", title, "**")),
      subtitle = subtitle
    ) |>

    # Column labels with better formatting
    gt::cols_label(.list = labels_list)

  # Add spanners if columns exist
  sample_size_cols <- intersect(c("n", "n1", "events"), col_names)
  if (length(sample_size_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Sample Size",
        columns = sample_size_cols
      )
  }

  survival_cols <- intersect(c("m1", "m0", "RMST"), col_names)
  if (length(survival_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Survival",
        columns = survival_cols
      )
  }

  hr_cols <- grep("HR", col_names, value = TRUE)
  if (length(hr_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Treatment Effect",
        columns = gt::starts_with("HR")
      )
  }

  # Add footnotes
  if (!is.na(hr_col) && length(hr_col) > 0) {
    # Use dplyr::all_of() which comes with gt, or just use the column directly
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md("**Unadjusted HR**: Standard Cox regression hazard ratio with robust standard errors"),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_col))
      )
  }

  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md(sprintf(
          "**Bias-corrected HR**: Bootstrap-adjusted estimate using infinitesimal jacknife method (%d iterations). Corrects for optimism in subgroup selection.",
          nb_boots
        )),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_adj_col))
      )
  }

  # Source note with technical details
  source_note_text <- "*Note*: Med = Median survival time (months). RMST<sub>diff</sub> = Restricted mean survival time difference."
  if (!is.null(boot_success_rate)) {
    source_note_text <- paste0(
      source_note_text,
      sprintf(" Subgroup identified in %.1f%% of bootstrap samples.", boot_success_rate * 100)
    )
  }

  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(source_note_text))

  # Styling
  tbl <- tbl |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f0f0f0"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_column_labels()
    )

  # Center align numeric columns if they exist
  numeric_cols <- intersect(c("n", "n1", "events", "m1", "m0", "RMST"), col_names)
  if (length(numeric_cols) > 0) {
    tbl <- tbl |>
      gt::tab_style(
        style = gt::cell_text(align = "center"),
        locations = gt::cells_body(columns = dplyr::all_of(numeric_cols))
      )
  }

  # Highlight the bias-corrected estimates if column exists
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    tbl <- tbl |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(columns = dplyr::all_of(hr_adj_col))
      )
  }

  # Add subtle borders
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
      table.font.size = gt::px(14)
    )

  return(tbl)
}



create_bootstrap_diagnostic_plots <- function(results, H_estimates, Hc_estimates) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }

  # Plot 1: Bootstrap distribution of bias-corrected estimates H
  p1 <- ggplot2::ggplot(results, ggplot2::aes(x = H_biasadj_2)) +
    ggplot2::geom_histogram(bins = 30, fill = "#4472C4", alpha = 0.7, color = "white") +
    ggplot2::geom_vline(xintercept = c(log(H_estimates$H2),log(H_estimates$H0)), color = c("green","red"),
                        linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "Bootstrap Distribution of Bias-Corrected
      log(hazard ratio): Subgroup H",
      subtitle = "Red line shows final bias-corrected estimate",
      x = "Log Hazard Ratio (bias-corrected)",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
    )

  # Plot 2: Bootstrap distribution of bias-corrected estimates Hc
  p2 <- ggplot2::ggplot(results, ggplot2::aes(x = Hc_biasadj_2)) +
    ggplot2::geom_histogram(bins = 30, fill = "#4472C4", alpha = 0.7, color = "white") +
    ggplot2::geom_vline(xintercept = c(log(Hc_estimates$H2),log(Hc_estimates$H0)), color = c("green","red"),
                        linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "Bootstrap Distribution of Bias-Corrected
      log(hazard ratio): Subgroup Hc",
      subtitle = "Red line shows final bias-corrected estimate",
      x = "Log Hazard Ratio (bias-corrected)",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
    )



  # Plot 3: Bias correction impact
  bias_data <- data.frame(
    iteration = 1:nrow(results),
    unadjusted = H_estimates$H0,
    adjusted = results$H_biasadj_2
  )

  p3 <- ggplot2::ggplot(results, ggplot2::aes(x = H_biasadj_1, y = H_biasadj_2)) +
    ggplot2::geom_point(alpha = 0.5, color = "#4472C4") +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Comparison of H Bias Correction Methods",
      x = "Method 1 (Simple optimism correction)",
      y = "Method 2 (Double bootstrap correction)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  p4 <- ggplot2::ggplot(results, ggplot2::aes(x = Hc_biasadj_1, y = Hc_biasadj_2)) +
    ggplot2::geom_point(alpha = 0.5, color = "#4472C4") +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Comparison of Hc Bias Correction Methods",
      x = "Method 1 (Simple optimism correction)",
      y = "Method 2 (Double bootstrap correction)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )



  list(
    H_distribution = p1,
    Hc_distribution = p2,
    H_methods_comparison = p3,
    Hc_methods_comparison = p4
  )
}


#' Simple Caption Generator for Bootstrap Table
#'
#' Generates an informative caption based on bootstrap results
#'
#' @param nb_boots Integer. Number of bootstrap iterations
#' @param boot_success_rate Numeric. Proportion of successful bootstraps
#' @param est.scale Character. "hr" or "1/hr"
#'
#' @return Character string with formatted caption
#' @export

generate_bootstrap_caption <- function(nb_boots, boot_success_rate, est.scale = "hr") {

  effect_name <- ifelse(est.scale == "hr", "hazard ratio", "inverse hazard ratio (1/HR)")

  caption <- sprintf(
    "Treatment effect estimates by identified subgroups. The unadjusted %s (HR) is from standard Cox regression with robust standard errors. The bias-corrected estimate (HR*) adjusts for optimism in subgroup selection using bootstrap methods (%d iterations, %.1f%% successful). Values less than 1.0 indicate benefit from treatment; bias correction typically shifts estimates toward the null, reflecting the optimism inherent in data-driven subgroup identification.",
    effect_name,
    nb_boots,
    boot_success_rate * 100
  )

  return(caption)
}


#' Summarize Bootstrap Event Counts
#'
#' Provides summary statistics for event counts across bootstrap iterations
#'
#' @param boot_results Data.table from bootstrap_results()
#' @param threshold Integer. Minimum event threshold (default: 5)
#' @return List with summary statistics
#' @export

summarize_bootstrap_events <- function(boot_results, threshold = 5) {

  results <- boot_results$results
  nb_boots <- nrow(results)

  # Check low events for ORIGINAL subgroup on BOOTSTRAP samples
  low_H_0 <- sum(results$events_H_0 < threshold, na.rm = TRUE)
  low_H_1 <- sum(results$events_H_1 < threshold, na.rm = TRUE)
  low_H_either <- sum(results$events_H_0 < threshold | results$events_H_1 < threshold, na.rm = TRUE)

  low_Hc_0 <- sum(results$events_Hc_0 < threshold, na.rm = TRUE)
  low_Hc_1 <- sum(results$events_Hc_1 < threshold, na.rm = TRUE)
  low_Hc_either <- sum(results$events_Hc_0 < threshold | results$events_Hc_1 < threshold, na.rm = TRUE)

  # Check low events for NEW subgroup on ORIGINAL data
  low_Hstar_0 <- sum(results$events_Hstar_0 < threshold, na.rm = TRUE)
  low_Hstar_1 <- sum(results$events_Hstar_1 < threshold, na.rm = TRUE)
  low_Hstar_either <- sum(results$events_Hstar_0 < threshold | results$events_Hstar_1 < threshold, na.rm = TRUE)

  low_Hcstar_0 <- sum(results$events_Hcstar_0 < threshold, na.rm = TRUE)
  low_Hcstar_1 <- sum(results$events_Hcstar_1 < threshold, na.rm = TRUE)
  low_Hcstar_either <- sum(results$events_Hcstar_0 < threshold | results$events_Hcstar_1 < threshold, na.rm = TRUE)

  # Summary statistics
  cat("\n=== Bootstrap Event Count Summary ===\n")
  cat(sprintf("Total bootstrap iterations: %d\n", nb_boots))
  cat(sprintf("Event threshold: <%d events\n\n", threshold))

  cat("ORIGINAL Subgroup H on BOOTSTRAP samples:\n")
  cat(sprintf("  Control arm <%d events: %d (%.1f%%)\n",
              threshold, low_H_0, 100*low_H_0/nb_boots))
  cat(sprintf("  Treatment arm <%d events: %d (%.1f%%)\n",
              threshold, low_H_1, 100*low_H_1/nb_boots))
  cat(sprintf("  Either arm <%d events: %d (%.1f%%)\n\n",
              threshold, low_H_either, 100*low_H_either/nb_boots))

  cat("ORIGINAL Subgroup Hc on BOOTSTRAP samples:\n")
  cat(sprintf("  Control arm <%d events: %d (%.1f%%)\n",
              threshold, low_Hc_0, 100*low_Hc_0/nb_boots))
  cat(sprintf("  Treatment arm <%d events: %d (%.1f%%)\n",
              threshold, low_Hc_1, 100*low_Hc_1/nb_boots))
  cat(sprintf("  Either arm <%d events: %d (%.1f%%)\n\n",
              threshold, low_Hc_either, 100*low_Hc_either/nb_boots))

  # Count successful bootstraps (found new subgroup)
  n_successful <- sum(!is.na(results$events_Hstar_0))

  if (n_successful > 0) {
    cat(sprintf("NEW Subgroups found: %d (%.1f%%)\n\n",
                n_successful, 100*n_successful/nb_boots))

    cat("NEW Subgroup H* on ORIGINAL data:\n")
    cat(sprintf("  Control arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hstar_0, 100*low_Hstar_0/n_successful))
    cat(sprintf("  Treatment arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hstar_1, 100*low_Hstar_1/n_successful))
    cat(sprintf("  Either arm <%d events: %d (%.1f%% of successful)\n\n",
                threshold, low_Hstar_either, 100*low_Hstar_either/n_successful))

    cat("NEW Subgroup Hc* on ORIGINAL data:\n")
    cat(sprintf("  Control arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_0, 100*low_Hcstar_0/n_successful))
    cat(sprintf("  Treatment arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_1, 100*low_Hcstar_1/n_successful))
    cat(sprintf("  Either arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_either, 100*low_Hcstar_either/n_successful))
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

#' Summarize Bootstrap Timing
#'
#' Creates a comprehensive timing summary for bootstrap results
#'
#' @param boot_results Data.table from bootstrap_results()
#' @return List with timing summaries
#' @export
summarize_bootstrap_timing <- function(boot_results) {

  # Overall timing from attributes
  overall <- attr(boot_results, "timing")

  # Per-iteration statistics
  iteration_stats <- list(
    mean = mean(boot_results$tmins_iteration, na.rm = TRUE),
    median = median(boot_results$tmins_iteration, na.rm = TRUE),
    sd = sd(boot_results$tmins_iteration, na.rm = TRUE),
    min = min(boot_results$tmins_iteration, na.rm = TRUE),
    max = max(boot_results$tmins_iteration, na.rm = TRUE),
    q25 = quantile(boot_results$tmins_iteration, 0.25, na.rm = TRUE),
    q75 = quantile(boot_results$tmins_iteration, 0.75, na.rm = TRUE)
  )

  # ForestSearch timing (subset where forestsearch ran)
  fs_times <- boot_results$tmins_search[!is.na(boot_results$tmins_search)]
  if (length(fs_times) > 0) {
    fs_stats <- list(
      n_runs = length(fs_times),
      mean = mean(fs_times),
      median = median(fs_times),
      sd = sd(fs_times),
      min = min(fs_times),
      max = max(fs_times)
    )
  } else {
    fs_stats <- NULL
  }

  # Success rate
  n_success <- sum(!is.na(boot_results$H_biasadj_2))
  success_rate <- n_success / nrow(boot_results)

  # Print summary
  cat("\n=== BOOTSTRAP TIMING SUMMARY ===\n\n")
  cat("Overall:\n")
  cat(sprintf("  Total time: %.2f minutes (%.2f hours)\n",
              overall$total_minutes, overall$total_hours))
  cat(sprintf("  Iterations: %d\n", overall$n_boots))
  cat(sprintf("  Success rate: %.1f%% (%d/%d)\n",
              success_rate * 100, n_success, overall$n_boots))

  cat("\nPer-iteration timing:\n")
  cat(sprintf("  Mean: %.2f minutes (%.1f seconds)\n",
              iteration_stats$mean, iteration_stats$mean * 60))
  cat(sprintf("  Median: %.2f minutes (%.1f seconds)\n",
              iteration_stats$median, iteration_stats$median * 60))
  cat(sprintf("  Range: [%.2f, %.2f] minutes\n",
              iteration_stats$min, iteration_stats$max))
  cat(sprintf("  IQR: [%.2f, %.2f] minutes\n",
              iteration_stats$q25, iteration_stats$q75))

  if (!is.null(fs_stats)) {
    cat("\nForestSearch timing (successful iterations only):\n")
    cat(sprintf("  Runs: %d (%.1f%%)\n",
                fs_stats$n_runs, fs_stats$n_runs / overall$n_boots * 100))
    cat(sprintf("  Mean: %.2f minutes (%.1f seconds)\n",
                fs_stats$mean, fs_stats$mean * 60))
    cat(sprintf("  Median: %.2f minutes (%.1f seconds)\n",
                fs_stats$median, fs_stats$median * 60))
  }

  cat("\n")

  # Return invisible
  invisible(list(
    overall = overall,
    iteration_stats = iteration_stats,
    fs_stats = fs_stats,
    success_rate = success_rate
  ))
}



#' Create Bootstrap Timing Plots
#'
#' Generates timing visualization plots for bootstrap analysis
#'
#' @param results Data.table of bootstrap results with timing columns
#' @return List of ggplot objects
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_point
#' @keywords internal

create_bootstrap_timing_plots <- function(results) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }

  plots <- list()

  # Check which timing columns are available
  has_iteration_timing <- "tmins_iteration" %in% names(results)
  has_search_timing <- "tmins_search" %in% names(results)

  # Plot 1: Iteration timing distribution
  if (has_iteration_timing) {
    p1 <- ggplot2::ggplot(results, ggplot2::aes(x = tmins_iteration)) +
      ggplot2::geom_histogram(bins = 30, fill = "#2E86AB", alpha = 0.7, color = "white") +
      ggplot2::geom_vline(xintercept = median(results$tmins_iteration, na.rm = TRUE),
                          color = "red", linetype = "dashed", linewidth = 1) +
      ggplot2::labs(
        title = "Bootstrap Iteration Timing Distribution",
        subtitle = sprintf("Median: %.2f min, Mean: %.2f min",
                           median(results$tmins_iteration, na.rm = TRUE),
                           mean(results$tmins_iteration, na.rm = TRUE)),
        x = "Minutes per iteration",
        y = "Frequency"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
      )

    plots$iteration_time <- p1
  }

  # Plot 2: ForestSearch timing (for successful iterations)
  if (has_search_timing) {
    fs_data <- results[!is.na(results$tmins_search), ]

    if (nrow(fs_data) > 0) {
      p2 <- ggplot2::ggplot(fs_data, ggplot2::aes(x = tmins_search)) +
        ggplot2::geom_histogram(bins = 30, fill = "#A23B72", alpha = 0.7, color = "white") +
        ggplot2::geom_vline(xintercept = median(fs_data$tmins_search, na.rm = TRUE),
                            color = "red", linetype = "dashed", linewidth = 1) +
        ggplot2::labs(
          title = "ForestSearch Timing Distribution",
          subtitle = sprintf("Successful iterations only (n=%d, %.1f%%)",
                             nrow(fs_data),
                             nrow(fs_data) / nrow(results) * 100),
          x = "Minutes for ForestSearch",
          y = "Frequency"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
        )

      plots$search_time <- p2
    }
  }

  # Plot 3: Timing over iterations (to detect trends)
  if (has_iteration_timing) {
    p3 <- ggplot2::ggplot(results, ggplot2::aes(x = boot_id, y = tmins_iteration)) +
      ggplot2::geom_point(alpha = 0.5, color = "#2E86AB") +
      ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE) +
      ggplot2::labs(
        title = "Bootstrap Timing Trend",
        subtitle = "Timing across iterations (with LOESS smoother)",
        x = "Bootstrap Iteration",
        y = "Minutes per iteration"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
      )

    plots$timing_trend <- p3
  }

  # Plot 4: Timing breakdown (if both iteration and search available)
  if (has_iteration_timing && has_search_timing) {
    # Calculate overhead
    timing_data <- data.frame(
      boot_id = results$boot_id,
      forestsearch = results$tmins_search,
      overhead = results$tmins_iteration - results$tmins_search
    )

    # Remove rows with NA
    timing_data <- timing_data[complete.cases(timing_data), ]

    if (nrow(timing_data) > 0) {
      # Reshape for stacked bar chart
      timing_long <- data.frame(
        boot_id = rep(timing_data$boot_id, 2),
        component = rep(c("ForestSearch", "Overhead"), each = nrow(timing_data)),
        minutes = c(timing_data$forestsearch, timing_data$overhead)
      )

      p4 <- ggplot2::ggplot(timing_long, ggplot2::aes(x = boot_id, y = minutes, fill = component)) +
        ggplot2::geom_col(position = "stack", alpha = 0.7) +
        ggplot2::scale_fill_manual(values = c("ForestSearch" = "#A23B72", "Overhead" = "#2E86AB")) +
        ggplot2::labs(
          title = "Bootstrap Timing Breakdown",
          subtitle = "ForestSearch vs Overhead (Cox models, bias correction)",
          x = "Bootstrap Iteration (successful only)",
          y = "Minutes",
          fill = "Component"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(size = 11, color = "gray40"),
          legend.position = "bottom"
        )

      plots$timing_breakdown <- p4
    }
  }

  # Plot 5: Box plot comparison
  if (has_iteration_timing && has_search_timing) {
    timing_data <- data.frame(
      boot_id = results$boot_id,
      Total = results$tmins_iteration,
      ForestSearch = results$tmins_search,
      Overhead = results$tmins_iteration - results$tmins_search
    )

    # Reshape for box plot
    timing_long <- tidyr::pivot_longer(
      timing_data,
      cols = c(Total, ForestSearch, Overhead),
      names_to = "Component",
      values_to = "Minutes"
    )

    # Remove NA
    timing_long <- timing_long[!is.na(timing_long$Minutes), ]

    if (nrow(timing_long) > 0) {
      p5 <- ggplot2::ggplot(timing_long, ggplot2::aes(x = Component, y = Minutes, fill = Component)) +
        ggplot2::geom_boxplot(alpha = 0.7) +
        ggplot2::scale_fill_manual(values = c(
          "Total" = "#2E86AB",
          "ForestSearch" = "#A23B72",
          "Overhead" = "#F18F01"
        )) +
        ggplot2::labs(
          title = "Bootstrap Timing Components",
          subtitle = "Distribution comparison",
          x = NULL,
          y = "Minutes"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(size = 11, color = "gray40"),
          legend.position = "none"
        )

      plots$timing_boxplot <- p5
    }
  }

  return(plots)
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
#' @importFrom gt gt tab_header tab_spanner tab_footnote tab_source_note md cols_label tab_style cell_fill cell_text
#' @export

format_bootstrap_timing_table <- function(timing_list, nb_boots, boot_success_rate) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table formatting. Install with: install.packages('gt')")
  }

  if (is.null(timing_list)) {
    stop("timing_list cannot be NULL")
  }

  # Extract components
  overall <- timing_list$overall
  iteration_stats <- timing_list$iteration_stats
  fs_stats <- timing_list$fs_stats
  overhead_stats <- timing_list$overhead_stats

  # Build the data frame for the table
  timing_rows <- list()

  # =========================================================================
  # SECTION 1: OVERALL TIMING
  # =========================================================================
  if (!is.null(overall)) {
    timing_rows$overall <- data.frame(
      Category = c("Overall", "", "", ""),
      Metric = c(
        "Total time",
        "Average per boot",
        "Successful boots",
        "Projected for 1000 boots"
      ),
      Value = c(
        sprintf("%.2f min (%.2f hrs)", overall$total_minutes, overall$total_hours),
        sprintf("%.2f min (%.1f sec)", overall$avg_minutes_per_boot, overall$avg_seconds_per_boot),
        sprintf("%d (%.1f%%)",
                round(boot_success_rate * nb_boots),
                boot_success_rate * 100),
        if (nb_boots < 1000) {
          projected_min <- overall$avg_minutes_per_boot * 1000
          sprintf("%.2f min (%.2f hrs)", projected_min, projected_min / 60)
        } else {
          "—"
        }
      ),
      stringsAsFactors = FALSE
    )
  }

  # =========================================================================
  # SECTION 2: PER-ITERATION STATISTICS
  # =========================================================================
  if (!is.null(iteration_stats)) {
    timing_rows$iteration <- data.frame(
      Category = c("Per-Iteration", "", "", "", "", "", ""),
      Metric = c(
        "Mean",
        "Median",
        "Std Dev",
        "Minimum",
        "Maximum",
        "25th percentile",
        "75th percentile"
      ),
      Value = c(
        sprintf("%.2f min (%.1f sec)",
                iteration_stats$mean, iteration_stats$mean * 60),
        sprintf("%.2f min (%.1f sec)",
                iteration_stats$median, iteration_stats$median * 60),
        sprintf("%.2f min", iteration_stats$sd),
        sprintf("%.2f min", iteration_stats$min),
        sprintf("%.2f min", iteration_stats$max),
        sprintf("%.2f min", iteration_stats$q25),
        sprintf("%.2f min", iteration_stats$q75)
      ),
      stringsAsFactors = FALSE
    )
  }

  # =========================================================================
  # SECTION 3: FORESTSEARCH TIMING
  # =========================================================================
  if (!is.null(fs_stats)) {
    timing_rows$forestsearch <- data.frame(
      Category = c("ForestSearch", "", "", "", ""),
      Metric = c(
        "Iterations with FS",
        "Mean FS time",
        "Median FS time",
        "Total FS time",
        "% of total time"
      ),
      Value = c(
        sprintf("%d (%.1f%%)", fs_stats$n_runs, fs_stats$pct_runs),
        sprintf("%.2f min (%.1f sec)", fs_stats$mean, fs_stats$mean * 60),
        sprintf("%.2f min (%.1f sec)", fs_stats$median, fs_stats$median * 60),
        sprintf("%.2f min", fs_stats$total),
        sprintf("%.1f%%", fs_stats$total / overall$total_minutes * 100)
      ),
      stringsAsFactors = FALSE
    )
  }

  # =========================================================================
  # SECTION 4: OVERHEAD TIMING
  # =========================================================================
  if (!is.null(overhead_stats)) {
    timing_rows$overhead <- data.frame(
      Category = c("Overhead", "", "", ""),
      Metric = c(
        "Mean overhead",
        "Median overhead",
        "Total overhead",
        "% of total time"
      ),
      Value = c(
        sprintf("%.2f min (%.1f sec)",
                overhead_stats$mean, overhead_stats$mean * 60),
        sprintf("%.2f min (%.1f sec)",
                overhead_stats$median, overhead_stats$median * 60),
        sprintf("%.2f min", overhead_stats$total),
        sprintf("%.1f%%", overhead_stats$pct_of_total)
      ),
      stringsAsFactors = FALSE
    )
  }

  # =========================================================================
  # SECTION 5: PERFORMANCE ASSESSMENT
  # =========================================================================
  if (!is.null(iteration_stats)) {
    avg_sec <- overall$avg_seconds_per_boot

    if (avg_sec < 5) {
      performance <- "Excellent ✓✓✓"
      perf_color <- "#d4edda"
    } else if (avg_sec < 15) {
      performance <- "Good ✓✓"
      perf_color <- "#d1ecf1"
    } else if (avg_sec < 30) {
      performance <- "Acceptable ✓"
      perf_color <- "#fff3cd"
    } else if (avg_sec < 60) {
      performance <- "Slow ⚠"
      perf_color <- "#f8d7da"
    } else {
      performance <- "Very Slow ⚠⚠"
      perf_color <- "#f8d7da"
    }

    timing_rows$performance <- data.frame(
      Category = c("Performance", ""),
      Metric = c(
        "Rating",
        "Speed"
      ),
      Value = c(
        performance,
        sprintf("%.1f sec/iteration", avg_sec)
      ),
      stringsAsFactors = FALSE
    )
  }

  # =========================================================================
  # SECTION 6: COMBINE ALL ROWS AND CREATE GT TABLE
  # =========================================================================

  timing_df <- do.call(rbind, timing_rows)
  rownames(timing_df) <- NULL

  # Create the gt table
  tbl <- timing_df |>
    gt::gt() |>

    # Title and subtitle
    gt::tab_header(
      title = gt::md("**Bootstrap Timing Analysis**"),
      subtitle = sprintf("%d iterations (%.1f%% successful)",
                         nb_boots, boot_success_rate * 100)
    ) |>

    # Column labels
    gt::cols_label(
      Category = "Category",
      Metric = "Metric",
      Value = "Value"
    ) |>

    # Style the category column (bold and slight background)
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f8f9fa"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = Category,
        rows = Category %in% c("Overall", "Per-Iteration", "ForestSearch",
                               "Overhead", "Performance")
      )
    ) |>

    # Highlight performance row if it exists
    gt::tab_style(
      style = list(
        gt::cell_fill(color = if (!is.null(iteration_stats)) perf_color else "#ffffff"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c(Metric, Value),
        rows = Category == "Performance"
      )
    )

  # =========================================================================
  # SECTION 7: ADD EXPLANATORY FOOTNOTES
  # =========================================================================

  tbl <- tbl |>
    gt::tab_footnote(
      footnote = gt::md("**Overall**: Total time includes all bootstrap iterations plus Ystar matrix generation"),
      locations = gt::cells_body(columns = Category, rows = Category == "Overall")
    )

  if (!is.null(fs_stats)) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md("**ForestSearch**: Time spent running subgroup search (successful iterations only)"),
        locations = gt::cells_body(columns = Category, rows = Category == "ForestSearch")
      )
  }

  if (!is.null(overhead_stats)) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md("**Overhead**: Time for Cox models, bias correction, and data management"),
        locations = gt::cells_body(columns = Category, rows = Category == "Overhead")
      )
  }

  # =========================================================================
  # SECTION 8: ADD RECOMMENDATIONS IF PERFORMANCE IS SLOW
  # =========================================================================

  if (!is.null(iteration_stats) && overall$avg_seconds_per_boot > 30) {
    recommendations <- "**Recommendations for improvement:**\n"

    if (!is.null(fs_stats) && fs_stats$mean > 0.5) {
      recommendations <- paste0(recommendations,
                                "• Consider reducing `max.minutes` in forestsearch\n",
                                "• Consider reducing `maxk` if currently > 2\n"
      )
    }

    if (nb_boots > 500) {
      recommendations <- paste0(recommendations,
                                "• Consider reducing `nb_boots` for initial testing\n"
      )
    }

    recommendations <- paste0(recommendations,
                              "• Ensure sufficient parallel workers are allocated"
    )

    tbl <- tbl |>
      gt::tab_source_note(
        source_note = gt::md(recommendations)
      )
  }

  # =========================================================================
  # SECTION 9: FINAL STYLING
  # =========================================================================

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


#' Format Bootstrap Diagnostics Table with gt
#'
#' Creates a publication-ready diagnostics table from bootstrap results,
#' showing success rates and quality metrics.
#'
#' @param diagnostics List. Diagnostics information from summarize_bootstrap_results()
#' @param nb_boots Integer. Number of bootstrap iterations
#' @param results Data.table. Bootstrap results with bias-corrected estimates
#' @param H_estimates List. H subgroup estimates
#' @param Hc_estimates List. Hc subgroup estimates
#'
#' @return A gt table object
#' @importFrom gt gt tab_header tab_spanner tab_footnote tab_source_note md cols_label tab_style cell_fill cell_text cells_body cells_column_labels
#' @export

format_bootstrap_diagnostics_table <- function(diagnostics, nb_boots, results,
                                               H_estimates = NULL, Hc_estimates = NULL) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table formatting. Install with: install.packages('gt')")
  }

  if (is.null(diagnostics)) {
    stop("diagnostics cannot be NULL")
  }

  # =========================================================================
  # SECTION 1: BOOTSTRAP SUCCESS METRICS
  # =========================================================================

  success_rate <- diagnostics$success_rate
  n_successful <- diagnostics$n_successful
  n_failed <- diagnostics$n_failed

  # Determine success rate color coding
  if (success_rate >= 0.90) {
    success_color <- "#d4edda"  # Green
    success_rating <- "Excellent ✓✓✓"
  } else if (success_rate >= 0.75) {
    success_color <- "#d1ecf1"  # Blue
    success_rating <- "Good ✓✓"
  } else if (success_rate >= 0.50) {
    success_color <- "#fff3cd"  # Yellow
    success_rating <- "Acceptable ✓"
  } else if (success_rate >= 0.25) {
    success_color <- "#f8d7da"  # Light red
    success_rating <- "Poor ⚠"
  } else {
    success_color <- "#f8d7da"  # Red
    success_rating <- "Very Poor ⚠⚠"
  }

  diagnostics_rows <- list()

  # Success metrics
  diagnostics_rows$success <- data.frame(
    Category = c("Success Rate", "", "", ""),
    Metric = c(
      "Total iterations",
      "Successful subgroup ID",
      "Failed to find subgroup",
      "Success rating"
    ),
    Value = c(
      sprintf("%d", nb_boots),
      sprintf("%d (%.1f%%)", n_successful, success_rate * 100),
      sprintf("%d (%.1f%%)", n_failed, (1 - success_rate) * 100),
      success_rating
    ),
    stringsAsFactors = FALSE
  )

  # =========================================================================
  # SECTION 2: BIAS-CORRECTED ESTIMATES QUALITY
  # =========================================================================

  if (!is.null(H_estimates) && !is.null(Hc_estimates)) {

    # Calculate bias correction impact
    H_bias_impact <- abs(H_estimates$H2 - H_estimates$H0) / H_estimates$H0 * 100
    Hc_bias_impact <- abs(Hc_estimates$H2 - Hc_estimates$H0) / Hc_estimates$H0 * 100

    # Calculate CI widths
    H_ci_width_raw <- H_estimates$H0_upper - H_estimates$H0_lower
    H_ci_width_bc <- H_estimates$H2_upper - H_estimates$H2_lower
    Hc_ci_width_raw <- Hc_estimates$H0_upper - Hc_estimates$H0_lower
    Hc_ci_width_bc <- Hc_estimates$H2_upper - Hc_estimates$H2_lower

    diagnostics_rows$estimates_H <- data.frame(
      Category = c("Subgroup H (Questionable)", "", "", ""),
      Metric = c(
        "Unadjusted estimate",
        "Bias-corrected estimate",
        "Bias correction impact",
        "CI width change"
      ),
      Value = c(
        sprintf("%.2f (%.2f, %.2f)",
                H_estimates$H0, H_estimates$H0_lower, H_estimates$H0_upper),
        sprintf("%.2f (%.2f, %.2f)",
                H_estimates$H2, H_estimates$H2_lower, H_estimates$H2_upper),
        sprintf("%.1f%%", H_bias_impact),
        sprintf("%.2f → %.2f", H_ci_width_raw, H_ci_width_bc)
      ),
      stringsAsFactors = FALSE
    )

    diagnostics_rows$estimates_Hc <- data.frame(
      Category = c("Subgroup Hc (Recommend)", "", "", ""),
      Metric = c(
        "Unadjusted estimate",
        "Bias-corrected estimate",
        "Bias correction impact",
        "CI width change"
      ),
      Value = c(
        sprintf("%.2f (%.2f, %.2f)",
                Hc_estimates$H0, Hc_estimates$H0_lower, Hc_estimates$H0_upper),
        sprintf("%.2f (%.2f, %.2f)",
                Hc_estimates$H2, Hc_estimates$H2_lower, Hc_estimates$H2_upper),
        sprintf("%.1f%%", Hc_bias_impact),
        sprintf("%.2f → %.2f", Hc_ci_width_raw, Hc_ci_width_bc)
      ),
      stringsAsFactors = FALSE
    )
  }

  # =========================================================================
  # SECTION 3: BOOTSTRAP DISTRIBUTION QUALITY METRICS
  # =========================================================================

  # Calculate quality metrics from bootstrap distribution
  if (!is.null(results)) {

    # For H subgroup
    H_valid <- results$H_biasadj_2[!is.na(results$H_biasadj_2)]
    if (length(H_valid) > 0) {
      H_mean <- mean(H_valid)
      H_sd <- sd(H_valid)
      H_cv <- (H_sd / abs(H_mean)) * 100
      H_skew <- calculate_skewness(H_valid)

      diagnostics_rows$quality_H <- data.frame(
        Category = c("Bootstrap Quality: H", "", "", ""),
        Metric = c(
          "Valid iterations",
          "Mean (SD)",
          "Coefficient of variation",
          "Skewness"
        ),
        Value = c(
          sprintf("%d", length(H_valid)),
          sprintf("%.2f (%.2f)", H_mean, H_sd),
          sprintf("%.1f%%", H_cv),
          sprintf("%.2f", H_skew)
        ),
        stringsAsFactors = FALSE
      )
    }

    # For Hc subgroup
    Hc_valid <- results$Hc_biasadj_2[!is.na(results$Hc_biasadj_2)]
    if (length(Hc_valid) > 0) {
      Hc_mean <- mean(Hc_valid)
      Hc_sd <- sd(Hc_valid)
      Hc_cv <- (Hc_sd / abs(Hc_mean)) * 100
      Hc_skew <- calculate_skewness(Hc_valid)

      diagnostics_rows$quality_Hc <- data.frame(
        Category = c("Bootstrap Quality: Hc", "", "", ""),
        Metric = c(
          "Valid iterations",
          "Mean (SD)",
          "Coefficient of variation",
          "Skewness"
        ),
        Value = c(
          sprintf("%d", length(Hc_valid)),
          sprintf("%.2f (%.2f)", Hc_mean, Hc_sd),
          sprintf("%.1f%%", Hc_cv),
          sprintf("%.2f", Hc_skew)
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  # =========================================================================
  # SECTION 4: SEARCH PERFORMANCE METRICS (if available)
  # =========================================================================

  if (!is.null(results) && "max_sg_est" %in% names(results)) {

    max_sg_valid <- results$max_sg_est[!is.na(results$max_sg_est)]
    L_valid <- results$L[!is.na(results$L)]
    max_count_valid <- results$max_count[!is.na(results$max_count)]

    if (length(max_sg_valid) > 0) {
      diagnostics_rows$search <- data.frame(
        Category = c("Search Performance", "", "", ""),
        Metric = c(
          "Mean max HR found",
          "Mean factors evaluated",
          "Mean combinations tried",
          "Proportion at maxk"
        ),
        Value = c(
          sprintf("%.2f (%.2f)", mean(max_sg_valid), sd(max_sg_valid)),
          sprintf("%.1f", mean(L_valid)),
          sprintf("%.0f", mean(max_count_valid)),
          if ("prop_maxk" %in% names(results)) {
            sprintf("%.1f%%", mean(results$prop_maxk, na.rm = TRUE) * 100)
          } else {
            "—"
          }
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  # =========================================================================
  # SECTION 5: COMBINE AND CREATE GT TABLE
  # =========================================================================

  diagnostics_df <- do.call(rbind, diagnostics_rows)
  rownames(diagnostics_df) <- NULL

  # Create the gt table
  tbl <- diagnostics_df |>
    gt::gt() |>

    # Title and subtitle
    gt::tab_header(
      title = gt::md("**Bootstrap Diagnostics Summary**"),
      subtitle = sprintf("Analysis of %d bootstrap iterations", nb_boots)
    ) |>

    # Column labels
    gt::cols_label(
      Category = "Category",
      Metric = "Metric",
      Value = "Value"
    ) |>

    # Style category rows (bold and background)
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f8f9fa"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = Category,
        rows = Category %in% c(
          "Success Rate",
          "Subgroup H (Questionable)",
          "Subgroup Hc (Recommend)",
          "Bootstrap Quality: H",
          "Bootstrap Quality: Hc",
          "Search Performance"
        )
      )
    ) |>

    # Highlight success rating row
    gt::tab_style(
      style = list(
        gt::cell_fill(color = success_color),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c(Metric, Value),
        rows = Metric == "Success rating"
      )
    )

  # =========================================================================
  # SECTION 6: ADD EXPLANATORY FOOTNOTES
  # =========================================================================

  tbl <- tbl |>
    gt::tab_footnote(
      footnote = gt::md("**Success Rate**: Proportion of bootstrap samples where ForestSearch identified a valid subgroup"),
      locations = gt::cells_body(columns = Category, rows = Category == "Success Rate")
    )

  if (!is.null(H_estimates)) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md("**Bias Correction Impact**: Percentage change from unadjusted to bias-corrected estimate"),
        locations = gt::cells_body(columns = Metric, rows = Metric == "Bias correction impact")
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**CI Width Change**: Confidence interval width before → after bias correction"),
        locations = gt::cells_body(columns = Metric, rows = Metric == "CI width change")
      )
  }

  if (!is.null(results) && length(H_valid) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md("**Coefficient of Variation**: Standard deviation as % of mean (lower is better)"),
        locations = gt::cells_body(columns = Metric, rows = Metric == "Coefficient of variation")
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**Skewness**: Measure of asymmetry (0 = symmetric, |skew| < 1 is generally good)"),
        locations = gt::cells_body(columns = Metric, rows = Metric == "Skewness")
      )
  }

  # =========================================================================
  # SECTION 7: ADD INTERPRETATION GUIDANCE
  # =========================================================================

  interpretation <- "**Interpretation Guide:**\n\n"

  if (success_rate >= 0.90) {
    interpretation <- paste0(interpretation,
                             "✓ **Excellent stability**: Subgroup is consistently identified across bootstrap samples.\n"
    )
  } else if (success_rate >= 0.75) {
    interpretation <- paste0(interpretation,
                             "✓ **Good stability**: Subgroup is reliably identified in most bootstrap samples.\n"
    )
  } else if (success_rate >= 0.50) {
    interpretation <- paste0(interpretation,
                             "⚠ **Moderate stability**: Subgroup identification is somewhat unstable. Consider:\n",
                             "  • Increasing sample size\n",
                             "  • Adjusting consistency threshold\n",
                             "  • Examining failed iterations for patterns\n"
    )
  } else {
    interpretation <- paste0(interpretation,
                             "⚠ **Poor stability**: Subgroup is rarely identified. Consider:\n",
                             "  • Reviewing subgroup criteria (n.min, hr.threshold)\n",
                             "  • Increasing sample size significantly\n",
                             "  • Simplifying search (reduce maxk)\n",
                             "  • Examining if subgroup is real or spurious\n"
    )
  }

  # Add CV interpretation
  if (!is.null(results) && length(H_valid) > 0) {
    if (H_cv < 10) {
      interpretation <- paste0(interpretation,
                               "\n✓ **Low variability**: Bootstrap estimates are precise (CV < 10%).\n"
      )
    } else if (H_cv < 25) {
      interpretation <- paste0(interpretation,
                               "\n✓ **Moderate variability**: Bootstrap estimates show acceptable precision.\n"
      )
    } else {
      interpretation <- paste0(interpretation,
                               "\n⚠ **High variability**: Bootstrap estimates are imprecise (CV ≥ 25%). ",
                               "Consider increasing nb_boots or sample size.\n"
      )
    }
  }

  tbl <- tbl |>
    gt::tab_source_note(
      source_note = gt::md(interpretation)
    )

  # =========================================================================
  # SECTION 8: FINAL STYLING
  # =========================================================================

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

#' Calculate Skewness
#'
#' Helper function to calculate sample skewness
#'
#' @param x Numeric vector
#' @return Numeric skewness value
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



#' Create Subgroup Characteristics Summary Table
#'
#' Analyzes the characteristics of identified subgroups across bootstrap iterations.
#' Includes agreement rates, consistency metrics, and factor frequencies.
#'
#' @param results Data.table. Bootstrap results with subgroup characteristics
#' @param nb_boots Integer. Total number of bootstrap iterations
#' @param original_sg Character vector. Original subgroup definition (M.1, M.2 from fs.est)
#' @param maxk Integer. Maximum number of factors allowed
#'
#' @return List with multiple summary components
#' @importFrom data.table data.table .N
#' @export

summarize_bootstrap_subgroups <- function(results, nb_boots,
                                          original_sg = NULL,
                                          maxk = 2) {

  # =========================================================================
  # SECTION 1: BASIC STATISTICS
  # =========================================================================

  # Filter to successful iterations (where subgroup was found)
  sg_found <- results[!is.na(Pcons)]
  n_found <- nrow(sg_found)
  pct_found <- 100 * n_found / nb_boots

  if (n_found == 0) {
    warning("No subgroups identified in any bootstrap iteration")
    return(list(
      basic_stats = NULL,
      consistency_dist = NULL,
      size_dist = NULL,
      factor_freq = NULL,
      agreement = NULL
    ))
  }

  # =========================================================================
  # SECTION 2: BASIC STATISTICS TABLE
  # =========================================================================

  basic_stats <- data.table::data.table(
    Metric = c(
      "Total bootstrap iterations",
      "Subgroups identified",
      "Success rate (%)",
      "",
      "Consistency (Pcons)",
      "  Mean",
      "  Median",
      "  SD",
      "  Min",
      "  Max",
      "  Q25",
      "  Q75",
      "",
      "Hazard Ratio (hr_sg)",
      "  Mean",
      "  Median",
      "  SD",
      "  Min",
      "  Max",
      "  Q25",
      "  Q75",
      "",
      "Subgroup Size (N_sg)",
      "  Mean",
      "  Median",
      "  SD",
      "  Min",
      "  Max",
      "  Q25",
      "  Q75",
      "",
      "Number of Factors (K_sg)",
      "  Mean",
      "  Median",
      "  Mode"
    ),
    Value = c(
      as.character(nb_boots),
      as.character(n_found),
      sprintf("%.1f%%", pct_found),
      "",
      "",
      sprintf("%.3f", mean(sg_found$Pcons, na.rm = TRUE)),
      sprintf("%.3f", median(sg_found$Pcons, na.rm = TRUE)),
      sprintf("%.3f", sd(sg_found$Pcons, na.rm = TRUE)),
      sprintf("%.3f", min(sg_found$Pcons, na.rm = TRUE)),
      sprintf("%.3f", max(sg_found$Pcons, na.rm = TRUE)),
      sprintf("%.3f", quantile(sg_found$Pcons, 0.25, na.rm = TRUE)),
      sprintf("%.3f", quantile(sg_found$Pcons, 0.75, na.rm = TRUE)),
      "",
      "",
      sprintf("%.2f", mean(sg_found$hr_sg, na.rm = TRUE)),
      sprintf("%.2f", median(sg_found$hr_sg, na.rm = TRUE)),
      sprintf("%.2f", sd(sg_found$hr_sg, na.rm = TRUE)),
      sprintf("%.2f", min(sg_found$hr_sg, na.rm = TRUE)),
      sprintf("%.2f", max(sg_found$hr_sg, na.rm = TRUE)),
      sprintf("%.2f", quantile(sg_found$hr_sg, 0.25, na.rm = TRUE)),
      sprintf("%.2f", quantile(sg_found$hr_sg, 0.75, na.rm = TRUE)),
      "",
      "",
      sprintf("%.0f", mean(sg_found$N_sg, na.rm = TRUE)),
      sprintf("%.0f", median(sg_found$N_sg, na.rm = TRUE)),
      sprintf("%.0f", sd(sg_found$N_sg, na.rm = TRUE)),
      sprintf("%.0f", min(sg_found$N_sg, na.rm = TRUE)),
      sprintf("%.0f", max(sg_found$N_sg, na.rm = TRUE)),
      sprintf("%.0f", quantile(sg_found$N_sg, 0.25, na.rm = TRUE)),
      sprintf("%.0f", quantile(sg_found$N_sg, 0.75, na.rm = TRUE)),
      "",
      "",
      sprintf("%.1f", mean(sg_found$K_sg, na.rm = TRUE)),
      sprintf("%.0f", median(sg_found$K_sg, na.rm = TRUE)),
      sprintf("%.0f", as.numeric(names(sort(table(sg_found$K_sg), decreasing = TRUE)[1])))
    )
  )

  # =========================================================================
  # SECTION 3: FACTOR FREQUENCY TABLE
  # =========================================================================

  # Count frequency of each factor appearing
  factor_cols <- paste0("M.", 1:maxk)
  factor_freq_list <- list()

  for (i in 1:maxk) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      freq <- sg_found[!is.na(get(col)) & get(col) != "", .N, by = col]
      setnames(freq, col, "Factor")
      freq[, Position := paste0("M.", i)]
      freq[, Percent := 100 * N / n_found]
      freq <- freq[order(-N)]
      factor_freq_list[[i]] <- freq
    }
  }

  factor_freq <- rbindlist(factor_freq_list, fill = TRUE)
  setcolorder(factor_freq, c("Position", "Factor", "N", "Percent"))

  # =========================================================================
  # SECTION 4: SUBGROUP DEFINITION AGREEMENT
  # =========================================================================

  # For maxk=1: just M.1
  # For maxk=2: M.1 & M.2 combination
  # For maxk=3: M.1 & M.2 & M.3 combination

  if (maxk == 1) {
    agreement <- sg_found[!is.na(M.1) & M.1 != "", .N, by = .(M.1)]
    setnames(agreement, "M.1", "Subgroup")
  } else if (maxk == 2) {
    # Count unique combinations
    agreement <- sg_found[!is.na(M.1) & M.1 != "", .N,
                          by = .(M.1, M.2, K_sg)]
    # Create readable label
    agreement[, Subgroup := ifelse(
      is.na(M.2) | M.2 == "",
      M.1,
      paste(M.1, "&", M.2)
    )]
    agreement <- agreement[, .(Subgroup, K_sg, N)]
  } else if (maxk == 3) {
    agreement <- sg_found[!is.na(M.1) & M.1 != "", .N,
                          by = .(M.1, M.2, M.3, K_sg)]
    agreement[, Subgroup := paste(
      ifelse(!is.na(M.1) & M.1 != "", M.1, ""),
      ifelse(!is.na(M.2) & M.2 != "", paste("&", M.2), ""),
      ifelse(!is.na(M.3) & M.3 != "", paste("&", M.3), "")
    )]
    agreement <- agreement[, .(Subgroup, K_sg, N)]
  }

  # # Add percentage and sort
  # agreement[, Percent := 100 * N / n_found]
  # agreement[, Percent_of_successful := 100 * N / n_found]
  # agreement <- agreement[order(-N)]
  # # Add rank
  # agreement[, Rank := .I]

  # Add percentage and sort
  agreement[, Percent_of_successful := 100 * N / n_found]
  agreement <- agreement[order(-N)]
  # Add rank
  agreement[, Rank := .I]


  # FILTER: Keep only subgroups that appear >= 20% of the time
  agreement_common <- agreement[Percent_of_successful >= 20]

  # ALSO: Create individual factor summary (ignoring position)
  #factor_presence <- summarize_factor_presence(sg_found, maxk = maxk)

  factor_presence_results <- summarize_factor_presence(sg_found, maxk = maxk, threshold = 20)


  # =========================================================================
  # SECTION 5: AGREEMENT WITH ORIGINAL SUBGROUP
  # =========================================================================

  original_agreement <- NULL
  if (!is.null(original_sg) && length(original_sg) > 0) {

    # Create original subgroup string for comparison
    if (maxk == 1) {
      original_string <- original_sg[1]
      match_col <- "M.1"
    } else if (maxk == 2 && length(original_sg) >= 2) {
      original_string <- paste(original_sg[1], "&", original_sg[2])
      sg_found[, match_string := ifelse(
        is.na(M.2) | M.2 == "",
        M.1,
        paste(M.1, "&", M.2)
      )]
      match_col <- "match_string"
    } else if (maxk == 3 && length(original_sg) >= 3) {
      original_string <- paste(original_sg[1], "&", original_sg[2], "&", original_sg[3])
      sg_found[, match_string := paste(
        ifelse(!is.na(M.1) & M.1 != "", M.1, ""),
        ifelse(!is.na(M.2) & M.2 != "", paste("&", M.2), ""),
        ifelse(!is.na(M.3) & M.3 != "", paste("&", M.3), "")
      )]
      match_col <- "match_string"
    } else {
      original_string <- original_sg[1]
      match_col <- "M.1"
    }

    # Calculate exact match rate
    n_exact_match <- sum(sg_found[[match_col]] == original_string, na.rm = TRUE)
    pct_exact_match <- 100 * n_exact_match / n_found

    # Calculate any-factor match rate (at least one factor matches)
    any_match <- 0
    for (i in 1:min(length(original_sg), maxk)) {
      col <- paste0("M.", i)
      if (col %in% names(sg_found)) {
        any_match <- any_match + sum(sg_found[[col]] == original_sg[i], na.rm = TRUE)
      }
    }
    pct_any_match <- 100 * (any_match > 0) / n_found

    original_agreement <- data.table::data.table(
      Metric = c(
        "Original subgroup",
        "Exact match count",
        "Exact match rate",
        "Any factor match rate"
      ),
      Value = c(
        original_string,
        as.character(n_exact_match),
        sprintf("%.1f%%", pct_exact_match),
        sprintf("%.1f%%", pct_any_match)
      )
    )
  }

  # =========================================================================
  # SECTION 6: CONSISTENCY DISTRIBUTION
  # =========================================================================

  # Categorize Pcons into bins
  sg_found[, Pcons_bin := cut(
    Pcons,
    breaks = c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0),
    labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.95", "≥0.95"),
    include.lowest = TRUE
  )]

  consistency_dist <- sg_found[, .N, by = Pcons_bin]
  consistency_dist[, Percent := 100 * N / n_found]
  consistency_dist <- consistency_dist[order(Pcons_bin)]
  setnames(consistency_dist, c("Consistency Range", "Count", "Percent"))

  # =========================================================================
  # SECTION 7: SUBGROUP SIZE DISTRIBUTION
  # =========================================================================

  # Categorize N_sg into meaningful bins
  size_breaks <- c(0, 50, 100, 150, 200, 300, Inf)
  size_labels <- c("<50", "50-99", "100-149", "150-199", "200-299", "≥300")

  sg_found[, N_sg_bin := cut(
    N_sg,
    breaks = size_breaks,
    labels = size_labels,
    include.lowest = TRUE
  )]

  size_dist <- sg_found[, .N, by = N_sg_bin]
  size_dist[, Percent := 100 * N / n_found]
  size_dist <- size_dist[order(N_sg_bin)]
  setnames(size_dist, c("Size Range", "Count", "Percent"))

  # =========================================================================
  # RETURN COMPILED RESULTS
  # =========================================================================

  list(
    basic_stats = basic_stats,
    consistency_dist = consistency_dist,
    size_dist = size_dist,
    factor_freq = factor_freq,
    agreement = agreement,
    factor_presence = factor_presence_results$base_factors,  # Base factors
    factor_presence_specific = factor_presence_results$specific_factors,  # Specific definitions >= 20%
    original_agreement = original_agreement,
    n_found = n_found,
    pct_found = pct_found
  )
}

#' Format Subgroup Summary Tables with gt
#'
#' Creates publication-ready gt tables for bootstrap subgroup analysis
#'
#' @param subgroup_summary List from summarize_bootstrap_subgroups()
#' @param nb_boots Integer. Number of bootstrap iterations
#'
#' @return List of gt table objects
#' @importFrom gt gt tab_header tab_spanner tab_footnote md
#' @export

format_subgroup_summary_tables <- function(subgroup_summary, nb_boots) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' required")
  }

  if (is.null(subgroup_summary)) {
    return(NULL)
  }

  tables <- list()

  # =========================================================================
  # TABLE 1: BASIC STATISTICS
  # =========================================================================

  if (!is.null(subgroup_summary$basic_stats)) {
    tables$basic_stats <- subgroup_summary$basic_stats |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Bootstrap Subgroup Characteristics**"),
        subtitle = sprintf("Summary statistics across %d bootstrap iterations", nb_boots)
      ) |>
      gt::cols_label(
        Metric = "",
        Value = "Value"
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#f0f0f0"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(
          rows = Metric %in% c(
            "Consistency (Pcons)",
            "Hazard Ratio (hr_sg)",
            "Subgroup Size (N_sg)",
            "Number of Factors (K_sg)"
          )
        )
      ) |>
      gt::tab_options(
        table.font.size = gt::px(13)
      )
  }

  # =========================================================================
  # TABLE 2: SUBGROUP DEFINITION AGREEMENT
  # =========================================================================

  if (!is.null(subgroup_summary$agreement)) {
    tables$agreement <- subgroup_summary$agreement |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Subgroup Definition Agreement**"),
        subtitle = sprintf("Top subgroup definitions across %d successful iterations",
                           subgroup_summary$n_found)
      ) |>
      gt::cols_label(
        Rank = "Rank",
        Subgroup = "Subgroup Definition",
        K_sg = "K",
        N = "Count",
        Percent_of_successful = "% of Successful"  # CHANGED: was Percent
      ) |>
      gt::fmt_number(
        columns = Percent_of_successful,  # CHANGED: was Percent
        decimals = 1
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(rows = Rank == 1)
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**K** = number of factors in subgroup definition"),
        locations = gt::cells_column_labels(columns = K_sg)
      )
  }



  # =========================================================================
  # TABLE 2C: INDIVIDUAL FACTOR PRESENCE (BASE)
  # =========================================================================

  if (!is.null(subgroup_summary$factor_presence)) {
    tables$factor_presence <- subgroup_summary$factor_presence |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Individual Factor Presence (Base)**"),
        subtitle = "How often each base factor appears (any threshold, any position)"
      ) |>
      gt::cols_label(
        Rank = "Rank",
        Factor = "Factor",
        Count = "Count",
        Percent = "% of Successful"
      ) |>
      gt::fmt_number(
        columns = Percent,
        decimals = 1
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#fff3cd")
        ),
        locations = gt::cells_body(
          rows = Percent >= 20
        )
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**Highlighted rows** appear in ≥20% of successful iterations"),
        locations = gt::cells_column_labels(columns = Percent)
      )
  }

  # =========================================================================
  # TABLE 2D: COMMON SPECIFIC FACTOR DEFINITIONS (>= 20%)
  # =========================================================================

  if (!is.null(subgroup_summary$factor_presence_specific) &&
      nrow(subgroup_summary$factor_presence_specific) > 0) {
    tables$factor_presence_specific <- subgroup_summary$factor_presence_specific |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Common Specific Factor Definitions**"),
        subtitle = "Specific factor definitions appearing in ≥20% of successful iterations"
      ) |>
      gt::cols_label(
        Rank = "Rank",
        Base_Factor = "Base Factor",
        Factor_Definition = "Specific Definition",
        Count = "Count",
        Percent = "% of Successful"
      ) |>
      gt::fmt_number(
        columns = Percent,
        decimals = 1
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8")
        ),
        locations = gt::cells_body(
          columns = Base_Factor
        )
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**Base Factor**: Variable name extracted from full definition"),
        locations = gt::cells_column_labels(columns = Base_Factor)
      )
  }


  # =========================================================================
  # TABLE 3: FACTOR FREQUENCY
  # =========================================================================

  if (!is.null(subgroup_summary$factor_freq)) {
    tables$factor_freq <- subgroup_summary$factor_freq |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Individual Factor Frequencies**"),
        subtitle = "How often each factor appears in identified subgroups"
      ) |>
      gt::cols_label(
        Position = "Position",
        Factor = "Factor",
        N = "Count",
        Percent = "% of Successful"
      ) |>
      gt::fmt_number(
        columns = Percent,
        decimals = 1
      ) |>
      gt::tab_style(
        style = gt::cell_fill(color = "#f0f0f0"),
        locations = gt::cells_body(
          rows = Position == "M.1"
        )
      )
  }

  # =========================================================================
  # TABLE 4: CONSISTENCY DISTRIBUTION
  # =========================================================================

  if (!is.null(subgroup_summary$consistency_dist)) {
    tables$consistency_dist <- subgroup_summary$consistency_dist |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Consistency (P<sub>cons</sub>) Distribution**"),
        subtitle = "Distribution of consistency scores across successful iterations"
      ) |>
      gt::fmt_number(
        columns = Percent,
        decimals = 1
      ) |>
      gt::cols_label(
        `Consistency Range` = gt::md("P<sub>cons</sub> Range"),
        Count = "Count",
        Percent = "%"
      )
  }

  # =========================================================================
  # TABLE 5: ORIGINAL AGREEMENT (if available)
  # =========================================================================

  if (!is.null(subgroup_summary$original_agreement)) {
    tables$original_agreement <- subgroup_summary$original_agreement |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Agreement with Original Subgroup**"),
        subtitle = "How often bootstrap identifies the same subgroup as original analysis"
      ) |>
      gt::cols_label(
        Metric = "",
        Value = "Value"
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(rows = 2:4)
      )
  }

  return(tables)
}


#' Summarize Factor Presence Across Bootstrap Subgroups
#'
#' Analyzes how often each individual factor appears in identified subgroups,
#' extracting base factor names from full definitions and identifying common
#' specific definitions.
#'
#' @param results Data.table. Bootstrap results with subgroup characteristics
#' @param maxk Integer. Maximum number of factors allowed
#' @param threshold Numeric. Percentage threshold for including specific definitions (default: 20)
#' @return List with base_factors and specific_factors data.tables
#' @keywords internal

summarize_factor_presence <- function(results, maxk = 2, threshold = 20) {

  # Filter to successful iterations
  sg_found <- results[!is.na(Pcons)]
  n_found <- nrow(sg_found)

  if (n_found == 0) return(NULL)

  # Extract all unique factors across all positions
  all_factors <- character()
  for (i in 1:maxk) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      factors <- sg_found[[col]][!is.na(sg_found[[col]]) & sg_found[[col]] != ""]
      all_factors <- c(all_factors, factors)
    }
  }

  # =========================================================================
  # PART 1: BASE FACTOR NAMES
  # =========================================================================

  # Function to extract base factor name
  extract_base_factor <- function(factor_string) {
    # Remove leading/trailing whitespace
    factor_string <- trimws(factor_string)

    # Remove negation operator !{...}
    factor_string <- gsub("^!\\{", "{", factor_string)

    # Remove curly braces
    factor_string <- gsub("[{}]", "", factor_string)

    # Extract the variable name (everything before <=, >=, <, >, ==, !=, or space)
    base_name <- sub("^([a-zA-Z_][a-zA-Z0-9_.]*).*", "\\1", factor_string)

    return(base_name)
  }

  # Apply extraction to all factors
  base_factors <- sapply(all_factors, extract_base_factor, USE.NAMES = FALSE)

  # Count each unique base factor
  base_factor_counts <- table(base_factors)

  # Create base factor summary table
  base_factor_summary <- data.table::data.table(
    Factor = names(base_factor_counts),
    Count = as.integer(base_factor_counts),
    Percent = 100 * as.integer(base_factor_counts) / n_found
  )

  # Sort by frequency
  base_factor_summary <- base_factor_summary[order(-Count)]
  base_factor_summary[, Rank := .I]

  # Reorder columns
  data.table::setcolorder(base_factor_summary, c("Rank", "Factor", "Count", "Percent"))

  # =========================================================================
  # PART 2: SPECIFIC FACTOR DEFINITIONS (>= threshold%)
  # =========================================================================

  # Count each specific factor definition
  specific_factor_counts <- table(all_factors)

  # Create specific factor summary table
  specific_factor_summary <- data.table::data.table(
    Factor_Definition = names(specific_factor_counts),
    Count = as.integer(specific_factor_counts),
    Percent = 100 * as.integer(specific_factor_counts) / n_found
  )

  # Filter to only those appearing >= threshold%
  specific_factor_summary <- specific_factor_summary[Percent >= threshold]

  # Add base factor name for reference
  specific_factor_summary[, Base_Factor := sapply(Factor_Definition, extract_base_factor)]

  # Sort by base factor, then by frequency
  specific_factor_summary <- specific_factor_summary[order(Base_Factor, -Count)]
  specific_factor_summary[, Rank := .I]

  # Reorder columns
  data.table::setcolorder(specific_factor_summary,
                          c("Rank", "Base_Factor", "Factor_Definition", "Count", "Percent"))

  # =========================================================================
  # RETURN BOTH SUMMARIES
  # =========================================================================

  return(list(
    base_factors = base_factor_summary,
    specific_factors = specific_factor_summary,
    threshold = threshold
  ))
}
