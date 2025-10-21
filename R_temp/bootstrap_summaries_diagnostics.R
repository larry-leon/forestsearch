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


#' Enhanced Bootstrap Results Summary (WITH TIMING)
#'
#' Creates comprehensive output including formatted table, diagnostic plots,
#' bootstrap quality metrics, and detailed timing analysis.
#'
#' @param boot_results List. Output from forestsearch_bootstrap_dofuture()
#' @param create_plots Logical. Generate diagnostic plots (default: TRUE)
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#'
#' @return List with formatted table, diagnostics, timing analysis, and plots
#' @export
summarize_bootstrap_results <- function(boot_results, create_plots = FALSE,
                                        est.scale = "hr") {

  # =========================================================================
  # SECTION 1: EXTRACT COMPONENTS
  # =========================================================================

  # Extract components
  FSsg_tab <- boot_results$FSsg_tab
  results <- boot_results$results
  H_estimates <- boot_results$H_estimates
  Hc_estimates <- boot_results$Hc_estimates

  # =========================================================================
  # SECTION 2: CALCULATE BOOTSTRAP SUCCESS METRICS
  # =========================================================================

  # Calculate bootstrap success rate
  boot_success_rate <- mean(!is.na(results$H_biasadj_2))
  nb_boots <- nrow(results)

  # Create formatted table
  formatted_table <- format_bootstrap_table(
    FSsg_tab = FSsg_tab,
    nb_boots = nb_boots,
    est.scale = est.scale,
    boot_success_rate = boot_success_rate
  )

  # =========================================================================
  # SECTION 3: EXTRACT AND ANALYZE TIMING INFORMATION
  # =========================================================================

  # Check if timing information is available
  has_timing <- !is.null(attr(results, "timing"))

  if (has_timing) {
    # Overall timing from attributes
    overall_timing <- attr(results, "timing")

    # Per-iteration statistics (if columns exist)
    has_iteration_timing <- "tmins_iteration" %in% names(results)
    has_search_timing <- "tmins_search" %in% names(results)

    if (has_iteration_timing) {
      iteration_stats <- list(
        mean = mean(results$tmins_iteration, na.rm = TRUE),
        median = median(results$tmins_iteration, na.rm = TRUE),
        sd = sd(results$tmins_iteration, na.rm = TRUE),
        min = min(results$tmins_iteration, na.rm = TRUE),
        max = max(results$tmins_iteration, na.rm = TRUE),
        q25 = quantile(results$tmins_iteration, 0.25, na.rm = TRUE),
        q75 = quantile(results$tmins_iteration, 0.75, na.rm = TRUE)
      )
    } else {
      iteration_stats <- NULL
    }

    # ForestSearch timing (subset where forestsearch ran)
    if (has_search_timing) {
      fs_times <- results$tmins_search[!is.na(results$tmins_search)]
      if (length(fs_times) > 0) {
        fs_stats <- list(
          n_runs = length(fs_times),
          pct_runs = length(fs_times) / nb_boots * 100,
          mean = mean(fs_times),
          median = median(fs_times),
          sd = sd(fs_times),
          min = min(fs_times),
          max = max(fs_times),
          total = sum(fs_times)
        )
      } else {
        fs_stats <- NULL
      }
    } else {
      fs_stats <- NULL
    }

    # Overhead timing (time not in forestsearch)
    if (has_iteration_timing && has_search_timing) {
      overhead_times <- results$tmins_iteration - results$tmins_search
      overhead_times <- overhead_times[!is.na(overhead_times)]

      if (length(overhead_times) > 0) {
        overhead_stats <- list(
          mean = mean(overhead_times),
          median = median(overhead_times),
          total = sum(overhead_times),
          pct_of_total = sum(overhead_times) / overall_timing$total_minutes * 100
        )
      } else {
        overhead_stats <- NULL
      }
    } else {
      overhead_stats <- NULL
    }

  } else {
    # No timing information available
    overall_timing <- NULL
    iteration_stats <- NULL
    fs_stats <- NULL
    overhead_stats <- NULL
  }

  # =========================================================================
  # SECTION 4: BOOTSTRAP DIAGNOSTICS
  # =========================================================================

  diagnostics <- list(
    n_boots = nb_boots,
    success_rate = boot_success_rate,
    n_successful = sum(!is.na(results$H_biasadj_2)),
    n_failed = sum(is.na(results$H_biasadj_2))
  )

  # Add timing to diagnostics if available
  if (has_timing && has_search_timing) {
    diagnostics$median_search_time = median(results$tmins_search, na.rm = TRUE)
    diagnostics$total_search_time = sum(results$tmins_search, na.rm = TRUE)
  }

  # =========================================================================
  # SECTION 5: PRINT COMPREHENSIVE SUMMARY (WITH FIXED sprintf)
  # =========================================================================

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("           BOOTSTRAP ANALYSIS SUMMARY                          \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  # Success metrics
  cat("BOOTSTRAP SUCCESS METRICS:\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat(sprintf("  Total iterations:              %d\n", diagnostics$n_boots))
  cat(sprintf("  Successful subgroup ID:        %d (%.1f%%)\n",
              diagnostics$n_successful, diagnostics$success_rate * 100))
  cat(sprintf("  Failed to find subgroup:       %d (%.1f%%)\n",
              diagnostics$n_failed, (1 - diagnostics$success_rate) * 100))
  cat("\n")

  # Timing summary
  if (has_timing) {
    cat("TIMING ANALYSIS:\n")
    cat("─────────────────────────────────────────────────────────────\n")

    # Overall timing
    cat("Overall:\n")
    cat(sprintf("  Total bootstrap time:          %.2f minutes (%.2f hours)\n",
                overall_timing$total_minutes, overall_timing$total_hours))
    cat(sprintf("  Average per iteration:         %.2f min (%.1f sec)\n",
                overall_timing$avg_minutes_per_boot,
                overall_timing$avg_seconds_per_boot))

    # Projected times
    if (nb_boots < 1000) {
      projected_1000 <- overall_timing$avg_minutes_per_boot * 1000
      cat(sprintf("  Projected for 1000 boots:      %.2f minutes (%.2f hours)\n",
                  projected_1000, projected_1000 / 60))
    }
    cat("\n")

    # Per-iteration timing
    if (!is.null(iteration_stats)) {
      cat("Per-iteration timing:\n")
      cat(sprintf("  Mean:                          %.2f min (%.1f sec)\n",
                  iteration_stats$mean, iteration_stats$mean * 60))
      cat(sprintf("  Median:                        %.2f min (%.1f sec)\n",
                  iteration_stats$median, iteration_stats$median * 60))
      cat(sprintf("  Std Dev:                       %.2f minutes\n",
                  iteration_stats$sd))
      cat(sprintf("  Range:                         [%.2f, %.2f] minutes\n",
                  iteration_stats$min, iteration_stats$max))
      cat(sprintf("  IQR:                           [%.2f, %.2f] minutes\n",
                  iteration_stats$q25, iteration_stats$q75))
      cat("\n")
    }

    # ForestSearch timing (FIXED sprintf statements)
    if (!is.null(fs_stats)) {
      cat("ForestSearch timing (successful iterations only):\n")
      cat(sprintf("  Iterations with FS:            %d (%.1f%%)\n",
                  fs_stats$n_runs, fs_stats$pct_runs))
      cat(sprintf("  Mean FS time:                  %.2f min (%.1f sec)\n",
                  fs_stats$mean, fs_stats$mean * 60))
      cat(sprintf("  Median FS time:                %.2f min (%.1f sec)\n",
                  fs_stats$median, fs_stats$median * 60))
      cat(sprintf("  Total FS time:                 %.2f minutes\n",
                  fs_stats$total))
      # FIXED: Changed from %.1f%% to %.1f and added separate %
      cat(sprintf("  FS time %s of total:            %.1f%s\n",
                  "%", fs_stats$total / overall_timing$total_minutes * 100, "%"))
      cat("\n")
    }

    # Overhead timing (FIXED sprintf statement)
    if (!is.null(overhead_stats)) {
      cat("Overhead timing (Cox models, bias correction, etc.):\n")
      cat(sprintf("  Mean overhead:                 %.2f min (%.1f sec)\n",
                  overhead_stats$mean, overhead_stats$mean * 60))
      cat(sprintf("  Median overhead:               %.2f min (%.1f sec)\n",
                  overhead_stats$median, overhead_stats$median * 60))
      cat(sprintf("  Total overhead:                %.2f minutes\n",
                  overhead_stats$total))
      # FIXED: Changed from %.1f%% to %.1f and added separate %
      cat(sprintf("  Overhead %s of total:           %.1f%s\n",
                  "%", overhead_stats$pct_of_total, "%"))
      cat("\n")
    }

    # Performance assessment
    cat("PERFORMANCE ASSESSMENT:\n")
    cat("─────────────────────────────────────────────────────────────\n")

    if (!is.null(iteration_stats)) {
      avg_sec <- overall_timing$avg_seconds_per_boot

      if (avg_sec < 5) {
        performance <- "Excellent"
        emoji <- "✓✓✓"
      } else if (avg_sec < 15) {
        performance <- "Good"
        emoji <- "✓✓"
      } else if (avg_sec < 30) {
        performance <- "Acceptable"
        emoji <- "✓"
      } else if (avg_sec < 60) {
        performance <- "Slow"
        emoji <- "⚠"
      } else {
        performance <- "Very Slow"
        emoji <- "⚠⚠"
      }

      cat(sprintf("  Performance rating:            %s %s\n", emoji, performance))
      cat(sprintf("  Average iteration speed:       %.1f seconds\n", avg_sec))

      # Recommendations based on performance
      if (avg_sec > 30) {
        cat("\n  Recommendations for improvement:\n")
        if (!is.null(fs_stats) && fs_stats$mean > 0.5) {
          cat("    • Consider reducing max.minutes in forestsearch\n")
          cat("    • Consider reducing maxk if currently > 2\n")
        }
        if (nb_boots > 500) {
          cat("    • Consider reducing nb_boots for initial testing\n")
        }
        cat("    • Ensure sufficient parallel workers are allocated\n")
      }
    }
    cat("\n")
  }

  cat("═══════════════════════════════════════════════════════════════\n\n")

  # =========================================================================
  # SECTION 6: CREATE DIAGNOSTIC PLOTS (INCLUDING TIMING PLOTS)
  # =========================================================================

  plots <- NULL
  if (create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- create_bootstrap_diagnostic_plots(results, H_estimates, Hc_estimates)

    # Add timing plots if data available
    if (has_timing && has_iteration_timing) {
      plots$timing <- create_bootstrap_timing_plots(results)
    }

    # Add combined plot if patchwork is available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      plots$combined <- (plots$H_distribution | plots$Hc_distribution) +
        patchwork::plot_annotation(
          title = "Bootstrap Distributions",
          subtitle = sprintf("%d iterations, %.1f%s successful",
                             nrow(results), boot_success_rate * 100, "%"),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold"))
        )

      # Add timing plots to combined if available
      if (!is.null(plots$timing)) {
        plots$combined_with_timing <- (plots$H_distribution | plots$Hc_distribution) /
          (plots$timing$iteration_time | plots$timing$search_time) +
          patchwork::plot_annotation(
            title = "Bootstrap Analysis: Distributions and Timing",
            subtitle = sprintf("%d iterations, %.1f%s successful, %.1f min total",
                               nrow(results), boot_success_rate * 100, "%",
                               overall_timing$total_minutes),
            theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold"))
          )
      }
    }
  }

  # =========================================================================
  # SECTION 7: COMPILE AND RETURN OUTPUT
  # =========================================================================

  output <- list(
    table = formatted_table,
    diagnostics = diagnostics,
    plots = plots
  )

  # Add timing information if available
  if (has_timing) {
    output$timing <- list(
      overall = overall_timing,
      iteration_stats = iteration_stats,
      fs_stats = fs_stats,
      overhead_stats = overhead_stats
    )
  }

  invisible(output)
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

