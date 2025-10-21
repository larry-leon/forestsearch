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
summarize_bootstrap_results <- function(sgharm, boot_results, create_plots = FALSE,
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
  # CREATE TIMING SUMMARY TABLE
  # =========================================================================

  timing_summary_table <- NULL

  if (has_timing) {
    timing_rows <- list()

    # Overall timing
    if (!is.null(overall_timing)) {
      timing_rows$overall <- data.frame(
        Category = "Overall",
        Metric = c("Total time (min)", "Total time (hours)", "Avg per boot (min)", "Avg per boot (sec)"),
        Value = c(
          sprintf("%.2f", overall_timing$total_minutes),
          sprintf("%.2f", overall_timing$total_hours),
          sprintf("%.2f", overall_timing$avg_minutes_per_boot),
          sprintf("%.1f", overall_timing$avg_seconds_per_boot)
        ),
        stringsAsFactors = FALSE
      )
    }

    # Per-iteration statistics
    if (!is.null(iteration_stats)) {
      timing_rows$iteration <- data.frame(
        Category = "Per-Iteration",
        Metric = c("Mean (min)", "Median (min)", "SD (min)", "Min (min)", "Max (min)", "Q25 (min)", "Q75 (min)"),
        Value = c(
          sprintf("%.2f", iteration_stats$mean),
          sprintf("%.2f", iteration_stats$median),
          sprintf("%.2f", iteration_stats$sd),
          sprintf("%.2f", iteration_stats$min),
          sprintf("%.2f", iteration_stats$max),
          sprintf("%.2f", iteration_stats$q25),
          sprintf("%.2f", iteration_stats$q75)
        ),
        stringsAsFactors = FALSE
      )
    }

    # ForestSearch timing
    if (!is.null(fs_stats)) {
      timing_rows$forestsearch <- data.frame(
        Category = "ForestSearch",
        Metric = c("Runs", "% of iterations", "Mean (min)", "Median (min)", "Total (min)", "% of total time"),
        Value = c(
          sprintf("%d", fs_stats$n_runs),
          sprintf("%.1f%%", fs_stats$pct_runs),
          sprintf("%.2f", fs_stats$mean),
          sprintf("%.2f", fs_stats$median),
          sprintf("%.2f", fs_stats$total),
          sprintf("%.1f%%", fs_stats$total / overall_timing$total_minutes * 100)
        ),
        stringsAsFactors = FALSE
      )
    }

    # Overhead timing
    if (!is.null(overhead_stats)) {
      timing_rows$overhead <- data.frame(
        Category = "Overhead",
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
    }
  }

  # =========================================================================
  # CREATE FORMATTED GT TIMING TABLE
  # =========================================================================

  timing_table_gt <- NULL
  if (requireNamespace("gt", quietly = TRUE)) {
    timing_table_gt <- tryCatch({
      format_bootstrap_timing_table(
        timing_list = list(
          overall = overall_timing,
          iteration_stats = iteration_stats,
          fs_stats = fs_stats,
          overhead_stats = overhead_stats
        ),
        nb_boots = nb_boots,
        boot_success_rate = boot_success_rate
      )
    }, error = function(e) {
      warning("Could not create gt timing table: ", e$message)
      NULL
    })
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
  # CREATE FORMATTED GT DIAGNOSTICS TABLE
  # =========================================================================

  diagnostics_table_gt <- NULL
  if (requireNamespace("gt", quietly = TRUE)) {
    diagnostics_table_gt <- tryCatch({
      format_bootstrap_diagnostics_table(
        diagnostics = diagnostics,
        nb_boots = nb_boots,
        results = results,
        H_estimates = H_estimates,
        Hc_estimates = Hc_estimates
      )
    }, error = function(e) {
      warning("Could not create gt diagnostics table: ", e$message)
      NULL
    })
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
    diagnostics_table_gt = diagnostics_table_gt,
    plots = plots
  )

  # Add timing information if available
  if (has_timing) {
    output$timing <- list(
      overall = overall_timing,
      iteration_stats = iteration_stats,
      fs_stats = fs_stats,
      overhead_stats = overhead_stats,
      time_table = timing_summary_table,
      time_table_gt = timing_table_gt
    )
  }

  subgroup_summary <- NULL
  if (!is.null(results) && "Pcons" %in% names(results)) {

    # Extract original subgroup definition if available
    original_sg <- NULL
    if (!is.null(boot_results$FSsg_tab)) {
      # Try to extract from FSsg_tab
      # This assumes the subgroup definition is available
    }

    subgroup_summary <- summarize_bootstrap_subgroups(
      results = results,
      nb_boots = nb_boots,
      original_sg = sgharm,
      maxk = 2  # Or extract from boot_results if stored
    )
  }

  # Add to output
  output$subgroup_summary <- subgroup_summary


  invisible(output)
}
