#' Enhanced Bootstrap Results Summary
#'
#' Creates comprehensive output including formatted table with subgroup footnote,
#' diagnostic plots, bootstrap quality metrics, and detailed timing analysis.
#'
#' @param sgharm The selected subgroup object from forestsearch results. Can be:
#'   \itemize{
#'     \item Character vector of factor definitions (e.g., c("\{age>=50\}", "\{nodes>=3\}"))
#'     \item List with `sgharm` element containing factor definitions
#'     \item List with `sg.harm_label` element (human-readable labels)
#'   }
#' @param boot_results List. Output from forestsearch_bootstrap_dofuture()
#' @param create_plots Logical. Generate diagnostic plots (default: FALSE)
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#'
#' @return List with components:
#'   \describe{
#'     \item{table}{gt table with treatment effects and subgroup footnote}
#'     \item{diagnostics}{List of bootstrap quality metrics}
#'     \item{diagnostics_table_gt}{gt table of diagnostics}
#'     \item{plots}{List of ggplot2 diagnostic plots (if create_plots=TRUE)}
#'     \item{timing}{List of timing analysis (if timing data available)}
#'     \item{subgroup_summary}{List from summarize_bootstrap_subgroups()}
#'   }
#'
#' @details
#' The `table` output includes a footnote displaying the identified subgroup
#' definition, analogous to the `tab_estimates` table from \code{\link{sg_tables}}.
#' This is achieved by extracting the subgroup definition from `sgharm` and
#' passing it to \code{\link{format_bootstrap_table}}.
#'
#' @seealso
#' \code{\link{format_bootstrap_table}} for table creation
#' \code{\link{sg_tables}} for analogous main analysis tables
#' \code{\link{summarize_bootstrap_subgroups}} for subgroup stability analysis
#'
#' @export
summarize_bootstrap_results <- function(sgharm,
                                        boot_results,
                                        create_plots = FALSE,
                                        est.scale = "hr") {

  # ===========================================================================
  # SECTION 1: EXTRACT COMPONENTS
  # ===========================================================================

  FSsg_tab <- boot_results$FSsg_tab
  results <- boot_results$results
  H_estimates <- boot_results$H_estimates
  Hc_estimates <- boot_results$Hc_estimates

  # ===========================================================================
  # SECTION 2: CALCULATE BOOTSTRAP SUCCESS METRICS
  # ===========================================================================

  boot_success_rate <- mean(!is.na(results$H_biasadj_2))
  nb_boots <- nrow(results)

  # ===========================================================================
  # SECTION 3: EXTRACT SUBGROUP DEFINITION FOR FOOTNOTE
 # ===========================================================================

  # Extract subgroup definition string from sgharm
  # This mirrors the approach in sg_tables() for tab_estimates
  sg_definition <- NULL

  if (!is.null(sgharm)) {
    if (is.character(sgharm)) {
      # Direct character vector: c("{age>=50}", "{nodes>=3}")
      sg_valid <- sgharm[!is.na(sgharm) & sgharm != ""]
      if (length(sg_valid) > 0) {
        sg_definition <- paste(sg_valid, collapse = " & ")
      }

    } else if (is.factor(sgharm)) {
      # Factor: convert to character first
      sg_valid <- as.character(sgharm)
      sg_valid <- sg_valid[!is.na(sg_valid) & sg_valid != ""]
      if (length(sg_valid) > 0) {
        sg_definition <- paste(sg_valid, collapse = " & ")
      }

    } else if (is.list(sgharm)) {
      # List structure: check for common element names
      # Priority: sg.harm_label > sgharm > sg.harm

      if (!is.null(sgharm$sg.harm_label)) {
        # Human-readable labels (preferred)
        sg_valid <- sgharm$sg.harm_label
        sg_valid <- sg_valid[!is.na(sg_valid) & sg_valid != ""]
        if (length(sg_valid) > 0) {
          sg_definition <- paste(sg_valid, collapse = " & ")
        }

      } else if (!is.null(sgharm$sgharm)) {
        # Technical factor names
        sg_valid <- sgharm$sgharm
        sg_valid <- sg_valid[!is.na(sg_valid) & sg_valid != ""]
        if (length(sg_valid) > 0) {
          sg_definition <- paste(sg_valid, collapse = " & ")
        }

      } else if (!is.null(sgharm$sg.harm)) {
        # Alternative naming
        sg_valid <- sgharm$sg.harm
        sg_valid <- sg_valid[!is.na(sg_valid) & sg_valid != ""]
        if (length(sg_valid) > 0) {
          sg_definition <- paste(sg_valid, collapse = " & ")
        }
      }
    }
  }

  # ===========================================================================
  # SECTION 4: CREATE FORMATTED TABLE WITH SUBGROUP FOOTNOTE
  # ===========================================================================

  formatted_table <- format_bootstrap_table(
    FSsg_tab = FSsg_tab,
    nb_boots = nb_boots,
    est.scale = est.scale,
    boot_success_rate = boot_success_rate,
    sg_definition = sg_definition  # NEW: pass subgroup definition
  )

  # ===========================================================================
  # SECTION 5: EXTRACT AND ANALYZE TIMING INFORMATION
  # ===========================================================================

  has_timing <- !is.null(attr(results, "timing"))
  overall_timing <- NULL
  iteration_stats <- NULL
  fs_stats <- NULL
  overhead_stats <- NULL
  timing_summary_table <- NULL
  timing_table_gt <- NULL

  if (has_timing) {
    overall_timing <- attr(results, "timing")

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
    }

    if (has_search_timing) {
      fs_times <- results$tmins_search[!is.na(results$tmins_search)]
      if (length(fs_times) > 0) {
        fs_stats <- list(
          n_runs = length(fs_times),
          pct_runs = 100 * length(fs_times) / nb_boots,
          mean = mean(fs_times),
          median = median(fs_times),
          sd = sd(fs_times),
          min = min(fs_times),
          max = max(fs_times),
          total = sum(fs_times)
        )
      }
    }

    if (has_iteration_timing && has_search_timing) {
      overhead <- results$tmins_iteration - results$tmins_search
      overhead <- overhead[!is.na(overhead)]
      if (length(overhead) > 0) {
        overhead_stats <- list(
          mean = mean(overhead),
          median = median(overhead),
          total = sum(overhead),
          pct_of_total = 100 * sum(overhead) / overall_timing$total_minutes
        )
      }
    }

    # Create timing summary table
    timing_summary_table <- tryCatch({
      create_timing_summary_table(
        overall_timing = overall_timing,
        iteration_stats = iteration_stats,
        fs_stats = fs_stats,
        overhead_stats = overhead_stats,
        nb_boots = nb_boots,
        boot_success_rate = boot_success_rate
      )
    }, error = function(e) {
      warning("Could not create timing summary table: ", e$message)
      NULL
    })

    # Create gt timing table
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

  # ===========================================================================
  # SECTION 6: CREATE DIAGNOSTICS
  # ===========================================================================

  diagnostics <- list(
    n_boots = nb_boots,
    n_successful = sum(!is.na(results$H_biasadj_2)),
    n_failed = sum(is.na(results$H_biasadj_2)),
    success_rate = boot_success_rate
  )

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

  # ===========================================================================
  # SECTION 7: CREATE PLOTS (if requested)
  # ===========================================================================

  plots <- NULL
  if (create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- tryCatch({
      create_bootstrap_diagnostic_plots(
        results = results,
        H_estimates = H_estimates,
        Hc_estimates = Hc_estimates,
        overall_timing = overall_timing
      )
    }, error = function(e) {
      warning("Could not create diagnostic plots: ", e$message)
      NULL
    })
  }

  # ===========================================================================
  # SECTION 8: PRINT SUMMARY
  # ===========================================================================

  cat("\n")
  cat("===============================================================\n")
  cat("           BOOTSTRAP ANALYSIS SUMMARY                          \n")
  cat("===============================================================\n\n")

  # Subgroup definition (NEW)
  if (!is.null(sg_definition)) {
    cat("IDENTIFIED SUBGROUP:\n")
    cat("-------------------------------------------------------------\n")
    cat(sprintf("  H: %s\n\n", sg_definition))
  }

  # Success metrics
  cat("BOOTSTRAP SUCCESS METRICS:\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  Total iterations:              %d\n", diagnostics$n_boots))
  cat(sprintf("  Successful subgroup ID:        %d (%.1f%%)\n",
              diagnostics$n_successful, diagnostics$success_rate * 100))
  cat(sprintf("  Failed to find subgroup:       %d (%.1f%%)\n",
              diagnostics$n_failed, (1 - diagnostics$success_rate) * 100))
  cat("\n")

  # Timing summary
  if (has_timing) {
    cat("TIMING ANALYSIS:\n")
    cat("-------------------------------------------------------------\n")
    cat("Overall:\n")
    cat(sprintf("  Total bootstrap time:          %.2f minutes (%.2f hours)\n",
                overall_timing$total_minutes, overall_timing$total_hours))
    cat(sprintf("  Average per iteration:         %.2f min (%.1f sec)\n",
                overall_timing$avg_minutes_per_boot,
                overall_timing$avg_seconds_per_boot))

    if (nb_boots < 1000) {
      projected_min <- overall_timing$avg_minutes_per_boot * 1000
      cat(sprintf("  Projected for 1000 boots:      %.2f min (%.2f hrs)\n",
                  projected_min, projected_min / 60))
    }
    cat("\n")
  }

  # ===========================================================================
  # SECTION 9: CREATE SUBGROUP SUMMARY
  # ===========================================================================

  subgroup_summary <- NULL

  if (!is.null(results) && "Pcons" %in% names(results)) {

    # Prepare original_sg for summarize_bootstrap_subgroups
    original_sg <- NULL
    if (!is.null(sgharm)) {
      if (is.character(sgharm) || is.factor(sgharm)) {
        original_sg <- as.character(sgharm)
      } else if (is.list(sgharm)) {
        if (!is.null(sgharm$sgharm)) {
          original_sg <- sgharm$sgharm
        } else if (!is.null(sgharm$sg.harm)) {
          original_sg <- sgharm$sg.harm
        } else if (!is.null(sgharm$sg.harm_label)) {
          original_sg <- sgharm$sg.harm_label
        }
      }
    }

    # Extract maxk from boot_results attributes or detect from columns
    maxk <- attr(boot_results, "maxk")
    if (is.null(maxk)) {
      maxk <- sum(c("M.1", "M.2", "M.3", "M.4", "M.5") %in% names(results))
      if (maxk == 0) maxk <- 2
    }

    subgroup_summary <- tryCatch({
      summarize_bootstrap_subgroups(
        results = results,
        nb_boots = nb_boots,
        original_sg = original_sg,
        maxk = maxk
      )
    }, error = function(e) {
      warning("Failed to create subgroup summary: ", e$message)
      NULL
    })
  }

  # ===========================================================================
  # SECTION 10: COMPILE AND RETURN OUTPUT
  # ===========================================================================

  output <- list(
    table = formatted_table,
    diagnostics = diagnostics,
    diagnostics_table_gt = diagnostics_table_gt,
    plots = plots,
    sg_definition = sg_definition,  # NEW: include for reference
    subgroup_summary = subgroup_summary
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

  invisible(output)
}
