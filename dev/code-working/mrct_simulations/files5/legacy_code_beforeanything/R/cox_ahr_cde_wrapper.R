#' Comprehensive Wrapper for Cox Spline Analysis with AHR and CDE Plotting
#'
#' This wrapper function combines Cox spline fitting with comprehensive visualization
#' of Average Hazard Ratios (AHRs) and Controlled Direct Effects (CDEs) as described
#' in the MRCT subgroups analysis documentation.
#'
#' @param df Data frame containing survival data with potential outcomes
#' @param tte_name Character string specifying time-to-event variable name. Default: "os_time"
#' @param event_name Character string specifying event indicator variable name. Default: "os_event"
#' @param treat_name Character string specifying treatment variable name. Default: "treat"
#' @param z_name Character string specifying continuous covariate/biomarker name. Default: "biomarker"
#' @param loghr_po_name Character string specifying potential outcome log HR variable. Default: "loghr_po"
#' @param theta1_name Optional: variable name for theta_1 (treated potential outcome). Default: "theta_1"
#' @param theta0_name Optional: variable name for theta_0 (control potential outcome). Default: "theta_0"
#' @param spline_df Integer degrees of freedom for spline fitting. Default: 3
#' @param alpha Numeric significance level for confidence intervals. Default: 0.20
#' @param hr_threshold Numeric hazard ratio threshold for subgroup identification. Default: 0.7
#' @param plot_style Character: "combined", "separate", or "grid" for plot layout. Default: "combined"
#' @param save_plots Logical whether to save plots to file. Default: FALSE
#' @param output_dir Character directory for saving plots. Default: "plots"
#' @param verbose Logical for diagnostic output. Default: TRUE
#'
#' @return List containing:
#' \describe{
#'   \item{cox_fit}{Results from cox_cs_fit function}
#'   \item{ahr_results}{AHR calculations for different subgroup definitions}
#'   \item{cde_results}{CDE calculations if theta variables available}
#'   \item{optimal_cutpoint}{Optimal biomarker cutpoint for subgroup}
#'   \item{subgroup_stats}{Statistics for recommended and questionable subgroups}
#'   \item{plots}{List of generated plot objects (if using ggplot2)}
#' }
#'
#' @examples
#' \dontrun{
#' # Load your ForestSearch output data
#' results <- cox_ahr_cde_analysis(
#'   df = df_nonAP,
#'   z_name = "biomarker",
#'   hr_threshold = 0.7,
#'   plot_style = "grid"
#' )
#'
#' # Access specific results
#' print(results$subgroup_stats)
#' print(results$optimal_cutpoint)
#' }
#'
#' @importFrom survival coxph Surv
#' @importFrom splines ns
#' @importFrom graphics par plot points lines polygon abline legend hist rug axis text
#' @importFrom grDevices rgb dev.copy dev.off pdf
#' @importFrom stats density quantile sd median
#' @export

cox_ahr_cde_analysis <- function(
    df,
    tte_name = "os_time",
    event_name = "os_event",
    treat_name = "treat",
    z_name = "biomarker",
    loghr_po_name = "loghr_po",
    theta1_name = "theta_1",
    theta0_name = "theta_0",
    spline_df = 3,
    alpha = 0.20,
    hr_threshold = 0.7,
    plot_style = c("combined", "separate", "grid"),
    save_plots = FALSE,
    output_dir = "plots",
    verbose = TRUE
) {

  plot_style <- match.arg(plot_style)

  # Validate input data
  required_cols <- c(tte_name, event_name, treat_name, z_name, loghr_po_name)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check for CDE variables
  has_cde <- all(c(theta1_name, theta0_name) %in% names(df))

  if (verbose) {
    cat("\n=====================================\n")
    cat("Cox AHR/CDE Analysis\n")
    cat("=====================================\n")
    cat("Sample size:", nrow(df), "\n")
    cat("Events:", sum(df[[event_name]]), "\n")
    cat("Treatment allocation:", table(df[[treat_name]]), "\n")
    cat("Biomarker range:", round(range(df[[z_name]], na.rm = TRUE), 2), "\n")
    cat("CDE variables available:", has_cde, "\n")
    cat("-------------------------------------\n\n")
  }

  # ===========================================================================
  # Step 1: Fit Cox Spline Model
  # ===========================================================================

  if (verbose) cat("Step 1: Fitting Cox spline model...\n")

  cox_results <- cox_cs_fit(
    df = df,
    tte_name = tte_name,
    event_name = event_name,
    treat_name = treat_name,
    z_name = z_name,
    alpha = alpha,
    spline_df = spline_df,
    show_plot = FALSE,
    verbose = verbose
  )

  # ===========================================================================
  # Step 2: Calculate AHR and CDE Metrics
  # ===========================================================================

  if (verbose) cat("\nStep 2: Calculating AHR and CDE metrics...\n")

  # Sort data by biomarker
  df_sorted <- df[order(df[[z_name]]), ]
  z_values <- df_sorted[[z_name]]
  loghr_po <- df_sorted[[loghr_po_name]]

  # Find optimal cutpoint (where treatment effect becomes beneficial)
  log_hr_threshold <- log(hr_threshold)
  beneficial_idx <- which(loghr_po < log_hr_threshold)

  if (length(beneficial_idx) > 0) {
    optimal_cutpoint <- min(z_values[beneficial_idx])
  } else {
    optimal_cutpoint <- max(z_values)
    warning("No subgroup with HR < ", hr_threshold, " found. Using max biomarker value.")
  }

  # Define subgroups
  df_sorted$subgroup <- ifelse(
    z_values >= optimal_cutpoint,
    "Recommend",
    "Questionable"
  )

  # Calculate AHR for different subgroup definitions
  ahr_results <- list()

  # Overall AHR (empirical)
  ahr_results$overall <- exp(mean(loghr_po))

  # AHR for recommended subgroup
  recommend_idx <- df_sorted$subgroup == "Recommend"
  ahr_results$recommend <- exp(mean(loghr_po[recommend_idx]))

  # AHR for questionable subgroup
  question_idx <- df_sorted$subgroup == "Questionable"
  ahr_results$questionable <- exp(mean(loghr_po[question_idx]))

  # Calculate AHR curves across biomarker range
  z_grid <- seq(min(z_values), max(z_values), length.out = 100)
  ahr_curve_pos <- numeric(length(z_grid))
  ahr_curve_neg <- numeric(length(z_grid))

  for (i in seq_along(z_grid)) {
    z_cut <- z_grid[i]
    idx_pos <- z_values >= z_cut
    idx_neg <- z_values <= z_cut

    if (sum(idx_pos) > 0) {
      ahr_curve_pos[i] <- exp(mean(loghr_po[idx_pos]))
    } else {
      ahr_curve_pos[i] <- NA
    }

    if (sum(idx_neg) > 0) {
      ahr_curve_neg[i] <- exp(mean(loghr_po[idx_neg]))
    } else {
      ahr_curve_neg[i] <- NA
    }
  }

  ahr_results$z_grid <- z_grid
  ahr_results$ahr_curve_pos <- ahr_curve_pos
  ahr_results$ahr_curve_neg <- ahr_curve_neg

  # Calculate CDE if available
  cde_results <- NULL
  if (has_cde) {
    if (verbose) cat("Calculating Controlled Direct Effects (CDEs)...\n")

    theta1 <- df_sorted[[theta1_name]]
    theta0 <- df_sorted[[theta0_name]]

    cde_results <- list()

    # Overall CDE
    cde_results$overall <- mean(exp(theta1)) / mean(exp(theta0))

    # CDE for subgroups
    cde_results$recommend <- mean(exp(theta1[recommend_idx])) / mean(exp(theta0[recommend_idx]))
    cde_results$questionable <- mean(exp(theta1[question_idx])) / mean(exp(theta0[question_idx]))

    # CDE curves
    cde_curve_pos <- numeric(length(z_grid))
    cde_curve_neg <- numeric(length(z_grid))

    for (i in seq_along(z_grid)) {
      z_cut <- z_grid[i]
      idx_pos <- z_values >= z_cut
      idx_neg <- z_values <= z_cut

      if (sum(idx_pos) > 0) {
        cde_curve_pos[i] <- mean(exp(theta1[idx_pos])) / mean(exp(theta0[idx_pos]))
      } else {
        cde_curve_pos[i] <- NA
      }

      if (sum(idx_neg) > 0) {
        cde_curve_neg[i] <- mean(exp(theta1[idx_neg])) / mean(exp(theta0[idx_neg]))
      } else {
        cde_curve_neg[i] <- NA
      }
    }

    cde_results$cde_curve_pos <- cde_curve_pos
    cde_results$cde_curve_neg <- cde_curve_neg
  }

  # ===========================================================================
  # Step 3: Calculate Subgroup Statistics
  # ===========================================================================

  if (verbose) cat("\nStep 3: Calculating subgroup statistics...\n")

  subgroup_stats <- data.frame(
    Subgroup = c("Questionable", "Recommend", "Overall"),
    n = c(
      sum(question_idx),
      sum(recommend_idx),
      nrow(df_sorted)
    ),
    n_pct = c(
      round(100 * sum(question_idx) / nrow(df_sorted), 1),
      round(100 * sum(recommend_idx) / nrow(df_sorted), 1),
      100
    ),
    events = c(
      sum(df_sorted[[event_name]][question_idx]),
      sum(df_sorted[[event_name]][recommend_idx]),
      sum(df_sorted[[event_name]])
    ),
    AHR_po = round(c(
      ahr_results$questionable,
      ahr_results$recommend,
      ahr_results$overall
    ), 3)
  )

  if (has_cde) {
    subgroup_stats$CDE <- round(c(
      cde_results$questionable,
      cde_results$recommend,
      cde_results$overall
    ), 3)
  }

  # Add Cox model estimates for each subgroup if desired
  # This would require fitting separate Cox models per subgroup

  if (verbose) {
    cat("\nSubgroup Summary:\n")
    print(subgroup_stats)
    cat("\nOptimal cutpoint (", z_name, " >= ", round(optimal_cutpoint, 2),
        ") defines 'Recommend' subgroup\n", sep = "")
  }

  # ===========================================================================
  # Step 4: Create Comprehensive Plots
  # ===========================================================================

  if (verbose) cat("\nStep 4: Creating visualizations...\n")

  # Create output directory if saving
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  # Set up plot layout based on style
  if (plot_style == "combined") {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  } else if (plot_style == "grid") {
    if (has_cde) {
      par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))
    } else {
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    }
  }

  # Plot 1: Treatment Effect Profile (Spline Fit)
  plot(cox_results$z_profile, cox_results$loghr_est,
       type = "l", lwd = 2, col = "darkblue",
       xlab = paste(z_name),
       ylab = "log(HR)",
       main = "Cox Spline: Treatment Effect Profile",
       ylim = range(c(cox_results$loghr_lower, cox_results$loghr_upper, loghr_po),
                    na.rm = TRUE))

  # Add confidence bands
  polygon(c(cox_results$z_profile, rev(cox_results$z_profile)),
          c(cox_results$loghr_lower, rev(cox_results$loghr_upper)),
          col = rgb(0, 0, 1, 0.2), border = NA)

  # Add true values if available
  points(z_values, loghr_po, pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3))

  # Add reference lines
  abline(h = 0, lty = 2, col = "gray")
  abline(h = log(hr_threshold), lty = 2, col = "red", lwd = 2)
  abline(v = optimal_cutpoint, lty = 2, col = "green", lwd = 2)

  legend("topright",
         legend = c("Spline estimate", "95% CI", "True values",
                    paste("HR =", hr_threshold), "Optimal cut"),
         lty = c(1, NA, NA, 2, 2),
         pch = c(NA, NA, 16, NA, NA),
         col = c("darkblue", rgb(0,0,1,0.2), rgb(0,0,0,0.3), "red", "green"),
         lwd = c(2, NA, NA, 2, 2),
         cex = 0.7, bty = "n")

  # Plot 2: Average Hazard Ratio (z >= threshold)
  plot(z_grid, ahr_curve_pos,
       type = "l", lwd = 2, col = "darkgreen",
       xlab = paste(z_name, "threshold"),
       ylab = "AHR",
       main = paste("AHR for", z_name, ">= threshold"),
       ylim = range(c(0.4, 2, ahr_curve_pos), na.rm = TRUE))

  abline(h = 1, lty = 2, col = "gray")
  abline(h = hr_threshold, lty = 2, col = "red", lwd = 2)
  abline(v = optimal_cutpoint, lty = 2, col = "green", lwd = 2)

  # Add CDE curve if available
  if (has_cde) {
    lines(z_grid, cde_curve_pos, lwd = 2, col = "purple", lty = 2)
    legend("topright",
           legend = c("AHR(po)", "CDE", paste("HR =", hr_threshold), "Optimal cut"),
           lty = c(1, 2, 2, 2),
           col = c("darkgreen", "purple", "red", "green"),
           lwd = c(2, 2, 2, 2),
           cex = 0.7, bty = "n")
  } else {
    legend("topright",
           legend = c("AHR(po)", paste("HR =", hr_threshold), "Optimal cut"),
           lty = c(1, 2, 2),
           col = c("darkgreen", "red", "green"),
           lwd = c(2, 2, 2),
           cex = 0.7, bty = "n")
  }

  # Plot 3: Average Hazard Ratio (z <= threshold)
  plot(z_grid, ahr_curve_neg,
       type = "l", lwd = 2, col = "darkorange",
       xlab = paste(z_name, "threshold"),
       ylab = "AHR",
       main = paste("AHR for", z_name, "<= threshold"),
       ylim = range(c(0.4, 2, ahr_curve_neg), na.rm = TRUE))

  abline(h = 1, lty = 2, col = "gray")
  abline(h = hr_threshold, lty = 2, col = "red", lwd = 2)
  abline(v = optimal_cutpoint, lty = 2, col = "green", lwd = 2)

  if (has_cde) {
    lines(z_grid, cde_curve_neg, lwd = 2, col = "purple", lty = 2)
  }

  # Plot 4: Comparison of AHR methods
  plot(z_grid, ahr_curve_pos,
       type = "l", lwd = 2, col = "darkgreen",
       xlab = paste(z_name, "threshold"),
       ylab = "Hazard Ratio",
       main = "AHR vs CDE Comparison",
       ylim = range(c(0.4, 2, ahr_curve_pos, ahr_curve_neg), na.rm = TRUE))

  lines(z_grid, ahr_curve_neg, lwd = 2, col = "darkorange")

  if (has_cde) {
    lines(z_grid, cde_curve_pos, lwd = 2, col = "purple", lty = 2)
    lines(z_grid, cde_curve_neg, lwd = 2, col = "brown", lty = 2)
  }

  abline(h = 1, lty = 2, col = "gray")
  abline(h = hr_threshold, lty = 2, col = "red", lwd = 1)
  abline(v = optimal_cutpoint, lty = 2, col = "green", lwd = 2)

  legend_items <- c("AHR (z>=)", "AHR (z<=)")
  legend_cols <- c("darkgreen", "darkorange")
  legend_ltys <- c(1, 1)

  if (has_cde) {
    legend_items <- c(legend_items, "CDE (z>=)", "CDE (z<=)")
    legend_cols <- c(legend_cols, "purple", "brown")
    legend_ltys <- c(legend_ltys, 2, 2)
  }

  legend("topright",
         legend = legend_items,
         lty = legend_ltys,
         col = legend_cols,
         lwd = 2,
         cex = 0.7, bty = "n")

  # Additional plots if grid layout and CDE available
  if (plot_style == "grid" && has_cde) {
    # Plot 5: Forest plot style comparison
    plot_df <- subgroup_stats[1:2, ]  # Exclude overall

    plot(1:2, plot_df$AHR_po,
         xlim = c(0.5, 2.5),
         ylim = c(0.3, 3),
         pch = 19, cex = 2,
         xaxt = "n",
         xlab = "",
         ylab = "Hazard Ratio",
         main = "Subgroup Treatment Effects",
         log = "y")

    axis(1, at = 1:2, labels = plot_df$Subgroup)

    if (has_cde) {
      points(1:2 + 0.1, plot_df$CDE, pch = 17, cex = 2, col = "purple")
    }

    abline(h = 1, lty = 2, col = "gray")
    abline(h = hr_threshold, lty = 2, col = "red")

    # Add error bars (would need CI data)
    for (i in 1:2) {
      text(i, plot_df$AHR_po[i] * 1.2,
           paste("n=", plot_df$n[i], sep = ""),
           cex = 0.8)
    }

    # Plot 6: Distribution overlay
    hist(z_values, breaks = 30, probability = TRUE,
         main = "Biomarker Distribution with Subgroups",
         xlab = z_name,
         col = rgb(0.5, 0.5, 0.5, 0.5))

    # Overlay density for each subgroup
    dens_rec <- density(z_values[recommend_idx])
    dens_que <- density(z_values[question_idx])

    lines(dens_rec, col = "darkgreen", lwd = 2)
    lines(dens_que, col = "darkorange", lwd = 2)
    abline(v = optimal_cutpoint, lty = 2, col = "green", lwd = 2)

    legend("topright",
           legend = c("Recommend", "Questionable", "Cutpoint"),
           col = c("darkgreen", "darkorange", "green"),
           lty = c(1, 1, 2),
           lwd = 2,
           cex = 0.7, bty = "n")
  }

  # Save plots if requested
  if (save_plots) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- file.path(output_dir, paste0("ahr_cde_analysis_", timestamp, ".pdf"))
    dev.copy(pdf, filename, width = 12, height = 8)
    dev.off()
    if (verbose) cat("\nPlots saved to:", filename, "\n")
  }

  # ===========================================================================
  # Step 5: Compile and Return Results
  # ===========================================================================

  results <- list(
    cox_fit = cox_results,
    ahr_results = ahr_results,
    cde_results = cde_results,
    optimal_cutpoint = optimal_cutpoint,
    subgroup_stats = subgroup_stats,
    data = list(
      z_values = z_values,
      loghr_po = loghr_po,
      subgroup = df_sorted$subgroup
    )
  )

  class(results) <- c("cox_ahr_cde", "list")

  if (verbose) {
    cat("\n=====================================\n")
    cat("Analysis Complete\n")
    cat("=====================================\n")
    cat("Optimal cutpoint:", round(optimal_cutpoint, 2), "\n")
    cat("Recommend subgroup size:", sum(recommend_idx),
        "(", round(100 * sum(recommend_idx) / nrow(df_sorted), 1), "%)\n")
    cat("AHR (Recommend):", round(ahr_results$recommend, 3), "\n")
    if (has_cde) {
      cat("CDE (Recommend):", round(cde_results$recommend, 3), "\n")
    }
    cat("=====================================\n\n")
  }

  return(results)
}

#' Print method for cox_ahr_cde objects
#'
#' @param x A cox_ahr_cde object from cox_ahr_cde_analysis()
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the input object
#' @export
#' @method print cox_ahr_cde
print.cox_ahr_cde <- function(x, ...) {
  cat("\n===== Cox AHR/CDE Analysis Results =====\n\n")

  cat("Optimal Cutpoint:", round(x$optimal_cutpoint, 3), "\n\n")

  cat("Subgroup Statistics:\n")
  print(x$subgroup_stats, row.names = FALSE)

  cat("\n")

  if (!is.null(x$cde_results)) {
    cat("Method Comparison:\n")
    cat("  Overall:      AHR =", round(x$ahr_results$overall, 3),
        "  CDE =", round(x$cde_results$overall, 3), "\n")
    cat("  Recommend:    AHR =", round(x$ahr_results$recommend, 3),
        "  CDE =", round(x$cde_results$recommend, 3), "\n")
    cat("  Questionable: AHR =", round(x$ahr_results$questionable, 3),
        "  CDE =", round(x$cde_results$questionable, 3), "\n")
  } else {
    cat("Average Hazard Ratios (Potential Outcomes):\n")
    cat("  Overall:     ", round(x$ahr_results$overall, 3), "\n")
    cat("  Recommend:   ", round(x$ahr_results$recommend, 3), "\n")
    cat("  Questionable:", round(x$ahr_results$questionable, 3), "\n")
  }

  cat("\n=========================================\n")

  invisible(x)
}

#' Summary method for cox_ahr_cde objects
#'
#' @param object A cox_ahr_cde object from cox_ahr_cde_analysis()
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the input object
#' @export
#' @method summary cox_ahr_cde
summary.cox_ahr_cde <- function(object, ...) {
  cat("\n===== Cox AHR/CDE Analysis Summary =====\n\n")

  # Cox model summary
  cat("Cox Spline Model:\n")
  cat("  Primary effect (no interaction):", round(object$cox_fit$cox_primary, 3), "\n")
  cat("  Spline df:", length(object$cox_fit$spline_basis$knots) + 1, "\n\n")

  # Subgroup summary
  print(object)

  # Additional diagnostics
  cat("\nDiagnostics:\n")
  cat("  Total observations:", sum(object$subgroup_stats$n[1:2]), "\n")
  cat("  Total events:", sum(object$subgroup_stats$events[1:2]), "\n")
  cat("  Biomarker range:",
      round(range(object$data$z_values, na.rm = TRUE), 2), "\n")

  invisible(object)
}
