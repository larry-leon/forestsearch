# =============================================================================
# plot_effect_distribution.R -- Three-panel effect distribution violin
# =============================================================================
#
# Companion to SGplot_estimates() in mrct_simulation.R.  Where SGplot_estimates
# is designed for the MRCT simulator (hr_itt / hr_train / hr_test / hr_sg
# column schema), plot_effect_distribution() operates on
# run_simulation_analysis() output (hr.itt / hr.H.hat / hr.Hc.hat) and
# supports both GLM (OR, RD, IRR) and survival (HR) endpoints.
# =============================================================================


#' Three-Panel Effect-Distribution Violin from Simulation Results
#'
#' Builds a three-panel violin plot (ITT / identified \eqn{\hat H} /
#' identified \eqn{\hat H^c}) of treatment-effect estimates across
#' simulation replicates produced by
#' \code{\link{run_simulation_analysis}()}.  Handles both GLM endpoints
#' (OR, RD, IRR) and survival (HR), optionally inverts effects to a
#' benefit scale (\code{subgroup_notation = "benefit"}), and overlays a
#' dotted reference line at the DGM-implied true effect in \eqn{\hat H}.
#'
#' The MRCT analogue is \code{\link{SGplot_estimates}}, which works on
#' the four-panel \code{mrct_region_sims()} schema (ITT / training /
#' testing / subgroup).  This function targets the three-panel
#' \code{run_simulation_analysis()} schema that does not involve a
#' train/test split.
#'
#' @param results \code{data.frame} returned by
#'   \code{\link{run_simulation_analysis}()} (optionally filtered).
#'   Must contain the columns \code{analysis}, \code{any.H},
#'   \code{hr.itt}, \code{hr.H.hat}, and \code{hr.Hc.hat}.  The column
#'   names end in \code{.hr} regardless of the underlying effect
#'   measure (OR / HR / RD / IRR); this follows the package's
#'   established naming convention.
#' @param analysis_method Character.  Method label in the
#'   \code{analysis} column to plot (e.g. \code{"FS"}, \code{"GRF"}).
#'   Default: \code{"FS"}.
#' @param effect_measure Character.  One of \code{"OR"}, \code{"HR"},
#'   \code{"RD"}, \code{"IRR"}, \code{"MD"}.  Used only for axis
#'   labelling and to decide the y-axis reference line (\eqn{1} for
#'   ratio measures, \eqn{0} for risk / mean differences).  \code{"MD"}
#'   (mean difference) is treated additively like \code{"RD"}.
#'   Default: \code{"OR"}.
#' @param subgroup_notation Character.  \code{"harm"} leaves effects
#'   on the original (switched) scale; \code{"benefit"} inverts
#'   multiplicative effects via \eqn{1/\hat\theta} and flips the sign
#'   of additive effects so that larger values favour the benefit
#'   subgroup.  Default: \code{"harm"}.
#' @param reference_effect Numeric or \code{NA}.  Optional overlay
#'   showing the DGM-implied true effect in \eqn{\hat H}
#'   (e.g. \code{1 / dgm$hazard_ratios$harm_subgroup} on a benefit
#'   scale).  Drawn as a dotted darkgreen line.  Default: \code{NA}.
#' @param title Character.  Plot title.  Default: auto-constructed.
#' @param subtitle Character or \code{NULL}.  Plot subtitle.
#'   Default: \code{NULL}.
#' @param panel_labels Named character vector with names \code{"itt"},
#'   \code{"H"}, \code{"Hc"} and the display labels for each violin.
#'   Default: \code{c(itt = "ITT", H = "Identified H-hat",
#'   Hc = "Identified Hc-hat")}.  When \code{subgroup_notation =
#'   "benefit"}, the default \code{Hc} label is adjusted to
#'   \code{"Identified Hc-hat (benefit)"}.
#' @param drop_undetected Logical.  Restrict the \eqn{\hat H} and
#'   \eqn{\hat H^c} panels to replicates where a subgroup was
#'   identified (\code{any.H == 1}).  Default: \code{TRUE}.  The ITT
#'   panel always uses all replicates.
#' @param trim_threshold Numeric or \code{NULL}.  Trimming is applied
#'   only if any of the three groups has a raw mean exceeding this
#'   absolute value -- i.e., extreme outliers that distort the violin
#'   or push annotated mean/SD into scientific notation.  Default:
#'   \code{1000}, matching \code{\link{build_estimation_table}} so the
#'   plot and the companion gt table auto-trim on the same condition.
#'   Set \code{NULL} to disable trimming entirely; set a smaller value
#'   (e.g., \code{100}) to be more aggressive on plots where reasonable
#'   OR / HR / IRR estimates sit well below 100 and even moderate
#'   outliers visibly distort the violin axis.
#' @param trim_fraction Numeric in (0, 0.5).  Fraction of observations
#'   to trim from each tail of each group when trimming triggers.
#'   Default: \code{0.01} (1\% from each tail; i.e., the central 98\%
#'   of values per group is kept).  Has no effect when
#'   \code{trim_threshold = NULL}.
#'
#' @return A \code{ggplot2} object.  Has \code{attr(p, "panel_data")}
#'   set to the long-format \code{data.table} used for plotting (for
#'   downstream summaries or diagnostic inspection).  When trimming is
#'   active, also has \code{attr(p, "trim_info")} containing per-group
#'   diagnostics (\code{n_total}, \code{n_trimmed}, \code{n_flagged},
#'   \code{raw_mean}, \code{raw_sd}, \code{trimmed_mean},
#'   \code{trimmed_sd}, \code{lower_bound}, \code{upper_bound}).
#'
#' @seealso \code{\link{run_simulation_analysis}} for the simulation
#'   pipeline, \code{\link{SGplot_estimates}} for the MRCT analogue,
#'   \code{\link{build_estimation_table}} for the tabular counterpart.
#'
#' @examples
#' \dontrun{
#' # Binary endpoint, benefit search, OR scale
#' true_or_benefit <- 1 / dgm_calibrated$hazard_ratios$harm_subgroup
#' p <- plot_effect_distribution(
#'   results_alt,
#'   analysis_method   = "FS",
#'   effect_measure    = "OR",
#'   subgroup_notation = "benefit",
#'   reference_effect  = true_or_benefit,
#'   title             = "ACTG175 Binary (HTE): OR Estimates Across Simulations",
#'   subtitle          = sprintf("n = 1000, 50 sims | truth OR(Q) = %.2f",
#'                               true_or_benefit)
#' )
#' print(p)
#'
#' # Survival endpoint, harm search, HR scale (no inversion)
#' p_surv <- plot_effect_distribution(
#'   results_alt,
#'   effect_measure    = "HR",
#'   subgroup_notation = "harm",
#'   reference_effect  = dgm_alt$hazard_ratios$harm_subgroup
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_hline
#'   scale_x_discrete scale_fill_brewer labs theme_minimal theme
#'   element_text
#' @export
plot_effect_distribution <- function(
    results,
    analysis_method   = "FS",
    effect_measure    = c("OR", "HR", "RD", "IRR", "MD"),
    subgroup_notation = c("harm", "benefit"),
    reference_effect  = NA_real_,
    title             = NULL,
    subtitle          = NULL,
    panel_labels      = NULL,
    drop_undetected   = TRUE,
    trim_threshold    = 1000,
    trim_fraction     = 0.01
) {

  if (!requireNamespace("ggplot2",    quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' required.", call. = FALSE)
  }

  effect_measure    <- match.arg(effect_measure)
  subgroup_notation <- match.arg(subgroup_notation)

  required_cols <- c("analysis", "any.H", "hr.itt", "hr.H.hat", "hr.Hc.hat")
  missing_cols  <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0L) {
    stop("Required column(s) missing from 'results': ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Slice to the requested method
  # ---------------------------------------------------------------------------
  rows <- results[results$analysis == analysis_method, , drop = FALSE]
  if (nrow(rows) == 0L) {
    stop(sprintf("No rows with analysis == '%s' in 'results'.",
                 analysis_method), call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Scale transform: invert multiplicative effects for "benefit";
  # flip sign of additive effects
  # ---------------------------------------------------------------------------
  is_ratio <- effect_measure %in% c("OR", "HR", "IRR")
  flip <- function(x) {
    if (subgroup_notation == "harm") return(x)
    if (is_ratio) {
      ifelse(is.finite(x) & x > 0, 1 / x, NA_real_)
    } else {
      -x
    }
  }

  any_H <- rows$any.H == 1L
  subset_H  <- if (drop_undetected) any_H else rep(TRUE, nrow(rows))
  subset_Hc <- subset_H

  panel_df <- data.table::rbindlist(list(
    data.table::data.table(est = flip(rows$hr.itt),                    group = "itt"),
    data.table::data.table(est = flip(rows$hr.H.hat[subset_H]),        group = "H"),
    data.table::data.table(est = flip(rows$hr.Hc.hat[subset_Hc]),      group = "Hc")
  ))

  # ---------------------------------------------------------------------------
  # Optional per-group trimming (activated only when trim_threshold is set
  # AND at least one group has |raw_mean| exceeding it).  Mirrors the pattern
  # in SGplot_estimates() but adds threshold-gating so backward-compatible
  # default (trim_threshold = NULL) leaves data untouched.
  # ---------------------------------------------------------------------------
  trim_info  <- NULL
  trim_active <- FALSE

  if (!is.null(trim_threshold)) {
    if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L ||
        trim_threshold <= 0) {
      stop("'trim_threshold' must be a positive numeric scalar or NULL.",
           call. = FALSE)
    }
    if (!is.numeric(trim_fraction) || length(trim_fraction) != 1L ||
        trim_fraction <= 0 || trim_fraction >= 0.5) {
      stop("'trim_fraction' must be in (0, 0.5).", call. = FALSE)
    }

    # Decide whether trimming triggers: any group's raw |mean| exceeds threshold
    raw_means <- vapply(c("itt", "H", "Hc"), function(g) {
      v <- panel_df$est[panel_df$group == g]
      v <- v[!is.na(v)]
      if (length(v) == 0L) NA_real_ else mean(v)
    }, numeric(1))

    trim_active <- any(!is.na(raw_means) &
                       abs(raw_means) > trim_threshold)

    if (trim_active) {
      panel_df[, trimmed := FALSE]
      trim_info     <- list()
      total_flagged <- 0L

      for (g in c("itt", "H", "Hc")) {
        idx        <- which(panel_df$group == g)
        vals       <- panel_df$est[idx]
        vals_clean <- vals[!is.na(vals)]

        if (length(vals_clean) < 5L) {
          trim_info[[g]] <- list(
            n_total = length(vals), n_trimmed = NA_integer_,
            n_flagged = 0L,
            raw_mean = if (length(vals_clean) > 0L) mean(vals_clean) else NA_real_,
            raw_sd   = if (length(vals_clean) > 1L) stats::sd(vals_clean) else NA_real_,
            trimmed_mean = NA_real_, trimmed_sd = NA_real_,
            lower_bound  = NA_real_, upper_bound = NA_real_
          )
          next
        }

        lo <- stats::quantile(vals_clean, trim_fraction,     na.rm = TRUE)
        hi <- stats::quantile(vals_clean, 1 - trim_fraction, na.rm = TRUE)

        is_extreme <- !is.na(panel_df$est[idx]) &
          (panel_df$est[idx] < lo | panel_df$est[idx] > hi)
        panel_df$trimmed[idx[is_extreme]] <- TRUE

        vals_trimmed   <- vals_clean[vals_clean >= lo & vals_clean <= hi]
        n_flagged      <- sum(is_extreme)
        total_flagged  <- total_flagged + n_flagged

        trim_info[[g]] <- list(
          n_total      = length(vals),
          n_trimmed    = length(vals_trimmed),
          n_flagged    = n_flagged,
          raw_mean     = mean(vals_clean),
          raw_sd       = stats::sd(vals_clean),
          trimmed_mean = mean(vals_trimmed),
          trimmed_sd   = stats::sd(vals_trimmed),
          lower_bound  = as.numeric(lo),
          upper_bound  = as.numeric(hi)
        )
      }

      # Auto-extend subtitle when caller did not supply one
      if (is.null(subtitle)) {
        subtitle <- sprintf(
          "Trimmed mean/SD shown (raw |mean| exceeded %s; %.0f%% symmetric trim per group; %d obs flagged)",
          format(trim_threshold, big.mark = ","),
          100 * trim_fraction, total_flagged
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Panel labels
  # ---------------------------------------------------------------------------
  default_labels <- c(
    itt = "ITT",
    H   = "Identified H-hat",
    Hc  = if (subgroup_notation == "benefit")
      "Identified Hc-hat (benefit)" else "Identified Hc-hat"
  )
  if (!is.null(panel_labels)) {
    for (nm in names(panel_labels)) default_labels[nm] <- panel_labels[nm]
  }
  panel_df$group <- factor(
    panel_df$group,
    levels = c("itt", "H", "Hc"),
    labels = default_labels[c("itt", "H", "Hc")]
  )

  # ---------------------------------------------------------------------------
  # x-axis annotations: mean (sd), N per panel
  # ---------------------------------------------------------------------------
  lvls <- levels(panel_df$group)
  # Map factor levels back to internal keys for trim_info lookup
  internal_key <- c("itt", "H", "Hc")
  names(internal_key) <- lvls

  annotated <- vapply(lvls, function(g) {
    key <- internal_key[g]

    if (trim_active && !is.null(trim_info[[key]]) &&
        !is.na(trim_info[[key]]$trimmed_mean)) {
      ti <- trim_info[[key]]
      sprintf("%s\nMean*: %.3f (%.3f),  N = %d  (%d flagged)",
              g, ti$trimmed_mean, ti$trimmed_sd, ti$n_trimmed,
              ti$n_flagged)
    } else {
      v <- panel_df$est[panel_df$group == g]
      v <- v[!is.na(v)]
      if (length(v) < 2L) g
      else sprintf("%s\nMean: %.3f (%.3f),  N = %d",
                   g, mean(v), stats::sd(v), length(v))
    }
  }, character(1))

  # ---------------------------------------------------------------------------
  # Axis label + null reference line
  # ---------------------------------------------------------------------------
  y_lab <- switch(
    effect_measure,
    "OR"  = if (subgroup_notation == "benefit")
             "Odds Ratio (benefit scale)" else "Odds Ratio",
    "HR"  = if (subgroup_notation == "benefit")
             "Hazard Ratio (benefit scale)" else "Hazard Ratio",
    "IRR" = if (subgroup_notation == "benefit")
             "Incidence Rate Ratio (benefit scale)" else "Incidence Rate Ratio",
    "RD"  = if (subgroup_notation == "benefit")
             "Risk Difference (benefit scale)" else "Risk Difference",
    "MD"  = if (subgroup_notation == "benefit")
             "Mean Difference (benefit scale)" else "Mean Difference"
  )
  null_ref <- if (is_ratio) 1 else 0

  if (is.null(title)) {
    title <- sprintf("Effect Estimates Across Simulations (%s)", analysis_method)
  }

  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  # When trimming is active, plot data excludes flagged rows so the violins
  # are not distorted by extremes; full panel_df is still returned as an
  # attribute for downstream inspection.
  plot_df <- if (trim_active) {
    panel_df[panel_df$trimmed == FALSE, ]
  } else {
    panel_df
  }

  p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$group, y = .data$est, fill = .data$group)
    ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.6) +
    ggplot2::geom_hline(yintercept = null_ref, colour = "red",
                        linetype = "dashed") +
    ggplot2::scale_x_discrete(labels = annotated) +
    ggplot2::scale_fill_brewer(palette = "Pastel1") +
    ggplot2::labs(x = "Analysis Population", y = y_lab,
                  title = title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(size = 10))

  if (!is.na(reference_effect)) {
    p <- p + ggplot2::geom_hline(yintercept = reference_effect,
                                 colour = "darkgreen", linetype = "dotted")
  }

  attr(p, "panel_data") <- panel_df
  if (trim_active) attr(p, "trim_info") <- trim_info
  p
}
