# =============================================================================
# plot_sg_glm_outcomes.R
#
# Subgroup outcome visualization for GLM (binary and continuous) outcomes.
# Analogous to plot_sg_weighted_km() for survival endpoints.
#
# For binary outcomes: grouped bar chart of event rates by treatment arm.
# For continuous outcomes: point-and-errorbar chart of means by treatment arm.
# Both annotated with effect estimates (OR, RR, RD, IRR, IRD, MD) and
# optional bias-corrected estimates from the bootstrap.
# =============================================================================


#' Plot GLM Subgroup Outcomes
#'
#' Creates a publication-ready visualization of identified subgroup outcomes
#' for binary or continuous endpoints.  This is the GLM counterpart to
#' \code{\link{plot_sg_weighted_km}} for survival endpoints.
#'
#' For binary outcomes, displays grouped bar charts of event rates by
#' treatment arm within each subgroup.  For continuous outcomes, displays
#' point-and-errorbar charts of means (+/- SE) by treatment arm.  Both
#' are annotated with the effect estimate and optional bias-corrected
#' estimate from the bootstrap.
#'
#' @param fs.est A forestsearch object containing \code{df.est} with
#'   \code{treat.recommend} column.
#' @param fs_bc Optional. Bootstrap results from
#'   \code{\link{forestsearch_bootstrap_dofuture}} for bias-corrected
#'   annotations.
#' @param outcome.name Character. Name of outcome column.
#' @param treat.name Character. Name of treatment column.
#' @param offset.name Character or \code{NULL}. Name of follow-up time
#'   column for rate-based measures (IRR, IRD).
#' @param effect_measure Character or \code{NULL}. Effect measure.
#'   Auto-detected from \code{fs.est$effect_measure} if \code{NULL}.
#' @param outcome_type Character or \code{NULL}. One of \code{"binary"},
#'   \code{"continuous"}.  Auto-detected from \code{fs.est$outcome_type}
#'   if \code{NULL}.
#' @param sg0_name Character. Label for H subgroup (treat.recommend == 0).
#'   Default: \code{"Questionable"}.
#' @param sg1_name Character. Label for Hc subgroup (treat.recommend == 1).
#'   Default: \code{"Recommend"}.
#' @param E.name Character. Label for experimental/treatment arm.
#'   Default: \code{"Treatment"}.
#' @param C.name Character. Label for control arm.
#'   Default: \code{"Control"}.
#' @param bar_colors Character vector of length 2. Colors for control and
#'   treatment bars.  Default: \code{c("steelblue", "coral")}.
#' @param show_effect Logical. Annotate with effect estimate.  Default
#'   \code{TRUE}.
#' @param show_bc Logical. Annotate with bias-corrected estimate when
#'   \code{fs_bc} is provided.  Default \code{TRUE}.
#' @param show_itt Logical. Include an ITT (full population) panel.
#'   Default \code{TRUE}.
#' @param conf.level Numeric. Confidence level for intervals.
#'   Default \code{0.95}.
#' @param title Character or \code{NULL}. Plot title.
#' @param subtitle Character or \code{NULL}. Plot subtitle.
#' @param verbose Logical. Print diagnostic messages. Default \code{FALSE}.
#'
#' @return A list of class \code{"fs_binary_plot"} containing:
#'   \describe{
#'     \item{plot}{A ggplot2 object.}
#'     \item{data}{Data frame of per-arm summaries used for the plot.}
#'     \item{effects}{Data frame of effect estimates per subgroup.}
#'     \item{sg_definition}{Character. Subgroup definition label.}
#'     \item{figure_note}{Character. Auto-generated figure note.}
#'   }
#'
#' @examples
#' \dontrun{
#' fs <- forestsearch(df, outcome_type = "binary", effect_measure = "OR",
#'                    ...)
#' result <- plot_sg_glm_outcomes(fs, outcome.name = "status",
#'                                 treat.name = "hormon")
#' print(result)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_segment
#'   geom_text geom_hline labs theme theme_minimal element_text
#'   element_blank scale_fill_manual scale_colour_manual
#'   coord_cartesian facet_wrap position_dodge
#' @export
plot_sg_glm_outcomes <- function(
    fs.est,
    fs_bc = NULL,
    outcome.name,
    treat.name,
    offset.name = NULL,
    effect_measure = NULL,
    outcome_type = NULL,
    sg0_name = "Questionable",
    sg1_name = "Recommend",
    E.name = "Treatment",
    C.name = "Control",
    bar_colors = c("steelblue", "coral"),
    show_effect = TRUE,
    show_bc = TRUE,
    show_itt = TRUE,
    conf.level = 0.95,
    title = NULL,
    subtitle = NULL,
    verbose = FALSE
) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION AND AUTO-DETECTION
  # ===========================================================================

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  if (is.null(fs.est) || is.null(fs.est$df.est)) {
    stop("fs.est must contain a non-empty df.est data frame.", call. = FALSE)
  }

  df <- fs.est$df.est

  # Auto-detect outcome_type first -- reject survival before other checks
  if (is.null(outcome_type)) {
    outcome_type <- fs.est$outcome_type
  }
  if (is.null(outcome_type) || outcome_type == "survival") {
    stop("plot_sg_glm_outcomes is for binary/continuous outcomes. ",
         "Use plot_sg_weighted_km for survival.", call. = FALSE)
  }

  if (!"treat.recommend" %in% names(df)) {
    stop("df.est must contain a 'treat.recommend' column.", call. = FALSE)
  }

  # Auto-detect remaining parameters
  if (is.null(effect_measure)) {
    effect_measure <- fs.est$effect_measure
  }
  if (is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      binary = "RD", continuous = "MD")
  }

  if (is.null(offset.name) && !is.null(fs.est$args_call_all$offset.name)) {
    offset.name <- fs.est$args_call_all$offset.name
  }

  is_log_scale <- effect_measure %in% c("OR", "RR", "IRR")
  z_alpha <- qnorm(1 - (1 - conf.level) / 2)

  # Subgroup definition label
  sg_definition <- if (!is.null(fs.est$sg.harm)) {
    paste(fs.est$sg.harm, collapse = " & ")
  } else {
    "Identified subgroup"
  }

  # ===========================================================================
  # SECTION 2: BUILD PER-ARM SUMMARIES
  # ===========================================================================

  # Split data into subgroups
  df_H  <- df[df$treat.recommend == 0, ]   # H: questionable
  df_Hc <- df[df$treat.recommend == 1, ]   # Hc: recommend

  build_arm_summary <- function(dfa, sg_label) {
    Y <- dfa[[outcome.name]]
    Treat <- dfa[[treat.name]]

    if (outcome_type == "binary") {
      # Event rates per arm
      rate_c <- mean(Y[Treat == 0], na.rm = TRUE)
      rate_t <- mean(Y[Treat == 1], na.rm = TRUE)
      n_c <- sum(Treat == 0)
      n_t <- sum(Treat == 1)
      # Wilson CI for proportions
      se_c <- sqrt(rate_c * (1 - rate_c) / n_c)
      se_t <- sqrt(rate_t * (1 - rate_t) / n_t)
      y_label <- "Event Rate"
    } else {
      # Means per arm
      rate_c <- mean(Y[Treat == 0], na.rm = TRUE)
      rate_t <- mean(Y[Treat == 1], na.rm = TRUE)
      n_c <- sum(Treat == 0)
      n_t <- sum(Treat == 1)
      se_c <- sd(Y[Treat == 0], na.rm = TRUE) / sqrt(n_c)
      se_t <- sd(Y[Treat == 1], na.rm = TRUE) / sqrt(n_t)
      y_label <- "Mean"
    }

    data.frame(
      Subgroup = sg_label,
      Arm = c(C.name, E.name),
      y = c(rate_c, rate_t),
      se = c(se_c, se_t),
      lo = c(rate_c - z_alpha * se_c, rate_t - z_alpha * se_t),
      hi = c(rate_c + z_alpha * se_c, rate_t + z_alpha * se_t),
      n = c(n_c, n_t),
      stringsAsFactors = FALSE
    )
  }

  # Build summaries
  summ_list <- list(
    build_arm_summary(df_H, sg0_name),
    build_arm_summary(df_Hc, sg1_name)
  )
  if (show_itt) {
    summ_list <- c(list(build_arm_summary(df, "ITT")), summ_list)
  }
  summ_df <- do.call(rbind, summ_list)

  # Order subgroups for display
  sg_levels <- if (show_itt) c("ITT", sg0_name, sg1_name) else c(sg0_name, sg1_name)
  summ_df$Subgroup <- factor(summ_df$Subgroup, levels = sg_levels)
  summ_df$Arm <- factor(summ_df$Arm, levels = c(C.name, E.name))

  # ===========================================================================
  # SECTION 3: COMPUTE EFFECT ESTIMATES
  # ===========================================================================

  effects_df <- NULL

  if (show_effect) {
    fn <- make_effect_estimator(
      outcome_type = outcome_type,
      treat.name = treat.name,
      outcome.name = outcome.name,
      offset.name = offset.name,
      effect_measure = effect_measure,
      ps_adjust_method = "none"
    )

    compute_effect <- function(dfa, sg_label) {
      res <- tryCatch(fn(dfa), error = function(e) NULL)
      if (is.null(res) || is.na(res$estimate)) return(NULL)

      est <- res$estimate
      se  <- res$se

      if (is_log_scale) {
        est_disp <- exp(est)
        lo_disp  <- exp(est - z_alpha * se)
        hi_disp  <- exp(est + z_alpha * se)
        label <- sprintf("%s = %.2f (%.2f, %.2f)",
                         effect_measure, est_disp, lo_disp, hi_disp)
      } else {
        lo_disp <- est - z_alpha * se
        hi_disp <- est + z_alpha * se
        label <- sprintf("%s = %.3f (%.3f, %.3f)",
                         effect_measure, est, lo_disp, hi_disp)
      }

      data.frame(
        Subgroup = sg_label,
        effect_label = label,
        stringsAsFactors = FALSE
      )
    }

    eff_list <- list(
      compute_effect(df_H, sg0_name),
      compute_effect(df_Hc, sg1_name)
    )
    if (show_itt) {
      eff_list <- c(list(compute_effect(df, "ITT")), eff_list)
    }
    effects_df <- do.call(rbind, Filter(Negate(is.null), eff_list))

    if (!is.null(effects_df)) {
      effects_df$Subgroup <- factor(effects_df$Subgroup, levels = sg_levels)
    }
  }

  # ===========================================================================
  # SECTION 4: BIAS-CORRECTED ANNOTATIONS
  # ===========================================================================

  bc_df <- NULL

  if (show_bc && !is.null(fs_bc)) {
    H_bc  <- fs_bc$H_estimates
    Hc_bc <- fs_bc$Hc_estimates

    build_bc_label <- function(bc_est, sg_label) {
      if (is.null(bc_est)) return(NULL)
      hr_bc <- bc_est$H2
      hr_lo <- bc_est$H2_lower
      hr_hi <- bc_est$H2_upper
      if (is.na(hr_bc)) return(NULL)

      label <- sprintf("%s(bc) = %.2f (%.2f, %.2f)",
                       effect_measure, hr_bc, hr_lo, hr_hi)
      data.frame(Subgroup = sg_label, bc_label = label,
                 stringsAsFactors = FALSE)
    }

    bc_list <- list(
      build_bc_label(H_bc, sg0_name),
      build_bc_label(Hc_bc, sg1_name)
    )
    bc_df <- do.call(rbind, Filter(Negate(is.null), bc_list))

    if (!is.null(bc_df)) {
      bc_df$Subgroup <- factor(bc_df$Subgroup, levels = sg_levels)
    }
  }

  # ===========================================================================
  # SECTION 5: BUILD GGPLOT
  # ===========================================================================

  y_label <- if (outcome_type == "binary") "Event Rate" else "Mean Outcome"

  if (is.null(title)) {
    title <- paste0("Subgroup Outcomes: ", sg_definition)
  }
  if (is.null(subtitle)) {
    subtitle <- paste0(
      effect_measure, " | ",
      outcome_type, " outcome | N = ", nrow(df)
    )
  }

  dodge_width <- 0.7

  if (outcome_type == "binary") {
    # Grouped bar chart for event rates
    p <- ggplot2::ggplot(
      summ_df,
      ggplot2::aes(x = Subgroup, y = y, fill = Arm)
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = dodge_width),
        width = 0.6, colour = "grey30", linewidth = 0.3
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(lo, 0), ymax = pmin(hi, 1)),
        position = ggplot2::position_dodge(width = dodge_width),
        width = 0.15, linewidth = 0.5
      ) +
      # N labels above bars
      ggplot2::geom_text(
        ggplot2::aes(
          y = pmin(hi, 1) + 0.02,
          label = paste0("n=", n)
        ),
        position = ggplot2::position_dodge(width = dodge_width),
        size = 3, vjust = 0
      ) +
      ggplot2::coord_cartesian(ylim = c(0, min(1, max(summ_df$hi, na.rm = TRUE) + 0.10)))

  } else {
    # Point-and-errorbar chart for continuous means
    p <- ggplot2::ggplot(
      summ_df,
      ggplot2::aes(x = Subgroup, y = y, colour = Arm)
    ) +
      ggplot2::geom_point(
        position = ggplot2::position_dodge(width = dodge_width),
        size = 3
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lo, ymax = hi),
        position = ggplot2::position_dodge(width = dodge_width),
        width = 0.2, linewidth = 0.6
      ) +
      # N labels
      ggplot2::geom_text(
        ggplot2::aes(
          y = hi + (max(summ_df$hi) - min(summ_df$lo)) * 0.03,
          label = paste0("n=", n)
        ),
        position = ggplot2::position_dodge(width = dodge_width),
        size = 3, vjust = 0
      )
  }

  # Common aesthetics
  p <- p +
    ggplot2::scale_fill_manual(
      values = stats::setNames(bar_colors, c(C.name, E.name)),
      name = NULL
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(bar_colors, c(C.name, E.name)),
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey40"),
      axis.text.x = ggplot2::element_text(face = "bold", size = 11),
      legend.position = "top",
      panel.grid.major.x = ggplot2::element_blank()
    )

  # ===========================================================================
  # SECTION 6: ADD EFFECT ESTIMATE ANNOTATIONS
  # ===========================================================================

  if (!is.null(effects_df) && nrow(effects_df) > 0) {
    # Compute y position for annotation (below bars)
    if (outcome_type == "binary") {
      annot_y <- -0.06
    } else {
      y_range <- max(summ_df$hi, na.rm = TRUE) - min(summ_df$lo, na.rm = TRUE)
      annot_y <- min(summ_df$lo, na.rm = TRUE) - y_range * 0.12
    }

    p <- p + ggplot2::geom_text(
      data = effects_df,
      ggplot2::aes(x = Subgroup, y = annot_y, label = effect_label),
      inherit.aes = FALSE,
      size = 3.2, fontface = "italic", colour = "grey20"
    )

    # Add BC annotation below the effect estimate
    if (!is.null(bc_df) && nrow(bc_df) > 0) {
      bc_y <- if (outcome_type == "binary") {
        annot_y - 0.05
      } else {
        annot_y - y_range * 0.08
      }

      p <- p + ggplot2::geom_text(
        data = bc_df,
        ggplot2::aes(x = Subgroup, y = bc_y, label = bc_label),
        inherit.aes = FALSE,
        size = 3.0, fontface = "bold.italic", colour = "darkblue"
      )
    }

    # Extend y-axis to accommodate annotations
    if (outcome_type == "binary") {
      bottom <- if (!is.null(bc_df)) annot_y - 0.08 else annot_y - 0.04
      p <- p + ggplot2::coord_cartesian(
        ylim = c(bottom, min(1, max(summ_df$hi, na.rm = TRUE) + 0.10)),
        clip = "off"
      )
    }
  }

  # ===========================================================================
  # SECTION 7: FIGURE NOTE
  # ===========================================================================

  note_parts <- character(0)
  note_parts <- c(note_parts,
    paste0("Identified subgroup (H): ", sg_definition, "."))

  if (outcome_type == "binary") {
    note_parts <- c(note_parts,
      "Bars show event rates per arm; error bars show 95% CI.")
  } else {
    note_parts <- c(note_parts,
      "Points show arm means; error bars show +/- 1.96 SE.")
  }

  if (!is.null(bc_df) && nrow(bc_df) > 0) {
    note_parts <- c(note_parts,
      paste0(effect_measure, "(bc) = bootstrap bias-corrected estimate."))
  }

  figure_note <- paste(note_parts, collapse = " ")

  if (verbose) {
    cat("Figure note:", figure_note, "\n")
  }

  # ===========================================================================
  # SECTION 8: RETURN
  # ===========================================================================

  result <- list(
    plot = p,
    data = summ_df,
    effects = effects_df,
    bc = bc_df,
    sg_definition = sg_definition,
    sg0_name = sg0_name,
    sg1_name = sg1_name,
    outcome_type = outcome_type,
    effect_measure = effect_measure,
    figure_note = figure_note
  )

  class(result) <- c("fs_binary_plot", "list")
  return(result)
}


#' @rdname plot_sg_glm_outcomes
#' @param x An fs_binary_plot object.
#' @param ... Additional arguments (unused).
#' @export
print.fs_binary_plot <- function(x, ...) {
  print(x$plot)
  invisible(x)
}

#' @rdname plot_sg_glm_outcomes
#' @export
plot.fs_binary_plot <- function(x, ...) {
  print(x$plot)
  invisible(x)
}
