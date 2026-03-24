#' Plot Subgroup Results Forest Plot
#'
#' Creates a publication-ready forest plot displaying identified subgroups with
#' hazard ratios, bias-corrected estimates, and cross-validation metrics.
#' This wrapper integrates ForestSearch results with the forestploter package.
#'
#' @description
#' Generates a comprehensive forest plot showing:
#' - ITT (Intent-to-Treat) population estimate
#' - Reference subgroups (e.g., by biomarker levels)
#' - Post-hoc identified subgroups with bias-corrected estimates
#' - Cross-validation agreement metrics as annotations
#'
#' @param fs_results List. A list containing ForestSearch analysis results with elements:
#'   \itemize{
#'     \item \code{fs.est}: ForestSearch estimation object from \code{\link{forestsearch}}
#'     \item \code{fs_bc}: Bootstrap bias-corrected results from
#'       \code{\link{forestsearch_bootstrap_dofuture}}
#'     \item \code{fs_kfold}: K-fold cross-validation results from
#'       \code{\link{forestsearch_Kfold}} or \code{\link{forestsearch_tenfold}} (optional)
#'     \item \code{fs_OOB}: Out-of-bag cross-validation results (optional, alternative to fs_kfold)
#'   }
#' @param df_analysis Data frame. The analysis dataset with outcome, event, and treatment
#'   variables.
#' @param subgroup_list List. Named list of subgroup definitions to include in the plot.
#'   Each element should be a list with:
#'   \itemize{
#'     \item \code{subset_expr}: Character string for subsetting (e.g., "BM> 1")
#'     \item \code{name}: Display name for the subgroup
#'     \item \code{type}: Either "reference" for pre-specified or "posthoc" for identified
#'   }
#' @param outcome.name Character. Name of the survival time variable.
#' @param event.name Character. Name of the event indicator variable.
#' @param treat.name Character. Name of the treatment variable.
#' @param E.name Character. Label for experimental arm (default: "Experimental").
#' @param C.name Character. Label for control arm (default: "Control").
#' @param est.scale Character. Estimate scale: "hr" or "1/hr" (default: "hr").
#' @param xlog Logical. If TRUE (default), the x-axis is plotted on a
#'   logarithmic scale. This is standard for hazard ratio forest plots
#'   where equal distances represent equal relative effects.
#' @param title_text Character. Plot title (default: NULL).
#' @param arrow_text Character vector of length 2. Arrow labels for forest plot
#'   (default: c("Favors Experimental", "Favors Control")).
#' @param footnote_text Character vector. Footnote text for the plot explaining CV metrics
#'   (default provides CV interpretation guidance; set to NULL to omit).
#' @param xlim Numeric vector of length 2. X-axis limits (default: c(0.25, 1.5)).
#' @param ticks_at Numeric vector. X-axis tick positions
#'   (default: c(0.25, 0.70, 1.0, 1.5)).
#' @param show_cv_metrics Logical. Whether to show cross-validation metrics
#'   (default: TRUE if fs_kfold or fs_OOB available).
#' @param cv_source Character. Source for CV metrics:
#'   "auto" (default, uses both if available, otherwise whichever is present),
#'   "kfold" (use fs_kfold only),
#'   "oob" (use fs_OOB only), or
#'   "both" (explicitly use both fs_kfold and fs_OOB, with K-fold first then OOB).
#' @param posthoc_colors Character vector. Colors for post-hoc subgroup rows
#'   (default: c("powderblue", "beige")).
#' @param reference_colors Character vector. Colors for reference subgroup rows
#'   (default: c("yellow", "powderblue")).
#' @param ci_column_spaces Integer. Number of spaces for the CI plot column width.
#'   More spaces = wider CI column (default: 20).
#' @param conf.level Numeric. Confidence level for intervals (default: 0.95 for 95% CI).
#'   Used to calculate the z-multiplier as qnorm(1 - (1 - conf.level)/2).
#' @param theme An fs_forest_theme object from \code{create_forest_theme()}.
#'   Use this to control plot sizing (fonts, row height, CI appearance,
#'   CV annotation font size). Default: NULL (uses default theme).
#'
#' @return A list containing:
#'   \describe{
#'     \item{plot}{The forestploter grob object (can be rendered with plot())}
#'     \item{data}{The data frame used for the forest plot}
#'     \item{row_types}{Character vector of row types for styling reference}
#'     \item{cv_metrics}{Cross-validation metrics text (if available)}
#'   }
#'
#' @details
#' ## ForestSearch Labeling Convention
#'
#' ForestSearch identifies subgroups based on hazard ratio thresholds:
#' \itemize{
#'   \item \code{sg.harm}: Contains the definition of the "harm" or "questionable"
#'     subgroup (H)
#'   \item \code{treat.recommend == 0}: Patient is IN the harm subgroup (H)
#'   \item \code{treat.recommend == 1}: Patient is in the COMPLEMENT subgroup
#'     (Hc, typically benefit)
#' }
#'
#' For \code{est.scale = "hr"} (searching for harm):
#' \itemize{
#'   \item H (treat.recommend=0): Subgroup defined by sg.harm with elevated HR
#'     (harm/questionable)
#'   \item Hc (treat.recommend=1): Complement of sg.harm (potential benefit)
#' }
#'
#' For \code{est.scale = "1/hr"} (searching for benefit):
#' \itemize{
#'   \item Roles are reversed: H becomes the benefit group
#' }
#'
#' @examples
#' \dontrun{
#' # Load ForestSearch results
#' load("fs_analysis_results.Rdata")  # Contains fs.est, fs_bc, fs_kfold
#'
#' # Define subgroups to display
#' subgroups <- list(
#'   bm_gt1 = list(
#'     subset_expr = "BM > 1",
#'     name = "BM > 1",
#'     type = "reference"
#'   ),
#'   bm_gt1_size_gt19 = list(
#'     subset_expr = "BM > 1 & tmrsize > 19",
#'     name = "BM > 1 & Tumor Size > 19",
#'     type = "posthoc"
#'   )
#' )
#'
#' # Create the forest plot with default theme
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs.est, fs_bc = fs_bc, fs_kfold = fs_kfold),
#'   df_analysis = df_itt,
#'   subgroup_list = subgroups,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "combo",
#'   E.name = "Experimental+CT",
#'   C.name = "CT"
#' )
#'
#' # Create with custom theme for larger plot
#' large_theme <- create_forest_theme(
#'   base_size = 14,
#'   row_padding = c(6, 4),
#'   cv_fontsize = 11
#' )
#'
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs.est, fs_bc = fs_bc, fs_kfold = fs_kfold),
#'   df_analysis = df_itt,
#'   subgroup_list = subgroups,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "combo",
#'   theme = large_theme
#' )
#'
#' # Display the plot
#' render_forestplot(result)
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for running the subgroup analysis
#' \code{\link{forestsearch_bootstrap_dofuture}} for bootstrap bias correction
#' \code{\link{forestsearch_Kfold}} for cross-validation
#' \code{\link{create_forest_theme}} for customizing plot appearance
#' \code{\link{render_forestplot}} for rendering the plot
#'
#' @importFrom survival coxph Surv
#' @importFrom grid gpar
#' @export

plot_subgroup_results_forestplot <- function(
    fs_results,
    df_analysis,
    subgroup_list = NULL,
    outcome.name,
    event.name,
    treat.name,
    E.name = "Experimental",
    C.name = "Control",
    est.scale = "hr",
    xlog = TRUE,
    title_text = NULL,
    arrow_text = c("Favors Experimental", "Favors Control"),
    footnote_text = c("Eg 80% of training found SG: 70% of B (+) also B in CV testing"),
    xlim = c(0.25, 1.5),
    ticks_at = c(0.25, 0.70, 1.0, 1.5),
    show_cv_metrics = TRUE,
    cv_source = c("auto", "kfold", "oob", "both"),
    posthoc_colors = c("powderblue", "beige"),
    reference_colors = c("yellow", "powderblue"),
    ci_column_spaces = 20,
    conf.level = 0.95,
    theme = NULL
) {

  # ==========================================================================
  # Input Validation
  # ==========================================================================

  if (!requireNamespace("forestploter", quietly = TRUE)) {
    stop("Package 'forestploter' is required but not installed.", call. = FALSE)
  }

  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required.")
  }

  # Validate fs_results structure
  fs.est <- fs_results$fs.est
  fs_bc <- fs_results$fs_bc
  fs_kfold <- fs_results$fs_kfold
  fs_OOB <- fs_results$fs_OOB

  # Match cv_source argument

  cv_source <- match.arg(cv_source)

 # Determine which CV sources to use based on cv_source parameter
  cv_data_kfold <- NULL
  cv_data_oob <- NULL

  if (cv_source == "auto" || cv_source == "both") {
    # When "auto" or "both", use all available sources
    if (!is.null(fs_kfold)) {
      cv_data_kfold <- fs_kfold
    }
    if (!is.null(fs_OOB)) {
      cv_data_oob <- fs_OOB
    }
  } else if (cv_source == "kfold") {
    if (!is.null(fs_kfold)) {
      cv_data_kfold <- fs_kfold
    }
  } else if (cv_source == "oob") {
    if (!is.null(fs_OOB)) {
      cv_data_oob <- fs_OOB
    }
  }

  # Determine if we have any CV data to show
  has_cv_data <- !is.null(cv_data_kfold) || !is.null(cv_data_oob)

  # Allow NULL fs.est if subgroup_list is provided
  if (is.null(fs.est) && is.null(subgroup_list)) {
    stop("Either fs_results$fs.est or subgroup_list must be provided")
  }

  # Calculate z-multiplier for confidence intervals
  z_alpha <- qnorm(1 - (1 - conf.level) / 2)

  # Create dynamic CI column label
  ci_label <- sprintf("HR (%d%% CI)", round(conf.level * 100))

  # ==========================================================================
  # Helper Functions
  # ==========================================================================

  # Create HR table row for a subgroup
  # Uses robust (sandwich) standard errors for consistency with cox_summary()
  create_hr_row <- function(dfa, sg_name, outcome.name, event.name, treat.name,
                            E.name, C.name) {

    # Check minimum sample size
    if (nrow(dfa) < 10) {
      warning(paste("Subgroup", sg_name, "has fewer than 10 observations"))
      return(NULL)
    }

    # Extract vectors for validation
    Y <- dfa[[outcome.name]]
    E <- dfa[[event.name]]
    Treat <- dfa[[treat.name]]

    # Check for sufficient events
    n_events <- sum(E)
    if (n_events < 2) {
      warning(paste("Subgroup", sg_name, "has fewer than 2 events"))
      return(NULL)
    }

    # Check treatment variation
    if (length(unique(Treat)) < 2) {
      warning(paste("Subgroup", sg_name, "has no variation in treatment"))
      return(NULL)
    }

    sf <- paste0("Surv(", outcome.name, ",", event.name, ") ~ ", treat.name)
    cox.formula <- as.formula(sf)

    fit <- tryCatch(
      survival::coxph(
        cox.formula,
        data = dfa,
        robust = TRUE,
        model = FALSE,
        x = FALSE,
        y = FALSE
      ),
      error = function(e) {
        warning(paste("Cox model failed for", sg_name, ":", e$message))
        NULL
      }
    )

    if (is.null(fit)) return(NULL)

    # Use conf.level for confidence intervals
    fit_summary <- summary(fit, conf.int = conf.level)
    conf_int <- fit_summary$conf.int

    # Handle edge case where conf.int might not exist
    if (is.null(conf_int) || nrow(conf_int) == 0) {
      warning(paste("Subgroup", sg_name, "has no confidence interval available"))
      return(NULL)
    }

    hr <- conf_int[1, c(1, 3, 4)]

    ntreat <- sum(Treat)
    ncontrol <- sum(1 - Treat)

    est <- hr[1]
    low <- hr[2]
    hi <- hr[3]
    se <- (hi - est) / z_alpha

    row <- data.frame(
      Subgroup = sg_name,
      stringsAsFactors = FALSE
    )
    row[[E.name]] <- ntreat
    row[[C.name]] <- ncontrol
    row$est <- est
    row$low <- low
    row$hi <- hi
    row$se <- se

    return(row)
  }

  # Create placeholder/header row
  create_header_row <- function(header_text, E.name, C.name) {
    row <- data.frame(
      Subgroup = header_text,
      stringsAsFactors = FALSE
    )
    row[[E.name]] <- ""
    row[[C.name]] <- ""
    row$est <- NA
    row$low <- NA
    row$hi <- NA
    row$se <- NA
    return(row)
  }

  # Create bias-corrected row from bootstrap results
  # REVISED: sg_name is used directly without appending "(bias-corrected)"
  create_bc_row <- function(bc_estimates, sg_name, ntreat, ncontrol, E.name, C.name) {

    if (is.null(bc_estimates)) {
      return(NULL)
    }

    row <- data.frame(
      Subgroup = paste0("  ", sg_name),
      stringsAsFactors = FALSE
    )
    row[[E.name]] <- ntreat
    row[[C.name]] <- ncontrol
    row$est <- bc_estimates$H2
    row$low <- bc_estimates$H2_lower
    row$hi <- bc_estimates$H2_upper
    row$se <- (bc_estimates$H2_upper - bc_estimates$H2) / z_alpha

    return(row)
  }

  # Generate CV sensitivity text
 # Works with either fs_kfold or fs_OOB (both have same structure)
  generate_sens_text <- function(cv_data, est.scale, cv_label = "CV") {
    if (is.null(cv_data)) return(NULL)

    cv <- cv_data$find_summary["Any"]
    if (est.scale == "hr") {
      Q <- cv_data$sens_summary["sens_H"]
      B <- cv_data$sens_summary["sens_Hc"]
    } else {
      Q <- cv_data$sens_summary["sens_Hc"]
      B <- cv_data$sens_summary["sens_H"]
    }

    cv_text <- paste0(cv_label, " found = ", round(100 * cv, 0), "%")
    aa <- paste0(round(100 * B, 0), "%,")
    bb <- paste0(round(100 * Q, 0), "%")
    sense_text <- paste("Agree(+,-) = ", aa, bb, collapse = ",")
    sg_text <- paste(cv_text, sense_text, sep = ", ")

    return(sg_text)
  }

  # ==========================================================================
  # Build the Forest Plot Data
  # ==========================================================================

  # Start with ITT row
  dt <- create_hr_row(df_analysis, "ITT", outcome.name, event.name, treat.name,
                      E.name, C.name)
  row_types <- c("itt")
  cv_texts <- list()
  cv_row_positions <- list()

  # Process subgroup_list if provided
  if (!is.null(subgroup_list)) {

    # -----------------------------------------------------------------
    # SECTION 1: Add reference subgroups
    # -----------------------------------------------------------------
    ref_sgs <- Filter(function(x) x$type == "reference", subgroup_list)
    if (length(ref_sgs) > 0) {

      for (sg in ref_sgs) {

        df_sg <- safe_subset(df_analysis, sg$subset_expr)


        if (!is.null(df_sg) && nrow(df_sg) > 10) {
          sg_row <- create_hr_row(df_sg, sg$name, outcome.name, event.name,
                                  treat.name, E.name, C.name)
          if (!is.null(sg_row)) {
            dt <- rbind(dt, sg_row)
            row_types <- c(row_types, "reference")
          }
        }
      }
    }

    # -----------------------------------------------------------------
    # SECTION 2: Add post-hoc subgroups from subgroup_list
    # -----------------------------------------------------------------
    posthoc_sgs <- Filter(function(x) x$type == "posthoc", subgroup_list)
    if (length(posthoc_sgs) > 0) {
      # Add separator header
      separator_row <- create_header_row("Post-hoc subgroups", E.name, C.name)
      dt <- rbind(dt, separator_row)
      row_types <- c(row_types, "separator")

      # Loop through and add each posthoc subgroup from subgroup_list
      for (sg in posthoc_sgs) {

        df_sg <- safe_subset(df_analysis, sg$subset_expr)

        if (!is.null(df_sg) && nrow(df_sg) > 10) {
          sg_row <- create_hr_row(df_sg, paste0("  ", sg$name), outcome.name,
                                  event.name, treat.name, E.name, C.name)
          if (!is.null(sg_row)) {
            dt <- rbind(dt, sg_row)
            row_types <- c(row_types, "posthoc_list")
          }
        }
      }

      # Add blank row after posthoc list subgroups
      blank_row <- create_header_row(" ", E.name, C.name)
      dt <- rbind(dt, blank_row)
      row_types <- c(row_types, "blank")
    }
  }

  # -----------------------------------------------------------------
  # SECTION 3: Add identified subgroups from ForestSearch
  # -----------------------------------------------------------------
  if (!is.null(fs.est) && !is.null(fs.est$df.est)) {

    df_fs <- fs.est$df.est

    # Handle est.scale for treatment assignment
    if (est.scale == "1/hr") {
      df_fs$treat <- 1 - df_fs[, treat.name]
    } else {
      df_fs$treat <- df_fs[, treat.name]
    }

    # Get subgroup definition from sg.harm
    # IMPORTANT: sg.harm defines the HARM subgroup (H), NOT benefit
    sg_harm <- fs.est$sg.harm
    if (!is.null(sg_harm)) {
      harm_label <- paste(sg_harm, collapse = " & ")
    } else {
      harm_label <- "Identified Subgroup"
    }

    # -----------------------------------------------------------------
    # Define benefit/questionable subgroups based on est.scale
    #
    # ForestSearch convention:
    #   - sg.harm defines H (the harm/questionable subgroup)
    #   - treat.recommend == 0 means patient is IN H (harm)
    #   - treat.recommend == 1 means patient is IN Hc (complement/benefit)
    #
    # For est.scale = "hr":
    #   - H (treat.recommend=0) = HARM subgroup = defined by sg.harm
    #   - Hc (treat.recommend=1) = BENEFIT subgroup = NOT(sg.harm)
    #
    # For est.scale = "1/hr":
    #   - Roles reversed: H becomes benefit, Hc becomes questionable
    # -----------------------------------------------------------------

    if (est.scale == "hr") {
      # For hr scale:
      #   treat.recommend==1 is Hc (complement = BENEFIT)
      #   treat.recommend==0 is H  (harm = QUESTIONABLE)
      df_benefit <- subset(df_fs, treat.recommend == 1)   # Hc
      df_question <- subset(df_fs, treat.recommend == 0)  # H

      # REVISED: Improved label format - "Not {def}: Benefit" instead of "NOT ({def}) (Benefit)"
      benefit_label <- paste0("Not ", harm_label)         # Hc = complement
      question_label <- harm_label                         # H = sg.harm

      # Bias-corrected estimates mapping
      bc_benefit_key <- "Hc_estimates"
      bc_question_key <- "H_estimates"

    } else {
      # For 1/hr scale: roles are reversed
      #   treat.recommend==0 is H (now BENEFIT)
      #   treat.recommend==1 is Hc (now QUESTIONABLE)
      df_benefit <- subset(df_fs, treat.recommend == 0)   # H is benefit
      df_question <- subset(df_fs, treat.recommend == 1)  # Hc is questionable

      benefit_label <- harm_label                          # H = sg.harm
      question_label <- paste0("Not ", harm_label)         # Hc = complement

      # Bias-corrected estimates mapping (reversed)
      bc_benefit_key <- "H_estimates"
      bc_question_key <- "Hc_estimates"
    }

    # -----------------------------------------------------------------
    # Add BENEFIT subgroup rows
    # REVISED: Label format changed to "Label: Benefit" instead of "Label (Benefit)"
    # -----------------------------------------------------------------
    benefit_header <- create_header_row(
      paste0("                   ", benefit_label, ":"), E.name, C.name)
    dt <- rbind(dt, benefit_header)
    row_types <- c(row_types, "posthoc_header")

    if (nrow(df_benefit) > 10) {
      benefit_row <- create_hr_row(df_benefit, "  Full-analysis",
                                   outcome.name, event.name, "treat",
                                   E.name, C.name)
      if (!is.null(benefit_row)) {
        dt <- rbind(dt, benefit_row)
        row_types <- c(row_types, "posthoc")
      }

      # Add bias-corrected if available
      # REVISED: Pass "Bias-corrected" as the label directly
      if (!is.null(fs_bc)) {
        bc_est <- fs_bc[[bc_benefit_key]]
        if (!is.null(bc_est)) {
          ntreat_bc <- nrow(subset(df_benefit, treat == 1))
          ncontrol_bc <- nrow(subset(df_benefit, treat == 0))
          bc_row <- create_bc_row(bc_est, "Bias-corrected", ntreat_bc,
                                  ncontrol_bc, E.name, C.name)
          if (!is.null(bc_row)) {
            dt <- rbind(dt, bc_row)
            row_types <- c(row_types, "posthoc_bc")
          }
        }
      }
    }

    # -----------------------------------------------------------------
    # Add QUESTIONABLE/HARM subgroup rows
    # REVISED: Label format changed to "Label: Questionable" instead of "Label (Questionable)"
    # -----------------------------------------------------------------
    question_header <- create_header_row(
      paste0(question_label, ":"), E.name, C.name)
    dt <- rbind(dt, question_header)
    row_types <- c(row_types, "posthoc_complement_header")

    if (nrow(df_question) > 10) {
      question_row <- create_hr_row(df_question, "  Full-analysis",
                                    outcome.name, event.name, "treat",
                                    E.name, C.name)
      if (!is.null(question_row)) {
        dt <- rbind(dt, question_row)
        row_types <- c(row_types, "posthoc_complement")
      }

      # Add bias-corrected if available
      # REVISED: Pass "Bias-corrected" as the label directly
      if (!is.null(fs_bc)) {
        bc_est <- fs_bc[[bc_question_key]]
        if (!is.null(bc_est)) {
          ntreat_bc <- nrow(subset(df_question, treat == 1))
          ncontrol_bc <- nrow(subset(df_question, treat == 0))
          bc_row <- create_bc_row(bc_est, "Bias-corrected", ntreat_bc,
                                  ncontrol_bc, E.name, C.name)
          if (!is.null(bc_row)) {
            dt <- rbind(dt, bc_row)
            row_types <- c(row_types, "posthoc_complement_bc")
          }
        }
      }

      # CV metrics for questionable subgroup (add blank rows for annotation space)
      # Supports both K-fold and OOB annotations when both are available
      if (has_cv_data && show_cv_metrics) {

        # First: K-fold annotation (if available)
        if (!is.null(cv_data_kfold)) {
          cv_text_kfold <- generate_sens_text(cv_data_kfold, est.scale, "K-fold")
          if (!is.null(cv_text_kfold)) {
            # Create blank row - text will be added via insert_text
            cv_row <- create_header_row("", E.name, C.name)
            dt <- rbind(dt, cv_row)
            row_types <- c(row_types, "cv_annotation")
            cv_texts[[length(cv_texts) + 1]] <- cv_text_kfold
            cv_row_positions[[length(cv_row_positions) + 1]] <- nrow(dt)
          }
        }

        # Second: OOB annotation (if available)
        if (!is.null(cv_data_oob)) {
          cv_text_oob <- generate_sens_text(cv_data_oob, est.scale, "OOB")
          if (!is.null(cv_text_oob)) {
            # Create blank row - text will be added via insert_text
            cv_row <- create_header_row("", E.name, C.name)
            dt <- rbind(dt, cv_row)
            row_types <- c(row_types, "cv_annotation")
            cv_texts[[length(cv_texts) + 1]] <- cv_text_oob
            cv_row_positions[[length(cv_row_positions) + 1]] <- nrow(dt)
          }
        }
      }
    }
  }

  # ==========================================================================
  # Create Color Scheme
  # ==========================================================================

  sg_colors <- sapply(row_types, function(type) {
    switch(type,
           "itt" = "yellow",
           "reference" = reference_colors[1],
           "separator" = "yellowgreen",
           "blank" = "white",
           "posthoc_list" = "lightyellow",
           "posthoc_header" = posthoc_colors[1],
           "posthoc" = posthoc_colors[1],
           "posthoc_bc" = posthoc_colors[1],
           "cv_annotation" = "white",
           "posthoc_complement_header" = posthoc_colors[2],
           "posthoc_complement" = posthoc_colors[2],
           "posthoc_complement_bc" = posthoc_colors[2],
           "white"
    )
  })

  # ==========================================================================
  # Create Forest Plot
  # ==========================================================================

  # Extract theme parameters or use defaults
  if (!is.null(theme) && inherits(theme, "fs_forest_theme")) {
    # Use custom theme parameters
    base_size <- theme$base_size
    row_padding <- theme$row_padding
    body_fontsize <- theme$body_fontsize
    header_fontsize <- theme$header_fontsize
    footnote_fontsize <- theme$footnote_fontsize
    footnote_col <- theme$footnote_col
    cv_fontsize <- theme$cv_fontsize
    cv_col <- theme$cv_col
    ci_pch <- theme$ci_pch
    ci_lwd <- theme$ci_lwd
    ci_Theight <- theme$ci_Theight
    ci_col <- theme$ci_col
    refline_lwd <- theme$refline_lwd
    refline_lty <- theme$refline_lty
    refline_col <- theme$refline_col
  } else {
    # Default values
    base_size <- 10
    row_padding <- c(4, 4)
    body_fontsize <- 10
    header_fontsize <- 11
    footnote_fontsize <- 9
    footnote_col <- "darkcyan"
    cv_fontsize <- 10
    cv_col <- "gray30"
    ci_pch <- 15
    ci_lwd <- 1.5
    ci_Theight <- 0.2
    ci_col <- "black"
    refline_lwd <- 1
    refline_lty <- "dashed"
    refline_col <- "gray30"
  }

  # Create forestploter theme with row-specific colors
  tm <- forestploter::forest_theme(
    base_size = base_size,
    core = list(
      fg_params = list(
        hjust = 1,
        x = 0.9,
        fontsize = body_fontsize
      ),
      bg_params = list(
        fill = sg_colors
      ),
      padding = grid::unit(row_padding, "mm")
    ),
    colhead = list(
      fg_params = list(
        hjust = 0.5,
        x = 0.5,
        fontsize = header_fontsize,
        fontface = "bold"
      )
    ),
    ci_pch = ci_pch,
    ci_col = ci_col,
    ci_lwd = ci_lwd,
    ci_Theight = ci_Theight,
    refline_lwd = refline_lwd,
    refline_lty = refline_lty,
    refline_col = refline_col,
    footnote_gp = grid::gpar(
      fontsize = footnote_fontsize,
      fontface = "italic",
      col = footnote_col
    )
  )

  # Add spacing column for CI plot (width controlled by ci_column_spaces)
  dt$` ` <- paste(rep(" ", ci_column_spaces), collapse = " ")

  # Create HR (xx% CI) text column with dynamic label
  dt[[ci_label]] <- ifelse(is.na(dt$se), "",
                           sprintf("%.2f (%.2f to %.2f)", dt$est, dt$low, dt$hi))

  # Suppress footnote if no CV data
  if (!has_cv_data) {
    footnote_text <- NULL
  }

  # Generate the forest plot
  p <- forestploter::forest(
    dt[, c("Subgroup", E.name, C.name, " ", ci_label)],
    title = title_text,
    est = dt$est,
    lower = dt$low,
    upper = dt$hi,
    sizes = 0.4,
    ci_column = 4,
    ref_line = 1,
    arrow_lab = arrow_text,
    xlim = xlim,
    ticks_at = ticks_at,
    xlog = xlog,
    footnote = footnote_text,
    theme = tm
  )

  # Add CV annotation text using insert_text (spans across first 3 columns)
  # Uses cv_fontsize from theme for proper sizing
  if (length(cv_row_positions) > 0) {
    for (i in seq_along(cv_texts)) {
      p <- forestploter::insert_text(
        p,
        text = paste0("  ", cv_texts[[i]]),
        row = cv_row_positions[[i]],
        col = 1:3,
        part = "body",
        just = "center",
        gp = grid::gpar(fontsize = cv_fontsize, fontface = "italic", col = cv_col)
      )
    }
  }

  # ==========================================================================
  # Return Results
  # ==========================================================================

  result <- list(
    plot = p,
    data = dt,
    row_types = row_types,
    cv_metrics = cv_texts
  )

  class(result) <- c("fs_forestplot", "list")

  return(result)
}


# ==============================================================================
# ADDITIONAL HELPER FUNCTIONS
# ==============================================================================

#' Create Subgroup Summary Data Frame for Forest Plot
#'
#' Creates a data frame suitable for forestploter from multiple subgroup analyses.
#' This is a more flexible alternative for complex subgroup configurations.
#'
#' @param df_analysis Data frame. The analysis dataset.
#' @param subgroups Named list of subgroup definitions.
#' @param outcome.name Character. Name of survival time variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @param E.name Character. Label for experimental arm.
#' @param C.name Character. Label for control arm.
#' @param fs_bc_list List. Named list of bootstrap results for each subgroup.
#' @param fs_kfold_list List. Named list of k-fold results for each subgroup.
#' @param conf.level Numeric. Confidence level for intervals (default: 0.95).
#'
#' @return Data frame with HR estimates for all subgroups.
#' @keywords internal

create_subgroup_summary_df <- function(
    df_analysis,
    subgroups,
    outcome.name,
    event.name,
    treat.name,
    E.name = "E",
    C.name = "C",
    fs_bc_list = NULL,
    fs_kfold_list = NULL,
    conf.level = 0.95
) {

  # Calculate z-multiplier for confidence intervals
  z_alpha <- qnorm(1 - (1 - conf.level) / 2)

  results <- list()

  # ITT
  results[["ITT"]] <- compute_sg_hr(
    df_analysis, "ITT", outcome.name, event.name, treat.name, E.name, C.name,
    z_alpha, conf.level
  )

  # Each subgroup
  for (sg_name in names(subgroups)) {
    sg <- subgroups[[sg_name]]

    df_sg <- safe_subset(df_analysis, sg$subset_expr)

    if (nrow(df_sg) > 10) {
      results[[sg_name]] <- compute_sg_hr(
        df_sg, sg$name, outcome.name, event.name, treat.name, E.name, C.name,
        z_alpha, conf.level
      )

      # Add bootstrap results if available
      if (!is.null(fs_bc_list[[sg_name]])) {
        bc <- fs_bc_list[[sg_name]]
        results[[paste0(sg_name, "_bc")]] <- data.frame(
          Subgroup = paste0(sg$name, " (bias-corrected)"),
          n_treat = results[[sg_name]]$n_treat,
          n_control = results[[sg_name]]$n_control,
          est = bc$H2,
          low = bc$H2_lower,
          hi = bc$H2_upper,
          se = (bc$H2_upper - bc$H2) / z_alpha
        )
      }
    }
  }

  # Combine into single data frame
  dt <- do.call(rbind, results)
  rownames(dt) <- NULL

  # Rename columns
  names(dt)[names(dt) == "n_treat"] <- E.name
  names(dt)[names(dt) == "n_control"] <- C.name

  return(dt)
}


#' Compute Hazard Ratio for a Single Subgroup
#'
#' Internal helper function to compute HR and CI for a subgroup.
#' Uses robust (sandwich) standard errors for consistency with cox_summary().
#'
#' @param df Data frame for the subgroup.
#' @param sg_name Character. Name of the subgroup.
#' @param outcome.name Character. Name of survival time variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @param E.name Character. Label for experimental arm.
#' @param C.name Character. Label for control arm.
#' @param z_alpha Numeric. Z-multiplier for CI (default: qnorm(0.975) for 95% CI).
#' @param conf.level Numeric. Confidence level for intervals (default: 0.95).
#'
#' @return Data frame with single row of HR estimates, or NULL if model fails.
#' @keywords internal

compute_sg_hr <- function(df, sg_name, outcome.name, event.name, treat.name,
                          E.name, C.name, z_alpha = qnorm(0.975),
                          conf.level = 0.95) {

 # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------

  Y <- df[[outcome.name]]
  E <- df[[event.name]]
  Treat <- df[[treat.name]]

  # Check for sufficient events
 n_events <- sum(E)
  if (n_events < 2) {
    warning("Subgroup '", sg_name, "': Fewer than 2 events; returning NA")
    ntreat <- sum(Treat)
    ncontrol <- sum(1 - Treat)
    return(data.frame(
      Subgroup = sg_name,
      n_treat = ntreat,
      n_control = ncontrol,
      est = NA_real_,
      low = NA_real_,
      hi = NA_real_,
      se = NA_real_
    ))
  }

  # Check treatment variation
  if (length(unique(Treat)) < 2) {
    warning("Subgroup '", sg_name, "': No variation in treatment; returning NA")
    ntreat <- sum(Treat)
    ncontrol <- sum(1 - Treat)
    return(data.frame(
      Subgroup = sg_name,
      n_treat = ntreat,
      n_control = ncontrol,
      est = NA_real_,
      low = NA_real_,
      hi = NA_real_,
      se = NA_real_
    ))
  }

  # -------------------------------------------------------------------------
  # Fit Cox model with robust SE and memory optimization
  # -------------------------------------------------------------------------

  sf <- paste0("Surv(", outcome.name, ",", event.name, ") ~ ", treat.name)
  cox.formula <- as.formula(sf)

  fit <- tryCatch({
    survival::coxph(
      cox.formula,
      data = df,
      robust = TRUE,
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
  }, error = function(e) {
    warning("Subgroup '", sg_name, "': Cox model failed: ", e$message)
    return(NULL)
  })

  # Handle fitting errors
  if (is.null(fit)) {
    ntreat <- sum(Treat)
    ncontrol <- sum(1 - Treat)
    return(data.frame(
      Subgroup = sg_name,
      n_treat = ntreat,
      n_control = ncontrol,
      est = NA_real_,
      low = NA_real_,
      hi = NA_real_,
      se = NA_real_
    ))
  }

  # -------------------------------------------------------------------------
  # Extract results
  # -------------------------------------------------------------------------

  # Use conf.level for confidence intervals
  fit_summary <- summary(fit, conf.int = conf.level)
  conf_int <- fit_summary$conf.int

  # Handle edge case where conf.int might not exist
  if (is.null(conf_int) || nrow(conf_int) == 0) {
    warning("Subgroup '", sg_name, "': No confidence interval available")
    ntreat <- sum(Treat)
    ncontrol <- sum(1 - Treat)
    return(data.frame(
      Subgroup = sg_name,
      n_treat = ntreat,
      n_control = ncontrol,
      est = NA_real_,
      low = NA_real_,
      hi = NA_real_,
      se = NA_real_
    ))
  }

  hr <- conf_int[1, c(1, 3, 4)]

  ntreat <- sum(Treat)
  ncontrol <- sum(1 - Treat)

  data.frame(
    Subgroup = sg_name,
    n_treat = ntreat,
    n_control = ncontrol,
    est = hr[1],
    low = hr[2],
    hi = hr[3],
    se = (hr[3] - hr[1]) / z_alpha
  )
}


#' Generate Cross-Validation Sensitivity Text
#'
#' Creates formatted text summarizing cross-validation agreement metrics.
#'
#' @param fs_kfold K-fold cross-validation results from forestsearch_Kfold.
#' @param est.scale Character. "hr" or "1/hr".
#'
#' @return Character string with formatted CV metrics.
#' @keywords internal

sens_text <- function(fs_kfold, est.scale = "hr") {
  if (is.null(fs_kfold)) return(NULL)

  cv <- fs_kfold$find_summary["Any"]
  if (est.scale == "hr") {
    Q <- fs_kfold$sens_summary["sens_H"]
    B <- fs_kfold$sens_summary["sens_Hc"]
  } else {
    Q <- fs_kfold$sens_summary["sens_Hc"]
    B <- fs_kfold$sens_summary["sens_H"]
  }

  cv_text <- paste0("CV found = ", round(100 * cv, 0), "%")
  aa <- paste0(round(100 * B, 0), "%,")
  bb <- paste0(round(100 * Q, 0), "%")
  sense_text <- paste("Agree(+,-) = ", aa, bb, collapse = ",")
  sg_text <- paste(cv_text, sense_text, sep = ", ")

  return(sg_text)
}
