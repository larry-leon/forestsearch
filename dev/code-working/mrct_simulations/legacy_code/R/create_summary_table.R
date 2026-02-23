#' Create Enhanced Summary Table for Baseline Characteristics
#'
#' @param data Data frame containing the analysis data
#' @param treat_var Character. Name of treatment variable (must have 2 levels)
#' @param vars_continuous Character vector. Names of continuous variables
#' @param vars_categorical Character vector. Names of categorical variables
#' @param vars_binary Character vector. Names of binary (0/1) variables
#' @param var_labels Named list. Custom labels for variables (optional)
#' @param digits Integer. Number of decimal places for continuous variables
#' @param show_pvalue Logical. Include p-values column
#' @param show_smd Logical. Include SMD (effect size) column
#' @param show_missing Logical. Include missing data rows
#' @param table_title Character. Main title for the table
#' @param table_subtitle Character. Subtitle for the table (optional)
#' @param source_note Character. Source note at bottom (optional)
#' @param font_size Numeric. Base font size in pixels (default: 12)
#' @param header_font_size Numeric. Header font size in pixels (default: 14)
#' @param footnote_font_size Numeric. Footnote font size in pixels (default: 10)
#' @param use_alternating_rows Logical. Apply zebra striping (default: TRUE)
#' @param stripe_color Character. Color for alternating rows (default: "#f9f9f9")
#' @param indent_size Numeric. Indentation for sub-levels in pixels (default: 20)
#' @param highlight_pval Numeric. Highlight p-values below this threshold (default: 0.05)
#' @param highlight_smd Numeric. Highlight SMD values above this threshold (default: 0.2)
#' @param highlight_color Character. Color for highlighting (default: "#fff3cd")
#' @param compact_mode Logical. Reduce spacing for compact display (default: FALSE)
#' @param column_width_var Numeric. Width for Variable column in pixels (default: 200)
#' @param column_width_stats Numeric. Width for stat columns in pixels (default: 120)
#' @param show_column_borders Logical. Show vertical column borders (default: FALSE)
#' @param custom_css Character. Additional custom CSS styling (optional)
#'
#' @return A gt table object (or data frame if gt not available)
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' create_summary_table(
#'   data = trial_data,
#'   treat_var = "treatment",
#'   vars_continuous = c("age", "bmi"),
#'   vars_categorical = c("sex", "stage")
#' )
#'
#' # Customized appearance
#' create_summary_table(
#'   data = trial_data,
#'   treat_var = "treatment",
#'   vars_continuous = c("age", "bmi"),
#'   vars_categorical = c("sex", "stage"),
#'   font_size = 11,
#'   header_font_size = 13,
#'   use_alternating_rows = TRUE,
#'   highlight_pval = 0.05,
#'   compact_mode = TRUE
#' )
#' }

create_summary_table <- function(data,
                                 treat_var = "treat",
                                 vars_continuous = NULL,
                                 vars_categorical = NULL,
                                 vars_binary = NULL,
                                 var_labels = NULL,
                                 digits = 1,
                                 show_pvalue = TRUE,
                                 show_smd = TRUE,
                                 show_missing = TRUE,
                                 table_title = "Baseline Characteristics by Treatment Arm",
                                 table_subtitle = NULL,
                                 source_note = NULL,
                                 font_size = 12,
                                 header_font_size = 14,
                                 footnote_font_size = 10,
                                 use_alternating_rows = TRUE,
                                 stripe_color = "#f9f9f9",
                                 indent_size = 20,
                                 highlight_pval = 0.05,
                                 highlight_smd = 0.2,
                                 highlight_color = "#fff3cd",
                                 compact_mode = FALSE,
                                 column_width_var = 200,
                                 column_width_stats = 120,
                                 show_column_borders = FALSE,
                                 custom_css = NULL) {

  # Check if gt package is available
  use_gt <- requireNamespace("gt", quietly = TRUE)

  # Get treatment groups
  treat_levels <- sort(unique(data[[treat_var]]))
  n_treat <- length(treat_levels)

  if(n_treat != 2) {
    stop("Currently supports only 2 treatment arms. Found ", n_treat, " levels in ", treat_var)
  }

  # Split data by treatment
  data_0 <- data[data[[treat_var]] == treat_levels[1], ]
  data_1 <- data[data[[treat_var]] == treat_levels[2], ]

  # Initialize results table
  results <- data.frame(
    Variable = character(),
    Level = character(),
    Control = character(),
    Treatment = character(),
    P_value = numeric(),
    SMD = numeric(),
    is_header = logical(),  # New: to identify header rows
    indent_level = numeric(), # New: for indentation control
    stringsAsFactors = FALSE
  )

  # =====================================
  # Process continuous variables
  # =====================================
  if(!is.null(vars_continuous)) {
    for(var in vars_continuous) {

      # Check if variable exists
      if(!var %in% names(data)) {
        warning("Variable '", var, "' not found in data")
        next
      }

      # Calculate statistics
      mean_0 <- mean(data_0[[var]], na.rm = TRUE)
      sd_0 <- stats::sd(data_0[[var]], na.rm = TRUE)
      n_0 <- sum(!is.na(data_0[[var]]))
      miss_0 <- sum(is.na(data_0[[var]]))

      mean_1 <- mean(data_1[[var]], na.rm = TRUE)
      sd_1 <- stats::sd(data_1[[var]], na.rm = TRUE)
      n_1 <- sum(!is.na(data_1[[var]]))
      miss_1 <- sum(is.na(data_1[[var]]))

      # Format as mean (SD)
      val_0 <- sprintf(paste0("%.", digits, "f (", "%.", digits, "f)"),
                       mean_0, sd_0)
      val_1 <- sprintf(paste0("%.", digits, "f (", "%.", digits, "f)"),
                       mean_1, sd_1)

      # Handle NaN values
      if(is.na(mean_0) || is.na(sd_0)) val_0 <- "N/A"
      if(is.na(mean_1) || is.na(sd_1)) val_1 <- "N/A"

      # Calculate p-value (t-test)
      if(show_pvalue && n_0 > 1 && n_1 > 1) {
        p_val <- tryCatch(
          stats::t.test(data_0[[var]], data_1[[var]])$p.value,
          error = function(e) NA
        )
      } else {
        p_val <- NA
      }

      # Calculate SMD (Cohen's d)
      if(show_smd && n_0 > 1 && n_1 > 1 && sd_0 > 0 && sd_1 > 0) {
        pooled_sd <- sqrt(((n_0 - 1) * sd_0^2 + (n_1 - 1) * sd_1^2) /
                            (n_0 + n_1 - 2))
        if(pooled_sd > 0) {
          smd <- abs(mean_1 - mean_0) / pooled_sd
        } else {
          smd <- NA
        }
      } else {
        smd <- NA
      }

      # Variable label
      var_label <- ifelse(!is.null(var_labels) && var %in% names(var_labels),
                          var_labels[[var]], var)

      # Add main row
      results <- rbind(results, data.frame(
        Variable = var_label,
        Level = "Mean (SD)",
        Control = val_0,
        Treatment = val_1,
        P_value = p_val,
        SMD = smd,
        is_header = FALSE,
        indent_level = 0,
        stringsAsFactors = FALSE
      ))

      # Add missing data row if needed
      if(show_missing && (miss_0 > 0 || miss_1 > 0)) {
        results <- rbind(results, data.frame(
          Variable = "",
          Level = "Missing",
          Control = paste0(miss_0, " (",
                           sprintf("%.1f", 100 * miss_0/nrow(data_0)), "%)"),
          Treatment = paste0(miss_1, " (",
                             sprintf("%.1f", 100 * miss_1/nrow(data_1)), "%)"),
          P_value = NA,
          SMD = NA,
          is_header = FALSE,
          indent_level = 1,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # =====================================
  # Process categorical variables
  # =====================================
  if(!is.null(vars_categorical)) {
    for(var in vars_categorical) {

      # Check if variable exists
      if(!var %in% names(data)) {
        warning("Variable '", var, "' not found in data")
        next
      }

      # Get levels
      levels_var <- sort(unique(stats::na.omit(data[[var]])))

      # Variable label
      var_label <- ifelse(!is.null(var_labels) && var %in% names(var_labels),
                          var_labels[[var]], var)

      # Calculate p-value (chi-square)
      if(show_pvalue && length(levels_var) > 1) {
        tab <- table(data[[var]], data[[treat_var]])
        if(all(dim(tab) > 1) && min(tab) >= 5) {
          p_val <- tryCatch(
            stats::chisq.test(tab)$p.value,
            error = function(e) NA
          )
        } else if(all(dim(tab) > 1)) {
          # Use Fisher's exact test for small cell counts
          p_val <- tryCatch(
            stats::fisher.test(tab, simulate.p.value = TRUE)$p.value,
            error = function(e) NA
          )
        } else {
          p_val <- NA
        }
      } else {
        p_val <- NA
      }

      # Calculate Cramer's V for effect size
      if(show_smd && length(levels_var) > 1) {
        tab <- table(data[[var]], data[[treat_var]])
        if(all(dim(tab) > 1) && sum(tab) > 0) {
          chi2 <- tryCatch(
            stats::chisq.test(tab)$statistic,
            error = function(e) NA
          )
          if(!is.na(chi2)) {
            n <- sum(tab)
            k <- min(nrow(tab), ncol(tab))
            cramers_v <- sqrt(chi2 / (n * (k - 1)))
          } else {
            cramers_v <- NA
          }
        } else {
          cramers_v <- NA
        }
      } else {
        cramers_v <- NA
      }

      # Add header row
      results <- rbind(results, data.frame(
        Variable = var_label,
        Level = "",
        Control = "",
        Treatment = "",
        P_value = p_val,
        SMD = cramers_v,
        is_header = TRUE,
        indent_level = 0,
        stringsAsFactors = FALSE
      ))

      # Add level rows
      for(lev in levels_var) {
        n_0 <- sum(data_0[[var]] == lev, na.rm = TRUE)
        n_1 <- sum(data_1[[var]] == lev, na.rm = TRUE)

        total_0 <- sum(!is.na(data_0[[var]]))
        total_1 <- sum(!is.na(data_1[[var]]))

        pct_0 <- if(total_0 > 0) 100 * n_0 / total_0 else 0
        pct_1 <- if(total_1 > 0) 100 * n_1 / total_1 else 0

        val_0 <- sprintf("%d (%.1f%%)", n_0, pct_0)
        val_1 <- sprintf("%d (%.1f%%)", n_1, pct_1)

        results <- rbind(results, data.frame(
          Variable = "",
          Level = as.character(lev),
          Control = val_0,
          Treatment = val_1,
          P_value = NA,
          SMD = NA,
          is_header = FALSE,
          indent_level = 1,
          stringsAsFactors = FALSE
        ))
      }

      # Add missing row if needed
      miss_0 <- sum(is.na(data_0[[var]]))
      miss_1 <- sum(is.na(data_1[[var]]))

      if(show_missing && (miss_0 > 0 || miss_1 > 0)) {
        results <- rbind(results, data.frame(
          Variable = "",
          Level = "Missing",
          Control = paste0(miss_0, " (",
                           sprintf("%.1f", 100 * miss_0/nrow(data_0)), "%)"),
          Treatment = paste0(miss_1, " (",
                             sprintf("%.1f", 100 * miss_1/nrow(data_1)), "%)"),
          P_value = NA,
          SMD = NA,
          is_header = FALSE,
          indent_level = 1,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # =====================================
  # Process binary variables
  # =====================================
  if(!is.null(vars_binary)) {
    for(var in vars_binary) {

      # Check if variable exists
      if(!var %in% names(data)) {
        warning("Variable '", var, "' not found in data")
        next
      }

      # Calculate proportions
      n_0_yes <- sum(data_0[[var]] == 1, na.rm = TRUE)
      n_1_yes <- sum(data_1[[var]] == 1, na.rm = TRUE)

      n_0_total <- sum(!is.na(data_0[[var]]))
      n_1_total <- sum(!is.na(data_1[[var]]))

      if(n_0_total > 0 && n_1_total > 0) {
        pct_0 <- 100 * n_0_yes / n_0_total
        pct_1 <- 100 * n_1_yes / n_1_total

        val_0 <- sprintf("%d (%.1f%%)", n_0_yes, pct_0)
        val_1 <- sprintf("%d (%.1f%%)", n_1_yes, pct_1)

        # Calculate p-value (Fisher's exact or chi-square)
        if(show_pvalue) {
          tab <- table(factor(data[[var]], levels = c(0, 1)),
                       data[[treat_var]])
          if(all(dim(tab) == c(2, 2))) {
            if(min(tab) < 5) {
              p_val <- tryCatch(
                stats::fisher.test(tab)$p.value,
                error = function(e) NA
              )
            } else {
              p_val <- tryCatch(
                stats::chisq.test(tab)$p.value,
                error = function(e) NA
              )
            }
          } else {
            p_val <- NA
          }
        } else {
          p_val <- NA
        }

        # Calculate SMD for binary
        if(show_smd) {
          p0 <- n_0_yes / n_0_total
          p1 <- n_1_yes / n_1_total
          denom <- sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 2)
          if(denom > 0) {
            smd <- abs(p1 - p0) / denom
          } else {
            smd <- NA
          }
        } else {
          smd <- NA
        }
      } else {
        val_0 <- "0 (0.0%)"
        val_1 <- "0 (0.0%)"
        p_val <- NA
        smd <- NA
      }

      # Variable label
      var_label <- ifelse(!is.null(var_labels) && var %in% names(var_labels),
                          var_labels[[var]], var)

      results <- rbind(results, data.frame(
        Variable = var_label,
        Level = "",
        Control = val_0,
        Treatment = val_1,
        P_value = p_val,
        SMD = smd,
        is_header = FALSE,
        indent_level = 0,
        stringsAsFactors = FALSE
      ))

      # Add missing row if needed
      if(show_missing) {
        miss_0 <- sum(is.na(data_0[[var]]))
        miss_1 <- sum(is.na(data_1[[var]]))

        if(miss_0 > 0 || miss_1 > 0) {
          results <- rbind(results, data.frame(
            Variable = "",
            Level = "Missing",
            Control = paste0(miss_0, " (",
                             sprintf("%.1f", 100 * miss_0/nrow(data_0)), "%)"),
            Treatment = paste0(miss_1, " (",
                               sprintf("%.1f", 100 * miss_1/nrow(data_1)), "%)"),
            P_value = NA,
            SMD = NA,
            is_header = FALSE,
            indent_level = 1,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  # Format p-values
  results$P_value_formatted <- ifelse(
    is.na(results$P_value), "",
    ifelse(results$P_value < 0.001, "<0.001",
           sprintf("%.3f", results$P_value))
  )

  # Format SMD
  results$SMD_formatted <- ifelse(
    is.na(results$SMD), "",
    sprintf("%.2f", results$SMD)
  )

  # =====================================
  # Create gt table if package is available
  # =====================================
  if(use_gt) {
    # Select columns for display
    display_data <- results[, c("Variable", "Level", "Control", "Treatment",
                                "P_value_formatted", "SMD_formatted")]

    gt_table <- gt::gt(display_data)

    # ===== Column Labels =====
    gt_table <- gt::cols_label(
      gt_table,
      Variable = "Characteristic",
      Level = "",
      Control = paste0("Control\n(n=", nrow(data_0), ")"),
      Treatment = paste0("Treatment\n(n=", nrow(data_1), ")"),
      P_value_formatted = "P-value",
      SMD_formatted = "SMD"
    )

    # ===== Header =====
    gt_table <- gt::tab_header(
      gt_table,
      title = table_title,
      subtitle = table_subtitle
    )

    # ===== Font Sizes =====
    gt_table <- gt::tab_options(
      gt_table,
      table.font.size = gt::px(font_size),
      heading.title.font.size = gt::px(header_font_size),
      heading.subtitle.font.size = gt::px(font_size),
      column_labels.font.size = gt::px(font_size),
      footnotes.font.size = gt::px(footnote_font_size),
      source_notes.font.size = gt::px(footnote_font_size)
    )

    # ===== Spacing (Compact Mode) =====
    if(compact_mode) {
      gt_table <- gt::tab_options(
        gt_table,
        data_row.padding = gt::px(2),
        heading.padding = gt::px(4),
        column_labels.padding = gt::px(4)
      )
    } else {
      gt_table <- gt::tab_options(
        gt_table,
        data_row.padding = gt::px(5),
        heading.padding = gt::px(8),
        column_labels.padding = gt::px(6)
      )
    }

    # ===== Column Widths =====
    # Only attempt to set column widths if explicitly requested AND not NULL
    if (!is.null(column_width_var) || !is.null(column_width_stats)) {
      # Set defaults if one is NULL
      if (is.null(column_width_var)) column_width_var <- 200
      if (is.null(column_width_stats)) column_width_stats <- 120

      # Try to apply column widths
      width_applied <- FALSE

      # Attempt 1: Standard gt::px() syntax (most common)
      if (!width_applied) {
        gt_table <- tryCatch({
          result <- gt::cols_width(
            gt_table,
            Variable ~ gt::px(column_width_var),
            Level ~ gt::px(100),
            Control ~ gt::px(column_width_stats),
            Treatment ~ gt::px(column_width_stats),
            P_value_formatted ~ gt::px(80),
            SMD_formatted ~ gt::px(80)
          )
          width_applied <- TRUE
          result
        }, error = function(e) {
          gt_table  # Return unchanged if error
        })
      }

      # Attempt 2: String syntax (newer gt versions)
      if (!width_applied) {
        gt_table <- tryCatch({
          result <- gt::cols_width(
            gt_table,
            Variable ~ paste0(column_width_var, "px"),
            Level ~ "100px",
            Control ~ paste0(column_width_stats, "px"),
            Treatment ~ paste0(column_width_stats, "px"),
            P_value_formatted ~ "80px",
            SMD_formatted ~ "80px"
          )
          width_applied <- TRUE
          result
        }, error = function(e) {
          gt_table  # Return unchanged if error
        })
      }

    }
    # If column width parameters are NULL, skip entirely (no message)

    # ===== Column Alignment =====
    gt_table <- gt::cols_align(
      gt_table,
      align = "left",
      columns = c("Variable", "Level")
    )

    gt_table <- gt::cols_align(
      gt_table,
      align = "center",
      columns = c("Control", "Treatment", "P_value_formatted", "SMD_formatted")
    )

    # ===== Handle Missing Values =====
    if (utils::packageVersion("gt") >= "0.6.0") {
      gt_table <- gt::sub_missing(
        gt_table,
        columns = gt::everything(),
        missing_text = ""
      )
    } else {
      gt_table <- gt::fmt_missing(
        gt_table,
        columns = gt::everything(),
        missing_text = ""
      )
    }

    # ===== Style: Header Rows (Bold) =====
    gt_table <- gt::tab_style(
      gt_table,
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        columns = "Variable",
        rows = results$is_header | (results$Variable != "" & results$Level != "Missing")
      )
    )

    # ===== Style: Indentation =====
    if(any(results$indent_level > 0)) {
      for(level in unique(results$indent_level[results$indent_level > 0])) {
        gt_table <- gt::tab_style(
          gt_table,
          style = gt::cell_text(indent = gt::px(indent_size * level)),
          locations = gt::cells_body(
            columns = "Level",
            rows = results$indent_level == level
          )
        )
      }
    }

    # ===== Style: Alternating Rows =====
    if(use_alternating_rows) {
      even_rows <- seq(2, nrow(display_data), by = 2)
      if(length(even_rows) > 0) {
        gt_table <- gt::tab_style(
          gt_table,
          style = gt::cell_fill(color = stripe_color),
          locations = gt::cells_body(rows = even_rows)
        )
      }
    }

    # ===== Style: Column Headers Border =====
    gt_table <- gt::tab_style(
      gt_table,
      style = gt::cell_borders(
        sides = "bottom",
        color = "black",
        weight = gt::px(2)
      ),
      locations = gt::cells_column_labels()
    )

    # ===== Style: Column Borders (Optional) =====
    if(show_column_borders) {
      gt_table <- gt::tab_style(
        gt_table,
        style = gt::cell_borders(
          sides = "right",
          color = "#d3d3d3",
          weight = gt::px(1)
        ),
        locations = gt::cells_body(
          columns = c("Variable", "Level", "Control", "Treatment")
        )
      )
    }

    # ===== Style: Highlight Significant P-values =====
    if(show_pvalue && !is.null(highlight_pval)) {
      sig_rows <- which(!is.na(results$P_value) & results$P_value < highlight_pval)
      if(length(sig_rows) > 0) {
        gt_table <- gt::tab_style(
          gt_table,
          style = gt::cell_fill(color = highlight_color),
          locations = gt::cells_body(
            columns = "P_value_formatted",
            rows = sig_rows
          )
        )
      }
    }

    # ===== Style: Highlight Large SMD =====
    if(show_smd && !is.null(highlight_smd)) {
      large_smd_rows <- which(!is.na(results$SMD) & results$SMD > highlight_smd)
      if(length(large_smd_rows) > 0) {
        gt_table <- gt::tab_style(
          gt_table,
          style = gt::cell_fill(color = highlight_color),
          locations = gt::cells_body(
            columns = "SMD_formatted",
            rows = large_smd_rows
          )
        )
      }
    }

    # ===== Footnotes =====
    if(show_smd) {
      gt_table <- gt::tab_footnote(
        gt_table,
        footnote = "SMD = Standardized mean difference (Cohen's d for continuous, Cramer's V for categorical)",
        locations = gt::cells_column_labels(columns = "SMD_formatted")
      )
    }

    if(show_pvalue) {
      gt_table <- gt::tab_footnote(
        gt_table,
        footnote = "P-values: t-test for continuous, chi-square/Fisher's exact for categorical/binary variables",
        locations = gt::cells_column_labels(columns = "P_value_formatted")
      )
    }

    # ===== Source Note =====
    if(!is.null(source_note)) {
      gt_table <- gt::tab_source_note(
        gt_table,
        source_note = source_note
      )
    }

    # ===== Custom CSS =====
    if(!is.null(custom_css)) {
      gt_table <- gt::tab_options(
        gt_table,
        table.additional_css = custom_css
      )
    }

    # ===== Hide Columns if Needed =====
    if(!show_pvalue) {
      gt_table <- gt::cols_hide(gt_table, columns = "P_value_formatted")
    }

    if(!show_smd) {
      gt_table <- gt::cols_hide(gt_table, columns = "SMD_formatted")
    }

    return(gt_table)

  } else {
    # Return data frame if gt not available
    message("Note: gt package not available. Returning data frame instead of formatted table.")

    # Clean up column names
    names(results)[names(results) == "Control"] <- paste0("Control (n=", nrow(data_0), ")")
    names(results)[names(results) == "Treatment"] <- paste0("Treatment (n=", nrow(data_1), ")")

    # Select columns based on options
    cols_to_keep <- c("Variable", "Level",
                      paste0("Control (n=", nrow(data_0), ")"),
                      paste0("Treatment (n=", nrow(data_1), ")"))

    if(show_pvalue) {
      cols_to_keep <- c(cols_to_keep, "P_value_formatted")
    }

    if(show_smd) {
      cols_to_keep <- c(cols_to_keep, "SMD_formatted")
    }

    return(results[, cols_to_keep])
  }
}


# ============================================================================
# Helper Functions - Preset Configurations
# ============================================================================

#' Preset: Compact Table
#'
#' @param ... Arguments passed to create_summary_table()
#' @keywords internal

create_summary_table_compact <- function(...) {
  create_summary_table(
    ...,
    font_size = 10,
    header_font_size = 12,
    footnote_font_size = 9,
    compact_mode = TRUE,
    indent_size = 15
    # Note: No column_width_* for maximum compatibility
  )
}


#' Preset: Publication-Ready Table
#'
#' @param ... Arguments passed to create_summary_table()
#' @keywords internal

create_summary_table_publication <- function(...) {
  create_summary_table(
    ...,
    font_size = 11,
    header_font_size = 13,
    footnote_font_size = 9,
    use_alternating_rows = FALSE,
    show_column_borders = FALSE,
    compact_mode = FALSE,
    highlight_pval = NULL,  # No highlighting for publications
    highlight_smd = NULL
    # Note: No column_width_* for maximum compatibility
  )
}


#' Preset: Presentation Table (Large Fonts)
#'
#' @param ... Arguments passed to create_summary_table()
#' @keywords internal

create_summary_table_presentation <- function(...) {
  create_summary_table(
    ...,
    font_size = 14,
    header_font_size = 16,
    footnote_font_size = 12,
    use_alternating_rows = TRUE,
    compact_mode = FALSE
    # Note: No column_width_* for maximum compatibility
    # Users can add manually if needed for their gt version
  )
}


#' Preset: Minimal Table (No Highlighting, No Alternating)
#'
#' @param ... Arguments passed to create_summary_table()
#' @keywords internal

create_summary_table_minimal <- function(...) {
  create_summary_table(
    ...,
    use_alternating_rows = FALSE,
    highlight_pval = NULL,
    highlight_smd = NULL,
    show_column_borders = FALSE
    # Note: No column_width_* for maximum compatibility
  )
}
