# =============================================================================
# plot_km_band_forestsearch.R - KM Band Plots for ForestSearch Subgroups
# =============================================================================
#
# Wrapper for plotKM.band_subgroups() that automatically configures subgroup
# indicators and labels from ForestSearch results.
#
# =============================================================================

#' Plot Kaplan-Meier Survival Difference Bands for ForestSearch Subgroups
#'
#' Creates Kaplan-Meier survival difference band plots comparing the identified
#' ForestSearch subgroup (sg.harm) and its complement against the ITT population.
#' This function wraps \code{plotKM.band_subgroups()} from the weightedsurv
#' package, automatically extracting subgroup definitions from ForestSearch
#' results.
#'
#' @param df Data frame. The analysis dataset containing all required variables
#'   including subgroup indicator columns.
#' @param fs.est A forestsearch object containing the identified subgroup,
#'   or \code{NULL} if using pre-defined subgroup indicators.
#' @param sg_cols Character vector. Names of columns in \code{df} containing
#'   subgroup indicators (0/1). These columns must already exist in \code{df}.
#'   If \code{NULL} and \code{fs.est} is provided, columns will be created
#'   automatically. Default: \code{NULL}
#' @param sg_labels Character vector. Subsetting expressions for each subgroup,
#'   corresponding to \code{sg_cols}. These are passed to
#'   \code{plotKM.band_subgroups()} which evaluates them as R expressions
#'   (e.g., \code{"age < 65"}, \code{"er <= 0"}, \code{"Qrecommend == 1"}).
#'   Must be same length as \code{sg_cols}.
#'   Default: \code{NULL} (auto-generated as \code{"colname == 1"}).
#' @param sg_colors Character vector. Colors for each subgroup curve,
#'   corresponding to \code{sg_cols}. Must be same length as \code{sg_cols}.
#'   Default: \code{NULL} (uses default color palette).
#' @param itt_color Character. Color for ITT population band.
#'   Default: \code{"azure3"}.
#' @param outcome.name Character. Name of time-to-event column.
#'   Default: \code{"tte"}.
#' @param event.name Character. Name of event indicator column.
#'   Default: \code{"event"}.
#' @param treat.name Character. Name of treatment column.
#'   Default: \code{"treat"}.
#' @param xlabel Character. X-axis label. Default: \code{"Time"}.
#' @param ylabel Character. Y-axis label. Default: \code{"Survival differences"}.
#' @param yseq_length Integer. Number of y-axis tick marks.
#'   Default: \code{5}.
#' @param draws_band Integer. Number of bootstrap draws for confidence band.
#'   Default: \code{1000}.
#' @param tau_add Numeric. Time horizon for the plot. If \code{NULL},
#'   auto-calculated from data. Default: \code{NULL}.
#' @param by_risk Numeric. Interval for risk table. Default: \code{6}.
#' @param risk_cex Numeric. Character expansion for risk table text.
#'   Default: \code{0.75}.
#' @param risk_delta Numeric. Vertical spacing for risk table.
#'   Default: \code{0.035}.
#' @param risk_pad Numeric. Padding for risk table. Default: \code{0.015}.
#' @param ymax_pad Numeric. Y-axis maximum padding. Default: \code{0.11}.
#' @param show_legend Logical. Whether to display the legend.
#'   Default: \code{TRUE}.
#' @param legend_pos Character. Legend position (e.g., "topleft", "bottomright").
#'   Default: \code{"topleft"}.
#' @param legend_cex Numeric. Character expansion for legend text.
#'   Default: \code{0.75}.
#' @param ref_subgroups Named list. Optional additional reference subgroups to
#'   include. Each element should be a list with:
#'   \describe{
#'     \item{subset_expr}{Character. R expression to define subgroup (e.g., "age < 65")}
#'     \item{label}{Character. Display label (optional, defaults to subset_expr)}
#'     \item{color}{Character. Color for the curve}
#'   }
#'   The function automatically creates indicator columns from the expressions.
#'   Default: \code{NULL}.
#' @param verbose Logical. Print diagnostic messages. Default: \code{FALSE}.
#'
#' @return Invisibly returns a list containing:
#'   \describe{
#'     \item{df}{The modified data frame with subgroup indicators}
#'     \item{sg_cols}{Character vector of subgroup column names used}
#'     \item{sg_labels}{Character vector of subgroup labels used}
#'     \item{sg_colors}{Character vector of colors used}
#'     \item{sg_harm_definition}{The subgroup definition extracted from fs.est}
#'     \item{ref_subgroups}{The reference subgroups list (if provided)}
#'   }
#'
#' @details
#' This function simplifies the workflow of creating KM survival difference
#' band plots for ForestSearch-identified subgroups. It can work in two modes:
#'
#' \strong{Mode 1: With ForestSearch result (\code{fs.est} provided)}
#' \enumerate{
#'   \item Extracts the subgroup definition from the ForestSearch result
#'   \item Creates binary indicator columns (Qrecommend, Brecommend) in \code{df}
#'   \item Generates appropriate labels from the subgroup definition
#'   \item Calls \code{plotKM.band_subgroups()} with configured parameters
#' }
#'
#' \strong{Mode 2: With pre-defined columns (\code{sg_cols} provided)}
#' \enumerate{
#'   \item Uses existing indicator columns in \code{df}
#'   \item Requires \code{sg_labels} and \code{sg_colors} to match \code{sg_cols}
#' }
#'
#' The sg.harm subgroup (Qrecommend) represents patients with questionable
#' treatment benefit (where \code{treat.recommend == 0} in ForestSearch output).
#' The complement (Brecommend) represents patients recommended for treatment.
#'
#' @section Subgroup Extraction:
#' When \code{fs.est} is provided, the subgroup definition is extracted from:
#' \itemize{
#'   \item \code{fs.est$grp.consistency$out_sg$sg.harm_label} - Human-readable labels
#'   \item \code{fs.est$sg.harm} - Technical factor names (fallback)
#'   \item \code{fs.est$df.est$treat.recommend} - Subgroup membership indicator
#' }
#'
#' @examples
#' \dontrun{
#' # Mode 1: Using ForestSearch result (auto-creates Qrecommend/Brecommend)
#' # This will use labels "Qrecommend == 1" and "Brecommend == 1"
#' plot_km_band_forestsearch(
#'   df = df.analysis,
#'   fs.est = fs_result,
#'   outcome.name = "os_months",
#'   event.name = "os_event",
#'   treat.name = "treatment"
#' )
#'
#' # Mode 1 with additional reference subgroups (auto-creates columns)
#' ref_sgs <- list(
#'   age_young = list(subset_expr = "age < 65", color = "brown"),
#'   age_old = list(subset_expr = "age >= 65", color = "grey")
#' )
#'
#' plot_km_band_forestsearch(
#'   df = df.analysis,
#'   fs.est = fs_result,
#'   ref_subgroups = ref_sgs,
#'   outcome.name = "os_months",
#'   event.name = "os_event",
#'   treat.name = "treatment",
#'   tau_add = 60,
#'   by_risk = 6
#' )
#'
#' # Mode 2: Using pre-defined subgroup columns with expression labels
#' # Note: sg_labels must be valid R expressions that evaluate against df
#' df.analysis$age_lt65 <- ifelse(df.analysis$age < 65, 1, 0)
#' df.analysis$age_ge65 <- ifelse(df.analysis$age >= 65, 1, 0)
#' df.analysis$Qrecommend <- ifelse(df.analysis$er <= 0, 1, 0)
#' df.analysis$Brecommend <- ifelse(df.analysis$er > 0, 1, 0)
#'
#' plot_km_band_forestsearch(
#'   df = df.analysis,
#'   sg_cols = c("age_lt65", "age_ge65", "Qrecommend", "Brecommend"),
#'   sg_labels = c("age < 65", "age >= 65", "er <= 0", "er > 0"),
#'   sg_colors = c("brown", "grey", "blue", "red"),
#'   outcome.name = "os_months",
#'   event.name = "os_event",
#'   treat.name = "treatment",
#'   tau_add = 60,
#'   by_risk = 6
#' )
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for running the subgroup analysis
#' \code{\link{plot_sg_weighted_km}} for weighted KM plots
#' \code{\link{plot_sg_results}} for comprehensive subgroup visualization
#'
#' @note
#' This function requires the \pkg{weightedsurv} package, which can be
#' installed from GitHub: \code{devtools::install_github("larry-leon/weightedsurv")}
#'
#' @importFrom graphics legend
#' @export
plot_km_band_forestsearch <- function(
    df,
    fs.est = NULL,
    sg_cols = NULL,
    sg_labels = NULL,
    sg_colors = NULL,
    itt_color = "azure3",
    outcome.name = "tte",
    event.name = "event",
    treat.name = "treat",
    xlabel = "Time",
    ylabel = "Survival differences",
    yseq_length = 5,
    draws_band = 1000,
    tau_add = NULL,
    by_risk = 6,
    risk_cex = 0.75,
    risk_delta = 0.035,
    risk_pad = 0.015,
    ymax_pad = 0.11,
    show_legend = TRUE,
    legend_pos = "topleft",
    legend_cex = 0.75,
    ref_subgroups = NULL,
    verbose = FALSE
) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION
  # ===========================================================================

  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }

  # Validate required columns exist
  required_cols <- c(outcome.name, event.name, treat.name)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ===========================================================================
  # SECTION 2: DETERMINE SUBGROUPS TO PLOT
  # ===========================================================================

  sg_harm_definition <- NULL

  if (!is.null(fs.est)) {
    # Extract subgroup definition from ForestSearch
    sg_harm_definition <- extract_fs_subgroup_definition(fs.est, verbose = verbose)

    # Create Qrecommend and Brecommend columns if not provided
    if (is.null(sg_cols)) {
      df <- create_fs_subgroup_indicators(
        df = df,
        fs.est = fs.est,
        col_names = c("Qrecommend", "Brecommend"),
        verbose = verbose
      )
      sg_cols <- c("Qrecommend", "Brecommend")

      # Auto-generate labels - try to get human-readable expressions
      # that plotKM.band_subgroups can evaluate
      if (is.null(sg_labels)) {
        sg_labels <- generate_readable_sg_labels(fs.est, verbose = verbose)
      }

      # Default colors for FS subgroups
      if (is.null(sg_colors)) {
        sg_colors <- c("blue", "red")
      }
    }
  }

  # ===========================================================================
  # SECTION 2b: PROCESS REFERENCE SUBGROUPS
 # ===========================================================================

  if (!is.null(ref_subgroups) && length(ref_subgroups) > 0) {
    ref_result <- create_reference_subgroup_columns(
      df = df,
      ref_subgroups = ref_subgroups,
      verbose = verbose
    )

    df <- ref_result$df
    ref_cols <- ref_result$cols
    ref_labels <- ref_result$labels
    ref_colors <- ref_result$colors

    # Append to existing subgroups
    if (!is.null(sg_cols)) {
      sg_cols <- c(sg_cols, ref_cols)
      sg_labels <- c(sg_labels, ref_labels)
      sg_colors <- c(sg_colors, ref_colors)
    } else {
      sg_cols <- ref_cols
      sg_labels <- ref_labels
      sg_colors <- ref_colors
    }
  }

  # Validate sg_cols
  if (is.null(sg_cols)) {
    stop("Either 'fs.est' or 'sg_cols' must be provided")
  }

  # Check all sg_cols exist in df
  missing_sg_cols <- setdiff(sg_cols, names(df))
  if (length(missing_sg_cols) > 0) {
    stop("Subgroup columns not found in data frame: ",
         paste(missing_sg_cols, collapse = ", "))
  }

  # Set default labels if not provided
  if (is.null(sg_labels)) {
    sg_labels <- sg_cols
  }

  # Set default colors if not provided
  if (is.null(sg_colors)) {
    default_palette <- c("blue", "red", "brown", "grey", "darkgreen",
                         "orange", "purple", "cyan", "magenta", "gold")
    sg_colors <- default_palette[seq_along(sg_cols)]
  }

  # Validate lengths match

  if (length(sg_labels) != length(sg_cols)) {
    stop("sg_labels must have same length as sg_cols")
  }
  if (length(sg_colors) != length(sg_cols)) {
    stop("sg_colors must have same length as sg_cols")
  }

  # ===========================================================================
  # SECTION 3: CALCULATE TAU IF NOT PROVIDED
  # ===========================================================================

  if (is.null(tau_add)) {
    tau_add <- max(df[[outcome.name]], na.rm = TRUE)
    if (verbose) {
      message(sprintf("Auto-calculated tau_add: %.1f", tau_add))
    }
  }

  # ===========================================================================
  # SECTION 4: DIAGNOSTIC OUTPUT
  # ===========================================================================

  if (verbose) {
    message("Subgroups to plot:")
    for (i in seq_along(sg_cols)) {
      n_in_sg <- sum(df[[sg_cols[i]]] == 1, na.rm = TRUE)
      message(sprintf("  %s (%s): n = %d, color = %s",
                      sg_labels[i], sg_cols[i], n_in_sg, sg_colors[i]))
    }
  }

  # ===========================================================================
  # SECTION 5: CALL plotKM.band_subgroups
  # ===========================================================================

  # Check that weightedsurv is available
  if (!requireNamespace("weightedsurv", quietly = TRUE)) {
    stop("Package 'weightedsurv' is required. ",
         "Install with: devtools::install_github('larry-leon/weightedsurv')")
  }

  # The function expects subgroup columns to exist in df with specific names
  # and sg_labels to be the display labels
  # Note: We explicitly list all parameters instead of using ... to avoid
  # passing wrapper-specific arguments (like ref_subgroups) to the underlying function
  result <- weightedsurv::plotKM.band_subgroups(
    df = df,
    sg_labels = sg_labels,
    sg_colors = sg_colors,
    color = itt_color,
    xlabel = xlabel,
    ylabel = ylabel,
    yseq_length = yseq_length,
    draws.band = draws_band,
    tau_add = tau_add,
    by.risk = by_risk,
    risk_cex = risk_cex,
    risk_delta = risk_delta,
    risk.pad = risk_pad,
    ymax.pad = ymax_pad,
    tte.name = outcome.name,
    treat.name = treat.name,
    event.name = event.name
  )

  # ===========================================================================
  # SECTION 6: ADD LEGEND
  # ===========================================================================

  if (show_legend) {
    legend_labels <- c("ITT", sg_labels)
    legend_colors <- c("black", sg_colors)

    legend(
      legend_pos,
      legend = legend_labels,
      lwd = 2,
      col = legend_colors,
      bty = "n",
      cex = legend_cex
    )
  }

  # ===========================================================================
  # SECTION 7: RETURN INVISIBLY
  # ===========================================================================

  invisible(list(
    df = df,
    sg_cols = sg_cols,
    sg_labels = sg_labels,
    sg_colors = sg_colors,
    sg_harm_definition = sg_harm_definition,
    ref_subgroups = ref_subgroups,
    plot_result = result
  ))
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Extract Subgroup Definition from ForestSearch Object
#'
#' Internal helper to extract human-readable subgroup definition.
#'
#' @param fs.est A forestsearch object.
#' @param verbose Logical. Print diagnostic messages.
#'
#' @return Character string describing the subgroup definition.
#' @keywords internal
extract_fs_subgroup_definition <- function(fs.est, verbose = FALSE) {

  definition <- NULL

  # Try grp.consistency first (most readable)
  if (!is.null(fs.est$grp.consistency$out_sg$sg.harm_label)) {
    definition <- paste(fs.est$grp.consistency$out_sg$sg.harm_label,
                        collapse = " & ")
  } else if (!is.null(fs.est$sg.harm)) {
    # Fall back to sg.harm factor names
    definition <- paste(fs.est$sg.harm, collapse = " & ")
  } else if (!is.null(fs.est$subgroup_definition)) {
    # Check for direct subgroup_definition
    definition <- fs.est$subgroup_definition
  }

  if (verbose && !is.null(definition)) {
    message("Extracted subgroup definition: ", definition)
  }

  return(definition)
}


#' Generate Readable Subgroup Labels from ForestSearch Object
#'
#' Extracts human-readable subgroup labels that are also valid R expressions
#' for use with \code{plotKM.band_subgroups()}. Attempts to extract the actual
#' subgroup definition (e.g., "er <= 0") rather than column references.
#'
#' @param fs.est A forestsearch object.
#' @param verbose Logical. Print diagnostic messages.
#'
#' @return Character vector of length 2: c(harm_label, benefit_label)
#' @keywords internal
generate_readable_sg_labels <- function(fs.est, verbose = FALSE) {

 # Try to extract the actual subgroup definition expression
  harm_expr <- NULL
  benefit_expr <- NULL

 # Method 1: Try grp.consistency$out_sg$sg.harm_label (most readable)
  if (!is.null(fs.est$grp.consistency$out_sg$sg.harm_label)) {
    labels <- fs.est$grp.consistency$out_sg$sg.harm_label
    # These are typically like "er <= 0" - valid R expressions
    harm_expr <- paste(labels, collapse = " & ")

    # Generate complement by negating
    benefit_expr <- generate_complement_expression(labels)

    if (verbose) {
      message("Extracted labels from grp.consistency$out_sg$sg.harm_label")
      message(sprintf("  Harm: %s", harm_expr))
      message(sprintf("  Benefit: %s", benefit_expr))
    }
  }

  # Method 2: Try subgroup_definition
  if (is.null(harm_expr) && !is.null(fs.est$subgroup_definition)) {
    harm_expr <- fs.est$subgroup_definition
    benefit_expr <- paste0("!(", harm_expr, ")")

    if (verbose) {
      message("Extracted labels from subgroup_definition")
    }
  }

  # Method 3: Try to parse sg.harm factor names into expressions
  if (is.null(harm_expr) && !is.null(fs.est$sg.harm)) {
    harm_expr <- parse_sg_harm_to_expression(fs.est$sg.harm, fs.est)
    if (!is.null(harm_expr)) {
      benefit_expr <- generate_complement_expression(harm_expr)
    }

    if (verbose && !is.null(harm_expr)) {
      message("Parsed labels from sg.harm factors")
    }
  }

  # Fallback: use column-based expressions
  if (is.null(harm_expr)) {
    harm_expr <- "Qrecommend == 1"
    benefit_expr <- "Brecommend == 1"

    if (verbose) {
      message("Using fallback column-based labels")
    }
  }

  return(c(harm_expr, benefit_expr))
}


#' Generate Complement Expression
#'
#' Creates the logical complement of a subgroup expression.
#' Handles common patterns like "var <= x" -> "var > x".
#'
#' @param expr Character vector of expressions to negate.
#'
#' @return Character string with negated expression.
#' @keywords internal
generate_complement_expression <- function(expr) {
  if (length(expr) == 0) return("Brecommend == 1")

  # If single expression, try smart negation
 if (length(expr) == 1) {
    e <- expr[1]

    # Pattern matching for common comparison operators
    # "var <= x" -> "var > x"
    if (grepl("<=", e)) {
      return(gsub("<=", ">", e))
    }
    # "var < x" -> "var >= x"
    if (grepl("<", e) && !grepl("<=", e)) {
      return(gsub("<", ">=", e))
    }
    # "var >= x" -> "var < x"
    if (grepl(">=", e)) {
      return(gsub(">=", "<", e))
    }
    # "var > x" -> "var <= x"
    if (grepl(">", e) && !grepl(">=", e)) {
      return(gsub(">", "<=", e))
    }
    # "var == x" -> "var != x"
    if (grepl("==", e)) {
      return(gsub("==", "!=", e))
    }
    # "var != x" -> "var == x"
    if (grepl("!=", e)) {
      return(gsub("!=", "==", e))
    }

    # Fallback: wrap in negation
    return(paste0("!(", e, ")"))
  }

  # Multiple expressions: negate the conjunction
  # "a & b" -> "!(a & b)" or could do "!a | !b" but former is cleaner
  combined <- paste(expr, collapse = " & ")
  return(paste0("!(", combined, ")"))
}


#' Parse sg.harm Factor Names to Expression
#'
#' Converts ForestSearch factor names (e.g., "er.0", "grade3.1") into
#' human-readable R expressions (e.g., "er <= 0", "grade3 == 1").
#'
#' @param sg_harm Character vector of factor names from fs.est$sg.harm.
#' @param fs.est ForestSearch object (for accessing confs_labels if available).
#'
#' @return Character string expression or NULL if parsing fails.
#' @keywords internal
parse_sg_harm_to_expression <- function(sg_harm, fs.est = NULL) {
  if (is.null(sg_harm) || length(sg_harm) == 0) return(NULL)

  # Check if we have confs_labels mapping
  confs_labels <- NULL
  if (!is.null(fs.est$parameters$confs_labels)) {
    confs_labels <- fs.est$parameters$confs_labels
  }

  expressions <- character(length(sg_harm))

  for (i in seq_along(sg_harm)) {
    factor_name <- sg_harm[i]

    # Try to use confs_labels if available
    if (!is.null(confs_labels) && factor_name %in% names(confs_labels)) {
      expressions[i] <- confs_labels[factor_name]
    } else {
      # Parse factor name pattern like "varname.level" -> "varname == level"
      # or use as-is if it looks like an expression already
      if (grepl("[<>=!]", factor_name)) {
        # Already an expression
        expressions[i] <- factor_name
      } else if (grepl("\\.", factor_name)) {
        # Pattern like "er.0" - convert to "er == 0" or similar
        parts <- strsplit(factor_name, "\\.")[[1]]
        if (length(parts) >= 2) {
          var_name <- parts[1]
          level <- parts[2]
          # This is a simplified conversion - may need adjustment
          expressions[i] <- paste0(var_name, " == ", level)
        } else {
          expressions[i] <- paste0(factor_name, " == 1")
        }
      } else {
        # Use as column indicator
        expressions[i] <- paste0(factor_name, " == 1")
      }
    }
  }

  if (length(expressions) == 1) {
    return(expressions[1])
  } else {
    return(paste(expressions, collapse = " & "))
  }
}


#' Create Reference Subgroup Indicator Columns
#'
#' Creates indicator columns (0/1) in the data frame for each reference
#' subgroup based on the provided subset expressions.
#'
#' @param df Data frame to modify.
#' @param ref_subgroups Named list of reference subgroup definitions.
#'   Each element should have \code{subset_expr} and optionally
#'   \code{label} and \code{color}.
#' @param verbose Logical. Print diagnostic messages.
#'
#' @return List with modified df, cols, labels, and colors vectors.
#' @keywords internal
create_reference_subgroup_columns <- function(df, ref_subgroups, verbose = FALSE) {

  cols <- character(0)
  labels <- character(0)
  colors <- character(0)

  # Default color palette for reference subgroups
  default_ref_colors <- c("brown", "grey", "darkgreen", "orange", "purple",
                          "cyan", "magenta", "gold", "pink", "navy")
  color_idx <- 1

  for (sg_name in names(ref_subgroups)) {
    sg <- ref_subgroups[[sg_name]]

    # Get subset expression
    subset_expr <- sg$subset_expr
    if (is.null(subset_expr)) {
      warning("Reference subgroup '", sg_name, "' missing 'subset_expr', skipping")
      next
    }

    # Create indicator column
    col_name <- paste0("ref_", sg_name)

    indicator <- as.integer(safe_eval_expr(df, subset_expr))

    if (is.null(indicator)) next

    df[[col_name]] <- indicator
    n_in_sg <- sum(indicator, na.rm = TRUE)

    # Get label (default to subset_expr)
    label <- if (!is.null(sg$label)) sg$label else subset_expr

    # Get color (default from palette)
    color <- if (!is.null(sg$color)) {
      sg$color
    } else {
      col <- default_ref_colors[color_idx]
      color_idx <- (color_idx %% length(default_ref_colors)) + 1
      col
    }

    cols <- c(cols, col_name)
    labels <- c(labels, label)
    colors <- c(colors, color)

    if (verbose) {
      message(sprintf("Created reference subgroup '%s': expr = '%s', n = %d, color = %s",
                      sg_name, subset_expr, n_in_sg, color))
    }
  }

  list(
    df = df,
    cols = cols,
    labels = labels,
    colors = colors
  )
}


#' Create Subgroup Indicator Columns from ForestSearch
#'
#' Internal helper to create Qrecommend and Brecommend indicator columns.
#'
#' @param df Data frame to modify.
#' @param fs.est A forestsearch object.
#' @param col_names Character vector of length 2. Names for the indicator
#'   columns: first for harm/questionable (treat.recommend == 0),
#'   second for benefit/recommend (treat.recommend == 1).
#'   Default: c("Qrecommend", "Brecommend")
#' @param verbose Logical. Print diagnostic messages.
#'
#' @return Modified data frame with indicator columns.
#' @keywords internal
create_fs_subgroup_indicators <- function(df, fs.est,
                                          col_names = c("Qrecommend", "Brecommend"),
                                          verbose = FALSE) {

  harm_col <- col_names[1]
  benefit_col <- col_names[2]

  # Method 1: Use treat.recommend from df.est
  if (!is.null(fs.est$df.est) && "treat.recommend" %in% names(fs.est$df.est)) {
    df_est <- fs.est$df.est

    # Match by ID if available
    if ("id" %in% names(df) && "id" %in% names(df_est)) {
      df[[harm_col]] <- as.integer(
        df$id %in% df_est$id[df_est$treat.recommend == 0]
      )
      df[[benefit_col]] <- as.integer(
        df$id %in% df_est$id[df_est$treat.recommend == 1]
      )
    } else {
      # Assume same row order
      if (nrow(df) == nrow(df_est)) {
        df[[harm_col]] <- as.integer(df_est$treat.recommend == 0)
        df[[benefit_col]] <- as.integer(df_est$treat.recommend == 1)
      } else {
        warning("Cannot match df to fs.est$df.est: different row counts")
      }
    }

    if (verbose) {
      message(sprintf("Created indicators: %s = %d, %s = %d",
                      harm_col, sum(df[[harm_col]], na.rm = TRUE),
                      benefit_col, sum(df[[benefit_col]], na.rm = TRUE)))
    }

    return(df)
  }

  # Method 2: Parse sg.harm factor names
  if (!is.null(fs.est$sg.harm)) {
    sg_factors <- fs.est$sg.harm

    # Create indicator by checking all factors
    df[[harm_col]] <- rep(1L, nrow(df))
    for (factor_name in sg_factors) {
      if (factor_name %in% names(df)) {
        df[[harm_col]] <- df[[harm_col]] * as.integer(df[[factor_name]] == 1)
      } else {
        warning("Factor '", factor_name, "' not found in data frame")
      }
    }
    df[[benefit_col]] <- as.integer(df[[harm_col]] == 0)

    if (verbose) {
      message("Created indicators from sg.harm factors: ",
              paste(sg_factors, collapse = ", "))
      message(sprintf("  %s = %d, %s = %d",
                      harm_col, sum(df[[harm_col]], na.rm = TRUE),
                      benefit_col, sum(df[[benefit_col]], na.rm = TRUE)))
    }

    return(df)
  }

  # Method 3: Check for grp.consistency
  if (!is.null(fs.est$grp.consistency$out_sg$sg.harm.id)) {
    df[[harm_col]] <- fs.est$grp.consistency$out_sg$sg.harm.id
    df[[benefit_col]] <- as.integer(df[[harm_col]] == 0)
    return(df)
  }

  warning("Could not create subgroup indicators from ForestSearch object")
  return(df)
}


#' Quick Plot KM Bands from ForestSearch
#'
#' Convenience wrapper with sensible defaults for quick visualization.
#'
#' @param df Data frame with analysis data.
#' @param fs.est ForestSearch result object.
#' @param outcome.name Character. Time-to-event column name.
#' @param event.name Character. Event indicator column name.
#' @param treat.name Character. Treatment column name.
#' @param ... Additional arguments passed to \code{plot_km_band_forestsearch()},
#'   such as \code{ref_subgroups}, \code{tau_add}, \code{by_risk}, etc.
#'
#' @return Invisibly returns the plot result.
#' @export
#'
#' @examples
#' \dontrun{
#' # Quick plot with minimal configuration
#' quick_km_band_plot(df_analysis, fs_result, "os_months", "os_event", "treat")
#'
#' # With reference subgroups
#' ref_sgs <- list(
#'   age_young = list(subset_expr = "age < 65", color = "brown"),
#'   age_old = list(subset_expr = "age >= 65", color = "grey")
#' )
#' quick_km_band_plot(df_analysis, fs_result, "os_months", "os_event", "treat",
#'                    ref_subgroups = ref_sgs, tau_add = 60)
#' }
quick_km_band_plot <- function(df, fs.est, outcome.name, event.name,
                               treat.name, ...) {
  plot_km_band_forestsearch(
    df = df,
    fs.est = fs.est,
    outcome.name = outcome.name,
    event.name = event.name,
    treat.name = treat.name,
    verbose = FALSE,
    ...
  )
}
