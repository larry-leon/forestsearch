# =============================================================================
# CROSS-VALIDATION SUMMARY TABLE FUNCTIONS
# =============================================================================


#' Create Summary Tables from forestsearch_KfoldOut Results
#'
#' Formats the detailed output from \code{\link{forestsearch_KfoldOut}(outall=TRUE)}
#' into publication-ready gt tables. This includes ITT estimates, original subgroup
#' estimates, and K-fold subgroup estimates.
#'
#' @param kfold_out List. Result from \code{forestsearch_KfoldOut(res, outall = TRUE)}.
#'   Must contain \code{itt_tab}, \code{SG_tab_original}, \code{SG_tab_Kfold},
#'   and optionally \code{tab_all}.
#' @param title Character. Main title for combined table.
#'   Default: "Cross-Validation Summary".
#' @param subtitle Character. Subtitle for table. Default: NULL (auto-generated).
#' @param show_metrics Logical. Include agreement and finding metrics in output.
#'   Default: TRUE.
#' @param digits Integer. Decimal places for numeric formatting. Default: 3.
#' @param font_size Integer. Font size in pixels. Default: 12.
#' @param use_gt Logical. Return gt table if TRUE, data.frame if FALSE.
#'   Default: TRUE.
#'
#' @return If \code{use_gt = TRUE}, returns a list with gt table objects:
#'   \itemize{
#'     \item \code{combined_table}: Combined ITT and subgroup estimates
#'     \item \code{itt_table}: ITT estimates only
#'     \item \code{original_table}: Original full-data subgroup estimates
#'     \item \code{kfold_table}: K-fold subgroup estimates
#'     \item \code{metrics_table}: Agreement and finding metrics (if \code{show_metrics = TRUE})
#'   }
#'   If \code{use_gt = FALSE}, returns equivalent data.frames.
#'
#' @examples
#' \dontrun{
#' # Run K-fold CV
#' cv_results <- forestsearch_Kfold(fs.est = fs_result, Kfolds = 10)
#'
#' # Get detailed output
#' kfold_out <- forestsearch_KfoldOut(cv_results, outall = TRUE)
#'
#' # Create summary tables
#' cv_tables <- cv_summary_tables(kfold_out)
#' cv_tables$combined_table
#' cv_tables$metrics_table
#' }
#'
#' @seealso \code{\link{cv_metrics_tables}} for formatting \code{forestsearch_tenfold()} results
#'
#' @importFrom gt gt tab_header tab_spanner cols_label fmt_number tab_style
#'   cell_text cells_column_labels tab_footnote cells_body tab_options md px
#' @export
cv_summary_tables <- function(
    kfold_out,
    title = "Cross-Validation Summary",
    subtitle = NULL,
    show_metrics = TRUE,
    digits = 3,
    font_size = 12,
    use_gt = TRUE
) {

  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================

  if (!is.list(kfold_out)) {
    stop("kfold_out must be a list from forestsearch_KfoldOut(outall = TRUE)")
  }

  # Check for required elements (from outall = TRUE)
  required <- c("itt_tab", "SG_tab_original")
  missing <- setdiff(required, names(kfold_out))

  if (length(missing) > 0) {
    stop("kfold_out missing required elements: ", paste(missing, collapse = ", "),
         "\nDid you use forestsearch_KfoldOut(res, outall = TRUE)?")
  }

  # ===========================================================================
  # EXTRACT TABLES
  # ===========================================================================

  itt_tab <- kfold_out$itt_tab
  SG_tab_original <- kfold_out$SG_tab_original
  SG_tab_Kfold <- kfold_out$SG_tab_Kfold
  tab_all <- kfold_out$tab_all

  # Extract metrics if available
  sens_metrics <- kfold_out$sens_metrics_original
  find_metrics <- kfold_out$find_metrics

  # ===========================================================================
  # CREATE COMBINED TABLE (tab_all)
  # ===========================================================================

  if (!use_gt) {
    # Return data.frames
    result <- list(
      combined_table = tab_all,
      itt_table = itt_tab,
      original_table = SG_tab_original,
      kfold_table = SG_tab_Kfold
    )

    if (show_metrics && !is.null(sens_metrics)) {
      result$metrics <- list(
        sens_metrics = sens_metrics,
        find_metrics = find_metrics
      )
    }

    return(result)
  }

  # ===========================================================================
  # CREATE GT TABLES
  # ===========================================================================

  result <- list()

  # --- Combined Table ---
  if (!is.null(tab_all)) {
    # Auto-generate subtitle if not provided
    if (is.null(subtitle)) {
      subtitle <- "ITT and Subgroup Estimates: Full Analysis (FA) vs K-fold (Kf)"
    }

    combined_gt <- tab_all |>
      as.data.frame() |>
      gt::gt(rownames_to_stub = TRUE) |>
      gt::tab_header(
        title = gt::md(paste0("**", title, "**")),
        subtitle = subtitle
      ) |>
      gt::tab_stubhead(label = "Analysis") |>
      gt::tab_options(
        table.font.size = gt::px(font_size),
        heading.title.font.size = gt::px(font_size + 2),
        heading.subtitle.font.size = gt::px(font_size)
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_column_labels()
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_stub()
      )

    # Add row group labels
    combined_gt <- combined_gt |>
      gt::tab_footnote(
        footnote = "FA = Full Analysis (original data), Kf = K-fold CV estimate"
      )

    result$combined_table <- combined_gt
  }

  # --- ITT Table ---
  if (!is.null(itt_tab)) {
    itt_gt <- itt_tab |>
      as.data.frame() |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Intent-to-Treat Estimates**")
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size)
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_column_labels()
      )

    result$itt_table <- itt_gt
  }

  # --- Original Subgroup Table ---
  if (!is.null(SG_tab_original)) {
    original_gt <- SG_tab_original |>
      as.data.frame() |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Subgroup Estimates (Full Analysis)**"),
        subtitle = "Subgroups identified using all data"
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size)
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_column_labels()
      )

    result$original_table <- original_gt
  }

  # --- K-fold Subgroup Table ---
  if (!is.null(SG_tab_Kfold)) {
    kfold_gt <- SG_tab_Kfold |>
      as.data.frame() |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Subgroup Estimates (K-fold CV)**"),
        subtitle = "Subgroups identified via cross-validation"
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size)
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_column_labels()
      )

    result$kfold_table <- kfold_gt
  }

  # --- Metrics Table ---
  if (show_metrics && !is.null(sens_metrics) && !is.null(find_metrics)) {
    metrics_df <- data.frame(
      Category = c(
        rep("Agreement", 4),
        rep("Finding", 8)
      ),
      Metric = c(
        "Sensitivity (H)", "Sensitivity (Hc)", "PPV (H)", "PPV (Hc)",
        "Any Found", "Exact Match", "At Least 1", "Cov1 Any",
        "Cov2 Any", "Cov1 & Cov2", "Cov1 Exact", "Cov2 Exact"
      ),
      Value = c(
        sens_metrics["sens_H"], sens_metrics["sens_Hc"],
        sens_metrics["ppv_H"], sens_metrics["ppv_Hc"],
        find_metrics["Any"], find_metrics["Exact"],
        find_metrics["At least 1"], find_metrics["Cov1"],
        find_metrics["Cov2"], find_metrics["Cov 1 & 2"],
        find_metrics["Cov1 exact"], find_metrics["Cov2 exact"]
      ) * 100,
      stringsAsFactors = FALSE
    )
    rownames(metrics_df) <- NULL

    metrics_gt <- metrics_df |>
      gt::gt(groupname_col = "Category") |>
      gt::tab_header(
        title = gt::md("**Cross-Validation Metrics**"),
        subtitle = "Agreement and subgroup finding rates (%)"
      ) |>
      gt::fmt_number(columns = "Value", decimals = 1) |>
      gt::cols_label(
        Metric = "Metric",
        Value = "Value (%)"
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size)
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_column_labels()
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_row_groups()
      )

    result$metrics_table <- metrics_gt
  }

  return(result)
}


# =============================================================================
# CV METRICS TABLES (for forestsearch_tenfold / forestsearch_Kfold results)
# =============================================================================


#' Create Metrics Tables for Cross-Validation Results
#'
#' Formats the \code{find_summary} and \code{sens_summary} outputs from
#' \code{\link{forestsearch_tenfold}} or \code{\link{forestsearch_Kfold}}
#' into publication-ready gt tables.
#'
#' @param cv_result List. Result from \code{forestsearch_tenfold()} or
#'   \code{forestsearch_Kfold()}. Must contain \code{find_summary} and
#'   \code{sens_summary} elements.
#' @param sg_definition Character vector. Subgroup factor definitions for
#'   labeling (optional). If NULL, extracted from \code{cv_result$sg_analysis}.
#' @param title Character. Main title for combined table. Default: "Cross-Validation Metrics".
#' @param show_percentages Logical. Display metrics as percentages (0-100) instead
#'   of proportions (0-1). Default: TRUE.
#' @param digits Integer. Decimal places for formatting. Default: 1.
#' @param include_raw Logical. Include raw matrices (\code{sens_out}, \code{find_out})
#'   in the output for detailed analysis. Default: FALSE.
#' @param table_style Character. One of "combined", "separate", or "minimal".
#'   \itemize{
#'     \item "combined": Single table with both agreement and finding metrics
#'     \item "separate": Two separate gt tables
#'     \item "minimal": Compact single-row summary
#'   }
#'   Default: "combined".
#' @param use_gt Logical. Return gt table(s) if TRUE, data.frame(s) if FALSE.
#'   Default: TRUE.
#'
#' @return Depending on \code{table_style}:
#'   \itemize{
#'     \item "combined": A single gt table (or data.frame)
#'     \item "separate": A list with \code{agreement_table} and \code{finding_table}
#'     \item "minimal": A single-row gt table (or data.frame)
#'   }
#'   If \code{include_raw = TRUE}, also includes \code{sens_out} and \code{find_out}
#'   matrices in the returned list.
#'
#' @examples
#' \dontrun{
#' # After running forestsearch_tenfold
#' tenfold_results <- forestsearch_tenfold(
#'   fs.est = fs_result,
#'   sims = 100,
#'   Kfolds = 10
#' )
#'
#' # Create combined metrics table
#' cv_tables <- cv_metrics_tables(tenfold_results)
#' cv_tables
#'
#' # Create separate tables
#' cv_tables <- cv_metrics_tables(tenfold_results, table_style = "separate")
#' cv_tables$agreement_table
#' cv_tables$finding_table
#'
#' # Minimal one-row summary
#' cv_metrics_tables(tenfold_results, table_style = "minimal")
#' }
#'
#' @seealso \code{\link{cv_summary_tables}} for formatting \code{forestsearch_KfoldOut(outall=TRUE)} results
#'
#' @importFrom gt gt tab_header tab_spanner cols_label fmt_number tab_style
#'   cell_text cells_column_labels tab_footnote cells_body tab_options
#' @export
cv_metrics_tables <- function(
    cv_result,
    sg_definition = NULL,
    title = "Cross-Validation Metrics",
    show_percentages = TRUE,
    digits = 1,
    include_raw = FALSE,
    table_style = c("combined", "separate", "minimal"),
    use_gt = TRUE
) {

 # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================

  table_style <- match.arg(table_style)

  if (!is.list(cv_result)) {
    stop("cv_result must be a list from forestsearch_tenfold() or forestsearch_Kfold()")
  }

  # Check for required elements
  if (is.null(cv_result$sens_summary) || is.null(cv_result$find_summary)) {
    stop("cv_result must contain 'sens_summary' and 'find_summary' elements")
  }

  sens_summary <- cv_result$sens_summary
  find_summary <- cv_result$find_summary

  # Extract metadata
  sims <- if (!is.null(cv_result$sims)) cv_result$sims else 1
  Kfolds <- if (!is.null(cv_result$Kfolds)) cv_result$Kfolds else NA

  # Get subgroup definition
  if (is.null(sg_definition) && !is.null(cv_result$sg_analysis)) {
    sg_definition <- cv_result$sg_analysis
  }
  sg_label <- if (!is.null(sg_definition)) {
    paste(sg_definition, collapse = " & ")
  } else {
    "Identified Subgroup"
  }

  # ===========================================================================
  # PREPARE DATA
  # ===========================================================================

  # Multiplier for percentages
  mult <- if (show_percentages) 100 else 1
  suffix <- if (show_percentages) "%" else ""

  # Agreement metrics (sensitivity/PPV)
  agreement_df <- data.frame(
    Metric = c("Sensitivity (H)", "Sensitivity (Hc)",
               "PPV (H)", "PPV (Hc)"),
    Description = c(
      "Agreement rate for subgroup H",
      "Agreement rate for complement Hc",
      "Positive predictive value for H",
      "Positive predictive value for Hc"
    ),
    Value = c(
      sens_summary["sens_H"],
      sens_summary["sens_Hc"],
      sens_summary["ppv_H"],
      sens_summary["ppv_Hc"]
    ) * mult,
    stringsAsFactors = FALSE
  )
  rownames(agreement_df) <- NULL

  # Finding metrics
  finding_df <- data.frame(
    Metric = c("Any Found", "Exact Match", "At Least 1",
               "Cov1 Any", "Cov2 Any", "Cov1 & Cov2",
               "Cov1 Exact", "Cov2 Exact"),
    Description = c(
      "Any subgroup identified",
      "Exact match on all factors",
      "At least one factor matches",
      "First covariate found (any cut)",
      "Second covariate found (any cut)",
      "Both covariates found",
      "First covariate exact match",
      "Second covariate exact match"
    ),
    Value = c(
      find_summary["Any"],
      find_summary["Exact"],
      find_summary["At least 1"],
      find_summary["Cov1"],
      find_summary["Cov2"],
      find_summary["Cov 1 & 2"],
      find_summary["Cov1 exact"],
      find_summary["Cov2 exact"]
    ) * mult,
    stringsAsFactors = FALSE
  )
  rownames(finding_df) <- NULL

  # ===========================================================================
  # CREATE TABLES BASED ON STYLE
  # ===========================================================================

  if (table_style == "combined") {

    # Combine into single table
    combined_df <- data.frame(
      Category = c(rep("Agreement", 4), rep("Subgroup Finding", 8)),
      Metric = c(agreement_df$Metric, finding_df$Metric),
      Description = c(agreement_df$Description, finding_df$Description),
      Value = c(agreement_df$Value, finding_df$Value),
      stringsAsFactors = FALSE
    )

    if (!use_gt) {
      result <- combined_df
    } else {
      # Create gt table
      gt_table <- combined_df |>
        gt::gt(groupname_col = "Category") |>
        gt::tab_header(
          title = gt::md(paste0("**", title, "**")),
          subtitle = paste0("Subgroup: ", sg_label)
        ) |>
        gt::fmt_number(
          columns = "Value",
          decimals = digits
        ) |>
        gt::cols_label(
          Metric = "Metric",
          Description = "Description",
          Value = if (show_percentages) "Value (%)" else "Value"
        ) |>
        gt::tab_style(
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_column_labels()
        ) |>
        gt::tab_style(
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_row_groups()
        )

      # Add footnote with metadata
      footnote_text <- paste0(
        "Based on ", sims, " simulation(s)",
        if (!is.na(Kfolds)) paste0(" with ", Kfolds, "-fold CV") else "",
        ". Values are ", if (sims > 1) "medians" else "proportions",
        if (show_percentages) " shown as percentages." else "."
      )

      gt_table <- gt_table |>
        gt::tab_footnote(
          footnote = footnote_text
        )

      result <- gt_table
    }

  } else if (table_style == "separate") {

    if (!use_gt) {
      result <- list(
        agreement_table = agreement_df,
        finding_table = finding_df
      )
    } else {
      # Agreement table
      agreement_gt <- agreement_df |>
        gt::gt() |>
        gt::tab_header(
          title = gt::md("**Agreement Metrics**"),
          subtitle = paste0("Subgroup: ", sg_label)
        ) |>
        gt::fmt_number(columns = "Value", decimals = digits) |>
        gt::cols_label(
          Metric = "Metric",
          Description = "Description",
          Value = if (show_percentages) "Value (%)" else "Value"
        ) |>
        gt::tab_style(
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_column_labels()
        )

      # Finding table
      finding_gt <- finding_df |>
        gt::gt() |>
        gt::tab_header(
          title = gt::md("**Subgroup Finding Metrics**"),
          subtitle = paste0("Subgroup: ", sg_label)
        ) |>
        gt::fmt_number(columns = "Value", decimals = digits) |>
        gt::cols_label(
          Metric = "Metric",
          Description = "Description",
          Value = if (show_percentages) "Value (%)" else "Value"
        ) |>
        gt::tab_style(
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_column_labels()
        )

      result <- list(
        agreement_table = agreement_gt,
        finding_table = finding_gt
      )
    }

  } else if (table_style == "minimal") {

    # Single-row compact summary
    minimal_df <- data.frame(
      Simulations = sims,
      Kfolds = Kfolds,
      `Any Found` = find_summary["Any"] * mult,
      `Exact Match` = find_summary["Exact"] * mult,
      `Sens H` = sens_summary["sens_H"] * mult,
      `Sens Hc` = sens_summary["sens_Hc"] * mult,
      `PPV H` = sens_summary["ppv_H"] * mult,
      `PPV Hc` = sens_summary["ppv_Hc"] * mult,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    if (!use_gt) {
      result <- minimal_df
    } else {
      gt_table <- minimal_df |>
        gt::gt() |>
        gt::tab_header(
          title = gt::md(paste0("**", title, "**")),
          subtitle = paste0("Subgroup: ", sg_label)
        ) |>
        gt::fmt_number(
          columns = c("Any Found", "Exact Match", "Sens H", "Sens Hc", "PPV H", "PPV Hc"),
          decimals = digits
        ) |>
        gt::fmt_integer(columns = c("Simulations", "Kfolds")) |>
        gt::tab_spanner(
          label = "Finding",
          columns = c("Any Found", "Exact Match")
        ) |>
        gt::tab_spanner(
          label = "Agreement",
          columns = c("Sens H", "Sens Hc", "PPV H", "PPV Hc")
        ) |>
        gt::tab_style(
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_column_labels()
        )

      result <- gt_table
    }
  }

  # ===========================================================================
  # INCLUDE RAW DATA IF REQUESTED
  # ===========================================================================

  if (include_raw && table_style != "minimal") {
    if (is.list(result)) {
      result$sens_out <- cv_result$sens_out
      result$find_out <- cv_result$find_out
      result$metadata <- list(
        sims = sims,
        Kfolds = Kfolds,
        sg_definition = sg_definition
      )
    } else {
      result <- list(
        table = result,
        sens_out = cv_result$sens_out,
        find_out = cv_result$find_out,
        metadata = list(
          sims = sims,
          Kfolds = Kfolds,
          sg_definition = sg_definition
        )
      )
    }
  }

  return(result)
}


#' Create Compact CV Summary Text
#'
#' Generates a compact text string summarizing CV results, suitable for
#' annotations in plots or reports.
#'
#' @param cv_result List. Result from \code{forestsearch_tenfold()} or
#'   \code{forestsearch_Kfold()}.
#' @param est.scale Character. "hr" or "1/hr" to determine label orientation.
#'   Default: "hr".
#' @param include_finding Logical. Include subgroup finding rate. Default: TRUE.
#' @param include_agreement Logical
#'
#' @return Character string with formatted CV metrics.
#'
#' @examples
#' \dontrun{
#' cv_text <- cv_summary_text(tenfold_results)
#' # Returns: "CV found = 95%, Agree(+,-) = 88%, 92%"
#' }
#'
#' @export
cv_summary_text <- function(
    cv_result,
    est.scale = "hr",
    include_finding = TRUE,
    include_agreement = TRUE
) {

  if (is.null(cv_result)) return(NULL)

  parts <- character(0)

  # Finding rate
  if (include_finding && !is.null(cv_result$find_summary)) {
    cv_found <- cv_result$find_summary["Any"]
    parts <- c(parts, paste0("CV found = ", round(100 * cv_found, 0), "%"))
  }

  # Agreement rates
  if (include_agreement && !is.null(cv_result$sens_summary)) {
    if (est.scale == "hr") {
      sens_pos <- cv_result$sens_summary["sens_H"]
      sens_neg <- cv_result$sens_summary["sens_Hc"]
    } else {
      sens_pos <- cv_result$sens_summary["sens_Hc"]
      sens_neg <- cv_result$sens_summary["sens_H"]
    }

    agree_text <- paste0(
      "Agree(+,-) = ",
      round(100 * sens_neg, 0), "%, ",
      round(100 * sens_pos, 0), "%"
    )
    parts <- c(parts, agree_text)
  }

  paste(parts, collapse = ", ")
}


#' Compare Multiple CV Results
#'
#' Creates a comparison table from multiple cross-validation runs with
#' different configurations.
#'
#' @param cv_list Named list of cv_result objects from \code{forestsearch_tenfold()}
#'   or \code{forestsearch_Kfold()}.
#' @param metrics Character vector. Which metrics to include. Options:
#'   "finding", "agreement", "all". Default: "all".
#' @param show_percentages Logical. Display as percentages. Default: TRUE.
#' @param digits Integer. Decimal places. Default: 1.
#' @param use_gt Logical. Return gt table if TRUE. Default: TRUE.
#'
#' @return A gt table or data.frame comparing CV results across configurations.
#'
#' @examples
#' \dontrun{
#' # Compare CV results from different configurations
#' cv_comparison <- cv_compare_results(
#'   cv_list = list(
#'     "maxk=1" = cv_maxk1,
#'     "maxk=2" = cv_maxk2,
#'     "maxk=3" = cv_maxk3
#'   )
#' )
#' }
#'
#' @export
cv_compare_results <- function(
    cv_list,
    metrics = c("all", "finding", "agreement"),
    show_percentages = TRUE,
    digits = 1,
    use_gt = TRUE
) {

  metrics <- match.arg(metrics)

  if (!is.list(cv_list) || length(cv_list) == 0) {
    stop("cv_list must be a non-empty named list of cv_result objects")
  }

  # Ensure names exist
  if (is.null(names(cv_list))) {
    names(cv_list) <- paste0("Config_", seq_along(cv_list))
  }

  mult <- if (show_percentages) 100 else 1

  # Build comparison data frame
  comparison_rows <- lapply(names(cv_list), function(config_name) {
    cv <- cv_list[[config_name]]

    row_data <- data.frame(
      Configuration = config_name,
      Sims = if (!is.null(cv$sims)) cv$sims else 1,
      Kfolds = if (!is.null(cv$Kfolds)) cv$Kfolds else NA,
      stringsAsFactors = FALSE
    )

    # Finding metrics
    if (metrics %in% c("all", "finding") && !is.null(cv$find_summary)) {
      row_data$Any_Found <- cv$find_summary["Any"] * mult
      row_data$Exact_Match <- cv$find_summary["Exact"] * mult
      row_data$At_Least_1 <- cv$find_summary["At least 1"] * mult
    }

    # Agreement metrics
    if (metrics %in% c("all", "agreement") && !is.null(cv$sens_summary)) {
      row_data$Sens_H <- cv$sens_summary["sens_H"] * mult
      row_data$Sens_Hc <- cv$sens_summary["sens_Hc"] * mult
      row_data$PPV_H <- cv$sens_summary["ppv_H"] * mult
      row_data$PPV_Hc <- cv$sens_summary["ppv_Hc"] * mult
    }

    row_data
  })

  comparison_df <- do.call(rbind, comparison_rows)
  rownames(comparison_df) <- NULL

  if (!use_gt) {
    return(comparison_df)
  }

  # Create gt table
  gt_table <- comparison_df |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Cross-Validation Comparison**"),
      subtitle = paste0(length(cv_list), " configurations compared")
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    )

  # Format numeric columns
  numeric_cols <- setdiff(names(comparison_df), c("Configuration", "Sims", "Kfolds"))

  if (length(numeric_cols) > 0) {
    gt_table <- gt_table |>
      gt::fmt_number(columns = gt::all_of(numeric_cols), decimals = digits)
  }

  # Format integer columns
  gt_table <- gt_table |>
    gt::fmt_integer(columns = c("Sims", "Kfolds"))

  # Add column spanners
  finding_cols <- intersect(c("Any_Found", "Exact_Match", "At_Least_1"), names(comparison_df))
  agreement_cols <- intersect(c("Sens_H", "Sens_Hc", "PPV_H", "PPV_Hc"), names(comparison_df))

  if (length(finding_cols) > 0) {
    gt_table <- gt_table |>
      gt::tab_spanner(label = "Finding", columns = gt::all_of(finding_cols))
  }

  if (length(agreement_cols) > 0) {
    gt_table <- gt_table |>
      gt::tab_spanner(label = "Agreement", columns = gt::all_of(agreement_cols))
  }

  # Rename columns for display
  gt_table <- gt_table |>
    gt::cols_label(
      Configuration = "Config",
      .fn = ~ gsub("_", " ", .x)
    )

  gt_table
}
