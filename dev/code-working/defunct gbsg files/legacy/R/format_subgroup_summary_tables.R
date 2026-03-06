#' Format Subgroup Summary Tables with gt
#'
#' Creates publication-ready gt tables for bootstrap subgroup analysis
#'
#' @param subgroup_summary List from summarize_bootstrap_subgroups()
#' @param nb_boots Integer. Number of bootstrap iterations
#'
#' @return List of gt table objects
#' @importFrom gt gt tab_header tab_spanner tab_footnote md
#' @keywords internal

format_subgroup_summary_tables <- function(subgroup_summary, nb_boots) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' required")
  }

  if (is.null(subgroup_summary)) {
    return(NULL)
  }

  tables <- list()

  # =========================================================================
  # TABLE 1: BASIC STATISTICS
  # =========================================================================

  if (!is.null(subgroup_summary$basic_stats)) {
    tables$basic_stats <- subgroup_summary$basic_stats |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Bootstrap Subgroup Characteristics**"),
        subtitle = sprintf("Summary statistics across %d bootstrap iterations", nb_boots)
      ) |>
      gt::cols_label(
        Metric = "",
        Value = "Value"
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#f0f0f0"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(
          rows = Metric %in% c(
            "Consistency (Pcons)",
            "Hazard Ratio (hr_sg)",
            "Subgroup Size (N_sg)",
            "Number of Factors (K_sg)"
          )
        )
      ) |>
      gt::tab_options(
        table.font.size = gt::px(13)
      )
  }

  # =========================================================================
  # TABLE 2: SUBGROUP DEFINITION AGREEMENT
  # =========================================================================

  if (!is.null(subgroup_summary$agreement)) {
    tables$agreement <- subgroup_summary$agreement |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Subgroup Definition Agreement**"),
        subtitle = sprintf("Top subgroup definitions across %d successful iterations",
                           subgroup_summary$n_found)
      ) |>
      gt::cols_label(
        Rank = "Rank",
        Subgroup = "Subgroup Definition",
        K_sg = "K",
        N = "Count",
        Percent_of_successful = "% of Successful"  # CHANGED: was Percent
      ) |>
      gt::fmt_number(
        columns = Percent_of_successful,  # CHANGED: was Percent
        decimals = 1
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(rows = Rank == 1)
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**K** = number of factors in subgroup definition"),
        locations = gt::cells_column_labels(columns = K_sg)
      )
  }



  # =========================================================================
  # TABLE 2C: INDIVIDUAL FACTOR PRESENCE (BASE)
  # =========================================================================

  if (!is.null(subgroup_summary$factor_presence)) {
    tables$factor_presence <- subgroup_summary$factor_presence |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Individual Factor Presence (Base)**"),
        subtitle = "How often each base factor appears (any threshold, any position)"
      ) |>
      gt::cols_label(
        Rank = "Rank",
        Factor = "Factor",
        Count = "Count",
        Percent = "% of Successful"
      ) |>
      gt::fmt_number(
        columns = Percent,
        decimals = 1
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#fff3cd")
        ),
        locations = gt::cells_body(
          rows = Percent >= 20
        )
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**Highlighted rows** appear in >=20% of successful iterations"),
        locations = gt::cells_column_labels(columns = Percent)
      )
  }

  # =========================================================================
  # TABLE 2D: COMMON SPECIFIC FACTOR DEFINITIONS (>= 20%)
  # =========================================================================

  if (!is.null(subgroup_summary$factor_presence_specific) &&
      nrow(subgroup_summary$factor_presence_specific) > 0) {
    tables$factor_presence_specific <- subgroup_summary$factor_presence_specific |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Common Specific Factor Definitions**"),
        subtitle = "Specific factor definitions appearing in >=20% of successful iterations"
      ) |>
      gt::cols_label(
        Rank = "Rank",
        Base_Factor = "Base Factor",
        Factor_Definition = "Specific Definition",
        Count = "Count",
        Percent = "% of Successful"
      ) |>
      gt::fmt_number(
        columns = Percent,
        decimals = 1
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8")
        ),
        locations = gt::cells_body(
          columns = Base_Factor
        )
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**Base Factor**: Variable name extracted from full definition"),
        locations = gt::cells_column_labels(columns = Base_Factor)
      )
  }


  if (!is.null(subgroup_summary$factor_freq)) {
    factor_tables <- create_factor_summary_tables(
      factor_freq = subgroup_summary$factor_freq,
      n_found = subgroup_summary$n_found,
      min_percent = 2
    )

    # Add to output
    tables$factor_freq <- factor_tables$by_position
    tables$factor_overall <- factor_tables$overall
  }

  # =========================================================================
  # TABLE 4: CONSISTENCY DISTRIBUTION
  # =========================================================================


  if (!is.null(subgroup_summary$consistency_dist)) {

    # Check if table has data
    if (nrow(subgroup_summary$consistency_dist) > 0) {

      tables$consistency_dist <- subgroup_summary$consistency_dist |>
        gt::gt() |>
        gt::tab_header(
          title = gt::md("**Consistency (P<sub>cons</sub>) Distribution**"),
          subtitle = "Distribution of consistency scores across successful iterations"
        ) |>

        # Format all numeric columns appropriately
        gt::fmt_number(
          columns = c(Count),  # Integer columns
          decimals = 0
        ) |>
        gt::fmt_number(
          columns = c(Percent),  # Percentage columns
          decimals = 2
        ) |>

        # Add column labels for clarity
        gt::cols_label(
          `Consistency Range` = "Consistency Range",
          Count = "Count",
          Percent = "Percentage (%)"
        ) |>

        # Color coding: High consistency (green)
        gt::tab_style(
          style = list(
            gt::cell_fill(color = "#d4edda"),
            gt::cell_text(weight = "bold")
          ),
          locations = gt::cells_body(
            rows = `Consistency Range` %in% c("0.85-0.95", ">= 0.95")
          )
        ) |>

        # Color coding: Medium consistency (yellow)
        gt::tab_style(
          style = list(
            gt::cell_fill(color = "#fff3cd")
          ),
          locations = gt::cells_body(
            rows = `Consistency Range` %in% c("0.7-0.8", "0.8-0.85")
          )
        ) |>

        # Color coding: Low consistency (red)
        gt::tab_style(
          style = list(
            gt::cell_fill(color = "#f8d7da")
          ),
          locations = gt::cells_body(
            rows = `Consistency Range` %in% c("<0.5", "0.5-0.7")
          )
        ) |>

        # Add helpful footnote
        gt::tab_footnote(
          footnote = "Color coding: Green (>=0.85), Yellow (0.7-0.85), Red (<0.7)",
          locations = gt::cells_column_labels(columns = `Consistency Range`)
        )
    }
  }



  # =========================================================================
  # TABLE 5: ORIGINAL AGREEMENT (if available)
  # =========================================================================

  if (!is.null(subgroup_summary$original_agreement)) {
    tables$original_agreement <- subgroup_summary$original_agreement |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Agreement with Original Subgroup**"),
        subtitle = "How often bootstrap identifies the same subgroup as original analysis"
      ) |>
      gt::cols_label(
        Metric = "",
        Value = "Value"
      ) |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(rows = 2:4)
      )
  }


  return(tables)
}


#' Create Factor Summary Tables from Bootstrap Results
#'
#' @description
#' Generates formatted GT tables summarizing factor frequencies from bootstrap
#' subgroup analysis. Creates two complementary tables: one showing factor
#' selection frequencies within each position (M.1, M.2, etc.), and another
#' showing overall factor frequencies across all positions.
#'
#' @param factor_freq Data.frame or data.table. Factor frequency table from
#'   \code{\link{summarize_bootstrap_subgroups}}, containing columns:
#'   \itemize{
#'     \item \code{Position}: Position identifier (e.g., "M.1", "M.2")
#'     \item \code{Factor}: Factor definition string
#'     \item \code{N}: Count of times factor appeared in this position
#'     \item \code{Percent}: Percentage relative to times position was populated
#'   }
#' @param n_found Integer. Number of successful bootstrap iterations (where a
#'   subgroup was identified). Used to calculate overall percentages.
#' @param min_percent Numeric. Minimum percentage threshold for including factors
#'   in the tables. Factors with selection frequencies below this threshold are
#'   excluded. Default is 2 (i.e., 2%).
#'
#' @return
#' A list with up to two GT table objects:
#' \describe{
#'   \item{\code{by_position}}{GT table showing factor frequencies within each
#'     position. Percentages represent conditional probability of factor selection
#'     given that the position was populated. Within each position, percentages
#'     sum to approximately 100% (may not sum exactly to 100% after filtering).}
#'   \item{\code{overall}}{GT table showing total factor frequencies across all
#'     positions. Includes additional columns indicating which positions each
#'     factor appeared in and how many unique positions used the factor.
#'     Percentages represent proportion of successful iterations where the factor
#'     appeared in any position.}
#' }
#'
#' If no factors meet the minimum threshold, the corresponding table element
#' will be NULL.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{summarize_bootstrap_subgroups}} for generating the factor_freq input
#'   \item \code{\link{format_subgroup_summary_tables}} for creating all subgroup summary tables
#'   \item \code{\link{summarize_bootstrap_results}} for complete bootstrap analysis workflow
#'   \item \code{\link{forestsearch_bootstrap_dofuture}} for running bootstrap analysis
#' }
#'
#' @note
#' This function requires the \pkg{gt} package for table creation. The overall
#' table also requires \pkg{dplyr} for data aggregation. If \pkg{dplyr} is not
#' available, only the position-specific table will be created and the overall
#' element will be NULL.
#'
#' Always check for NULL before using the returned tables:
#' \preformatted{
#' if (!is.null(factor_tables$by_position)) {
#'   print(factor_tables$by_position)
#' }
#' }
#'
#' If all factors have percentages below \code{min_percent}, both table elements
#' will be NULL.
#'
#' @importFrom gt gt tab_header tab_spanner cols_label fmt_number tab_style
#'   cell_fill cell_text cells_body cells_column_labels tab_footnote
#'   tab_source_note md
#' @importFrom dplyr group_by summarise mutate arrange n row_number select
#'
#' @keywords internal
#'
#' @family bootstrap summary functions
#' @family table formatting functions
#'
#' @keywords bootstrap subgroups tables

create_factor_summary_tables <- function(factor_freq, n_found, min_percent = 2) {

  tables <- list()

  if (is.null(factor_freq) || nrow(factor_freq) == 0) {
    return(tables)
  }

  # ========================================================================
  # TABLE 1: By Position (existing)
  # ========================================================================

  factor_freq_filtered <- factor_freq[factor_freq$Percent > min_percent, ]

  if (nrow(factor_freq_filtered) > 0) {

    footnote_text <- paste0(
      "Percentage of times this factor was selected ",
      "when this position was populated"
    )

    tables$by_position <- factor_freq_filtered |>
      gt::gt() |>
      gt::tab_header(
        title = gt::md("**Factor Frequencies by Position**"),
        subtitle = sprintf("Factors with >%.0f%% selection frequency", min_percent)
      ) |>
      gt::cols_label(
        Position = "Position",
        Factor = "Factor",
        N = "Count",
        Percent = "% When Used"
      ) |>
      gt::fmt_number(columns = Percent, decimals = 1) |>
      gt::tab_style(
        style = gt::cell_fill(color = "#f8f9fa"),
        locations = gt::cells_body(rows = Position == "M.1")
      ) |>
      gt::tab_footnote(
        footnote = footnote_text,
        locations = gt::cells_column_labels(columns = Percent)
      )
  }

  # ========================================================================
  # TABLE 2: Overall (NEW)
  # ========================================================================

  if (requireNamespace("dplyr", quietly = TRUE)) {

    factor_overall <- factor_freq |>
      dplyr::group_by(Factor) |>
      dplyr::summarise(
        N_total = sum(N),
        Positions = paste(sort(unique(Position)), collapse = ", "),
        N_positions = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        Percent_of_successful = 100 * N_total / n_found
      ) |>
      dplyr::arrange(desc(N_total)) |>
      dplyr::mutate(Rank = row_number()) |>
      dplyr::select(Rank, Factor, N_total, Percent_of_successful, N_positions, Positions)

    # Filter overall
    factor_overall_filtered <- factor_overall[
      factor_overall$Percent_of_successful > min_percent,
    ]

    if (nrow(factor_overall_filtered) > 0) {

      overall_footnote <- paste0(
        "Total appearances across all positions. ",
        "A factor may appear in multiple positions in the same iteration."
      )

      tables$overall <- factor_overall_filtered |>
        gt::gt() |>
        gt::tab_header(
          title = gt::md("**Overall Factor Frequencies**"),
          subtitle = sprintf(
            "Factors appearing in >%.0f%% of %d successful iterations",
            min_percent,
            n_found
          )
        ) |>
        gt::cols_label(
          Rank = "Rank",
          Factor = "Factor",
          N_total = "Total Count",
          Percent_of_successful = "% of Successful",
          N_positions = "# Positions",
          Positions = "Found In"
        ) |>
        gt::fmt_number(columns = Percent_of_successful, decimals = 1) |>
        gt::tab_style(
          style = list(
            gt::cell_text(weight = "bold"),
            gt::cell_fill(color = "#d4edda")
          ),
          locations = gt::cells_body(
            columns = c(Factor, N_total),
            rows = Percent_of_successful > 50
          )
        ) |>
        gt::tab_footnote(
          footnote = overall_footnote,
          locations = gt::cells_column_labels(columns = N_total)
        ) |>
        gt::tab_footnote(
          footnote = "Number of different positions where this factor appeared",
          locations = gt::cells_column_labels(columns = N_positions)
        )
    }
  }

  return(tables)
}

