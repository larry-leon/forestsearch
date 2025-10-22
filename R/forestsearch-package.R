#' @keywords internal
## usethis namespace: start
#' @keywords internal
"_PACKAGE"

## usethis namespace: end
NULL

#' @importFrom graphics par
#' @importFrom stats as.formula coef na.exclude na.omit predict setNames
#' @importFrom utils hasName head tail
#' @importFrom data.table := .SD .GRP .I .N rbindlist setnames setcolorder is.data.table
#' @importFrom future multisession multicore sequential
NULL

# Note: install.packages is intentionally NOT imported as it's generally not
# appropriate for package code. Consider revising ensure_packages() function
# to use requireNamespace() instead.

#' Suppress R CMD check notes for NSE (Non-Standard Evaluation) variables
#'
#' These variables are used in data.table and dplyr operations where column
#' names are referenced directly without quotes.
#' R CMD check cannot detect these are valid column references, so we declare
#' them here to suppress false positive warnings.
#'
#' @noRd
utils::globalVariables(c(
  # Bootstrap and treatment variables
  "boot",
  "treat.recommend",

  # Cross-validation indices
  "cv_index",
  "cvindex",

  # Estimation and scaling
  "est.scale",
  "insplit1",

  # Confounder information
  "confounders.name",

  # Column name helpers
  "..col_names",

  # Summary table columns
  "Category",
  "Metric",
  "Value",

  # Hazard ratio and sample size metrics
  "HR",
  "hr",
  "K",
  "K_sg",
  "N",
  "N_sg",
  "N_sg_bin",

  # Percentage and ranking metrics
  "Percent",
  "Percent_of_successful",
  "Rank",

  # Consistency metrics
  "Pcons",
  "Pcons_bin",

  # Factor summaries
  "Base_Factor",
  "Factor_Definition",
  "Count",

  # Event and outcome names
  "event.name",
  "outcome.name",
  "treat.name",

  # Harm flag
  "flag.harm",

  # Duplicate detection
  "dup_key",
  "dup_group",
  "dup_count",
  "dup_rank",

  # Grouping variable
  "grp",

  # Bootstrap matching
  "M.1",
  "M.2",
  "M.3",

  # Subgroup identifiers
  "Subgroup",
  "match_string"
))
