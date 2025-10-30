#' Summarize Factor Presence Across Bootstrap Subgroups (Original)
#'
#' @description
#' Analyzes how often each individual factor appears in identified subgroups.
#' This is the original implementation that may encounter data.table ordering
#' errors when column names conflict with objects in the environment.
#'
#' @param results Data.table or data.frame. Bootstrap results with subgroup
#'   characteristics containing columns like Pcons, M.1, M.2, etc.
#' @param maxk Integer. Maximum number of factors allowed (default: 2).
#'   Determines how many M.* columns to examine.
#' @param threshold Numeric. Percentage threshold for including specific
#'   definitions in the output (default: 20). Only factor definitions appearing
#'   in at least this percentage of successful iterations are included.
#'
#' @return A list with three components:
#' \describe{
#'   \item{base_factors}{Data.table with columns: Rank, Factor, Count, Percent.
#'     Shows frequency of base factor names (e.g., "age", "bmi") across all positions.}
#'   \item{specific_factors}{Data.table with columns: Rank, Base_Factor,
#'     Factor_Definition, Count, Percent. Shows specific factor definitions
#'     (e.g., "age <= 50") that appear >= threshold%.}
#'   \item{threshold}{Numeric. The threshold value used for filtering.}
#' }
#'
#' @examples
#' \dontrun{
#' # May fail if Base_Factor or Count exist as objects
#' result <- summarize_factor_presence(bootstrap_results, maxk = 2)
#' }
#'
#' @seealso
#' \code{\link{summarize_factor_presence_safe}} for a robust version
#' \code{\link{summarize_factor_presence_fixed}} for a wrapper with options
#'
#' @keywords internal
summarize_factor_presence_dtable <- function(results, maxk = 2, threshold = 20) {

  # Ensure data.table with proper conversion
  if (!inherits(results, "data.table")) {
    if (is.data.frame(results)) {
      results <- data.table::as.data.table(results)
    } else {
      stop("results must be a data.frame or data.table")
    }
  }

  # Filter successful iterations
  if ("Pcons" %in% names(results)) {
    sg_found <- results[!is.na(results$Pcons), ]
  } else {
    sg_found <- results
  }

  n_found <- nrow(sg_found)
  if (n_found == 0) return(NULL)

  # Extract all factors
  all_factors <- character()
  for (i in 1:maxk) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      factors <- sg_found[[col]][!is.na(sg_found[[col]]) & sg_found[[col]] != ""]
      all_factors <- c(all_factors, factors)
    }
  }

  if (length(all_factors) == 0) return(NULL)

  # Function to extract base factor name
  extract_base_factor <- function(factor_string) {
    factor_string <- trimws(factor_string)
    factor_string <- gsub("^!\\{", "{", factor_string)
    factor_string <- gsub("[{}]", "", factor_string)
    base_name <- sub("^([a-zA-Z_][a-zA-Z0-9_.]*).*", "\\1", factor_string)
    return(base_name)
  }

  # =========================================================================
  # PART 1: BASE FACTOR SUMMARY
  # =========================================================================

  # Count base factors
  base_factors <- sapply(all_factors, extract_base_factor, USE.NAMES = FALSE)
  base_factor_counts <- table(base_factors)

  # Create base factor summary
  base_factor_summary <- data.table::data.table(
    Factor = names(base_factor_counts),
    Count = as.integer(base_factor_counts),
    Percent = 100 * as.integer(base_factor_counts) / n_found
  )

  # FIX: Use setorder instead of subsetting with order
  data.table::setorder(base_factor_summary, -Count)
  base_factor_summary[, Rank := .I]
  data.table::setcolorder(base_factor_summary, c("Rank", "Factor", "Count", "Percent"))

  # =========================================================================
  # PART 2: SPECIFIC FACTOR SUMMARY
  # =========================================================================

  # Count specific factors
  specific_factor_counts <- table(all_factors)

  # Create specific factor summary
  specific_factor_summary <- data.table::data.table(
    Factor_Definition = names(specific_factor_counts),
    Count = as.integer(specific_factor_counts),
    Percent = 100 * as.integer(specific_factor_counts) / n_found
  )

  # Filter by threshold
  specific_factor_summary <- specific_factor_summary[Percent >= threshold, ]

  if (nrow(specific_factor_summary) > 0) {
    # Add base factor column
    specific_factor_summary[, Base_Factor := sapply(Factor_Definition, extract_base_factor)]

    # FIX: Use setorder with proper column references
    # This was the problematic line - now fixed
    data.table::setorder(specific_factor_summary, Base_Factor, -Count)

    # Add rank
    specific_factor_summary[, Rank := .I]

    # Reorder columns
    data.table::setcolorder(specific_factor_summary,
                            c("Rank", "Base_Factor", "Factor_Definition", "Count", "Percent"))
  }

  return(list(
    base_factors = base_factor_summary,
    specific_factors = specific_factor_summary,
    threshold = threshold
  ))
}

#' Summarize Factor Presence - Safe Implementation
#'
#' @description
#' A robust version of \code{summarize_factor_presence} that avoids data.table
#' ordering issues by using base R operations. This function is guaranteed not
#' to fail due to column name conflicts.
#'
#' @param results Data.table or data.frame. Bootstrap results with subgroup
#'   characteristics. Will be converted to data.frame internally for processing.
#' @param maxk Integer. Maximum number of factors allowed (default: 2).
#' @param threshold Numeric. Percentage threshold (default: 20).
#'
#' @return Same as \code{summarize_factor_presence} but guaranteed to work:
#' \describe{
#'   \item{base_factors}{Data.table with base factor frequencies}
#'   \item{specific_factors}{Data.table with specific factor definitions}
#'   \item{threshold}{The threshold value used}
#' }
#'
#' @details
#' Implementation differences from original:
#' \itemize{
#'   \item Converts data.table to data.frame at start
#'   \item Uses base R \code{order()} with explicit column references
#'   \item Converts back to data.table only at the end
#'   \item No reliance on data.table's special evaluation
#' }
#'
#' Performance: Slightly slower than pure data.table operations but more reliable.
#' The performance difference is negligible for typical bootstrap results (<10k rows).
#'
#' @examples
#' # Always works, regardless of environment
#' result <- summarize_factor_presence_safe(bootstrap_results, maxk = 2)
#'
#' # Access results
#' print(result$base_factors)     # Frequency of base factors
#' print(result$specific_factors) # Specific definitions >= 20%
#'
#' @note
#' This is the recommended function for production use when stability is more
#' important than minor performance gains.
#'
#' @seealso
#' \code{\link{summarize_factor_presence}} for the original version
#' \code{\link{summarize_factor_presence_fixed}} for a flexible wrapper
#'
#' @export
summarize_factor_presence_safe <- function(results, maxk = 2, threshold = 20) {

  # Convert to data.frame for safer operations
  if (inherits(results, "data.table")) {
    results <- as.data.frame(results)
  }

  # Filter successful iterations
  if ("Pcons" %in% names(results)) {
    sg_found <- results[!is.na(results$Pcons), ]
  } else {
    sg_found <- results
  }

  n_found <- nrow(sg_found)
  if (n_found == 0) return(NULL)

  # Extract all factors
  all_factors <- character()
  for (i in 1:maxk) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      factors <- sg_found[[col]][!is.na(sg_found[[col]]) & sg_found[[col]] != ""]
      all_factors <- c(all_factors, factors)
    }
  }

  if (length(all_factors) == 0) return(NULL)

  # Function to extract base factor name
  extract_base_factor <- function(factor_string) {
    factor_string <- trimws(factor_string)
    factor_string <- gsub("^!\\{", "{", factor_string)
    factor_string <- gsub("[{}]", "", factor_string)
    base_name <- sub("^([a-zA-Z_][a-zA-Z0-9_.]*).*", "\\1", factor_string)
    return(base_name)
  }

  # =========================================================================
  # PART 1: BASE FACTOR SUMMARY
  # =========================================================================

  base_factors <- sapply(all_factors, extract_base_factor, USE.NAMES = FALSE)
  base_factor_counts <- table(base_factors)

  # Create as data.frame first
  base_factor_summary <- data.frame(
    Factor = names(base_factor_counts),
    Count = as.integer(base_factor_counts),
    Percent = 100 * as.integer(base_factor_counts) / n_found,
    stringsAsFactors = FALSE
  )

  # Sort using base R
  base_factor_summary <- base_factor_summary[order(-base_factor_summary$Count), ]
  base_factor_summary$Rank <- seq_len(nrow(base_factor_summary))
  base_factor_summary <- base_factor_summary[, c("Rank", "Factor", "Count", "Percent")]

  # Convert to data.table at the end
  base_factor_summary <- data.table::as.data.table(base_factor_summary)

  # =========================================================================
  # PART 2: SPECIFIC FACTOR SUMMARY
  # =========================================================================

  specific_factor_counts <- table(all_factors)

  # Create as data.frame first
  specific_factor_summary <- data.frame(
    Factor_Definition = names(specific_factor_counts),
    Count = as.integer(specific_factor_counts),
    Percent = 100 * as.integer(specific_factor_counts) / n_found,
    stringsAsFactors = FALSE
  )

  # Filter by threshold
  specific_factor_summary <- specific_factor_summary[specific_factor_summary$Percent >= threshold, ]

  if (nrow(specific_factor_summary) > 0) {
    # Add base factor
    specific_factor_summary$Base_Factor <- sapply(
      specific_factor_summary$Factor_Definition,
      extract_base_factor
    )

    # Sort using base R - THIS AVOIDS THE ERROR
    ord <- order(specific_factor_summary$Base_Factor, -specific_factor_summary$Count)
    specific_factor_summary <- specific_factor_summary[ord, ]

    # Add rank
    specific_factor_summary$Rank <- seq_len(nrow(specific_factor_summary))

    # Reorder columns
    specific_factor_summary <- specific_factor_summary[, c("Rank", "Base_Factor", "Factor_Definition", "Count", "Percent")]

    # Convert to data.table at the end
    specific_factor_summary <- data.table::as.data.table(specific_factor_summary)
  } else {
    # Empty data.table with correct structure
    specific_factor_summary <- data.table::data.table(
      Rank = integer(),
      Base_Factor = character(),
      Factor_Definition = character(),
      Count = integer(),
      Percent = numeric()
    )
  }

  return(list(
    base_factors = base_factor_summary,
    specific_factors = specific_factor_summary,
    threshold = threshold
  ))
}

#' Summarize Factor Presence - Fixed Wrapper
#'
#' @description
#' A flexible wrapper that can use either the safe (base R) or optimized
#' (data.table) implementation. Provides a migration path and testing capability.
#'
#' @param results Data.table or data.frame. Bootstrap results with subgroup
#'   characteristics.
#' @param maxk Integer. Maximum number of factors allowed (default: 2).
#' @param threshold Numeric. Percentage threshold (default: 20).
#' @param use_safe Logical. If TRUE (default), uses the safe base R implementation.
#'   If FALSE, attempts to use the data.table implementation with fixes.
#'
#' @return Same structure as other variants:
#' \describe{
#'   \item{base_factors}{Data.table with base factor frequencies}
#'   \item{specific_factors}{Data.table with specific factor definitions}
#'   \item{threshold}{The threshold value used}
#' }
#'
#' @details
#' This wrapper function provides:
#' \itemize{
#'   \item Backward compatibility with existing code
#'   \item Ability to test both implementations
#'   \item Gradual migration path from problematic to safe code
#'   \item Performance comparison capability
#' }
#'
#' When to use each mode:
#' \itemize{
#'   \item \code{use_safe = TRUE}: Production environments, automated pipelines
#'   \item \code{use_safe = FALSE}: Testing, performance-critical applications
#'     with controlled environments
#' }
#'
#' @examples
#' # Default: use safe implementation
#' result_safe <- summarize_factor_presence_fixed(bootstrap_results)
#'
#' # Explicitly use safe mode
#' result_safe <- summarize_factor_presence_fixed(
#'   bootstrap_results,
#'   maxk = 2,
#'   use_safe = TRUE
#' )
#'
#' # Try data.table mode (may fail)
#' result_dt <- tryCatch(
#'   summarize_factor_presence_fixed(bootstrap_results, use_safe = FALSE),
#'   error = function(e) {
#'     message("data.table mode failed, using safe mode")
#'     summarize_factor_presence_fixed(bootstrap_results, use_safe = TRUE)
#'   }
#' )
#'
#' @note
#' The default \code{use_safe = TRUE} is recommended unless you have specific
#' performance requirements and a controlled environment.
#'
#' @seealso
#' \code{\link{summarize_factor_presence}} for the original
#' \code{\link{summarize_factor_presence_safe}} for direct safe version
#' \code{\link{test_summarize_factor_presence}} for testing the implementations
#'
#' @export

summarize_factor_presence <- function(results, maxk = 2, threshold = 20, use_safe = TRUE) {
  if (use_safe) {
    # Use the safe implementation that avoids data.table ordering issues
    return(summarize_factor_presence_safe(results, maxk, threshold))
  } else {
    # Use the data.table implementation with fixes
    return(summarize_factor_presence_dtable(results, maxk, threshold))
  }
}
