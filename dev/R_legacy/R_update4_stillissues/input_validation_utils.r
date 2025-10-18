#' Input Validation Utility Functions
#'
#' A collection of standardized input validation functions to ensure 
#' consistent error handling across the ForestSearch package.
#'
#' @name validation_utils
NULL

#' Validate Data Frame Input
#'
#' Checks that input is a valid data frame with minimum required rows.
#'
#' @param df Input object to validate.
#' @param arg_name Character. Name of the argument for error messages.
#' @param min_rows Integer. Minimum number of rows required (default: 1).
#' @param required_cols Character vector. Required column names (optional).
#' @return Invisible TRUE if validation passes, stops with error otherwise.
#' @export
validate_dataframe <- function(df, arg_name = "df", min_rows = 1, required_cols = NULL) {
  # Check if data frame
  if (!is.data.frame(df)) {
    stop(sprintf("'%s' must be a data.frame, got %s", arg_name, class(df)[1]))
  }
  
  # Check minimum rows
  if (nrow(df) < min_rows) {
    stop(sprintf("'%s' must have at least %d rows, has %d", arg_name, min_rows, nrow(df)))
  }
  
  # Check required columns
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop(sprintf("'%s' missing required columns: %s", 
                   arg_name, paste(missing_cols, collapse = ", ")))
    }
  }
  
  invisible(TRUE)
}

#' Validate Numeric Input
#'
#' Checks that input is numeric with proper bounds.
#'
#' @param x Input object to validate.
#' @param arg_name Character. Name of the argument for error messages.
#' @param length Integer. Expected length (NULL for any length).
#' @param min Numeric. Minimum value (NULL for no minimum).
#' @param max Numeric. Maximum value (NULL for no maximum).
#' @param allow_na Logical. Allow NA values (default: FALSE).
#' @param allow_inf Logical. Allow infinite values (default: FALSE).
#' @return Invisible TRUE if validation passes, stops with error otherwise.
#' @export
validate_numeric <- function(x, arg_name = "x", length = NULL, min = NULL, 
                             max = NULL, allow_na = FALSE, allow_inf = FALSE) {
  # Check if numeric
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be numeric, got %s", arg_name, class(x)[1]))
  }
  
  # Check length
  if (!is.null(length) && length(x) != length) {
    stop(sprintf("'%s' must have length %d, has length %d", 
                 arg_name, length, length(x)))
  }
  
  # Check for NA
  if (!allow_na && any(is.na(x))) {
    stop(sprintf("'%s' contains NA values", arg_name))
  }
  
  # Check for infinite
  if (!allow_inf && any(is.infinite(x[!is.na(x)]))) {
    stop(sprintf("'%s' contains infinite values", arg_name))
  }
  
  # Check minimum
  if (!is.null(min) && any(x[!is.na(x)] < min)) {
    stop(sprintf("'%s' must be >= %s, found values < %s", arg_name, min, min))
  }
  
  # Check maximum
  if (!is.null(max) && any(x[!is.na(x)] > max)) {
    stop(sprintf("'%s' must be <= %s, found values > %s", arg_name, max, max))
  }
  
  invisible(TRUE)
}

#' Validate Character Input
#'
#' Checks that input is character with proper constraints.
#'
#' @param x Input object to validate.
#' @param arg_name Character. Name of the argument for error messages.
#' @param length Integer. Expected length (NULL for any length).
#' @param allowed Character vector. Allowed values (NULL for any value).
#' @param allow_na Logical. Allow NA values (default: FALSE).
#' @return Invisible TRUE if validation passes, stops with error otherwise.
#' @export
validate_character <- function(x, arg_name = "x", length = NULL, 
                               allowed = NULL, allow_na = FALSE) {
  # Check if character
  if (!is.character(x)) {
    stop(sprintf("'%s' must be character, got %s", arg_name, class(x)[1]))
  }
  
  # Check length
  if (!is.null(length) && length(x) != length) {
    stop(sprintf("'%s' must have length %d, has length %d", 
                 arg_name, length, length(x)))
  }
  
  # Check for NA
  if (!allow_na && any(is.na(x))) {
    stop(sprintf("'%s' contains NA values", arg_name))
  }
  
  # Check allowed values
  if (!is.null(allowed)) {
    invalid <- setdiff(x[!is.na(x)], allowed)
    if (length(invalid) > 0) {
      stop(sprintf("'%s' contains invalid values: %s. Allowed: %s", 
                   arg_name, paste(invalid, collapse = ", "),
                   paste(allowed, collapse = ", ")))
    }
  }
  
  invisible(TRUE)
}

#' Validate Logical Input
#'
#' Checks that input is logical.
#'
#' @param x Input object to validate.
#' @param arg_name Character. Name of the argument for error messages.
#' @param length Integer. Expected length (NULL for any length).
#' @return Invisible TRUE if validation passes, stops with error otherwise.
#' @export
validate_logical <- function(x, arg_name = "x", length = 1) {
  # Check if logical
  if (!is.logical(x)) {
    stop(sprintf("'%s' must be logical, got %s", arg_name, class(x)[1]))
  }
  
  # Check length
  if (!is.null(length) && length(x) != length) {
    stop(sprintf("'%s' must have length %d, has length %d", 
                 arg_name, length, length(x)))
  }
  
  # Check for NA
  if (any(is.na(x))) {
    stop(sprintf("'%s' contains NA values", arg_name))
  }
  
  invisible(TRUE)
}

#' Validate Column Names
#'
#' Checks that specified column names exist in a data frame.
#'
#' @param df Data frame.
#' @param col_names Character vector of column names to check.
#' @param df_name Character. Name of the data frame for error messages.
#' @return Invisible TRUE if validation passes, stops with error otherwise.
#' @export
validate_columns <- function(df, col_names, df_name = "df") {
  validate_dataframe(df, df_name)
  
  missing <- setdiff(col_names, names(df))
  if (length(missing) > 0) {
    stop(sprintf("'%s' missing required columns: %s", 
                 df_name, paste(missing, collapse = ", ")))
  }
  
  invisible(TRUE)
}

#' Validate Matrix Input
#'
#' Checks that input is a valid matrix with proper dimensions.
#'
#' @param x Input object to validate.
#' @param arg_name Character. Name of the argument for error messages.
#' @param nrow Integer. Expected number of rows (NULL for any).
#' @param ncol Integer. Expected number of columns (NULL for any).
#' @param min_nrow Integer. Minimum number of rows (NULL for no minimum).
#' @param min_ncol Integer. Minimum number of columns (NULL for no minimum).
#' @return Invisible TRUE if validation passes, stops with error otherwise.
#' @export
validate_matrix <- function(x, arg_name = "x", nrow = NULL, ncol = NULL, 
                            min_nrow = NULL, min_ncol = NULL) {
  # Check if matrix
  if (!is.matrix(x)) {
    stop(sprintf("'%s' must be a matrix, got %s", arg_name, class(x)[1]))
  }
  
  # Check exact dimensions
  if (!is.null(nrow) && nrow(x) != nrow) {
    stop(sprintf("'%s' must have %d rows, has %d", arg_name, nrow, nrow(x)))
  }
  
  if (!is.null(ncol) && ncol(x) != ncol) {
    stop(sprintf("'%s' must have %d columns, has %d", arg_name, ncol, ncol(x)))
  }
  
  # Check minimum dimensions
  if (!is.null(min_nrow) && nrow(x) < min_nrow) {
    stop(sprintf("'%s' must have at least %d rows, has %d", 
                 arg_name, min_nrow, nrow(x)))
  }
  
  if (!is.null(min_ncol) && ncol(x) < min_ncol) {
    stop(sprintf("'%s' must have at least %d columns, has %d", 
                 arg_name, min_ncol, ncol(x)))
  }
  
  invisible(TRUE)
}

#' Validate Formula
#'
#' Checks that input is a valid formula.
#'
#' @param formula Input object to validate.
#' @param arg_name Character. Name of the argument for error messages.
#' @return Invisible TRUE if validation passes, stops with error otherwise.
#' @export
validate_formula <- function(formula, arg_name = "formula") {
  if (!inherits(formula, "formula")) {
    stop(sprintf("'%s' must be a formula", arg_name))
  }
  invisible(TRUE)
}

#' Safe Division with Division-by-Zero Handling
#'
#' Performs division with proper handling of division by zero.
#'
#' @param numerator Numeric. Numerator value(s).
#' @param denominator Numeric. Denominator value(s).
#' @param replace_with Value to return when denominator is zero (default: NA).
#' @return Numeric result of division, or replace_with when denominator is zero.
#' @export
safe_divide <- function(numerator, denominator, replace_with = NA) {
  validate_numeric(numerator, "numerator")
  validate_numeric(denominator, "denominator")
  
  result <- ifelse(denominator == 0, replace_with, numerator / denominator)
  return(result)
}

#' Calculate Proportion with Division-by-Zero Handling
#'
#' Calculates proportion (count / total) with proper handling of edge cases.
#'
#' @param count Numeric. Count value(s).
#' @param total Numeric. Total value(s).
#' @param replace_with Value to return when total is zero (default: NA).
#' @return Numeric proportion, or replace_with when total is zero.
#' @export
safe_proportion <- function(count, total, replace_with = NA) {
  validate_numeric(count, "count", min = 0)
  validate_numeric(total, "total", min = 0)
  
  result <- ifelse(total == 0, replace_with, count / total)
  return(result)
}

#' Calculate Sensitivity with Division-by-Zero Handling
#'
#' Calculates sensitivity (TP / (TP + FN)) with proper handling of edge cases.
#'
#' @param true_positive Numeric. True positive count.
#' @param false_negative Numeric. False negative count.
#' @param replace_with Value to return when denominator is zero (default: NA).
#' @return Numeric sensitivity, or replace_with when denominator is zero.
#' @export
calculate_sensitivity <- function(true_positive, false_negative, replace_with = NA) {
  validate_numeric(true_positive, "true_positive", min = 0)
  validate_numeric(false_negative, "false_negative", min = 0)
  
  denominator <- true_positive + false_negative
  result <- ifelse(denominator == 0, replace_with, true_positive / denominator)
  return(result)
}

#' Calculate Specificity with Division-by-Zero Handling
#'
#' Calculates specificity (TN / (TN + FP)) with proper handling of edge cases.
#'
#' @param true_negative Numeric. True negative count.
#' @param false_positive Numeric. False positive count.
#' @param replace_with Value to return when denominator is zero (default: NA).
#' @return Numeric specificity, or replace_with when denominator is zero.
#' @export
calculate_specificity <- function(true_negative, false_positive, replace_with = NA) {
  validate_numeric(true_negative, "true_negative", min = 0)
  validate_numeric(false_positive, "false_positive", min = 0)
  
  denominator <- true_negative + false_positive
  result <- ifelse(denominator == 0, replace_with, true_negative / denominator)
  return(result)
}

#' Calculate PPV with Division-by-Zero Handling
#'
#' Calculates positive predictive value (TP / (TP + FP)) with proper handling of edge cases.
#'
#' @param true_positive Numeric. True positive count.
#' @param false_positive Numeric. False positive count.
#' @param replace_with Value to return when denominator is zero (default: NA).
#' @return Numeric PPV, or replace_with when denominator is zero.
#' @export
calculate_ppv <- function(true_positive, false_positive, replace_with = NA) {
  validate_numeric(true_positive, "true_positive", min = 0)
  validate_numeric(false_positive, "false_positive", min = 0)
  
  denominator <- true_positive + false_positive
  result <- ifelse(denominator == 0, replace_with, true_positive / denominator)
  return(result)
}

#' Calculate NPV with Division-by-Zero Handling
#'
#' Calculates negative predictive value (TN / (TN + FN)) with proper handling of edge cases.
#'
#' @param true_negative Numeric. True negative count.
#' @param false_negative Numeric. False negative count.
#' @param replace_with Value to return when denominator is zero (default: NA).
#' @return Numeric NPV, or replace_with when denominator is zero.
#' @export
calculate_npv <- function(true_negative, false_negative, replace_with = NA) {
  validate_numeric(true_negative, "true_negative", min = 0)
  validate_numeric(false_negative, "false_negative", min = 0)
  
  denominator <- true_negative + false_negative
  result <- ifelse(denominator == 0, replace_with, true_negative / denominator)
  return(result)
}

#' Standardize Variable Names
#'
#' Converts variable names to a standardized format (lowercase with underscores).
#'
#' @param names Character vector. Variable names to standardize.
#' @param style Character. Naming style: "snake_case" (default), "camelCase", or "original".
#' @return Character vector of standardized names.
#' @export
standardize_names <- function(names, style = "snake_case") {
  validate_character(names, "names")
  validate_character(style, "style", length = 1, 
                    allowed = c("snake_case", "camelCase", "original"))
  
  if (style == "original") return(names)
  
  if (style == "snake_case") {
    # Convert to snake_case
    names <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", names)  # camelCase to snake_case
    names <- gsub("([A-Z])([A-Z][a-z])", "\\1_\\2", names)  # Handle acronyms
    names <- gsub("\\.", "_", names)  # dots to underscores
    names <- gsub("-", "_", names)  # dashes to underscores
    names <- tolower(names)
    names <- gsub("_+", "_", names)  # Remove duplicate underscores
    names <- gsub("^_|_$", "", names)  # Remove leading/trailing underscores
  }
  
  if (style == "camelCase") {
    # Convert to camelCase
    names <- gsub("[._-]([a-z])", "\\U\\1", names, perl = TRUE)
    names <- gsub("^([A-Z])", "\\L\\1", names, perl = TRUE)
  }
  
  return(names)
}

#' Check Variable Name Consistency
#'
#' Checks if variable naming is consistent across a vector of names.
#'
#' @param names Character vector. Variable names to check.
#' @return List with style detected and consistency flag.
#' @export
check_name_consistency <- function(names) {
  validate_character(names, "names")
  
  has_underscore <- any(grepl("_", names))
  has_dot <- any(grepl("\\.", names))
  has_camelCase <- any(grepl("[a-z][A-Z]", names))
  
  styles <- c()
  if (has_underscore) styles <- c(styles, "snake_case")
  if (has_dot) styles <- c(styles, "dot.case")
  if (has_camelCase) styles <- c(styles, "camelCase")
  
  consistent <- length(styles) <= 1
  
  list(
    styles_detected = if (length(styles) == 0) "lowercase" else styles,
    is_consistent = consistent,
    recommendation = if (!consistent) "Consider using a single naming convention" else "Names are consistent"
  )
}