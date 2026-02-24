# =============================================================================
# forestsearch_helpers.R - Utility Functions for ForestSearch
# =============================================================================
#
# Helper functions used across the ForestSearch package:
#   - add_id_column()          : Ensure data frame has a unique ID column
#   - get_dfpred()             : Apply subgroup definition to new data
#   - evaluate_comparison()    : Safe operator-dispatch expression evaluator
#   - get_param()              : Extract parameter with default fallback
# =============================================================================


#' Add ID Column to Data Frame
#'
#' Ensures that a data frame has a unique ID column. If \code{id.name} is not
#' provided, a column named \code{"id"} is added. If \code{id.name} is provided
#' but does not exist in the data frame, it is created with unique integer
#' values.
#'
#' @param df.analysis Data frame to which the ID column will be added.
#' @param id.name Character. Name of the ID column to add (default is
#'   \code{NULL}, which uses \code{"id"}).
#'
#' @return Data frame with the ID column added if necessary.
#' @keywords internal
add_id_column <- function(df.analysis, id.name = NULL) {
  if (is.null(id.name)) {
    df.analysis$id <- seq_len(nrow(df.analysis))
    id.name <- "id"
  } else if (!(id.name %in% names(df.analysis))) {
    df.analysis[[id.name]] <- seq_len(nrow(df.analysis))
  }
  return(df.analysis)
}


# =============================================================================
# SUBGROUP APPLICATION
# =============================================================================

#' Generate Prediction Dataset with Subgroup Treatment Recommendation
#'
#' Creates a prediction dataset with a treatment recommendation flag based
#' on the subgroup definition. Supports both label expressions
#' (e.g., \code{"\{er <= 0\}"}) and bare column names (e.g., \code{"q3.1"}).
#'
#' Each element of \code{sg.harm} is processed as follows:
#' \enumerate{
#'   \item Outer braces and leading \code{!} are stripped.
#'   \item If the result matches \code{"var op value"} (where \code{op} is
#'     one of \code{<=}, \code{<}, \code{>=}, \code{>}, \code{==},
#'     \code{!=}), the comparison is executed directly on
#'     \code{df.predict[[var]]}.
#'   \item Otherwise the expression is treated as a column name and
#'     membership is \code{df.predict[[name]] == 1}.
#' }
#'
#' @param df.predict Data frame for prediction (test or validation set).
#' @param sg.harm Character vector of subgroup-defining labels. Values may
#'   be wrapped in braces and optionally negated, e.g. \code{"\{er <= 0\}"}
#'   or \code{"!\{size <= 35\}"}. Plain column names (e.g., \code{"q3.1"})
#'   are treated as binary indicators that must equal 1.
#' @param version Integer; encoding version (maintained for backward
#'   compatibility). Default: 1.
#'
#' @return Data frame with treatment recommendation flag
#'   (\code{treat.recommend}): 0 for harm subgroup, 1 for complement.
#'
#' @seealso \code{\link{evaluate_comparison}} for the operator-dispatch
#'   logic, \code{\link{forestsearch}} for the main analysis function.
#'
#' @examples
#' \dontrun{
#' # With brace-wrapped label expressions
#' sg <- c("{er <= 0}", "{size <= 35}")
#' df_out <- get_dfpred(df.predict = test_data, sg.harm = sg)
#'
#' # With negation
#' sg_neg <- c("{er <= 0}", "!{size <= 35}")
#' df_neg <- get_dfpred(df.predict = test_data, sg.harm = sg_neg)
#'
#' # With bare column names (binary indicators)
#' sg_col <- c("q1.1", "q3.1")
#' df_col <- get_dfpred(df.predict = encoded_data, sg.harm = sg_col)
#' }
#'
#' @export
get_dfpred <- function(df.predict, sg.harm, version = 1) {

  df.pred <- df.predict
  labels <- if (!is.null(names(sg.harm))) unname(sg.harm) else sg.harm

  # Build membership indicator for each factor
  in_harm <- rep(TRUE, nrow(df.pred))

  for (lab in labels) {
    # Strip braces: "{er <= 0}" -> "er <= 0"
    clean <- gsub("^!?\\{(.*)\\}$", "\\1", lab)
    is_negated <- grepl("^!", lab)

    member <- evaluate_comparison(clean, df.pred)

    # Apply negation if needed
    if (is_negated) member <- !member

    in_harm <- in_harm & member
  }

  df.pred$treat.recommend <- ifelse(in_harm, 0L, 1L)
  df.pred
}


#' Evaluate a Comparison Expression Without eval(parse())
#'
#' Parses a string of the form \code{"var op value"} and evaluates it
#' directly against a data frame column using operator dispatch. Falls back
#' to column-name lookup for bare names.
#'
#' @param expr Character. An expression like \code{"er <= 0"},
#'   \code{"size > 35"}, \code{"grade3 == 1"}, or a bare column name
#'   like \code{"q3.1"}.
#' @param df Data frame whose columns are referenced by \code{expr}.
#'
#' @return Logical vector of length \code{nrow(df)}.
#'
#' @details
#' Supported operators (matched longest-first to avoid partial-match
#' ambiguity): \code{<=}, \code{>=}, \code{!=}, \code{==}, \code{<},
#' \code{>}.
#'
#' If no operator is found, \code{expr} is treated as a column name and
#' the result is \code{df[[expr]] == 1}.
#'
#' The value on the right-hand side is coerced to numeric when possible,
#' otherwise kept as character for string comparisons.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(er = c(-1, 0, 1, 2), size = c(10, 20, 30, 40))
#' evaluate_comparison("er <= 0", df)
#' # [1]  TRUE  TRUE FALSE FALSE
#'
#' evaluate_comparison("size > 25", df)
#' # [1] FALSE FALSE  TRUE  TRUE
#' }
#'
#' @export
evaluate_comparison <- function(expr, df) {

  expr <- trimws(expr)

  # Operators ordered longest-first to avoid partial matching
  # (e.g., "<=" must be tried before "<")
  ops <- c("<=", ">=", "!=", "==", "<", ">")

  for (op in ops) {
    if (grepl(op, expr, fixed = TRUE)) {
      parts <- strsplit(expr, op, fixed = TRUE)[[1L]]
      if (length(parts) != 2L) next

      var_name <- trimws(parts[1L])
      value <- trimws(parts[2L])

      if (!var_name %in% names(df)) {
        warning("Column '", var_name, "' not found in data frame",
                call. = FALSE)
        return(rep(NA, nrow(df)))
      }

      col <- df[[var_name]]
      val <- suppressWarnings(as.numeric(value))
      if (is.na(val)) val <- value # keep as character for string comparisons

      result <- switch(
        op,
        "<=" = col <= val,
        ">=" = col >= val,
        "!=" = col != val,
        "==" = col == val,
        "<"  = col <  val,
        ">"  = col >  val
      )

      return(result)
    }
  }

  # No operator found - treat as a column name (binary indicator)
  if (expr %in% names(df)) {
    return(df[[expr]] == 1)
  }

  warning("Could not parse expression and column '", expr,
          "' not found in data frame", call. = FALSE)
  rep(NA, nrow(df))
}


# =============================================================================
# PARAMETER HELPERS
# =============================================================================

#' Get Parameter with Default Fallback
#'
#' Safely retrieves a named element from a list, returning a default value
#' if the element is missing or \code{NULL}.
#'
#' @param args_list List to extract from.
#' @param param_name Character. Name of the element to retrieve.
#' @param default_value Default value to return if element is missing or
#'   \code{NULL}.
#'
#' @return The value of \code{args_list[[param_name]]} if present and
#'   non-\code{NULL}, otherwise \code{default_value}.
#'
#' @keywords internal
get_param <- function(args_list, param_name, default_value) {
  if (hasName(args_list, param_name) && !is.null(args_list[[param_name]])) {
    return(args_list[[param_name]])
  }
  return(default_value)
}


# ─────────────────────────────────────────────────────────────────────
# Core helper: evaluates any expression string against a data frame
# ─────────────────────────────────────────────────────────────────────

#' Evaluate an expression string in a data-frame scope
#'
#' Parses and evaluates \code{expr} in a restricted environment
#' containing only the columns of \code{df} (parent: \code{baseenv()}).
#' This isolates evaluation from the global environment, reducing
#' scope for unintended side effects.
#'
#' @param df Data frame providing column names as variables.
#' @param expr Character. Expression to evaluate
#'   (e.g., \code{"BM > 1 & tmrsize > 19"}).
#'
#' @return Result of evaluating \code{expr}, or \code{NULL} on failure.
#'
#' @keywords internal
safe_eval_expr <- function(df, expr) {
  tryCatch({
    env <- list2env(as.list(df), parent = baseenv())
    eval(parse(text = expr), envir = env)
  }, error = function(e) {
    warning(
      "Failed to evaluate expression: '", expr, "' - ", e$message,
      call. = FALSE
    )
    NULL
  })
}


# ─────────────────────────────────────────────────────────────────────
# Convenience wrapper: subset rows using an expression string
# ─────────────────────────────────────────────────────────────────────

#' Subset a data frame using an expression string
#'
#' Thin wrapper around \code{\link{safe_eval_expr}} that uses the
#' logical result to subset rows.
#'
#' @param df Data frame.
#' @param expr Character. Subset expression
#'   (e.g., \code{"BM > 1 & tmrsize > 19"}).
#'
#' @return Subset of \code{df}, or \code{NULL} on failure.
#'
#' @keywords internal
safe_subset <- function(df, expr) {
  idx <- safe_eval_expr(df, expr)
  if (is.null(idx)) return(NULL)
  df[idx, , drop = FALSE]
}



