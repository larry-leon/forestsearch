# =============================================================================
# forestsearch_helpers.R - Utility Functions for ForestSearch
# =============================================================================
#
# Helper functions used across the ForestSearch package:
#   - add_id_column()          : Ensure data frame has a unique ID column
#   - get_dfpred()             : Apply subgroup definition to new data
#   - evaluate_comparison()    : Safe operator-dispatch expression evaluator
#   - get_param()              : Extract parameter with default fallback
#   - collect_results()        : Post-hoc error/result separation for foreach
#   - reset_workers()          : Reset parallel worker pool between sim phases
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

  # GRF union-of-leaves subgroups are carried as a single disjunctive
  # definition string, e.g. "(x1 > 0 & x2 <= 1) | (x1 <= 0 & x3 > 0)", which
  # is NOT an AND-composed factor vector.  Detect that form and evaluate it as
  # a disjunction of conjunctions (each conjunction an AND of "var op value"
  # comparisons) so membership matches the GRF policy tree.  All other inputs
  # fall through to the standard AND-composition below.
  if (length(labels) == 1L && grepl("|", labels, fixed = TRUE)) {
    conj_strings <- strsplit(labels, "\\s*\\|\\s*")[[1]]
    in_any <- rep(FALSE, nrow(df.pred))
    for (cs in conj_strings) {
      cs_clean <- gsub("^\\s*\\(|\\)\\s*$", "", cs)          # strip parens
      comps    <- strsplit(cs_clean, "\\s*&\\s*")[[1]]
      in_cj    <- rep(TRUE, nrow(df.pred))
      for (cmp in comps) {
        cmp_clean  <- gsub("^!?\\{(.*)\\}$", "\\1", trimws(cmp))
        is_negated <- grepl("^!", trimws(cmp))
        member     <- evaluate_comparison(cmp_clean, df.pred)
        if (is_negated) member <- !member
        in_cj <- in_cj & member
      }
      in_any <- in_any | in_cj
    }
    df.pred$treat.recommend <- ifelse(in_any, 0L, 1L)
    return(df.pred)
  }

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
      # Coerce factors with all-numeric levels to numeric so that
      # numeric comparisons (<=, >=, <, >) are meaningful.
      if (is.factor(col)) {
        lvls <- levels(col)
        if (!anyNA(suppressWarnings(as.numeric(lvls))))
          col <- as.numeric(as.character(col))
      }
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


#' Translate a Subgroup Label into an Evaluable Expression String
#'
#' Internal helper used by \code{compute_frontier_cis()} and
#' \code{compute_frontier_member_flags()} to turn the human-readable
#' subgroup-cut labels produced by \code{FS_labels()} into expression
#' strings safe to pass to \code{eval(parse(text = .))} against the
#' analysis data frame.
#'
#' Input is a single \code{c_label} of one of two recognised forms:
#' \itemize{
#'   \item \code{"\{<expr>\}"}      -- subgroup member  (action = 1)
#'   \item \code{"!\{<expr>\}"}     -- subgroup complement (action = 0)
#' }
#' The \code{<expr>} inside the braces is either a comparison
#' (\code{"wtkg <= 86.6"}) or a bare identifier naming a binary
#' indicator factor (\code{"symptom"}).
#'
#' Translation rules:
#' \itemize{
#'   \item \code{"\{<cmp>\}"}    -> \code{"<cmp>"}
#'   \item \code{"!\{<cmp>\}"}   -> \code{"!(<cmp>)"} (parenthesised so '!'
#'                                   binds to the whole comparison rather
#'                                   than to the left operand alone)
#'   \item \code{"\{<id>\}"}     -> \code{"<id> == 1"}
#'   \item \code{"!\{<id>\}"}    -> \code{"<id> == 0"} (the operationally
#'                                   correct complement of a binary 0/1
#'                                   indicator factor; \code{!factor} is
#'                                   not well-defined in base R and emits
#'                                   the \code{Ops.factor} warning while
#'                                   silently producing an all-NA mask)
#' }
#'
#' Malformed input -- non-character, NA, empty string, missing braces,
#' empty braces, or a bare-name path whose content is not a valid R
#' identifier -- raises a loud \code{stop()}.  Subgroup labels are
#' load-bearing for inferential output, so any deviation from the
#' expected forms must surface as an error rather than degrade silently
#' to a wrong-but-quiet result.
#'
#' @param c_label Character scalar.  A single cut label from a frontier
#'   row's \code{m_cols} entries.
#'
#' @return Character scalar suitable for \code{eval(parse(text = .))}
#'   against the analysis data frame.
#'
#' @keywords internal
#' @noRd
.label_to_expr <- function(c_label) {
  # Hard validation: helper input must be a single non-NA non-empty string.
  if (!is.character(c_label) || length(c_label) != 1L ||
      is.na(c_label) || nchar(c_label) == 0L) {
    stop(".label_to_expr(): 'c_label' must be a single non-empty character ",
         "string; got: ",
         paste(utils::capture.output(utils::str(c_label)), collapse = " "),
         call. = FALSE)
  }

  # Detect leading "!" (complement marker) and extract the inner content
  # of {...}.  The full label must match one of the two recognised forms;
  # anything else is a malformed label and we error loudly rather than
  # silently mangling it.
  negate <- startsWith(c_label, "!")
  if (negate) {
    m <- regmatches(c_label, regexec("^!\\{(.+)\\}$", c_label))[[1L]]
  } else {
    m <- regmatches(c_label, regexec("^\\{(.+)\\}$",  c_label))[[1L]]
  }
  if (length(m) != 2L || nchar(m[2L]) == 0L) {
    stop(".label_to_expr(): unrecognised subgroup label form: '", c_label,
         "'. Expected '{<expr>}' or '!{<expr>}' (FS_labels() output).",
         call. = FALSE)
  }
  inner <- m[2L]

  # Bare column name (no comparison operator) => factor-equality semantics.
  is_bare <- !grepl("[<>=!]", inner)
  if (is_bare) {
    # Inner must look like a single R identifier; reject anything weirder
    # (e.g. embedded spaces or operators we missed) so labels can't quietly
    # become wrong code.
    if (!grepl("^[.A-Za-z_][.A-Za-z0-9_]*$", inner)) {
      stop(".label_to_expr(): bare-name label '", c_label,
           "' does not contain a valid R identifier (got: '", inner, "').",
           call. = FALSE)
    }
    return(if (negate) paste0(inner, " == 0") else paste0(inner, " == 1"))
  }

  # Comparison form: wrap in parens when negating so '!' binds to the
  # whole comparison rather than to the left operand.
  if (negate) paste0("!(", inner, ")") else inner
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
#' @note
#'   \code{eval(parse())} is used intentionally here.
#'   \code{\link{evaluate_comparison}} handles only single comparisons
#'   (e.g., \code{"er <= 0"}); this function is needed for the compound
#'   logical expressions produced by the ForestSearch subgroup enumeration
#'   algorithm (e.g., \code{"er <= 0 & nodes > 3"}).  Evaluation is
#'   sandboxed: the environment contains only the columns of \code{df}
#'   with \code{baseenv()} as parent, so neither the global environment
#'   nor any package namespace is in scope.  No user-supplied strings are
#'   evaluated; only internally-constructed subgroup definition strings
#'   reach this function.
#'
#' @seealso \code{\link{evaluate_comparison}} for the single-comparison
#'   operator-dispatch alternative that avoids \code{eval(parse())}.
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


# =============================================================================
# ARGUMENT FILTERING
# =============================================================================

#' Filter a Parameter List to a Function's Formals
#'
#' Given a named list of user-supplied parameters and a target function,
#' returns only the parameter entries whose names appear in the target's
#' \code{\link[base]{formals}}.  Entries with names that do not match any
#' formal (typos, renamed parameters, parameters from a different function)
#' are either silently dropped (default) or reported via \code{warning()}.
#'
#' This replaces the older "hand-maintained \code{valid_pnames} whitelist"
#' pattern in \code{.run_fs_analysis_gen()} and \code{.run_grf_analysis_gen()}.
#' New \code{forestsearch()} or \code{grf.subg.harm.*()} arguments flow
#' through automatically; authors of new parameters no longer need to
#' update a whitelist in a separate file to make the new parameter
#' reachable from user-level \code{fs_params} / \code{grf_params} lists.
#'
#' @param params Named list of user-supplied parameter values.
#' @param target_fun A function (not its name).  Its \code{\link{formals}}
#'   determine the set of valid parameter names.
#' @param exclude Character vector of formal-argument names to explicitly
#'   \emph{exclude} from the allowlist, even though they appear in
#'   \code{formals(target_fun)}.  Useful for arguments the caller manages
#'   internally (e.g., \code{df.analysis}, \code{confounders.name},
#'   \code{details}, \code{plot.sg}, \code{quiet}) and does not want users
#'   to override via the params list.
#' @param warn_unknown Logical.  If \code{TRUE}, emit a \code{warning()}
#'   listing any \code{params} names that were dropped because they do not
#'   appear in \code{formals(target_fun)} (typo-detection aid).  Default
#'   \code{FALSE} (silent drop, matching legacy behavior).
#'
#' @return A named list, the subset of \code{params} whose names are in
#'   \code{formals(target_fun)} and not in \code{exclude}.
#'
#' @keywords internal
# Typical call sites inside the package:
#
#   # In .run_fs_analysis_gen(), params = user's fs_params list:
#   user_fs_args <- filter_valid_args(
#     params       = params,
#     target_fun   = forestsearch,
#     exclude      = c("df.analysis", "confounders.name",
#                      "details", "plot.sg", "quiet"),
#     warn_unknown = FALSE
#   )
#
# Return value for a params list with a typo:
#   params = list(hr.threshold = 2.0, sg_focus = "maxSG", my_typo = "x")
#   filter_valid_args(params, forestsearch) -> list(hr.threshold = 2.0,
#                                                    sg_focus = "maxSG")
filter_valid_args <- function(params,
                              target_fun,
                              exclude = character(0),
                              warn_unknown = FALSE) {
  if (is.null(params) || length(params) == 0L) return(list())
  if (!is.function(target_fun))
    stop("'target_fun' must be a function.", call. = FALSE)

  allowed <- setdiff(names(formals(target_fun)), exclude)
  # Exclude '...' — passing arbitrary extras via ... creates its own
  # silent-failure risks and is not how fs_params is intended to be used.
  allowed <- setdiff(allowed, "...")

  param_names <- names(params)
  if (is.null(param_names)) return(list())

  keep_mask <- param_names %in% allowed
  dropped   <- param_names[!keep_mask]

  if (warn_unknown && length(dropped) > 0L) {
    warning(sprintf(
      "filter_valid_args: dropped %d unknown parameter(s): %s",
      length(dropped),
      paste(shQuote(dropped), collapse = ", ")
    ), call. = FALSE)
  }

  params[keep_mask]
}


# =============================================================================
# PARALLEL RESULT HELPERS
# =============================================================================

#' Collect Results from foreach with Error Handling
#'
#' Separates successful results from errors in the output of a
#' \code{\link[foreach]{foreach}} loop run with
#' \code{.errorhandling = "pass"}.  Successful results are combined via
#' \code{\link[data.table]{rbindlist}} (with an \code{rbind} fallback
#' if \pkg{data.table} is not available) and failed tasks are counted,
#' reported, and optionally escalated to a hard error.
#'
#' This replaces the fragile pattern of calling \code{nrow()} or
#' \code{rbind()} directly on a raw \code{foreach} result list, which
#' fails silently when error objects are mixed in with data-frame rows.
#'
#' @details
#' Design rationale for the reporting strategy:
#' \itemize{
#'   \item A \code{\link{warning}} is emitted so that calling code
#'     (and \code{tryCatch} / \code{withCallingHandlers} blocks) can
#'     capture the failure condition programmatically.
#'   \item The error summary is also printed with \code{\link{cat}} to
#'     \code{stderr}, because \code{\link{warning}} and
#'     \code{\link{message}} output is sometimes buffered or
#'     suppressed inside future/foreach parallel workers or when
#'     \code{options(warn = -1)} is set by user code upstream.  The
#'     \code{cat()} output is the visible-no-matter-what fallback.
#'   \item When \code{stop_on_all_fail = TRUE} (the default), a
#'     completely empty result list is escalated to an error rather
#'     than silently returning an empty data frame.  This avoids
#'     downstream code proceeding with zero rows and producing
#'     misleading summaries.  Calibration code that may legitimately
#'     tolerate some failures (e.g. Monte Carlo FPR estimation) can
#'     set \code{stop_on_all_fail = FALSE}.
#' }
#'
#' The returned object is always a data frame -- never \code{NULL} --
#' because attributes (\code{n_failed}, \code{n_total}) attached to
#' \code{NULL} are silently discarded by R, whereas attributes on an
#' empty data frame are preserved and inspectable by the caller.
#'
#' @param raw_list List returned by \code{foreach(...,
#'   .errorhandling = "pass")}.  Each element is either a data frame,
#'   data.table, or row-bindable object (success) or a condition
#'   object (failure).
#' @param label Character.  Optional label included in the warning
#'   message to identify which simulation phase produced the errors
#'   (e.g. \code{"H1 alt"}, \code{"Bootstrap"}).  Default \code{""}.
#' @param stop_on_all_fail Logical.  If \code{TRUE} (default),
#'   raise an error when every task in \code{raw_list} failed.  If
#'   \code{FALSE}, return an empty data frame and only emit a
#'   warning.  Default \code{TRUE}.
#' @param max_error_messages Integer.  Maximum number of distinct
#'   error messages to print in the diagnostic output.  Default
#'   \code{3L}.
#' @param keep_diagnostics Logical.  If \code{TRUE}, scan each
#'   successful element for a \code{"diagnostics"} attribute and
#'   collect them into a list attached to the returned object as
#'   \code{attr(result, "diagnostics")}.  Used with
#'   \code{\link{run_simulation_analysis}}'s \code{keep} / 
#'   \code{keep_first_n} mechanism to retain per-replicate heavy
#'   objects (FS or GRF result objects, candidate tables) that
#'   \code{data.table::rbindlist()} would otherwise drop.
#'   Default \code{FALSE} (preserves minimal-memory behavior).
#'
#' @return A data frame of row-bound successful results.  Always a
#'   data frame; \strong{never} \code{NULL} -- an empty data frame is
#'   returned when every task failed and \code{stop_on_all_fail = FALSE}.
#'   Two attributes are always attached:
#'   \describe{
#'     \item{\code{n_failed}}{Integer.
#'       Number of tasks that returned errors.}
#'     \item{\code{n_total}}{Integer.
#'       Total number of tasks (successes + failures).}
#'   }
#'   When \code{keep_diagnostics = TRUE}, a third attribute
#'   \code{"diagnostics"} may be attached -- a named list of
#'   per-replicate diagnostic objects, keyed by \code{"sim<ID>"}
#'   when the source element carries an \code{"sim_id"} attribute
#'   (\code{\link{run_simulation_analysis}} sets this automatically).
#'
#' @section Preserving per-replicate diagnostics:
#' \code{data.table::rbindlist()} drops per-element attributes during
#' row-binding.  When the caller has populated each replicate with
#' heavy diagnostic objects (e.g. via
#' \code{\link{run_simulation_analysis}} with a non-empty
#' \code{keep} argument), set \code{keep_diagnostics = TRUE} to
#' retain them.  The diagnostics attribute on each input element
#' must carry \code{attr(..., "sim_id")} for stable keying; that is
#' set automatically by \code{run_simulation_analysis()}.
#'
#' @examples
#' \donttest{
#' library(foreach)
#' library(doFuture)
#' future::plan("sequential")
#'
#' raw <- foreach(i = 1:5, .errorhandling = "pass") %dofuture% {
#'   if (i == 3L) stop("simulated failure")
#'   data.frame(task = i, value = rnorm(1))
#' }
#' res <- collect_results(raw, label = "demo")
#' # Warning: demo: 1 of 5 tasks failed.
#' nrow(res)
#' # [1] 4
#' attr(res, "n_failed")
#' # [1] 1
#' }
#'
#' @seealso \code{\link{reset_workers}} for resetting the
#'   \pkg{future} plan and releasing worker memory between
#'   simulation phases.
#'
#' @export
collect_results <- function(raw_list,
                            label = "",
                            stop_on_all_fail = TRUE,
                            max_error_messages = 3L,
                            keep_diagnostics = FALSE) {

  if (!is.list(raw_list)) {
    stop("'raw_list' must be a list returned by foreach().", call. = FALSE)
  }
  n_total <- length(raw_list)

  if (n_total == 0L) {
    result <- data.frame()
    attr(result, "n_failed") <- 0L
    attr(result, "n_total")  <- 0L
    return(result)
  }

  is_err   <- vapply(raw_list, inherits, logical(1), "error")
  n_failed <- sum(is_err)

  prefix <- if (nzchar(label)) paste0(label, ": ") else ""

  # --------------------------------------------------------------------
  # Report failures via BOTH cat() (always visible) and warning()
  # (programmatically catchable).  warning() is sometimes buffered or
  # swallowed inside future/foreach workers, so cat() is the belt-and-
  # braces visible fallback.
  # --------------------------------------------------------------------
  if (n_failed > 0L) {
    err_msgs  <- vapply(raw_list[is_err], conditionMessage, character(1))
    uniq_msgs <- utils::head(unique(err_msgs),
                             max(1L, as.integer(max_error_messages)))

    cat(sprintf("*** %s%d of %d tasks failed ***\n",
                prefix, n_failed, n_total),
        file = stderr())
    cat("  First error message(s):\n", file = stderr())
    for (msg in uniq_msgs) cat("    ", msg, "\n", file = stderr())

    warning(sprintf("%s%d of %d tasks failed. First error: %s",
                    prefix, n_failed, n_total, uniq_msgs[1L]),
            call. = FALSE)
  }

  # --------------------------------------------------------------------
  # Harvest per-replicate diagnostics BEFORE rbindlist (which drops
  # element-level attributes).  Only non-NULL "diagnostics" attributes
  # are collected.  See run_simulation_analysis(keep = ..., keep_first_n = ...)
  # for the producer side.
  # --------------------------------------------------------------------
  diag_list <- NULL
  if (isTRUE(keep_diagnostics)) {
    diag_list <- list()
    good_idx  <- which(!is_err)
    for (i in good_idx) {
      d <- attr(raw_list[[i]], "diagnostics")
      if (!is.null(d)) {
        sim_id <- attr(d, "sim_id")
        nm <- if (!is.null(sim_id) && length(sim_id) == 1L)
          sprintf("sim%d", as.integer(sim_id)) else sprintf("i%d", i)
        diag_list[[nm]] <- d
      }
    }
  }

  # --------------------------------------------------------------------
  # Row-bind successful results.  Prefer data.table::rbindlist for
  # fill = TRUE semantics (column-union across heterogeneous rows);
  # fall back to do.call(rbind, ...) if data.table is unavailable.
  # --------------------------------------------------------------------
  good <- raw_list[!is_err]

  if (length(good) == 0L) {
    if (isTRUE(stop_on_all_fail)) {
      stop(sprintf(
        "%sAll %d tasks failed; see warning() for first error.",
        prefix, n_total),
        call. = FALSE)
    }
    result <- data.frame()
  } else if (requireNamespace("data.table", quietly = TRUE)) {
    result <- data.table::rbindlist(good, use.names = TRUE, fill = TRUE)
    result <- as.data.frame(result)
  } else {
    result <- tryCatch(
      do.call(rbind, good),
      error = function(e) {
        stop(sprintf(
          "%sCould not rbind() successful results; install 'data.table' ",
          "for heterogeneous-column support. Original error: %s",
          prefix, conditionMessage(e)),
          call. = FALSE)
      }
    )
    result <- as.data.frame(result)
  }

  attr(result, "n_failed") <- as.integer(n_failed)
  attr(result, "n_total")  <- as.integer(n_total)
  if (isTRUE(keep_diagnostics) && length(diag_list) > 0L) {
    attr(result, "diagnostics") <- diag_list
  }
  result
}


#' Reset a Parallel Worker Pool
#'
#' Tears down the current \pkg{future} plan, forces garbage
#' collection, optionally pauses, and spins up a fresh worker pool.
#' Used to release worker-side memory between simulation phases
#' (e.g. between alternative and null Monte Carlo runs) so that
#' accumulated R-session state in long-running workers does not
#' bloat memory or cause sporadic failures.
#'
#' @details
#' The worker reset sequence is:
#' \enumerate{
#'   \item \code{future::plan("sequential")} -- dismantles the
#'     existing worker pool; R processes terminate and their memory
#'     is returned to the OS.
#'   \item \code{gc(verbose = FALSE)} -- reclaims main-process memory
#'     holding references to now-defunct futures.
#'   \item \code{Sys.sleep(sleep)} -- brief pause to let the OS
#'     reclaim sockets / pipes before the new pool is requested.
#'     A small non-zero sleep (the default \code{1}) materially
#'     reduces \code{multisession} setup failures observed on busy
#'     Linux servers; set to \code{0} to skip.
#'   \item \code{future::plan(strategy, workers = workers)} -- starts
#'     a fresh pool.
#' }
#'
#' This function is a no-op if \code{future} is not installed or if
#' the requested strategy cannot be initialised; a warning is emitted
#' and the previous plan is left in place.
#'
#' @param workers Integer.  Number of workers for the new pool.
#'   Default \code{NULL} leaves the current plan's worker count
#'   unchanged (useful if a reset without resize is wanted; this
#'   requires \code{future} to be installed).
#' @param strategy Character.  The \pkg{future} strategy to
#'   reinstate.  Default \code{"multisession"}.  Other reasonable
#'   choices are \code{"multicore"} (Linux/macOS only) and
#'   \code{"sequential"} (to disable parallelism entirely).
#' @param sleep Numeric.  Seconds to pause between teardown and
#'   setup.  Default \code{1}.  Set to \code{0} to skip.
#' @param gc Logical.  Whether to force garbage collection between
#'   teardown and setup.  Default \code{TRUE}.
#' @param quiet Logical.  Suppress informational \code{message()}
#'   output.  Default \code{FALSE}.
#'
#' @return Invisibly, a list with the previous and current plan
#'   class names: \code{list(previous = <chr>, current = <chr>)}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("future", quietly = TRUE)) {
#'   future::plan("multisession", workers = 2L)
#'   # ... run alternative simulations ...
#'   reset_workers(workers = 2L)
#'   # ... run null simulations ...
#'   future::plan("sequential")
#' }
#' }
#'
#' @seealso \code{\link{collect_results}} for error-tolerant
#'   collection of \code{foreach} output.
#'
#' @export
reset_workers <- function(workers   = NULL,
                          strategy  = "multisession",
                          sleep     = 1,
                          gc        = TRUE,
                          quiet     = FALSE) {

  if (!requireNamespace("future", quietly = TRUE)) {
    warning("'future' is not installed; reset_workers() is a no-op.",
            call. = FALSE)
    return(invisible(list(previous = NA_character_, current = NA_character_)))
  }

  prev_class <- class(future::plan())[1L]

  # -- Teardown ---------------------------------------------------------
  tryCatch(
    future::plan("sequential"),
    error = function(e) {
      warning(sprintf("Could not tear down existing plan: %s",
                      conditionMessage(e)),
              call. = FALSE)
    }
  )

  if (isTRUE(gc)) {
    gc(verbose = FALSE)
  }

  if (is.numeric(sleep) && length(sleep) == 1L && sleep > 0) {
    Sys.sleep(sleep)
  }

  # -- Setup ------------------------------------------------------------
  # If the caller wants sequential, we're already there -- skip the setup.
  if (identical(strategy, "sequential")) {
    if (!isTRUE(quiet)) {
      message(sprintf("reset_workers: plan reset to sequential (was %s).",
                      prev_class))
    }
    return(invisible(list(previous = prev_class, current = "sequential")))
  }

  new_plan <- tryCatch({
    if (is.null(workers)) {
      future::plan(strategy)
    } else {
      future::plan(strategy, workers = as.integer(workers))
    }
    class(future::plan())[1L]
  }, error = function(e) {
    warning(sprintf("Could not start '%s' plan: %s. Leaving plan as sequential.",
                    strategy, conditionMessage(e)),
            call. = FALSE)
    "sequential"
  })

  if (!isTRUE(quiet)) {
    if (is.null(workers)) {
      message(sprintf("reset_workers: plan reset to %s (was %s).",
                      new_plan, prev_class))
    } else {
      message(sprintf("reset_workers: plan reset to %s with %d worker(s) (was %s).",
                      new_plan, as.integer(workers), prev_class))
    }
  }

  invisible(list(previous = prev_class, current = new_plan))
}


# ============================================================================
# .normalize_sg_focus() -- internal vocabulary alias
# ============================================================================
# The package's `sg_focus` argument originally used a Cox / hazard-ratio
# vocabulary ("hr", "hrMaxSG", "hrMinSG") that pre-dates the GLM extension
# (continuous MD, binary OR/RR/RD, count IRR).  The GLM-natural vocabulary
# uses "eff*" prefixes -- "effMaxSG", "effMinSG", "eff" -- which read
# correctly regardless of effect measure.
#
# Both vocabularies are accepted throughout the user-facing API; this
# helper translates the new names to the canonical internal "hr*" form
# at entry points so that downstream code (which is keyed on the
# canonical form) needs no changes.
#
# Aliases (case-sensitive, matches `match.arg` conventions):
#   "effMaxSG" -> "hrMaxSG"
#   "effMinSG" -> "hrMinSG"
#   "eff"      -> "hr"
#   "maxSG"    -> "maxSG"   (unchanged; pre-existing alias for hrMaxSG)
#   "minSG"    -> "minSG"   (unchanged; pre-existing alias for hrMinSG)
# Any other value is returned unchanged so the downstream whitelist
# check produces its usual error message.
# ----------------------------------------------------------------------------

#' Translate the new effect-vocabulary `sg_focus` aliases to canonical form
#'
#' Accepts the GLM-natural vocabulary (\code{"effMaxSG"},
#' \code{"effMinSG"}, \code{"eff"}) and translates to the canonical
#' Cox-flavored internal names (\code{"hrMaxSG"}, \code{"hrMinSG"},
#' \code{"hr"}).  Unrecognized values pass through unchanged so the
#' caller's own whitelist check fires.
#'
#' This is an internal helper.  User code does not call it directly;
#' \code{\link{forestsearch}} and \code{\link{subgroup.consistency}}
#' invoke it on entry.
#'
#' @param sg_focus Character scalar.
#'
#' @return Character scalar in canonical form.
#'
#' @keywords internal
#' @noRd
.normalize_sg_focus <- function(sg_focus) {
  if (!is.character(sg_focus) || length(sg_focus) != 1L) {
    return(sg_focus)  # downstream whitelist will reject
  }
  switch(sg_focus,
    effMaxSG = "hrMaxSG",
    effMinSG = "hrMinSG",
    eff      = "hr",
    sg_focus
  )
}


# =============================================================================
# DINA screening-stage helpers (Phase 2: use_dina)
# =============================================================================

#' Map a forestsearch outcome_type to a DINA family
#'
#' @param outcome_type One of "survival", "binary", "continuous", "count".
#' @return The corresponding DINA family string.
#' @noRd
.map_dina_family <- function(outcome_type) {
  switch(outcome_type,
    survival   = "cox",
    binary     = "binomial",
    continuous = "gaussian",
    count      = "poisson",
    stop(".map_dina_family(): unrecognised outcome_type '", outcome_type, "'.")
  )
}

#' Resolve and validate the \code{dina_args} list for \code{forestsearch()}
#'
#' Fills defaults (several inherited from the enclosing \code{forestsearch()}
#' call) and rejects unknown keys with a loud error.  Splits the result into
#' the keys destined for \code{dina()} (the fit) and those destined for
#' \code{dina_frontier()} (the extraction).
#'
#' @param dina_args User-supplied list (possibly empty).
#' @param outcome_type forestsearch outcome_type, for the family default.
#' @param n_min_default forestsearch \code{n.min}, the default frontier
#'   \code{n_min}.
#' @param seed_default forestsearch \code{seedit}, the default DINA seed.
#' @return A list with elements \code{fit} (named list of \code{dina()}
#'   arguments) and \code{frontier} (named list of \code{dina_frontier()}
#'   arguments).
#' @noRd
.resolve_dina_args <- function(dina_args, outcome_type, n_min_default,
                               seed_default) {
  if (is.null(dina_args)) dina_args <- list()
  if (!is.list(dina_args)) {
    stop("`dina_args` must be a list.", call. = FALSE)
  }
  fit_keys      <- c("family", "seed", "n_folds", "cens_type", "cens_params")
  frontier_keys <- c("scope", "m_diff", "n_min", "direction",
                      "max_per_covariate", "max_subgroups", "digits")
  # screening-behavior keys: control HOW use_dina screening turns the DINA
  # fit into candidate cuts (frontier vs. single selected cut) and how deep
  # the selected-cut search runs.  Not passed to dina() or dina_frontier();
  # max_depth / grid_probs are forwarded to dina_subgroup() at the selection
  # sites, selected_only is consumed directly by forestsearch().
  behavior_keys <- c("selected_only", "max_depth", "grid_probs",
                     "select_statistic")
  recognised    <- c(fit_keys, frontier_keys, behavior_keys)

  nms <- names(dina_args)
  if (length(dina_args) > 0L && (is.null(nms) || any(nms == ""))) {
    stop("`dina_args` must be a fully named list.", call. = FALSE)
  }
  unknown <- setdiff(nms, recognised)
  if (length(unknown) > 0L) {
    stop("`dina_args` has unrecognised key(s): ",
         paste(shQuote(unknown), collapse = ", "),
         ".  Recognised keys: ",
         paste(shQuote(recognised), collapse = ", "), ".",
         call. = FALSE)
  }

  get_arg <- function(key, default) if (key %in% nms) dina_args[[key]] else default

  # selected_only: when TRUE (the default, matching GRF's
  # return_selected_cuts_only), use_dina screening contributes the SINGLE cut
  # dina_subgroup() selects (using forestsearch's own sg_focus /
  # selection_rule / etc.); set FALSE to contribute the full dina_frontier()
  # candidate set instead.
  selected_only <- get_arg("selected_only", TRUE)
  if (!is.logical(selected_only) || length(selected_only) != 1L ||
      is.na(selected_only)) {
    stop("`dina_args$selected_only` must be a single TRUE/FALSE.",
         call. = FALSE)
  }

  # max_depth / grid_probs: forwarded to dina_subgroup() at the selected-cut
  # and subgroup_method = "dina" sites.  Default max_depth = 2L aligns DINA with
  # the consistency search (maxk = 2) and GRF (grf_depth = 2): the selected cut
  # may be an AND-conjunction of two covariates, whose component cuts are then
  # composed by the consistency search as usual.  Set 1L for depth-1
  # (single-covariate) candidates only.
  max_depth <- get_arg("max_depth", 2L)
  max_depth <- suppressWarnings(as.integer(max_depth))
  if (length(max_depth) != 1L || is.na(max_depth) ||
      !max_depth %in% c(1L, 2L)) {
    stop("`dina_args$max_depth` must be 1 or 2.", call. = FALSE)
  }
  grid_probs <- get_arg("grid_probs", seq(0.1, 0.9, by = 0.1))
  if (!is.numeric(grid_probs) || length(grid_probs) < 1L ||
      anyNA(grid_probs) || any(grid_probs < 0) || any(grid_probs > 1)) {
    stop("`dina_args$grid_probs` must be a numeric vector in [0, 1] ",
         "with no NAs.", call. = FALSE)
  }

  # select_statistic: which statistic ranks the qualifying candidate family.
  # "dina" (default) ranks on DINA's native subgroup-mean tau-hat -- the
  # legacy behaviour, unchanged.  "effect" ranks on the inferential effect
  # measure (Cox HR for survival; OR/MD/IRR under the GLM extension), computed
  # with the SAME per-candidate estimator the Tier-2 de-biased gate uses, so
  # the realized selection is the exact event the gate then de-biases.
  select_statistic <- get_arg("select_statistic", "dina")
  if (!is.character(select_statistic) || length(select_statistic) != 1L ||
      !select_statistic %in% c("dina", "effect")) {
    stop("`dina_args$select_statistic` must be \"dina\" or \"effect\".",
         call. = FALSE)
  }

  # Fit arguments: family + seed always set; the remaining fit keys are
  # forwarded to dina() ONLY when the user supplied them, so dina()'s own
  # defaults otherwise apply (no silently invented values).
  fit <- list(
    family = get_arg("family", .map_dina_family(outcome_type)),
    seed   = get_arg("seed", seed_default)
  )
  for (k in c("n_folds", "cens_type", "cens_params")) {
    if (k %in% nms) fit[[k]] <- dina_args[[k]]
  }

  # Frontier arguments: explicit defaults, n_min inherits forestsearch n.min.
  frontier <- list(
    scope             = get_arg("scope", "wide"),
    m_diff            = get_arg("m_diff", NULL),
    n_min             = get_arg("n_min", n_min_default),
    direction         = get_arg("direction", "both"),
    max_per_covariate = get_arg("max_per_covariate", 3L),
    max_subgroups     = get_arg("max_subgroups", 10L),
    digits            = get_arg("digits", 3L)
  )

  list(fit = fit, frontier = frontier, selected_only = selected_only,
       select = list(max_depth = max_depth, grid_probs = grid_probs,
                     select_statistic = select_statistic))
}


#' Re-rank DINA's qualifying family on the inferential effect measure
#'
#' Replaces DINA's native tau-hat winner with the candidate that maximises the
#' effect measure the Tier-2 gate de-biases (Cox HR for survival; the resolved
#' GLM effect otherwise), scored with the gate's own per-candidate estimator
#' (`.fs_dg_pieces`) so the realized selection is the exact event the gate
#' corrects.  Ranking honours the same `sg_focus` / `selection_rule` /
#' `effect_neighborhood` band logic as `dina_subgroup()`.  Returns `sg` with the
#' winner fields overwritten and a `sel_effect` (link-scale) column attached to
#' `sg$candidates`; if no candidate is scorable it returns `sg` unchanged.
#'
#' @keywords internal
#' @noRd
.dina_reselect_on_effect <- function(sg, df, outcome_type, effect_measure,
                                     treat.name, outcome.name, event.name,
                                     offset.name, adjust_covariates,
                                     adverse_outcome, sg_focus,
                                     selection_rule, effect_neighborhood) {
  tab <- sg$candidates
  if (is.null(tab) || !nrow(tab)) return(sg)

  # Effect measure / spec, mirroring the DINA Tier-2 gate branch: survival uses
  # "HR" on the log scale; GLM uses the resolved effect_measure.
  em   <- if (identical(outcome_type, "survival")) "HR" else effect_measure
  adv  <- if (identical(outcome_type, "survival")) TRUE else adverse_outcome
  spec <- .fs_dg_spec(outcome_type, em, treat.name, outcome.name, event.name,
                      offset.name, adjust_covariates, adverse_outcome = adv,
                      df = df)

  # Per-candidate effect on the SAME statistic the gate uses (.fs_dg_pieces),
  # row-aligned to `tab`; skip rows with < 6 members exactly as the gate does.
  nr        <- nrow(tab)
  eff_link  <- rep(NA_real_, nr)
  sz        <- rep(NA_integer_, nr)
  log_scale <- TRUE
  for (i in seq_len(nr)) {
    op1   <- if (identical(as.character(tab$d1[i]), "left")) "<=" else ">="
    comps <- list(data.frame(variable = as.character(tab$v1[i]), op = op1,
                             value = as.numeric(tab$c1[i]),
                             stringsAsFactors = FALSE))
    if (!is.na(tab$v2[i])) {
      op2 <- if (identical(as.character(tab$d2[i]), "left")) "<=" else ">="
      comps[[2L]] <- data.frame(variable = as.character(tab$v2[i]), op = op2,
                                value = as.numeric(tab$c2[i]),
                                stringsAsFactors = FALSE)
    }
    cj  <- do.call(rbind, comps)
    mem <- tryCatch(.fs_dg_members_from_conj(df, cj),
                    error = function(e) integer(0))
    if (length(mem) < 6L) next
    pc <- tryCatch(.fs_dg_pieces(df[mem, , drop = FALSE], spec),
                   error = function(e) NULL)
    if (is.null(pc) || !is.finite(pc$beta_hat)) next
    eff_link[i] <- pc$beta_hat
    sz[i]       <- length(mem)
    log_scale   <- isTRUE(pc$log_scale)
  }
  sg$candidates$sel_effect <- eff_link  # link scale; consumed by the gate branch

  ok <- which(is.finite(eff_link) & !is.na(sz))
  if (!length(ok)) return(sg)           # nothing scorable -> keep native winner

  # Natural-scale effect for ordering / band (ratio families exponentiated; the
  # monotone transform leaves non-band foci unaffected).  Match dina_subgroup()'s
  # canonical sg_focus and its order() switch, with row index as a stable
  # insertion-order tiebreak.
  sgf <- .normalize_sg_focus(sg_focus)
  eff <- if (log_scale) exp(eff_link) else eff_link
  ord <- switch(
    sgf,
    maxSG   = ok[order(-sz[ok], -eff[ok], ok)],
    minSG   = ok[order( sz[ok], -eff[ok], ok)],
    hr      = ok[order(-eff[ok], ok)],
    hrMaxSG = {
      band <- .compute_inclusion_band(hr_vec = eff[ok], n_vec = sz[ok],
                                      selection_rule = selection_rule,
                                      effect_neighborhood = effect_neighborhood)
      ok[order(-band, -sz[ok], -eff[ok], ok)]
    },
    hrMinSG = {
      band <- .compute_inclusion_band(hr_vec = eff[ok], n_vec = sz[ok],
                                      selection_rule = selection_rule,
                                      effect_neighborhood = effect_neighborhood)
      ok[order(-band, sz[ok], -eff[ok], ok)]
    },
    ok[order(-eff[ok], ok)]             # fallback: plain effect argmax
  )
  w <- ord[1L]

  # Overwrite the winner fields from the winning candidate row, rebuild the AND
  # mask so the stored out_sg is internally consistent.  The caller builds
  # sg.harm / membership from sg$covariate/direction/threshold.
  win_v <- c(as.character(tab$v1[w]), as.character(tab$v2[w]))
  win_d <- c(as.character(tab$d1[w]), as.character(tab$d2[w]))
  win_c <- c(as.numeric(tab$c1[w]),  as.numeric(tab$c2[w]))
  keep  <- !is.na(win_v)
  sg$covariate  <- win_v[keep]
  sg$direction  <- win_d[keep]
  sg$threshold  <- win_c[keep]
  sg$depth      <- sum(keep)
  sg$n_subgroup <- sz[w]
  if ("tau_hat" %in% names(tab)) sg$mean_tau_hat <- as.numeric(tab$tau_hat[w])

  mask <- rep(TRUE, nrow(df))
  for (t in seq_len(sg$depth)) {
    x_t <- df[[sg$covariate[t]]]
    mask <- mask & (if (identical(sg$direction[t], "left"))
      x_t <= sg$threshold[t] else x_t >= sg$threshold[t])
  }
  sg$mask <- mask
  sg$select_statistic <- "effect"
  sg
}


#' DINA-selection mode for forestsearch (subgroup_method = "dina")
#'
#' Fits a DINA model and delegates subgroup selection to
#' \code{dina_subgroup()}, bypassing GRF/LASSO/consistency.  Returns the
#' selected subgroup as a forestsearch cut label plus the membership
#' (\code{treat.recommend}) tables the estimation/bootstrap machinery
#' consumes.  Selection criteria come from the enclosing forestsearch()
#' call; \code{dina_args} supplies only the DINA fit tuning.
#'
#' @return A list with \code{found}, \code{sg.harm} (label vector or NULL),
#'   \code{grp.consistency} (a consistency-shaped list), \code{dina_res},
#'   and \code{df.est}/\code{df.predict}/\code{df.test} carrying
#'   \code{treat.recommend} (0 = selected subgroup).
#' @noRd
.forestsearch_dina_select <- function(df, df.predict, df.test,
                                      confounders.name, outcome.name, event.name,
                                      treat.name, id.name, outcome_type,
                                      hr.threshold, n.min, sg_focus,
                                      selection_rule, effect_neighborhood,
                                      dina_args, dina_res, seedit, details,
                                      effect_measure = NULL, offset.name = NULL,
                                      adjust_covariates = NULL,
                                      adverse_outcome = TRUE) {
  da <- .resolve_dina_args(dina_args, outcome_type,
                           n_min_default = n.min, seed_default = seedit)

  # Fit DINA unless a fit was supplied.
  if (is.null(dina_res)) {
    status_arg <- if (identical(da$fit$family, "cox")) event.name else NULL
    fit_call <- c(list(df = df, outcome = outcome.name, treatment = treat.name,
                       covariates = confounders.name, status = status_arg),
                  da$fit)
    dina_res <- do.call(dina, fit_call)
  }

  # Harm floor on the link scale: log(hr.threshold) for ratio families
  # (cox/binomial/poisson), identity (mean difference) for gaussian.
  m_diff <- if (identical(da$fit$family, "gaussian")) hr.threshold
            else log(hr.threshold)

  if (isTRUE(details)) {
    lines <- c(
      "",
      "[forestsearch] DINA selection (subgroup_method = \"dina\")",
      paste0("  Family:              ", da$fit$family),
      paste0("  sg_focus:            ", sg_focus),
      paste0("  selection_rule:      ", selection_rule),
      paste0("  effect_neighborhood: ", effect_neighborhood),
      paste0("  Harm floor:          ", sprintf("m_diff = %.4f", m_diff),
             if (!identical(da$fit$family, "gaussian"))
               sprintf("  (hr.threshold = %.4g)", hr.threshold) else ""),
      paste0("  n.min:               ", n.min)
    )

    # DINA frontier: the per-covariate candidate cuts under consideration.
    fr <- tryCatch(
      dina_frontier(fit = dina_res, df = df, covariates = confounders.name,
                    scope = "wide", n_min = n.min),
      error = function(e) {
        lines <<- c(lines, paste0("  (frontier unavailable: ",
                                  conditionMessage(e), ")"))
        NULL
      }
    )
    if (!is.null(fr) && nrow(fr) > 0L) {
      lines <- c(lines, "  DINA frontier candidates (per-covariate non-dominated):")
      fr_show <- fr[, intersect(c("covariate", "direction", "threshold",
                                  "n_subgroup", "effect", "cut_expr"),
                                names(fr)), drop = FALSE]
      lines <- c(lines,
                 paste0("    ",
                        utils::capture.output(print(fr_show, row.names = FALSE))))
    } else if (!is.null(fr)) {
      lines <- c(lines, "  DINA frontier: no candidates met the size constraint.")
    }
    # Single message() call: stderr is forwarded by doFuture workers, whereas
    # cat() (stdout) is captured by the future and never reaches the document.
    message(paste(lines, collapse = "\n"))
  }

  sg <- dina_subgroup(
    fit                 = dina_res,
    df                  = df,
    covariates          = confounders.name,
    m_diff              = m_diff,
    n_min               = n.min,
    max_depth           = da$select$max_depth,
    grid_probs          = da$select$grid_probs,
    sg_focus            = sg_focus,
    selection_rule      = selection_rule,
    effect_neighborhood = effect_neighborhood
  )

  # Effect-based re-selection (dina_args$select_statistic = "effect"): re-rank
  # DINA's qualifying family on the inferential effect measure the Tier-2 gate
  # de-biases, overriding the native tau-hat winner.  The default "dina" path
  # leaves `sg` untouched.  Any failure falls back to the native winner.
  if (identical(da$select$select_statistic, "effect") && isTRUE(sg$found) &&
      !is.null(sg$candidates) && nrow(sg$candidates) > 0L) {
    sg <- tryCatch(
      .dina_reselect_on_effect(
        sg = sg, df = df, outcome_type = outcome_type,
        effect_measure = effect_measure, treat.name = treat.name,
        outcome.name = outcome.name, event.name = event.name,
        offset.name = offset.name, adjust_covariates = adjust_covariates,
        adverse_outcome = adverse_outcome, sg_focus = sg_focus,
        selection_rule = selection_rule,
        effect_neighborhood = effect_neighborhood),
      error = function(e) sg)
  }

  if (isTRUE(details)) {
    lines2 <- c(
      paste0("  Candidates searched:  ",
             if (!is.null(sg$n_candidates_searched)) sg$n_candidates_searched else NA),
      paste0("  Candidates qualifying (>= floor, >= n.min): ",
             if (!is.null(sg$n_candidates_qualifying)) sg$n_candidates_qualifying else NA)
    )
    if (isTRUE(sg$found)) {
      op_d    <- ifelse(sg$direction == "left", "<=", ">=")
      sel_lbl <- paste(sprintf("{%s %s %s}", sg$covariate, op_d,
                               format(sg$threshold)),
                       collapse = " & ")
      lines2 <- c(lines2,
                  paste0("  SELECTED: ", sel_lbl,
                         sprintf("  (n = %d, mean tau-hat = %.4f)",
                                 sg$n_subgroup, sg$mean_tau_hat)))
    } else {
      lines2 <- c(lines2,
                  "  SELECTED: none -- no candidate met the harm floor and size constraint.")
    }
    message(paste(lines2, collapse = "\n"))
  }

  if (!isTRUE(sg$found)) {
    return(list(found = FALSE, sg.harm = NULL, grp.consistency = NULL,
                dina_res = dina_res, df.est = df,
                df.predict = df.predict, df.test = df.test))
  }

  # Cut label(s): direction "left" => {x <= q}, "right" => {x >= q}.  For a
  # depth-2 conjunction sg$covariate/direction/threshold are length-2, so
  # sg.harm becomes a two-element label vector that get_dfpred() AND-composes
  # -- the same representation a 2-factor consistency subgroup uses.
  op      <- ifelse(sg$direction == "left", "<=", ">=")
  sg.harm <- sprintf("{%s %s %s}", sg$covariate, op,
                     formatC(sg$threshold, format = "g", digits = 17, width = 1L))

  # Membership via the standard label evaluator (treat.recommend == 0 = harm).
  df.est         <- get_dfpred(df, sg.harm, version = 2)
  df.predict_out <- if (!is.null(df.predict)) get_dfpred(df.predict, sg.harm, version = 2) else NULL
  df.test_out    <- if (!is.null(df.test))    get_dfpred(df.test,    sg.harm, version = 2) else NULL

  id_vals <- if (id.name %in% names(df)) df[[id.name]]
             else if ("id" %in% names(df)) df[["id"]]
             else seq_len(nrow(df))
  df_flag <- data.frame(id = id_vals,
                        treat.recommend = df.est$treat.recommend)

  grp.consistency <- list(
    out_sg     = sg,
    sg.harm    = sg.harm,
    sg.harm.id = as.integer(df.est$treat.recommend == 0L),
    df_flag    = df_flag,
    algorithm  = "dina"
  )

  list(found = TRUE, sg.harm = sg.harm, grp.consistency = grp.consistency,
       dina_res = dina_res, df.est = df.est,
       df.predict = df.predict_out, df.test = df.test_out,
       candidates = sg$candidates,    # qualifying family for the Tier-2 gate
       select_statistic = da$select$select_statistic)
}


#' Re-rank GRF's DR-candidate frontier on the inferential effect measure
#'
#' The GRF analogue of \code{.dina_reselect_on_effect()}, for
#' \code{grf_selection = "frontier"}.  Re-scores the DR-candidate family
#' (\code{grf_res$candidates}) on the effect measure the Tier-2 gate de-biases
#' (Cox HR for survival; the resolved GLM effect otherwise), using the gate's
#' own per-candidate estimator (\code{.fs_dg_pieces}), then re-selects the winner
#' with GRF's own \code{.grf_frontier_select()} logic on those Cox-HR scores
#' (harm floor HR >= 1).  The winning row is turned back into the standard
#' \code{sg_def} via \code{.grf_sg_def_from_candidate()}, so all downstream
#' membership/labelling is unchanged.  Attaches a \code{sel_effect} (link-scale)
#' column to \code{grf_res$candidates}; on any failure to score/select it returns
#' \code{grf_res} unchanged.  Tree-mode selection has no enumerated family to
#' rank, so the caller applies this only in frontier mode.
#'
#' @keywords internal
#' @noRd
.grf_reselect_on_effect <- function(grf_res, df, outcome_type, effect_measure,
                                    treat.name, outcome.name, event.name,
                                    offset.name, adjust_covariates,
                                    adverse_outcome, frontier_rule,
                                    effect_neighborhood) {
  cand <- grf_res$candidates
  if (is.null(cand) || !nrow(cand)) return(grf_res)

  em   <- if (identical(outcome_type, "survival")) "HR" else effect_measure
  adv  <- if (identical(outcome_type, "survival")) TRUE else adverse_outcome
  spec <- .fs_dg_spec(outcome_type, em, treat.name, outcome.name, event.name,
                      offset.name, adjust_covariates, adverse_outcome = adv,
                      df = df)

  # Per-candidate effect on the SAME statistic the gate uses, with GRF-native
  # membership (sg_def -> .grf_evaluate_subgroup); skip < 6-member candidates.
  nr        <- nrow(cand)
  eff_link  <- rep(NA_real_, nr)
  log_scale <- TRUE
  for (i in seq_len(nr)) {
    sgd_i <- tryCatch(.grf_sg_def_from_candidate(cand[i, , drop = FALSE]),
                      error = function(e) NULL)
    if (is.null(sgd_i)) next
    mem <- tryCatch(which(.grf_evaluate_subgroup(sgd_i, df) == 0L),
                    error = function(e) integer(0))
    if (length(mem) < 6L) next
    pc <- tryCatch(.fs_dg_pieces(df[mem, , drop = FALSE], spec),
                   error = function(e) NULL)
    if (is.null(pc) || !is.finite(pc$beta_hat)) next
    eff_link[i] <- pc$beta_hat
    log_scale   <- isTRUE(pc$log_scale)
  }
  grf_res$candidates$sel_effect <- eff_link  # link scale; consumed by gate branch
  if (!any(is.finite(eff_link))) return(grf_res)  # nothing scorable -> keep native

  # Re-select via GRF's own frontier logic, scoring on the natural-scale effect
  # with an HR-harm floor (>= 1) -- the ratio analogue of dmin.grf on the
  # additive DR scale.
  cand_hr <- grf_res$candidates
  cand_hr$effect <- if (log_scale) exp(eff_link) else eff_link
  cand_hr <- cand_hr[is.finite(cand_hr$effect), , drop = FALSE]
  win <- tryCatch(
    .grf_frontier_select(cand_hr, dmin = 1, rule = frontier_rule,
                         nbhd = effect_neighborhood),
    error = function(e) NULL)
  if (is.null(win) || !nrow(win)) return(grf_res)  # no HR-harm winner -> keep native

  grf_res$sg_def <- .grf_sg_def_from_candidate(win)
  grf_res$select_statistic <- "effect"
  grf_res
}


#' GRF-selection mode for forestsearch (subgroup_method = "grf")
#'
#' Delegates subgroup *selection* to the GRF subgroup-identification routine
#' (\code{grf.subg.harm.survival()} for survival, \code{grf.subg.harm.glm()}
#' for GLM outcomes), bypassing the LASSO/DINA screening and the consistency
#' search.  This is the GRF analogue of \code{\link{.forestsearch_dina_select}}:
#' the policy-tree harm leaf is returned as a forestsearch cut-label vector
#' plus the membership (\code{treat.recommend}) tables that the
#' estimation/bootstrap machinery consumes.
#'
#' The GRF run is configured from the enclosing \code{forestsearch()} call's
#' existing GRF arguments (\code{frac.tau}, \code{dmin.grf}, \code{grf_depth},
#' \code{is.RCT}, \code{seedit}, and the GLM extras), built through the same
#' \code{.build_grf_survival_args()} / \code{.build_grf_glm_args()} helpers the
#' screening path uses, so a GRF-selection run is identical to the GRF screening
#' fit -- only the consumption differs.  \code{return_selected_cuts_only} is
#' forced to \code{TRUE} here: the selected harm leaf is the subgroup, so the
#' single root-to-leaf cut set (not the full per-depth cut pool) is what defines
#' it.
#'
#' GRF emits each split as a \code{"var <= value"} expression, and a depth-d
#' harm leaf is the AND-intersection of its d ancestor splits, so
#' \code{grf_res$sg.harm.id} is a length-d character vector that maps directly
#' onto a d-element forestsearch label vector (wrapped in braces for
#' \code{get_dfpred()}).  Membership is then re-derived from those labels via
#' \code{get_dfpred()}, exactly as the DINA path does, which makes it apply
#' cleanly to \code{df.predict}/\code{df.test} and is verified to reproduce
#' GRF's own tree-node assignment.
#'
#' @return A list with \code{found}, \code{sg.harm} (label vector or NULL),
#'   \code{grp.consistency} (a consistency-shaped list), \code{grf_res} (the
#'   raw GRF result object), and \code{df.est}/\code{df.predict}/\code{df.test}
#'   carrying \code{treat.recommend} (0 = selected subgroup).
#' @noRd
.forestsearch_grf_select <- function(df, df.predict, df.test,
                                     confounders.name, outcome.name,
                                     event.name, treat.name, id.name,
                                     outcome_type, n.min,
                                     frac.tau, dmin.grf, grf_depth, is.RCT,
                                     adverse_outcome = FALSE,
                                     offset.name = NULL,
                                     overdispersion = "none",
                                     grf_count_transform = "log",
                                     grf_res = NULL, seedit = 8316951L,
                                     grf_selection = "tree",
                                     frontier_rule = "effMaxSG",
                                     effect_neighborhood = 0.10,
                                     details = FALSE,
                                     grf_select_statistic = "dr",
                                     effect_measure = NULL,
                                     adjust_covariates = NULL) {

  # Fit GRF unless a fit was supplied.  Build args through the same helpers the
  # screening path uses so the GRF-selection fit cannot drift from GRF
  # screening.  return_selected_cuts_only = TRUE: the harm leaf IS the subgroup.
  if (is.null(grf_res)) {
    if (identical(outcome_type, "survival")) {
      surv_args <- .build_grf_survival_args(
        data                      = df,
        confounders.name          = confounders.name,
        outcome.name              = outcome.name,
        event.name                = event.name,
        id.name                   = id.name,
        treat.name                = treat.name,
        frac.tau                  = frac.tau,
        n.min                     = n.min,
        dmin.grf                  = dmin.grf,
        is.RCT                    = is.RCT,
        grf_depth                 = grf_depth,
        seedit                    = seedit,
        return_selected_cuts_only = TRUE,
        grf_selection             = grf_selection,
        frontier_rule             = frontier_rule,
        effect_neighborhood       = effect_neighborhood
      )
      grf_res <- do.call(grf.subg.harm.survival, surv_args)
    } else {
      glm_args <- .build_grf_glm_args(
        data                      = df,
        confounders.name          = confounders.name,
        outcome.name              = outcome.name,
        treat.name                = treat.name,
        id.name                   = id.name,
        outcome_type              = outcome_type,
        n.min                     = n.min,
        dmin.grf                  = dmin.grf,
        is.RCT                    = is.RCT,
        grf_depth                 = grf_depth,
        seedit                    = seedit,
        return_selected_cuts_only = TRUE,
        adverse_outcome           = adverse_outcome,
        offset.name               = offset.name,
        overdispersion            = overdispersion,
        grf_count_transform       = grf_count_transform,
        grf_selection             = grf_selection,
        frontier_rule             = frontier_rule,
        effect_neighborhood       = effect_neighborhood
      )
      grf_res <- do.call(grf.subg.harm.glm, glm_args)
    }
  }

  # Effect-based re-selection (grf_select_statistic = "effect"): re-rank GRF's
  # DR-candidate frontier on the inferential effect the Tier-2 gate de-biases,
  # overriding the native DR-score winner.  Frontier-only: tree-mode selection
  # is the policy-tree leaf, which has no enumerated family to rank, so the leaf
  # stands.  Default "dr" leaves grf_res untouched.  Failures fall back to native.
  if (identical(grf_select_statistic, "effect") &&
      identical(grf_selection, "frontier") &&
      !is.null(grf_res) && !inherits(grf_res, "try-error") &&
      !is.null(grf_res$candidates) && nrow(grf_res$candidates) > 0L) {
    grf_res <- tryCatch(
      .grf_reselect_on_effect(
        grf_res = grf_res, df = df, outcome_type = outcome_type,
        effect_measure = effect_measure, treat.name = treat.name,
        outcome.name = outcome.name, event.name = event.name,
        offset.name = offset.name, adjust_covariates = adjust_covariates,
        adverse_outcome = adverse_outcome, frontier_rule = frontier_rule,
        effect_neighborhood = effect_neighborhood),
      error = function(e) grf_res)
  } else if (identical(grf_select_statistic, "effect") &&
             identical(grf_selection, "tree") && isTRUE(details)) {
    message("[forestsearch] grf_select_statistic = \"effect\" is ignored for ",
            "grf_selection = \"tree\" (the policy-tree leaf is the selection); ",
            "use grf_selection = \"frontier\" for effect-based re-ranking.")
  }

  # Structured, path-based subgroup definition (correct directions; a
  # disjunction of conjunctions when the GRF subgroup spans multiple leaves).
  # Built inside grf.subg.harm.* and returned as $sg_def.  Membership is
  # derived from this definition via .grf_evaluate_subgroup(), which reproduces
  # predict(tree, X) exactly and applies cleanly to df.predict / df.test.
  sg_def <- grf_res$sg_def
  found  <- !is.null(grf_res) && !inherits(grf_res, "try-error") &&
            !is.null(sg_def) && length(sg_def$conjunctions) > 0L

  if (isTRUE(details)) {
    lines <- c(
      "",
      "[forestsearch] GRF selection (subgroup_method = \"grf\")",
      paste0("  Outcome type:   ", outcome_type),
      paste0("  grf_depth:      ", grf_depth),
      paste0("  dmin.grf:       ", dmin.grf),
      paste0("  frac.tau:       ", frac.tau),
      paste0("  n.min:          ", n.min)
    )
    if (found) {
      lines <- c(lines,
                 paste0("  SELECTED (depth ",
                        if (!is.null(grf_res$selected_depth))
                          grf_res$selected_depth
                        else if (!is.null(grf_res$best_depth))
                          grf_res$best_depth else NA,
                        if (isTRUE(sg_def$is_disjunction))
                          ", union of leaves" else "",
                        "): ", sg_def$definition))
    } else {
      lines <- c(lines, "  SELECTED: none -- GRF identified no harm subgroup.")
    }
    message(paste(lines, collapse = "\n"))
  }

  if (!found) {
    return(list(found = FALSE, sg.harm = NULL, grp.consistency = NULL,
                grf_res = grf_res, df.est = df,
                df.predict = df.predict, df.test = df.test))
  }

  # sg.harm representation:
  #   * single conjunction -> brace-label vector ("{var op val}"), the same
  #     AND-composed representation a DINA/consistency subgroup uses;
  #   * union of leaves     -> a single disjunctive definition string
  #     "(a & b) | (c & d)", which get_dfpred() cannot evaluate, so membership
  #     for ALL cases is computed from the structured definition instead.
  sg.harm <- if (!is.null(sg_def$labels)) sg_def$labels else sg_def$definition

  # Membership from the structured definition (reproduces predict(tree, X);
  # valid on df / df.predict / df.test).
  tr_est  <- .grf_evaluate_subgroup(sg_def, df)
  df.est  <- df; df.est$treat.recommend <- tr_est
  df.predict_out <- if (!is.null(df.predict)) {
    o <- df.predict; o$treat.recommend <- .grf_evaluate_subgroup(sg_def, df.predict); o
  } else NULL
  df.test_out <- if (!is.null(df.test)) {
    o <- df.test; o$treat.recommend <- .grf_evaluate_subgroup(sg_def, df.test); o
  } else NULL

  id_vals <- if (id.name %in% names(df)) df[[id.name]]
             else if ("id" %in% names(df)) df[["id"]]
             else seq_len(nrow(df))
  df_flag <- data.frame(id = id_vals,
                        treat.recommend = df.est$treat.recommend)

  grp.consistency <- list(
    out_sg     = grf_res,
    sg.harm    = sg.harm,
    sg.harm.id = as.integer(df.est$treat.recommend == 0L),
    df_flag    = df_flag,
    algorithm  = "grf"
  )

  list(found = TRUE, sg.harm = sg.harm, grp.consistency = grp.consistency,
       grf_res = grf_res, df.est = df.est,
       df.predict = df.predict_out, df.test = df.test_out,
       candidates = grf_res$candidates,   # DR-candidate family for the Tier-2 gate
       select_statistic = grf_select_statistic)
}
