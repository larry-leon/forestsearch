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
