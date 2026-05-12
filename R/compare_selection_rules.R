# ============================================================================
# compare_selection_rules.R
# ============================================================================
# Wrapper to run forestsearch() across multiple (sg_focus, selection_rule)
# combinations with all other parameters held fixed, capture each run's
# in-flight diagnostic output, build the Pareto-frontier plots, and
# return a structured comparison object.  Pairs naturally with
# `pareto_frontier_table()`, `explain_pareto_selection()`, and the dual
# diagnostic printed by `forestsearch(show_candidate_summary = TRUE)`.
# ----------------------------------------------------------------------------

# Silence R CMD check NOTE for the inline data.table column reference.
utils::globalVariables(c("m"))

#' Compare forestsearch Runs Across Selection-Rule Combinations
#'
#' Runs \code{\link{forestsearch}} once per combination of
#' \code{sg_focus[i]} + \code{selection_rule[i]} (tuple semantics; the
#' two vectors are paired element-by-element and must have the same
#' length).  All other arguments are held fixed across runs and passed
#' through via \code{...}.
#'
#' For each combination the wrapper:
#' \enumerate{
#'   \item Captures the in-flight console output from forestsearch
#'     (including the pre-consistency \strong{Candidate Evaluation
#'     Preview} and post-consistency \strong{Candidate Evaluation
#'     Summary} when \code{show_candidate_summary = TRUE} is passed).
#'   \item Optionally computes the frontier CIs via
#'     \code{\link{compute_frontier_cis}}.
#'   \item Builds the Pareto-frontier plot via
#'     \code{\link{plot_pareto_frontier}}.
#' }
#'
#' The captured stdout is returned per combination so the caller can
#' \code{cat()} it later for side-by-side inspection in a Quarto chunk.
#' Plot objects are returned individually and (when
#' \pkg{patchwork} is installed) composed into a single side-by-side
#' figure.  The full \code{fs} object for each combination is returned
#' so downstream diagnostics (\code{\link{pareto_frontier_table}},
#' \code{\link{explain_pareto_selection}}) can be run without re-fitting.
#'
#' @section Tuple semantics:
#' With \code{sg_focus = c("hrMaxSG", "hrMaxSG")} and
#' \code{selection_rule = c("pareto", "both")}, two runs are
#' performed: one with \code{(hrMaxSG, pareto)} and one with
#' \code{(hrMaxSG, both)}.  The length-mismatch case raises an error.
#'
#' @section Console capture:
#' All stdout from each \code{forestsearch()} call is captured (not
#' just the diagnostic tables) via \code{capture.output()}.  The
#' PREVIEW / SUMMARY blocks are clearly delineated by banner separators
#' within the captured text.  Use \code{cat(out$console[[i]])} to
#' replay the output in a Quarto chunk; \code{results = "asis"} works
#' but is not required since the captured text already includes
#' newlines.
#'
#' @param df.analysis Analysis-ready data frame; passed to
#'   \code{forestsearch()}.
#' @param sg_focus Character vector of \code{sg_focus} values to try.
#' @param selection_rule Character vector of \code{selection_rule}
#'   values to try.  Must have the same length as \code{sg_focus}.
#' @param compute_cis Logical.  If \code{TRUE} (default), runs
#'   \code{compute_frontier_cis()} for each fit so the plot includes
#'   CI bars.  Set to \code{FALSE} to skip the (often expensive) CI
#'   step.
#' @param n_splits Integer.  Passed to
#'   \code{compute_frontier_cis(n_splits = )}.  Default \code{1000L}.
#' @param ci_seed Integer.  Passed to
#'   \code{compute_frontier_cis(seed = )}.  Default \code{1L}.
#' @param plot_xlim Numeric vector of length 2 (or \code{NULL}).
#'   Common x-axis range to apply to all plots so they're directly
#'   comparable.  Default \code{NULL} = let each plot auto-range.
#' @param show_band Logical.  Whether to draw the effect-band
#'   shading on each plot.  Default \code{TRUE}.
#' @param combo_labels Character vector of length \code{length(sg_focus)},
#'   used as the title for each plot and key for the returned lists.
#'   Default \code{NULL} = auto, e.g.\ \code{"hrMaxSG / pareto"}.
#' @param verbose Logical.  Print one line per combination during
#'   execution.  Default \code{TRUE}.
#' @param ... Additional arguments passed through to every
#'   \code{forestsearch()} call (data names, thresholds, parallel
#'   args, etc.).  Do NOT pass \code{sg_focus} or
#'   \code{selection_rule} here -- they're handled by the wrapper.
#'
#' @return A list with class \code{"forestsearch_comparison"} and
#'   components:
#'   \describe{
#'     \item{\code{combos}}{\code{data.frame} with columns
#'       \code{sg_focus}, \code{selection_rule}, \code{label}.}
#'     \item{\code{fs}}{Named list of \code{forestsearch} result
#'       objects, one per combination.  NULL entries mark fits that
#'       errored.}
#'     \item{\code{ci_tab}}{Named list of frontier-CI tables (or all
#'       NULL if \code{compute_cis = FALSE}).}
#'     \item{\code{plots}}{Named list of \code{ggplot} objects, one per
#'       combination.}
#'     \item{\code{plot_grid}}{A single \pkg{patchwork} object placing
#'       the plots side by side, or \code{NULL} if patchwork isn't
#'       installed.}
#'     \item{\code{console}}{Named list of character vectors -- the
#'       captured stdout from each forestsearch run.}
#'     \item{\code{errors}}{Named list of error messages (or NULL) for
#'       combos that failed.}
#'   }
#'
#' @examples
#' \dontrun{
#' out <- compare_selection_rules(
#'   df.analysis      = actg_df,
#'   sg_focus         = c("hrMaxSG", "hrMaxSG"),
#'   selection_rule   = c("pareto",  "both"),
#'   # All other forestsearch args:
#'   confounders.name = confounders.name,
#'   outcome.name     = adverse_outcome,
#'   treat.name       = treat.name,
#'   id.name          = id.name,
#'   outcome_type     = "continuous",
#'   effect_measure   = "MD",
#'   adverse_outcome  = TRUE,
#'   hr.threshold     = 10,
#'   hr.consistency   = 5,
#'   pconsistency.threshold = 0.90,
#'   show_candidate_summary = TRUE,
#'   details          = TRUE
#' )
#'
#' # Inspect captured console output for combo 1
#' cat(out$console[[1]])
#'
#' # Side-by-side plot
#' print(out$plot_grid)
#'
#' # Use fs object for downstream diagnostics
#' pareto_frontier_table(out$fs[[1]], ci_table = out$ci_tab[[1]])
#' explain_pareto_selection(out$fs[[2]], ci_table = out$ci_tab[[2]])
#' }
#'
#' @seealso \code{\link{forestsearch}},
#'   \code{\link{plot_pareto_frontier}},
#'   \code{\link{pareto_frontier_table}},
#'   \code{\link{explain_pareto_selection}},
#'   \code{\link{compute_frontier_cis}}.
#' @export
compare_selection_rules <- function(df.analysis,
                                    sg_focus,
                                    selection_rule,
                                    compute_cis  = TRUE,
                                    n_splits     = 1000L,
                                    ci_seed      = 1L,
                                    plot_xlim    = NULL,
                                    show_band    = TRUE,
                                    combo_labels = NULL,
                                    verbose      = TRUE,
                                    ...) {
  # --- 1. Validate inputs --------------------------------------------------
  if (!is.character(sg_focus) || !is.character(selection_rule)) {
    stop("`sg_focus` and `selection_rule` must be character vectors.",
         call. = FALSE)
  }
  if (length(sg_focus) != length(selection_rule)) {
    stop(sprintf(
      "Tuple-semantics mismatch: length(sg_focus) = %d but length(selection_rule) = %d.  The two vectors are paired element-by-element and must have the same length.  See ?compare_selection_rules for grid alternatives.",
      length(sg_focus), length(selection_rule)),
      call. = FALSE)
  }
  n_combos <- length(sg_focus)
  if (n_combos == 0L) {
    stop("`sg_focus` and `selection_rule` must have at least one element.",
         call. = FALSE)
  }

  # Auto-generate labels if not provided
  if (is.null(combo_labels)) {
    combo_labels <- paste(sg_focus, selection_rule, sep = " / ")
  } else {
    if (length(combo_labels) != n_combos) {
      stop("`combo_labels` length must match `sg_focus` and `selection_rule`.",
           call. = FALSE)
    }
  }
  # Ensure label uniqueness (used as list names)
  if (anyDuplicated(combo_labels) > 0L) {
    combo_labels <- make.unique(combo_labels, sep = "_")
  }

  # Guard against the user accidentally passing sg_focus/selection_rule in ...
  dots <- list(...)
  if (any(c("sg_focus", "selection_rule") %in% names(dots))) {
    stop("`sg_focus` and `selection_rule` are wrapper-controlled axes; do not pass them via `...`.",
         call. = FALSE)
  }

  # --- 2. Loop combinations -----------------------------------------------
  fs_list      <- vector("list", n_combos); names(fs_list)      <- combo_labels
  ci_list      <- vector("list", n_combos); names(ci_list)      <- combo_labels
  plot_list    <- vector("list", n_combos); names(plot_list)    <- combo_labels
  console_list <- vector("list", n_combos); names(console_list) <- combo_labels
  error_list   <- vector("list", n_combos); names(error_list)   <- combo_labels

  for (k in seq_len(n_combos)) {
    label_k <- combo_labels[k]
    focus_k <- sg_focus[k]
    rule_k  <- selection_rule[k]

    if (verbose) {
      cat(sprintf("[%d/%d] Running forestsearch: sg_focus = %s, selection_rule = %s\n",
                  k, n_combos, focus_k, rule_k))
    }

    # Capture stdout from the forestsearch call, including the
    # show_candidate_summary diagnostic banners.  Use an explicit
    # holder environment so the assignment survives the
    # capture.output / tryCatch nesting (a bare `<<-` can leak to the
    # global env depending on where the lookup finds a match).
    holder <- new.env(parent = emptyenv())
    holder$fs    <- NULL
    holder$error <- NULL
    captured <- tryCatch(
      utils::capture.output({
        holder$fs <- tryCatch(
          do.call(forestsearch, c(
            list(df.analysis    = df.analysis,
                 sg_focus       = focus_k,
                 selection_rule = rule_k),
            dots
          )),
          error = function(e) {
            holder$error <- conditionMessage(e)
            cat(sprintf("ERROR: %s\n", conditionMessage(e)))
            NULL
          }
        )
      }),
      error = function(e) {
        holder$error <- conditionMessage(e)
        character(0)
      }
    )
    fs_k  <- holder$fs
    err_k <- holder$error

    fs_list[[k]]      <- fs_k
    console_list[[k]] <- captured
    error_list[[k]]   <- err_k

    if (is.null(fs_k)) {
      if (verbose) cat(sprintf("    failed: %s\n",
                               if (is.null(err_k)) "no error reported" else err_k))
      next
    }

    # Compute frontier CIs (optional)
    if (isTRUE(compute_cis)) {
      ci_k <- tryCatch(
        compute_frontier_cis(fs_k, n_splits = n_splits, seed = ci_seed),
        error = function(e) {
          if (verbose) cat(sprintf("    compute_frontier_cis() failed: %s\n",
                                   conditionMessage(e)))
          NULL
        }
      )
      ci_list[[k]] <- ci_k
    }

    # Build the plot
    plot_k <- tryCatch(
      plot_pareto_frontier(
        fs        = fs_k,
        ci_table  = ci_list[[k]],
        show_band = show_band,
        xlim      = plot_xlim
      ),
      error = function(e) {
        if (verbose) cat(sprintf("    plot_pareto_frontier() failed: %s\n",
                                 conditionMessage(e)))
        NULL
      }
    )

    # Annotate each plot with its label (so a side-by-side composition reads cleanly)
    if (!is.null(plot_k) && requireNamespace("ggplot2", quietly = TRUE)) {
      plot_k <- plot_k + ggplot2::labs(subtitle = label_k)
    }
    plot_list[[k]] <- plot_k
  }

  # --- 3. Side-by-side composition (patchwork if available) ---------------
  plot_grid <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    valid_plots <- Filter(Negate(is.null), plot_list)
    if (length(valid_plots) > 0L) {
      plot_grid <- patchwork::wrap_plots(valid_plots, ncol = length(valid_plots))
    }
  } else if (verbose) {
    cat("Note: package 'patchwork' not installed; returning individual plots without side-by-side composition.\n")
  }

  # --- 4. Assemble return -------------------------------------------------
  combos_df <- data.frame(
    sg_focus       = sg_focus,
    selection_rule = selection_rule,
    label          = combo_labels,
    stringsAsFactors = FALSE
  )

  out <- list(
    combos     = combos_df,
    fs         = fs_list,
    ci_tab     = ci_list,
    plots      = plot_list,
    plot_grid  = plot_grid,
    console    = console_list,
    errors     = error_list
  )
  class(out) <- c("forestsearch_comparison", "list")
  out
}


#' Print Method for forestsearch_comparison
#'
#' Short summary of the comparison.  Use \code{print(x$plot_grid)} for
#' the side-by-side plot and \code{cat(x$console[[i]])} for each
#' combo's diagnostic output.
#'
#' @param x A \code{forestsearch_comparison} object.
#' @param ... Ignored.
#' @export
print.forestsearch_comparison <- function(x, ...) {
  cat("forestsearch comparison across selection-rule combinations\n")
  cat("---\n")
  print(x$combos)
  cat("---\n")
  n_ok    <- sum(!vapply(x$fs,     is.null, logical(1)))
  n_plots <- sum(!vapply(x$plots,  is.null, logical(1)))
  n_ci    <- sum(!vapply(x$ci_tab, is.null, logical(1)))
  cat(sprintf("Successful fits: %d / %d\n", n_ok, nrow(x$combos)))
  cat(sprintf("Plots built:     %d / %d\n", n_plots, nrow(x$combos)))
  cat(sprintf("CIs computed:    %d / %d\n", n_ci, nrow(x$combos)))
  cat("\nAccess the contents via x$fs, x$ci_tab, x$plots, x$plot_grid, x$console.\n")
  cat("Replay diagnostic output for combo i:  cat(x$console[[i]])\n")
  cat("View side-by-side plot:                print(x$plot_grid)\n")
  invisible(x)
}
