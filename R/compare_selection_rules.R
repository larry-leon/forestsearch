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


# ============================================================================
# extract_candidate_diagnostics() -- internal helper
# ============================================================================
# Slice forestsearch()'s captured stdout into the PREVIEW and SUMMARY
# blocks so the comparison qmd (or any downstream caller) can present
# just the diagnostic tables as the main display, with the full
# captured output kept accessible as a collapsible / on-demand view.
#
# Banner format (defined by print_candidate_preview / print_candidate_summary):
#
#   ====================== ... (line of '=')
#   CANDIDATE EVALUATION PREVIEW (pre-consistency) (...)
#   ====================== ... (line of '=')
#   ... body ...
#   ====================== ... (line of '=')
#
#   ====================== ... (line of '=')
#   CANDIDATE EVALUATION SUMMARY  (...)
#   ====================== ... (line of '=')
#   ... body ...
#   ====================== ... (line of '=')
#
# A block is delimited by THREE '=' lines: the top, the post-banner
# divider, and the bottom.  We slice from the FIRST '=' line of the
# header to the THIRD '=' line, inclusive, so the user sees the full
# framed block.
# ----------------------------------------------------------------------------

#' Slice the PREVIEW and SUMMARY blocks out of captured forestsearch stdout
#'
#' Internal helper used by \code{\link{compare_selection_rules}} to
#' populate the per-combo \code{diagnostics} slot.  Identifies the
#' \code{CANDIDATE EVALUATION PREVIEW (pre-consistency)} and
#' \code{CANDIDATE EVALUATION SUMMARY} blocks in a captured-output
#' character vector and returns them as separate strings, alongside
#' the full original.
#'
#' If a block is not found (e.g.\ \code{show_candidate_summary = FALSE}
#' was passed, or the run errored before reaching that block) the
#' corresponding element is \code{NA_character_}.
#'
#' @param captured Character vector as returned by
#'   \code{\link{capture.output}}.  One element per line.
#'
#' @return A list with character-scalar fields \code{preview},
#'   \code{summary}, and \code{full}.  The first two are NA when the
#'   block is absent; \code{full} always has the original text.
#'
#' @keywords internal
#' @noRd
extract_candidate_diagnostics <- function(captured) {
  if (is.null(captured) || length(captured) == 0L) {
    return(list(preview = NA_character_,
                summary = NA_character_,
                full    = ""))
  }

  # Find the header lines.  Use fixed = TRUE so the banner phrase is
  # treated literally (not regex), and trimws() to be tolerant of any
  # leading/trailing whitespace that might creep in via different cat()
  # call sites.
  ix_preview_hdr <- which(grepl(
    "CANDIDATE EVALUATION PREVIEW", captured, fixed = TRUE))
  ix_summary_hdr <- which(grepl(
    "CANDIDATE EVALUATION SUMMARY", captured, fixed = TRUE))

  # Bar separator: a line of one or more '=' characters (no other
  # printable content).
  is_bar <- grepl("^=+\\s*$", captured)
  bar_ix <- which(is_bar)

  slice_block <- function(hdr_line) {
    # The first bar above the header is the block-open bar; we want
    # the three-bar window: open / under-header / close.
    if (length(hdr_line) == 0L) return(NA_character_)
    hdr_line <- hdr_line[1L]
    bars_before <- bar_ix[bar_ix < hdr_line]
    bars_after  <- bar_ix[bar_ix > hdr_line]
    if (length(bars_before) < 1L || length(bars_after) < 2L) {
      # Not enough bars to bound the block; fall back to including
      # everything from the header to the next bar (or end).
      start_ix <- hdr_line
      end_ix   <- if (length(bars_after) >= 1L) bars_after[1L] else length(captured)
    } else {
      start_ix <- max(bars_before)        # bar immediately above header
      end_ix   <- bars_after[2L]          # second bar after header = closing bar
    }
    paste(captured[start_ix:end_ix], collapse = "\n")
  }

  list(
    preview = slice_block(ix_preview_hdr),
    summary = slice_block(ix_summary_hdr),
    full    = paste(captured, collapse = "\n")
  )
}


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
#'     \item{\code{plot_combined}}{A single \code{ggplot} composing all
#'       configurations onto one Pareto plot, with each winner labeled
#'       \code{S1: <combo_label>}, \code{S2: <combo_label>}, etc.
#'       Returned only when ALL successful configurations share an
#'       identical passing set (see \code{\link{plot_pareto_combined}}
#'       for the equality criterion); \code{NULL} otherwise.  See
#'       \code{plot_combined_subsets} below for the per-subset
#'       fallback when some -- but not all -- configurations match.}
#'     \item{\code{plot_combined_subsets}}{Named list of \code{ggplot}
#'       objects, one per equivalence class of configurations sharing
#'       an identical passing set.  Only groups of size \eqn{\ge 2}
#'       are included (a singleton group has no peer to combine
#'       with).  Names are derived from the shared \code{sg_focus}
#'       when constant within a group (e.g., \code{"effMaxSG"}); from
#'       the concatenated combo labels otherwise.  When all valid
#'       combos share one passing set this list contains exactly one
#'       element, equal to \code{plot_combined}.  The list is empty
#'       when no two combos match.}
#'     \item{\code{combined_skip_reason}}{Character scalar or
#'       \code{NULL}.  When \code{plot_combined} is \code{NULL}
#'       because the equality precondition failed, this captures the
#'       specific reason -- size mismatch, definition-set mismatch,
#'       or value drift on a named column -- as a single string.
#'       \code{NULL} when the combined plot was built or when fewer
#'       than two fits succeeded.}
#'     \item{\code{console}}{Named list of character vectors -- the
#'       captured stdout from each forestsearch run (full output).}
#'     \item{\code{diagnostics}}{Named list of per-combo slices of the
#'       captured output: each element is a list with character-scalar
#'       fields \code{preview} (the \code{CANDIDATE EVALUATION
#'       PREVIEW} block), \code{summary} (the \code{CANDIDATE
#'       EVALUATION SUMMARY} block), and \code{full} (the entire
#'       captured stdout, same content as \code{console[[i]]}
#'       collapsed to a single string).  \code{preview} and
#'       \code{summary} are \code{NA_character_} when the corresponding
#'       banner is absent (e.g.\ \code{show_candidate_summary} was
#'       \code{FALSE}, or the run errored before reaching it).  This
#'       slice is what a comparison document should display as the
#'       primary diagnostic; \code{full} (or \code{console[[i]]})
#'       remains available for an expandable / on-demand view of the
#'       complete run output.}
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
#'   \code{\link{plot_pareto_combined}},
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
  fs_list          <- vector("list", n_combos); names(fs_list)          <- combo_labels
  ci_list          <- vector("list", n_combos); names(ci_list)          <- combo_labels
  plot_list        <- vector("list", n_combos); names(plot_list)        <- combo_labels
  console_list     <- vector("list", n_combos); names(console_list)     <- combo_labels
  diagnostics_list <- vector("list", n_combos); names(diagnostics_list) <- combo_labels
  error_list       <- vector("list", n_combos); names(error_list)       <- combo_labels

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

    fs_list[[k]]          <- fs_k
    console_list[[k]]     <- captured
    diagnostics_list[[k]] <- extract_candidate_diagnostics(captured)
    error_list[[k]]       <- err_k

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

  # --- 3b. Combined plot(s) by matching passing-set group ---------------
  #
  # plot_pareto_combined() requires combos to share an identical passing
  # set (same subgroup definitions, same hr/N/E/K values).  The first
  # version of this wrapper passed ALL valid combos at once: if even one
  # disagreed, no combined plot was built.  That all-or-nothing behavior
  # missed cases where natural subsets DO agree.  For example, the
  # comparison vignette with
  #
  #     sg_focus       = c("effMaxSG", "effMaxSG", "effMinSG", "effMinSG")
  #     selection_rule = c("pareto",   "both",     "pareto",   "both")
  #
  # produces matching passing sets within each sg_focus value (combos
  # 1+2 agree; combos 3+4 agree) but not across them (the preview-sort
  # order under effMaxSG vs effMinSG drives different random-split
  # seeds, which can push borderline candidates over or under the
  # pconsistency threshold).
  #
  # Strategy: compute a passing-set signature for each valid fs object,
  # group combos by signature, then build one combined plot per group
  # of size >= 2.  Populate the backward-compatible $plot_combined slot
  # iff there's a single group covering all valid combos.
  plot_combined         <- NULL
  plot_combined_subsets <- list()
  combined_warn         <- NULL
  valid_idx <- which(!vapply(fs_list, is.null, logical(1)))

  if (length(valid_idx) >= 2L && exists("plot_pareto_combined", mode = "function")) {

    # Compute a signature per valid fs: the sorted-collapsed
    # subgroup-definition set.  Two fs objects with identical
    # signatures share an identical passing set (modulo value drift,
    # which plot_pareto_combined() also checks but is typically zero
    # across rules on the same data).
    sig_for_fs <- function(fs) {
      out_sg <- tryCatch(fs$grp.consistency$out_sg, error = function(e) NULL)
      if (is.null(out_sg) || is.null(out_sg$result) ||
          nrow(out_sg$result) == 0L) return(NA_character_)
      res <- out_sg$result
      m_cols <- grep("^M\\.", names(res), value = TRUE)
      keys <- vapply(seq_len(nrow(res)), function(i) {
        v <- unlist(res[i, m_cols, with = FALSE], use.names = FALSE)
        v <- v[!is.na(v) & nzchar(v)]
        if (length(v) == 0L) "(empty)" else paste(sort(v), collapse = " & ")
      }, character(1))
      paste(sort(keys), collapse = " | ")
    }
    sigs <- vapply(fs_list[valid_idx], sig_for_fs, character(1))

    # Group valid combo indices by signature.
    groups <- split(valid_idx, sigs)
    # Drop singleton groups (combined plot needs >= 2 combos).
    groups <- groups[lengths(groups) >= 2L]

    # Build a combined plot per group.  Use the first combo's
    # sg_focus + selection_rule as the group's natural label,
    # falling back to combo labels otherwise.
    for (g in seq_along(groups)) {
      idx <- groups[[g]]
      group_warn <- NULL
      group_plot <- tryCatch(
        withCallingHandlers(
          plot_pareto_combined(
            fs_list       = fs_list[idx],
            combo_labels  = combo_labels[idx],
            ci_table_list = ci_list[idx],
            show_band     = show_band,
            xlim          = plot_xlim,
            verbose       = TRUE
          ),
          warning = function(w) {
            group_warn <<- conditionMessage(w)
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) NULL
      )
      if (!is.null(group_plot)) {
        # Label this group by the shared sg_focus, when constant; else
        # by the concatenated combo labels.
        sg_in_group <- unique(sg_focus[idx])
        group_label <- if (length(sg_in_group) == 1L) {
          sg_in_group
        } else {
          paste(combo_labels[idx], collapse = " + ")
        }
        plot_combined_subsets[[group_label]] <- group_plot
      } else if (!is.null(group_warn)) {
        # A group that signature-matched but plot_pareto_combined()
        # still rejected (e.g. value drift inside matching definitions).
        # Surface the most recent such reason via combined_warn.
        combined_warn <- group_warn
      }
    }

    # Backward-compatible: $plot_combined is set when ALL valid combos
    # share one signature (i.e., exactly one group, covering everything).
    if (length(groups) == 1L && length(groups[[1L]]) == length(valid_idx)) {
      plot_combined <- plot_combined_subsets[[1L]]
    } else if (length(groups) == 0L) {
      # No group of size >= 2; report a useful reason.
      if (is.null(combined_warn)) {
        combined_warn <- sprintf(
          "All %d valid configurations have distinct passing sets; no two combos can be combined.",
          length(valid_idx))
      }
    } else {
      # Multiple matching subsets, but no single group covers all combos.
      sizes <- vapply(groups, length, integer(1))
      combined_warn <- sprintf(
        paste0("Passing sets fall into %d matching group(s) of sizes %s; ",
               "no single combined plot covers all %d configurations.  ",
               "See $plot_combined_subsets for per-group plots."),
        length(groups), paste(sizes, collapse = ", "), length(valid_idx))
    }

    if (verbose) {
      n_subsets <- length(plot_combined_subsets)
      if (!is.null(plot_combined)) {
        cat("Combined plot built (all", length(valid_idx),
            "configurations share one passing set).\n")
      } else if (n_subsets > 0L) {
        cat(sprintf(
          "Combined plot built per matching subset (%d subset%s: %s).\n",
          n_subsets, if (n_subsets == 1L) "" else "s",
          paste(names(plot_combined_subsets), collapse = ", ")))
        if (!is.null(combined_warn)) {
          cat("  ", combined_warn, "\n", sep = "")
        }
      } else if (!is.null(combined_warn)) {
        cat("Combined plot skipped:\n  ", combined_warn, "\n", sep = "")
      } else {
        cat("Combined plot skipped (no diagnostic captured).\n")
      }
    }
  }

  # --- 4. Assemble return -------------------------------------------------
  combos_df <- data.frame(
    sg_focus       = sg_focus,
    selection_rule = selection_rule,
    label          = combo_labels,
    stringsAsFactors = FALSE
  )

  out <- list(
    combos        = combos_df,
    fs            = fs_list,
    ci_tab        = ci_list,
    plots         = plot_list,
    plot_grid     = plot_grid,
    plot_combined = plot_combined,
    plot_combined_subsets = plot_combined_subsets,
    combined_skip_reason = combined_warn,  # NULL if combined built or never attempted
    console       = console_list,
    diagnostics   = diagnostics_list,
    errors        = error_list
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
  n_prev  <- sum(vapply(x$diagnostics,
                        function(d) !is.null(d) && !is.na(d$preview),
                        logical(1)))
  n_sumb  <- sum(vapply(x$diagnostics,
                        function(d) !is.null(d) && !is.na(d$summary),
                        logical(1)))
  cat(sprintf("Successful fits: %d / %d\n", n_ok, nrow(x$combos)))
  cat(sprintf("Plots built:     %d / %d\n", n_plots, nrow(x$combos)))
  cat(sprintf("CIs computed:    %d / %d\n", n_ci, nrow(x$combos)))
  cat(sprintf("Diagnostics:     PREVIEW found %d / %d; SUMMARY found %d / %d\n",
              n_prev, nrow(x$combos), n_sumb, nrow(x$combos)))
  n_subsets <- if (is.null(x$plot_combined_subsets)) 0L
               else length(x$plot_combined_subsets)
  cat(sprintf("Combined plot:   %s\n",
              if (!is.null(x$plot_combined)) "yes (all configurations share one passing set)"
              else if (n_subsets > 0L)
                sprintf("per-subset only (%d matching subset%s: %s)",
                        n_subsets,
                        if (n_subsets == 1L) "" else "s",
                        paste(names(x$plot_combined_subsets), collapse = ", "))
              else "no (see combined_skip_reason)"))
  if (!is.null(x$combined_skip_reason)) {
    msg <- x$combined_skip_reason
    cat("  note:   ", msg, "\n", sep = "")
  }
  cat("\nAccess the contents via x$fs, x$ci_tab, x$plots, x$plot_grid,\n",
      "  x$plot_combined, x$plot_combined_subsets, x$diagnostics, x$console.\n", sep = "")
  cat("Show PREVIEW + SUMMARY for combo i:    cat(x$diagnostics[[i]]$preview,\n",
      "                                          x$diagnostics[[i]]$summary,\n",
      "                                          sep = \"\\n\\n\")\n", sep = "")
  cat("Full captured output for combo i:      cat(x$console[[i]], sep = \"\\n\")\n")
  cat("View side-by-side plot:                print(x$plot_grid)\n")
  cat("View combined plot (if available):     print(x$plot_combined)\n")
  invisible(x)
}
