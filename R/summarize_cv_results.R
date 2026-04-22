# =============================================================================
# summarize_cv_results.R
#
# Post-hoc diagnostic summary for forestsearch_tenfold() output.  Analogous to
# summarize_bootstrap_results() but operating on the sim x fold grid captured
# in `fold_summary`.
#
# Exports:
#   - summarize_cv_results()      : main user-facing function
#   - print.fs_cv_summary()       : print method for returned object
#
# Internal helpers (dot-prefixed):
#   - .cv_build_sg_full()         : combine sg1, sg2 into canonical label
#   - .cv_clean_grf()             : normalise NA / empty grf_cuts strings
#   - .cv_canonicalize_sg()       : canonical form for exact-match comparison
#   - .cv_tabulate_subgroups()    : frequency table of identified subgroups
#   - .cv_tabulate_grf_cuts()     : frequency table of raw GRF outputs
#   - .cv_xtab_cut_vs_sg()        : GRF cut x identified subgroup cross-tab
#   - .cv_decompose_no_sg()       : breakdown of no-subgroup folds
#   - .cv_original_agreement()    : agreement with full-data analysis
#   - .cv_metrics_table_safe()    : delegate to cv_metrics_tables()
#   - .cv_diagnostic_plots()      : optional ggplot2 diagnostics
#   - .format_cv_*()              : gt-table formatters for each component
# =============================================================================


#' Summarise ForestSearch Cross-Validation Diagnostics
#'
#' Post-hoc diagnostic summary of \code{\link{forestsearch_tenfold}} output,
#' aggregating the sim x fold grid captured in \code{fold_summary} into
#' publication-ready \code{gt} tables.  Parallel to
#' \code{\link{summarize_bootstrap_results}} but focused on cross-validation
#' stability and GRF candidate-discovery behaviour.
#'
#' @section Diagnostic components:
#' The returned object includes (as \code{gt} tables by default):
#' \itemize{
#'   \item \code{identification_summary}: top identified subgroups across
#'     all sim x fold pairs, with counts and percentages.
#'   \item \code{grf_cut_summary}: top raw \code{grf_cuts} strings across
#'     training sets -- how often GRF returned each output (and how often
#'     it returned none).  \code{NULL} if \code{fold_summary$grf_cuts} is
#'     not available (e.g., objects from \pkg{forestsearch} prior to the
#'     GLM-extension branch).
#'   \item \code{cut_vs_subgroup_xtab}: cross-tabulation of unique
#'     \code{grf_cuts} strings against identified subgroups.  Reveals
#'     whether GRF's discovery aligns with the pipeline's final selection.
#'   \item \code{no_subgroup_decomposition}: tabulation of folds that
#'     identified no subgroup, broken into "GRF returned no cut" vs
#'     "GRF returned a cut but consistency rejected it".
#'   \item \code{pconsistency_distribution}: binned distribution of the
#'     consistency probability (\code{Pcons}) achieved by the identified
#'     subgroup on each fold, restricted to folds that identified a
#'     subgroup; includes median and IQR.  Reveals whether identifying
#'     folds are scraping \code{pconsistency.threshold} (e.g., Pcons just
#'     above 0.90) or are well above it.  \code{NULL} if
#'     \code{fold_summary$pconsistency} is not available.  Rejected-
#'     candidate Pcons values are not surfaced -- only the selected
#'     candidate's Pcons is captured in \code{fold_summary}.
#'   \item \code{fold_numeric_summary}: tidy median / IQR / range table
#'     with one row per numeric diagnostic column present in
#'     \code{fold_summary} (\code{n_test}, \code{pconsistency},
#'     \code{training_fs_hr}, \code{n_candidates_evaluated}).  Absent
#'     columns simply don't produce rows, so the table adapts to
#'     pre-feature-branch objects without error.  Training-fold
#'     subgroup effects are on the effect measure's natural scale
#'     (HR for survival; OR, RR, IRR, RD, IRD, or MD for GLM,
#'     exponentiated at capture for ratio measures so the column is
#'     always on natural scale regardless of outcome type); they are
#'     in-sample estimates (optimistically biased) and are surfaced
#'     for diagnostic comparison across folds only.
#'   \item \code{original_agreement}: if \code{original_sg} and/or
#'     \code{original_grf_cuts} are supplied, fraction of folds matching
#'     the original full-data analysis (exact subgroup, partial shared
#'     factor, GRF-cut match, both).
#'   \item \code{metrics_table}: formatted version of
#'     \code{cv_output$sens_summary} and \code{cv_output$find_summary}
#'     (delegates to \code{\link{cv_metrics_tables}}).
#'   \item \code{plots}: if \code{create_plots = TRUE}, bar charts of
#'     top identified subgroups and top GRF cuts.
#' }
#' Raw \code{data.frame} versions of each tabulation are available at
#' \code{$data$*} for custom post-processing.
#'
#' @section Backward compatibility:
#' The GRF-dependent slots (\code{grf_cut_summary},
#' \code{cut_vs_subgroup_xtab}, \code{no_subgroup_decomposition}) require
#' a \code{grf_cuts} column in \code{fold_summary}; the
#' \code{pconsistency_distribution} slot requires a \code{pconsistency}
#' column.  Both were added on the \code{feature/glm-extension} branch
#' (\code{grf_cuts} first, \code{pconsistency} later).  Older objects
#' are detected and the affected slots return \code{NULL} without error.
#' Metadata fields \code{has_grf_cuts} and \code{has_pconsistency}
#' record which columns were present.
#'
#' @param cv_output List of class \code{"fs_tenfold"}, as returned by
#'   \code{\link{forestsearch_tenfold}}.  Must contain a non-empty
#'   \code{fold_summary} data frame.
#' @param original_sg Character vector or \code{NULL}.  Original full-data
#'   subgroup definition (e.g., \code{fs$sg.harm}).  When supplied,
#'   \code{original_agreement} reports exact and partial-factor match
#'   rates against this definition.
#' @param original_grf_cuts Character vector or \code{NULL}.  Original
#'   full-data \code{grf_cuts} (e.g., \code{fs$grf_cuts}).  When supplied
#'   alongside a \code{grf_cuts} column in \code{fold_summary},
#'   \code{original_agreement} also reports GRF-cut match rate.
#' @param create_plots Logical.  If \code{TRUE} and \pkg{ggplot2} is
#'   installed, produce bar-chart diagnostics (default: \code{FALSE}).
#' @param top_n Integer.  Maximum number of detail rows to retain in
#'   frequency tables (default: \code{15L}).  If the true number of
#'   distinct values exceeds \code{top_n}, a single
#'   \code{"(N others)"} row aggregates the tail so that \code{n_folds}
#'   totals remain honest.  Set to a large value (e.g., \code{Inf}) to
#'   disable truncation entirely.  Plots use \code{min(top_n, 10L)}.
#'
#' @return An object of class \code{"fs_cv_summary"}, a list with
#'   components:
#'   \describe{
#'     \item{\code{identification_summary}, \code{grf_cut_summary},
#'       \code{cut_vs_subgroup_xtab}, \code{no_subgroup_decomposition},
#'       \code{pconsistency_distribution}, \code{fold_numeric_summary},
#'       \code{original_agreement}, \code{metrics_table}}{\code{gt}
#'       tables (or \code{data.frame} fallback if \pkg{gt} is
#'       unavailable).  Some may be \code{NULL} depending on input
#'       columns and optional arguments.}
#'     \item{\code{plots}}{\code{NULL} or a list of \code{ggplot} objects.}
#'     \item{\code{data}}{List of raw \code{data.frame} tabulations
#'       underlying the tables (\code{identification}, \code{grf_cuts},
#'       \code{cut_vs_subgroup}, \code{no_subgroup}, \code{pconsistency},
#'       \code{fold_numeric_summary}, \code{original_agreement}).}
#'     \item{\code{n_sims}, \code{n_folds}, \code{total_pairs}}{Grid
#'       dimensions.}
#'     \item{\code{has_grf_cuts}, \code{has_pconsistency}}{Logical.
#'       Whether the corresponding \code{fold_summary} columns were
#'       available.  The \code{fold_numeric_summary} slot adapts
#'       automatically to whichever numeric diagnostic columns are
#'       present and does not have a dedicated flag.}
#'   }
#'
#' @seealso \code{\link{forestsearch_tenfold}} for running repeated
#'   K-fold CV; \code{\link{cv_metrics_tables}} for the underlying
#'   metrics formatter; \code{\link{summarize_bootstrap_results}} for
#'   the bootstrap analogue.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' df <- survival::gbsg
#' df$grade3      <- as.integer(df$grade == "3")
#' df$time_months <- df$rfstime / 30.4375
#'
#' fs <- forestsearch(
#'   df.analysis      = df,
#'   confounders.name = c("age", "meno", "size", "grade3",
#'                        "nodes", "pgr", "er"),
#'   outcome.name     = "time_months",
#'   event.name       = "status",
#'   treat.name       = "hormon"
#' )
#'
#' tf <- forestsearch_tenfold(fs.est = fs, sims = 100, Kfolds = 10)
#'
#' # Basic diagnostic summary
#' cv_diag <- summarize_cv_results(tf)
#' cv_diag$identification_summary
#' cv_diag$grf_cut_summary
#' cv_diag$no_subgroup_decomposition
#'
#' # With original-analysis comparison and plots
#' cv_diag_full <- summarize_cv_results(
#'   tf,
#'   original_sg        = fs$sg.harm,
#'   original_grf_cuts  = fs$grf_cuts,
#'   create_plots       = TRUE
#' )
#' cv_diag_full$original_agreement
#' cv_diag_full$plots$identification
#' }
#'
#' @importFrom utils head
#' @export
summarize_cv_results <- function(cv_output,
                                 original_sg       = NULL,
                                 original_grf_cuts = NULL,
                                 create_plots      = FALSE,
                                 top_n             = 15L) {

  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(cv_output)) {
    stop("`cv_output` must be a list returned by forestsearch_tenfold().",
         call. = FALSE)
  }

  if (inherits(cv_output, "fs_kfold") && !inherits(cv_output, "fs_tenfold")) {
    stop(
      "`summarize_cv_results()` operates on forestsearch_tenfold() output ",
      "(class 'fs_tenfold'), which exposes a `fold_summary` data frame.\n",
      "For single K-fold results (class 'fs_kfold'), use ",
      "`forestsearch_KfoldOut(res, outall = TRUE)` together with ",
      "`cv_summary_tables()` / `cv_metrics_tables()`.",
      call. = FALSE
    )
  }

  if (!inherits(cv_output, "fs_tenfold")) {
    warning(
      "`cv_output` does not inherit from class 'fs_tenfold'; ",
      "attempting to proceed using duck-typed `fold_summary` slot.",
      call. = FALSE
    )
  }

  fsum <- cv_output$fold_summary
  if (is.null(fsum) || !is.data.frame(fsum) || nrow(fsum) == 0L) {
    stop(
      "`cv_output$fold_summary` is missing or empty. ",
      "summarize_cv_results() requires the per-fold summary produced by ",
      "forestsearch_tenfold() on the feature/glm-extension branch.",
      call. = FALSE
    )
  }

  required_cols <- c("sim", "fold", "sg1", "sg2", "any_found")
  missing_cols  <- setdiff(required_cols, names(fsum))
  if (length(missing_cols) > 0L) {
    stop(
      "`fold_summary` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
      top_n < 1L) {
    stop("`top_n` must be a positive integer (or `Inf` to disable ",
         "truncation).", call. = FALSE)
  }
  top_n <- if (is.finite(top_n)) as.integer(top_n) else Inf

  # ---------------------------------------------------------------------------
  # FEATURE DETECTION AND NORMALISATION
  # ---------------------------------------------------------------------------
  has_grf_cuts <- "grf_cuts" %in% names(fsum)

  n_sims      <- length(unique(fsum$sim))
  n_folds     <- length(unique(fsum$fold))
  total_pairs <- nrow(fsum)

  fsum$sg_full <- .cv_build_sg_full(fsum$sg1, fsum$sg2)

  if (has_grf_cuts) {
    fsum$grf_clean <- vapply(fsum$grf_cuts, .cv_clean_grf, character(1))
  } else {
    fsum$grf_clean <- NA_character_
  }

  # ---------------------------------------------------------------------------
  # COMPONENT 1: IDENTIFICATION SUMMARY
  # ---------------------------------------------------------------------------
  identification_df <- .cv_tabulate_subgroups(
    fsum$sg_full, total_pairs, top_n
  )

  # ---------------------------------------------------------------------------
  # COMPONENT 2: GRF CUT SUMMARY
  # ---------------------------------------------------------------------------
  grf_cuts_df <- if (has_grf_cuts) {
    .cv_tabulate_grf_cuts(fsum$grf_clean, total_pairs, top_n)
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # COMPONENT 3: CUT vs SUBGROUP CROSS-TAB
  # ---------------------------------------------------------------------------
  cut_vs_sg_df <- if (has_grf_cuts) {
    .cv_xtab_cut_vs_sg(fsum$grf_clean, fsum$sg_full, top_n)
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # COMPONENT 4: NO-SUBGROUP DECOMPOSITION
  # ---------------------------------------------------------------------------
  no_sg_df <- if (has_grf_cuts) {
    .cv_decompose_no_sg(fsum$sg_full, fsum$grf_clean, total_pairs)
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # COMPONENT 4b: PCONSISTENCY DISTRIBUTION (Phase B)
  # ---------------------------------------------------------------------------
  has_pconsistency <- "pconsistency" %in% names(fsum)
  pconsistency_df  <- if (has_pconsistency) {
    .cv_pconsistency_distribution(fsum$pconsistency, fsum$sg_full,
                                   total_pairs)
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # COMPONENT 4c: FOLD NUMERIC SUMMARY (Phase B)
  #
  # Tidy table of median / IQR / range / n for every numeric diagnostic
  # column present in fold_summary.  Backward-compatible by row: absent
  # columns simply don't produce rows.  Preferred over multiplying
  # single-column distribution slots whenever a new diagnostic is added.
  # ---------------------------------------------------------------------------
  fold_numeric_df <- .cv_fold_numeric_summary(fsum)

  # ---------------------------------------------------------------------------
  # COMPONENT 5: ORIGINAL-ANALYSIS AGREEMENT
  # ---------------------------------------------------------------------------
  original_agreement_df <- .cv_original_agreement(
    fsum              = fsum,
    original_sg       = original_sg,
    original_grf_cuts = original_grf_cuts,
    has_grf_cuts      = has_grf_cuts
  )

  # ---------------------------------------------------------------------------
  # COMPONENT 6: METRICS TABLE (delegation)
  # ---------------------------------------------------------------------------
  metrics_table <- .cv_metrics_table_safe(cv_output, original_sg)

  # ---------------------------------------------------------------------------
  # COMPONENT 7: PLOTS (optional)
  # ---------------------------------------------------------------------------
  plots <- if (isTRUE(create_plots)) {
    .cv_diagnostic_plots(fsum, has_grf_cuts, top_n = min(top_n, 10L))
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # FORMAT gt TABLES (graceful fallback to data.frame if gt missing)
  # ---------------------------------------------------------------------------
  gt_available <- requireNamespace("gt", quietly = TRUE)
  if (!gt_available) {
    warning(
      "Package 'gt' not installed; returning data.frame tables. ",
      "Install gt for publication-ready formatted output.",
      call. = FALSE
    )
  }

  identification_summary <- if (gt_available) {
    .format_cv_identification_gt(
      identification_df, total_pairs, n_sims, n_folds
    )
  } else {
    identification_df
  }

  grf_cut_summary <- if (!is.null(grf_cuts_df)) {
    if (gt_available) {
      .format_cv_grf_cut_gt(grf_cuts_df, total_pairs, n_sims, n_folds)
    } else {
      grf_cuts_df
    }
  } else {
    NULL
  }

  cut_vs_subgroup_xtab <- if (!is.null(cut_vs_sg_df)) {
    if (gt_available) {
      .format_cv_xtab_gt(cut_vs_sg_df, total_pairs)
    } else {
      cut_vs_sg_df
    }
  } else {
    NULL
  }

  no_subgroup_decomposition <- if (!is.null(no_sg_df)) {
    if (gt_available) {
      .format_cv_no_sg_gt(no_sg_df, total_pairs)
    } else {
      no_sg_df
    }
  } else {
    NULL
  }

  pconsistency_distribution <- if (!is.null(pconsistency_df)) {
    if (gt_available) {
      .format_cv_pconsistency_gt(pconsistency_df, total_pairs)
    } else {
      pconsistency_df
    }
  } else {
    NULL
  }

  fold_numeric_summary <- if (!is.null(fold_numeric_df) &&
                              nrow(fold_numeric_df) > 0L) {
    if (gt_available) {
      .format_cv_fold_numeric_summary_gt(fold_numeric_df, total_pairs)
    } else {
      fold_numeric_df
    }
  } else {
    NULL
  }

  original_agreement <- if (!is.null(original_agreement_df)) {
    if (gt_available) {
      .format_cv_original_agreement_gt(original_agreement_df)
    } else {
      original_agreement_df
    }
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # ASSEMBLE AND RETURN
  # ---------------------------------------------------------------------------
  structure(
    list(
      # gt tables (or data.frame fallback)
      identification_summary    = identification_summary,
      grf_cut_summary           = grf_cut_summary,
      cut_vs_subgroup_xtab      = cut_vs_subgroup_xtab,
      no_subgroup_decomposition = no_subgroup_decomposition,
      pconsistency_distribution = pconsistency_distribution,
      fold_numeric_summary      = fold_numeric_summary,
      original_agreement        = original_agreement,
      metrics_table             = metrics_table,
      plots                     = plots,

      # Raw data frames for custom post-processing
      data = list(
        identification       = identification_df,
        grf_cuts             = grf_cuts_df,
        cut_vs_subgroup      = cut_vs_sg_df,
        no_subgroup          = no_sg_df,
        pconsistency         = pconsistency_df,
        fold_numeric_summary = fold_numeric_df,
        original_agreement   = original_agreement_df
      ),

      # Metadata
      n_sims           = n_sims,
      n_folds          = n_folds,
      total_pairs      = total_pairs,
      has_grf_cuts     = has_grf_cuts,
      has_pconsistency = has_pconsistency
    ),
    class = c("fs_cv_summary", "list")
  )
}


#' Print Method for fs_cv_summary Objects
#'
#' @param x An \code{fs_cv_summary} object from
#'   \code{\link{summarize_cv_results}}.
#' @param n Integer.  Maximum rows to print from each raw data frame
#'   (default: \code{5L}).
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @export
print.fs_cv_summary <- function(x, n = 5L, ...) {
  cat("ForestSearch CV Diagnostic Summary\n")
  cat("==================================\n")
  cat(sprintf("Simulations:            %d\n", x$n_sims))
  cat(sprintf("Folds per simulation:   %d\n", x$n_folds))
  cat(sprintf("Total sim x fold pairs: %d\n", x$total_pairs))
  cat(sprintf("grf_cuts available:     %s\n\n",
              if (x$has_grf_cuts) "yes" else "no"))

  cat("Top identified subgroups:\n")
  print(utils::head(x$data$identification, n))

  if (!is.null(x$data$grf_cuts)) {
    cat("\nTop GRF cuts:\n")
    print(utils::head(x$data$grf_cuts, n))
  }

  if (!is.null(x$data$no_subgroup)) {
    cat("\nNo-subgroup decomposition:\n")
    print(x$data$no_subgroup)
  }

  cat("\nAccess gt-formatted tables via:\n")
  cat("  $identification_summary, $grf_cut_summary, $cut_vs_subgroup_xtab,\n")
  cat("  $no_subgroup_decomposition, $original_agreement, $metrics_table\n")
  if (!is.null(x$plots)) cat("Plots available at $plots.\n")

  invisible(x)
}


# =============================================================================
# INTERNAL: NORMALISATION HELPERS
# =============================================================================

#' Build canonical subgroup label from sg1, sg2
#' @param sg1 Character vector.
#' @param sg2 Character vector.
#' @return Character vector same length as inputs.
#' @keywords internal
#' @noRd
.cv_build_sg_full <- function(sg1, sg2) {
  sg1 <- as.character(sg1)
  sg2 <- as.character(sg2)
  ifelse(
    is.na(sg1) | sg1 == "",
    "(no subgroup)",
    ifelse(
      is.na(sg2) | sg2 == "",
      sg1,
      paste(sg1, sg2, sep = " & ")
    )
  )
}


#' Normalise a single grf_cuts entry for display
#' @param x Character scalar.
#' @return Character scalar with empty/NA replaced by "(no GRF cut)".
#' @keywords internal
#' @noRd
.cv_clean_grf <- function(x) {
  if (is.null(x) || length(x) == 0L) return("(no GRF cut)")
  if (is.na(x) || !nzchar(x)) return("(no GRF cut)")
  x
}


#' Canonical form for subgroup exact-match comparison
#'
#' Sorts factor components and collapses with " & " so that
#' \code{c("a", "b")} and \code{c("b", "a")} both canonicalise to
#' \code{"a & b"}.
#'
#' @param sg1 Character vector.
#' @param sg2 Character vector.
#' @return Character vector of canonical labels.
#' @keywords internal
#' @noRd
.cv_canonicalize_sg <- function(sg1, sg2) {
  sg1 <- as.character(sg1)
  sg2 <- as.character(sg2)
  out <- character(length(sg1))
  for (i in seq_along(sg1)) {
    factors <- c(sg1[i], sg2[i])
    factors <- factors[!is.na(factors) & factors != ""]
    if (length(factors) == 0L) {
      out[i] <- "(no subgroup)"
    } else {
      out[i] <- paste(sort(factors), collapse = " & ")
    }
  }
  out
}


# =============================================================================
# INTERNAL: TABULATION HELPERS
# =============================================================================

#' Aggregate tail rows of a top-n truncated frequency table
#'
#' Replaces rows beyond \code{top_n} with a single "(N others)" row
#' carrying the tail sum.  Preserves the \code{n_folds} total across
#' the table -- users see an honest accounting instead of silent
#' truncation.  Applied to \code{.cv_tabulate_subgroups()} and
#' \code{.cv_tabulate_grf_cuts()} output.
#'
#' @param tab data.frame with columns \code{Rank}, \code{value_col},
#'   \code{n_folds}, \code{pct_folds} (in that order) already sorted
#'   by descending \code{n_folds}.
#' @param top_n Integer.  Maximum number of detail rows to keep.
#' @param value_col Character scalar.  Name of the column holding the
#'   frequency-key (\code{"Subgroup"} or \code{"grf_cut"}).
#' @param total_pairs Integer.  Denominator for the percentage
#'   recomputation of the "(N others)" row.
#' @return data.frame with at most \code{top_n + 1L} rows; the last
#'   row has \code{Rank = top_n + 1L} and the key column labelled
#'   \code{"(N others)"}.
#' @keywords internal
#' @noRd
.cv_aggregate_others <- function(tab, top_n, value_col, total_pairs) {
  head_df <- utils::head(tab, top_n)
  tail_df <- tab[seq.int(top_n + 1L, nrow(tab)), , drop = FALSE]
  n_tail  <- sum(tail_df$n_folds)
  k_tail  <- nrow(tail_df)

  others_row <- data.frame(
    Rank      = top_n + 1L,
    placeholder = sprintf("(%d others)", k_tail),
    n_folds   = n_tail,
    pct_folds = round(100 * n_tail / total_pairs, 2),
    stringsAsFactors = FALSE
  )
  names(others_row)[names(others_row) == "placeholder"] <- value_col

  rbind(head_df, others_row)
}


#' Frequency table of identified subgroups
#' @keywords internal
#' @noRd
.cv_tabulate_subgroups <- function(sg_full, total_pairs, top_n = 15L) {
  tab <- as.data.frame(
    table(Subgroup = sg_full),
    stringsAsFactors = FALSE
  )
  names(tab)[names(tab) == "Freq"] <- "n_folds"
  tab <- tab[order(-tab$n_folds), , drop = FALSE]
  tab$pct_folds <- round(100 * tab$n_folds / total_pairs, 2)
  tab$Rank      <- seq_len(nrow(tab))
  tab <- tab[, c("Rank", "Subgroup", "n_folds", "pct_folds"), drop = FALSE]
  rownames(tab) <- NULL
  if (!is.null(top_n) && nrow(tab) > top_n) {
    tab <- .cv_aggregate_others(tab, top_n, value_col = "Subgroup",
                                total_pairs = total_pairs)
  }
  tab
}


#' Frequency table of raw GRF cut strings
#' @keywords internal
#' @noRd
.cv_tabulate_grf_cuts <- function(grf_clean, total_pairs, top_n = 15L) {
  tab <- as.data.frame(
    table(grf_cut = grf_clean),
    stringsAsFactors = FALSE
  )
  names(tab)[names(tab) == "Freq"] <- "n_folds"
  tab <- tab[order(-tab$n_folds), , drop = FALSE]
  tab$pct_folds <- round(100 * tab$n_folds / total_pairs, 2)
  tab$Rank      <- seq_len(nrow(tab))
  tab <- tab[, c("Rank", "grf_cut", "n_folds", "pct_folds"), drop = FALSE]
  rownames(tab) <- NULL
  if (!is.null(top_n) && nrow(tab) > top_n) {
    tab <- .cv_aggregate_others(tab, top_n, value_col = "grf_cut",
                                total_pairs = total_pairs)
  }
  tab
}


#' Cross-tab of GRF cut x identified subgroup
#' @keywords internal
#' @noRd
.cv_xtab_cut_vs_sg <- function(grf_clean, sg_full, top_n = 15L) {
  tab <- as.data.frame(
    table(grf_cut = grf_clean, subgroup = sg_full),
    stringsAsFactors = FALSE
  )
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  names(tab)[names(tab) == "Freq"] <- "n_folds"
  total <- sum(tab$n_folds)
  tab$pct_folds <- round(100 * tab$n_folds / total, 2)
  tab <- tab[order(-tab$n_folds), , drop = FALSE]
  rownames(tab) <- NULL
  if (!is.null(top_n) && nrow(tab) > top_n) {
    # Two-key xtab: aggregate the tail into a single synthetic row that
    # labels both key columns as "(N others)" so downstream gt formatting
    # stays column-consistent.
    head_df <- utils::head(tab, top_n)
    tail_df <- tab[seq.int(top_n + 1L, nrow(tab)), , drop = FALSE]
    n_tail  <- sum(tail_df$n_folds)
    k_tail  <- nrow(tail_df)
    others_row <- data.frame(
      grf_cut   = sprintf("(%d others)", k_tail),
      subgroup  = sprintf("(%d others)", k_tail),
      n_folds   = n_tail,
      pct_folds = round(100 * n_tail / total, 2),
      stringsAsFactors = FALSE
    )
    tab <- rbind(head_df, others_row)
  }
  tab
}


#' Decompose no-subgroup folds into GRF-none vs consistency-rejected
#' @keywords internal
#' @noRd
.cv_decompose_no_sg <- function(sg_full, grf_clean, total_pairs) {
  is_no_sg <- sg_full == "(no subgroup)"
  n_no_sg  <- sum(is_no_sg)

  if (n_no_sg == 0L) {
    return(data.frame(
      Category     = c(
        "Folds with no subgroup identified",
        "  GRF returned no cut",
        "  GRF returned a cut; consistency rejected"
      ),
      n_folds      = c(0L, 0L, 0L),
      pct_of_total = c(0, 0, 0),
      pct_of_no_sg = c(NA_real_, NA_real_, NA_real_),
      stringsAsFactors = FALSE
    ))
  }

  grf_none         <- grf_clean[is_no_sg] == "(no GRF cut)"
  n_grf_none       <- sum(grf_none, na.rm = TRUE)
  n_consist_reject <- n_no_sg - n_grf_none

  data.frame(
    Category     = c(
      "Folds with no subgroup identified",
      "  GRF returned no cut",
      "  GRF returned a cut; consistency rejected"
    ),
    n_folds      = c(n_no_sg, n_grf_none, n_consist_reject),
    pct_of_total = round(
      100 * c(n_no_sg, n_grf_none, n_consist_reject) / total_pairs, 2
    ),
    pct_of_no_sg = c(
      100,
      round(100 * n_grf_none / n_no_sg, 2),
      round(100 * n_consist_reject / n_no_sg, 2)
    ),
    stringsAsFactors = FALSE
  )
}


#' Distribution of Pcons among identifying folds
#'
#' Summarises the per-fold consistency probability (\code{Pcons}) captured
#' in \code{fold_summary$pconsistency}, restricted to folds that
#' identified a subgroup.  Reveals whether identifying folds are
#' scraping the \code{pconsistency.threshold} (e.g., Pcons values just
#' above 0.90) or are well above it.  Rejected-candidate Pcons values
#' are \emph{not} surfaced here -- the current \code{forestsearch()}
#' return object only carries the selected candidate's Pcons.
#'
#' @param pcons Numeric vector.  \code{fold_summary$pconsistency}.
#' @param sg_full Character vector.  Identified subgroup label per fold.
#' @param total_pairs Integer.  Total sim x fold pairs (denominator for
#'   percentages).
#' @return data.frame with columns \code{Pcons_range}, \code{n_folds},
#'   \code{pct_of_identifying}, plus two rows at the end with median
#'   and IQR of Pcons among identifying folds.  Returns \code{NULL} if
#'   no folds identified a subgroup.
#' @keywords internal
#' @noRd
.cv_pconsistency_distribution <- function(pcons, sg_full, total_pairs) {
  keep <- !is.na(pcons) & sg_full != "(no subgroup)"
  n_identifying <- sum(keep)

  if (n_identifying == 0L) return(NULL)

  p <- pcons[keep]

  # Bin breaks mirror summarize_bootstrap_subgroups() consistency_dist.
  breaks <- c(0, 0.5, 0.7, 0.8, 0.85, 0.95, 1.0)
  labels <- c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85",
              "0.85-0.95", ">=0.95")

  bins <- cut(p, breaks = breaks, labels = labels,
              include.lowest = TRUE, right = FALSE)

  # Table preserves factor level ordering
  tab <- table(bins)
  dist_df <- data.frame(
    Pcons_range       = labels,
    n_folds           = as.integer(tab),
    pct_of_identifying = round(100 * as.integer(tab) / n_identifying, 2),
    stringsAsFactors = FALSE
  )

  # Append summary rows (median + IQR) so the full picture lives in one df
  q <- stats::quantile(p, probs = c(0.25, 0.5, 0.75), na.rm = TRUE,
                        names = FALSE)
  summary_rows <- data.frame(
    Pcons_range = c("Median (identifying)", "IQR (identifying)",
                    "Identifying folds total"),
    n_folds     = c(NA_integer_, NA_integer_, n_identifying),
    pct_of_identifying = c(
      NA_real_, NA_real_,
      round(100 * n_identifying / total_pairs, 2)
    ),
    stringsAsFactors = FALSE
  )
  # Use the Pcons_range column to carry the numeric formatting so the
  # data.frame stays rectangular (character column accepts both labels
  # and formatted numbers).
  summary_rows$Pcons_range[1] <- sprintf("Median (identifying): %.3f",
                                         q[2])
  summary_rows$Pcons_range[2] <- sprintf("IQR (identifying): %.3f - %.3f",
                                         q[1], q[3])

  rbind(dist_df, summary_rows)
}


#' Tidy numeric summary of fold_summary diagnostic columns
#'
#' Produces a tidy data.frame with rows for each numeric diagnostic
#' column present in \code{fold_summary} and columns for n / median /
#' IQR / min / max.  Backward-compatible: absent columns simply produce
#' no rows.  Preferred over proliferating single-column distribution
#' slots as additional diagnostics are surfaced in future phases.
#'
#' @param fsum Data frame.  The \code{fold_summary} data frame as
#'   captured by \code{\link{forestsearch_tenfold}}.
#' @return data.frame with columns \code{Metric}, \code{n}, \code{median},
#'   \code{Q1}, \code{Q3}, \code{min}, \code{max}, \code{n_na},
#'   \code{Context}.  One row per numeric diagnostic column that is
#'   present in \code{fsum}.  Empty data frame if none are present.
#' @keywords internal
#' @noRd
.cv_fold_numeric_summary <- function(fsum) {

  # Metric specifications: name, display label, subset rule, and a
  # human-readable context noting the denominator population.  Rules
  # apply in order and are skipped when the column is absent.
  metric_specs <- list(
    list(
      col     = "n_test",
      label   = "Test-fold size (n_test)",
      filter  = NULL,
      context = "All sim x fold pairs"
    ),
    list(
      col     = "pconsistency",
      label   = "Pcons (identifying folds)",
      filter  = function(fsum) !is.na(fsum$pconsistency),
      context = "Folds that identified a subgroup"
    ),
    list(
      col     = "training_fs_hr",
      label   = "Training-fold subgroup effect (identifying folds)",
      filter  = function(fsum) !is.na(fsum$training_fs_hr),
      context = "Natural scale; in-sample; optimistically biased"
    ),
    list(
      col     = "n_candidates_evaluated",
      label   = "Candidates evaluated per fold",
      filter  = function(fsum) !is.na(fsum$n_candidates_evaluated),
      context = "Folds where consistency stage ran"
    )
  )

  rows <- list()
  for (spec in metric_specs) {
    if (!spec$col %in% names(fsum)) next

    idx <- if (is.null(spec$filter)) seq_len(nrow(fsum))
           else which(spec$filter(fsum))
    values <- fsum[[spec$col]][idx]
    values <- values[!is.na(values)]  # belt-and-braces

    if (length(values) == 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        Metric  = spec$label,
        n       = 0L,
        median  = NA_real_,
        Q1      = NA_real_,
        Q3      = NA_real_,
        min     = NA_real_,
        max     = NA_real_,
        n_na    = nrow(fsum) - if (is.null(spec$filter)) 0L
                               else sum(spec$filter(fsum)),
        Context = spec$context,
        stringsAsFactors = FALSE
      )
      next
    }

    q <- stats::quantile(values, probs = c(0.25, 0.5, 0.75),
                         na.rm = TRUE, names = FALSE)

    n_total <- if (is.null(spec$filter)) nrow(fsum)
               else sum(spec$filter(fsum))
    n_na    <- nrow(fsum) - n_total

    rows[[length(rows) + 1L]] <- data.frame(
      Metric  = spec$label,
      n       = length(values),
      median  = q[2],
      Q1      = q[1],
      Q3      = q[3],
      min     = min(values),
      max     = max(values),
      n_na    = n_na,
      Context = spec$context,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    return(data.frame(
      Metric  = character(0),
      n       = integer(0),
      median  = numeric(0),
      Q1      = numeric(0),
      Q3      = numeric(0),
      min     = numeric(0),
      max     = numeric(0),
      n_na    = integer(0),
      Context = character(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
}


#' Agreement with original full-data analysis
#' @keywords internal
#' @noRd
.cv_original_agreement <- function(fsum,
                                   original_sg       = NULL,
                                   original_grf_cuts = NULL,
                                   has_grf_cuts      = FALSE) {

  if (is.null(original_sg) && is.null(original_grf_cuts)) return(NULL)

  total <- nrow(fsum)
  rows  <- list(
    data.frame(
      Metric = "Total sim x fold pairs",
      Value  = as.character(total),
      stringsAsFactors = FALSE
    )
  )

  # Subgroup-definition agreement
  if (!is.null(original_sg)) {
    orig_sg_char <- as.character(original_sg)
    orig_sg_char <- orig_sg_char[!is.na(orig_sg_char) & orig_sg_char != ""]

    if (length(orig_sg_char) > 0L) {
      orig_canon <- paste(sort(orig_sg_char), collapse = " & ")
      fold_canon <- .cv_canonicalize_sg(fsum$sg1, fsum$sg2)

      exact <- sum(fold_canon == orig_canon, na.rm = TRUE)

      partial <- sum(vapply(seq_len(total), function(i) {
        fold_factors <- c(fsum$sg1[i], fsum$sg2[i])
        fold_factors <- fold_factors[!is.na(fold_factors) &
                                       fold_factors != ""]
        if (length(fold_factors) == 0L) return(FALSE)
        any(fold_factors %in% orig_sg_char)
      }, logical(1)), na.rm = TRUE)

      rows <- c(rows, list(
        data.frame(
          Metric = "Original subgroup definition",
          Value  = orig_canon,
          stringsAsFactors = FALSE
        ),
        data.frame(
          Metric = "Exact subgroup match",
          Value  = sprintf("%d (%.1f%%)", exact, 100 * exact / total),
          stringsAsFactors = FALSE
        ),
        data.frame(
          Metric = "Partial match (>=1 shared factor)",
          Value  = sprintf("%d (%.1f%%)", partial, 100 * partial / total),
          stringsAsFactors = FALSE
        )
      ))
    }
  }

  # GRF-cut agreement
  if (!is.null(original_grf_cuts) && has_grf_cuts) {
    orig_grf_char <- as.character(original_grf_cuts)
    orig_grf_char <- orig_grf_char[!is.na(orig_grf_char) &
                                     orig_grf_char != ""]

    if (length(orig_grf_char) > 0L) {
      orig_grf_str <- paste(orig_grf_char, collapse = " | ")
      grf_exact    <- sum(fsum$grf_cuts == orig_grf_str, na.rm = TRUE)

      rows <- c(rows, list(
        data.frame(
          Metric = "Original GRF cut(s)",
          Value  = orig_grf_str,
          stringsAsFactors = FALSE
        ),
        data.frame(
          Metric = "Exact GRF-cut match",
          Value  = sprintf("%d (%.1f%%)", grf_exact,
                           100 * grf_exact / total),
          stringsAsFactors = FALSE
        )
      ))

      # Both-match (subgroup AND GRF cut)
      if (!is.null(original_sg)) {
        orig_sg_char <- as.character(original_sg)
        orig_sg_char <- orig_sg_char[!is.na(orig_sg_char) &
                                       orig_sg_char != ""]
        if (length(orig_sg_char) > 0L) {
          orig_canon <- paste(sort(orig_sg_char), collapse = " & ")
          fold_canon <- .cv_canonicalize_sg(fsum$sg1, fsum$sg2)
          both       <- sum(
            fold_canon == orig_canon & fsum$grf_cuts == orig_grf_str,
            na.rm = TRUE
          )
          rows <- c(rows, list(data.frame(
            Metric = "Both subgroup and GRF cut match",
            Value  = sprintf("%d (%.1f%%)", both, 100 * both / total),
            stringsAsFactors = FALSE
          )))
        }
      }
    }
  }

  do.call(rbind, rows)
}


#' Safe delegation to cv_metrics_tables()
#' @keywords internal
#' @noRd
.cv_metrics_table_safe <- function(cv_output, original_sg = NULL) {
  tryCatch(
    cv_metrics_tables(
      cv_result     = cv_output,
      sg_definition = original_sg,
      table_style   = "combined",
      use_gt        = TRUE
    ),
    error = function(e) {
      warning("Could not build metrics_table: ", e$message, call. = FALSE)
      NULL
    }
  )
}


# =============================================================================
# INTERNAL: gt FORMATTERS
# =============================================================================

#' @keywords internal
#' @noRd
.format_cv_identification_gt <- function(df, total_pairs, n_sims, n_folds) {
  gt::gt(df) |>
    gt::tab_header(
      title    = gt::md("**Identified subgroups across sim x fold pairs**"),
      subtitle = sprintf(
        "%d simulations x %d folds = %d total pairs",
        n_sims, n_folds, total_pairs
      )
    ) |>
    gt::cols_label(
      Rank      = "Rank",
      Subgroup  = "Identified subgroup",
      n_folds   = "n folds",
      pct_folds = "% folds"
    ) |>
    gt::fmt_number(columns = "pct_folds", decimals = 2) |>
    gt::tab_options(table.font.size = gt::px(13))
}


#' @keywords internal
#' @noRd
.format_cv_grf_cut_gt <- function(df, total_pairs, n_sims, n_folds) {
  gt::gt(df) |>
    gt::tab_header(
      title    = gt::md("**GRF cuts across training sets**"),
      subtitle = sprintf(
        "Raw grf_cuts strings; %d simulations x %d folds = %d total pairs",
        n_sims, n_folds, total_pairs
      )
    ) |>
    gt::cols_label(
      Rank      = "Rank",
      grf_cut   = "GRF cuts",
      n_folds   = "n folds",
      pct_folds = "% folds"
    ) |>
    gt::fmt_number(columns = "pct_folds", decimals = 2) |>
    gt::tab_options(table.font.size = gt::px(13))
}


#' @keywords internal
#' @noRd
.format_cv_xtab_gt <- function(df, total_pairs) {
  gt::gt(df) |>
    gt::tab_header(
      title    = gt::md("**GRF cut x identified subgroup**"),
      subtitle = sprintf(
        "Cross-tabulation across %d sim x fold pairs", total_pairs
      )
    ) |>
    gt::cols_label(
      grf_cut   = "GRF cuts",
      subgroup  = "Identified subgroup",
      n_folds   = "n folds",
      pct_folds = "% folds"
    ) |>
    gt::fmt_number(columns = "pct_folds", decimals = 2) |>
    gt::tab_source_note(
      paste0(
        "Rows where the identified subgroup does not align with GRF's ",
        "output reveal divergence between candidate discovery and ",
        "consistency evaluation."
      )
    ) |>
    gt::tab_options(table.font.size = gt::px(13))
}


#' @keywords internal
#' @noRd
.format_cv_no_sg_gt <- function(df, total_pairs) {
  gt::gt(df) |>
    gt::tab_header(
      title    = gt::md("**Decomposition of no-subgroup folds**"),
      subtitle = sprintf(
        "Out of %d sim x fold pairs", total_pairs
      )
    ) |>
    gt::cols_label(
      Category     = "Category",
      n_folds      = "n folds",
      pct_of_total = "% of all pairs",
      pct_of_no_sg = "% of no-sg folds"
    ) |>
    gt::fmt_number(columns = c("pct_of_total", "pct_of_no_sg"), decimals = 2) |>
    gt::sub_missing(columns = "pct_of_no_sg", missing_text = "--") |>
    gt::tab_source_note(
      paste0(
        "Guidance: a high 'GRF returned no cut' share argues for relaxing ",
        "dmin.grf; a high 'consistency rejected' share argues for relaxing ",
        "pconsistency.threshold."
      )
    ) |>
    gt::tab_options(table.font.size = gt::px(13))
}


#' @keywords internal
#' @noRd
.format_cv_pconsistency_gt <- function(df, total_pairs) {
  # The last three rows carry summary statistics (median, IQR, total);
  # the first six carry the binned distribution.  Highlight the summary
  # rows so readers orient quickly.
  is_summary <- grepl("^(Median|IQR|Identifying)", df$Pcons_range)

  gt::gt(df) |>
    gt::tab_header(
      title    = gt::md(
        "**Pcons distribution among identifying folds**"
      ),
      subtitle = sprintf(
        "Out of %d sim x fold pairs", total_pairs
      )
    ) |>
    gt::cols_label(
      Pcons_range        = "Pcons range",
      n_folds            = "n folds",
      pct_of_identifying = "% of identifying folds"
    ) |>
    gt::fmt_number(columns = "pct_of_identifying", decimals = 2) |>
    gt::sub_missing(
      columns     = c("n_folds", "pct_of_identifying"),
      missing_text = "--"
    ) |>
    gt::tab_style(
      style    = list(
        gt::cell_fill(color = "#f0f0f0"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(rows = is_summary)
    ) |>
    gt::tab_source_note(
      paste0(
        "Restricted to folds that identified a subgroup.  A concentration ",
        "near the pconsistency.threshold (typically 0.90) suggests many ",
        "identifications are marginal; concentration above 0.95 suggests ",
        "robust identifications.  Values for rejected candidates are not ",
        "surfaced -- only the selected candidate's Pcons is captured."
      )
    ) |>
    gt::tab_options(table.font.size = gt::px(13))
}


#' @keywords internal
#' @noRd
.format_cv_fold_numeric_summary_gt <- function(df, total_pairs) {
  gt::gt(df) |>
    gt::tab_header(
      title    = gt::md("**Fold-level numeric diagnostics**"),
      subtitle = sprintf(
        "Across %d sim x fold pairs; per-metric denominator shown in Context",
        total_pairs
      )
    ) |>
    gt::cols_label(
      Metric  = "Metric",
      n       = "n",
      median  = "Median",
      Q1      = "Q1",
      Q3      = "Q3",
      min     = "Min",
      max     = "Max",
      n_na    = "n NA",
      Context = "Context"
    ) |>
    gt::fmt_number(
      columns  = c("median", "Q1", "Q3", "min", "max"),
      decimals = 3
    ) |>
    gt::sub_missing(
      columns      = c("median", "Q1", "Q3", "min", "max"),
      missing_text = "--"
    ) |>
    gt::tab_source_note(
      paste0(
        "Training-fold subgroup effects are in-sample estimates on the ",
        "effect measure's natural scale (HR for survival; OR, RR, IRR, ",
        "RD, IRD, or MD for GLM).  They are optimistically biased and ",
        "shown here for diagnostic comparison across folds only; do ",
        "not use for inference."
      )
    ) |>
    gt::tab_options(table.font.size = gt::px(13))
}


#' @keywords internal
#' @noRd
.format_cv_original_agreement_gt <- function(df) {
  gt::gt(df) |>
    gt::tab_header(
      title    = gt::md("**Agreement with original full-data analysis**"),
      subtitle = "Fraction of sim x fold pairs matching the original result"
    ) |>
    gt::cols_label(
      Metric = "",
      Value  = "Value"
    ) |>
    gt::tab_style(
      style    = list(
        gt::cell_fill(color = "#e8f4f8"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = "Metric",
        rows    = grepl("^(Exact|Partial|Both)", df$Metric)
      )
    ) |>
    gt::tab_options(table.font.size = gt::px(13))
}


# =============================================================================
# INTERNAL: PLOT HELPERS
# =============================================================================

#' Construct ggplot diagnostics for CV summary
#' @keywords internal
#' @noRd
.cv_diagnostic_plots <- function(fsum, has_grf_cuts, top_n = 10L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' required for CV diagnostic plots; skipping.")
    return(NULL)
  }

  plots <- list()

  # ---- identified-subgroup bar chart ----
  sg_tab <- as.data.frame(
    table(Subgroup = fsum$sg_full),
    stringsAsFactors = FALSE
  )
  names(sg_tab)[names(sg_tab) == "Freq"] <- "n_folds"
  sg_tab   <- sg_tab[order(-sg_tab$n_folds), , drop = FALSE]
  sg_tab   <- utils::head(sg_tab, top_n)
  sg_tab$Subgroup <- factor(sg_tab$Subgroup, levels = rev(sg_tab$Subgroup))

  plots$identification <- ggplot2::ggplot(
    sg_tab,
    ggplot2::aes(x = .data$Subgroup, y = .data$n_folds)
  ) +
    ggplot2::geom_col(fill = "#2E86AB", alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = sprintf(
        "Top %d identified subgroups (sim x fold)",
        nrow(sg_tab)
      ),
      x = NULL,
      y = "n folds"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )

  # ---- GRF-cut bar chart ----
  if (has_grf_cuts) {
    grf_tab <- as.data.frame(
      table(grf_cut = fsum$grf_clean),
      stringsAsFactors = FALSE
    )
    names(grf_tab)[names(grf_tab) == "Freq"] <- "n_folds"
    grf_tab <- grf_tab[order(-grf_tab$n_folds), , drop = FALSE]
    grf_tab <- utils::head(grf_tab, top_n)
    grf_tab$grf_cut <- factor(grf_tab$grf_cut, levels = rev(grf_tab$grf_cut))

    plots$grf_cut <- ggplot2::ggplot(
      grf_tab,
      ggplot2::aes(x = .data$grf_cut, y = .data$n_folds)
    ) +
      ggplot2::geom_col(fill = "#A23B72", alpha = 0.8) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = sprintf(
          "Top %d GRF cuts across training sets",
          nrow(grf_tab)
        ),
        x = NULL,
        y = "n folds"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12)
      )
  }

  plots
}
