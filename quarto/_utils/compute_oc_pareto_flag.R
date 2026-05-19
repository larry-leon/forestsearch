# =============================================================================
# compute_oc_pareto_flag.R --- Pareto frontier flag for OC summary tables
# =============================================================================
# Location: quarto/_utils/compute_oc_pareto_flag.R
# Dependencies: none beyond base R.
#
# Adds a logical column to a data.frame / data.table of OC summary rows
# indicating which rows lie on the Pareto frontier of a user-specified
# set of metrics.  By default the frontier is computed over
# {Detection, Sen, Spec}, all higher-is-better; the metrics argument
# is configurable and `higher_better` can be set per-metric for
# mixed-direction frontiers.
#
# Dominance is defined in the standard weak-dominance form:
#
#   p_j dominates p_i  iff  p_j[k] >= p_i[k] for all metrics k
#                           AND p_j[k] >  p_i[k] for at least one k
#
# A row is on the Pareto frontier when no other row dominates it.
# Tied configurations are both on the frontier.
#
# Sourcing from a .qmd setup chunk:
#
#   source("../../../_utils/compute_oc_pareto_flag.R")   # adjust depth
#   display_dt <- compute_oc_pareto_flag(display_dt)
#   table(display_dt$pareto_frontier)
#
# This file is NOT part of the forestsearch package; it is an
# exploratory utility shared across simulation evaluation documents.
#
# Note: this function is structurally identical to the Pareto routine
# used inside forestsearch (compute_pareto_frontier() in
# subgroup_consistency_helpers.R), but operates on user-chosen OC
# columns rather than the package-internal (hr, N) coordinates.  The
# two functions are deliberately separate: this one supports an
# arbitrary metric set with per-axis direction; the package version
# is specialised to subgroup-selection coordinates.
# =============================================================================


# compute_oc_pareto_flag() --------------------------------------------------
#
# Arguments
#   data           data.frame / data.table of OC summary rows.
#   metrics        character vector (length >= 2) of column names that
#                  define the Pareto axes.  Default
#                  c("Detection", "Sen", "Spec").
#   flag_col       character scalar; name of the logical column to add.
#                  Default "pareto_frontier".
#   higher_better  logical scalar or vector of length(metrics).  Per-
#                  metric direction-of-better flag.  TRUE means
#                  larger values are preferred (default); FALSE means
#                  smaller values are preferred (the function flips
#                  the sign internally).  Recycled to length(metrics)
#                  if scalar.  Default TRUE.
#
# Value
#   The input `data` with a new logical column (named by `flag_col`)
#   appended.  TRUE = on the Pareto frontier; FALSE = dominated;
#   FALSE also when the row has NA on any frontier metric (since
#   dominance cannot be evaluated).  The column carries an attribute
#   `pareto_meta` recording the metrics and direction flags used.
#
# Complexity
#   O(n^2) in the number of rows.  For OC summary tables with ~20
#   rows this is trivially fast and the code stays transparent.
#
# Edge cases
#   * Identical configurations (matching on every metric) are all
#     on the frontier under the weak-dominance definition.
#   * NA on any frontier metric -> row flagged FALSE (cannot evaluate).
#   * A single row is trivially on its own frontier (flagged TRUE).
#
# Examples (display_dt is built by an OC summary article):
#
#   # Default: Pareto on {Detection, Sen, Spec}, all higher-is-better
#   display_dt <- compute_oc_pareto_flag(display_dt)
#   sum(display_dt$pareto_frontier)   # how many frontier configs?
#
#   # Frontier over a different metric set
#   display_dt <- compute_oc_pareto_flag(
#     display_dt,
#     metrics  = c("Detection", "Sen", "Spec", "PPV"),
#     flag_col = "pareto_frontier_4d"
#   )
#
#   # Mixed direction (e.g. include a "lower is better" metric)
#   display_dt <- compute_oc_pareto_flag(
#     display_dt,
#     metrics       = c("Detection", "Sen", "Spec", "Trim_count"),
#     higher_better = c(TRUE, TRUE, TRUE, FALSE)
#   )
# ---------------------------------------------------------------------------

compute_oc_pareto_flag <- function(data,
                                   metrics       = c("Detection", "Sen", "Spec"),
                                   flag_col      = "pareto_frontier",
                                   higher_better = TRUE) {

  # ---- Input validation --------------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'data' has zero rows; nothing to evaluate.", call. = FALSE)
  }
  if (!is.character(metrics) || length(metrics) < 2L ||
      any(is.na(metrics)) || any(nchar(metrics) == 0L)) {
    stop("'metrics' must be a character vector of length >= 2 with no NA / empty entries.",
         call. = FALSE)
  }
  if (!is.character(flag_col) || length(flag_col) != 1L ||
      is.na(flag_col) || nchar(flag_col) == 0L) {
    stop("'flag_col' must be a single non-empty character string.",
         call. = FALSE)
  }
  if (!is.logical(higher_better) || anyNA(higher_better)) {
    stop("'higher_better' must be logical (TRUE/FALSE), no NAs.",
         call. = FALSE)
  }
  if (length(higher_better) == 1L) {
    higher_better <- rep(higher_better, length(metrics))
  } else if (length(higher_better) != length(metrics)) {
    stop(sprintf(
      "'higher_better' must be length 1 (recycled) or length(metrics) = %d; got %d.",
      length(metrics), length(higher_better)), call. = FALSE)
  }

  missing_cols <- setdiff(metrics, names(data))
  if (length(missing_cols)) {
    stop(sprintf(
      "Metric column(s) not found in 'data': %s.  Available columns: %s.",
      paste(shQuote(missing_cols), collapse = ", "),
      paste(shQuote(names(data)),  collapse = ", ")),
      call. = FALSE)
  }
  for (col in metrics) {
    if (!is.numeric(data[[col]])) {
      stop(sprintf("Metric column '%s' must be numeric; got %s.",
                   col, paste(class(data[[col]]), collapse = "/")),
           call. = FALSE)
    }
  }

  # ---- Build the metric matrix in canonical "higher is better" form ------
  # Build column-by-column via `[[` so the extraction works identically
  # on data.frame and data.table inputs (the `data[, metrics]` form
  # triggers data.table's NSE and fails when metrics is a variable).
  M <- do.call(cbind, lapply(metrics, function(col) data[[col]]))
  colnames(M) <- metrics
  for (k in seq_along(metrics)) {
    if (!higher_better[k]) M[, k] <- -M[, k]
  }

  # ---- Identify rows that cannot be evaluated (any NA on a metric) -------
  has_na <- apply(M, 1, anyNA)

  # ---- Pairwise dominance scan -------------------------------------------
  # O(n^2).  For the small n typical of OC tables (~20 rows) this is
  # fine and keeps the code transparent.  Each row starts assumed-on-
  # frontier; finding any dominator flips it to FALSE.
  n  <- nrow(M)
  is_frontier <- !has_na   # NA rows pre-emptively FALSE

  for (i in seq_len(n)) {
    if (!is_frontier[i]) next        # already excluded
    for (j in seq_len(n)) {
      if (i == j || has_na[j]) next  # skip self and NA-row comparisons
      diff_ji <- M[j, ] - M[i, ]
      # j weakly dominates i: all diff_ji >= 0, and at least one > 0
      if (all(diff_ji >= 0) && any(diff_ji > 0)) {
        is_frontier[i] <- FALSE
        break
      }
    }
  }

  # Attach provenance and return
  attr(is_frontier, "pareto_meta") <- list(
    metrics       = metrics,
    higher_better = higher_better
  )

  result <- data
  result[[flag_col]] <- is_frontier
  result
}
