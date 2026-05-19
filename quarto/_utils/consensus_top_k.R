# =============================================================================
# consensus_top_k.R --- Flag configurations that are in the top-k across
#                        all supplied score columns
# =============================================================================
# Location: quarto/_utils/consensus_top_k.R
# Dependencies: none beyond base R.
#
# Adds a logical column to a data.frame / data.table of OC summary rows
# indicating which rows rank in the top-k under EVERY supplied score
# column.  Intended use: after running compute_oc_score_grid() to add
# multiple weighting-based scores, this utility identifies the
# configurations that survive a "robustness across weightings" filter --
# configurations whose preference does not depend on the specific
# weighting choice.
#
# Sourcing from a .qmd setup chunk:
#
#   source("../../../_utils/consensus_top_k.R")
#   display_dt <- consensus_top_k(
#     display_dt,
#     score_cols = c("S_equal", "S_bp_default", "S_spec_emphasis",
#                    "S_detect_emphasis", "S_sen_emphasis"),
#     k = 5
#   )
#   # adds display_dt$consensus_topK = TRUE for rows top-5 under ALL cols
#
# This file is NOT part of the forestsearch package; it is an
# exploratory utility shared across simulation evaluation documents.
# =============================================================================


# consensus_top_k() ---------------------------------------------------------
#
# Arguments
#   data           data.frame / data.table of OC summary rows.
#   score_cols     character vector (length >= 1) of numeric score
#                  columns to evaluate.  Each is ranked descending
#                  (higher is better); a row is in the consensus
#                  top-k if its rank is <= k under every column.
#   k              positive integer.  Top-k cutoff.  Default 5.
#   flag_col       character scalar: name of the logical column to
#                  add.  Default "consensus_topK".
#   ties           tie-breaking method passed to base::rank().  One
#                  of "min", "average", "first", "last", "random",
#                  "max".  Default "min" (best-case for tied scores
#                  -- conservative against false exclusion from the
#                  consensus).
#   add_rank_cols  logical.  If TRUE, also adds one rank column per
#                  score column (named paste0("rank_", score_col)).
#                  Default FALSE.
#
# Value
#   The input `data` with the consensus flag column appended (and
#   optionally rank columns).  The flag column carries a `consensus_meta`
#   attribute recording the score columns, k, and tie-breaking method.
#
# Edge cases
#   * If a score column contains NA, that row's rank is NA, and it is
#     flagged FALSE (a row that cannot be ranked under one weighting
#     cannot be in the consensus).
#   * If k >= nrow(data), every fully-ranked row is in the consensus
#     under every column trivially.
#
# Examples (display_dt has been augmented with grid scores):
#
#   # Top-5 under all five default weightings:
#   display_dt <- consensus_top_k(
#     display_dt,
#     score_cols = paste0("S_", c("equal","bp_default","spec_emphasis",
#                                 "detect_emphasis","sen_emphasis")),
#     k = 5
#   )
#   sum(display_dt$consensus_topK)
#
#   # Stricter: top-3, and include the rank columns for inspection
#   display_dt <- consensus_top_k(
#     display_dt,
#     score_cols    = paste0("S_", c("equal","bp_default","spec_emphasis")),
#     k             = 3,
#     flag_col      = "consensus_top3",
#     add_rank_cols = TRUE
#   )
# ---------------------------------------------------------------------------

consensus_top_k <- function(data,
                            score_cols,
                            k             = 5L,
                            flag_col      = "consensus_topK",
                            ties          = c("min", "average", "first",
                                              "last", "random", "max"),
                            add_rank_cols = FALSE) {

  # ---- Input validation --------------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'data' has zero rows; nothing to rank.", call. = FALSE)
  }
  if (!is.character(score_cols) || length(score_cols) < 1L ||
      anyNA(score_cols) || any(nchar(score_cols) == 0L)) {
    stop("'score_cols' must be a non-empty character vector with no NA / empty entries.",
         call. = FALSE)
  }
  missing_cols <- setdiff(score_cols, names(data))
  if (length(missing_cols)) {
    stop(sprintf(
      "Score column(s) not found in 'data': %s.  Available columns: %s.",
      paste(shQuote(missing_cols), collapse = ", "),
      paste(shQuote(names(data)),  collapse = ", ")),
      call. = FALSE)
  }
  for (col in score_cols) {
    if (!is.numeric(data[[col]])) {
      stop(sprintf("Score column '%s' must be numeric; got %s.",
                   col, paste(class(data[[col]]), collapse = "/")),
           call. = FALSE)
    }
  }

  if (!is.numeric(k) || length(k) != 1L || is.na(k) ||
      k < 1 || k != as.integer(k)) {
    stop("'k' must be a positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  if (!is.character(flag_col) || length(flag_col) != 1L ||
      is.na(flag_col) || nchar(flag_col) == 0L) {
    stop("'flag_col' must be a single non-empty character string.",
         call. = FALSE)
  }
  ties <- match.arg(ties)

  if (!is.logical(add_rank_cols) || length(add_rank_cols) != 1L ||
      is.na(add_rank_cols)) {
    stop("'add_rank_cols' must be TRUE or FALSE.", call. = FALSE)
  }

  # ---- Compute ranks (higher score -> better rank, so negate) ------------
  # Per-column rank, NA-aware.  rank(-x, na.last = "keep") gives NA
  # rank for NA entries, which then fail the "rank <= k" test.
  rank_mat <- vapply(score_cols, function(col) {
    base::rank(-data[[col]], ties.method = ties, na.last = "keep")
  }, numeric(nrow(data)))

  # Worst (largest) rank across columns; NA propagates.
  worst_rank <- apply(rank_mat, 1, max)

  is_consensus <- !is.na(worst_rank) & worst_rank <= k

  # ---- Attach provenance and assemble result -----------------------------
  attr(is_consensus, "consensus_meta") <- list(
    score_cols = score_cols,
    k          = k,
    ties       = ties
  )

  result <- data
  result[[flag_col]] <- is_consensus

  if (add_rank_cols) {
    for (j in seq_along(score_cols)) {
      result[[paste0("rank_", score_cols[j])]] <- rank_mat[, j]
    }
  }
  result
}
