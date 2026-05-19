# =============================================================================
# compute_oc_score_grid.R --- Score a configuration table under a grid
#                              of weight vectors
# =============================================================================
# Location: quarto/_utils/compute_oc_score_grid.R
# Dependencies: compute_oc_score.R (sourced separately).
#
# Computes the composite OC score under several named weightings in
# one call and appends one column per weighting.  Intended for
# weight-sensitivity / robustness analysis: a configuration that
# ranks near the top across every plausible weighting is robustly
# preferred; a configuration whose rank swings sharply between
# weightings is an artefact of the weight choice.
#
# Sourcing from a .qmd setup chunk:
#
#   source("../../../_utils/compute_oc_score.R")
#   source("../../../_utils/compute_oc_score_grid.R")
#
#   display_dt <- compute_oc_score_grid(display_dt)        # default grid
#   # Adds S_equal, S_bp_default, S_spec_emphasis, ...
#
# This file is NOT part of the forestsearch package; it is an
# exploratory utility shared across simulation evaluation documents.
# =============================================================================


# Default sensitivity-analysis weight grid for the canonical metric set
# (Detection, Sen, Spec).  Each named entry is a length-3 non-negative
# numeric vector; weights are auto-normalised inside compute_oc_score().
# - "equal"           : neutral benchmark
# - "bp_default"      : the benefit-protection default (Spec-heaviest)
# - "spec_emphasis"   : pushes further toward Spec
# - "detect_emphasis" : pushes toward Detection
# - "sen_emphasis"    : pushes toward Sen (stress-test the bp_default)
.OC_WEIGHT_GRID_DEFAULT <- list(
  equal           = c(1.00, 1.00, 1.00),
  bp_default      = c(1.50, 1.00, 1.75),
  spec_emphasis   = c(1.00, 0.50, 2.50),
  detect_emphasis = c(2.00, 1.00, 1.00),
  sen_emphasis    = c(1.00, 2.00, 1.00)
)


# compute_oc_score_grid() ---------------------------------------------------
#
# Arguments
#   data           data.frame / data.table of OC summary rows.
#   weight_grid    named list of non-negative numeric vectors, each of
#                  length(metrics).  Each list element is one weighting
#                  scheme; the entry name becomes the score-column
#                  suffix.  Default: see .OC_WEIGHT_GRID_DEFAULT
#                  above, which contains five weightings spanning
#                  equal weighting through directional emphases on
#                  each of {Detection, Sen, Spec}.
#   method         "linear" or "geometric".  Passed to
#                  compute_oc_score() for each weighting in the grid.
#                  Default "linear".
#   metrics        character vector of OC metric column names.
#                  Default c("Detection","Sen","Spec").
#   score_prefix   prefix for the added score columns.  The full
#                  column name is paste0(score_prefix, names(weight_grid)).
#                  Default "S_".
#
# Value
#   The input `data` with length(weight_grid) new score columns
#   appended.  Each column carries a `score_meta` attribute (set by
#   compute_oc_score()) recording method, metrics, and normalised
#   weights for that weighting.
#
# Edge cases
#   * If a column with the same target name already exists, it is
#     overwritten with a message.
#   * If any individual scoring call fails, the function stops with
#     a clear identification of which weighting failed.
#
# Examples (display_dt is built by an OC summary article):
#
#   # Default grid: five weightings spanning equal -> directional
#   display_dt <- compute_oc_score_grid(display_dt)
#   names(display_dt)   # ... S_equal, S_bp_default, ...
#
#   # Custom two-element grid
#   display_dt <- compute_oc_score_grid(
#     display_dt,
#     weight_grid = list(
#       conservative = c(1, 1, 3),       # very Spec-heavy
#       aggressive   = c(3, 1, 1)        # very Detection-heavy
#     )
#   )
# ---------------------------------------------------------------------------

compute_oc_score_grid <- function(data,
                                  weight_grid  = NULL,
                                  method       = c("linear", "geometric"),
                                  metrics      = c("Detection", "Sen", "Spec"),
                                  score_prefix = "S_") {

  # ---- Resolve & validate arguments --------------------------------------
  method <- match.arg(method)

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'data' has zero rows; nothing to score.", call. = FALSE)
  }

  if (!exists("compute_oc_score", mode = "function")) {
    stop("compute_oc_score() not found.  Source ",
         "compute_oc_score.R before this file.", call. = FALSE)
  }

  if (is.null(weight_grid)) {
    weight_grid <- .OC_WEIGHT_GRID_DEFAULT
  }
  if (!is.list(weight_grid) || length(weight_grid) == 0L) {
    stop("'weight_grid' must be a non-empty named list of numeric vectors.",
         call. = FALSE)
  }
  if (is.null(names(weight_grid)) || any(names(weight_grid) == "") ||
      anyDuplicated(names(weight_grid))) {
    stop("'weight_grid' must have unique, non-empty names ",
         "(used as score-column suffixes).", call. = FALSE)
  }
  bad <- vapply(weight_grid, function(w) {
    !is.numeric(w) || length(w) != length(metrics) || anyNA(w) || any(w < 0)
  }, logical(1))
  if (any(bad)) {
    stop(sprintf(
      "weight_grid entries must be non-negative numeric vectors of length(metrics) = %d.  ",
      length(metrics)),
      "Offending entries: ",
      paste(shQuote(names(weight_grid)[bad]), collapse = ", "),
      ".", call. = FALSE)
  }

  if (!is.character(score_prefix) || length(score_prefix) != 1L ||
      is.na(score_prefix)) {
    stop("'score_prefix' must be a single character string.", call. = FALSE)
  }

  # ---- Loop over the grid ------------------------------------------------
  out <- data
  for (nm in names(weight_grid)) {
    col_name <- paste0(score_prefix, nm)
    if (col_name %in% names(out)) {
      message(sprintf("Overwriting existing column '%s'.", col_name))
    }
    out <- tryCatch(
      compute_oc_score(
        out,
        method    = method,
        metrics   = metrics,
        weights   = weight_grid[[nm]],
        score_col = col_name
      ),
      error = function(e) {
        stop(sprintf(
          "Scoring failed for weighting '%s' (weights = %s): %s",
          nm, paste(weight_grid[[nm]], collapse = ", "),
          conditionMessage(e)), call. = FALSE)
      }
    )
  }
  out
}
