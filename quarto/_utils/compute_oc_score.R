# =============================================================================
# compute_oc_score.R --- Composite scalar scores for OC summary tables
# =============================================================================
# Location: quarto/_utils/compute_oc_score.R
# Dependencies: none beyond base R.
#
# Adds a scalar composite-score column to a data.frame / data.table of
# OC summary rows.  Supports three scoring rules:
#
#   * "linear"    : S = sum_k w_k * x_k                 (weighted sum)
#   * "geometric" : S = prod_k x_k^{w_k}                (weighted geometric mean)
#   * "youden"    : S = x_1 + x_2 - 1                   (Youden's J)
#
# Linear and geometric methods accept any number of metric columns
# (length >= 2) and a corresponding weight vector; weights are
# auto-normalised to sum to 1.  Youden's J is a fixed two-component
# formula (default Sen + Spec - 1) and ignores any supplied weights.
#
# Default weights ("benefit-protection" weighting)
# --------------------------------------------------
# When both 'metrics' and 'weights' are left at their defaults (i.e.
# the canonical {Detection, Sen, Spec} under method = "linear" or
# "geometric"), weights default to (1.5, 1, 1.75) -- the largest
# weight on Spec.  The rationale: in the harm-subgroup framing of
# forestsearch, high Spec on the harm class corresponds to low
# false-positive harm calls, which is equivalent to high *Sen* for
# the complementary benefit class -- i.e. capturing as many true
# beneficiaries as possible in the treatment-recommendation set
# Hc = complement of the identified harm subgroup.  Weighting Spec
# above Sen therefore reflects a "protect the beneficiaries"
# decision rule that is the operational goal of harm-subgroup
# discovery in clinical-trial subgroup analysis.  The direction
# also matches standard clinical-trial conservatism around false
# positives.  If 'metrics' is changed, weights default to equal
# weights (1/k each) unless the user supplies their own.
#
# Sourcing from a .qmd setup chunk:
#
#   source("../../../_utils/compute_oc_score.R")   # adjust depth
#
#   # Default: benefit-protection weighting on {Detection, Sen, Spec}
#   display_dt <- compute_oc_score(display_dt)
#
#   # Custom weights:
#   display_dt <- compute_oc_score(display_dt,
#                                  weights = c(0.4, 0.3, 0.3))
#
# This file is NOT part of the forestsearch package; it is an
# exploratory utility shared across simulation evaluation documents.
# =============================================================================


# compute_oc_score() --------------------------------------------------------
#
# Arguments
#   data       data.frame / data.table of OC summary rows.
#   method     character scalar, one of:
#              "linear"    -- weighted linear sum
#              "geometric" -- weighted geometric mean
#              "youden"    -- Youden's J statistic (Sen + Spec - 1)
#              Default "linear".
#   metrics    character vector of column names in `data` whose values
#              are combined into the score.  Default depends on method:
#              c("Detection","Sen","Spec") for linear/geometric;
#              c("Sen","Spec") fixed for youden.
#   weights   non-negative numeric vector, same length as `metrics`,
#              giving the relative weight of each metric.  NULL means
#              the benefit-protection default c(1.5, 1, 1.75) when
#              metrics is also defaulted (canonical setting); equal
#              weights otherwise.  See file header for rationale.
#              Auto-normalised to sum to 1.  Ignored for method =
#              "youden".
#   score_col character scalar, name of the column to add to `data`.
#              Default "S".
#
# Value
#   The input `data` with the score column appended.  The score
#   column carries an attribute `score_meta` listing the method, the
#   metrics used, and the normalised weights for transparency.
#
# Edge cases
#   * If any metric value is NA, the score for that row is NA.
#   * For "geometric", any zero among the metric values forces S = 0
#     for that row (standard geometric-mean behaviour).  Negative
#     metric values are rejected (geometric mean undefined).
#   * For "youden", `metrics` must have length 2; defaults to
#     c("Sen","Spec").  `weights` is ignored (set to NA in attributes).
#
# Examples (display_dt is built by an OC summary article):
#
#   # Default: equal-weight linear of Detection, Sen, Spec
#   display_dt <- compute_oc_score(display_dt)
#
#   # Weighted linear emphasising detection
#   display_dt <- compute_oc_score(display_dt,
#                                  weights = c(0.5, 0.25, 0.25))
#
#   # Geometric mean (penalises imbalance)
#   display_dt <- compute_oc_score(display_dt,
#                                  method = "geometric",
#                                  score_col = "S_geom")
#
#   # Youden's J
#   display_dt <- compute_oc_score(display_dt,
#                                  method = "youden",
#                                  score_col = "J")
# ---------------------------------------------------------------------------

compute_oc_score <- function(data,
                             method    = c("linear", "geometric", "youden"),
                             metrics   = NULL,
                             weights   = NULL,
                             score_col = "S") {

  # ---- Argument validation -----------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'data' has zero rows; nothing to score.", call. = FALSE)
  }
  method <- match.arg(method)

  if (!is.character(score_col) || length(score_col) != 1L ||
      is.na(score_col) || nchar(score_col) == 0L) {
    stop("'score_col' must be a single non-empty character string.",
         call. = FALSE)
  }

  # Default metrics per method.  Track whether the user supplied
  # 'metrics' explicitly; the canonical "benefit-protection" default
  # weights (see below) only apply when BOTH metrics and weights
  # are left at their defaults.
  metrics_defaulted <- is.null(metrics)
  if (metrics_defaulted) {
    metrics <- switch(method,
                      linear    = c("Detection", "Sen", "Spec"),
                      geometric = c("Detection", "Sen", "Spec"),
                      youden    = c("Sen", "Spec"))
  }
  if (!is.character(metrics) || length(metrics) < 2L ||
      any(is.na(metrics)) || any(nchar(metrics) == 0L)) {
    stop("'metrics' must be a character vector of length >= 2 with no NA / empty entries.",
         call. = FALSE)
  }

  # Method-specific arity checks
  if (method == "youden" && length(metrics) != 2L) {
    stop(sprintf(
      "method = 'youden' requires exactly 2 metrics (got %d). ",
      length(metrics)),
      "Youden's J is the fixed formula x_1 + x_2 - 1; the default ",
      "is c(\"Sen\", \"Spec\").", call. = FALSE)
  }

  # Columns present?
  missing_cols <- setdiff(metrics, names(data))
  if (length(missing_cols)) {
    stop(sprintf(
      "Metric column(s) not found in 'data': %s.  Available columns: %s.",
      paste(shQuote(missing_cols), collapse = ", "),
      paste(shQuote(names(data)),  collapse = ", ")),
      call. = FALSE)
  }

  # All metric columns numeric?
  for (col in metrics) {
    if (!is.numeric(data[[col]])) {
      stop(sprintf("Metric column '%s' must be numeric; got %s.",
                   col, paste(class(data[[col]]), collapse = "/")),
           call. = FALSE)
    }
  }

  # Weights handling
  if (method == "youden") {
    if (!is.null(weights)) {
      message("'weights' is ignored for method = 'youden'.")
    }
    weights_used <- NA_real_
  } else {
    if (is.null(weights)) {
      # Canonical default: when metrics is also defaulted (i.e. the
      # user is using {Detection, Sen, Spec} under linear/geometric)
      # use the "benefit-protection" weighting (1.5, 1, 1.75).  This
      # gives Spec the largest weight, which in the harm-subgroup
      # frame minimises false-positive harm calls -- equivalently,
      # minimises wrongful exclusion of true beneficiaries from the
      # treatment recommendation Hc.  See file header for rationale.
      # Any other metric set falls back to equal weights.
      if (metrics_defaulted &&
          identical(metrics, c("Detection", "Sen", "Spec"))) {
        weights <- c(1.5, 1, 1.75)
      } else {
        weights <- rep(1 / length(metrics), length(metrics))
      }
    }
    if (!is.numeric(weights) || length(weights) != length(metrics) ||
        anyNA(weights)) {
      stop("'weights' must be a numeric vector with no NAs, same length as 'metrics'.",
           call. = FALSE)
    }
    if (any(weights < 0)) {
      stop("'weights' must be non-negative.", call. = FALSE)
    }
    if (sum(weights) <= 0) {
      stop("'weights' must sum to a positive value.", call. = FALSE)
    }
    weights_used <- weights / sum(weights)   # auto-normalise to sum to 1
  }

  # ---- Compute the score -------------------------------------------------
  # Build column-by-column via `[[` so the extraction works identically
  # on data.frame and data.table inputs (the `data[, metrics]` form
  # triggers data.table's NSE and fails when metrics is a variable).
  M <- do.call(cbind, lapply(metrics, function(col) data[[col]]))
  colnames(M) <- metrics

  if (method == "geometric" && any(M < 0, na.rm = TRUE)) {
    stop("method = 'geometric' requires all metric values to be non-negative.",
         call. = FALSE)
  }

  S <- switch(
    method,

    linear = as.numeric(M %*% weights_used),

    geometric = {
      # Per row: prod_k x_k^{w_k}.  Zero in any component -> 0 in product.
      # Use log-sum-exp safe path: when any value is exactly 0, result is 0.
      apply(M, 1, function(row) {
        if (anyNA(row))       return(NA_real_)
        if (any(row == 0))    return(0)
        exp(sum(weights_used * log(row)))
      })
    },

    youden = M[, 1L] + M[, 2L] - 1
  )

  # Attach provenance for transparency.
  attr(S, "score_meta") <- list(
    method  = method,
    metrics = metrics,
    weights = weights_used
  )

  result <- data
  result[[score_col]] <- S
  result
}
