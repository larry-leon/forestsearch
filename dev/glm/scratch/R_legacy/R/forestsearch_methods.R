# =============================================================================
# forestsearch_methods.R - S3 Methods for forestsearch Objects
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helpers (not exported)
# -----------------------------------------------------------------------------

#' Safely extract a value from nested list slots
#'
#' Tries multiple accessor paths and returns the first non-NULL value.
#' Eliminates fragile slot-specific code scattered across methods.
#'
#' @param x A forestsearch object.
#' @param ... Character vectors, each a path of nested names to try.
#' @return The first non-NULL value found, or NULL.
#' @keywords internal
#' @noRd
.fs_get <- function(x, ...) {
  paths <- list(...)
  for (path in paths) {
    val <- x
    ok <- TRUE
    for (key in path) {
      if (is.null(val) || !is.list(val) || !key %in% names(val)) {
        ok <- FALSE
        break
      }
      val <- val[[key]]
    }
    if (ok && !is.null(val)) return(val)
  }
  NULL
}


#' Extract human-readable subgroup labels
#'
#' Prioritises \code{grp.consistency$out_sg$sg.harm_label} (readable labels
#' like \code{"er <= 0"}) over bare factor names in \code{sg.harm}.
#'
#' @param x A forestsearch object.
#' @return Character vector of labels, or NULL.
#' @keywords internal
#' @noRd
.fs_sg_labels <- function(x) {
  labels <- .fs_get(
    x,
    c("grp.consistency", "out_sg", "sg.harm_label")
  )
  if (!is.null(labels)) {
    labels <- labels[!is.na(labels) & labels != ""]
    if (length(labels) > 0) return(labels)
  }
  x$sg.harm
}


# =============================================================================
# print.forestsearch
# =============================================================================

#' Print Method for forestsearch Objects
#'
#' Displays a concise summary of ForestSearch results including the
#' identified subgroup definition, consistency metrics, algorithm details,
#' and computation time.
#'
#' @param x A \code{forestsearch} object returned by
#'   \code{\link{forestsearch}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \dontrun{
#' fs <- forestsearch(df.analysis = mydata, ...)
#' print(fs)
#' # or simply:
#' fs
#' }
#'
#' @seealso \code{\link{summary.forestsearch}} for detailed output,
#'   \code{\link{plot.forestsearch}} for visualization.
#' @export
print.forestsearch <- function(x, ...) {
  cat("ForestSearch Results\n")
  cat("====================\n\n")

  if (is.null(x$sg.harm)) {
    cat("No subgroup identified.\n")
    return(invisible(x))
  }

  # --- Subgroup definition ---
  labels <- .fs_sg_labels(x)
  cat("Selected Subgroup:\n")
  cat("  Definition:", paste(labels, collapse = " & "), "\n")

  # --- sg_focus (try multiple locations) ---
  sg_focus <- .fs_get(
    x,
    c("grp.consistency", "sg_focus"),
    c("args_call_all", "sg_focus"),
    c("sg_focus")
  )
  if (!is.null(sg_focus)) {
    cat("  sg_focus:", sg_focus, "\n")
  }

  # --- Top-ranked consistency result ---
  top_result <- .fs_get(x, c("grp.consistency", "out_sg", "result"))

  if (!is.null(top_result) && nrow(top_result) > 0) {
    top <- top_result[1, ]
    if ("N" %in% names(top))
      cat("  N:", top$N, "\n")
    if ("hr" %in% names(top))
      cat("  HR:", round(as.numeric(top$hr), 3), "\n")
    if ("Pcons" %in% names(top))
      cat("  Pcons:", round(as.numeric(top$Pcons), 3), "\n")
  }

  # --- Algorithm info ---
  algorithm <- .fs_get(x, c("grp.consistency", "algorithm"))
  if (!is.null(algorithm)) {
    cat("  Algorithm:", algorithm, "\n")
  }

  n_eval <- .fs_get(x, c("grp.consistency", "n_candidates_evaluated"))
  n_pass <- .fs_get(x, c("grp.consistency", "n_passed"))
  if (!is.null(n_eval)) {
    cat("  Candidates evaluated:", n_eval, "\n")
    cat("  Candidates passed:", n_pass, "\n")
  }

  # --- Timing ---
  if (!is.null(x$minutes_all)) {
    cat("\nComputation time:", round(x$minutes_all, 2), "minutes\n")
  }

  invisible(x)
}


# =============================================================================
# summary.forestsearch
# =============================================================================

#' Summary Method for forestsearch Objects
#'
#' Provides a detailed summary of a ForestSearch analysis including input
#' parameters, variable selection results, consistency evaluation, and
#' the selected subgroup with key metrics.
#'
#' @param object A \code{forestsearch} object returned by
#'   \code{\link{forestsearch}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{object}.
#'
#' @examples
#' \dontrun{
#' fs <- forestsearch(df.analysis = mydata, ...)
#' summary(fs)
#' }
#'
#' @export
summary.forestsearch <- function(object, ...) {
  cat("ForestSearch Summary\n")
  cat("====================\n\n")

  # -------------------------------------------------------------------------
  # Analysis Parameters
  # -------------------------------------------------------------------------
  params <- object$args_call_all
  if (!is.null(params)) {
    cat("Analysis Parameters:\n")
    .print_param <- function(label, value) {
      if (!is.null(value)) cat("  ", label, ": ", value, "\n", sep = "")
    }
    .print_param("sg_focus",                params$sg_focus)
    .print_param("hr.threshold",            params$hr.threshold)
    .print_param("hr.consistency",          params$hr.consistency)
    .print_param("pconsistency.threshold",  params$pconsistency.threshold)
    .print_param("n.min",                   params$n.min)
    .print_param("fs.splits",               params$fs.splits)
    .print_param("maxk",                    params$maxk)
    .print_param("use_twostage",            params$use_twostage)
    .print_param("use_lasso",               params$use_lasso)
    .print_param("use_grf",                 params$use_grf)
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Variable Selection
  # -------------------------------------------------------------------------
  n_candidate <- length(object$confounders.candidate)
  n_evaluated <- length(object$confounders.evaluated)
  if (n_candidate > 0 || n_evaluated > 0) {
    cat("Variable Selection:\n")
    cat("  Candidate confounders:", n_candidate, "\n")
    cat("  Confounders evaluated:", n_evaluated, "\n")

    # GRF screening info
    if (!is.null(object$grf_res) && !inherits(object$grf_res, "try-error")) {
      n_grf_cuts <- length(object$grf_cuts)
      if (n_grf_cuts > 0) {
        cat("  GRF cuts applied:", n_grf_cuts, "\n")
      }
    }
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Search Space
  # -------------------------------------------------------------------------
  if (!is.null(object$prop_maxk) || !is.null(object$max_sg_est)) {
    cat("Search Space:\n")
    if (!is.null(object$prop_maxk)) {
      cat("  Proportion of max combinations searched:",
          round(object$prop_maxk, 4), "\n")
    }
    if (!is.null(object$max_sg_est)) {
      cat("  Maximum subgroup HR estimate:",
          round(object$max_sg_est, 3), "\n")
    }
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Consistency Evaluation
  # -------------------------------------------------------------------------
  gc <- object$grp.consistency
  if (!is.null(gc)) {
    cat("Consistency Evaluation:\n")
    algorithm <- gc$algorithm %||% "fixed"
    cat("  Algorithm:", algorithm, "\n")

    if (!is.null(gc$n_candidates_evaluated)) {
      cat("  Candidates evaluated:", gc$n_candidates_evaluated, "\n")
      cat("  Candidates passed:", gc$n_passed, "\n")
    }

    # Two-stage specific info
    if (identical(algorithm, "twostage")) {
      ts <- params$twostage_args
      if (!is.null(ts)) {
        cat("  Stage 1 screening splits:", ts$n.splits.screen, "\n")
        cat("  Stage 2 batch size:", ts$batch.size, "\n")
      }
    }
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Selected Subgroup
  # -------------------------------------------------------------------------
  if (!is.null(object$sg.harm)) {
    labels <- .fs_sg_labels(object)
    cat("Selected Subgroup:\n")
    cat("  Definition:", paste(labels, collapse = " & "), "\n")

    # Factor-level names (technical)
    if (!identical(labels, object$sg.harm)) {
      cat("  Factor names:", paste(object$sg.harm, collapse = " & "), "\n")
    }

    top_result <- .fs_get(object, c("grp.consistency", "out_sg", "result"))
    if (!is.null(top_result) && nrow(top_result) > 0) {
      top <- top_result[1, ]
      if ("N" %in% names(top))
        cat("  Sample size:", top$N, "\n")
      if ("hr" %in% names(top))
        cat("  Hazard ratio:", round(as.numeric(top$hr), 3), "\n")
      if ("Pcons" %in% names(top))
        cat("  Consistency:",
            round(as.numeric(top$Pcons) * 100, 1), "%\n")
    }

    # Dataset sizes
    if (!is.null(object$df.est)) {
      n_est <- nrow(object$df.est)
      n_harm <- sum(object$df.est$treat.recommend == 0, na.rm = TRUE)
      cat("  Estimation data: n =", n_est,
          "(harm =", n_harm,
          ", complement =", n_est - n_harm, ")\n")
    }
  } else {
    cat("No subgroup identified.\n")
  }

  # -------------------------------------------------------------------------
  # Timing
  # -------------------------------------------------------------------------
  if (!is.null(object$minutes_all)) {
    cat("\nComputation time:", round(object$minutes_all, 2), "minutes\n")
  }

  invisible(object)
}


# =============================================================================
# plot.forestsearch
# =============================================================================

#' Plot ForestSearch Results
#'
#' Dispatches to \code{\link{plot_sg_results}} for Kaplan-Meier curves,
#' hazard-ratio forest plots, or combined panels.
#'
#' @param x A \code{forestsearch} object returned by
#'   \code{\link{forestsearch}}.
#' @param type Character. Type of plot:
#'   \describe{
#'     \item{\code{"combined"}}{KM curves + forest plot (default)}
#'     \item{\code{"km"}}{Kaplan-Meier survival curves only}
#'     \item{\code{"forest"}}{Hazard-ratio forest plot only}
#'     \item{\code{"summary"}}{Summary statistics panel}
#'   }
#' @param outcome.name Character. Name of time-to-event column.
#'   Default: \code{"Y"}.
#' @param event.name Character. Name of event indicator column.
#'   Default: \code{"Event"}.
#' @param treat.name Character. Name of treatment column.
#'   Default: \code{"Treat"}.
#' @param ... Additional arguments passed to \code{\link{plot_sg_results}},
#'   such as \code{by.risk}, \code{conf.level}, \code{est.scale},
#'   \code{sg0_name}, \code{sg1_name}, \code{treat_labels}, \code{colors},
#'   \code{title}, \code{show_events}, \code{show_ci}, \code{show_logrank},
#'   \code{show_hr}.
#'
#' @return Invisibly returns the plot result from
#'   \code{\link{plot_sg_results}}.
#'
#' @seealso \code{\link{plot_sg_results}} for full control over appearance,
#'   \code{\link{plot_sg_weighted_km}} for weighted KM curves,
#'   \code{\link{plot_subgroup_results_forestplot}} for publication-ready
#'   forest plots.
#'
#' @examples
#' \dontrun{
#' fs <- forestsearch(df.analysis = mydata, ...)
#'
#' # Combined KM + forest plot (default)
#' plot(fs)
#'
#' # KM curves only
#' plot(fs, type = "km")
#'
#' # Forest plot only
#' plot(fs, type = "forest")
#'
#' # With non-standard column names
#' plot(fs, type = "km",
#'      outcome.name = "os_months",
#'      event.name = "os_event",
#'      treat.name = "treatment")
#'
#' # With custom labels
#' plot(fs, sg0_name = "High Risk", sg1_name = "Standard Risk",
#'      treat_labels = c("0" = "Placebo", "1" = "Active Drug"))
#' }
#'
#' @export
plot.forestsearch <- function(x,
                              type = c("combined", "km",
                                       "forest", "summary"),
                              outcome.name = "Y",
                              event.name = "Event",
                              treat.name = "Treat",
                              ...) {

  type <- match.arg(type)

  if (is.null(x$df.est)) {
    message("No subgroup identified -- nothing to plot.")
    return(invisible(x))
  }

  result <- plot_sg_results(
    fs.est       = x,
    plot_type    = type,
    outcome.name = outcome.name,
    event.name   = event.name,
    treat.name   = treat.name,
    ...
  )

  invisible(result)
}
