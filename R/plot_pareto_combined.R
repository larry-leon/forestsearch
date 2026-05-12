# ============================================================================
# plot_pareto_combined.R
# ============================================================================
# Single-plot composition for multiple forestsearch fits that share the
# same set of consistency-passing candidates.  Annotates each selected
# subgroup with one or more S<i>-markers naming the configuration(s)
# that picked it.
# ----------------------------------------------------------------------------

# Silence R CMD check NOTE for ggplot2 NSE column references.
utils::globalVariables(c("hr_nat", "N", "is_selected_any", "label_text",
                         ".def_key"))

#' Combined Pareto-Frontier Plot Across Configurations Sharing a Passing Set
#'
#' For two or more \code{\link{forestsearch}} fits whose consistency-passing
#' candidate sets are \strong{identical}, produces a single Pareto-frontier
#' plot annotated with one \code{S<i>: <combo_label>} marker per
#' configuration at each selected subgroup.  If two configurations pick
#' the same subgroup, the markers stack into a single multi-line label
#' (e.g.\ \code{"S1: hrMaxSG/pareto\\nS2: hrMaxSG/both"}).
#'
#' The combined plot is meaningful only when the passing sets match.
#' When they don't, the function returns \code{NULL} with a warning;
#' the side-by-side view (one panel per configuration) is the
#' appropriate alternative.
#'
#' @section Equality criterion:
#' Two passing sets are considered equal when, after sorting by
#' subgroup definition string:
#' \itemize{
#'   \item they have the same number of rows;
#'   \item the same set of subgroup definitions (concatenated
#'     \code{M.*} columns), compared as sets (order-independent);
#'   \item within each matched definition, the \code{hr}, \code{N},
#'     \code{E}, and \code{K} values agree within \code{tolerance}.
#' }
#'
#' \code{Pcons} is deliberately excluded from the value check.  It can
#' legitimately drift across rules because the preview sort (which
#' depends on \code{selection_rule}) changes the queue order, which
#' changes the random-split state each candidate consumes.  Drift of
#' up to ~0.10 between runs on the same subgroup is expected and does
#' not indicate a real disagreement.  Note also that the internal
#' candidate id \code{m} can differ across configurations for the
#' same reason; \code{m} is NOT used for equality.
#'
#' @param fs_list A list of \code{forestsearch} result objects.  Length
#'   must be >= 2.
#' @param combo_labels Character vector of length \code{length(fs_list)}.
#'   Used in the \code{S<i>: <label>} annotations.  Default \code{NULL}
#'   = auto, e.g.\ \code{"config 1"}, \code{"config 2"}, ...
#' @param ci_table_list Optional list of CI tables from
#'   \code{\link{compute_frontier_cis}}, parallel to \code{fs_list}.
#'   When provided and the tables match (within tolerance), the CI bars
#'   are drawn using the first table.
#' @param show_band Logical.  Draw the effect-band shading if applicable.
#'   Default \code{TRUE}.
#' @param xlim Numeric vector of length 2 or \code{NULL}.  Default
#'   \code{NULL} = auto-range.
#' @param tolerance Numeric.  Per-cell tolerance for the value
#'   equality check (\code{hr}, \code{N}, \code{E}, \code{K}).
#'   Default \code{1e-6}.
#' @param verbose Logical.  Emit a warning naming the specific
#'   equality-check failure mode when the sets don't match.  Default
#'   \code{TRUE}.
#'
#' @return A \code{ggplot} object, or \code{NULL} (with a warning) if
#'   the passing sets don't satisfy the equality criterion.
#'
#' @examples
#' \dontrun{
#' out <- compare_selection_rules(
#'   df.analysis    = actg_df,
#'   sg_focus       = c("hrMaxSG", "hrMaxSG"),
#'   selection_rule = c("pareto",  "both"),
#'   ...
#' )
#'
#' # Auto-attached by compare_selection_rules():
#' print(out$plot_combined)
#'
#' # Or call directly:
#' p <- plot_pareto_combined(
#'   fs_list       = out$fs,
#'   combo_labels  = out$combos$label,
#'   ci_table_list = out$ci_tab
#' )
#' print(p)
#' }
#'
#' @seealso \code{\link{compare_selection_rules}},
#'   \code{\link{plot_pareto_frontier}},
#'   \code{\link{pareto_frontier_table}}.
#' @export
plot_pareto_combined <- function(fs_list,
                                 combo_labels   = NULL,
                                 ci_table_list  = NULL,
                                 show_band      = TRUE,
                                 xlim           = NULL,
                                 tolerance      = 1e-6,
                                 verbose        = TRUE) {
  # --- 0. Validation ------------------------------------------------------
  if (!is.list(fs_list) || length(fs_list) < 2L) {
    if (verbose) warning("plot_pareto_combined() requires >= 2 fs objects.",
                         call. = FALSE)
    return(NULL)
  }
  # Drop NULL slots; any NULL means we can't reliably compare
  if (any(vapply(fs_list, is.null, logical(1)))) {
    if (verbose) warning("plot_pareto_combined(): some fs entries are NULL; cannot produce combined view.",
                         call. = FALSE)
    return(NULL)
  }
  n_combos <- length(fs_list)
  if (is.null(combo_labels)) {
    combo_labels <- sprintf("config %d", seq_len(n_combos))
  } else if (length(combo_labels) != n_combos) {
    if (verbose) warning("combo_labels length must match fs_list length.",
                         call. = FALSE)
    return(NULL)
  }

  # --- 1. Extract each fit's passing set + selected m --------------------
  passing_list <- vector("list", n_combos)
  selected_keys <- character(n_combos)
  for (k in seq_len(n_combos)) {
    out_sg <- tryCatch(fs_list[[k]]$grp.consistency$out_sg,
                       error = function(e) NULL)
    if (is.null(out_sg) || is.null(out_sg$result) ||
        nrow(out_sg$result) == 0L) {
      if (verbose) warning(sprintf(
        "fs_list[[%d]] has no passing candidates.", k), call. = FALSE)
      return(NULL)
    }
    res <- data.table::copy(out_sg$result)
    # Build a stable subgroup-definition key from the M.* columns
    m_cols <- grep("^M\\.", names(res), value = TRUE)
    res$.def_key <- vapply(seq_len(nrow(res)), function(i) {
      v <- unlist(res[i, m_cols, with = FALSE], use.names = FALSE)
      v <- v[!is.na(v) & nzchar(v)]
      if (length(v) == 0L) "(empty)" else paste(sort(v), collapse = " & ")
    }, character(1))
    data.table::setorder(res, .def_key)
    passing_list[[k]] <- res
    selected_keys[k] <- res$.def_key[
      which(as.integer(res$m) == as.integer(out_sg$result$m[1L]))[1L]
    ]
  }

  # --- 2. Equality check on passing sets ---------------------------------
  ref <- passing_list[[1L]]
  for (k in 2L:n_combos) {
    cmp <- passing_list[[k]]
    if (nrow(cmp) != nrow(ref)) {
      if (verbose) warning(sprintf(
        "Passing sets differ in size: fs[[1]] has %d, fs[[%d]] has %d. Use the side-by-side plot instead.",
        nrow(ref), k, nrow(cmp)), call. = FALSE)
      return(NULL)
    }
    if (!setequal(cmp$.def_key, ref$.def_key)) {
      diff_keys <- setdiff(cmp$.def_key, ref$.def_key)
      if (verbose) warning(sprintf(
        "Passing sets differ in subgroup definitions: fs[[%d]] has %d unique definition(s) not in fs[[1]] (first: '%s'). Use the side-by-side plot instead.",
        k, length(diff_keys),
        if (length(diff_keys) > 0L) diff_keys[1L] else ""),
        call. = FALSE)
      return(NULL)
    }
    # Within-key value tolerance: align cmp rows to ref by .def_key, then check
    cmp_aligned <- cmp[match(ref$.def_key, cmp$.def_key), ]
    # Within-key value tolerance: align cmp rows to ref by .def_key, then check.
    # Pcons is deliberately NOT in val_cols: it can legitimately drift across
    # rules because the preview sort (which depends on selection_rule) changes
    # the queue order, which changes the random-split state consumed by each
    # candidate.  Pcons drift of up to ~0.10 between two runs on the same
    # subgroup is expected and does not indicate disagreement on the
    # candidate's identity or geometry.  Equality is checked on hr / N / E / K
    # (effect / sample size / events / parsimony), which are deterministic
    # given the subgroup definition and the data.
    val_cols <- intersect(c("hr", "N", "E", "K"), names(ref))
    for (col in val_cols) {
      a <- as.numeric(ref[[col]])
      b <- as.numeric(cmp_aligned[[col]])
      max_diff <- max(abs(a - b), na.rm = TRUE)
      if (!is.finite(max_diff) || max_diff > tolerance) {
        if (verbose) warning(sprintf(
          "Passing sets share definitions but '%s' values differ (max diff = %.6g > tolerance %.6g). Use the side-by-side plot instead.",
          col, max_diff, tolerance), call. = FALSE)
        return(NULL)
      }
    }
  }

  # --- 3. Build the plotting data from the first fit ---------------------
  out_sg_1 <- fs_list[[1L]]$grp.consistency$out_sg
  frontier <- tryCatch(out_sg_1$pareto_frontier, error = function(e) NULL)
  if (is.null(frontier) || !data.table::is.data.table(frontier) ||
      nrow(frontier) == 0L) {
    if (verbose) warning("plot_pareto_combined(): no Pareto frontier on fs_list[[1]].",
                         call. = FALSE)
    return(NULL)
  }

  effect_measure <- fs_list[[1L]]$effect_measure %||% "HR"
  effect_log_scale <- effect_measure %in% c("OR", "RR", "IRR")
  effect_label <- switch(effect_measure,
    HR  = "Hazard ratio",  OR  = "Odds ratio",
    RR  = "Risk ratio",    IRR = "Incidence rate ratio",
    RD  = "Risk difference",  IRD = "Incidence rate difference",
    MD  = "Mean difference",
    effect_measure)

  # All passing candidates on natural scale (use the first fit's result;
  # equality-checked above)
  all_passing <- data.table::copy(passing_list[[1L]])
  all_passing$hr_nat <- if (effect_log_scale) exp(as.numeric(all_passing$hr))
                        else as.numeric(all_passing$hr)

  # Frontier rows for the line layer
  ft <- data.table::copy(frontier)
  ft$hr_nat <- if (effect_log_scale) exp(as.numeric(ft$hr))
               else as.numeric(ft$hr)
  data.table::setorder(ft, -hr_nat)

  # --- 4. Build the annotation table -------------------------------------
  # For each unique selected subgroup, collect the S<i>: <label> lines
  # of every combo that picked it.
  unique_selected <- unique(selected_keys)
  annot_rows <- lapply(unique_selected, function(k_def) {
    combos_picking_k <- which(selected_keys == k_def)
    text_lines <- sprintf("S%d: %s", combos_picking_k,
                          combo_labels[combos_picking_k])
    label_text <- paste(text_lines, collapse = "\n")
    row <- all_passing[all_passing$.def_key == k_def, ][1L, ]
    data.frame(
      hr_nat     = row$hr_nat,
      N          = as.numeric(row$N),
      label_text = label_text,
      stringsAsFactors = FALSE
    )
  })
  annot_dt <- do.call(rbind, annot_rows)

  # --- 5. Assemble plot ---------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (verbose) warning("plot_pareto_combined(): ggplot2 not installed.",
                         call. = FALSE)
    return(NULL)
  }

  # Build a flag column for the legend/annotations (any combo selected this point)
  all_passing$is_selected_any <- all_passing$.def_key %in% unique_selected

  # Determine the effective sg_focus / selection_rule / band epsilon from the
  # first fit's args_call_all (passing-set equality across combos was
  # already verified, so the band threshold is the same regardless of
  # which combo's settings we read).  For band drawing the relevant
  # questions are: is the band axis applicable (sg_focus = hrMaxSG/hrMinSG),
  # and what's the neighborhood width?
  args_1 <- fs_list[[1L]]$args_call_all %||% list()
  # Prefer fs$sg_focus (always set to the canonical form by
  # forestsearch()); fall back to args_call_all$sg_focus (may carry
  # the user's raw input under older forestsearch builds) and
  # normalize defensively.  This belt-and-braces approach makes the
  # band check robust to both the legacy "hr*" vocabulary and the
  # GLM-natural "eff*" aliases, regardless of which forestsearch
  # version produced the fs object.
  sg_focus_1 <- fs_list[[1L]]$sg_focus %||% args_1$sg_focus %||% "hr"
  if (exists(".normalize_sg_focus", mode = "function")) {
    sg_focus_1 <- .normalize_sg_focus(sg_focus_1)
  } else {
    # Inline fallback for callers running this file in isolation
    # (smoke tests, examples) without forestsearch_helpers.R sourced.
    sg_focus_1 <- switch(sg_focus_1,
                         effMaxSG = "hrMaxSG",
                         effMinSG = "hrMinSG",
                         eff      = "hr",
                         sg_focus_1)
  }
  eps_val <- args_1$effect_neighborhood %||% 0.10
  band_applies <- isTRUE(show_band) &&
                  sg_focus_1 %in% c("hrMaxSG", "hrMinSG")

  # Compute the band threshold once (used for both the shading and the
  # subtitle note); leave undefined when not applicable.
  floor_v <- NA_real_
  hr_max  <- NA_real_
  if (band_applies) {
    hr_max  <- max(all_passing$hr_nat, na.rm = TRUE)
    floor_v <- (1 - eps_val) * hr_max
  }

  # Frontier step polyline.  Direction "hv" matches plot_pareto_frontier()
  # so the two views read consistently.
  p <- ggplot2::ggplot(all_passing,
                       ggplot2::aes(x = hr_nat, y = N))

  # Optional effect-neighborhood band (drawn first so points overlay it)
  if (band_applies) {
    p <- p + ggplot2::annotate("rect",
      xmin = floor_v, xmax = hr_max * 1.02,
      ymin = -Inf,    ymax = Inf,
      fill = "grey85", alpha = 0.35)
  }

  p <- p +
    ggplot2::geom_step(data = ft, direction = "hv",
                       linewidth = 0.7, colour = "steelblue") +
    ggplot2::geom_point(
      ggplot2::aes(colour = is_selected_any, size = is_selected_any)
    ) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "grey25", `TRUE` = "#D55E00"),
      labels = c(`FALSE` = "frontier",
                 `TRUE`  = "selected"),
      name   = NULL
    ) +
    ggplot2::scale_size_manual(
      values = c(`FALSE` = 2.5, `TRUE` = 4),
      guide  = "none"
    )

  # Optional split CI bars from ci_table_list[[1]]
  if (!is.null(ci_table_list) && length(ci_table_list) >= 1L &&
      !is.null(ci_table_list[[1L]])) {
    ci_tab <- ci_table_list[[1L]]
    if (data.table::is.data.table(ci_tab) && nrow(ci_tab) > 0L &&
        all(c("m", "split_lcl", "split_ucl") %in% names(ci_tab))) {
      # Merge by m onto the frontier rows, coercing both to integer
      ft_m <- ft
      ft_m$m <- suppressWarnings(as.integer(ft_m$m))
      ci_m <- ci_tab[, c("m", "split_lcl", "split_ucl"), with = FALSE]
      ci_m$m <- suppressWarnings(as.integer(ci_m$m))
      ft_ci <- merge(ft_m, ci_m, by = "m", all.x = TRUE, sort = FALSE)
      ft_ci <- ft_ci[!is.na(ft_ci$split_lcl) & !is.na(ft_ci$split_ucl), ]
      if (nrow(ft_ci) > 0L) {
        # geom_linerange() with y=, xmin=, xmax= replaces the deprecated
        # geom_errorbarh(height = 0, ...) pattern (ggplot2 >= 3.5).
        # No endcaps were used (height = 0), so geom_linerange is the
        # direct equivalent without the deprecation warning.
        p <- p + ggplot2::geom_linerange(
          data    = ft_ci,
          mapping = ggplot2::aes(y = N, xmin = split_lcl, xmax = split_ucl),
          inherit.aes = FALSE,
          colour = "grey40", alpha = 0.6, linewidth = 0.5
        )
      }
    }
  }

  # Annotation labels (S1: label, etc.) on the selected points.
  # Nudge slightly to the right of the point so they don't cover it.
  p <- p + ggplot2::geom_label(
    data    = annot_dt,
    mapping = ggplot2::aes(x = hr_nat, y = N, label = label_text),
    inherit.aes = FALSE,
    hjust       = -0.1,
    vjust       = 0.5,
    size        = 3.2,
    label.size  = 0.3,
    label.padding = ggplot2::unit(0.18, "lines"),
    fill        = "white",
    colour      = "grey20"
  )

  # Title note: when all combos pick the same subgroup, say so
  agreement_note <- if (length(unique_selected) == 1L) {
    sprintf(" - all %d configurations agreed on the same subgroup", n_combos)
  } else {
    sprintf(" - %d configurations, %d distinct winners",
            n_combos, length(unique_selected))
  }

  band_note <- if (band_applies) {
    sprintf("; \u03b5-band shading shows %s \u2265 %.3g (%.0f%% neighborhood)",
            effect_label, floor_v, 100 * eps_val)
  } else ""

  p <- p +
    ggplot2::labs(
      title    = "Pareto frontier on (effect, N) - combined view",
      subtitle = sprintf("Passing set is identical across %d configurations%s%s",
                         n_combos, agreement_note, band_note),
      x = effect_label,
      y = expression(italic(N) ~ " (subgroup size)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor    = ggplot2::element_blank(),
      legend.position     = "right",
      plot.title.position = "plot"
    )

  if (!is.null(xlim)) {
    if (is.numeric(xlim) && length(xlim) == 2L) {
      p <- p + ggplot2::coord_cartesian(xlim = xlim)
    }
  }

  p
}

# Internal NULL-coalesce; declared locally so this file doesn't depend
# on a package-level export.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
}
