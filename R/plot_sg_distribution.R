#' Plot Distribution of Identified Subgroups
#'
#' Bar chart of subgroups identified across simulations, filtered to
#' those appearing in at least \code{min_pct} of the found simulations.
#' Supports two column schemas (MRCT and \code{run_simulation_analysis}),
#' top-\eqn{K} capping with a pooled "Other" bar, automatic threshold
#' halving when no label clears the initial cut, and an \code{n_bars}
#' attribute for adaptive figure-height computation.
#'
#' @param results \code{data.frame} or \code{data.table} with at least
#'   two columns: a per-replicate detection flag (0/1) and an
#'   identified-subgroup label string.  Column names are configurable
#'   via \code{any_col} and \code{label_col}.  The defaults
#'   (\code{"any_found"}, \code{"sg_found"}) match
#'   \code{\link{mrct_region_sims}()}; pass
#'   \code{any_col = "any.H"}, \code{label_col = "sg.def"} to plot
#'   \code{\link{run_simulation_analysis}()} output.
#' @param min_pct Numeric.  Minimum percentage threshold for display
#'   (0-100).  Default: \code{5}.
#' @param title Character.  Plot title.
#' @param wrap_width Integer.  Character width for wrapping long
#'   subgroup labels.  Default: \code{25}.
#' @param any_col Character.  Column name of the per-replicate
#'   detection flag.  Default: \code{"any_found"}.
#' @param label_col Character.  Column name of the identified-subgroup
#'   label string.  Default: \code{"sg_found"}.
#' @param top_k Integer or \code{NULL}.  Maximum number of bars to
#'   display; additional passing labels are pooled into a single
#'   "Other (\eqn{k} labels)" bar (see \code{show_other}).  Prevents
#'   unreadable charts under diffuse null hypotheses where finds are
#'   fragmented across many singleton labels.  Default: \code{NULL}
#'   (no cap, i.e. the original behaviour).
#' @param show_other Logical.  When \code{top_k} is set and extra
#'   labels are pooled, attach the "Other" bar.  Default: \code{TRUE}.
#'   Set to \code{FALSE} to suppress the pooled bar entirely; its
#'   count and percent are still reported in the subtitle for honesty.
#' @param min_floor_pct Numeric.  Lower bound on automatic threshold
#'   halving.  Default: \code{0.5}.
#' @param max_halvings Integer.  Maximum number of halvings applied to
#'   \code{min_pct} when no label clears the initial threshold.
#'   Default: \code{0} (no halving, i.e. the original behaviour).
#'   Set to \code{2} to match the extended-diagnostic behaviour used
#'   in the benefit-search vignettes.
#' @param bar_label_inside Logical.  Render bar labels inside bars
#'   (\code{hjust = 1.05}) instead of outside (\code{hjust = -0.1}).
#'   Prevents small-count labels from clipping off the right edge.
#'   Default: \code{FALSE} (original outside-of-bar placement).
#' @param placeholder_on_empty Logical.  When no subgroups clear the
#'   threshold (even after halving), return a minimal
#'   "no data" placeholder plot with \code{attr(p, "n_bars") = 0L}
#'   instead of \code{invisible(NULL)}.  Default: \code{FALSE} (the
#'   original behaviour).  Set to \code{TRUE} to simplify adaptive-
#'   figure-height rendering in vignettes.
#'
#' @return A \code{ggplot2} object, or \code{invisible(NULL)} when no
#'   subgroups are found and \code{placeholder_on_empty = FALSE}.
#'   When a plot is returned it carries two attributes:
#'   \describe{
#'     \item{\code{n_bars}}{Integer.  Number of bars in the plot
#'       (0 in the placeholder case).  Use
#'       \code{\link{sgdist_fig_height}()} to compute an adaptive
#'       figure height for knitr.}
#'     \item{\code{effective_min_pct}}{Numeric.  The threshold in
#'       effect after any automatic halving.}
#'   }
#'
#' @seealso \code{\link{sgdist_fig_height}} for computing an adaptive
#'   figure height from \code{n_bars}, and
#'   \code{\link{run_simulation_analysis}} /
#'   \code{\link{mrct_region_sims}} for the upstream simulation
#'   engines.
#'
#' @examples
#' \dontrun{
#' # MRCT / mrct_region_sims() output (default column schema)
#' fs <- forestsearch(gbsg,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
#' plot_sg_distribution(fs$grp.consistency$result_new)
#'
#' # run_simulation_analysis() output with top-K cap and auto-halving
#' fs_rows <- results_alt[results_alt$analysis == "FS", ]
#' p <- plot_sg_distribution(
#'   fs_rows,
#'   any_col              = "any.H",
#'   label_col            = "sg.def",
#'   top_k                = 15L,
#'   show_other           = FALSE,
#'   max_halvings         = 2L,
#'   bar_label_inside     = TRUE,
#'   placeholder_on_empty = TRUE,
#'   title                = "HTE: Subgroups Identified (H-hat)"
#' )
#' fig_h <- sgdist_fig_height(attr(p, "n_bars"))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text coord_flip scale_y_continuous expansion labs theme_minimal theme element_text element_blank annotate xlim ylim theme_void
#' @export
plot_sg_distribution <- function(results,
                                 min_pct              = 5,
                                 title                = "Distribution of Identified Subgroups",
                                 wrap_width           = 25,
                                 any_col              = "any_found",
                                 label_col            = "sg_found",
                                 top_k                = NULL,
                                 show_other           = TRUE,
                                 min_floor_pct        = 0.5,
                                 max_halvings         = 0L,
                                 bar_label_inside     = FALSE,
                                 placeholder_on_empty = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Please install it.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Schema validation
  # ---------------------------------------------------------------------------
  if (!(any_col   %in% names(results))) {
    stop(sprintf("Column '%s' not found in 'results'.", any_col), call. = FALSE)
  }
  if (!(label_col %in% names(results))) {
    stop(sprintf("Column '%s' not found in 'results'.", label_col), call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Restrict to detecting replicates; normalise NA / "" labels
  # ---------------------------------------------------------------------------
  any_flag  <- as.integer(results[[any_col]])
  labels    <- results[[label_col]]
  labels[is.na(labels) | !nzchar(labels)] <- "(none identified)"

  n_total  <- length(any_flag)
  found    <- which(any_flag == 1L)
  n_found  <- length(found)

  if (n_found == 0L) {
    if (placeholder_on_empty) {
      return(.empty_sgdist_plot(
        title = title,
        msg   = sprintf("No subgroups identified\n(0 of %d sims)", n_total)
      ))
    }
    message("No subgroups identified in any simulation.")
    return(invisible(NULL))
  }

  tab <- sort(table(labels[found]), decreasing = TRUE)
  counts <- data.frame(
    sg_found         = names(tab),
    count            = as.integer(tab),
    stringsAsFactors = FALSE
  )
  counts$pct <- 100 * counts$count / n_found

  # ---------------------------------------------------------------------------
  # Auto-halving: lower the threshold (within floor) until at least one
  # label clears it, up to `max_halvings` attempts.
  # ---------------------------------------------------------------------------
  current    <- min_pct
  attempts   <- c(current)
  n_halvings <- 0L
  keep       <- counts[counts$pct >= current, , drop = FALSE]

  while (nrow(keep) == 0L &&
         n_halvings < as.integer(max_halvings) &&
         (current / 2) >= min_floor_pct) {
    current    <- current / 2
    attempts   <- c(attempts, current)
    n_halvings <- n_halvings + 1L
    keep       <- counts[counts$pct >= current, , drop = FALSE]
  }

  if (nrow(keep) == 0L) {
    msg <- if (n_halvings > 0L) {
      sprintf(
        "No label \u2265 %g%% even after halving to %g%%\n(%d of %d sims found)",
        attempts[1], current, n_found, n_total
      )
    } else {
      sprintf("No subgroup appears in \u2265 %g%% of found simulations.",
              min_pct)
    }
    if (placeholder_on_empty) {
      return(.empty_sgdist_plot(title = title, msg = msg))
    }
    message(msg)
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Top-K cap: always compute pool metadata (for the subtitle) but only
  # attach the "Other" bar when show_other = TRUE.
  # ---------------------------------------------------------------------------
  n_distinct_passing <- nrow(keep)
  n_distinct_total   <- nrow(counts)
  pool_info          <- NULL

  if (!is.null(top_k) && n_distinct_passing > as.integer(top_k)) {
    top  <- keep[seq_len(top_k), , drop = FALSE]
    rest <- keep[seq.int(top_k + 1L, n_distinct_passing), , drop = FALSE]
    pool_info <- list(
      n_labels = nrow(rest),
      count    = sum(rest$count),
      pct      = sum(rest$pct)
    )
    if (show_other) {
      pooled_row <- data.frame(
        sg_found = sprintf("Other (%d labels)", pool_info$n_labels),
        count    = pool_info$count,
        pct      = pool_info$pct,
        stringsAsFactors = FALSE
      )
      keep <- rbind(top, pooled_row)
    } else {
      keep <- top
    }
  }

  n_other <- n_found - sum(counts$count[counts$pct >= current])

  # ---------------------------------------------------------------------------
  # Display: wrap long labels, set factor order, build bar text
  # ---------------------------------------------------------------------------
  keep$sg_label <- vapply(keep$sg_found, function(x)
    paste(strwrap(x, width = wrap_width), collapse = "\n"), character(1))
  keep$sg_label <- factor(keep$sg_label, levels = rev(keep$sg_label))
  keep$bar_label <- sprintf("%d (%.1f%%)", keep$count, keep$pct)

  # ---------------------------------------------------------------------------
  # Subtitle: find rate + visible / passing / pooled counts + halving note
  # ---------------------------------------------------------------------------
  sub <- sprintf(
    "Found: %d / %d sims | Showing %s %d of %d labels \u2265 %g%%",
    n_found, n_total,
    if (!is.null(top_k)) "top" else "subgroups",
    if (!is.null(top_k)) min(as.integer(top_k), n_distinct_passing)
    else n_distinct_passing,
    if (!is.null(top_k)) n_distinct_passing else n_distinct_total,
    current
  )
  if (!is.null(pool_info)) {
    verb <- if (show_other) "Other" else "Suppressed"
    sub  <- paste0(sub, sprintf(" | %s: %d (%.1f%%)",
                                verb, pool_info$count, pool_info$pct))
  } else if (is.null(top_k) && n_other > 0L) {
    sub <- paste0(sub, sprintf(" | Other: %d (%.1f%%)",
                               n_other, 100 * n_other / n_found))
  }

  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  hjust_val <- if (bar_label_inside) 1.05 else -0.1
  size_val  <- if (bar_label_inside) 3.2  else 3.5
  expand_r  <- if (bar_label_inside) c(0, 0.02) else c(0, 0.15)

  p <- ggplot2::ggplot(
    keep,
    ggplot2::aes(x = .data$sg_label, y = .data$count, fill = .data$sg_label)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = .data$bar_label),
                       hjust = hjust_val, size = size_val) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = expand_r)) +
    ggplot2::labs(x = NULL, y = "Count",
                  title = title, subtitle = sub) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position   = "none",
      plot.subtitle     = ggplot2::element_text(size = 10, colour = "grey40"),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y       = ggplot2::element_text(size = 10)
    )

  if (n_halvings > 0L) {
    p <- p + ggplot2::labs(
      caption = sprintf(
        "Threshold auto-halved %d \u00d7 from %g%% to %g%% (no label cleared the initial value).",
        n_halvings, attempts[1], current
      )
    )
  }

  attr(p, "n_bars")            <- nrow(keep)
  attr(p, "effective_min_pct") <- current
  p
}


#' Placeholder panel when no subgroup bars can be drawn
#'
#' Internal helper used by \code{plot_sg_distribution()} when
#' \code{placeholder_on_empty = TRUE}.
#'
#' @param title Plot title.
#' @param msg Message rendered inside the panel.
#' @return A ggplot2 object with \code{attr(p, "n_bars") = 0L}.
#' @keywords internal
#' @noRd
.empty_sgdist_plot <- function(title, msg) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5,
                      label = msg, size = 5, hjust = 0.5) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  attr(p, "n_bars")            <- 0L
  attr(p, "effective_min_pct") <- NA_real_
  p
}


#' Adaptive Figure Height for a Subgroup-Distribution Plot
#'
#' Computes a reasonable knitr \code{fig.height} (in inches) for a
#' \code{\link{plot_sg_distribution}()} panel based on the number of
#' bars, keeping the visual quality consistent across scenarios with
#' very different numbers of distinct subgroup labels.
#'
#' @param n_bars Integer.  Number of bars in the plot, typically
#'   retrieved via \code{attr(p, "n_bars")} on a plot returned by
#'   \code{\link{plot_sg_distribution}()}.
#' @param floor_in Numeric.  Minimum height in inches.
#'   Default: \code{4.0}.
#' @param ceil_in Numeric.  Maximum height in inches.
#'   Default: \code{10.0}.
#' @param per_bar Numeric.  Height contribution per bar in inches.
#'   Default: \code{0.35}.
#' @param base_in Numeric.  Base height added to \code{per_bar * n_bars}
#'   before applying the floor / ceiling.  Default: \code{3.0}.
#'
#' @return Numeric scalar, a height in inches in \code{[floor_in,
#'   ceil_in]}.  Suitable for the \code{fig.height} chunk option.
#'
#' @examples
#' \dontrun{
#' p <- plot_sg_distribution(results, placeholder_on_empty = TRUE)
#' h <- sgdist_fig_height(attr(p, "n_bars"))
#' # Use `h` as the fig.height chunk option.
#' }
#'
#' @seealso \code{\link{plot_sg_distribution}}
#'
#' @export
sgdist_fig_height <- function(n_bars,
                              floor_in = 4.0,
                              ceil_in  = 10.0,
                              per_bar  = 0.35,
                              base_in  = 3.0) {
  n_bars <- as.integer(n_bars)
  if (length(n_bars) != 1L || is.na(n_bars)) n_bars <- 1L
  max(floor_in, min(ceil_in, base_in + per_bar * max(1L, n_bars)))
}
