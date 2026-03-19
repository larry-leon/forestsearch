#' Plot Distribution of Identified Subgroups
#'
#' Bar chart of subgroups identified across simulations, filtered to
#' those appearing in at least \code{min_pct} of the found simulations.
#'
#' @param results data.table from \code{mrct_region_sims()}
#' @param min_pct Numeric. Minimum percentage threshold for display (default: 5)
#' @param title Character. Plot title. Default: "Distribution of Identified Subgroups"
#' @param wrap_width Integer. Character width for wrapping long subgroup labels.
#'   Default: 25
#'
#' @keywords internal
#' @return A ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_col geom_text coord_flip
#'   scale_y_continuous expansion labs theme_minimal theme element_text

plot_sg_distribution <- function(results,
                                 min_pct = 5,
                                 title = "Distribution of Identified Subgroups",
                                 wrap_width = 25) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Please install it.")
  }

  # Filter to simulations where a subgroup was found
  found <- results[results$any_found == 1, ]
  n_found <- nrow(found)
  n_total <- nrow(results)

  if (n_found == 0) {
    message("No subgroups identified in any simulation.")
    return(invisible(NULL))
  }

  # Count each subgroup and compute percentage
  sg_counts <- as.data.frame(table(sg_found = found$sg_found), stringsAsFactors = FALSE)
  names(sg_counts) <- c("sg_found", "count")
  sg_counts$pct <- 100 * sg_counts$count / n_found

  # Filter to subgroups meeting min_pct threshold
  sg_keep <- sg_counts[sg_counts$pct >= min_pct, ]
  n_other <- n_found - sum(sg_keep$count)

  if (nrow(sg_keep) == 0) {
    message(sprintf("No subgroup appears in >= %.0f%% of found simulations.", min_pct))
    return(invisible(NULL))
  }

  # Sort by count descending
  sg_keep <- sg_keep[order(-sg_keep$count), ]

  # Wrap long labels
  sg_keep$sg_label <- sapply(sg_keep$sg_found, function(x) {
    paste(strwrap(x, width = wrap_width), collapse = "\n")
  })

  # Set factor order for plot
  sg_keep$sg_label <- factor(sg_keep$sg_label,
                             levels = rev(sg_keep$sg_label))

  # Build subtitle
  subtitle_text <- sprintf(
    "Found: %d / %d sims | Showing subgroups >= %.0f%% (%d of %d distinct)",
    n_found, n_total, min_pct, nrow(sg_keep), nrow(sg_counts)
  )

  if (n_other > 0) {
    subtitle_text <- paste0(subtitle_text,
                            sprintf(" | Other: %d (%.1f%%)", n_other, 100 * n_other / n_found))
  }

  # Bar labels: count (pct%)
  sg_keep$bar_label <- sprintf("%d (%.1f%%)", sg_keep$count, sg_keep$pct)

  # Plot
  p <- ggplot2::ggplot(sg_keep, ggplot2::aes(x = sg_label, y = count, fill = sg_label)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = bar_label), hjust = -0.1, size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      x = NULL,
      y = "Count",
      title = title,
      subtitle = subtitle_text
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40"),
      panel.grid.major.y = ggplot2::element_blank()
    )

  p
}
