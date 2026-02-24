# =============================================================================
# gg_forest.R
# ggplot2 + patchwork forest plot for forestsearch vignettes
#
# Architecture
# ─────────────────────────────────────────────────────────────────────────────
# Three separate ggplots stitched together with patchwork::wrap_plots():
#
#   [  label_panel  |      ci_panel      | annot_panel_1 | annot_panel_2 | … ]
#
#   label_panel   — subgroup names; geom_text on a blank canvas
#   ci_panel      — confidence intervals; geom_errorbarh + geom_point, log x
#   annot_panel_k — one per annotation column; geom_text on blank canvas
#
# Because all panels share the same discrete y factor, patchwork aligns rows
# perfectly.  fig.height controls row density directly:
#
#     row_height_in = fig.height / n_rows   (no hidden scaling)
#
# Usage
# ─────────────────────────────────────────────────────────────────────────────
# p <- gg_forest(
#   subgroups  = names_ok,
#   est        = hr_q[ok, 2],
#   lo         = hr_q[ok, 1],
#   hi         = hr_q[ok, 3],
#   cat_vec    = cat_ok,
#   cat_colours = cat_cols,
#   annot      = list("Mean N" = mean_n_chr, "Pr(UB<1)" = pr_chr),
#   ref_line   = 0.70,
#   vert_lines = 1.00,
#   xlim       = c(0.15, 3.5),
#   ticks_at   = c(0.20, 0.50, 0.70, 1.00, 1.50, 2.50),
#   xlab       = "Hazard Ratio",
#   base_size  = 11,
#   widths     = c(3, 5, 1, 1)   # label | CI | annot1 | annot2
# )
# print(p)   # or plot(p)
#
# =============================================================================

#' ggplot2 / patchwork forest plot
#'
#' Creates a publication-quality forest plot using ggplot2 for the CI panel
#' and patchwork to assemble label and annotation columns alongside it.
#' Unlike forestploter, `fig.height` maps directly to row density —
#' `row_height = fig.height / n_rows` with no hidden scaling.
#'
#' @param subgroups Character vector of subgroup names (displayed top to bottom).
#' @param est Numeric vector of point estimates (median HR or similar).
#' @param lo  Numeric vector of lower bounds (e.g. 1st percentile ECI).
#' @param hi  Numeric vector of upper bounds (e.g. 99th percentile ECI).
#' @param cat_vec Optional character vector of category labels (one per row).
#'   Used to colour CI lines and label text.
#' @param cat_colours Optional named character vector mapping category labels
#'   to colours. Defaults to grey for all rows.
#' @param annot Optional named list of character vectors, one per annotation
#'   column. Names become column headers. Each vector must match `length(subgroups)`.
#' @param ref_line Numeric. X position of the primary reference line (default 1).
#'   Drawn as a dashed red line.
#' @param vert_lines Numeric vector. X positions of secondary vertical lines
#'   (default NULL). Drawn as dotted grey lines.
#' @param ref_col  Colour of the primary reference line (default "firebrick").
#' @param ref_lty  Line type of the primary reference line (default "dashed").
#' @param vert_col Colour of secondary vertical lines (default "grey50").
#' @param vert_lty Line type of secondary vertical lines (default "dotted").
#' @param xlim Numeric vector length 2. X-axis limits for the CI panel.
#' @param ticks_at Numeric vector. X-axis tick positions.
#' @param tick_labels Character vector. Custom tick labels (default: as.character(ticks_at)).
#' @param xlog Logical. If TRUE (default), x-axis on log scale.
#' @param xlab Character. X-axis label (default "Hazard Ratio").
#' @param title   Character. Overall plot title (default NULL).
#' @param subtitle Character. Plot subtitle (default NULL).
#' @param footnote Character. Footnote appended below the CI panel (default NULL).
#' @param point_size Numeric. Size of point estimate symbol (default 2.5).
#' @param line_size  Numeric. Line width of CI segments (default 0.8).
#' @param point_shape Integer. pch for point estimates (default 21, filled circle).
#' @param base_size  Numeric. ggplot2 base font size in pt (default 11).
#'   Controls all text — increase to make the plot larger; no other knob needed.
#' @param widths Numeric vector. Relative patchwork column widths:
#'   c(label, ci, annot_1, annot_2, …). Default: c(3.5, 5, rep(1, n_annot)).
#' @param row_expand Numeric. Extra space above and below row range on y-axis,
#'   in row units (default 0.6).
#'
#' @return A patchwork object. Render with `print()` or `plot()`.
#'   Control dimensions entirely via knitr chunk options `fig.width` /
#'   `fig.height`:  row height = fig.height / n_rows.
#'
#' @examples
#' \dontrun{
#' # Recommended fig.height: n_rows * 0.45 + 1.5 (for title/axis overhead)
#' # e.g. 20 rows -> fig.height = 20 * 0.45 + 1.5 = 10.5
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_text
#'   geom_vline scale_x_log10 scale_x_continuous scale_y_discrete
#'   scale_colour_identity scale_fill_identity labs theme theme_minimal
#'   element_blank element_text element_line margin coord_cartesian
#'   annotation_logticks
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
#' @export
gg_forest <- function(
    subgroups,
    est,
    lo,
    hi,
    cat_vec       = NULL,
    cat_colours   = NULL,
    annot         = NULL,
    ref_line      = 1,
    vert_lines    = NULL,
    ref_col       = "firebrick",
    ref_lty       = "dashed",
    vert_col      = "grey50",
    vert_lty      = "dotted",
    xlim          = NULL,
    ticks_at      = NULL,
    tick_labels   = NULL,
    xlog          = TRUE,
    xlab          = "Hazard Ratio",
    title         = NULL,
    subtitle      = NULL,
    footnote      = NULL,
    point_size    = 2.5,
    line_size     = 0.8,
    point_shape   = 21,
    base_size     = 11,
    widths        = NULL,
    row_expand    = 0.6
) {

  # ── Input validation ───────────────────────────────────────────────────────
  n <- length(subgroups)
  stopifnot(length(est) == n, length(lo) == n, length(hi) == n)

  if (!is.null(cat_vec))    stopifnot(length(cat_vec) == n)
  if (!is.null(annot))      lapply(annot, function(v) stopifnot(length(v) == n))

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' required.  Install with: install.packages('patchwork')")

  # ── Colour vector (one colour per row) ────────────────────────────────────
  default_col <- "grey30"
  if (!is.null(cat_vec) && !is.null(cat_colours)) {
    row_cols <- dplyr::coalesce(cat_colours[cat_vec], default_col)
  } else {
    row_cols <- rep(default_col, n)
  }

  # ── Shared y factor (rows displayed top-to-bottom) ────────────────────────
  # Factor levels: last level = bottom of plot in ggplot's coord system.
  # We want subgroups[1] at the TOP, so reverse the factor.
  y_levels <- rev(subgroups)   # ggplot draws level 1 at y = 1 (bottom)

  df <- data.frame(
    subgroup = factor(subgroups, levels = y_levels),
    est      = est,
    lo       = lo,
    hi       = hi,
    row_col  = row_cols,
    stringsAsFactors = FALSE
  )

  # ── Axis limits and breaks ─────────────────────────────────────────────────
  if (is.null(xlim)) {
    padding <- if (xlog) 1.15 else diff(range(c(lo, hi), na.rm = TRUE)) * 0.1
    xlim <- if (xlog)
      c(min(lo, na.rm = TRUE) / padding, max(hi, na.rm = TRUE) * padding)
    else
      c(min(lo, na.rm = TRUE) - padding, max(hi, na.rm = TRUE) + padding)
  }

  if (is.null(ticks_at)) {
    ticks_at <- if (xlog)
      c(0.25, 0.5, 1.0, 2.0, 4.0)
    else
      pretty(xlim, n = 5)
  }
  if (is.null(tick_labels)) tick_labels <- as.character(ticks_at)

  # ── Y-axis range (with row_expand breathing room) ─────────────────────────
  y_lo <- 1 - row_expand
  y_hi <- n + row_expand

  # ── Shared minimal theme for blank side panels ────────────────────────────
  theme_blank <- function(bs = base_size) {
    theme_minimal(base_size = bs) +
    theme(
      axis.title   = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      plot.margin  = margin(2, 2, 2, 2, "pt")
    )
  }

  # ── 1. LABEL PANEL ──────────────────────────────────────────────────────────
  df_lbl <- data.frame(
    subgroup = factor(subgroups, levels = y_levels),
    col      = row_cols,
    stringsAsFactors = FALSE
  )

  p_label <- ggplot2::ggplot(df_lbl,
      ggplot2::aes(x = 0, y = subgroup, label = subgroup, colour = col)) +
    ggplot2::geom_text(
      hjust     = 1,
      size      = base_size / ggplot2::.pt,   # convert pt to ggplot size units
      fontface  = ifelse(df_lbl$col == (if (!is.null(cat_colours)) cat_colours["ITT"]
                                        else default_col), "bold", "plain")
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(limits = c(-1, 0.05), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi), clip = "off") +
    ggplot2::labs(title = "Subgroup") +
    theme_blank(base_size) +
    theme(
      plot.title = ggplot2::element_text(
        size = base_size, face = "bold", hjust = 1,
        margin = margin(b = 4, unit = "pt")
      ),
      plot.margin = margin(2, 4, 2, 2, "pt")
    )

  # ── 2. CI PANEL ─────────────────────────────────────────────────────────────
  # Build scale
  if (xlog) {
    x_scale <- ggplot2::scale_x_log10(
      limits = xlim,
      breaks = ticks_at,
      labels = tick_labels,
      expand = c(0, 0)
    )
  } else {
    x_scale <- ggplot2::scale_x_continuous(
      limits = xlim,
      breaks = ticks_at,
      labels = tick_labels,
      expand = c(0, 0)
    )
  }

  p_ci <- ggplot2::ggplot(df, ggplot2::aes(y = subgroup)) +
    # Reference lines
    ggplot2::geom_vline(
      xintercept = ref_line,
      colour     = ref_col,
      linetype   = ref_lty,
      linewidth  = 0.7
    )

  if (!is.null(vert_lines)) {
    p_ci <- p_ci +
      ggplot2::geom_vline(
        xintercept = vert_lines,
        colour     = vert_col,
        linetype   = vert_lty,
        linewidth  = 0.55
      )
  }

  p_ci <- p_ci +
    # CI lines (whiskers)
    ggplot2::geom_segment(
      ggplot2::aes(x = lo, xend = hi, yend = subgroup, colour = row_col),
      linewidth = line_size,
      lineend   = "round"
    ) +
    # Point estimates
    ggplot2::geom_point(
      ggplot2::aes(x = est, fill = row_col),
      shape  = point_shape,
      colour = "white",
      size   = point_size,
      stroke = 0.7
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    x_scale +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi), clip = "on") +
    ggplot2::labs(
      x        = xlab,
      title    = if (xlog) paste0(xlab, " (log scale)") else xlab,
      subtitle = subtitle,
      caption  = footnote
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.title.y      = element_blank(),
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.text.x       = ggplot2::element_text(size = base_size * 0.85),
      axis.title.x      = ggplot2::element_text(size = base_size * 0.9,
                                                 margin = margin(t = 4, unit = "pt")),
      panel.grid.major.y = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey88", linewidth = 0.3),
      plot.title   = ggplot2::element_text(size = base_size, face = "bold",
                                            margin = margin(b = 4, unit = "pt")),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.8, colour = "grey45"),
      plot.caption  = ggplot2::element_text(size = base_size * 0.75, colour = "grey45",
                                             hjust = 0, margin = margin(t = 6, unit = "pt")),
      plot.margin  = margin(2, 8, 2, 4, "pt")
    )

  # ── 3. ANNOTATION PANELS ─────────────────────────────────────────────────
  p_annots <- list()

  if (!is.null(annot)) {
    for (col_name in names(annot)) {
      df_a <- data.frame(
        subgroup = factor(subgroups, levels = y_levels),
        value    = as.character(annot[[col_name]]),
        stringsAsFactors = FALSE
      )

      pa <- ggplot2::ggplot(df_a,
          ggplot2::aes(x = 0.5, y = subgroup, label = value)) +
        ggplot2::geom_text(
          hjust   = 0.5,
          size    = base_size / ggplot2::.pt * 0.9,
          colour  = "grey20"
        ) +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ggplot2::scale_y_discrete(drop = FALSE) +
        ggplot2::coord_cartesian(ylim = c(y_lo, y_hi), clip = "off") +
        ggplot2::labs(title = col_name) +
        theme_blank(base_size) +
        theme(
          plot.title = ggplot2::element_text(
            size = base_size, face = "bold", hjust = 0.5,
            margin = margin(b = 4, unit = "pt")
          ),
          plot.margin = margin(2, 4, 2, 4, "pt")
        )

      p_annots[[col_name]] <- pa
    }
  }

  # ── 4. ASSEMBLE WITH PATCHWORK ─────────────────────────────────────────────
  n_annot <- length(p_annots)

  if (is.null(widths)) {
    widths <- c(3.5, 5, rep(1.1, n_annot))
  }

  panels <- c(list(p_label, p_ci), p_annots)

  assembled <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_layout(widths = widths)

  if (!is.null(title)) {
    assembled <- assembled +
      patchwork::plot_annotation(
        title = title,
        theme = theme(
          plot.title = ggplot2::element_text(
            size = base_size * 1.2, face = "bold",
            margin = margin(b = 6, unit = "pt")
          )
        )
      )
  }

  assembled
}
