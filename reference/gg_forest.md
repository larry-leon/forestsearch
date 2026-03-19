# ggplot2 / patchwork forest plot

Creates a publication-quality forest plot using ggplot2 for the CI panel
and patchwork to assemble label and annotation columns alongside it.
Unlike forestploter, `fig.height` maps directly to row density —
`row_height = fig.height / n_rows` with no hidden scaling.

## Usage

``` r
gg_forest(
  subgroups,
  est,
  lo,
  hi,
  cat_vec = NULL,
  cat_colours = NULL,
  annot = NULL,
  ref_line = 1,
  vert_lines = NULL,
  ref_col = "firebrick",
  ref_lty = "dashed",
  vert_col = "grey50",
  vert_lty = "dotted",
  xlim = NULL,
  ticks_at = NULL,
  tick_labels = NULL,
  xlog = TRUE,
  xlab = "Hazard Ratio",
  title = NULL,
  subtitle = NULL,
  footnote = NULL,
  point_size = 2.5,
  line_size = 0.8,
  point_shape = 21,
  base_size = 11,
  widths = NULL,
  row_expand = 0.6
)
```

## Arguments

- subgroups:

  Character vector of subgroup names (displayed top to bottom).

- est:

  Numeric vector of point estimates (median HR or similar).

- lo:

  Numeric vector of lower bounds (e.g. 1st percentile ECI).

- hi:

  Numeric vector of upper bounds (e.g. 99th percentile ECI).

- cat_vec:

  Optional character vector of category labels (one per row). Used to
  colour CI lines and label text.

- cat_colours:

  Optional named character vector mapping category labels to colours.
  Defaults to grey for all rows.

- annot:

  Optional named list of character vectors, one per annotation column.
  Names become column headers. Each vector must match
  `length(subgroups)`.

- ref_line:

  Numeric. X position of the primary reference line (default 1). Drawn
  as a dashed red line.

- vert_lines:

  Numeric vector. X positions of secondary vertical lines (default
  NULL). Drawn as dotted grey lines.

- ref_col:

  Colour of the primary reference line (default "firebrick").

- ref_lty:

  Line type of the primary reference line (default "dashed").

- vert_col:

  Colour of secondary vertical lines (default "grey50").

- vert_lty:

  Line type of secondary vertical lines (default "dotted").

- xlim:

  Numeric vector length 2. X-axis limits for the CI panel.

- ticks_at:

  Numeric vector. X-axis tick positions.

- tick_labels:

  Character vector. Custom tick labels (default:
  as.character(ticks_at)).

- xlog:

  Logical. If TRUE (default), x-axis on log scale.

- xlab:

  Character. X-axis label (default "Hazard Ratio").

- title:

  Character. Overall plot title (default NULL).

- subtitle:

  Character. Plot subtitle (default NULL).

- footnote:

  Character. Footnote appended below the CI panel (default NULL).

- point_size:

  Numeric. Size of point estimate symbol (default 2.5).

- line_size:

  Numeric. Line width of CI segments (default 0.8).

- point_shape:

  Integer. pch for point estimates (default 21, filled circle).

- base_size:

  Numeric. ggplot2 base font size in pt (default 11). Controls all text
  — increase to make the plot larger; no other knob needed.

- widths:

  Numeric vector. Relative patchwork column widths: c(label, ci,
  annot_1, annot_2, …). Default: c(3.5, 5, rep(1, n_annot)).

- row_expand:

  Numeric. Extra space above and below row range on y-axis, in row units
  (default 0.6).

## Value

A patchwork object. Render with
[`print()`](https://rdrr.io/r/base/print.html) or
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). Control
dimensions entirely via knitr chunk options `fig.width` / `fig.height`:
row height = fig.height / n_rows.
