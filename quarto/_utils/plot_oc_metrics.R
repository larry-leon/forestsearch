# =============================================================================
# plot_oc_metrics.R --- Graphical display of OC summary tables
# =============================================================================
# Location: quarto/_utils/plot_oc_metrics.R
# Dependencies: ggplot2 (load in the consuming .qmd's setup chunk).
#
# Renders a 2D bubble chart of an operating-characteristic summary table
# of the kind produced by the cross-cut simulation summary articles
# (e.g. actg175_binary_m1_harm_fs1_summary.qmd).  The x and y axes are
# spatial; the third metric (z) is encoded as point size.  Categorical
# columns can be mapped to color and shape, and an optional column can
# be used to label individual points.
#
# Sourcing from a .qmd setup chunk:
#
#   library(ggplot2)
#   source("../../_utils/plot_oc_metrics.R")    # adjust relative depth
#
#   p <- plot_oc_metrics(display_dt)
#   print(p)
#
# This file is NOT part of the forestsearch package; it is an exploratory
# utility shared across simulation evaluation documents.  If/when this
# utility is promoted to the package, restore the roxygen header and
# place under R/ with the @export and @keywords internal tags.
# =============================================================================


# Default mapping from internal column name -> display label.  Mirrors
# the gt::cols_label() mapping used in the OC summary articles.
# Unicode: U+0124 LATIN CAPITAL LETTER H WITH CIRCUMFLEX,
# U+00D4 LATIN CAPITAL LETTER O WITH CIRCUMFLEX,
# U+1D9C MODIFIER LETTER SMALL C.
.OC_DEFAULT_LABELS <- c(
  Detection   = "Detect",
  Sen         = "Sen",
  Spec        = "Spec",
  PPV         = "PPV",
  NPV         = "NPV",
  OR_H_hat    = "\u00d4R(\u0124)",          # ÔR(Ĥ)
  OR_Hc_hat   = "\u00d4R(\u0124\u1d9c)",    # ÔR(Ĥᶜ)
  OR_ITT      = "\u00d4R(ITT)",             # ÔR(ITT)
  Size_H_mean = "|\u0124| mean",            # |Ĥ| mean
  Size_H_min  = "|\u0124| min",
  Size_H_max  = "|\u0124| max",
  N_sims      = "N"
)


# plot_oc_metrics() ---------------------------------------------------------
#
# Arguments
#   data         data.frame / data.table of OC summary rows
#                (one row per configuration x analysis method).
#   x            column in `data` for x-axis (numeric).  Default "Spec".
#   y            column in `data` for y-axis (numeric).  Default "NPV".
#   z            column in `data` encoded as point size (numeric).
#                Default "Detection".
#   color_by     column used for point color, or NULL.  Default "group".
#   shape_by     column used for point shape, or NULL.  Default "Analysis".
#   label_by     column used for inline point labels, or NULL.  Default NULL.
#   size_limits  length-2 numeric vector or NULL.  Pass c(0, 1) when `z`
#                is a proportion-scaled metric (Detection, Sen, Spec,
#                PPV, NPV) so bubble size is absolutely interpretable.
#                Default NULL (data range).
#   size_range   length-2 numeric vector: min/max point size (mm).
#                Default c(2, 10).
#   axis_labels  named character vector overriding the default display
#                labels.  Names = column names, values = labels.  Default
#                NULL.  Falls back to the built-in default map, then to
#                the bare column name.
#   pareto_col   character scalar naming a logical column in `data`
#                that indicates Pareto-frontier membership, or NULL
#                to disable.  When non-NULL, points where the flag
#                is TRUE are drawn with full opacity and a thicker
#                stroke; points where the flag is FALSE (or NA) are
#                dimmed.  See compute_oc_pareto_flag() for building
#                this column.  Default NULL.
#   title        plot title (default NULL).
#   subtitle     plot subtitle (default NULL).
#
# Returns: a ggplot object.
#
# Default display labels (column -> label):
#   Detection   -> "Detect"
#   Sen, Spec, PPV, NPV -> unchanged
#   OR_H_hat    -> ÔR(Ĥ)
#   OR_Hc_hat   -> ÔR(Ĥᶜ)
#   OR_ITT      -> ÔR(ITT)
#   Size_H_*    -> |Ĥ| mean/min/max
#   anything else -> column name unchanged
#
# Why bubble chart, not 3D scatter: a 2D bubble renders reliably in
# static HTML reports without perspective ambiguity or occlusion, and
# avoids a plotly dependency.  To get an interactive view, install
# plotly and call plotly::ggplotly(plot_oc_metrics(...)).
#
# Color/shape with many configurations: if the OC summary spans many
# focus x fs cells, color_by = "group" gives a crowded legend.  Two
# recipes: (i) collapse the group label (e.g. data$focus_only <-
# gsub(" .*", "", data$group); color_by = "focus_only"), or
# (ii) facet the returned plot via + facet_wrap(~ fs).
#
# Examples (assuming display_dt is the table built by an OC summary
# article's unified-oc-data chunk):
#
#   p_default <- plot_oc_metrics(display_dt)
#   print(p_default)
#
#   # Sensitivity vs. PPV with effect estimate as bubble size:
#   p_sen <- plot_oc_metrics(
#     display_dt,
#     x = "Sen", y = "PPV", z = "OR_H_hat",
#     size_limits = NULL                  # OR is on the ratio scale
#   )
#
#   # With per-point labels (good for small grids):
#   p_lab <- plot_oc_metrics(
#     display_dt,
#     label_by = "group",
#     size_limits = c(0, 1)
#   )
#
#   # Custom axis labels:
#   p_relabel <- plot_oc_metrics(
#     display_dt,
#     axis_labels = c(Spec      = "Specificity (1 - FPR)",
#                     NPV       = "Negative predictive value",
#                     Detection = "Detection rate")
#   )
#
#   # Convert to interactive plotly view:
#   if (requireNamespace("plotly", quietly = TRUE)) {
#     plotly::ggplotly(plot_oc_metrics(display_dt, label_by = "group"))
#   }
# ---------------------------------------------------------------------------

plot_oc_metrics <- function(data,
                            x           = "Spec",
                            y           = "NPV",
                            z           = "Detection",
                            color_by    = "group",
                            shape_by    = "Analysis",
                            label_by    = NULL,
                            size_limits = NULL,
                            size_range  = c(2, 10),
                            axis_labels = NULL,
                            pareto_col  = NULL,
                            title       = NULL,
                            subtitle    = NULL) {

  # ---- Input validation ---------------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'data' has zero rows; nothing to plot.", call. = FALSE)
  }

  # Coerce data.table -> data.frame for ggplot2 (it accepts both, but
  # this normalizes column-access behaviour).
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  .require_scalar_char <- function(arg, name, allow_null = FALSE) {
    if (is.null(arg)) {
      if (allow_null) return(invisible(TRUE))
      stop(sprintf("'%s' must not be NULL.", name), call. = FALSE)
    }
    if (!is.character(arg) || length(arg) != 1L || is.na(arg)) {
      stop(sprintf("'%s' must be a single non-NA character string.", name),
           call. = FALSE)
    }
    invisible(TRUE)
  }
  .require_scalar_char(x, "x")
  .require_scalar_char(y, "y")
  .require_scalar_char(z, "z")
  .require_scalar_char(color_by,   "color_by",   allow_null = TRUE)
  .require_scalar_char(shape_by,   "shape_by",   allow_null = TRUE)
  .require_scalar_char(label_by,   "label_by",   allow_null = TRUE)
  .require_scalar_char(pareto_col, "pareto_col", allow_null = TRUE)

  # Columns present?  c() drops NULL args automatically.
  required_cols <- c(x, y, z, color_by, shape_by, label_by, pareto_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols)) {
    stop(sprintf(
      "Column(s) not found in 'data': %s.  Available columns: %s.",
      paste(shQuote(missing_cols), collapse = ", "),
      paste(shQuote(names(data)),  collapse = ", ")),
      call. = FALSE)
  }

  # x / y / z must be numeric.
  for (axis in c("x", "y", "z")) {
    col <- get(axis)
    if (!is.numeric(data[[col]])) {
      stop(sprintf(
        "Column '%s' (mapped to %s-axis) must be numeric; got %s.",
        col, axis, paste(class(data[[col]]), collapse = "/")),
        call. = FALSE)
    }
  }

  # pareto_col must be logical when supplied.
  if (!is.null(pareto_col) && !is.logical(data[[pareto_col]])) {
    stop(sprintf(
      "Column '%s' (pareto_col) must be logical; got %s.",
      pareto_col, paste(class(data[[pareto_col]]), collapse = "/")),
      call. = FALSE)
  }

  # size_limits, size_range shape checks
  if (!is.null(size_limits)) {
    if (!is.numeric(size_limits) || length(size_limits) != 2L ||
        anyNA(size_limits) || size_limits[1L] >= size_limits[2L]) {
      stop("'size_limits' must be a length-2 numeric vector with limits[1] < limits[2].",
           call. = FALSE)
    }
  }
  if (!is.numeric(size_range) || length(size_range) != 2L ||
      anyNA(size_range) || size_range[1L] >= size_range[2L]) {
    stop("'size_range' must be a length-2 numeric vector with range[1] < range[2].",
         call. = FALSE)
  }

  # axis_labels: named character or NULL
  if (!is.null(axis_labels)) {
    if (!is.character(axis_labels) || is.null(names(axis_labels)) ||
        any(names(axis_labels) == "" | is.na(names(axis_labels)))) {
      stop("'axis_labels' must be a named character vector (names = column names, values = display labels).",
           call. = FALSE)
    }
  }

  # ---- Resolve display labels --------------------------------------------
  .lookup_label <- function(col) {
    if (!is.null(axis_labels) && col %in% names(axis_labels)) {
      return(unname(axis_labels[col]))
    }
    if (col %in% names(.OC_DEFAULT_LABELS)) {
      return(unname(.OC_DEFAULT_LABELS[col]))
    }
    col
  }
  lab_x     <- .lookup_label(x)
  lab_y     <- .lookup_label(y)
  lab_z     <- .lookup_label(z)
  lab_color <- if (!is.null(color_by)) .lookup_label(color_by) else NULL
  lab_shape <- if (!is.null(shape_by)) .lookup_label(shape_by) else NULL

  # ---- Build the plot ----------------------------------------------------
  # Each aes() captures the function's environment, so the local
  # character vars x, y, z, color_by, shape_by, label_by are visible
  # at render time -- the standard ggplot2 idiom for dynamic columns.
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[x]],
      y = .data[[y]]
    )
  )

  # Point layer: branch on which optional aesthetics are active so
  # inactive ones are omitted entirely (avoids NULL-aes legend quirks).
  point_mapping <- if (!is.null(color_by) && !is.null(shape_by)) {
    ggplot2::aes(size   = .data[[z]],
                 colour = .data[[color_by]],
                 shape  = .data[[shape_by]])
  } else if (!is.null(color_by)) {
    ggplot2::aes(size   = .data[[z]],
                 colour = .data[[color_by]])
  } else if (!is.null(shape_by)) {
    ggplot2::aes(size  = .data[[z]],
                 shape = .data[[shape_by]])
  } else {
    ggplot2::aes(size = .data[[z]])
  }

  # Point layer(s).  When pareto_col is NULL: single layer (original
  # behaviour).  When non-NULL: split into a dimmed dominated-points
  # layer (alpha 0.30, stroke 0.4) and an emphasised frontier layer
  # (alpha 0.95, stroke 1.2).  NA entries in the flag column are
  # treated as FALSE (dimmed) so non-evaluable rows don't quietly
  # appear on the frontier.
  if (is.null(pareto_col)) {
    p <- p + ggplot2::geom_point(mapping = point_mapping,
                                 alpha = 0.85, stroke = 0.4)
  } else {
    is_front <- data[[pareto_col]]
    is_front[is.na(is_front)] <- FALSE
    frontier_data  <- data[is_front,  , drop = FALSE]
    dominated_data <- data[!is_front, , drop = FALSE]

    if (nrow(dominated_data) > 0L) {
      p <- p + ggplot2::geom_point(
        data    = dominated_data,
        mapping = point_mapping,
        alpha   = 0.30,
        stroke  = 0.4
      )
    }
    if (nrow(frontier_data) > 0L) {
      p <- p + ggplot2::geom_point(
        data    = frontier_data,
        mapping = point_mapping,
        alpha   = 0.95,
        stroke  = 1.2
      )
    }
  }

  # Size scale
  size_scale_args <- list(
    name  = lab_z,
    range = size_range
  )
  if (!is.null(size_limits)) size_scale_args$limits <- size_limits
  p <- p + do.call(ggplot2::scale_size_continuous, size_scale_args)

  # Optional labels
  if (!is.null(label_by)) {
    label_mapping <- ggplot2::aes(label = .data[[label_by]])
    p <- p + ggplot2::geom_text(
      mapping       = label_mapping,
      size          = 3,
      hjust         = -0.15,
      vjust         = -0.15,
      check_overlap = TRUE,
      show.legend   = FALSE
    )
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(
      x        = lab_x,
      y        = lab_y,
      colour   = lab_color,
      shape    = lab_shape,
      title    = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      legend.box      = "vertical",
      plot.title      = ggplot2::element_text(face = "bold"),
      plot.subtitle   = ggplot2::element_text(size = 10)
    ) +
    ggplot2::guides(
      colour = if (!is.null(color_by))
        ggplot2::guide_legend(order = 1, override.aes = list(size = 4))
        else "none",
      shape  = if (!is.null(shape_by))
        ggplot2::guide_legend(order = 2, override.aes = list(size = 4))
        else "none",
      size   = ggplot2::guide_legend(order = 3)
    )

  p
}
