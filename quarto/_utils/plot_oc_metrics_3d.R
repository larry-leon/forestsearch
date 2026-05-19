# =============================================================================
# plot_oc_metrics_3d.R --- Interactive 3D scatter of OC summary tables
# =============================================================================
# Location: quarto/_utils/plot_oc_metrics_3d.R
# Dependencies: plotly (load in the consuming .qmd's setup chunk, or rely
#               on the requireNamespace() runtime gate below).
#
# Sibling of plot_oc_metrics.R that renders three OC metrics as a true 3D
# scatter via plotly.  Each of x, y, z is a spatial axis; categorical
# columns can be mapped to color and symbol; an optional 4th metric can
# be encoded as marker size; an optional column drives hover text.
#
# HTML-only.  The returned plotly object renders inline in HTML Quarto
# output (format: html) but will not render in PDF.  For PDF, use the
# 2D bubble version in plot_oc_metrics.R.
#
# Why 3D scatter at all: with rotation, all four metric values (three
# spatial + optional size) are simultaneously legible.  Without
# rotation, static screenshots can be misleading because perspective
# alone does not resolve depth; the first thing a reader should do
# with the rendered plot is rotate it.
#
# Sourcing from a .qmd setup chunk:
#
#   source("../../../_utils/plot_oc_metrics_3d.R")   # adjust depth
#
#   plot_oc_metrics_3d(display_dt)
#
# This file is NOT part of the forestsearch package; it is an
# exploratory utility shared across simulation evaluation documents.
# =============================================================================


# Default mapping from internal column name -> display label.  Duplicated
# from plot_oc_metrics.R so each utility file is self-contained.  Keep
# the two copies synchronized when adding new label entries.
.OC_DEFAULT_LABELS_3D <- c(
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


# plot_oc_metrics_3d() ------------------------------------------------------
#
# Arguments
#   data         data.frame / data.table of OC summary rows.
#   x, y, z      column names for spatial axes (numeric).  Defaults
#                "Spec", "NPV", "Detection".
#   color_by     column used for point color, or NULL.  Default "group".
#   symbol_by    column used for point symbol (plotly's term for
#                shape), or NULL.  Default "Analysis".  plotly's
#                scatter3d supports ~8 distinct symbols; with more
#                levels the palette cycles.
#   hover_by     column used for hover-text labels, or NULL.  Default
#                "group" (hover shows the configuration label + axes).
#   size_by      column used as a 4th dimension encoded in marker size,
#                or NULL.  Default NULL (uniform marker size).
#   size_range   length-2 numeric vector giving min/max marker size in
#                pixels when size_by is non-NULL, or a length-1
#                numeric scalar giving the uniform size when size_by
#                is NULL.  Default c(6, 18).
#   axis_labels  named character vector overriding default display
#                labels (names = column names, values = labels).
#                Default NULL.
#   title        plot title (default NULL).
#
# Returns: a plotly object.  Print the object to render it (or let it
#   be the last expression in a Quarto chunk).
#
# Limitations:
#   * HTML rendering only -- will not appear in PDF output.
#   * Best for small grids (~5-30 points).  At higher counts,
#     occlusion dominates and the 2D bubble version is more readable.
#   * Static screenshots are misleading.  Rotation is essential.
#   * symbol palette is limited (~8 distinct symbols).
#
# Examples (assuming display_dt is the table built by an OC summary
# article's unified-oc-data chunk):
#
#   # Defaults: Spec x NPV x Detection, color = group, symbol = Analysis,
#   # hover shows group label.
#   plot_oc_metrics_3d(display_dt)
#
#   # Alt metric on z-axis:
#   plot_oc_metrics_3d(display_dt, z = "OR_H_hat",
#                      title = "Spec x NPV x OR(H)")
#
#   # 4D view: spatial = Spec/NPV/Detect; marker size = sensitivity
#   plot_oc_metrics_3d(display_dt, size_by = "Sen")
#
#   # Custom labels:
#   plot_oc_metrics_3d(
#     display_dt,
#     axis_labels = c(Spec = "Specificity (1 - FPR)",
#                     NPV  = "Negative predictive value")
#   )
# ---------------------------------------------------------------------------

plot_oc_metrics_3d <- function(data,
                               x           = "Spec",
                               y           = "NPV",
                               z           = "Detection",
                               color_by    = "group",
                               symbol_by   = "Analysis",
                               hover_by    = "group",
                               size_by     = NULL,
                               size_range  = c(6, 18),
                               axis_labels = NULL,
                               title       = NULL) {

  # ---- Dependency gate ---------------------------------------------------
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for plot_oc_metrics_3d().\n",
         "  Install with: install.packages(\"plotly\")",
         call. = FALSE)
  }

  # ---- Input validation --------------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'data' has zero rows; nothing to plot.", call. = FALSE)
  }
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
  .require_scalar_char(x,         "x")
  .require_scalar_char(y,         "y")
  .require_scalar_char(z,         "z")
  .require_scalar_char(color_by,  "color_by",  allow_null = TRUE)
  .require_scalar_char(symbol_by, "symbol_by", allow_null = TRUE)
  .require_scalar_char(hover_by,  "hover_by",  allow_null = TRUE)
  .require_scalar_char(size_by,   "size_by",   allow_null = TRUE)

  # Columns present?  c() drops NULL args automatically.
  required_cols <- c(x, y, z, color_by, symbol_by, hover_by, size_by)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols)) {
    stop(sprintf(
      "Column(s) not found in 'data': %s.  Available columns: %s.",
      paste(shQuote(missing_cols), collapse = ", "),
      paste(shQuote(names(data)),  collapse = ", ")),
      call. = FALSE)
  }

  # Spatial axes and size column must be numeric.
  for (axis in c("x", "y", "z")) {
    col <- get(axis)
    if (!is.numeric(data[[col]])) {
      stop(sprintf(
        "Column '%s' (mapped to %s-axis) must be numeric; got %s.",
        col, axis, paste(class(data[[col]]), collapse = "/")),
        call. = FALSE)
    }
  }
  if (!is.null(size_by) && !is.numeric(data[[size_by]])) {
    stop(sprintf(
      "Column '%s' (size_by) must be numeric; got %s.",
      size_by, paste(class(data[[size_by]]), collapse = "/")),
      call. = FALSE)
  }

  # size_range shape: length-2 when size_by is non-NULL, length-1 OK otherwise.
  if (!is.numeric(size_range) || anyNA(size_range)) {
    stop("'size_range' must be numeric and non-NA.", call. = FALSE)
  }
  if (!is.null(size_by)) {
    if (length(size_range) != 2L || size_range[1L] >= size_range[2L]) {
      stop("When 'size_by' is non-NULL, 'size_range' must be a length-2 ",
           "numeric vector with size_range[1] < size_range[2].",
           call. = FALSE)
    }
  } else {
    if (!length(size_range) %in% c(1L, 2L)) {
      stop("'size_range' must be length 1 or 2 when 'size_by' is NULL.",
           call. = FALSE)
    }
  }

  # axis_labels: named character or NULL
  if (!is.null(axis_labels)) {
    if (!is.character(axis_labels) || is.null(names(axis_labels)) ||
        any(names(axis_labels) == "" | is.na(names(axis_labels)))) {
      stop("'axis_labels' must be a named character vector ",
           "(names = column names, values = display labels).",
           call. = FALSE)
    }
  }

  # ---- Resolve display labels --------------------------------------------
  .lookup_label <- function(col) {
    if (!is.null(axis_labels) && col %in% names(axis_labels)) {
      return(unname(axis_labels[col]))
    }
    if (col %in% names(.OC_DEFAULT_LABELS_3D)) {
      return(unname(.OC_DEFAULT_LABELS_3D[col]))
    }
    col
  }
  lab_x <- .lookup_label(x)
  lab_y <- .lookup_label(y)
  lab_z <- .lookup_label(z)

  # ---- Build the plotly call --------------------------------------------
  # Use formula syntax (~`col`) for plotly NSE.  Backticks let arbitrary
  # column names through safely.
  .f <- function(col) stats::as.formula(paste0("~`", col, "`"))

  # Marker spec: uniform size or rescaled vector
  marker_spec <- if (!is.null(size_by)) {
    sv  <- data[[size_by]]
    rng <- range(sv, na.rm = TRUE)
    rescaled <- if (diff(rng) > 0) {
      size_range[1L] + (sv - rng[1L]) / diff(rng) * diff(size_range)
    } else {
      rep(mean(size_range), length(sv))
    }
    list(size = rescaled, sizemode = "diameter",
         line = list(width = 0.5, color = "#444444"))
  } else {
    uniform_size <- if (length(size_range) == 1L) size_range else mean(size_range)
    list(size = uniform_size,
         line = list(width = 0.5, color = "#444444"))
  }

  args <- list(
    data   = data,
    type   = "scatter3d",
    mode   = "markers",
    x      = .f(x),
    y      = .f(y),
    z      = .f(z),
    marker = marker_spec
  )
  if (!is.null(color_by))  args$color  <- .f(color_by)
  if (!is.null(symbol_by)) args$symbol <- .f(symbol_by)
  if (!is.null(hover_by)) {
    args$text      <- .f(hover_by)
    args$hoverinfo <- "text+x+y+z"
  } else {
    args$hoverinfo <- "x+y+z"
  }

  p <- do.call(plotly::plot_ly, args)

  # ---- Layout: axis titles, plot title, legend orientation ---------------
  p <- plotly::layout(
    p,
    title = title,
    scene = list(
      xaxis = list(title = lab_x),
      yaxis = list(title = lab_y),
      zaxis = list(title = lab_z),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.2))
    ),
    legend = list(orientation = "v", x = 1.02, y = 1)
  )

  p
}
