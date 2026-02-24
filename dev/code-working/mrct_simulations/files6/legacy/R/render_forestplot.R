# =============================================================================
# render_forestplot.R - Helper for Rendering ForestSearch Forest Plots
# =============================================================================
#
# Provides functions to create custom themes and render forestploter plots
# with size control in Quarto/RMarkdown documents.
#
# =============================================================================

#' Create Forest Plot Theme with Size Controls
#'
#' Creates a forestploter theme with parameters that control overall plot
#' sizing and appearance. This is the primary way to control how large
#' the forest plot renders.
#'
#' @param base_size Numeric. Base font size in points. This is the primary
#'   scaling parameter - increasing it will proportionally scale all fonts,
#'   row padding, and line widths. Default: 10.
#' @param scale Numeric. Additional scaling multiplier applied on top of
#'   base_size. Use for quick overall scaling. Default: 1.0.
#' @param row_padding Numeric vector of length 2. Padding around row content
#'   in mm as c(vertical, horizontal). If NULL, auto-calculated from base_size.
#'   Default: NULL.
#' @param ci_pch Integer. Point character for CI. 15=square, 16=circle,
#'   18=diamond. Default: 15.
#' @param ci_lwd Numeric. Line width for CI lines. If NULL, auto-calculated
#'   from base_size. Default: NULL.
#' @param ci_Theight Numeric. Height of T-bar ends on CI. If NULL, auto-calculated
#'   from base_size. Default: NULL.
#' @param ci_col Character. Color for CI lines and points. Default: "black".
#' @param header_fontsize Numeric. Font size for column headers.
#'   If NULL, auto-calculated as base_size * scale + 1. Default: NULL.
#' @param body_fontsize Numeric. Font size for body text.
#'   If NULL, auto-calculated as base_size * scale. Default: NULL.
#' @param footnote_fontsize Numeric. Font size for footnotes.
#'   If NULL, auto-calculated as base_size * scale - 1. Default: NULL.
#' @param footnote_col Character. Color for footnote text. Default: "darkcyan".
#' @param title_fontsize Numeric. Font size for title.
#'   If NULL, auto-calculated as base_size * scale + 4. Default: NULL.
#' @param cv_fontsize Numeric. Font size for CV annotation text.
#'   If NULL, auto-calculated as base_size * scale. Default: NULL.
#' @param cv_col Character. Color for CV annotation text. Default: "gray30".
#' @param refline_lwd Numeric. Reference line width. If NULL, auto-calculated.
#'   Default: NULL.
#' @param refline_lty Character. Reference line type. Default: "dashed".
#' @param refline_col Character. Reference line color. Default: "gray30".
#' @param vertline_lwd Numeric. Vertical line width. If NULL, auto-calculated.
#'   Default: NULL.
#' @param vertline_lty Character. Vertical line type. Default: "dashed".
#' @param vertline_col Character. Vertical line color. Default: "gray20".
#' @param arrow_type Character. Arrow type: "open" or "closed". Default: "closed".
#' @param arrow_col Character. Arrow color. Default: "black".
#' @param summary_fill Character. Fill color for summary diamonds. Default: "black".
#' @param summary_col Character. Border color for summary diamonds. Default: "black".
#'
#' @return A list of class "fs_forest_theme" containing all theme parameters.
#'
#' @details
#' The \code{base_size} parameter is the primary way to control plot size.
#' When you change \code{base_size}, the following are automatically scaled:
#' \itemize{
#'   \item All font sizes (body, header, footnote, CV, title)
#'   \item Row padding (vertical and horizontal)
#'   \item CI line width and T-bar height
#'   \item Reference and vertical line widths
#' }
#'
#' The scaling formula uses base_size = 10 as the reference point:
#' \itemize{
#'   \item base_size = 10: Default sizing
#'   \item base_size = 12: 20% larger
#'   \item base_size = 14: 40% larger
#'   \item base_size = 16: 60% larger
#' }
#'
#' You can override any individual parameter by specifying it explicitly.
#'
#' The theme does NOT set row background colors - those are determined
#' automatically by \code{plot_subgroup_results_forestplot()} based on
#' row types (ITT, reference, posthoc, etc.).
#'
#' @examples
#' \dontrun{
#' # Simple: just increase base_size for larger plot
#' large_theme <- create_forest_theme(base_size = 14)
#'
#' # Or use scale for quick adjustment
#' large_theme <- create_forest_theme(base_size = 10, scale = 1.4)
#'
#' # Fine-tune specific elements
#' custom_theme <- create_forest_theme(
#'   base_size = 14,
#'   cv_fontsize = 12,  # Override auto-calculated CV font size
#'   ci_lwd = 2.5       # Override auto-calculated CI line width
#' )
#'
#' # Use with plot_subgroup_results_forestplot
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   theme = large_theme
#' )
#'
#' render_forestplot(result)
#' }
#'
#' @seealso \code{\link{plot_subgroup_results_forestplot}}, \code{\link{render_forestplot}}
#'
#' @importFrom grid unit gpar
#' @export
create_forest_theme <- function(
    base_size = 10,
    scale = 1.0,
    row_padding = NULL,
    ci_pch = 15,
    ci_lwd = NULL,
    ci_Theight = NULL,
    ci_col = "black",
    header_fontsize = NULL,
    body_fontsize = NULL,
    footnote_fontsize = NULL,
    footnote_col = "darkcyan",
    title_fontsize = NULL,
    cv_fontsize = NULL,
    cv_col = "gray30",
    refline_lwd = NULL,
    refline_lty = "dashed",
    refline_col = "gray30",
    vertline_lwd = NULL,
    vertline_lty = "dashed",
    vertline_col = "gray20",
    arrow_type = "closed",
    arrow_col = "black",
    summary_fill = "black",
    summary_col = "black"
) {

  if (!requireNamespace("forestploter", quietly = TRUE)) {
    stop("Package 'forestploter' is required but not installed.", call. = FALSE)
  }

  # Calculate scaling factor relative to base_size = 10
  size_factor <- (base_size / 10) * scale

  # Auto-calculate font sizes if not provided
  if (is.null(body_fontsize)) body_fontsize <- base_size * scale
  if (is.null(header_fontsize)) header_fontsize <- base_size * scale + 1
  if (is.null(footnote_fontsize)) footnote_fontsize <- max(base_size * scale - 1, 8)
  if (is.null(title_fontsize)) title_fontsize <- base_size * scale + 4
  if (is.null(cv_fontsize)) cv_fontsize <- base_size * scale

 # Auto-calculate row padding if not provided
  # Base padding is c(4, 4) at base_size = 10
  if (is.null(row_padding)) {
    row_padding <- c(4 * size_factor, 4 * size_factor)
  }

  # Auto-calculate line widths if not provided
  # Base ci_lwd is 1.5 at base_size = 10
  if (is.null(ci_lwd)) ci_lwd <- 1.5 * size_factor
  if (is.null(ci_Theight)) ci_Theight <- 0.2 * size_factor
  if (is.null(refline_lwd)) refline_lwd <- 1.0 * size_factor
  if (is.null(vertline_lwd)) vertline_lwd <- 1.0 * size_factor

  # Create result object with parameters
  result <- list(
    base_size = base_size,
    scale = scale,
    size_factor = size_factor,
    row_padding = row_padding,
    body_fontsize = body_fontsize,
    header_fontsize = header_fontsize,
    footnote_fontsize = footnote_fontsize,
    footnote_col = footnote_col,
    title_fontsize = title_fontsize,
    cv_fontsize = cv_fontsize,
    cv_col = cv_col,
    ci_pch = ci_pch,
    ci_lwd = ci_lwd,
    ci_Theight = ci_Theight,
    ci_col = ci_col,
    refline_lwd = refline_lwd,
    refline_lty = refline_lty,
    refline_col = refline_col,
    vertline_lwd = vertline_lwd,
    vertline_lty = vertline_lty,
    vertline_col = vertline_col,
    arrow_type = arrow_type,
    arrow_col = arrow_col,
    summary_fill = summary_fill,
    summary_col = summary_col
  )

  class(result) <- c("fs_forest_theme", "list")

  return(result)
}


#' Print Method for ForestSearch Forest Theme
#'
#' @param x An fs_forest_theme object
#' @param ... Additional arguments (ignored)
#' @export
print.fs_forest_theme <- function(x, ...) {
  cat("ForestSearch Forest Plot Theme\n")
  cat("==============================\n")
  cat("Base size:       ", x$base_size, "\n")
  cat("Scale:           ", x$scale, "\n")
  cat("Size factor:     ", round(x$size_factor, 2), "x (relative to base_size=10)\n")
  cat("\nCalculated values:\n")
  cat("  Body font:     ", round(x$body_fontsize, 1), "\n")
  cat("  Header font:   ", round(x$header_fontsize, 1), "\n")
  cat("  CV font:       ", round(x$cv_fontsize, 1), "\n")
  cat("  Footnote font: ", round(x$footnote_fontsize, 1), "\n")
  cat("  Row padding:   ", paste(round(x$row_padding, 1), collapse = ", "), " mm\n")
  cat("  CI line width: ", round(x$ci_lwd, 2), "\n")
  cat("\nUse with: plot_subgroup_results_forestplot(..., theme = x)\n")
  invisible(x)
}


#' Render ForestSearch Forest Plot
#'
#' Renders a forest plot from \code{plot_subgroup_results_forestplot()}.
#'
#' @param x An fs_forestplot object from \code{plot_subgroup_results_forestplot()}.
#' @param newpage Logical. Call grid.newpage() before drawing. Default: TRUE.
#'
#' @return Invisibly returns the grob object.
#'
#' @details
#' To control plot sizing, create a custom theme using \code{create_forest_theme()}
#' and pass it to \code{plot_subgroup_results_forestplot()}:
#'
#' \code{my_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))}
#'
#' \code{result <- plot_subgroup_results_forestplot(..., theme = my_theme)}
#'
#' \code{render_forestplot(result)}
#'
#' @examples
#' \dontrun{
#' # For larger plot, use theme parameter in plot_subgroup_results_forestplot:
#' large_theme <- create_forest_theme(
#'   base_size = 14,
#'   row_padding = c(6, 4),
#'   cv_fontsize = 11
#' )
#'
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   theme = large_theme
#' )
#'
#' render_forestplot(result)
#' }
#'
#' @importFrom grid grid.newpage grid.draw
#' @export
render_forestplot <- function(x, newpage = TRUE) {

  if (!inherits(x, "fs_forestplot")) {
    stop("x must be an fs_forestplot object from plot_subgroup_results_forestplot()")
  }

  if (newpage) grid::grid.newpage()
  grid::grid.draw(x$plot)

  invisible(x$plot)
}


#' Save ForestSearch Forest Plot to File
#'
#' Saves a forest plot to a file (PDF, PNG, etc.) with explicit dimensions.
#'
#' @param x An fs_forestplot object.
#' @param filename Character. Output filename. Extension determines format.
#' @param width Numeric. Plot width in inches. Default: 12.
#' @param height Numeric. Plot height in inches. Default: 10.
#' @param dpi Numeric. Resolution for raster formats. Default: 300.
#' @param bg Character. Background color. Default: "white".
#'
#' @return Invisibly returns the filename.
#'
#' @examples
#' \dontrun{
#' # Create plot with custom theme
#' large_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))
#'
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   theme = large_theme
#' )
#'
#' # Save to file
#' save_forestplot(result, "forest_plot.pdf", width = 14, height = 12)
#' }
#'
#' @importFrom grDevices pdf png svg jpeg dev.off
#' @export
save_forestplot <- function(
    x,
    filename,
    width = 12,
    height = 10,
    dpi = 300,
    bg = "white"
) {

  if (!inherits(x, "fs_forestplot")) {
    stop("x must be an fs_forestplot object")
  }

  # Determine format from extension
  ext <- tolower(tools::file_ext(filename))

  # Open appropriate device
  if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = height, bg = bg)
  } else if (ext == "png") {
    grDevices::png(filename, width = width, height = height, units = "in",
                   res = dpi, bg = bg)
  } else if (ext == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(filename, width = width, height = height, bg = bg)
    } else {
      grDevices::svg(filename, width = width, height = height, bg = bg)
    }
  } else if (ext %in% c("jpg", "jpeg")) {
    grDevices::jpeg(filename, width = width, height = height, units = "in",
                    res = dpi, bg = bg, quality = 95)
  } else {
    stop("Unsupported file format: ", ext, ". Use .pdf, .png, .svg, or .jpg")
  }

  # Draw the plot
  grid::grid.newpage()
  grid::grid.draw(x$plot)

  # Close device
  grDevices::dev.off()

  message("Saved forest plot to: ", filename)
  invisible(filename)
}


#' Print Method for ForestSearch Forest Plot
#'
#' @param x An fs_forestplot object
#' @param ... Additional arguments (ignored)
#' @export
print.fs_forestplot <- function(x, ...) {
  cat("ForestSearch Subgroup Results Forest Plot\n")
  cat("=========================================\n")
  cat("Number of rows:", nrow(x$data), "\n")
  cat("Row types:", paste(unique(x$row_types), collapse = ", "), "\n")

  if (length(x$cv_metrics) > 0) {
    cat("\nCross-validation metrics:\n")
    for (i in seq_along(x$cv_metrics)) {
      cat("  -", x$cv_metrics[[i]], "\n")
    }
  }

  cat("\nTo display: render_forestplot(x) or plot(x)\n")
  cat("To save:    save_forestplot(x, 'plot.pdf', width = 14, height = 12)\n")
  cat("\nFor larger plot or bigger CV text, recreate with custom theme:\n")
  cat("  theme <- create_forest_theme(base_size = 14, cv_fontsize = 11)\n")
  cat("  result <- plot_subgroup_results_forestplot(..., theme = theme)\n")

  invisible(x)
}


#' Plot Method for ForestSearch Forest Plot
#'
#' @param x An fs_forestplot object
#' @param ... Additional arguments (ignored)
#' @export
plot.fs_forestplot <- function(x, ...) {
  render_forestplot(x)
  invisible(x)
}
