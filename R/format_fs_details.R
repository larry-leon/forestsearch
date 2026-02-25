#' Format ForestSearch Details Output for Beamer Two-Column Display
#'
#' Captures \code{forestsearch(details = TRUE)} console output and splits it
#' into two columns for readable beamer slides. Left column shows variable
#' selection (GRF, LASSO, candidate factors); right column shows subgroup
#' search, consistency evaluation, and results.
#'
#' @param fs_output Character vector of captured output lines from
#'   \code{capture.output(forestsearch(..., details = TRUE))}.
#' @param split_after Character string (regex). The output is split after the
#'   block matching this pattern. Default: \code{"Candidate factors"} which
#'   splits after the candidate factor list.
#' @param fontsize Character. LaTeX font size for the output text.
#'   One of \code{"tiny"}, \code{"scriptsize"}, \code{"footnotesize"},
#'   \code{"small"}. Default: \code{"scriptsize"}.
#' @param col_widths Numeric vector of length 2. Column widths as fractions
#'   of \code{\\textwidth}. Default: \code{c(0.48, 0.52)}.
#'
#' @return Invisibly returns a list with \code{left} and \code{right}
#'   character vectors. Side effect: emits LaTeX via \code{cat()} for use
#'   in a chunk with \code{results='asis'}.
#'
#' @section Quarto Setup:
#' No special LaTeX packages required. Works in any beamer frame
#' without the \code{fragile} option.
#'
#' @section Usage:
#' In a Quarto beamer chunk with \code{results='asis'} and
#' \code{echo=FALSE}, first capture the forestsearch output with
#' \code{capture.output()}, then call \code{format_fs_details(fs_output)}.
#'
#' @examples
#' \dontrun{
#' fs_output <- capture.output(
#'   fs_res <- forestsearch(
#'     df = dat, outcome.name = "time", event.name = "status",
#'     treat.name = "trt", id.name = "id",
#'     confounders.name = confs, details = TRUE
#'   )
#' )
#' # Default: split after candidate factors
#' format_fs_details(fs_output, fontsize = "tiny")
#'
#' # Custom split: after filtering summary
#' format_fs_details(fs_output, split_after = "Found.*candidate")
#' }
#'
#' @export
format_fs_details <- function(fs_output,
                              split_after = "Candidate factors",
                              fontsize = "scriptsize",
                              col_widths = c(0.48, 0.52)) {

  stopifnot(is.character(fs_output))

  # --- Clean lines ---
  lines <- trimws(fs_output, which = "right")
  while (length(lines) > 0 && lines[1] == "") lines <- lines[-1]
  while (length(lines) > 0 && lines[length(lines)] == "") {
    lines <- lines[-length(lines)]
  }

  if (length(lines) == 0) {
    cat("% format_fs_details: no output to display\n")
    return(invisible(list(left = character(0), right = character(0))))
  }

  # --- Find split point ---
  anchor <- grep(split_after, lines)

  if (length(anchor) == 0) {
    # No match: split roughly in half
    warning("Pattern '", split_after, "' not found. Splitting at midpoint.")
    split_idx <- ceiling(length(lines) / 2)
  } else {
    split_idx <- anchor[1]

    # For "Candidate factors", extend past the printed vector (lines with "[")
    if (grepl("Candidate factors", split_after, fixed = TRUE)) {
      for (i in seq(split_idx + 1, min(split_idx + 20, length(lines)))) {
        if (grepl("^\\s*\\[", lines[i])) {
          split_idx <- i
        } else {
          break
        }
      }
    }
  }

  left_lines <- lines[seq_len(split_idx)]
  right_lines <- if (split_idx < length(lines)) {
    lines[seq(split_idx + 1, length(lines))]
  } else {
    character(0)
  }

  # Trim leading/trailing blank lines from each column
  trim_blank <- function(x) {
    while (length(x) > 0 && x[1] == "") x <- x[-1]
    while (length(x) > 0 && x[length(x)] == "") x <- x[-length(x)]
    x
  }
  left_lines <- trim_blank(left_lines)
  right_lines <- trim_blank(right_lines)

  # --- Validate fontsize ---
  valid_sizes <- c("tiny", "scriptsize", "footnotesize", "small", "normalsize")
  if (!fontsize %in% valid_sizes) {
    warning("Invalid fontsize '", fontsize, "'. Using 'scriptsize'.")
    fontsize <- "scriptsize"
  }

  # --- Escape LaTeX special characters for \texttt{} ---
  escape_latex <- function(x) {
    x <- gsub("\\", "\\textbackslash ", x, fixed = TRUE)
    x <- gsub("{", "\\{", x, fixed = TRUE)
    x <- gsub("}", "\\}", x, fixed = TRUE)
    x <- gsub("#", "\\#", x, fixed = TRUE)
    x <- gsub("$", "\\$", x, fixed = TRUE)
    x <- gsub("%", "\\%", x, fixed = TRUE)
    x <- gsub("&", "\\&", x, fixed = TRUE)
    x <- gsub("_", "\\_", x, fixed = TRUE)
    x <- gsub("~", "\\textasciitilde ", x, fixed = TRUE)
    x <- gsub("^", "\\textasciicircum ", x, fixed = TRUE)
    # Preserve spaces with \phantom{x}-free approach
    x <- gsub("  ", "~~ ", x, fixed = TRUE)
    x
  }

  emit_column <- function(lines_vec) {
    for (ln in lines_vec) {
      if (ln == "") {
        cat("\\vspace{2pt}\n")
      } else {
        cat("\\texttt{", escape_latex(ln), "}\\\\[-2pt]\n", sep = "")
      }
    }
  }

  # --- Emit LaTeX: texttt lines inside beamer columns ---
  # No fragile frame or fancyvrb needed
  cat("\\begin{columns}[T]\n")

  # Left column: variable selection
  cat("\\begin{column}{", col_widths[1], "\\textwidth}\n", sep = "")
  cat("{\\", fontsize, "\\raggedright\n", sep = "")
  emit_column(left_lines)
  cat("}\n")
  cat("\\end{column}\n")

  # Right column: search + consistency
  cat("\\begin{column}{", col_widths[2], "\\textwidth}\n", sep = "")
  cat("{\\", fontsize, "\\raggedright\n", sep = "")
  emit_column(right_lines)
  cat("}\n")
  cat("\\end{column}\n")

  cat("\\end{columns}\n")

  invisible(list(left = left_lines, right = right_lines))
}
