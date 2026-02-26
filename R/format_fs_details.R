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
#' @param max_width Integer. Maximum character width per line before wrapping.
#'   Long lines are wrapped at comma or space boundaries with a 4-space
#'   continuation indent. Default: 48 (suitable for half-slide columns at
#'   scriptsize).
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
                              col_widths = c(0.48, 0.52),
                              max_width = 48) {

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

  # --- Wrap long lines at comma or space boundaries ---
  wrap_line <- function(line, width, indent = 4) {
    if (nchar(line) <= width) return(line)

    result <- character(0)
    remaining <- line
    # Detect leading whitespace of original line for continuation
    lead <- nchar(sub("\\S.*", "", remaining))
    cont_prefix <- paste(rep(" ", lead + indent), collapse = "")

    first <- TRUE
    while (nchar(remaining) > width) {
      chunk <- substr(remaining, 1, width)

      # Prefer breaking at ", " (comma-space) within the chunk
      comma_pos <- gregexpr(",\\s?", chunk)[[1]]
      comma_pos <- comma_pos[comma_pos > 0]

      if (length(comma_pos) > 0) {
        # Break after the last comma+space within width
        last_comma <- max(comma_pos)
        match_len <- attr(gregexpr(",\\s?", chunk)[[1]], "match.length")
        idx <- which(comma_pos == last_comma)
        break_at <- last_comma + match_len[idx] - 1
      } else {
        # Fall back to last space
        space_pos <- gregexpr(" ", chunk)[[1]]
        space_pos <- space_pos[space_pos > 0]
        if (length(space_pos) > 0) {
          break_at <- max(space_pos)
        } else {
          # No good break point — hard break at width
          break_at <- width
        }
      }

      result <- c(result, substr(remaining, 1, break_at))
      remaining <- paste0(cont_prefix, trimws(substr(remaining, break_at + 1,
                                                      nchar(remaining)),
                                               which = "left"))
    }
    if (nchar(remaining) > 0) {
      result <- c(result, remaining)
    }
    result
  }

  # --- Escape LaTeX special characters for \texttt{} ---
  # Order: braces and specials first, then space → ~ last.
  # Note: backslashes never appear in R console output, so not escaped.
  escape_latex <- function(x) {
    x <- gsub("{", "\\{", x, fixed = TRUE)
    x <- gsub("}", "\\}", x, fixed = TRUE)
    x <- gsub("#", "\\#", x, fixed = TRUE)
    x <- gsub("$", "\\$", x, fixed = TRUE)
    x <- gsub("%", "\\%", x, fixed = TRUE)
    x <- gsub("&", "\\&", x, fixed = TRUE)
    x <- gsub("_", "\\_", x, fixed = TRUE)
    x <- gsub("^", "\\textasciicircum{}", x, fixed = TRUE)
    x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
    # Last: replace EVERY space with ~ (non-breaking space).
    # In \texttt{}, regular spaces collapse; ~ preserves exact
    # character positions for column-aligned tabular output.
    x <- gsub(" ", "~", x, fixed = TRUE)
    x
  }

  emit_column <- function(lines_vec) {
    # Detect table blocks: lines between === or --- separators should
    # NOT be wrapped, as they have precise column alignment.
    # Also skip wrapping for any line with 3+ consecutive spaces
    # (indicates columnar alignment, e.g. table headers/data rows).
    in_table <- FALSE

    for (ln in lines_vec) {
      if (ln == "") {
        cat("\\vspace{2pt}\n")
        next
      }

      # Detect separator lines (=== or --- patterns)
      is_separator <- grepl("^[=\\-]{10,}", ln)

      if (is_separator) {
        in_table <- !in_table
        cat("\\texttt{", escape_latex(ln), "}\\\\[-2pt]\n", sep = "")
        next
      }

      # Detect columnar alignment: 3+ consecutive spaces mid-line
      is_columnar <- grepl("\\S {3,}\\S", ln)

      if (in_table || is_columnar) {
        # Table block or columnar line: emit as-is, no wrapping
        cat("\\texttt{", escape_latex(ln), "}\\\\[-2pt]\n", sep = "")
      } else {
        # Non-tabular line: wrap if too long
        wrapped <- wrap_line(ln, max_width)
        for (wl in wrapped) {
          cat("\\texttt{", escape_latex(wl), "}\\\\[-2pt]\n", sep = "")
        }
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
