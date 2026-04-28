#' Format Pareto Frontier of Candidate Subgroups
#'
#' Renders the post-hoc Pareto frontier on (effect, N) -- both maximized
#' -- as a formatted \pkg{gt} table or returns it as a
#' \code{data.table} for programmatic use.  Works uniformly across
#' survival (Cox PH) and GLM outcome types: the effect-column label and
#' scale handling are derived from the \code{forestsearch} object's
#' \code{effect_measure}.
#'
#' The frontier is a diagnostic: it lists candidate subgroups that are
#' \emph{not dominated} on (effect, N) simultaneously.  It is computed
#' inside \code{\link{sg_consistency_out}} (see
#' \code{\link{compute_pareto_frontier}}) and attached to the result
#' object as \code{fs$grp.consistency$out_sg$pareto_frontier}.  It is
#' \strong{not} used for selection -- the selected subgroup is chosen
#' by \code{sg_focus} and may or may not appear on the frontier.
#'
#' @param fs A \code{forestsearch} object returned by
#'   \code{\link{forestsearch}}.
#' @param format Character.  Either \code{"gt"} (default) for a
#'   formatted \pkg{gt} table or \code{"data.table"} for a plain
#'   \code{data.table} with effect on the natural scale and an
#'   \code{is_selected} flag.
#' @param digits_effect Integer.  Decimal places for the effect column.
#'   Default \code{3}.
#' @param digits_pcons Integer.  Decimal places for \code{Pcons}.
#'   Default \code{3}.
#' @param include_factor_columns Logical.  If \code{TRUE} (default),
#'   include the \code{M.<factor>} columns identifying each subgroup's
#'   defining cuts.  If \code{FALSE}, show only summary metrics.
#' @param highlight_selected Logical.  If \code{TRUE} (default), the
#'   selected subgroup row is highlighted in the \pkg{gt} table.
#'   Ignored when \code{format = "data.table"}.
#'
#' @return Depending on \code{format}:
#'   \describe{
#'     \item{\code{"gt"}}{A \code{gt_tbl} object.  Returns \code{NULL}
#'       invisibly (with a message) if the frontier is unavailable or
#'       empty.}
#'     \item{\code{"data.table"}}{A \code{data.table} of the frontier
#'       with effect on the natural scale, an \code{is_selected}
#'       logical column, and effect column renamed to the
#'       \code{effect_measure} label (e.g., \code{"HR"}, \code{"OR"}).
#'       Returns an empty \code{data.table} if the frontier is
#'       unavailable.}
#'   }
#'
#' @details
#' \strong{Scale handling.}  For ratio measures stored on the log
#' scale internally (\code{"OR"}, \code{"RR"}, \code{"IRR"}), the
#' effect column is exponentiated for display.  For \code{"HR"} (Cox
#' PH, natural scale) and identity-scale measures (\code{"RD"},
#' \code{"IRD"}, \code{"MD"}), values pass through unchanged.
#'
#' \strong{Selected-row identification.}  The selected subgroup is
#' identified by matching the original-table row index \code{m}
#' against the top row of \code{fs$grp.consistency$out_sg$result}.
#' This is robust to sorting and outcome-type differences.
#'
#' \strong{For \code{sg_focus = "hrMinSG"}}, the selected subgroup may
#' not appear on the frontier -- that focus deliberately prefers small
#' subgroups, which are typically N-dominated.
#'
#' @examples
#' \dontrun{
#' # Survival example
#' data(gbsg, package = "survival")
#' fs <- forestsearch(gbsg, ...)
#' pareto_frontier_table(fs)                       # gt table
#' pareto_frontier_table(fs, format = "data.table") # raw frontier
#'
#' # Programmatic use: how many alternatives are on the frontier?
#' nrow(pareto_frontier_table(fs, format = "data.table"))
#' }
#'
#' @seealso \code{\link{compute_pareto_frontier}},
#'   \code{\link{sort_subgroups}}, \code{\link{forestsearch}}.
#' @importFrom data.table is.data.table copy setnames data.table
#' @export
pareto_frontier_table <- function(fs,
                                  format = c("gt", "data.table"),
                                  digits_effect = 3L,
                                  digits_pcons  = 3L,
                                  include_factor_columns = TRUE,
                                  highlight_selected = TRUE) {

  format <- match.arg(format)

  # --- 1. Locate the frontier on the fs object -------------------------
  out_sg   <- tryCatch(fs$grp.consistency$out_sg, error = function(e) NULL)
  frontier <- tryCatch(out_sg$pareto_frontier,    error = function(e) NULL)

  if (is.null(out_sg) || is.null(frontier)) {
    message("No Pareto frontier available on this forestsearch object.")
    return(.empty_pareto_return(format))
  }
  if (!data.table::is.data.table(frontier) || nrow(frontier) == 0L) {
    message("Pareto frontier is empty (no subgroups passed consistency).")
    return(.empty_pareto_return(format))
  }

  # --- 2. Resolve effect-measure label and scale -----------------------
  effect_measure <- if (!is.null(fs$effect_measure)) fs$effect_measure
                    else "HR"
  effect_log_scale <- isTRUE(effect_measure %in% c("OR", "RR", "IRR"))

  # --- 3. Identify the selected row ------------------------------------
  selected_m <- tryCatch(
    as.integer(out_sg$result[1, ]$m),
    error = function(e) NA_integer_
  )

  ft <- data.table::copy(frontier)
  ft_m <- if ("m" %in% names(ft)) as.integer(ft[["m"]]) else NA_integer_
  is_selected_vec <- !is.na(selected_m) & !is.na(ft_m) & ft_m == selected_m
  ft[["is_selected"]] <- is_selected_vec

  # --- 4. Convert effect to natural scale, rename column ---------------
  if (effect_log_scale && "hr" %in% names(ft)) {
    ft[["hr"]] <- exp(as.numeric(ft[["hr"]]))
  }
  if ("hr" %in% names(ft)) {
    data.table::setnames(ft, "hr", effect_measure)
  }

  # --- 5. Choose display columns ---------------------------------------
  m_cols <- grep("^M\\.", names(ft), value = TRUE)
  if (!include_factor_columns) m_cols <- character(0)

  metric_cols <- intersect(c("N", "E", effect_measure, "Pcons", "K"),
                           names(ft))
  display_cols <- c("is_selected", m_cols, metric_cols)
  display_cols <- intersect(display_cols, names(ft))

  ft <- ft[, display_cols, with = FALSE]

  # --- 6. Branch on format ---------------------------------------------
  if (format == "data.table") return(ft[])

  .render_pareto_gt(
    ft                  = ft,
    is_selected_vec     = is_selected_vec,
    effect_measure      = effect_measure,
    effect_log_scale    = effect_log_scale,
    digits_effect       = digits_effect,
    digits_pcons        = digits_pcons,
    sg_focus            = .resolve_sg_focus(fs),
    highlight_selected  = highlight_selected
  )
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.empty_pareto_return <- function(format) {
  if (format == "data.table") {
    return(data.table::data.table())
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.resolve_sg_focus <- function(fs) {
  sgf <- tryCatch(fs$grp.consistency$sg_focus, error = function(e) NULL)
  if (is.null(sgf)) sgf <- tryCatch(fs$args_call_all$sg_focus,
                                    error = function(e) NULL)
  if (is.null(sgf)) "(unknown)" else as.character(sgf)
}

#' @keywords internal
#' @noRd
.render_pareto_gt <- function(ft, is_selected_vec,
                              effect_measure, effect_log_scale,
                              digits_effect, digits_pcons, sg_focus,
                              highlight_selected) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    warning("Package 'gt' not installed; returning data.table instead.",
            call. = FALSE)
    return(ft[])
  }

  # Build "Selected" column with star marker; drop is_selected from body
  selected_marker <- ifelse(is_selected_vec, "\u2605", "")
  ft_display <- data.table::copy(ft)
  if ("is_selected" %in% names(ft_display)) {
    ft_display[["is_selected"]] <- NULL
  }
  ft_display <- cbind(Selected = selected_marker,
                      as.data.frame(ft_display))

  scale_note <- if (effect_log_scale)
    sprintf("%s shown on natural scale (exponentiated from log-scale storage)",
            effect_measure)
  else
    sprintf("%s shown on natural scale", effect_measure)

  subtitle <- sprintf(
    "%d non-dominated candidate(s) | sg_focus = '%s' | %s | diagnostic only",
    nrow(ft_display), sg_focus, scale_note
  )

  tbl <- gt::gt(ft_display) |>
    gt::tab_header(
      title    = sprintf("Pareto frontier on (%s, N)", effect_measure),
      subtitle = subtitle
    ) |>
    gt::cols_label(Selected = "")

  # Numeric formatting (only for columns that exist)
  if (effect_measure %in% names(ft_display)) {
    tbl <- gt::fmt_number(tbl, columns = effect_measure,
                          decimals = digits_effect)
  }
  if ("Pcons" %in% names(ft_display)) {
    tbl <- gt::fmt_number(tbl, columns = "Pcons",
                          decimals = digits_pcons)
  }
  int_cols <- intersect(c("N", "E", "K"), names(ft_display))
  if (length(int_cols) > 0L) {
    tbl <- gt::fmt_number(tbl, columns = int_cols, decimals = 0)
  }

  # Right-align the marker column
  tbl <- gt::cols_align(tbl, align = "center", columns = "Selected")

  # Highlight the selected row
  if (isTRUE(highlight_selected) && any(is_selected_vec)) {
    sel_rows <- which(is_selected_vec)
    tbl <- gt::tab_style(
      tbl,
      style = list(
        gt::cell_fill(color = "#EEF6FB"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(rows = sel_rows)
    )
  }

  # Footnote anchored to the Selected column header
  tbl <- gt::tab_footnote(
    tbl,
    footnote  = sprintf(
      "%s marks the subgroup chosen by sg_focus = '%s'.",
      "\u2605", sg_focus
    ),
    locations = gt::cells_column_labels(columns = "Selected")
  )

  tbl
}
