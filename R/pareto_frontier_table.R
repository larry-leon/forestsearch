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
#' \strong{Optional confidence intervals.}  When \code{add_cis = TRUE},
#' three CI columns are added to the table.  All three are computed
#' post-hoc by \code{\link{compute_frontier_cis}} and do not modify the
#' \code{fs} object.
#' \itemize{
#'   \item \code{Naive 95\% CI} -- full-sample Wald CI from a Cox or
#'     GLM refit on each frontier member's data.  Ignores
#'     subgroup-search selection; anti-conservative by construction.
#'   \item \code{Split ~ 95\% CI} -- subsample-derived approximation,
#'     computed from the empirical SD of averaged half-sample effects
#'     across \code{n_splits} random 50/50 splits.
#'   \item \code{FSBC ~ 95\% CI} -- bias-corrected interval following
#'     the bootstrap algorithm of \cite{Leon2024fs} (eq 7, 9) but
#'     treating the selected subgroup as fixed across splits.  The
#'     cell shows \code{est (lcl, ucl)} where \code{est} is the
#'     bias-corrected effect estimate \eqn{2\hat\beta - \bar{\tilde\beta}}.
#'     A trailing dagger marks rows where the bias-corrected variance
#'     was non-positive (CI computed with \eqn{V} clipped to 0).
#' }
#' See \code{\link{compute_frontier_cis}} for the algebra.
#'
#' @param add_cis Logical.  If \code{TRUE}, compute and display naive
#'   and split-derived 95\% CIs for each frontier member.  Default
#'   \code{FALSE} (backward-compatible).
#' @param n_splits Integer.  Number of 50/50 splits per frontier member
#'   for the split-derived CI.  Used only when \code{add_cis = TRUE}.
#'   Default \code{1000L}.
#' @param ci_seed Integer or \code{NULL}.  Seed for split-derived CI
#'   reproducibility.  Default \code{NULL}.
#'
#' @examples
#' \dontrun{
#' # Survival example
#' data(gbsg, package = "survival")
#' fs <- forestsearch(gbsg, ...)
#' pareto_frontier_table(fs)                       # gt table
#' pareto_frontier_table(fs, format = "data.table") # raw frontier
#' pareto_frontier_table(fs, add_cis = TRUE)       # with CIs
#'
#' # Programmatic use: how many alternatives are on the frontier?
#' nrow(pareto_frontier_table(fs, format = "data.table"))
#' }
#'
#' @seealso \code{\link{compute_pareto_frontier}},
#'   \code{\link{compute_frontier_cis}},
#'   \code{\link{sort_subgroups}}, \code{\link{forestsearch}}.
#' @importFrom data.table is.data.table copy setnames data.table
#' @export
pareto_frontier_table <- function(fs,
                                  format = c("gt", "data.table"),
                                  digits_effect = 3L,
                                  digits_pcons  = 3L,
                                  include_factor_columns = TRUE,
                                  highlight_selected = TRUE,
                                  add_cis  = FALSE,
                                  n_splits = 1000L,
                                  ci_seed  = NULL) {

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

  # --- 4b. Optional confidence intervals -------------------------------
  ci_cols <- character(0)
  if (isTRUE(add_cis)) {
    ci_dt <- tryCatch(
      compute_frontier_cis(fs, n_splits = n_splits, seed = ci_seed),
      error = function(e) {
        warning("compute_frontier_cis() failed: ", conditionMessage(e),
                call. = FALSE)
        NULL
      })
    if (!is.null(ci_dt) && data.table::is.data.table(ci_dt) &&
        nrow(ci_dt) > 0L) {
      # Defensive: not every ci_dt build will contain the FSBC columns
      # (e.g. if a future release strips them).  Only request what exists.
      ci_cols_want <- intersect(
        c("m", "naive_lcl", "naive_ucl",
          "split_lcl", "split_ucl",
          "fsbc_estimate", "fsbc_lcl", "fsbc_ucl", "fsbc_var_pos"),
        names(ci_dt))
      ft <- merge(
        ft,
        ci_dt[, ci_cols_want, with = FALSE],
        by = "m", all.x = TRUE, sort = FALSE
      )
      # Pre-format the CI strings on the natural scale.  Width is fixed
      # at digits_effect for both bounds so columns align visually.
      fmt <- function(x) {
        if (is.na(x)) "NA" else
          formatC(x, format = "f", digits = digits_effect)
      }
      ft[["Naive 95% CI"]] <- vapply(seq_len(nrow(ft)), function(k) {
        if (is.na(ft$naive_lcl[k]) || is.na(ft$naive_ucl[k])) "NA"
        else sprintf("(%s, %s)",
                     fmt(ft$naive_lcl[k]), fmt(ft$naive_ucl[k]))
      }, character(1))
      ft[["Split ~ 95% CI"]] <- vapply(seq_len(nrow(ft)), function(k) {
        if (is.na(ft$split_lcl[k]) || is.na(ft$split_ucl[k])) "NA"
        else sprintf("(%s, %s)",
                     fmt(ft$split_lcl[k]), fmt(ft$split_ucl[k]))
      }, character(1))
      ci_cols <- c("Naive 95% CI", "Split ~ 95% CI")

      # FSBC-mimic column: show "est (lcl, ucl)" so the bias-corrected
      # estimate appears alongside its interval.  Tilde signals
      # "approximation; not full FSBC".  Flag degenerate-variance rows
      # with a trailing dagger.
      if (all(c("fsbc_estimate", "fsbc_lcl",
                "fsbc_ucl", "fsbc_var_pos") %in% names(ft))) {
        ft[["FSBC ~ 95% CI"]] <- vapply(seq_len(nrow(ft)), function(k) {
          if (is.na(ft$fsbc_estimate[k]) ||
              is.na(ft$fsbc_lcl[k]) ||
              is.na(ft$fsbc_ucl[k])) {
            return("NA")
          }
          flag <- if (isFALSE(ft$fsbc_var_pos[k])) " \u2020" else ""
          sprintf("%s (%s, %s)%s",
                  fmt(ft$fsbc_estimate[k]),
                  fmt(ft$fsbc_lcl[k]),
                  fmt(ft$fsbc_ucl[k]),
                  flag)
        }, character(1))
        ci_cols <- c(ci_cols, "FSBC ~ 95% CI")
      }
    }
  }

  # --- 5. Choose display columns ---------------------------------------
  m_cols <- grep("^M\\.", names(ft), value = TRUE)
  if (!include_factor_columns) m_cols <- character(0)

  metric_cols <- intersect(c("N", "E", effect_measure, "Pcons", "K"),
                           names(ft))
  display_cols <- c("is_selected", m_cols, metric_cols, ci_cols)
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

  # --- CI spanner: group the three CI columns under "95% CIs" --------------
  # Use only the CI columns that exist in ft_display; relabel each to a
  # short approach tag and group them under a single spanner header.
  ci_label_map <- c(
    "Naive 95% CI"   = "Naive",
    "Split ~ 95% CI" = "Split ~",
    "FSBC ~ 95% CI"  = "FSBC ~"
  )
  ci_cols_present <- intersect(names(ci_label_map), names(ft_display))
  if (length(ci_cols_present) > 0L) {
    short_labels <- as.list(ci_label_map[ci_cols_present])
    names(short_labels) <- ci_cols_present
    tbl <- gt::cols_label(tbl, .list = short_labels)
    tbl <- gt::cols_align(tbl, align = "center",
                          columns = ci_cols_present)
    tbl <- gt::tab_spanner(tbl,
                           label   = "95% CIs",
                           columns = ci_cols_present)
    # Dagger footnote anchored to FSBC subcolumn (only if present)
    if ("FSBC ~ 95% CI" %in% ci_cols_present) {
      tbl <- gt::tab_footnote(
        tbl,
        footnote  = paste0(
          "\u2020 marks rows where the bias-corrected variance was ",
          "non-positive (Monte-Carlo correction in eq. 9 of Le\u00f3n ",
          "et al. 2024 exceeded the IJ term).  The CI is computed with ",
          "V clipped to 0 and should be interpreted as degenerate."
        ),
        locations = gt::cells_column_labels(columns = "FSBC ~ 95% CI")
      )
    }
  }

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
