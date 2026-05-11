# Internal NULL-coalesce; declared locally so this file doesn't depend
# on a package-level export.  (data.table has no equivalent; rlang has
# one but we don't import rlang.)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

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
#' @param digits_ci Integer.  Decimal places used for the three CI
#'   columns (Naive, Split, FSBC) when \code{ci_table} is supplied.
#'   Default \code{2}.  CIs typically read more cleanly at 2 dp; raise
#'   to 3 if you need tighter precision.
#' @param include_dominated Logical.  If \code{FALSE} (default), the
#'   table shows only the Pareto-non-dominated subgroups (the
#'   frontier).  If \code{TRUE}, the table shows \strong{all} passing
#'   candidates -- frontier members AND dominated candidates -- with
#'   additional columns indicating frontier membership and band
#'   eligibility (when the selection rule uses one).  This is the
#'   "rich post-consistency view" that complements the in-flight
#'   summary printed by \code{forestsearch(show_candidate_summary = TRUE)}.
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
#' \strong{Optional confidence intervals.}  When \code{ci_table} is
#' supplied (the data.table returned by
#' \code{\link{compute_frontier_cis}}), three CI columns are added to
#' the table.  All three CIs are computed by
#' \code{\link{compute_frontier_cis}}; this function only displays
#' what is passed in.  Pass the same \code{ci_table} to
#' \code{\link{plot_pareto_frontier}} for consistent display across
#' table and plot.
#' \itemize{
#'   \item \code{Naive 95\% CI} -- full-sample Wald CI from a Cox or
#'     GLM refit on each frontier member's data.  Ignores
#'     subgroup-search selection; anti-conservative by construction.
#'   \item \code{Split ~ 95\% CI} -- subsample-derived approximation,
#'     computed from the empirical SD of averaged half-sample effects
#'     across the splits used by \code{compute_frontier_cis}.
#'   \item \code{FSBC ~ 95\% CI} -- bias-corrected interval following
#'     the bootstrap algorithm of \cite{Leon2024fs} (eq 7) but
#'     treating the selected subgroup as fixed across half-jackknife
#'     replicates.  The cell shows \code{est (lcl, ucl)} where
#'     \code{est} is the bias-corrected effect estimate
#'     \eqn{2\hat\beta - \overline{\hat\beta^{(h)}}}.
#' }
#' See \code{\link{compute_frontier_cis}} for the algebra.
#'
#' @param ci_table Optional data.table of frontier CIs from
#'   \code{\link{compute_frontier_cis}}.  When supplied, three CI
#'   columns are merged into the table; when \code{NULL} (default),
#'   no CI columns are shown.  This function does not run any
#'   computation on \code{ci_table} -- it simply displays the bounds
#'   already computed.
#'
#' @examples
#' \dontrun{
#' # Survival example
#' data(gbsg, package = "survival")
#' fs <- forestsearch(gbsg, ...)
#'
#' # Basic table (no CIs)
#' pareto_frontier_table(fs)
#' pareto_frontier_table(fs, format = "data.table")
#'
#' # With CIs (compute once, display anywhere)
#' ci_tab <- compute_frontier_cis(fs, n_splits = 1000, seed = 1)
#' pareto_frontier_table(fs, ci_table = ci_tab)
#' plot_pareto_frontier(fs,   ci_table = ci_tab)  # same intervals
#' }
#'
#' @seealso \code{\link{compute_pareto_frontier}},
#'   \code{\link{compute_frontier_cis}},
#'   \code{\link{plot_pareto_frontier}},
#'   \code{\link{sort_subgroups}}, \code{\link{forestsearch}}.
#' @importFrom data.table is.data.table copy setnames data.table
#' @export
pareto_frontier_table <- function(fs,
                                  format = c("gt", "data.table"),
                                  digits_effect = 3L,
                                  digits_pcons  = 3L,
                                  digits_ci     = 2L,
                                  include_dominated = FALSE,
                                  include_factor_columns = TRUE,
                                  highlight_selected = TRUE,
                                  ci_table = NULL) {

  format <- match.arg(format)

  # --- 1. Locate the frontier (or full passing set) on the fs object ---
  out_sg <- tryCatch(fs$grp.consistency$out_sg, error = function(e) NULL)
  if (is.null(out_sg)) {
    message("No subgroup-consistency output available on this forestsearch object.")
    return(.empty_pareto_return(format))
  }

  # Source table: full passing set if include_dominated, else just the frontier
  if (isTRUE(include_dominated)) {
    src <- tryCatch(out_sg$result, error = function(e) NULL)
    src_label <- "passing candidates"
  } else {
    src <- tryCatch(out_sg$pareto_frontier, error = function(e) NULL)
    src_label <- "Pareto frontier"
  }

  if (is.null(src)) {
    message(sprintf("No %s available on this forestsearch object.", src_label))
    return(.empty_pareto_return(format))
  }
  if (!data.table::is.data.table(src) || nrow(src) == 0L) {
    message(sprintf("%s is empty (no subgroups passed consistency).",
                    tools::toTitleCase(src_label)))
    return(.empty_pareto_return(format))
  }

  # Frontier-membership lookup (only meaningful when include_dominated = TRUE)
  frontier <- tryCatch(out_sg$pareto_frontier, error = function(e) NULL)
  frontier_m <- if (!is.null(frontier) && data.table::is.data.table(frontier) &&
                    nrow(frontier) > 0L && "m" %in% names(frontier)) {
    as.integer(frontier$m)
  } else {
    integer(0)
  }

  # Selection-rule + band parameters (for in_band flag in include_dominated mode)
  args_call <- fs$args_call_all %||% list()
  sg_focus_v       <- args_call$sg_focus %||% "hr"
  selection_rule_v <- args_call$selection_rule %||% "neighborhood"
  effect_nbhd_v    <- args_call$effect_neighborhood %||% 0.10
  band_used <- sg_focus_v %in% c("hrMaxSG", "hrMinSG") &&
               selection_rule_v %in% c("neighborhood", "both")

  # --- 2. Resolve effect-measure label and scale -----------------------
  effect_measure <- if (!is.null(fs$effect_measure)) fs$effect_measure
                    else "HR"
  effect_log_scale <- isTRUE(effect_measure %in% c("OR", "RR", "IRR"))

  # --- 3. Identify the selected row ------------------------------------
  selected_m <- tryCatch(
    as.integer(out_sg$result[1, ]$m),
    error = function(e) NA_integer_
  )

  ft <- data.table::copy(src)
  ft_m <- if ("m" %in% names(ft)) as.integer(ft[["m"]]) else NA_integer_
  is_selected_vec <- !is.na(selected_m) & !is.na(ft_m) & ft_m == selected_m
  ft[["is_selected"]] <- is_selected_vec

  # Add on_frontier and in_band flags for the include_dominated view
  if (isTRUE(include_dominated)) {
    ft[["on_frontier"]] <- !is.na(ft_m) & ft_m %in% frontier_m
    if (band_used) {
      # Use natural-scale effect for the band threshold
      hr_nat_for_band <- if (effect_log_scale && "hr" %in% names(ft)) {
        exp(as.numeric(ft$hr))
      } else if ("hr" %in% names(ft)) {
        as.numeric(ft$hr)
      } else {
        rep(NA_real_, nrow(ft))
      }
      floor_v <- (1 - effect_nbhd_v) * max(hr_nat_for_band, na.rm = TRUE)
      ft[["in_band"]] <- hr_nat_for_band >= floor_v
    }
  }

  # --- 4. Convert effect to natural scale, rename column ---------------
  if (effect_log_scale && "hr" %in% names(ft)) {
    ft[["hr"]] <- exp(as.numeric(ft[["hr"]]))
  }
  if ("hr" %in% names(ft)) {
    data.table::setnames(ft, "hr", effect_measure)
  }

  # --- 4b. Optional confidence intervals (merge user-supplied ci_table) ----
  # No computation here -- the user has already called
  # compute_frontier_cis() and passes the result.  This keeps the table and
  # plot in sync when they share the same ci_table object.
  ci_cols <- character(0)
  if (!is.null(ci_table)) {
    if (!data.table::is.data.table(ci_table) || nrow(ci_table) == 0L) {
      warning("ci_table is not a non-empty data.table; ignoring.",
              call. = FALSE)
    } else if (!"m" %in% names(ci_table)) {
      warning("ci_table has no 'm' column; cannot merge with frontier. ",
              "Pass the data.table returned by compute_frontier_cis(fs).",
              call. = FALSE)
    } else {
      # Defensive: only request columns that exist on the supplied object.
      # This lets a caller pass a hand-trimmed ci_table without crashing.
      ci_cols_want <- intersect(
        c("m", "naive_lcl", "naive_ucl",
          "split_lcl", "split_ucl",
          "fsbc_estimate", "fsbc_lcl", "fsbc_ucl"),
        names(ci_table))
      # Coerce both join keys to integer.  frontier$m on some fs objects
      # is character (legacy artifact of do.call(rbind, ...) coercion in
      # subgroup.consistency); ci_table$m is always integer.  Without
      # this coercion the join errors with a type-mismatch.
      ft[["m"]] <- suppressWarnings(as.integer(ft[["m"]]))
      ci_local <- ci_table[, ci_cols_want, with = FALSE]
      ci_local[["m"]] <- suppressWarnings(as.integer(ci_local[["m"]]))
      ft <- merge(
        ft,
        ci_local,
        by = "m", all.x = TRUE, sort = FALSE
      )
      # Pre-format the CI strings on the natural scale.  Width is fixed
      # at digits_ci for both bounds so columns align visually.
      fmt <- function(x) {
        if (is.na(x)) "NA" else
          formatC(x, format = "f", digits = digits_ci)
      }
      if (all(c("naive_lcl", "naive_ucl") %in% names(ft))) {
        ft[["Naive 95% CI"]] <- vapply(seq_len(nrow(ft)), function(k) {
          if (is.na(ft$naive_lcl[k]) || is.na(ft$naive_ucl[k])) "NA"
          else sprintf("(%s, %s)",
                       fmt(ft$naive_lcl[k]), fmt(ft$naive_ucl[k]))
        }, character(1))
        ci_cols <- c(ci_cols, "Naive 95% CI")
      }
      if (all(c("split_lcl", "split_ucl") %in% names(ft))) {
        ft[["Split ~ 95% CI"]] <- vapply(seq_len(nrow(ft)), function(k) {
          if (is.na(ft$split_lcl[k]) || is.na(ft$split_ucl[k])) "NA"
          else sprintf("(%s, %s)",
                       fmt(ft$split_lcl[k]), fmt(ft$split_ucl[k]))
        }, character(1))
        ci_cols <- c(ci_cols, "Split ~ 95% CI")
      }
      # FSBC-mimic column: show "est (lcl, ucl)" so the bias-corrected
      # estimate appears alongside its interval.  Tilde signals
      # "approximation; not full FSBC".  The point estimate uses
      # digits_effect (matching the HR column), while the bounds use
      # digits_ci (matching the other CI cells).
      if (all(c("fsbc_estimate", "fsbc_lcl",
                "fsbc_ucl") %in% names(ft))) {
        fmt_pt <- function(x) {
          if (is.na(x)) "NA" else
            formatC(x, format = "f", digits = digits_effect)
        }
        ft[["FSBC ~ 95% CI"]] <- vapply(seq_len(nrow(ft)), function(k) {
          if (is.na(ft$fsbc_estimate[k]) ||
              is.na(ft$fsbc_lcl[k]) ||
              is.na(ft$fsbc_ucl[k])) {
            return("NA")
          }
          sprintf("%s (%s, %s)",
                  fmt_pt(ft$fsbc_estimate[k]),
                  fmt(ft$fsbc_lcl[k]),
                  fmt(ft$fsbc_ucl[k]))
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
  # Flag columns are added when include_dominated = TRUE (they exist on ft
  # only in that branch)
  flag_cols <- intersect(c("on_frontier", "in_band"), names(ft))
  display_cols <- c("is_selected", m_cols, metric_cols, flag_cols, ci_cols)
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

  # Convert on_frontier / in_band logical columns to star markers for gt.
  # Logical columns are kept as-is in the data.table format (handled
  # upstream); only the gt rendering path converts them.
  if ("on_frontier" %in% names(ft_display)) {
    ft_display[["on_frontier"]] <- ifelse(
      isTRUE(ft_display[["on_frontier"]]) | ft_display[["on_frontier"]],
      "\u2605", "")
  }
  if ("in_band" %in% names(ft_display)) {
    ft_display[["in_band"]] <- ifelse(
      is.na(ft_display[["in_band"]]), "",
      ifelse(ft_display[["in_band"]], "\u2605", ""))
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

  # Relabel flag columns and center-align them
  if ("on_frontier" %in% names(ft_display)) {
    tbl <- tbl |>
      gt::cols_label(on_frontier = "Frontier") |>
      gt::cols_align(align = "center", columns = "on_frontier")
  }
  if ("in_band" %in% names(ft_display)) {
    tbl <- tbl |>
      gt::cols_label(in_band = "InBand") |>
      gt::cols_align(align = "center", columns = "in_band")
  }

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
