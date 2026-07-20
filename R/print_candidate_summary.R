# ============================================================================
# Candidate summary printer
# ============================================================================
# Renders a plain-text summary table of all passing candidates after
# the consistency stage, with frontier / band / selected flags.  Called
# from subgroup.consistency() when show_candidate_summary = TRUE.
# ----------------------------------------------------------------------------

#' Print a candidate summary table (internal)
#'
#' Plain-text summary of all candidates that passed the consistency
#' threshold, with flags for Pareto-frontier membership,
#' effect-neighborhood band eligibility (when applicable), and the
#' selected subgroup.  Columns are computed from
#' \code{out_sg$result} and \code{out_sg$pareto_frontier}; no
#' additional data structures are required.
#'
#' Pcons is shown only for passing candidates (Path 1 / light-touch
#' design): the consistency evaluator returns NULL for non-passing
#' candidates, so their Pcons values are not retained.  The footer
#' line reports the evaluated vs.\ passed counts so the reader can
#' see how the filter operated.
#'
#' @param out_sg The \code{out_sg} list returned by
#'   \code{sg_consistency_out()}; must contain \code{result} (passing
#'   candidates, lex-key sorted) and \code{pareto_frontier} (subset).
#' @param n_evaluated Integer; number of candidates that went through
#'   consistency evaluation (= length of \code{results_list} = top
#'   \code{max_subgroups_search} candidates).
#' @param sg_focus Character; selection focus
#'   (\code{"hr"}/\code{"maxSG"}/\code{"minSG"}/\code{"hrMaxSG"}/\code{"hrMinSG"}).
#' @param selection_rule Character; rule for choice
#'   (\code{"neighborhood"}/\code{"pareto"}/\code{"both"}).
#' @param effect_neighborhood Numeric; band tolerance, used to compute
#'   the in-band indicator when applicable.
#' @param effect_log_scale Logical; whether \code{out_sg$result$hr} is
#'   stored on the log scale (TRUE for OR/RR/IRR).
#' @param effect_label Character; label for the effect column (e.g.\ "HR").
#' @param max_width Integer; total table width in characters.  Default
#'   \code{NULL} = auto.
#'
#' @keywords internal
#' @noRd
print_candidate_summary <- function(out_sg,
                                    n_evaluated,
                                    sg_focus,
                                    selection_rule,
                                    effect_neighborhood,
                                    effect_log_scale = FALSE,
                                    effect_label     = "HR",
                                    max_width        = NULL,
                                    max_print        = Inf) {
  # ---------------------------------------------------------------------------
  # Always emit the SUMMARY banner, even in the zero-passed case.
  # Downstream tooling (extract_candidate_diagnostics() in
  # compare_selection_rules.R) keys on the banner phrase to slice the
  # captured stdout into the "post-consistency summary" block.  Without
  # the banner, a legitimate "no candidates passed" outcome renders as
  # "(not available)" in the comparison vignette, which mis-suggests
  # that the call never reached the summary step.
  # ---------------------------------------------------------------------------
  width_lim <- if (!is.null(max_width)) max_width else 110L
  bar       <- paste(rep("=", width_lim), collapse = "")
  rule_str  <- sprintf("sg_focus = %s, selection_rule = %s",
                       shQuote(sg_focus,        type = "cmd"),
                       shQuote(selection_rule,  type = "cmd"))
  cat("\n", bar, "\n", sep = "")
  cat("CANDIDATE EVALUATION SUMMARY  (", rule_str, ")\n", sep = "")
  cat(bar, "\n", sep = "")

  # Empty-result branch: emit a compact diagnostic body and return.
  if (is.null(out_sg) || is.null(out_sg$result) ||
      nrow(out_sg$result) == 0L) {
    cat(sprintf("Evaluated: %d   Passed: 0\n", n_evaluated))
    cat("No candidates met the consistency threshold.\n")
    cat(bar, "\n", sep = "")
    return(invisible(NULL))
  }

  res <- data.table::copy(out_sg$result)
  fr  <- out_sg$pareto_frontier

  # Natural-scale effect column
  hr_nat <- if (isTRUE(effect_log_scale)) exp(as.numeric(res$hr))
            else as.numeric(res$hr)

  # Build M.* cut strings
  m_cols <- grep("^M\\.", names(res), value = TRUE)
  cuts <- if (length(m_cols) == 0L) {
    rep("(no cuts)", nrow(res))
  } else {
    vapply(seq_len(nrow(res)), function(k) {
      v <- unlist(res[k, m_cols, with = FALSE], use.names = FALSE)
      v <- v[!is.na(v) & nzchar(v)]
      if (length(v) == 0L) "(empty)" else paste(v, collapse = " & ")
    }, character(1))
  }

  # Frontier flag
  is_frontier <- if (is.null(fr) || nrow(fr) == 0L) rep(FALSE, nrow(res))
                 else as.integer(res$m) %in% as.integer(fr$m)

  # Selected flag (row 1 of out_sg$result is the winner per sort_subgroups)
  is_selected <- c(TRUE, rep(FALSE, nrow(res) - 1L))

  # Under sg_focus = "maxeff", consistency is a property of the WINNER only
  # (per-candidate Pcons is not computed -- see the maxeff fast path in
  # subgroup.consistency()).  Suppress the per-candidate P column, which would be
  # NA for every non-winner, and report the selected subgroup's Pcons in the
  # footer instead.  All other foci keep the P column unchanged.
  show_pcons <- !identical(sg_focus, "maxeff")

  # Band flag (only when the rule uses one)
  band_used <- sg_focus %in% c("hrMaxSG", "hrMinSG") &&
               selection_rule %in% c("neighborhood", "both")
  in_band <- if (band_used) {
    floor_v <- (1 - effect_neighborhood) * max(hr_nat, na.rm = TRUE)
    hr_nat >= floor_v
  } else {
    rep(NA, nrow(res))
  }

  # Column formatting -------------------------------------------------------
  n   <- nrow(res)
  rank_strs    <- formatC(seq_len(n),       format = "d", width = 4)
  hr_strs      <- formatC(hr_nat,           format = "f", digits = 3, width = 6)
  N_strs       <- formatC(as.integer(res$N), format = "d", width = 5)
  E_strs       <- formatC(as.integer(res$E), format = "d", width = 4)
  K_strs       <- formatC(as.integer(res$K), format = "d", width = 3)
  flag <- function(b, width = 3L) {
    ch <- ifelse(is.na(b), ".",
                 ifelse(isTRUE(b) | b == TRUE, "*", "-"))
    formatC(ch, width = width, flag = " ", big.mark = "")
  }
  # Header widths for the flag/scalar columns.  Headers abbreviated:
  #   Pcons -> "P",  Frontier -> "OF",  InBand -> "IB",  Selected -> "S".
  # See legend printed below the table for definitions.
  pcons_hdr_w <- 5L   # "P" header, data needs 5 chars at 3-dp Pcons
  fr_w        <- 2L   # "OF"
  band_w      <- 2L   # "IB"
  sel_w       <- 1L   # "S"
  # Pad the single character to the column width with leading/trailing spaces
  center_flag <- function(b, width) {
    out <- character(length(b))
    for (i in seq_along(b)) {
      ch <- if (is.na(b[i])) "." else if (isTRUE(b[i]) | b[i] == TRUE) "*" else "-"
      left  <- floor((width - 1) / 2)
      right <- width - 1 - left
      out[i] <- paste0(strrep(" ", left), ch, strrep(" ", right))
    }
    out
  }
  # Pcons strings: center the values within pcons_hdr_w
  pcons_strs <- formatC(as.numeric(res$Pcons), format = "f", digits = 3,
                        width = pcons_hdr_w)
  fr_strs    <- center_flag(is_frontier, fr_w)
  sel_strs   <- center_flag(is_selected, sel_w)
  band_strs  <- if (band_used) center_flag(in_band, band_w) else NULL

  # Header ------------------------------------------------------------------
  hdr_cols <- c("Rank", effect_label, "N", "E", "K",
                if (show_pcons) "P", "OF",
                if (band_used) "IB", "S", "Subgroup")
  widths   <- c(4, max(6L, nchar(effect_label)), 5, 4, 3,
                if (show_pcons) pcons_hdr_w, fr_w,
                if (band_used) band_w, sel_w)
  # Width of "Subgroup" col absorbs leftover space; compute total
  fixed_total <- sum(widths) + 2 * length(widths)   # padding
  width_lim   <- if (is.null(max_width)) 110L else max_width
  cuts_w <- max(20L, width_lim - fixed_total - 4L)
  # Truncate cuts to fit
  cuts_disp <- vapply(cuts, function(s) {
    if (nchar(s) > cuts_w) paste0(substr(s, 1L, cuts_w - 3L), "...") else s
  }, character(1))

  # Table separators (banner already emitted at top of function)
  thin <- paste(rep("-", width_lim), collapse = "")

  # Header row
  hdr_parts <- c(
    formatC("Rank",         width = 4,          flag = "-"),
    formatC(effect_label,   width = max(6L, nchar(effect_label)), flag = "-"),
    formatC("N",            width = 5,          flag = "-"),
    formatC("E",            width = 4,          flag = "-"),
    formatC("K",            width = 3,          flag = "-"),
    if (show_pcons) formatC("P", width = pcons_hdr_w, flag = "-"),
    formatC("OF",           width = fr_w,        flag = "-"),
    if (band_used) formatC("IB", width = band_w, flag = "-"),
    formatC("S",            width = sel_w,       flag = "-"),
    "Subgroup"
  )
  cat(paste(hdr_parts, collapse = "  "), "\n", sep = "")
  cat(thin, "\n", sep = "")

  # Body rows (capped at max_print; candidates are sorted selection-first, so
  # the selected subgroup is always shown).
  n_show <- min(n, max_print)
  for (k in seq_len(n_show)) {
    row_parts <- c(
      rank_strs[k], hr_strs[k], N_strs[k], E_strs[k], K_strs[k],
      if (show_pcons) pcons_strs[k], fr_strs[k],
      if (band_used) band_strs[k],
      sel_strs[k], cuts_disp[k]
    )
    cat(paste(row_parts, collapse = "  "), "\n", sep = "")
  }
  if (n_show < n) {
    cat(sprintf("... %d more candidate%s not shown (max_print = %s).\n",
                n - n_show, if (n - n_show == 1L) "" else "s",
                format(max_print)))
  }
  cat(thin, "\n", sep = "")

  # Footer
  n_front <- sum(is_frontier)
  n_band  <- if (band_used) sum(in_band, na.rm = TRUE) else NA_integer_
  footer  <- sprintf(
    "Evaluated: %d   Passed: %d   On frontier: %d%s   Selected: m=%s",
    n_evaluated, n, n_front,
    if (band_used) sprintf("   In band: %d", n_band) else "",
    as.character(res$m[1])
  )
  cat(footer, "\n", sep = "")
  # maxeff: report the selected subgroup's consistency here (it is the only one
  # computed) rather than as a per-candidate column.
  if (!show_pcons) {
    sel_pcons <- suppressWarnings(as.numeric(res$Pcons[1]))
    cat(sprintf("Selected subgroup Pcons (consistency probability): %s\n",
                if (is.na(sel_pcons)) "NA (inestimable)" else formatC(sel_pcons, format = "f", digits = 3)))
  }
  legend <- sprintf(
    "Legend: %sOF = on Pareto frontier;%s S = selected.",
    if (show_pcons) "P = Pcons (consistency probability); " else "",
    if (band_used) "  IB = in effect-size band;" else ""
  )
  cat(legend, "\n", sep = "")
  cat(bar, "\n\n", sep = "")
  invisible(NULL)
}


# ============================================================================
# Pre-consistency preview printer
# ============================================================================
# Prints all candidates that are about to be evaluated for consistency,
# in their preview-sort order (the same order used to apply the
# stop_Kgroups cutoff).  Called from subgroup.consistency() when
# show_candidate_summary = TRUE, BEFORE consistency runs.
# ----------------------------------------------------------------------------

#' Print a pre-consistency candidate preview table (internal)
#'
#' Plain-text preview of the candidates that will be evaluated for
#' consistency, in their preview-sort order.  Same column conventions
#' as \code{\link{print_candidate_summary}} but without Pcons or
#' Selected (which are not yet known).  Frontier and InBand flags are
#' computed pre-consistency from the candidate HR / N values.
#'
#' Use this together with \code{print_candidate_summary()} to give the
#' reader two views: which candidates went IN to consistency (this
#' function) and which came OUT (the summary).  The two views together
#' make the rule's filter visible end-to-end.
#'
#' @param found.hrs A data.table of candidates after preview sorting
#'   and \code{stop_Kgroups} trimming.  Must contain columns
#'   \code{HR}, \code{n}, \code{E}, \code{K}.
#' @param index.Z Matrix or data.table of factor indicators (rows
#'   align with \code{found.hrs}).
#' @param names.Z Character vector of factor column names.
#' @param confs_labels Named character vector for human-readable cut
#'   labels (passed to \code{FS_labels()}).
#' @param sg_focus,selection_rule,effect_neighborhood,effect_log_scale
#'   See \code{\link{print_candidate_summary}}.
#' @param effect_label Character; column label (e.g. \code{"HR"}).
#' @param max_width Integer; total table width.  Default \code{NULL} = auto.
#'
#' @keywords internal
#' @noRd
print_candidate_preview <- function(found.hrs,
                                    index.Z,
                                    names.Z,
                                    confs_labels,
                                    sg_focus,
                                    selection_rule,
                                    effect_neighborhood,
                                    effect_log_scale = FALSE,
                                    effect_label     = "HR",
                                    max_width        = NULL,
                                    max_print        = Inf) {
  if (is.null(found.hrs) || nrow(found.hrs) == 0L) {
    cat("\nNo candidates to preview.\n")
    return(invisible(NULL))
  }

  n <- nrow(found.hrs)

  # Natural-scale effect
  hr_raw <- as.numeric(found.hrs$HR)
  hr_nat <- if (isTRUE(effect_log_scale)) exp(hr_raw) else hr_raw

  N_vec <- as.integer(found.hrs$n)
  E_vec <- as.integer(found.hrs$E)
  K_vec <- as.integer(found.hrs$K)

  # Build subgroup-definition strings from index.Z and names.Z
  cuts_str <- character(n)
  for (i in seq_len(n)) {
    idx_i <- as.numeric(unlist(index.Z[i, ]))
    factors_i <- names.Z[idx_i == 1]
    labels_i <- vapply(factors_i, FS_labels, character(1),
                       confs_labels = confs_labels,
                       USE.NAMES = FALSE)
    cuts_str[i] <- if (length(labels_i) == 0L) "(empty)"
                   else paste(labels_i, collapse = " & ")
  }

  # In-band flag (pre-consistency, computed from HR vector of preview set)
  band_used <- sg_focus %in% c("hrMaxSG", "hrMinSG") &&
               selection_rule %in% c("neighborhood", "both")
  in_band <- if (band_used) {
    hr_max  <- max(hr_nat, na.rm = TRUE)
    floor_v <- (1 - effect_neighborhood) * hr_max
    hr_nat >= floor_v
  } else {
    rep(NA, n)
  }

  # Frontier flag (pre-consistency, computed from HR-N joint dominance)
  is_frontier <- !.pareto_dominated_xy(hr_nat, as.numeric(N_vec))

  # Column-width plumbing (mirrors print_candidate_summary)
  rank_strs  <- formatC(seq_len(n), format = "d", width = 4)
  hr_strs    <- formatC(hr_nat,     format = "f", digits = 3, width = 6)
  N_strs     <- formatC(N_vec,      format = "d", width = 5)
  E_strs     <- formatC(E_vec,      format = "d", width = 4)
  K_strs     <- formatC(K_vec,      format = "d", width = 3)
  center_flag <- function(b, width) {
    out <- character(length(b))
    for (i in seq_along(b)) {
      ch <- if (is.na(b[i])) "." else if (isTRUE(b[i]) | b[i] == TRUE) "*" else "-"
      left  <- floor((width - 1) / 2)
      right <- width - 1 - left
      out[i] <- paste0(strrep(" ", left), ch, strrep(" ", right))
    }
    out
  }
  fr_w <- 2L; band_w <- 2L
  fr_strs   <- center_flag(is_frontier, fr_w)
  band_strs <- if (band_used) center_flag(in_band, band_w) else NULL

  width_lim <- if (is.null(max_width)) 110L else max_width
  fixed_total <- 4 + max(6L, nchar(effect_label)) + 5 + 4 + 3 + fr_w +
                 (if (band_used) band_w else 0) +
                 2 * (6 + if (band_used) 1 else 0)
  cuts_w <- max(20L, width_lim - fixed_total - 4L)
  cuts_disp <- vapply(cuts_str, function(s) {
    if (nchar(s) > cuts_w) paste0(substr(s, 1L, cuts_w - 3L), "...") else s
  }, character(1))

  # Banner + header + body + footer
  rule_str <- sprintf("sg_focus = %s, selection_rule = %s",
                      shQuote(sg_focus, type = "cmd"),
                      shQuote(selection_rule, type = "cmd"))
  bar  <- paste(rep("=", width_lim), collapse = "")
  thin <- paste(rep("-", width_lim), collapse = "")
  cat("\n", bar, "\n", sep = "")
  cat("CANDIDATE EVALUATION PREVIEW (pre-consistency) (", rule_str, ")\n",
      sep = "")
  cat(bar, "\n", sep = "")
  hdr_parts <- c(
    formatC("Rank",         width = 4,    flag = "-"),
    formatC(effect_label,   width = max(6L, nchar(effect_label)), flag = "-"),
    formatC("N",            width = 5,    flag = "-"),
    formatC("E",            width = 4,    flag = "-"),
    formatC("K",            width = 3,    flag = "-"),
    formatC("OF",           width = fr_w, flag = "-"),
    if (band_used) formatC("IB", width = band_w, flag = "-"),
    "Subgroup"
  )
  cat(paste(hdr_parts, collapse = "  "), "\n", sep = "")
  cat(thin, "\n", sep = "")
  # Body rows (capped at max_print; sorted selection-first).
  n_show <- min(n, max_print)
  for (k in seq_len(n_show)) {
    row_parts <- c(
      rank_strs[k], hr_strs[k], N_strs[k], E_strs[k], K_strs[k],
      fr_strs[k],
      if (band_used) band_strs[k],
      cuts_disp[k]
    )
    cat(paste(row_parts, collapse = "  "), "\n", sep = "")
  }
  if (n_show < n) {
    cat(sprintf("... %d more candidate%s not shown (max_print = %s).\n",
                n - n_show, if (n - n_show == 1L) "" else "s",
                format(max_print)))
  }
  cat(thin, "\n", sep = "")
  n_front <- sum(is_frontier)
  n_band  <- if (band_used) sum(in_band, na.rm = TRUE) else NA_integer_
  footer <- sprintf(
    "To evaluate: %d   On frontier: %d%s",
    n, n_front,
    if (band_used) sprintf("   In band: %d", n_band) else ""
  )
  cat(footer, "\n", sep = "")
  legend <- sprintf(
    "Legend: OF = on Pareto frontier%s.",
    if (band_used) ";  IB = in effect-size band" else ""
  )
  cat(legend, "\n", sep = "")
  cat(bar, "\n\n", sep = "")
  invisible(NULL)
}
