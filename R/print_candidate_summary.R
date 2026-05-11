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
                                    max_width        = NULL) {
  if (is.null(out_sg) || is.null(out_sg$result) ||
      nrow(out_sg$result) == 0L) {
    cat("\nNo candidates passed the consistency threshold.\n")
    cat(sprintf("Evaluated: %d  Passed: 0\n", n_evaluated))
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
  pcons_strs   <- formatC(as.numeric(res$Pcons), format = "f", digits = 3,
                          width = 6)
  flag <- function(b, width = 3L) {
    ch <- ifelse(is.na(b), ".",
                 ifelse(isTRUE(b) | b == TRUE, "*", "-"))
    formatC(ch, width = width, flag = " ", big.mark = "")
  }
  # Header widths for the flag columns (centered)
  fr_w   <- 5L   # "Front"
  band_w <- 6L   # "InBand"
  sel_w  <- 3L   # "Sel"
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
  fr_strs   <- center_flag(is_frontier, fr_w)
  sel_strs  <- center_flag(is_selected, sel_w)
  band_strs <- if (band_used) center_flag(in_band, band_w) else NULL

  # Header ------------------------------------------------------------------
  hdr_cols <- c("Rank", effect_label, "N", "E", "K", "Pcons", "Front",
                if (band_used) "InBand", "Sel", "Subgroup")
  widths   <- c(4, max(6L, nchar(effect_label)), 5, 4, 3, 6, 5,
                if (band_used) 6, 3)
  # Width of "Subgroup" col absorbs leftover space; compute total
  fixed_total <- sum(widths) + 2 * length(widths)   # padding
  width_lim   <- if (is.null(max_width)) 110L else max_width
  cuts_w <- max(20L, width_lim - fixed_total - 4L)
  # Truncate cuts to fit
  cuts_disp <- vapply(cuts, function(s) {
    if (nchar(s) > cuts_w) paste0(substr(s, 1L, cuts_w - 3L), "...") else s
  }, character(1))

  # Banner ------------------------------------------------------------------
  rule_str <- sprintf("sg_focus = %s, selection_rule = %s",
                      shQuote(sg_focus, type = "cmd"),
                      shQuote(selection_rule, type = "cmd"))
  bar <- paste(rep("=", width_lim), collapse = "")
  thin <- paste(rep("-", width_lim), collapse = "")
  cat("\n", bar, "\n", sep = "")
  cat("CANDIDATE EVALUATION SUMMARY  (", rule_str, ")\n", sep = "")
  cat(bar, "\n", sep = "")

  # Header row
  hdr_parts <- c(
    formatC("Rank",         width = 4,  flag = "-"),
    formatC(effect_label,   width = max(6L, nchar(effect_label)), flag = "-"),
    formatC("N",            width = 5,  flag = "-"),
    formatC("E",            width = 4,  flag = "-"),
    formatC("K",            width = 3,  flag = "-"),
    formatC("Pcons",        width = 6,  flag = "-"),
    formatC("Front",        width = 5,  flag = "-"),
    if (band_used) formatC("InBand", width = 6, flag = "-"),
    formatC("Sel",          width = 3,  flag = "-"),
    "Subgroup"
  )
  cat(paste(hdr_parts, collapse = "  "), "\n", sep = "")
  cat(thin, "\n", sep = "")

  # Body rows
  for (k in seq_len(n)) {
    row_parts <- c(
      rank_strs[k], hr_strs[k], N_strs[k], E_strs[k], K_strs[k],
      pcons_strs[k], fr_strs[k],
      if (band_used) band_strs[k],
      sel_strs[k], cuts_disp[k]
    )
    cat(paste(row_parts, collapse = "  "), "\n", sep = "")
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
  cat(bar, "\n\n", sep = "")
  invisible(NULL)
}
