# Silence R CMD check NOTE for bare data.table column reference.
utils::globalVariables(c("m"))

#' Describe the Pareto-Frontier Selection in Words
#'
#' Wraps \code{\link{pareto_frontier_table}} with a verbal account of
#' the selection: why the chosen subgroup satisfies the configured
#' \code{sg_focus} + \code{selection_rule} criterion, and, for each
#' non-dominated alternative, on which lexicographic axis it lost to
#' the winner.
#'
#' The wrapper is \strong{descriptive}, not algorithmic: it reads the
#' frontier and the configuration off the returned \code{fs} object and
#' renders prose.  It does not re-run the selection.  A future change
#' to the lexicographic key in \code{\link{forestsearch}} would require
#' updating this wrapper to remain accurate.
#'
#' @section Verbosity:
#' \describe{
#'   \item{\code{"compact"}}{Criterion sentence + selected line + a
#'     one-line bullet per loser.}
#'   \item{\code{"medium"} (default)}{Criterion sentence + selected
#'     paragraph (including the Naive CI if \code{ci_table} is given)
#'     + one short paragraph per loser, each naming the lex-key axis
#'     on which the loser was beaten.}
#'   \item{\code{"detailed"}}{Medium + a closing paragraph that
#'     reproduces the formal lex key, the resolved indicator values,
#'     and a note about any ties broken downstream.}
#' }
#'
#' @section Loser ordering:
#' Losers are sorted by "closeness-to-winning", which is operationalised
#' as: \strong{later loss-axis comes first}, then \strong{smaller gap on
#' that axis}.  Specifically, a loser that lost on \code{K} (the last
#' tiebreaker) is reported before one that lost on \code{N} (an early
#' axis), and a loser that beat the winner on every axis until the
#' final one is reported before a loser that fell off the band entirely.
#' This puts the boundary cases first, which is the most informative
#' order for understanding how the rule made its choice.
#'
#' @param fs A \code{forestsearch} result object.
#' @param ci_table Optional data.table from
#'   \code{\link{compute_frontier_cis}}.  When supplied, the Naive
#'   95\% CI is woven into the selected-candidate paragraph (only).
#'   Loser paragraphs do not reference CIs by design -- the
#'   selection rule does not consider CIs.
#' @param effect_label Optional character override for the effect-
#'   measure label in prose, e.g. \code{"Hazard ratio"}.  Defaults to
#'   a sensible string for the outcome type.
#' @param digits Integer.  Decimal places for effect estimates and
#'   bound values in the prose.  Default \code{3L}.
#' @param verbosity One of \code{"compact"}, \code{"medium"} (default),
#'   or \code{"detailed"}.
#' @param format One of \code{"markdown"} (default), \code{"character"},
#'   or \code{"list"}.
#'
#' @return Depending on \code{format}:
#'   \describe{
#'     \item{\code{"markdown"}}{A single string with embedded line
#'       breaks and bullet markers, ready for \code{cat()} in a Quarto
#'       chunk with \code{results = "asis"}.}
#'     \item{\code{"character"}}{A character vector with one element
#'       per logical paragraph or bullet (no markdown markers).}
#'     \item{\code{"list"}}{A structured list with named elements
#'       \code{criterion}, \code{selected}, \code{losers},
#'       \code{closing}, and \code{meta}.}
#'   }
#'
#' @examples
#' \dontrun{
#' ci_tab <- compute_frontier_cis(fs, n_splits = 1000, seed = 1)
#'
#' # Inline in a Quarto chunk
#' cat(explain_pareto_selection(fs, ci_table = ci_tab))
#'
#' # Or get the structured object
#' explain_pareto_selection(fs, ci_table = ci_tab, format = "list")
#' }
#'
#' @seealso \code{\link{pareto_frontier_table}},
#'   \code{\link{plot_pareto_frontier}},
#'   \code{\link{compute_frontier_cis}}.
#' @export
explain_pareto_selection <- function(fs,
                                     ci_table     = NULL,
                                     effect_label = NULL,
                                     digits       = 3L,
                                     verbosity    = c("medium", "compact", "detailed"),
                                     format       = c("markdown", "character", "list")) {
  verbosity <- match.arg(verbosity)
  format    <- match.arg(format)

  # --- 1. Config + frontier + winner --------------------------------------
  cfg <- .epx_parse_config(fs)
  out_sg <- tryCatch(fs$grp.consistency$out_sg, error = function(e) NULL)
  if (is.null(out_sg) || is.null(out_sg$pareto_frontier)) {
    return(.epx_empty_return(format, "No Pareto frontier on this fs object."))
  }
  effect_measure <- fs$effect_measure %||% "HR"
  if (is.null(effect_label)) effect_label <- .epx_effect_label(effect_measure)

  ft <- .epx_tidy_frontier(out_sg$pareto_frontier, effect_measure)
  if (is.null(ft) || nrow(ft) == 0L) {
    return(.epx_empty_return(format, "Frontier is empty."))
  }

  selected_m <- tryCatch(as.integer(out_sg$result[1L, ]$m),
                         error = function(e) NA_integer_)
  if (is.na(selected_m)) {
    return(.epx_empty_return(format, "Could not identify the selected subgroup."))
  }
  selected_on_frontier <- selected_m %in% ft$m

  # --- 2. Compute the band indicator for hrMaxSG / hrMinSG ---------------
  band_relevant <- cfg$sg_focus %in% c("hrMaxSG", "hrMinSG") &&
                   cfg$selection_rule %in% c("neighborhood", "both")
  if (band_relevant) {
    hr_max <- max(ft$hr_nat, na.rm = TRUE)
    floor_v <- (1 - cfg$effect_neighborhood) * hr_max
    ft$in_band <- ft$hr_nat >= floor_v
  } else {
    ft$in_band <- rep(TRUE, nrow(ft))
    floor_v <- NA_real_
  }

  # --- 3. Identify winner + loser rows -----------------------------------
  if (selected_on_frontier) {
    winner <- ft[ft$m == selected_m, , drop = FALSE]
    losers <- ft[ft$m != selected_m, , drop = FALSE]
  } else {
    winner <- NULL
    losers <- ft
  }

  # --- 4. Classify each loser's loss-axis + closeness sort --------------
  loser_list <- if (is.null(winner) || nrow(losers) == 0L) {
    list()
  } else {
    Ls <- lapply(seq_len(nrow(losers)), function(k) {
      L <- losers[k, ]
      cls <- .epx_classify_loss(L, winner, cfg$sg_focus)
      list(row = L, loss = cls)
    })
    .epx_order_losers(Ls)
  }

  # --- 5. Build sentences -----------------------------------------------
  parts <- list(
    criterion = .epx_criterion(cfg, floor_v, band_relevant,
                               effect_label, digits),
    selected  = .epx_selected(winner, cfg, ci_table, effect_label,
                              digits, selected_on_frontier),
    losers    = lapply(loser_list, function(L)
      .epx_loser(L, winner, cfg, effect_label, digits, verbosity,
                 band_relevant = band_relevant)),
    closing   = if (verbosity == "detailed")
                  .epx_closing(cfg, ft, selected_on_frontier)
                else NULL,
    meta      = list(
      n_frontier           = nrow(ft),
      selected_on_frontier = selected_on_frontier,
      band_relevant        = band_relevant,
      floor_v              = floor_v,
      has_ci               = !is.null(ci_table)
    )
  )

  # --- 6. Format dispatch -----------------------------------------------
  .epx_render(parts, verbosity, format)
}


# =========================================================================
# Internal helpers
# =========================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Capitalize first letter (for sentence-start formatting).
.epx_capitalize_first <- function(s) {
  if (is.null(s) || !nzchar(s)) return(s)
  paste0(toupper(substring(s, 1L, 1L)), substring(s, 2L))
}

# Parse selection configuration from fs$args_call_all
.epx_parse_config <- function(fs) {
  args <- fs$args_call_all %||% list()
  list(
    sg_focus            = args$sg_focus            %||% "hr",
    selection_rule      = args$selection_rule      %||% "neighborhood",
    effect_neighborhood = args$effect_neighborhood %||% 0.10,
    hr_threshold        = args$hr.threshold        %||% NA_real_,
    p_threshold         = args$pconsistency.threshold %||% NA_real_
  )
}

# Label for prose
.epx_effect_label <- function(em) {
  switch(em,
    HR  = "hazard ratio",
    OR  = "odds ratio",
    RR  = "risk ratio",
    IRR = "incidence rate ratio",
    RD  = "risk difference",
    IRD = "incidence rate difference",
    MD  = "mean difference",
    em)
}

# Tidy the frontier into a data.frame with cuts_str and hr_nat
.epx_tidy_frontier <- function(frontier, effect_measure) {
  if (!data.table::is.data.table(frontier) && !is.data.frame(frontier)) {
    return(NULL)
  }
  ft <- as.data.frame(frontier)
  if (!nrow(ft)) return(ft)
  ft$m <- suppressWarnings(as.integer(ft$m))
  m_cols <- grep("^M\\.", names(ft), value = TRUE)
  ft$cuts_str <- vapply(seq_len(nrow(ft)), function(k) {
    cuts <- unlist(ft[k, m_cols, drop = TRUE], use.names = FALSE)
    cuts <- cuts[!is.na(cuts) & nzchar(cuts)]
    if (length(cuts) == 0L) "(empty subgroup)"
    else paste(cuts, collapse = " AND ")
  }, character(1))
  is_log <- effect_measure %in% c("OR", "RR", "IRR")
  ft$hr_nat <- if (is_log) exp(ft$hr) else ft$hr
  ft
}

# Classify a loser's loss-axis vs the winner under the given sg_focus.
# Returns list(axis, gap) where axis is in c("band","N","Pcons","effect","K","tie")
# and gap is the magnitude of the loss on that axis (NA for "band" and "tie").
.epx_classify_loss <- function(loser, winner, sg_focus) {
  # 1. Band axis (only relevant for hrMaxSG / hrMinSG with band-using rule)
  if (isTRUE(winner$in_band) && !isTRUE(loser$in_band)) {
    return(list(axis = "band", gap = NA_real_,
                gap_label = "outside the effect-size band"))
  }
  # 2. N axis (direction depends on focus)
  N_loses <- switch(sg_focus,
    maxSG    = loser$N <  winner$N,
    hrMaxSG  = loser$N <  winner$N,
    minSG    = loser$N >  winner$N,
    hrMinSG  = loser$N >  winner$N,
    FALSE)
  if (N_loses) {
    gap <- abs(loser$N - winner$N)
    direction_label <- if (sg_focus %in% c("maxSG","hrMaxSG"))
      sprintf("smaller subgroup (N=%d vs N=%d, %d fewer subjects)",
              loser$N, winner$N, gap)
    else
      sprintf("larger subgroup (N=%d vs N=%d, %d more subjects)",
              loser$N, winner$N, gap)
    return(list(axis = "N", gap = gap, gap_label = direction_label))
  }
  if (loser$N != winner$N) {
    # Loser's N is BETTER than winner's; cannot lose on N. Fall through.
  }
  # 3. Pcons axis
  if (loser$Pcons < winner$Pcons) {
    gap <- winner$Pcons - loser$Pcons
    return(list(axis = "Pcons", gap = gap,
                gap_label = sprintf("lower consistency (Pcons=%.3f vs %.3f)",
                                    loser$Pcons, winner$Pcons)))
  }
  # 4. Effect axis (lower effect is worse under hrMaxSG / hrMinSG / hr)
  if (loser$hr_nat < winner$hr_nat) {
    gap <- winner$hr_nat - loser$hr_nat
    return(list(axis = "effect", gap = gap,
                gap_label = sprintf("weaker effect (%.3f vs %.3f)",
                                    loser$hr_nat, winner$hr_nat)))
  }
  # 5. K axis (larger K is worse)
  if (loser$K > winner$K) {
    gap <- loser$K - winner$K
    return(list(axis = "K", gap = gap,
                gap_label = sprintf("less parsimonious (K=%d vs K=%d)",
                                    loser$K, winner$K)))
  }
  list(axis = "tie", gap = NA_real_,
       gap_label = "tied on all lex-key axes (downstream tiebreaker)")
}

# Order losers by closeness-to-winning: later lex axis first, smaller gap first.
.epx_order_losers <- function(losers) {
  if (length(losers) == 0L) return(losers)
  axis_rank <- c(band = 1L, N = 2L, Pcons = 3L, effect = 4L, K = 5L, tie = 6L)
  ranks <- vapply(losers, function(L) axis_rank[L$loss$axis], integer(1))
  gaps  <- vapply(losers, function(L) {
    g <- L$loss$gap
    if (is.null(g) || is.na(g)) Inf else g
  }, numeric(1))
  ord <- order(-ranks, gaps)  # higher rank first, smaller gap first
  losers[ord]
}

# Build criterion sentence
.epx_criterion <- function(cfg, floor_v, band_relevant, effect_label, digits) {
  # Focus description depends on whether the rule actually enforces a band
  focus_desc <- if (band_relevant) {
    switch(cfg$sg_focus,
      hrMaxSG = sprintf("largest-N within the effect-size band (sg_focus = \"hrMaxSG\")"),
      hrMinSG = sprintf("smallest-N within the effect-size band (sg_focus = \"hrMinSG\")"),
      cfg$sg_focus)
  } else {
    switch(cfg$sg_focus,
      hr      = sprintf("highest-consistency (sg_focus = \"hr\")"),
      maxSG   = sprintf("largest-N (sg_focus = \"maxSG\")"),
      minSG   = sprintf("smallest-N (sg_focus = \"minSG\")"),
      hrMaxSG = sprintf("largest-N among non-dominated candidates (sg_focus = \"hrMaxSG\")"),
      hrMinSG = sprintf("smallest-N among non-dominated candidates (sg_focus = \"hrMinSG\")"),
      cfg$sg_focus)
  }

  rule_desc <- switch(cfg$selection_rule,
    neighborhood = "the effect-size neighborhood rule",
    pareto       = "the Pareto-non-dominance rule",
    both         = "the intersection of effect-band and Pareto-non-dominance rules",
    cfg$selection_rule)

  band_clause <- if (band_relevant) {
    sprintf(" The effect-size band requires %s within %d%% of the maximum on the frontier (i.e., %s >= %s).",
            effect_label, round(100 * cfg$effect_neighborhood),
            effect_label, formatC(floor_v, format = "f", digits = digits))
  } else ""

  tie_clause <- if (cfg$sg_focus %in% c("hrMaxSG","hrMinSG")) {
    " Within the eligible set, ties on N are broken by Pcons, then by effect strength, then by parsimony (smaller K)."
  } else if (cfg$sg_focus %in% c("maxSG","minSG")) {
    " Ties on N are broken by Pcons, then by parsimony (smaller K)."
  } else ""

  plain <- sprintf("**Selection criterion.** The configured criterion is %s under %s.%s%s",
                   focus_desc, rule_desc, band_clause, tie_clause)
  list(plain = plain)
}

# Build the selected paragraph
.epx_selected <- function(winner, cfg, ci_table, effect_label, digits,
                          selected_on_frontier) {
  if (!selected_on_frontier) {
    return(list(plain = paste0(
      "**Selected subgroup.** The selected subgroup does not appear on the ",
      "Pareto frontier shown above. This is allowed under `hrMinSG`, which ",
      "may pick a small-N candidate that is itself dominated by a larger-N ",
      "candidate with similar effect. The frontier rows below describe the ",
      "non-dominated alternatives that were evaluated but not chosen.")))
  }

  fmt_e <- function(x) formatC(x, format = "f", digits = digits)
  fmt_p <- function(x) formatC(x, format = "f", digits = digits)
  cuts <- winner$cuts_str
  hr   <- winner$hr_nat
  N    <- winner$N
  Pc   <- winner$Pcons
  K    <- winner$K

  # CI snippet (Naive only, per design)
  ci_snippet <- ""
  if (!is.null(ci_table) && data.table::is.data.table(ci_table)) {
    ci_row <- ci_table[ci_table$m == as.integer(winner$m), ]
    if (nrow(ci_row) == 1L &&
        !is.na(ci_row$naive_lcl) && !is.na(ci_row$naive_ucl)) {
      ci_snippet <- sprintf(
        " The Naive 95%% CI for the %s on this subgroup is (%s, %s).",
        effect_label,
        formatC(ci_row$naive_lcl, format = "f", digits = digits),
        formatC(ci_row$naive_ucl, format = "f", digits = digits))
    }
  }

  reason <- switch(cfg$sg_focus,
    hrMaxSG = "the largest-N candidate in the eligible set",
    hrMinSG = "the smallest-N candidate in the eligible set",
    maxSG   = "the largest-N candidate",
    minSG   = "the smallest-N candidate",
    hr      = "the candidate with the highest Pcons",
    "the candidate satisfying the selection rule")

  plain <- sprintf(
    paste0("**Selected subgroup.** %s (%s = %s, N = %d, Pcons = %s, K = %d). ",
           "This is %s among non-dominated frontier candidates.%s"),
    cuts, effect_label, fmt_e(hr), N, fmt_p(Pc), K, reason, ci_snippet)
  list(plain = plain)
}

# Build one loser paragraph
.epx_loser <- function(L, winner, cfg, effect_label, digits, verbosity,
                       band_relevant = TRUE) {
  fmt_e <- function(x) formatC(x, format = "f", digits = digits)
  fmt_p <- function(x) formatC(x, format = "f", digits = digits)
  row <- L$row
  loss <- L$loss

  header <- sprintf("**%s** (%s = %s, N = %d, Pcons = %s, K = %d)",
                    row$cuts_str, effect_label, fmt_e(row$hr_nat),
                    row$N, fmt_p(row$Pcons), row$K)

  # Build phrasing that depends on whether a band was enforced
  n_tiebreaker_clause <- if (band_relevant)
    "N is the primary tiebreaker once the band is satisfied"
  else
    "N is the primary tiebreaker among non-dominated candidates"

  why <- switch(loss$axis,
    band   = sprintf(
      "fell outside the effect-size band (%s = %s below the floor); excluded before the N tiebreaker even applies.",
      effect_label, fmt_e(row$hr_nat)),
    N      = sprintf(
      "%s. %s, so this candidate's %s could not overcome the size loss.",
      .epx_capitalize_first(loss$gap_label),
      n_tiebreaker_clause,
      if (row$hr_nat > winner$hr_nat) "stronger effect" else "comparable effect"),
    Pcons  = sprintf("%s; Pcons is the secondary tiebreaker after N.",
                     .epx_capitalize_first(loss$gap_label)),
    effect = sprintf("%s; effect strength is the tertiary tiebreaker.",
                     .epx_capitalize_first(loss$gap_label)),
    K      = sprintf("%s; K is the final tiebreaker.",
                     .epx_capitalize_first(loss$gap_label)),
    tie    = sprintf("%s.", .epx_capitalize_first(loss$gap_label)))

  if (verbosity == "compact") {
    plain <- sprintf("- %s -- %s", header, why)
  } else {
    plain <- sprintf("- %s. %s", header, why)
  }
  list(plain = plain, loss_axis = loss$axis)
}

# Build closing paragraph for detailed verbosity
.epx_closing <- function(cfg, ft, selected_on_frontier) {
  key_text <- switch(cfg$sg_focus,
    hr      = "(-Pcons, -effect, +K)",
    maxSG   = "(-N, -Pcons, +K)",
    minSG   = "(+N, -Pcons, +K)",
    hrMaxSG = "(-in_band, -N, -Pcons, -effect, +K)",
    hrMinSG = "(-in_band, +N, -Pcons, -effect, +K)",
    "(see forestsearch documentation)")

  detail <- sprintf(
    paste0("**Lexicographic key.** The ranking key for sg_focus = \"%s\" is ",
           "%s, ordered ascending. The winner is the lex-minimum across all ",
           "%d frontier candidates."),
    cfg$sg_focus, key_text, nrow(ft))

  symmetry_note <- if (cfg$sg_focus == "hrMinSG" &&
                       cfg$selection_rule %in% c("pareto", "both")) {
    " Note: under hrMinSG, the `pareto` and `both` rules always agree on the winner (see methodology paper, Proposition 3.4)."
  } else ""

  list(plain = paste0(detail, symmetry_note))
}

# Format dispatch
.epx_render <- function(parts, verbosity, format) {
  if (format == "list") return(parts)

  pieces <- character(0)
  pieces <- c(pieces, parts$criterion$plain)
  pieces <- c(pieces, parts$selected$plain)

  if (length(parts$losers) > 0L) {
    header <- if (verbosity == "compact") "**Non-selected frontier candidates:**"
              else "**Non-selected frontier candidates** (ordered by closeness to winning):"
    pieces <- c(pieces, header)
    for (L in parts$losers) pieces <- c(pieces, L$plain)
  } else {
    pieces <- c(pieces,
      "**Non-selected candidates.** The frontier has a single candidate, so no comparison applies.")
  }

  if (!is.null(parts$closing)) pieces <- c(pieces, parts$closing$plain)

  if (format == "character") return(pieces)

  # markdown: join with double newline
  paste(pieces, collapse = "\n\n")
}

# Empty/abort return
.epx_empty_return <- function(format, msg) {
  switch(format,
    markdown  = msg,
    character = msg,
    list      = list(error = msg)
  )
}
