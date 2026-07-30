#' Tier-2 de-biased estimates table (bootstrap-table format)
#'
#' Renders the selected subgroup's fast de-biased gate estimate
#' ([forestsearch()] with `mr_inference = TRUE`) in the same `gt` layout as the
#' bias-corrected table produced by [summarize_bootstrap_results()].  The
#' descriptive columns (sample size, medians/RMST or rates) and the naive effect
#' are tier-independent, so they are inherited from a bootstrap result used as a
#' template; the bias-corrected column is replaced with the Tier-2 values from
#' `fs$mr_inference` -- the de-biased estimate for the selected subgroup and, when
#' the gate was run with `include_complement = TRUE`, for the complement as well
#' (otherwise the complement cell is marked, since the gate de-biases the
#' selected subgroup only).
#'
#' @param fs A [forestsearch()] result with a non-NULL `mr_inference` element.
#' @param boot_results A [forestsearch_bootstrap_dofuture()] result used as the
#'   template for the descriptive columns, the naive effect, the complement row,
#'   and the subgroup footnote.  Must contain `FSsg_tab` and `SG_CIs`.
#' @param est.scale Character, `"hr"` or `"1/hr"`; passed to
#'   [format_bootstrap_table()].  Should match the scale of `boot_results`.
#' @param digits Integer display precision (default `2`), applied with the same
#'   reformatter [summarize_bootstrap_results()] uses so the two tables align.
#'
#' @return A `gt` table.
#'
#' @details
#' The selected-subgroup row is matched by its bias-corrected string
#' (`boot_results$SG_CIs$H_bc`) rather than by position, so the substitution is
#' robust to row ordering.  The Tier-2 de-biased confidence interval uses the
#' subgroup robust SE and is an approximation of the bootstrap
#' infinitesimal-jackknife interval, so the bias-corrected cell may differ
#' slightly from the full bootstrap even when the naive effect matches exactly.
#'
#' @seealso [summarize_bootstrap_results()], [format_bootstrap_table()]
#' @export
mr_estimates_table <- function(fs, boot_results, est.scale = "hr",
                               digits = 2) {

  g <- fs$mr_inference
  if (is.null(g)) {
    stop("fs$mr_inference is NULL; run forestsearch(mr_inference = TRUE).",
         call. = FALSE)
  }
  if (is.null(boot_results$FSsg_tab) || is.null(boot_results$SG_CIs)) {
    stop("boot_results must contain FSsg_tab and SG_CIs to serve as a template.",
         call. = FALSE)
  }

  tab <- as.data.frame(boot_results$FSsg_tab, stringsAsFactors = FALSE)
  cn  <- colnames(tab)

  # Bias-corrected column ("HR*", "OR*", "IRR*", ...).
  bc_col <- grep("(HR|IRR|OR|RR|RD|IRD|MD)\\*", cn, value = TRUE)[1]
  if (is.na(bc_col)) {
    stop("Could not locate a bias-corrected column in boot_results$FSsg_tab.",
         call. = FALSE)
  }

  # Match the selected (H) row by the template's H bias-corrected string; mark
  # the complement (Hc) row, which the gate does not de-bias.
  h_row  <- which(tab[[bc_col]] == boot_results$SG_CIs$H_bc)
  hc_row <- which(tab[[bc_col]] == boot_results$SG_CIs$Hc_bc)
  if (length(h_row) != 1L) {
    stop("Could not uniquely match the selected-subgroup row in FSsg_tab.",
         call. = FALSE)
  }

  # Tier-2 de-biased CI at construction precision; the shared reformatter below
  # renders it at `digits`, matching the bootstrap table.
  ci_str <- function(e, l, u) {
    paste0(round(e, 4), " (", round(l, 4), ",", round(u, 4), ")")
  }
  tab[[bc_col]][h_row] <- ci_str(g$debiased$est, g$debiased$lower,
                                 g$debiased$upper)

  # Complement: use the gate's complement de-biased estimate when available
  # (include_complement = TRUE); otherwise mark it, since the gate de-biases the
  # selected subgroup only.
  cc <- g$complement
  has_cc <- !is.null(cc) && !is.null(cc$debiased) && is.finite(cc$debiased$est)
  if (length(hc_row) == 1L && has_cc) {
    tab[[bc_col]][hc_row] <- ci_str(cc$debiased$est, cc$debiased$lower,
                                    cc$debiased$upper)
  }

  reformat <- tryCatch(
    get(".reformat_FSsg_tab_digits", envir = asNamespace("forestsearch")),
    error = function(e) NULL)
  if (!is.null(reformat) && !is.null(digits)) {
    tab <- tryCatch(reformat(tab, digits = digits), error = function(e) tab)
  }

  # If the complement was not de-biased, mark it after reformatting
  # (an em dash is not a CI string).
  if (length(hc_row) == 1L && !has_cc) tab[[bc_col]][hc_row] <- "\u2014"

  draws <- if (!is.null(g$settings$draws)) g$settings$draws else NA_integer_
  cc_note <- if (has_cc) "" else "; complement not de-biased"
  format_bootstrap_table(
    FSsg_tab          = tab,
    nb_boots          = draws,
    est.scale         = est.scale,
    boot_success_rate = NULL,
    sg_definition     = paste(fs$sg.harm, collapse = " & "),
    title             = "Treatment Effect by Subgroup",
    subtitle          = sprintf(
      "Tier-2 de-biased estimate (multiplier bootstrap, %s draws)%s",
      if (is.na(draws)) "?" else format(draws), cc_note))
}
