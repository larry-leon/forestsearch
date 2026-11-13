#' Multiplier-resampling (MR) estimates table (bootstrap-table format)
#'
#' Renders the selected subgroup's MR de-biased estimate ([forestsearch()] with
#' `mr_inference = TRUE`) in the `gt` layout of the bias-corrected table
#' produced by [summarize_bootstrap_results()].
#'
#' **With a bootstrap template** (`boot_results` supplied): the descriptive
#' columns (sample size, medians/RMST or rates) and the naive effect do not
#' depend on which de-biasing method was used, so they are inherited from a
#' full-bootstrap (FB) result used as a template; the bias-corrected column is
#' replaced with the MR values from `fs$mr_inference` -- the de-biased
#' estimate for the selected subgroup and, when MR was run with
#' `include_complement = TRUE`, for the complement as well (otherwise the
#' complement cell is marked, since MR de-biases the selected subgroup only).
#'
#' **Without a template** (`boot_results = NULL`, the default): the table is
#' built from `fs` alone. The descriptive columns are computed from
#' `fs$df.est` (the analysis data with the declared rule's `treat.recommend`)
#' by the builder the bootstrap uses ([SG_tab_estimates()]; for GLM outcomes
#' its GLM counterpart), one row for the selected subgroup and one for its
#' complement in the analysis data. The naive effect is `fs$mr_inference$naive`
#' (and `$complement$naive`) with its robust (sandwich) limits, and its header
#' says so. The MR columns follow: the de-biased estimate with its IJ interval
#' (`debiased`); when `fs$mr_inference$field` exists, the field one-sided 95%
#' bound (`field$lower_1s` for the selected subgroup, `field$complement$upper_1s`
#' for the complement); and when the joint block exists, the Bonferroni pair
#' (`field$joint$bonf_lower_H`, `field$joint$bonf_upper_Hc`). The field and
#' Bonferroni columns appear only when those quantities exist. Complement
#' cells carry MR values only when MR was run with `include_complement = TRUE`
#' (and the complement field with `field_complement = TRUE`); otherwise they
#' are marked with an em dash, as with a template.
#'
#' @param fs A [forestsearch()] result with a non-NULL `mr_inference` element.
#' @param boot_results `NULL` (default) to build the table from `fs` alone, or
#'   a [forestsearch_bootstrap_dofuture()] result used as the template for the
#'   descriptive columns, the naive effect, the complement row, and the
#'   subgroup footnote.  A template must contain `FSsg_tab` and `SG_CIs`.
#' @param est.scale Character, `"hr"` or `"1/hr"`; passed to
#'   [format_bootstrap_table()].  With a template it should match the scale of
#'   `boot_results`; without one, `"1/hr"` inverts every ratio (limits swap,
#'   and a one-sided lower bound becomes an upper bound) and is accepted for
#'   ratio measures only.
#' @param digits Integer display precision (default `2`), applied with the same
#'   reformatter [summarize_bootstrap_results()] uses so the two tables align.
#'
#' @return A `gt` table.  Without a template it carries the attribute
#'   `"mr_values"`: a data frame, one row per subgroup (`row` `"H"`, `"Hc"`),
#'   of the numbers the table displays, on the effect scale of
#'   `fs$mr_inference` (before any `"1/hr"` inversion): `naive_est`,
#'   `naive_lower`, `naive_upper`, `debiased_est`, `ij_lower`, `ij_upper`, and,
#'   when present, `field_bound` and `bonferroni_bound` (`NA` where a
#'   complement quantity was not computed; the complement's naive columns are
#'   then the descriptive table's own robust Cox fit, not stored here).
#'
#' @details
#' With a template, the selected-subgroup row is matched by its bias-corrected string
#' (`boot_results$SG_CIs$H_bc`) rather than by position, so the substitution is
#' robust to row ordering.  The MR de-biased confidence interval is a
#' leading-order approximation of the full-bootstrap infinitesimal-jackknife
#' interval -- it comes from the IJ variance computed on the multiplier draws
#' under `ci_method = "ij"` or `"field"` (the forest-search path defaults to
#' `"field"`, which adds the field bounds; the DINA and GRF paths default to
#' `"ij"`), and from the subgroup robust SE only under `ci_method = "wald"`.
#' Either way the bias-corrected cell may differ
#' slightly from the full bootstrap even when the naive effect matches exactly.
#'
#' For GLM outcomes the effect columns follow the fit's orientation: the
#' naive and MR values are those of `fs$mr_inference`, on the scale of the
#' fit's own estimates; the descriptive table's estimator (which gives the
#' complement's naive cell when MR did not cover the complement) is built
#' with the fit's stored `adverse_outcome` or, when none is stored, with
#' [forestsearch()]'s rule (`TRUE` for binary and count outcomes, `FALSE`
#' otherwise); and a template's effect column comes from
#' [forestsearch_bootstrap_dofuture()], which builds its estimator with the
#' same rule.  With the default `est.scale = "hr"`, the descriptive columns
#' `Rate(C)` or `Mean(C)`, `Rate(T)` or `Mean(T)`, and their raw difference
#' `Diff` (treated minus control) are the raw outcome and do not follow it.
#' So for an unadjusted fit where a
#' lower outcome is adverse (`adverse_outcome = FALSE`, the default for
#' continuous outcomes) the effect and `Diff` have opposite signs: an MD is
#' `-Diff`, an RD is `-Diff / 100`, and the log of an OR or RR has the sign
#' opposite to `Diff`.  Count outcomes and the rate measures (IRR, IRD) are
#' not reoriented (see [make_effect_estimator()]).
#'
#' @seealso [forestsearch()] for the `mr_inference` switch and the vocabulary
#'   section; [fs_mr_inference()] for the MR computation itself;
#'   [forestsearch_bootstrap_dofuture()] for the full bootstrap that supplies
#'   the template; [fs_fdr_report()] for harm-confirmation threshold sweeps;
#'   [summarize_bootstrap_results()], [format_bootstrap_table()]
#' @export
mr_estimates_table <- function(fs, boot_results = NULL, est.scale = "hr",
                               digits = 2) {

  g <- fs$mr_inference
  if (is.null(g)) {
    stop("fs$mr_inference is NULL; run forestsearch(mr_inference = TRUE).",
         call. = FALSE)
  }
  if (is.null(boot_results)) {
    return(.mr_table_from_fs(fs, g, est.scale = est.scale, digits = digits))
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
  # the complement (Hc) row, which MR does not de-bias.
  h_row  <- which(tab[[bc_col]] == boot_results$SG_CIs$H_bc)
  hc_row <- which(tab[[bc_col]] == boot_results$SG_CIs$Hc_bc)
  if (length(h_row) != 1L) {
    stop("Could not uniquely match the selected-subgroup row in FSsg_tab.",
         call. = FALSE)
  }

  # MR de-biased CI at construction precision; the shared reformatter below
  # renders it at `digits`, matching the bootstrap table.
  ci_str <- function(e, l, u) {
    paste0(round(e, 4), " (", round(l, 4), ",", round(u, 4), ")")
  }
  tab[[bc_col]][h_row] <- ci_str(g$debiased$est, g$debiased$lower,
                                 g$debiased$upper)

  # Complement: use MR's complement de-biased estimate when available
  # (include_complement = TRUE); otherwise mark it, since MR de-biases the
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
  # The bias-corrected footnote names the method used here, not the
  # bootstrap's: the MR de-biased estimate, its interval from the IJ variance
  # on the multiplier draws (the robust SE under ci_method = "wald").
  bc_measure <- sub("\\*$", "", bc_col)
  se_txt <- if (identical(g$ci_method, "wald")) {
    "the subgroup robust (sandwich) standard error"
  } else {
    "the infinitesimal-jackknife (IJ) variance computed on the multiplier draws"
  }
  bc_footnote <- sprintf(paste0(
    "**Bias-corrected %s**: multiplier-resampling (MR) de-biased estimate ",
    "(%s draws), with its 95%% interval from %s. Corrects for optimism in ",
    "subgroup selection."),
    bc_measure, if (is.na(draws)) "?" else format(draws), se_txt)
  format_bootstrap_table(
    FSsg_tab          = tab,
    nb_boots          = draws,
    est.scale         = est.scale,
    boot_success_rate = NULL,
    sg_definition     = paste(fs$sg.harm, collapse = " & "),
    title             = "Treatment Effect by Subgroup",
    subtitle          = sprintf(
      "MR de-biased estimate (multiplier resampling, %s draws)%s",
      if (is.na(draws)) "?" else format(draws), cc_note),
    bc_footnote       = bc_footnote)
}


#' The MR estimates table built from the forestsearch() result alone
#'
#' The `boot_results = NULL` route of [mr_estimates_table()]: the descriptive
#' columns from `fs$df.est` by the builder the bootstrap uses
#' ([SG_tab_estimates()] or `SG_tab_estimates_glm()`), the naive effect and
#' the MR columns from `fs$mr_inference`.
#' @noRd
.mr_table_from_fs <- function(fs, g, est.scale = "hr", digits = 2) {

  if (is.null(g$debiased) || is.null(g$naive)) {
    stop("fs$mr_inference carries no estimates (",
         if (is.null(g$note)) "short return" else paste(g$note, collapse = "; "),
         "); nothing to tabulate.", call. = FALSE)
  }
  df <- fs$df.est
  if (is.null(df) || !"treat.recommend" %in% names(df)) {
    stop("fs$df.est with a treat.recommend column is required.", call. = FALSE)
  }
  est.scale <- match.arg(est.scale, c("hr", "1/hr"))
  log_scale <- isTRUE(g$log_scale)
  if (est.scale == "1/hr" && !log_scale) {
    stop("est.scale = \"1/hr\" applies to ratio measures only.", call. = FALSE)
  }
  a <- fs$args_call_all
  outcome_type <- if (is.null(a$outcome_type)) "survival" else a$outcome_type[1L]

  # Descriptive columns, laid out as the bootstrap's FSsg_tab (rows: Qstnbl =
  # H, treat.recommend == 0; Recmnd = the complement), without the
  # bias-corrected column.
  if (outcome_type == "survival") {
    tab <- SG_tab_estimates(
      df = df, SG_flag = "treat.recommend", draws = 0, details = FALSE,
      outcome.name = a$outcome.name, event.name = a$event.name,
      treat.name = a$treat.name, strata.name = NULL,
      potentialOutcome.name = a$potentialOutcome.name,
      est.scale = est.scale, sg0_name = "Qstnbl", sg1_name = "Recmnd")
  } else {
    # Orientation: the fit's resolved adverse_outcome, so the effect cells are
    # the fit's own estimates, sign included.  NULL only on an object not
    # built by forestsearch(); resolved by forestsearch()'s rule.
    adverse_outcome <- a$adverse_outcome
    if (is.null(adverse_outcome)) {
      adverse_outcome <- outcome_type %in% c("binary", "count")
    }
    est_fn <- make_effect_estimator(
      outcome_type = outcome_type, treat.name = a$treat.name,
      outcome.name = a$outcome.name, event.name = NULL,
      offset.name = a$offset.name, effect_measure = a$effect_measure,
      adverse_outcome = adverse_outcome,
      adjust_covariates = a$adjust_covariates)
    tab <- SG_tab_estimates_glm(
      df = df, SG_flag = "treat.recommend", outcome.name = a$outcome.name,
      treat.name = a$treat.name, estimator_fn = est_fn,
      effect_measure = a$effect_measure, outcome_type = outcome_type,
      sg1_name = "Recmnd", sg0_name = "Qstnbl", est.scale = est.scale,
      digits = 4)
  }
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  eff_col <- grep("^(HR|IRR|OR|RR|RD|IRD|MD)\\s*\\(.*CI\\)$", colnames(tab),
                  value = TRUE)[1]
  if (is.na(eff_col)) {
    stop("Could not locate the effect column of the descriptive table.",
         call. = FALSE)
  }
  h_row  <- which(tab$Subgroup == "Qstnbl")
  hc_row <- which(tab$Subgroup == "Recmnd")
  if (length(h_row) != 1L || length(hc_row) != 1L) {
    stop("Could not match the selected-subgroup and complement rows.",
         call. = FALSE)
  }

  # On est.scale = "1/hr" every ratio is inverted and its limits swap; a
  # one-sided lower bound becomes an upper bound.
  inv <- est.scale == "1/hr"
  tr <- function(x) if (inv) 1 / x else x
  num <- function(x) if (is.null(x) || length(x) != 1L) NA_real_ else as.numeric(x)

  cc <- g$complement
  has_cc <- !is.null(cc) && !is.null(cc$debiased) && is.finite(num(cc$debiased$est))
  fl <- g$field
  has_fl <- !is.null(fl) && is.finite(num(fl$lower_1s))
  fc <- if (has_fl) fl$complement else NULL
  has_fc <- has_cc && !is.null(fc) && is.finite(num(fc$upper_1s))
  jt <- if (has_fl) fl$joint else NULL
  has_jt <- !is.null(jt) && is.finite(num(jt$bonf_lower_H))

  # The numbers of the table, on the effect scale of fs$mr_inference (before
  # any 1/hr inversion); returned as the "mr_values" attribute.
  vals <- data.frame(
    row = c("H", "Hc"), Subgroup = c("Qstnbl", "Recmnd"),
    naive_est   = c(num(g$naive$est),   if (has_cc) num(cc$naive$est)   else NA_real_),
    naive_lower = c(num(g$naive$lower), if (has_cc) num(cc$naive$lower) else NA_real_),
    naive_upper = c(num(g$naive$upper), if (has_cc) num(cc$naive$upper) else NA_real_),
    debiased_est = c(num(g$debiased$est),   if (has_cc) num(cc$debiased$est)   else NA_real_),
    ij_lower     = c(num(g$debiased$lower), if (has_cc) num(cc$debiased$lower) else NA_real_),
    ij_upper     = c(num(g$debiased$upper), if (has_cc) num(cc$debiased$upper) else NA_real_),
    stringsAsFactors = FALSE)
  if (has_fl) {
    vals$field_bound <- c(num(fl$lower_1s), if (has_fc) num(fc$upper_1s) else NA_real_)
  }
  if (has_jt) {
    vals$bonferroni_bound <- c(num(jt$bonf_lower_H),
                               if (has_fc) num(jt$bonf_upper_Hc) else NA_real_)
  }

  f <- function(x) formatC(tr(x), format = "f", digits = digits)
  ci <- function(e, l, u) {
    lo <- if (inv) u else l; up <- if (inv) l else u
    paste0(f(e), " (", f(lo), ", ", f(up), ")")
  }
  dash <- "\u2014"
  # A one-sided bound, with its side: H's lower bound, the complement's upper
  # bound on the HR scale; reversed on 1/hr.
  bnd <- function(x, lower) {
    if (!is.finite(x)) return(dash)
    side <- if (xor(lower, inv)) "\u2265 " else "\u2264 "
    paste0(side, f(x))
  }

  # The naive effect: fs$mr_inference's estimate with its robust limits. For
  # the complement without include_complement, the descriptive table's own
  # robust Cox fit (the same estimate) stays.
  tab[[eff_col]] <- as.character(tab[[eff_col]])
  tab[[eff_col]][h_row] <- ci(vals$naive_est[1], vals$naive_lower[1], vals$naive_upper[1])
  if (has_cc) {
    tab[[eff_col]][hc_row] <- ci(vals$naive_est[2], vals$naive_lower[2], vals$naive_upper[2])
  }
  reformat <- .reformat_FSsg_tab_digits
  if (!is.null(digits)) tab <- tryCatch(reformat(tab, digits = digits), error = function(e) tab)

  col_by_row <- function(h, hc) { v <- rep(dash, nrow(tab)); v[h_row] <- h; v[hc_row] <- hc; v }
  tab$MR_debiased <- col_by_row(
    ci(vals$debiased_est[1], vals$ij_lower[1], vals$ij_upper[1]),
    if (has_cc) ci(vals$debiased_est[2], vals$ij_lower[2], vals$ij_upper[2]) else dash)
  if (has_fl) {
    tab$MR_field <- col_by_row(bnd(vals$field_bound[1], TRUE), bnd(vals$field_bound[2], FALSE))
  }
  if (has_jt) {
    tab$MR_bonferroni <- col_by_row(bnd(vals$bonferroni_bound[1], TRUE),
                                    bnd(vals$bonferroni_bound[2], FALSE))
  }

  draws <- if (!is.null(g$settings$draws)) g$settings$draws else NA_integer_
  measure <- sub("\\s*\\(.*", "", eff_col)
  tbl <- format_bootstrap_table(
    FSsg_tab          = tab,
    nb_boots          = draws,
    est.scale         = est.scale,
    boot_success_rate = NULL,
    sg_definition     = paste(fs$sg.harm, collapse = " & "),
    title             = "Treatment Effect by Subgroup",
    subtitle          = sprintf(
      "Naive and multiplier-resampling (MR) estimates (%s draws)%s",
      if (is.na(draws)) "?" else format(draws),
      if (has_cc) "" else "; complement not de-biased"))

  mr_cols <- intersect(c("MR_debiased", "MR_field", "MR_bonferroni"), colnames(tab))
  lab <- list()
  lab[[eff_col]] <- gt::md(paste0(measure, "<br/>(95% CI, robust)<sup>\u2020</sup>"))
  lab$MR_debiased <- gt::md("De-biased<br/>(95% IJ CI)")
  if (has_fl) lab$MR_field <- gt::md("Field<br/>one-sided 95%")
  if (has_jt) lab$MR_bonferroni <- gt::md("Bonferroni pair<br/>(joint 95%)")
  tbl <- gt::cols_label(tbl, .list = lab)
  tbl <- gt::tab_spanner(tbl, label = "Multiplier resampling (MR)", columns = mr_cols)
  tbl <- gt::tab_footnote(
    tbl,
    footnote = gt::md(paste0(
      "**MR**: the selected estimate de-biased for its selection over the ",
      "re-selections of multiplier draws within the fixed candidate family; ",
      "two-sided infinitesimal-jackknife (IJ) interval",
      if (has_fl) paste0("; field: one-sided 95% bound from the nested Monte Carlo ",
                         "of the re-selection (lower for the selected subgroup, upper ",
                         "for its complement on the HR scale)") else "",
      if (has_jt) "; Bonferroni pair: the two bounds each at 97.5%, jointly 95%" else "",
      if (has_cc) "." else "; the complement is not de-biased (include_complement = FALSE).")),
    locations = gt::cells_column_spanners(spanners = "Multiplier resampling (MR)"))

  attr(tbl, "mr_values") <- vals
  tbl
}
