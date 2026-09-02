# plot_subgroup_sims.R ----------------------------------------------------
# Figure construction for the extreme-subgroups study: a plot() method
# assembling the gg_forest() calls of the committed vignettes (commit
# 002f4f37) parameter-for-parameter, and forest_height() converting
# panel row counts to the fig.height values the vignette chunks use.
#
# Every default below -- axis limits, tick positions, reference-line
# styling, annotation labels and rounding, footnote wording -- is lifted
# from the corresponding vignette chunk for that (metric, panel) pair,
# so the ported documents render identical figures.  All of it is
# overridable through the arguments or `...` (which is merged over the
# defaults before the gg_forest() call).

#' Recommended figure height for a forest panel
#'
#' Converts a panel's row count to the `fig.height` (inches) used by the
#' vignette chunks: `rows * row_height + overhead`, rounded to 0.1.
#'
#' @param x A `"subgroup_sims_summary"` object.
#' @param panel `"single"`, `"combo"`, or `"highrisk"`.
#' @param row_height Inches per row (default 0.45).
#' @param overhead Fixed inches for title, x-axis, and footnote
#'   (default 1.5).
#'
#' @return A single numeric, rounded to one decimal.
#' @seealso [plot.subgroup_sims_summary()]
#' @export
#' @examples
#' \dontrun{
#' fh_single <- forest_height(S, "single")
#' }
forest_height <- function(x, panel = c("single", "combo", "highrisk"),
                          row_height = 0.45, overhead = 1.5) {
  stopifnot(inherits(x, "subgroup_sims_summary"))
  panel <- match.arg(panel)
  rows <- switch(panel,
                 single   = x$n_single,
                 combo    = x$n_combo,
                 highrisk = x$highrisk$n)
  round(rows * row_height + overhead, 1)
}

# Internal: footnote text per (metric, panel), exactly as in the
# vignette chunks.  `np` (per-trial N) appears only in the hr/single
# footnote; when unknown (NA) that clause is dropped.
.subgroup_sims_footnote <- function(metric, panel, n_sims, hr_true,
                                    np, thr) {
  ht <- sprintf("%.2f", hr_true)
  if (metric == "hr") {
    if (panel == "single" && !is.na(np)) {
      paste0("Median HR (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials, N = ", np, " per trial  |  ",
             "True HR = ", ht, " (red dashed)")
    } else if (panel %in% c("single", "combo")) {
      paste0("Median HR (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials  |  ",
             "True HR = ", ht, " (red dashed)")
    } else {
      paste0("Same subgroups and ordering as the UB(HR) panel above  |  ",
             n_sims, " simulated trials  |  ",
             "Arrowheads mark whiskers extending beyond the axis  |  ",
             "Dashed: true HR = ", ht, ";  Dotted: HR = 1.0")
    }
  } else {
    if (panel %in% c("single", "combo")) {
      paste0("Median UB(HR) (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials  |  ",
             "Dotted line: UB = 1.0;  Dashed line: true HR = ", ht)
    } else {
      paste0("Filter: Pr(UB(HR) >= 2.0) >= ", round(100 * thr), "%  |  ",
             "ITT included for comparison  |  ",
             n_sims, " simulated trials  |  ",
             "Arrowheads mark whiskers extending beyond the axis  |  ",
             "Dotted: UB = 1.0;  Dashed: true HR = ", ht)
    }
  }
}

#' Forest plot of an extreme-subgroups summary
#'
#' Draws one of the vignettes' six forest panels from a
#' [summary.subgroup_sims()] object via [gg_forest()]: the empirical-CI
#' distribution of the HR point estimate or of the upper 95% bound
#' UB(HR), for the single-variable panel, the combination + random-
#' benchmark panel, or the high-risk panel (Pr(UB>=2) filter with the
#' ITT row anchored first).
#'
#' Defaults reproduce the committed vignette figures exactly:
#' metric-specific axis limits, tick positions, reference/vertical line
#' styling, the four annotation columns (N plus tail probabilities and
#' the unconditional median), and the panel-specific footnotes. Any
#' default can be overridden through the named arguments or `...`,
#' which is merged over the assembled argument list before the
#' [gg_forest()] call (so e.g. `ref_col = "black"` or
#' `widths = ...` in `...` win).
#'
#' @param x A `"subgroup_sims_summary"` object.
#' @param metric `"hr"` (HR point estimate) or `"ub"` (upper 95% bound).
#' @param panel `"single"`, `"combo"`, or `"highrisk"`.
#' @param hr_true True hazard ratio for the reference/vertical line and
#'   footnotes; defaults to the value stored on `x` and must be
#'   non-`NULL`.
#' @param base_size,point_size,line_size Passed to [gg_forest()];
#'   defaults match the vignettes' `fp_*` values.
#' @param widths Column width vector passed to [gg_forest()].
#' @param clip_marker Passed to [gg_forest()]; `"arrow"` marks whiskers
#'   clipped at the axis limits.
#' @param xlim,ticks_at,xlab,footnote `NULL` (default) selects the
#'   vignette value for the chosen `metric`/`panel`.
#' @param ... Further arguments merged over the defaults and passed to
#'   [gg_forest()].
#'
#' @return The [gg_forest()] plot object.
#' @seealso [forest_height()], [summary.subgroup_sims()]
#' @export
#' @examples
#' \dontrun{
#' S <- summary(sims_uniform, hr_true = 0.70)
#' plot(S, metric = "ub", panel = "combo")
#' }
plot.subgroup_sims_summary <- function(x,
                                       metric = c("hr", "ub"),
                                       panel = c("single", "combo",
                                                 "highrisk"),
                                       hr_true = x$hr_true,
                                       base_size = 14,
                                       point_size = 2.2,
                                       line_size = 0.75,
                                       widths = c(3.5, 5, 1.1, 1.2, 1.2,
                                                  1.1),
                                       clip_marker = "arrow",
                                       xlim = NULL,
                                       ticks_at = NULL,
                                       xlab = NULL,
                                       footnote = NULL,
                                       ...) {
  metric <- match.arg(metric)
  panel  <- match.arg(panel)
  if (is.null(hr_true)) {
    stop("`hr_true` is required for the reference lines and footnotes; ",
         "pass it here or to summary().")
  }

  q <- if (metric == "hr") x$hr_q else x$ub_q

  # Row selection per panel, exactly as in the vignette chunks
  if (panel == "highrisk") {
    sel <- x$highrisk$include
    ord <- x$highrisk$ord
  } else {
    sel <- if (metric == "hr") x$ok else x$ok_ub
    ord <- if (metric == "hr") {
      if (panel == "single") x$idx_hr_single else x$idx_hr_combo
    } else {
      if (panel == "single") x$idx_ub_single else x$idx_ub_combo
    }
  }
  pick <- function(v) v[sel][ord]

  # Annotation block per metric (labels, rounding as in the chunks)
  annot <- if (metric == "hr") {
    list(
      "N"      = as.character(pick(x$mean_n)),
      "HR<0.5" = paste0(round(100 * pick(x$pr_hr_lt050), 1), "%"),
      "HR>1.0" = paste0(round(100 * pick(x$pr_hr_gt1),   1), "%"),
      "mHR"    = pick(x$mhr_uncond)
    )
  } else {
    list(
      "N"     = as.character(pick(x$mean_n)),
      "UB>=2" = paste0(round(100 * pick(x$pr_ub_ge2), 1), "%"),
      "UB>=3" = paste0(round(100 * pick(x$pr_ub_ge3), 1), "%"),
      "mUB"   = pick(x$mub_uncond)
    )
  }

  args <- list(
    subgroups   = pick(x$sg_names),
    est         = q[sel, 2][ord],
    lo          = q[sel, 1][ord],
    hi          = q[sel, 3][ord],
    cat_vec     = pick(x$cat),
    cat_colours = x$cat_cols,
    annot       = annot,
    point_size  = point_size,
    line_size   = line_size,
    base_size   = base_size,
    widths      = widths,
    clip_marker = clip_marker,
    xlog        = TRUE
  )

  if (metric == "hr") {
    args$ref_line   <- hr_true
    args$vert_lines <- 1.00
    args$ref_col    <- "firebrick"
    args$vert_col   <- "grey55"
    args$vert_lty   <- "dotted"
    if (panel == "highrisk") args$ref_lty <- "dashed"
    args$xlim     <- if (is.null(xlim)) c(0.15, 3.5) else xlim
    args$ticks_at <- if (is.null(ticks_at)) {
      c(0.20, 0.35, 0.50, 0.70, 1.00, 1.50, 2.50)
    } else ticks_at
    args$xlab <- if (is.null(xlab)) {
      "Hazard Ratio (Cox, stratified by grade)"
    } else xlab
  } else {
    args$ref_line   <- 1.00
    args$vert_lines <- hr_true
    args$ref_col    <- "grey55"
    args$ref_lty    <- "dotted"
    args$vert_col   <- "firebrick"
    args$vert_lty   <- "dashed"
    args$xlim     <- if (is.null(xlim)) c(0.30, 9.0) else xlim
    args$ticks_at <- if (is.null(ticks_at)) {
      c(0.40, 0.70, 1.00, 1.50, 2.50, 4.00, 8.00)
    } else ticks_at
    args$xlab <- if (is.null(xlab)) {
      "Upper Bound of 95% CI  [UB(HR)]"
    } else xlab
  }

  args$footnote <- if (is.null(footnote)) {
    .subgroup_sims_footnote(metric, panel, x$n_sims, hr_true,
                            np  = if (is.null(x$n_per_trial)) NA_integer_
                                  else x$n_per_trial,
                            thr = x$highrisk$threshold)
  } else footnote

  dots <- list(...)
  if (length(dots)) args[names(dots)] <- dots

  do.call(gg_forest, args)
}
