#' Plot the Pareto Frontier of Candidate Subgroups
#'
#' Produces a 2D scatter of candidate subgroups in (effect, N) space,
#' with the Pareto frontier drawn as a step polyline and the selected
#' subgroup highlighted.  When \code{ci_table} is supplied (the return
#' value of \code{\link{compute_frontier_cis}}), horizontal error bars
#' for the split-derived 95\% CI are overlaid on each frontier member.
#'
#' @param fs A \code{forestsearch} result object.
#' @param ci_table Optional \code{data.table} from
#'   \code{\link{compute_frontier_cis}}.  When supplied, horizontal
#'   error bars depict the split-derived 95\% CI on the effect axis.
#'   Default \code{NULL} (no error bars).
#' @param show_band Logical.  If \code{TRUE} and \code{sg_focus} is
#'   \code{"hrMaxSG"} or \code{"hrMinSG"}, shade the effect-size
#'   neighborhood band on the right of the plot.  Default \code{FALSE}.
#' @param effect_neighborhood Numeric.  Override of the band width;
#'   defaults to the value used in the original \code{forestsearch}
#'   call (or \code{0.10} if not recoverable).
#' @param label_members Logical.  If \code{TRUE}, label each frontier
#'   member with its short cut description.  Default \code{TRUE}.
#' @param point_size Numeric.  Size of frontier-member points.
#'   Default \code{3}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The frontier polyline is drawn as a downward step from large-effect
#' / small-N candidates to small-effect / large-N candidates.  Points
#' off the frontier (if any are present in the \code{fs$grp.consistency$out_sg$result}
#' table beyond those on the frontier) are not currently plotted; this
#' function only displays the frontier and its selected member.
#'
#' For \code{sg_focus = "hrMinSG"}, the selected subgroup may not lie
#' on the Pareto frontier (the focus deliberately prefers small
#' subgroups, which are typically N-dominated).  In that case the
#' selected marker appears off the polyline.
#'
#' @examples
#' \dontrun{
#' p <- plot_pareto_frontier(fs)
#' print(p)
#'
#' # With split-derived 95% CIs:
#' ci_tab <- compute_frontier_cis(fs, n_splits = 1000)
#' p2 <- plot_pareto_frontier(fs, ci_table = ci_tab)
#' print(p2)
#' }
#'
#' @seealso \code{\link{pareto_frontier_table}},
#'   \code{\link{compute_frontier_cis}}, \code{\link{frontier_member_flags}}.
#'
#' @importFrom data.table is.data.table copy setorder
#' @importFrom ggplot2 ggplot aes geom_step geom_point geom_text
#'             geom_errorbarh annotate labs theme_bw theme
#'             scale_colour_manual element_blank
#' @param xlim Optional numeric vector of length 2 controlling the
#'   x-axis range, e.g.\ \code{c(0.5, 3)}.  By default the axis
#'   auto-expands to include split-CI bars when supplied, which can
#'   make individual points look crowded.  Pass an explicit
#'   \code{xlim} (or use \code{xlim_trim = TRUE}) to zoom in.
#' @param xlim_trim Logical.  If \code{TRUE} (and \code{xlim} is
#'   \code{NULL}), the x-axis is trimmed to the range of frontier
#'   point estimates with a small padding, ignoring CI bars.
#'   Default \code{FALSE}.
#'
#' @export
plot_pareto_frontier <- function(fs,
                                 ci_table          = NULL,
                                 show_band         = FALSE,
                                 effect_neighborhood = NULL,
                                 label_members     = TRUE,
                                 point_size        = 3,
                                 xlim              = NULL,
                                 xlim_trim         = FALSE) {

  out_sg   <- tryCatch(fs$grp.consistency$out_sg, error = function(e) NULL)
  frontier <- tryCatch(out_sg$pareto_frontier,    error = function(e) NULL)
  if (is.null(frontier) || !data.table::is.data.table(frontier) ||
      nrow(frontier) == 0L) {
    message("No Pareto frontier available; cannot plot.")
    return(invisible(NULL))
  }

  effect_measure <- fs$effect_measure %||% "HR"
  is_log_scale   <- effect_measure %in% c("OR", "RR", "IRR")
  sg_focus       <- fs$sg_focus %||% "hr"

  ft <- data.table::copy(frontier)
  if (is_log_scale && "hr" %in% names(ft)) {
    ft[["hr"]] <- exp(as.numeric(ft[["hr"]]))
  }

  selected_m <- tryCatch(
    as.integer(out_sg$result[1, ]$m),
    error = function(e) NA_integer_
  )
  ft[["is_selected"]] <-
    !is.na(selected_m) & !is.na(ft$m) & as.integer(ft$m) == selected_m

  # Build short labels for each frontier member from M.* columns
  m_cols <- grep("^M\\.", names(ft), value = TRUE)
  if (length(m_cols) > 0L) {
    ft[["label"]] <- vapply(seq_len(nrow(ft)), function(k) {
      cuts <- unlist(ft[k, m_cols, with = FALSE], use.names = FALSE)
      cuts <- cuts[!is.na(cuts) & nzchar(cuts)]
      if (length(cuts) == 0L) "(empty)" else
        paste(cuts, collapse = " & ")
    }, character(1))
  } else {
    ft[["label"]] <- as.character(ft$m)
  }

  # Step polyline points (sort by hr descending)
  step_dt <- data.table::copy(ft[, c("hr", "N"), with = FALSE])
  data.table::setorder(step_dt, -hr)

  # Optional split CI merge.
  # NB: frontier$m on some fs objects is character (legacy artifact of
  # do.call(rbind, ...) coercion in subgroup.consistency); ci_table$m is
  # always integer.  Coerce both to integer here to avoid a
  # bmerge() join-type-mismatch error.
  if (!is.null(ci_table) && data.table::is.data.table(ci_table) &&
      nrow(ci_table) > 0L) {
    ft_local <- ft
    ft_local[["m"]] <- suppressWarnings(as.integer(ft_local[["m"]]))
    ci_local <- ci_table[, c("m", "split_lcl", "split_ucl"),
                         with = FALSE]
    ci_local[["m"]] <- suppressWarnings(as.integer(ci_local[["m"]]))
    ft <- merge(ft_local, ci_local,
                by = "m", all.x = TRUE, sort = FALSE)
  } else {
    ft[["split_lcl"]] <- NA_real_
    ft[["split_ucl"]] <- NA_real_
  }

  effect_label <- switch(effect_measure,
    HR  = "Hazard ratio",
    OR  = "Odds ratio",
    RR  = "Risk ratio",
    IRR = "Incidence rate ratio",
    RD  = "Risk difference",
    IRD = "Incidence rate difference",
    MD  = "Mean difference",
    effect_measure
  )

  p <- ggplot2::ggplot(ft, ggplot2::aes(x = hr, y = N)) +
    ggplot2::geom_step(data = step_dt,
                       ggplot2::aes(x = hr, y = N),
                       direction = "hv",
                       linewidth = 0.7, colour = "steelblue") +
    ggplot2::geom_point(
      ggplot2::aes(colour = is_selected),
      size = point_size
    ) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "grey25", `TRUE` = "#D55E00"),
      labels = c("frontier", "selected"),
      name   = NULL
    ) +
    ggplot2::labs(x = effect_label,
                  y = expression(italic(N) ~ " (subgroup size)")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position  = "right")

  # Optional split CI error bars
  if (any(!is.na(ft$split_lcl))) {
    p <- p + ggplot2::geom_errorbarh(
      data = ft[!is.na(split_lcl) & !is.na(split_ucl)],
      ggplot2::aes(xmin = split_lcl, xmax = split_ucl, y = N),
      height = 0, alpha = 0.5, colour = "grey40",
      inherit.aes = FALSE
    )
  }

  if (label_members) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label),
      nudge_y = max(ft$N) * 0.03,
      size = 3.2, colour = "grey20"
    )
  }

  # Resolve x-axis range.
  #   - If user supplied xlim, honor it.
  #   - Else if xlim_trim, zoom to the range of point estimates
  #     (ignoring CI bars) plus a small padding.
  #   - Else: let ggplot auto-expand (default).
  if (!is.null(xlim)) {
    if (!is.numeric(xlim) || length(xlim) != 2L) {
      warning("xlim must be a numeric vector of length 2; ignoring.",
              call. = FALSE)
    } else {
      p <- p + ggplot2::coord_cartesian(xlim = xlim)
    }
  } else if (isTRUE(xlim_trim)) {
    rng <- range(ft$hr, na.rm = TRUE)
    pad <- 0.05 * diff(rng)
    p <- p + ggplot2::coord_cartesian(xlim = c(rng[1] - pad, rng[2] + pad))
  }

  # Optional effect-neighborhood band
  if (show_band && sg_focus %in% c("hrMaxSG", "hrMinSG")) {
    eps_val <- effect_neighborhood %||%
      fs$args_call_all$effect_neighborhood %||% 0.10
    hr_max  <- max(ft$hr, na.rm = TRUE)
    floor_v <- (1 - eps_val) * hr_max
    p <- p + ggplot2::annotate("rect",
      xmin = floor_v, xmax = hr_max * 1.02,
      ymin = -Inf,    ymax = Inf,
      fill = "grey85", alpha = 0.35)
  }

  p
}


`%||%` <- function(a, b) if (is.null(a)) b else a
