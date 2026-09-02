# summary_subgroup_sims.R -------------------------------------------------
# Single source of truth for the extreme-subgroups summary statistics.
#
# Definitions are lifted verbatim from the committed vignettes (commit
# 002f4f37): chunk `uniform-results-table` (quantiles, tail
# probabilities, convergence counts, the formatted results table) plus
# the panel machinery of `forest-theme-params` and
# `ub-forest-highrisk-setup` (category mapping, colour map, ok masks,
# single/combo index vectors, high-risk inclusion and ordering).  The
# design-comparison memo's per-design statistics are these same
# definitions; consolidating them here removes the three-way duplication
# that previously had to be kept in sync by hand.

#' Summarize an extreme-subgroups simulation study
#'
#' Computes, from the raw matrices of a [run_subgroup_sims()] result (or
#' a loaded RDS payload coerced with `class<-`), every statistic used by
#' the vignettes' results table and forest panels: 1st/50th/99th
#' percentile matrices for HR and UB(HR), convergence counts, mean
#' subgroup N, tail probabilities `Pr(HR<0.5)`, `Pr(HR>1.0)`,
#' `Pr(UB>=2)`, `Pr(UB>=3)` (denominated over converged fits,
#' `na.rm = TRUE`), unconditional median display strings, the formatted
#' `results_tbl`, category labels and colours, validity masks, the
#' single/combination panel index vectors, and the high-risk panel
#' (threshold filter, ITT anchor, ordering).
#'
#' Structurally empty subgroups (all-`NA` columns) yield `NA` medians
#' and `NaN` means/proportions, exactly as in the vignettes; the
#' formatted table renders them as `"-"`.
#'
#' @param object A `"subgroup_sims"` object.
#' @param hr_true True hazard ratio for reference lines; defaults to the
#'   value stored on `object`.
#' @param probs Quantile probabilities for the ECI matrices; the default
#'   `c(0.01, 0.50, 0.99)` matches the vignettes and downstream code
#'   indexes column 2 as the median.
#' @param single_cats,combo_cats Category labels routed to the
#'   single-variable and combination forest panels.
#' @param pr_thresh_highrisk High-risk panel inclusion threshold on
#'   `Pr(UB>=2)`, default `0.10`.
#' @param itt_name Name of the ITT row, always included in and anchored
#'   at the top of the high-risk panel.
#' @param cat_cols Named colour map for categories.
#' @param ... Unused.
#'
#' @return An object of class `"subgroup_sims_summary"`: a list with
#'   fields named as in the vignettes (`hr_q`, `ub_q`, `n_valid`,
#'   `pct_valid`, `mean_n`, `pr_hr_lt050`, `pr_hr_gt1`, `pr_ub_ge2`,
#'   `pr_ub_ge3`, `mhr_uncond`, `mub_uncond`, `results_tbl`,
#'   `sg_names`, `cat`, `cat_cols`, `ok`, `ok_ub`, `idx_hr_single`,
#'   `idx_hr_combo`, `idx_ub_single`, `idx_ub_combo`, `n_single`,
#'   `n_combo`, `highrisk`, `n_sims`, `design`, `hr_true`).
#' @seealso [run_subgroup_sims()]
#' @export
summary.subgroup_sims <- function(object,
                                  hr_true = object$hr_true,
                                  probs = c(0.01, 0.50, 0.99),
                                  single_cats = c("ITT", "Clinical",
                                                  "Continuous"),
                                  combo_cats = c("Interaction", "3-way",
                                                 "Random"),
                                  pr_thresh_highrisk = 0.10,
                                  itt_name = "All Patients",
                                  cat_cols = c(
                                    "ITT"         = "black",
                                    "Clinical"    = "steelblue",
                                    "Continuous"  = "darkgreen",
                                    "Interaction" = "purple4",
                                    "3-way"       = "violetred3",
                                    "Random"      = "darkorange"),
                                  ...) {

  sim_hrs <- object$sim_hrs
  sim_ubs <- object$sim_ubs
  sim_ns  <- object$sim_ns
  n_sims  <- object$n_sims
  sg_names <- colnames(sim_hrs)

  # Median and 1st-99th percentile ECI across simulations
  hr_q <- t(apply(sim_hrs, 2, stats::quantile, probs = probs, na.rm = TRUE))
  ub_q <- t(apply(sim_ubs, 2, stats::quantile, probs = probs, na.rm = TRUE))

  # Convergence: number of simulations where estimation succeeded
  n_valid   <- apply(sim_hrs, 2, function(x) sum(!is.na(x)))
  pct_valid <- round(100 * n_valid / n_sims, 1)

  mean_n <- round(apply(sim_ns, 2, mean, na.rm = TRUE))

  # Tail probabilities, denominated over converged fits (na.rm = TRUE)
  pr_hr_lt050 <- apply(sim_hrs, 2, function(x) mean(x < 0.5,  na.rm = TRUE))
  pr_hr_gt1   <- apply(sim_hrs, 2, function(x) mean(x > 1.0,  na.rm = TRUE))
  pr_ub_ge2   <- apply(sim_ubs, 2, function(x) mean(x >= 2.0, na.rm = TRUE))
  pr_ub_ge3   <- apply(sim_ubs, 2, function(x) mean(x >= 3.0, na.rm = TRUE))

  # Unconditional median displays (character, "-" when empty)
  mhr_uncond <- ifelse(is.na(hr_q[, 2]), "-",
                       as.character(round(hr_q[, 2], 2)))
  mub_uncond <- ifelse(is.na(ub_q[, 2]), "-",
                       as.character(round(ub_q[, 2], 2)))

  fmt_pct <- function(x) ifelse(is.nan(x), "-",
                                paste0(round(100 * x, 1), "%"))

  results_tbl <- data.frame(
    Subgroup    = sg_names,
    Mean_N      = ifelse(is.nan(mean_n), "-", as.character(mean_n)),
    N_valid     = paste0(n_valid, " (", pct_valid, "%)"),
    Pr_HR_lt050 = fmt_pct(pr_hr_lt050),
    Pr_HR_gt1   = fmt_pct(pr_hr_gt1),
    mHR         = mhr_uncond,
    Pr_UB_ge2   = fmt_pct(pr_ub_ge2),
    Pr_UB_ge3   = fmt_pct(pr_ub_ge3),
    mUB         = mub_uncond,
    stringsAsFactors = FALSE
  )
  colnames(results_tbl) <- c(
    "Subgroup", "N", "#(% converged)",
    "Pr(HR<0.5)", "Pr(HR>1.0)", "mHR",
    "Pr(UB>=2)", "Pr(UB>=3)", "mUB"
  )

  # Category mapping -- first-match order identical to the vignettes'
  # dplyr::case_when block, in base R to avoid the dependency.
  grp <- vapply(object$subgroups, `[[`, character(1L), "grp")
  cat_lab <- ifelse(grp == "ITT", "ITT",
             ifelse(grp == "Clinical", "Clinical",
             ifelse(grp == "Continuous", "Continuous",
             ifelse(grepl("^Interaction", grp), "Interaction",
             ifelse(grp == "3-way", "3-way", "Random")))))

  ok    <- !is.na(hr_q[, 2])
  ok_ub <- !is.na(ub_q[, 2])

  idx_hr_single <- which(cat_lab[ok]    %in% single_cats)
  idx_hr_combo  <- which(cat_lab[ok]    %in% combo_cats)
  idx_ub_single <- which(cat_lab[ok_ub] %in% single_cats)
  idx_ub_combo  <- which(cat_lab[ok_ub] %in% combo_cats)

  # High-risk panel: meets the Pr(UB>=2) threshold OR is the ITT anchor;
  # ordered ITT first, then descending Pr(UB>=2).
  include_highrisk <- ok_ub &
    (pr_ub_ge2 >= pr_thresh_highrisk | sg_names == itt_name)
  ord_highrisk <- order(
    sg_names[include_highrisk] != itt_name,
    -pr_ub_ge2[include_highrisk]
  )

  out <- list(
    hr_q = hr_q, ub_q = ub_q,
    n_valid = n_valid, pct_valid = pct_valid, mean_n = mean_n,
    pr_hr_lt050 = pr_hr_lt050, pr_hr_gt1 = pr_hr_gt1,
    pr_ub_ge2 = pr_ub_ge2, pr_ub_ge3 = pr_ub_ge3,
    mhr_uncond = mhr_uncond, mub_uncond = mub_uncond,
    results_tbl = results_tbl,
    sg_names = sg_names, cat = cat_lab, cat_cols = cat_cols,
    ok = ok, ok_ub = ok_ub,
    idx_hr_single = idx_hr_single, idx_hr_combo = idx_hr_combo,
    idx_ub_single = idx_ub_single, idx_ub_combo = idx_ub_combo,
    n_single = length(idx_hr_single), n_combo = length(idx_hr_combo),
    highrisk = list(threshold = pr_thresh_highrisk, itt_name = itt_name,
                    include = include_highrisk, ord = ord_highrisk,
                    n = sum(include_highrisk)),
    single_cats = single_cats, combo_cats = combo_cats,
    n_sims = n_sims, design = object$design, hr_true = hr_true
  )
  class(out) <- c("subgroup_sims_summary", "list")
  out
}

#' @export
print.subgroup_sims_summary <- function(x, ...) {
  cat(sprintf(
    "<subgroup_sims_summary>  design: %s | n_sims: %s | subgroups: %d%s\n\n",
    x$design, format(x$n_sims, big.mark = ","), length(x$sg_names),
    if (!is.null(x$hr_true)) sprintf(" | true HR: %.2f", x$hr_true) else ""))
  print(x$results_tbl, row.names = FALSE)
  invisible(x)
}
