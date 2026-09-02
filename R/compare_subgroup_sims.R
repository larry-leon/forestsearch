# compare_subgroup_sims.R -------------------------------------------------
# Two-design comparison for extreme-subgroups studies: the single
# package-level implementation of the design-comparison memo's former
# inline machinery (quarto/extreme_subgroups/fixed_random/
# extreme_subgroups_design_comparison.qmd).  The per-design statistics
# (.col_stats) and the interleaved comparison frame are lifted verbatim
# from that document, so live-mode memo renders are column-identical by
# construction (verified by dev/accept_phase3_memo_compare.R).

#' Compare two extreme-subgroups simulation studies
#'
#' Builds the per-subgroup comparison frame between two
#' [run_subgroup_sims()] results (or RDS payloads loaded from the
#' vignettes -- any list carrying `sim_hrs` / `sim_ubs` / `sim_ns`
#' matrices with identical column names). For each input the per-design
#' statistics replicate the vignettes' Section 6.6.1 definitions
#' exactly: tail probabilities are percentages denominated over
#' converged fits (`na.rm = TRUE`), medians are unconditional, `N` is
#' the across-trial mean subgroup size, and structurally empty
#' subgroups (all-`NA` columns) are mapped from `NaN` to `NA`.
#'
#' Two validations guard the row-wise alignment (formerly the memo's
#' Guard 1 and Guard 2): the subgroup panels must be identical, and,
#' when `expect_designs` is supplied, each input's `design` label must
#' match it exactly.
#'
#' The medians here use `median()`, matching the memo; they are equal
#' to [summary.subgroup_sims()]'s type-7 `quantile(..., 0.50)` values
#' up to floating-point formulation (`(a + b) / 2` versus
#' `a + 0.5 * (b - a)`), which is why the memo's numbers and the
#' vignettes' tables have always agreed to the printed digit.
#'
#' @param x,y The two studies to compare, in display order (the memo
#'   passes random-X then fixed-X).
#' @param expect_designs Optional length-2 character: required `design`
#'   labels for `x` and `y` respectively (e.g.
#'   `c("resample", "fixed")`); `NULL` skips the check.
#' @param suffixes Length-2 character appended to each statistic's
#'   column for `x` and `y`; the default `c("_r", "_f")` reproduces the
#'   memo's column names (`N_r`, `N_f`, `ub2_r`, ...).
#'
#' @return A plain `data.frame` (one row per subgroup, matching the
#'   memo's static-mode frame type) with columns `subgroup`, then for
#'   each of `N`, `ub2`, `ub3`, `mUB`, `hr05`, `hr1`, `mHR` the `x` and
#'   `y` columns interleaved. The inputs' `n_sims` and `design` values
#'   are attached as attributes `"n_sims"` and `"designs"`
#'   (informational; dropped by most transformations).
#' @seealso [run_subgroup_sims()], [summary.subgroup_sims()]
#' @export
#' @examples
#' \dontrun{
#' pl_r <- readRDS("results/extreme_sims_resample_10000_payload.rds")
#' pl_f <- readRDS("results/extreme_sims_fixed_10000_payload.rds")
#' cmp <- compare_subgroup_sims(pl_r, pl_f,
#'                              expect_designs = c("resample", "fixed"))
#' }
compare_subgroup_sims <- function(x, y,
                                  expect_designs = NULL,
                                  suffixes = c("_r", "_f")) {
  for (o in list(x, y)) {
    if (!is.list(o) ||
        !all(c("sim_hrs", "sim_ubs", "sim_ns") %in% names(o)) ||
        is.null(colnames(o$sim_hrs))) {
      stop("each input must carry sim_hrs / sim_ubs / sim_ns matrices ",
           "with subgroup column names (a run_subgroup_sims() result ",
           "or a vignette payload).")
    }
  }
  if (!identical(colnames(x$sim_hrs), colnames(y$sim_hrs))) {
    stop("Subgroup names differ between the two result sets; re-run ",
         "both studies from the same subgroup definitions.")
  }
  if (!is.null(expect_designs)) {
    stopifnot(is.character(expect_designs), length(expect_designs) == 2L)
    if (!identical(x$design, expect_designs[1]) ||
        !identical(y$design, expect_designs[2])) {
      stop("Design labels are not ('", expect_designs[1], "', '",
           expect_designs[2], "'): got ('", x$design, "', '", y$design,
           "').  Check the inputs.")
    }
  }
  stopifnot(is.character(suffixes), length(suffixes) == 2L,
            !anyDuplicated(suffixes))

  # Per-design statistics -- verbatim from the memo's .col_stats().
  nan_to_na <- function(v) { v[is.nan(v)] <- NA_real_; v }
  .cs <- function(pl) {
    data.frame(
      subgroup = colnames(pl$sim_hrs),
      N    = nan_to_na(round(colMeans(pl$sim_ns, na.rm = TRUE))),
      ub2  = nan_to_na(100 * apply(pl$sim_ubs, 2,
                                   function(v) mean(v >= 2.0, na.rm = TRUE))),
      ub3  = nan_to_na(100 * apply(pl$sim_ubs, 2,
                                   function(v) mean(v >= 3.0, na.rm = TRUE))),
      mUB  = suppressWarnings(apply(pl$sim_ubs, 2, stats::median,
                                    na.rm = TRUE)),
      hr05 = nan_to_na(100 * apply(pl$sim_hrs, 2,
                                   function(v) mean(v < 0.5, na.rm = TRUE))),
      hr1  = nan_to_na(100 * apply(pl$sim_hrs, 2,
                                   function(v) mean(v > 1.0, na.rm = TRUE))),
      mHR  = suppressWarnings(apply(pl$sim_hrs, 2, stats::median,
                                    na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
  a <- .cs(x)
  b <- .cs(y)

  out <- data.frame(subgroup = a$subgroup, stringsAsFactors = FALSE)
  for (nm in c("N", "ub2", "ub3", "mUB", "hr05", "hr1", "mHR")) {
    out[[paste0(nm, suffixes[1])]] <- a[[nm]]
    out[[paste0(nm, suffixes[2])]] <- b[[nm]]
  }
  attr(out, "n_sims")  <- c(x = if (is.null(x$n_sims)) NA_integer_ else
                                  as.integer(x$n_sims),
                            y = if (is.null(y$n_sims)) NA_integer_ else
                                  as.integer(y$n_sims))
  attr(out, "designs") <- c(x = if (is.null(x$design)) NA_character_ else
                                  x$design,
                            y = if (is.null(y$design)) NA_character_ else
                                  y$design)
  out
}
