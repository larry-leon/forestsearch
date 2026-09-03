# compare_subgroup_sims.R -------------------------------------------------
# Two-design comparison for extreme-subgroups studies: the single
# package-level implementation of the design-comparison memo's former
# inline machinery (quarto/extreme_subgroups/fixed_random/
# extreme_subgroups_design_comparison.qmd).  The per-design statistics
# (.col_stats) and the interleaved comparison frame are lifted verbatim
# from that document, so live-mode memo renders are column-identical by
# construction (verified by dev/accept_phase3_memo_compare.R).
#
# Phase 4.2 (effect-aware summaries): tail thresholds resolve explicit
# arguments -> the inputs' `effect` metadata -> the HR legacy values;
# with legacy inputs and no overrides the computed frame and attributes
# are identical to the pre-4.2 output (same threshold values flow
# through the same expressions).  Guard order (4.2r2): metadata
# compatibility runs before panel alignment, so the categorical
# incompatibility (MD vs HR) is reported first.

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
#' Column names are structural and retained across outcome types
#' (`hr05`, `hr1`, `ub2`, `ub3`, `mHR`, `mUB`), exactly as `sim_hrs`
#' serves generic duty on GLM results: for [subgroup_glm()] fits they
#' hold estimate-scale statistics at the resolved thresholds. Threshold
#' resolution: explicit arguments win; otherwise the inputs' `effect`
#' metadata supplies them; otherwise the HR legacy `c(0.5, 1.0)` /
#' `c(2, 3)`. An `NA` threshold yields an all-`NA` column. The two
#' inputs must carry compatible metadata (both none, or the same
#' effect measure) -- comparing an MD study against an HR study is a
#' hard error, raised before the panel-alignment guard so the
#' categorical incompatibility is reported first.
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
#' @param est_thresholds Length-2 numeric `c(low, high)` for the
#'   `hr05`-style (`est < low`) and `hr1`-style (`est > high`)
#'   percentages; `NULL` resolves via `effect` metadata, then the HR
#'   legacy `c(0.5, 1.0)`.
#' @param ub_thresholds Length-2 numeric for the `ub2`-style and
#'   `ub3`-style (`UB >= t`) percentages; resolution as above (HR
#'   legacy `c(2, 3)`).
#'
#' @return A plain `data.frame` (one row per subgroup, matching the
#'   memo's static-mode frame type) with columns `subgroup`, then for
#'   each of `N`, `ub2`, `ub3`, `mUB`, `hr05`, `hr1`, `mHR` the `x` and
#'   `y` columns interleaved. The inputs' `n_sims` and `design` values
#'   are attached as attributes `"n_sims"` and `"designs"`; when the
#'   inputs carry effect metadata it is attached as attribute
#'   `"effect"` (informational; attributes are dropped by most
#'   transformations).
#' @seealso [run_subgroup_sims()], [summary.subgroup_sims()],
#'   [subgroup_glm()]
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
                                  suffixes = c("_r", "_f"),
                                  est_thresholds = NULL,
                                  ub_thresholds = NULL) {
  for (o in list(x, y)) {
    if (!is.list(o) ||
        !all(c("sim_hrs", "sim_ubs", "sim_ns") %in% names(o)) ||
        is.null(colnames(o$sim_hrs))) {
      stop("each input must carry sim_hrs / sim_ubs / sim_ns matrices ",
           "with subgroup column names (a run_subgroup_sims() result ",
           "or a vignette payload).")
    }
  }
  # Effect-metadata compatibility guard -- runs BEFORE the panel guard:
  # comparing an MD study against an HR study is categorically wrong,
  # and re-running with aligned subgroups (the panel guard's advice)
  # could not fix it. Legacy inputs (both NULL) pass trivially.
  ex <- x$effect; ey <- y$effect
  if (xor(is.null(ex), is.null(ey)) ||
      (!is.null(ex) && !identical(ex$measure, ey$measure))) {
    stop("The two results carry incompatible effect metadata (",
         if (is.null(ex)) "none" else ex$measure, " vs ",
         if (is.null(ey)) "none" else ey$measure,
         "); compare like with like.")
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

  # Threshold resolution (explicit -> metadata -> HR legacy); ex/ey were
  # validated by the metadata guard above.
  if (is.null(est_thresholds)) {
    est_thresholds <- if (!is.null(ex$est_thresholds)) {
      ex$est_thresholds
    } else c(0.5, 1.0)
  }
  if (is.null(ub_thresholds)) {
    ub_thresholds <- if (!is.null(ex$ub_thresholds)) {
      ex$ub_thresholds
    } else c(2, 3)
  }
  stopifnot(is.numeric(est_thresholds), length(est_thresholds) == 2L,
            is.numeric(ub_thresholds),  length(ub_thresholds)  == 2L)

  # Per-design statistics -- verbatim from the memo's .col_stats(), with
  # the thresholds flowing in as values (legacy values reproduce the
  # memo's numbers identically; an NA threshold propagates to an all-NA
  # column via the NaN -> NA mapping).
  nan_to_na <- function(v) { v[is.nan(v)] <- NA_real_; v }
  .cs <- function(pl) {
    data.frame(
      subgroup = colnames(pl$sim_hrs),
      N    = nan_to_na(round(colMeans(pl$sim_ns, na.rm = TRUE))),
      ub2  = nan_to_na(100 * apply(pl$sim_ubs, 2,
                                   function(v) mean(v >= ub_thresholds[1],
                                                    na.rm = TRUE))),
      ub3  = nan_to_na(100 * apply(pl$sim_ubs, 2,
                                   function(v) mean(v >= ub_thresholds[2],
                                                    na.rm = TRUE))),
      mUB  = suppressWarnings(apply(pl$sim_ubs, 2, stats::median,
                                    na.rm = TRUE)),
      hr05 = nan_to_na(100 * apply(pl$sim_hrs, 2,
                                   function(v) mean(v < est_thresholds[1],
                                                    na.rm = TRUE))),
      hr1  = nan_to_na(100 * apply(pl$sim_hrs, 2,
                                   function(v) mean(v > est_thresholds[2],
                                                    na.rm = TRUE))),
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
  if (!is.null(ex)) attr(out, "effect") <- ex
  out
}
