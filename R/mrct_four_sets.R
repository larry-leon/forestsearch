#' The four patient sets of a regional-consistency analysis
#'
#' Fits the unadjusted two-sample Cox model on each set of the regional-consistency
#' analysis, from columns you name: all patients, the region, region(\eqn{G^*}) (the
#' region without the true harmed subgroup, when `truth.name` is given) and
#' region(\eqn{\hat G}) (the region without the subgroup \eqn{\hat H}). The sets and
#' their labels are those of [mrct_case_study()] and [mrct_trial_analysis()], and
#' `print()` and `plot()` draw on the same code.
#'
#' @param df A data frame with one row per patient.
#' @param outcome.name,event.name,treat.name Columns: the time, the event indicator
#'   (1 = event) and the treatment (1 = experimental).
#' @param region.name The region column (1 = in the region, 0 = not).
#' @param in_Hhat Logical, one value per row: membership of \eqn{\hat H}, for instance
#'   `get_dfpred(df, sg.harm = fit$sg.harm, version = 2)$treat.recommend == 0`. `NULL`
#'   (the default) means nothing is declared: region(\eqn{\hat G}) is the whole region.
#' @param truth.name Optional column of the true harmed subgroup (1 = in \eqn{H^*}).
#'   It adds the region(\eqn{G^*}) row and the 2 x 2 tables of \eqn{\hat H} against
#'   \eqn{H^*}; without it neither appears.
#'
#' @return An object of class `mrct_four_sets`, a list with
#'   - `sets`: the table, also returned by `as.data.frame()`: `set`, `label`, `n`,
#'     `n_exp`, `n_ctl`, `events_exp`, `events_ctl`, and the Cox fit on the set: `hr`,
#'     `lo`, `hi` (Wald 95%), `se_loghr`. The rows are `itt_all`, `itt_region`,
#'     `region_Gstar` (with `truth.name`) and `region_Ghat`.
#'   - `confusion`: with `truth.name`, a list of two integer vectors
#'     `c(TP, FP, FN, TN)` of \eqn{\hat H} against \eqn{H^*}, `nonregion` and `region`;
#'     else `NULL`.
#'   - `trial`, `in_Hhat`, and `call_args`, the column names.
#'
#' @seealso [mrct_trial_analysis()] for the whole pipeline on an observed trial;
#'   [mrct_case_study()] for a simulated trial.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- data.frame(time = rexp(n, 0.1), event = rbinom(n, 1, 0.8),
#'                 treat = rbinom(n, 1, 0.5), region = rbinom(n, 1, 0.25),
#'                 x = rbinom(n, 1, 0.3))
#' fs4 <- mrct_four_sets(d, "time", "event", "treat", region.name = "region",
#'                       in_Hhat = d$x == 1)
#' fs4
#' as.data.frame(fs4)[, c("set", "n", "hr", "lo", "hi")]
#' plot(fs4, which = "forest")
#'
#' @export
mrct_four_sets <- function(df, outcome.name, event.name, treat.name, region.name,
                           in_Hhat = NULL, truth.name = NULL) {
  .mrct_check_trial(df, outcome.name, event.name, treat.name, region.name, truth.name)
  if (!is.null(in_Hhat) &&
      (!is.logical(in_Hhat) || length(in_Hhat) != nrow(df) || anyNA(in_Hhat)))
    stop("in_Hhat must be NULL or a logical vector with one non-missing value per row of df.",
         call. = FALSE)
  trial <- as.data.frame(df)
  if (is.null(in_Hhat)) in_Hhat <- rep(FALSE, nrow(trial))
  nm <- list(outcome = outcome.name, treat = treat.name, event = event.name)
  built <- .mrct_four_set_table(trial, region.name, truth.name, in_Hhat, nm)
  structure(list(
    sets = built$sets, confusion = built$confusion, trial = trial, in_Hhat = in_Hhat,
    call_args = list(outcome.name = outcome.name, event.name = event.name,
                     treat.name = treat.name, region.name = region.name,
                     truth.name = truth.name)),
    class = "mrct_four_sets")
}


#' @rdname mrct_four_sets
#' @param x An `mrct_four_sets` object.
#' @param digits Digits for hazard ratios. Default 3.
#' @param ... Unused.
#' @export
print.mrct_four_sets <- function(x, digits = 3, ...) {
  a <- x$call_args
  cat(sprintf("MRCT four sets: n = %d, region '%s'\n\n", nrow(x$trial), a$region.name))
  .mrct_print_sets(x$sets, .mrct_fm(digits))
  .mrct_print_confusion(x$confusion)
  invisible(x)
}


#' @rdname mrct_four_sets
#' @param row.names,optional Unused, as in `as.data.frame()`.
#' @export
as.data.frame.mrct_four_sets <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$sets
}


#' @rdname mrct_four_sets
#' @param which Any of `"forest"`, `"km"`. Default both.
#' @param xlim Hazard-ratio range of the forest plot. Default `c(0.3, 2)`; intervals are
#'   clipped to it.
#' @param time_unit Character. The label of the Kaplan-Meier time axis; `NULL` (the
#'   default) labels it `"time"`.
#' @param analysis_time The Kaplan-Meier axis's upper limit. `NULL` (the default) uses
#'   the longest follow-up in the trial.
#' @export
plot.mrct_four_sets <- function(x, which = c("forest", "km"), xlim = c(0.3, 2),
                                time_unit = NULL, analysis_time = NULL, ...) {
  which <- match.arg(which, several.ok = TRUE)
  a <- x$call_args
  .mrct_plot_sets(x$sets, x$trial, a$region.name, a$truth.name, x$in_Hhat, which, xlim, time_unit,
                  if (is.null(analysis_time)) max(x$trial[[a$outcome.name]]) else analysis_time,
                  list(outcome = a$outcome.name, treat = a$treat.name, event = a$event.name),
                  "plot.mrct_four_sets")
  invisible(x)
}


#' Check the columns of an observed trial
#'
#' The analysis columns must exist and be complete; the treatment, region and truth
#' columns must be 0/1.
#' @noRd
.mrct_check_trial <- function(df, outcome.name, event.name, treat.name, region.name,
                              truth.name = NULL, id.name = NULL) {
  if (!is.data.frame(df)) stop("df must be a data frame.", call. = FALSE)
  args <- list(outcome.name = outcome.name, event.name = event.name, treat.name = treat.name,
               region.name = region.name, truth.name = truth.name, id.name = id.name)
  args <- args[!vapply(args, is.null, NA)]
  for (a in names(args)) {
    v <- args[[a]]
    if (!is.character(v) || length(v) != 1L)
      stop(a, " must be one column name.", call. = FALSE)
    if (!v %in% names(df)) stop(a, " '", v, "' not found in df.", call. = FALSE)
    if (a != "id.name" && anyNA(df[[v]]))
      stop(a, " '", v, "' has missing values.", call. = FALSE)
  }
  for (a in intersect(names(args), c("treat.name", "region.name", "truth.name")))
    if (!all(df[[args[[a]]]] %in% c(0, 1)))
      stop(a, " '", args[[a]], "' must be 0/1.", call. = FALSE)
  invisible(TRUE)
}
