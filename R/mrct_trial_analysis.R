#' Regional-consistency analysis of an observed multiregional trial
#'
#' The observed-trial counterpart of [mrct_case_study()], running the same analysis
#' code: the all-patient and region Cox fits and the trigger; [forestsearch()] on the
#' non-region patients with the region's patients as the test set; the declared
#' subgroup \eqn{\hat H} applied to the whole trial; the four-set Cox summary; and
#' multiplier-resampling inference on \eqn{\hat H}.
#'
#' @section Stages:
#' The result reports how far the pipeline went in `$stage`; none of these is an error.
#' The all-patient and region fits are always computed (the `itt_all` and `itt_region`
#' rows of `$sets`).
#' - `"not_positive"`: the all-patient upper 95% limit is not below 1.0.
#' - `"no_trigger"`: the trial is positive and the region's hazard ratio is not above
#'   1.0.
#' - `"no_declaration"`: the identifier ran and declared nothing. As in
#'   [mrct_case_study()], `region_Ghat` is then the whole region.
#' - `"complete"`: the identifier ran and declared \eqn{\hat H}.
#'
#' The identifier does not run at the first two stages, so `fs`, `mr`, `in_Hhat` and
#' `confusion` are `NULL` there, `Hhat` is `NA`, and `$sets` has no `region_Ghat` row.
#' With `proceed_without_trigger = TRUE` it runs regardless: the stage is then
#' `"no_declaration"` or `"complete"`, and `$trigger` still reports the trigger.
#'
#' @param df A data frame with one row per patient.
#' @param outcome.name,event.name,treat.name Columns: the time, the event indicator
#'   (1 = event) and the treatment (1 = experimental).
#' @param id.name The patient-id column; [forestsearch()] needs one.
#' @param region.name The observed region column (1 = in the region, 0 = not).
#' @param truth.name Optional column of a known harmed subgroup (1 = in \eqn{H^*}), as
#'   in simulation. It adds the region(\eqn{G^*}) row and the 2 x 2 tables of \eqn{\hat
#'   H} against \eqn{H^*}; without it, as in a real trial, neither appears.
#' @param fs_args Named list of [forestsearch()] arguments, as in [mrct_case_study()].
#'   `confounders.name` is required. `sg_focus`, `maxk`, `hr.threshold`,
#'   `hr.consistency` and `pconsistency.threshold` not given take
#'   [mrct_case_study()]'s defaults, and `conf_force` defaults to `character(0)`. The
#'   data, the column names, the display flags and the (sequential) parallel plan are
#'   set here.
#' @param mr_inference Logical. Run multiplier resampling on \eqn{\hat H}
#'   ([fs_mr_inference()]). Default `TRUE`.
#' @param mr_inference_args List passed to `forestsearch(mr_inference_args = )`; give
#'   `seed` here for a reproducible result. Default: the field method with the
#'   complement de-biased and given its own field.
#' @param proceed_without_trigger Logical. Run the identifier on a trial without the
#'   trigger. Default `FALSE`.
#'
#' @return An object of class `mrct_trial_analysis`, a list with
#'   - `stage`: see the Stages section.
#'   - `trial`: the trial (data.frame).
#'   - `fs`: the [forestsearch()] result; [mr_values()] or [mr_estimates_table()] read
#'     its multiplier-resampling results.
#'   - `sets`: the set table, as in [mrct_four_sets()]. The unadjusted Cox fits are its
#'     rows.
#'   - `in_Hhat`: logical, per trial row: membership of \eqn{\hat H}.
#'   - `confusion`: with `truth.name`, a list of two integer vectors
#'     `c(TP, FP, FN, TN)` of \eqn{\hat H} against \eqn{H^*}, `nonregion` and `region`.
#'   - `Hhat`: the definition (`"none"` when nothing is declared).
#'   - `mr`: `fs$mr_inference`.
#'   - `trigger`: logical: all-patient upper limit < 1.0 and region hazard ratio > 1.0.
#'   - `call_args`: the column names and the [forestsearch()] arguments (without data).
#'   - `warnings`: warnings raised by [forestsearch()], muffled and kept here.
#'
#' @seealso [mrct_four_sets()] for the sets alone; [mrct_case_study()] for a simulated
#'   trial; [mr_values()] for the multiplier-resampling quantities.
#'
#' @examples
#' set.seed(1)
#' n <- 800
#' d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
#' # a harmed subgroup in every region, and more of it in the region
#' d$region <- rbinom(n, 1, plogis(-1.5 - 1.5 * d$x1))
#' d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
#' d$event <- as.integer(d$time < 30)
#' d$time <- pmin(d$time, 30)
#' fit <- mrct_trial_analysis(
#'   d, "time", "event", "treat", "id", region.name = "region",
#'   fs_args = list(sg_focus = "hrMaxSG", confounders.name = c("x1", "x2"),
#'                  use_grf = FALSE, use_lasso = FALSE, fs.splits = 200),
#'   mr_inference = FALSE)
#' fit$stage
#' fit
#' plot(fit, which = "forest")
#'
#' @export
mrct_trial_analysis <- function(df, outcome.name, event.name, treat.name, id.name, region.name,
                                truth.name = NULL, fs_args = list(),
                                mr_inference = TRUE,
                                mr_inference_args = list(ci_method = "field",
                                                         include_complement = TRUE,
                                                         field_complement = TRUE),
                                proceed_without_trigger = FALSE) {
  .fs_preserve_rng()
  .mrct_check_trial(df, outcome.name, event.name, treat.name, region.name, truth.name, id.name)
  if (!is.list(fs_args) || !is.list(mr_inference_args))
    stop("fs_args and mr_inference_args must be lists.", call. = FALSE)
  if (is.null(fs_args$confounders.name))
    stop("fs_args$confounders.name is required: the covariates the identifier may use.",
         call. = FALSE)

  dfs <- data.table::as.data.table(df)
  trial <- as.data.frame(dfs)
  nm <- list(outcome = outcome.name, treat = treat.name, event = event.name, id = id.name)
  columns <- list(outcome.name = outcome.name, event.name = event.name, treat.name = treat.name,
                  region.name = region.name, id.name = id.name, truth.name = truth.name)

  fits <- .mrct_four_set_table(trial, region.name, truth.name, NULL, nm)
  trigger <- .mrct_trigger(fits$sets)
  if (!trigger && !isTRUE(proceed_without_trigger)) {
    hi <- fits$sets$hi[fits$sets$set == "itt_all"]
    positive <- !is.na(hi) && hi < 1
    return(structure(list(
      stage = if (positive) "no_trigger" else "not_positive",
      trial = trial, fs = NULL, sets = fits$sets, in_Hhat = NULL, confusion = NULL,
      Hhat = NA_character_, mr = NULL, trigger = trigger,
      call_args = c(columns, list(fs_args = fs_args)), warnings = character(0)),
      class = "mrct_trial_analysis"))
  }

  defaults <- list(sg_focus = "minSG", maxk = 1, hr.threshold = 0.90, hr.consistency = 0.80,
                   pconsistency.threshold = 0.90, conf_force = character(0))
  core <- .mrct_analyze_trial(dfs, region.name, truth.name, utils::modifyList(defaults, fs_args),
                              mr_inference, mr_inference_args, nm, caller = "mrct_trial_analysis")
  structure(list(
    stage = if (identical(core$Hhat, "none")) "no_declaration" else "complete",
    trial = core$trial, fs = core$fs, sets = core$sets, in_Hhat = core$in_Hhat,
    confusion = core$confusion, Hhat = core$Hhat, mr = core$mr, trigger = core$trigger,
    call_args = c(columns, list(fs_args = core$fs_call_args[setdiff(names(core$fs_call_args),
                                                                    c("df.analysis", "df.test"))])),
    warnings = core$warnings),
    class = "mrct_trial_analysis")
}


#' @rdname mrct_trial_analysis
#' @param x An `mrct_trial_analysis` object.
#' @param digits Digits for hazard ratios. Default 3.
#' @param ... Unused.
#' @export
print.mrct_trial_analysis <- function(x, digits = 3, ...) {
  fm <- .mrct_fm(digits)
  a <- x$call_args
  cat(sprintf("MRCT trial analysis: n = %d, region '%s'\n", nrow(x$trial), a$region.name))
  ran <- !is.null(x$fs)
  cat(sprintf("Stage: %s%s\n", x$stage, if (ran) "" else " (the identifier was not run)"))
  cat(sprintf("H-hat: %s\n", if (ran) x$Hhat else "not evaluated"))
  cat(sprintf("Trigger (ITT ub < 1.0 and region ITT > 1.0): %s\n\n", x$trigger))
  .mrct_print_sets(x$sets, fm)
  .mrct_print_confusion(x$confusion)
  .mrct_print_mr(x$mr, fm, evaluated = TRUE)
  invisible(x)
}


#' @rdname mrct_trial_analysis
#' @param which Any of `"forest"`, `"km"`. Default both.
#' @param xlim Hazard-ratio range of the forest plot. Default `c(0.3, 2)`; intervals are
#'   clipped to it.
#' @param time_unit Character. The label of the Kaplan-Meier time axis; `NULL` (the
#'   default) labels it `"time"`.
#' @param analysis_time The Kaplan-Meier axis's upper limit. `NULL` (the default) uses
#'   the longest follow-up in the trial.
#' @export
plot.mrct_trial_analysis <- function(x, which = c("forest", "km"), xlim = c(0.3, 2),
                                     time_unit = NULL, analysis_time = NULL, ...) {
  which <- match.arg(which, several.ok = TRUE)
  a <- x$call_args
  .mrct_plot_sets(x$sets, x$trial, a$region.name, a$truth.name, x$in_Hhat, which, xlim, time_unit,
                  if (is.null(analysis_time)) max(x$trial[[a$outcome.name]]) else analysis_time,
                  list(outcome = a$outcome.name, treat = a$treat.name, event = a$event.name),
                  "plot.mrct_trial_analysis")
  invisible(x)
}
