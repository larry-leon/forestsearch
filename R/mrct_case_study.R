# =============================================================================
# MRCT case study: find a simulated trial and examine it
# =============================================================================
#
#   - mrct_select_case(): the trials of a mrct_region_sims(record = "full")
#     result meeting the trigger, a declaration and region(G-hat) < 1.0
#   - mrct_case_study(): redraw one trial exactly as mrct_region_sims() drew
#     it, run the identifier (with multiplier resampling) on its non-region
#     patients and analyse the four patient sets
#   - print() / plot() methods for mrct_case_study
#
# The analysis and the print / plot workhorses are in R/mrct_trial_helpers.R,
# shared with mrct_trial_analysis() and mrct_four_sets().
#
# =============================================================================


#' Select a case-study trial from an MRCT simulation
#'
#' Reads a \code{\link{mrct_region_sims}} result run with
#' \code{record = "full"} and returns the trials that meet the conditions of
#' a regional-consistency case study: the \strong{trigger} (the upper 95\%
#' limit of the all-patient hazard ratio below \code{itt_ub} and the region's
#' hazard ratio above \code{region_itt}), a \strong{declared} subgroup
#' \eqn{\hat H}, and a Cox hazard ratio below 1.0 in \strong{region(\eqn{\hat
#' G})}, the region without \eqn{\hat H}.
#'
#' The flags are evaluated per trial from the \code{"estimators"} attribute:
#' the trigger from the \code{itt_all} row's \code{hi} and the
#' \code{itt_region} row's \code{hr}; \code{declared} from the trial's
#' \code{any_found == 1}; region(\eqn{\hat G}) < 1.0 from the
#' \code{region_Ghat} row's \code{hr}. A missing hazard ratio never satisfies
#' a condition. With \code{which = "first"} the result is the first trial by
#' index meeting every condition in \code{conditions}: with the three default
#' conditions, the rule that chose the case example of the mrct-consistency
#' report (replicate 5 of its S2 cell).
#'
#' @param res A \code{\link{mrct_region_sims}} result with
#'   \code{record = "full"} (it carries the \code{"estimators"} attribute;
#'   take the attribute's owner, not a subset of it).
#' @param conditions Character. The conditions a trial must meet, any of
#'   \code{"trigger"}, \code{"declared"}, \code{"region_Ghat_lt1"}. Default:
#'   all three.
#' @param which \code{"first"} (default): the first trial by index meeting
#'   the conditions; \code{"all"}: every such trial.
#' @param itt_ub Numeric. The trigger's bound on the all-patient upper 95\%
#'   limit (strictly below); \code{Inf} drops the bound. Default 1.0.
#' @param region_itt Numeric. The trigger's bound on the region's hazard
#'   ratio (strictly above); \code{-Inf} drops the bound. Default 1.0.
#'
#' @return A data.frame with one row per selected trial (zero rows, with a
#'   message, when none qualifies): \code{sim}; the flags \code{trigger},
#'   \code{declared}, \code{region_Ghat_lt1} (logical); and the trial's hazard
#'   ratios \code{hr_itt_all}, \code{itt_ub} (its upper 95\% limit),
#'   \code{hr_itt_region}, \code{hr_region_Ghat} and, when the record has the
#'   \code{region_Gstar} row, \code{hr_region_Gstar}; and \code{sg_found}.
#'   The attribute \code{"n_meeting"} is the number of trials meeting the
#'   conditions.
#'
#' @seealso \code{\link{mrct_case_study}} to redraw and examine the selected
#'   trial; \code{\link{mrct_region_sims}}.
#'
#' @examples
#' \dontrun{
#' res <- mrct_region_sims(dgm, n_sims = 20, n_sample = 1000,
#'                         region_var = "z_region", seed = 613800000,
#'                         record = "full")
#' mrct_select_case(res)                    # the first qualifying trial
#' mrct_select_case(res, which = "all")     # every qualifying trial
#' }
#' @export
mrct_select_case <- function(res,
                             conditions = c("trigger", "declared", "region_Ghat_lt1"),
                             which = c("first", "all"),
                             itt_ub = 1.0,
                             region_itt = 1.0) {
  which <- match.arg(which)
  allowed <- c("trigger", "declared", "region_Ghat_lt1")
  if (!is.character(conditions) || !length(conditions) || !all(conditions %in% allowed))
    stop("conditions must be one or more of ", paste(sprintf('"%s"', allowed), collapse = ", "),
         ".", call. = FALSE)
  est <- attr(res, "estimators")
  if (is.null(est))
    stop("res has no \"estimators\" attribute: mrct_select_case() needs a ",
         "mrct_region_sims(record = \"full\") result (record = \"point\" has no per-set ",
         "estimates), and a data.table subset of it drops the attribute.", call. = FALSE)
  if (!all(c("sim", "any_found") %in% names(res)))
    stop("res must have the columns sim and any_found.", call. = FALSE)
  for (b in list(itt_ub, region_itt))
    if (!is.numeric(b) || length(b) != 1L || is.na(b))
      stop("itt_ub and region_itt must be single numbers (Inf / -Inf drop the bound).",
           call. = FALSE)

  est <- as.data.frame(est)
  sims <- sort(unique(res$sim))
  row_of <- function(k, col) {
    y <- est[est$estimator == k, , drop = FALSE]
    y[[col]][match(sims, y$sim)]
  }
  out <- data.frame(sim = sims)
  hr_itt <- row_of("itt_all", "hr"); ub <- row_of("itt_all", "hi")
  hr_reg <- row_of("itt_region", "hr"); hr_gh <- row_of("region_Ghat", "hr")
  af <- res$any_found[match(sims, res$sim)]
  lt <- function(v, b) !is.na(v) & v < b
  gt <- function(v, b) !is.na(v) & v > b
  out$trigger <- lt(ub, itt_ub) & gt(hr_reg, region_itt)
  out$declared <- !is.na(af) & af == 1
  out$region_Ghat_lt1 <- lt(hr_gh, 1)
  out$hr_itt_all <- hr_itt
  out$itt_ub <- ub
  out$hr_itt_region <- hr_reg
  out$hr_region_Ghat <- hr_gh
  if (any(est$estimator == "region_Gstar")) out$hr_region_Gstar <- row_of("region_Gstar", "hr")
  if ("sg_found" %in% names(res)) out$sg_found <- res$sg_found[match(sims, res$sim)]

  meet <- Reduce(`&`, lapply(conditions, function(k) out[[k]]))
  n_meeting <- sum(meet)
  sel <- which(meet)
  if (which == "first") sel <- utils::head(sel, 1L)
  if (!length(sel))
    message("mrct_select_case(): no trial of the ", length(sims), " meets ",
            paste(conditions, collapse = " & "), ".")
  out <- out[sel, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n_meeting") <- n_meeting
  attr(out, "conditions") <- conditions
  out
}


#' Redraw and examine one trial of an MRCT simulation
#'
#' Redraws trial \code{k} of a \code{\link{mrct_region_sims}} run exactly as
#' the simulation drew it, runs \code{\link{forestsearch}} on its non-region
#' patients (the region's patients as the test set) with the run's arguments
#' and multiplier resampling, and analyses four patient sets: all patients,
#' the region, region(\eqn{G^*}) (the region without the true harmed subgroup,
#' when the trial has \code{flag_harm}) and region(\eqn{\hat G}) (the region
#' without the declared subgroup \eqn{\hat H}).
#'
#' @section Reproducing a simulated trial:
#' \code{mrct_region_sims()} draws trial \code{k} in a worker with
#' \code{simulate_from_dgm(seed = seed + k)}, under the L'Ecuyer-CMRG
#' generator that the workers run; \code{mrct_case_study()} switches to that
#' generator for the draw (the same \code{simulate_from_dgm()} arguments: the
#' \code{sim_args} defaults overridden by \code{n}, \code{analysis_time},
#' \code{cens_adjust} and \code{seed + k}, and \code{entry_var =
#' "entrytime"} when \code{dgm$df_super} has that column) and restores the
#' caller's generator, and the caller's random-number state, on exit. The
#' identifier seeds its own random steps (data preparation and the
#' consistency splits), so it returns the simulation's \eqn{\hat H} under any
#' generator. Multiplier resampling runs after the draw under the caller's
#' generator with \code{seed = seed + k} unless \code{mr_inference_args}
#' gives a seed. (Inside \code{mrct_region_sims()} the resampling, if
#' requested through \code{fs_args}, runs under the workers' generator, so
#' its MR columns need not equal this object's.)
#'
#' @param dgm The DGM of the run (\code{\link{generate_aft_dgm_flex}},
#'   \code{\link{inject_mrct_structure}}).
#' @param seed The run's base seed (\code{mrct_region_sims(seed = )}).
#' @param k The trial's index (\code{sim}).
#' @param n The trial size (\code{mrct_region_sims(n_sample = )}); \code{NULL}
#'   (default) uses \code{nrow(dgm$df_super)}, as \code{mrct_region_sims()}
#'   does.
#' @param region_var The region indicator column. Default \code{"z_regA"}.
#' @param fs_args Named list of \code{forestsearch()} arguments. It may also
#'   hold the arguments that \code{mrct_region_sims()} takes directly
#'   (\code{sg_focus}, \code{maxk}, \code{hr.threshold},
#'   \code{hr.consistency}, \code{pconsistency.threshold},
#'   \code{confounders.name}, \code{conf_force}); those not given take
#'   \code{mrct_region_sims()}'s defaults, so a run's direct arguments and
#'   its \code{fs_args} together reproduce its identifier. The data, the
#'   outcome, treatment, event, id and potential-outcome names, the display
#'   flags and the (sequential) parallel plan are set here, as in
#'   \code{mrct_region_sims()}.
#' @param sim_args Named list of \code{simulate_from_dgm()} arguments, as in
#'   \code{mrct_region_sims()}. Default \code{list(rand_ratio = 1,
#'   draw_treatment = TRUE)}.
#' @param analysis_time Calendar time of the analysis. Default 60.
#' @param cens_adjust Censoring adjustment on the log scale. Default 0.
#' @param mr_inference Logical. Run multiplier resampling on \eqn{\hat H}
#'   (\code{\link{fs_mr_inference}}). Default \code{TRUE}.
#' @param mr_inference_args List passed to \code{forestsearch(mr_inference_args
#'   = )}; \code{seed} defaults to \code{seed + k}. Default: the field method
#'   with the complement de-biased and given its own field.
#' @param entry_delay \code{NULL} (default) or \code{list(var = <column>,
#'   delay = a)}, passed to \code{simulate_from_dgm()} (overrides
#'   \code{sim_args$entry_delay}).
#'
#' @return An object of class \code{mrct_case_study}, a list with
#' \describe{
#'   \item{\code{trial}}{The redrawn trial (data.frame).}
#'   \item{\code{fs}}{The \code{forestsearch()} result (its
#'     \code{mr_inference} element holds the MR result);
#'     \code{\link{mr_estimates_table}(x$fs)} gives the inference table.}
#'   \item{\code{sets}}{The four-set table: \code{set}, \code{label},
#'     \code{n}, \code{n_exp}, \code{n_ctl}, \code{events_exp},
#'     \code{events_ctl}, and the unadjusted two-sample Cox fit
#'     \code{Surv(y_sim, event_sim) ~ treat_sim} on the set: \code{hr},
#'     \code{lo}, \code{hi} (Wald 95\%), \code{se_loghr}. The rows are
#'     \code{itt_all}, \code{itt_region}, \code{region_Gstar} (when
#'     \code{flag_harm} exists) and \code{region_Ghat} (the whole region when
#'     nothing is declared).}
#'   \item{\code{in_Hhat}}{Logical, per trial row: membership of \eqn{\hat H}
#'     by \code{\link{get_dfpred}} on the selected definition.}
#'   \item{\code{confusion}}{When \code{flag_harm} exists, a list of two
#'     integer vectors \code{c(TP, FP, FN, TN)} of \eqn{\hat H} against
#'     \eqn{H^*} (\code{flag_harm == 1}), \code{nonregion} and \code{region};
#'     else \code{NULL}.}
#'   \item{\code{Hhat}}{The definition as \code{mrct_region_sims()} records
#'     it (\code{"none"} when nothing is declared).}
#'   \item{\code{mr}}{\code{fs$mr_inference}.}
#'   \item{\code{trigger}}{Logical: all-patient upper limit < 1.0 and region
#'     hazard ratio > 1.0.}
#'   \item{\code{call_args}}{\code{seed}, \code{k}, \code{n},
#'     \code{region_var}, \code{analysis_time}, the \code{simulate_from_dgm()}
#'     and \code{forestsearch()} argument lists (without data).}
#'   \item{\code{warnings}}{Warnings raised by \code{forestsearch()}, muffled
#'     and kept here.}
#' }
#'
#' @seealso \code{\link{mrct_select_case}}, \code{\link{mrct_region_sims}},
#'   \code{\link{mr_estimates_table}}.
#'
#' @examples
#' \dontrun{
#' k <- mrct_select_case(res)$sim
#' cs <- mrct_case_study(dgm, seed = 613800000, k = k, n = 1000,
#'                       region_var = "z_region",
#'                       fs_args = list(sg_focus = "hrMaxSG", maxk = 2,
#'                                      conf_force = character(0)),
#'                       sim_args = list(rand_ratio = 1, draw_treatment = TRUE,
#'                                       max_entry = 11),
#'                       analysis_time = 40, cens_adjust = 50)
#' cs
#' plot(cs)
#' mr_estimates_table(cs$fs)
#' }
#' @importFrom survival coxph Surv
#' @export
mrct_case_study <- function(dgm, seed, k, n = NULL, region_var = "z_regA",
                            fs_args = list(),
                            sim_args = list(rand_ratio = 1, draw_treatment = TRUE),
                            analysis_time = 60, cens_adjust = 0,
                            mr_inference = TRUE,
                            mr_inference_args = list(ci_method = "field",
                                                     include_complement = TRUE,
                                                     field_complement = TRUE),
                            entry_delay = NULL) {
  .fs_preserve_rng()
  if (!inherits(dgm, c("aft_dgm_flex", "aft_dgm")))
    stop("dgm must be an object created by generate_aft_dgm_flex()", call. = FALSE)
  if (missing(seed) || !is.numeric(seed) || length(seed) != 1L || !is.finite(seed))
    stop("seed must be the run's base seed, one finite number.", call. = FALSE)
  if (missing(k) || !is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1 || k != round(k))
    stop("k must be one positive integer, the trial's index.", call. = FALSE)
  if (!is.character(region_var) || length(region_var) != 1L || !region_var %in% names(dgm$df_super))
    stop("region_var '", region_var, "' not found in dgm$df_super", call. = FALSE)
  if (!is.list(fs_args) || !is.list(sim_args) || !is.list(mr_inference_args))
    stop("fs_args, sim_args and mr_inference_args must be lists.", call. = FALSE)
  if (is.null(n)) n <- nrow(dgm$df_super)

  # The caller's generator and stream are restored on exit.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)

  # ---- the draw: mrct_region_sims() section 4, under the workers' generator
  if (!is.null(entry_delay)) sim_args$entry_delay <- entry_delay
  sim_call_args <- utils::modifyList(sim_args, list(
    dgm = dgm, n = n,
    entry_var = if ("entrytime" %in% names(dgm$df_super)) "entrytime" else NULL,
    analysis_time = analysis_time, cens_adjust = cens_adjust,
    seed = seed + k))
  RNGkind("L'Ecuyer-CMRG")
  dfs <- data.table::as.data.table(do.call(simulate_from_dgm, sim_call_args))
  do.call(RNGkind, as.list(old_kind))

  # ---- the identifier: mrct_region_sims()'s construction of the call
  defaults <- list(
    sg_focus = "minSG", maxk = 1, hr.threshold = 0.90, hr.consistency = 0.80,
    pconsistency.threshold = 0.90,
    confounders.name = setdiff(
      grep("^z_|^ecog|^strat", names(dgm$df_super), value = TRUE),
      c("id", "treat", "tte", "event", "y_sim", "treat_sim", "event_sim", "t_true",
        "c_time", "lin_pred_0", "lin_pred_1", "loghr_po", "entrytime", "flag_harm")),
    conf_force = c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5"))
  fs_spec <- utils::modifyList(defaults, fs_args)
  if (!is.null(mr_inference_args$seed)) {
    mr_args <- mr_inference_args
  } else {
    mr_args <- c(mr_inference_args, list(seed = seed + k))
  }

  # ---- the analysis: the core shared with mrct_trial_analysis()
  core <- .mrct_analyze_trial(
    dfs, region_var,
    truth_var = if ("flag_harm" %in% names(dfs)) "flag_harm" else NULL,
    fs_spec, mr_inference, mr_args,
    nm = list(outcome = "y_sim", treat = "treat_sim", event = "event_sim",
              id = if ("id" %in% names(dfs)) "id" else NULL,
              po = if ("loghr_po" %in% names(dfs)) "loghr_po" else NULL),
    caller = "mrct_case_study")
  if (is.null(core$confusion))
    message("mrct_case_study(): the trial has no flag_harm column; the region(G*) set ",
            "and the 2 x 2 tables of H-hat against H* are omitted.")

  fs_call_args <- core$fs_call_args
  structure(list(
    trial = core$trial, fs = core$fs, sets = core$sets, in_Hhat = core$in_Hhat,
    confusion = core$confusion, Hhat = core$Hhat, mr = core$mr, trigger = core$trigger,
    call_args = list(seed = seed, k = k, n = n, region_var = region_var,
                     analysis_time = analysis_time,
                     sim_args = sim_call_args[setdiff(names(sim_call_args), "dgm")],
                     fs_args = fs_call_args[setdiff(names(fs_call_args), c("df.analysis", "df.test"))]),
    warnings = core$warnings),
    class = "mrct_case_study")
}


#' Print an MRCT case study
#'
#' The four-set table and, when multiplier resampling ran on a declared
#' subgroup, the MR estimates and bounds for \eqn{\hat H} and its non-region
#' complement \eqn{\hat H^c}.
#'
#' @param x An \code{mrct_case_study} object.
#' @param digits Digits for hazard ratios. Default 3.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.mrct_case_study <- function(x, digits = 3, ...) {
  fm <- .mrct_fm(digits)
  a <- x$call_args
  cat(sprintf("MRCT case study: trial %d (seed %s + %d), n = %d, region '%s'\n",
              as.integer(a$k), format(a$seed, scientific = FALSE), as.integer(a$k),
              nrow(x$trial), a$region_var))
  cat(sprintf("H-hat: %s\n", x$Hhat))
  cat(sprintf("Trigger (ITT ub < 1.0 and region ITT > 1.0): %s\n\n", x$trigger))
  .mrct_print_sets(x$sets, fm)
  .mrct_print_confusion(x$confusion)
  .mrct_print_mr(x$mr, fm)
  invisible(x)
}


#' Plot an MRCT case study
#'
#' The forest plot of the four sets' hazard ratios with Wald 95\% intervals
#' (log scale) and the Kaplan-Meier curves by arm in each set (one panel per
#' set, 2 x 2), the latter with \pkg{weightedsurv}.
#'
#' @param x An \code{mrct_case_study} object.
#' @param which Any of \code{"forest"}, \code{"km"}. Default both.
#' @param xlim Hazard-ratio range of the forest plot. Default
#'   \code{c(0.3, 2)}; intervals are clipped to it.
#' @param time_unit Character. The unit of the Kaplan-Meier time axis, used
#'   as its label (e.g. \code{"months"}); \code{NULL} (default) labels it
#'   \code{"time"}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @importFrom graphics par plot segments abline axis mtext text title
#' @export
plot.mrct_case_study <- function(x, which = c("forest", "km"), xlim = c(0.3, 2),
                                 time_unit = NULL, ...) {
  which <- match.arg(which, several.ok = TRUE)
  .mrct_plot_sets(x$sets, x$trial, x$call_args$region_var, "flag_harm", x$in_Hhat, which, xlim,
                  time_unit, x$call_args$analysis_time,
                  list(outcome = "y_sim", treat = "treat_sim", event = "event_sim"),
                  "plot.mrct_case_study")
  invisible(x)
}
