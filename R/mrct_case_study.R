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

  df_non <- dfs[dfs[[region_var]] == 0]
  df_reg <- dfs[dfs[[region_var]] == 1]

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
  fs_base <- list(
    df.analysis = as.data.frame(df_non), df.test = as.data.frame(df_reg),
    outcome.name = "y_sim", treat.name = "treat_sim", event.name = "event_sim",
    id.name = if ("id" %in% names(df_non)) "id" else NULL,
    potentialOutcome.name = if ("loghr_po" %in% names(df_non)) "loghr_po" else NULL,
    show_candidate_summary = FALSE, details = FALSE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential"), quiet = TRUE,
    mr_inference = isTRUE(mr_inference),
    mr_inference_args = if (isTRUE(mr_inference)) mr_args else list())
  fs_call_args <- utils::modifyList(fs_spec, fs_base)
  fs_warn <- character(0)
  fs <- withCallingHandlers(
    do.call(forestsearch, fs_call_args),
    warning = function(w) {
      if (!grepl("parallel_args missing required elements", conditionMessage(w)))
        fs_warn <<- c(fs_warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (length(fs_warn))
    message("mrct_case_study(): forestsearch() raised ", length(fs_warn),
            " warning(s), kept in $warnings.")

  # ---- the sets
  trial <- as.data.frame(dfs)
  declared <- !is.null(fs$sg.harm)
  in_Hhat <- if (declared) {
    get_dfpred(trial, sg.harm = fs$sg.harm, version = 2)$treat.recommend == 0
  } else {
    rep(FALSE, nrow(trial))
  }
  reg <- trial[[region_var]] == 1
  has_band <- "flag_harm" %in% names(trial)
  band <- if (has_band) trial$flag_harm == 1 else NULL
  set_idx <- list(itt_all = rep(TRUE, nrow(trial)), itt_region = reg)
  if (has_band) set_idx$region_Gstar <- reg & !band
  set_idx$region_Ghat <- reg & !in_Hhat
  labels <- c(itt_all = "All patients, ITT", itt_region = "Region, ITT",
              region_Gstar = "Region minus H* (G*)", region_Ghat = "Region minus H-hat (G-hat)")
  dt <- data.table::as.data.table(trial)
  sets <- do.call(rbind, lapply(names(set_idx), function(s) {
    f <- as.data.frame(.mrct_record_fit(dt, set_idx[[s]]))
    data.frame(set = s, label = unname(labels[s]), n = f$n, n_exp = f$n_exp,
               n_ctl = f$n - f$n_exp, events_exp = f$events_exp, events_ctl = f$events_ctl,
               hr = f$hr, lo = f$lo, hi = f$hi, se_loghr = f$se_loghr,
               stringsAsFactors = FALSE)
  }))
  rownames(sets) <- NULL

  confusion <- NULL
  if (has_band) {
    cc <- function(w) c(TP = sum(in_Hhat[w] & band[w]), FP = sum(in_Hhat[w] & !band[w]),
                        FN = sum(!in_Hhat[w] & band[w]), TN = sum(!in_Hhat[w] & !band[w]))
    confusion <- list(nonregion = cc(!reg), region = cc(reg))
  } else {
    message("mrct_case_study(): the trial has no flag_harm column; the region(G*) set ",
            "and the 2 x 2 tables of H-hat against H* are omitted.")
  }

  itt <- sets[sets$set == "itt_all", ]; rgn <- sets[sets$set == "itt_region", ]
  trigger <- !is.na(itt$hi) && itt$hi < 1 && !is.na(rgn$hr) && rgn$hr > 1

  structure(list(
    trial = trial, fs = fs, sets = sets, in_Hhat = in_Hhat, confusion = confusion,
    Hhat = if (declared) paste(fs$sg.harm, collapse = " & ") else "none",
    mr = fs$mr_inference, trigger = trigger,
    call_args = list(seed = seed, k = k, n = n, region_var = region_var,
                     analysis_time = analysis_time,
                     sim_args = sim_call_args[setdiff(names(sim_call_args), "dgm")],
                     fs_args = fs_call_args[setdiff(names(fs_call_args), c("df.analysis", "df.test"))]),
    warnings = fs_warn),
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
  fm <- function(v) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = digits))
  a <- x$call_args
  cat(sprintf("MRCT case study: trial %d (seed %s + %d), n = %d, region '%s'\n",
              as.integer(a$k), format(a$seed, scientific = FALSE), as.integer(a$k),
              nrow(x$trial), a$region_var))
  cat(sprintf("H-hat: %s\n", x$Hhat))
  cat(sprintf("Trigger (ITT ub < 1.0 and region ITT > 1.0): %s\n\n", x$trigger))
  s <- x$sets
  tab <- data.frame(Set = s$label, n = s$n,
                    `exp n/events` = sprintf("%d/%d", s$n_exp, s$events_exp),
                    `ctl n/events` = sprintf("%d/%d", s$n_ctl, s$events_ctl),
                    HR = fm(s$hr), `95% CI` = sprintf("%s, %s", fm(s$lo), fm(s$hi)),
                    check.names = FALSE)
  print(tab, row.names = FALSE, right = FALSE)
  if (!is.null(x$confusion)) {
    cat("\nH-hat against H* (TP, FP, FN, TN):\n")
    for (p in names(x$confusion))
      cat(sprintf("  %-10s %s\n", p, paste(x$confusion[[p]], collapse = ", ")))
  }
  mr <- x$mr
  if (is.null(mr) || is.null(mr$debiased)) {
    cat("\nMultiplier resampling: not run, or nothing declared.\n")
    return(invisible(x))
  }
  cat(sprintf("\nMultiplier resampling (%s, %s draws):\n", mr$ci_method,
              if (!is.null(mr$settings$draws)) mr$settings$draws else "?"))
  g <- function(...) { v <- mr; for (k in c(...)) { if (!is.list(v) || is.null(v[[k]])) return(NA_real_); v <- v[[k]] }; v }
  rows <- data.frame(
    set = c("H-hat", "H-hat complement (non-region)"),
    naive = fm(c(g("naive", "est"), g("complement", "naive", "est"))),
    debiased = fm(c(g("debiased", "est"), g("complement", "debiased", "est"))),
    `IJ 95% CI` = sprintf("%s, %s", fm(c(g("debiased", "lower"), g("complement", "debiased", "lower"))),
                          fm(c(g("debiased", "upper"), g("complement", "debiased", "upper")))),
    `field 1-sided 95%` = c(paste(">=", fm(g("field", "lower_1s"))),
                            paste("<=", fm(g("field", "complement", "upper_1s")))),
    Bonferroni = c(paste(">=", fm(g("field", "joint", "bonf_lower_H"))),
                   paste("<=", fm(g("field", "joint", "bonf_upper_Hc")))),
    check.names = FALSE)
  print(rows, row.names = FALSE, right = FALSE)
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
  if (is.null(time_unit)) time_unit <- "time"
  s <- x$sets
  if ("forest" %in% which) {
    ord <- intersect(c("itt_all", "itt_region", "region_Ghat", "region_Gstar"), s$set)
    lab <- c(itt_all = "All patients, ITT", itt_region = "Region, ITT",
             region_Ghat = "Region minus the identified subgroup",
             region_Gstar = "Region minus the true harmed subgroup")
    cols <- c(itt_all = "black", itt_region = "firebrick", region_Ghat = "steelblue",
              region_Gstar = "grey40")[ord]
    f <- s[match(ord, s$set), ]; y <- rev(seq_along(ord))
    op <- graphics::par(mar = c(4, 17, 1, 9)); on.exit(graphics::par(op), add = TRUE)
    graphics::plot(pmin(pmax(f$hr, xlim[1]), xlim[2]), y, log = "x", xlim = xlim,
                   ylim = c(0.5, length(ord) + 0.5), pch = 15, cex = 1.3, yaxt = "n", xaxt = "n",
                   xlab = "hazard ratio, experimental vs control (log scale, 95% CI)", ylab = "",
                   col = cols)
    graphics::segments(pmax(f$lo, xlim[1]), y, pmin(f$hi, xlim[2]), y, lwd = 2, col = cols)
    graphics::abline(v = 1, lty = 2)
    at <- c(0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5)
    at <- at[at >= xlim[1] & at <= xlim[2]]
    graphics::axis(1, at = at, labels = format(at))
    graphics::axis(2, at = y, labels = lab[ord], las = 1, tick = FALSE)
    graphics::mtext(sprintf("%.3f (%.3f, %.3f)", f$hr, f$lo, f$hi), side = 4, at = y,
                    las = 1, line = 0.5)
    graphics::par(op)
  }
  if ("km" %in% which) {
    if (!requireNamespace("weightedsurv", quietly = TRUE)) {
      message("plot.mrct_case_study(): weightedsurv is not installed; the Kaplan-Meier ",
              "figure is skipped.")
      return(invisible(x))
    }
    ta <- x$call_args$analysis_time
    br <- ta / 4
    op <- graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.2, 3, 1))
    on.exit(graphics::par(op), add = TRUE)
    for (i in seq_len(nrow(s))) {
      d <- x$trial[.mrct_case_set(x, s$set[i]), ]
      dfc <- suppressWarnings(weightedsurv::df_counting(
        df = d, by.risk = br, tte.name = "y_sim", event.name = "event_sim",
        treat.name = "treat_sim", arms = c("experimental", "control")))
      weightedsurv::plot_weighted_km(
        dfc, conf.int = FALSE, show.cox = FALSE, show.logrank = FALSE, show.med = FALSE,
        show_arm_legend = TRUE, arms = c("experimental", "control"), put.legend.arms = "left",
        by.risk = br, xmin = 0, xmax = ta, ymin = 0, ymax = 1.12, Xlab = time_unit,
        Ylab = "proportion event-free", col.0 = "black", col.1 = "blue", lwds = c(1.5, 1.5),
        censor.mark.all = FALSE, risk.cex = 0.8, arm.cex = 0.85)
      graphics::axis(1, at = seq(0, ta, br))
      graphics::text(ta, 1.12, adj = c(1, 1), cex = 0.85, labels = paste(
        sprintf("HR %.3f (95%% CI %.3f, %.3f)", s$hr[i], s$lo[i], s$hi[i]),
        sprintf("experimental n = %d, events %d", s$n_exp[i], s$events_exp[i]),
        sprintf("control n = %d, events %d", s$n_ctl[i], s$events_ctl[i]), sep = "\n"))
      graphics::title(sprintf("%d. %s", i, s$label[i]), cex.main = 1.05)
    }
  }
  invisible(x)
}

#' Rows of the trial in one of the case study's sets
#' @noRd
.mrct_case_set <- function(x, set) {
  reg <- x$trial[[x$call_args$region_var]] == 1
  switch(set,
         itt_all = rep(TRUE, nrow(x$trial)),
         itt_region = reg,
         region_Gstar = reg & x$trial$flag_harm != 1,
         region_Ghat = reg & !x$in_Hhat)
}
